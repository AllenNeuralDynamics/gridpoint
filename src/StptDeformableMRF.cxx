/*=========================================================================
 *
 *  StptDeformableMRF.cxx
 *
 *  Discrete symmetric deformable registration, derived from
 *  AllenInstitute/stpt_registration (see StptDeformableMRF.h for the exact upstream
 *  files). Copyright (c) Allen Institute for Brain Science.
 *  Licensed under the Allen Institute Software License; see LICENSE.
 *
 *  The numeric routines below are the upstream implementations, lifted out of the
 *  LIMS-bound dispatcher class into free functions. The algorithm is unchanged; what
 *  changed is:
 *    - class members (image size, histogram bins, sample count, intensity maxima,
 *      metric outputs) became an explicit MrfContext passed by reference;
 *    - `delete x;` on array allocations became `delete[] x;` (the originals were
 *      undefined behaviour);
 *    - the orchestrator no longer reads an atlas or a model directory.
 *
 *=========================================================================*/
#include "StptDeformableMRF.h"

#include "itkMultiThreaderBase.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

namespace stpt
{
namespace
{

using std::max;
using std::min;
using std::pop_heap;
using std::push_heap;

/** Index type for anything that can address a full-resolution volume, or a cost
 *  volume (control-grid nodes x label space).
 *
 *  Upstream used plain ``int`` throughout. That caps a volume at 2^31-1 voxels, but
 *  the real ceiling is lower and less obvious: the cost volumes are indexed as
 *  ``node + label * nodeCount``, which for a CCF 10 um volume at the finest level is
 *  about 2.35e9 -- over the limit while the volume itself (1.20e9 voxels) is still
 *  under it. Overflow there wraps to a negative index and corrupts silently rather
 *  than failing, so every expression that can reach either scale is 64-bit. Grid
 *  vertex ids, histogram bin numbers and label counts stay ``int``: they are bounded
 *  by the control grid and the search neighbourhood, which are orders of magnitude
 *  smaller. */
using vidx = std::int64_t;

/** Clamp an index into [0, hi]. Replaces nested min/max whose literal 0 would not
 *  deduce against a 64-bit index. */
inline vidx
clampIdx(vidx v, vidx hi)
{
  return v < 0 ? 0 : (v > hi ? hi : v);
}

/** What the upstream routines used to read off their enclosing class. */
struct MrfContext
{
  vidx         sizeX{ 0 };
  vidx         sizeY{ 0 };
  vidx         sizeZ{ 0 };
  unsigned int numBin{ 32 };
  int          sampleNum{ 1024 };
  float        maxFixed{ 255.0f };
  float        maxMoving{ 255.0f };
  float        metricBefore{ 0.0f };
  float        metricAfter{ 0.0f };
};

/** Edge in the control-grid graph; ordered so std::push_heap yields a MIN-heap.
 *  Upstream: Common/idpRegistrationUtilities.{h,cxx}. */
class Relation
{
public:
  double weight;
  int    vert1;
  int    vert2;
  Relation(double w = 0, int v1 = 0, int v2 = 0)
    : weight(w)
    , vert1(v1)
    , vert2(v2)
  {}
  bool
  operator<(const Relation & b) const
  {
    return (this->weight > b.weight);
  }
};

/** Accept an edge only if it joins the tree to a vertex not yet in it. */
int
newRelation(Relation edge1, Relation & edgeout, bool * vertices)
{
  bool new1 = vertices[edge1.vert1];
  bool new2 = vertices[edge1.vert2];
  int  out1;
  if (new1 ^ new2)
  {
    if (new1)
    {
      out1 = edge1.vert2;
      edgeout = Relation(edge1.weight, edge1.vert1, edge1.vert2);
    }
    else
    {
      out1 = edge1.vert1;
      edgeout = Relation(edge1.weight, edge1.vert2, edge1.vert1);
    }
  }
  else
  {
    out1 = -1;
  }
  return out1;
}

// ---------------------------------------------------------------------------
// Parallel execution.
//
// The deformable stage was entirely single-threaded, which made it 93% of wall
// time on a 24-core box and left larger machines almost completely idle. The
// loops below are data-parallel: each iteration writes its own disjoint slice
// of the output, so splitting the index range into contiguous chunks preserves
// the result exactly, including its floating-point rounding -- no iteration
// contributes to a shared accumulator.
//
// ITK thread pooling is used rather than OpenMP or raw std::thread: ITK is
// already the only dependency, the pool avoids paying thread-creation cost on
// the many short parallel regions inside one pyramid level, and it honours the
// existing --threads flag with no extra plumbing.
// ---------------------------------------------------------------------------
template <typename Body>
void
parallelFor(vidx count, Body body)
{
  if (count <= 0) { return; }
  const vidx maxThreads =
    std::max<vidx>(1, static_cast<vidx>(itk::MultiThreaderBase::GetGlobalDefaultNumberOfThreads()));
  // More chunks than threads: per-chunk cost is uneven (boundary nodes take a
  // slower path), and finer chunks let the pool even that out.
  const vidx chunks = std::min<vidx>(count, maxThreads * 4);
  if (chunks <= 1)
  {
    body(static_cast<vidx>(0), count);
    return;
  }
  itk::MultiThreaderBase::Pointer mt = itk::MultiThreaderBase::New();
  mt->ParallelizeArray(
    0,
    static_cast<itk::SizeValueType>(chunks),
    [&](itk::SizeValueType c) {
      const vidx begin = static_cast<vidx>(c) * count / chunks;
      const vidx end = (static_cast<vidx>(c) + 1) * count / chunks;
      body(begin, end);
    },
    nullptr);
}

// ---------------------------------------------------------------------------
// Upstream numeric core, verbatim apart from the changes noted in the file header.
// ---------------------------------------------------------------------------
        void
        minimumInd(float *numbers, float &value, int &index, vidx length) {
            value = numbers[0];
            vidx i;
            index = 0;
            for (i = 1; i < length; i++) {
                if (numbers[i] < value) {
                    index = static_cast<int>(i);
                    value = numbers[i];
                }
            }
        }

        void
        DfmfldModulation(MrfContext & ctx, float* dfmfld, int* ordered, int* parents, int stepsize) {
            vidx xsize = ctx.sizeX;
            vidx ysize = ctx.sizeY;
            vidx zsize = ctx.sizeZ;
            int xgridsize = static_cast<int>(xsize / stepsize);
            int ygridsize = static_cast<int>(ysize / stepsize);
            int zgridsize = static_cast<int>(zsize / stepsize);

            int num_vertices = xgridsize * ygridsize*zgridsize;

            int num_neighbours = 6;
            float* edgecost = new float[num_vertices * num_neighbours];
            int* index_neighbours = new int[num_vertices * num_neighbours];
            for (int i = 0; i < num_vertices * num_neighbours; i++) {
                edgecost[i] = 0.0;
                index_neighbours[i] = -1;
            }

            int xdelta[6] = {-1, 1, 0, 0, 0, 0};
            int ydelta[6] = {0, 0, -1, 1, 0, 0};
            int zdelta[6] = {0, 0, 0, 0, -1, 1};

            //calculate weights 
            parallelFor(zgridsize, [&](vidx kBegin, vidx kEnd) {
            for (int k = static_cast<int>(kBegin); k < static_cast<int>(kEnd); k++) {
                for (int j = 0; j < ygridsize; j++) {
                    for (int i = 0; i < xgridsize; i++) {
                        for (int nb = 0; nb < num_neighbours; nb++) {
                            if ((i + ydelta[nb]) >= 0 & (i + ydelta[nb]) < xgridsize & (j + xdelta[nb]) >= 0 & (j + xdelta[nb]) < ygridsize & (k + zdelta[nb]) >= 0 & (k + zdelta[nb]) < zgridsize) {
                                index_neighbours[i + j * xgridsize + k * xgridsize * ygridsize + nb * num_vertices] = i + ydelta[nb]+(j + xdelta[nb]) * xgridsize + (k + zdelta[nb]) * xgridsize*ygridsize;
                                for (int k1 = 0; k1 < stepsize; k1++) {
                                    for (int j1 = 0; j1 < stepsize; j1++) {
                                        for (int i1 = 0; i1 < stepsize; i1++) {
                                            const vidx x = j * stepsize + j1;
                                            const vidx y = i * stepsize + i1;
                                            const vidx z = k * stepsize + k1;
                                            const vidx xnew = (j + xdelta[nb]) * stepsize + j1;
                                            const vidx ynew = (i + ydelta[nb]) * stepsize + i1;
                                            const vidx znew = (k + zdelta[nb]) * stepsize + k1;
                                            edgecost[i + j * xgridsize + k * xgridsize * ygridsize + nb * num_vertices] += fabs(dfmfld[y + x * xsize + z * xsize * ysize] - dfmfld[ynew + xnew * xsize + znew * xsize * ysize]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
            });

            int root = xgridsize / 2 + ygridsize / 2 * xgridsize + zgridsize / 2 * xgridsize*ygridsize;

            std::vector<Relation> sorted;
            bool* vertices = new bool[num_vertices];
            int* level = new int[num_vertices];
            for (int i = 0; i < num_vertices; i++) {
                vertices[i] = false;
                parents[i] = -1;
            }

            level[root] = 0;
            int last = root;
            vertices[root] = true;
            Relation edgeout = Relation(0.0, -1, -1);
            Relation minedge = Relation(0.0, -1, -1);
            float cost = 0.0;
            
            for (int i = 0; i < num_vertices - 1; i++) { // have all grids added
                for (int j = 0; j < num_neighbours; j++) {
                    int n = index_neighbours[last + j * num_vertices];
                    if (n >= 0) {
                        sorted.push_back(Relation(edgecost[last + j * num_vertices], last, n));
                        push_heap(sorted.begin(), sorted.end());
                    }
                }
                last = -1;
                //find lowest
                while (last == -1) {
                    minedge = sorted.front();
                    pop_heap(sorted.begin(), sorted.end());
                    sorted.pop_back();
                    last = newRelation(minedge, edgeout, vertices); //return next valid 
                }
                cost += edgeout.weight;
                vertices[last] = true;
                level[edgeout.vert2] = level[edgeout.vert1] + 1;
                parents[edgeout.vert2] = edgeout.vert1;
            }

            //re-caculate
            int maxlevel = 0;
            for (int i = 0; i < num_vertices; i++) {
                if (level[i] > maxlevel)
                    maxlevel = level[i];
            }
            maxlevel++;
            int* leveloffset = new int[maxlevel];
            int* levelcount = new int[maxlevel];
            for (int i = 0; i < maxlevel; i++) {
                leveloffset[i] = 0;
                levelcount[i] = 0;
            }
            for (int i = 0; i < num_vertices; i++) {
                if (level[i] < maxlevel - 1)
                    leveloffset[level[i] + 1]++; //counting number 
            }
            for (int i = 1; i < maxlevel; i++) {
                leveloffset[i] += leveloffset[i - 1];
            }
            for (int i = 0; i < num_vertices; i++) {
                int num = leveloffset[level[i]] + levelcount[level[i]];
                levelcount[level[i]]++;
                ordered[num] = i;
            }

            sorted.clear();

            delete[] edgecost;
            delete[] index_neighbours;
            delete[] levelcount;
            delete[] leveloffset;
            delete[] vertices;
            delete[] level;

        }

        void
        energy(float *val, int* index, int length, float offset, int k, int* v, float* z, float* f, int* indexinv) {
            float INF = 1e10;

            int j = 0;
            z[0] = -INF;
            z[1] = INF;
            v[0] = 0;
            for (int q = 1; q < length; q++) {
                float s = ((val[q * k] + powf((float) q + offset, 2.0))-(val[v[j] * k] + powf((float) v[j] + offset, 2.0))) / (2.0 * (float) (q - v[j]));
                while (s <= z[j]) {
                    j--;
                    s = ((val[q * k] + powf((float) q + offset, 2.0))-(val[v[j] * k] + powf((float) v[j] + offset, 2.0))) / (2.0 * (float) (q - v[j]));
                }
                j++;
                v[j] = q;
                z[j] = s;
                z[j + 1] = INF;
            }
            j = 0;
            for (int q = 0; q < length; q++) {
                f[q] = val[q * k];
                indexinv[q] = index[q * k];
            }
            for (int q = 0; q < length; q++) {
                while (z[j + 1] < q) {
                    j++;
                }
                index[q * k] = indexinv[v[j]];
                val[q * k] = powf((float) q - ((float) v[j] + offset), 2.0) + f[v[j]];
            }

        }

        void
        regulate(float* r, int* indr, int rl, float dx, float dy, float dz) {

            for (int i = 0; i < rl * rl * rl; i++) {
                indr[i] = i;
            }
            int* v = new int[rl];
            float* z = new float[rl + 1];
            float* f = new float[rl];
            int* i1 = new int[rl];

            for (int k = 0; k < rl; k++) {
                for (int i = 0; i < rl; i++) {
                    energy(r + i + k * rl*rl, indr + i + k * rl*rl, rl, -dx, rl, v, z, f, i1);
                }
            }
            for (int k = 0; k < rl; k++) {
                for (int j = 0; j < rl; j++) {
                    energy(r + j * rl + k * rl*rl, indr + j * rl + k * rl*rl, rl, -dy, 1, v, z, f, i1);
                }
            }

            for (int j = 0; j < rl; j++) {
                for (int i = 0; i < rl; i++) {
                    energy(r + i + j*rl, indr + i + j*rl, rl, -dz, rl*rl, v, z, f, i1);
                }
            }
            delete []i1;
            delete []f;

            delete []v;
            delete []z;


        }

        void
        interpolate3d(float* interp, float* input, float* x1, float* y1, float* z1, vidx xsize, vidx ysize, vidx zsize, vidx xsizenew, vidx ysizenew, vidx zsizenew, bool flag) {
            using namespace std;

            // Output rows are independent, and there are ysize*zsize of them -- ample
            // parallel width even on a coarse control grid.
            const vidx rowCount = ysize * zsize;
            parallelFor(rowCount, [&](vidx rowBegin, vidx rowEnd) {
                for (vidx r = rowBegin; r < rowEnd; r++) {
                    const vidx k = r / ysize;
                    const vidx j = r - k * ysize;
                    for (vidx i = 0; i < xsize; i++) {
                        vidx x = (vidx) floor(x1[i + j * xsize + k * xsize * ysize]);
                        vidx y = (vidx) floor(y1[i + j * xsize + k * xsize * ysize]);
                        vidx z = (vidx) floor(z1[i + j * xsize + k * xsize * ysize]);
                        float deltax = x1[i + j * xsize + k * xsize * ysize] - x;
                        float deltay = y1[i + j * xsize + k * xsize * ysize] - y;
                        float deltaz = z1[i + j * xsize + k * xsize * ysize] - z;

                        if (flag) {
                            x += j;
                            y += i;
                            z += k;
                        }
                        interp[i + j * xsize + k * xsize * ysize] = (1.0 - deltax)*(1.0 - deltay)*(1.0 - deltaz) * input[clampIdx(y, xsizenew - 1) + clampIdx(x, ysizenew - 1) * xsizenew + clampIdx(z, zsizenew - 1) * xsizenew * ysizenew]+
                                (1.0 - deltax) * deltay * (1.0 - deltaz) * input[clampIdx(y + 1, xsizenew - 1) + clampIdx(x, ysizenew - 1) * xsizenew + clampIdx(z, zsizenew - 1) * xsizenew * ysizenew] +
                                deltax * (1.0 - deltay)*(1.0 - deltaz) * input[clampIdx(y, xsizenew - 1) + clampIdx(x + 1, ysizenew - 1) * xsizenew + clampIdx(z, zsizenew - 1) * xsizenew * ysizenew]+
                                (1.0 - deltax)*(1.0 - deltay) * deltaz * input[clampIdx(y, xsizenew - 1) + clampIdx(x, ysizenew - 1) * xsizenew + clampIdx(z + 1, zsizenew - 1) * xsizenew * ysizenew] +
                                deltax * deltay * (1.0 - deltaz) * input[clampIdx(y + 1, xsizenew - 1) + clampIdx(x + 1, ysizenew - 1) * xsizenew + clampIdx(z, zsizenew - 1) * xsizenew * ysizenew]+
                                (1.0 - deltax) * deltay * deltaz * input[clampIdx(y + 1, xsizenew - 1) + clampIdx(x, ysizenew - 1) * xsizenew + clampIdx(z + 1, zsizenew - 1) * xsizenew * ysizenew] +
                                deltax * (1.0 - deltay) * deltaz * input[clampIdx(y, xsizenew - 1) + clampIdx(x + 1, ysizenew - 1) * xsizenew + clampIdx(z + 1, zsizenew - 1) * xsizenew * ysizenew] +
                                deltax * deltay * deltaz * input[clampIdx(y + 1, xsizenew - 1) + clampIdx(x + 1, ysizenew - 1) * xsizenew + clampIdx(z + 1, zsizenew - 1) * xsizenew * ysizenew];
                    }
                }
            });
        }

        void
        InterpDeformationFld(float* fldx_in, float* fldy_in, float* fldz_in, float* fldx, float* fldy, float* fldz, vidx xsize_in, vidx ysize_in, vidx zsize_in, vidx xsize, vidx ysize, vidx zsize) {


            float scale_m = (float) xsize_in / (float) xsize;
            float scale_n = (float) ysize_in / (float) ysize;
            float scale_o = (float) zsize_in / (float) zsize;

            float* x1 = new float[xsize_in * ysize_in * zsize_in];
            float* y1 = new float[xsize_in * ysize_in * zsize_in];
            float* z1 = new float[xsize_in * ysize_in * zsize_in];
            const vidx fillRows = ysize_in * zsize_in;
            parallelFor(fillRows, [&](vidx rowBegin, vidx rowEnd) {
                for (vidx r = rowBegin; r < rowEnd; r++) {
                    const vidx k = r / ysize_in;
                    const vidx j = r - k * ysize_in;
                    for (vidx i = 0; i < xsize_in; i++) {
                        x1[i + j * xsize_in + k * xsize_in * ysize_in] = j / scale_n;
                        y1[i + j * xsize_in + k * xsize_in * ysize_in] = i / scale_m;
                        z1[i + j * xsize_in + k * xsize_in * ysize_in] = k / scale_o;
                    }
                }
            });

            interpolate3d(fldx_in, fldx, x1, y1, z1, xsize_in, ysize_in, zsize_in, xsize, ysize, zsize, false);
            interpolate3d(fldy_in, fldy, x1, y1, z1, xsize_in, ysize_in, zsize_in, xsize, ysize, zsize, false);
            interpolate3d(fldz_in, fldz, x1, y1, z1, xsize_in, ysize_in, zsize_in, xsize, ysize, zsize, false);

            delete []x1;
            delete []y1;
            delete []z1;

        }

        void
        fldinv(float* fldx_in, float* fldy_in, float* fldz_in, float* fldx, float* fldy, float* fldz, vidx m, vidx n, vidx o) {
            float* un = new float[m * n * o];
            float* vn = new float[m * n * o];
            float* wn = new float[m * n * o];
            float* uin = new float[m * n * o];
            float* vin = new float[m * n * o];
            float* win = new float[m * n * o];

            for (vidx i = 0; i < m * n * o; i++) {
                uin[i] = -fldx[i];
                vin[i] = -fldy[i];
                win[i] = -fldz[i];
                fldx_in[i] = 0.0;
                fldy_in[i] = 0.0;
                fldz_in[i] = 0.0;
            }
            for (int it = 0; it < 10; it++) {
                interpolate3d(un, uin, fldx_in, fldy_in, fldz_in, m, n, o, m, n, o, true);
                interpolate3d(vn, vin, fldx_in, fldy_in, fldz_in, m, n, o, m, n, o, true);
                interpolate3d(wn, win, fldx_in, fldy_in, fldz_in, m, n, o, m, n, o, true);
                for (vidx i = 0; i < m * n * o; i++) {
                    fldx_in[i] = un[i];
                    fldy_in[i] = vn[i];
                    fldz_in[i] = wn[i];

                }
            }
            delete []uin;
            delete []vin;
            delete []win;
            delete []un;
            delete []vn;
            delete []wn;
        }

        void
        combineDeformation(float* fldx, float* fldy, float* fldz, float* fldx1, float* fldy1, float* fldz1, float* fldx2, float* fldy2, float* fldz2, vidx m, vidx n, vidx o) {
            float* uc = new float[m * n * o];
            float* vc = new float[m * n * o];
            float* wc = new float[m * n * o];

            interpolate3d(uc, fldx1, fldx2, fldy2, fldz2, m, n, o, m, n, o, true);
            interpolate3d(vc, fldy1, fldx2, fldy2, fldz2, m, n, o, m, n, o, true);
            interpolate3d(wc, fldz1, fldx2, fldy2, fldz2, m, n, o, m, n, o, true);

            for (vidx i = 0; i < m * n * o; i++) {
                fldx[i] = uc[i] + fldx2[i];
                fldy[i] = vc[i] + fldy2[i];
                fldz[i] = wc[i] + fldz2[i];

            }
            delete []uc;
            delete []vc;
            delete []wc;

        }

        void
        smoothDfm(float* fldx, float* fldy, float *fldz, vidx m, vidx n, vidx o, int expsteps, int factor) {
            float* fldx1 = new float[m * n * o];
            float* fldy1 = new float[m * n * o];
            float* fldz1 = new float[m * n * o];

            float* fldx2 = new float[m * n * o];
            float* fldy2 = new float[m * n * o];
            float* fldz2 = new float[m * n * o];

            float scale = 1.0 / (float) factor;
            float coeff = 1.0 / (float) powf(2.0, expsteps);
            for (vidx i = 0; i < m * n * o; i++) {
                fldx1[i] = coeff * fldx[i] * scale;
                fldy1[i] = coeff * fldy[i] * scale;
                fldz1[i] = coeff * fldz[i] * scale;

            }
            for (int it = 0; it < expsteps; it++) {
                combineDeformation(fldx2, fldy2, fldz2, fldx1, fldy1, fldz1, fldx1, fldy1, fldz1, m, n, o);
                for (vidx i = 0; i < m * n * o; i++) {
                    fldx1[i] = fldx2[i];
                    fldy1[i] = fldy2[i];
                    fldz1[i] = fldz2[i];
                }
            }
            for (vidx i = 0; i < m * n * o; i++) {
                fldx[i] = fldx2[i]*(float) factor;
                fldy[i] = fldy2[i]*(float) factor;
                fldz[i] = fldz2[i]*(float) factor;

            }
            delete []fldx2;
            delete []fldy2;
            delete []fldz2;
            delete []fldx1;
            delete []fldy1;
            delete []fldz1;

        }

        void
        symmetricMapping(float* xforward, float* yforward, float* zforward, float* xbackward, float* ybackward, float* zbackward, vidx m, vidx n, vidx o, int factor) {
            float factor1 = 1.0 / (float) factor;
            float* usym = new float[m* n * o];
            float* vsym = new float[m* n * o];
            float* wsym = new float[m* n * o];
            float* usym2 = new float[m* n * o];
            float* vsym2 = new float[m* n * o];
            float* wsym2 = new float[m* n * o];

            for (vidx i = 0; i < m * n * o; i++) {
                usym[i] = xforward[i]*(factor1 * 0.5);
                vsym[i] = yforward[i]*(factor1 * 0.5);
                wsym[i] = zforward[i]*(factor1 * 0.5);
                usym2[i] = xbackward[i]*(factor1 * 0.5);
                vsym2[i] = ybackward[i]*(factor1 * 0.5);
                wsym2[i] = zbackward[i]*(factor1 * 0.5);

            }
            float* ui = new float[m * n * o];
            float* vi = new float[m * n * o];
            float* wi = new float[m * n * o];
            float* u2i = new float[m * n * o];
            float* v2i = new float[m * n * o];
            float* w2i = new float[m * n * o];
            fldinv(ui, vi, wi, usym, vsym, wsym, m, n, o);
            fldinv(u2i, v2i, w2i, usym2, vsym2, wsym2, m, n, o);

            combineDeformation(xforward, yforward, zforward, u2i, v2i, w2i, usym, vsym, wsym, m, n, o);
            combineDeformation(xbackward, ybackward, zbackward, ui, vi, wi, usym2, vsym2, wsym2, m, n, o);

            for (vidx i = 0; i < m * n * o; i++) {
                xforward[i] *= (float) factor;
                yforward[i] *= (float) factor;
                zforward[i] *= (float) factor;
                xbackward[i] *= (float) factor;
                ybackward[i] *= (float) factor;
                zbackward[i] *= (float) factor;
            }

            delete []ui;
            delete []vi;
            delete []wi;
            delete []u2i;
            delete []v2i;
            delete []w2i;
            delete []usym;
            delete []vsym;
            delete []wsym;
            delete []usym2;
            delete []vsym2;
            delete []wsym2;

        }

        void
        genDfm(MrfContext & ctx, float* fldx, float* fldy, float* fldz, float* fldx0, float* fldy0, float* fldz0, float* costMatrix, float lambda, int hw, int stepsize, float gridmin, int* sorted, int* relations) {

            using namespace std;

            vidx sizex = ctx.sizeX;
            vidx sizey = ctx.sizeY;
            vidx sizez = ctx.sizeZ;

            int gridsizex = static_cast<int>(sizex / stepsize);
            int gridsizey = static_cast<int>(sizey / stepsize);
            int gridsizez = static_cast<int>(sizez / stepsize);


            // node count x label count is the cost-volume stride: 64-bit.
            vidx gridsizevol = static_cast<vidx>(gridsizex) * gridsizey * gridsizez;

            // deformation field
            int searchsize = hw * 2 + 1;
            float* xs = new float[searchsize * searchsize * searchsize];
            float* ys = new float[searchsize * searchsize * searchsize];
            float* zs = new float[searchsize * searchsize * searchsize];

            for (int i = 0; i < searchsize; i++) {
                for (int j = 0; j < searchsize; j++) {
                    for (int k = 0; k < searchsize; k++) {
                        xs[i + j * searchsize + k * searchsize * searchsize] = (j - hw) * gridmin;
                        ys[i + j * searchsize + k * searchsize * searchsize] = (i - hw) * gridmin;
                        zs[i + j * searchsize + k * searchsize * searchsize] = (k - hw) * gridmin;
                    }
                }
            }
            vidx searchvol = static_cast<vidx>(searchsize) * searchsize * searchsize;

            int *labels = new int[gridsizevol];
            short *indexall = new short[gridsizevol * searchvol];
            float *energy = new float[searchvol];
            float *values = new float[searchvol];
            int *index = new int[searchvol];


            for (vidx i = 0; i < searchvol; i++) {
                energy[i] = 0;
            }

            int frac = (int) (gridsizevol / 25);  // progress ticks only
            if (frac < 1) { frac = 1; }  // guard: i % frac below

            //calculate cost
            for (vidx i = (gridsizevol - 1); i > 0; i--) { //do for each control point
                if ((i % frac) == 0) {
                    cout << "x" << flush;
                }
                int next = sorted[i];
                int last = relations[sorted[i]];
                int znext = next / (gridsizex * gridsizey);
                int xnext = (next - znext * gridsizex * gridsizey) / gridsizex;
                int ynext = next - znext * gridsizex * gridsizey - xnext*gridsizex;
                
                int zlast = last / (gridsizex * gridsizey);
                int xlast = (last - zlast * gridsizex * gridsizey) / gridsizex;
                int ylast = last - zlast * gridsizex * gridsizey - xlast*gridsizex;

                // Contiguous: one node's labels now sit next to each other.
                const vidx nextBase = static_cast<vidx>(next) * searchvol;
                const vidx lastBase = static_cast<vidx>(last) * searchvol;
                for (vidx l = 0; l < searchvol; l++) {
                    energy[l] = costMatrix[nextBase + l];
                }

                float deltax = (fldx0[ylast + xlast * gridsizex + zlast * gridsizex * gridsizey] - fldx0[ynext + xnext * gridsizex + znext * gridsizex * gridsizey]) / (float) gridmin;
                float deltay = (fldy0[ylast + xlast * gridsizex + zlast * gridsizex * gridsizey] - fldy0[ynext + xnext * gridsizex + znext * gridsizex * gridsizey]) / (float) gridmin;
                float deltaz = (fldz0[ylast + xlast * gridsizex + zlast * gridsizex * gridsizey] - fldz0[ynext + xnext * gridsizex + znext * gridsizex * gridsizey]) / (float) gridmin;

                regulate(energy, index, searchsize, deltax, deltay, deltaz);

                for (vidx l = 0; l < searchvol; l++) {
                    indexall[nextBase + l] = index[l];
                }

                for (vidx l = 0; l < searchvol; l++) {
                    costMatrix[lastBase + l] += energy[l];
                }

            }

            //compuate label
            int i = 0;
            int oroot = sorted[i];
            int z1 = oroot / (gridsizex * gridsizey);
            int x1 = (oroot - z1 * gridsizex * gridsizey) / gridsizex;
            int y1 = oroot - z1 * gridsizex * gridsizey - x1*gridsizex;
            const vidx orootBase = static_cast<vidx>(oroot) * searchvol;
            for (vidx l = 0; l < searchvol; l++) {
                energy[l] = costMatrix[orootBase + l];
            }
            float value;
            int ind;
            minimumInd(energy, value, ind, searchvol);
            for (vidx l = 0; l < searchvol; l++) {
                indexall[orootBase + l] = l;
            }
            labels[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = ind;
            fldx[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = xs[ind] + fldx0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];
            fldy[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = ys[ind] + fldy0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];
            fldz[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = zs[ind] + fldz0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];


            //compute displacements 
            for (vidx i = 1; i < gridsizevol; i++) {
                int next = sorted[i];
                int last = relations[sorted[i]];
                int z1 = next / (gridsizex * gridsizey);
                int x1 = (next - z1 * gridsizex * gridsizey) / gridsizex;
                int y1 = next - z1 * gridsizex * gridsizey - x1*gridsizex;
                int z2 = last / (gridsizex * gridsizey);
                int x2 = (last - z2 * gridsizex * gridsizey) / gridsizex;
                int y2 = last - z2 * gridsizex * gridsizey - x2*gridsizex;

                ind = indexall[static_cast<vidx>(next) * searchvol
                               + labels[y2 + x2 * gridsizex + z2 * gridsizex * gridsizey]];
                labels[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = ind;
                fldx[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = xs[ind] + fldx0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];
                fldy[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = ys[ind] + fldy0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];
                fldz[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey] = zs[ind] + fldz0[y1 + x1 * gridsizex + z1 * gridsizex * gridsizey];

            }

            delete[] energy;
            delete[] values;
            delete[] index;
            delete[] indexall;
            delete[] labels;

            //return 0;
        }

        void
        optimization(MrfContext & ctx, float* imgFix, float* imgMov, float* energy, float lambda, int hw, float gridstep, float mingrid) {
            using namespace std;

            float* fixed;
            float* moving = imgMov;

            bool subpixel = false;
            if (mingrid == 0.5)
                subpixel = true;

            float gamma = (float) gridstep / (lambda * (float) mingrid);

            vidx m = ctx.sizeX;
            vidx n = ctx.sizeY;
            vidx o = ctx.sizeZ;
            vidx sz = m*n*o;

            int gridx = static_cast<int>(m / gridstep);
            int gridy = static_cast<int>(n / gridstep);
            int gridz = static_cast<int>(o / gridstep);
            // node count; multiplied by the label count below as a cost-volume stride.
            vidx sz1 = static_cast<vidx>(gridx)*gridy*gridz;

            int len = hw * 2 + 1;
            int lengthfull = (int) pow(len, 3.0);
            float* xs = new float[lengthfull];
            float* ys = new float[lengthfull];
            float* zs = new float[lengthfull];

            for (int i = 0; i < len; i++) {
                for (int j = 0; j < len; j++) {
                    for (int k = 0; k < len; k++) {
                        xs[i + j * len + k * len * len] = (float) ((j - hw) * mingrid);
                        ys[i + j * len + k * len * len] = (float) ((i - hw) * mingrid);
                        zs[i + j * len + k * len * len] = (float) ((k - hw) * mingrid);

                    }
                }
            }

            int hw_local;
            if (subpixel)
                hw_local = hw;
            else
                hw_local = hw * (int) mingrid;

            float* movingi;
            int* movingiInd; //array of the index of bin for the histgram computing
            int* fixedInd;

            const float binSizeFixed = ctx.maxFixed / ctx.numBin;
            const float binSizeMoving = ctx.maxMoving / ctx.numBin;

            vidx mi = m;
            vidx ni = n;
            vidx oi = o;

            if (subpixel) {
                //interpolation 
                mi *= 2;
                ni *= 2;
                oi *= 2;
                vidx szi = mi*ni*oi;
                float* x1 = new float[szi];
                float* y1 = new float[szi];
                float* z1 = new float[szi];
                movingi = new float[szi];
                movingiInd = new int[szi];
                for (int k = 0; k < oi; k++) {
                    for (int j = 0; j < ni; j++) {
                        for (int i = 0; i < mi; i++) {
                            x1[i + j * mi + k * mi * ni] = 0.5 * (float) j;
                            y1[i + j * mi + k * mi * ni] = 0.5 * (float) i;
                            z1[i + j * mi + k * mi * ni] = 0.5 * (float) k;
                        }
                    }
                }
                interpolate3d(movingi, moving, x1, y1, z1, mi, ni, oi, m, n, o, false);
                delete []x1;
                delete []y1;
                delete []z1;
                for (int i = 0; i < lengthfull; i++) {
                    xs[i] *= 2.0;
                    ys[i] *= 2.0;
                    zs[i] *= 2.0;
                }


            } else {
                movingi = new float[sz];
                movingiInd = new int[sz];
                for (vidx i = 0; i < sz; i++) {
                    movingi[i] = moving[i];
                }


            }

            fixed = new float[m * n * o];
            fixedInd = new int[m * n * o];
            for (vidx i = 0; i < m * n * o; i++) {

                fixed[i] = imgFix[i];
                fixedInd[i] = (int) floor(fixed[i] / binSizeFixed);
            }

            // The label offsets xs/ys/zs are fixed for the whole call, and in subpixel mode
            // already carry the x2 scaling applied above. Fold each label's three offsets into
            // a single flat buffer offset so the interior case below costs one add per label
            // instead of three float-to-int conversions and two multiplies.
            std::vector<int>  xsi(lengthfull), ysi(lengthfull), zsi(lengthfull);
            std::vector<vidx> labelOffset(lengthfull);
            for (int l = 0; l < lengthfull; l++) {
                xsi[l] = (int) xs[l];
                ysi[l] = (int) ys[l];
                zsi[l] = (int) zs[l];
                labelOffset[l] = static_cast<vidx>(ysi[l]) + static_cast<vidx>(xsi[l]) * mi
                               + static_cast<vidx>(zsi[l]) * mi * ni;
            }

            int samples = ctx.sampleNum;
            bool randommode = samples < pow(gridstep, 3);
            int maxsamp;
            if (randommode) {
                maxsamp = samples;
            } else {
                maxsamp = (int) pow(gridstep, 3);
            }
            float epsilon = gamma / (float) maxsamp;

            {
                // The sample indices come off the global rand() stream, which is sequential
                // state: drawing them inside a parallel loop would both race and reorder them.
                // Drawing them all up front, in the original order, lets the parallel loop
                // reproduce the serial result exactly rather than merely equivalently.
                std::vector<int> sampleIdx;
                if (randommode) {
                    const int cellCount = (int) pow(gridstep, 3);
                    sampleIdx.resize(static_cast<size_t>(sz1) * maxsamp);
                    for (size_t sIdx = 0; sIdx < sampleIdx.size(); sIdx++) {
                        int i1 = (int) (rand() * pow(gridstep, 3) / float(RAND_MAX));
                        if (i1 >= cellCount) { i1 = cellCount - 1; }
                        if (i1 < 0) { i1 = 0; }
                        sampleIdx[sIdx] = i1;
                    }
                }

                // Each node writes only energy[i + l*sz1], so nodes are independent and the
                // chunks below are disjoint in both their reads and their writes.
                parallelFor(sz1, [&](vidx nodeBegin, vidx nodeEnd) {

                // Scratch the serial version shared has to become per-chunk.
                std::vector<float> cost1Buf(lengthfull), normBBuf(lengthfull);
                float* cost1 = cost1Buf.data();
                float* normalizerB = normBBuf.data();
                vidx xind, yind, zind;

                for (vidx i = nodeBegin; i < nodeEnd; i++) {
                    vidx z1 = i / (static_cast<vidx>(gridx) * gridy);
                    vidx x1 = (i - z1 * gridx * gridy) / gridx;
                    vidx y1 = i - z1 * gridx * gridy - x1*gridx;

                    z1 *= gridstep;
                    x1 *= gridstep;
                    y1 *= gridstep;

                    bool boundaries = true; //check boundaries 
                    if (subpixel) {
                        if (x1 * 2 + (gridstep - 1)*2 + hw_local >= ni | y1 * 2 + (gridstep - 1)*2 + hw_local >= mi | z1 * 2 + (gridstep - 1)*2 + hw_local >= oi)
                            boundaries = false;
                        if (x1 * 2 - hw_local < 0 | y1 * 2 - hw_local < 0 | z1 * 2 - hw_local < 0)
                            boundaries = false;
                    }
                    else {
                        if (x1 + (gridstep - 1) + hw_local >= ni | y1 + (gridstep - 1) + hw_local >= mi | z1 + (gridstep - 1) + hw_local >= oi)
                            boundaries = false;
                        if (x1 - hw_local < 0 | y1 - hw_local < 0 | z1 - hw_local < 0)
                            boundaries = false;
                    }


                    // normalizerA was an array indexed by label, but what it accumulates --
                    // the squared fixed-image intensity -- has no label dependence whatsoever,
                    // so every entry held the same number. One scalar replaces the array, and
                    // the accumulation moves out of the label loop entirely.
                    float normA = 0.0f;
                    for (int l = 0; l < lengthfull; l++) {
                        cost1[l] = 0.0;
                        normalizerB[l] = 0.0;
                    }


                    for (int j1 = 0; j1 < maxsamp; j1++) {
                        const int i1 = randommode
                                         ? sampleIdx[static_cast<size_t>(i) * maxsamp + j1]
                                         : j1;
                        vidx zz = i1 / (gridstep * gridstep);
                        vidx xx = (i1 - zz * gridstep * gridstep) / gridstep;
                        vidx yy = i1 - zz * gridstep * gridstep - xx*gridstep;

                        xx += x1;
                        yy += y1;
                        zz += z1;

                        // One load and one square of the fixed sample per sample, rather than
                        // repeating both for every one of the (up to 29791) labels.
                        const float fv = fixed[yy + xx * m + zz * m * n];
                        normA += fv * fv;

                        const vidx bx = subpixel ? xx * 2 : xx;
                        const vidx by = subpixel ? yy * 2 : yy;
                        const vidx bz = subpixel ? zz * 2 : zz;

                        if (boundaries) {
                            // Interior node: every label lands inside the volume, so both the
                            // clamp and the per-label index arithmetic fall away and the moving
                            // sample is one offset from a fixed base.
                            const vidx base = by + bx * mi + bz * mi * ni;
                            for (int l = 0; l < lengthfull; l++) {
                                const float mv = movingi[base + labelOffset[l]];
                                cost1[l] += fv * mv;
                                normalizerB[l] += mv * mv;
                            }
                        } else {
                            for (int l = 0; l < lengthfull; l++) {
                                xind = clampIdx(bx + xsi[l], ni - 1);
                                yind = clampIdx(by + ysi[l], mi - 1);
                                zind = clampIdx(bz + zsi[l], oi - 1);
                                const float mv = movingi[yind + xind * mi + zind * mi * ni];
                                cost1[l] += fv * mv;
                                normalizerB[l] += mv * mv;
                            }
                        }
                    }

                    // Contiguous run of labelCount floats for this node, rather than
                    // lengthfull writes each one node-count apart.
                    const vidx energyBase = i * lengthfull;
                    for (int l = 0; l < lengthfull; l++) {
                        // Squaring, not a general pow. The ratio stays in float exactly as
                        // upstream computed it, so the rounding sequence is unchanged.
                        const float q = cost1[l] * cost1[l] / normA / normalizerB[l];
                        energy[energyBase + l] = epsilon * (1 - double(q) * double(q));
                    }
                }
                });
            }

            delete []movingi;
            delete []fixed;
            delete []movingiInd;
            delete []fixedInd;

            delete []xs;
            delete []ys;
            delete []zs;

            //return 0;
        }

        void
        updateImage(MrfContext & ctx, float* warped, float* imgFix, float* imgMov, float* u1, float* v1, float* w1) {
            vidx m = ctx.sizeX;
            vidx n = ctx.sizeY;
            vidx o = ctx.sizeZ;


            float metric = 0;
            float metric_start = 0;


            interpolate3d(warped, imgFix, u1, v1, w1, m, n, o, m, n, o, true);

            for (vidx i = 0; i < m; i++) {
                for (vidx j = 0; j < n; j++) {
                    for (vidx k = 0; k < o; k++) {
                        metric += pow(imgMov[i + j * m + k * m * n] - warped[i + j * m + k * m * n], 2);
                        metric_start += pow(imgMov[i + j * m + k * m * n] - imgFix[i + j * m + k * m * n], 2);
                    }
                }
            }

            metric /= m * n*o;
            metric_start /= m * n*o;

            ctx.metricBefore = metric_start;
            ctx.metricAfter = metric;
        }

} // anonymous namespace

// ---------------------------------------------------------------------------
// Orchestrator. Replaces the upstream deformableReg(), which read an atlas and a
// model directory off the LIMS dispatcher; the coarse-to-fine schedule, the
// forward/backward pairing and the symmetric composition are unchanged.
// ---------------------------------------------------------------------------
long long
EstimateDeformablePeakBytes(int sizeX, int sizeY, int sizeZ, const DeformableParameters & params)
{
  const long long vox = static_cast<long long>(sizeX) * sizeY * sizeZ;

  // Mirror the schedule RunDeformableRegistration derives: ceil-halve the grid each
  // level, and stop once the spacing bottoms out and would repeat.
  long long worstLevel = 0;
  int       step = params.gridStart > 0 ? params.gridStart : 1;
  int       previous = 0;
  for (int level = 0; level < params.numLevels; ++level)
  {
    if (level > 0)
    {
      previous = step;
      step = static_cast<int>(std::ceil(step / 2.0));
      if (step == previous) { break; }
    }
    const long long hw = static_cast<long long>(std::ceil(step / 2.0));
    const long long side = 2 * hw + 1;
    const long long labels = side * side * side;
    const long long nodes = std::max(1LL, static_cast<long long>(sizeX / step)) *
                            std::max(1LL, static_cast<long long>(sizeY / step)) *
                            std::max(1LL, static_cast<long long>(sizeZ / step));
    // forward + backward cost volume, float each
    worstLevel = std::max(worstLevel, 2LL * 4LL * nodes * labels);
  }

  return 190LL * vox + worstLevel + 20LL * 1024 * 1024;
}

// ---------------------------------------------------------------------------
DeformableResult
RunDeformableRegistration(const float *                imgFixIn,
                          const float *                imgMovIn,
                          int                          sizeX,
                          int                          sizeY,
                          int                          sizeZ,
                          const DeformableParameters & params)
{
  const vidx volsize = static_cast<vidx>(sizeX) * sizeY * sizeZ;

  MrfContext ctx;
  ctx.sizeX = sizeX;
  ctx.sizeY = sizeY;
  ctx.sizeZ = sizeZ;
  ctx.numBin = params.numBin;
  ctx.sampleNum = params.sampleNum;

  // Upstream assumes 8-bit data and hard-codes 255 here. Derive the real maxima so the
  // similarity histogram bins correctly whatever the input range is.
  ctx.maxFixed = 0.0f;
  ctx.maxMoving = 0.0f;
  for (vidx i = 0; i < volsize; ++i)
  {
    ctx.maxFixed = max(ctx.maxFixed, imgFixIn[i]);
    ctx.maxMoving = max(ctx.maxMoving, imgMovIn[i]);
  }
  if (ctx.maxFixed <= 0.0f) { ctx.maxFixed = 1.0f; }
  if (ctx.maxMoving <= 0.0f) { ctx.maxMoving = 1.0f; }

  // Coarse-to-fine control-grid spacings: 30, 15, 8, 4 upstream (ceil-halving).
  // Cap the starting spacing so the coarsest grid still has several nodes per axis;
  // otherwise small images collapse to a 1-node grid and the optimization is vacuous.
  const int shortestAxis = min(min(sizeX, sizeY), sizeZ);
  int       gridStart = params.gridStart;
  const int kMinNodesPerAxis = 4;
  if (gridStart > shortestAxis / kMinNodesPerAxis)
  {
    gridStart = max(2, shortestAxis / kMinNodesPerAxis);
    std::cerr << "warning: reducing deformable grid spacing " << params.gridStart << " -> "
              << gridStart << " so the coarsest control grid keeps at least "
              << kMinNodesPerAxis << " nodes per axis (shortest axis is " << shortestAxis
              << ")." << std::endl;
  }

  int numLevels = max(1, params.numLevels);
  std::vector<int> gridSize(numLevels);
  gridSize[0] = gridStart;
  for (int i = 1; i < numLevels; ++i)
  {
    gridSize[i] = static_cast<int>(std::ceil(gridSize[i - 1] / 2.0));
  }
  // Drop levels whose spacing has bottomed out at 1 and would repeat.
  while (numLevels > 1 && gridSize[numLevels - 1] == gridSize[numLevels - 2])
  {
    --numLevels;
    gridSize.resize(numLevels);
  }

  std::vector<int> gridHw(numLevels);
  for (int i = 0; i < numLevels; ++i)
  {
    gridHw[i] = static_cast<int>(std::ceil(gridSize[i] / 2.0));
  }
  std::vector<float> labelSizeMin(numLevels, 1.0f);
  labelSizeMin[numLevels - 1] = 0.5f;

  // Working copies: the routines take non-const buffers.
  std::vector<float> imgFix(imgFixIn, imgFixIn + volsize);
  std::vector<float> imgMov(imgMovIn, imgMovIn + volsize);
  std::vector<float> warpedImgFix(volsize, 0.0f);
  std::vector<float> warpedImgMov(volsize, 0.0f);

  std::vector<float> xfld(volsize, 0.0f), yfld(volsize, 0.0f), zfld(volsize, 0.0f);
  std::vector<float> xfldinv(volsize, 0.0f), yfldinv(volsize, 0.0f), zfldinv(volsize, 0.0f);

  int xgridsize = sizeX / gridSize[0];
  int ygridsize = sizeY / gridSize[0];
  int zgridsize = sizeZ / gridSize[0];
  std::vector<float> xlabel(xgridsize * ygridsize * zgridsize, 0.0f);
  std::vector<float> ylabel(xlabel.size(), 0.0f);
  std::vector<float> zlabel(xlabel.size(), 0.0f);
  std::vector<float> xlabelinv(xlabel.size(), 0.0f);
  std::vector<float> ylabelinv(xlabel.size(), 0.0f);
  std::vector<float> zlabelinv(xlabel.size(), 0.0f);

  for (int level = 0; level < numLevels; ++level)
  {
    if (params.verbose) { std::cerr << "      [step] warp fixed/moving" << std::endl; }
    updateImage(ctx, warpedImgFix.data(), imgMov.data(), imgFix.data(), xfld.data(), yfld.data(), zfld.data());
    updateImage(ctx, warpedImgMov.data(), imgFix.data(), imgMov.data(), xfldinv.data(), yfldinv.data(), zfldinv.data());

    const float minsize = labelSizeMin[level];
    const int   stepsize = gridSize[level];
    const int   hwidth = gridHw[level];
    const int   labelCount = static_cast<int>(std::pow(hwidth * 2 + 1, 3.0));

    const int xsizethis = sizeX / stepsize;
    const int ysizethis = sizeY / stepsize;
    const int zsizethis = sizeZ / stepsize;
    const int volsizethis = xsizethis * ysizethis * zsizethis;
    if (volsizethis <= 0)
    {
      break;
    }

    // Carry the previous level's solution onto this level's finer grid.
    if (params.verbose) { std::cerr << "      [step] interp previous level" << std::endl; }
    std::vector<float> xfldstart(volsizethis), yfldstart(volsizethis), zfldstart(volsizethis);
    std::vector<float> xfldstartinv(volsizethis), yfldstartinv(volsizethis), zfldstartinv(volsizethis);
    InterpDeformationFld(xfldstart.data(), yfldstart.data(), zfldstart.data(),
                         xlabel.data(), ylabel.data(), zlabel.data(),
                         xsizethis, ysizethis, zsizethis, xgridsize, ygridsize, zgridsize);
    InterpDeformationFld(xfldstartinv.data(), yfldstartinv.data(), zfldstartinv.data(),
                         xlabelinv.data(), ylabelinv.data(), zlabelinv.data(),
                         xsizethis, ysizethis, zsizethis, xgridsize, ygridsize, zgridsize);

    if (params.verbose)
    {
      std::cout << "  deformable level " << level << ": grid spacing " << stepsize
                << ", nodes " << xsizethis << "x" << ysizethis << "x" << zsizethis
                << ", labels " << labelCount << ", samples " << ctx.sampleNum << std::endl;
    }

    xlabel.assign(volsizethis, 0.0f);
    ylabel.assign(volsizethis, 0.0f);
    zlabel.assign(volsizethis, 0.0f);
    xlabelinv.assign(volsizethis, 0.0f);
    ylabelinv.assign(volsizethis, 0.0f);
    zlabelinv.assign(volsizethis, 0.0f);

    // Minimum spanning tree over the control grid, in both directions.
    if (params.verbose) { std::cerr << "      [step] mst (DfmfldModulation)" << std::endl; }
    std::vector<int> optimizedForward(volsizethis), relationForward(volsizethis);
    std::vector<int> optimizedBackward(volsizethis), relationBackward(volsizethis);
    DfmfldModulation(ctx, imgFix.data(), optimizedForward.data(), relationForward.data(), stepsize);
    DfmfldModulation(ctx, imgMov.data(), optimizedBackward.data(), relationBackward.data(), stepsize);

    // Discrete label costs, then message passing along the tree.
    if (params.verbose) { std::cerr << "      [step] cost (optimization)" << std::endl; }
    std::vector<float> costforward(static_cast<size_t>(volsizethis) * labelCount);
    std::vector<float> costbackward(static_cast<size_t>(volsizethis) * labelCount);
    optimization(ctx, imgFix.data(), warpedImgFix.data(), costforward.data(), params.lambda, hwidth, stepsize, minsize);
    optimization(ctx, imgMov.data(), warpedImgMov.data(), costbackward.data(), params.lambda, hwidth, stepsize, minsize);

    if (params.verbose) { std::cerr << "      [step] genDfm" << std::endl; }
    genDfm(ctx, xlabel.data(), ylabel.data(), zlabel.data(),
           xfldstart.data(), yfldstart.data(), zfldstart.data(),
           costforward.data(), params.lambda, hwidth, stepsize, minsize,
           optimizedForward.data(), relationForward.data());
    genDfm(ctx, xlabelinv.data(), ylabelinv.data(), zlabelinv.data(),
           xfldstartinv.data(), yfldstartinv.data(), zfldstartinv.data(),
           costbackward.data(), params.lambda, hwidth, stepsize, minsize,
           optimizedBackward.data(), relationBackward.data());

    if (params.verbose) { std::cerr << "      [step] smoothDfm" << std::endl; }
    smoothDfm(xlabel.data(), ylabel.data(), zlabel.data(), xsizethis, ysizethis, zsizethis, 4, stepsize);
    smoothDfm(xlabelinv.data(), ylabelinv.data(), zlabelinv.data(), xsizethis, ysizethis, zsizethis, 4, stepsize);
    // The step that makes the result symmetric and invertible.
    if (params.verbose) { std::cerr << "      [step] symmetricMapping" << std::endl; }
    symmetricMapping(xlabel.data(), ylabel.data(), zlabel.data(),
                     xlabelinv.data(), ylabelinv.data(), zlabelinv.data(),
                     xsizethis, ysizethis, zsizethis, stepsize);

    if (params.verbose) { std::cerr << "      [step] interp to full res" << std::endl; }
    InterpDeformationFld(xfld.data(), yfld.data(), zfld.data(),
                         xlabel.data(), ylabel.data(), zlabel.data(),
                         sizeX, sizeY, sizeZ, xsizethis, ysizethis, zsizethis);
    InterpDeformationFld(xfldinv.data(), yfldinv.data(), zfldinv.data(),
                         xlabelinv.data(), ylabelinv.data(), zlabelinv.data(),
                         sizeX, sizeY, sizeZ, xsizethis, ysizethis, zsizethis);

    updateImage(ctx, warpedImgFix.data(), imgMov.data(), imgFix.data(), xfld.data(), yfld.data(), zfld.data());
    if (params.verbose)
    {
      std::cout << "    mean squared difference " << ctx.metricBefore << " -> " << ctx.metricAfter
                << std::endl;
    }

    xgridsize = xsizethis;
    ygridsize = ysizethis;
    zgridsize = zsizethis;

    ctx.sampleNum = ctx.sampleNum / 4;
    if (ctx.sampleNum < 64) { ctx.sampleNum = 64; }
  }

  DeformableResult result;
  result.warped.swap(warpedImgFix);
  // Upstream names the two in-plane field components the other way round from the buffer
  // layout: interpolate3d adds the fastest loop index to the component it calls "y".
  // Export in x-fastest (ITK) order.
  result.dispX.swap(yfld);
  result.dispY.swap(xfld);
  result.dispZ.swap(zfld);
  // Same component transposition for the backward field -- it comes off the identical
  // interpolate3d convention, so it must get the identical swap or it exports transposed.
  result.invX.swap(yfldinv);
  result.invY.swap(xfldinv);
  result.invZ.swap(zfldinv);
  result.metricBefore = ctx.metricBefore;
  result.metricAfter = ctx.metricAfter;
  result.levelsUsed = numLevels;
  result.gridStartUsed = gridStart;
  return result;
}

} // namespace stpt
