//  Copyright (C) 2018-2021  Pedro Gomes
//  See full notice in NOTICE.md

#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <omp.h>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <x86intrin.h>

#include "definitions.hpp"

using namespace std;

/*!
 * Read input files (point coordinates and displacements)
 */
int readInput(size_t nPt, size_t nCP, size_t nDim, type* Xcp, type* Ucp, type* X);

/*!
 * Write output file (modified coordinates)
 */
int writeOutput(size_t nPt, const type* X);

/*!
 * Apply the RBF deformation defined by control point coordinates
 * Xcp and displacements Ucp to the points with coordinates X.
 * nDim is the number of coordinates being moved (i.e. number of
 * columns of Ucp, 2 or 3), Xcp and X are assumed to be 3D with the
 * xyz coordinates in column major storage (xxx,yyy,zzz). All need
 * to be 32B aligned with the columns padded to a SIMD len multiple.
 */
template<size_t nDim>
void deform(size_t nPt, size_t nCP, const double* Xcp, const double* Ucp, double* X)
{
    static_assert(nDim==2 || nDim==3, "Invalid number of dimensions.");

    double t0 = omp_get_wtime();

    if(nPt > 10000) {
        stringstream ss;
        ss << "   Applying " << nCP << "*" << nDim << " RBF to " << nPt << " points";
        cout << setfill('.') << setw(44) << left << ss.str() << flush;
    }

    // tile ("cache block") sizes
    const size_t TILEH = 512;
    const size_t TILEW = 256;

    // fill two registers with the constants we need
    __m256d vone  = _mm256_set1_pd(1.0),
            vfour = _mm256_set1_pd(4.0);

    // explicitly vectorized (multiple i's at a time) and blocked
    // outer most loop over sets of #TILEH points
    #pragma omp parallel for schedule(static,1)
    for(size_t i0=0; i0<nPt; i0+=TILEH)
    {
        // displacement buffer, due to tiling, coordinates cannot be modified in-place
        alignas(64) double U[nDim*TILEH*sizeof(double)];

        // zero the tile displacements
        for(size_t k=0; k<nDim*TILEH; k+=SIMDLEN)
            _mm256_store_pd(&U[k], _mm256_setzero_pd());

        // stop point for inner i loop
        size_t iend = min(i0+TILEH,nPt);

        // second loop over sets of #TILEW control points
        for(size_t j0=0; j0<nCP; j0+=TILEW)
        {
            // stop point for inner j loop
            size_t jend = min(j0+TILEW,nCP);

            // inner i loop, over #TILEH points
            // vectorized, operate on #SIMDLEN points at a time
            for(size_t i=i0; i<iend; i+=SIMDLEN)
            {
                // coordinates and displacements of points i
                __m256d wi,
                xi = _mm256_load_pd(&X[   i   ]),
                yi = _mm256_load_pd(&X[ i+nPt ]),
                zi = _mm256_load_pd(&X[i+2*nPt]),
                ui = _mm256_load_pd(&U[    i-i0    ]),
                vi = _mm256_load_pd(&U[ i-i0+TILEH ]);
                if(nDim==3)
                wi = _mm256_load_pd(&U[i-i0+2*TILEH]);

                // inner j loop, over #TILEW control points, vectorized loads
                for(size_t j=j0; j<jend; j+=SIMDLEN)
                {
                    // coordinates of points j, and an aux var
                    __m256d t,
                    xj = _mm256_load_pd(&Xcp[   j   ]),
                    yj = _mm256_load_pd(&Xcp[ j+nCP ]),
                    zj = _mm256_load_pd(&Xcp[j+2*nCP]);

                    // compute the possible 4 distances from i to j...
                    #define COMPUTE_DIST(D) __m256d                         \
                    D = _mm256_sub_pd(xi,xj);  D = _mm256_mul_pd(D,D);      \
                    t = _mm256_sub_pd(yi,yj);  D = _mm256_fmadd_pd(t,t,D);  \
                    t = _mm256_sub_pd(zi,zj);  D = _mm256_fmadd_pd(t,t,D);  \
                    D = _mm256_sqrt_pd(D)

                    // ...by going through the different permutations
                    #define SHUFFLE(FUN,IMM8)   \
                    xj = FUN(xj,xj,IMM8);       \
                    yj = FUN(yj,yj,IMM8);       \
                    zj = FUN(zj,zj,IMM8)

                    COMPUTE_DIST(d0);

                    SHUFFLE(_mm256_shuffle_pd,0b0101);
                    COMPUTE_DIST(d1);

                    SHUFFLE(_mm256_permute2f128_pd,1);
                    COMPUTE_DIST(d2);

                    SHUFFLE(_mm256_shuffle_pd,0b0101);
                    COMPUTE_DIST(d3);

                    // coordinate registers now hold the displacements
                    xj = _mm256_load_pd(&Ucp[   j   ]),
                    yj = _mm256_load_pd(&Ucp[ j+nCP ]);
                    if(nDim==3)
                    zj = _mm256_load_pd(&Ucp[j+2*nCP]);

                    // coefficients for each set of distances...
                    #define COMPUTE_COEFF(C)                                \
                    t = _mm256_min_pd(vone,C);  t = _mm256_sub_pd(vone,t);  \
                    t = _mm256_mul_pd(t,t);     t = _mm256_mul_pd(t,t);     \
                    C = _mm256_fmadd_pd(vfour,C,vone);                      \
                    C = _mm256_mul_pd(t,C)

                    // ...+ update i point displacements
                    #define UPDATE_DISP(C)          \
                    COMPUTE_COEFF(C);               \
                    ui = _mm256_fmadd_pd(C,xj,ui);  \
                    vi = _mm256_fmadd_pd(C,yj,vi);  \
                    if(nDim==3) wi = _mm256_fmadd_pd(C,zj,wi)

                    UPDATE_DISP(d0);

                    SHUFFLE(_mm256_shuffle_pd,0b0101);
                    UPDATE_DISP(d1);

                    SHUFFLE(_mm256_permute2f128_pd,1);
                    UPDATE_DISP(d2);

                    SHUFFLE(_mm256_shuffle_pd,0b0101);
                    UPDATE_DISP(d3);
                }

                // store updated displacements
                _mm256_store_pd(&U[    i-i0    ], ui);
                _mm256_store_pd(&U[ i-i0+TILEH ], vi);
                if(nDim==3)
                _mm256_store_pd(&U[i-i0+2*TILEH], wi);
            }
        }

        // add tile displacements to the coordinates
        for(size_t k=0; k<nDim; ++k)
        {
            for(size_t i=i0; i<iend; i+=SIMDLEN)
            {
                __m256d
                x = _mm256_load_pd(&X[i+k*nPt]),
                u = _mm256_load_pd(&U[i-i0+k*TILEH]);
                x = _mm256_add_pd(x,u);
                _mm256_stream_pd(&X[i+k*nPt], x);
            }
        }
    }

    if(nPt > 10000) cout << " " << omp_get_wtime()-t0 << "s\n";
}

#endif // FUNCTIONS_HPP
