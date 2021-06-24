//  Copyright (C) 2018-2021  Pedro Gomes
//  See full notice in NOTICE.md

#include "functions.hpp"

#include <fstream>
#include <iostream>

#include <cmath>
#include <algorithm>

#if defined HAVE_MKL
#include <mkl.h>
#else
#include <lapacke.h>
#endif
#include <omp.h>
#include <x86intrin.h>

using namespace std;


int readInput(size_t nPt,
              size_t nRhs,
              vector<type> & x,
              vector<type> & y,
              vector<type> & z,
              vector<type> & b)
{
    size_t bytes = nPt*sizeof(type);

    // read point coordinates
    fstream file;
    file.open(coordsFile, ios::in | ios::binary);
    if(file.bad()) return -1;

    x.resize(nPt);  y.resize(nPt);  z.resize(nPt);

    file.read(reinterpret_cast<char*>(x.data()), bytes);
    file.read(reinterpret_cast<char*>(y.data()), bytes);
    file.read(reinterpret_cast<char*>(z.data()), bytes);
    if(file.bad()) return -1;

    // read point displacements (the RBF RHS)
    file.close();
    file.open(displsFile, ios::in | ios::binary);
    if(file.bad()) return -1;

    b.resize(nPt*nRhs);

    file.read(reinterpret_cast<char*>(b.data()), bytes*nRhs);
    if(file.bad()) return -1;

    return 0;
}

int writeOutput(std::vector<type> const & sol)
{
    fstream file;
    file.open(ctrlptFile, ios::out | ios::binary);
    if(file.bad()) return -1;

    file.write(reinterpret_cast<const char*>(sol.data()), sol.size()*sizeof(type));
    if(file.bad()) return -1;

    return 0;
}

void buildKernel(vector<type> const & x,
                 vector<type> const & y,
                 vector<type> const & z,
                 vector<type> & M)
{
    const size_t n = x.size();

    // info about the size needed for the RFP format
    const size_t odd = n&1;
    const size_t even = 1-odd;
    const size_t k = n/2;

    M.resize(n*(n+1)/2);

    #define WENDLAND(C) {               \
    type d = sqrt(pow(x[i]-x[j], 2)+    \
                  pow(y[i]-y[j], 2)+    \
                  pow(z[i]-z[j], 2));   \
    type t = pow(max(0.0,1.0-d),2);     \
    C = t*t*(1.0+4.0*d);}

    //  This is how lower triangular indices ij map to a
    //  column major RFP matrix's 1D storage array:
    //    if(j<k+odd) idx = j*(n+even) + i+even;
    //    else        idx = (i-k)*(n+even) + j-k-odd;
    //  M is built in two loops with swaped inner/outer
    //  indices to keep indexing linearly into it.

    #pragma omp parallel for schedule(dynamic,512)
    for(size_t j=0; j<k+odd; ++j) // all columns
    {
        // index of diagonal entry in the 1D RFP array
        size_t idx = j*(n+even+1)+even;

        M[idx] = 1.0; // by definition

        // lower rows, offset by j since we START at the diagonal
        for(size_t i=j+1; i<n; ++i) WENDLAND(M[idx+i-j]);
    }

    #pragma omp parallel for schedule(dynamic,512)
    for(size_t i=k+odd; i<n; ++i) // bottom rows
    {
        // index of diagonal entry in the 1D RFP array
        size_t idx = (i-k)*(n+even+1)-odd;

        // left columns, offset by i since we END at the diagonal
        for(size_t j=k+odd; j<i; ++j) WENDLAND(M[idx+j-i]);

        M[idx] = 1.0; // by definition
    }

    #undef WENDLAND
}

void kernelProduct(size_t m,
                   size_t n,
                   double const * x,
                   double const * y,
                   double const * z,
                   double const * A,
                   double * B)
{
    // AVX register size
    const size_t SIZE = 4;

    // fill two registers with the constants we need
    __m256d vone  = _mm256_set1_pd(1.0),
            vfour = _mm256_set1_pd(4.0);

    // vectorizable loop size (multiple of AVX register size)
    size_t m_simd = (m/SIZE)*SIZE;

    // scalar (last) part of the loop (easier to read, so shown first)
    for(size_t i=m_simd; i<m; ++i)
    {
        // initialize products to 0
        for(size_t k=0; k<n; ++k) B[i+k*m] = 0.0;

        for(size_t j=0; j<m; ++j)
        {
            // compute distance from i to j
            double d = sqrt(pow(x[i]-x[j],2)+
                            pow(y[i]-y[j],2)+
                            pow(z[i]-z[j],2));

            // compute the RBF kernel coefficient
            double t = pow(max(0.0,1.0-d),2);
            t *= t*(1.0+4.0*d);

            // update products
            for(size_t k=0; k<n; ++k) B[i+k*m] += t*A[j+k*m];
        }
    }

    // as above but explicitly vectorized, multiple i's at a time
    #pragma omp parallel for schedule(static,512)
    for(size_t i=0; i<m_simd; i+=SIZE)
    {
        // coordinates of point(s) i
        __m256d
        xi = _mm256_loadu_pd(&x[i]),
        yi = _mm256_loadu_pd(&y[i]),
        zi = _mm256_loadu_pd(&z[i]);

        // initialize products to 0
        for(size_t k=0; k<n; ++k)
            _mm256_storeu_pd(&B[i+k*m], _mm256_setzero_pd());

        // compute products
        for(size_t j=0; j<m; ++j)
        {
            // compute distance from point(s) i to point j
            __m256d
            a = _mm256_set1_pd(x[j]),
            t = _mm256_sub_pd(xi,a),
            d = _mm256_mul_pd(t,t);     // d^2  = (xi-xj)^2

            a = _mm256_set1_pd(y[j]);
            t = _mm256_sub_pd(yi,a);
            d = _mm256_fmadd_pd(t,t,d); // d^2 += (yi-yj)^2

            a = _mm256_set1_pd(z[j]);
            t = _mm256_sub_pd(zi,a);
            d = _mm256_fmadd_pd(t,t,d); // d^2 += (zi-zj)^2

            d = _mm256_sqrt_pd(d);      // d = sqrt(d^2)

            // compute the RBF kernel coefficient, registers
            // are re-used, final result is in "t"
            t = _mm256_min_pd(vone,d);  // t = min(1,d)
            t = _mm256_sub_pd(vone,t);  // t = 1-t <=> max(0,1-d)
            t = _mm256_mul_pd(t,t);     // t*= t   <=> t^2
            t = _mm256_mul_pd(t,t);     // t*= t   <=> t^4
            d = _mm256_fmadd_pd(vfour,d,vone); // d = 1+4*d
            t = _mm256_mul_pd(t,d);     // t*= d   <=> max(0,1-d)^4*(1+4*d)

            // update products e.g. B += t*A
            for(size_t k=0; k<n; ++k)
            {
                d = _mm256_loadu_pd(&B[i+k*m]); // read
                a = _mm256_set1_pd(A[j+k*m]);
                d = _mm256_fmadd_pd(t,a,d);     // modify
                _mm256_storeu_pd(&B[i+k*m], d); // write
            }
        }
    }
}

int solveRBF(size_t& iters,
             type& error,
             vector<type> const & x,
             vector<type> const & y,
             vector<type> const & z,
             vector<type> & res,
             vector<type> & lhs)
{
    #define TIC t0 = omp_get_wtime()
    #define TOC omp_get_wtime() - t0

    const size_t maxIters = iters;
    const type tol = error;
    double t0;

    size_t m = x.size(), n = res.size()/m;

    cout << "   RBF solver started\n";

    // init lhs and compute rhs norm
    lhs.clear();
    lhs.resize(m*n,0.0);

    type bNorm = 0.0;
    for(auto r : res) bNorm += r*r;
    bNorm = sqrt(bNorm);


    cout << "    Building RBF kernel...  " << flush;
    vector<type> ap;
    TIC;
    buildKernel(x,y,z,ap);
    cout << TOC << "s\n";


    cout << "    Computing LLT factr...  " << flush;
    TIC;
    int info = LAPACKE_dpftrf(LAPACK_COL_MAJOR,'N','L',m,ap.data());
    if(info!=0) return -1;
    cout << TOC << "s\n";


    cout << "    Iterative refinement... " << flush;
    TIC;
    for(iters=0; iters<=maxIters; ++iters)
    {
        // compute and apply correction
        vector<type> corr(m*n);
        #pragma omp simd
        for(size_t i=0; i<m*n; ++i) corr[i] = res[i];

        info = LAPACKE_dpftrs(LAPACK_COL_MAJOR,'N','L',m,n,ap.data(),corr.data(),m);
        if(info!=0) return -1;

        #pragma omp simd
        for(size_t i=0; i<m*n; ++i) lhs[i] += corr[i];

        // update and compute residual
        vector<type> delta(m*n);
        kernelProduct(m, n, x.data(), y.data(), z.data(), corr.data(), delta.data());

        error = 0.0;
        #pragma omp simd
        for(size_t i=0; i<m*n; ++i)
        {
            res[i] -= delta[i];
            error += pow(res[i],2);
        }
        error = sqrt(error)/bNorm;

        if(error < tol) break;
    }
    cout << TOC << "s\n";

    cout << "   RBF solver finished\n";
    cout << "    Refinement iters:       " << iters << "\n";
    cout << "    Norm of residuals:      " << error << "\n";

    return 0;
}
