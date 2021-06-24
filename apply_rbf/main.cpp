//  Copyright (C) 2018-2021  Pedro Gomes
//  See full notice in NOTICE.md

#include <cstdlib>
#include <omp.h>

#include "definitions.hpp"
#include "functions.hpp"

#define CLEAN_AND_EXIT(RETCODE) {free(X); free(Xcp); free(Ucp); return RETCODE;}

int main(int argc, char* argv[])
{
    if(argc<5) return -1;

    const size_t nPt = atoi(argv[1]);
    const size_t nCP = atoi(argv[2]);
    const size_t nDim = atoi(argv[3]);
    const size_t nThrd = atoi(argv[4]);

    omp_set_num_threads(nThrd);

    // aligned and padded allocation for points so that each
    // coordinate can be accessed using aligned load/stores
    const size_t nPt_simd = ROUNDSIZE(nPt);
    auto X = static_cast<type*>(aligned_alloc(64,3*nPt_simd*sizeof(type)));

    const size_t nCP_simd = ROUNDSIZE(nCP);
    auto Xcp = static_cast<type*>(aligned_alloc(64,   3*nCP_simd*sizeof(type)));
    auto Ucp = static_cast<type*>(aligned_alloc(64,nDim*nCP_simd*sizeof(type)));

    if(readInput(nPt, nCP, nDim, Xcp, Ucp, X)!=0) CLEAN_AND_EXIT(-1)

    switch(nDim) {
      case 2: deform<2>(nPt_simd, nCP_simd, Xcp, Ucp, X); break;
      case 3: deform<3>(nPt_simd, nCP_simd, Xcp, Ucp, X); break;
      default: CLEAN_AND_EXIT(-1)
    }

    if(writeOutput(nPt,X)!=0) CLEAN_AND_EXIT(-1)

    CLEAN_AND_EXIT(0)
}
