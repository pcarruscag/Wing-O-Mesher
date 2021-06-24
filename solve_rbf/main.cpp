//  Copyright (C) 2018-2021  Pedro Gomes
//  See full notice in NOTICE.md

#include <cstdlib>

#include "definitions.hpp"
#include "functions.hpp"

#include <vector>
#include <omp.h>


int main(int argc, char* argv[])
{
    if(argc<4) return -1;

    const size_t nPt = atoi(argv[1]);
    const size_t nRhs = atoi(argv[2]);
    const size_t nThrd = atoi(argv[3]);

    omp_set_num_threads(nThrd);

    size_t iters = 10;
    type   error = 1e-21;
    std::vector<type> x, y, z, B, X;

    if(readInput(nPt,nRhs,x,y,z,B)!=0) return -1;

    if(solveRBF(iters,error,x,y,z,B,X)!=0) return -1;

    if(writeOutput(X)!=0) return -1;

    return 0;
}
