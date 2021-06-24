//  Copyright (C) 2018-2021  Pedro Gomes
//  See full notice in NOTICE.md

#include "functions.hpp"

#include <fstream>

using namespace std;

bool readIntoPaddedArray(const char fname[], size_t m, size_t n, type* x)
{
    fstream file;
    file.open(fname, ios::in | ios::binary);
    if(file.bad()) return true;

    size_t mAlign = ROUNDSIZE(m);
    for(size_t j=0; j<n; ++j)
    {
        auto xj = &x[j*mAlign];
        file.read(reinterpret_cast<char*>(xj), m*sizeof(type));
        // fill padding with zeros
        for(size_t i=mAlign; i<m; ++i) xj[i] = 0.0;
    }
    if(file.bad()) return true;

    return false;
}

int readInput(size_t nPt, size_t nCP, size_t nDim, type* Xcp, type* Ucp, type* X)
{
    bool error = false;

    #pragma omp parallel reduction(||:error)
    {
      #pragma omp single
      {
        // read control point coordinates
        #pragma omp task
        error = readIntoPaddedArray(coordsFile,nCP,3,Xcp);

        // read control point displacements
        #pragma omp task
        error = readIntoPaddedArray(ctrlptFile,nCP,nDim,Ucp);

        // read point coordinates
        #pragma omp task
        error = readIntoPaddedArray(pointsFile,nPt,3,X);
      }
    }

    return error? -1 : 0;
}

int writeOutput(size_t nPt, const type* X)
{
    fstream file;
    file.open(pointsFile, ios::out | ios::binary);
    if(file.bad()) return -1;

    size_t m = ROUNDSIZE(nPt);
    for(size_t k=0; k<3; ++k)
    {
        file.write(reinterpret_cast<const char*>(&X[k*m]), nPt*sizeof(type));
    }
    if(file.bad()) return -1;

    return 0;
}
