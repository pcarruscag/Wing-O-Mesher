//  Copyright (C) 2018-2021  Pedro Gomes
//  See full notice in NOTICE.md

#ifndef DEFINITIONS_HPP
#define DEFINITIONS_HPP

typedef double type;

enum : size_t { SIMDLEN = 4 }; // for double and 256bit AVX

#define ROUNDSIZE(N) ((N+SIMDLEN-1)/SIMDLEN)*SIMDLEN

const char coordsFile[] = "__coords.dat";
const char ctrlptFile[] = "__ctrlpt.dat";
const char pointsFile[] = "__points.dat";

#endif // DEFINITIONS_HPP
