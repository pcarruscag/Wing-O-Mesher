//  Copyright (C) 2018-2021  Pedro Gomes
//  See full notice in NOTICE.md

#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP

#include <stdlib.h>
#include <vector>
#include "definitions.hpp"


/*!
 * Read input files (point coordinates and displacements)
 */
int readInput(size_t nPt,
              size_t nRhs,
              std::vector<type> & x,
              std::vector<type> & y,
              std::vector<type> & z,
              std::vector<type> & b);

/*!
 * Write output file (control point weights)
 */
int writeOutput(std::vector<type> const & sol);


/*!
 * Build Wendland C2 RBF kernel assuming reference radius = 1
 * Rectangular full packed storage (RFP), lower triangular, column-major.
 */
void buildKernel(std::vector<type> const & x,
                 std::vector<type> const & y,
                 std::vector<type> const & z,
                 std::vector<type> & M);

/*!
 * Compute the product (B) of a column-major matrix (A[m*n])
 * with the RBF kernel without holding the kernel in memory
 * (hard-coded type as avx intrinsics are used).
 */
void kernelProduct(size_t m,
                   size_t n,
                   double const * x,
                   double const * y,
                   double const * z,
                   double const * A,
                   double * B);

/*!
 * Solve RBF system by iterative refinement with
 * the LLT factorization using minimum memory.
 * On entry res is the rhs and on exit the residual vector.
 */
int solveRBF(size_t& iters,
             type& error,
             std::vector<type> const & x,
             std::vector<type> const & y,
             std::vector<type> const & z,
             std::vector<type> & res,
             std::vector<type> & lhs);


#endif // FUNCTIONS_HPP
