/**
 * Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

// A bunch of operation that benefits from OMP hyperthreading

#include "omp_la.h"

namespace omp_la {

/// @brief Reproduces 'SUBROUTINE OMPMULS (nNo, r, U)'.
//
void omp_mul_s(const int nNo, const double r, Vector<double>& U)
{
  for (int i = 0; i < nNo; i++) {
    U(i) = r * U(i);
  }
}

/// @brief Reproduces 'SUBROUTINE OMPMULV (dof, nNo, r, U)'.
//
void omp_mul_v(const int dof, const int nNo, const double r, Array<double>& U)
{
  switch (dof) {
    case 1:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = r * U(0,i);
      }
    break;

    case 2:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = r*U(0,i);
        U(1,i) = r*U(1,i);
      }
    break;

    case 3:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = r*U(0,i);
        U(1,i) = r*U(1,i);
        U(2,i) = r*U(2,i);
      }
    break;

    case 4:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = r*U(0,i);
        U(1,i) = r*U(1,i);
        U(2,i) = r*U(2,i);
        U(3,i) = r*U(3,i);
      }
    break;

    default: 
      for (int i = 0; i < nNo; i++) {
        for (int j = 0; j < U.nrows(); j++) {
          U(j,i) = r*U(j,i);
        }
        //U(:,i) = r*U(:,i)
      }
  }
}

/// @brief Reproduces 'SUBROUTINE OMPSUMS (nNo, r, U, V)'.
//
void omp_sum_s(const int nNo, const double r, Vector<double>& U, const Vector<double>& V)
{
  for (int i = 0; i < nNo; i++) {
    U(i) = U(i) + r*V(i);
  }
}

/// @brief Reproduces 'SUBROUTINE OMPSUMV (dof, nNo, r, U, V)'.
//
void omp_sum_v(const int dof, const int nNo, const double r, Array<double>& U, const Array<double>& V)
{
  switch (dof) {

    case 1:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = U(0,i) + r*V(0,i);
      }
    break;

    case 2:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = U(0,i) + r*V(0,i);
        U(1,i) = U(1,i) + r*V(1,i);
      }
    break;

    case 3:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = U(0,i) + r*V(0,i);
        U(1,i) = U(1,i) + r*V(1,i);
        U(2,i) = U(2,i) + r*V(2,i);
      }
    break;

    case 4:
      for (int i = 0; i < nNo; i++) {
        U(0,i) = U(0,i) + r*V(0,i);
        U(1,i) = U(1,i) + r*V(1,i);
        U(2,i) = U(2,i) + r*V(2,i);
        U(3,i) = U(3,i) + r*V(3,i);
      }
    break;

    default:
      for (int j = 0; j < U.nrows(); j++) {
        for (int i = 0; i < nNo; i++) {
          U(j,i) = U(j,i) + r*V(j,i);
        }
        //U(:,i) = U(:,i) + r*V(:,i)
      } 
  }
}


};


