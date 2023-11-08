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

// The dot product between two scaler-vector container vectors are
// calculated here.
//
// Only the part of U and V which are owned by this processor is
// included in dot product calculation
// In order to have the correct answer it is needed that COMMU has
// been done before calling this function (or the ansesters of U and
// V are passed through COMMU)

#include "dot.h"

#include "fils_struct.hpp"

namespace dot {

/// @brief Reproduces 'FUNCTION FSILS_DOTS(nNo, commu, U, V)'. 
//
double fsils_dot_s(const int nNo, FSILS_commuType& commu, const Vector<double>& U, const Vector<double>& V)
{
  double result = 0.0; 

  for (int i = 0; i < nNo; i++) {
    result = result + U(i)*V(i);
  }

  if (commu.nTasks == 1) {
    return result;
  }

  double tmp{0.0};

  MPI_Allreduce(&result, &tmp, 1, cm_mod::mpreal, MPI_SUM, commu.comm);

  return tmp;
}

/// @brief Reproduces 'FUNCTION FSILS_DOTV(dof, nNo, commu, U, V)'.
//
double fsils_dot_v(const int dof, const int nNo, FSILS_commuType& commu, const Array<double>& U, const Array<double>& V)
{
  double result = 0.0; 

  switch (dof) {
    case 1:
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*V(0,i);
      }
    break;

    case 2:
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i);
      }
    break;

    case 3:
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i) +  U(2,i)*V(2,i);
      }
    break;

    case 4:
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i) + U(2,i)*V(2,i) + U(3,i)*V(3,i);
      }
    break;

    defualt: 
      for (int i = 0; i < nNo; i++) {
        double sum{0.0};
        for (int j = 0; j < U.nrows(); j++) {
          sum += U(j,i) * V(j,i);
        }
        result = result + sum; 
      }
  }

  if (commu.nTasks == 1) {
    return result;
  }

  double tmp{0.0};

  MPI_Allreduce(&result, &tmp, 1, cm_mod::mpreal, MPI_SUM, commu.comm);

  return tmp;
}

/// @brief Reproduces Fortran 'FSILS_NCDOTS(nNo, , U, V)'.
//
double fsils_nc_dot_s(const int nNo, const Vector<double>& U, const Vector<double>& V)
{
  double result{0.0};

  for (int i = 0; i < nNo; i++) { 
    result = result + U(i)*V(i);
  } 

  return result;
}

/// @brief Reproduces 'FUNCTION FSILS_NCDOTV(dof, nNo, U, V) RESULT(FSILS_DOTV)'.
//
double fsils_nc_dot_v(const int dof, const int nNo, const Array<double>& U, const Array<double>& V)
{
  double result = 0.0;

  switch (dof) {
    case 1: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U(0,i)*V(0,i);
      }
    } break;

    case 2: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i);
      }
    } break;

    case 3: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i) + U(2,i)*V(2,i);
      }
    } break;

    case 4: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U(0,i)*V(0,i) + U(1,i)*V(1,i) +  U(2,i)*V(2,i) + U(3,i)*V(3,i);
      }
    } break;

    default: {
      for (int i = 0; i < nNo; i++) { 
        result = result + U.col(i) * V.col(i);
      }
    } break;
  } 

  return result;
}

};


