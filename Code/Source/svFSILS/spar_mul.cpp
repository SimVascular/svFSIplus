/* Copyright (c) Stanford University, The Regents of the University of California, and others.
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

//--------------------------------------------------------------------
// Product of a sparse matrix and a vector. The matrix might be
// vector in neither, one or both dimensions.
//--------------------------------------------------------------------
//
// Reproduces code in SPARMUL.f.

#include "spar_mul.h"

#include "fsils_api.hpp"

namespace spar_mul {

/// @brief Reproduces 'SUBROUTINE FSILS_SPARMULSS(lhs, rowPtr, colPtr, K, U, KU)'
//
void fsils_spar_mul_ss(FSILS_lhsType& lhs, const Array<int>& rowPtr, const Vector<int>& colPtr, 
    const Vector<double>& K, const Vector<double>& U, Vector<double>& KU)
{
  int nNo = lhs.nNo;
  KU = 0.0;

  for (int i = 0; i < nNo; i++) {
    for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
      KU(i) = KU(i) + K(j) * U(colPtr(j));
    } 
  }

  fsils_commus(lhs, KU);
}

/// @brief Reproduces 'SUBROUTINE FSILS_SPARMULSV(lhs, rowPtr, colPtr, dof, K, U, KU)'. 
//
void fsils_spar_mul_sv(FSILS_lhsType& lhs, const Array<int>& rowPtr, const Vector<int>& colPtr, 
    const int dof, const Array<double>& K, const Vector<double>& U, Array<double>& KU)
{
  int nNo = lhs.nNo;
  KU = 0.0;

  switch (dof) {

    case 1:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          KU(0,i) = KU(0,i) + K(0,j)*U(colPtr(j));
        }
      }
    break; 

    case 2: {
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU(0,i) = KU(0,i) + K(0,j)*U(col);
          KU(1,i) = KU(1,i) + K(1,j)*U(col);
        }
      }

    } break; 

    case 3: {
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU(0,i) += K(0,j) * U(col);
          KU(1,i) += K(1,j) * U(col);
          KU(2,i) += K(2,j) * U(col);
        }
      }
    } break; 

    case 4:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU(0,i) = KU(0,i) + K(0,j)*U(col);
          KU(1,i) = KU(1,i) + K(1,j)*U(col);
          KU(2,i) = KU(2,i) + K(2,j)*U(col);
          KU(3,i) = KU(3,i) + K(3,j)*U(col);
        }
      }
    break; 

    default: 
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          for (int m = 0; m < KU.nrows(); m++) {
            KU(m,i) = KU(m,i) + K(m,j) * U(col);
          }
        }
      }
  } 

  fsils_commuv(lhs, dof, KU);
}

/// @brief Reproduces 'SUBROUTINE FSILS_SPARMULVS(lhs, rowPtr, colPtr, dof, K, U, KU)'.
//
void fsils_spar_mul_vs(FSILS_lhsType& lhs, const Array<int>& rowPtr, const Vector<int>& colPtr, 
    const int dof, const Array<double>& K, const Array<double>& U, Vector<double>& KU)
{
  int nNo = lhs.nNo;
  KU = 0.0;

  switch (dof) {

    case 1:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          KU(i) = KU(i) + K(0,j) * U(0,colPtr(j));
        }
      }
    break; 

    case 2:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU(i) = KU(i) + K(0,j)*U(0,col) + K(1,j)*U(1,col);
        }
      }
    break; 

    case 3:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU(i) = KU(i) + K(0,j)*U(0,col) + K(1,j)*U(1,col) + K(2,j)*U(2,col);
        }
      }
    break; 

    case 4:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU(i) = KU(i) + K(0,j)*U(0,col) + K(1,j)*U(1,col) + K(2,j)*U(2,col) + K(3,j)*U(3,col);
        }
      }
    break; 

    default: 
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          double sum = 0.0;
          for (int m = 0; m < K.nrows(); m++) {
            sum += K(m,j) * U(m,col);
          }
          KU(i) = KU(i) + sum; 
          //KU(i) = KU(i) + SUM(K(:,j)*U(:,colPtr(j)))
        }
     }
  } 

  fsils_commus(lhs, KU);
}

/// @brief Reproduces 'SUBROUTINE FSILS_SPARMULVV(lhs, rowPtr, colPtr, dof, K, U, KU)'. 
//
void fsils_spar_mul_vv(FSILS_lhsType& lhs, const Array<int>& rowPtr, const Vector<int>& colPtr, 
    const int dof, const Array<double>& K, const Array<double>& U, Array<double>& KU)
{
  int nNo = lhs.nNo;
  KU = 0.0;

  switch (dof) {

    case 1:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          KU(0,i) = KU(0,i) + K(0,j)*U(0,colPtr(j));
        }
      }
    break;

    case 2:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU(0,i) = KU(0,i) + K(0,j)*U(0,col) + K(1,j)*U(1,col);
          KU(1,i) = KU(1,i) + K(2,j)*U(0,col) + K(3,j)*U(1,col);
        }
      }
    break;

    case 3:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU(0,i) = KU(0,i) + K(0,j)*U(0,col) + K(1,j)*U(1,col) + K(2,j)*U(2,col);
          KU(1,i) = KU(1,i) + K(3,j)*U(0,col) + K(4,j)*U(1,col) + K(5,j)*U(2,col);
          KU(2,i) = KU(2,i) + K(6,j)*U(0,col) + K(7,j)*U(1,col) + K(8,j)*U(2,col);
        }
      }
    break;

    case 4:
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU(0,i) = KU(0,i) + K(0 ,j)*U(0,col) + K(1 ,j)*U(1,col) + K(2 ,j)*U(2,col) + K(3 ,j)*U(3,col);
          KU(1,i) = KU(1,i) + K(4 ,j)*U(0,col) + K(5 ,j)*U(1,col) + K(6 ,j)*U(2,col) + K(7 ,j)*U(3,col);
          KU(2,i) = KU(2,i) + K(8 ,j)*U(0,col) + K(9,j)*U(1,col) + K(10,j)*U(2,col) + K(11,j)*U(3,col);
          KU(3,i) = KU(3,i) + K(12,j)*U(0,col) + K(13,j)*U(1,col) + K(14,j)*U(2,col) + K(15,j)*U(3,col);
        }
      }
    break;

    default: 
      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          for (int l = 0; l < dof; l++) {
            int e = l*dof;
            int s = e - dof + 1;;
            double sum = 0.0;
            for (int k = 0; k < dof; k++) {
              sum += K(k+s,j) * U(k,col);
            }
            KU(l,i) = KU(l,i) + sum;
          }
        }
     }
  } 

  fsils_commuv(lhs, dof, KU);
}

};


