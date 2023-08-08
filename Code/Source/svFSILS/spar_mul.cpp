
//--------------------------------------------------------------------
// Product of a sparse matrix and a vector. The matrix might be
// vector in neither, one or both dimensions.
//--------------------------------------------------------------------
//
// Reproduces code in SPARMUL.f.

#include "spar_mul.h"

#include "fsils_api.hpp"

namespace spar_mul {

//-------------------
// fsils_spar_mul_ss
//-------------------
// Reproduces 'SUBROUTINE FSILS_SPARMULSS(lhs, rowPtr, colPtr, K, U, KU)'
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

//-------------------
// fsils_spar_mul_sv
//-------------------
// Reproduces 'SUBROUTINE FSILS_SPARMULSV(lhs, rowPtr, colPtr, dof, K, U, KU)'. 
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
      auto KU_data = KU.data();
      auto nr_ku = KU.nrows();
      auto nc_ku = KU.ncols();
      //data_[row + col*nr_ku]

      auto K_data = K.data();
      auto nr_k = K.nrows();
      auto nc_k = K.ncols();
      //K_data[row + col*nr_k]

      auto U_data = U.data();

      for (int i = 0; i < nNo; i++) {
        for (int j = rowPtr(0,i); j <= rowPtr(1,i); j++) {
          int col = colPtr(j);
          KU_data[0] += K_data[0+j*nr_k] * U_data[col];
          KU_data[1] += K_data[1+j*nr_k] * U_data[col];
          KU_data[2] += K_data[2+j*nr_k] * U_data[col];
          //KU(0,i) += K(0,j) * U_data[col];
          //KU(1,i) += K(1,j) * U_data[col];
          //KU(2,i) += K(2,j) * U_data[col];
        }
        KU_data += nr_ku;
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

//-------------------
// fsils_spar_mul_vs
//-------------------
// Reproduces 'SUBROUTINE FSILS_SPARMULVS(lhs, rowPtr, colPtr, dof, K, U, KU)'.
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

//--------------------
// fsi_ls_spar_mul_vv
//--------------------
// Reproduces 'SUBROUTINE FSILS_SPARMULVV(lhs, rowPtr, colPtr, dof, K, U, KU)'. 
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


