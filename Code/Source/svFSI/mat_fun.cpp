
#include "mat_fun.h"

/*
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
*/

#include "consts.h"
#include "utils.h"

#include <math.h>

#include "lapack_defs.h"

namespace mat_fun {

Array<int> t_ind;
//static Array<int> t_ind;

//----------
// mat_ddot
//----------
// Double dot product of 2 square matrices
//
double mat_ddot(const Array<double>& A, const Array<double>& B, const int nd)
{
  double s = 0.0;

  for (int j = 0; j < nd; j++) { 
    for (int i = 0; i < nd; i++) { 
      s = s + A(i,j) * B(i,j);
    }
  }

  return s;
}

//---------
// mat_det
//---------
//
double 
mat_det(const Array<double>& A, const int nd)
{
  double D = 0.0;

  if (nd == 2) {
    D = A(0,0)*A(1,1) - A(0,1)*A(1,0);

  } else { 
    Array<double> Am(nd-1, nd-1);

    for (int i = 0; i < nd; i++) { 
      int n = 0;

      for (int j = 0; j < nd; j++) { 
        if (i == j) {
          continue; 
        } else { 

          for (int k = 0; k < nd-1; k++) { 
            Am(k,n) = A(k+1,j);
          } 

          n = n + 1;
        }
      }
      D = D + pow(-1.0, static_cast<double>(2+i)) * A(0,i) * mat_det(Am,nd-1);
    }
  }

  return D;
}

//---------
// mat_dev
//---------
//
Array<double> 
mat_dev(const Array<double>& A, const int nd)
{
  Array<double> result(nd,nd);

  double trA = mat_trace(A,nd);

  result = A - (trA / static_cast<double>(nd)) * mat_id(nd);

  return result;
}

//---------------
// mat_dyad_prod
//---------------
// Create a matrix from outer product of two vectors.
//
Array<double> 
mat_dyad_prod(const Vector<double>& u, const Vector<double>& v, const int nd)
{
  Array<double> result(nd,nd);

  for (int j = 0; j < nd; j++) {
    for (int i = 0; i < nd; i++) {
      result(i,j) = u(i) * v(j);
    }
  }

  return result;
}

//--------
// mat_id
//--------
//
Array<double> 
mat_id(const int nd)
{
  Array<double> A(nd,nd);

  for (int i = 0; i < nd; i++) {
    A(i,i) = 1.0;
  }

  return A;
}

//---------
// mat_inv
//---------
// This function computes inverse of a square matrix
//
Array<double> 
mat_inv(const Array<double>& A, const int nd, bool debug)
{
  int iok = 0;
  Array<double> Ainv(nd,nd);

  if (nd == 2) {
    double d = mat_det(A, nd);
    if (utils::is_zero(fabs(d))) {
      iok = -1;
    }

    // [NOTE] Divide by zero is possible here?
    Ainv(0,0) =  A(1,1) / d;
    Ainv(0,1) = -A(0,1) / d;

    Ainv(1,0) = -A(1,0) / d;
    Ainv(1,1) =  A(0,0) / d;

  } else if (nd == 3) {
    double d = mat_det(A, nd);
    if (utils::is_zero(fabs(d))) {
      iok = -1;
    }

    Ainv(0,0) = (A(1,1)*A(2,2) - A(1,2)*A(2,1)) / d;
    Ainv(0,1) = (A(0,2)*A(2,1) - A(0,1)*A(2,2)) / d;
    Ainv(0,2) = (A(0,1)*A(1,2) - A(0,2)*A(1,1)) / d;

    Ainv(1,0) = (A(1,2)*A(2,0) - A(1,0)*A(2,2)) / d;
    Ainv(1,1) = (A(0,0)*A(2,2) - A(0,2)*A(2,0)) / d;
    Ainv(1,2) = (A(0,2)*A(1,0) - A(0,0)*A(1,2)) / d;

    Ainv(2,0) = (A(1,0)*A(2,1) - A(1,1)*A(2,0)) / d;
    Ainv(2,1) = (A(0,1)*A(2,0) - A(0,0)*A(2,1)) / d;
    Ainv(2,2) = (A(0,0)*A(1,1) - A(0,1)*A(1,0)) / d;

  } else if ((nd > 3) && (nd < 10)) {
    double d = mat_det(A, nd);
    if (utils::is_zero(fabs(d))) {
      iok = -1;
    }
    Ainv = mat_inv_ge(A, nd, debug);

  } else {
    Ainv = mat_inv_lp(A, nd);
  } 

  if (iok != 0) {
     throw std::runtime_error("Singular matrix detected to compute inverse");
  }

  return Ainv;
}

//------------
// mat_inv_ge
//------------
// This function computes inverse of a square matrix using Gauss Elimination method
//
Array<double>
mat_inv_ge(const Array<double>& Ain, const int n, bool debug)
{
  Array<double> A(n,n);
  Array<double> B(n,n);
  A = Ain;

  if (debug) {
    std::cout << "[mat_inv_ge] ========== mat_inv_ge =========" << std::endl;
    if (std::numeric_limits<double>::is_iec559) {
      std::cout << "[mat_inv_ge] is_iec559 " << std::endl;
    }
  }

  // Auxillary matrix
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      B(i,j) = 0.0;
    }
    B(i,i) = 1.0;
  }

  // Loop over columns of A to find a good pivot. 
  //
  double max_val = 0.0; 
  int irow = 0;
  double d = 0.0; 

  if (debug) {
    A.print("A");
  }

  for (int i = 0; i < n; i++) {
    if (debug) {
      std::cout << "[mat_inv_ge] " << std::endl;
      std::cout << "[mat_inv_ge] ---------- i: " << i+1 << std::endl;
      A.print("A");
    }
    double max_val = fabs(A(i,i));
    irow = i;

    for (int j = i; j < n; j++) {
      //if (debug) {
       // std::cout << "[mat_inv_ge] A(j,i): " << A(j,i) << std::endl;
        //std::cout << "[mat_inv_ge] max_val: " << max_val << std::endl;
      //}
      if (fabs(A(j,i)) > max_val) {
        max_val = fabs(A(j,i));
        irow = j;
      }
    }

    if (debug) {
      std::cout << "[mat_inv_ge] max_val: " << max_val << std::endl;
      std::cout << "[mat_inv_ge] irow: " << irow+1 << std::endl;
    }

    // Interchange rows.
    //
    if (max_val > fabs(A(i,i))) {
      if (debug) {
        std::cout << "[mat_inv_ge] " << std::endl;
        std::cout << "[mat_inv_ge] Interchange rows " << std::endl;
      }

      for (int k = 0; k < n; k++) {
        d = A(i,k);
        A(i,k) = A(irow,k);
        A(irow,k) = d;

        d = B(i,k) ;
        B(i,k) = B(irow,k);
        B(irow,k) = d;
      }
    }

    d = A(i,i);

    if (debug) {
      std::cout << "[mat_inv_ge]  " << std::endl;
      std::cout << "[mat_inv_ge]  Scale ..." << std::endl;
      std::cout << "[mat_inv_ge] d: " << d << std::endl;
    }

    for (int j = 0; j < n; j++) {
      A(i,j) = A(i,j) / d;
      B(i,j) = B(i,j) / d;
    }

    if (debug) {
      std::cout << "[mat_inv_ge]  " << std::endl;
      std::cout << "[mat_inv_ge]  Reduce ..." << std::endl;
    }

    for (int j = i+1; j < n; j++) {
      d = A(j,i);

      if (debug) {
        std::cout << "[mat_inv_ge]  " << std::endl;
        std::cout << "[mat_inv_ge]  j: " << j+1 << std::endl;
        std::cout << "[mat_inv_ge]  d: " << d << std::endl;
      }

      for (int k = 0; k < n; k++) {
        A(j,k) = A(j,k) - d*A(i,k);
        B(j,k) = B(j,k) - d*B(i,k);
      }
    }
  }

  if (debug) {
    std::cout << "[mat_inv_ge]  " << std::endl;
    std::cout << "[mat_inv_ge] Final reduce ..." << std::endl;
  }

  for (int i = 0; i < n-1; i++) {
    for (int j = i+1; j < n; j++) {
      d = A(i,j);
      if (debug) {
        std::cout << "[mat_inv_ge] i j " << i+1 << " " << j+1 << std::endl;
        std::cout << "[mat_inv_ge] d: " << d << std::endl;
      }
      for (int k = 0; k < n; k++) {
        A(i,k) = A(i,k) - d*A(j,k);
        B(i,k) = B(i,k) - d*B(j,k);
        if (debug) {
          std::cout << "[mat_inv_ge] B(i,k): " << B(i,k) << std::endl;
        }
      }
    }
  }

  return B;
}

//------------
// mat_inv_ge
//------------
// This function computes inverse of a square matrix using Gauss Elimination method
//
// [TODO:DaveP] The original version sometimes produced NaNs.
//
Array<double> 
mat_inv_ge_orig(const Array<double>& A, const int nd, bool debug)
{
  Array<double> B(nd,2*nd);
  Array<double> Ainv(nd,nd);

  // Auxillary matrix
  for (int i = 0; i < nd; i++) { 
    for (int j = 0; j < nd; j++) { 
      B(i,j) = A(i,j);
    }
    B(i,nd+i) = 1.0;
  }

  // Pivoting
  for (int i = nd-1; i > 0; i--) { 
    if (B(i,0) > B(i-1,0)) {
      for (int j = 0; j < 2*nd; j++) { 
        double d = B(i,j);
        B(i,j) = B(i-1,j);
        B(i-1,j) = d;
      }
    }
  }

  if (debug) {
    std::cout << "[mat_inv]  " << std::endl;
    std::cout << "[mat_inv] B: " << B << std::endl;
    std::cout << "[mat_inv]  " << std::endl;
  }

  // Do row-column operations and reduce to diagonal
  double d;

  for (int i = 0; i < nd; i++) { 
    if (debug) {
      std::cout << "[mat_inv] ========== i " << i+1 << " ==========" << std::endl;
    }
    for (int j = 0; j < nd; j++) { 
      if (debug) {
        std::cout << "[mat_inv] ########## j " << j+1 << " ##########" << std::endl;
      }
      if (j != i) {
        d = B(j,i) / B(i,i);
        if (debug) {
          std::cout << "[mat_inv] B(j,i): " << B(j,i) << std::endl;
          std::cout << "[mat_inv] B(i,i): " << B(i,i) << std::endl;
          std::cout << "[mat_inv] d: " << d << std::endl;
          std::cout << "[mat_inv]  " << std::endl;
        }
        for (int k = 0; k < 2*nd; k++) { 
          if (debug) {
            std::cout << "[mat_inv] ------- k " << k+1 << " -------" << std::endl;
          }
          if (debug) {
            std::cout << "[mat_inv] B(j,k): " << B(j,k) << std::endl;
            std::cout << "[mat_inv] B(i,k): " << B(i,k) << std::endl;
            std::cout << "[mat_inv] d: " << d << std::endl;
            std::cout << "[mat_inv] d*B(i,k): " << d*B(i,k) << std::endl;
          }
          B(j,k) = B(j,k) - d*B(i,k);
          if (debug) {
            std::cout << "[mat_inv] B(j,k): " << B(j,k) << std::endl;
          }
        }
      }
    }
  }

  // Unit matrix
  for (int i = 0; i < nd; i++) { 
    double d = B(i,i);
    for (int j = 0; j < 2*nd; j++) { 
      B(i,j) = B(i,j) / d;
    }
  }

  // Inverse
  for (int i = 0; i < nd; i++) { 
    for (int j = 0; j < nd; j++) { 
      Ainv(i,j) = B(i,j+nd);
    }
  }

  return Ainv;
}

//------------
// mat_inv_lp
//------------
// This function computes inverse of a square matrix using Lapack functions (DGETRF + DGETRI)
//
// Replaces 'FUNCTION MAT_INV_LP(A, nd) RESULT(Ainv)' defined in MATFUN.f.
//
Array<double> 
mat_inv_lp(const Array<double>& A, const int nd)
{
  Vector<int> IPIV(nd);
  int iok;
  int n = nd;
  auto Ad = A;

  dgetrf_(&n, &n, Ad.data(), &n, IPIV.data(), &iok);

  if (iok != 0) {
    throw std::runtime_error("Singular matrix detected to compute inverse");
  }

  Vector<double> WORK(2*nd);
  int nd_2 = 2*nd;

  dgetri_(&n, Ad.data(), &n, IPIV.data(), WORK.data(), &nd_2, &iok);

  if (iok != 0) {
    throw std::runtime_error("ERROR: Matrix inversion failed (LAPACK)");
  }

  Array<double> Ainv(nd,nd);

  for (int i = 0; i < nd; i++) {
    for (int j = 0; j < nd; j++) {
      Ainv(i,j) = Ad(i,j);
    }
  }

  return Ainv;
}

//------------------
// mat_inv_lp_eigen
//------------------
// not used, just a test.
//
Array<double> 
mat_inv_lp_eigen(const Array<double>& A, const int nd)
{
#ifdef use_eigen
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixType;

  MatrixType Ad(nd, nd);

  for (int i = 0; i < nd; i++) { 
    for (int j = 0; j < nd; j++) { 
      Ad(i,j) = A(i,j);
    }
  }

  Eigen::FullPivLU<MatrixType> lu(Ad);

  lu.inverse();
  Array<double> Ainv(nd,nd);

  for (int i = 0; i < nd; i++) { 
    for (int j = 0; j < nd; j++) { 
      //Ainv(i,j) = Ad_inverse(i,j);
    }
  }

  return Ainv;

#endif
}

//---------
// mat_mul
//---------
// Multiply a matrix by a vector.
//
// Reproduces Fortran MATMUL.
//
Vector<double> 
mat_mul(const Array<double>& A, const Vector<double>& v)
{
  int num_rows = A.nrows();
  int num_cols = A.ncols();

  if (num_cols != v.size()) {
    throw std::runtime_error("[mat_mul] The number of columns of A (" + std::to_string(num_cols) + ") does not equal the size of v (" + 
        std::to_string(v.size()) + ").");
  }

  Vector<double> result(num_rows);

  for (int i = 0; i < num_rows; i++) {
    double sum = 0.0;

    for (int j = 0; j < num_cols; j++) {
      sum += A(i,j) * v(j);
    }

    result(i) = sum;
  }

  return result;
}

//---------
// mat_mul
//---------
// Multiply a matrix by a matrix.
//
// Reproduces Fortran MATMUL.
//
Array<double> 
mat_mul(const Array<double>& A, const Array<double>& B)
{
  int A_num_rows = A.nrows();
  int A_num_cols = A.ncols();
  int B_num_rows = B.nrows();
  int B_num_cols = B.ncols();

  if (A_num_cols != B_num_rows) {
    throw std::runtime_error("[mat_mul] The number of columns of A (" + std::to_string(A_num_cols) + ") does not equal " +
        " the number of rows of B (" + std::to_string(B_num_rows) + ").");
  }

  Array<double> result(A_num_rows, B_num_cols);

  for (int i = 0; i < A_num_rows; i++) {
    for (int j = 0; j < B_num_cols; j++) {
      double sum = 0.0;
      for (int k = 0; k < A_num_cols; k++) {
        sum += A(i,k) * B(k,j);
      }
    result(i,j) = sum;
    }
  }

  return result;
}

//---------
// mat_mul
//---------
// Multiply a matrix by a matrix.
//
// Compute result directly into the passed argument.
//
void mat_mul(const Array<double>& A, const Array<double>& B, Array<double>& result)
{ 
  int A_num_rows = A.nrows();
  int A_num_cols = A.ncols();
  int B_num_rows = B.nrows();
  int B_num_cols = B.ncols();
  
  if (A_num_cols != B_num_rows) {
    throw std::runtime_error("[mat_mul] The number of columns of A (" + std::to_string(A_num_cols) + ") does not equal " +
        " the number of rows of B (" + std::to_string(B_num_rows) + ").");
  }
  
  for (int i = 0; i < A_num_rows; i++) { 
    for (int j = 0; j < B_num_cols; j++) {
      double sum = 0.0; 

      for (int k = 0; k < A_num_cols; k++) {
        sum += A(i,k) * B(k,j);
      }

      result(i,j) = sum;
    }
  }
}

//---------
// mat_sym
//---------
// Symmetric part of a matrix, S = (A + A.T)/2
//
Array<double> 
mat_symm(const Array<double>& A, const int nd)
{
  Array<double> S(nd, nd);

  for (int i = 0; i < nd; i++) { 
    for (int j = 0; j < nd; j++) { 
      S(i,j) = 0.5* (A(i,j) + A(j,i));
    }
  }

  return S;
}


//---------------
// mat_symm_prod
//---------------
// Create a matrix from symmetric product of two vectors
//
Array<double> 
mat_symm_prod(const Vector<double>& u, const Vector<double>& v, const int nd)
{
  Array<double> result(nd, nd);

  for (int i = 0; i < nd; i++) { 
    for (int j = 0; j < nd; j++) { 
      result(i,j) = 0.5 * (u(i)*v(j) + u(j)*v(i));
    }
  }

  return result;
}

//-----------
// mat_trace
//-----------
// Trace of second order matrix of rank nd
//
double mat_trace(const Array<double>& A, const int nd)
{
  double result = 0.0;

  for (int i = 0; i < nd; i++) {
    result += A(i,i);
  }

  return result;
}

//-----------------
// ten_asym_prod12
//-----------------
// Create a 4th order tensor from antisymmetric outer product of
// two matrices
//
//   Cijkl = Aij*Bkl-Ail*Bjk
//
Tensor4<double> 
ten_asym_prod12(const Array<double>& A, const Array<double>& B, const int nd)
{
  Tensor4<double> C(nd,nd,nd,nd);

  int nn = pow(nd,4);

  for (int ii = 0; ii < nn; ii) { 
    int i = t_ind(0,ii);
    int j = t_ind(1,ii);
    int k = t_ind(2,ii);
    int l = t_ind(3,ii);
    C(i,j,k,l) = 0.5 * ( A(i,j)*B(k,l) - A(i,l)*B(j,k) );
  }

  return C;
}

//----------
// ten_ddot
//----------
//Double dot product of 2 4th order tensors T_ijkl = A_ijmn * B_klmn
//
// Reproduces 'FUNCTION TEN_DDOT_3434(A, B, nd) RESULT(C)'.
//
Tensor4<double>
ten_ddot(const Tensor4<double>& A, const Tensor4<double>& B, const int nd)
{
  int nn = pow(nd,4);
  Tensor4<double> C(nd,nd,nd,nd);

  if (nd == 2) {
    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C(i,j,k,l) = C(i,j,k,l) + A(i,j,0,0)*B(k,l,0,0) + 
                                A(i,j,0,1)*B(k,l,0,1) + 
                                A(i,j,1,0)*B(k,l,1,0) +                      
                                A(i,j,1,1)*B(k,l,1,1);
    }

  } else {

    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C(i,j,k,l) = C(i,j,k,l) + A(i,j,0,0)*B(k,l,0,0)
                              + A(i,j,0,1)*B(k,l,0,1)
                              + A(i,j,0,2)*B(k,l,0,2)
                              + A(i,j,1,0)*B(k,l,1,0)
                              + A(i,j,1,1)*B(k,l,1,1)
                              + A(i,j,1,2)*B(k,l,1,2)
                              + A(i,j,2,0)*B(k,l,2,0)
                              + A(i,j,2,1)*B(k,l,2,1)
                              + A(i,j,2,2)*B(k,l,2,2);

    }
  }

  return C;
}

//---------------
// ten_ddot_2412
//---------------
// T_ijkl = A_imjn * B_mnkl
//
Tensor4<double>
ten_ddot_2412(const Tensor4<double>& A, const Tensor4<double>& B, const int nd)
{
  int nn = pow(nd,4);
  Tensor4<double> C(nd,nd,nd,nd);

  if (nd == 2) {
    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C(i,j,k,l) = C(i,j,k,l) + A(i,0,j,0)*B(0,0,k,l)
                              + A(i,0,j,1)*B(0,1,k,l)
                              + A(i,1,j,0)*B(1,0,k,l)
                              + A(i,1,j,1)*B(1,1,k,l);

    }

  } else {

    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C(i,j,k,l) = C(i,j,k,l) + A(i,0,j,0)*B(0,0,k,l)
                              + A(i,0,j,1)*B(0,1,k,l)
                              + A(i,0,j,2)*B(0,2,k,l)
                              + A(i,1,j,0)*B(1,0,k,l)
                              + A(i,1,j,1)*B(1,1,k,l)
                              + A(i,1,j,2)*B(1,2,k,l)
                              + A(i,2,j,0)*B(2,0,k,l)
                              + A(i,2,j,1)*B(2,1,k,l)
                              + A(i,2,j,2)*B(2,2,k,l);
    }
  }

  return C;
}


//---------------
// ten_ddot_3424
//---------------
//
Tensor4<double>
ten_ddot_3424(const Tensor4<double>& A, const Tensor4<double>& B, const int nd)
{
  int nn = pow(nd,4);
  Tensor4<double> C(nd,nd,nd,nd);

  if (nd == 2) {
    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C(i,j,k,l) = C(i,j,k,l) + A(i,j,0,0)*B(k,0,l,0) +
                                A(i,j,0,1)*B(k,0,l,1) +
                                A(i,j,1,0)*B(k,1,l,0) +
                                A(i,j,1,1)*B(k,1,l,1);
    }

  } else {

    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C(i,j,k,l) = C(i,j,k,l) + A(i,j,0,0)*B(k,0,l,0)
                              + A(i,j,0,1)*B(k,0,l,1)
                              + A(i,j,0,2)*B(k,0,l,2)
                              + A(i,j,1,0)*B(k,1,l,0)
                              + A(i,j,1,1)*B(k,1,l,1)
                              + A(i,j,1,2)*B(k,1,l,2)
                              + A(i,j,2,0)*B(k,2,l,0)
                              + A(i,j,2,1)*B(k,2,l,1)
                              + A(i,j,2,2)*B(k,2,l,2);

    }
  }

  return C;
}


//----------
// ten_init
//----------
// Initialize tensor index pointer
//
void ten_init(const int nd)
{
  int nn = pow(nd, 4);
  t_ind.resize(4, nn);

  int ii = 0;
  for (int l = 0; l < nd; l++) {
    for (int k = 0; k < nd; k++) {
      for (int j = 0; j < nd; j++) {
        for (int i = 0; i < nd; i++) {
          t_ind(0,ii) = i;
          t_ind(1,ii) = j;
          t_ind(2,ii) = k;
          t_ind(3,ii) = l;
          ii = ii + 1;
        }
      }
    }
  }
}

//---------------
// ten_dyad_prod
//---------------
// Create a 4th order tensor from outer product of two matrices.
//
Tensor4<double>  
ten_dyad_prod(const Array<double>& A, const Array<double>& B, const int nd)
{   
  int nn = pow(nd,4);
  Tensor4<double> C(nd,nd,nd,nd);
  
  for (int ii = 0; ii < nn; ii++) {
    int i = t_ind(0,ii);
    int j = t_ind(1,ii);
    int k = t_ind(2,ii);
    int l = t_ind(3,ii);
    C(i,j,k,l) = A(i,j) * B(k,l);
  }

  return C;
}

//---------
// ten_ids
//---------
// Create a 4th order order symmetric identity tensor
//
Tensor4<double>
ten_ids(const int nd)
{
  Tensor4<double> A(nd,nd,nd,nd);

  for (int i = 0; i < nd; i++) {
    for (int j = 0; j < nd; j++) {
      A(i,j,i,j) = A(i,j,i,j) + 0.5;
      A(i,j,j,i) = A(i,j,j,i) + 0.5;
    }      
  }      

  return A;
}

//-----------
// ten_mddot
//-----------
// Double dot product of a 4th order tensor and a 2nd order tensor
//
//   C_ij = (A_ijkl * B_kl)
//
Array<double> 
ten_mddot(const Tensor4<double>& A, const Array<double>& B, const int nd) 
{
  Array<double> C(nd,nd);

  if (nd == 2) {
    for (int i = 0; i < nd; i++) {
      for (int j = 0; j < nd; j++) {
        C(i,j) = A(i,j,0,0)*B(0,0) + A(i,j,0,1)*B(0,1) + A(i,j,1,0)*B(1,0) + A(i,j,1,1)*B(1,1);
      }
    }

  } else { 
    for (int i = 0; i < nd; i++) {
      for (int j = 0; j < nd; j++) {
        C(i,j) = A(i,j,0,0)*B(0,0) + A(i,j,0,1)*B(0,1) + A(i,j,0,2)*B(0,2) + A(i,j,1,0)*B(1,0) + 
            A(i,j,1,1)*B(1,1) + A(i,j,1,2)*B(1,2) + A(i,j,2,0)*B(2,0) + A(i,j,2,1)*B(2,1) + 
            A(i,j,2,2)*B(2,2);
        }
     }
  }

  return C;
}


//---------------
// ten_symm_prod
//---------------
// Create a 4th order tensor from symmetric outer product of two matrices.
//
// Reproduces 'FUNCTION TEN_SYMMPROD(A, B, nd) RESULT(C)'.
//
Tensor4<double> 
ten_symm_prod(const Array<double>& A, const Array<double>& B, const int nd)
{
  int nn = pow(nd,4);
  Tensor4<double> C(nd,nd,nd,nd);
  
  for (int ii = 0; ii < nn; ii++) {
    int i = t_ind(0,ii);
    int j = t_ind(1,ii);
    int k = t_ind(2,ii);
    int l = t_ind(3,ii);
    C(i,j,k,l) = 0.5* ( A(i,k)*B(j,l) + A(i,l)*B(j,k) );
  } 

  return C;
}

//---------------
// ten_transpose
//---------------
//
Tensor4<double>
ten_transpose(const Tensor4<double>& A, const int nd)
{ 
  int nn = pow(nd,4);
  Tensor4<double> result(nd,nd,nd,nd);
  
  for (int ii = 0; ii < nn; ii++) {
    int i = t_ind(0,ii);
    int j = t_ind(1,ii);
    int k = t_ind(2,ii);
    int l = t_ind(3,ii);
    result(i,j,k,l) = A(k,l,i,j);
  }
  
  return result;
}

//-----------
// transpose
//-----------
// Reproduces Fortran TRANSPOSE.
//
Array<double> 
transpose(const Array<double>& A)
{
  int num_rows = A.nrows();
  int num_cols = A.ncols();
  Array<double> result(num_cols, num_rows);

  for (int i = 0; i < num_rows; i++) {
    for (int j = 0; j < num_cols; j++) {
      result(j,i) = A(i,j);
    }
  }

  return result;
}

};


