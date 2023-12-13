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

// These are template functions reproducing some mat_fun functions 
// using C++ arrays rather than Array objects. 

#ifndef MAT_FUN_CARRAY_H 
#define MAT_FUN_CARRAY_H 

#include "mat_fun.h"
#include "utils.h"

namespace mat_fun_carray {

extern Array<int> t_ind;

void ten_init(const int nd);

template <size_t N>
void mat_zero(double A[N][N])
{
  memset(A, 0, N*sizeof(A[0]));
}

template <size_t N>
void mat_zero(double A[N][N][N])
{
  memset(A, 0, N*sizeof(A[0]));
}

/// @brief Print a 2D array. 
//
template <size_t N>
void print(const std::string& name, const double A[N][N])
{
  std::cout << name << ": ";
  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      std::cout << A[i][j] << " ";
    }
  }
  std::cout << std::endl;
}

/// @brief Print a 3D array. 
//
template <size_t N>
void print(const std::string& name, const double A[N][N][N][N])
{
  std::cout << name << ": ";
  for (int l = 0; l < N; l++) {
    for (int k = 0; k < N; k++) {
      for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
          std::cout << A[i][j][k][l] << " ";
        }
      }
    }
  }
  std::cout << std::endl;
}


template <size_t N>
double mat_det(const double A[N][N])
{
  double D = 0.0;

  if (N == 2) {
    D = A[0][0]*A[1][1] - A[0][1]*A[1][0];

  } else {
    double Am[2][2];

    for (int i = 0; i < N; i++) {
      int n = 0;

      for (int j = 0; j < N; j++) {
        if (i == j) {
          continue;
        } else {

          for (int k = 0; k < N-1; k++) {
            Am[k][n] = A[k+1][j];
          }

          n = n + 1;
        }
      }
      D = D + pow(-1.0, static_cast<double>(2+i)) * A[0][i] * mat_det<2>(Am);
    }
  }

  return D;
  }

template <size_t N>
void mat_id(double A[N][N])
{
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      if (i == j) {
        A[i][j] = 1.0;
      } else {
        A[i][j] = 0.0;
      }
    }
  }
}

template <size_t N>
void transpose(const double A[N][N],  double result[N][N])
{
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result[j][i] = A[i][j];
    }
  }
}

//---------
// mat_mul
//---------
//
template <size_t N>
void mat_mul(const double A[N][N], const double B[N][N], double result[N][N])
{
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      double sum = 0.0;

      for (int k = 0; k < N; k++) {
        sum += A[i][k] * B[k][j];
      }

    result[i][j] = sum;
    }
  }
}

template <size_t N>
void mat_mul(const double A[N][N], const Vector<double>& v, double result[N])
{
  for (int i = 0; i < N; i++) {
    double sum = 0.0;

    for (int j = 0; j < N; j++) {
      sum += A[i][j] * v(j);
    }

    result[i] = sum;
  }
}

template <size_t N>
void mat_mul6x3(const double A[2*N][2*N], const Array<double>& B, Array<double>& C)
{
  int B_num_rows = B.nrows();
  int B_num_cols = B.ncols();

  if (B_num_rows != 2*N) {
    throw std::runtime_error("[mat_mul] The number of columns of A (" + std::to_string(2*N) + ") does not equal " +
        " the number of rows of B (" + std::to_string(B_num_rows) + ").");
  }
 
  for (int i = 0; i < 2*N; i++) {
    for (int j = 0; j < B_num_cols; j++) {
      double sum = 0.0;

      for (int k = 0; k < 2*N; k++) {
        sum += A[i][k] * B(k,j);
      }

      C(i,j) = sum;
    }
  }
}

template <size_t N>
void mat_inv(const double A[N][N], double Ainv[N][N])
{
  int iok = 0;

  if (N == 2) {
    double d = mat_det<N>(A);
    if (utils::is_zero(fabs(d))) {
      iok = -1;
    }

    Ainv[0][0] =  A[1][1] / d;
    Ainv[0][1] = -A[0][1] / d;

    Ainv[1][0] = -A[1][0] / d;
    Ainv[1][1] =  A[0][0] / d;

  } else if (N == 3) {
    double d = mat_det<N>(A);
    if (utils::is_zero(fabs(d))) {
      iok = -1;
    }

    Ainv[0][0] = (A[1][1]*A[2][2] - A[1][2]*A[2][1]) / d;
    Ainv[0][1] = (A[0][2]*A[2][1] - A[0][1]*A[2][2]) / d;
    Ainv[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1]) / d;

    Ainv[1][0] = (A[1][2]*A[2][0] - A[1][0]*A[2][2]) / d;
    Ainv[1][1] = (A[0][0]*A[2][2] - A[0][2]*A[2][0]) / d;
    Ainv[1][2] = (A[0][2]*A[1][0] - A[0][0]*A[1][2]) / d;

    Ainv[2][0] = (A[1][0]*A[2][1] - A[1][1]*A[2][0]) / d;
    Ainv[2][1] = (A[0][1]*A[2][0] - A[0][0]*A[2][1]) / d;
    Ainv[2][2] = (A[0][0]*A[1][1] - A[0][1]*A[1][0]) / d;
  }
}

template <size_t N>
double mat_trace(const double A[N][N])
{
  double result = 0.0;

  for (int i = 0; i < N; i++) {
    result += A[i][i];
  }

  return result;
}


template <size_t N>
void ten_zero(double A[N][N][N][N])
{
  memset(A, 0, N*sizeof(A[0]));
}

template <size_t N>
void ten_dyad_prod(const double A[N][N], const double B[N][N], double C[N][N][N][N])
{  
  int nn = pow(N,4);
 
  for (int ii = 0; ii < nn; ii++) {
    int i = t_ind(0,ii);
    int j = t_ind(1,ii);
    int k = t_ind(2,ii);
    int l = t_ind(3,ii);
    C[i][j][k][l] = A[i][j] * B[k][l];
  }
}

template <size_t N>
void ten_ids(double A[N][N][N][N])
{
  ten_zero(A);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      A[i][j][i][j] += 0.5;
      A[i][j][j][i] += 0.5;
    }
  }
}

template <size_t N>
void mat_dyad_prod(const Vector<double>& u, const Vector<double>& v, double result[N][N])
{
  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      result[i][j] = u(i) * v(j);
    }
  }
}

template <size_t N>
void mat_symm_prod(const Vector<double>& u, const Vector<double>& v, double result[N][N])
{
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      result[i][j] = 0.5 * (u(i)*v(j) + u(j)*v(i));
    }
  }
}

template <size_t N>
void ten_dyad_prod(const Array<double>& A, const Array<double>& B, double C[N][N][N][N])
{   
  int nn = pow(N,4);

  for (int ii = 0; ii < nn; ii++) {
    int i = t_ind(0,ii);
    int j = t_ind(1,ii);
    int k = t_ind(2,ii);
    int l = t_ind(3,ii);
    C[i][j][k][l] = A(i,j) * B(k,l);
  }
}

template <size_t N>
void ten_symm_prod(const double A[N][N], const double B[N][N], double C[N][N][N][N])
{
  int nn = pow(N,4);
 
  for (int ii = 0; ii < nn; ii++) {
    int i = t_ind(0,ii);
    int j = t_ind(1,ii);
    int k = t_ind(2,ii);
    int l = t_ind(3,ii);
    C[i][j][k][l] = 0.5 * (A[i][k] * B[j][l]  +  A[i][l] * B[j][k] ) ;
    //C(i,j,k,l) = 0.5* ( A(i,k)*B(j,l) + A(i,l)*B(j,k) );
  }

}

template <size_t N>
double mat_ddot(const double A[N][N], const double B[N][N])
{
  double s = 0.0;

  for (int j = 0; j < N; j++) {
    for (int i = 0; i < N; i++) {
      s = s + A[i][j] * B[i][j];
    }
  }

  return s;
}

template <size_t N>
void ten_transpose(const double A[N][N][N][N], double result[N][N][N][N])
{
  int nn = pow(N,4);
 
  for (int ii = 0; ii < nn; ii++) {
    int i = t_ind(0,ii);
    int j = t_ind(1,ii);
    int k = t_ind(2,ii);
    int l = t_ind(3,ii);
    result[i][j][k][l] = A[k][l][i][j];
  }
}

template <size_t N>
void ten_ddot(const double A[N][N][N][N], const double B[N][N][N][N], double C[N][N][N][N])
{
  int nn = pow(N,4);
  ten_zero(C);

  if (N == 2) {
    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C[i][j][k][l] += A[i][j][0][0] * B[k][l][0][0] +
                       A[i][j][0][1] * B[k][l][0][1] +
                       A[i][j][1][0] * B[k][l][1][0] +
                       A[i][j][1][1] * B[k][l][1][1];
    }

  } else {

    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C[i][j][k][l] += A[i][j][0][0]*B[k][l][0][0] + 
                       A[i][j][0][1]*B[k][l][0][1] + 
                       A[i][j][0][2]*B[k][l][0][2] + 
                       A[i][j][1][0]*B[k][l][1][0] + 
                       A[i][j][1][1]*B[k][l][1][1] + 
                       A[i][j][1][2]*B[k][l][1][2] + 
                       A[i][j][2][0]*B[k][l][2][0] + 
                       A[i][j][2][1]*B[k][l][2][1] + 
                       A[i][j][2][2]*B[k][l][2][2];

    }
  }

}

template <size_t N>
double norm(const double u[N], const double v[N])
{
  double norm = 0.0;
  for (int i = 0; i < N; i++) {
    norm += u[i]*v[i];
  }
  return norm;
}

template <size_t N>
double norm(const Vector<double>& u, const double v[N])
{
  double norm = 0.0;
  for (int i = 0; i < N; i++) {
    norm += u[i]*v[i];
  }
  return norm;
}


template <size_t N>
void mat_symm(const double A[N][N], double S[N][N])
{
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      S[i][j] = 0.5 * (A[i][j] + A[j][i]);
    }
  }
}

template <size_t N>
void mat_dev(const double A[N][N], double dev[N][N])
{
  double trA = mat_trace<N>(A);
  double id[N][N];
  mat_id<N>(id);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      dev[i][j] = A[i][j] - (trA / static_cast<double>(N)) * id[i][j];
    }
  }
}


template <size_t N>
void ten_ddot_3424(const double A[N][N][N][N], const double B[N][N][N][N], double C[N][N][N][N])
{
  int nn = pow(N,4);
  ten_zero(C);

  if (N == 2) {
    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C[i][j][k][l] += A[i][j][0][0]*B[k][0][l][0] +
                       A[i][j][0][1]*B[k][0][l][1] +
                       A[i][j][1][0]*B[k][1][l][0] +
                       A[i][j][1][1]*B[k][1][l][1];
    }

  } else {

    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C[i][j][k][l] += A[i][j][0][0]*B[k][0][l][0] + 
                       A[i][j][0][1]*B[k][0][l][1] + 
                       A[i][j][0][2]*B[k][0][l][2] + 
                       A[i][j][1][0]*B[k][1][l][0] + 
                       A[i][j][1][1]*B[k][1][l][1] + 
                       A[i][j][1][2]*B[k][1][l][2] + 
                       A[i][j][2][0]*B[k][2][l][0] + 
                       A[i][j][2][1]*B[k][2][l][1] + 
                       A[i][j][2][2]*B[k][2][l][2];

    }
  }
}

template <size_t N>
void ten_ddot_2412(const double A[N][N][N][N], const double B[N][N][N][N], double C[N][N][N][N])
{
  int nn = pow(N,4);
  ten_zero(C);

  if (N == 2) {
    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C[i][j][k][l] += A[i][0][j][0]*B[0][0][k][l] + 
                       A[i][0][j][1]*B[0][1][k][l] + 
                       A[i][1][j][0]*B[1][0][k][l] + 
                       A[i][1][j][1]*B[1][1][k][l];

    }

  } else {

    for (int ii = 0; ii < nn; ii++) {
      int i = t_ind(0,ii);
      int j = t_ind(1,ii);
      int k = t_ind(2,ii);
      int l = t_ind(3,ii);

      C[i][j][k][l] += A[i][0][j][0]*B[0][0][k][l] + 
                       A[i][0][j][1]*B[0][1][k][l] + 
                       A[i][0][j][2]*B[0][2][k][l] + 
                       A[i][1][j][0]*B[1][0][k][l] + 
                       A[i][1][j][1]*B[1][1][k][l] + 
                       A[i][1][j][2]*B[1][2][k][l] + 
                       A[i][2][j][0]*B[2][0][k][l] + 
                       A[i][2][j][1]*B[2][1][k][l] + 
                       A[i][2][j][2]*B[2][2][k][l];
    }
  }
}


};

#endif

