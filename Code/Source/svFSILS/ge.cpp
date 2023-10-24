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

// This subroutine solve a linear system of equations AX=B using
// Gauss elimination and replaces the B with X

#include "ge.h"

#include <math.h>

namespace ge {

bool ge(const int nV, const int N, const Array<double>& A, Vector<double>& B)
{
  // Constructing a preconditioner. This is to prevent latter problem
  // with singular matrix

  Vector<double> W(N);
  double tol = std::numeric_limits<double>::denorm_min();
  double eps = std::numeric_limits<double>::epsilon();

  // Constructing a preconditioner.
  //
  for (int i = 0; i < N; i++) {
    if (fabs(A(i,i)) < tol) { 
      B = 0.0;
      return false;
    }

    W(i) = 1.0 / sqrt(fabs(A(i,i)));
  }

  Array<double> C(N,N+1);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      C(i,j) = W(i)*W(j)*A(i,j);
    }
    C(i,N) = W(i)*B(i);
  }

  if (N <= 0) {
    return false;

  } else if (N == 1) {
    B(0) = C(0,1) / C(0,0);
    B(0) = B(0)*W(0);
    return true;

  } else if (N == 2) {
    double pivot = C(0,0)*C(1,1) - C(1,0)*C(0,1);

    // Singular matrix
    if (fabs(pivot) < eps) {
      B  = 0.0;
      return false;
    }

    B(0) = (C(0,2)*C(1,1) - C(1,2)*C(0,1)) / pivot;
    B(1) = (C(1,2)*C(0,0) - C(0,2)*C(1,0)) / pivot;
    B(0) = W(0)*B(0);
    B(1) = W(1)*B(1);
    return true;
  }

  for (int m = 0; m < N-1; m++) {
    int ipv = m;
    double pivot = fabs(C(m,m));

    for (int i = m+1; i < N; i++) {
      if (fabs(C(i,m)) > pivot) {
        ipv = i;
        pivot = fabs(C(i,m));
      }
    }

    // Singular matrix
    if (fabs(pivot) < eps) {
      B = 0.0;
      return false; 
    }

    if (ipv != m) {
      for (int j = m; j < N+1; j++) {
        double saveEl = C(m,j);
        C(m,j) = C(ipv,j);
        C(ipv,j) = saveEl;
      }
    }

    for (int i = m+1; i < N; i++) {
      double saveEl = C(i,m) / C(m,m);
      C(i,m) = 0.0;
      for (int j = m+1; j < N+1; j++) {
        C(i,j) = C(i,j) - saveEl*C(m,j);
      }
    }
  }

  for (int j = N-1; j >= 0; j--) { 
    for (int i = j+1; i < N; i++) { 
      C(j,N) = C(j,N) - C(j,i)*C(i,N);
    }
    C(j,N) = C(j,N) / C(j,j);
  }

  for (int i = 0; i < N; i++) {
    B(i) = W(i)*C(i,N);
  }

  return true;
}

};


