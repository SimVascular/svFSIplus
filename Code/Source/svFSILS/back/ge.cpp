
// This subroutine solve a linear system of equations AX=B using
// Gauss elimination and replaces the B with X

#include "ge.h"

#include <math.h>

namespace ge {

//----
// ge
//----
//
bool ge(const int nV, const int N, const Array<double>& A, Vector<double>& B)
{
  int tid = 0;
  auto msg_prefix = std::string("[ge:") + std::to_string(tid) + "] ";
  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "========== ge ==========" << std::endl;
  //std::cout << msg_prefix << "nV: " << nV << std::endl;
  //std::cout << msg_prefix << "N: " << N << std::endl;

  // Constructing a preconditioner. This is to prevent latter problem
  // with singular matrix

  Vector<double> W(N);
  double tol = std::numeric_limits<double>::denorm_min();
  double eps = std::numeric_limits<double>::epsilon();
  //std::cout << msg_prefix << "tol: " << tol << std::endl;

  // Constructing a preconditioner.
  //
  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "Constructing a preconditioner ..." << std::endl;

  for (int i = 0; i < N; i++) {
    if (fabs(A(i,i)) < tol) { 
      B = 0.0;
      //std::cout << msg_prefix << "Done" << std::endl;
      return false;
    }

    W(i) = 1.0 / sqrt(fabs(A(i,i)));
    //std::cout << msg_prefix << "W(" << i << "): " << W(i) << std::endl;
  }

  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "Constructing C  ..." << std::endl;
  Array<double> C(N,N+1);

  for (int i = 0; i < N; i++) {
    for (int j = 0; j < N; j++) {
      C(i,j) = W(i)*W(j)*A(i,j);
      //std::cout << msg_prefix << "C(" << i+1 << "," << j+1 << "): " << C(i,j) << std::endl;
    }
    C(i,N) = W(i)*B(i);
  }

  if (N <= 0) {
    //std::cout << msg_prefix << "Done" << std::endl;
    return false;

  } else if (N == 1) {
    B(0) = C(0,1) / C(0,0);
    B(0) = B(0)*W(0);
    //std::cout << msg_prefix << "Done" << std::endl;
    return true;

  } else if (N == 2) {
    double pivot = C(0,0)*C(1,1) - C(1,0)*C(0,1);

    // Singular matrix
    if (fabs(pivot) < eps) {
      B  = 0.0;
      //std::cout << msg_prefix << "Done" << std::endl;
      return false;
    }

    B(0) = (C(0,2)*C(1,1) - C(1,2)*C(0,1)) / pivot;
    B(1) = (C(1,2)*C(0,0) - C(0,2)*C(1,0)) / pivot;
    B(0) = W(0)*B(0);
    B(1) = W(1)*B(1);
    //std::cout << msg_prefix << "Done" << std::endl;
    return true;
  }

  //std::cout << msg_prefix << "General with pivots .. " << std::endl;
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
      //std::cout << msg_prefix << "Done" << std::endl;
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
  //for (int j=N, 1, -1
    for (int i = j+1; i < N; i++) { 
      C(j,N) = C(j,N) - C(j,i)*C(i,N);
      //C(j,N+1) = C(j,N+1) - C(j,i)*C(i,N+1);
    }
    C(j,N) = C(j,N) / C(j,j);
    //C(j,N+1) = C(j,N+1) / C(j,j);
  }

  for (int i = 0; i < N; i++) {
    B(i) = W(i)*C(i,N);
    //B(i) = W(i)*C(i,N+1);
  }

  //std::cout << msg_prefix << "Done" << std::endl;
  return true;
}

};


