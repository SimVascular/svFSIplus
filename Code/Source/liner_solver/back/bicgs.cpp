
//------------------------------------------------------------------
// Biconjugate-gradient algorithm, available for scaler and vectors.
//-------------------------------------------------------------------

#include "bicgs.h"

#include "fsils_api.hpp"

#include "add_bc_mul.h"
#include "bcast.h"
#include "dot.h"
#include "norm.h"
#include "omp_la.h"
#include "spar_mul.h"

#include "Array3.h"

#include <math.h>

namespace bicgs {

//--------
// bicgsv
//--------
//
void bicgsv (fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof, 
    const Array<double>& K, Array<double>& R)
{
  #define n_debug_bicgsv
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[bicgsv:") + std::to_string(tid) + "] ";
  #ifdef debug_bicgsv
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== bicgsv ==========" << std::endl;
  #endif

  using namespace fsi_linear_solver;

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_bicgsv
  std::cout << msg_prefix << "ls.mItr: " << ls.mItr << std::endl;
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "mynNo: " << mynNo << std::endl;
  //K.write(msg_prefix+"K");
  //R.write(msg_prefix+"R");
  //exit(0);
  #endif

  Array<double> P(dof,nNo), Rh(dof,nNo), X(dof,nNo), V(dof,nNo),   
      S(dof,nNo), T(dof,nNo);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  double err = norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
  double errO = err;
  ls.iNorm = err;
  double eps = std::max(ls.absTol,ls.relTol*err);
  double rho = err*err;
  double beta = rho;
  X = 0.0;
  P = R;
  Rh = R;
  int i_itr = 1;
  #ifdef debug_bicgsv
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "err: " << err << std::endl;
  std::cout << msg_prefix << "eps: " << eps << std::endl;
  #endif

  for (int i = 0; i < ls.mItr; i++) {
    #ifdef debug_bicgsv
    std::cout << msg_prefix << std::endl;
    std::cout << msg_prefix << "----- i " << i+1 << " -----" << std::endl;
    std::cout << msg_prefix << "err: " << err << std::endl;
    std::cout << msg_prefix << "eps: " << eps << std::endl;
    #endif
    if (err < eps) { 
      ls.suc = true;
      break;
    }

    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof, K, P, V);
    double alpha = rho / dot::fsils_dot_v(dof, mynNo, lhs.commu, Rh, V);
    S = R - alpha*V;

    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof, K, S, T);
    double omega = norm::fsi_ls_normv(dof, mynNo, lhs.commu, T);
    omega = dot::fsils_dot_v(dof, mynNo, lhs.commu, T, S) / (omega * omega);

    X = X + alpha*P + omega*S;
    R = S - omega*T;

    errO = err;
    err =  norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
    double rhoO  = rho;
    rho = dot::fsils_dot_v(dof, mynNo, lhs.commu, R, Rh);
    beta = rho*alpha / (rhoO*omega);

    #ifdef debug_bicgsv
    std::cout << msg_prefix << "alpha: " << alpha << std::endl;
    std::cout << msg_prefix << "omega: " << omega << std::endl;
    std::cout << msg_prefix << "rho: " << rho << std::endl;
    std::cout << msg_prefix << "beta: " << beta << std::endl;
    #endif

    P = R + beta * (P - omega*V);
    //P.write(msg_prefix+"P");
    //exit(0);
    i_itr += 1;
  } 

  R = X;
  ls.itr = i_itr - 1;
  ls.fNorm = err;
  ls.callD =  fsi_linear_solver::fsils_cpu_t() - ls.callD;
  #ifdef debug_bicgsv
  std::cout << msg_prefix << "ls.itr: " << ls.itr << std::endl;
  #endif

  if (errO < std::numeric_limits<double>::epsilon()) { 
     ls.dB = 0.0;
  } else { 
     ls.dB = 10.0 * log(err / errO);
  }

  //R.write(msg_prefix+"R");
  //exit(0);

}


};



