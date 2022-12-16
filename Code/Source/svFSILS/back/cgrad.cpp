
// Conjugate-gradient algorithm for scaler, vector and Schur
// complement cases.

#include "cgrad.h"

#include "fsils_api.hpp"
#include "add_bc_mul.h"
#include "dot.h"
#include "omp_la.h"
#include "norm.h"
#include "spar_mul.h"

#include <math.h>

namespace cgrad {

//-------
// schur
//-------
// Reproduces 'SUBROUTINE CGRAD_SCHUR(lhs, ls, dof, D, G, L, R)'
//
// Note [DaveP] This is not producing matching results from the Fortran version.
// The MPI calls seem to be fine.
//
void schur(FSILS_lhsType& lhs, FSILS_subLsType& ls, const int dof, const Array<double>& D, 
    const Array<double>& G, const Vector<double>& L, Vector<double>& R)
{
  #define n_debug_schur
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[schur:") + std::to_string(tid) + "] ";
  #ifdef debug_schur
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== schur ==========" << std::endl;
  #endif

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;

  Vector<double> X(nNo), P(nNo), SP(nNo), DGP(nNo); 
  Array<double> GP(dof,nNo), unCondU(dof,nNo);

  double time = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  ls.iNorm = norm::fsi_ls_norms(mynNo, lhs.commu, R);
  double eps = pow(std::max(ls.absTol,ls.relTol*ls.iNorm),2.0);
  double errO = ls.iNorm*ls.iNorm;
  double err = errO;
  #ifdef debug_schur
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "mynNo: " << mynNo << std::endl;
  std::cout << msg_prefix << "ls.iNorm: " << ls.iNorm << std::endl;
  std::cout << msg_prefix << "eps: " << eps << std::endl;
  std::cout << msg_prefix << "errO: " << errO << std::endl;
  #endif

  X = 0.0;
  P = R;
  int last_i = 0;

  #ifdef debug_schur
  D.write(msg_prefix+"D");
  G.write(msg_prefix+"G");
  L.write(msg_prefix+"L");
  R.write(msg_prefix+"R");
  #endif

  for (int i = 0; i < ls.mItr; i++) {
    #ifdef debug_schur
    std::cout << msg_prefix << std::endl;
    std::cout << msg_prefix << "----- i " << i+1 << " -----" << std::endl;
    std::cout << msg_prefix << "err: " << err << std::endl;
    auto istr = "_" + std::to_string(i+1);
    P.write(msg_prefix+"P"+istr);
    #endif
    last_i = i;

    if (err < eps) {
      ls.suc = true;
      break;
    }

    errO = err;

    // GP = G * P
    spar_mul::fsils_spar_mul_sv(lhs, lhs.rowPtr, lhs.colPtr, dof, G, P, GP);
    //CALL FSILS_SPARMULSV(lhs, lhs.rowPtr, lhs.colPtr, dof, G, P,GP)
    #ifdef debug_schur
    GP.write(msg_prefix+"GP"+istr);
    #endif

    //if (ANY(lhs.face.coupledFlag)) {
    for (auto& face : lhs.face) {
      if (face.coupledFlag) {
        auto unCondU = GP;
        add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, GP);
        //CALL ADDBCMUL(lhs, BCOP_TYPE_PRE, dof, unCondU, GP)
        //std::cout << msg_prefix << "#### Face is couple " << std::endl;
        //exit(0);
        break;
      }
    }

    // DGP = K * GP
    spar_mul::fsils_spar_mul_vs(lhs, lhs.rowPtr, lhs.colPtr, dof, D, GP, DGP);
    //CALL FSILS_SPARMULVS(lhs, lhs.rowPtr, lhs.colPtr, dof, D,GP,DGP)
    #ifdef debug_schur
    DGP.write(msg_prefix+"DGP"+istr);
    #endif

    // SP = L * P
    spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, L, P, SP);
    //CALL FSILS_SPARMULSS(lhs, lhs.rowPtr, lhs.colPtr, L, P, SP)
    #ifdef debug_schur
    SP.write(msg_prefix+"SP_a"+istr);
    #endif

    // SP = SP - DGP
    omp_la::omp_sum_s(nNo, -1.0, SP, DGP);
    //CALL OMPSUMS(nNo, -1._LSRP, SP, DGP)
    //!SP = SP - DGP
    #ifdef debug_schur
    SP.write(msg_prefix+"SP_b"+istr);
    #endif

    // Here 
    double alpha = errO / dot::fsils_dot_s(mynNo, lhs.commu, P, SP);
    #ifdef debug_schur
    //alpha = errO/FSILS_DOTS(mynNo, lhs.commu, P, SP)
    std::cout << msg_prefix << "alpha: " << alpha << std::endl;
    #endif

    // X = X + alpha * P
    omp_la::omp_sum_s(nNo, alpha, X, P);
    //CALL OMPSUMS(nNo, alpha, X, P)
    //!X = X + alpha*P
    #ifdef debug_schur
    X.write(msg_prefix+"X"+istr);
    #endif

    // R = R - alpha * SP
    omp_la::omp_sum_s(nNo, -alpha, R, SP);
    //CALL OMPSUMS(nNo, -alpha, R, SP);
    //!R = R - alpha*SP
    #ifdef debug_schur
    R.write(msg_prefix+"R"+istr);
    #endif

    err = norm::fsi_ls_norms(mynNo, lhs.commu, R);
    err = err * err;
    #ifdef debug_schur
    std::cout << msg_prefix << "err: " << err << std::endl;
    std::cout << msg_prefix << "errO/err: " << errO/err << std::endl;
    std::cout << msg_prefix << "err/errO: " << err/errO << std::endl;
    #endif

    // P = P + errO/err * R 
    double c1 = errO / err;
    omp_la::omp_sum_s(nNo, c1, P, R);
    //omp_la::omp_sum_s(nNo, errO/err, P, R);
    //CALL OMPSUMS(nNo, errO/err, P, R)          

    // P = err/errO * P
    double c2 = err / errO;
    omp_la::omp_mul_s(nNo, c2, P);
    //omp_la::omp_mul_s(nNo, err/errO, P);
    //CALL OMPMULS(nNo, err/errO, P)
    //P = R + err/errO * P;
  }

  R = X;
  ls.fNorm = sqrt(err);
  ls.callD = fsi_linear_solver::fsils_cpu_t() - time + ls.callD;
  ls.itr = ls.itr + last_i;
  #ifdef debug_schur
  std::cout << msg_prefix << "errO: " << errO << std::endl;
  std::cout << msg_prefix << "ls.fNorm: " << ls.fNorm << std::endl;
  std::cout << msg_prefix << "ls.itr: " << ls.itr << std::endl;
  #endif

  if (errO < std::numeric_limits<double>::epsilon()) {
    ls.dB = 0.0;
  } else {
    ls.dB = 5.0 * log(err/errO);
  }

  //R.write(msg_prefix+"R");
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);
}

//---------
// cgrad_v
//---------
//
void cgrad_v(FSILS_lhsType& lhs, FSILS_subLsType& ls, const int dof, const Array<double>& K, Array<double>& R)
{
  #define n_debug_cgrad_v 
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[cgrad_v:") + std::to_string(tid) + "] ";
  #ifdef debug_cgrad_v
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== cgrad_v ==========" << std::endl;
  #endif

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_cgrad_v
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "mynNo: " << mynNo << std::endl;
  std::cout << msg_prefix << "ls.mItr: " << ls.mItr << std::endl;
  #endif

  Array<double> P(dof,nNo), KP(dof,nNo), X(dof,nNo);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  ls.iNorm = norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
  double eps = pow(std::max(ls.absTol, ls.relTol* ls.iNorm), 2.0);

  double errO = ls.iNorm * ls.iNorm;
  double err  = errO;
  X = 0.0;
  P = R;
  int last_i = 0;
  #ifdef debug_cgrad_v
  std::cout << msg_prefix << "ls.iNorm: " << ls.iNorm << std::endl;
  std::cout << msg_prefix << "eps: " << eps << std::endl;
  std::cout << msg_prefix << "err: " << eps << std::endl;
  #endif

  for (int i = 0; i < ls.mItr; i++) {
    #ifdef debug_cgrad_v
    std::cout << msg_prefix << std::endl;
    std::cout << msg_prefix << "----- i " << i+1 << " -----" << std::endl;
    std::cout << msg_prefix << "err: " << err << std::endl;
    auto istr = "_" + std::to_string(i+1);
    P.write(msg_prefix+"P"+istr);
    #endif
    last_i = i;

    if (err < eps) {
      ls.suc = true;
      break;
    }

    errO = err;

    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof, K, P, KP);

    double alpha = errO / dot::fsils_dot_v(dof, mynNo, lhs.commu, P, KP);
    omp_la::omp_sum_v(dof, nNo, alpha, X, P);
    omp_la::omp_sum_v(dof, nNo, -alpha, R, KP);
    #ifdef debug_cgrad_v
    std::cout << msg_prefix << "alpha: " << alpha << std::endl;
    #endif

    err = norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
    err = err * err;

    omp_la::omp_sum_v(dof, nNo, errO/err, P, R);
    omp_la::omp_mul_v(dof, nNo, err/errO, P);
  }

  R = X;
  ls.itr = last_i;
  ls.fNorm = sqrt(err);
  ls.callD = fsi_linear_solver::fsils_cpu_t() - ls.callD;

  if (errO < std::numeric_limits<double>::epsilon()) {
    ls.dB = 0.0;
  } else {
    ls.dB = 5.0 * log(err/errO);
  }
}


};


