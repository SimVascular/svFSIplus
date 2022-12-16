
//-------------------------------------------------------------------------
// Graduate [sic] minimum residual algorithm is implemented here for vector
// and scaler problems.
//-------------------------------------------------------------------------

#include "gmres.h"

#include "fsils_api.hpp"

#include "add_bc_mul.h"
#include "bcast.h"
#include "dot.h"
#include "norm.h"
#include "omp_la.h"
#include "spar_mul.h"

#include "Array3.h"

#include <math.h>

namespace gmres {

//--------
// bc_pre
//--------
//
void bc_pre(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof, 
    const int mynNo, const int nNo)
{
  int nsd = dof - 1;
  Array<double> v(nsd,nNo);

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto &face = lhs.face[faIn];
    if (face.coupledFlag) {
      if (face.sharedFlag) {
        v = 0.0;

        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            v(i,Ac) = face.valM(i,a);
          }
        }

        face.nS = pow(norm::fsi_ls_normv(nsd, mynNo, lhs.commu, v), 2.0);

      } else { 
        face.nS = 0.0;
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            face.nS = face.nS + pow(face.valM(i,a), 2.0);
          }
        }
      }
    }
  }
}

//-------
// gmres
//-------
// Solver the system Val * X = R.
//
// Reproduces the Fortran 'GMRES' subroutine.
//
void gmres(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof, 
    const Array<double>& Val, const Array<double>& R, Array<double>& X)
{
  #define n_debug_gmres
  #ifdef debug_gmres
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[gmres:") + std::to_string(tid) + "] ";
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== gmres ==========" << std::endl;
  std::cout << msg_prefix << "lhs.debug_active: " << lhs.debug_active << std::endl;
  #endif

  using namespace fsi_linear_solver;

  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_gmres
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "mynNo: " << mynNo << std::endl;
  std::cout << msg_prefix << "ls.sD: " << ls.sD << std::endl;
  std::cout << msg_prefix << "ls.mItr: " << ls.mItr << std::endl;
  std::cout << msg_prefix << "ls.absTol: " << ls.absTol << std::endl;
  std::cout << msg_prefix << "ls.relTol: " << ls.relTol << std::endl;
  //R.write(msg_prefix+"R");
  //Val.write(msg_prefix+"Val");
  #endif

  Array<double> h(ls.sD+1,ls.sD); 
  Array3<double> u(dof,nNo,ls.sD+1); 
  Array<double> unCondU(dof,nNo);
  Vector<double> y(ls.sD), c(ls.sD), s(ls.sD), err(ls.sD+1);

  double time = fsi_linear_solver::fsils_cpu_t(); 
  ls.suc = false;
  double eps = 0.0;
  int last_i = 0;
  X = 0.0;

  //lhs.rowPtr.write(msg_prefix+"rowPtr");
  //lhs.colPtr.write(msg_prefix+"colPtr");

  //R.write(msg_prefix+"R");
  //Val.write(msg_prefix+"Val");
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);

  for (int l = 0; l < ls.mItr; l++) {
    #ifdef debug_gmres
    std::cout << msg_prefix << std::endl;
    std::cout << msg_prefix << "---------- l " << l+1 << " ----------" << std::endl;
    #endif

    if (l == 0) {
      u.set_slice(0, R);
      //u(:,:,1) = R
      //u.write(msg_prefix+"u_0");
    } else {
      auto u_slice = u.slice(0);
      spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof,  Val, X, u_slice);
      //CALL FSILS_SPARMULVV(lhs, lhs.rowPtr, lhs.colPtr, dof,  Val, X, u(:,:,1))

      add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, dof, X, u_slice);
      //CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, dof, X, u(:,:,1))

      //u_slice.write(msg_prefix+"u_slice");
      //exit(0);

      ls.itr = ls.itr + 1;
      u.set_slice(0, R - u_slice);
      //u(:,:,1) = R - u(:,:,1)
    }

    //if (l == 0) {
    //u.write(msg_prefix+"u");
    //}
    //MPI_Barrier(lhs.commu.comm);
    //exit(0);

    for (auto& face : lhs.face) {
    //if (ANY(lhs.face.coupledFlag)) {
      if (face.coupledFlag) {
        auto u_slice = u.slice(0);
        auto unCondU = u.slice(0);
        add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, u_slice);
        //unCondU = u(:,:,1)
        //CALL ADDBCMUL(lhs, BCOP_TYPE_PRE, dof, unCondU, u(:,:,1))
        u.set_slice(0, u_slice);
        //std::cout << msg_prefix << "###### Has coupled face" << std::endl;
        break; 
      }
    }

    err[0] = norm::fsi_ls_normv(dof, mynNo, lhs.commu, u.slice(0));
    #ifdef debug_gmres
    std::cout << msg_prefix << "err(1): " << err[0] << std::endl;
    #endif

    if (l == 0) {
      eps = err[0];

      if (eps <= ls.absTol) {
        ls.callD = std::numeric_limits<double>::epsilon();
        ls.dB = 0.0;
        return; 
      }

      ls.iNorm = eps;
      ls.fNorm = eps;
      eps = std::max(ls.absTol, ls.relTol*eps);
    }
    #ifdef debug_gmres
    std::cout << msg_prefix << "eps: " << eps << std::endl;
    #endif

    ls.dB = ls.fNorm;
    auto u_slice = u.slice(0) / err(0);
    u.set_slice(0, u_slice);
    //u(:,:,1) = u(:,:,1) / err(1)

    //Val.write(msg_prefix+"Val");
    //exit(0);

    for (int i = 0; i < ls.sD; i++) {
      #ifdef debug_gmres
      std::cout << msg_prefix << std::endl;
      std::cout << msg_prefix << "----- i " << i+1 << " -----" << std::endl;
      #endif
      last_i = i;
      auto u_slice = u.slice(i);
      auto u_slice_1 = u.slice(i+1);
      spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof,  Val, u_slice, u_slice_1);
      //CALL FSILS_SPARMULVV(lhs, lhs.rowPtr, lhs.colPtr, dof, Val, u(:,:,i), u(:,:,i+1))

      //if (i+1 == 3) {
      //  lhs.debug_active = false;
      //}

      add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, dof, u_slice, u_slice_1);
      u.set_slice(i+1, u_slice_1);
      //CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, dof, u(:,:,i), u(:,:,i+1))

      //u.write(msg_prefix+"u");
      //MPI_Barrier(lhs.commu.comm);
      //exit(0);

      ls.itr = ls.itr + 1;

      //if (ANY(lhs.face.coupledFlag)) {
      for (auto& face : lhs.face) {
        if (face.coupledFlag) {
          auto u_slice_1 = u.slice(i+1);
          auto unCondU = u.slice(i+1);
          add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, u_slice_1);
          u.set_slice(i+1, u_slice_1);
          //unCondU = u(:,:,i+1)
          //CALL ADDBCMUL(lhs,BCOP_TYPE_PRE,dof, unCondU, u(:,:,i+1))
          break;
        }
      }

      if (i+1 == 1 && lhs.debug_active) {
        //u.write(msg_prefix+"u");
        //MPI_Barrier(lhs.commu.comm);
        //exit(0);
      }

      // [TODO:DaveP] careful here, should j <= i+1 ?   
      for (int j = 0; j <= i+1; j++) {
      //DO j=1, i+1
        h(j,i) = dot::fsils_nc_dot_v(dof, mynNo, u.slice(j), u.slice(i+1));
        //h(j,i) = FSILS_NCDOTV(dof, mynno, u(:,:,j), u(:,:,i+1))
        //std::cout << msg_prefix << "1: h(" << j+1 << "," << i+1 << "): " << h(j,i) << std::endl;
      }

      auto h_col = h.col(i);
      bcast::fsils_bcast_v(i+2, h_col, lhs.commu);
      //bcast::fsils_bcast_v(i+1, h_col, lhs.commu);
      h.set_col(i, h_col);
      //CALL FSILS_BCASTV(i+1, h(:,i), lhs.commu)
      //h.write(msg_prefix+"h");
      //MPI_Barrier(lhs.commu.comm);
      //exit(0);

      // [TODO:DaveP] careful here, should j <= i ?   
      for (int j = 0; j <= i; j++) {
        auto u_slice_1 = u.slice(i+1);
        omp_la::omp_sum_v(dof, nNo, -h(j,i), u_slice_1, u.slice(j));
        u.set_slice(i+1, u_slice_1);
        //CALL OMPSUMV(dof, nNo, -h(j,i), u(:,:,i+1), u(:,:,j))
        h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i);
        //std::cout << msg_prefix << "h(" << i+2 << "," << i+1 << "): " << h(i+1,i) << std::endl;
      }

      h(i+1,i) = sqrt(fabs(h(i+1,i)));
      //std::cout << msg_prefix << "h(" << i+2 << "," << i+1 << "): " << h(i+1,i) << std::endl;

      u_slice_1 = u.slice(i+1);
      omp_la::omp_mul_v(dof, nNo, 1.0/h(i+1,i), u_slice_1);
      u.set_slice(i+1, u_slice_1);
      //CALL OMPMULV(dof, nNo, 1.0/h(i+1,i), u(:,:,i+1))

      // [TODO:DaveP] careful here, should j <= i-1 ?   
      for (int j = 0; j <= i-1; j++) {
      //DO j=1, i-1
        double tmp = c(j)*h(j,i) + s(j)*h(j+1,i);
        h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i);
        h(j,i) = tmp;
        //std::cout << msg_prefix << "### h(j,i): " << h(j,i) << std::endl;
      }

      //std::cout << msg_prefix << "h(i,i): " << h(i,i) << std::endl;
      //std::cout << msg_prefix << "h(i+1,i): " << h(i+1,i) << std::endl;

      double tmp = sqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i));
      c(i) = h(i,i) / tmp;
      s(i) = h(i+1,i) / tmp;
      h(i,i) = tmp;
      h(i+1,i) = 0.0;
      err(i+1) = -s(i)*err(i);
      err(i) = c(i)*err(i);
      #ifdef debug_gmres
      std::cout << msg_prefix << std::endl;
      std::cout << msg_prefix << "tmp: " << tmp << std::endl;
      std::cout << msg_prefix << "err(i): " << err(i) << std::endl;
      std::cout << msg_prefix << "err(i+1): " << err(i+1) << std::endl;
      std::cout << msg_prefix << "eps: " << eps << std::endl;
      #endif

      //u.write(msg_prefix+"u");
      //h.write(msg_prefix+"h");
      //MPI_Barrier(lhs.commu.comm);
      //exit(0);

      if (fabs(err(i+1)) < eps) {
        ls.suc = true;
        break;
      }
    } // for int i = 0; i < ls.sD
    #ifdef debug_gmres
    std::cout << msg_prefix << "End of i loop" << std::endl;
    #endif


    //h.write(msg_prefix+"h");
    //exit(0);
    //std::cout << msg_prefix << std::endl;
    //std::cout << msg_prefix << "last_i:  " << last_i << std::endl;

    // [TODO:DaveP] careful here, last_i is used to index arrays.
    //
    if (last_i >= ls.sD) {
      last_i = ls.sD - 1;
    }

    for (int i = 0; i <= last_i; i++) {
      y(i) = err(i);
    }
    //y = err(1:i)

    //std::cout << msg_prefix << std::endl;
    //std::cout << msg_prefix << "y: " << std::endl;
    for (int j = last_i; j >= 0; j--) { 
      for (int k = j+1; k <= last_i; k++) {
        y(j) = y(j) - h(j,k)*y(k);
      }
      y(j) = y(j) / h(j,j);
      //std::cout << msg_prefix << j+1 << " " << y(j) << std::endl;
    }

    for (int j = 0; j <= last_i; j++) {
      omp_la::omp_sum_v(dof, nNo, y(j), X, u.slice(j));
      //CALL OMPSUMV(dof, nNo, y(j), X, u(:,:,j))
    }

    ls.fNorm = fabs(err(last_i+1));
    if (ls.suc) {
      break;
    }

  } // for l = 0; l < ls.mItr
  #ifdef debug_gmres
  std::cout << msg_prefix << "End of l loop" << std::endl;
  #endif
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);

  ls.callD = fsi_linear_solver::fsils_cpu_t() - time + ls.callD;
  ls.dB  = 10.0 * log(ls.fNorm / ls.dB);

  #ifdef debug_gmres
  std::cout << msg_prefix << "Done" << std::endl;
  #endif
}

//---------
// gmres_s
//---------
// Reproduces the Fortran 'GMRESS' subroutine.
//
void gmres_s(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof,
    const Vector<double>& Val, Vector<double>& R)
{
  #define n_debug_gmres_s
  #ifdef debug_gmres_s
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[gmres_s:") + std::to_string(tid) + "] ";
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== gmres_s ==========" << std::endl;
  #endif

  using namespace fsi_linear_solver;

  bool flag = false;
  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_gmres_s
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "mynNo: " << mynNo << std::endl;
  std::cout << msg_prefix << "ls.sD: " << ls.sD << std::endl;
  std::cout << msg_prefix << "ls.mItr: " << ls.mItr << std::endl;
  std::cout << msg_prefix << "ls.absTol: " << ls.absTol << std::endl;
  std::cout << msg_prefix << "ls.relTol: " << ls.relTol << std::endl;
  #endif

  Array<double> h(ls.sD+1,ls.sD);
  Array<double> u(nNo,ls.sD+1);
  Vector<double> X(nNo), y(ls.sD), c(ls.sD), s(ls.sD), err(ls.sD+1);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  double eps = norm::fsi_ls_norms(mynNo, lhs.commu, R);
  ls.iNorm = eps;
  ls.fNorm = eps;
  eps = std::max(ls.absTol, ls.relTol*eps);
  ls.itr = 0;
  int last_i = 0;
  #ifdef debug_gmres_s
  std::cout << msg_prefix << "ls.iNorm: " << ls.iNorm << std::endl;
  std::cout << msg_prefix << "eps: " << eps << std::endl;
  #endif

  if (ls.iNorm <= ls.absTol) {
    ls.callD = std::numeric_limits<double>::epsilon();
    ls.dB = 0.0;
    return; 
  }

  //Vector<double>::write_disabled = false;
  //R.write(msg_prefix+"R");
  //Val.write(msg_prefix+"Val");

  for (int l = 0; l < ls.mItr; l++) {
    #ifdef debug_gmres_s
    std::cout << msg_prefix << std::endl;
    std::cout << msg_prefix << "======== l " << l+1 << " ======== " << std::endl;
    #endif
    ls.dB = ls.fNorm;
    ls.itr = ls.itr + 1;
    auto u_col = u.col(0);
    spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, Val, X, u_col);
    //CALL FSILS_SPARMULSS(lhs, lhs%rowPtr, lhs%colPtr, Val,X,u(:,1))

    u.set_col(0, R - u_col);
    //u(:,:,1) = R - u(:,:,1)

    err[0] = norm::fsi_ls_norms(mynNo, lhs.commu, u.col(0));
    u_col = u.col(0) / err[0];
    u.set_col(0, u_col);
    //u(:,:,1) = u(:,:,1)/err(1)
    #ifdef debug_gmres_s
    std::cout << msg_prefix << "err(1): " << err[0] << std::endl;
    #endif


    for (int i = 0; i < ls.sD; i++) {
      #ifdef debug_gmres_s
      std::cout << msg_prefix << std::endl;
      std::cout << msg_prefix << "----- i " << i+1 << " ----- " << std::endl;
      #endif
      ls.itr = ls.itr + 1;
      last_i = i;
      auto u_col = u.col(i);
      auto u_col_1 = u.col(i+1);
      spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, Val, u_col, u_col_1);
      //CALL FSILS_SPARMULSS(lhs, lhs%rowPtr, lhs%colPtr, Val, u(:,i), u(:,i+1))
      u.set_col(i+1, u_col_1);

      for (int j = 0; j <= i+1; j++) {
        h(j,i) = dot::fsils_nc_dot_s(mynNo, u.col(j), u.col(i+1));
        //std::cout << msg_prefix << "h(j,i): " << h(j,i) << std::endl;
      }

#if 0
      if (i+1 == 10) { 
        Array<double>::write_disabled = false;
        u.write(msg_prefix+"u");
        MPI_Barrier(lhs.commu.comm);
        exit(0);
      }
#endif

      auto h_col = h.col(i);
      bcast::fsils_bcast_v(i+2, h_col, lhs.commu);
      //CALL FSILS_BCASTV(i+1, h(:,i), lhs.commu)
      h.set_col(i, h_col);

      for (int j = 0; j <= i; j++) {
        auto u_col_1 = u.col(i+1);
        omp_la::omp_sum_s(nNo, -h(j,i), u_col_1, u.col(j));
        u.set_col(i+1, u_col_1);
        h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i);
      }
      h(i+1,i) = sqrt(fabs(h(i+1,i)));

      u_col_1 = u.col(i+1);
      omp_la::omp_mul_s(nNo, 1.0/h(i+1,i), u_col_1);
      u.set_col(i+1, u_col_1);
      //CALL OMPMULV(dof, nNo, 1._LSRP/h(i+1,i), u(:,:,i+1))

      for (int j = 0; j <= i-1; j++) {
        double tmp = c(j)*h(j,i) + s(j)*h(j+1,i);
        h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i);
        h(j,i) = tmp;
      }

      double tmp = sqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i));
      c(i) = h(i,i) / tmp;
      s(i) = h(i+1,i) / tmp;
      h(i,i) = tmp;
      h(i+1,i) = 0.0;
      err(i+1) = -s(i)*err(i);
      err(i) = c(i)*err(i);
      #ifdef debug_gmres_s
      std::cout << msg_prefix << "err(i+1): " << err(i+1) << std::endl;
      std::cout << msg_prefix << "tmp: " << tmp << std::endl;
      #endif

      if (fabs(err(i+1)) < eps) {
        ls.suc = true;
        break;
      }
    } // for int i = 0; i < ls.sD

    if (last_i >= ls.sD) {
      last_i = ls.sD - 1;
    }

    for (int i = 0; i <= last_i; i++) {
      y(i) = err(i);
    }
    //y = err(1:i)

    for (int j = last_i; j >= 0; j--) { 
      for (int k = j+1; k <= last_i; k++) {
        y(j) = y(j) - h(j,k)*y(k);
      }
      y(j) = y(j) / h(j,j);
      //std::cout << msg_prefix << j+1 << " " << y(j) << std::endl;
    }

    for (int j = 0; j <= last_i; j++) {
      omp_la::omp_sum_s(nNo, y(j), X, u.col(j));
    }

    ls.fNorm = fabs(err(last_i+1));
    if (ls.suc) {
      break;
    }
  }

  R = X;
  ls.callD = fsi_linear_solver::fsils_cpu_t() - ls.callD;
  ls.dB  = 10.0 * log(ls.fNorm / ls.dB);

  //exit(0);
}

//---------
// gmres_v
//---------
// Reproduces the Fortran 'GMRESV' subroutine.
//
void gmres_v(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof,
    const Array<double>& Val, Array<double>& R)
{
  #define n_debug_gmres_v
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[gmres_v:") + std::to_string(tid) + "] ";
  #ifdef debug_gmres_v
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== gmres_v ==========" << std::endl;
  #endif
  double time = fsi_linear_solver::fsils_cpu_t();

  using namespace fsi_linear_solver;

  bool flag = false;
  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  #ifdef debug_gmres_v
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "mynNo: " << mynNo << std::endl;
  std::cout << msg_prefix << "ls.sD: " << ls.sD << std::endl;
  std::cout << msg_prefix << "ls.mItr: " << ls.mItr << std::endl;
  std::cout << msg_prefix << "ls.absTol: " << ls.absTol << std::endl;
  std::cout << msg_prefix << "ls.relTol: " << ls.relTol << std::endl;
  //Array<double>::write_disabled = false;
  R.write(msg_prefix+"R");
  Val.write(msg_prefix+"Val");
  //exit(0);
  #endif

  Array<double> h(ls.sD+1,ls.sD), X(dof,nNo);
  Array3<double> u(dof,nNo,ls.sD+1);
  Array<double> unCondU(dof,nNo);
  Vector<double> y(ls.sD), c(ls.sD), s(ls.sD), err(ls.sD+1);

  ls.callD = fsi_linear_solver::fsils_cpu_t();
  ls.suc = false;
  double eps = norm::fsi_ls_normv(dof, mynNo, lhs.commu, R);
  ls.iNorm = eps;
  ls.fNorm = eps;
  eps = std::max(ls.absTol, ls.relTol*eps);
  ls.itr = 0;
  int last_i = 0;
  #ifdef debug_gmres_v
  std::cout << msg_prefix << "ls.iNorm: " << ls.iNorm << std::endl;
  std::cout << msg_prefix << "eps: " << eps << std::endl;
  #endif

  bc_pre(lhs, ls, dof, mynNo, nNo);

  if (ls.iNorm <= ls.absTol) {
    ls.callD = std::numeric_limits<double>::epsilon();
    ls.dB = 0.0;
    return; 
  }

  for (int l = 0; l < ls.mItr; l++) {
    //std::cout << msg_prefix << "[GMRES] ===== l " << l+1 << " ===== " << std::endl;
    #ifdef debug_gmres_v
    std::cout << msg_prefix << std::endl;
    std::cout << msg_prefix << "===== l " << l+1 << " ===== " << std::endl;
    #endif
    ls.dB = ls.fNorm;
    ls.itr = ls.itr + 1;
    auto u_slice = u.slice(0);
    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof,  Val, X, u_slice);
    //CALL FSILS_SPARMULVV(lhs, lhs.rowPtr, lhs.colPtr, dof, Val, X, u(:,:,1))

    add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, dof, X, u_slice);
    //CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, dof, X, u(:,:,1))

    u.set_slice(0, R - u_slice);
    //u(:,:,1) = R - u(:,:,1)

    for (auto& face : lhs.face) {
    //if (ANY(lhs.face.coupledFlag)) {
      if (face.coupledFlag && flag) {
        auto u_slice = u.slice(0);
        auto unCondU = u.slice(0);
        add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, u_slice);
        u.set_slice(0, u_slice);
        std::cout << msg_prefix << "face.coupledFlag " << std::endl;
        break;
      }
    }

    err[0] = norm::fsi_ls_normv(dof, mynNo, lhs.commu, u.slice(0));
    u_slice = u.slice(0) / err[0];
    u.set_slice(0, u_slice);
    //u(:,:,1) = u(:,:,1)/err(1)
    #ifdef debug_gmres_v
    std::cout << msg_prefix << "err(1): " << err[0] << std::endl;
    #endif

    for (int i = 0; i < ls.sD; i++) {
      #ifdef debug_gmres_v
      std::cout << msg_prefix << std::endl;
      std::cout << msg_prefix << "----- i " << i+1 << " ----- " << std::endl;
      #endif
      ls.itr = ls.itr + 1;
      last_i = i;
      auto u_slice = u.slice(i);
      auto u_slice_1 = u.slice(i+1);
      spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, dof,  Val, u_slice, u_slice_1);
      //CALL FSILS_SPARMULVV(lhs, lhs.rowPtr, lhs.colPtr, dof, Val, u(:,:,i), u(:,:,i+1))

      add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, dof, u_slice, u_slice_1);
      u.set_slice(i+1, u_slice_1);
      //CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, dof, u(:,:,i), u(:,:,i+1))

      for (auto& face : lhs.face) {
        if (face.coupledFlag && flag) {
          auto u_slice_1 = u.slice(i+1);
          auto unCondU = u.slice(i+1);
          add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_PRE, dof, unCondU, u_slice_1);
          u.set_slice(i+1, u_slice_1);
          //unCondU = u(:,:,i+1)
          //CALL ADDBCMUL(lhs,BCOP_TYPE_PRE,dof, unCondU, u(:,:,i+1))
          //std::cout << msg_prefix << "face.coupledFlag " << std::endl;
          break;
        }
      }
      //if (ANY(lhs.face.coupledFlag).AND.flag) {
      //  unCondU = u(:,:,i+1)
      //  CALL ADDBCMUL(lhs,BCOP_TYPE_PRE,dof, unCondU, u(:,:,i+1))
      //}

      for (int j = 0; j <= i+1; j++) {
        h(j,i) = dot::fsils_nc_dot_v(dof, mynNo, u.slice(j), u.slice(i+1));
        //h(j,i) = FSILS_NCDOTV(dof, mynno, u(:,:,j), u(:,:,i+1))
        #ifdef debug_gmres_v
        std::cout << msg_prefix << "h(j,i): " << h(j,i) << std::endl;
        #endif
      }

      auto h_col = h.col(i);
      bcast::fsils_bcast_v(i+2, h_col, lhs.commu);
      //CALL FSILS_BCASTV(i+1, h(:,i), lhs.commu)
      h.set_col(i, h_col);

      for (int j = 0; j <= i; j++) {
        auto u_slice_1 = u.slice(i+1);
        omp_la::omp_sum_v(dof, nNo, -h(j,i), u_slice_1, u.slice(j));
        u.set_slice(i+1, u_slice_1);
        //CALL OMPSUMV(dof, nNo, -h(j,i), u(:,:,i+1), u(:,:,j))
        h(i+1,i) = h(i+1,i) - h(j,i)*h(j,i);
        //std::cout << msg_prefix << "h(" << i+2 << "," << i+1 << "): " << h(i+1,i) << std::endl;
      }
      h(i+1,i) = sqrt(fabs(h(i+1,i)));

      u_slice_1 = u.slice(i+1);
      omp_la::omp_mul_v(dof, nNo, 1.0/h(i+1,i), u_slice_1);
      u.set_slice(i+1, u_slice_1);
      //CALL OMPMULV(dof, nNo, 1._LSRP/h(i+1,i), u(:,:,i+1))

      for (int j = 0; j <= i-1; j++) {
        double tmp = c(j)*h(j,i) + s(j)*h(j+1,i);
        h(j+1,i) = -s(j)*h(j,i) + c(j)*h(j+1,i);
        h(j,i) = tmp;
      }

      double tmp = sqrt(h(i,i)*h(i,i) + h(i+1,i)*h(i+1,i));
      c(i) = h(i,i) / tmp;
      s(i) = h(i+1,i) / tmp;
      h(i,i) = tmp;
      h(i+1,i) = 0.0;
      err(i+1) = -s(i)*err(i);
      err(i) = c(i)*err(i);
      #ifdef debug_gmres_v
      std::cout << msg_prefix << "err(i+1): " << err(i+1) << std::endl;
      std::cout << msg_prefix << "tmp: " << tmp << std::endl;
      #endif

      if (fabs(err(i+1)) < eps) {
        ls.suc = true;
        break;
      }
    } // for int i = 0; i < ls.sD

    if (last_i >= ls.sD) {
      last_i = ls.sD - 1;
    }

    for (int i = 0; i <= last_i; i++) {
      y(i) = err(i);
    }
    //y = err(1:i)

    for (int j = last_i; j >= 0; j--) { 
      for (int k = j+1; k <= last_i; k++) {
        y(j) = y(j) - h(j,k)*y(k);
      }
      y(j) = y(j) / h(j,j);
      //std::cout << msg_prefix << j+1 << " " << y(j) << std::endl;
    }

    for (int j = 0; j <= last_i; j++) {
      omp_la::omp_sum_v(dof, nNo, y(j), X, u.slice(j));
      //CALL OMPSUMV(dof, nNo, y(j), X, u(:,:,j))
    }

    ls.fNorm = fabs(err(last_i+1));
    if (ls.suc) {
      break;
    }
  }

  R = X;
  ls.callD = fsi_linear_solver::fsils_cpu_t() - ls.callD;
  ls.dB  = 10.0 * log(ls.fNorm / ls.dB);

  #ifdef debug_gmres_v
  double exec_time = fsi_linear_solver::fsils_cpu_t() - time;
  std::cout << msg_prefix << "Execution time: " << exec_time << std::endl;
  std::cout << msg_prefix << "Done" << std::endl;
  #endif
  double exec_time = fsi_linear_solver::fsils_cpu_t() - time;
  std::cout << msg_prefix << "Execution time: " << exec_time << std::endl;
  std::cout << msg_prefix << "Done" << std::endl;

  //Array<double>::write_disabled = false;
  //R.write(msg_prefix+"R");
  //exit(0);
}


};


