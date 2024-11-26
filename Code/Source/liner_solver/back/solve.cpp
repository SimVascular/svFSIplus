
// In this routine, the appropriate LS algorithm is called and
// the solution is returned.

#include "lhs.h"
#include "CmMod.h"
#include "bicgs.h"
#include "cgrad.h"
#include "gmres.h"
#include "ns_solver.h"
#include "precond.h"

namespace fsi_linear_solver {

//-------------
// fsils_solve
//-------------
//
// Modifies: Val, Ri
//
// Ri(dof,lhs.nNo): Residual
// Val(dof*dof,lhs.nnz): LHS
//
// Reproduces 'SUBROUTINE FSILS_SOLVE (lhs, ls, dof, Ri, Val, prec, incL, res)'.
//
void fsils_solve(FSILS_lhsType& lhs, FSILS_lsType& ls, const int dof, Array<double>& Ri, Array<double>& Val, 
    const consts::PreconditionerType prec, const Vector<int>& incL, const Vector<double>& res)
{
  using namespace consts;

  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[fsils_solve:") + std::to_string(tid) + "] ";
  #define n_debug_fsils_solve
  #ifdef debug_fsils_solve
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== fsils_solve ==========" << std::endl;
  #endif

  const int nNo = lhs.nNo;
  const int nnz = lhs.nnz;
  const int nFaces = lhs.nFaces;
  #ifdef debug_fsils_solve
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "nnz: " << nnz << std::endl;
  std::cout << msg_prefix << "nFaces: " << nFaces << std::endl;
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "ls.LS_type: " << ls.LS_type << std::endl;
  #endif

  if (lhs.nFaces != 0) {
    for (auto& face : lhs.face) { 
      face.incFlag = true;
    }

    if (incL.size() != 0) {
      #ifdef debug_fsils_solve
      std::cout << msg_prefix << "incL is present" << std::endl;
      #endif
      for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
        #ifdef debug_fsils_solve
        std::cout << msg_prefix << "incL[" << faIn << "]: " << incL[faIn] << std::endl;
        #endif
        if (incL(faIn) == 0) {
          lhs.face[faIn].incFlag = false;
        }
      }
    }

    bool flag = false;
    for (auto& face : lhs.face) { 
      if (face.bGrp == BcType::BC_TYPE_Neu) {
        flag = true;
        break;
      }
    }
    //flag = ANY(lhs%face%bGrp.EQ.BC_TYPE_Neu)
    //std::cout << msg_prefix << "flag: " << flag << std::endl;

    if (res.size() == 0 && flag) {
      //PRINT *, "FSILS: res is required for Neu surfaces"
      //STOP "FSILS: FATAL ERROR"
      throw std::runtime_error("[fsils_solve] res is required for Neu surfaces");
    }

    for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
      auto& face = lhs.face[faIn];
      face.coupledFlag = false;
      if (!face.incFlag) {
        continue;
      }
      bool flag = (face.bGrp == BcType::BC_TYPE_Neu);
      if (flag && res(faIn) != 0.0) {
        face.res = res(faIn);
        face.coupledFlag = true;
      }
    }
  }

  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "Copy Ri into R ..." << std::endl;
  Array<double> R(dof,nNo), Wr(dof,nNo), Wc(dof,nNo);

  for (int a = 0; a < nNo; a++) {
    //std::cout << msg_prefix << "a: " << a+1 << "  lhs.map(a): " << lhs.map(a) << std::endl;
    for (int i = 0; i < dof; i++) {
      R(i,lhs.map(a)) = Ri(i,a);
    }
    //R(:,lhs.map(a)) = Ri(:,a)
  }

  //Ri.write(msg_prefix+"Ri");
  //R.write(msg_prefix+"R");
  //Val.write(msg_prefix+"Val");
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);

  // Apply preconditioner.
  //
  // Modifies Val and R.
  //
  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "Compute preconditioner ..." << std::endl;

  if (prec == PreconditionerType::PREC_FSILS) {
    #ifdef debug_fsils_solve
    std::cout << msg_prefix << "preconditioner PREC_FSILS " << std::endl;
    #endif
    precond::precond_diag(lhs, lhs.rowPtr, lhs.colPtr, lhs.diagPtr, dof, Val, R, Wc);
    //CALL PRECONDDIAG(lhs, lhs.rowPtr, lhs.colPtr, lhs.diagPtr, dof, Val, R, Wc)
  } else if (prec == PreconditionerType::PREC_RCS) {
    #ifdef debug_fsils_solve
    std::cout << msg_prefix << "preconditioner PREC_RCS" << std::endl;
    #endif
    precond::precond_rcs(lhs, lhs.rowPtr, lhs.colPtr, lhs.diagPtr, dof, Val, R, Wr, Wc);
    //CALL PRECONDRCS(lhs, lhs.rowPtr, lhs.colPtr, lhs.diagPtr, dof, Val, R, Wr, Wc)
  } else {
    //PRINT *, "This linear solver and preconditioner combination is not supported."
  }

  //Wc.write(msg_prefix+"Wc");
  //R.write(msg_prefix+"R");
  //Val.write(msg_prefix+"Val");
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);

  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "Solve ..." << std::endl;

  // Solve for 'R'.
  //
  switch (ls.LS_type) {
    case LinearSolverType::LS_TYPE_NS:
      //std::cout << msg_prefix << "solve using LS_TYPE_NS " << std::endl;
      ns_solver::ns_solver(lhs, ls, dof, Val, R);
      //CALL NSSOLVER(lhs, ls, dof, Val, R)
    break;

    case LinearSolverType::LS_TYPE_GMRES:
      if (dof == 1) {
        auto Valv = Val.row(0);
        auto Rv = R.row(0);
        gmres::gmres_s(lhs, ls.RI, dof, Valv, Rv);
        Val.set_row(0,Valv);
        R.set_row(0,Rv);
        //CALL GMRESS(lhs, ls.RI, Val, R)
      } else {
        gmres::gmres_v(lhs, ls.RI, dof, Val, R);
        //CALL GMRESV(lhs, ls.RI, dof, Val, R)
      }
    break;

    case LinearSolverType::LS_TYPE_CG:
      if (dof == 1) {
         throw std::runtime_error("FSILS: CGRADS is not implemented");
        //CALL CGRADS(lhs, ls.RI, Val, R)
      } else {
        cgrad::cgrad_v(lhs, ls.RI, dof, Val, R);
        //CALL CGRADV(lhs, ls.RI, dof, Val, R)
      }
    break;

    case LinearSolverType::LS_TYPE_BICGS:
      if (dof == 1) {
         throw std::runtime_error("FSILS: BICGSS is not implemented");
        //CALL BICGSS(lhs, ls.RI, Val, R)
      } else {
        bicgs::bicgsv(lhs, ls.RI, dof, Val, R);
        //CALL BICGSV(lhs, ls.RI, dof, Val, R)
      }
    break;

    default:
      throw std::runtime_error("FSILS: LS_type not defined");
  }

  // Element-wise multiplication.
  //
  for (int i = 0; i < Wc.size(); i++) {
    R(i) = Wc(i) * R(i);
  }
  //R = Wc*R;

  for (int a = 0; a < nNo; a++) {
    for (int i = 0; i < R.num_rows(); i++) {
      Ri(i,a) = R(i,lhs.map(a));
    }
    //Ri(:,a) = R(:,lhs.map(a))
  }
}

};
