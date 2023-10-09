
// Subroutines related to initializing linear solver arrays and
// function calls to svFSILS and Trilinos solver library

#include "ls.h"

#include <math.h>

#include "consts.h"
#include "fsils_api.hpp"

#ifdef WITH_TRILINOS
#include "trilinos_linear_solver.h"
#endif

namespace ls_ns {

/// @brief Reproduces Fortran 'SUBROUTINE INIT_DIR_AND_COUPNEU_BC(incL, res)'.
//
void init_dir_and_coup_neu(ComMod& com_mod, const Vector<int>& incL,
                           const Vector<double>& res)
{
  using namespace consts;
  using namespace fsi_linear_solver;

  int dof = com_mod.dof;
  int gtnNo = com_mod.gtnNo;
  int tnNo = com_mod.tnNo;
  auto& lhs = com_mod.lhs;

  if (lhs.nFaces != 0) {
    for (auto& face : lhs.face) {
      face.incFlag = true;
    }

    for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
      if (incL(faIn) == 0) {
        lhs.face[faIn].incFlag = false;
      }
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

  auto& tls = com_mod.tls;
  tls.W = 1.0;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    if (!face.incFlag) {
      continue;
    }

    int faDof = std::min(face.dof, dof);

    if (face.bGrp == BcType::BC_TYPE_Dir) {
      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < faDof; i++) {
          tls.W(i, Ac) = tls.W(i, Ac) * face.val(i, a);
        }
      }
    }
  }

  Array<double> v(dof, tnNo);
  bool isCoupledBC = false;

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    if (face.coupledFlag) {
      isCoupledBC = true;
      int faDof = std::min(face.dof, dof);

      for (int a = 0; a < face.nNo; a++) {
        int Ac = face.glob(a);
        for (int i = 0; i < faDof; i++) {
          v(i, Ac) = v(i, Ac) + sqrt(fabs(res(faIn))) * face.val(i, a);
        }
      }
    }
  }

#ifdef WITH_TRILINOS
  trilinos_bc_create_(v.data(), isCoupledBC);
#endif
}

/// @brief Allocate com_mod.R and com_mod.Val arrays.
///
/// Modifies:
///    com_mod.R - Residual vector
///    com_mod.Val - LHS matrix
///
/// Reproduces 'SUBROUTINE LSALLOC(lEq)'.
//
void ls_alloc(ComMod& com_mod, eqType& lEq)
{
  int dof = com_mod.dof;
  int tnNo = com_mod.tnNo;
  int gtnNo = com_mod.gtnNo;
  auto& lhs = com_mod.lhs;

  com_mod.R.resize(dof, tnNo);

  if (!lEq.assmTLS) {
    com_mod.Val.resize(dof * dof, com_mod.lhs.nnz);
  }

#ifdef WITH_TRILINOS

  auto& tls = com_mod.tls;

  if (lEq.useTLS) {
    if (tls.W.size() != 0) {
      tls.W.clear();
      tls.R.clear();
      trilinos_lhs_free_();
    }

    tls.W.resize(dof, tnNo);
    tls.R.resize(dof, tnNo);

    int cpp_index = 1;
    int task_id = com_mod.cm.idcm();

    trilinos_lhs_create_(gtnNo, lhs.mynNo, tnNo, lhs.nnz, tls.ltg.data(),
                         com_mod.ltg.data(), com_mod.rowPtr.data(),
                         com_mod.colPtr.data(), dof, cpp_index, task_id);
  }

#endif
}

/// @brief Modifies:
///  com_mod.R      // Residual vector
///  com_mod.Val    // LHS matrix
///
/// Reproduces ' SUBROUTINE LSSOLVE(lEq, incL, res)'.
//
void ls_solve(ComMod& com_mod, eqType& lEq, const Vector<int>& incL,
              const Vector<double>& res)
{
#define n_debug_ls_solve
#ifdef debug_ls_solve
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lEq.sym: " << lEq.sym;
  dmsg << "lEq.useTLS: " << lEq.useTLS;
  dmsg << "lEq.assmTLS: " << lEq.assmTLS;
#endif

#ifdef WITH_TRILINOS

  if (lEq.useTLS) {
    init_dir_and_coup_neu(com_mod, incL, res);
  }

  if (lEq.assmTLS) {
    lEq.FSILS.RI.suc = false;
    auto& tls = com_mod.tls;
    int solver_type = static_cast<int>(lEq.ls.LS_type);
    int prec_type = static_cast<int>(lEq.ls.PREC_Type);

    trilinos_solve_(tls.R.data(), tls.W.data(), lEq.FSILS.RI.fNorm,
                    lEq.FSILS.RI.iNorm, lEq.FSILS.RI.itr, lEq.FSILS.RI.callD,
                    lEq.FSILS.RI.dB, lEq.FSILS.RI.suc, solver_type,
                    lEq.FSILS.RI.relTol, lEq.FSILS.RI.mItr, lEq.FSILS.RI.sD,
                    prec_type, lEq.assmTLS);

  } else if (lEq.useTLS) {
    auto& Val = com_mod.Val;
    auto& R = com_mod.R;
    auto& tls = com_mod.tls;
    int solver_type = static_cast<int>(lEq.ls.LS_type);
    int prec_type = static_cast<int>(lEq.ls.PREC_Type);

    trilinos_global_solve_(
        Val.data(), R.data(), tls.R.data(), tls.W.data(), lEq.FSILS.RI.fNorm,
        lEq.FSILS.RI.iNorm, lEq.FSILS.RI.itr, lEq.FSILS.RI.callD,
        lEq.FSILS.RI.dB, lEq.FSILS.RI.suc, solver_type, lEq.FSILS.RI.relTol,
        lEq.FSILS.RI.mItr, lEq.FSILS.RI.sD, prec_type);

  } else {
#endif

    auto& lhs = com_mod.lhs;
    int dof = com_mod.dof;
    auto& R = com_mod.R;      // Residual vector
    auto& Val = com_mod.Val;  // LHS matrix

    fsi_linear_solver::fsils_solve(lhs, lEq.FSILS, dof, R, Val,
                                   lEq.ls.PREC_Type, incL, res);

#ifdef WITH_TRILINOS
  }

  if (lEq.useTLS) {
    for (int a = 0; a < com_mod.tnNo; a++) {
      for (int i = 0; i < com_mod.R.nrows(); i++) {
        com_mod.R(i, a) = com_mod.tls.R(i, com_mod.lhs.map(a));
      }
    }
  }
#endif
}

};  // namespace ls_ns
