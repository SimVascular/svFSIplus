
#include "fsils_api.hpp"
#include "fils_struct.hpp"

namespace fsi_linear_solver {

//-----------------
// fsils_bc_create
//-----------------
//
// Modifies:
//  lhs.face[faIn].nNo 
//  lhs.face[faIn].dof
//  lhs.face[faIn].bGrp 
//  lhs.face[faIn].glob
//  lhs.face[faIn].val
//  lhs.face[faIn].valM
//
void fsils_bc_create(FSILS_lhsType& lhs, int faIn, int nNo, int dof, BcType BC_type, const Vector<int>& gNodes, 
    const Array<double>& Val)
{
  using namespace consts;

  #define n_debug_fsils_bc_create
  #ifdef debug_fsils_bc_create
  DebugMsg dmsg(__func__,  lhs.commu.task);
  dmsg.banner();
  dmsg << "lhs.nNo: " << lhs.nNo;
  dmsg << "lhs.nFaces: " << lhs.nFaces;
  dmsg << "faIn: " << faIn;
  dmsg << "nNo: " << nNo;
  dmsg << "dof: " << dof;
  dmsg << "BC_type: " << BC_type;
  dmsg << "Val.size(): " << Val.size();
  dmsg << "Val.nrows: " << Val.nrows_;
  dmsg << "Val.ncols: " << Val.ncols_;
  #endif

  if (faIn >= lhs.nFaces) {
    throw std::runtime_error("FSILS: faIn is exceeding lhs structure maximum number of faces (" 
        + std::to_string(lhs.nFaces) + ") is less than " + std::to_string(faIn) + ".");
  }

  if (faIn <= -1) {
    throw std::runtime_error("FSILS: faIn is smaller than zero");
  }

  lhs.face[faIn].nNo = nNo;
  lhs.face[faIn].dof = dof;
  lhs.face[faIn].bGrp = BC_type;

  lhs.face[faIn].glob.resize(nNo); 
  lhs.face[faIn].val.resize(dof,nNo);   
  lhs.face[faIn].valM.resize(dof,nNo);

  for (int a = 0; a < nNo; a++) {
    int Ac = lhs.map(gNodes(a));
    lhs.face[faIn].glob(a) = Ac;
  }

  if (Val.size() != 0) {
    for (int a = 0; a < nNo; a++) {
      for (int i = 0; i < Val.nrows(); i++) {
        lhs.face[faIn].val(i,a) = Val(i,a);
      }
    }
  } else { 
     lhs.face[faIn].val = 0.0;
  }

  if (lhs.commu.nTasks > 1) {
    int a = 0;

    if (lhs.face[faIn].nNo != 0) {
      a = 1;
    }

    int Ac;

    MPI_Allreduce(&a, &Ac, 1, cm_mod::mpint, MPI_SUM, lhs.commu.comm);

    if (Ac > 1) {
      lhs.face[faIn].sharedFlag = true;
      Array<double> v(dof,lhs.nNo);

      for (int a = 0; a < nNo; a++) {
        int Ac = lhs.face[faIn].glob(a);
        for (int i = 0; i < dof; i++) {
          v(i,Ac) = lhs.face[faIn].val(i,a);
        }
      }

      fsils_commuv(lhs, dof, v); 

      for (int a = 0; a < nNo; a++) {
        int Ac = lhs.face[faIn].glob(a);
        for (int i = 0; i < dof; i++) {
          lhs.face[faIn].val(i,a) = v(i,Ac);
        }
      }
    }
  }
}

//-----------------
// fsils_bc_create
//-----------------
// fsils_bc_create() without optional 'Val' parameter.
//
void fsils_bc_create(FSILS_lhsType& lhs, int faIn, int nNo, int dof, BcType BC_type, const Vector<int>& gNodes)
{
  Array<double> Val;
  fsils_bc_create(lhs, faIn, nNo, dof, BC_type, gNodes, Val);
}

};

