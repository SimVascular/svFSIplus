
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
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[fsils_bc_create:") + std::to_string(tid) + "] ";
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== fsils_bc_create ==========" << std::endl;
  std::cout << msg_prefix << "lhs.nNo: " << lhs.nNo << std::endl;
  std::cout << msg_prefix << "lhs.nFaces: " << lhs.nFaces << std::endl;
  std::cout << msg_prefix << "faIn: " << faIn << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "BC_type: " << BC_type << std::endl;
  std::cout << msg_prefix << "Val.size(): " << Val.size() << std::endl;
  std::cout << msg_prefix << "Val.nrows: " << Val.nrows_ << std::endl;
  std::cout << msg_prefix << "Val.ncols: " << Val.ncols_ << std::endl;
  #endif

  if (faIn >= lhs.nFaces) {
    throw std::runtime_error("FSILS: faIn is exceeding lhs structure maximum number of faces (" 
        + std::to_string(lhs.nFaces) + ") is less than " + std::to_string(faIn) + ".");
  }

  if (faIn <= -1) {
  //if (faIn <= 0) {
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
      for (int i = 0; i < Val.num_rows(); i++) {
        lhs.face[faIn].val(i,a) = Val(i,a);
        //std::cout << msg_prefix << "Val(i,a): " << a+1 << " " << i+1 << " " << Val(i,a) << std::endl;
      }
      //lhs.face[faIn].val(:,a) = Val(:,a);
    }
  } else { 
     lhs.face[faIn].val = 0.0;
  }

  if (lhs.commu.nTasks > 1) {
    int a = 0;

    //std::cout << msg_prefix << "lhs.face[faIn].nNo: " << lhs.face[faIn].nNo << std::endl;
    if (lhs.face[faIn].nNo != 0) {
      a = 1;
    }

    int Ac;

    MPI_Allreduce(&a, &Ac, 1, cm_mod::mpint, MPI_SUM, lhs.commu.comm);
    //CALL MPI_ALLREDUCE(a, Ac, 1, mpint, MPI_SUM, lhs.commu.comm, i)
    //std::cout << msg_prefix << "Ac: " << Ac << std::endl;

    if (Ac > 1) {
      lhs.face[faIn].sharedFlag = true;
      Array<double> v(dof,lhs.nNo);

      for (int a = 0; a < nNo; a++) {
        int Ac = lhs.face[faIn].glob(a);
        for (int i = 0; i < dof; i++) {
          v(i,Ac) = lhs.face[faIn].val(i,a);
        }
        //v(:,Ac) = lhs.face[faIn].val(:,a)
      }

      //std::cout << msg_prefix << "fsils_commuv ..." << std::endl;
      //std::cout << msg_prefix << "v.size(): " << v.size() << std::endl;
      fsils_commuv(lhs, dof, v); 
      // CALL FSILS_COMMUV(lhs, dof, v)

      for (int a = 0; a < nNo; a++) {
        int Ac = lhs.face[faIn].glob(a);
        for (int i = 0; i < dof; i++) {
          lhs.face[faIn].val(i,a) = v(i,Ac);
        }
        //lhs.face(faIn).val(:,a) = v(:,Ac)
      }
    }
  }
  //std::cout << msg_prefix << "Done" << std::endl;
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);
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

