
// The functions here reproduce the subroutines defined in svFSILS/INCOMMU.f.

// To syncronize the data on the boundaries between the processors.
//
// This a both way communication with three main part:
// 1 - rTmp {in master} = R          {from slave}
// 2 - R    {in master} = R + rTmp   {both from master}
// 3 - rTmp {in master} = R          {from master}
// 4 - R    {in slave}  = rTmp       {from master}

#include "fsils.hpp"
#include "CmMod.h"
#include "Array3.h"

#include "fsils_std.h"

namespace fsi_linear_solver {

//--------------
// fsils_commus
//--------------
//
void fsils_commus(const FSILS_lhsType& lhs, Vector<double>& R)
{
  if (lhs.commu.nTasks == 1) {
    return;
  }

  #define n_debug_fsils_commus
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[fsils_commus:") + std::to_string(tid) + "] ";
  #ifdef debug_fsils_commus
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== fsils_commus ==========" << std::endl;
  std::cout << msg_prefix << "lhs.cS.size(): " << lhs.cS.size() << std::endl;
  #endif

  // [TODO:DaveP] nReq can be 0 so lhs.cS.size() = 0, just return.
  //
  int nReq = lhs.nReq;
  if (lhs.cS.size() == 0) {
    return;
  }

  int nmax = std::max_element(lhs.cS.begin(), lhs.cS.end(), 
    [](const FSILS_cSType& a, const FSILS_cSType& b){return a.n < b.n;})->n; 
  //i = MAXVAL(lhs.cS.n)

  Array<double> sB(nmax,nReq); 
  Array<double> rB(nmax,nReq); 
  std::vector<MPI_Request> rReq(nReq); 
  std::vector<MPI_Request> sReq(nReq);
  //ALLOCATE(sB(i,nReq), rB(i,nReq), rReq(nReq), sReq(nReq))
  #ifdef debug_fsils_commus
  std::cout << msg_prefix << "nReq: " << nReq << std::endl;
  std::cout << msg_prefix << "i: " << nmax << std::endl;
  #endif

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      //sB(j,i) = j+1;
      sB(j,i) = R(k);
      //std::cout <<  msg_prefix << "i,j,k: " << i+1 << " " << j+1 << " " << k+1 << std::endl;
      //std::cout << "[fsils_commus] sB(j,i): " << j+1 << " " << i+1 << std::endl;
    }
  }

  //sB.set_debug(true, "sB " + msg_prefix);

  // [TODO:DaveP] possible problems here. 
  //
  int mpi_tag = 1;

  for (int i = 0; i < nReq; i++) {
    #ifdef debug_fsils_commus
    std::cout << msg_prefix << ">>> i: " << i << std::endl;
    std::cout << msg_prefix << "lhs.cS[i].n: " << lhs.cS[i].n << std::endl;
    std::cout << msg_prefix << "Data to recv/send lhs.cS[i].n: " << lhs.cS[i].n << std::endl;
    std::cout << msg_prefix << "Recv rB from lhs.cS[i].iP: " << lhs.cS[i].iP << std::endl;
    #endif
    auto rec_err = MPI_Irecv(rB.col_data(i), lhs.cS[i].n, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &rReq[i]);
    //CALL MPI_IRECV(rB(:,i), lhs%cS(i)%n, mpreal, lhs%cS(i)%iP-1, 1, lhs%commu%comm, rReq(i), ierr)

    #ifdef debug_fsils_commus
    std::cout << msg_prefix << "Send sB to lhs.cS[i].iP: " << lhs.cS[i].iP << std::endl;
    #endif
    auto send_err = MPI_Isend(sB.col_data(i), lhs.cS[i].n, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &sReq[i]);
    //CALL MPI_ISEND(sB(:,i), lhs%cS(i)%n, mpreal, lhs%cS(i)%iP-1, 1, lhs%commu%comm, sReq(i), ierr)
  }

  // Wait for the MPI receive to complete.
  //
  #ifdef debug_fsils_commus
  std::cout << msg_prefix << "Waiting for rB receive complete ..." << std::endl;
  #endif
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&rReq[i], &stat);
    #ifdef debug_fsils_commus
    std::cout << msg_prefix << "Done Waiting for rB receive" << std::endl;
    //std::cout << msg_prefix << "err: " << err << std::endl; 
    #endif
    //CALL MPI_Wait(rReq(i), stat, ierr)
  }

  #ifdef debug_fsils_commus
  std::cout << msg_prefix << "Set R ..." << std::endl;
  #endif
  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      R(k) = R(k) + rB(j,i);
    }
  }

  // Wait for the MPI send to complete.
  //
  #ifdef debug_fsils_commus
  std::cout << msg_prefix << "Waiting for sB send complete ..." << std::endl;
  //std::cout << std::resetiosflags( std::cout.flags() ); 
  #endif
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&sReq[i], &stat);
    #ifdef debug_fsils_commus
    std::cout << msg_prefix << "Done Waiting for sB send" << std::endl;
    //std::cout << msg_prefix << "err: " << err << std::endl; 
    #endif
    //CALL MPI_WAIT(sReq(i), stat, ierr)
  }

  /*
  if (lhs.debug_active) {  
    sB.write(msg_prefix+"sB");
    rB.write(msg_prefix+"rB");
  }
  */

  #ifdef debug_fsils_commus
  rB.write(msg_prefix+"rB");
  #endif

  #ifdef debug_fsils_commus
  std::cout << msg_prefix << "Done" << std::endl;
  #endif
}

//--------------
// fsils_commuv
//--------------
// This a both way communication with three main part:
//
// 1 - rTmp {in master} = R          {from slave}
// 2 - R    {in master} = R + rTmp   {both from master}
// 3 - rTmp {in master} = R          {from master}
// 4 - R    {in slave}  = rTmp       {from master}
//
void fsils_commuv(const FSILS_lhsType& lhs, int dof, Array<double>& R)
{
  if (lhs.commu.nTasks == 1) {
    return;
  }

  #define n_debug_fsils_commuv
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[fsils_commuv:") + std::to_string(tid) + "] ";
  #ifdef debug_fsils_commuv
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== fsils_commuv ==========" << std::endl;
  #endif

  // [TODO:DaveP] nReq can be 0 so lhs.cS.size() = 0, just return.
  if (lhs.cS.size() == 0) {
    return;
  }

  int nReq = lhs.nReq;
  int nmax = std::max_element(lhs.cS.begin(), lhs.cS.end(), 
      [](const FSILS_cSType& a, const FSILS_cSType& b){return a.n < b.n;})->n; 
  //i = MAXVAL(lhs.cS.n)

  #ifdef debug_fsils_commuv
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "nReq: " << nReq << std::endl;
  std::cout << msg_prefix << "nmax: " << nmax << std::endl;
  #endif

  Array3<double> sB(dof,nmax,nReq); 
  Array3<double> rB(dof,nmax,nReq); 
  std::vector<MPI_Request> rReq(nReq); 
  std::vector<MPI_Request> sReq(nReq);

  int n = 1;

  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      //std::cout << msg_prefix << "j: " << j+1 << " " << k+1 << std::endl;
      for (int l = 0; l < dof; l++) { 
        //sB(l,j,i) = n++;
        sB(l,j,i) = R(l,k);
      }
      // sB(:,j,i) = R(:,k)
    }
  }

  //int tid = lhs.commu.task;
  //auto msg_prefix = std::string("[fsils_commuv:") + std::to_string(tid) + "] ";
  //sB.set_debug(true, "sB " + msg_prefix);

  // [TODO:DaveP] possible problems here. 
  //
  int mpi_tag = 1;
  //int mpi_tag = 2;

  for (int i = 0; i < nReq; i++) {
    #ifdef debug_fsils_commuv
    std::cout << msg_prefix << std::endl;
    std::cout << msg_prefix << "lhs.cS[i].n: " << lhs.cS[i].n << std::endl;
    std::cout << msg_prefix << "Data to recv/send lhs.cS[i].n*dof: " << lhs.cS[i].n*dof << std::endl;
    std::cout << msg_prefix << "Receive rB slice from lhs.cS[i].iP: " << lhs.cS[i].iP << std::endl;
    #endif

    auto rec_err = MPI_Irecv(rB.slice_data(i), lhs.cS[i].n*dof, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &rReq[i]);
    //CALL MPI_IRECV(rB(:,:,i), lhs%cS(i)%n*dof, mpreal, lhs%cS(i)%iP-1, 1, lhs%commu%comm, rReq(i), ierr)

    #ifdef debug_fsils_commuv
    std::cout << msg_prefix << "Send sB slice to lhs.cS[i].iP: " << lhs.cS[i].iP << std::endl;
    #endif

    auto send_err = MPI_Isend(sB.slice_data(i), lhs.cS[i].n*dof, mpreal, lhs.cS[i].iP, mpi_tag, lhs.commu.comm, &sReq[i]);
    //CALL MPI_ISEND(sB(:,:,i), lhs%cS(i)%n*dof, mpreal, lhs%cS(i)%iP-1, 1, lhs%commu%comm, sReq(i), ierr)
  }

  // Wait for the MPI receive to complete.
  //
  #ifdef debug_fsils_commuv
  std::cout << msg_prefix << "Waiting for rB receive complete ..." << std::endl;
  #endif
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&rReq[i], &stat);
    //std::cout << msg_prefix << "err: " << err << std::endl;
    #ifdef debug_fsils_commuv
    std::cout << msg_prefix << "Done Waiting for rB receive" << std::endl;
    #endif
    //CALL MPI_Wait(rReq(i), stat, ierr)
  }

  #ifdef debug_fsils_commuv
  std::cout << msg_prefix << "Set R ..." << std::endl;
  #endif
  for (int i = 0; i < nReq; i++) {
    for (int j = 0; j < lhs.cS[i].n; j++) { 
      int k = lhs.cS[i].ptr(j);
      for (int l = 0; l < dof; l++) { 
        R(l,k) = R(l,k) + rB(l,j,i);
      }
      // R(:,k) = R(:,k) + rB(:,j,i)
    }
  }

  // Wait for the MPI send to complete.
  //
  #ifdef debug_fsils_commuv
  std::cout << msg_prefix << "Waiting for sB send complete ..." << std::endl;
  #endif
  for (int i = 0; i < nReq; i++) {
    MPI_Status stat;
    auto err = MPI_Wait(&sReq[i], &stat);
    #ifdef debug_fsils_commuv
    std::cout << msg_prefix << "Done Waiting for sB send" << std::endl;
    //std::cout << msg_prefix << "err: " << err << std::endl;
    #endif
    //CALL MPI_Wait(sReq(i), stat, ierr)
  }

  #ifdef debug_fsils_commuv
  if (lhs.debug_active) {  
    sB.write(msg_prefix+"sB");
    rB.write(msg_prefix+"rB");
  }
  #endif

  //sB.write(msg_prefix+"sB");
  //rB.write(msg_prefix+"rB");
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);

  #ifdef debug_fsils_commuv
  std::cout << msg_prefix << "Done" << std::endl;
  #endif
}


};


