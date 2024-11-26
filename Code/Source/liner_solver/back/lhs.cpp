
// The functions here reproduce the subroutines defined in svFSILS/LHS.f.

#include "lhs.h"
#include "CmMod.h"

#include "mpi.h"

namespace fsi_linear_solver {

//------------------
// fsils_lhs_create 
//------------------
//
// Modifies:
//
//  lhs.foC 
//  lhs.gnNo 
//  lhs.nNo 
//  lhs.nnz 
//  lhs.commu 
//  lhs.nFaces 
//  lhs.mynNo 
//
//  lhs.colPtr
//  lhs.rowPtr
//  lhs.diagPtr
//  lhs.map
//  lhs.face
//
void fsils_lhs_create(FSILS_lhsType& lhs, FSILS_commuType& commu, int gnNo, int nNo, int nnz, Vector<int>& gNodes,  
       Vector<int> &rowPtr, Vector<int>& colPtr, int nFaces)
{
  #define n_debug_fsils_lhs_create
  #ifdef debug_fsils_lhs_create
  auto msg_prefix = std::string("[fsils_lhs_create:") + std::to_string(commu.task) + "] ";
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== fsils_lhs_create ==========" << std::endl;
  #endif

  lhs.foC = true; 
  lhs.gnNo = gnNo;
  lhs.nNo = nNo;
  lhs.nnz = nnz;
  lhs.commu = commu;
  lhs.nFaces = nFaces;
  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << "gnNo: " << gnNo << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "nnz: " << nnz << std::endl;
  std::cout << msg_prefix << "nFaces: " << nFaces << std::endl;
  #endif

  int nTasks = commu.nTasks;
  auto comm = commu.comm;
  auto tF = commu.tF;
  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << "nTasks: " << nTasks << std::endl;
  std::cout << msg_prefix << "tF: " << tF << std::endl;
  #endif
  //rowPtr.write(msg_prefix+"rowPtr",1);

  lhs.colPtr.resize(nnz); 
  lhs.rowPtr.resize(2,nNo); 
  lhs.diagPtr.resize(nNo);
  lhs.map.resize(nNo); 

  // [TODO:DaveP] lhs.face[] is indexed by faIn which seems to be 1-based
  // so allocate lhs.face[] to be nFaces+1.
  lhs.face.resize(nFaces);
  //lhs.face.resize(nFaces+1);
  // ALLOCATE (lhs%colPtr(nnz), lhs%rowPtr(2,nNo), lhs%diagPtr(nNo), lhs%map(nNo), lhs%face(nFaces))

  // For a sequential simulation. 
  //
  if (nTasks == 1) {

    for (int i = 0; i < nnz; i++) {
      lhs.colPtr(i) = colPtr(i);
    }

    for (int Ac = 0; Ac < nNo; Ac++) {
      int s = rowPtr(Ac);
      int e = rowPtr(Ac+1) - 1;

      for (int i = s; i <= e; i++) {
        int a = colPtr(i);
        if (Ac == a) {
          lhs.diagPtr(Ac) = i;
          break; 
        }
      }

      lhs.rowPtr(0,Ac) = s;
      lhs.rowPtr(1,Ac) = e;
      lhs.map(Ac) = Ac;
    }

    /*
    std::cout << msg_prefix << "lhs.diagPtr: " << std::endl;
    for (int i = 0; i < nNo; i++) {
      std::cout << msg_prefix << i+1 << ": " << lhs.diagPtr(i) << std::endl;
    }

    std::cout << msg_prefix << "lhs.rowPtr: " << std::endl;
    for (int i = 0; i < nNo; i++) {
      std::cout << msg_prefix << i+1 << ": " << lhs.rowPtr(0,i) << " " << lhs.rowPtr(1,i) << std::endl;
    }
    */

    lhs.mynNo = nNo;
    return; 
  }

  // Get the number of nodes for this process and the max number of 
  // nodes for all processes.
  //
  // Combines values from all processes and distributes the result back to all processes.
  //
  int maxnNo;
  MPI_Allreduce(&nNo, &maxnNo, 1, cm_mod::mpint, MPI_MAX, comm);

  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "maxnNo: " << maxnNo << std::endl;
  #endif
  //ALLOCATE(aNodes(maxnNo,nTasks), part(maxnNo), sCount(nTasks), disp(nTasks), gtlPtr(gnNo), ltg(nNo))

  // Initialize some data structures.
  //
  Vector<int> part(maxnNo); 
  part = -1;

  for (int i = 0; i < nNo; i++) {
    part(i) = gNodes(i);
  }
  //part = 0
  //part(1:nNo) = gNodes

  Vector<int> sCount(nTasks); 
  Vector<int> disp(nTasks); 

  for (int i = 0; i < nTasks; i++) {
     disp(i) = i*maxnNo;
     //disp(i) = (i-1)*maxnNo
     sCount(i) = maxnNo;
  }

  Array<int> aNodes(maxnNo,nTasks); 
  Vector<int> ltg(nNo, "ltg");

  // Gather data from all processes and deliver it to all. Each process may contribute a different amount of data.
  MPI_Allgatherv(part.data(), maxnNo, cm_mod::mpint, aNodes.data(), sCount.data(), disp.data(), cm_mod::mpint, comm);
  // CALL MPI_ALLGATHERV(part, maxnNo, mpint, aNodes, sCount, disp, mpint, comm, ierr)

  Vector<int> gtlPtr(gnNo, "gtlPtr"); 
  gtlPtr = -1;
  //gtlPtr = 0

  for (int a = 0; a < nNo; a++) {
     int Ac = gNodes(a);
     gtlPtr(Ac) = a;
  }

  // Including the nodes shared by processors with higher ID at the end,
  // and including the nodes shared by lower processors IDs at the front.
  // shnNo is counter for lower ID and mynNo is counter for higher ID
  //
  // [TODO:DaveP] lhs.mynNo indexes into an array so must be 0-based?
  //
  lhs.mynNo = nNo;
  //lhs.mynNo = nNo - 1;
  lhs.shnNo = 0;
  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << "Including the nodes shared ..." << std::endl;
  std::cout << msg_prefix << "lhs.mynNo: " << lhs.mynNo << std::endl;
  std::cout << msg_prefix << "tF: " << tF << std::endl;
  #endif

  for (int i = nTasks-1; i >= 0; i--) { 
    //std::cout << msg_prefix << "---- i " << i << " ----" << std::endl;
    // Will include local nodes later
    if (i == tF) {
      continue; 
    }

    for (int a = 0; a < maxnNo; a++) { 
      //std::cout << msg_prefix << "  ---- a " << a << " ----" << std::endl;
      // Global node number in processor i at location a
      int Ac = aNodes(a,i);
      //std::cout << msg_prefix << "  Ac: " << Ac << std::endl;
      // Exit if this is the last node
      if (Ac == -1) {
        break; 
      }

      // Corresponding local node in current processor.
      int ai = gtlPtr(Ac);
      if (ai != -1) {
        // If this node has not been included already
        if (aNodes(ai,tF) != -1) {
          // If the processor ID is lower, it is appended to the beginning
          if (i < tF) {
            //std::cout << msg_prefix << "  beg Ac: " << Ac << std::endl;
            ltg(lhs.shnNo) = Ac;
            lhs.shnNo = lhs.shnNo + 1;
          // If the processor ID is higher, it is appended to the end
          } else {
            //std::cout << msg_prefix << "  end Ac: " << Ac << std::endl;
            ltg(lhs.mynNo-1) = Ac;
            lhs.mynNo = lhs.mynNo - 1;
          }
          aNodes(ai,tF) = -1;
        }
      }
    }
  }

  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "Including local nodes ..." << std::endl;
  std::cout << msg_prefix << "lhs.shnNo: " << lhs.shnNo << std::endl;
  std::cout << msg_prefix << "lhs.mynNo: " << lhs.mynNo << std::endl;
  #endif
  // Now including the local nodes that are left behind
  int j = lhs.shnNo + 1;
  for (int a = 0; a < nNo; a++) {
    int Ac = aNodes(a,tF);
    //std::cout << msg_prefix << "a: " << a+1 << " Ac: " << Ac+1 << std::endl;
    //  If this node has not been included already
    if (Ac != -1) {
      ltg(j-1) = Ac;
      j = j + 1;
      //std::cout << msg_prefix << "j: " << j << std::endl;
    }
  }

  if (j != lhs.mynNo+1) {
    //PRINT *, "FSILS: Unexpected behavior", j, lhs.mynNo
    throw std::runtime_error("FSILS: Unexpected behavior: j=" + std::to_string(j) + " lhs.mynNo: " + std::to_string(lhs.mynNo) + ".");
    MPI_Finalize();
    //STOP "FSILS: FATAL ERROR"
  }

  // Having the new ltg pointer, map is constructed
  //
  gtlPtr = -1;
  //gtlPtr = 0;
  for (int a = 0; a < nNo; a++) {
     int Ac = ltg(a);
     gtlPtr(Ac) = a;
  }

  for (int a = 0; a < nNo; a++) {
     int Ac = gNodes(a);
     lhs.map(a) = gtlPtr(Ac);
  }

  // Based on the new ordering of the nodes, rowPtr and colPtr are constructed
  //
  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "Construct rowPtr and colPtr ..." << std::endl;
  #endif
  for (int a = 0; a < nNo; a++) {
    int Ac = lhs.map(a);
    //std::cout << msg_prefix << a+1 << "   Ac: " << Ac+1 << " " << rowPtr(a)+1 << " " << rowPtr(a+1) << std::endl;
    lhs.rowPtr(0,Ac) = rowPtr(a);
    lhs.rowPtr(1,Ac) = rowPtr(a+1) - 1;
  }
  //lhs.rowPtr.write(msg_prefix+"rowPtr",true,1);

  for (int i = 0; i < nnz; i++) {
    lhs.colPtr(i) = lhs.map(colPtr(i));
  }

  // diagPtr points to the diagonal entries of LHS
  for (int Ac = 0; Ac < nNo; Ac++) {
    for (int i = lhs.rowPtr(0,Ac); i <= lhs.rowPtr(1,Ac); i++) {
      int a = lhs.colPtr(i);
      if (Ac == a) {
        lhs.diagPtr(Ac) = i;
        break; 
      }
    }
  }

  // Constructing the communication data structure based on the ltg
  //
  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "Constructing the communication data ..." << std::endl;
  #endif
  for (int i = 0; i < nNo; i++) {
    part(i) = ltg(i);
  }
  // part(1:nNo) = ltg

  MPI_Allgatherv(part.data(), maxnNo, cm_mod::mpint, aNodes.data(), sCount.data(), disp.data(), cm_mod::mpint, comm);

  // This variable keeps track of number of shared nodes
  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "Process shared nodes ..." << std::endl;
  #endif
  disp = 0;
  //disp = 0;
  lhs.nReq = 0;

  for (int i = 0; i < nTasks; i++) {
    //std::cout << msg_prefix << "----- i " << i << " -----" << std::endl;
    if (i == tF) {
      continue; 
    }
    for (int a = 0; a <  maxnNo; a++) {
      // Global node number in processor i at location a
      int Ac = aNodes(a,i);
      // Exit if this is the last node
      if (Ac == -1) {
        break;
      }
      // Corresponding local node in current processor
      int ai = gtlPtr(Ac);
      if (ai != -1) {
        disp(i) = disp(i) + 1;
      }
    }

    if (disp(i) != 0) {
      lhs.nReq = lhs.nReq + 1;
    }
  }
  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << "lhs.nReq: " << lhs.nReq << std::endl;
  #endif

  lhs.cS.resize(lhs.nReq);
  //ALLOCATE(lhs.cS(lhs.nReq))

  // Now that we know which processor is communicating to which, we can
  // setup the handles and structures
  //
  // lhs.cS[j].iP is the processor to communicate with.
  //
  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "Setup the handles ..." << std::endl;
  #endif
  j = 0;
  for (int i = 0; i < nTasks; i++) {
    int a = disp(i);
    if (a != 0) {
      lhs.cS[j].iP = i;
      lhs.cS[j].n = a;
      //std::cout << msg_prefix << "j: " << j << std::endl;
      //std::cout << msg_prefix << "a: " << a << std::endl;
      //std::cout << msg_prefix << "lhs.cS[j].iP: " << lhs.cS[j].iP << std::endl;
      //std::cout << msg_prefix << "lhs.cS[j].n: " << lhs.cS[j].n << std::endl;
      lhs.cS[j].ptr.resize(a);
      j = j + 1;
      //ALLOCATE(lhs.cS(j).ptr(a))
    }
  }

  // Order of nodes in ptr is based on the node order in processor
  // with higher ID. ptr is calculated for tF+1:nTasks and will be
  // sent over.
  //
  #ifdef debug_fsils_lhs_create
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "Order of nodes ..." << std::endl;
  #endif
  MPI_Status status;

  for (int i = 0; i < lhs.nReq; i++) {
    int iP = lhs.cS[i].iP;

    if (iP < tF) {
      //std::cout << msg_prefix << "MPI_Recv from iP: " << iP << std::endl;
      //std::cout << msg_prefix << "lhs.cS[i].n: " << lhs.cS[i].n << std::endl;
      MPI_Recv(lhs.cS[i].ptr.data(), lhs.cS[i].n, cm_mod::mpint, iP, 1,  comm, &status);
      // MPI_Recv(lhs.cS[i].ptr, lhs.cS[i].n, cm_mod::mpint, iP-1, 1,  comm);

      for (int j = 0; j < lhs.cS[i].n; j++) { 
        lhs.cS[i].ptr[j] = gtlPtr(lhs.cS[i].ptr[j]);
      }
    } else {
      // This is a counter for the shared nodes
      j = 0;
      for (int a = 0; a < maxnNo; a++) {
        // Global node number in processor i at location a
        int Ac = aNodes(a,iP);
        // Exit if this is the last node
        if (Ac == -1) {
          break; 
        }
        // Corresponding local node in current processor
        int ai = gtlPtr(Ac);
        if (ai != -1) {
          // Just for now global node ID is used. Later on this will be changed
          // to make sure nodes corresponds to each other on both processors
          // and then will be transformed to local node IDs
          lhs.cS[i].ptr[j] = Ac;
          j = j + 1;
        }
      }

      //std::cout << msg_prefix << "MPI_Send to iP: " << iP << std::endl;
      //std::cout << msg_prefix << "lhs.cS[i].n: " << lhs.cS[i].n << std::endl;
      MPI_Send(lhs.cS[i].ptr.data(), lhs.cS[i].n, cm_mod::mpint, iP, 1, comm);
      //MPI_Send(lhs.cS[i].ptr, lhs.cS[i].n, cm_mod::mpint, iP-1, 1, comm);

      for (int j = 0; j < lhs.cS[i].n; j++) {
        lhs.cS[i].ptr[j] = gtlPtr(lhs.cS[i].ptr[j]);
      }
    }
  }

  /*
  for (int i = 0; i < lhs.nReq; i++) {
    lhs.cS[i].ptr.write(msg_prefix+"cs_ptr_"+std::to_string(i+1),1);
  }
  */

  //std::cout << msg_prefix << "Done" << std::endl;
}


};


