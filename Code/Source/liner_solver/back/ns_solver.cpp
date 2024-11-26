
// This routine is mainley intended for solving incompressible NS or
// FSI equations with a form of AU=R, in which A = [K D;-G L] and
// G = -D^t

#include "ns_solver.h"

#include "fsils_api.hpp"
#include "fils_struct.hpp"

#include "add_bc_mul.h"
#include "cgrad.h"
#include "dot.h"
#include "ge.h"
#include "gmres.h"
#include "norm.h"
#include "spar_mul.h"

#include "Array3.h"

#include <math.h>

namespace ns_solver {

//--------
// bc_pre
//--------
//
// Modifies: lhs.face[].nS
//
void bc_pre(fsi_linear_solver::FSILS_lhsType& lhs, const int nsd, const int dof, const int nNo, const int mynNo) 
{
  /*
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[bc_pre:") + std::to_string(tid) + "] ";
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== bc_pre ==========" << std::endl;
  */

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];

    if (face.coupledFlag) {
      if (face.sharedFlag) {
        Array<double> v(nsd,nNo);

        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            v(i,Ac) = face.valM(i,a);
          }
        }

        face.nS = pow(norm::fsi_ls_normv(nsd, mynNo, lhs.commu, v),2.0);

      } else {
        face.nS = 0.0;
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            face.nS = face.nS + pow(face.valM(i,a),2.0);
          }
        }
      }
    }
  }
  //std::cout << msg_prefix << "Done" << std::endl;
}

//--------
// depart
//--------
// Store sections of the 'Val' into separate arrays: 'Gt', 'mK', etc.
//
// Modifies: no globals
//
void depart(fsi_linear_solver::FSILS_lhsType& lhs, const int nsd, const int dof, const int nNo, const int nnz, 
    const Array<double>& Val, Array<double>& Gt, Array<double>& mK, Array<double>& mG, Array<double>& mD, Vector<double>& mL)
{
  /*
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[depart:") + std::to_string(tid) + "] ";
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== depart ==========" << std::endl;
  std::cout << msg_prefix << "nsd: " << nsd << std::endl;
  */
  Vector<double> tmp((nsd+1)*(nsd+1));

  if (nsd == 2) {
    for (int i = 0; i < nnz; i++) {
      auto tmp = Val.col(i);

      mK(0,i) = tmp(0);
      mK(1,i) = tmp(1);
      mK(2,i) = tmp(3);
      mK(3,i) = tmp(4);

      mG(0,i) = tmp(2);
      mG(1,i) = tmp(5);

      mD(0,i) = tmp(6);
      mD(1,i) = tmp(7);

      mL(i) = tmp(8);
    }

  } else if (nsd == 3) {

    for (int i = 0; i < nnz; i++) {
      auto tmp = Val.col(i);

      mK(0,i) = tmp(0);
      mK(1,i) = tmp(1);
      mK(2,i) = tmp(2);
      mK(3,i) = tmp(4);
      mK(4,i) = tmp(5);
      mK(5,i) = tmp(6);
      mK(6,i) = tmp(8);
      mK(7,i) = tmp(9);
      mK(8,i) = tmp(10);

      mG(0,i) = tmp(3);
      mG(1,i) = tmp(7);
      mG(2,i) = tmp(11);

      mD(0,i) = tmp(12);
      mD(1,i) = tmp(13);
      mD(2,i) = tmp(14);

      mL(i)   = tmp(15);
    }

  } else {
    //PRINT *, "FSILS: Not defined nsd for DEPART", nsd
    //STOP "FSILS: FATAL ERROR"
  }

  for (int i = 0; i < nNo; i++) {
    for (int j = lhs.rowPtr(0,i); j <= lhs.rowPtr(1,i); j++) {
      int k = lhs.colPtr(j);

      for (int l = lhs.rowPtr(0,k); l <= lhs.rowPtr(1,k); l++) {
        if (lhs.colPtr(l) == i) {
          for (int m = 0; m < Gt.num_rows(); m++) {
            Gt(m,l) = -mG(m,j);
          }
          //Gt(:,l) = -mG(:,j)
          break;
        }
      }
    }
  }
}

//-----------
// ns_solver
//-----------
// This routine is mainley intended for solving incompressible NS or
// FSI equations with a form of AU=R, in which A = [K D;-G L] and
// G = -D^t
//
// Ri (dof, lhs.nNo )
//
void ns_solver(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_lsType& ls, const int dof, const Array<double>& Val, Array<double>& Ri)
{
  using namespace consts;
  using namespace fsi_linear_solver;

  #define n_debug_ns_solver
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[ns_solver:") + std::to_string(tid) + "] ";
  #ifdef debug_ns_solver
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== ns_solver ==========" << std::endl;
  double time = fsi_linear_solver::fsils_cpu_t();
  #endif

  const int nNo = lhs.nNo;
  const int nnz = lhs.nnz;
  const int mynNo = lhs.mynNo;
  const int nsd = dof - 1;
  const int iB = ls.RI.mItr;
  const int nB = 2*iB;

  #ifdef debug_ns_solver
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "nsd: " << nsd << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "nnz: " << nnz << std::endl;
  std::cout << msg_prefix << "mynNo: " << mynNo << std::endl;
  std::cout << msg_prefix << "iB: " << iB << std::endl;
  std::cout << msg_prefix << "nB: " << nB << std::endl;
  Ri.write(msg_prefix+"Ri");
  #endif

  Vector<double> Rc(nNo), Rci(nNo), tmp(nB*nB+nB), tmpG(nB*nB+nB), B(nB), xB(nB), oldxB(nB); 
  Array<double> Rm(nsd,nNo), Rmi(nsd,nNo), A(nB,nB), P(nNo,iB), MP(nNo,nB);
  Array3<double> U(nsd,nNo,iB), MU(nsd,nNo,nB);

  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < Ri.num_cols(); j++) {
      Rmi(i,j) = Ri(i,j);
    }
  }
  //Rmi = Ri(1:nsd,:)

  for (int i = 0; i < Ri.num_cols(); i++) {
    Rci(i) = Ri(dof-1,i);
  }
  //Rci = Ri(dof,:)

  Rm = Rmi;
  Rc = Rci;

  double eps = sqrt(pow(norm::fsi_ls_normv(nsd,mynNo,lhs.commu,Rm),2.0) + 
                    pow(norm::fsi_ls_norms(mynNo,lhs.commu,Rc),2.0));

  #ifdef debug_ns_solver
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "eps (Rm/Rc): " << eps << std::endl;
  Rm.write(msg_prefix+"Rm");
  Rc.write(msg_prefix+"Rc");
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);
  #endif

  ls.RI.iNorm = eps;
  ls.RI.fNorm = eps*eps;

  // Calling duration 
  ls.CG.callD = 0.0; 
  ls.GM.callD = 0.0;
  ls.RI.callD = fsi_linear_solver::fsils_cpu_t(); 

  ls.CG.itr   = 0;
  ls.GM.itr   = 0;
  ls.RI.suc   = false;
  eps = std::max(ls.RI.absTol, ls.RI.relTol*eps);
  #ifdef debug_ns_solver
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "eps: " << eps << std::endl;
  std::cout << msg_prefix << "ls.RI.iNorm: " << ls.RI.iNorm << std::endl;
  std::cout << msg_prefix << "ls.RI.fNorm: " << ls.RI.fNorm << std::endl;
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);
  #endif

  Array<double> Gt(nsd,nnz), mK(nsd*nsd,nnz), mG(nsd,nnz), mD(nsd,nnz);
  Vector<double> mL(nnz); 

  // Store sections of the 'Val' array into separate arrays: 'Gt', 'mK', etc.
  //
  // Modfies: Gt, mK, mG, mD, and mL.
  //
  //std::cout << msg_prefix << std::endl;
  #ifdef debug_ns_solver
  std::cout << msg_prefix << "depart ... " << std::endl;
  #endif
  depart(lhs, nsd, dof, nNo, nnz, Val, Gt, mK, mG, mD, mL);
  //CALL DEPART

  #ifdef debug_ns_solver
  Gt.write(msg_prefix+"Gt");
  mD.write(msg_prefix+"mD");
  mG.write(msg_prefix+"mG");
  mK.write(msg_prefix+"mK");
  mL.write(msg_prefix+"mL");
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);
  #endif

  // Computes lhs.face[].nS for each face.
  //
  //std::cout << msg_prefix << std::endl;
  #ifdef debug_ns_solver
  std::cout << msg_prefix << "bc_pre ... " << std::endl;
  #endif
  bc_pre(lhs, nsd, dof, nNo, mynNo);
  //CALL BCPRE
  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    auto& face = lhs.face[faIn];
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "faIn: " << faIn << "  face.nS: " << face.nS << std::endl;
    #endif
   }

  #ifdef debug_ns_solver
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "Loop i on ls.RI.mItr ... " << std::endl;
  #endif
  int iBB{0}; 
  int i_count{0};
 
  // Note: iB and iBB appear to index into arrays.
  //
  for (int i = 0; i < ls.RI.mItr; i++) {
  //for (int i = 0; i < 1; i++) {
    #ifdef debug_ns_solver
    auto istr = "_" + std::to_string(i+1);
    std::cout << msg_prefix << std::endl;
    std::cout << msg_prefix << "---------- i " << i+1 << " ----------" << std::endl;
    #endif
   
    int iB = 2*i;
    iBB = 2*i + 1;
    //iB  = 2*i - 1
    //iBB = 2*i
    ls.RI.dB = ls.RI.fNorm;
    i_count = i;
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "iB: " << iB << std::endl;
    std::cout << msg_prefix << "iBB: " << iBB << std::endl;
    std::cout << msg_prefix << "ls.RI.fNorm: " << ls.RI.fNorm << std::endl;
    #endif

    // Solve for U = inv(mK) * Rm
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "Call gmres  U = K^-1*Rm ..." << std::endl;
    #endif
    auto U_slice = U.slice(i);
    gmres::gmres(lhs, ls.GM, nsd, mK, Rm, U_slice);
    U.set_slice(i, U_slice);
    //CALL GMRES(lhs, ls.GM, nsd, mK, Rm, U(:,:,i))
    #ifdef debug_ns_solver_write
    U.write(msg_prefix+"U_a"+istr);
    #endif

    // P = D*U
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "P = D * U ..." << std::endl;
    #endif
    auto P_col = P.col(i);
    spar_mul::fsils_spar_mul_vs(lhs, lhs.rowPtr, lhs.colPtr, nsd, mD, U.slice(i), P_col);
    P.set_col(i, P_col);
    //CALL FSILS_SPARMULVS(lhs, lhs.rowPtr, lhs.colPtr, nsd, mD, U(:,:,i), P(:,i))
    #ifdef debug_ns_solver_write
    P.write(msg_prefix+"P_a"+istr);
    #endif

    // P = Rc - P
    //
    P.set_col(i, Rc - P.col(i));
    //P(:,i) = Rc - P(:,i)
    #ifdef debug_ns_solver_write
    P.write(msg_prefix+"P_b"+istr);
    #endif

    // P = [L + G^t*G]^-1*P
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "P = inv[L + G^t*G] * P  ..." << std::endl;
    #endif

    P_col = P.col(i);
    #ifdef debug_ns_solver
    P_col.write(msg_prefix+"P_col"+istr);
    #endif
    cgrad::schur(lhs, ls.CG, nsd, Gt, mG, mL, P_col);
    P.set_col(i, P_col);
    //CALL CGRAD_SCHUR(lhs, ls.CG, nsd, Gt, mG, mL, P(:,i))
    #ifdef debug_ns_solver_write
    P.write(msg_prefix+"P_c"+istr);
    #endif
    //MPI_Barrier(lhs.commu.comm);
    //exit(0);

    // MU1 = G*P
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "MU1 = G * P  ..." << std::endl;
    std::cout << msg_prefix << "i: " << i+1 << std::endl;
    std::cout << msg_prefix << "iB: " << iB+1 << std::endl;
    #endif
    P_col = P.col(i);
    auto MU_iB = MU.slice(iB);
    spar_mul::fsils_spar_mul_sv(lhs, lhs.rowPtr, lhs.colPtr, nsd, mG, P_col, MU_iB);
    MU.set_slice(iB, MU_iB); 
    //CALL FSILS_SPARMULSV(lhs, lhs.rowPtr, lhs.colPtr, nsd, mG, P(:,i), MU(:,:,iB))
    #ifdef debug_ns_solver_write
    MU.write(msg_prefix+"MU_a"+istr);
    #endif
    //MPI_Barrier(lhs.commu.comm);
    //exit(0);

    // MU2 = Rm - G*P
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "MU2 = Rm - G*P ..." << std::endl;
    #endif
    MU.set_slice(iBB, Rm - MU_iB);
    //MU(:,:,iBB) = Rm - MU(:,:,iB)
    #ifdef debug_ns_solver_write
    MU.write(msg_prefix+"MU_b"+istr);
    #endif

    // U = inv(K) * [Rm - G*P]
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "U = inv(K) * [Rm - G*P]  ..." << std::endl;
    std::cout << msg_prefix << "iBB: " << iBB << std::endl;
    #endif
    lhs.debug_active = true;
    auto U_i = U.slice(i);
    gmres::gmres(lhs, ls.GM, nsd, mK, MU.slice(iBB), U_i);
    U.set_slice(i, U_i);
    // CALL GMRES(lhs, ls.GM, nsd, mK, MU(:,:,iBB), U(:,:,i))
    #ifdef debug_ns_solver_write
    U.write(msg_prefix+"U_b"+istr);
    #endif
    //MPI_Barrier(lhs.commu.comm);
    //exit(0);

    // MU2 = K*U
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "MU2 = K*U ..." << std::endl;
    #endif
    auto MU_iBB = MU.slice(iBB);
    spar_mul::fsils_spar_mul_vv(lhs, lhs.rowPtr, lhs.colPtr, nsd, mK, U.slice(i), MU_iBB);
    //CALL FSILS_SPARMULVV(lhs, lhs.rowPtr, lhs.colPtr, nsd, mK, U(:,:,i), MU(:,:,iBB))
    MU.set_slice(iBB, MU_iBB);
    #ifdef debug_ns_solver_write
    MU.write(msg_prefix+"MU_c"+istr);
    #endif

    add_bc_mul::add_bc_mul(lhs, BcopType::BCOP_TYPE_ADD, nsd, U.slice(i), MU_iBB);
    MU.set_slice(iBB, MU_iBB);
    //CALL ADDBCMUL(lhs, BCOP_TYPE_ADD, nsd, U(:,:,i), MU(:,:,iBB))
    #ifdef debug_ns_solver_write
    MU.write(msg_prefix+"MU_d"+istr);
    #endif

    // MP1 = L*P
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "MP1 = L*P ..." << std::endl;
    #endif
    auto MP_iB = MP.col(iB);
    spar_mul::fsils_spar_mul_ss(lhs, lhs.rowPtr, lhs.colPtr, mL, P.col(i), MP_iB);
    MP.set_col(iB, MP_iB);
    //CALL FSILS_SPARMULSS(lhs, lhs.rowPtr, lhs.colPtr, mL, P(:,i), MP(:,iB))
    #ifdef debug_ns_solver_write
    MP.write(msg_prefix+"MP_a"+istr);
    #endif

    // MP2 = D*U
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "MP2 = D*U ... " << std::endl;
    #endif
    auto MP_iBB = MP.col(iBB);
    spar_mul::fsils_spar_mul_vs(lhs, lhs.rowPtr, lhs.colPtr, nsd, mD, U.slice(i), MP_iBB);

    MP.set_col(iBB, MP_iBB);
    //CALL FSILS_SPARMULVS(lhs, lhs.rowPtr, lhs.colPtr, nsd, mD, U(:,:,i), MP(:,iBB))
    #ifdef debug_ns_solver_write
    MP.write(msg_prefix+"MP_b"+istr);
    #endif

    int c = 0;

    for (int k = iB; k <= iBB; k++) {
      for (int j = 0; j <= k; j++) {
        tmp(c) = dot::fsils_nc_dot_v(nsd, mynNo, MU.slice(j), MU.slice(k))  +  
                 dot::fsils_nc_dot_s(mynNo, MP.col(j), MP.col(k));
        //tmp(c) = FSILS_NCDOTV(nsd, mynNo, MU(:,:,j), MU(:,:,k))  + FSILS_NCDOTS(     mynNo, MP(:,j),   MP(:,k))
        #ifdef debug_ns_solver
        std::cout << msg_prefix << "tmp(" << c+1 << "): " << tmp(c) << std::endl;
        #endif
        c = c + 1;
      }

      tmp(c) = dot::fsils_nc_dot_v(nsd, mynNo, MU.slice(k), Rmi)  +  
               dot::fsils_nc_dot_s(mynNo, MP.col(k), Rci);
      //tmp(c) = FSILS_NCDOTV(nsd, mynNo, MU(:,:,k), Rmi) + FSILS_NCDOTS(     mynNo,  MP(:,k),  Rci)
      #ifdef debug_ns_solver
      std::cout << msg_prefix << "tmp(" << c+1 << "): " << tmp(c) << std::endl;
      #endif
      c = c + 1;
    }

    #ifdef debug_ns_solver
    std::cout << msg_prefix << "c: " << c << std::endl;
    std::cout << msg_prefix << "Call MPI_ALLREDUCE .." << std::endl;
    #endif
    if (lhs.commu.nTasks > 1) {
      MPI_Allreduce(tmp.data(), tmpG.data(), c, cm_mod::mpreal, MPI_SUM, lhs.commu.comm);
      //CALL MPI_ALLREDUCE(tmp, tmpG, c, mpreal, MPI_SUM, lhs.commu.comm, j)
      tmp = tmpG;
    }

    // Set arrays for Gauss elimination
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << std::endl;
    std::cout << msg_prefix << "Set arrays for Gauss elimination ... " << std::endl;
    #endif
    c = 0;

    for (int k = iB; k <= iBB; k++) {
      for (int j = 0; j <= k; j++) {
        A(j,k) = tmp(c);
        A(k,j) = tmp(c);
        c = c + 1;
      }

      B(k) = tmp(c);
      c  = c + 1;
    }
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "c: " << c << std::endl;
    #endif

    xB = B;

    // Perform Gauss elimination.
    //
    // [TODO:DaveP] careful with iBB, seems to be an index and a counter.
    //
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "Perform Gauss elimination ... " << std::endl;
    #endif
    if (ge::ge(nB, iBB+1, A, xB)) {
      oldxB = xB;

    } else { 
      if (lhs.commu.masF) {
        std::cout << "[FSILS] Singular matrix detected" << std::endl;
        //throw std::runtime_error("FSILS: Singular matrix detected");
        //PRINT *,"FSILS: Singular matrix detected"
      }

      xB = oldxB;

      // [TODO:DaveP] careful here, shift iB and iBB by -1?
      //
      if (i > 0) {
        iB = iB - 2;
        iBB = iBB - 2;
        #ifdef debug_ns_solver
        std::cout << msg_prefix << "New iB: " << iB << std::endl;
        std::cout << msg_prefix << "New iBB: " << iBB << std::endl;
        std::cout << msg_prefix << "New iB: " << iB << std::endl;
        #endif
      }
      break;
    }

    //std::cout << msg_prefix << "Compute sum ... "  << std::endl;
    double sum = 0.0;
    for (int i = 0; i <= iBB; i++) { 
      sum += xB(i) * B(i);
    }
    ls.RI.fNorm = pow(ls.RI.iNorm,2.0) - sum;
    //ls.RI.fNorm = ls.RI.iNorm**2.0 - SUM(xB(1:iBB)*B(1:iBB))
    #ifdef debug_ns_solver
    std::cout << msg_prefix << "sum: " << sum << std::endl;
    std::cout << msg_prefix << "ls.RI.fNorm: " << ls.RI.fNorm << std::endl;
    #endif

    //MPI_Barrier(lhs.commu.comm);
    //exit(0);

    if (ls.RI.fNorm < eps*eps) {
      ls.RI.suc = true;
      break;
    }

    Rm = Rmi - xB(0)*MU.slice(0);
    //Rm = Rmi - xB(1)*MU(:,:,1)

    Rc = Rci - xB(0)*MP.col(0);
    //Rc = Rci - xB(1)*MP(:,1)

    for (int j = 1; j <= iBB; j++) {
      Rm = Rm - xB(j)*MU.slice(j);
      //Rm = Rm - xB(j)*MU(:,:,j)

      Rc = Rc - xB(j)*MP.col(j);
      //Rc = Rc - xB(j)*MP(:,j)
    }

  //std::cout << msg_prefix << "Exit at i=" << i+1 << std::endl;
  //MPI_Barrier(lhs.commu.comm);
  //exit(0);
  } // for i = 0; i < ls.RI.mItr

  #ifdef debug_ns_solver
  std::cout << msg_prefix << "Finished iteration" << std::endl;
  std::cout << msg_prefix << std::endl;
  #endif

  // [TODO:DaveP] careful here with i_count.
  if (i_count >= ls.RI.mItr) {
  //if (i > ls.RI.mItr) {
    ls.RI.itr = ls.RI.mItr;
  } else { 
    ls.RI.itr = i_count;
    //ls.RI.itr = i;

    Rc = Rci - xB(0)*MP.col(0);
    //Rc = Rci - xB(1)*MP(:,1)

    for (int j = 1; j <= iBB; j++) {
      Rc = Rc - xB(j)*MP.col(j);
      //Rc = Rc - xB(j)*MP(:,j)
    }
  }

  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "i_count: " << i_count << std::endl;

  ls.Resc = static_cast<int>(100.0 * pow(norm::fsi_ls_norms(mynNo, lhs.commu, Rc),2.0) / ls.RI.fNorm);
  //ls.Resc = NINT(100.0*FSILS_NORMS(mynNo, lhs.commu, Rc)**2.0 / ls.RI.fNorm, KIND=LSIP)
  ls.Resm = 100 - ls.Resc;

  //std::cout << msg_prefix << std::endl;
  #ifdef debug_ns_solver
  std::cout << msg_prefix << "ls.Resc: " << ls.Resc << std::endl;
  std::cout << msg_prefix << "ls.Resm: " << ls.Resm << std::endl;
  std::cout << msg_prefix << "ls.RI.itr: " << ls.RI.itr << std::endl;
  #endif

  Rmi = xB(1) * U.slice(0);
  //Rmi = xB(2)*U(:,:,1)

  Rci = xB(0) * P.col(0);
  //Rci = xB(1)*P(:,1)

  // [TODO:DaveP] carfeful with ls.RI.itr.
  //
  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "Compute Rmi ... " << std::endl;

  for (int i = 1; i <= ls.RI.itr; i++) {
    //std::cout << msg_prefix << std::endl;
    //std::cout << msg_prefix << "i: " << i << std::endl;
    int iB = 2*i;
    int iBB = 2*i + 1;
    //int iB = 2*i - 1;
    //int iBB = 2*i;
    //std::cout << msg_prefix << "iB: " << iB << std::endl;
    //std::cout << msg_prefix << "iBB: " << iBB << std::endl;
    //std::cout << msg_prefix << "xB(iBB): " << xB(iBB) << std::endl;

    Rmi = Rmi + xB(iBB) * U.slice(i);
    //Rmi = Rmi + xB(iBB)*U(:,:,i)

    Rci = Rci + xB(iB) * P.col(i);
    //Rci = Rci + xB(iB)*P(:,i)
  }

  /*
  U.write(msg_prefix+"U");
  P.write(msg_prefix+"P");
  xB.write(msg_prefix+"xB");
  Rc.write(msg_prefix+"Rc");
  MP.write(msg_prefix+"MP");
  Rmi.write(msg_prefix+"Rmi");
  Rci.write(msg_prefix+"Rci");
  */

  // Set Calling duration.
  ls.RI.callD = fsi_linear_solver::fsils_cpu_t() - ls.RI.callD;
  //ls.RI.callD = FSILS_CPUT() - ls.RI.callD

  ls.RI.dB = 5.0 * log(ls.RI.fNorm / ls.RI.dB);
  //ls.RI.dB    = 5.0*LOG(ls.RI.fNorm/ls.RI.dB)

  if (ls.Resc < 0.0 || ls.Resm < 0.0) {
    std::cout << msg_prefix << "ls.Resc: " << ls.Resc << std::endl;
    std::cout << msg_prefix << "ls.Resm: " << ls.Resm << std::endl;
    ls.Resc = 0;
    ls.Resm = 0;
    ls.RI.dB = 0;
    ls.RI.fNorm = 0.0;

    if (lhs.commu.masF) {
      std::cout << "[WARNING:ns_solve] unexpected behavior in FSILS (likely due to the ill-conditioned LHS matrix)." << std::endl; 
      //throw std::runtime_error("FSILS: unexpected behavior in FSILS (likely due to the ill-conditioned LHS matrix)"); 
      //PRINT "(A)", "Warning: unexpected behavior in FSILS (likely due to the ill-conditioned LHS matrix)"
    }
  }

  ls.RI.fNorm = sqrt(ls.RI.fNorm);
  #ifdef debug_ns_solver
  std::cout << msg_prefix << "ls.RI.callD: " << ls.RI.callD << std::endl;
  std::cout << msg_prefix << "ls.RI.dB: " << ls.RI.dB << std::endl;
  std::cout << msg_prefix << "ls.RI.fNorm: " << ls.RI.fNorm << std::endl;
  #endif

  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < Rmi.num_cols(); j++) {
      Ri(i, j) = Rmi(i, j);
    }
  }
  //Ri(1:nsd,:) = Rmi

  Ri.set_row(dof-1, Rci);
  //Ri(dof,:) = Rci
  //Ri.write(msg_prefix+"Ri");

  if (lhs.commu.masF) {
    //CALL LOGFILE
  }

  #ifdef debug_ns_solver
  double exec_time = fsi_linear_solver::fsils_cpu_t() - time;
  std::cout << msg_prefix << "Execution time: " << exec_time << std::endl;
  std::cout << msg_prefix << "Done" << std::endl;
  #endif

  //MPI_Barrier(lhs.commu.comm);
  //exit(0);
}


};

