
#include "remesh.h"

#include "all_fun.h"
#include "mat_fun.h"
#include "nn.h"
#include "output.h"
#include "read_msh.h"
#include "remeshTet.h"
#include "vtk_xml.h"

#include <array>
#include<iostream>
#include <filesystem>
#include<fstream>

namespace remesh {

//--------
// distre
//--------
// Distribute the new mesh elements to all processors
//
void distre(ComMod& com_mod, CmMod& cm_mod, mshType& lM, int& nEl, Vector<int>& gE)
{
  auto& cm = com_mod.cm;
  int gnEl = lM.gnEl;
  int eNoN = lM.eNoN;

  #define debug_distre
  #ifdef debug_distre
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  dmsg << "eNoN: " << eNoN;
  #endif

  Vector<int> part(gnEl);
  nEl = 0;

  for (int e = 0; e < gnEl; e++) {
    int i = 0;
    for (int a = 0; a < eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      if (lM.gpN(Ac) == cm.tF(cm_mod)) {
        i = i + 1;
      }
    }

    if (i == eNoN) { 
      nEl = nEl + 1;
      part(e) = 1;
    }
  }

  //DEALLOCATE(lM.gpN)

  //std::cout << "[distrn] ----- gE -----" << std::endl;
  gE.resize(nEl);
  nEl = 0;

  for (int e = 0; e < gnEl; e++) {
    if (part(e) > 0) {
      gE(nEl) = e;
      //std::cout << "[distrn] " << nEl+1 << "  " << gE(nEl)+1 << std::endl;
      nEl = nEl + 1;
    }
  }

  dmsg << "Done " << "";
  //exit(0);
}

//--------------
// dist_msh_srf
//--------------
// Reproduces Fortran 'SUBROUTINE DISTMSHSRF(lFa, lM, iOpt)'
//
void dist_msh_srf(ComMod& com_mod, ChnlMod& chnl_mod, faceType& lFa, mshType& lM, const int iOpt)
{
  const int nsd = com_mod.nsd;
  auto& rmsh = com_mod.rmsh;
  std::string sTmp = chnl_mod.appPath + ".remesh_tmp.dir";
  /*
  sTmp = TRIM(appPath)//".remesh_tmp.dir"
  INQUIRE(FILE=TRIM(sTmp)//"/.", EXIST=flag)
  if (.NOT. flag) {
     CALL SYSTEM("mkdir  -p  "//TRIM(sTmp))
  }
  */

  for (int e = 0; e < lFa.nEl; e++) {
    for (int a = 0; a < lFa.eNoN; a++) {
      int Ac = lFa.IEN(a,e);
      Ac = lFa.gN(Ac);
      lFa.IEN(a,e) = Ac;
    }
  }

  for (int iFa = 0; iFa < lM.nFa; iFa++) {
    lM.fa[iFa].nEl = lM.fa[iFa].gnEl;
    //if (ALLOCATED(lM.fa(iFa).IEN)) DEALLOCATE(lM.fa(iFa).IEN)
    //if (ALLOCATED(lM.fa(iFa).gE)) DEALLOCATE(lM.fa(iFa).gE)

    lM.fa[iFa].IEN.resize(lM.fa[iFa].eNoN,lM.fa[iFa].nEl);
    //ALLOCATE(lM.fa(iFa).IEN(lM.fa(iFa).eNoN,lM.fa(iFa).nEl))

    lM.fa[iFa].gE.resize(lM.fa[iFa].nEl);
    //ALLOCATE(lM.fa(iFa).gE(lM.fa(iFa).nEl))

    int eoff = 0;
    if (iFa > 0) {
      for (int i = 0; i < iFa-1; i++) {
        eoff = lM.fa[i].gnEl;
      }
      //eoff = SUM(lM.fa(1:iFa-1).gnEl)
    }

    for (int e = 0; e < lM.fa[iFa].gnEl; e++) {
      lM.fa[iFa].gE(e) = lFa.gE(eoff+e);
      for (int i = 0; i < lM.fa[iFa].eNoN; i++) {
        lM.fa[iFa].IEN(i,e) = lFa.IEN(i,eoff+e);
      }
      //lM.fa(iFa).IEN(:,e) = lFa.IEN(:,eoff+e);
    }

    read_msh_ns::calc_nbc(lM, lM.fa[iFa]);
    //CALL CALCNBC(lM, lM.fa(iFa))

    lM.fa[iFa].x.resize(nsd, lM.fa[iFa].nNo);
    //ALLOCATE(lM.fa(iFa).x(nsd, lM.fa(iFa).nNo))

    for (int a = 0; a < lM.fa[iFa].nNo; a++) {
      int Ac = lM.fa[iFa].gN(a);
      for (int i = 0; i < nsd; i++) {
        lM.fa[iFa].x(i,a) = lM.x(i,Ac);
      }
      //lM.fa(iFa).x(:,a) = lM.x(:,Ac)
    }

    for (int i = 0; i < lM.fa[iFa].gebc.ncols(); i++) {
      lM.fa[iFa].gebc(0,i) = lM.fa[iFa].gE(i);
    }
    //lM.fa(iFa).gebc(1,:) = lM.fa(iFa).gE(:)

    // [NOTE] I'm not sure about this.
    for (int i = 0; i < lM.fa[iFa].eNoN; i++) {
      for (int j = 0; i < lM.fa[iFa].nEl; i++) {
        lM.fa[iFa].gebc(i+1,j) = lM.fa[iFa].IEN(i,j);
      }
    }
    //lM.fa(iFa).gebc(2:1+lM.fa(iFa).eNoN,:) = lM.fa(iFa).IEN(:,:)

    if (iOpt == 1) {
      std::string fTmp = sTmp + lM.fa[iFa].name + "_" + std::to_string(rmsh.rTS) + "_cpp.vtp";
      //fTmp = TRIM(sTmp)//"/"//TRIM(lM.fa(iFa).name)//"_"// STR(rmsh.rTS)//".vtp"
      Vector<int> incNd(lM.gnNo);
      //ALLOCATE(incNd(lM.gnNo))

      for (int a = 0; a < lM.fa[iFa].nNo; a++) {
        int Ac = lM.fa[iFa].gN(a);
        incNd(Ac) = a;
      }

      for (int e = 0; e < lM.fa[iFa].nEl; e++) {
        for (int a = 0; lM.fa[iFa].eNoN; a++) {
          int Ac = lM.fa[iFa].IEN(a,e);
          lM.fa[iFa].IEN(a,e) = incNd(Ac) - 1;
        }
      }

      vtk_xml::write_vtp(com_mod, lM.fa[iFa], fTmp);
      //CALL WRITEVTP(lM.fa(iFa), fTmp)

      for (int e = 0; e < lM.fa[iFa].nEl; e++) {
        for (int a = 0; lM.fa[iFa].eNoN; a++) {
          int Ac = lM.fa[iFa].IEN(a,e) + 1;
          Ac = lM.fa[iFa].gN(Ac);
          lM.fa[iFa].IEN(a,e) = Ac;
        }
      }

      //DEALLOCATE(incNd)
    }
  }

}

//--------
// distrn
//--------
//
// Modifies
//   gN(nNo) - list of node indices 0  
//   lM.gpN - processor ID (1, 2, 3, ...) for each node
//
void distrn(ComMod& com_mod, CmMod& cm_mod, const int iM, mshType& lM, Array<double>& Dg, int& nNo, Vector<int>& gN)
{
  const int nsd = com_mod.nsd;
  auto& rmsh = com_mod.rmsh;
  auto& cm = com_mod.cm;

  #define debug_distrn
  #ifdef debug_distrn
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  int gnNo = lM.gnNo;
  int i = 0;
  double f = 2.5e-2;
  nNo = 0;

  Vector<int> part(gnNo), tmpI(gnNo);

  dmsg << "iM: " << iM;
  dmsg << "lM.gnNo: " << lM.gnNo;
  dmsg << "lM.nNo: " << lM.nNo;
  dmsg << "gnNo: " << gnNo;
  dmsg << "msh(iM).nNo: " << com_mod.msh[iM].nNo;
  dmsg << "cm.tF(cm_mod): " << cm.tF(cm_mod);

  while (true) {
    part = 0;
    tmpI = 0;
    nNo = 0;
    f = 2.0 * f;
    double tol = (1.0 + f) * rmsh.maxEdgeSize(iM);
    i = i+1;
    dmsg << "---------- pass " << std::to_string(i) + " ----------";
    dmsg << "tol: " << tol;

    for (int a = 0; a < gnNo; a++) {
      //dmsg << "---------- a " << std::to_string(a) + " ----------";
      if (part(a) != 0) {
        continue; 
      }
      double minS = std::numeric_limits<double>::max();

      for (int b = 0; b < com_mod.msh[iM].nNo; b++) {
        int Ac = com_mod.msh[iM].gN(b);
        double dS = 0.0;

        for (int i = 0; i < nsd; i++) {
          double diff = com_mod.x(i,Ac) + Dg(i,Ac) - lM.x(i,a);
          dS += diff * diff;
        }
        dS = sqrt(dS);
        //dS = SQRT(SUM((x(:,Ac)+Dg(:,Ac)-lM.x(:,a))**2._RKIND))

        if (minS > dS) minS = dS;
      }

      if (minS < tol) {
        nNo = nNo + 1;
        part(a) = cm.tF(cm_mod);
      }
    }
    dmsg << "nNo: " << nNo;

    MPI_Allreduce(part.data(), tmpI.data(), gnNo, cm_mod::mpint, MPI_MAX, cm.com());

    int b = 0;

    for (int a = 0; a < gnNo; a++) {
      if (tmpI(a) > 0) b = b + 1;
    }
    dmsg << "b: " << b;

    if (b == gnNo) {
      break;
    } else { 
      std::cout << "[distrn] Found only " + std::to_string(b) + " nodes in pass " + std::to_string(i) + 
        " out of " + std::to_string(gnNo) + " nodes." << std::endl;

      if (i > 5) {
        throw std::runtime_error("[distrn] Could not distribute all nodes in " + std::to_string(i) + " passes.");
      }
    }
  }

  gN.resize(nNo);
  lM.gpN.resize(gnNo);
  nNo = 0;

  //std::cout << "[distrn] ----- gN -----" << std::endl;

  for (int a = 0; a < gnNo; a++) {
    if (part(a) == cm.tF(cm_mod)) {
      gN(nNo) = a;
      //std::cout << "[distrn] " << nNo+1 << "  " << gN(nNo) << std::endl;
      nNo = nNo + 1;
    }
  }

  lM.gpN = part;
  dmsg << "Done" << " ";
}

//--------
// find_n
//--------
// Find element in the old mesh for each node in the new mesh
//
// Reproduces Fortran 'SUBROUTINE FINDN(xp, iM, Dg, eList, Ec, Nsf)'
//
void find_n(ComMod& com_mod, const Vector<double>& Xp, const int iM, const Array<double>& Dg, 
    const Vector<int>& eList, int& Ec, Vector<double>& Nsf, int pcount)
{
  const int nsd = com_mod.nsd;
  auto& msh = com_mod.msh[iM];
  int ne = eList.size();
  Array<double> Amat(nsd+1,msh.eNoN);

  Nsf = 0.0;
  bool debug = false;

#ifdef debug_find_c
  if (pcount == 19195+1) {
    std::cout << "[find_n] ========== find_n ==========  " << std::endl;
    std::cout << "[find_n] ne: " << ne << std::endl;
    std::cout << "[find_n] Xp: " << Xp << std::endl;
    debug = true;
  }
#endif

  for (int e = 0; e < ne; e++) {
    Ec = eList(e);

    if (Ec == -1) {
      break;
    }

#ifdef debug_find_c
    if (pcount == 19195+1) {
      std::cout << "[find_n] ---------- e " << e+1 << " ----------" << std::endl;
      std::cout << "[find_n] Ec: " << Ec+1 << std::endl;
    }
#endif

    Amat = 1.0;

    for (int a = 0; a <  msh.eNoN; a++) {
      int Ac = msh.IEN(a,Ec);
#ifdef debug_find_c
      if (pcount == 19195+1) {
        std::cout << "[find_n] Ac: " << Ac+1 << std::endl;
        std::cout << "[find_n] x: " << com_mod.x.col(Ac) << std::endl;
        std::cout << "[find_n] Dg: " << Dg.col(Ac) << std::endl;
      }
#endif
      for (int i = 0; i < nsd; i++) {
        Amat(i,a) = com_mod.x(i,Ac) + Dg(i,Ac);
      }
    }

#ifdef debug_find_c
    if (pcount == 19195+1) {
      std::cout << "[find_n]  " << std::endl;
      std::cout << "[find_n] Amat: " << Amat << std::endl;
    }
#endif

    // [NOTE] Sometimes the inverse computed by mat_inv() is not 
    // correct, A * Amt != I.  mat_inv_lp() seems to work all
    // the time.
    //
    //Amat = mat_fun::mat_inv(Amat, msh.eNoN, debug);
    Amat = mat_fun::mat_inv_lp(Amat, msh.eNoN);

    int a = 0;

#ifdef debug_find_c
    if (pcount == 19195+1) {
      std::cout << "[find_n]  " << std::endl;
      std::cout << "[find_n] inv Amat: " << Amat << std::endl;
      //if (e == 69) exit(0);
    }
#endif

    for (int i = 0; i < nsd+1; i++) {
      Nsf(i) = 0.0;

      for (int j = 0; j < msh.eNoN; j++) {
        Nsf(i) = Nsf(i) + Amat(i,j)*Xp(j);
      }

#ifdef debug_find_c
      if (pcount == 19195+1) {
        std::cout << "[find_n] i Nsf(i): " << i+1 << "  " << Nsf(i) << std::endl;
      }
#endif

      if (Nsf(i) > -1.0e-14 && Nsf(i) < (1.0+1.0e-14)) {
        a = a + 1;
#ifdef debug_find_c
        if (pcount == 19195+1) {
          std::cout << "[find_n] a: " << a << std::endl;
        }
#endif
      }
    }

    if (a == nsd+1) {
      return;
    }
  }

  Ec = -1;
  Nsf = 0.0;
}

//--------------
// get_adj_esrc
//--------------
// Create list of connected/adjacent elements for old/source mesh
//
// Reproduces Fortran 'SUBROUTINE GETADJESRC(lM, kneList)'
//
void get_adj_esrc(ComMod& com_mod, mshType& lM, Array<int>& kneList)
{
  #define debug_get_adj_esrc
  #ifdef debug_get_adj_esrc
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& rmsh = com_mod.rmsh;

  // First get elements around all nodes
  //
  Vector<int> nL(lM.nNo);

  for (int e = 0; e < lM.nEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.IEN(a,e);
      Ac = lM.lN(Ac);
      nL(Ac) = nL(Ac) + 1;
    }
  }

  int maxKNE = std::max(nL.max(), 0);
  dmsg << "maxKNE: " << maxKNE;

  Array<int> tmpList(maxKNE, lM.nNo);
  tmpList = -1;  // Stores element IDs
  nL = 0;        // Stores element counts

  for (int e = 0; e < lM.nEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.IEN(a,e);
      Ac = lM.lN(Ac);
      tmpList(nL(Ac), Ac) = e;
      nL(Ac) = nL(Ac) + 1;
    }
  }

  int b = 2 * maxKNE;

  // Now get elements around each element and avoid redundancy
  //
  get_elems:
    b = b + maxKNE;
    nL.resize(lM.nEl);
    kneList.resize(b, lM.nEl);
    kneList = -1;
    bool flag = false;

    for (int e = 0; e < lM.nEl; e++) {
      for (int a = 0; a < lM.eNoN; a++) {
        int Ac = lM.IEN(a,e);
        Ac = lM.lN(Ac);

        for (int i = 0; i < maxKNE; i++) {
          if (tmpList(i,Ac) != -1) {
            flag = true;

            for (int j = 0; j < nL(e); j++) {
              if (kneList(j,e) == tmpList(i,Ac)) {
                flag = false;
                break;
              }
            }

            if (flag) {
              nL(e) = nL(e) + 1;
              if (nL(e) >= b) goto get_elems;
              kneList(nL(e)-1, e) = tmpList(i,Ac);
            }
          } else {
            break;
          }

        } // for i

      } // for a

    } // for e

  maxKNE = std::max(nL.max(), 0);
  tmpList.resize(b,lM.nEl);
  tmpList = kneList;
  kneList.resize(maxKNE, lM.nEl);
  kneList = -1;

  for (int e = 0; e < lM.nEl; e++) {
    //std::cout << "[get_adj_esrc] " << e+1 << std::endl;
    for (int i = 0; i < maxKNE; i++) {
      kneList(i, e) = tmpList(i, e);
      //std::cout << kneList(i, e) << " ";
    }
    //std::cout << std::endl;
  }

  dmsg << "Done " << "";
}

//--------------
// get_adj_ntgt
//--------------
// Create list of connected/adjacent nodes for new/target mesh
//
// Reproduces Fortran 'SUBROUTINE GETADJNTGT(lM, nNo, nEl, gN, gE, knnList)'
//
void get_adj_ntgt(ComMod& com_mod, mshType& lM, const int nNo, const int nEl, const Vector<int>& gN, 
    const Vector<int>& gE, Array<int>& knnList)
{
  Vector<int> lN(lM.gnNo);

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    lN(Ac) = a;
  }

  Vector<int> nL(nNo);

  for (int e = 0; e < nEl; e++) {
    int Ec = gE(e);

    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,Ec);
      int i = lN(Ac);

      for (int b = 0; b < lM.eNoN; b++) {
        if (b == a) {
          continue;
        }
        int Bc = lM.gIEN(b,Ec);
        nL(i) = nL(i) + 1;
      }
    }
  }

  int maxKNN = std::max(nL.max(), 0);

  Array<int> tmpList(maxKNN, nNo);

  tmpList = -1;
  nL = 0;

  for (int e = 0; e < nEl; e++) {
    int Ec = gE(e);

    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,Ec);
      int i = lN(Ac);

      for (int b = 0; b < lM.eNoN; b++) {
        if (a == b) {
          continue;
        }

        int Bc = lM.gIEN(b,Ec);
        int j = lN(Bc);
        bool flag = true;

        for (int n = 0; n < nL(i); n++) {
          if (tmpList(n,i) ==  j) {
            flag = false;
            break;
          }
        }

        if (flag) {
          tmpList(nL(i), i) = j;
          nL(i) = nL(i) + 1;
        }
      }
    }
  }

  maxKNN = std::max(nL.max(), 0);
  knnList.resize(maxKNN, nNo);
  knnList = -1;

  for (int a = 0; a < nNo; a++) {
    //std::cout << "[get_adj_ntgt] " << a+1 << std::endl;
    for (int i = 0; i < maxKNN; i++) {
      knnList(i,a) = tmpList(i,a);
      //std::cout << knnList(i,a)+1 << " ";
    }
    //std::cout << std::endl;
  }

}

//---------
// interp
//---------
// Interpolation of data variables from source mesh to target mesh
//
void interp(ComMod& com_mod, CmMod& cm_mod, const int lDof, const int iM, mshType& tMsh, Array<double>& sD, Array<double>& tgD)
{
  const int nsd = com_mod.nsd;
  const int tnNo = com_mod.tnNo;
  auto& msh = com_mod.msh;
  auto& rmsh = com_mod.rmsh;
  auto& cm = com_mod.cm;
  #define debug_interp
  #ifdef debug_interp
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  if (nsd+1 != msh[iM].eNoN) {
    throw std::runtime_error("[interp] Inconsistent element type for interpolation. Can support 2D Tri or 3D Tet elements only.");
  }

  int gnNo = tMsh.gnNo;
  int eNoN = tMsh.eNoN;
  int gnEl = tMsh.gnEl;
  int nNo = 0;
  #ifdef debug_interp
  dmsg << "tMsh.gnNo: " << tMsh.gnNo;
  dmsg << "tMsh.eNoN: " << tMsh.eNoN;
  dmsg << "tMsh.gnEl: " << tMsh.gnEl;
  #endif

  // Distribute the new mesh nodes among all the processors
  // Need to transfer mesh displacement
  //
  Vector<int> gN;
  Array<double> Dg(nsd,tnNo);
  int i = nsd+1;

  for (int j = 0; j < nsd; j++) { 
    for (int k = 0; k < tnNo; k++) { 
      Dg(j,k) = rmsh.D0(j+nsd+1,k);
    }
  }
  //Dg(:,:) = rmsh.D0(i+1:i+nsd,:)

  distrn(com_mod, cm_mod, iM, tMsh, Dg, nNo, gN);
  //CALL DISTRN(iM, tMsh, Dg, nNo, gN)

  // Distribute elements of the new mesh to all processors
  int nEl = 0;
  Vector<int> gE;
  distre(com_mod, cm_mod, tMsh, nEl, gE);
  //CALL DISTRE(tMsh, nEl, gE)

  // Setup data structures for octree search
  // Get adjacent cells for source (old) mesh
  //
  Array<int> srcAdjEl;
  get_adj_esrc(com_mod, msh[iM], srcAdjEl);
  //CALL GETADJESRC(msh(iM), srcAdjEl)
  int maxKNE = srcAdjEl.nrows();
  #ifdef debug_interp
  dmsg << "maxKNE: " << maxKNE;
  #endif

  // Get adjacent nodes for each node on the new mesh
  //
  // tgtAdjNd stores node IDs.
  //
  #ifdef debug_interp
  dmsg << "get_adj_ntgt " << " ...";
  #endif
  Array<int> tgtAdjNd;
  get_adj_ntgt(com_mod, tMsh, nNo, nEl, gN, gE, tgtAdjNd);
  //CALL GETADJNTGT(tMsh, nNo, nEl, gN, gE, tgtAdjNd)

  int maxKNN = tgtAdjNd.nrows();
  //DEALLOCATE(gE)
  #ifdef debug_interp
  dmsg << "maxKNN: " << maxKNN;
  dmsg << "nNo: " << nNo;
  #endif

  Vector<double> Xp(nsd+1), Nsf(eNoN); 
  Array<double> gNsf(eNoN,nNo); 
  Vector<int> tagNd(gnNo), rootEl(nNo), masEList(msh[iM].nEl);
  //ALLOCATE(Xp(nsd+1), Nsf(eNoN), gNsf(eNoN,nNo), tagNd(gnNo), gE(nNo), rootEl(nNo), chckNp(nNo), masEList(msh(iM).nEl))
  rootEl = -1;

  std::vector<bool> chckNp(nNo);
  std::fill(chckNp.begin(), chckNp.end(), false);

  for (int e = 0; e < msh[iM].nEl; e++) {
    masEList(e) = e;
  }

  // Determine boundary nodes on the new mesh, where interpolation is
  // not needed, or boundary search is performed
  //
  Vector<int> tmpL(gnNo);
  //tmpL = -1;

  for (int e = 0; e < tMsh.fa[0].nEl; e++) {
    for (int a = 0; a < tMsh.fa[0].eNoN; a++) {
      int Ac = tMsh.fa[0].IEN(a,e);
      if (tmpL(Ac) == 0) {
        tmpL(Ac) = 1;
      }
    }
  }

  // srfNds is a really a bool array.
  //
  Vector<int>srfNds(nNo);

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    if (tmpL(Ac) > 0) {
      srfNds(a) = 1;
    }
  }

  //DEALLOCATE(tmpL)

  // tagNd stores procesors IDs ?
  int bTag = 2*cm.np();

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);

    if (srfNds(a) > 0) {
      chckNp[a] = true;
      tagNd(Ac) = bTag;
      gE(a) = -1;

      for (int i = 0; i < eNoN; i++) {
        gNsf(i,a) = 0.0;
      }
    }
  }
  
  // Find a starting element through brute search and the starting node
  // is used to initialize the queue (local numbering)
  //
  #ifdef debug_interp
  dmsg << "Find a starting element ... " << "";
  #endif
  utils::queueType rootNdQ;

  for (int a = 0; a < nNo; a++) {
    #ifdef debug_interp_1
    dmsg << "---------- a " << std::to_string(a+1) + " ----------";
    #endif
    if (chckNp[a]) {
      continue; 
    }
    int Ac = gN(a);
    Xp = 1.0;
    int Ec = -1;
    for (int i = 0; i < nsd; i++) {
      Xp(i) = tMsh.x(i,Ac);
    }

    find_n(com_mod, Xp, iM, Dg, masEList, Ec, Nsf, 0);
    #ifdef debug_interp_1
    dmsg << "Ec: " << Ec+1; 
    #endif
    //CALL FINDN(Xp, iM, Dg, masEList, Ec, Nsf)
    chckNp[a] = true;

    // Once an element in the source mesh is found, set 'this' element as
    // 'root' element for all the neighbor nodes of the target mesh. Then
    // add these neighbor nodes to the queue to form the 'front'
    //
    if (Ec > -1) {
      gE(a) = Ec;
      tagNd(Ac) = cm.tF(cm_mod);
      #ifdef debug_interp_s
      if (Ac == 8864-1) {
        dmsg << "#### 1 Ac: " << Ac+1; 
      }
      #endif

      for (int i = 0; i < eNoN; i++) {
        gNsf(i,a) = Nsf(i);
      }
      //gNsf(:,a) = Nsf(:)

      for (int nn = 0; nn < maxKNN; nn++) {
        int b = tgtAdjNd(nn,a);

        if (b > -1) {
          #ifdef debug_interp_1
          dmsg << "nn: " << std::to_string(nn+1) +  "  b: " + std::to_string(b+1); 
          #endif
          utils::enqueue(rootNdQ, b);
          rootEl(b) = Ec;
        }
      }

      if (rootNdQ.n > 1) {
        break;
      }
    }
  }

  //exit(0);

  // Node-Cell search begins here. Uses fastest grid-to-grid algorithm
  // based on advancing front in the nearest vicinity
  //
  //dmsg << " " << "";
  //dmsg << "Node-Cell search begins ... " << "";
  Vector<int> eList(maxKNE); 
  int probe;
  int pcount = 1;
 
  while (dequeue(rootNdQ, probe)) {
    #ifdef debug_interp_1
    dmsg << "---------- probe " << std::to_string(probe+1) + " ----------";
    dmsg << "pcount: " << pcount;
    #endif
    pcount += 1;

     if (std::all_of(chckNp.begin(), chckNp.end(), [](bool v) { return v; })) {
       //dmsg << "break, chckNp all true " << "";
       break;
     }
    //if (ALL(chckNp(:))) EXIT

    if (chckNp[probe]) {
      //dmsg << "continue" << "";
      continue;
    }

    int Ac = gN(probe);
    Xp = 1.0;
    for (int i = 0; i < nsd; i++) {
      Xp(i) = tMsh.x(i,Ac);
    }
    eList = -1;
    bool flag = false;
    #ifdef debug_interp_1
    dmsg << "Ac: " << Ac+1;
    dmsg << "rootEl(probe): " << rootEl(probe)+1;
    #endif

    for (int e = 0; e < maxKNE; e++) {
      int Ec = srcAdjEl(e,rootEl(probe));
      if (Ec > -1) {
        eList(e) = Ec;
        #ifdef debug_interp_1
        dmsg << "srcAdjEl:  e Ec: " << std::to_string(e+1) + " " + std::to_string(Ec+1);
        #endif
      }
    }

    for (int itry = 1; itry <= 2; itry++) {
      #ifdef debug_interp_1
      dmsg << ">>>>> try: " << itry; 
      dmsg << "Xp: " << Xp; 
      #endif
      int Ec;
      find_n(com_mod, Xp, iM, Dg, eList, Ec, Nsf, pcount);
      #ifdef debug_interp_1
      dmsg << "Ec: " << Ec+1; 
      #endif

      #ifdef debug_interp_s
      if (pcount == 19195+1) {
        dmsg << "eList: " << eList; 
        exit(0);
      }
      #endif

      // If an element is found, the neighbors of the node are initialized
      // with this element as root element and are added to the queue
      // (front propagation)

      if (Ec > -1) {
        gE(probe) = Ec;
        tagNd(Ac) = cm.tF(cm_mod);
        #ifdef debug_interp_s
        if (Ac == 8864-1) {
          dmsg << "#### 2 Ac: " << Ac+1; 
          dmsg << "     Ec: " << Ec+1; 
          dmsg << "     probe: " << probe+1; 
          dmsg << "     rootEl(probe): " << rootEl(probe)+1;
          dmsg << "     pcount: " << pcount; 
          dmsg << "     try: " << itry; 
          dmsg << "     elist: " << eList; 
          dmsg << "     Xp: " << Xp; 
        }
        #endif

        for (int i = 0; i < Nsf.size(); i++) {
          gNsf(i,probe) = Nsf(i);
        }

        for (int nn = 0; nn < maxKNN; nn++) {
          int a = tgtAdjNd(nn,probe);

          if (a > -1) {
            #ifdef debug_interp_1
            dmsg << "tgtAdjNd:   nn a Ec: " << std::to_string(nn+1) + " " + std::to_string(a+1) + " " + std::to_string(Ec+1); 
            #endif
            utils::enqueue(rootNdQ, a);
            rootEl(a) = Ec;
            #ifdef debug_interp_s
            if (a == 8864-1) {
              dmsg << "set rootEl(a): " << rootEl(a)+1; 
            }
            #endif
          }
        }

        eList.resize(maxKNE); 
        eList = -1;
        break;
      }

      // If failed in first attempt, try again with increased search area
      // but add neighbors uniquely without repetition
      //
      if (itry == 1) {
        tmpL.resize(nEl);
        tmpL = 0;
        int ne = 0;
        flag = false; 

        for (int e = 0; e < maxKNE; e++) {
          int Ec = eList(e);

          if (Ec > -1) {
            for (int a = 0; a < maxKNE; a++) {
              int Bc = srcAdjEl(a,Ec);

              if (Bc != -1 && Bc != Ec) {
                flag = true; 

                for (int i = 0; i < ne; i++) {
                  if (Bc == tmpL(i)) {
                    flag = false; 
                    break;
                  }
                }

                if (flag) {
                  tmpL(ne) = Bc;
                  ne = ne + 1;
                }
              }
            }
          }
        }

        eList.resize(ne);
        eList = -1;
        for (int i = 0; i < ne; i++) {
          eList(i) = tmpL(i);
        }

      } else if (itry == 2) {
        eList.resize(maxKNE); 
        eList = -1;
      }
    } // itry

    chckNp[probe] = true;
  }

  /*
  dmsg << "tagNd: " << "";
  for (int a = 0; a < gnNo; a++) {
    dmsg << "a tagNd: " << std::to_string(a+1) + " " + std::to_string(tagNd(a));
  }
  */

  /*
  Vector<int>::write_disabled = false;
  tagNd.write("tagNd");

  Array<int>::write_disabled = false;
  tgtAdjNd.write("tgtAdjNd");
  Array<int>::write_disabled = false;
  srcAdjEl.write("srcAdjEl");
  */

  tmpL.resize(gnNo);
  tmpL = 0;

  #ifdef debug_interp
  dmsg << "MPI_Allreduce ... " << "";
  #endif

  MPI_Allreduce(tagNd.data(), tmpL.data(), gnNo, cm_mod::mpint, MPI_MAX, cm.com());
  //CALL MPI_ALLREDUCE(tagNd, tmpL, gnNo, mpint, MPI_MAX, cm.com(), ierr)

  // tmpL was previously a list of element IDs but it now holds
  // processor IDs.
  //
  int nn = 0;
  int nbnd = 0;

  for (int a = 0; a < gnNo; a++) {
    if (tmpL(a) > 0) {
      nn = nn + 1;
    }
    if (tmpL(a) == 2*cm.np()) {
      nbnd = nbnd + 1;
    }
  }
  #ifdef debug_interp
  dmsg << "nn: " << nn;
  //exit(0);
  #endif

  // Assign tag for nodes that were interpolated in other processors
  //
  #ifdef debug_interp
  dmsg << "Assign tag for nodes ... " << "";
  #endif

  tagNd = 0;

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    tagNd(Ac) = tmpL(Ac);
  }

  // Now use brute force to find left over non-interpolated nodes
  //
  #ifdef debug_interp
  dmsg << "Find left over non-interpolated nodes ... " << "";
  #endif
  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);

    if (tagNd(Ac) == 0) {
      Xp = 1.0;
      for (int i = 0; i < nsd; i++) {
        Xp(i) = tMsh.x(i,Ac);
      }

      int Ec = -1;
      find_n(com_mod, Xp, iM, Dg, masEList, Ec, Nsf, 0);

      if (Ec > -1) {
        gE(a) = Ec;
        tagNd(Ac) = cm.tF(cm_mod);
        for (int i = 0; i < Nsf.size(); i++) {
          gNsf(i,a) = Nsf(i);
        }
      }
    }
  }

  // Nodes belonging to other procs are reassigned 0
  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    if (tagNd(Ac) != cm.tF(cm_mod) && tagNd(Ac) != bTag) {
      gE(a) = 0;
      tagNd(Ac) = 0;
      for (int i = 0; i < gNsf.nrows(); i++) {
        gNsf(i,a) = 0.0;
      }
    }
  }
  //DEALLOCATE(tmpL)

  // Now that all the elements have been found, data is interpolated
  // from the source to the target mesh
  //
  #ifdef debug_interp
  dmsg << "Interpolate from the source to the target mesh ... " << "";
  #endif
  Array<double> tmpX(lDof,nNo);

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);

    if (tagNd(Ac) == cm.tF(cm_mod)) {
      int Ec = gE(a);
      for (int i = 0; i < gNsf.nrows(); i++) {
        Nsf(i) = gNsf(i,a);
      }

      for (int i = 0; i < eNoN; i++) {
        int Bc = msh[iM].IEN(i,Ec);
        Bc = msh[iM].lN(Bc);
        for (int j = 0; j < tmpX.nrows(); j++) {
          tmpX(j,a) = tmpX(j,a) + Nsf(i)*sD(j,Bc);
        }
      }
    }
  }

  // Since there is no direct mapping for face data, we use L2 norm
  // to find the nearest face node and copy its solution. This requires
  // face node/IEN structure to NOT be changed during remeshing.
  //
  tmpL.resize(nNo);
  bool flag = false; 

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);

    if (srfNds(a) > 0) {
      flag = false; 

      for (int iFa = 0; iFa < msh[iM].nFa; iFa++) {
        auto& fa = msh[iM].fa[iFa];

        for (int b = 0; b <fa.nNo; b++) {
          int Bc = fa.gN(b);
          double dS = 0.0; 
          for (int i = 0; i < nsd; i++) {
            double sum = com_mod.x(i,Bc) + Dg(i,Bc) - tMsh.x(i,Ac);
            dS += dS * dS;
          }

          dS = sqrt(dS);

          if (dS < 1.E-12) {
            tmpL(a) = Bc;
            flag = true;
            break;
          }
        }

        if (flag) break;
      }

      if (flag) {
        int Bc = msh[iM].lN(tmpL(a));
        for (int i = 0; i < tmpX.nrows(); i++) { 
          tmpX(i,a) = sD(i,Bc);
        }
      } else { 
        tagNd(Ac) = 0;
      }
    }
  }
  //DEALLOCATE(tmpL)

  // Map the tagged nodes and solution to local vector within a proc,
  // including boundary nodes. Since the boundary nodes can be overlapping
  // across different procs, these are repeated. But this will not cause
  // problem as the solution is simply overwritten depending on the face pointer.
  //
  #ifdef debug_interp
  dmsg << "Map the tagged nodes and solutio ... " << "";
  #endif
  nn = 0;

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    if (tagNd(Ac) > 0) nn = nn + 1;
  }
  tmpL.resize(nn);
  nn = 0;

  for (int a = 0; a < nNo; a++) {
    int Ac = gN(a);
    if (tagNd(Ac) > 0) {
      tmpL(nn) = a;
      nn = nn + 1;
    }
  }

  gNsf.resize(lDof,nn);

  for (int i = 0; i < nn; i++) {
    int a = tmpL(i);
    int Ac = gN(a);
    for (int j = 0; j < tmpX.nrows(); j++) { 
      gNsf(j,i) = tmpX(j,a);
    }
    tmpL(i) = Ac;
  }

  Vector<int> disp;
  if (cm.mas(cm_mod)) {
    disp.resize(cm.np());
  } else { 
    //ALLOCATE(disp(0))
  }

  #ifdef debug_interp
  dmsg << "MPI_Gather ... " << "";
  #endif
  MPI_Gather(&nn, 1, cm_mod::mpint, disp.data(), 1, cm_mod::mpint, cm_mod.master, cm.com());
  //CALL MPI_GATHER(nn, 1, mpint, disp, 1, mpint, master, cm.com(), ierr)

  Vector<int> sCount, gvec;

  if (cm.mas(cm_mod)) {
    int i = disp.sum();
    i = i * (1 + lDof);
    gvec.resize(i);
    sCount.resize(cm.np());
    for (int i = 0; i < cm.np(); i++) {
      sCount(i) = disp(i)*(1+lDof);
    }
    disp(0) = 0;
    for (int i = 1; i < cm.np(); i++) {
      disp(i) = disp(i-1) + sCount(i-1);
    }
  } else { 
    //ALLOCATE(gvec(0), sCount(0))
  }

  Vector<int> vec((lDof+1)*nn);
  int e = 0;

  for (int a = 0; a < nn; a++) {
    vec(e) = static_cast<double>(tmpL(a));
    e = e + 1;

    for (int b = 0; b < lDof; b++) {
      vec(e) = gNsf(b,a);
      e = e + 1;
    }
  }

  #ifdef debug_interp
  dmsg << "" << "";
  dmsg << "MPI_Gatherv ... " << "";
  dmsg << "nn: " << nn;
  dmsg << "(1+lDof)*nn: " << (1+lDof)*nn;
  dmsg << "vec.size(): " << vec.size();
  dmsg << "gvec.size(): " << gvec.size();
  dmsg << "sCount.size(): " << sCount.size();
  dmsg << "sCount(0): " << sCount(0);
  dmsg << "disp.size(): " << disp.size();
  dmsg << "disp(0): " << disp(0);
  int csize;
  MPI_Comm_size( cm.com(), &csize);
  dmsg << "csize: " << csize;
  #endif

  if (disp.size() > 1) {
  MPI_Gatherv(vec.data(), (1+lDof)*nn, cm_mod::mpreal, gvec.data(), sCount.data(), disp.data(), 
      cm_mod::mpreal, cm_mod.master, cm.com());
  //CALL MPI_GATHERV(vec, (1+lDof)*nn, mpreal, gvec, sCount, disp, mpreal, master, cm.com(), ierr)
  } else {
    gvec = vec;
  }

  #ifdef debug_interp
  dmsg << "done MPI_Gatherv " << "";
  //exit(0);
  #endif

  if (cm.mas(cm_mod)) {
    #ifdef debug_interp
    dmsg << "" << "";
    dmsg << "Set tgD ... " << "";
    #endif
    tgD = 0.0;
    nn = sCount.sum();
    i = 0;
    dmsg << "nn: " << nn;

    while (true) { 
      int Ac = round(gvec(i));
      //dmsg << "  Ac: " << std::to_string(i+1) + " " + std::to_string(Ac);
      i = i + 1;

      for (int b = 0; b < lDof; b++) {
        tgD(b,Ac) = gvec(i);
        i = i + 1;
      }
      if (i == nn) break;
    }
  }

  #ifdef debug_interp
  dmsg << "i: " << i;
  //exit(0);
  #endif
}

//-------------
// int_msh_srf
//-------------
// Reproduces 'SUBROUTINE INTMSHSRF(lM, lFa)'
//
void int_msh_srf(ComMod& com_mod, CmMod& cm_mod, mshType& lM, faceType& lFa)
{
  auto& cm = com_mod.cm;
  #define n_debug_int_msh_srf
  #ifdef debug_int_msh_srf 
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  auto& rmsh = com_mod.rmsh;

  int eNoNb = lM.fa[0].eNoN;
  lFa.eNoN = eNoNb;

  lFa.gnEl = 0;
  for (auto& face : lM.fa) {
    lFa.gnEl += face.gnEl;
  }
  //lFa.gnEl = SUM(lM.fa(:).gnEl)

  lFa.nNo  = 0;
  lFa.name = "old_mesh_surface";
  lFa.nEl = lFa.gnEl;

  lFa.IEN.resize(lFa.eNoN,lFa.gnEl); 
  lFa.gE.resize(lFa.gnEl);
  //ALLOCATE(lFa.IEN(lFa.eNoN,lFa.gnEl), lFa.gE(lFa.gnEl))
  //dmsg << "lFa.gnEle: " << lFa.gnEl;

  if (eNoNb != lM.fa[lM.nFa-1].eNoN) {
    //wrn = "    Remeshing not formulated for different face types"
    return;
  }

  if (cm.mas(cm_mod)) {
    for (int iFa = 0; iFa < lM.nFa; iFa++) {
      //dmsg << "---------- iFa " << std::to_string(iFa+1) + " ----------";
      int eoff = 0;
      if (iFa > 0) {
        for (int i = 0; i < iFa; i++) {
          eoff += lM.fa[i].gnEl;
        }
        //eoff = SUM(lM.fa(1:iFa-1).gnEl)
      }
      //dmsg << "eoff: " << eoff;

      for (int e = eoff; e < eoff+lM.fa[iFa].gnEl; e++) {
        lFa.gE(e) = lM.fa[iFa].gebc(0,e-eoff);
        //lFa.gE(e) = lM.fa(iFa).gebc(1,e-eoff)
        //dmsg << "e: " << e+1;
        //dmsg << "lFa.gE(e): " << lFa.gE(e);

        for (int i = 0; i < lFa.eNoN; i++) {
          lFa.IEN(i,e) = lM.fa[iFa].gebc(i+1,e-eoff);
        }
        //dmsg << "lFa.IEN(:,e):  " << lFa.IEN.col(e);
        //lFa.IEN(:,e) = lM.fa(iFa).gebc(2:1+eNoNb,e-eoff)
      }
    }

    read_msh_ns::calc_nbc(lM, lFa);
    //CALL CALCNBC(lM, lFa)

    Vector<int> incNd(lM.gnNo);

    for (int a = 0; a < lFa.nNo; a++) {
      int Ac = lFa.gN(a);
      incNd(Ac) = a;
    }

    for (int e = 0; e < lFa.gnEl; e++) {
      for (int a = 0; a < lFa.eNoN; a++) {
        int Ac = lFa.IEN(a,e);
        lFa.IEN(a,e) = incNd(Ac);
      }
    }
  }

  cm.bcast(cm_mod, &lFa.nNo);

  if (cm.slv(cm_mod)) {
    lFa.gN.resize(lFa.nNo);
    //ALLOCATE(lFa.gN(lFa.nNo))
   }
}

//-------------
// remesher_3d
//-------------
//
void remesher_3d(ComMod& com_mod, CmMod& cm_mod, int iM, faceType& lFa, mshType& lM)
{
  using namespace consts;

  #define debug_remesher_3d 
  #ifdef debug_remesher_3d
  auto& cm = com_mod.cm;
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  dmsg << "iM: " << iM;
  #endif

  auto& rmsh = com_mod.rmsh;
  dmsg << "rmsh.maxEdgeSize(iM): " << rmsh.maxEdgeSize(iM);

  std::array<double,3> rparams = {
    rmsh.maxRadRatio,
    rmsh.minDihedAng,
    rmsh.maxEdgeSize(iM)
  };

  int iOK = 0;

  if (rmsh.method == MeshGeneratorType::RMSH_TETGEN) {
     remesh3d_tetgen(lFa.nNo, lFa.nEl, lFa.x.data(), lFa.IEN.data(), rparams, &iOK);
     //CALL remesh3d_tetgen(lFa.nNo, lFa.nEl, lFa.x, lFa.IEN, rparams, iOK)
     //if (iOK .LT. 0)
     //2  err = "Fatal! TetGen returned with error. Check log"
  } else { 
     //err = "Unknown remesher choice."
  }

  std::string elem_file_name  = "new-vol-mesh-cpp.ele";
  std::ifstream new_elem_mesh;
  new_elem_mesh.open(elem_file_name);

  // Get the number of remeshed elements.
  //
  lM.gnEl = 0;
  int value;
  std::string line;

  while (std::getline(new_elem_mesh, line)) {
    lM.gnEl += 1;
  }
/*
  fileType fTmp;
  fTmp.fname  = "new-vol-mesh.ele"
  fid = fTmp.open()
  lM.gnEl = 0
  for (int 
     READ (fid,*,END=111)
     lM.gnEl = lM.gnEl + 1
  }
 111  REWIND(fid)
*/

  new_elem_mesh.clear();
  new_elem_mesh.seekg(0);

  if (rmsh.method == MeshGeneratorType::RMSH_TETGEN) {
    lM.gnEl = lM.gnEl - 1;
  }
  std::cout << "[remesher_3d] Number of elements after remesh " << lM.gnEl << std::endl;
  lM.gIEN.resize(lM.eNoN,lM.gnEl);

  // Read element connectivity.
  //
  bool first_list = true;
  int e = 0;

  while (std::getline(new_elem_mesh, line)) {
    if (line == "" || first_list) {
      first_list = false;
      continue;
    }
    std::istringstream line_input(line);
    std::vector<int> values;
    //std::cout << "[remesher_3d] conn ... "  << std::endl;
    //std::cout << "[remesher_3d] " << n << ": ";
    int i = 0;
    line_input >> value;
    while (line_input >> value) {
      //std::cout << value << "  ";
      //values.push_back(value);
      lM.gIEN(i,e) = value;
      i += 1;
    }

    e += 1;
    //std::cout << std::endl;
  }

  /*
  if (ALLOCATED(lM.gIEN)) DEALLOCATE(lM.gIEN)
  ALLOCATE(lM.gIEN(lM.eNoN,lM.gnEl))
  if (rmsh.method .EQ. RMSH_TETGEN) READ(fid,*)
  for (int  e=1, lM.gnEl
     READ(fid,*) i, lM.gIEN(:,e)
  }
  CLOSE(fid, STATUS='DELETE')
  */

  // Get the number of remeshed nodes.
  //
  std::string node_file_name  = "new-vol-mesh-cpp.node";
  std::ifstream new_node_mesh;
  new_node_mesh.open(node_file_name);

  lM.gnNo = 0;
  while (std::getline(new_node_mesh, line)) {
    //std::cout << line << std::endl;
    lM.gnNo += 1;
  }

  if (rmsh.method == MeshGeneratorType::RMSH_TETGEN) {
    lM.gnNo = lM.gnNo - 1;
  }
 
  std::cout << "[remesher_3d] Number of nodes after remesh " << lM.gnNo << std::endl;
  new_node_mesh.clear();
  new_node_mesh.seekg(0);
  lM.x.resize(com_mod.nsd,lM.gnNo);

  // Read node coordinates.
  //
  first_list = true;
  int n = 0;

  while (std::getline(new_node_mesh, line)) {
    if (line == "" || first_list) {
      first_list = false;
      continue;
    }
    std::istringstream line_input(line);
    int i = 0;
    double value;
    line_input >> value;
    while (line_input >> value) {
      lM.x(i,n) = value;
      i += 1;
    }

    n += 1;
  }

  /*
  fTmp.fname  = "new-vol-mesh.node"
  fid = fTmp.open()
  lM.gnNo = 0
  for (int 
     READ (fid,*,END=112)
     lM.gnNo = lM.gnNo + 1
  }
 112  REWIND(fid)
  if (rmsh.method .EQ. RMSH_TETGEN) lM.gnNo = lM.gnNo - 1
  std = "    Number of vertices after remesh "//STR(lM.gnNo)

  if (ALLOCATED(lM.x)) DEALLOCATE(lM.x)
  ALLOCATE(lM.x(nsd,lM.gnNo))
  if (rmsh.method .EQ. RMSH_TETGEN) READ(fid,*)
  for (int  Ac=1, lM.gnNo
     READ(fid,*) i, lM.x(:,Ac)
  }
  CLOSE(fid, STATUS='DELETE')
  */

  nn::select_ele(com_mod, lM);
}

//----------------
// remesh_restart
//----------------
// Reproduces Fortran 'SUBROUTINE REMESHRESTART(timeP)'
//
void remesh_restart(Simulation* simulation)
{
  using namespace consts;

  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;
  auto& cm = com_mod.cm;
  auto& chnl_mod = simulation->chnl_mod;

  #define debug_remesh_restart
  #ifdef debug_remesh_restart 
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  auto& stFileName = com_mod.stFileName;
  auto& rmsh = com_mod.rmsh;

  auto sTmp = stFileName + "_last_cpp.bin";
  #ifdef debug_remesh_restart 
  dmsg << "sTmp: " << sTmp;
  #endif
  std::filesystem::remove(sTmp);

  // Create the file.
  //
  if (cm.mas(cm_mod)) {
    std::ofstream restart_file(sTmp, std::ios::out | std::ios::binary);
    restart_file.close();
  }

  // This call is to block all processors
  cm.bcast(cm_mod, &rmsh.rTS);

  // Write something.
  //
  auto fTmp = stFileName + "_" + std::to_string(rmsh.rTS) + "_cpp.bin";
  auto const recLn = com_mod.recLn;
  const bool dFlag = com_mod.dFlag;
  auto& timeP = com_mod.timeP;
 
  #ifdef debug_remesh_restart 
  dmsg << "fTmp: " << fTmp;
  dmsg << "cm.tF: " << cm.tF(cm_mod);
  dmsg << "dFlag: " << dFlag;
  #endif

  std::ofstream restart_file(fTmp, std::ios::out | std::ios::binary | std::ios::app);
  std::streampos write_pos = (cm.tF(cm_mod) - 1) * recLn;
  restart_file.seekp(write_pos);
  output::write_restart_header(com_mod, timeP, restart_file);

  auto& cplBC = com_mod.cplBC;
  restart_file.write((char*)cplBC.xn.data(), cplBC.xn.msize());
  restart_file.write((char*)rmsh.Y0.data(), rmsh.Y0.msize());
  restart_file.write((char*)rmsh.A0.data(), rmsh.A0.msize());

  if (dFlag) {
    restart_file.write((char*)rmsh.D0.data(), rmsh.D0.msize());
  }

  restart_file.close();

  if (cm.mas(cm_mod)) {
    std::string link_cmd = "ln -f " + fTmp + " " + sTmp;
    std::system(link_cmd.c_str());
  }

  auto& x = com_mod.x;
  auto& gtnNo = com_mod.gtnNo;
  const int nsd = com_mod.nsd;
  const int nMsh = com_mod.nMsh;
  const int tDof= com_mod.tDof;
  gtnNo = 0;
  int lDof = 3 * com_mod.tDof;
  mshType tMsh;

  Array<double> gtX, gtD, gX, gD;
  Array<double> tempD, tempX, gnD;

  for (int iM = 0; iM < nMsh; iM++) {
    auto& msh = com_mod.msh[iM];

    if (rmsh.flag[iM]) {
      if (rmsh.method == MeshGeneratorType::RMSH_TETGEN) {
         //std = " Remeshing <"//CLR(TRIM(msh(iM).name))// "> using <"//CLR("Tetgen")//"> library at time "// STR(rmsh.rTS)
      } else { 
         //err = "Unexpected behavior in Remesher"
      }

      tempX.resize(nsd,msh.nNo);
      gD.resize(lDof,msh.nNo);
      //ALLOCATE(tempX(nsd,msh(iM).nNo))
      //ALLOCATE(gD(lDof,msh(iM).nNo))

      for (int a = 0; a < msh.nNo; a++) {
        int Ac = msh.gN(a);
        for (int i = 0; i < nsd; i++) {
          tempX(i,a) = x(i,Ac) + rmsh.D0(i+nsd+1);
        }
        //tempX(:,a) = x(:,Ac) + rmsh.D0(nsd+2:2*nsd+1,Ac)
   
        for (int i = 0; i < tDof; i++) {
          gD(i,a) = rmsh.A0(i,Ac);
          gD(i+tDof-1,a) = rmsh.Y0(i,Ac);
          gD(2*tDof+i,a) = rmsh.D0(i,Ac);
        }
        //gD(1:tDof,a) = rmsh.A0(:,Ac)
        //gD(tDof+1:2*tDof,a) = rmsh.Y0(:,Ac)
        //gD(2*tDof+1:3*tDof,a) = rmsh.D0(:,Ac)
      }

      tMsh.nFa = 1;
      tMsh.fa.resize(tMsh.nFa);
      //ALLOCATE(tMsh.fa(tMsh.nFa))

      if (cm.mas(cm_mod)) {
        tMsh.gnNo = msh.gnNo;
        tMsh.gnEl = msh.gnEl;
        tMsh.eNoN = msh.eNoN;
      } else {
        tMsh.gnNo = 0;
        tMsh.gnEl = 0;
        tMsh.eNoN = 0;
      }

      gX.resize(nsd,tMsh.gnNo);
      //ALLOCATE(gX(nsd,tMsh.gnNo));
      gX = all_fun::global(com_mod, cm_mod, msh, tempX);
      //gX = GLOBAL(msh(iM), tempX);
      //DEALLOCATE(tempX);

      if (cm.mas(cm_mod)) {
        tMsh.gIEN.resize(tMsh.eNoN,tMsh.gnEl);
        //ALLOCATE(tMsh.gIEN(tMsh.eNoN,tMsh.gnEl))
        tMsh.gIEN = msh.gIEN;

        for (int e = 0; e < tMsh.gnEl; e++) {
          int Ec = msh.otnIEN(e);
          for (int a = 0; a < tMsh.eNoN; a++) {
            msh.gIEN(a,e) = tMsh.gIEN(a,Ec);
          }
        }
      }

      int_msh_srf(com_mod, cm_mod, msh, tMsh.fa[0]);
      //CALL INTMSHSRF(msh(iM), tMsh.fa(1))

      if (cm.mas(cm_mod)) {
        tMsh.fa[0].x.resize(nsd,tMsh.fa[0].nNo);
        //ALLOCATE(tMsh.fa(1).x(nsd,tMsh.fa(1).nNo))
        for (int a = 0; a < tMsh.fa[0].nNo; a++) {
          int Ac = tMsh.fa[0].gN(a);
          for (int i = 0; i < nsd; i++) {
            tMsh.fa[0].x(i,a) = gX(i,Ac);
          }
          //tMsh.fa(1).x(:,a) = gX(:,Ac)
        }

        if (nsd == 2) {
          //err = "Remesher not yet developed for 2D objects"
        } else {
          remesher_3d(com_mod, cm_mod, iM, tMsh.fa[0], tMsh);
          //CALL REMESHER_3D(iM, tMsh.fa(1), tMsh)
        }

        gnD.resize(lDof,tMsh.gnNo);
        //ALLOCATE(gnD(lDof,tMsh.gnNo))
      } else {
        //ALLOCATE(tMsh.fa(1).x(0,0), gnD(0,0))
      }

      cm.bcast(cm_mod, &tMsh.gnNo);
      cm.bcast(cm_mod, &tMsh.eNoN);
      cm.bcast(cm_mod, &tMsh.gnEl);
      //CALL cm.bcast(tMsh.gnNo)
      //CALL cm.bcast(tMsh.eNoN)
      //CALL cm.bcast(tMsh.gnEl)
      dmsg << "tMsh.gnNo: " << tMsh.gnNo;
      dmsg << "tMsh.eNoN: " << tMsh.eNoN;

      if (cm.slv(cm_mod)) {
        tMsh.x.resize(nsd,tMsh.gnNo);
        tMsh.gIEN.resize(tMsh.eNoN,tMsh.gnEl);
        //ALLOCATE(tMsh.x(nsd,tMsh.gnNo))
        //ALLOCATE(tMsh.gIEN(tMsh.eNoN,tMsh.gnEl))
      }
      dmsg << "1 size tMsh.x: " << tMsh.x.size();

      MPI_Bcast(tMsh.x.data(), nsd*tMsh.gnNo, cm_mod::mpreal, cm_mod.master, cm.com());
      //CALL MPI_BCAST(tMsh.x, nsd*tMsh.gnNo, mpreal, master, cm.com(), ierr)

      MPI_Bcast(tMsh.gIEN.data(), tMsh.eNoN*tMsh.gnEl, cm_mod::mpint, cm_mod.master, cm.com());
      //CALL MPI_BCAST(tMsh.gIEN, tMsh.eNoN*tMsh.gnEl, mpint, master, cm.com(), ierr)

      MPI_Bcast(tMsh.fa[0].IEN.data(), tMsh.fa[0].eNoN * tMsh.fa[0].nEl, cm_mod::mpint, cm_mod.master, cm.com());
      //CALL MPI_BCAST(tMsh.fa(1).IEN, tMsh.fa(1).eNoN * tMsh.fa(1).nEl, mpint, master, cm.com(), ierr)

      MPI_Bcast(tMsh.fa[0].gN.data(), tMsh.fa[0].nNo, cm_mod::mpint, cm_mod.master, cm.com());
      //CALL MPI_BCAST(tMsh.fa(1).gN, tMsh.fa(1).nNo, mpint, master, cm.com(), ierr)

      dmsg << "2 size tMsh.x: " << tMsh.x.size();
      dmsg << "fTmp: " << fTmp;

      interp(com_mod, cm_mod, lDof, iM, tMsh, gD, gnD);
      //CALL INTERP(lDof, iM, tMsh, gD, gnD)

      /*
      for (int i = 0; i < lDof; i++) {
        for (int j = 0; j < msh.nNo; j++) {
          std::cout << "i j gD: " << i+1 << " " << j+1 << " " << gD(i,j) << std::endl;
        }
      }
      */

      if (cm.mas(cm_mod)) {
        msh.gnNo = tMsh.gnNo;
        int a = gtnNo + msh.gnNo;

        if (iM > 0) {
          tempX.resize(nsd,gtnNo); 
          tempD.resize(lDof,gtnNo);
          //ALLOCATE(tempX(nsd,gtnNo), tempD(lDof,gtnNo))
          tempX = gtX;
          tempD = gtD;

          gtX.resize(nsd,a); 
          gtD.resize(lDof,a);
          //ALLOCATE(gtX(nsd,a), gtD(lDof,a))

          for (int i = 0; i < nsd; i++) {
            for (int j = 0; j < gtnNo; j++) {
              gtX(i,j) = tempX(i,j);
              gtD(i,j) = tempD(i,j);
            }
          }
          //gtX(:,1:gtnNo) = tempX(:,:)
          //gtD(:,1:gtnNo) = tempD(:,:)
          //DEALLOCATE(tempX, tempD)
        } else {
          gtX.resize(nsd,a);
          //ALLOCATE(gtX(nsd,a))

          gtD.resize(lDof,a);
          //ALLOCATE(gtD(lDof,a))
        }

        msh.x.resize(nsd, msh.gnNo);
        //ALLOCATE(msh(iM).x(nsd,msh(iM).gnNo))

        for (int i = 0; i < nsd; i++) {
          for (int j = gtnNo; j < a; j++) {
            gtX(i,j) = tMsh.x(i,j) - gnD(i + 2*tDof+nsd+1, j);
            gtD(i,j) = gnD(i,j);
            msh.x(i,j) = gtX(i,j);
          }
        }
        //gtX(:,gtnNo+1:a) = tMsh.x(:,:) - gnD(2*tDof+nsd+2:2*tDof+2*nsd+1,:)
        //gtD(:,gtnNo+1:a) = gnD(:,:) msh(iM).x(:,:) = gtX(:,gtnNo+1:a)

        gtnNo = a;

        set_face_ebc(com_mod, cm_mod, tMsh.fa[0], tMsh);
        //CALL SETFACEEBC(tMsh.fa(1), tMsh)

        msh.eNoN = tMsh.eNoN;
        msh.gnEl = tMsh.gnEl;


        msh.gIEN.resize(msh.eNoN,msh.gnEl);
        for (int i = 0; i < msh.eNoN; i++) {
          for (int j = 0; i < msh.gnEl; i++) {
            msh.gIEN(i,j) = tMsh.gIEN(i,j);
          }
        }
        //msh(iM).gIEN(:,:) = tMsh.gIEN(:,:)

        dist_msh_srf(com_mod, chnl_mod, tMsh.fa[0], msh, 1);
        //CALL DISTMSHSRF(tMsh.fa(1), msh(iM), 1)

        // HERE DaveP
        exit(0);

        sTmp = chnl_mod.appPath + ".remesh_tmp.dir";
        fTmp = sTmp + msh.name +  "_" + std::to_string(rmsh.rTS) + ".vtu";
        msh.gIEN -= 1;

        //sTmp = TRIM(appPath)//".remesh_tmp.dir"
        //fTmp = TRIM(sTmp)//"/"//TRIM(msh(iM).name)//  "_"//STR(rmsh.rTS)//".vtu"
        //msh(iM).gIEN(:,:) = msh(iM).gIEN(:,:) - 1
        //CALL WRITEVTU(msh(iM), fTmp)
        //msh(iM).gIEN(:,:) = msh(iM).gIEN(:,:) + 1

      } else {
        //if (.NOT.ALLOCATED(gtX)) }
          //ALLOCATE(gtX(0,0), gtD(0,0))
         //}
      }

      //CALL MPI_BARRIER(cm.com(), ierr)

      //CALL DESTROY(tMsh)
      //DEALLOCATE(gX, gD, gnD)

    } else {
     //ALLOCATE(tempX(nsd,msh(iM).nNo))
     //ALLOCATE(tempD(lDof,msh(iM).nNo))

     for (int a = 0; a < msh.nNo; a++) {
       int Ac = msh.gN(a);
       //tempX(1:nsd,a) = x(:,Ac)
       //tempD(1:tDof,a) = rmsh.A0(:,Ac)
       //tempD(tDof+1:2*tDof,a) = rmsh.Y0(:,Ac)
       //tempD(2*tDof+1:3*tDof,a) = rmsh.D0(:,Ac)
      }

      int a;

      if (cm.mas(cm_mod)) {
        a = msh.gnNo;
      } else {
        a = 0;
      }

      //ALLOCATE(gX(nsd,a), gD(lDof,a))
      //gX = GLOBAL(msh(iM), tempX)
      //gD = GLOBAL(msh(iM), tempD)

      if (cm.mas(cm_mod)) {
        //ALLOCATE(tMsh.gIEN(msh(iM).eNoN,msh(iM).gnEl))
        tMsh.gIEN = msh.gIEN;

        for (int e = 0; e < msh.gnEl; e++) { 
          int Ec = msh.otnIEN(e);
          //msh.gIEN(:,e) = tMsh.gIEN(:,Ec)
        }

        //DEALLOCATE(tMsh.gIEN)
      }

      //DEALLOCATE(tempX, tempD)
      //if (ALLOCATED(msh(iM).x)) DEALLOCATE(msh(iM).x)
      //ALLOCATE(msh(iM).x(nsd,msh(iM).gnNo))

      tMsh.nFa = 1;
      //ALLOCATE(tMsh.fa(tMsh.nFa))
      //CALL INTMSHSRF(msh(iM), tMsh.fa(1))

      if (cm.mas(cm_mod)) {
        //ALLOCATE(tMsh.fa(1).x(nsd,tMsh.fa(1).nNo))
        for (int i = 0; i < tMsh.fa[0].nNo; i++) {
          int Ac = tMsh.fa[0].gN(i);
          //tMsh.fa(1).x(:,i) = gX(:,Ac)
        }
        msh.x = gX;

        //CALL DISTMSHSRF(tMsh.fa(1), msh(iM), 2)

        a = gtnNo + a;

        if (iM > 0) {
          //ALLOCATE(tempX(nsd,gtnNo))
          //ALLOCATE(tempD(lDof,gtnNo))
          tempX = gtX;
          tempD = gtD;
          //if (ALLOCATED(gtX)) DEALLOCATE(gtX)
          //if (ALLOCATED(gtD)) DEALLOCATE(gtD)
          //ALLOCATE(gtX(nsd,a), gtD(lDof,a))
          //gtX(:,1:gtnNo) = tempX(:,:)
          //gtD(:,1:gtnNo) = tempD(:,:)
          //DEALLOCATE(tempX, tempD)
        } else {
          //ALLOCATE(gtX(nsd,a), gtD(lDof,a))
        }

        //gtX(:,gtnNo+1:a) = gX(:,:)
        //gtD(:,gtnNo+1:a) = gD(:,:)
        gtnNo = a;

      } else {
        //ALLOCATE(tMsh.fa(1).x(0,0))
        //if (ALLOCATED(gtX)) DEALLOCATE(gtX)
        //if (ALLOCATED(gtD)) DEALLOCATE(gtD)
        //ALLOCATE(gtX(0,0), gtD(0,0))
      }
      //DEALLOCATE(gX, gD)
      //CALL DESTROY(tMsh)
    } // reMesh flag
  }

  //CALL cm.bcast(gtnNo)
  //DEALLOCATE(x, rmsh.A0, rmsh.Y0, rmsh.D0)

  if (cm.mas(cm_mod)) {
    //ALLOCATE(x(nsd,gtnNo))
    //ALLOCATE(rmsh.A0(tDof,gtnNo))
    //ALLOCATE(rmsh.Y0(tDof,gtnNo))
    //ALLOCATE(rmsh.D0(tDof,gtnNo))
    x = gtX;

    for (int a = 0; a < gtnNo; a++) { 
      //rmsh.A0(:,a) = gtD(1:tDof,a)
      //rmsh.Y0(:,a) = gtD(tDof+1:2*tDof,a)
      //rmsh.D0(:,a) = gtD(2*tDof+1:3*tDof,a)
    }

  } else {
    //ALLOCATE(rmsh.A0(0,0),rmsh.Y0(0,0),rmsh.D0(0,0))
  }

  //DEALLOCATE(gtX, gtD)

  for (int iM = 0; iM < nMsh; iM++) {
    auto& msh = com_mod.msh[iM];

    if (cm.mas(cm_mod)) {
      //if (ALLOCATED(msh(iM).eDist))  DEALLOCATE(msh(iM).eDist)
      //if (ALLOCATED(msh(iM).eId))    DEALLOCATE(msh(iM).eId)
      //if (ALLOCATED(msh(iM).gN))     DEALLOCATE(msh(iM).gN)
      //if (ALLOCATED(msh(iM).gpN))    DEALLOCATE(msh(iM).gpN)
      //if (ALLOCATED(msh(iM).IEN))    DEALLOCATE(msh(iM).IEN)
      //if (ALLOCATED(msh(iM).otnIEN)) DEALLOCATE(msh(iM).otnIEN)
      //if (ALLOCATED(msh(iM).INN))    DEALLOCATE(msh(iM).INN)
      //if (ALLOCATED(msh(iM).lN))     DEALLOCATE(msh(iM).lN)
      //if (ALLOCATED(msh(iM).eIEN))   DEALLOCATE(msh(iM).eIEN)
      //if (ALLOCATED(msh(iM).sbc))    DEALLOCATE(msh(iM).sbc)
      //if (ALLOCATED(msh(iM).iGC))    DEALLOCATE(msh(iM).iGC)
      //if (ALLOCATED(msh(iM).nW))     DEALLOCATE(msh(iM).nW)
      //if (ALLOCATED(msh(iM).w))      DEALLOCATE(msh(iM).w)
      //if (ALLOCATED(msh(iM).xi))     DEALLOCATE(msh(iM).xi)
      //if (ALLOCATED(msh(iM).xib))    DEALLOCATE(msh(iM).xib)
      //if (ALLOCATED(msh(iM).x))      DEALLOCATE(msh(iM).x)
      //if (ALLOCATED(msh(iM).N))      DEALLOCATE(msh(iM).N)
      //if (ALLOCATED(msh(iM).Nb))     DEALLOCATE(msh(iM).Nb)
      //if (ALLOCATED(msh(iM).nV))     DEALLOCATE(msh(iM).nV)
      //if (ALLOCATED(msh(iM).fN))     DEALLOCATE(msh(iM).fN)
      //if (ALLOCATED(msh(iM).Nx))     DEALLOCATE(msh(iM).Nx)
      //if (ALLOCATED(msh(iM).Nxx))    DEALLOCATE(msh(iM).Nxx)

      //CALL DESTROY(msh(iM).nAdj)
      //CALL DESTROY(msh(iM).eAdj)

      for (int i = 0; i < msh.nFs; i++) {
        //CALL DESTROY(msh(iM).fs(i))
      }
      //if (ALLOCATED(msh(iM).fs))     DEALLOCATE(msh(iM).fs)
    } else {
      //CALL DESTROY(msh(iM))
    }
  }

  //if (cm.slv()) DEALLOCATE(msh)

/*
  if (ALLOCATED(eq)) {
    for (int  iEq=1, nEq
      CALL DESTROY(eq(iEq))
    }
    DEALLOCATE(eq)
  }
*/

  /*
  if (ALLOCATED(colPtr))   DEALLOCATE(colPtr)
  if (ALLOCATED(dmnID))    DEALLOCATE(dmnID)
  if (ALLOCATED(ltg))      DEALLOCATE(ltg)
  if (ALLOCATED(rowPtr))   DEALLOCATE(rowPtr)
  if (ALLOCATED(idMap))    DEALLOCATE(idMap)
  if (ALLOCATED(cmmBdry))  DEALLOCATE(cmmBdry)
  if (ALLOCATED(iblank))   DEALLOCATE(iblank)
  if (ALLOCATED(Ao))       DEALLOCATE(Ao)
  if (ALLOCATED(An))       DEALLOCATE(An)
  if (ALLOCATED(Do))       DEALLOCATE(Do)
  if (ALLOCATED(Dn))       DEALLOCATE(Dn)
  if (ALLOCATED(R))        DEALLOCATE(R)
  if (ALLOCATED(Val))      DEALLOCATE(Val)
  if (ALLOCATED(Yo))       DEALLOCATE(Yo)
  if (ALLOCATED(Yn))       DEALLOCATE(Yn)
  if (ALLOCATED(Bf))       DEALLOCATE(Bf)
  */
  cplBC.nFa = 0;

  // Additional physics based variables to be deallocated
/*
  if (ALLOCATED(Ad))       DEALLOCATE(Ad)
  if (ALLOCATED(Rd))       DEALLOCATE(Rd)
  if (ALLOCATED(Kd))       DEALLOCATE(Kd)
  if (ALLOCATED(pS0))      DEALLOCATE(pS0)
  if (ALLOCATED(pSn))      DEALLOCATE(pSn)
  if (ALLOCATED(pSa))      DEALLOCATE(pSa)
*/




  exit(0);

}

//--------------
// set_face_ebc
//--------------
//
void set_face_ebc(ComMod& com_mod, CmMod& cm_mod, faceType& lFa, mshType& lM)
{
  const int nsd = com_mod.nsd;

  Vector<int> nAssocEl(lM.gnNo);

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      nAssocEl(Ac) = nAssocEl(Ac) + 1;
    }
  }

  int maxAssocEl = nAssocEl.max();
  Array<int> assocEl(maxAssocEl, lM.gnNo);
  nAssocEl = 0;

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      assocEl(nAssocEl(Ac),Ac) = e;
      nAssocEl(Ac) = nAssocEl(Ac) + 1;
    }
  }

  lFa.gE.resize(lFa.nEl);
  Array<int> bin(maxAssocEl,2);

  for (int e = 0; e < lFa.gnEl; e++) {
    bin = 0;

    for (int a = 0; a < lFa.eNoN; a++) {
      int Ac = lFa.IEN(a,e);

      if (a == 0) {
        for (int i = 0; i < nAssocEl(Ac); i++) {
          bin(i,0) = assocEl(i,Ac);
          bin(i,1) = 1;
        }
      } else {
        for (int i = 0; i < nAssocEl(Ac); i++) {
          for (int j = 0; j < maxAssocEl; j++) {
            if (bin(j,0) == 0) {
              bin(j,0) = assocEl(i,Ac);
              bin(j,1) = 1;
              break;
            } else if (bin(j,0) == assocEl(i,Ac)) {
              bin(j,1) = bin(j,1) + 1;
            }
          } 
        } 
      }
    } 

    for (int j = 0; j < maxAssocEl; j++) {
      if (bin(j,1) == lFa.eNoN) {
        lFa.gE(e) = bin(j,0);
        break;
      }
    } 
  } 

  Array<double> xl(nsd,lM.eNoN), v(nsd,lM.eNoN);

  // Check mesh quality and reset IEN if necessary
  //
  int e = lFa.gE(0);
  for (int i = 0; i < nsd; i++) {
    for (int j = 0; j < lM.eNoN; j++) {
      xl(i,j) = lM.x(i,lM.gIEN(j,e));
    }
  }
  //xl = lM.x(:,lM.gIEN(:,e))

  for (int i = 0; i < nsd; i++) {
    v(i,0) = xl(i,1) - xl(i,0);
    v(i,1) = xl(i,2) - xl(i,1);
    v(i,2) = xl(i,3) - xl(i,2);
  }

  Array<double> v12(3,2);
  for (int i = 0; i < nsd; i++) {
    v12(i,0) = v(i,0); 
    v12(i,1) = v(i,1); 
  }

  v.set_col(3, utils::cross(v12));

  double sum = 0.0;
  for (int i = 0; i < nsd; i++) {
    sum += v(i,2) * v(i,3);
  }
  int sn = utils::sign(sum);
  //sn = SGN(SUM(v(:,3)*v(:,4)))

  if (sn == 1) {
    int a=0, b=1;

    for (int e = 0; e < lM.gnEl; e++) {
      int Ac = lM.gIEN(a,e);
      lM.gIEN(a,e) = lM.gIEN(b,e);
      lM.gIEN(b,e) = Ac;
    }

  } else if (sn == 0) {
    throw std::runtime_error("[set_face_ebc] Surface element " + std::to_string(e+1) + " is distorted.");
  }

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      for (int i = 0; i < nsd; i++) {
        xl(i,a) = lM.x(i,Ac);
      }
    }

    double Jac = all_fun::jacobian(com_mod, nsd, lM.eNoN, xl, lM.Nx.rslice(0));
    //std::cout << "[set_face_ebc] e Jac: " << e+1 << " " << Jac << std::endl;

    if (Jac < 0.0) {
      throw std::runtime_error("[set_face_ebc] Remeshing didn't improve mesh quality.");
    }
  }

}


};
