
#include "all_fun.h"

#include "fsils_api.hpp"
#include "mat_fun.h"
#include "nn.h"
#include "utils.h"

#include <bitset>
#include <math.h>

namespace all_fun {

//--------------
// aspect_ratio
//--------------
//
double aspect_ratio(ComMod& com_mod, const int nDim, const int eNoN, const Array<double>& x)
{
  Array<int> rowM(eNoN,eNoN-1); 
  Array<int> colM(nDim,nDim-1);

  Vector<double> s(eNoN);
  Array<double> Dsub(nDim,nDim); 
  Vector<double> detD(nDim);

  if (nDim == 2) {
    for (int a = 0; a < eNoN; a++) {
      int ap = a + 1;
      if (a == eNoN-1) {
        ap = 0;
      }
      auto x_diff = x.col(a) - x.col(ap);
      s(a) = sqrt( x_diff * x_diff );
    }

  // This is only for tri and tets so for nDim=3 eNoN=4.
  //
  } else if (nDim == 3) {
    Array<int> rowM{ {0, 1, 2}, {0, 1, 3}, {0, 2, 3}, {1, 2, 3} }; 
    Array<int> colM{ {0, 1}, {1, 2}, {2, 0} };

    for (int a = 0; a < eNoN; a++) {
      for (int b = 0; b < nDim; b++) {
        Dsub = 1.0;

        for (int i = 0; i < eNoN-1; i++) {
          int irow = rowM(a,i);
          for (int j = 0; j < nDim-1; j++) {
            int icol = colM(b,j);
            Dsub(i,j) = x(icol,irow);
          } 
        } 
        detD(b) = mat_fun::mat_det(Dsub,nDim);
      } 
      s(a) = 0.5 * sqrt(detD * detD);
    } 
  }

  return s.max() / s.min();
}

//-------
// commu
//-------
//
void commu(const ComMod& com_mod, Vector<double>& U)
{
  if (com_mod.cm.seq()) {
    return;
  }

  if (U.size() != com_mod.lhs.nNo) {
    throw std::runtime_error("COMMU is only specified for vector with size nNo");
  }

  U = mkc(com_mod, U);

  fsi_linear_solver::fsils_commus(com_mod.lhs, U);
  //CALL FSILS_COMMUS(lhs, U)

  mkci(com_mod, U);
  //CALL MKCI(U)
}

//-------
// commu
//-------
//
void commu(const ComMod& com_mod, Array<double>& U)
{
  if (com_mod.cm.seq()) {
    return;
  }

  int m = U.nrows();
  if (U.ncols() != com_mod.lhs.nNo) {
    throw std::runtime_error("COMMU is only specified for vector with size nNo");
  }

  U = mkc(com_mod, U);

  fsi_linear_solver::fsils_commuv(com_mod.lhs, m, U);

  mkci(com_mod, U);

  //U = MKC(U)
  //CALL FSILS_COMMUV(lhs, m, U)
  //CALL MKCI(U)

}

//--------
// domain
//--------
// This function returns the domain that an element of a mesh belongs to
//
// Reproduces 'FUNCTION DOMAIN(lM, iEq, e)' defined in ALLFUN.f.
//
int domain(const ComMod& com_mod, const mshType& lM, const int iEq, const int e)
{
  int domain_id = -1;
  auto& eq = com_mod.eq[iEq];

  // Domain Id of -1 counts for the entire domain
  //
  for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
    domain_id = iDmn;
    if (eq.dmn[iDmn].Id == -1) {
      return domain_id;
    }
  }

  if (lM.eId.size() == 0) { 
    throw std::runtime_error("eId is not allocated");
  }

  for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
    domain_id = iDmn;
    if (utils::btest(lM.eId[e], eq.dmn[iDmn].Id)) {
      return domain_id;
    }
  }

  return domain_id;
}

//-----------
// find_face
//-----------
// Find the face ID and mesh ID based on the face name.
//
void find_face(const std::vector<mshType>& mesh_list, const std::string& faceName, int& iM, int& iFa)
{
  iFa = -1;
  iM = -1;

  for (int mi = 0; mi < mesh_list.size(); mi++) { 
    auto& mesh = mesh_list[mi];
    for (int fi = 0; fi < mesh.nFa; fi++) {
      auto& face = mesh.fa[fi];
      if (face.name == faceName) {
        iM = mi;
        iFa = fi;
        return;
      }
    }
  }

  if (iM == -1) {
    throw std::runtime_error("Can't find face named '" + faceName + "' from defined mesh names."); 
  }
}

//----------
// find_msh
//----------
// Find the mesh ID based on the mesh name.
//
void find_msh(const std::vector<mshType>& mesh_list, const std::string& mesh_name, int& iM)
{
  iM = -1;
  
  for (int i = 0; i < mesh_list.size(); i++) {
    if (mesh_list[i].name == mesh_name) { 
      iM = i;
      break;
    }
  }
}

//--------
// global
//--------
// Reproduces 'FUNCTION GLOBALRV(lM, U)' defined in ALLFUN.f.
//
Array<double> 
global(const ComMod& com_mod, const CmMod& cm_mod, const mshType& lM, const Array<double>& U)
{
  auto& cm = com_mod.cm;

  #define n_debug_global_rv
  #ifdef debug_global_rv 
  DebugMsg dmsg(__func__, cm.idcm());
  dmsg.banner();
  #endif

  int m = U.nrows();
  #ifdef debug_global_rv 
  dmsg << "m: " << m;
  dmsg << "U.ncols(): " << U.ncols();
  dmsg << "lM.gnNo: " << lM.gnNo;
  #endif 

  if (U.ncols() != lM.nNo) {
    throw std::runtime_error("GLOBAL is only specified for array with columns size nNo");
  }

  // [TODO:DaveP] what's going on here?
  //
  if (cm.seq()) {
    // ALLOCATE(GLOBALRV(m,lM.gnNo))
    return U;
  }

  Array<double> result;
  Vector<int> sCount(cm.np()); 
  Vector<int> disp(cm.np());
  Array<double> ienU(m*lM.eNoN, lM.nEl);
  Array<double> gienU;

  if (cm.mas(cm_mod)) {
    gienU.resize(m*lM.eNoN, lM.gnEl); 
    result.resize(m,lM.gnNo);
   } else {
    //ALLOCATE(gienU(0,0), GLOBALRV(0,0))
  }

  for (int e = 0; e < lM.nEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.IEN(a,e);
      Ac = lM.lN(Ac);
      for (int i = 0; i < m; i++) {
        ienU(m*a+i,e) = U(i,Ac);
      }
    }
  }

  int a = lM.eNoN*m;

  for (int i = 0; i < cm.np(); i++) {
    disp(i) = lM.eDist(i)*a;
    sCount(i) = lM.eDist(i+1)*a - disp(i);
  }

  MPI_Gatherv(ienU.data(), lM.nEl*a, cm_mod::mpreal, gienU.data(), sCount.data(), disp.data(), cm_mod::mpreal, cm_mod.master, cm.com());

  // If a slave process return an empty result.
  if (cm.slv(cm_mod)) {
    return result;
  }

  for (int e = 0; e < lM.gnEl; e++) {
    for (int a = 0; a < lM.eNoN; a++) {
      int Ac = lM.gIEN(a,e);
      for (int i = 0; i < m; i++) {
        result(i,Ac) = gienU(m*a+i,e);
      }
    }
  }

  return result;
}

//-------
// integ
//-------
// This routine integrate an equation over a particular domain
//
// Replicates 'FUNCTION vInteg(dId, s, l, u, pFlag)' defined in ALLFUN.f.
//
double integ(const ComMod& com_mod, const CmMod& cm_mod, int dId, const Array<double>& s, int l, int u, bool pFlag)
{
  using namespace consts;

  #define n_debug_integ_v
  #ifdef debug_integ_v
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "dId: " << dId;
  dmsg << "l: " << l;
  dmsg << "pFlag: " << pFlag;
  #endif

  int nNo = s.ncols();
  int tnNo = com_mod.tnNo;
  bool ibFlag = com_mod.ibFlag;

  if (nNo != tnNo) {
    if (ibFlag) {
      if (nNo != com_mod.ib.tnNo) {
          throw std::runtime_error("Incompatible vector size in vInteg");
      }
    } else { 
      throw std::runtime_error("Incompatible vector size in vInteg");
    } 
  }

  #ifdef debug_integ_v
  dmsg << "tnNo: " << tnNo;
  dmsg << "nNo: " << nNo;
  #endif

  bool flag = pFlag; 
  if (l != u) {
    throw std::runtime_error("Incompatible spatial output setting and element type");
  }

  bool isIB = false;
  if (ibFlag) {
    if (nNo == com_mod.ib.tnNo) {
      isIB = true;
    }
  } 

  double result = 0.0; 
  int nsd = com_mod.nsd;
  auto& ib = com_mod.ib;
  fsType fs;

  if (!isIB) {
    for (int iM = 0; iM < com_mod.nMsh; iM++) {
      auto& msh = com_mod.msh[iM];
      int insd = nsd;
      if (msh.lShl) insd = nsd-1;
      if (msh.lFib) insd = 1;

      // Update pressure function space for Taylor-Hood type element

      if (flag && (com_mod.msh[iM].nFs == 2)) {
        fs.nG    = msh.fs[1].nG;
        fs.eType = msh.fs[1].eType;
        fs.lShpF = msh.fs[1].lShpF;
        fs.eNoN  = msh.fs[1].eNoN;

        fs.w.resize(fs.nG); 
        fs.N.resize(fs.eNoN,fs.nG); 
        fs.Nx.resize(nsd,fs.eNoN,fs.nG);

        if (fs.eType != ElementType::NRB) {
          fs.w  = msh.fs[1].w;
          fs.N  = msh.fs[1].N;
          fs.Nx = msh.fs[1].Nx;
        }

      } else { 
        fs.nG    = msh.fs[0].nG;
        fs.eType = msh.fs[0].eType;
        fs.lShpF = msh.fs[0].lShpF;
        fs.eNoN  = msh.fs[0].eNoN;

        fs.w.resize(fs.nG); 
        fs.N.resize(fs.eNoN,fs.nG); 
        fs.Nx.resize(nsd,fs.eNoN,fs.nG);

        if (fs.eType != ElementType::NRB) {
          fs.w  = msh.fs[0].w;
          fs.N  = msh.fs[0].N;
          fs.Nx = msh.fs[0].Nx;
        }
      }
      int eNoN = fs.eNoN;

      Array<double> xl(nsd,eNoN); 
      Array<double> Nxi(insd,eNoN); 
      Array<double> Nx(insd,eNoN); 
      Vector<double> sl(eNoN); 
      Array<double> tmps(nsd,insd);
      Array<double> tmp(nsd,nsd);

      for (int e = 0; e < msh.nEl; e++) {
        if (dId > 0 &&  msh.eId.size() != 0) {
          if (!utils::btest(msh.eId(e),dId)) {
            continue;
          }
        }

        // Updating the shape functions, if this is a NURB
        //
        // [TODO:DaveP] not implemented. 
        //
        if (msh.eType == ElementType::NRB) {
          //CALL NRBNNX(msh(iM), e)
          fs.w  = msh.w;
          fs.N  = msh.N;
          fs.Nx = msh.Nx;
        }

        int ibl = 0;
        for (int a = 0; a < eNoN; a++) { 
          int Ac = msh.IEN(a,e);
          xl.set_col(a, com_mod.x.col(Ac));

          if (com_mod.mvMsh) {
            for (int i = 0; i < nsd; i++) { 
              xl(i,a) += com_mod.Do(i+nsd+1,Ac);
            }
          }

          if (l == u) {
            sl(a) = s(l,Ac);
          } else { 
            auto rows = s.col(Ac, {l,u});
            sl(a) = sqrt(utils::norm(rows));
          }
          ibl = ibl + com_mod.iblank(Ac);
        }

        if (ibl == eNoN) {
          continue;
        }

        double Jac = 0.0;

        for (int g = 0; g < fs.nG; g++) {
          Nxi = fs.Nx.slice(g);

          if (g == 0 || !fs.lShpF) {
            if (msh.lShl) {
              Vector<double> nV(nsd);
              nn::gnns(nsd, eNoN, Nxi, xl, nV, tmps, tmps);
              Jac  = sqrt(utils::norm(nV));
            } else { 
              nn::gnn(eNoN, nsd, insd, Nxi, xl, Nx, Jac, tmp);
            }
          }

          if (utils::is_zero(Jac)) {
            throw std::runtime_error("Jac < 0 for element: " + std::to_string(e) + ")");
          }
          double sHat = 0.0;

          for (int a = 0; a < eNoN; a++) {
            int Ac = msh.IEN(a,e);
            sHat = sHat + sl(a)*fs.N(a,g);
          }
          result += fs.w(g) * Jac * sHat;
        }
      }
    }

  } else { 

    for (int iM = 0; iM < ib.nMsh; iM++) {
      int eNoN = ib.msh[iM].eNoN;
      int insd = nsd;

      Array<double> xl(nsd,eNoN); 
      Array<double> Nxi(insd,eNoN); 
      Array<double> Nx(insd,eNoN); 
      Vector<double> sl(eNoN); 
      Array<double> tmps(nsd,insd);
      Array<double> tmp(nsd,nsd);

      for (int e = 0; e < ib.msh[iM].nEl; e++) {
        if (dId > 0 && ib.msh[iM].eId.size() != 0) {
          if (!utils::btest(ib.msh[iM].eId(e),dId)) {
            continue;
          }
        }

        // Updating the shape functions, if this is a NURB
        // [TODO:DaveP] not implemented. 
        if (ib.msh[iM].eType == ElementType::NRB) {
          //CALL NRBNNX(ib.msh(iM), e)
        }

        for (int a = 0; a < eNoN; a++) {
          int Ac = ib.msh[iM].IEN(a,e);
          for (int i = 0; i < nsd; i++) { 
            xl(i,a) = ib.x(i,Ac) + ib.Ubo(i,Ac);
          }

          if (l == u) {
            sl(a) = s(l,Ac);
          } else { 
            auto rows = s.col(Ac, {l,u});
            sl(a) = sqrt(utils::norm(rows));
          }
        }

        for (int g = 0; g < ib.msh[iM].nG; g++) {
          double Jac = 0.0;
          Nxi = ib.msh[iM].Nx.slice(g);
          if (g == 0 ||  !ib.msh[iM].lShpF) {
            nn::gnn(eNoN, nsd, insd, Nxi, xl, Nx, Jac, tmp);
          }
          double sHat = 0.0;

          for (int a = 0; a < eNoN; a++) {
            int Ac = ib.msh[iM].IEN(a,e);
            sHat = sHat + sl(a)*ib.msh[iM].N(a,g);
          }
          result += ib.msh[iM].w(g)*Jac*sHat;
        }
      }
    }
  }

  if (com_mod.cm.seq() || isIB) {
    return result;
  }

  result = com_mod.cm.reduce(cm_mod, result);
  return result;
}

//-------
// integ
//-------
// This routine integrate s over the surface faId.
//
// Reproduces 'FUNCTION IntegS(lFa, s, pflag)'.
//
double integ(const ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const Vector<double>& s, bool pFlag)
{
  using namespace consts;
  #define n_debug_integ_s
  #ifdef debug_integ_s
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "lFa.iM: " << lFa.iM+1;
  dmsg << "lFa.name: " << lFa.name;
  dmsg << "lFa.eType: " << lFa.eType;
  #endif 

  bool flag = pFlag; 
  int nsd = com_mod.nsd;
  int insd = nsd - 1;

  if (com_mod.msh[lFa.iM].lShl) {
    insd = insd - 1;
  }

  if (com_mod.msh[lFa.iM].lFib) {
    insd = 0;
  }

  int nNo = s.size();
  #ifdef debug_integ_s
  dmsg << "nNo: " << nNo;
  dmsg << "insd: " << insd;
  dmsg << "flag: " << flag;
  #endif

  if (nNo != com_mod.tnNo) {
    if (com_mod.ibFlag) {
      if (nNo != com_mod.ib.tnNo) {
        throw std::runtime_error("Incompatible vector size in Integ");
      }
    } else {
      throw std::runtime_error("Incompatible vector size in vInteg");
    }
  }

  bool isIB = false;

  if (com_mod.ibFlag) {
    if (nNo == com_mod.ib.tnNo) {
      isIB = true;
    }
  }

  // Update pressure function space for Taylor-Hood element
  //
  fsType fs;

  if (flag) {
    if (lFa.nFs != 2) {
      throw std::runtime_error("Incompatible boundary integral function call and face element type");
    }

    fs.nG    = lFa.fs[1].nG;
    fs.eType = lFa.fs[1].eType;
    fs.lShpF = lFa.fs[1].lShpF;
    fs.eNoN  = lFa.fs[1].eNoN;

    fs.w.resize(fs.nG); 
    fs.N.resize(fs.eNoN,fs.nG); 
    fs.Nx.resize(insd,fs.eNoN,fs.nG);

    if (fs.eType != ElementType::NRB) {
      fs.w  = lFa.fs[1].w;
      fs.N  = lFa.fs[1].N;
      fs.Nx = lFa.fs[1].Nx;
    }
  } else {
    fs.nG    = lFa.fs[0].nG;
    fs.eType = lFa.fs[0].eType;
    fs.lShpF = lFa.fs[0].lShpF;
    fs.eNoN  = lFa.fs[0].eNoN;

    fs.w.resize(fs.nG); 
    fs.N.resize(fs.eNoN,fs.nG); 
    fs.Nx.resize(insd,fs.eNoN,fs.nG);

    if (fs.eType != ElementType::NRB) {
      fs.w  = lFa.fs[0].w;
      fs.N  = lFa.fs[0].N;
      fs.Nx = lFa.fs[0].Nx;
    }
  }

  #ifdef debug_integ_s
  dmsg << "fs.nG: " << fs.nG;
  dmsg << "fs.eType: " << fs.eType;
  dmsg << "fs.w: " << fs.w;
  #endif
  double result = 0.0;

  for (int e = 0; e < lFa.nEl; e++) {
    // [TODO:DaveP] not implemented.
    if (lFa.eType == ElementType::NRB) {
      if (!isIB) {
        //CALL NRBNNXB(msh(lFa.iM), lFa, e)
      } else {
        //CALL NRBNNXB(ib.msh(lFa.iM), lFa, e)
      }
      fs.w  = lFa.w;
      fs.N  = lFa.N;
      fs.Nx = lFa.Nx;
    }

    for (int g = 0; g < fs.nG; g++) {
      Vector<double> n(nsd);
      if (!isIB) {
        auto Nx = fs.Nx.slice(g);
        nn::gnnb(com_mod, lFa, e, g, nsd, insd, fs.eNoN, Nx, n);
      }

      double Jac = sqrt(utils::norm(n));

      // Calculating the function value
      double sHat = 0.0;
      for (int a = 0; a < fs.eNoN; a++) {
        int Ac = lFa.IEN(a,e);
        sHat = sHat + s(Ac)*fs.N(a,g);
      }

      // Now integrating
      result = result + Jac*fs.w(g)*sHat;
     }
  }

  if (com_mod.cm.seq() || isIB) {
    return result; 
  }

  result = com_mod.cm.reduce(cm_mod, result);
  return result; 
}

//-------
// integ
//-------
// This routine integrate s over the surface faId. 
//
// Reproduces 'FUNCTION IntegV(lFa, s)'
//
double integ(const ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const Array<double>& s)
{
  using namespace consts;

  #define n_debug_integ_V
  #ifdef debug_integ_V
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  #endif

  auto& cm = com_mod.cm;
  int nsd = com_mod.nsd;
  int insd = nsd - 1;
  int tnNo = com_mod.tnNo;

  if (s.nrows() != nsd) {
    throw std::runtime_error("Incompatible vector size in integ");
  }

  int nNo = s.ncols();
  #ifdef debug_integ_V
  dmsg << "nNo: " << nNo;
  dmsg << "nsd: " << nsd;
  #endif

  if (nNo != tnNo) {
    if (com_mod.ibFlag) {
      if (nNo != com_mod.ib.tnNo) {
        throw std::runtime_error("Incompatible vector size in integ");
      }
    } else {
      throw std::runtime_error("Incompatible vector size in integ");
    }
  }

  bool isIB = false;
  if (com_mod.ibFlag) {
    if (nNo ==  com_mod.ib.tnNo) {
      isIB = true;
    }
  }

  double result =  0.0;

  for (int e = 0; e < lFa.nEl; e++) {
    //dmsg << "----- e " << e+1 << " -----";
    //  Updating the shape functions, if this is a NURB
    if (lFa.eType == ElementType::NRB) {
      if (!isIB) {
         //CALL NRBNNXB(msh(lFa.iM), lFa, e)
      } else {
         //CALL NRBNNXB(ib.msh(lFa.iM), lFa, e)
      }
    }

    for (int g = 0; g < lFa.nG; g++) {
      //dmsg << ">>> g: " << g+1;
      Vector<double> n(nsd);
      if (!isIB) {
        auto Nx = lFa.Nx.slice(g);
        nn::gnnb(com_mod, lFa, e, g, nsd, nsd-1, lFa.eNoN, Nx, n);
        //CALL GNNB(lFa, e, g, nsd-1, lFa.eNoN, lFa.Nx(:,:,g), n)
      } else {
        //CALL GNNIB(lFa, e, g, n)
      }

      //  Calculating the function value
      //
      double sHat = 0.0;

      for (int a = 0; a < lFa.eNoN; a++) {
        int Ac = lFa.IEN(a,e);
        for (int i = 0; i < nsd; i++) {
          sHat = sHat + lFa.N(a,g) * s(i,Ac) * n(i);
          //dmsg << "s(i,Ac): " << s(i,Ac);
        }
      }

      //  Now integrating
      //dmsg << "sHat: " << sHat;
      //dmsg << "lFa.w(g): " << lFa.w(g);
      result = result + lFa.w(g) * sHat;
    }
  }

  if (cm.seq() || isIB) {
    return result; 
  }

  result = cm.reduce(cm_mod, result);

  return result; 
}

//-------
// integ
//-------
// This routine integrate s(l:u,:) over the surface faId.
//
// Note that 'l' seems to be a length and 'uo' an offset. 'l' should never be 0.
//
// Reproduces 'FUNCTION IntegG(lFa, s, l, uo, THflag)'.
//
double integ(const ComMod& com_mod, const CmMod& cm_mod, const faceType& lFa, const Array<double>& s, const int l, int uo, bool THflag)
{
  using namespace consts;

  #define n_debug_integ_g
  #ifdef debug_integ_g
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "l: " << l;
  dmsg << "uo: " << uo;
  #endif

  auto& cm = com_mod.cm;
  int nsd = com_mod.nsd;
  int insd = nsd - 1;
  int tnNo = com_mod.tnNo;

  int u = l;
  if (uo != -1) { 
    u = uo;
  }

  bool flag = THflag; 
  int nNo = s.ncols();

  if (nNo != tnNo) {
    if (com_mod.ibFlag) {
      if (nNo != com_mod.tnNo) { 
        throw std::runtime_error("Incompatible vector size in integ");
      }
    } else {
      throw std::runtime_error("Incompatible vector size in integ");
    }
  }

  double result = 0.0; 

  if (u-l+1 == nsd) {
     Array<double> vec(nsd,nNo);
     for (int a = 0; a < nNo; a++) {
       for (int i = 0; i < nsd; i++) {
         vec(i,a) = s(i+l-1,a);                 
       }
     }
     result = integ(com_mod, cm_mod, lFa, vec);

  } else if (l == u) {
     Vector<double> sclr(nNo);
     for (int a = 0; a < nNo; a++) {
        sclr(a) = s(l-1,a);
     }
     result = integ(com_mod, cm_mod, lFa, sclr, flag);
  } else {
    throw std::runtime_error("Unexpected dof in integ");
  }

  return result; 
}

//-----------
// is_domain
//-----------
//
bool is_domain(const ComMod& com_mod, const eqType& eq, const int node, const consts::EquationType phys)
{
  bool result = false;

  // Single domain is assumed and we only need to check that
  //
  if (eq.nDmn == 1) {
    if (eq.dmn[0].phys == phys) {
      result = true;
    } 

  // Domain partition is expected
  //
  } else { 
    if (com_mod.dmnId.size() == 0) {
      throw std::runtime_error("Domain partitioning info is not provided.");
    }
    //  IF (.NOT.ALLOCATED(dmnId)) err = "Domain partitioning info "// 2      "is not provided"

    for (int iDmn = 0; iDmn < eq.nDmn; iDmn++) {
      if (eq.dmn[iDmn].phys == phys) {
        if (utils::btest(com_mod.dmnId(node),eq.dmn[iDmn].Id)) {
          result = true;
          break;
        }
      }
    }
  }

  return result;
}

//----------
// jacobian
//----------
// Computes the JACOBIAN of an element.
//
double jacobian(ComMod& com_mod, const int nDim, const int eNoN, const Array<double>& x, const Array<double>&Nxi)
{
  double Jac = 0.0;
  Array<double> xXi(nDim,nDim);

  for (int a = 0; a < eNoN; a++) {
    xXi.set_col(0, xXi.col(0) + x.col(a) * Nxi(0,a));
    xXi.set_col(1, xXi.col(1) + x.col(a) * Nxi(1,a));
    // xXi(:,1) = xXi(:,1) + (x(:,a) * Nxi(1,a))
    // xXi(:,2) = xXi(:,2) + (x(:,a) * Nxi(2,a))

    if (nDim == 3) {
      xXi.set_col(2, xXi.col(2) + x.col(a) * Nxi(2,a));
      //xXi(:,3) = xXi(:,3) + (x(:,a) * Nxi(3,a))
    }
  }

  Jac = mat_fun::mat_det(xXi,nDim);
  //Jac = MAT_DET(xXi,nDim)

  return Jac;
}

//-------
// local
//-------
//
Vector<int> 
local(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, Vector<int>& U)
{
  Vector<int> local_vector;

  if (com_mod.ltg.size() == 0) {
    throw std::runtime_error("ltg is not set yet");
  }

  if (cm.mas(cm_mod)) {
    if (U.size() != com_mod.gtnNo) {
      throw std::runtime_error("local_is is only specified for vector with size gtnNo");
    }
  }

  if (cm.seq()) {
    local_vector.resize(com_mod.gtnNo);
    local_vector = U;
    return local_vector; 
  }

  local_vector.resize(com_mod.tnNo); 
  Vector<int> tmpU(com_mod.gtnNo);

  if (cm.mas(cm_mod)) {
    tmpU = U;
  }

  cm.bcast(cm_mod, tmpU);

  for (int a = 0; a < com_mod.tnNo; a++) {
    int Ac = com_mod.ltg[a];
    local_vector[a] = tmpU[Ac];
  }

  return local_vector;
}

//-------
// local
//-------
//
Array<double> 
local(const ComMod& com_mod, const CmMod& cm_mod, const cmType& cm, Array<double>& U)
{
  //int task_id = cm.idcm();
  //std::string msg_prefix = std::string("[local_rv:") + std::to_string(task_id) + "] ";
  //dmsg;
  //dmsg << "========== local_rv ==========";

  if (com_mod.ltg.size() == 0) {
    throw std::runtime_error("ltg is not set yet");
  }

  Array<double> local_array; 
  int m;

  if (cm.mas(cm_mod)) {
    m = U.nrows();
    if (U.ncols() != com_mod.gtnNo) {
      throw std::runtime_error("local_rv is only specified for vector with size gtnNo");
    }
  }

  if (cm.seq()) {
    local_array.resize(m, com_mod.gtnNo);
    local_array = U;
    return local_array; 
  }

  cm.bcast(cm_mod, &m);

  local_array.resize(m, com_mod.tnNo); 
  Vector<double> tmpU(m * com_mod.gtnNo);
  //ALLOCATE(LOCALRV(m,tnNo), tmpU(m*gtnNo))

  if (cm.mas(cm_mod)) {
    for (int a = 0; a < com_mod.gtnNo; a++) {
      int s = m * a;
      int e = m * (a + 1);
      for (int i = 0; i < U.nrows(); i++) {
        tmpU(i+s) = U(i, a);
      }
      // tmpU(s:e) = U(:,a)
    }
  }

  cm.bcast(cm_mod, tmpU);

  //dmsg << "Compute local_rv ";

  for (int a = 0; a < com_mod.tnNo; a++) {
    int Ac = com_mod.ltg[a];
    int s = m * Ac;
    int e = m * (Ac+1);
    //if (a < 5) dmsg << "localrv Ac: " << Ac;
    // int s  = m * (Ac-1) + 1;
    // int e  = m * Ac;
    // LOCALRV(:,a) = tmpU(s:e)
    // if (a < 5) dmsg << "localrv s and e: " << s << " " << e;
    //if (a < 5) dmsg << "localrv a: " << a;
    for (int i = 0; i < m; i++) {
      //if (a < 5) dmsg << "localrv tmpU(i+s): " << tmpU(i+s);
      local_array(i,a) = tmpU(i+s);
      //if (a < 5) dmsg << "localrv (i,a): " << local_array(i,a);
    }
  }

  return local_array;
}

//-----
// mkc
//-----
//
Vector<double> 
mkc(const ComMod& com_mod, Vector<double>& U)
{
  if (U.size() != com_mod.lhs.nNo) {
    throw std::runtime_error("MKC is only specified for vector with size nNo");
  }

  Vector<double> result(com_mod.lhs.nNo);

  if (com_mod.cm.seq()) {
    result = U;
  } else {  
    for (int a = 0; a < com_mod.lhs.nNo; a++) {
      result(com_mod.lhs.map(a)) = U(a);
    }
  }

  return result;
}

Array<double> 
mkc(const ComMod& com_mod, Array<double>& U)
{
  int m = U.nrows();
  if (U.ncols() != com_mod.lhs.nNo) {
    throw std::runtime_error("MKC is only specified for vector with size nNo");
  }

  Array<double> result(m, com_mod.lhs.nNo);

  if (com_mod.cm.seq()) {
    result = U;
  } else {
    for (int a = 0; a < com_mod.lhs.nNo; a++) {
      for (int i = 0; i < m; i++) {
        result(i, com_mod.lhs.map(a)) = U(i, a);
      }
    }
  }

  return result;
}

//------
// mkci
//------
//
void mkci(const ComMod& com_mod, Vector<double>& U)
{
  if (com_mod.cm.seq()) {
    return;
  }

  if (U.size() != com_mod.lhs.nNo) {
    throw std::runtime_error("MKC is only specified for vector with size nNo");
  }

  Vector<double> tmp(com_mod.lhs.nNo);
  tmp = U;

  for (int a = 0; a < com_mod.lhs.nNo; a++) {
    U(a) = tmp(com_mod.lhs.map(a));
  }
}

void mkci(const ComMod& com_mod, Array<double>& U)
{
  if (com_mod.cm.seq()) {
    return;
  }

  auto& lhs = com_mod.lhs;
  int nNo = lhs.nNo;
  int m = U.nrows();

  if (U.ncols() != nNo) {
    throw std::runtime_error("MKC is only specified for vector with size nNo");
  }

  Array<double> tmp(m, nNo);
  tmp = U;

  for (int a = 0; a < nNo; a++) {
    for (int i = 0; i < m; i++) {
      U(i,a) = tmp(i,lhs.map(a));
    } 
  }
}

//------------
// set_dmn_id
//------------
// Set domain ID to a given number for the entire or a range of elements in a mesh.
//
void set_dmn_id(mshType& mesh, const int iDmn, const int ifirst, const int ilast)
{
  int first = 0;
  if (ifirst != consts::int_inf) {
    first = ifirst;
  }
  if ((first < 0) || (first > mesh.gnEl-1)) {
    throw std::runtime_error("Setting domain ID with a start element range " + std::to_string(first) + ".");
  }

  int last = mesh.gnEl - 1;
  if (ilast != consts::int_inf) {
    last = ilast;
  }
  if ((last < first) || (last > mesh.gnEl-1)) {
    throw std::runtime_error("Setting domain ID with a end element range " + std::to_string(last) + ".");
  }

  if (mesh.eId.size() == 0) {
    mesh.eId = Vector<int>(mesh.gnEl);
  }

  // Set the iDimn'th bit for each element ID.
  for (int e = first; e <= last; e++) {
    mesh.eId[e] |= 1UL << iDmn;
  }
}

//----------
// skewness
//----------
// Computes the Skewness of an element.
//
double skewness(ComMod& com_mod, const int nDim, const int eNoN, const Array<double>& x)
{
  Vector<double> coeff;

  if (nDim == 2) {
    coeff.set_values({1.0, -1.0, 1.0, -1.0});
  } else {
    coeff.set_values({1.0, 1.0, -1.0, 1.0, 1.0});
  }

  Array<double> Dmat(eNoN,nDim+2);
  Dmat = 1.0;
  Array<double> Dsub(eNoN,nDim+1);

  for (int a = 0; a < eNoN; a++) {
    auto col = x.col(a);
    Dmat(a,0) = col * col;

    for (int i = 0; i < nDim; i++) {
      Dmat(a, i+1) = x(i,a);
    }
  }

  Vector<double> detD(nDim+2);

  for (int j = 0; j < nDim+2; j++) {
    int cnt = 0;
    for (int i = 0; i < nDim+2; i++) {
      if (i == j) {
        continue;
      }
      Dsub.set_col(cnt, Dmat.col(i));
      cnt = cnt + 1;
    }
    detD(j) = coeff(j) * mat_fun::mat_det(Dsub,nDim+1);
  }

  double circumRad = 0.0; 
  for (int i = 0; i < nDim; i++) {
    circumRad += detD(i+1)*detD(i+1);
  }
  circumRad = sqrt(circumRad - 4.0*detD(0)*detD(nDim+1)) / (2.0*fabs(detD(0)));

  double integ_eq, integ_el;

  if (nDim == 2) {
    integ_eq = 0.25 * sqrt(27.0) * circumRad*circumRad;
    integ_el = 0.5 * fabs(detD(0));
  } else if (nDim == 3) {
    integ_eq = 80.0 * pow(circumRad,3.0) / sqrt(243.0);
    integ_el = fabs(detD(0)) / 6.0;
  }

  return fabs(integ_eq - integ_el) / integ_eq;
}

//------------
// split_jobs
//------------
// Spliting "m" jobs between "n" workers. "b" contains amount of jobs
// and "A" will store the distribution of jobs
//
// A(nMsh, num_proc)
//
// b(nMsh)
//
// n = nMsh
//
// m = num_proc
//
// Replicates 'RECURSIVE SUBROUTINE SPLITJOBS(m, n, A, b)' in ALLFUN.f.
//
void split_jobs(int tid, int m, int n, Array<double>& A, Vector<double>& b)
{
  #define ndebug_split_jobs
  #ifdef debug_split_jobs
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "m: " << m;
  dmsg << "n: " << n;
  #endif

  if ((m <= 0) || (n <= 0)) {
    return;
  }

  // Single worker, but multiple jobs.
  //
  if (n == 1) {
    #ifdef debug_split_jobs
    dmsg << "n == 1  ";
    #endif
    int j = 0;
    for (int i = 0; i < m; i++) {
      A(i,j) = b[i];
      #ifdef debug_split_jobs
      dmsg << "n=1 A(i,0): " << A(i,0);
      #endif
    }
    return;
  }

  // Multiple workers, but a single job.
  //
  if (m == 1) {
    #ifdef debug_split_jobs
    dmsg << "m == 1  ";
    #endif
    int i = 0;
    for (int j = 0; j < n; j++) {
      A(i,j) = b[i] / static_cast<double>(n); 
      #ifdef debug_split_jobs
      dmsg << "m=1 A(0,j): " << A(i,j);
      #endif
    }
    return;
  }

  // Multiple workers and multiple jobs
  // This is the initial guess for nl, nr
  int nl  = n / 2;
  int nr  = n - nl;

  // This is the total amount of work
  #ifdef debug_split_jobs
  dmsg << "b: ";
  for (int i = 0; i < m; i++) {
    dmsg << "b( " << i << "): " << b[i] ;
  }
  #endif
  double sb  = b.sum();

  // The work that suppose to be done by "l"
  double sbl = sb * static_cast<double>(nl) / static_cast<double>(n); 
  #ifdef debug_split_jobs
  dmsg << "nl: " << nl;
  dmsg << "nr: " << nr;
  dmsg << "sb: " << sb;
  dmsg << "sbl: " << sbl;
  #endif

  double sl = 0.0;
  int ival = 0;
  for (int i = 0; i < m; i++) {
    if (sl + b[i] > sbl) {
      ival = i;
      break;
    }
    sl = sl + b[i];
  }

  double optsl = fabs(sl - sbl);
  int ml = ival;
  #ifdef debug_split_jobs
  dmsg;
  dmsg << "optsl: " << optsl;
  dmsg << "ml: " << ml;
  dmsg;
  dmsg << "Set j ... ";
  #endif

  int j = -1;

  for (int i = ml; i < m; i++) {
    #ifdef debug_split_jobs
    dmsg << "---- i " << i << " ----";
    #endif
    if (fabs(sl + b[i] - sbl) < optsl) {
      j = i;
      optsl = fabs(sl + b[i] - sbl);
      #ifdef debug_split_jobs
      dmsg << "j: " << j;
      dmsg << "optsl: " << optsl;
      #endif
    }
  }

  if (j !=  -1) {
    ml = ml + 1;
  }

  int mr = m - ml;
  Vector<double> bl(ml); 
  Vector<double> br(mr);

  #ifdef debug_split_jobs
  dmsg;
  dmsg << "j: " << j;
  dmsg << "mr: " << mr;
  dmsg << "ml: " << ml;
  #endif

  if (j != -1) {
    for (int i = 0; i < ml; i++) {
      bl[i] = b[i];
    }
    bl[ml-1] = b[j];
    int k = 0;
    for (int i = ml-1; i < m; i++) {
      if (i == j) {
        continue; 
      }
      br[k] = b[i];
      k = k + 1;
    }
  } else { 
    for (int i = 0; i < ml; i++) {
      bl[i] = b[i];
    }
    for (int i = 0; i < mr; i++) {
      br[i] = b[i+ml-1];
    }
  }

  nl = round(static_cast<double>(n) * bl.sum() / sb);
  if (nl == 0) {
    nl = 1;
  }
  if (nl == n) { 
    nl = n - 1;
  }
  nr = n - nl;

  #ifdef debug_split_jobs
  dmsg << "nl: " << nl;
  dmsg << "nr: " << nr;
  #endif

  Array<double> Al(ml,nl); 
  split_jobs(tid, ml, nl, Al, bl);

  Array<double> Ar(mr,nr);
  split_jobs(tid, mr, nr, Ar, br);

  #ifdef debug_split_jobs
  dmsg << "Allocate Al: " << ml << " x " << nl;
  dmsg << "Allocate Ar: " << mr << " x " << nr;
  dmsg << "Allocate A: " << A.nrows() << " x " << A.ncols();
  #endif

  A = 0.0;

  if (j != -1) {
    #ifdef debug_split_jobs
    dmsg << "";
    dmsg << "set A from Al ... ml: " << ml;
    #endif
    for (int i = 0; i < ml-1; i++) {
      auto Al_row = Al.row(i);
      A.set_row(i, Al_row);
      #ifdef debug_split_jobs
      for (int ii = 0; ii < n; ii++) {
        dmsg << " Al set A(i, " << i << ", " << ii << "): " << A(i,ii);
      }
      #endif
    }

    A.set_row(j, Al.row(ml-1));
    #ifdef debug_split_jobs
    dmsg << "----------";
    for (int ii = 0; ii < n; ii++) {
      dmsg << " set from Al:  A(i, " << j << ", " << ii << "): " << A(j,ii);
    }
    #endif

    int k = 0;
    for (int i = ml-1; i < m; i++) {
      if (i == j) {
        continue; 
      }
      A.set_row(i, nl, Ar.row(k));
      k = k + 1;
    }
  } else { 
    for (int i = 0; i < ml; i++) {
      for (int j = 0; j < nl; j++) {
        A(i, j) = Al(i, j);
        A(i+ml, j+nl) = Al(i, j);
      }
    }
  }

}


};
