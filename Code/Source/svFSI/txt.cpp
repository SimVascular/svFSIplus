
// This is to calculate average/flux of variables at the boundaries
// and print them into the "B_" files.

#include "txt.h"

#include "all_fun.h"
#include "consts.h"
#include "post.h"
#include "set_bc.h"
#include "utils.h"
#include <math.h>

namespace txt_ns {

//-------
// cctxt
//-------
// This is to check/create the txt file
//
void cctxt(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const std::array<std::string,2>& fName, const std::vector<bool>& wtn)
{
  int fid = 1;

  if (com_mod.cm.slv(cm_mod)) {
    return;
  }

  // i=1 are the boundary values and i=2 are volume values
  //
  for (int i = 0; i < 2; i++) {
    if (!wtn[i]) {
      continue;
    }

    bool flag = false;
    auto file_name = fName[i];
    if (FILE *fp = fopen(file_name.c_str(), "r")) {
      flag = true;
      fclose(fp);
    }

    // [NOTE] What is this all about?
    if (com_mod.cTS != 0 && flag) {
      //CALL TRIMFILE(cTS+3,fName(i))
      continue;
    }

    auto fp = fopen(file_name.c_str(), "w");

    if (i == 0) {
      for (int iM = 0; iM < com_mod.nMsh; iM++) {
        auto& msh = com_mod.msh[iM];

        for (int iFa = 0; iFa < msh.nFa; iFa++) {
          auto& fa = msh.fa[iFa];
          auto stmp = fa.name;
          fprintf(fp, " %s ", stmp.c_str());
        }
      }

      fprintf(fp, "\n");

      for (int iM = 0; iM < com_mod.nMsh; iM++) {
        auto& msh = com_mod.msh[iM];
        for (int iFa = 0; iFa < msh.nFa; iFa++) {
          auto& fa = msh.fa[iFa];
          fprintf(fp, " %4.3e ", fa.area);
        }
      }

    } else {
      for (int iDmn = 0; iDmn < lEq.nDmn; iDmn++) {
        std::string stmp;
        if (lEq.dmn[iDmn].Id == -1) {
          stmp = "ENTIRE";
        }
      }

      for (int iDmn = 0; iDmn < lEq.nDmn; iDmn++) {
        std::string stmp;
      }
    }

    fprintf(fp, "\n");
    fprintf(fp, "\n");
    fclose(fp);
  }

}

//-----
// txt
//-----
// [NOTE] This is not fully implemented.
//
void txt(Simulation* simulation, const bool flag) 
{
  using namespace consts;
  using namespace utils;

  auto& com_mod = simulation->com_mod;
  auto& cm_mod = simulation->cm_mod;
  auto& chnl_mod = simulation->chnl_mod;
  int nsd = com_mod.nsd;
  int tnNo = com_mod.tnNo;
  auto& cplBC = com_mod.cplBC;

  #define n_debug_txt
  #ifdef debug_txt
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "maxNSD: " << maxNSD;
  dmsg << "cplBC.coupled: " << cplBC.coupled;
  #endif

  int fid = 1;

  // Writing data related to cplBC
  //
  if (!com_mod.resetSim) {
    if (cplBC.coupled) {
      bool ltmp = false;

      if (!flag) {
        if (cplBC.useGenBC) {
          set_bc::genBC_Integ_X(com_mod, cm_mod, "L");

        } else {

          for (auto& bc : com_mod.eq[0].bc) {
            if (utils::btest(bc.bType, iBC_RCR)) {
              ltmp = true;
              break;
            }
          }

          set_bc::cplBC_Integ_X(com_mod, cm_mod, ltmp);
        }
      }

      if (com_mod.cm.mas(cm_mod) && !cplBC.useGenBC) {
        if (flag) {
          if (com_mod.cTS == 0 || !ltmp) {
            //OPEN(fid, FILE=cplBC.saveName)
            //CLOSE(fid, STATUS='DELETE')
          } else {
            //CALL TRIMFILE(cTS,cplBC.saveName)
          }
        } else {
          //OPEN(fid, FILE=cplBC.saveName, POSITION='APPEND')
          //WRITE(fid,'(ES14.6E2)',ADVANCE='NO') cplBC.xp(1)

          for (int i = 0; i < cplBC.nX; i++) {
            // WRITE(fid,'(2(1X,ES14.6E2))',ADVANCE='NO') cplBC.xn(i), cplBC.fa(i).y
          }

          for (int i = 1; i < cplBC.nXp; i++) {
            //WRITE(fid,'(ES14.6E2)',ADVANCE='NO') cplBC.xp(i)
          }

          //WRITE(fid,*)
          //CLOSE(fid)
        }
      }
    }
  } // resetSim

  // Process eq outputs
  //
  //dmsg << "Process eq outputs ... ";
  //dmsg << "com_mod.nEq: " << com_mod.nEq;

  for (int iEq = 0; iEq < com_mod.nEq; iEq++) {
    #ifdef debug_txt
    dmsg << "---------- iEq " << iEq+1 << " ----------";
    #endif
    Array<double> tmpV(maxNSD,tnNo);
    auto& eq = com_mod.eq[iEq];

    for (int iOut = 0; iOut < eq.nOutput; iOut++) {
      auto& output = eq.output[iOut];
      // Don't write it when it doesn't suppose to be written
      #ifdef debug_txt
      dmsg;
      dmsg << "----- iOut " << iOut+1 << " -----";
      #endif

      std::vector<bool> wtn;
      for (int i = 1; i < 3; i++) {
        wtn.push_back( eq.output[iOut].wtn[i] );
      }

      if (std::count(wtn.begin(), wtn.end(), false) == wtn.size()) {
        continue;
      } 

      // eq.s - Pointer to start of unknown Yo(:,s:e)
      // eq.output[iOut].l - Length of the outputed variable
      // eq.output[iOut].o - Offset from the first index
      //
      int l = eq.output[iOut].l;
      int s = eq.s + eq.output[iOut].o;  
      int e = s + l - 1;
      #ifdef debug_txt
      dmsg << "l: " << l;
      dmsg << "s: " << s;
      dmsg << "e: " << e;
      #endif

      auto oGrp = eq.output[iOut].grp;
      bool div = true;
      bool pflag = false;
      #ifdef debug_txt
      dmsg << oGrp;
      dmsg << "oGrp: " << oGrp;
      dmsg << "name: " << eq.output[iOut].name;
      #endif

      switch (oGrp) {
        case OutputType::outGrp_NA:
          throw std::runtime_error("NA outGrp in TXT");
        break;

        case OutputType::outGrp_A:
          for (int j = 0; j < tnNo; j++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,j) = com_mod.An(i+s,j);
            }
          }
        break;

        case OutputType::outGrp_Y:
          for (int j = 0; j < tnNo; j++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,j) = com_mod.Yn(i+s,j);
            }
          }

          if (l == 1) {
            pflag = true;
          }
        break;

        case OutputType::outGrp_D:
          for (int j = 0; j < tnNo; j++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,j) = com_mod.Dn(i+s,j);
            }
          }
        break;

        case OutputType::outGrp_WSS: 
        case OutputType::outGrp_vort: 
        case OutputType::outGrp_trac:
          post::all_post(simulation, tmpV, com_mod.Yn, com_mod.Dn, oGrp, iEq);
          for (int a = 0; a < tnNo; a++) {
            auto vec = tmpV.col(a, {0,nsd-1});
            tmpV(0,a) = sqrt(norm(vec));
          }
          l = 1;
        break;

        case OutputType::outGrp_eFlx: 
        case OutputType::outGrp_hFlx: 
        case OutputType::outGrp_divV: 
        case OutputType::outGrp_J: 
        case OutputType::outGrp_mises:
          post::all_post(simulation, tmpV, com_mod.Yn, com_mod.Dn, oGrp, iEq);
        break;

        case OutputType::outGrp_absV:
          for (int a = 0; a < tnNo; a++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,a) = com_mod.Yn(i,a) - com_mod.Yn(i+nsd+1,a);
            }
          }
        break;

        case OutputType::outGrp_I:
          div = false;
          for (int a = 0; a < tnNo; a++) {
            for (int i = 0; i < l; i++) {
              tmpV(i,a) = 1.0;
            }
          }
        break;

        default:
          throw std::runtime_error("Undefined output");
      } 

      // Write average/flux of variables at the boundaries
      //
      #ifdef debug_txt
      dmsg << "Write average/flux of variables ...";
      #endif
      const auto& appPath = chnl_mod.appPath;
      std::array<std::string,2> fName; 
      fName[0] = eq.sym + "_" + output.name;
      fName[1] = fName[0];

      if (l == nsd) {
        fName[0] = appPath + "B_" + fName[0] + "_flux_cpp.txt";
      } else {
        fName[0] = appPath + "B_" + fName[0] + "_average_cpp.txt";
      }

      fName[1] = appPath + "V_" + fName[1] + "_flux_cpp.txt";
      #ifdef debug_txt
      dmsg << "flag: " << flag;
      dmsg << "fName[0]: " << fName[0];
      dmsg << "fName[1]: " << fName[1];
      #endif

      if (flag) {
        cctxt(com_mod, cm_mod, eq, fName, wtn);
      } else {
        wtxt(com_mod, cm_mod, eq, l, fName, tmpV, wtn, div, pflag);
      }
    }

    // IB outputs

    // TODO:DaveP] not implemented.
    //
/* 
    if (.NOT.ibFlag) continue;
    if (ALLOCATED(tmpV)) DEALLOCATE(tmpV)
    ALLOCATE (tmpV(maxnsd,ib.tnNo))
    for (int iOut=1, eq(iEq).nOutIB
      // Don't write it when it doesn't suppose to be written
      wtn = eq(iEq).outIB(iOut).wtn(2:3)
      if (ALL(.NOT.wtn)) continue;
      l = eq(iEq).outIB(iOut).l
      s = eq(iEq).s + eq(iEq).outIB(iOut).o
      e = s + l - 1

      oGrp = eq(iEq).outIB(iOut).grp
      div  = .TRUE.

      SELECT case (oGrp)
      case (outGrp_NA)
        err = "NA outGrp in TXT"
        case (outGrp_Y)
           for (int a=1, ib.tnNo
              tmpV(1:l,a) = ib.Yb(s:e,a)/REAL(cm.np(), KIND=RKIND)
           {
        case (outGrp_D)
           for (int a=1, ib.tnNo
              tmpV(1:l,a) = ib.Ubo(s:e,a)/REAL(cm.np(), KIND=RKIND)
           {
        case (outGrp_I)
           div = .FALSE.
           for (int a=1, ib.tnNo
              tmpV(1:l,a) = 1.
           {
        case DEFAULT
           err = "Undefined output"
        END SELECT

        fName = eq(iEq).sym//"_"//TRIM(eq(iEq).outIB(iOut).name)
        if (l .EQ. nsd) {
           fName(1) = TRIM(appPath)//"IB_B_"//TRIM(fName(1))//
     2        "_flux.txt"
        } else {
           fName(1) = TRIM(appPath)//"IB_B_"//TRIM(fName(1))//
     2        "_average.txt"
        }
        fName(2) = TRIM(appPath)//"IB_V_"//TRIM(fName(2))//
     2     "_average.txt"

        if (flag) {
           CALL IB_CCTXT(eq(iEq), fName, wtn)
        } else {
           CALL IB_WTXT(eq(iEq), l, fName, tmpV, wtn, div)
        }
     }
*/
  }

  #ifdef debug_txt
  dmsg << "Done txt";
  #endif
}

//-----
// wtx
//-----
// This is to write to txt file
//
// TODO:DaveP] not fully implemented.
//
// WARNING: There is an indexing problem here because 'm' is a length and not an index.
//
void wtxt(const ComMod& com_mod, CmMod& cm_mod, const eqType& lEq, const int m, const std::array<std::string,2>& fName, 
    const Array<double>& tmpV, const std::vector<bool>& wtn, const bool div, const bool pFlag)
{
  return;

  #ifdef debug_wtxt
  DebugMsg dmsg(__func__, com_mod.cm.idcm());
  dmsg.banner();
  dmsg << "m: " << m;
  dmsg << "div: " << div;
  dmsg << "pFlag: " << pFlag;
  #endif

  int fid = 1;
  int nsd = com_mod.nsd;
  FILE* fp = nullptr;

  for (int i = 0; i < 2; i++) {
    if (!wtn[i]) {
      continue;
    }

    if (com_mod.cm.mas(cm_mod)) {
      fp = fopen(fName[i].c_str(), "a");
    }

    if (i == 0) {
      for (int iM = 0; iM < com_mod.nMsh; iM++) {
        auto& msh = com_mod.msh[iM];
        bool lTH = false; 

        if (msh.nFs == 2) {
          lTH = true;
        }

        for (int iFa = 0; iFa < msh.nFa; iFa++) {
          auto& fa = msh.fa[iFa];
          double tmp = 0.0;

          if (m == 1) {
            if (div) {
              tmp = fa.area;
              tmp = all_fun::integ(com_mod, cm_mod, fa, tmpV, 1) / tmp;
            } else {
              if (pFlag && lTH) {
                tmp = all_fun::integ(com_mod, cm_mod, fa, tmpV, 1, true);
              } else {
                tmp = all_fun::integ(com_mod, cm_mod, fa, tmpV, 1);
              }
            }

          } else if (m == nsd) {
            tmp = all_fun::integ(com_mod, cm_mod, fa, tmpV, 1, m);
          } else {
            throw std::runtime_error("WTXT only accepts 1 and nsd");
          }

          if (com_mod.cm.mas(cm_mod)) {
            fprintf(fp, " %.10e ", tmp);
          }
        }
      }

    } else {
      for (int iDmn = 0; iDmn < lEq.nDmn; iDmn++) {
        auto& dmn = lEq.dmn[iDmn];
        double tmp = 0.0;

        if (div) {
          tmp = dmn.v;
          tmp = all_fun::integ(com_mod, cm_mod, dmn.Id, tmpV, 1, m) / tmp;
        } else {
          tmp = all_fun::integ(com_mod, cm_mod, dmn.Id, tmpV, 1, m, pFlag);
        }
        if (com_mod.cm.mas(cm_mod)) {
          fprintf(fp, " %.10e ", tmp);
        }
      }
    }

    if (com_mod.cm.mas(cm_mod)) {
      fprintf(fp, "\n");
      fclose(fp);
    }
  }
}

};
