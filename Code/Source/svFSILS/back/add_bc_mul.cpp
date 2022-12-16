
//--------------------------------------------------------------------
// The contribution of coupled BCs is added to the matrix-vector
// product operation. Depending on the type of operation (adding the
// contribution or compution the PC contribution) different
// coefficients are used.
//--------------------------------------------------------------------

// Reproduces code in ADDBCMUL.f.

#include "add_bc_mul.h"

#include "dot.h"

namespace add_bc_mul {

//------------
// add_bc_mul
//------------
//
void add_bc_mul(FSILS_lhsType& lhs, const BcopType op_Type, const int dof, const Array<double>& X, Array<double>& Y)
{
  //int tid = lhs.commu.task;
  //auto msg_prefix = std::string("[add_bc_mul:") + std::to_string(tid) + "] ";
  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "========== add_bc_mul ==========" << std::endl;

  Vector<double> coef(lhs.nFaces); 
  Array<double> v(dof,lhs.nNo);

  if (op_Type == BcopType::BCOP_TYPE_ADD) {
    //std::cout << msg_prefix << "BCOP_TYPE_ADD" << std::endl;
    for (int i = 0; i < lhs.nFaces; i++) {
      coef(i) = lhs.face[i].res;
      //std::cout << msg_prefix << "coef(" << i+1 << "): " << coef(i) << std::endl;
    }
  } else if (op_Type == BcopType::BCOP_TYPE_PRE) {
    //std::cout << msg_prefix << "BCOP_TYPE_PRE" << std::endl;
    for (int i = 0; i < lhs.nFaces; i++) {
      coef(i) = -lhs.face[i].res / (1.0 + (lhs.face[i].res*lhs.face[i].nS));
      //std::cout << msg_prefix << "coef(" << i+1 << "): " << coef(i) << std::endl;
    }
  } else { 
    //PRINT *, "FSILS: op_Type is not defined"
    //STOP "FSILS: FATAL ERROR"
  }

  for (int faIn = 0; faIn < lhs.nFaces; faIn++) {
    //std::cout << msg_prefix << "----- faIn " << faIn+1 << " -----" << std::endl;
    auto& face = lhs.face[faIn];
    int nsd = std::min(face.dof, dof);
    //std::cout << msg_prefix << "nsd: " << nsd << std::endl;

    if (face.coupledFlag) {
      //std::cout << msg_prefix << "coupledFlag faIn: " << faIn+1 << std::endl;
      if (face.sharedFlag) {
        //std::cout << msg_prefix << "  sharedFlag " << std::endl;
        v = 0.0;

        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            v(i,Ac) = face.valM(i,a);
          }
        }

        double S = coef(faIn) * dot::fsils_dot_v(dof, lhs.mynNo, lhs.commu, v, X);

        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            Y(i,Ac) = Y(i,Ac) + v(i,Ac)*S;
          }
        }

      } else  {
        //std::cout << msg_prefix << "  not sharedFlag " << std::endl;
        double S = 0.0;
        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          for (int i = 0; i < nsd; i++) {
            S = S + face.valM(i,a)*X(i,Ac);
            //std::cout << msg_prefix << "  valM(i,a): " << a+1 << " " << i+1 << " " << face.valM(i,a) << std::endl;
          }
        }

        S = coef(faIn) * S;
        //std::cout << msg_prefix << "  S: " << S << std::endl;

        for (int a = 0; a < face.nNo; a++) {
          int Ac = face.glob(a);
          //std::cout << msg_prefix << "  Ac: " << Ac+1 << std::endl;
          for (int i = 0; i < nsd; i++) {
            Y(i,Ac) = Y(i,Ac) + face.valM(i,a)*S;
            //std::cout << msg_prefix << "  Y(i,Ac): " << Y(i,Ac) << std::endl;
          }
        }
      }
    }
  }

}

};
