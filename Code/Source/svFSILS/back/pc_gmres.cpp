
#include "pc_gmres.h"

#include "fsils_api.hpp"

#include "add_bc_mul.h"
#include "bcast.h"
#include "dot.h"
#include "norm.h"
#include "omp_la.h"
#include "spar_mul.h"

#include "Array3.h"

#include <math.h>

namespace pc_gmres {

//----------
// pc_gmres
//----------
//
void pc_gmres(fsi_linear_solver::FSILS_lhsType& lhs, fsi_linear_solver::FSILS_subLsType& ls, const int dof, 
    const Array<double>& Val, const Array<double>& R)
{
  #define n_debug_pc_gmres
  #ifdef debug_pc_gmres
  int tid = lhs.commu.task;
  auto msg_prefix = std::string("[pc_gmres:") + std::to_string(tid) + "] ";
  std::cout << msg_prefix << std::endl;
  std::cout << msg_prefix << "========== pc_gmres ==========" << std::endl;
  #endif

  using namespace fsi_linear_solver;

/*
  int nNo = lhs.nNo;
  int mynNo = lhs.mynNo;
  int sD = ls.GM.sD;
  #ifdef debug_gmres
  std::cout << msg_prefix << "dof: " << dof << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  std::cout << msg_prefix << "mynNo: " << mynNo << std::endl;
  std::cout << msg_prefix << "sD: " << sD << std::endl;
  #endif
*/

}

};
