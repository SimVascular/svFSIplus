
// The functions here reproduce the subroutines defined in svFSILS/COMMU.f.

#include "fsils.hpp"
#include "CmMod.h"

#include "mpi.h"

namespace fsi_linear_solver {

//--------------------
// fsils_commu_create
//--------------------
//
// Modifies:
//   commu.task
//   commu.nTasks
//   
void fsils_commu_create(FSILS_commuType& commu, cm_mod::MpiCommWorldType commi)
{
  // Some of these parameters are set for sequential version
  //
  commu.foC = true;
  commu.comm = commi;
  commu.nTasks = 1;
  commu.task   = 0;
  commu.master = 0;

  MPI_Comm_rank(commi, &commu.task);
  MPI_Comm_size(commi, &commu.nTasks);

  //auto msg_prefix = std::string("[fsils_commu_create:") + std::to_string(commu.task) + "] ";
  //std::cout << msg_prefix << std::endl;
  //std::cout << msg_prefix << "========== fsils_commu_create ==========" << std::endl;
  //std::cout << msg_prefix << "commu.task: " << commu.task << std::endl;
  //std::cout << msg_prefix << "commu.nTasks: " << commu.nTasks << std::endl;

  MPI_Allreduce(&commu.task, &commu.master, 1, cm_mod::mpint, MPI_MIN, commi);

  commu.masF = false;
  commu.tF = commu.task;
  //commu.tF = commu.task + 1;    // For Fortran task number
  if (commu.task == commu.master) {
    commu.masF = true; 
  }

  //std::cout << msg_prefix << "commu.masF: " << commu.masF << std::endl;
  //std::cout << msg_prefix << "commu.tF: " << commu.tF << std::endl;
  //std::cout << msg_prefix << "Done." << std::endl;
}


};


