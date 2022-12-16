
// For calculating the norm of a scaler or vector based vector.
//
//
// Only the part of U which is owned by this processor is included in
// norm calculation, i.e. U(:,1:cS(tF)%ptr+cS(tF)%n-1)
// In order to have the correct answer it is needed that COMMU has
// been done before calling this function (or the ansesters of U
// are passed through COMMU)

#include "norm.h"

#include "CmMod.h"

#include "mpi.h"

#include <math.h>

namespace norm {

//--------------
// fsi_ls_norms
//--------------
//
double fsi_ls_norms(const int nNo, FSILS_commuType& commu, const Vector<double>& U)
{
  #define n_debug_fsi_ls_norms 
  #ifdef debug_fsi_ls_norms
  int tid = commu.task;
  auto msg_prefix = std::string("[fsi_ls_norms:") + std::to_string(tid) + "] ";
  std::cout << msg_prefix << "========== fsi_ls_norms ==========" << std::endl;
  std::cout << msg_prefix << "nNo: " << nNo << std::endl;
  #endif

  double result = 0.0;

  for (int i = 0; i < nNo; i++) {
    result = result + U(i)*U(i);
  } 

  #ifdef debug_fsi_ls_norms
  std::cout << msg_prefix << "result: " << result << std::endl; 
  #endif

  if (commu.nTasks != 1) {
    double tmp;
    MPI_Allreduce(&result, &tmp, 1, cm_mod::mpreal, MPI_SUM, commu.comm);
    //CALL MPI_ALLREDUCE(FSILS_NORMS, tmp, 1, mpreal, MPI_SUM, commu%comm, ierr)
    result = tmp;
    #ifdef debug_fsi_ls_norms
    std::cout << msg_prefix << "tmp: " << tmp << std::endl; 
    std::cout << msg_prefix << "sqrt(tmp): " << sqrt(tmp) << std::endl; 
    #endif
  }

  return sqrt(result);
}

//--------------
// fsi_ls_normv 
//--------------
//
double fsi_ls_normv(const int dof, const int nNo, FSILS_commuType& commu, const Array<double>& U)
{
  double result = 0.0;

  switch (dof) {
    case 1: {
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*U(0,i);
      }
    } break; 

    case 2: {
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*U(0,i) + U(1,i)*U(1,i);
      }
    } break; 

    case 3: {
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*U(0,i) + U(1,i)*U(1,i) + U(2,i)*U(2,i);
      }
    } break; 

    case 4: {
      for (int i = 0; i < nNo; i++) {
        result = result + U(0,i)*U(0,i) + U(1,i)*U(1,i) + U(2,i)*U(2,i) + U(3,i)*U(3,i);
      }
    } break; 

    default: { 
      for (int i = 0; i < nNo; i++) {
        for (int j = 0; j < U.num_rows(); j++) {
          result = result + U(j,i)*U(j,i);
        }
        //result = result + U(:,i)*U(:,i);
      }
    } break; 
  } 

  if (commu.nTasks != 1) {
    double tmp;
    MPI_Allreduce(&result, &tmp, 1, cm_mod::mpreal, MPI_SUM, commu.comm);
    //CALL MPI_ALLREDUCE(result, tmp, 1, mpreal, MPI_SUM, commu.comm, ierr)
    result = tmp;
  }

  return sqrt(result);
}

};
