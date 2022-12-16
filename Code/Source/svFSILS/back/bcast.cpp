
// To broadcast a variable to all processors
//
// Reproduces code in BCAST.f.

#include "bcast.h"

#include "mpi.h"

namespace bcast {

//-------------
// fsils_bcast
//-------------
//
void fsils_bcast(double& u, FSILS_commuType& commu)
{
  if (commu.nTasks > 1) { 
    double uG;
    MPI_Allreduce(&u, &uG, 1, cm_mod::mpreal, MPI_SUM, commu.comm);
    //CALL MPI_ALLREDUCE(u, uG, 1, mpreal, MPI_SUM, commu%comm, ierr)
    u = uG;
  } 
}

//---------------
// fsils_bcast_v
//---------------
//
void fsils_bcast_v(const int n, Vector<double>& u, FSILS_commuType& commu)
{
  //std::cout << "[fsils_bcast_v] n: " << n << std::endl;
  if (commu.nTasks > 1) { 
    Vector<double> uG(n);
    MPI_Allreduce(u.data(), uG.data(), n, cm_mod::mpreal, MPI_SUM, commu.comm);
    u = uG;
  } 
}

};


