
#ifndef LOAD_MSH_H 
#define LOAD_MSH_H 

#include "ComMod.h"
#include "Parameters.h"
#include "Simulation.h"

#include <string>

namespace load_msh {

  void read_ccne(Simulation* simulation, mshType& mesh, const MeshParameters* mesh_param);

  void read_ndnlff(const std::string& file_name, faceType& face);

  void read_sv(Simulation* simulation, mshType& mesh, const MeshParameters* param);

};

#endif

