
#include "Simulation.h"

#include "all_fun.h"
#include "bf.h"
#include "contact.h"
#include "distribute.h"
#include "eq_assem.h"
#include "fs.h"
#include "initialize.h"
#include "ls.h"
#include "output.h"
#include "pic.h"
#include "read_files.h"
#include "read_msh.h"
#include "set_bc.h"
#include "txt.h"
#include "ustruct.h"
#include "vtk_xml.h"

#include <stdlib.h>
#include <iomanip>
#include <iostream>

//------------
// read_files
//------------
// Read in XML file and all mesh and BC data.  
//
void read_files(Simulation* simulation, const std::string& file_name)
{
  simulation->com_mod.timer.set_time();

  if (simulation->com_mod.cm.slv(simulation->cm_mod)) {
    return;
  }

  try {
    read_files_ns::read_files(simulation, file_name);

  } catch (const std::exception& exception) {
    std::cout << "[svFSI] ERROR The svFSI program has failed." << std::endl;
    std::cout << "[svFSI] ERROR " << exception.what() << std::endl;
    exit(1);
  }
  
}

//------
// main
//------
//
int main(int argc, char *argv[])
{
  if (argc != 2) {
    std::cout << "[svFSI:ERROR] The svFSI program requires the solver input file name as an argument." << std::endl;
    exit(1);
  }

  std::string file_name(argv[1]);

  MPI_Init(&argc, &argv);

  auto simulation = new Simulation();

  std::cout << std::scientific << std::setprecision(16);

  read_files(simulation, file_name);
}

