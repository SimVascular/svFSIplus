/* unit test for get_pk2cc
*
*
*/

// include GoogleTest
#include "gtest/gtest.h"

// include header files from svFSIplus
#include "mat_fun.h"
#include "mat_fun_carray.h"
#include "mat_models.h"
#include "mat_models_carray.h"
#include "Simulation.h"
#include "sv_struct.h"

// from main:
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
#include "remesh.h"
#include "set_bc.h"
#include "txt.h"
#include "ustruct.h"
#include "vtk_xml.h"


// include libraries
#include <stdlib.h>
#include <iomanip>
#include <iostream>

// Using global variables
Simulation *simulation;
ComMod com_mod;
CmMod cm_mod_;
cmType cm;
CepMod cep_mod;

// S, Dm: target variables 
Array<double> S(3,3), Dm(6,6);
// I: identity matrix
Array<double> I = mat_fun::mat_id(3);
// Z: zero matrix
Array<double> Z(3,3);



void read_files(Simulation* simulation, const std::string& file_name)
{
  simulation->com_mod.timer.set_time();

  if (simulation->com_mod.cm.slv(simulation->cm_mod)) {
    return;
  }

  read_files_ns::read_files(simulation, file_name);

/*
  try {
    read_files_ns::read_files(simulation, file_name);

  } catch (const std::exception& exception) {
    std::cout << "[svFSIplus] ERROR The svFSIplus program has failed." << std::endl;
    std::cout << "[svFSIplus] ERROR " << exception.what() << std::endl;
    exit(1);
  }
*/
  
}

TEST(Dummy, OneEqualsToOne) {
    EXPECT_EQ(1, 1);
}

TEST(GetPK2ccTest, Identity) {
    for (short i = 0; i < 3; i++){
        for (short j = 0; j < 3; j++){
            EXPECT_EQ(S(i,j), Z(i,j));   
        }
    }
}

// Main function running the tests
int main(int argc, char **argv) {
    // std::cout << "Current working directory: " 
    //           << std::filesystem::current_path() << std::endl;

    // Initialize MPI.
    //
    int mpi_rank, mpi_size;
    // MPI_Init(&argc, &argv);
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
    // Original: initiate in main():

    simulation = new Simulation();

    // filename for the unit element .xml: (To be complete !!)
    // current_path: /path_to_svFSIplus/svFSIplus/build/svFSI-build/Source/svFSI
    std::string file_name = std::filesystem::current_path().string() + "/../../../../tests/unitTests/unitTest.xml";
    // std::string file_name = "/Users/yuechengyu/Work/Cardiac/svFSIplus_unitTest/svFSIplus/tests/unitTests/unitTest.xml"; 
    // std::string file_name(argv[1]);


    read_files(simulation, file_name);

    std::cout << "Finish reading files" << std::endl;

    distribute(simulation);
    Vector<double> init_time(3);
    initialize(simulation, init_time);

    // Original: define from struct_3d/struct_3d_carray(...) <- construct_dsolid(...)

    com_mod = simulation->com_mod;
    cm_mod_ = simulation->cm_mod;
    cm = com_mod.cm;
    cep_mod = simulation->get_cep_mod();

    const int nsd  = com_mod.nsd;

    std::cout << "The number of meshes: " << com_mod.nMsh << std::endl;

    auto iM = com_mod.nMsh;   // nMsh = 1 from unit mesh 
    auto& lM = com_mod.msh[iM-1]; 
    auto e = lM.gnEl;   // global number of elements

    std::cout << "The number of elements: " << e << std::endl;

    int eNoN = lM.eNoN;
    int nFn = lM.nFn;
    if (nFn == 0) {
      nFn = 1;
    }
    std::cout << "The variable eNoN is : " << eNoN << std::endl;
    std::cout << "The variable nFn is : " << nFn << std::endl;

    Array<double> fN(nsd,nFn);
    fN  = 0.0;
    if (lM.fN.size() != 0) {
      for (int iFn = 0; iFn < nFn; iFn++) {
        for (int i = 0; i < nsd; i++) {
          fN(i,iFn) = lM.fN(i+nsd*iFn,e);
        }
      }
    }

    int cEq = com_mod.cEq;
    auto& eq = com_mod.eq[cEq];
    int cDmn = com_mod.cDmn;
    auto& dmn = eq.dmn[cDmn];



    // ya_g or ya: a constant related to electromechanics ??
    double ya_g = 0.0;
    


    std::cout << "Before calling get_pk2cc: input argument F" << std::endl;
    std::cout << "=================== INPUT ===================" << std::endl;
    // F: deformation gradient
    Array<double> F(3,3);
    F = I;
    // F(0, 0) = 2.0;
    // const auto& stM = dmn.stM;
    // std::cout << "Material type: " << stM.isoType << std::endl;
    F.print("Input Deformation ");

    // Original: call in sv_struct.cpp
    mat_models::get_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);

    std::cout <<"=================== OUTPUT ===================" << std::endl;
    std::cout << "After calling get_pk2cc: output S and Dm" << std::endl;
    S.print("Output stress tensor:");

    std::cout << "After calling get_pk2cc: output S and Dm" << std::endl;
    Dm.print("Output stiffness matrix:");

    // mat_models_carray::get_pk2cc<3>(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
    // from struct_3d(...): com_mod, cep_mod, dmn <- eq <- com_mod, nFn, fN 
    // from construct_dsolid(...): com_mod, cep_mod, nFn <- lM, fN <- lM
    // from global_eq_assem(...): com_mod, cep_mod, lM
    // from iterate_simulation(simulation): com_mod <- simulation, cep_mod <- simulation, lM <- com_mod.msh[iM] (iM < com_mod.nMsh)
    // iterate_simulation(simulation) <- run_simulation(simulation) <- main(...)

    std::cout << "This is svFSIplus/tests/unitTests/test.cpp now!" << std::endl;

    std::cout << "Let's try EXPECT_EQ(S, actual_value)" << std::endl;


    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

}
