/* unit test for get_pk2cc
*  this program generates a simple svFSI simulation to test function get_pk2cc.
*  it allows you to create input deformation gradients to check the output 
*  stress and tangent modulus. it usese Google Test for test cases.
*/

#include <stdlib.h>
#include <iomanip>
#include <iostream>
#include "gtest/gtest.h"   // include GoogleTest
#include "mat_fun.h"
#include "mat_fun_carray.h"
#include "mat_models.h"
#include "mat_models_carray.h"
#include "Simulation.h"
#include "sv_struct.h"
#include "distribute.h"
#include "initialize.h"
#include "read_files.h"


// Using global variables: 
// for the output: stress tensor and tangent tensor

Array<double> S_0(3,3), Dm_0(6,6);      // for identity test
Array<double> S_1(3,3), Dm_1(6,6);      // for non-identity test


// Google Test test cases

TEST(ExampleTest, OneEqualsToOne) {
  EXPECT_EQ(1, 1);
}

TEST(GetPK2ccTest, Identity) {
  // input: 
  // 1 0 0
  // 0 1 0
  // 0 0 1
  for (short i = 0; i < 3; i++){
      for (short j = 0; j < 3; j++){
        EXPECT_EQ(S_0(i,j), 0);   
      }
  }
}

TEST(GetPK2ccTest, Non_Identity) {
  // input: 
  // 2 0 0
  // 0 1 0
  // 0 0 1
  for (short i = 0; i < 3; i++){
      for (short j = 0; j < 3; j++){
        if (i == j) {
          EXPECT_NE(S_1(i,j), 0);   
        } else {
          EXPECT_EQ(S_1(i,j), 0); 
        }   
      }
  }
}


void read_files(Simulation* simulation, const std::string& file_name)
{
  simulation->com_mod.timer.set_time();

  if (simulation->com_mod.cm.slv(simulation->cm_mod)) {
    return;
  }

  read_files_ns::read_files(simulation, file_name);
}


// Main function running the tests
int main(int argc, char **argv) {
    // define deformation gradient as input for get_pk2cc
    Array<double> F_0 = mat_fun::mat_id(3);
    Array<double> F_1 = mat_fun::mat_id(3);
    F_1(0,0) = 2;

    // Initialize MPI.
    //
    int mpi_rank, mpi_size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    auto simulation = new Simulation();

    std::string file_path = "/../../../../tests/unitTests/unitTest.xml";
    std::string file_name = std::filesystem::current_path().string() + file_path; 

    read_files(simulation, file_name);

    std::cout << "Finish reading files" << std::endl;

    distribute(simulation);
    Vector<double> init_time(3);
    initialize(simulation, init_time);

    // preparing input arguments for get_pk2cc
    auto com_mod = simulation->com_mod;
    auto cm_mod_ = simulation->cm_mod;
    auto cm = com_mod.cm;
    auto cep_mod = simulation->get_cep_mod();
    const int nsd  = com_mod.nsd;
    auto iM = com_mod.nMsh;   // nMsh = 1 from unit mesh 
    auto& lM = com_mod.msh[iM-1]; 
    auto e = lM.gnEl;   // global number of elements
    int eNoN = lM.eNoN;
    int nFn = lM.nFn;
    if (nFn == 0) {nFn = 1;}
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
    double ya_g = 0.0;   // ya_g or ya: a constant related to electromechanics ?
    
    // Call get_pk2cc, which is originally called from sv_struct.cpp
    mat_models::get_pk2cc(com_mod, cep_mod, dmn, F_0, nFn, fN, ya_g, S_0, Dm_0);
    mat_models::get_pk2cc(com_mod, cep_mod, dmn, F_1, nFn, fN, ya_g, S_1, Dm_1);

    // results printing
    // std::cout <<"==================== INPUT ===================" << std::endl;
    // F.print("Input Deformation ");
    // std::cout <<"=================== OUTPUT ===================" << std::endl;
    // S.print("Output stress tensor:");
    // Dm.print("Output stiffness matrix:");

    // the following is where each argument of get_pk2cc from: 
    // mat_models_carray::get_pk2cc<3>(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
    // from struct_3d(...): com_mod, cep_mod, dmn <- eq <- com_mod, nFn, fN 
    // from construct_dsolid(...): com_mod, cep_mod, nFn <- lM, fN <- lM
    // from global_eq_assem(...): com_mod, cep_mod, lM
    // from iterate_simulation(simulation): com_mod <- simulation, cep_mod <- simulation, 
    // lM <- com_mod.msh[iM] (iM < com_mod.nMsh)
    // iterate_simulation(simulation) <- run_simulation(simulation) <- main(...)

    // run Google Test
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();

}
