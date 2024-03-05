/* Copyright (c) Stanford University, The Regents of the University of California, and others.
 *
 * All Rights Reserved.
 *
 * See Copyright-SimVascular.txt for additional details.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

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


int argc; 
char **argv;

class MockSimulation : public Simulation {
public:
    MockSimulation() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockComMod : public ComMod {
public:
    MockComMod() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockCepMod : public CepMod {
public:
    MockCepMod() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockcmType : public cmType {
public:
    MockcmType() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockstModelType : public stModelType {
public:
    MockstModelType() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockdmnType : public dmnType {
public:
    MockdmnType() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockcemModelType : public cemModelType {
public:
    MockcemModelType() {
        // initialize if needed 
    }
    // Mock methods if needed
};

TEST(TestMockObject, sampleTest) {
    MockSimulation *simulation;

    MockComMod com_mod;
    MockCepMod cep_mod;
    MockcmType cm;
    MockstModelType stM;
    MockdmnType dmn;
    MockcemModelType cem;

    // preparing input arguments for get_pk2cc
    // auto com_mod = simulation->com_mod;

    std::cout << "test brreak point 1" << std::endl;

    // auto cm_mod_ = simulation->cm_mod;
    // auto cm = com_mod.cm;
    // auto cep_mod = simulation->get_cep_mod();
    com_mod.nsd  = 3;   // setup by hand
    int nsd = com_mod.nsd;
    dmn.stM = stM;
    int iM = 1;   // nMsh = 1 from unit mesh 
    // auto& lM = com_mod.msh[iM-1]; 
    // auto e = lM.gnEl;   // global number of elements
    // int eNoN = lM.eNoN;
    // int nFn = lM.nFn;
    int nFn = 0; 
    if (nFn == 0) {nFn = 1;}

    // for potential anisotropic material:
    // Array<double> fN(nsd,nFn);
    // fN  = 0.0;
    // if (lM.fN.size() != 0) {
    //   for (int iFn = 0; iFn < nFn; iFn++) {
    //     for (int i = 0; i < nsd; i++) {
    //       fN(i,iFn) = lM.fN(i+nsd*iFn,e);
    //     }
    //   }
    // }

    std::cout << "test brreak point 2" << std::endl;

    // int cEq = com_mod.cEq;
    // auto& eq = com_mod.eq[cEq];
    // int cDmn = com_mod.cDmn;
    // auto dmn = eq.dmn[cDmn];
    double ya_g = 0.0;   // ya_g or ya: a constant related to electromechanics ?

    std::cout << "test brreak point 3" << std::endl;


    // Step 1: define material properties
    // set stM.isoType as NeoHookean
    // auto& stM = dmn.stM;
    dmn.stM.isoType = consts::ConstitutiveModelType::stIso_nHook;
    // set material parameters
    double E = 1e6;   // Elasticity_modulus
    double nu = 0.5;   // Poisson_ratio
    double mu  = 0.5 * E / (1.0 + nu);   // Shear_modulus
    dmn.stM.C10 = 0.5 * mu;   // set_material_props.h 
    dmn.stM.volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
    dmn.stM.Kpen = 4e9;   // Penalty_parameter

    std::cout << "test brreak point 4" << std::endl;

    // Step 2: define fiber orientation (anisotropic only)
    Array<double> fN(nsd,nFn);
    fN  = 0.0;

    // Step 3: define the input 
    Array<double> F = mat_fun::mat_id(3);
    // Step 4: define the reference output
    Array<double> S_ref(3,3), Dm_ref(6,6);
    S_ref = 0.0; Dm_ref = 0.0;
    // Step 5: compare
    Array<double> S(3,3), Dm(6,6);

    std::cout << "test brreak point 5" << std::endl;

    cep_mod.cem = cem;
    cep_mod.cem.aStress = false;
    cep_mod.cem.aStrain = false;

    mat_models::get_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);

    std::cout << "test brreak point 6" << std::endl;

    for (int i = 0; i < 3; i++){
      for (int j = 0; j < 3; j++){
        EXPECT_EQ(S(i,j), S_ref(i,j));   
      }
    }

    // more tests ...
    F = mat_fun::mat_id(3);
    F(0,0) = 2.0;
    mat_models::get_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
    for (int i = 0; i < 3; i++){
        for (int j = 0; j < 3; j++){
          if (i == j) {
            EXPECT_NE(S(i,j), S_ref(i,j));   
          } else {
            EXPECT_EQ(S(i,j), S_ref(i,j)); 
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


// class UnitTest_getpk2cc : public ::testing::Test {
// protected:
//   // parameters for simulation setup 
//   Simulation *simulation;
//   ComMod com_mod;
//   CepMod cep_mod; 
//   dmnType dmn;
//   int nFn;
//   double ya_g;
//   int nsd;

//   static void SetUpTestSuite() {
//     // Initialize MPI.
//     //
//     int mpi_rank, mpi_size;
//     MPI_Init(&argc, &argv);
//     MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
//   }

//   static void TearDownTestSuite() {
//     // Finalize MPI once after all tests are done
//     MPI_Finalize();
//   }

//   void setUpSimulation (std::string file_name) {
//     // this function can be used if the user has designed input file
//   }

//   void SetUp() override {
    
//     std::string file_name = "../../../../tests/unitTests/unitTest.xml";
//     simulation = new Simulation();
    
//     read_files(simulation, file_name);

//     distribute(simulation);
//     Vector<double> init_time(3);
//     initialize(simulation, init_time);

//     // preparing input arguments for get_pk2cc
//     com_mod = simulation->com_mod;
//     auto cm_mod_ = simulation->cm_mod;
//     auto cm = com_mod.cm;
//     cep_mod = simulation->get_cep_mod();
//     nsd  = com_mod.nsd;
//     auto iM = com_mod.nMsh;   // nMsh = 1 from unit mesh 
//     auto& lM = com_mod.msh[iM-1]; 
//     auto e = lM.gnEl;   // global number of elements
//     int eNoN = lM.eNoN;
//     nFn = lM.nFn;
//     if (nFn == 0) {nFn = 1;}

//     // for potential anisotropic material:
//     // Array<double> fN(nsd,nFn);
//     // fN  = 0.0;
//     // if (lM.fN.size() != 0) {
//     //   for (int iFn = 0; iFn < nFn; iFn++) {
//     //     for (int i = 0; i < nsd; i++) {
//     //       fN(i,iFn) = lM.fN(i+nsd*iFn,e);
//     //     }
//     //   }
//     // }

//     int cEq = com_mod.cEq;
//     auto& eq = com_mod.eq[cEq];
//     int cDmn = com_mod.cDmn;
//     dmn = eq.dmn[cDmn];
//     ya_g = 0.0;   // ya_g or ya: a constant related to electromechanics ?
//   }

// };


// TEST_F(UnitTest_getpk2cc, NeoHookean) {
//   // Step 1: define material properties
//   // set stM.isoType as NeoHookean
//   auto& stM = dmn.stM;
//   stM.isoType = consts::ConstitutiveModelType::stIso_nHook;
//   // set material parameters
//   double E = 1e6;   // Elasticity_modulus
//   double nu = 0.5;   // Poisson_ratio
//   double mu  = 0.5 * E / (1.0 + nu);   // Shear_modulus
//   stM.C10 = 0.5 * mu;   // set_material_props.h 
//   stM.volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
//   stM.Kpen = 4e9;   // Penalty_parameter

//   // Step 2: define fiber orientation (anisotropic only)
//   Array<double> fN(nsd,nFn);
//   fN  = 0.0;

//   // Step 3: define the input 
//   Array<double> F = mat_fun::mat_id(3);
//   // Step 4: define the reference output
//   Array<double> S_ref(3,3), Dm_ref(6,6);
//   S_ref = 0.0; Dm_ref = 0.0;
//   // Step 5: compare
//   Array<double> S(3,3), Dm(6,6);
//   mat_models::get_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
//   for (int i = 0; i < 3; i++){
//     for (int j = 0; j < 3; j++){
//       EXPECT_EQ(S(i,j), S_ref(i,j));   
//     }
//   }

//   // more tests ...
//   F = mat_fun::mat_id(3);
//   F(0,0) = 2.0;
//   mat_models::get_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
//   for (int i = 0; i < 3; i++){
//       for (int j = 0; j < 3; j++){
//         if (i == j) {
//           EXPECT_NE(S(i,j), S_ref(i,j));   
//         } else {
//           EXPECT_EQ(S(i,j), S_ref(i,j)); 
//         }   
//       }
//   }

// }


// TEST_F(UnitTest_getpk2cc, NeoHookean_2) {
//   // Step 1: define material properties
//   // set stM.isoType as NeoHookean
//   auto& stM = dmn.stM;
//   stM.isoType = consts::ConstitutiveModelType::stIso_nHook;
//   // set material parameters
//   double E = 2e6;   // Elasticity_modulus
//   double nu = 0.5;   // Poisson_ratio
//   double mu  = 0.5 * E / (1.0 + nu);   // Shear_modulus
//   stM.C10 = 0.5 * mu;   // set_material_props.h 
//   stM.volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
//   stM.Kpen = 4e9;   // Penalty_parameter

//   // Step 2: define fiber orientation (anisotropic only)
//   Array<double> fN(nsd,nFn);
//   fN  = 0.0;

//   // Step 3: define the input 
//   Array<double> F = mat_fun::mat_id(3);
//   // Step 4: define the reference output
//   Array<double> S_ref(3,3), Dm_ref(6,6);
//   S_ref = 0.0; Dm_ref = 0.0;
//   // Step 5: compare
//   Array<double> S(3,3), Dm(6,6);
//   mat_models::get_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
//   for (int i = 0; i < 3; i++){
//     for (int j = 0; j < 3; j++){
//       EXPECT_EQ(S(i,j), S_ref(i,j));   
//     }
//   }

//   // more tests ...
//   F = mat_fun::mat_id(3);
//   F(0,0) = 2.0;
//   mat_models::get_pk2cc(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
//   for (int i = 0; i < 3; i++){
//       for (int j = 0; j < 3; j++){
//         if (i == j) {
//           EXPECT_NE(S(i,j), S_ref(i,j));   
//         } else {
//           EXPECT_EQ(S(i,j), S_ref(i,j)); 
//         }   
//       }
//   }
// }


int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}

/* ++++++++++++++++++++++++ why do we need an input file for unit test? ++++++++++++++++++++++++++
  the following is where each argument of get_pk2cc from: 
  mat_models_carray::get_pk2cc<3>(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
  from struct_3d(...): com_mod, cep_mod, dmn <- eq <- com_mod, nFn, fN 
  from construct_dsolid(...): com_mod, cep_mod, nFn <- lM, fN <- lM
  from global_eq_assem(...): com_mod, cep_mod, lM
  from iterate_simulation(simulation): com_mod <- simulation, cep_mod <- simulation, 
  lM <- com_mod.msh[iM] (iM < com_mod.nMsh)
  iterate_simulation(simulation) <- run_simulation(simulation) <- main(...) 
*/ 