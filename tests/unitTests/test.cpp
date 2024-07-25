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

#include "test.h"

using namespace mat_fun;
using namespace std;

// ----------------------------------------------------------------------------
// --------------------------- Neo-Hookean Material ---------------------------
// ----------------------------------------------------------------------------

// Test fixture class for Neo-Hookean material model
class NeoHookeanTest : public ::testing::Test {
protected:
    // Variables common across tests
    NeoHookeanParams params;
    double F[3][3] = {}; // Deformation gradient
    int n_iter = 10;       // Number of random perturbations to test
    double rel_tol = 1e-3; // relative tolerance for comparing values
    double abs_tol = 1e-11; // absolute tolerance for comparing values
    double delta = 1e-7; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi

    // Add the test object
    TestNeoHookean* TestNH;

    // Setup method to initialize variables before each test
    void SetUp() override {

        // Set random values for the Neo-Hookean parameters between 1000 and 10000
        params.C10 = getRandomDouble(1000.0, 10000.0);

        // Initialize the test object
        TestNH = new TestNeoHookean(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestNH;
        TestNH = nullptr;
    }
};

// Test PK2 stress zero for F = I
TEST_F(NeoHookeanTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    double F[3][3] = {{1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0},
                       {0.0, 0.0, 1.0}};
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestNH->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress with finite difference for random F
TEST_F(NeoHookeanTest, TestPK2StressRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Compute reference PK2 stress with finite difference for random F
    create_random_F(F);
    double S_ref[3][3]; // PK2 stress
    TestNH->calcPK2StressFiniteDifference(F, delta, S_ref);

    // Check PK2 stress against reference value
    TestNH->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress consistent with strain energy for random F
TEST_F(NeoHookeanTest, TestPK2StressConsistentRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Check random F produces consistent PK2 stress
    create_random_F(F);
    TestNH->testPK2StressConsistentWithStrainEnergy(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// Test material elasticity consistent with PK2 stress
TEST_F(NeoHookeanTest, TestMaterialElasticityConsistentRandomF) {
    //verbose = true; // Show values of CC, dE, CCdE and dS

    // Check with random F
    create_random_F(F);
    TestNH->testMaterialElasticityConsistentWithPK2Stress(F, n_iter, rel_tol, abs_tol, delta, verbose);
}


// ----------------------------------------------------------------------------
// --------------------------- Mooney-Rivlin Material -------------------------
// ----------------------------------------------------------------------------

// Test fixture class for Mooney-Rivlin material model
class MooneyRivlinTest : public ::testing::Test {
protected:
    // Variables common across tests
    MooneyRivlinParams params;
    double F[3][3] = {}; // Deformation gradient
    int n_iter = 10;       // Number of random perturbations to test
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double abs_tol = 1e-11; // absolute tolerance for comparing values
    double delta = 1e-7; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi


    // Add the test object
    TestMooneyRivlin* TestMR;

    // Setup method to initialize variables before each test
    void SetUp() override {

        // Set random values for the Mooney-Rivlin parameters between 1000 and 10000
        params.C01 = getRandomDouble(1000.0, 10000.0);
        params.C10 = getRandomDouble(1000.0, 10000.0);

        // Initialize the test object
        TestMR = new TestMooneyRivlin(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestMR;
        TestMR = nullptr;
    }
};

// Test PK2 stress zero for F = I
TEST_F(MooneyRivlinTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    double F[3][3] = {1.0, 0.0, 0.0,
                      0.0, 1.0, 0.0,
                      0.0, 0.0, 1.0};
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestMR->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress with finite difference for random F
TEST_F(MooneyRivlinTest, TestPK2StressRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Compute reference PK2 stress with finite difference for random F
    create_random_F(F);
    double S_ref[3][3]; // PK2 stress
    TestMR->calcPK2StressFiniteDifference(F, delta, S_ref);

    // Check PK2 stress against reference value
    TestMR->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress consistent with strain energy for random F
TEST_F(MooneyRivlinTest, TestPK2StressConsistentRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Check random F produces consistent PK2 stress
    create_random_F(F);
    TestMR->testPK2StressConsistentWithStrainEnergy(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// Test material elasticity consistent with PK2 stress
TEST_F(MooneyRivlinTest, TestMaterialElasticityConsistentRandomF) {
    //verbose = true; // Show values of CC, dE, CCdE and dS

    // Check with random F
    create_random_F(F);
    TestMR->testMaterialElasticityConsistentWithPK2Stress(F, n_iter, rel_tol, abs_tol, delta, verbose);
}


// ----------------------------------------------------------------------------
// ------------------- Holzapfel-Gasser-Ogden Material -------------------------
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// ----------------------- Holzapfel-Ogden Material ---------------------------
// ----------------------------------------------------------------------------

// Test fixture class for Holzapfel-Ogden material model
class HolzapfelOgdenTest : public ::testing::Test {
protected:
    // Variables common across tests
    HolzapfelOgdenParams params;
    double F[3][3] = {}; // Deformation gradient
    int n_iter = 10;       // Number of random perturbations to test
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double abs_tol = 1e-11; // absolute tolerance for comparing values
    double delta = 1e-7; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi

    // Add the test object
    TestHolzapfelOgden* TestHO;

    // Setup method to initialize variables before each test
    void SetUp() override {

        // Set Holzapfel-Ogden parameters from cardiac benchmark paper
        params.a = 59.0; // Pa
        params.a_f = 18472.0; // Pa
        params.a_s = 2481.0; // Pa
        params.a_fs = 216.0; // Pa
        params.b = 8.023; // Pa
        params.b_f = 16.026; // Pa
        params.b_s = 11.12; // Pa
        params.b_fs = 11.436; // Pa
        params.kappa = 1.0e6; // Pa
        params.k = 100000.0; // Pa 

        // Set random values for f between 0 and 1 and normalize
        params.f[0] = getRandomDouble(0.0, 1.0);
        params.f[1] = getRandomDouble(0.0, 1.0);
        params.f[2] = getRandomDouble(0.0, 1.0);
        double norm_f = sqrt(params.f[0]*params.f[0] + params.f[1]*params.f[1] + params.f[2]*params.f[2]);
        params.f[0] /= norm_f; params.f[1] /= norm_f; params.f[2] /= norm_f;

        // Create s orthogonal to f
        if (fabs(params.f[0]) < 0.9) { // Check if f[0] is not the dominant component
            params.s[0] = 0;
            params.s[1] = params.f[2];
            params.s[2] = -params.f[1];
        } else { // If f[0] is the dominant component, use another approach
            params.s[0] = -params.f[2];
            params.s[1] = 0;
            params.s[2] = params.f[0];
        }

        // Normalize s
        double norm_s = sqrt(params.s[0]*params.s[0] + params.s[1]*params.s[1] + params.s[2]*params.s[2]);
        params.s[0] /= norm_s; params.s[1] /= norm_s; params.s[2] /= norm_s;

        // Check f.s = 0
        double dot_fs = params.f[0]*params.s[0] + params.f[1]*params.s[1] + params.f[2]*params.s[2];
        if (fabs(dot_fs) > 1e-6) {
            cout << "f.s = " << dot_fs << endl;
            cout << "f = [" << params.f[0] << ", " << params.f[1] << ", " << params.f[2] << "]" << endl;
            cout << "s = [" << params.s[0] << ", " << params.s[1] << ", " << params.s[2] << "]" << endl;
            throw runtime_error("f and s are not orthogonal");
        }

        // Flag to use full anisotropic invariants for strain energy computation
        params.full_anisotropic_invariants = false;


        // Initialize the test object
        TestHO = new TestHolzapfelOgden(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestHO;
        TestHO = nullptr;
    }
};

// Test PK2 stress zero for F = I
TEST_F(HolzapfelOgdenTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    double F[3][3] = {{1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0},
                       {0.0, 0.0, 1.0}};
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress for triaxial stretch
TEST_F(HolzapfelOgdenTest, TestPK2StressTriaxialStretch) {
    //verbose = true; // Show values of S and S_ref

    // Check triaxial stretch produces PK2 stress consistent with svFSI
    double F[3][3] = {{1.1, 0.0, 0.0},
                       {0.0, 1.2, 0.0},
                       {0.0, 0.0, 1.3}};
    
    // Compute reference PK2 stress with finite difference for triaxial stretch
    double S_ref[3][3]; // PK2 stress
    TestHO->calcPK2StressFiniteDifference(F, delta, S_ref);
    
    // Check PK2 stress against reference value
    TestHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress for triaxial compression
TEST_F(HolzapfelOgdenTest, TestPK2StressTriaxialCompression) {
    //verbose = true; // Show values of S and S_ref

    // Check triaxial compression produces PK2 stress consistent with svFSI
    double F[3][3] = {{0.9, 0.0, 0.0},
                       {0.0, 0.8, 0.0},
                       {0.0, 0.0, 0.7}};
    
    // Compute reference PK2 stress with finite difference for triaxial compression
    double S_ref[3][3]; // PK2 stress
    TestHO->calcPK2StressFiniteDifference(F, delta, S_ref);
    
    // Check PK2 stress against reference value
    TestHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress with finite difference for random F
TEST_F(HolzapfelOgdenTest, TestPK2StressRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Compute reference PK2 stress with finite difference for F = I + random perturbations
    create_random_perturbed_identity_F(F, 0.5);
    double S_ref[3][3]; // PK2 stress
    TestHO->calcPK2StressFiniteDifference(F, delta, S_ref);

    // Check PK2 stress against reference value
    TestHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress consistent with strain energy for random F
TEST_F(HolzapfelOgdenTest, TestPK2StressConsistentRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Check F = I + random perturbations produces consistent PK2 stress
    create_random_perturbed_identity_F(F, 0.5);
    TestHO->testPK2StressConsistentWithStrainEnergy(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// Test material elasticity consistent with PK2 stress
TEST_F(HolzapfelOgdenTest, TestMaterialElasticityConsistentRandomF) {
    //verbose = true; // Show values of CC, dE, CCdE and dS

    // Check F = I + random perturbations produces consistent PK2 stress
    create_random_perturbed_identity_F(F, 0.5);
    TestHO->testMaterialElasticityConsistentWithPK2Stress(F, n_iter, rel_tol, abs_tol, delta, verbose);
}






// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------
// ----------------------- Volumetric penalty models -------------------------
// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// --------------------------- Quadratic Volumetric Penalty Model ------------------------
// ----------------------------------------------------------------------------
// Test fixture class for Quadratic penalty model
class QuadraticVolumetricPenaltyTest : public ::testing::Test {
protected:
    // Variables common across tests
    VolumetricPenaltyParams params;
    double F[3][3] = {}; // Deformation gradient
    int n_iter = 10;       // Number of random perturbations to test
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double abs_tol = 1e-11; // absolute tolerance for comparing values
    double delta = 1e-7; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi

    // Add the test object
    TestQuadraticVolumetricPenalty* TestQVP;

    // Setup method to initialize variables before each test
    void SetUp() override {

        // Set random values for the Quadratic penalty parameters between 1000 and 10000
        params.kappa = getRandomDouble(1000.0, 10000.0);

        // Initialize the test object
        TestQVP = new TestQuadraticVolumetricPenalty(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestQVP;
        TestQVP = nullptr;
    }
};

// Test PK2 stress zero for F = I
TEST_F(QuadraticVolumetricPenaltyTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    double F[3][3] = {{1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0},
                       {0.0, 0.0, 1.0}};
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestQVP->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress zero for prescribed isochoric deformation
TEST_F(QuadraticVolumetricPenaltyTest, TestPK2StressPrescribedIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Check isochoric deformation produces zero PK2 stress
    double F[3][3] = {{1.1, 0.0, 0.0},
                       {0.0, 1.2, 0.0},
                       {0.0, 0.0, 1.0/(1.1*1.2)}};
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestQVP->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress zero for random isochoric deformation
TEST_F(QuadraticVolumetricPenaltyTest, TestPK2StressRandomIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Create random deformation gradient tensor
    create_random_F(F);

    // Make det(F) = 1
    double J = mat_fun_carray::mat_det(F);
    double J13 = pow(J, 1.0/3.0);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            F[i][j] /= J13;
        }
    }

    // Check random isochoric deformation produces zero PK2 stress
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestQVP->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress with finite difference for random F
TEST_F(QuadraticVolumetricPenaltyTest, TestPK2StressRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Compute reference PK2 stress with finite difference for random F
    create_random_F(F);
    double S_ref[3][3]; // PK2 stress
    TestQVP->calcPK2StressFiniteDifference(F, delta, S_ref);

    // Check PK2 stress against reference value
    TestQVP->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress consistent with strain energy for random F
TEST_F(QuadraticVolumetricPenaltyTest, TestPK2StressConsistentRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Check random F produces consistent PK2 stress
    create_random_F(F);
    TestQVP->testPK2StressConsistentWithStrainEnergy(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// Test material elasticity consistent with PK2 stress
TEST_F(QuadraticVolumetricPenaltyTest, TestMaterialElasticityConsistentRandomF) {
    //verbose = true; // Show values of CC, dE, CCdE and dS

    // Check with random F
    create_random_F(F);
    TestQVP->testMaterialElasticityConsistentWithPK2Stress(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// ----------------------------------------------------------------------------
// --------------------------- Simo-Taylor91 Volumetric Penalty Model ---------
// ----------------------------------------------------------------------------

// Test fixture class for Simo-Taylor91 penalty model
class SimoTaylor91VolumetricPenaltyTest : public ::testing::Test {
protected:
    // Variables common across tests
    VolumetricPenaltyParams params;
    double F[3][3] = {}; // Deformation gradient
    int n_iter = 10;       // Number of random perturbations to test
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double abs_tol = 1e-9; // absolute tolerance for comparing values
    double delta = 1e-7; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi

    // Add the test object
    TestSimoTaylor91VolumetricPenalty* TestST91;

    // Setup method to initialize variables before each test
    void SetUp() override {

        // Set random values for the Simo-Taylor91 penalty parameters between 1000 and 10000
        params.kappa = getRandomDouble(1000.0, 10000.0);

        // Initialize the test object
        TestST91 = new TestSimoTaylor91VolumetricPenalty(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestST91;
        TestST91 = nullptr;
    }
};

// Test PK2 stress zero for F = I
TEST_F(SimoTaylor91VolumetricPenaltyTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    double F[3][3] = {{1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0},
                       {0.0, 0.0, 1.0}};
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestST91->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress zero for prescribed isochoric deformation
TEST_F(SimoTaylor91VolumetricPenaltyTest, TestPK2StressPrescribedIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Check isochoric deformation produces zero PK2 stress
    double F[3][3] = {{1.1, 0.0, 0.0},
                       {0.0, 1.2, 0.0},
                       {0.0, 0.0, 1.0/(1.1*1.2)}};
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestST91->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress zero for random isochoric deformation
TEST_F(SimoTaylor91VolumetricPenaltyTest, TestPK2StressRandomIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Create random deformation gradient tensor
    create_random_F(F);

    // Make det(F) = 1
    double J = mat_fun_carray::mat_det(F);
    double J13 = pow(J, 1.0/3.0);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            F[i][j] /= J13;
        }
    }

    // Check random isochoric deformation produces zero PK2 stress
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestST91->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress with finite difference for random F
TEST_F(SimoTaylor91VolumetricPenaltyTest, TestPK2StressRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Compute reference PK2 stress with finite difference for random F
    create_random_F(F);
    double S_ref[3][3]; // PK2 stress
    TestST91->calcPK2StressFiniteDifference(F, delta, S_ref);

    // Check PK2 stress against reference value
    TestST91->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress consistent with strain energy for random F
TEST_F(SimoTaylor91VolumetricPenaltyTest, TestPK2StressConsistentRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Check random F produces consistent PK2 stress
    create_random_F(F);
    TestST91->testPK2StressConsistentWithStrainEnergy(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// Test material elasticity consistent with PK2 stress
TEST_F(SimoTaylor91VolumetricPenaltyTest, TestMaterialElasticityConsistentRandomF) {
    //verbose = true; // Show values of CC, dE, CCdE and dS

    // Check with random F
    create_random_F(F);
    TestST91->testMaterialElasticityConsistentWithPK2Stress(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// ----------------------------------------------------------------------------
// --------------------------- Miehe94 Volumetric Penalty Model ---------------
// ----------------------------------------------------------------------------
// Test fixture class for Miehe94 penalty model
class Miehe94VolumetricPenaltyTest : public ::testing::Test {
protected:
    // Variables common across tests
    VolumetricPenaltyParams params;
    double F[3][3] = {}; // Deformation gradient
    int n_iter = 10;       // Number of random perturbations to test
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double abs_tol = 1e-9; // absolute tolerance for comparing values
    double delta = 1e-7; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi

    // Add the test object
    TestMiehe94VolumetricPenalty* TestM94;

    // Setup method to initialize variables before each test
    void SetUp() override {

        // Set random values for the Miehe94 penalty parameters between 1000 and 10000
        params.kappa = getRandomDouble(1000.0, 10000.0);

        // Initialize the test object
        TestM94 = new TestMiehe94VolumetricPenalty(params);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the test object
        delete TestM94;
        TestM94 = nullptr;
    }
};

// Test PK2 stress zero for F = I
TEST_F(Miehe94VolumetricPenaltyTest, TestPK2StressIdentityF) {
    //verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    double F[3][3] = {{1.0, 0.0, 0.0},
                       {0.0, 1.0, 0.0},
                       {0.0, 0.0, 1.0}};
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestM94->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress zero for prescribed isochoric deformation
TEST_F(Miehe94VolumetricPenaltyTest, TestPK2StressPrescribedIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Check isochoric deformation produces zero PK2 stress
    double F[3][3] = {{1.1, 0.0, 0.0},
                       {0.0, 1.2, 0.0},
                       {0.0, 0.0, 1.0/(1.1*1.2)}};
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestM94->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress zero for random isochoric deformation
TEST_F(Miehe94VolumetricPenaltyTest, TestPK2StressRandomIsochoricDeformation) {
    //verbose = true; // Show values of S and S_ref

    // Create random deformation gradient tensor
    create_random_F(F);

    // Make det(F) = 1
    double J = mat_fun_carray::mat_det(F);
    double J13 = pow(J, 1.0/3.0);
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            F[i][j] /= J13;
        }
    }

    // Check random isochoric deformation produces zero PK2 stress
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestM94->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress with finite difference for random F
TEST_F(Miehe94VolumetricPenaltyTest, TestPK2StressRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Compute reference PK2 stress with finite difference for random F
    create_random_F(F);
    double S_ref[3][3]; // PK2 stress
    TestM94->calcPK2StressFiniteDifference(F, delta, S_ref);

    // Check PK2 stress against reference value
    TestM94->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress consistent with strain energy for random F
TEST_F(Miehe94VolumetricPenaltyTest, TestPK2StressConsistentRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Check random F produces consistent PK2 stress
    create_random_F(F);
    TestM94->testPK2StressConsistentWithStrainEnergy(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// Test material elasticity consistent with PK2 stress
TEST_F(Miehe94VolumetricPenaltyTest, TestMaterialElasticityConsistentRandomF) {
    //verbose = true; // Show values of CC, dE, CCdE and dS

    // Check with random F
    create_random_F(F);
    TestM94->testMaterialElasticityConsistentWithPK2Stress(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// ----------------------------------------------------------------------------
// ----------------------------------------------------------------------------


