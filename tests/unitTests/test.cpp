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
    double delta = 1e-6; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi

    // Add the test object
    TestNeoHookean* TestNH;

    // Setup method to initialize variables before each test
    void SetUp() override {

        // Seed random number generator
        //srand(static_cast<unsigned int>(time(0)));

        // Set random values for the Neo-Hookean parameter between 1000 and 10000
        params.C10 = 1000 + 9000 * (double)rand() / RAND_MAX;

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
    create_identity_F(F);
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
    double delta = 1e-6; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi


    // Add the test object
    TestMooneyRivlin* TestMR;



    // Setup method to initialize variables before each test
    void SetUp() override {

        // Seed random number generator
        //srand(static_cast<unsigned int>(time(0)));

        // Set random values for the Mooney-Rivlin parameters between 1000 and 10000
        params.C01 = 1000 + 9000 * (double)rand() / RAND_MAX;
        params.C10 = 1000 + 9000 * (double)rand() / RAND_MAX;

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
    create_identity_F(F);
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
    double delta = 1e-6; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi

    // Add the test object
    TestHolzapfelOgden* TestHO;

    // Setup method to initialize variables before each test
    void SetUp() override {

        // Seed random number generator
        //srand(static_cast<unsigned int>(time(0)));

        // Set random values for the Holzapfel-Ogden 'a' parameters between 1000 and 10000
        params.a = 1000 + 9000 * (double)rand() / RAND_MAX;
        params.a_f = 1000 + 9000 * (double)rand() / RAND_MAX;
        params.a_s = 1000 + 9000 * (double)rand() / RAND_MAX;
        params.a_fs = 1000 + 9000 * (double)rand() / RAND_MAX;

        // Set random values for the Holzapfel-Ogden 'b' parameters between 1 and 20
        params.b = 1 + 19 * (double)rand() / RAND_MAX;
        params.b_f = 1 + 19 * (double)rand() / RAND_MAX;
        params.b_s = 1 + 19 * (double)rand() / RAND_MAX;
        params.b_fs = 1 + 19 * (double)rand() / RAND_MAX;

        // Set volumetric penalty parameter between 10^5 and 10^7
        params.kappa = 1e5 + 9e6 * (double)rand() / RAND_MAX;

        // Set smoothed Heaviside parameter between 10 and 200
        params.k = 10 + 190 * (double)rand() / RAND_MAX;

        // Set random values for f and s fiber directions and normalize them
        params.f[0] = (double)rand() / RAND_MAX;
        params.f[1] = (double)rand() / RAND_MAX;
        params.f[2] = (double)rand() / RAND_MAX;
        params.s[0] = (double)rand() / RAND_MAX;
        params.s[1] = (double)rand() / RAND_MAX;
        params.s[2] = (double)rand() / RAND_MAX;
        double norm_f = sqrt(params.f[0]*params.f[0] + params.f[1]*params.f[1] + params.f[2]*params.f[2]);
        double norm_s = sqrt(params.s[0]*params.s[0] + params.s[1]*params.s[1] + params.s[2]*params.s[2]);
        params.f[0] /= norm_f; params.f[1] /= norm_f; params.f[2] /= norm_f;
        params.s[0] /= norm_s; params.s[1] /= norm_s; params.s[2] /= norm_s;


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
    verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    create_identity_F(F);
    double S_ref[3][3] = {}; // PK2 stress initialized to zero
    TestHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress with finite difference for random F
TEST_F(HolzapfelOgdenTest, TestPK2StressRandomF) {
    verbose = true; // Show values of S, dE, SdE and dPsi

    // Compute reference PK2 stress with finite difference for random F
    create_random_F(F);
    double S_ref[3][3]; // PK2 stress
    TestHO->calcPK2StressFiniteDifference(F, delta, S_ref);

    // Check PK2 stress against reference value
    TestHO->testPK2StressAgainstReference(F, S_ref, rel_tol, abs_tol, verbose);
}

// Test PK2 stress consistent with strain energy for random F
TEST_F(HolzapfelOgdenTest, TestPK2StressConsistentRandomF) {
    //verbose = true; // Show values of S, dE, SdE and dPsi

    // Check random F produces consistent PK2 stress
    create_random_F(F);
    TestHO->testPK2StressConsistentWithStrainEnergy(F, n_iter, rel_tol, abs_tol, delta, verbose);
}

// Test material elasticity consistent with PK2 stress
TEST_F(HolzapfelOgdenTest, TestMaterialElasticityConsistentRandomF) {
    //verbose = true; // Show values of CC, dE, CCdE and dS

    // Check with random F
    create_random_F(F);
    TestHO->testMaterialElasticityConsistentWithPK2Stress(F, n_iter, rel_tol, abs_tol, delta, verbose);
}


