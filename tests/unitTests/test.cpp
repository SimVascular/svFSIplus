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
    rel_tol = 1e-6; // relative tolerance for comparing S with reference value
    bool verbose = true; // Show values of S and S_ref

    // Check identity F produces zero PK2 stress
    create_identity_F(F);
    double S_ref[3][3] = {0}; // PK2 stress initialized to zero
    TestMR->testPK2StressAgainstReference(F, S_ref, rel_tol, verbose);
}

// Test PK2 stress with finite difference for random F
TEST_F(MooneyRivlinTest, TestPK2StressRandomF) {
    rel_tol = 1e-3; // relative tolerance for comparing S with values from svFSI
    delta = 1e-6; // perturbation scaling factor
    verbose = true; // Show values of S, dE, SdE and dPsi

    // Compute reference PK2 stress with finite difference for random F
    create_random_F(F);
    double S_ref[3][3]; // PK2 stress
    TestMR->calcPK2StressFiniteDifference(F, delta, S_ref);

    // Check PK2 stress against reference value
    TestMR->testPK2StressAgainstReference(F, S_ref, rel_tol, verbose);
}

// Test PK2 stress consistent with strain energy for random F
TEST_F(MooneyRivlinTest, TestPK2StressConsistentRandomF) {
    n_iter = 10;      // Number of random perturbations to test
    rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    delta = 1e-6; // perturbation scaling factor
    verbose = true; // Show values of S, dE, SdE and dPsi

    // Check random F produces consistent PK2 stress
    create_random_F(F);
    TestMR->testPK2StressConsistentWithStrainEnergy(F, n_iter, rel_tol, delta, verbose);
}

// Test material elasticity consistent with PK2 stress
TEST_F(MooneyRivlinTest, TestMaterialElasticityConsistentRandomF) {
    n_iter = 10;       // Number of random perturbations to test
    rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    delta = 1e-6; // perturbation scaling factor
    verbose = false; // Show values of CC, dE, CCdE and dS

    // Check with random F
    create_random_F(F);
    TestMR->testMaterialElasticityConsistentWithPK2Stress(F, n_iter, rel_tol, delta, verbose);
}


// ----------------------------------------------------------------------------
// ------------------- Holzapfel-Gasser-Ogden Material -------------------------
// ----------------------------------------------------------------------------


// ----------------------------------------------------------------------------
// ----------------------- Holzapfel-Ogden Material ---------------------------
// ----------------------------------------------------------------------------


