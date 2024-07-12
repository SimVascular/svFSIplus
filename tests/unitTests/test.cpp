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


// Test fixture class for Mooney-Rivlin material model
class MooneyRivlinTest : public ::testing::Test {
protected:
    // Variables common across tests
    MooneyRivlinParams params;

    // Add the test object
    TestMooneyRivlin* TestMR;

    // Setup method to initialize variables before each test
    void SetUp() override {

        // Set random values for the Mooney-Rivlin parameters between 1000 and 10000
        srand(static_cast<unsigned int>(time(0))); // seed random number generator
        params.C01 = 1000 + 9000 * (double)rand() / RAND_MAX;
        params.C10 = 1000 + 9000 * (double)rand() / RAND_MAX;

        std::cout << "C01 = " << params.C01 << ", C10 = " << params.C10 << std::endl;

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


// Test PK2 stress consistent with strain energy
TEST_F(MooneyRivlinTest, TestPK2Stress) {
    int n_iter = 100;
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double delta = 1e-8; // perturbation scaling factor
    bool verbose = false; // Show values of S, dE, SdE and dPsi
    TestMR->testPK2StressConsistentWithStrainEnergy(n_iter, rel_tol, delta, verbose);
}

// Test material elasticity consistent with PK2 stress
TEST_F(MooneyRivlinTest, TestMaterialElasticity) {
    int n_iter = 100;
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double delta = 1e-8; // perturbation scaling factor
    bool verbose = false; // Show values of CC, dE, CCdE and dS
    TestMR->testMaterialElasticityConsistentWithPK2Stress(n_iter, rel_tol, delta, verbose);
}



