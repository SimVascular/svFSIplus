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


TEST(UnitTestIso_1, nHK) {
    // Step 1: define parameters
    auto matType = consts::ConstitutiveModelType::stIso_nHook;   // Material_model: options refer to consts.h 
    auto volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
    double E = 1e6;   // Elasticity_modulus
    double nu = 0.5;   // Poisson_ratio
    double pen = 4e9;   // Penalty_parameter
    double C01;   // additional parameter to C10 (optional)

    // Step 2: construct test object
    UnitTestIso nHK(matType, E, nu, volType, pen, C01);

    // Step 3: define the input 
    double F[3][3]; 
    memset(F, 0, sizeof(F));   // initialize
    F[0][0] = 1.0; F[1][1] = 1.0; F[2][2] = 1.0;   // set to Identity

    // Step 4: define the reference output 
    double S_ref[3][3], Dm_ref[6][6];
    memset(S_ref, 0, sizeof(S_ref));
    memset(Dm_ref, 0, sizeof(Dm_ref));

    // Step 5: run unit test
    nHK.runUnitTest(F, S_ref, Dm_ref);
      
}

TEST(UnitTestIso_2, MR) {
    // Step 1: define parameters
    auto matType = consts::ConstitutiveModelType::stIso_MR;   // Material_model: options refer to consts.h 
    auto volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
    double E = 1e6;   // Elasticity_modulus
    double nu = 0.495;   // Poisson_ratio
    double pen = 4e9;   // Penalty_parameter
    double C01 = 0.1;   // additional parameter to C10 (optional)

    // Step 2: construct test object
    UnitTestIso MR(matType, E, nu, volType, pen, C01);

    // Step 3: define the input 
    double F[3][3]; 
    memset(F, 0, sizeof(F));   // initialize
    F[0][0] = 1.0; F[1][1] = 1.0; F[2][2] = 1.0;   // set to Identity

    // Step 4: define the reference output 
    double S_ref[3][3], Dm_ref[6][6];
    memset(S_ref, 0, sizeof(S_ref));
    memset(Dm_ref, 0, sizeof(Dm_ref));

    // Step 5: run unit test
    MR.runUnitTest(F, S_ref, Dm_ref);
      
}

