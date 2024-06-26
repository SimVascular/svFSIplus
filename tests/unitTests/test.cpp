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


// Creates a random deformation gradient F with values between 0.1 and 10.0,
// and det(F) > 0
template<int N>
void create_random_F(double F[N][N]) {
    // Initialize random seed
    srand(static_cast<unsigned int>(time(0))); // seed random number generator

    // Create a random deformation gradient with values between 0.1 and 10.0, 
    // and det(F) > 0
    double J = -1.0;
    while (J < 0) {
        for (int i = 0; i < N; i++) {
            for (int J = 0; J < N; J++) {
                F[i][J] = 0.1 + (10.0 - 0.1) * rand() / RAND_MAX;
            }
        }
        J = mat_fun_carray::mat_det<N>(F);
    }
}

// Function to compute the strain energy density Psi for the Mooney-Rivlin material model
// given the material parameters and deformation gradient F
// Templated for the number of spatial dimensions N
template<int N>
double Psi_MR_ref(double C01, double C10, double F[N][N]) {          

    // ----------------------------------------------------------------
    // ------------------------ Setup ---------------------------------
    // ----------------------------------------------------------------

    // Cast N to a double
    double N_d = static_cast<double>(N);

    // Jacobian of F
    double J = mat_fun_carray::mat_det<N>(F);

    // Transpose of F
    double F_T[N][N];
    mat_fun_carray::transpose<N>(F, F_T);

    // Right Cauchy-Green deformation tensor
    double C[N][N];
    mat_fun_carray::mat_mul<N>(F_T, F, C);

    // Right Cauchy-Green deformation tensor squared
    double C2[N][N];
    mat_fun_carray::mat_mul<N>(C, C, C2);

    // Modified right Cauchy-Green deformation tensor
    double C_bar[N][N];
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C_bar[i][j] = pow(J, (-2.0/N_d)) * C[i][j];
        }
    }

    // Modified right Cauchy-Green deformation tensor squared
    double C_bar2[N][N];
    mat_fun_carray::mat_mul<N>(C_bar, C_bar, C_bar2);

    // Invariants of C
    double I1 = mat_fun_carray::mat_trace<N>(C);
    double I2 = 0.5 * (pow(I1, 2) - mat_fun_carray::mat_trace<N>(C2));

    // Invariants of C_bar
    double Ib1 = mat_fun_carray::mat_trace<N>(C_bar);
    double Ib2 = 0.5 * (pow(Ib1, 2) - mat_fun_carray::mat_trace<N>(C_bar2));

    // Check that invariants satisfy expected relationship
    EXPECT_NEAR( pow(J, (-2.0/3.0)) * I1, Ib1, 1e-9 * Ib1);
    EXPECT_NEAR( pow(J, (-4.0/3.0)) * I2, Ib2, 1e-9 * Ib2);

    // ----------------------------------------------------------------
    // --------------- Compute strain energy density ------------------
    // ----------------------------------------------------------------

    // Strain energy density for Mooney-Rivlin material model
    double Psi_iso = C10 * (Ib1 - 3) + C01 * (Ib2 - 3);

    return Psi_iso;

}

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
    double F[3][3] = {};
    F[0][0] = 1.0; F[1][1] = 1.0; F[2][2] = 1.0;   // set to Identity

    // Step 4: define the reference output 
    double S_ref[3][3] = {};
    double Dm_ref[6][6] = {};

    // Step 5: run unit test
    double rel_tol = 1e-6;
    nHK.compare_S_Dm(F, S_ref, Dm_ref, rel_tol);
      
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
    double F[3][3] = {};
    F[0][0] = 1.0; F[1][1] = 1.0; F[2][2] = 1.0;   // set to Identity

    // Step 4: define the reference output 
    double S_ref[3][3] = {};
    double Dm_ref[6][6] = {};

    // Step 5: run unit test
    double rel_tol = 1e-6;
    MR.compare_S_Dm(F, S_ref, Dm_ref, rel_tol);
      
}

TEST(TestMR_S, 3D) {
    // -------------------------
    // Step 1: define parameters
    // -------------------------

    auto matType = consts::ConstitutiveModelType::stIso_MR;   // Material_model: options refer to consts.h 
    auto volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
    double E_mod = 1e6;   // Elasticity_modulus
    double nu = 0.495;   // Poisson_ratio
    double pen = 0.0;   // Penalty_parameter
    double C01 = 0.1;   // additional parameter to C10 (optional)


    // Compute derived material parameters
    double mu  = 0.5 * E_mod / (1.0 + nu);                    // Shear_modulus
    double C10 = 0.5 * mu - C01;                 

    // -------------------------
    // Step 2: construct test object
    // -------------------------
    UnitTestIso MR(matType, E_mod, nu, volType, pen, C01);

    // -------------------------
    // Step 3: create random deformation gradient
    // -------------------------

    double F[3][3];
    create_random_F(F);

    // -------------------------
    // Step 4: compute the reference S_ij using finite difference
    // -------------------------

    // Compute strain energy density given F
    double Psi = Psi_MR_ref<3>(C01, C10, F);

    // Compute 1st PK stress P_ref_iJ = dPsi / dF[i][J] using finite difference, component by component
    double P_ref[3][3] = {};
    double delta = 1e-9; // perturbation size
    double F_tilde[3][3]; // perturbed deformation gradient
    for (int i = 0; i < 3; i++) {
        for (int J = 0; J < 3; J++) {
            // Perturb the iJ-th component of F by delta * rand()
            for (int k = 0; k < 3; k++) {
                for (int l = 0; l < 3; l++) {
                    F_tilde[k][l] = F[k][l];
                }
            }
            F_tilde[i][J] += delta;

            // Compute Psi_MR_ref for perturbed deformation gradient
            double Psi_tilde = Psi_MR_ref<3>(C01, C10, F_tilde);

            // Compute differences in Psi and E
            double dPsi = Psi_tilde - Psi;
            double dF_iJ = delta;

            // Compute P_ref[i][J] = dPsi / dF[i][J]
            P_ref[i][J] = dPsi / dF_iJ;
        }
    }

    // Compute S_ref = F^-1 * P_ref
    double S_ref[3][3] = {};
    double F_inv[3][3];
    mat_fun_carray::mat_inv<3>(F, F_inv);
    mat_fun_carray::mat_mul<3>(F_inv, P_ref, S_ref);


    // TODO: Compute Dm_ref
    double Dm_ref[6][6] = {};

    // Step 5: Compare S and Dm from svFSI with reference solutions
    double rel_tol = 1e-4;
    MR.compare_S_Dm(F, S_ref, Dm_ref, rel_tol);
}


TEST(TestMR_S_Dm, 3D) {
    // -------------------------
    // Step 1: define parameters
    // -------------------------
    auto matType = consts::ConstitutiveModelType::stIso_MR;   // Material_model: options refer to consts.h 
    auto volType = consts::ConstitutiveModelType::stVol_ST91;   // Dilational_penalty_model
    double E_mod = 1e6;   // Elasticity_modulus
    double nu = 0.495;   // Poisson_ratio
    double pen = 0.0;   // Penalty_parameter
    double C01 = 0.1;   // additional parameter to C10 (optional)
    // Compute derived material parameters
    double mu  = 0.5 * E_mod / (1.0 + nu);                    // Shear_modulus
    double C10 = 0.5 * mu - C01;    

    // -------------------------
    // Step 2: construct test object
    // -------------------------
    UnitTestIso MR(matType, E_mod, nu, volType, pen, C01);

    // -------------------------
    // Step 3: create random deformation gradient
    // -------------------------
    double F[3][3];
    create_random_F(F);

    // -------------------------
    // Step 4: compute Psi and S_ij for F
    // -------------------------
    double Psi = Psi_MR_ref<3>(C01, C10, F);

    double S[3][3], Dm[6][6];
    MR.get_pk2cc(F, S, Dm); // S from svFSI

    // -------------------------
    // Step 4: generate random dF and check that dPsi and dS are consistent
    // -------------------------
    int n_iter = 10;
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double delta = 1e-6; // perturbation scaling factor
    double dF[3][3]; // perturbation to deformation gradient
    double F_tilde[3][3]; // perturbed deformation gradient
    
    for (int i = 0; i < n_iter; i++) {
        // Perturb the deformation gradient
        double dF[3][3] = {};
        for (int i = 0; i < 3; i++) {
            for (int J = 0; J < 3; J++) {
                dF[i][J] = delta * (2.0 * rand() / RAND_MAX - 1.0); // (random number between -1 and 1) * delta
                F_tilde[i][J] = F[i][J] + dF[i][J]; // perturbed deformation gradient
            }
        }

        // Compute Psi and S for perturbed deformation gradient
        double Psi_tilde = Psi_MR_ref<3>(C01, C10, F_tilde);
        double S_tilde[3][3], Dm_tilde[6][6];
        MR.get_pk2cc(F_tilde, S_tilde, Dm_tilde);

        // Compute dPsi and dS
        double dPsi = Psi_tilde - Psi;

        double dS[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dS[i][j] = S_tilde[i][j] - S[i][j];
            }
        }

        // Check that S and Dm are consistent with the computed dPsi and dS
        MR.check_consistent_S_Dm(F, dF, dPsi, dS, rel_tol);
    }
}



