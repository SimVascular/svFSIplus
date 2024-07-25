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

// --------------------------------------------------------------
// To run the tests in test.cpp
// 1.  Navigate to <svFSIplus_root_directory>/build/svFSI-build/Source/svFSI
// 2.  Run `make` to build the tests
// 3.  Run `ctest --verbose` to run the tests


// --------------------------------------------------------------
// To add a new material model to test:
// In test.h:
// 1. Create a new material parameters class derived from MatParams (e.g. NeoHookeanParams)
// 2. Create a new material model test class derived from TestMaterialModel (e.g. TestNeoHookean)
// 3. Implement required functions in the material model test class:
//    - Constructor: Sets the material model type and parameters for svFSIplus (e.g. TestNeoHookean())
//    - printMaterialParameters(): Prints the material parameters
//    - computeStrainEnergy(): Computes the strain energy density function
// In test.cpp:
// 4. Create a new text fixture class derived from ::testing::Test (e.g. NeoHookeanTest)
//    - In this you set the values of the material parameters for testing
// 5. Add tests for your new material model (e.g. TEST_F(NeoHookeanTest, TestPK2StressIdentityF))
// --------------------------------------------------------------


#include <stdlib.h>
#include <iostream>
#include <random>
#include <chrono>
#include "gtest/gtest.h"   // include GoogleTest
#include "mat_fun.h"
#include "mat_fun_carray.h"
#include "mat_models.h"
#include "mat_models_carray.h"

// Class to contain material parameters
class MatParams {
public:
    virtual ~MatParams() {} // Virtual destructor for proper cleanup
};

// Class to contain Neo-Hookean material parameters
class NeoHookeanParams : public MatParams {
public:
    double C10;

    // Default constructor
    NeoHookeanParams() : C10(0.0) {}

    // Constructor with parameters
    NeoHookeanParams(double c10) : C10(c10) {}

};

// Class to contain Mooney-Rivlin material parameters
class MooneyRivlinParams : public MatParams {
public:
    double C01;
    double C10;

    // Default constructor
    MooneyRivlinParams() : C01(0.0), C10(0.0) {}

    // Constructor with parameters
    MooneyRivlinParams(double c01, double c10) : C01(c01), C10(c10) {}

};

// Class to contain Holzapfel-Ogden material parameters
class HolzapfelOgdenParams : public MatParams {
public:
    double a;    
    double b;
    double a_f;
    double b_f;
    double a_s;
    double b_s;
    double a_fs;
    double b_fs;
    double kappa;
    double f[3];    // Fiber direction
    double s[3];    // Sheet direction

    double k; // Smoothed Heaviside function parameter

    bool full_anisotropic_invariants; // Flag to use full anisotropic invariants in strain energy density function

    // Default constructor
    HolzapfelOgdenParams() : a(0.0), b(0.0), a_f(0.0), b_f(0.0), a_s(0.0), b_s(0.0), a_fs(0.0), b_fs(0.0), kappa(0.0), k(0.0) {
        for (int i = 0; i < 3; i++) {
            f[i] = 0.0;
            s[i] = 0.0;
        }
    }

    // Constructor with parameters
    HolzapfelOgdenParams(double a, double b, double a_f, double b_f, double a_s, double b_s, double a_fs, double b_fs, double kappa, double k, double f[3], double s[3]) : a(a), b(b), a_f(a_f), b_f(b_f), a_s(a_s), b_s(b_s), a_fs(a_fs), b_fs(b_fs), kappa(kappa), k(k) {
        for (int i = 0; i < 3; i++) {
            this->f[i] = f[i];
            this->s[i] = s[i];
        }
    }
};

// --------------------------------------------------------------
// ---------------------- Helper functions ----------------------
// --------------------------------------------------------------

// Creates an identity deformation gradient F
template<int N>
void create_identity_F(double F[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int J = 0; J < N; J++) {
            F[i][J] = (i == J);
        }
    }
}

// Inline getRandomDouble function
inline double getRandomDouble(double min, double max) {
    // Uncomment to use a random seed
    //unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
    unsigned int seed = 42;
    static std::default_random_engine engine(seed);
    std::uniform_real_distribution<double> distribution(min, max);
    return distribution(engine);
}

// Creates a random deformation gradient F with values between min and max,
// and det(F) > 0
template<int N>
void create_random_F(double F[N][N], double min=0.1, double max=10.0) {
    // Create a random deformation gradient with values between min and max, 
    // and det(F) > 0
    double J = -1.0;
    while (J < 0) {
        for (int i = 0; i < N; i++) {
            for (int J = 0; J < N; J++) {
                F[i][J] = getRandomDouble(min, max);
            }
        }
        J = mat_fun_carray::mat_det<N>(F);
    }
}

// Creates a deformation matrix F of random deviations from the identity matrix
template<int N>
void create_random_perturbed_identity_F(double F[N][N], double max_deviation) {
    // Create a random deformation gradient with values perturbed from the identity matrix, 
    // and det(F) > 0
    double J = -1.0;
    while (J < 0) {
        for (int i = 0; i < N; i++) {
            for (int J = 0; J < N; J++) {
                F[i][J] = (i == J) + max_deviation * getRandomDouble(-1.0, 1.0);
            }
        }
        J = mat_fun_carray::mat_det<N>(F);
    }
}

// Perturb the deformation gradient F by delta times a random number between -1 and 1
// and store the perturbed deformation gradient in F_tilde
template<int N>
void perturb_random_F(const double F[N][N], const double delta, double F_tilde[N][N]) {

    // Perturb the deformation gradient and store in F_tilde
    double dF_iJ;
    for (int i = 0; i < N; i++) {
        for (int J = 0; J < N; J++) {
            dF_iJ = delta * getRandomDouble(-1.0, 1.0);
            F_tilde[i][J] = F[i][J] + dF_iJ; // perturbed deformation gradient
        }
    }
}

// Compute Jacobian J, right Cauchy-Green deformation tensor C, and Green-Lagrange
// strain tensor E from the deformation gradient F
template<int N>
void calc_JCE(const double F[N][N], double &J, double C[N][N], double E[N][N]) {
    // Compute Jacobian of F
    J = mat_fun_carray::mat_det<N>(F);

    // Compute transpose of F
    double F_T[N][N];
    mat_fun_carray::transpose<N>(F, F_T);

    // Compute right Cauchy-Green deformation tensor
    mat_fun_carray::mat_mul<N>(F_T, F, C);

    // Compute Green-Lagrange strain tensor
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            E[i][j] = 0.5 * (C[i][j] - (i == j));
        }
    }
}

// Structure to store solid mechanics terms used to compute strain energy density
// functions.
template<int N>
struct solidMechanicsTerms {
    double J;
    double C[N][N];
    double E[N][N];
    double E2[N][N];
    double C_bar[N][N];
    double I1, I2, Ib1, Ib2;
};

// Function to compute the solid mechanics terms used to compute strain energy density
// functions 
template<int N>
solidMechanicsTerms<N> calcSolidMechanicsTerms(const double F[N][N]) {
    solidMechanicsTerms<N> out;

    const double N_d = static_cast<double>(N); // Convert N to double for calculations

    // Jacobian of F
    out.J = mat_fun_carray::mat_det<N>(F);

    // Transpose of F
    double F_T[N][N];
    mat_fun_carray::transpose<N>(F, F_T);

    // Right Cauchy-Green deformation tensor
    mat_fun_carray::mat_mul<N>(F_T, F, out.C);

    // Right Cauchy-Green deformation tensor squared
    double C2[N][N];
    mat_fun_carray::mat_mul<N>(out.C, out.C, C2);

    // Green-Lagrange strain tensor
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            out.E[i][j] = 0.5 * (out.C[i][j] - (i == j));
        }
    }

    // Green-Lagrange strain tensor squared
    mat_fun_carray::mat_mul<N>(out.E, out.E, out.E2);

    // Modified right Cauchy-Green deformation tensor
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            out.C_bar[i][j] = pow(out.J, (-2.0/N_d)) * out.C[i][j];
        }
    }

    // Modified right Cauchy-Green deformation tensor squared
    double C_bar2[N][N];
    mat_fun_carray::mat_mul<N>(out.C_bar, out.C_bar, C_bar2);

    // Invariants of C
    out.I1 = mat_fun_carray::mat_trace<N>(out.C);
    out.I2 = 0.5 * (pow(out.I1, 2) - mat_fun_carray::mat_trace<N>(C2));

    // Invariants of C_bar
    out.Ib1 = mat_fun_carray::mat_trace<N>(out.C_bar);
    out.Ib2 = 0.5 * (pow(out.Ib1, 2) - mat_fun_carray::mat_trace<N>(C_bar2));

    // Check that invariants satisfy expected relationship
    EXPECT_NEAR( pow(out.J, (-2.0/3.0)) * out.I1, out.Ib1, 1e-9 * out.Ib1);
    EXPECT_NEAR( pow(out.J, (-4.0/3.0)) * out.I2, out.Ib2, 1e-9 * out.Ib2);

    return out;
}


// --------------------------------------------------------------
// -------------------- Mock svFSIplus object -------------------
// --------------------------------------------------------------


class MockCepMod : public CepMod {
public:
    MockCepMod() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockdmnType : public dmnType {
public:
    MockdmnType() {
        // initialize if needed 
    }
    // MockstModelType mockStM;
    // Mock methods if needed
};
class MockmshType : public mshType {
public:
    MockmshType() {
        // initialize if needed 
    }
    // Mock methods if needed
};
class MockeqType : public eqType {
public:
    MockeqType() {
        // initialize if needed 
    }
    MockdmnType mockDmn;
    // Mock methods if needed
};
class MockComMod : public ComMod {
public:
    MockComMod() {
        // initialize if needed 
        nsd = 3;
    }
    MockeqType mockEq;
    MockmshType mockMsh;
    // Mock methods if needed
};


// --------------------------------------------------------------
// ------------------ Test Material Model Classes ---------------
// --------------------------------------------------------------

// Class for testing material models in svFSI
class TestMaterialModel {
public:
    MockComMod com_mod;
    MockCepMod cep_mod;
    int nFn;
    Array<double> fN;
    double ya_g;

    TestMaterialModel(const consts::ConstitutiveModelType matType, const consts::ConstitutiveModelType penType) {
        int nsd = com_mod.nsd;
        mat_fun_carray::ten_init(nsd);                        // initialize tensor index pointer

        // Set material and penalty models
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.isoType = matType;                            // Mat_model
        dmn.stM.volType = penType;                            // Dilational_penalty_model
        
        // Initialize fibers and other material parameters
        nFn = 2;                          // Number of fiber directions
        fN = Array<double>(nsd, nFn);     // Fiber directions array (initialized to zeros)
        ya_g = 0.0;                       // ?

    // Material parameters are set in each derived class
    }

    // Pure virtual method to print material parameters
    virtual void printMaterialParameters() = 0;

    // Pure virtual method for computing Strain Energy
    virtual double computeStrainEnergy(const double F[3][3]) = 0;

    // Function to get S and Dm for a given F, from the get_pk2cc function in mat_models_carray.h
    void get_pk2cc(double F[3][3], double S[3][3], double Dm[6][6]) {
        auto &dmn = com_mod.mockEq.mockDmn;

        mat_models_carray::get_pk2cc<3>(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
    }

    // Function to compute the PK2 stress tensor S(F) from the strain energy density Psi(F)
    // using finite differences
    // Analytically, we should have S = dPsi/dE. Since we have Psi(F), we cannot
    // directly compute S. Instead, we compute S = F^-1 * P, where P = dPsi/dF 
    // is computed using finite differences in each component of F.
    // ARGS:
    // - F: Deformation gradient
    // - delta: Perturbation scaling factor
    // - S: PK2 stress tensor
    // 
    // Pseudocode:
    // - Compute strain energy density Psi(F)
    // - For each component of F, F[i][J]
    //      - Perturb F[i][J] by delta to get F_tilde
    //      - Compute Psi(F_tilde)
    //      - Compute dPsi = Psi(F_tilde) - Psi(F)
    //      - Compute P[i][J] = dPsi / delta
    // - Compute S = F^-1 * P
    template<int N>
    void calcPK2StressFiniteDifference(const double F[N][N], double delta, double (&S)[N][N]) {
        // Compute strain energy density given F
        double Psi = computeStrainEnergy(F);

        // Compute 1st PK stress P_iJ = dPsi / dF[i][J] using finite difference, component by component
        double P[3][3] = {};
        double F_tilde[N][N]; // perturbed deformation gradient
        for (int i = 0; i < N; i++) {
            for (int J = 0; J < N; J++) {
                // Perturb the iJ-th component of F by delta
                for (int k = 0; k < N; k++) {
                    for (int l = 0; l < N; l++) {
                        F_tilde[k][l] = F[k][l];
                    }
                }
                F_tilde[i][J] += delta;

                // Compute Psi_MR for perturbed deformation gradient
                double Psi_tilde = computeStrainEnergy(F_tilde);

                // Compute differences in Psi
                double dPsi = Psi_tilde - Psi;

                // Compute P[i][J] = dPsi / dF[i][J]
                P[i][J] = dPsi / delta;
            }
        }

        // Compute S_ref = F^-1 * P_ref
        double F_inv[N][N];
        mat_fun_carray::mat_inv<N>(F, F_inv);
        mat_fun_carray::mat_mul<N>(F_inv, P, S);

        
    }

    // Test the consistency of the PK2 stress tensor S(F) from get_pk2cc() with the strain 
    // energy density Psi(F) provided by the user.
    // Analytically, we should have S = dPsi/dE. We are checking whether
    // S:dE = dPsi, where dE and dPsi are computed using finite differences in F.
    //
    // ARGS:
    // - F: Deformation gradient
    // - n_iter: Number of random perturbations to test
    // - rel_tol: Relative tolerance for comparing dPsi and S:dE
    // - delta: Perturbation scaling factor
    // - verbose: Show values of S, dE, SdE and dPsi
    //
    // Psuedocode:
    // - Compute Psi(F)
    // - Compute S(F) from get_pk2cc()
    // - For many random dF
    //      - Compute dPsi = Psi(F + dF) - Psi(F)
    //      - Compute dE from dF
    //      - Check that S:dE = dPsi
    void testPK2StressConsistentWithStrainEnergy(double F[3][3], int n_iter, double rel_tol, double abs_tol, double delta, bool verbose = false) {
        // Compute E from F
        double J, C[3][3], E[3][3];
        calc_JCE(F, J, C, E);

        // Compute Psi(F)
        double Psi = computeStrainEnergy(F);

        // Compute S(F) from svFSIplus
        double S[3][3], Dm[6][6];
        get_pk2cc(F, S, Dm);

        // Generate many random dF and check that S:dE = dPsi
        // S was obtained from get_pk2cc(), and dPsi = Psi(F + dF) - Psi(F)
        double dPsi, dE[3][3], SdE;
        double F_tilde[3][3], J_tilde, C_tilde[3][3], E_tilde[3][3];
        for (int i = 0; i < n_iter; i++) {
            // Perturb the deformation gradient
            perturb_random_F(F, delta, F_tilde);

            // Compute perturbed E_tilde from F_tilde
            calc_JCE(F_tilde, J_tilde, C_tilde, E_tilde);

            // Compute Psi(F_tilde)
            double Psi_tilde = computeStrainEnergy(F_tilde);

            // Compute dPsi = Psi(F_tilde) - Psi(F)
            dPsi = Psi_tilde - Psi;

            // Compute dE = E(F_tilde) - E(F)
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    dE[i][j] = E_tilde[i][j] - E[i][j];
                }
            }


            // Compute S:dE
            double SdE = mat_fun_carray::mat_ddot<3>(S, dE);

            // Check that S:dE = dPsi
            EXPECT_NEAR(SdE, dPsi, fmax(abs_tol, rel_tol * fabs(dPsi)));
            
            // Print results if verbose
            if (verbose) {
                std::cout << "Iteration " << i << ":" << std::endl;

                printMaterialParameters();

                std::cout << "F =" << std::endl;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        std::cout << F[i][j] << " ";
                    }
                    std::cout << std::endl;
                }

                std::cout << "S =" << std::endl;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        std::cout << S[i][j] << " ";
                    }
                    std::cout << std::endl;
                }

                std::cout << "dE =" << std::endl;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        std::cout << dE[i][j] << " ";
                    }
                    std::cout << std::endl;
                }

                std::cout << "SdE = " << SdE << ", dPsi = " << dPsi << std::endl;
                std::cout << std::endl;
            }
        }
    }

    // Test the consistency of the material elasticity tensor CC(F) from get_pk2cc() with the
    // PK2 stress tensor S(F) from get_pk2cc()
    // Analytically, we should have CC:dE = dS. We are checking whether
    // CC:dE = dS, where dE and dS are computed using finite differences in F.
    //
    // ARGS:
    // - F: Deformation gradient
    // - n_iter: Number of random perturbations to test
    // - rel_tol: Relative tolerance for comparing dS and CC:dE
    // - delta: Perturbation scaling factor
    // - verbose: Show values of CC, dE, CCdE and dS
    //
    // Psuedocode:
    // - Compute S(F) and CC(F) from get_pk2cc()
    // - For many random dF
    //      - Compute S(F + dF) from get_pk2cc()
    //      - Compute dS = S(F + dF) - S(F)
    //      - Compute dE from dF
    //      - Check that CC:dE = dS
    void testMaterialElasticityConsistentWithPK2Stress(double F[3][3], int n_iter, double rel_tol, double abs_tol, double delta, bool verbose = false) {

        // Compute E from F
        double J, C[3][3], E[3][3];
        calc_JCE(F, J, C, E);

        // Compute S_ij(F)
        // Compute CC_ijkl(F). 
        // CC is provided in Voigt notation as Dm, and we will convert it to CC
        double S[3][3], Dm[6][6];
        get_pk2cc(F, S, Dm); // S from svFSI

        // Calculate CC from Dm
        double CC[3][3][3][3];
        mat_models_carray::voigt_to_cc_carray<3>(Dm, CC);

        // ------- Ancillary test ---------
        // Calculate Dm_check from CC
        double Dm_check[6][6];
        mat_models_carray::cc_to_voigt_carray<3>(CC, Dm_check);

        // Check that Dm_check = Dm, for sanity
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 6; j++) {
                EXPECT_NEAR(Dm_check[i][j], Dm[i][j], abs_tol);
            }
        }
        // -------------------------------

        // Generate many random dF and check that CC:dE = dS
        // CC was obtained from get_pk2cc(), and dS = S(F + dF) - S(F), 
        // where S is also obtained from get_pk2cc()
        double dS[3][3], dE[3][3], CCdE[3][3];
        double F_tilde[3][3], J_tilde, C_tilde[3][3], E_tilde[3][3];
        double S_tilde[3][3], Dm_tilde[6][6];
        
        // Loop over many random perturbations to the deformation gradient
        for (int i = 0; i < n_iter; i++) {
            // Perturb the deformation gradient
            perturb_random_F<3>(F, delta, F_tilde);

            // Compute perturbed E_tilde from F_tilde
            calc_JCE(F_tilde, J_tilde, C_tilde, E_tilde);

            // Compute dE
            for (int i = 0; i < 3; i++) {
                for (int J = 0; J < 3; J++) {
                    dE[i][J] = E_tilde[i][J] - E[i][J];
                }
            }

            // Compute perturbed S_tilde with perturbed deformation gradient
            get_pk2cc(F_tilde, S_tilde, Dm_tilde);

            // Compute dS
            double dS[3][3];
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    dS[i][j] = S_tilde[i][j] - S[i][j];
                }
            }

            // Check that CC_ijkl dE_kl = dS_ij
            mat_fun_carray::ten_mat_ddot<3>(CC, dE, CCdE);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    EXPECT_NEAR(CCdE[i][j], dS[i][j], fmax(abs_tol, rel_tol * fabs(dS[i][j])));
                }
            }

            // Print results if verbose
            if (verbose) {
                std::cout << "Iteration " << i << ":" << std::endl;

                printMaterialParameters();

                std::cout << "F =" << std::endl;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        std::cout << F[i][j] << " ";
                    }
                    std::cout << std::endl;
                }

                std::cout << "CC =" << std::endl;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        for (int k = 0; k < 3; k++) {
                            for (int l = 0; l < 3; l++) {
                                std::cout << CC[i][j][k][l] << " ";
                            }
                            std::cout << std::endl;
                        }
                    }
                    std::cout << std::endl;
                }

                std::cout << "dE =" << std::endl;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        std::cout << dE[i][j] << " ";
                    }
                    std::cout << std::endl;
                }

                std::cout << "dS =" << std::endl;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        std::cout << dS[i][j] << " ";
                    }
                    std::cout << std::endl;
                }

                std::cout << "CCdE =" << std::endl;
                for (int i = 0; i < 3; i++) {
                    for (int j = 0; j < 3; j++) {
                        std::cout << CCdE[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }
        }
    }

    // Function to compare PK2 stress with reference solution
    // ARGS:
    // - F: Deformation gradient
    // - S_ref: Reference solution for PK2 stress
    // - rel_tol: Relative tolerance for comparing S with S_ref
    void testPK2StressAgainstReference(double F[3][3], double S_ref[3][3], double rel_tol, double abs_tol, bool verbose = false) {
        // Compute S(F) from get_pk2cc()
        double S[3][3], Dm[6][6];
        get_pk2cc(F, S, Dm);

        // Compare S with reference solution
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                EXPECT_NEAR(S[i][j], S_ref[i][j], fmax(abs_tol, rel_tol * fabs(S_ref[i][j])));
            }
        }

        // Print results if verbose
        if (verbose) {
            printMaterialParameters();

            std::cout << "F =" << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    std::cout << F[i][j] << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "S =" << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    std::cout << S[i][j] << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "S_ref =" << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    std::cout << S_ref[i][j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }

    // Function to compare material elasticity tensor with reference solution
    // ARGS:
    // - F: Deformation gradient
    // - CC_ref: Reference solution for material elasticity tensor
    // - rel_tol: Relative tolerance for comparing CC with CC_ref
    void testMaterialElasticityAgainstReference(double F[3][3], double CC_ref[3][3][3][3], double rel_tol, double abs_tol, bool verbose = false) {
        // Compute CC(F) from get_pk2cc()
        double S[3][3], Dm[6][6];
        get_pk2cc(F, S, Dm);

        // Calculate CC from Dm
        double CC[3][3][3][3];
        mat_models_carray::voigt_to_cc_carray<3>(Dm, CC);

        // Compare CC with reference solution
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                for (int k = 0; k < 3; k++){
                    for (int l = 0; l < 3; l++){
                        EXPECT_NEAR(CC[i][j][k][l], CC_ref[i][j][k][l], 
                        fmax(abs_tol, rel_tol * fabs(CC_ref[i][j][k][l])));   
                    }
                }
            }
        }

        // Print results if verbose
        if (verbose) {
            printMaterialParameters();

            std::cout << "F =" << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    std::cout << F[i][j] << " ";
                }
                std::cout << std::endl;
            }

            std::cout << "CC =" << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        for (int l = 0; l < 3; l++) {
                            std::cout << CC[i][j][k][l] << " ";
                        }
                        std::cout << std::endl;
                    }
                }
                std::cout << std::endl;
            }

            std::cout << "CC_ref =" << std::endl;
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    for (int k = 0; k < 3; k++) {
                        for (int l = 0; l < 3; l++) {
                            std::cout << CC_ref[i][j][k][l] << " ";
                        }
                        std::cout << std::endl;
                    }
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
    }
};

// --------------------------------------------------------------
// Class for test of Neo-Hookean material model
class TestNeoHookean : public TestMaterialModel {
public:

    NeoHookeanParams params;

    TestNeoHookean(const NeoHookeanParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_nHook, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {

        // Set Neo-Hookean material parameters for svFSIplus
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.C10 = params.C10;
        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter
    }

    // Print Neo-Hookean material parameters
    void printMaterialParameters() {
        std::cout << "C10 = " << params.C10 << std::endl;
    }

    // Compute strain energy for Neo-Hookean material model
    double computeStrainEnergy(const double F[3][3]) {

        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for Neo-Hookean material model
        // Psi_iso = C10 * (Ib1 - 3)
        double Psi_iso = params.C10 * (smTerms.Ib1 - 3.);

        return Psi_iso;
    }
};

// --------------------------------------------------------------
// Class for test of Mooney-Rivlin material model
class TestMooneyRivlin : public TestMaterialModel {
public:

    MooneyRivlinParams params;

    TestMooneyRivlin(const MooneyRivlinParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_MR, consts::ConstitutiveModelType::stVol_ST91),
        params(params_) 
        {

        // Set Mooney-Rivlin material parameters for svFSIplus
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.C01 = params.C01;
        dmn.stM.C10 = params.C10;
        dmn.stM.Kpen = 0.0;         // Zero volumetric penalty parameter
    }

    // Print Mooney-Rivlin material parameters
    void printMaterialParameters() {
        std::cout << "C01 = " << params.C01 << ", C10 = " << params.C10 << std::endl;
    }

    // Compute strain energy for Mooney-Rivlin material model
    double computeStrainEnergy(const double F[3][3]) {

        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for Mooney-Rivlin material model
        // Psi_iso = C10 * (Ib1 - 3) + C01 * (Ib2 - 3)
        double Psi_iso = params.C10 * (smTerms.Ib1 - 3.) + params.C01 * (smTerms.Ib2 - 3.);

        return Psi_iso;
    }
};

// ----------------------------------------------------------------------------
// Class for test of Holzapfel-Ogden material model
class TestHolzapfelOgden : public TestMaterialModel {
public:

    HolzapfelOgdenParams params;

    TestHolzapfelOgden(const HolzapfelOgdenParams &params_) : TestMaterialModel( consts::ConstitutiveModelType::stIso_HO, consts::ConstitutiveModelType::stVol_ST91),
        params(params_)
        {

        // Set Holzapfel-Ogden material parameters
        auto &dmn = com_mod.mockEq.mockDmn;
        dmn.stM.a = params.a;
        dmn.stM.b = params.b;
        dmn.stM.aff = params.a_f;
        dmn.stM.bff = params.b_f;
        dmn.stM.ass = params.a_s;
        dmn.stM.bss = params.b_s;
        dmn.stM.afs = params.a_fs;
        dmn.stM.bfs = params.b_fs;
        dmn.stM.Kpen = params.kappa;

        // Set number of fiber directions and fiber directions
        nFn = 2;
        Vector<double> f = {params.f[0], params.f[1], params.f[2]};
        Vector<double> s = {params.s[0], params.s[1], params.s[2]};
        fN.set_col(0, f);
        fN.set_col(1, s);
    }

    // Print Holzapfel-Ogden material parameters
    void printMaterialParameters() {
        std::cout << "a = " << params.a << std::endl;
        std::cout << "b = " << params.b << std::endl;
        std::cout << "a_f = " << params.a_f << std::endl;
        std::cout << "b_f = " << params.b_f << std::endl;
        std::cout << "a_s = " << params.a_s << std::endl;
        std::cout << "b_s = " << params.b_s << std::endl;
        std::cout << "a_fs = " << params.a_fs << std::endl;
        std::cout << "b_fs = " << params.b_fs << std::endl;
        std::cout << "kappa = " << params.kappa << std::endl;
        std::cout << "k = " << params.k << std::endl;
        std::cout << "f = " << "[" << params.f[0] << " " << params.f[1] << " " << params.f[2] << "]" << std::endl;
        std::cout << "s = " << "[" << params.s[0] << " " << params.s[1] << " " << params.s[2] << "]" << std::endl;
        std::cout << "full_anisotropic_invariants = " << params.full_anisotropic_invariants << std::endl;
    }

    // Smooth Heaviside function centered at 1
    double chi(const double x, const double k=100) const {
        return 1. / (1. + exp(-k * (x - 1.)));
    }

    // Compute strain energy for Holzapfel-Ogden material model
    double computeStrainEnergy(const double F[3][3]) {

        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for Holzapfel-Ogden material model
        // Psi = a/2b * exp{b(I1_bar - 3)} + Sum_{i=f,s} [a_i/2b_i * chi(I4_i) * (exp{b_i(I4_i - 1)^2} - 1)]
        //       + a_fs/2b_fs * (exp{b_fs*I8_fs^2} - 1) + kappa/4 * (J^2 - 1 - 2*ln(J))

        // Material parameters
        double a = params.a;
        double b = params.b;
        double a_f = params.a_f;
        double b_f = params.b_f;
        double a_s = params.a_s;
        double b_s = params.b_s;
        double a_fs = params.a_fs;
        double b_fs = params.b_fs;
        double kappa = params.kappa;

        // Smoothed Heaviside parameter
        double k = params.k;

        // Fiber and sheet directions
        double f[3] = {params.f[0], params.f[1], params.f[2]};
        double s[3] = {params.s[0], params.s[1], params.s[2]};

        // Formulation used by cardiac mechanics benchmark paper (Arostica et al., 2024)
        // Uses I1_bar (bar = isochoric), but I4_f, I4_s, I8_fs (not bar)
        // This corresponds to the HO-ma (modified anisotropy) implementation in svFSIplus
        if (params.full_anisotropic_invariants) {
            // Invariants
            double I1_bar = smTerms.Ib1;
            // I4_f = f . C . f
            double C_f[3]; mat_fun_carray::mat_mul<3>(smTerms.C, f, C_f);
            double I4_f = mat_fun_carray::norm<3>(f, C_f);
            // I4_s = s . C . s
            double C_s[3]; mat_fun_carray::mat_mul<3>(smTerms.C, s, C_s);
            double I4_s = mat_fun_carray::norm<3>(s, C_s);
            // I8_fs = f . C . s
            double I8_fs = mat_fun_carray::norm<3>(f, C_s);

            // Strain energy density for Holzapfel-Ogden material model with full anisotropic invariants
            double Psi = 0.0;
            Psi += a / (2.0 * b) * exp(b * (I1_bar - 3.0));                             // Isotropic term
            Psi += a_f / (2.0 * b_f) * chi(I4_f, k) * (exp(b_f * pow(I4_f - 1.0, 2)) - 1.0);   // Fiber term
            Psi += a_s / (2.0 * b_s) * chi(I4_s, k) * (exp(b_s * pow(I4_s - 1.0, 2)) - 1.0);   // Sheet term
            Psi += a_fs / (2.0 * b_fs) * (exp(b_fs * pow(I8_fs, 2)) - 1.0);                   // Cross-fiber term
            Psi += kappa / 4.0 * (pow(smTerms.J, 2) - 1.0 - 2.0 * log(smTerms.J));      // Simo-Taylor 91 volumetric penalty term

            return Psi;
        }
        // Formulation with fully decoupled isochoric-volumetric split
        // Uses I1_bar, I4_bar_f, I4_bar_s, I8_bar_fs (bar = isochoric)
        // This corresponds to the HO-d (decoupled) implementation in svFSIplus
        else {
            // Invariants
            double I1_bar = smTerms.Ib1;
            // I4_bar_f = f . C_bar . f
            double C_bar_f[3]; mat_fun_carray::mat_mul<3>(smTerms.C_bar, f, C_bar_f);
            double I4_bar_f = mat_fun_carray::norm<3>(f, C_bar_f);
            // I4_bar_s = s . C_bar . s
            double C_bar_s[3]; mat_fun_carray::mat_mul<3>(smTerms.C_bar, s, C_bar_s);
            double I4_bar_s = mat_fun_carray::norm<3>(s, C_bar_s);
            // I8_bar_fs = f . C_bar . s
            double I8_bar_fs = mat_fun_carray::norm<3>(f, C_bar_s);

            // Strain energy density for Holzapfel-Ogden material model with modified anisotropic invariants (bar quantities)
            double Psi = 0.0;
            Psi += a / (2.0 * b) * exp(b * (I1_bar - 3.0));                             // Isotropic term
            Psi += a_f / (2.0 * b_f) * chi(I4_bar_f, k) * (exp(b_f * pow(I4_bar_f - 1.0, 2)) - 1.0);   // Fiber term
            Psi += a_s / (2.0 * b_s) * chi(I4_bar_s, k) * (exp(b_s * pow(I4_bar_s - 1.0, 2)) - 1.0);   // Sheet term
            Psi += a_fs / (2.0 * b_fs) * (exp(b_fs * pow(I8_bar_fs, 2)) - 1.0);                   // Cross-fiber term
            Psi += kappa / 4.0 * (pow(smTerms.J, 2) - 1.0 - 2.0 * log(smTerms.J));      // Volumetric penalty term

            return Psi;
        }
    }
};
