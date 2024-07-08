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

// Perturb the deformation gradient F by delta times a random number between -1 and 1
// and store the perturbed deformation gradient in F_tilde
template<int N>
void perturb_random_F(const double F[N][N], const double delta, double F_tilde[N][N]) {

    // Perturb the deformation gradient and store in F_tilde
    double dF_iJ;
    for (int i = 0; i < N; i++) {
        for (int J = 0; J < N; J++) {
            dF_iJ = delta * (2.0 * std::rand() / RAND_MAX - 1.0); // (random number between -1 and 1) * delta
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

// Class to contain material parameters
class MatParams {
public:
    virtual ~MatParams() {} // Virtual destructor for proper cleanup
};

// Class to contain St. Venant-Kirchhoff material parameters
class StVKParams : public MatParams {
public:
    double C01;
    double C10;

    StVKParams(double c01, double c10) : C01(c01), C10(c10) {}
};

// Class to contain Neo-Hookean material parameters
class NeoHookeanParams : public MatParams {
public:
    double C10;

    NeoHookeanParams(double c10) : C10(c10) {}
};

// Class to contain Mooney-Rivlin material parameters
class MooneyRivlinParams : public MatParams {
public:
    double C01;
    double C10;

    MooneyRivlinParams(double c01, double c10) : C01(c01), C10(c10) {}
};





// Class to represent a strain energy density function
class StrainEnergy {
public:
    virtual ~StrainEnergy() {}
    virtual double compute(const double F[3][3], const MatParams& params) const = 0;
};

// Class to represent a St. Venant-Kirchhoff strain energy density function
class StVKStrainEnergy : public StrainEnergy {
public:
    double compute(const double F[3][3], const MatParams& params) const override {
        const StVKParams& stVKParams = dynamic_cast<const StVKParams&>(params);

        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for St. Venant-Kirchhoff material model
        // Psi_iso = C10/2 * tr(E)^2 + C01 * tr(E^2)
        double trE = mat_fun_carray::mat_trace<3>(smTerms.E);
        double trE2 = mat_fun_carray::mat_trace<3>(smTerms.E2);
        double Psi_iso = stVKParams.C10/2. * pow(trE, 2) + stVKParams.C01 * trE2;

        return Psi_iso;
    }
};

// Class to represent a Neo-Hookean strain energy density function
class NeoHookeanStrainEnergy : public StrainEnergy {
public:
    double compute(const double F[3][3], const MatParams& params) const override {
        const NeoHookeanParams& nhParams = dynamic_cast<const NeoHookeanParams&>(params);

        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for Neo-Hookean material model
        // Psi_iso = C10 * (Ib1 - 3)
        double Psi_iso = nhParams.C10 * (smTerms.Ib1 - 3.);

        return Psi_iso;
    }
};


// Class to represent a Mooney-Rivlin strain energy density function
class MooneyRivlinStrainEnergy : public StrainEnergy {
public:
    double compute(const double F[3][3], const MatParams& params) const override {
        const MooneyRivlinParams& mrParams = dynamic_cast<const MooneyRivlinParams&>(params);

        // Compute solid mechanics terms
        solidMechanicsTerms smTerms = calcSolidMechanicsTerms(F);

        // Strain energy density for Mooney-Rivlin material model
        // Psi_iso = C10 * (Ib1 - 3) + C01 * (Ib2 - 3)
        double Psi_iso = mrParams.C10 * (smTerms.Ib1 - 3.) + mrParams.C01 * (smTerms.Ib2 - 3.);

        return Psi_iso;

    }
};


// Function to compute the strain energy density Psi for the
// St. Venant-Kirchhoff material model given the material parameters and deformation gradient F
template<int N>
double Psi_StVK(const double C01, const double C10, const double F[N][N]) {

    // Compute solid mechanics terms
    solidMechanicsTerms<N> smTerms = calcSolidMechanicsTerms<N>(F);

    // Strain energy density for St. Venant-Kirchhoff material model
    // Psi_iso = C10/2 * tr(E)^2 + C01 * tr(E^2)
    std::cout << "smTerms.E = " << std::endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << smTerms.E[i][j] << " ";
        }
        std::cout << std::endl;
    }

    std::cout << "smTerms.E2 = " << std::endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            std::cout << smTerms.E2[i][j] << " ";
        }
        std::cout << std::endl;
    }

    double trE = mat_fun_carray::mat_trace<N>(smTerms.E);
    double trE2 = mat_fun_carray::mat_trace<N>(smTerms.E2);

    std::cout << "trE = " << trE << std::endl;
    std::cout << "trE2 = " << trE2 << std::endl;

    std::cout << "C10 = " << C10 << std::endl;
    std::cout << "C01 = " << C01 << std::endl;
    
    double Psi_iso = C10/2. * pow(trE, 2) + C01 * trE2;

}

// Function to compute the (isochoric) strain energy density Psi for the
// Neo-Hookean material model given the material parameters and deformation gradient F
template<int N>
double Psi_nHook(const double C10, const double F[N][N]) {

    // Compute solid mechanics terms
    solidMechanicsTerms<N> smTerms = calcSolidMechanicsTerms<N>(F);

    // Strain energy density for Neo-Hookean material model
    // Psi_iso = C10 * (Ib1 - 3)
    double Psi_iso = C10 * (smTerms.Ib1 - 3.);

    return Psi_iso;

}

// Function to compute the (isochoric) strain energy density Psi for the 
// Mooney-Rivlin material model given the material parameters and deformation gradient F
template<int N>
double Psi_MR(const double C01, const double C10, const double F[N][N]) {          

    // Compute solid mechanics terms
    solidMechanicsTerms<N> smTerms = calcSolidMechanicsTerms<N>(F);

    // Strain energy density for Mooney-Rivlin material model
    // Psi_iso = C10 * (Ib1 - 3) + C01 * (Ib2 - 3)
    double Psi_iso = C10 * (smTerms.Ib1 - 3.) + C01 * (smTerms.Ib2 - 3.);

    return Psi_iso;

}

// Assuming Psi is a template function that takes the size, mat_params, and a 3x3 matrix as arguments
template<int N>
void calcPK2StressFiniteDifference(const double F[N][N], const StrainEnergy &PsiFunc, const MatParams &matParams, double delta, double (&S_ref)[N][N]) {
    // Compute strain energy density given F
    double Psi = PsiFunc.compute(F, matParams);

    // Compute 1st PK stress P_ref_iJ = dPsi / dF[i][J] using finite difference, component by component
    double P_ref[3][3] = {};
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
            double Psi_tilde = PsiFunc.compute(F_tilde, matParams);

            // Compute differences in Psi
            double dPsi = Psi_tilde - Psi;

            // Compute P_ref[i][J] = dPsi / dF[i][J]
            P_ref[i][J] = dPsi / delta;
        }
    }

    // Compute S_ref = F^-1 * P_ref
    double F_inv[N][N];
    mat_fun_carray::mat_inv<N>(F, F_inv);
    mat_fun_carray::mat_mul<N>(F_inv, P_ref, S_ref);
}


// Test fixture class for St. Venant-Kirchhoff material model
class StVenantKirchhoffTest : public ::testing::Test {
protected:
    // Variables common across tests
    consts::ConstitutiveModelType matType; // Material model type
    consts::ConstitutiveModelType volType; // Dilational penalty model type
    double E_mod; // Elasticity modulus
    double nu; // Poisson ratio
    double pen; // Volumetric penalty parameter
    double C01; // St. Venant Kirckhoff parameter
    double mu;  // Shear modulus
    double C10; // St. Venant Kirckhoff parameter

    // Add the UnitTestIso object as a member variable
    UnitTestIso* StVK;

    // Setup method to initialize variables before each test
    void SetUp() override {
        matType = consts::ConstitutiveModelType::stIso_StVK;
        volType = consts::ConstitutiveModelType::stVol_ST91;
        E_mod = 1e6;
        nu = 0.495;
        pen = 0.0;
        C01 = 0.1;

        // Compute derived material parameters
        mu = 0.5 * E_mod / (1.0 + nu);
        C10 = 0.5 * mu - C01;

        // Initialize the UnitTestIso object
        StVK = new UnitTestIso(matType, E_mod, nu, volType, pen, C01);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the UnitTestIso object
        delete StVK;
        StVK = nullptr;
    }
};

// Test fixture class for Neo-Hookean material model
class NeoHookeanTest : public ::testing::Test {
protected:
    // Variables common across tests
    consts::ConstitutiveModelType matType; // Material model type
    consts::ConstitutiveModelType volType; // Dilational penalty model type
    double E_mod; // Elasticity modulus
    double nu;    // Poisson ratio
    double C01;   // NeoHookean parameter
    double pen;   // Volumetric penalty parameter
    double C10;   // NeoHookean parameter

    // Add the UnitTestIso object as a member variable
    UnitTestIso* nHook;

    // Setup method to initialize variables before each test
    void SetUp() override {
        matType = consts::ConstitutiveModelType::stIso_nHook;
        volType = consts::ConstitutiveModelType::stVol_ST91;
        E_mod = 1e6;
        nu = 0.495;
        pen = 0.0;
        C01 = 0.1;

        // Compute derived material parameters
        double mu = 0.5 * E_mod / (1.0 + nu);
        C10 = 0.5 * mu - C01;

        // Initialize the UnitTestIso object
        nHook = new UnitTestIso(matType, E_mod, nu, volType, pen, C01);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the UnitTestIso object
        delete nHook;
        nHook = nullptr;
    }
};

// Test fixture class for Mooney-Rivlin material model
class MooneyRivlinTest : public ::testing::Test {
protected:
    // Variables common across tests
    consts::ConstitutiveModelType matType; // Material model type
    consts::ConstitutiveModelType volType; // Dilational penalty model type
    double E_mod; // Elasticity modulus
    double nu;   // Poisson ratio
    double pen; // Volumetric penalty parameter
    double C01; // Mooney-Rivlin parameter
    double C10; // Mooney-Rivlin parameter

    // Add the UnitTestIso object as a member variable
    UnitTestIso* MR;

    // Setup method to initialize variables before each test
    void SetUp() override {
        matType = consts::ConstitutiveModelType::stIso_MR;
        volType = consts::ConstitutiveModelType::stVol_ST91;
        E_mod = 1e6;
        nu = 0.495;
        pen = 0.0;
        C01 = 0.1;

        // Compute derived material parameters
        double mu = 0.5 * E_mod / (1.0 + nu);
        C10 = 0.5 * mu - C01;

        // Initialize the UnitTestIso object
        MR = new UnitTestIso(matType, E_mod, nu, volType, pen, C01);
    }

    // TearDown method to clean up after each test, if needed
    void TearDown() override {
        // Clean up the UnitTestIso object
        delete MR;
        MR = nullptr;
    }
};


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

// Test the consistency of the PK2 stress tensor S(F) from get_pk2cc() with the strain 
// energy density Psi(F) provided by the user
// Psuedocode:
// - Generate random deformation gradient F
// - Compute Psi(F)
// - Compute S(F) from get_pk2cc()
// - For many random dF
//      - Compute dPsi = Psi(F + dF) - Psi(F)
//      - Compute dE from dF
//      - Check that S:dE = dPsi
TEST_F(StVenantKirchhoffTest, TestPK2Stress) {
    // Create random deformation gradient
    double F[3][3];
    create_random_F(F);

    // Compute E from F
    double J, C[3][3], E[3][3];
    calc_JCE(F, J, C, E);

    // Compute strain energy density Psi
    double Psi = Psi_StVK<3>(C01, C10, F);

    // Compute S(F)
    double S[3][3], Dm[6][6];
    StVK->get_pk2cc(F, S, Dm); // S from svFSI


    // Generate many random dF and check that S:dE = dPsi
    // S was obtained from get_pk2cc(), and dPsi = Psi(F + dF) - Psi(F)
    int n_iter = 10;
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double delta = 1e-6; // perturbation scaling factor

    double F_tilde[3][3]; // perturbed deformation gradient
    
    // Loop over many random perturbations to the deformation gradient
    for (int i = 0; i < n_iter; i++) {
        // Perturb the deformation gradient
        perturb_random_F<3>(F, delta, F_tilde);

        // Compute perturbed E
        double J_tilde, C_tilde[3][3], E_tilde[3][3];
        calc_JCE(F_tilde, J_tilde, C_tilde, E_tilde);

        // Compute dE
        double dE[3][3];
        for (int i = 0; i < 3; i++) {
            for (int J = 0; J < 3; J++) {
                dE[i][J] = E_tilde[i][J] - E[i][J];
            }
        }

        // Compute perturbed Psi with perturbed deformation gradient
        double Psi_tilde = Psi_StVK<3>(C01, C10, F_tilde);

        // Compute dPsi
        double dPsi = Psi_tilde - Psi;

        // Check that S_ij dE_ij = dPsi
        double SdE = mat_fun_carray::mat_ddot<3>(S, dE);
        EXPECT_NEAR(SdE, dPsi, rel_tol * fabs(dPsi));
        //std::cout << "SdE = " << SdE << ", dPsi = " << dPsi << std::endl;

    }
}

// Test the consistency of the material elasticity tensor CC(F) from get_pk2cc() with the
// PK2 stress tensor S(F) from get_pk2cc()
// Psuedocode:
// - Generate random deformation gradient F
// - Compute S(F) and CC(F) from get_pk2cc()
// - For many random dF
//      - Compute S(F + dF) from get_pk2cc()
//      - Compute dS = S(F + dF) - S(F)
//      - Compute dE from dF
//      - Check that CC:dE = dS
TEST_F(StVenantKirchhoffTest, TestMaterialElasticityTensor) {
    // Create random deformation gradient
    double F[3][3];
    create_random_F(F);

    // Compute E from F
    double J, C[3][3], E[3][3];
    calc_JCE(F, J, C, E);

    // Compute S_ij(F)
    // Compute CC_ijkl(F). 
    // CC is provided in Voigt notation as Dm, and we will convert it to CC
    double S[3][3], Dm[6][6];
    StVK->get_pk2cc(F, S, Dm); // S from svFSI

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
            EXPECT_NEAR(Dm_check[i][j], Dm[i][j], 1e-12 * fabs(Dm[i][j]));
        }
    }
    // -------------------------------

    // Generate many random dF and check that CC:dE = dS
    // CC was obtained from get_pk2cc(), and dS = S(F + dF) - S(F), 
    // where S is also obtained from get_pk2cc()
    int n_iter = 10;
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double delta = 1e-6; // perturbation scaling factor
    double F_tilde[3][3]; // perturbed deformation gradient
    
    // Loop over many random perturbations to the deformation gradient
    for (int i = 0; i < n_iter; i++) {
        // Perturb the deformation gradient
        perturb_random_F<3>(F, delta, F_tilde);

        // Compute perturbed E
        double J_tilde, C_tilde[3][3], E_tilde[3][3];
        calc_JCE(F_tilde, J_tilde, C_tilde, E_tilde);

        // Compute dE
        double dE[3][3];
        for (int i = 0; i < 3; i++) {
            for (int J = 0; J < 3; J++) {
                dE[i][J] = E_tilde[i][J] - E[i][J];
            }
        }

        // Compute perturbed S with perturbed deformation gradient
        double S_tilde[3][3], Dm_tilde[6][6];
        StVK->get_pk2cc(F_tilde, S_tilde, Dm_tilde);

        // Compute dS
        double dS[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dS[i][j] = S_tilde[i][j] - S[i][j];
            }
        }

        // Check that CC_ijkl dE_kl = dS_ij
        double CCdE[3][3];
        mat_fun_carray::ten_mat_ddot<3>(CC, dE, CCdE);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                EXPECT_NEAR(CCdE[i][j], dS[i][j], rel_tol * fabs(dS[i][j]));
                //std::cout << "CCdE[" << i << "][" << j << "] = " << CCdE[i][j] << ", dS[" << i << "][" << j << "] = " << dS[i][j] << std::endl;
            }
        }
    }
}

// Test the consistency of the PK2 stress tensor S(F) from get_pk2cc() with the strain 
// energy density Psi(F) provided by the user
// Psuedocode:
// - Generate random deformation gradient F
// - Compute Psi(F)
// - Compute S(F) from get_pk2cc()
// - For many random dF
//      - Compute dPsi = Psi(F + dF) - Psi(F)
//      - Compute dE from dF
//      - Check that S:dE = dPsi
TEST_F(NeoHookeanTest, TestPK2Stress) {
    // Create random deformation gradient
    double F[3][3];
    create_random_F(F);

    // Compute E from F
    double J, C[3][3], E[3][3];
    calc_JCE(F, J, C, E);

    // Compute Psi(F)
    double Psi = Psi_nHook<3>(C10, F);

    // Compute S(F)
    double S[3][3], Dm[6][6];
    nHook->get_pk2cc(F, S, Dm); // S from svFSI

    // Generate many random dF and check that S:dE = dPsi
    // S was obtained from get_pk2cc(), and dPsi = Psi(F + dF) - Psi(F)
    int n_iter = 10;
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double delta = 1e-6; // perturbation scaling factor

    double F_tilde[3][3]; // perturbed deformation gradient
    
    // Loop over many random perturbations to the deformation gradient
    for (int i = 0; i < n_iter; i++) {
        // Perturb the deformation gradient
        perturb_random_F<3>(F, delta, F_tilde);

        // Compute perturbed E
        double J_tilde, C_tilde[3][3], E_tilde[3][3];
        calc_JCE(F_tilde, J_tilde, C_tilde, E_tilde);

        // Compute dE
        double dE[3][3];
        for (int i = 0; i < 3; i++) {
            for (int J = 0; J < 3; J++) {
                dE[i][J] = E_tilde[i][J] - E[i][J];
            }
        }

        // Compute perturbed Psi with perturbed deformation gradient
        double Psi_tilde = Psi_nHook<3>(C10, F_tilde);

        // Compute dPsi
        double dPsi = Psi_tilde - Psi;

        // Check that S_ij dE_ij = dPsi
        double SdE = mat_fun_carray::mat_ddot<3>(S, dE);
        EXPECT_NEAR(SdE, dPsi, rel_tol * fabs(dPsi));
        //std::cout << "SdE = " << SdE << ", dPsi = " << dPsi << std::endl;

    }
}

// Test the consistency of the material elasticity tensor CC(F) from get_pk2cc() with the
// PK2 stress tensor S(F) from get_pk2cc()
// Psuedocode:
// - Generate random deformation gradient F
// - Compute S(F) and CC(F) from get_pk2cc()
// - For many random dF
//      - Compute S(F + dF) from get_pk2cc()
//      - Compute dS = S(F + dF) - S(F)
//      - Compute dE from dF
//      - Check that CC:dE = dS
TEST_F(NeoHookeanTest, TestMaterialElasticityTensor) {
    // Create random deformation gradient
    double F[3][3];
    create_random_F(F);

    // Compute E from F
    double J, C[3][3], E[3][3];
    calc_JCE(F, J, C, E);

    // Compute S_ij(F)
    // Compute CC_ijkl(F). 
    // CC is provided in Voigt notation as Dm, and we will convert it to CC
    double S[3][3], Dm[6][6];
    nHook->get_pk2cc(F, S, Dm); // S from svFSI

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
            EXPECT_NEAR(Dm_check[i][j], Dm[i][j], 1e-12 * fabs(Dm[i][j]));
        }
    }
    // -------------------------------

    // Generate many random dF and check that CC:dE = dS
    // CC was obtained from get_pk2cc(), and dS = S(F + dF) - S(F), 
    // where S is also obtained from get_pk2cc()
    int n_iter = 10;
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double delta = 1e-6; // perturbation scaling factor
    double F_tilde[3][3]; // perturbed deformation gradient
    
    // Loop over many random perturbations to the deformation gradient
    for (int i = 0; i < n_iter; i++) {
        // Perturb the deformation gradient
        perturb_random_F<3>(F, delta, F_tilde);

        // Compute perturbed E
        double J_tilde, C_tilde[3][3], E_tilde[3][3];
        calc_JCE(F_tilde, J_tilde, C_tilde, E_tilde);

        // Compute dE
        double dE[3][3];
        for (int i = 0; i < 3; i++) {
            for (int J = 0; J < 3; J++) {
                dE[i][J] = E_tilde[i][J] - E[i][J];
            }
        }

        // Compute perturbed S with perturbed deformation gradient
        double S_tilde[3][3], Dm_tilde[6][6];
        nHook->get_pk2cc(F_tilde, S_tilde, Dm_tilde);

        // Compute dS
        double dS[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dS[i][j] = S_tilde[i][j] - S[i][j];
            }
        }

        // Check that CC_ijkl dE_kl = dS_ij
        double CCdE[3][3];
        mat_fun_carray::ten_mat_ddot<3>(CC, dE, CCdE);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                EXPECT_NEAR(CCdE[i][j], dS[i][j], rel_tol * fabs(dS[i][j]));
                //std::cout << "CCdE[" << i << "][" << j << "] = " << CCdE[i][j] << ", dS[" << i << "][" << j << "] = " << dS[i][j] << std::endl;
            }
        }
    }
}

// Calculate S(F) using finite differences and compare with the S(F) from svFSI
// Psuedocode:
// - Generate random deformation gradient F
// - Compute P_ij(F) = dPsi / dF_ij using finite differences
//      - Perturb the iJ-th component of F by small delta
//      - dPsi = Psi(F + dF) - Psi(F)
//      - P_ij = dPsi / dF_ij
// - Compute S(F) = F^-1 * P(F)
TEST_F(MooneyRivlinTest, TestPK2StressDirect) {
    // Create random deformation gradient
    double F[3][3];
    create_random_F(F);

    // Compute the reference S_ij using finite difference
    MooneyRivlinStrainEnergy PsiFunc;
    MooneyRivlinParams matParams = {C01, C10};
    double delta = 1e-9; // perturbation size
    double S_ref[3][3];
    calcPK2StressFiniteDifference<3>(F, PsiFunc, matParams, delta, S_ref);

    // TODO: Compute Dm_ref
    double Dm_ref[6][6] = {};

    // Compare S and Dm from svFSI with reference solutions
    double rel_tol = 1e-4;
    MR->compare_S_Dm(F, S_ref, Dm_ref, rel_tol);
}

// Test the consistency of the PK2 stress tensor S(F) from get_pk2cc() with the strain 
// energy density Psi(F) provided by the user
// Psuedocode:
// - Generate random deformation gradient F
// - Compute Psi(F)
// - Compute S(F) from get_pk2cc()
// - For many random dF
//      - Compute dPsi = Psi(F + dF) - Psi(F)
//      - Compute dE from dF
//      - Check that S:dE = dPsi
TEST_F(MooneyRivlinTest, TestPK2Stress) {
    // Create random deformation gradient
    double F[3][3];
    create_random_F(F);

    // Compute E from F
    double J, C[3][3], E[3][3];
    calc_JCE(F, J, C, E);

    // Compute Psi(F)
    double Psi = Psi_MR<3>(C01, C10, F);

    // Compute S(F)
    double S[3][3], Dm[6][6];
    MR->get_pk2cc(F, S, Dm); // S from svFSI

    // Generate many random dF and check that S:dE = dPsi
    // S was obtained from get_pk2cc(), and dPsi = Psi(F + dF) - Psi(F)
    int n_iter = 10;
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double delta = 1e-6; // perturbation scaling factor

    double F_tilde[3][3]; // perturbed deformation gradient
    
    // Loop over many random perturbations to the deformation gradient
    for (int i = 0; i < n_iter; i++) {
        // Perturb the deformation gradient
        perturb_random_F<3>(F, delta, F_tilde);

        // Compute perturbed E
        double J_tilde, C_tilde[3][3], E_tilde[3][3];
        calc_JCE(F_tilde, J_tilde, C_tilde, E_tilde);

        // Compute dE
        double dE[3][3];
        for (int i = 0; i < 3; i++) {
            for (int J = 0; J < 3; J++) {
                dE[i][J] = E_tilde[i][J] - E[i][J];
            }
        }

        // Compute perturbed Psi with perturbed deformation gradient
        double Psi_tilde = Psi_MR<3>(C01, C10, F_tilde);

        // Compute dPsi
        double dPsi = Psi_tilde - Psi;

        // Check that S_ij dE_ij = dPsi
        double SdE = mat_fun_carray::mat_ddot<3>(S, dE);
        EXPECT_NEAR(SdE, dPsi, rel_tol * fabs(dPsi));
        //std::cout << "SdE = " << SdE << ", dPsi = " << dPsi << std::endl;

    }
}

// Test the consistency of the material elasticity tensor CC(F) from get_pk2cc() with the
// PK2 stress tensor S(F) from get_pk2cc()
// Psuedocode:
// - Generate random deformation gradient F
// - Compute S(F) and CC(F) from get_pk2cc()
// - For many random dF
//      - Compute S(F + dF) from get_pk2cc()
//      - Compute dS = S(F + dF) - S(F)
//      - Compute dE from dF
//      - Check that CC:dE = dS
TEST_F(MooneyRivlinTest, TestMaterialElasticityTensor) {
    // Create random deformation gradient
    double F[3][3];
    create_random_F(F);

    // Compute E from F
    double J, C[3][3], E[3][3];
    calc_JCE(F, J, C, E);

    // Compute S_ij(F)
    // Compute CC_ijkl(F). 
    // CC is provided in Voigt notation as Dm, and we will convert it to CC
    double S[3][3], Dm[6][6];
    MR->get_pk2cc(F, S, Dm); // S from svFSI

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
            EXPECT_NEAR(Dm_check[i][j], Dm[i][j], 1e-12 * fabs(Dm[i][j]));
        }
    }
    // -------------------------------

    // Generate many random dF and check that CC:dE = dS
    // CC was obtained from get_pk2cc(), and dS = S(F + dF) - S(F), 
    // where S is also obtained from get_pk2cc()
    int n_iter = 10;
    double rel_tol = 1e-3; // relative tolerance for comparing dPsi and dS with values from svFSI
    double delta = 1e-6; // perturbation scaling factor
    double F_tilde[3][3]; // perturbed deformation gradient
    
    // Loop over many random perturbations to the deformation gradient
    for (int i = 0; i < n_iter; i++) {
        // Perturb the deformation gradient
        perturb_random_F<3>(F, delta, F_tilde);

        // Compute perturbed E
        double J_tilde, C_tilde[3][3], E_tilde[3][3];
        calc_JCE(F_tilde, J_tilde, C_tilde, E_tilde);

        // Compute dE
        double dE[3][3];
        for (int i = 0; i < 3; i++) {
            for (int J = 0; J < 3; J++) {
                dE[i][J] = E_tilde[i][J] - E[i][J];
            }
        }

        // Compute perturbed S with perturbed deformation gradient
        double S_tilde[3][3], Dm_tilde[6][6];
        MR->get_pk2cc(F_tilde, S_tilde, Dm_tilde);

        // Compute dS
        double dS[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                dS[i][j] = S_tilde[i][j] - S[i][j];
            }
        }

        // Check that CC_ijkl dE_kl = dS_ij
        double CCdE[3][3];
        mat_fun_carray::ten_mat_ddot<3>(CC, dE, CCdE);
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                EXPECT_NEAR(CCdE[i][j], dS[i][j], rel_tol * fabs(dS[i][j]));
                //std::cout << "CCdE[" << i << "][" << j << "] = " << CCdE[i][j] << ", dS[" << i << "][" << j << "] = " << dS[i][j] << std::endl;
            }
        }
    }
}



