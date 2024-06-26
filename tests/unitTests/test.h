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

#include <stdlib.h>
#include <iostream>
#include "gtest/gtest.h"   // include GoogleTest
#include "mat_fun.h"
#include "mat_fun_carray.h"
#include "mat_models.h"
#include "mat_models_carray.h"


template<int N>
void calc_JCE(double F[N][N], double &J, double C[N][N], double E[N][N]) {
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

// Class for unit test of isotropic material model (currently only introduced two 
// parameters C10 and C01)
class UnitTestIso {
public:
    MockComMod com_mod;
    MockCepMod cep_mod;
    
    UnitTestIso(consts::ConstitutiveModelType matType, double E, double nu, 
                consts::ConstitutiveModelType penType, double pen, double C01 = 0.0) {
        int nsd = com_mod.nsd;
        auto &dmn = com_mod.mockEq.mockDmn;
        mat_fun_carray::ten_init(nsd);                        // initialize tensor index pointer
        dmn.stM.isoType = matType;                            // Mat_model
        double mu  = 0.5 * E / (1.0 + nu);                    // Shear_modulus
        dmn.stM.C10 = 0.5 * mu - C01;                         // set_material_props.h 
        dmn.stM.C01 = C01;
        dmn.stM.volType = penType;                            // Dilational_penalty_model
        dmn.stM.Kpen = pen;                                   // Penalty_parameter
    }

    void compare_S_Dm(double F[3][3], double S_ref[3][3], double Dm_ref[6][6], double rel_tol) {
        int nsd = com_mod.nsd;
        auto &dmn = com_mod.mockEq.mockDmn;
        // hard code for nHK
        int nFn = 1; 
        Array<double> fN(nsd, nFn);
        double ya_g = 0.0;   
        double S[3][3], Dm[6][6];

        mat_models_carray::get_pk2cc<3>(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);

        // Compare S with reference solution
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                EXPECT_NEAR(S[i][j], S_ref[i][j], rel_tol * fabs(S_ref[i][j]));   
            }
        }
    }

    // Function to get S and Dm from F
    void get_pk2cc(double F[3][3], double S[3][3], double Dm[6][6]) {
        int nsd = com_mod.nsd;
        auto &dmn = com_mod.mockEq.mockDmn;
        // hard code for nHK
        int nFn = 1; 
        Array<double> fN(nsd, nFn);
        double ya_g = 0.0;   

        mat_models_carray::get_pk2cc<3>(com_mod, cep_mod, dmn, F, nFn, fN, ya_g, S, Dm);
    }
};