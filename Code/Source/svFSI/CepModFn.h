/**
 * Copyright (c) Stanford University, The Regents of the University of California, and others.
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

#ifndef CEP_MOD_FN_H 
#define CEP_MOD_FN_H 

#include "Array.h"
#include "Vector.h"

/// @brief This module defines data structures for Fitzhugh-Nagumo cellular
/// activation model for cardiac electrophysiology.
///
/// The classes defined here duplicate the data structures in the Fortran FNMOD module defined in CEPMOD_FN.f 
/// and PARAMS_FN.f files. 
class CepModFn
{
  public:
    CepModFn();
    ~CepModFn();

    // Scaling factors
    /// Voltage scaling
    double Vscale  = 1.0;
    /// Time scaling
    double Tscale  = 1.0;
    /// Voltage offset parameter
    double Voffset = 0.0;

    /// Model parameters
    double alpha = -0.50;
    double a = 0.0;
    double b = -0.60;
    double c = 50.0;

    void getf(const int n, const Vector<double>& X, Vector<double>& f, const double fext);
    void getj(const int n, const Vector<double>& X, Array<double>& JAC);

    void init(const int nX, Vector<double> &X);
    void init(const int nX, Vector<double> &X, double X0);
    void init(const int nX, Vector<double> &X, Vector<double>& X0);

    void integ_cn2(const int nX, Vector<double>& Xn, const double Ts, const double Ti, const double Istim, 
        Vector<int>& IPAR, Vector<double>& RPAR);
    void integ_fe(const int nX, Vector<double>& X, const double Ts, const double Ti, const double Istim);
    void integ_rk(const int nX, Vector<double>& X, const double Ts, const double Ti, const double Istim);

};

#endif

