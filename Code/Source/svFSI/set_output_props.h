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

// The 'output_props_map' map defined here sets equation 
// output properties from the given OutputType.
//

#include <tuple>

/// @brief The 'OutputProps' tuple stores data for the 'outputType' object for
///
///   output.grp - The group that this belong to (one of outType_)
///   output.o - Offset from the first index
///   output.l - Length of the outputed variable 
///   output.name - The name to be used for the output and also in input file
//
using OutputProps = std::tuple<consts::OutputType, int, int, std::string>; 

/// @brief Reproduces Fortran READOUTPUTS.
//
std::map<consts::OutputType, OutputProps> output_props_map = 
{
  //                                             -----------------------------------------------------------------
  //                                                    output group        o   l                    name
  //                                             -----------------------------------------------------------------
  {OutputType::out_absVelocity,  std::make_tuple(OutputType::outGrp_absV,   0, nsd,           "Absolute_velocity") },
  {OutputType::out_acceleration, std::make_tuple(OutputType::outGrp_A,      0, nsd,           "Acceleration") },
  {OutputType::out_cauchy,       std::make_tuple(OutputType::outGrp_cauchy, 0, com_mod.nsymd, "Cauchy_stress") },

  {OutputType::out_CGInv1,       std::make_tuple(OutputType::out_CGInv1,   0,  1,             "CG_Strain_Trace") },
  {OutputType::out_CGstrain,     std::make_tuple(OutputType::outGrp_C,     0,  com_mod.nsymd, "CG_Strain") },

  {OutputType::out_defGrad,      std::make_tuple(OutputType::outGrp_F,      0, nsd*nsd,       "Def_grad") },
  {OutputType::out_displacement, std::make_tuple(OutputType::outGrp_D,      0, nsd,           "Displacement") },
  {OutputType::out_divergence,   std::make_tuple(OutputType::outGrp_divV,   0, 1,             "Divergence") },
  {OutputType::out_energyFlux,   std::make_tuple(OutputType::outGrp_eFlx,   0, nsd,           "Energy_flux") },

  {OutputType::out_fibAlign,     std::make_tuple(OutputType::outGrp_fA,     0, 1,             "Fiber_alignment") },
  {OutputType::out_fibDir,       std::make_tuple(OutputType::outGrp_fN,     0, nsd,           "Fiber_direction") },
  {OutputType::out_fibStrn,      std::make_tuple(OutputType::outGrp_fS,     0, 1,             "Fiber_shortening") },

  {OutputType::out_heatFlux,     std::make_tuple(OutputType::outGrp_hFlx,   0, nsd,           "Heat_flux") },
  {OutputType::out_integ,        std::make_tuple(OutputType::outGrp_I,      0,   1, nsd == 2 ?  "Area" : "Volume") },
  {OutputType::out_jacobian,     std::make_tuple(OutputType::outGrp_J,      0,   1,             "Jacobian") },
  {OutputType::out_mises,        std::make_tuple(OutputType::outGrp_mises,  0,   1,             "VonMises_stress") },
  {OutputType::out_pressure,     std::make_tuple(OutputType::outGrp_Y,      nsd, 1,           "Pressure") },
  {OutputType::out_strain,       std::make_tuple(OutputType::outGrp_strain, 0, com_mod.nsymd, "Strain") },
  {OutputType::out_strainInv,    std::make_tuple(OutputType::outGrp_stInv,  0, nsd,           "Strain_invariants") },
  {OutputType::out_stress,       std::make_tuple(OutputType::outGrp_stress, 0, com_mod.nsymd, "Stress") },
  {OutputType::out_temperature,  std::make_tuple(OutputType::outGrp_Y,      0, 1,             "Temperature") },
  {OutputType::out_traction,     std::make_tuple(OutputType::outGrp_trac,   0, nsd,           "Traction") },
  {OutputType::out_velocity,     std::make_tuple(OutputType::outGrp_Y,      0, nsd,           "Velocity") },
  {OutputType::out_viscosity,    std::make_tuple(OutputType::outGrp_Visc,   0, 1,             "Viscosity") },
  {OutputType::out_voltage,      std::make_tuple(OutputType::outGrp_Y,      0, 1,             "Action_potential") },
  {OutputType::out_vortex,       std::make_tuple(OutputType::outGrp_vortex, 0, 1,             "Vortex") },
  {OutputType::out_vorticity,    std::make_tuple(OutputType::outGrp_vort,   0, maxNSD,        "Vorticity") },
  {OutputType::out_WSS,          std::make_tuple(OutputType::outGrp_WSS,    0, maxNSD,        "WSS")},
  {OutputType::out_MBF,          std::make_tuple(OutputType::outGrp_MBF,    0, 1,             "MBF")}
};

