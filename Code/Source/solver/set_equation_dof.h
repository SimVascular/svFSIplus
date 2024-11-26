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

#include <map>
#include <tuple>

/// @brief The 'equation_dof_map' map defined here sets equation dof and sym data members. 
//
using EquationDofType = std::tuple<int, std::string>; 

std::map<consts::EquationType, EquationDofType> equation_dof_map =
{
  {EquationType::phys_fluid,    std::make_tuple(nsd+1, "NS") },
  {EquationType::phys_heatF,    std::make_tuple(1,     "HF") },
  {EquationType::phys_heatS,    std::make_tuple(1,     "HS") },
  {EquationType::phys_lElas,    std::make_tuple(nsd,   "LE") },
  {EquationType::phys_struct,   std::make_tuple(nsd,   "ST") },
  {EquationType::phys_ustruct,  std::make_tuple(nsd+1, "ST") },
  {EquationType::phys_CMM,      std::make_tuple(nsd+1, "CM") },
  {EquationType::phys_shell,    std::make_tuple(nsd,   "SH") },
  {EquationType::phys_FSI,      std::make_tuple(nsd+1, "FS") },
  {EquationType::phys_mesh,     std::make_tuple(nsd,   "MS") },
  {EquationType::phys_CEP,      std::make_tuple(1,     "EP") },
  {EquationType::phys_stokes,   std::make_tuple(nsd+1, "SS") }
};

