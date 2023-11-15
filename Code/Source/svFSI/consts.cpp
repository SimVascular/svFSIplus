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

#include "consts.h" 

namespace consts {


/// @brief Reproduces the 'SELECT CASE (lM%eType)' statements in the
/// Fortran 'PARTMSH' subroutine. 
//
const std::map<ElementType,int> element_type_to_elem_nonb = 
{
  {ElementType::LIN1,  1 },
  {ElementType::LIN2,  1 },
  {ElementType::TRI3,  2 },
  {ElementType::TRI6,  3 },
  {ElementType::QUD4,  2 },
  {ElementType::QUD8,  3 },
  {ElementType::QUD9,  3 },
  {ElementType::TET4,  3 },
  {ElementType::TET10, 6 },
  {ElementType::HEX8,  4 },
  {ElementType::HEX20, 8 },
  {ElementType::HEX27, 9 },
  {ElementType::WDG,   3 }
};

/// @brief Map for constitutive_model string to ConstitutiveModelType. 
///
/// Type of constitutive model (isochoric) for structure equation:
/// St.Venant-Kirchhoff, modified St.Venant-Kirchhoff, NeoHookean,
/// Mooney-Rivlin, modified Holzapfel-Gasser-Ogden with dispersion,
/// Linear model (S = mu*I), Guccione (1995), Holzapfel & Ogden model
//
const std::map<std::string,ConstitutiveModelType> constitutive_model_name_to_type = 
{
  {"lin", ConstitutiveModelType::stIso_lin}, 
  {"linear", ConstitutiveModelType::stIso_lin},

  {"stVK", ConstitutiveModelType::stIso_StVK},
  {"stVenantKirchhoff", ConstitutiveModelType::stIso_StVK},

  {"m-stVK", ConstitutiveModelType::stIso_mStVK},
  {"modified-stVK", ConstitutiveModelType::stIso_mStVK},
  {"modified-stVenantKirchhoff", ConstitutiveModelType::stIso_mStVK},

  {"nHK", ConstitutiveModelType::stIso_nHook}, 
  {"nHK91", ConstitutiveModelType::stIso_nHook},
  {"neoHookean", ConstitutiveModelType::stIso_nHook},
  {"neoHookeanSimo91", ConstitutiveModelType::stIso_nHook},

  {"MR", ConstitutiveModelType::stIso_MR}, 
  {"Mooney-Rivlin", ConstitutiveModelType::stIso_MR},

  {"HGO", ConstitutiveModelType::stIso_HGO}, 

  {"Guccione", ConstitutiveModelType::stIso_Gucci},
  {"Gucci", ConstitutiveModelType::stIso_Gucci},

  {"HO", ConstitutiveModelType::stIso_HO}, 
  {"Holzapfel", ConstitutiveModelType::stIso_HO},

  {"quad", ConstitutiveModelType::stVol_Quad},
  {"Quad", ConstitutiveModelType::stVol_Quad},
  {"quadratic", ConstitutiveModelType::stVol_Quad},
  {"Quadratic",ConstitutiveModelType::stVol_Quad},

  {"ST91", ConstitutiveModelType::stVol_ST91},
  {"Simo-Taylor91", ConstitutiveModelType::stVol_ST91},

  {"M94", ConstitutiveModelType::stVol_M94},
  {"Miehe94", ConstitutiveModelType::stVol_M94},

};

/// @brief Map for contact model string name to ContacteModelType
//
const std::map<std::string,ContactModelType> contact_model_name_to_type =
{
  {"penalty", ContactModelType::cntctM_penalty},
  {"potential", ContactModelType::cntctM_potential},
};


/// @brief Map for fluid viscosity model string to FluidViscosityModelType. 
//
const std::map<std::string,FluidViscosityModelType> fluid_viscosity_model_name_to_type
{
  {"const", FluidViscosityModelType::viscType_Const},
  {"constant", FluidViscosityModelType::viscType_Const},
  {"newtonian", FluidViscosityModelType::viscType_Const},

  {"cy", FluidViscosityModelType::viscType_CY},
  {"carreau-yasuda", FluidViscosityModelType::viscType_CY},

  {"cass", FluidViscosityModelType::viscType_Cass},
  {"cassons", FluidViscosityModelType::viscType_Cass},

};


/// @brief Map number of element nodes to element type.
///
/// Note that this is not an ambiguous mapping between, 
/// need to handle some overlap in number of element nodes.
//
static const std::map<int,ElementType> num_nodes_to_type = {
    {1, ElementType::PNT},
    {2, ElementType::LIN1}, 
    {2, ElementType::LIN2}, 
    {3, ElementType::TRI3},
    {6, ElementType::TRI6},
    {4, ElementType::QUD4},
    {8, ElementType::QUD8},
    {9, ElementType::QUD9},
    {4, ElementType::TET4},
    {10, ElementType::TET10}, 
    {8, ElementType::HEX8},
    {20, ElementType::HEX20},
    {27, ElementType::HEX27},
    {6, ElementType::WDG}, 
    {0, ElementType::NRB}
  };

const std::map<std::string,CplBCType> cplbc_name_to_type = {
    {"E", CplBCType::cplBC_E},
    {"I", CplBCType::cplBC_I},
    {"SI", CplBCType::cplBC_SI}
  };

/// @brief Map equation name to a type
//
const std::map<std::string,EquationType> equation_name_to_type = {
    {"advection_diffusion", EquationType::phys_heatF},
    {"AD", EquationType::phys_heatF},
    {"heatF", EquationType::phys_heatF},
    {"dyeTransport", EquationType::phys_heatF},
    {"scalarTransport", EquationType::phys_heatF},

    {"cardiac_electro_physiology", EquationType::phys_CEP},
    {"CEP", EquationType::phys_CEP},

    {"coupled_momentum", EquationType::phys_CMM},
    {"CMM", EquationType::phys_CMM},

    {"fluid", EquationType::phys_fluid},

    {"fluid-solid-interaction", EquationType::phys_FSI},
    {"FSI", EquationType::phys_FSI},

    {"linear_elasticity", EquationType::phys_lElas},
    {"lElas", EquationType::phys_lElas},

    {"mesh", EquationType::phys_mesh},

    {"shell", EquationType::phys_shell},

    {"solid_heat", EquationType::phys_heatS},
    {"heatS", EquationType::phys_heatS},
    {"laplace", EquationType::phys_heatS},
    {"poisson", EquationType::phys_heatS},

    {"stokes", EquationType::phys_stokes},

    {"structural", EquationType::phys_struct},
    {"struct", EquationType::phys_struct},

    {"structural_velocity_pressure", EquationType::phys_ustruct},
    {"ustruct", EquationType::phys_ustruct},

  };

const std::map<std::string,MeshGeneratorType> mesh_generator_name_to_type = {
    {"Tetgen", MeshGeneratorType::RMSH_TETGEN},
    {"Meshsim", MeshGeneratorType::RMSH_MESHSIM}
};

/// @brief Map for preconditioner type string to PreconditionerType enum
//
const std::map<std::string,PreconditionerMapType> preconditioner_name_to_type =
{
  {"fsils", std::make_pair(PreconditionerType::PREC_FSILS,false)},
  {"svfsi", std::make_pair(PreconditionerType::PREC_FSILS,false)},

  {"rcs", std::make_pair(PreconditionerType::PREC_RCS,false)},
  {"row-column-scaling", std::make_pair(PreconditionerType::PREC_RCS,false)},

  {"trilinos-diagonal", std::make_pair(PreconditionerType::PREC_TRILINOS_DIAGONAL,true)},

  {"trilinos-blockjacobi", std::make_pair(PreconditionerType::PREC_TRILINOS_BLOCK_JACOBI,true)},
  {"blockjacobi", std::make_pair(PreconditionerType::PREC_TRILINOS_BLOCK_JACOBI,true)},

  {"trilinos-ilu", std::make_pair(PreconditionerType::PREC_TRILINOS_ILU,true)},
  {"trilinos-ilut", std::make_pair(PreconditionerType::PREC_TRILINOS_ILUT,true)},

  {"trilinos-ic", std::make_pair(PreconditionerType::PREC_TRILINOS_IC,true)},
  {"trilinos-ict", std::make_pair(PreconditionerType::PREC_TRILINOS_ICT,true)},

  {"trilinos-ml", std::make_pair(PreconditionerType::PREC_TRILINOS_ML,true)},

  {"petsc", std::make_pair(PreconditionerType::PREC_PETSC,true)}
};

/// @brief Map solver type string to SolverType enum. 
//
const std::map<std::string,SolverType> solver_name_to_type 
{
  {"bi-partitioned", SolverType::lSolver_NS},
  {"ns", SolverType::lSolver_NS},
  {"bpn", SolverType::lSolver_NS},
  {"bipn", SolverType::lSolver_NS},

  {"gmres", SolverType::lSolver_GMRES},

  {"conjugate-gradient", SolverType::lSolver_CG},
  {"cg", SolverType::lSolver_CG},

  {"bi-conjugate-gradient", SolverType::lSolver_BICGS},
  {"bicg", SolverType::lSolver_BICGS},
  {"bicgs", SolverType::lSolver_BICGS}
};


};

