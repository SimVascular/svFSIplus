
# **Problem Description**

Solve the same problem as in [fluid/pipe_RCR_3d](../pipe_RCR_3D). Instead of using the RCR within the `svFSIplus` solver, this example demonstrates how to set up RCR boundary condition in a more generalized framework using sv0DSolver.

## Introduction

Both sv0DSolver and genBC (see [pipe_RCR_genBC](../pipe_RCR_genBC)) provide a framework to programmatically define custom inflow and outflow boundary conditions for a CFD simulation. The framework allows users to create an arbitrary lumped parameter network (LPN, or 0D model) layout suitable for their application. Some common examples include a lumped parameter heart model that models contraction of the heart chambers to use as an inlet boundary condition, sophisticated models of the downstream circulation for various areas of the body such as the legs and upper body, or a closed-loop formulation where all outflow of the SimVascular model returns back to the inflow after passing through the veins, heart, and pulmonary arteries.

**Essentially, sv0D and genBC are two different implementations of the same functionality, and sv0D is the preferred choice.**  genBC is a legacy code developed for [svSolver](https://github.com/SimVascular/svSolver) and requires direct input of a system of equations. sv0DSolver has pre-defined blocks with associated equations and is more user-friendly. Still, `svFSIplus` provides backward compatibility for genBC so that svSolver users can migrate to the new solver easily. 

## Configuration of genBC

There are excellent tutorials online that show the users how to set up genBC step-by-step.
SimVascular website: https://simvascular.github.io/docsGenBC.html
Youtube tutorial: https://www.youtube.com/watch?v=znfV0XLV79s&ab_channel=SimVascular

**We strongly encourage users to go through these tutorials first to become familiar with the concept and the workflow.** 

The input file [svFSI_genBC.inp](./svFSI_genBC.inp) follows the master input file as a template. Some specific input options are discussed below:

```
   <Couple_to_genBC type="SI">
      <ZeroD_code_file_path> genBC/genBC.exe </ZeroD_code_file_path>
   </Couple_to_genBC>
```

This tells the solver that the 0d models will be calculated through genBC. Options to couple 0D codes with svFSI are `N`: none; `I`: implicit; `SI`: semi-implicit; `E`: explicit.

```
   <Add_BC name="lumen_inlet" > 
      <Type> Dir </Type> 
      <Time_dependence> Unsteady </Time_dependence> 
      <Temporal_values_file_path> lumen_inlet.flw</Temporal_values_file_path> 
      <Zero_out_perimeter> true </Zero_out_perimeter> 
      <Impose_flux> true </Impose_flux> 
   </Add_BC> 

   <Add_BC name="lumen_outlet" > 
      <Type> Neu </Type> 
      <Time_dependence> Coupled </Time_dependence> 
   </Add_BC> 
```

In this example, we use the LPN for just the outlet RCR boundary condition and use a file to specify the inlet flow conditions.

### Matching faces between 3D and 0D

In genBC, user has to provide face tags in the USER.f file and make sure that these tags match the ones in solver.inp file for svSolver. In essence, when using genBC, it is the user's responsibility to match both the bcs provided in the svPresolver and svSolver files, and also in the USER.f file.

### Implementation of 0D model

The 0D model or LPN is defined in USER.f. Here are the boundary conditions implemented in genBC:

```fortran
!     BC
      Rp = 121D0
      C  = 1.5D-4
      Rd = 1212D0

      f(1) = (1D0/C) * (Q(1) - x(1)/Rd)
      offset(1) = Q(1)*Rp

!     Assign the additional parameters to be printed       
      Xprint(1)=t
      Xprint(2)= Q(1)
      Xprint(3) = offset(1)

```

Here, we specify the RCR boundary used at the outlet, with the corresponding formulation defined in `f(1)`.

In genBC, `Q(:)` and `P(:)` are defined as flow rates of Neumann faces and pressure of Dirichlet faces, respectively. It is the user's responsibility to carefully match the component `Q(i)` to the corresponding face, which can be error-prone when the number of faces are large. On the other hand, in cplBC, `Q(i)` simply represents the flow rate on the ith face.


### Outputs from 0D model

genBC writes out the results in `AllData` in the current folder, and it includes all the unknowns and the user-specified outputs. In this case, the order will be

```
outlet_pressure  time  outlet_flux  offset(1)
```

Here, the last three are user-defined outputs.

### Initial conditions

In genBC, the initial conditions are specified in USER.f through variable `tZeroX`. Hence, user needs to recompile genBC every time it changes. 
