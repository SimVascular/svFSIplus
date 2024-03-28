
# **Problem Description**

Solve dye transportation with fluid flow in a cylindrical tube with zero neumann boundary condition at the outlet and steady flow at the inlet. The dye is passively transported by the flow through advection and diffusion.

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

## Scalar transport equation
Unlike the test case in `04-fluid/02-dye_AD`, the advection-diffusion equation that governs the dye transportation is added to the input file without a proceeding fluid equation. 
Here we assume that we already have a velocity field that has been precomputed by some prior physics simulation (denoted `precomputed_velocity.vtu`). The specified velocity field solution
is then taken to advect a scalar field. 
This equation is built on top of the heat transfer equation, so some of the terminology in the input file follows those in the heat equation, e.g. "Conductivity", "Temperature".

Some noticeable setting in the input file are:

```
   <Use_precomputed_solution> true </Use_precomputed_solution>
   <Precomputed_solution_file_path> precomputed_velocity.vtu </Precomputed_solution_file_path>
   <Precomputed_solution_field_name> Velocity </Precomputed_solution_field_name>
```

This tell the solver that this is a one-way coupling study, and the dye is passively transported by the flow.

```
   <Output type="Alias" >
     <Temperature> Concentration </Temperature>
   </Output>
```

This renames  "Temperature" to "Concentration" for ease of interpretation.

```
   <Add_BC name="lumen_inlet" >
      <Type> Dirichlet </Type>
      <Time_dependence> Steady </Time_dependence>
      <Value> 1.0 </Value>
   </Add_BC>
```

The dye is constantly released at the inlet.
