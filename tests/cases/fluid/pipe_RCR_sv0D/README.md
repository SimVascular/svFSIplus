
# **Problem Description**

Solve the same problem as in [fluid/pipe_RCR_3d](../pipe_RCR_3D). Instead of using the RCR within the `svFSIplus` solver, this example demonstrates how to set up RCR boundary condition in a more generalized framework using sv0DSolver.

## Introduction

Both sv0DSolver and genBC (see [pipe_RCR_genBC](../pipe_RCR_genBC)) provide a framework to programmatically define custom inflow and outflow boundary conditions for a CFD simulation. The framework allows users to create an arbitrary lumped parameter network (LPN, or 0D model) layout suitable for their application. Some common examples include a lumped parameter heart model that models contraction of the heart chambers to use as an inlet boundary condition, sophisticated models of the downstream circulation for various areas of the body such as the legs and upper body, or a closed-loop formulation where all outflow of the SimVascular model returns back to the inflow after passing through the veins, heart, and pulmonary arteries.

**Essentially, sv0D and genBC are two different implementations of the same functionality, and sv0D is the preferred choice.**  genBC is a legacy code developed for [svSolver](https://github.com/SimVascular/svSolver) and requires direct input of a system of equations. sv0DSolver has pre-defined blocks with associated equations and is more user-friendly. Still, `svFSIplus` provides backward compatibility for genBC so that svSolver users can migrate to the new solver easily. 


## Configuration of sv0DSolver

The following files require user's attention: [svFSI.xml](./svFSI.xml), [svzerod_3Dcoupling.json](./svzerod_3Dcoupling.json) and [svZeroD_interface.dat](./svZeroD_interface.dat).

### svFSI.xml

The input file [svFSI_genBC.xml](./svFSI.xml) follows the master input file as a template. Some specific input options are discussed below:

```
   <Couple_to_svZeroD type="SI">
   </Couple_to_svZeroD>
```

This tells the solver that the 0d models will be calculated through sv0DSolver. Options to couple 0D codes with svFSI are `N`: none; `I`: implicit; `SI`: semi-implicit; `E`: explicit.

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

### svzerod_3Dcoupling.json

This is the configuration file for sv0DSolver and contains the elements of the 0D model being coupled to the 3D simulation. 

For more information on the available parameters and elements, documentation is available here: [svZeroDSolver](https://github.com/SimVascular/svZeroDSolver)

**The following are necessary in "simulation_parameters" for a coupled simulation:**
"coupled_simulation": true,
"steady_initial": false

The external coupling block is what connects the 3D element to the 0D model. sv0D allows you to create a name for this element and specify its type (in this case, we are interested in **flow** out of the pipe). It is connected at the **inlet** of the block with the name **RCR**. Values of **time** (t) are set to the beginning and end of a cardiac cycle (0.0 to 1.0 s) and the corresponding **flow values** (Q) are set to 1.0, as this flow will be received from the 3D simulation.

The RCR boundary condition block sets up the RCR element with the desired resistance and pressure values.

```
{
    "simulation_parameters": {
        "coupled_simulation": true,
        "number_of_time_pts": 100,
        "output_all_cycles": true,
        "steady_initial": false
    },
    "boundary_conditions": [
        {
            "bc_name": "RCR",
            "bc_type": "RCR",
            "bc_values": {
                "Rp": 121.0,
                "Rd": 1212.0,
                "C": 1.5e-4,
                "Pd": 0.0
            }
        }
    ],
    "external_solver_coupling_blocks": [
        {
            "name": "RCR_coupling",
            "type": "FLOW",
            "location": "inlet",
            "connected_block": "RCR",
            "periodic": false,
            "values": {
                "t": [0.0, 1.0],
                "Q": [1.0, 1.0]
            }
        }
    ],
    "junctions": [],
    "vessels": []
}
```

### svZeroD_interface.dat

This file sets up the interface between svFSIplus and sv0DSolver. It requires the path of the dynamic library for svZeroDSolver and the input file (svzerod_3Dcoupling.json) discussed above.

This file also matches the external coupling blocks in the 0D model to the coupled surfaces in svFSIplus:
The first element in each line should be the name of the block from the json file and the second element should be the index of the coupled surface in svFSIplus. In this case, there is only one coupled surface with index 0.

```
svZeroD external coupling block names to surface IDs (where surface IDs are from *.svpre file):
RCR_coupling 0
```

The next lines initialize the pressure and flow of these coupled surfaces in the 0D model:
0 indicates that the values will not be initialized, and 1 indicates that they will be initialized to the value provided afterwards.

```
Initialize external coupling block flows:
0

External coupling block initial flows (one number is provided, it is applied to all coupling blocks):
0.0

Initialize external coupling block pressures:
1

External coupling block initial pressures (one number is provided, it is applied to all coupling blocks):
0.0
```