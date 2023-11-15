
# **Problem Description**

This example simulates pulsatile flow in a pipe with a deformable wall using coupled momentum method (CMM) for continuous testing. Only a few time steps are tested and compared. Pusaltile inflow and RCR are used as the inflow and outflow boundary conditions (BCs), respectively. 

## Solution workflow

This problem is solved in three major steps:

1. Solve for fluid flow with a rigid wall and using the same inflow and outflow BCs. 

2. a) Inflate the wall under a diastolic pressure loading. or b) Prestress the wall under a distolic pressure loading.

3. Start CMM simulation with the intial pressure and velocity fields from the rigid wall simulation and the intial displacement field given by 2a or the presstress field by 2b.

## Rigid Wall Simulation (Step 1)

A python script (calcMeanPressTrac.py) is available in the rigid wall solution folder to obtain a mean pressure/traction on the wall used for **Setp 2**.

## Inflation (Step 2)

Four settings in the input file need closer attention when performing inflation step for CMM:

### (I) Load mesh as shell surface

The inflation step is performed on the wall set as a shell,

<Add_mesh name="wall"> 
  <Set_mesh_as_shell> true </Set_mesh_as_shell>
  <Mesh_file_path> ../../../fluid/pipe_RCR_3d/mesh/mesh-surfaces/lumen_wall.vtp </Mesh_file_path>
</Add_mesh>

### (II) Set Initialize keyword

Within <Add_equation type="CMM">, the keyword <Initialize> needs to be set to "inflate":

<Initialize> inflate </Initialize>

Alternately, one could also prestress the vessel wall by using the prestress option:

<Prestress> true </Prestress>
<Initialize> prestress </Initialize>

and output the stress field as well:

<Output type="Spatial">
  <Displacement> true </Displacement>
  <Stress> true </Stress>
</Output>


### (III) Diastolic loading

CMM is initialized via inflation or prestress approaches by applying a diastolic load on the wall. This is performed either using a constant pressure or traction, or from a spatially-varying pressure or traction computed from the rigid wall flow simulation. Examples are provided below for constant pressure,

<Add_BF mesh="wall">
  <Type> Neumann </Type>
  <Value> 8000.0 </Value>
</Add_BF>

and, spatially-varying traction field as,

<Add_BF mesh="wall">
  <Type> Trac </Type>
  <Time_dependence> Spatial </Time_dependence>
  <Temporal_and_spatial_values_file_path> rigid_wall_mean_traction.vtp </Temporal_and_spatial_values_file_path>
</Add_BF>

Note that because the wall is modeled as a shell surface, the above setting applies pressure/traction as a *body force* on the shell surface using `Add BF` keyword and **not** as a *surface traction* using the `Add BC` keyword.


## FSI using CMM (Step 3)
To run a CMM simulation, complete mesh files (instead of a wall surface mesh) are added into the solver input file. The parameters including <Prestress> and <Initialize> for **Step 2** need to be removed. 
 
### (I) Initialize flow field

For faster convergence, we initialize the flow for the FSI simulation using data from **Step 1** within `Add mesh` keyword as,

<Add_mesh name="msh" >
...
  <Initial_velocities_file_path> ../1-rigid-solution/result_800_cpp.vtu </Initial_velocities_file_path>
  <Initial_pressures_file_path> ../1-rigid-solution/result_800_cpp.vtu </Initial_pressures_file_path>
...
</Add_mesh>

In this case, the solution from the rigid wall simulation at time step 800 is used to initialize the velocity and pressure fields.

### (II) CMM loading initial displacement/prestress fields

Results from **Step 2** will be used for the wall BC. If an inflation procedure is used in **Step 2**, the following BC setting is used:

<Add_BC name="lumen_wall" > 
  <Type> CMM </Type> 
  <Initial_displacements_file_path> ../2a-inflate/result_003.vtu </Initial_displacements_file_path>
</Add_BC>


If the prestress procedure is employed in **Step 2**, replace <Initial_displacements_file_path> with <Prestress_file_path>:

<Add_BC name="lumen_wall" > 
  <Type> CMM </Type> 
  <Prestress_file_path> ../2b-prestress/result_003.vtu </Prestress_file_path>
</Add_BC> 




