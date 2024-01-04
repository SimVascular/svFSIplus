
# **Problem Description**

This is an example of using the coupled momentum method (CMM) with variable wall properties and a user-defined Dirichlet inflow condition. For continuous testing, march through 3 time steps only. Pulsatile inflow and resistance are used as the inflow and outflow boundary conditions (BCs), respectively.

## Solution workflow

The 3-step procedure is similar to the 3D pipe CMM example. To use variable wall properties, thickness and Young's modulus need to be prescribed for the deformable wall by providing a VTP file (e.g., svFSI_var_wall_props.vtp). In the solver input file, we define:

<Variable_wall_properties mesh_name="wall">
    <Wall_properties_file_path> ../svFSI_var_wall_props.vtp </Wall_properties_file_path>
</Variable_wall_properties>

The wall property file can be prepared by adding two point arrays, "Thickness" and "Elasticity modulus," to the walls_combined.vtp in the mesh-complete.

Similarly, a user-defined inflow velocity profile could be achieved by providing a bct.vtp file that prescribes time-velocity data for each inlet node.

<Add_BC name="inlet_aorta" > 
   <Type> Dir </Type>
   <Time_dependence> General </Time_dependence>
   <Bct_file_path> ../bct.vtp </Bct_file_path>
</Add_BC> 




