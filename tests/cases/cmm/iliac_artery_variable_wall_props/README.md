
# **Problem Description**

Simulate the fluid-structure interaction in an anatomic model of the abdominal aorta and common iliac arteries using the coupled momentum method (CMM).

## Material properties
The arterial wall material properties are defined using spatially varying values for elasticity modulus and wall thickness. These values are read in from a VTK-format **svFSI_var_wall_props.vtp** file containing the two *PointData DataArray*s named **Thickness** and **Elasticity modulus**.

The **svFSI_var_wall_props.vtp** file is specified in using the **Wall_properties_file_path** parameter
```
<Variable_wall_properties mesh_name="wall">
    <Wall_properties_file_path> ../svFSI_var_wall_props.vtp </Wall_properties_file_path>
</Variable_wall_properties>
```

## Boundary conditions

### Inflow boundary condition
The aorta inlet boundary condition uses user-defined pulsatile flow data read in from a VTK-format **bct.vtp** file containing a time-series of nodal velocities *PointData DataArray*s named **velocity_** suffixed by a time value

**bct.vtp**
```
<VTKFile type="PolyData" version="0.1" byte_order="LittleEndian" header_type="UInt32" compressor="vtkZLibDataCompressor">
  <PolyData>
    <Piece NumberOfPoints="108"                  NumberOfVerts="0"                    NumberOfLines="0"                    NumberOfStrips="0"                    NumberOfPolys="181"                 >
      <PointData Scalars="GlobalNodeID">
        <DataArray type="Int32" Name="GlobalNodeID" format="appended" RangeMin="95"                   RangeMax="6293"                 offset="0"                   />
        <DataArray type="Float64" Name="velocity_0.0000" NumberOfComponents="3" format="appended" RangeMin="0"                    RangeMax="2.1367791043"         offset="228"                 />
        <DataArray type="Float64" Name="velocity_0.0101" NumberOfComponents="3" format="appended" RangeMin="0"                    RangeMax="2.3399017892"         offset="1392"                />
        <DataArray type="Float64" Name="velocity_0.0202" NumberOfComponents="3" format="appended" RangeMin="0"                    RangeMax="2.5857762197"         offset="2564"                />
.
.
.
```
The  **bct.vtp** file is specified using the **Bct_file_path** parameter
```
<Add_BC name="inlet_aorta" > 
   <Type> Dir </Type>
   <Time_dependence> General </Time_dependence>
   <Bct_file_path> ../bct.vtp </Bct_file_path>
</Add_BC> 
```

### Outflow boundary condition
Resistance outflow boundary conditions are used for the common iliac arteries.
```
   <Add_BC name="outlet_aorta" >
      <Type> Neu </Type>
      <Time_dependence> Resistance </Time_dependence>
      <Value> 12000 </Value>
   </Add_BC>


  <Add_BC name="outlet_iliac" >
      <Type> Neu </Type>
      <Time_dependence> Resistance </Time_dependence>
      <Value> 12000 </Value>
   </Add_BC>
```

## Simulation Workflow

The simulation workflow is similar to the workflow used in the <a href="https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/cmm/pipe_3d"> 3D pipe CMM example</a>

1) Solve for fluid flow with a rigid wall and using the same inflow and outflow BCs

2) a) Inflate the wall under a diastolic pressure loading or b) Prestress the wall under a distolic pressure loading

3) Start CMM simulation with the intial pressure and velocity fields from the rigid wall simulation and the intial displacement field given by 2a or the presstress field by 2b











