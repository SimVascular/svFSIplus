
# **Problem Description**

Simulate the propagation of an electrical signal inside a 2D plane using the ten-Tusscher-Panfilov cell activation model.

## Pacemaker and Non-Pacemaker Cells

The Aliev-Panfilov cell activation model is used to represent non-pacemaker myocytes, cells which need to receive stimulus from neighboring cells in order to start the depolarization (action potential begins with the voltage becoming more positive). 

The pacemaker cells initiating a cardiac action potential are simulated by applying an external stimulus to a region representing a group of cells using an Aliev-Panfilov cell activation model. 

This pacemaker and non-pacemaker myocytes are defined for the simulation using two svFSIplus `domain`s. The `domain`s comprise separate sets of elements in the finite element mesh identified using integer values (1 or 2) read in from the `domain_info.dat` file using the solver input XML file `Domain_file_path` parameter

```
   <Domain_file_path> mesh/domain_info.dat </Domain_file_path> 
```

The integer values are used to identify `domain`s using
- 1 - non-pacemaker cells
- 2 - pacemaker cells.

A domain-specific electrophysiology model is then defined for the non-pacemaker cells 
```
<Domain id="1" >
  <Electrophysiology_model> AP </Electrophysiology_model>
  <Isotropic_conductivity> 1.0 </Isotropic_conductivity>
  <ODE_solver> Euler </ODE_solver>
 </Domain>
```
and the pacemaker cells
```
<Domain id="2" >
  <Electrophysiology_model> AP </Electrophysiology_model>
  <Isotropic_conductivity> 1.0 </Isotropic_conductivity>
  <ODE_solver> Euler </ODE_solver>
  <Stimulus type="Istim" >
    <Amplitude> 10.0 </Amplitude>
    <Start_time> 0.0 </Start_time>
    <Duration> 10.0 </Duration>
  </Stimulus>
</Domain>
```

The `Stimulus` parameter used in the pacemaker cells `domain 2` defines an external stimulus as a square wave.


## References
S. Göktepe and E. Kuhl. Computational modeling of cardiac electrophysiology: A novel finite
element approach. International Journal for Numerical Methods in Engineering, 79(2):156–178, Jul 2009.


