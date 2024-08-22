
# **Problem Description**

Simulate the propagation of an electrical signal inside a 2D plane using the ten-Tusscher-Panfilov cell activation model.

## Pacemaker and Non-Pacemaker Cells

The Aliev-Panfilov cell activation model is used to represent non-pacemaker myocytes, cells which need to receive stimulus from neighboring cells in order to start the depolarization (action potential begins with the voltage becoming more positive). 

The pacemaker cells initiating a cardiac action potential are simulated by applying an external stimulus to a region representing a group of cells using an Aliev-Panfilov cell activation model. 

This group of pacemaker cells is defined in the solver input XML file using a `Domain`

 `Domain` object. 

```
   Domain file path: mesh/h0.25/domain_info.dat
```

Here, a list of integers is provided for each element to serve as its domain ID. Elements with ID 1 are non-pacemaker cells and those with ID 2 will act as pacemaker cells.

```
   Domain: 1 {
      Electrophysiology model: AP
      Conductivity (iso): 1.0
      ODE solver: Euler
   }

   Domain: 2 {
      Electrophysiology model: AP
      Conductivity (iso): 1.0
      Stimulus: Istim {
         Amplitude: 10.0
         Start time: 0.0
         Duration: 10.0
      }
      ODE solver: Euler
   }
```

To define an external stimulus, users need to provide amplitude, start time and duration. The signal is essentially a square wave.

## References
S. Göktepe and E. Kuhl. Computational modeling of cardiac electrophysiology: A novel finite
element approach. International Journal for Numerical Methods in Engineering, 79(2):156–
178, jul 2009.


