
# **Problem Description**

Solve electrophysiology inside a 2D plane. The Aliev-Panfilov model is used to describe the cell activation. For more information regarding the model please refer to the following publication:

> S. Göktepe and E. Kuhl. Computational modeling of cardiac electrophysiology: A novel finite
> element approach. International Journal for Numerical Methods in Engineering, 79(2):156–
> 178, jul 2009.

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

## Pacemaker and Non-Pacemaker Cells

The Aliev-Panfilov model is developed to model the non-pacemaker myocytes, which need to receive stimulus from neighboring cells to start the depolarization. In this example, instead of using a different electrophysiological model for the pacemaker cells, we add external stimulus to the Aliev-Panfilov model in a small group of cells so that they can behave like pacemaker cells and initiate the wave propagation. This is achieved through `Domain` object. 

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

