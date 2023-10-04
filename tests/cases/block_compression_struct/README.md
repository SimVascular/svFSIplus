
# **Problem Description**

Solve block compression problem using displacement-based solid equation. The problem set-up is as follows [1]:

<p align="center">
   <img src="./configuration.png" width="600">
</p>

The final displacement is plotted below.

<p align="center">
   <img src="./displacement.png" width="600">
</p>

The input file `svFSI.inp` follows the master input file [`svFSI_master.inp`](./svFSI_master.inp) as a template. Some specific input options are discussed below:

## Directional Dirichlet BC

The problem configuration requires enforcing Dirichlet BC along a specific direction. For example, mesh points on the patch should have zero displacement along x and y directions. 

```
   Add BC: patch {
      Type: Dir
      Value: 0.0
      Effective direction: (1, 1, 0)
   }
```

This is achieved through `Effective direction` command. Currently, the software only supports specify Dirichlet BC along Cartesian directions. For example, `(1,0,0)` or `(1,0)` sets Dirichlet BC along x-axis in 3D or 2D.



## Reference
1. Liu, Ju, and Alison L. Marsden.  A Unified Continuum and Variational Multiscale Formulation for Fluids, Solids, and Fluid Structure Interaction.  *Computer Methods in Applied Mechanics and Engineering* 337 (August 2018): 549 97. https://doi.org/10.1016/j.cma.2018.03.045.
