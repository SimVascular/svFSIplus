
# **Problem Description**

Simulate pressure wave propagation in an arterial model using the Arbitrary Lagrangian-Eulerian method

The simulation differs from the <a href="https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/fsi/pipe_3d"> FSI in a 3D Pipe </a> test only in the linear algebra package it uses.

The simulation uses the Trilinos linear algebra package with the **trilinos-ml** preconditioner.
```
<Linear_algebra type="trilinos" >
  <Preconditioner> trilinos-ml </Preconditioner>
</Linear_algebra>
```
