
# **Problem Description**

Simulate unsteady fluid flow in a pipe.

The simulation differs from the <a href="https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/fluid/pipe_RCR_3d"> Fluid RCR 3D Pipe </a> test only in the linear algebra package it uses.

The simulation uses the Trilinos linear algebra package with a **trilinos-blockjacobi** preconditioner.
```
<LS type="GMRES" >
  <Linear_algebra type="trilinos" >
    <Preconditioner> trilinos-blockjacobi </Preconditioner>
  </Linear_algebra>
  <Max_iterations> 100 </Max_iterations>
  <Tolerance> 1e-12 </Tolerance>
</LS>
```

