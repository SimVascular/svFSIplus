
# **Problem Description**

Simulate unsteady fluid flow in a pipe.

The simulation differs from the <a href="https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/fluid/pipe_RCR_3d"> Fluid RCR 3D Pipe </a> test only in the linear algebra package it uses.

The simulation uses the PETSc linear algebra package with the **petsc-rcs** preconditioner.
```
<LS type="GMRES" >
  <Linear_algebra type="petsc" >
    <Preconditioner> petsc-rcs </Preconditioner>
  </Linear_algebra>
  <Max_iterations> 100 </Max_iterations>
  <Tolerance> 1e-12 </Tolerance>
</LS>
```

