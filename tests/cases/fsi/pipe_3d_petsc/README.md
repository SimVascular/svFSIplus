
# **Problem Description**

Simulate pressure wave propagation in an arterial model using the Arbitrary Lagrangian-Eulerian method

The simulation differs from the <a href="https://github.com/SimVascular/svFSIplus/tree/main/tests/cases/fsi/pipe_3d"> FSI in a 3D Pipe </a> test only in the linear algebra package it uses.

The simulation uses the PETSc linear algebra package with the **petsc-jacobi** preconditioner.
```
<LS type="GMRES" >
  <Linear_algebra type="petsc" >
    <Preconditioner> petsc-jacobi </Preconditioner>
  </Linear_algebra>
  <Max_iterations> 100 </Max_iterations>
  <Tolerance> 1e-12 </Tolerance>
  <Krylov_space_dimension> 50 </Krylov_space_dimension>
</LS>
```
