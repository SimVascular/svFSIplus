# svFSIplus 
## Table of Contents
[Introduction](#introduction)<br>
[Code Organization](#organization)

<h1 id="introduction"> Introduction </h1>

svFSIplus is a C++ implementation of the Fortran [svFSI](https://github.com/SimVascular/svFSI) multi-physics finite element solver designed for computational modeling of the cardiovascular system. It represents the <i>First Stage</i> in the development of a C++ multi-physics finite element solver and is essentailly a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. The *Second Stage* of the solver development will be an entirely new implementation encapsilating the *First Stage* code into C++ classes forming a clean and simple abstaction that can be enhanced by any developer.

The C++ implementation differs from the Fortran implementation in three fundamental ways
1) Custom C++ Array and Vector classes to reproduce Fortran dynamic arrays 
2) XML format file replaces the plain-text input file to specify simulation parameters
3) Direct calls to the VTK API replaces the custom code used to read/write VTK format files  

The following sections describe how the C++ implementation is organized and how it replicates the data structures and flow of control of the Fortran implementation. Some important details of the C++ implementation will also be discussed.


<h1 id="organization"> Organization </h1>

svFSIplus is essentailly a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. The C++ implementation attempts to replicate the data structures and flow of control of the Fortran implementation and maintains its organization.

Most of the Fortran code is replicated in C++ using the same file and subroutine names converted to lower case with underscores added to improve readability. For example
```
   ================================================================================================
                Fortran                       |                      C++ 
   ================================================================================================
           SUBROUTINE READFILES               |                  read_files()
   ------------------------------------------------------------------------------------------------
           SUBROUTINE READMSH                 |                  read_msh()
   ------------------------------------------------------------------------------------------------
               LOADMSH.f                      |                  load_msh.cpp
   ------------------------------------------------------------------------------------------------
               VTKXML.f                       |                  vtk_xml.cpp
   ------------------------------------------------------------------------------------------------
```

All Fortan subroutines located in a particular file will typically have a C++ implementation in a similarly named file. 

C++ functions are defined within a `namespace` defined for each Fortran file. For example, the functioins in  `load_msh.cpp` are defined within the `load_msh` `namespace`. Some `namespaces` are named with a `_ns` suffix to prevent conflicts with function names (e.g. `read_files_ns`).


