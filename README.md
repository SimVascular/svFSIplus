# svFSIplus 
## Table of Contents
[Introduction](#introduction)<br>
[Code Organization](#organization)<br>
[Translating Fortran into C++](#translate)<br>
[C++ Simulation Class](#simulation_class)<br>
[C++ Vector Class](#vector_class)<br>

[C++ Programming](#cpp_programming)<br>

<h1 id="introduction"> Introduction </h1>

svFSIplus is a C++ implementation of the Fortran [svFSI](https://github.com/SimVascular/svFSI) multi-physics finite element solver designed for computational modeling of the cardiovascular system. It represents the <i>First Stage</i> in the development of a C++ multi-physics finite element solver and is essentailly a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. This was done to maintain a simple mapping between the code of the two versions and to faciliate debugging by comparing intermediate data (i.e. element stiffness matrices) and simulation results.

The *Second Stage* of the solver development will be an entirely new implementation encapsilating portions of the *First Stage* code into C++ classes forming a clean and simple abstaction that can be enhanced by any developer.

The C++ implementation differs from the Fortran implementation in three fundamental ways
1) Custom C++ Array and Vector classes to reproduce Fortran dynamic arrays 
2) XML format file replaces the plain-text input file to specify simulation parameters
3) Direct calls to the VTK API replaces the custom code used to read/write VTK format files  

The following sections describe how the C++ implementation is organized and how it replicates the data structures and flow of control of the Fortran implementation. Some important details of the C++ implementation will also be discussed.


<h1 id="organization"> Code Organization </h1>

The C++ implementation attempts to replicate the data structures and flow of control of the Fortran implementation and to maintains its organization.

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

All Fortan subroutines located in a particular file will typically have a C++ implementation in a similarly named file. This was done to maintain a simple mapping between the locations of the C++ and Fortran code.

C++ functions are defined within a `namespace` defined for each Fortran file. For example, the functioins in  `load_msh.cpp` are defined within the `load_msh` `namespace`. Some `namespaces` are named with a `_ns` suffix to prevent conflicts with function names (e.g. `read_files_ns`).


<h1 id="translate"> Translating Fortran into C++ </h1>

This section provides some details about how the svFSI Fortran code was translated into C++ code. This will help to convert any new Fortran code developed in the Fortan svFSI code not included in svFSIplus. 

<h2 id="translate_vars"> Variable Names </h2>

svFSIplus is essentailly a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. The original Fortran variable names are typically small, contain no underscores for readability and are often ambiguous. However, the **same varaible names** are used in both the C++ and Fortran versions in order to maintain a clear correspondence between the variables used in the two codes. 


<h2 id="translate_modules"> Fortran Modules </h2>

Modules were introduced in Fortran to modualize a large code by splitting it into separate files containing procedures and data specific to a certain application. A module is like a C++ class because it can encapsulate both data and procedures. The svFSI Fortran code uses modules primarily to store and access global variables. 

C++ classes are used to implement Fortran modules. Fortran variable names are retained to prevent (or maintain) confusion. The C++ module name uses the same Fortan name in camel case. For example, several of the Fortan module names and the files that implements them are given below with the coressponding C++ class name and implementation files.

```
   ================================================================================================
                Fortran Module Name (file)        |                C++ Class Name (file)
   ================================================================================================
                CEPMOD (CEPMOD.f)                 |               CepMod (CepMod.h,cpp)
   ------------------------------------------------------------------------------------------------
                CHNLMOD (CHNL.f)                  |               ChnlMod (ChnlMod.h,cpp)
   ------------------------------------------------------------------------------------------------
                COMMOD (MOD.f)                    |               ComMod (ComMod.h,cpp)
   ------------------------------------------------------------------------------------------------
 
  ```         
             
The Fortan `USE` command provides access to all the variables defined in a module. Almost all of the svFSI Fortran procedures have a `USE COMMOD` command that provides access to all of the global varaibles (about 90) defined in the `COMMOD` module. For example
```
      SUBROUTINE CONSTRUCT_uSOLID(lM, Ag, Yg, Dg)

      USE COMMOD
      USE ALLFUN
```

**svFSIplus does not use any global variables.**  A C++ module object is passed to each procedure that needs to access its varaibles. For example, in C++ the `ComMod` object `com_mod` is explicity passed to the `construct_usolid` procedure. All C++ modules are stored in the [Simulation](#simulation_class) class.
```
void construct_usolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag,
    const Array<double>& Yg, const Array<double>& Dg)
```

<h2 id="translate_arrays"> Fortran Dynamic Arrays </h2>

Fortran dynamic arrays have been reproduced using custom [Vector](#vector_class), `Array` and `Array3` C++ classes. 




Note that the custom `Vector`, `Array` and `Array3` classes will most likey be replacesd by a more sophisticated implementation such as `Eigen`.

<h1 id="simulation_class"> C++ Simulation Class </h1>


<h1 id="vector_class"> C++ Vector Class </h1>


<h1 id="cpp_programming"> C++ Programming </h1>





