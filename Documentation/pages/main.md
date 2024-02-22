@mainpage svFSI

[TOC]

The following sections describe how to build and use svFSIplus. Implementation details can be found [here](@ref implementation). See the [testing guide](@ref testing) for instructions on running the provided examples.

# Solver Parameter Input XML File

The svFSI solver reads simulation parameters from a [custom text format](https://github.com/SimVascular/svFSI/blob/master/svFSI_master.inp) file.   Simulation parameters are set using keyword/value pairs separated by a colon. A keyword consists of series of alphanumeric characters separated by spaces.

Keyword/value example
```json
Coupled: 1
Min iterations: 1
Tolerance: 1e-6
```

Braces **{}** are used to define simulation parameters with sub-elements. For example the `Domain` keyword is associated with
several other sub-elements of keyword/values pairs.
```json
Domain: 0 {
  Equation: fluid
  Density: 1.0
  Viscosity: Constant {Value: 0.04}
  Backflow stabilization coefficient: 0.2
}
```

svFSIplus solver simulation parameters are stored in an XML-format file. The XML file organization and keyword names replicate the original text format using the following conversion rules 
- Spaces in keywords names are replaced by underscores
- An additional value after a colon is replaced by an XML attribute used to identify the value. For example: `Add equation: FSI` is replaced by `<Add_equation type="FSI" >`.
- Parameters with sub-elements are replaced by an XML element with sub-elements. For example `Add equation: FSI { sub-elements }` is replaced by `<Add_equation type="FSI" > xml-sub-elements </Add_equation>`.

Example: The original text 
```mfs
Add equation: FSI {
   Coupled: 1
   Min iterations: 1
   Max iterations: 10
   Tolerance: 1e-6

   Domain: 0 {
      Equation: fluid
      Density: 1.0
      Viscosity: Constant {Value: 0.04}
      Backflow stabilization coefficient: 0.2
   }
}
   
```
is replaced by new XML format with
```xml
<Add_equation type="FSI" >
   <Coupled> true </Coupled>
   <Min_iterations> 1 </Min_iterations>
   <Max_iterations> 10 </Max_iterations>
   <Tolerance> 1e-6 </Tolerance>

   <Domain id="0" >
      <Equation> fluid </Equation>
      <Density> 1.0 </Density>
      <Viscosity model="Constant" >
         <Value> 0.04 </Value>
      </Viscosity>
      <Backflow_stabilization_coefficient> 0.2 </Backflow_stabilization_coefficient>
   </Domain>
</Add_equation>
```


# Pre-built svFSIplus Binaries

Pre-built svFSIplus binaries can be downloaded from the [SimVascular SimTK downloads page](https://simtk.org/frs/?group_id=188).
Installers are provided for MacOS and Ubuntu operating systems. 

The compilers and other software used to build the binaries (listed on the downloads page) must be installed on the target 
platform in order to use the svFSIplus binary.


# Building svFSIplus from Source

The svFSIplus source code can be downloaded and built on computers running MacOS and other Linux operating systems. Builds are most 
commonly done on MacOS and Ubuntu platforms.

**The Windows operating system is currently not supported.**


## Build Process

Several [software packages](#building_packages) must be installed in order to build and use svFSIplus.

svFSIplus is built on a MacOS or Linux computer from a terminal using the following commands
```bash
mkdir svFSIplus-package
cd svFSIplus-package
git clone https://github.com/SimVascular/svFSIplus
mkdir build
cd build/
cmake ../svFSIplus/
make -j4
```

This creates the svFSIplus binary `svFSIplus-package/build/svFSI-build/bin/svFSIplus*`.

The binary can be installed in /usr/local using the following commands
- cd svFSI-build
- sudo make install

The `svFSIplus*` binary is now located in /usr/local/SV/bin.

It is convenient to update the PATH environment variable with the location of the `svFSIplus*` binary. Another option is to create
an alias for `svFSIplus`. For example: `alias svFsiPlus=/usr/local/SV/bin/svFSIplus`


## Required Software Packages

The following software packages are required to build svFSIplus

- CMake
- C++17 compiler 
- Visualization Toolkit
- Open MPI

Software packages can be installed using package management tools like [homebrew](https://brew.sh) and [apt](https://en.wikipedia.org/wiki/APT_(software)).


### CMake

svFSIplus is built using [CMake](https://cmake.org). CMake is used to control the software compilation process using the `CMakeLists.txt` files located within the svFSIplus source. CMake first attempts to discover the software dependencies (e.g. C++ compiler) required to build svFSIplus. It will then generate the native [makefiles](https://en.wikipedia.org/wiki/Make_(software)) used to compile it.


### Visualization Toolkit

svFSIplus uses the [The Visualization Toolkit (VTK)](https://vtk.org/) to read finite element mesh data stored in VTP and VTU [formats](https://kitware.github.io/vtk-examples/site/VTKFileFormats/). The [VTP](https://vtk.org/doc/nightly/html/classvtkPolyData.html#details) format is used to store 2D finite element meshes for inlet and outlet surfaces. The [VTU](https://vtk.org/doc/nightly/html/classvtkUnstructuredGrid.html#details) format is used to store finite element volume meshes. svFSIplus can also write simulation results to VTK format files. 


### Open MPI

svFSIplus uses the [Open MPI](https://www.open-mpi.org/) Message Passing Interface (MPI) libraries for parallel computing. 


## Using Trilinos with svFSIplus

By default svFSIplus uses custom linear algebra solvers (e.g. GMRES) in a simulation. It also supports several simple preconditioners. 

svFSIplus can also be built with the linear algebra solvers and preconditioners provided by the [Trilinos](https://trilinos.github.io/) scientific software package. Using Trilinos can often reduce the execution time of a simulation.


### Building Trilinos

The Trilinos libraries need to be built from the source downloaded from [here](https://github.com/trilinos/Trilinos). See the [Quick configure, build and install hints for Trilinos](https://github.com/trilinos/Trilinos/blob/master/INSTALL.rst) for the instructions about how to build Trilinos.

Build with the third-party libraries (TPLs) 
 - Boost
 - BLAS
 - HDF5
 - HYPRE
 - LAPACK
 - MUMPS
 
and with The Trilinos packages  
- Amesos
- AztecOO
- Epetra
- EpetraEXT
- Ifpack
- ML
- MueLU
- ROL
- Sacado
- Teuchos
- Zoltan


### Compiling with Trilinos Libraries

svFSIplus is compiled with the Trilinos libraries using the CMake `SV_USE_TRILINOS` option
```
cmake ../svFSIplus/ -DSV_USE_TRILINOS:BOOL=ON 
```
CMake will be able to find the Trilinos libraries once they have been built and installed.


# Running a Simulation

The svFSIplus binary is run from the command line with an input solver parameters XML file as an argument

```bash
svFSIplus solver_params.xml
```

A simulation can be run in parallel on four processors using 

```bash
mpiexec -np 4 svFSIplus solver_params.xml
```
