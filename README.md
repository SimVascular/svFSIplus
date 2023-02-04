# svFSIplus 

## Table of Contents
[Introduction](#introduction)<br>
[Solver Parameter Input XML File ](#xml_file)<br>
[Pre-built svFSIplus Binaries](#executables)<br>
[Building svFSIplus from Source](#building)<br>
[Running a Simulation](#simulation_run)


<!--- ====================================================================================================================== --->
<!--- ===================================================== Introduction  ================================================== --->
<!--- ====================================================================================================================== --->

<h1 id="introduction"> Introduction </h1>

svFSIplus is a C++ implementation of the Fortran [svFSI](https://github.com/SimVascular/svFSI) multi-physics finite element solver designed for computational modeling of the cardiovascular system. It represents the <i>First Stage</i> in the development of a C++ multi-physics finite element solver and is essentailly a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. This was done to maintain a simple mapping between the code of the two versions and to faciliate debugging by comparing intermediate data (i.e. element stiffness matrices) and simulation results.

The *Second Stage* of the solver development will be an entirely new implementation encapsilating portions of the *First Stage* code into C++ classes forming a clean and simple abstaction that can be enhanced by any developer.

The following sections describe how to build and use svFSIplus. Implementaton details can be found [here](https://github.com/ktbolt/svFSIplus/blob/update-README/Code/Source/README.md).


<!--- ====================================================================================================================== --->
<!--- ============================================== Solver Parameter Input XML File  ====================================== --->
<!--- ====================================================================================================================== --->

<h1 id="xml_file"> Solver Parameter Input XML File  </h1>

The svFSI solver reads simulation parametes from a [custom text format](https://github.com/SimVascular/svFSI/blob/master/svFSI_master.inp) file.   Simulation parameters are set using keyword/value pairs separated by a colon. A keyword consists of series of alphanumeric characters separated by spaces.

Keyword/value example
```
Coupled: 1
Min iterations: 1
Tolerance: 1e-6
```

Braces **{}** are used to define simulation parameters with sub-elements. For example the `Domain` keyword is associated with
several other sub-elements of keyword/values pairs.
```
Domain: 0 {
  Equation: fluid
  Density: 1.0
  Viscosity: Constant {Value: 0.04}
  Backflow stabilization coefficient: 0.2
}
```

svFSIplus solver simulation parameters are stored in an XML-format file. The XML file organization and keyword names replicate the original text format using the following conversion rules 
- Spaces in keywords names are replaced by underscores
- An additional value after a colon is replaced by an XML atttibute used to identify the value. For example: `Add equation: FSI` is replaced by `<Add_equation type="FSI" >`.
- Parameters with sub-elements are replaced by an XML element with sub-elements. For exmaple `Add equation: FSI { sub-elements }` is replaced by `<Add_equation type="FSI" > xml-sub-elements </Add_equation>`.

Example: The original text 
```
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
```
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

<!--- ====================================================================================================================== --->
<!--- =============================================  svFSIplus Binaries  =================================================== --->
<!--- ====================================================================================================================== --->

<h1 id="executables"> Pre-built svFSIplus Binaries </h1>

Pre-built svFSIplus binaries can be downloaded from the [SimVascular SimTK downloads page](https://simtk.org/frs/?group_id=188).
Installers are provided for MacOS and Ubuntu operating systems. 

The compilers and other software used to build the binaries (listed on the downloads page) must be installed on the target 
platform in order to use the svFSIplus binary.


<!--- ====================================================================================================================== --->
<!--- ============================================= Building svFSIplus  ==================================================== --->
<!--- ====================================================================================================================== --->

<h1 id="building"> Building svFSIplus from Source </h1>

The svFSIplus source code can be downloaed and built on computers running MacOS and other Linux operating systems. Builds are most 
commonly done on MacOS and Ubuntu platforms.

**The Windows operating system is currently not supported.**


<h2 id="building_build"> Build Process </h2>

Several [software packages](#building_packages) must be installed in order to build and use svFSIplus.

svFSIplus is built on a MacOS or Linux computer from a terminal using the following commands
```
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

It is convienient to update the PATH environment variable with the location of the `svFSIplus*` binary. Another option is to create
an alias for `svFSIplus`. For example: `alias svFsiPlus=/usr/local/SV/bin/svFSIplus`


<h2 id="building_packages"> Required Software Packages </h2>

The following software packages are required to build svFSIplus

- CMake
- C++17 compiler 
- Visualization Toolkit
- Open MPI

Software packages can be installed using package management tools like [homebrew](https://brew.sh) and [apt](https://en.wikipedia.org/wiki/APT_(software)).


<h3 id="building_cmake"> CMake </h3>

svFSIplus is built using [CMake](https://cmake.org). CMake is used to control the software compilation process using the `CMakeLists.txt` files located withn the svFSIplus source. CMake first attempts to discover the software dependencies (e.g. C++ compiler) required to build svFSIplus. It will then generate the native [makefiles](https://en.wikipedia.org/wiki/Make_(software)) used to compile it.


<h3 id="building_vtk"> Visualization Toolkit </h3>

svFSIplus uses the [The Visualization Toolkit (VTK)](https://vtk.org/) to read finite element mesh data stored in VTP and VTU [formats](https://kitware.github.io/vtk-examples/site/VTKFileFormats/). The [VTP](https://vtk.org/doc/nightly/html/classvtkPolyData.html#details) format is used to store 2D finite element meshes for inlet and outlet surfaces. The [VTU](https://vtk.org/doc/nightly/html/classvtkUnstructuredGrid.html#details) format is used to store finite element volume meshes. svFSIplus can also write simulation results to VTK format files. 


<h3 id="building_mpi"> Open MPI </h3>

svFSIplus uses the [Open MPI](https://www.open-mpi.org/) Message Passing Interface (MPI) libraries for parallel computing. 


<h2 id="building_build"> Using Trilinos with svFSIplus</h2>

By defualt svFSIplus uses its custom  linear algebra solvers (e.g. GMRES). It also supports several simple preconditioners. 

svFSIplus can be built with the linear algebra solvers and preconditioners provided by the [Trilinos](https://trilinos.github.io/) scientific software package. Using Trilinos can often reduce the execution time of a simulation.


<h3 id="building_Trilinos"> Building Trilinos </h3>

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


<h3 id="building_compile_with_Trilinos"> Compiling with Trilinos Libraries  </h3>

svFSIplus is compiled with the Trilinos libraries using the CMake `SV_USE_TRILINOS` option
```
cmake ../svFSIplus/ -DSV_USE_TRILINOS:BOOL=ON 
```
CMake will be able to find the Trilinos libraries once they have been built and installed.


<!--- ====================================================================================================================== --->
<!--- ============================================= Running a Simulation  ================================================== --->
<!--- ====================================================================================================================== --->

<h1 id="simulation_run"> Running a Simulation </h1>

The svFSIplus binary is run from the command line with an input solver parameters XML file as an argument

```
svFSIplus solver_params.xml
```

A simulation can be run in parallel on four processors using 

```
mpiexec -np 4 svFSIplus solver_params.xml
```








