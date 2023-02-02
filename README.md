# svFSIplus 

## Table of Contents
[Introduction](#introduction)<br>
[Solver Parameter Input XML File ](#xml_file)<br>
[Pre-built svFSIplus Binaries](#executables)<br>
[Building svFSIplus](#building)<br>

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

<h1 id="building"> Building svFSIplus </h1>

svFSIplus can be downloaed and built on computers running MacOS and various Linux operating systems. The Windows operating system is currently
not supported.

The following software packages are required to build svFSIplus

- CMake
- C++17compiler 
- Visualization Toolkit
- Open MPI

<h2 id="building_cmake"> CMake </h2>

svFSIplus is built using [CMake](https://cmake.org). CMake is used to control the software compilation process using the `CMakeLists.txt` files withn the svFSIplus source. CMake first attempts to discover the software dependencies (e.g. C++ compiler) required to build svFSIplus. It will then generate the native [makefiles](https://en.wikipedia.org/wiki/Make_(software)) used to compile it.


<h2 id="building_vtk"> Visualization Toolkit </h2>

[The Visualization Toolkit (VTK)](https://vtk.org/)


<h2 id="building_mpi"> Open MPI </h2>

[Open MPI](https://www.open-mpi.org/)




svFSIplus is built using CMake. It requires the following software 

The following packages are required to build and use svFSI.








