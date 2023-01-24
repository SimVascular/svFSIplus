# svFSIplus 

## Table of Contents
[Introduction](#introduction)<br>
[Solver Parameter Input XML File ](#xml_file)<br>
[Implementation Details](#cpp_programming)<br>

<!--- ====================================================================================================================== --->
<!--- ===================================================== Introduction  ================================================== --->
<!--- ====================================================================================================================== --->

<h1 id="introduction"> Introduction </h1>

svFSIplus is a C++ implementation of the Fortran [svFSI](https://github.com/SimVascular/svFSI) multi-physics finite element solver designed for computational modeling of the cardiovascular system. It represents the <i>First Stage</i> in the development of a C++ multi-physics finite element solver and is essentailly a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. This was done to maintain a simple mapping between the code of the two versions and to faciliate debugging by comparing intermediate data (i.e. element stiffness matrices) and simulation results.

The *Second Stage* of the solver development will be an entirely new implementation encapsilating portions of the *First Stage* code into C++ classes forming a clean and simple abstaction that can be enhanced by any developer.

The following sections describe how to build and use svFSIplus.


<!--- ====================================================================================================================== --->
<!--- ============================================== Solver Parameter Input XML File  ====================================== --->
<!--- ====================================================================================================================== --->

<h1 id="xml_file"> Solver Parameter Input XML File  </h1>

The svFSIplus solver parameters read in for a simulation are stored in an XML-format file. The XML file organization and parameter names replicate the old input text file except that parameter names have spaces replaced by underscores. 

Parameters using an additional value after the `:` have an XML atttibute added to identify the value: `Add equation: FSI` is converted to `<Add_equation type="FSI" >`.

Example: Old input text file
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
new XML format
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
<!--- ============================================= Implementation Details  ================================================ --->
<!--- ====================================================================================================================== --->

<h1 id="cpp_programming"> Implementation Details </h1>

This section covers some of the C++ implementation details that may be useful to developers adding new capabilities to svFSIplus.




