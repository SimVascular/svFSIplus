# svFSIplus Implementation 

This document descibes some of the implementation details of svFSIplus C++ code.


## Table of Contents
[Introduction](#introduction)<br>
[Code Organization](#organization)<br>
[Translating Fortran into C++](#translate)<br>
[Simulation Class](#simulation_class)<br>
[Array and Vector Classes](#array_vector_class)<br>
[Solver Parameter Input XML File ](#xml_file)<br>
[Implementation Details](#cpp_programming)<br>

<!--- ====================================================================================================================== --->
<!--- ===================================================== Introduction  ================================================== --->
<!--- ====================================================================================================================== --->

<h1 id="introduction"> Introduction </h1>

svFSIplus is essentailly a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code into C++. 
This provides a simple mapping between the code of the original Fortan and C++ versions and aided in debugging the C++ code. 

The C++ implementation differs from the Fortran implementation in four fundamental ways
1) Custom C++ Array and Vector classes to reproduce Fortran dynamic arrays 
2) XML format file replaces the plain-text input file to specify simulation parameters
3) Direct calls to the VTK API replaces the custom code used to read/write VTK format files
4) Uses 0-based indexing into arrays

The following sections describe how the C++ implementation is organized and how it replicates the data structures and flow of control of the Fortran implementation. Some important details of the C++ implementation will also be discussed.

<!--- ====================================================================================================================== --->
<!--- ===================================================== Code Organization  ============================================= --->
<!--- ====================================================================================================================== --->

<h1 id="organization"> Code Organization </h1>

The C++ implementation attempts to replicate the data structures and flow of control of the Fortran implementation and to maintains its organization. 

Most of the Fortran code is replicated in C++ using the same file and procedure names converted to lower case with underscores to improve readability. 
For example
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

All Fortan procedures located in a particular file will typically have a C++ implementation in a similarly named file. This was done to maintain 
a simple mapping between the locations of the C++ and Fortran code. 

The Fortan svFSI code is implemented using a procedural programming paradigm where data is passed to procedures to carry out a series of computational steps. It makes little or no use of the object orientation features providing by Fortran90. This organization is reproduced in the C++ implementation so there are essentailly no class methods used in the core simulation code.

C++ functions are defined within a `namespace` defined for each Fortran file. For example, the functions in `load_msh.cpp` are defined within the `load_msh` `namespace`. Some `namespaces` contain a `_ns` suffix to prevent conflicts with function names (e.g. `read_files_ns`). 

All simulation data is stored in the [Simulation](#simulation_class) class.

<!--- ====================================================================================================================== --->
<!--- =============================================== Translating Fortran into C++ ========================================= --->
<!--- ====================================================================================================================== --->

<h1 id="translate"> Translating Fortran into C++ </h1>

This section provides some details about how the svFSI Fortran code was translated into C++ code. This will help to convert any new 
Fortran code developed in the Fortan svFSI code not included in svFSIplus. 

<!--- ------------------------------- ---> 
<!--- -------  Variable Names ------- --->  
<!--- ------------------------------- --->  
    
<h2 id="translate_vars"> Variable Names </h2>

svFSIplus is essentailly a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. The original Fortran variable names are typically small, contain no underscores for readability and are often ambiguous. However, the **same varaible names** are used 
in both the C++ and Fortran codes in order to maintain a clear relationship between them.

For example, the following section of Fortran code
```
      p = 0._RKIND
      DO a=1, eNoNq
         p = p + Nq(a)*yl(4,a)
      END DO

      uh = 0._RKIND
      IF (mvMsh) THEN
         DO a=1, eNoNw
            uh(1) = uh(1) + Nw(a)*yl(5,a)
            uh(2) = uh(2) + Nw(a)*yl(6,a)
            uh(3) = uh(3) + Nw(a)*yl(7,a)
         END DO
      END IF
      un = (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2) + (u(3)-uh(3))*nV(3)
      un = (ABS(un) - un) * 0.5_RKIND

      u(:) = u(:) - ub(:)
      ubn  = u(1)*nV(1) + u(2)*nV(2) + u(3)*nV(3)
```
is replaced by the following section of C++ code
```
  double p = 0.0;
  for (int a = 0; a < eNoNq; a++) {
    p = p + Nq(a)*yl(3,a);
  }

  Vector<double> uh(3);

  if (com_mod.mvMsh) {
    for (int a = 0; a < eNoNw; a++) {
      uh(0) = uh(0) + Nw(a)*yl(4,a);
      uh(1) = uh(1) + Nw(a)*yl(5,a);
      uh(2) = uh(2) + Nw(a)*yl(6,a);
    }
  }

  double un = (u(0)-uh(0))*nV(0) + (u(1)-uh(1))*nV(1) + (u(2)-uh(2))*nV(2);
  un = (fabs(un) - un) * 0.50;

  u = u - ub;
  double ubn  = u(0)*nV(0) + u(1)*nV(1) + u(2)*nV(2);

  ```

In this example the Fortran `DO` loops are replaced by C++ `for` loops using C++ 0-based indexing. 
Array indexing is discussed in the [Fortran Dynamic Arrays](#translate_arrays) section below.

<!--- -------------------------------- ---> 
<!--- -------  Fortran Modules ------- --->  
<!--- -------------------------------- --->  

<h2 id="translate_modules"> Fortran Modules </h2>

Modules were introduced in Fortran to modualize a large code by splitting it into separate files containing procedures and 
data specific to a certain application. A module is like a C++ class because it can encapsilate both data and procedures. 
The svFSI Fortran code uses modules primarily to store and access global variables. 

C++ classes are used to implement Fortran modules. Fortran variable names are retained to prevent (or maintain) confusion. 
A C++ module name uses the same Fortan name converted to camel case. For example, several of the Fortan module names and the 
files that implements them are given below with the coressponding C++ class name and implementation files.

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
             
The Fortan `USE` command provides access to all the variables defined in a module. Almost all of the svFSI Fortran procedures have 
a `USE COMMOD` command that provides access to all of the global varaibles (about 90) defined in the `COMMOD` module. For example
```
      SUBROUTINE CONSTRUCT_uSOLID(lM, Ag, Yg, Dg)

      USE COMMOD
      USE ALLFUN
```

**svFSIplus does not use any global variables.**  A C++ module object is passed to each procedure that needs to access its variables. 
For example, in C++ the `ComMod` object `com_mod` is explicitly passed to the `construct_usolid` function. 
```
void construct_usolid(ComMod& com_mod, CepMod& cep_mod, const mshType& lM, const Array<double>& Ag,
    const Array<double>& Yg, const Array<double>& Dg)
```
All C++ modules are stored as member data in the [Simulation Class](#simulation_class).


<!--- -------------------------------- ---> 
<!---      Fortran Dynamic Arrays      --->  
<!--- -------------------------------- --->  

<h2 id="translate_arrays"> Fortran Dynamic Arrays </h2>

Fortran dynamic arrays have been reproduced using custom [Vector](#array_vector_class), [Array](array_vector_class) and
[Array3](array_vector_class) C++ class templates. Note that these user-defined classes will most likey be replaced by a 
more sophisticated matrix package such as `Eigen`.

Fortran dynamic arrays are declared using the `ALLOCATABLE` attribute. For example the `REAL, ALLOCATABLE :: A(:,:)` statement 
declares the two dimentional array of type `REAL` named  `A`. The Fortran `ALLOCATE A(3,10)` statement then dynamically creates 
storage for `A` as a 3x10 array. The `DEALLOCATE A` statement is used to return the memory used by `A`. Allocatable arrays are 
automatically deallocated when going out of scope.

C++ dynamic arrays are declared using the `Array<T>` template, where T is the array data type: double or int. The two dimentional 
array of type double named  `A` is declared and memory allocated using `Array<double> A(3,10);`. Memory is released when `A` goes 
out of scope or is explicitly freed using the `clear()` method.
   
C++ multidimensional arrays are referenced using 0-based indexing and are traversed in column-major order like Fortran. Array 
indexes use paranthesis `A(i,j)` not brackets `A[i][j]` to access array elements. 

For example, the following sections of Fortran code that declare and use dynamocs arrays
```
      INTEGER(KIND=IKIND), ALLOCATABLE :: ptr(:)

      REAL(KIND=RKIND), ALLOCATABLE :: xl(:,:), al(:,:), yl(:,:),
     2   bfl(:,:), lR(:,:), lK(:,:,:)

      ALLOCATE(ptr(eNoN), xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN),
     2   bfl(nsd,eNoN), lR(dof,eNoN), lK(dof*dof,eNoN,eNoN))

      ! Create local copies
      DO a=1, eNoN
         Ac = lM%IEN(a,e)
         ptr(a)   = Ac
         xl(:,a)  = x(:,Ac)
         al(:,a)  = Ag(:,Ac)
         yl(:,a)  = Yg(:,Ac)
         bfl(:,a) = Bf(:,Ac)
      END DO
```
are replaced by the following section of C++ code
```
  Vector<int> ptr(eNoN);
  Array<double> xl(nsd,eNoN), al(tDof,eNoN), yl(tDof,eNoN), bfl(nsd,eNoN), lR(dof,eNoN);
  Array3<double> lK(dof*dof,eNoN,eNoN);
  
  // Create local copies
  for (int a = 0; a < eNoN; a++) {
    int Ac = lM.IEN(a,e);
    ptr(a) = Ac;

    for (int i = 0; i < nsd; i++) {
      xl(i,a) = x(i,Ac);               
      bfl(i,a) = Bf(i,Ac);
    }
      
    for (int i = 0; i < tDof; i++) {
      al(i,a) = Ag(i,Ac);
      yl(i,a) = Yg(i,Ac);
    }
  }
```

Note that the `:` array operator used to copy a column of an array is part of the Fortran language. It was not 
always possible to efficiently (i.e. memory-to-memory copy) and cleanly replace Fortran array operators by C++ Array 
template methods. In the above example the Fortran `:` operator was replaced in C++ by an explicit `for` loop. 


<!--- ====================================================================================================================== --->
<!--- =============================================== Simulation Class  ==================================================== --->
<!--- ====================================================================================================================== --->

<h1 id="simulation_class"> Simulation Class </h1>

The C++ [Simulation](https://github.com/SimVascular/svFSIplus/blob/main/Code/Source/svFSI/Simulation.h) class encapsilates 
all of the objects (Fortran modules) used to store simulation data. It also contains a `Parameters` object used to store 
simulation parameters read in from an XML file.

The `Simulation` class does not contain any methods used in the core simulation code. Like the Fortan svFSI code it is only 
used to pass data to procedures to carry out a series of computational steps.


<!--- ====================================================================================================================== --->
<!--- ========================================== Array and Vector Class Templates  ========================================= --->
<!--- ====================================================================================================================== --->

<h1 id="array_vector_class"> Array and Vector Class Templates </h1>

Fortran dynamic arrays have been reproduced using custom `Vector`, `Array` and `Array3` C++ class templates. Note that these custom class templates will most likey be replaced by a more sophisticated matrix package such as `Eigen`.

The class templates are able to reproduce much of the functionality of Fortran arrays and array intrinsic functions (e.g. sum). The challenge is to  create class methods that are as efficient as the Fortan array operators. Because the operators are part of the Fortran language the compiler can
optimize them as a effient memory-to-memory copies. For example
```
A(:,n) = B(:,n)
A = B
```

The objects created from class templates are not part of the C++ language like arrays (i.e. double A[100]). They have the overhead associated with all C++ objects (construct/destroy). Object copy and assignment operators must also be handled efficiently so that intermediate objects are not created and extra data copys are avoided.

The `Vector`, `Array` and `Array3` class templates have a `data()` method that returns a point to the object's internal memory. This is neeed for
MPI calls that take raw C pointers as arguments. For example
```
MPI_Allreduce(part.data(), tmpI.data(), gnNo, cm_mod::mpint, MPI_MAX, cm.com());
```

The class templates are defined in the [Vector.h](https://github.com/SimVascular/svFSIplus/blob/main/Code/Source/svFSI/Vector.h), [Array.h](https://github.com/SimVascular/svFSIplus/blob/main/Code/Source/svFSI/Array.h) and [Array3.h](https://github.com/SimVascular/svFSIplus/blob/main/Code/Source/svFSI/Array3.h) files.


<!--- -------------------------------- ---> 
<!---  Allocating and Freeing Memoory  --->  
<!--- -------------------------------- --->  

<h2 id="array_vector_class"> Allocating and Freeing Memory </h2>

Objects can be created using its constructor
```
Array<double> A(2,2);
```
or defined and later resized
```
Array<double> A;

A.resize(2,2);
```

Object memory is initialized to 0.

An object's memory is freed using its `clear()` method
```
Array<double> A(2,2)

A.clear();
```
or when it goes out of scope.

<!--- -------------------------------- ---> 
<!---               Indexing           --->  
<!--- -------------------------------- --->  

<h2 id="array_vector_class"> Indexing and Memory Layout </h2>

C++ multidimensional arrays are referenced using 0-based indexing and are traversed in column-major order like Fortran. Array indexes use paranthesis `A(i,j)` not brackets `A[i][j]` to access array elements. This was done to make C++ code look more like Fortran and to simplify the conversion process. 
```
  Vector<double> u(2);
  Array<double> ux(2,2);
  Array3<double> uxx(2,2,2);

  for (int a = 0; a < eNoNw; a++) {
    ud(0) = ud(0) + Nw(a)*(al(0,a)-bfl(0,a));
    ud(1) = ud(1) + Nw(a)*(al(1,a)-bfl(1,a));

    u(0) = u(0) + Nw(a)*yl(0,a);
    u(1) = u(1) + Nw(a)*yl(1,a);

    ux(0,0) = ux(0,0) + Nwx(0,a)*yl(0,a);
    ux(1,0) = ux(1,0) + Nwx(1,a)*yl(0,a);
    ux(0,1) = ux(0,1) + Nwx(0,a)*yl(1,a);
    ux(1,1) = ux(1,1) + Nwx(1,a)*yl(1,a);

    uxx(0,0,0) = uxx(0,0,0) + Nwxx(0,a)*yl(0,a);
    uxx(1,0,1) = uxx(1,0,1) + Nwxx(1,a)*yl(0,a);
    uxx(1,0,0) = uxx(1,0,0) + Nwxx(2,a)*yl(0,a);

    uxx(0,1,0) = uxx(0,1,0) + Nwxx(0,a)*yl(1,a);
    uxx(1,1,1) = uxx(1,1,1) + Nwxx(1,a)*yl(1,a);
    uxx(1,1,0) = uxx(1,1,0) + Nwxx(2,a)*yl(1,a);
  }
```

Indexes can be checked by defining the `_check_enabled` directive within each template include file. An index out of bounds will throw an `std::runtime_error` exception. Note that index checking will substantially slow down a simulation so it should be disabled when not testing.


<!--- -------------------------------- ---> 
<!---             Operators            --->  
<!--- -------------------------------- --->  

<h2 id="array_vector_class"> Operators </h2>

Class templates support most mathematical operators: =,+,-,*,/,+=

Some Fortran array intrinsics (e.g. abs, sqrt) are also supported. 

Example
```
  Array<double> Wr(dof,nNo), Wc(dof,nNo);
  
  Wr = 1.0;
  
  Wr = Wr - 0.5;
  Wr = Wr / abs(Wr);
  Wr = (Wr + abs(Wr)) * 0.5;
  
  Wc = 1.0 / sqrt(Wc);

  W1 = W1 * Wr;
  W2 = W2 * Wc;

```

The Array `*` operator performs an element-by-element multiplicaton, not a matrix multiplicaton. This was done to replicate Fortran.

It is more effienct to use the `+=` operator `A += B` than `A = A + B` which performs a copy.

<!--- -------------------------------- ---> 
<!---     Getting an Array Column      --->  
<!--- -------------------------------- --->  

<h2 id="array_vector_class"> Getting an Array Column </h2>

A lot of Fortran code in svFSI operates on a column of a 2D array. For example
```
CALL FLUID3D_M(vmsStab, fs(1)%eNoN, fs(2)%eNoN, w, ksix,
     2            fs(1)%N(:,g), fs(2)%N(:,g), Nwx, Nqx, Nwxx, al, yl,
     3            bfl, lR, lK)
```
where `fs(1)%N(:,g)` gets the column `g` of the `fs(2)%N` array.


The operation of getting a column of data from an Array object is supported using two different methods
```
Vector<T> col(const int col, const std::array<int,2>& range={-1,-1}) const - Returns a new Vector<T> object containing a copy of the column data.

Vector<T> rcol(const int col) const - Returns a new Vector<T> object containing an address pointing into the Array internal data. Modifying the Vector<T> object's data modifies the orginal Array data.
```

Use the `col` method if the column data is not going to be modified
```
auto N0 = fs[0].N.col(g);
auto N1 = fs[1].N.col(g);
fluid_3d_m(com_mod, vmsStab, fs[0].eNoN, fs[1].eNoN, w, ksix, N0, N1, Nwx, Nqx, Nwxx, al, yl, bfl, lR, lK);
```

Use the `rcol` method if the column data is going to be modified; it might also help to speed up a procedure that 
is called a lot (e.g. in material models).

<!--- -------------------------------- ---> 
<!---    Getting an Array3 Slice       --->  
<!--- -------------------------------- --->  

<h2 id="array_vector_class"> Getting an Array3 Slice </h2>

A lot of Fortran code in svFSI operates on a slice (a 2D sub-array) of a 3D array. For example
```
CALL GNN(fs(1)%eNoN, nsd, fs(1)%Nx(:,:,g), xwl, Nwx, Jac, ksix)
```
where `fs(1)%Nx(:,:,g)` gets a slice (2D sub-array) `g` of the `fs(1)%Nx` array.

The operation of getting a slice of data from an Array3 object is supported using two different methods
```
Array<T> slice(const int slice) const - Returns a new Array<T> object containing a copy of the slice data.
Array<T> rslice(const int slice) const - Return an Array<T> object with data pointing into the Array3 internal data.
```

The `rslice()` method returns an `Array` object whose internal data is a pointer to the internal data of an `Array3` object. 
This was done to reduce the overhead of copying sub-arrays in some sections of the custom linear algebra code. The `Array` 
object will not free its data if it is a reference to the data of a `Array3` object. Use the `rslice` method if the slice data
is going to be modified. It can also speed up code that repeatedly extracts sub-arrarys used in computations but are not modified.


<!--- ====================================================================================================================== --->
<!--- ============================================== Solver Parameter Input XML File  ====================================== --->
<!--- ====================================================================================================================== --->

<h1 id="xml_file"> Solver Parameter Input XML File  </h1>

The Fortan svFSI solver read in simulation parameters in a custom text format. All parameters were read in at startup and stored as 
an array of characters (string) using custom code. Parameter values were then retrieved during various stages of the compuation. 
The string representation was converted when the value of a parameter was needed. An error in a parameter value could therefore not
be detected until later stages of the compuataion (e.g., when the mesh is being distributed over processors).

svFSIplus solver simulation parameters are stored in the Extensible Markup Language (XML) file format. XML is a simple text-based format for representing structured data as element tree. The XML tree starts at a root element and branches from the root to sub-elements. All elements can have sub-elements. An XML tag is a markup construct that begins with < and ends with >. There are three types of tag:
1) start-tag, such as `<section>`
2) end-tag, such as `</section>`
3) empty-element tag, such as `<line-break />`

XML _tags_ represent data structures and contain metadata. An XML _element_ is a logical component that either begins with a start-tag and ends with a matching end-tag or consists only of an empty-element tag. The characters between the start-tag and end-tag, if any, are the element's _content_, and may contain markup, including other elements, which are called child elements. An _attribute_ is a markup construct consisting of a name–value pair that exists within a start-tag or empty-element tag.

Exmaple: 
```
 <svFSIFile> |--------------------------------------------- root element
   <Add_equation type="FSI" >   |-------------------------- start-tag with an attribute named type
   <Coupled> true </Coupled>    |--------------------------- Coupled element with data true
   <Min_iterations> 1 </Min_iterations>
   <Max_iterations> 1 </Max_iterations>
  
   <Tolerance> 1e-6 </Tolerance>

   <Domain id="0" >    |----------------------------- Domain element with several child elements
      <Equation> fluid </Equation>
      <Density> 1.0 </Density>
      <Viscosity model="Constant" >
         <Value> 0.04 </Value>
      </Viscosity>
      <Backflow_stabilization_coefficient> 0.2 </Backflow_stabilization_coefficient>
   </Domain>

 </svFSIFile>
```

The elements in the svFSIplus simulation file are represented by _sections_ of related parameters. Sub-elements are refered to as _sub-sections_. 
The svFSIplus simulation file has four top-level sections
```
1) General
2) Mesh
3) Equation
4) Projection
```

<!--- -------------------------------- ---> 
<!---          Parameters Class        --->  
<!--- -------------------------------- --->  

<h2 id="xml_file"> Parameters class </h2>

The [Parameters](https://github.com/SimVascular/svFSIplus/blob/main/Code/Source/svFSI/Parameters.h) class is used to read and store simulation 
parameters parsed from an XML file using using [tinyxml2](https://github.com/leethomason/tinyxml2). Parameter types are checked as they are read so errors in parameter values are immediately detected.


The `Parameters` class contains objects for each of the top-level sections in the parameters file
```
    GeneralSimulationParameters general_simulation_parameters;
    std::vector<MeshParameters*> mesh_parameters;
    std::vector<EquationParameters*> equation_parameters;
    std::vector<ProjectionParameters*> projection_parameters;
```

Each section is represented as a class containing objects for the parameters defined for that section and objects 
representing any sub-sections. Objects representing parameters are named the same as the name used in the XML file except with a 
lower case first character. Each parameter has a name and a value with as a basic type (bool, double, int, etc.) using the 
`Parameter` template class. 

Example: MeshParameters class parameter objects
```
class MeshParameters : public ParameterLists
{
    Parameter<std::string> name;                       // <Add_mesh name=NAME >

    Parameter<bool> initialize_rcr_from_flow;          // <Initialize_RCR_from_flow> BOOL </Initialize_RCR_from_flow>

    Parameter<std::string> mesh_file_path;             // <Mesh_file_path> FILE_NAME </Mesh_file_path>
    
    Parameter<double> mesh_scale_factor;               // <Mesh_scale_factor> SCALE </Mesh_scale_factor>
};
```

All section classes inherit from the `ParameterLists` class which has methods to set parameter values and store them in
a map for processing (e.g. checking that all required parameters have been set). 

Parameter names and default values are set in each section object constructor using member data. The `ParameterLists::set_parameter()`
sets the name and default value for a parameter, and if a value for it is required to be given in the XML file.

Example: Setting parameter names and values in the MeshParameters constructor
```
MeshParameters::MeshParameters()
{
  bool required = true;

  set_parameter("Domain", 0,  !required, domain_id);
  set_parameter("Domain_file_path", "", !required, domain_file_path);

  set_parameter("Fiber_direction_file_path", {}, !required, fiber_direction_file_paths);

  set_parameter("Mesh_file_path", "", required, mesh_file_path);
  set_parameter("Mesh_scale_factor", 1.0, !required, mesh_scale_factor);
  set_parameter("Prestress_file_path", "", !required, prestress_file_path);

  set_parameter("Initial_displacements_file_path", "", !required, initial_displacements_file_path);
  set_parameter("Initial_pressures_file_path", "", !required, initial_pressures_file_path);
  set_parameter("Initial_velocities_file_path", "", !required, initial_velocities_file_path);

  set_parameter("Set_mesh_as_fibers", false, !required, set_mesh_as_fibers);
  set_parameter("Set_mesh_as_shell", false, !required, set_mesh_as_shell);
}
```

Parameter values are set using the 'set_values()' method which contains calls to tinyxml2 to parse parameter 
values from an XML file. The XML elements within a section are extracted in a while loop. Sub-sections or data 
will need to be checked and processed. The 'ParameterLists::set_parameter_value()' method is used to set the value of a 
parameter from a string. 

Example: Parsing XML and setting parameter values in MeshParameters::set_values()
```
void MeshParameters::set_values(tinyxml2::XMLElement* mesh_elem)
{
  using namespace tinyxml2;
  std::string error_msg = "Unknown " + xml_element_name_ + " XML element '";
  auto item = mesh_elem->FirstChildElement();

  while (item != nullptr) {
    auto name = std::string(item->Value());

    // Add_face sub-section.
    if (name == FaceParameters::xml_element_name_) {         // <Add_face name=NAME>
      auto face_params = new FaceParameters();
      face_params->set_values(item);
      face_parameters.push_back(face_params);

    // There may be multiple 'Fiber_direction' elements so store
    // them as a list of VectorParameter<double>. 
    //
    } else if (name == "Fiber_direction") {                   // <Fiber_direction> (x, y, z)  </Fiber_direction>
      auto value = item->GetText();
      VectorParameter<double> dir("Fiber_direction", {}, false, {});
      dir.set(value);
      fiber_directions.push_back(dir);

    // Just a simple element. 
    } else if (item->GetText() != nullptr) {                  // <Mesh_file_path>, <Mesh_scale_factor>, <Domain>, etc.
      auto value = item->GetText();
      try {
        set_parameter_value(name, value);
      } catch (const std::bad_function_call& exception) {
        throw std::runtime_error(error_msg + name + "'.");
      }
    } else {
      throw std::runtime_error(error_msg + name + "'.");
    }

    item = item->NextSiblingElement();
  }
}
```

Sections that contain simple elements (i.e., no sub-sections or special data processing) can be automatically parsed. 

Example: Automatically parsing XML and setting parameter values in `LinearSolverParameters::set_values(tinyxml2::XMLElement* xml_elem)`
```
    std::string error_msg = "Unknown " + xml_element_name + " XML element '";

  // Get the 'type' from the <LS type=TYPE> element.
  const char* stype;
  auto result = xml_elem->QueryStringAttribute("type", &stype);
  if (stype == nullptr) {
    throw std::runtime_error("No TYPE given in the XML <LStype=TYPE> element.");
  }
  type.set(std::string(stype));

  using std::placeholders::_1;
  using std::placeholders::_2;

  // Create a function pointer 'fptr' to 'LinearSolverParameters::set_parameter_value'.
  std::function<void(const std::string&, const std::string&)> ftpr =
      std::bind( &LinearSolverParameters::set_parameter_value, *this, _1, _2);

  // Parse XML and set parameter values.
  xml_util_set_parameters(ftpr, xml_elem, error_msg);
```

<!--- -------------------------------- ---> 
<!---      Accessing Parameters        --->  
<!--- -------------------------------- --->  

<h2 id="xml_file"> Accessing Parameters </h2>

Parameter values are accessed from the core simulation code using the `Simulation` object's `Parameters` object. 
The `Parameter` template class `()` operator or `value()` method is used to access the parameters's value, the `defined()` 
method is used to check if a parameter's value has been set.

Example: Accessing parameter values
```
auto& general_params = simulation->parameters.general_simulation_parameters'
const int nsd = general_params.number_of_spatial_dimensions();
auto file_path = simulation.parameters.mesh_parameters[0].mesh_file_path();

if (eq_params->variable_wall_properties.defined()) { }
```

<!--- ====================================================================================================================== --->
<!--- ============================================= Implementation Details  ================================================ --->
<!--- ====================================================================================================================== --->

<h1 id="cpp_programming"> Implementation Details </h1>

This section covers some of the C++ implementation details that may be useful to developers adding new capabilities to svFSIplus.





