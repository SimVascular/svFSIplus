svFSIplus is a C++ implementation of the Fortran [svFSI](https://github.com/SimVascular/svFSI) multi-physics finite element solver designed for computational modeling of the cardiovascular system. 

svFSIplus uses a procedural programming paradigm which it inherited from its Fortran predecessor. This means that the code is implemented as a set of functions that call each other. The solver program thus executes by passing data to functions to carry out a series of computational steps.

# Namespaces
svFSIplus uses namespaces to help organize the functions defined within each file. Namespace names are the same as the file name containing their functions.

# Classes
C++ classes are used to reproduce Fortran modules and user-defined data types. A module is a bit like a C++ class because it can encapsulate both data and procedures. Classes representing Fortan modules are stored in a file using the class name; the A class definition is in an A.h file and its implementation is in a A.cpp file.

A user-defined data type (e.g. mshType) is implemented as a C++ class and is used primarilally to store data. All member data is public. Some have a destroy method used to free memory. These classes are defined in the C++ module classes.

There is no significant class hierarchy.

