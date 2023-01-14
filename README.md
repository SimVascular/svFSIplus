# Introduction
*svFSIplus* is a C++ implementation of the Fortran [svFSI](https://github.com/SimVascular/svFSI) multi-physics finite element solver designed for computational modeling of the cardiovascular system. It represents the *First Stage* in the development of a C++ multi-physics finite element solver and is essentailly a direct line-by-line translation of the [svFSI](https://github.com/SimVascular/svFSI) Fortran code. The *Second Stage* of the solver development will be an entirely new implementation encapsilating the *First Stage* code into C++ classes forming a clean and simple abstaction that can be enhanced by any developer.

The C++ implementation differs from the Fortran implementation in three fundamental ways
1) Custom C++ Array and Vector classes to reproduce Fortran dynamic arrays 
2) XML format file replaces the plain-text input file to specify simulation parameters
3) Direct calls to the VTK API replaces the custom code used to read/write VTK format files  

The following sections describe how the C++ implementation is organized and how it replicates the data structures and flow of control of the Fortran implementation.
