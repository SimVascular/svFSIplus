# svFSIplus

[![Build Status](https://github.com/SimVascular/svFSIplus/actions/workflows/test.yml/badge.svg)](https://github.com/SimVascular/svFSIplus/actions)
[![codecov](https://codecov.io/github/SimVascular/svFSIplus/graph/badge.svg?token=I848DNIHSP)](https://codecov.io/github/SimVascular/svFSIplus)
![Latest Release](https://img.shields.io/github/v/release/SimVascular/svFSIplus?label=latest)
![Platform](https://img.shields.io/badge/platform-macOS%20|%20linux-blue)
[![DOI](https://img.shields.io/badge/DOI-10.21105%2Fjoss.04118-green)](https://doi.org/10.21105/joss.04118)

### Table of Contents
**[Introduction](#introduction)**<br>
**[Building the svFSIplus Program from Source](#building)**<br>
**[Building with External Linear Algebra Packages](#linear_algebra_packages)**<br>
**[Running the svFSIplus Program](#running_svfsiplus)**<br> 
**[Docker Container](#docker_container)**<br>
**[Testing](tests/README.md)**<br>


<!--- =================================================================================== -->
<!---                                 Introduction                                        -->
<!--- =================================================================================== -->

<h1 id="introduction"> Introduction </h1>

svFSIplus is an open-source, parallel, finite element multi-physics solver providing capabilities to simulate the partial differential equations (PDEs) governing solid and fluid mechanics, diffusion, and electrophysiology. Equations can be solved coupled to simulate the interaction between multiple regions representing different physical systems. For example, in a coupled fluid-solid simulation the motion of the fluid can deform a solid region while the changing geometry of the solid region changes the way the fluid flows. Equation coupling provides a framework for the computational modeling of whole heart dynamics.

svFSIplus is a C++ implementation of the Fortran [svFSI](https://github.com/SimVascular/svFSI) multi-physics finite element solver designed for computational modeling of the cardiovascular system. The C++ implementation is essentially a line-by-line translation of the svFSI Fortran code and therefore uses a procedural rather than an object oriented programming paradigm. The code will be incrementally refactored into an object oriented code. 

The [SimVascular svFSIplus Documentation](https://simvascular.github.io/documentation/svfsiplus.html) provides documentation describing how to use the svFSIplus solver. It also has developer guide describing the code organization and some implementation details. 

The [svFSIplus Internal Code Documentation](https://simvascular.github.io/svfsiplus/index.html) provides documentation of the svFSIplus source code. It is automatically generated using [Doxygen](https://www.doxygen.nl).

<!--- =================================================================================== -->
<!---                       Building the svFSIplus                                        -->
<!--- =================================================================================== -->

<h1 id="building"> Building the svFSIplus Program from Source </h1>
The svFSIplus program can be compiled and linked from the GitHub source using a CMake build process. The build process creates a binary executable file named **svfsiplus**.

## Supported Platforms

svFSIplus can be built on most Unix-like operating systems
- MacOS 
- Ubuntu
- CentOS

A native Windows build is currently not supported. However, svFSIplus can be built and run from an Ubuntu terminal environment on Windows using [Windows Subsystem for Linux (WSL)](https://ubuntu.com/desktop/wsl).

## Software Dependencies
The following software packages are required to be installed in order to build svFSIplus

- [git](https://git-scm.com/) - Distributed version control system used to interact with GitHub repositories
- [CMake](https://cmake.org/) - Used to the build the binary executable **svfsiplus**
- C++17 compiler - C++ compiler, linker and libraries 
- [Visualization Toolkit (VTK)](https://vtk.org/) - Used for reading and writing VTK-format VTP and VTU files 
- [Open MPI](https://www.open-mpi.org/) - Used for parallel processing
- [BLAS](https://www.netlib.org/blas/) - Used for performing basic vector and matrix operations (optional but may be needed for external linear algebra packages)
- [LAPACK](https://www.netlib.org/lapack/) - Used for solving systems of simultaneous linear equations (optional but may be needed for external linear algebra packages)



These software packages are installed using a package-management system
- Ubuntu - [apt](https://en.wikipedia.org/wiki/APT_(software))
- MacOS - [homebrew](https://en.wikipedia.org/wiki/Homebrew_(package_manager))
- high-performance computing (HPC) cluster - [module system](https://hpc-wiki.info/hpc/Modules#:~:text=From%20HPC%20Wiki,average%20user%20will%20ever%20use)


Installing VTK on a high-performance computing (HPC) cluster is typically not supported and may require building it from source. See [Building Visualization Toolkit (VTK) Libraries](#building_vtk).


## Building svFSIplus

svFSIplus is built using the following steps

1) Download the source from GitHub

   `git clone https://github.com/SimVascular/svFSIplus`

   Creates an `svFSIplus` directory.
   
3) Create a build directory and change directories to it

   ```
   mkdir fsiplus-build
   cd fsiplus-build
   ```
   
4) Execute the build

   ```
   cmake ../svFSIplus
   make
   ```
   This creates the **svfsiplus** binary executable located in
   ```
   fsiplus-build/svFSIplus-build/bin/svfsiplus
   ```


<h2 id="building_vtk"> Building Visualization Toolkit (VTK) Libraries </h2>

svFSIplus uses VTK to read finite element mesh data (created by the SimVascular mesh generation software), fiber geometry, initial conditions and write simulation results. Building the complete VTK library requires certain graphics libraries to be installed (e.g. OpenGL, X11) which make it difficult to build on an HPC cluster. 

However, a subset of the complete VTK library can be built to just include reading/writing functionality without graphics.

The following steps demonstrate how to build local VTK libraries in the `/user/shared/vtk` directory from the source downloaded from https://www.vtk.org/files/release/9.3/VTK-9.3.1.tar.gz. 


```
mkdir /user/shared/vtk
cd /user/shared/vtk
wget https://www.vtk.org/files/release/9.3/VTK-9.3.1.tar.gz
tar xvf VTK-9.3.1.tar.gz 
mkdir build
cd build
cmake -DBUILD_SHARED_LIBS:BOOL=OFF \
      -DCMAKE_BUILD_TYPE:STRING=RELEASE \
      -DBUILD_EXAMPLES=OFF \
      -DBUILD_TESTING=OFF \
      -DVTK_USE_SYSTEM_EXPAT:BOOL=ON \
      -DVTK_USE_SYSTEM_ZLIB:BOOL=ON \
      -DVTK_LEGACY_REMOVE=ON \
      -DVTK_Group_Rendering=OFF \
      -DVTK_Group_StandAlone=OFF \
      -DVTK_RENDERING_BACKEND=None \
      -DVTK_WRAP_PYTHON=OFF \
      -DModule_vtkChartsCore=ON \
      -DModule_vtkCommonCore=ON \
      -DModule_vtkCommonDataModel=ON \
      -DModule_vtkCommonExecutionModel=ON \
      -DModule_vtkFiltersCore=ON \
      -DModule_vtkFiltersFlowPaths=ON \
      -DModule_vtkFiltersModeling=ON \
      -DModule_vtkIOLegacy=ON \
      -DModule_vtkIOXML=ON \
      -DVTK_GROUP_ENABLE_Views=NO \
      -DVTK_GROUP_ENABLE_Web=NO \
      -DVTK_GROUP_ENABLE_Imaging=NO \
      -DVTK_GROUP_ENABLE_Qt=DONT_WANT \
      -DVTK_GROUP_ENABLE_Rendering=DONT_WANT \
      -DCMAKE_INSTALL_PREFIX=/user/shared/vtk/install \
      ../VTK-9.3.1
make -j4
make install
```

These commands will create the following directories under /user/shared/vtk/install
```
bin/ include/ lib/ share/
```

You can then build svFISplus with the local VTK libraries using the two CMake command-line arguments that sets the option to use a local VTK build and the location of the build
```
-DSV_USE_LOCAL_VTK=ON  -DSV_VTK_LOCAL_PATH=/user/shared/vtk/install 
```

<!--- =================================================================================== -->
<!---                      External Linear Algebra Packages                               -->
<!--- =================================================================================== -->

<h1 id="linear_algebra_packages"> Building with External Linear Algebra Packages </h1>

Numerical linear algebra uses computer algorithms to solve the linear system generated by finite element method. Linear algebra libraries provide access to specialized or 
general purpose routines implementing a significant number of computer algorithms useful when solving linear systems.

svFSIplus supports interfaces to the following numerical linear algebra packages
- Custom numerical linear algebra routines included in the svFSIplus implementation
- [Portable, Extensible Toolkit for Scientific Computation (PETSc)](https://petsc.org/release/overview/)
- [Trilinos](https://trilinos.github.io/)

These packages may be not be available for the latest version using a package-management system or be available on an HPC cluster so they must be built from source.

## Building svFSIplus with External Linear Algebra Packages

The following CMake command-line arguments are used to build svFSIplus with external linear algebra packages
```
-DSV_USE_TRILINOS:BOOL=ON - Sets an option for CMake to look for the Trilinos package
-DSV_PETSC_DIR:PathToPetsInstall - Sets the location where the PETSc package is installed
```

When built with an external linear algebra package svFSIplus can be run using that package by setting the **Linear_algebra** parameter in the solver input file.
For example: PETSc
```
<Linear_algebra type="petsc" >
  <Preconditioner> petsc-rcs </Preconditioner>
</Linear_algebra>
```
For example: Trilinos
```
<Linear_algebra type="trilinos" >
  <Preconditioner> trilinos-ilut </Preconditioner>
</Linear_algebra>
```


## Building Trilinos 

To build Trilinos first download the source 
```
git clone https://github.com/trilinos/Trilinos
```

Then follow the [Trilinos build instructions](https://github.com/trilinos/Trilinos/blob/master/INSTALL.rst).

When building select the following third-party libraries (TPLs)
```
Boost
BLAS
HDF5
HYPRE
LAPACK
MUMPS
```
and the following Trilinos packages
```
mesos
AztecOO
Epetra
EpetraEXT
Ifpack
ML
MueLU
ROL
Sacado
Teuchos
Zoltan
```

## Building PETSc

PETSc libraries can be [installed](https://petsc.org/release/install/) using package managers.

See the [Quick Start Tutorial](https://petsc.org/release/install/install_tutorial/) for building from source.


<!--- =================================================================================== -->
<!---                      Running the svFSIplus Program                                  -->
<!--- =================================================================================== -->

<h1 id="running_svfsiplus"> Running the svFSIplus Program  </h1>

Once the **svfsiplus** binary is available add its location to your environment variable (`PATH` on Unix-like systems).
The **svfsiplus** binary is run from the command line with the name of a input solver parameters XML file as an argument

```
svfsiplus fluid3d.xml
```

The simulation progress for each time step will be printed showing
various solver convergence measures
```
---------------------------------------------------------------------
 Eq     N-i     T       dB  Ri/R1   Ri/R0    R/Ri     lsIt   dB  %t
---------------------------------------------------------------------
 NS 1-1  2.790e-01  [0 1.000e+00 1.000e+00 6.522e-12]  [157 -258 82]
 NS 1-2  8.890e-01  [-57 1.373e-03 1.373e-03 3.007e-11]  [253 -99 92]
 NS 1-3  1.067e+00  [-125 5.082e-07 5.082e-07 1.680e-05]  [117 -110 74]
 NS 1-4  1.114e+00  [-197 1.304e-10 1.304e-10 6.799e-02]  [7 -27 2]
 NS 1-5  1.170e+00  [-221 8.863e-12 8.863e-12 1.000e+00]  !0 0 0!
 NS 2-1  2.341e+00  [0 1.000e+00 2.856e+01 7.250e-14]  [504 -124 96]
 NS 2-2  2.580e+00  [-75 1.586e-04 4.529e-03 2.002e-09]  [143 -200 81]
 NS 2-3  2.771e+00  [-146 4.871e-08 1.391e-06 6.319e-06]  [123 -120 75]
 NS 2-4  2.820e+00  [-216 1.483e-11 4.235e-10 1.781e-02]  [11 -40 6]
 NS 2-5s 2.894e+00  [-251 2.642e-13 7.547e-12 1.000e+00]  !0 0 0!
```
A directory named `1-procs` containing the simulation results output will be created
```
B_NS_Pressure_average.txt       histor.dat                      stFile_last.bin
B_NS_Velocity_flux.txt          result_002.vtu
B_NS_WSS_average.txt            stFile_002.bin
```

A simulation can be run in parallel on four processors using
```
mpiexec -np 4 svfsiplus fluid3.xml
```
In this case a directory named `4-procs` containing the simulation results output will be created. Results from different processors will be combined into a single file for a given time step.

<!--- =================================================================================== -->
<!---                             Docker container                                        -->
<!--- =================================================================================== -->

<h1 id="docker_container">  Docker </h1>
An alternative way with respect to building from source, is the Docker option. To use this option Docker must be installed first. Please refer to [Docker webpage](https://www.docker.com/products/docker-desktop/) to know more about Docker and how to install it on your machine. The following steps describe how to build a Docker image or pull an existent one from DockerHub, and how to run a Docker container. The last section is a brief guide to perform the same steps but in Singularity, since HPC systems usually use Singularity to handle containers.

## Docker image
A Docker image is a read-only template that may contain dependencies, libraries, and everything needed to run a program. It is like a snapshot of a particular environment. 
A Docker image can be created directly from a [dockerfile](https://docs.docker.com/reference/dockerfile/#:~:text=A%20Dockerfile%20is%20a%20text,can%20use%20in%20a%20Dockerfile.) or an existent image can be pulled from [DockerHub](https://hub.docker.com), if available. For this repository, both options are available.
The latest version of svFSIplus program is pre-compiled in a Docker image, built from a dockerfile provided in Docker/solver. The Docker image includes two different type of builds, one where the solver is compiled with Trilinos and the other one where the solver is compiled with PETSc. 
This Docker image can be downloaded (pulled) from the dockerhub simvascular repository [simvascular/solver](https://registry.hub.docker.com/u/simvascular). To pull an image, run the command:
```
docker pull simvascular/solver:latest
```
Note that this image was built for AMD64 (x86) architecture, and it will not work on other architectures such as ARM64 (AArch64). In this case, the image has to be built from the provided dockerfile. 
To create the image from the dockerfile provided in Docker/solver, follow the steps below:
1) build an Ubuntu-based image containing the whole environment in which svFSIplus program can be compiled. The provided dockerfiles are based on Ubuntu-20.04 and Ubuntu-22.04, but they can be easily adapted to use the latest version of Ubuntu, by changing the following line in Docker/ubuntu20/dockerfile or Docker/ubuntu22/dockerfile: 
```
FROM ubuntu:20.04 AS base / FROM ubuntu:22.04 AS base
```
to 
```
FROM ubuntu:latest AS base
```
Build the environmnet Docker image: 
```
cd Docker/ubuntu20 or cd Docker/ubuntu22
```
```
docker build -t RepositoryName:tagImage .
```
where -t allows the user to set the name for the image created. For example:
```
docker build -t libraries:latest .
```
2) build the image containing the compiled svFSIplus program. This image will be based on the environment created in the previous step (libraries:latest). In order to do this, open the Docker/solver/dockerfile and modify the following lines:
```
FROM simvascular/libraries:ubuntu22 AS builder 
```
to
```
FROM libraries:latest AS builder
```
and 
```
FROM simvascular/libraries:ubuntu22 AS final 
```
to
```
FROM libraries:latest AS final
```
Build the solver image:
```
cd Docker/solver
```
```
docker build -t solver:latest .
```
The image include the PETSc-based svFSIplus executable in:
```
/build-petsc/svFSIplus-build/bin/svfsiplus
```
and the Trilinos-based svFSIplus executable in:
```
/build-trilinos/svFSIplus-build/bin/svfsiplus
```

## Docker container
A Docker container is a running instance of a Docker image. It is a lightweight, isolated, and executable unit. 
Once the image is created, it can be run interactively by running the following command:
```
docker run -it -v FolderToUpload:/NameOfFolder solver:latest
```
In this command:
- -it: means run interactively Docker image
- -v: mounts a directory 'FolderToUpload' from the host machine in the container where the directory has the name '/NameOfFolder'. For example the folder containing the mesh and the input file necessary to run a simulation should be mounted.
A command can also be run in the container directly without running the container interactively. For example we want to run the solver using the mpirun, we may want to run the following command:
```
docker run -v FolderToUpload:/NameOfFolder solver:latest mpirun -n 4 /build-trilinos/svFSIplus-build/bin/svfsiplus svFSIplus.xml
```
The previous command will run the solver on 4 processors using the input file svFSIplus.xml and the mesh in the folder 'FolderToUpload' mounted inside the container.

## Containers on HPC: Singularity
Most of the HPC systems (if not all) are based on AMD64 architecture and the solver image can be directly pulled from [simvascular/solver](https://hub.docker.com/r/simvascular/solver). First of all, make sure the singularity module is loaded on the HPC system. Then, pull the solver image (it is recommended to run the following command on the compute node for example through an interactive job):
```
singularity pull docker://simvascular/solver:latest
``` 
After the pull is complete, you should have a file with extension .sif (solver image). This image contains the two executables of the svFSIplus program build with PETSc and Trilinos support, respectively. 
In the following, we provide two example of job submission's scripts that can be used as a reference to run a simulation using the svFSIplus solver on an HPC cluster. 
1) single-node job script:
```
#!/bin/bash
#SBATCH --job-name
#SBATCH --output
#SBATCH --partition
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=
#SBATCH --mem=0 
#SBATCH -t 48:00:00

# For single node, no modules should be loaded to avoid incongruences between HPC and containers environments
module purge

singularity run --bind #mount-the-necessary-folders-the-container-should-have-access-to \
/PathToTheImage/NameOfImagePulled.sif \
mpirun -n #TotalNumberOfTasks /build-trilinos/svFSIplus-build/bin/svfsiplus svFSIplus.xml
```

2) multi-node job script
```
#!/bin/bash
#SBATCH --job-name
#SBATCH --output
#SBATCH --partition
#SBATCH --nodes
#SBATCH --ntasks-per-node
#SBATCH --mem
#SBATCH -t 00:00:00

# The following 'export' may not work on all the HPC systems
export UCX_TLS=ib
export PMIX_MCA_gds=hash
export OMPI_MCA_btl_tcp_if_include=ib0

module purge
# Load here all the modules necessary to use the HPC MPI, for example: module load openmpi

mpirun -n #TotalNumberOfTasks singularity run --bind #mount-the-necessary-folders-the-container-should-have-access-to \
/PathToTheImage/NameOfImagePulled.sif \
/build-trilinos/svFSIplus-build/bin/svfsiplus svFSIplus.xml
```
Since the multi-node relies on both MPI, the one on the HPC and the one inside the container, there may be some problems. In the following, we give a solution (workaround) for two common problems:
- if the HPC OpenMPI was built with cuda support, then it may happen that it is expecting that OpenMPI inside the container to be built with cuda support too, which is not the case. Possible solution is to add --mca mpi_cuda_support 0:
```
mpirun --mca mpi_cuda_support 0 -n #TotalNumberOfTasks ... 
```
- if for some reason, it is complaining about not finding 'munge' then add --mca psec ^munge:
```
mpirun --mca psec ^munge -n #TotalNumberOfTasks ... 
```