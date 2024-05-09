# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/dcodoni/research/svFSIplus-package/svFSIplus/Code"
  "/Users/dcodoni/research/svFSIplus-package/svFSIplus/build-petsc/svFSI-build"
  "/Users/dcodoni/research/svFSIplus-package/svFSIplus/build-petsc/svFSI-prefix"
  "/Users/dcodoni/research/svFSIplus-package/svFSIplus/build-petsc/svFSI-prefix/tmp"
  "/Users/dcodoni/research/svFSIplus-package/svFSIplus/build-petsc/svFSI-prefix/src/svFSI-stamp"
  "/Users/dcodoni/research/svFSIplus-package/svFSIplus/build-petsc/svFSI-prefix/src"
  "/Users/dcodoni/research/svFSIplus-package/svFSIplus/build-petsc/svFSI-prefix/src/svFSI-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/dcodoni/research/svFSIplus-package/svFSIplus/build-petsc/svFSI-prefix/src/svFSI-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/dcodoni/research/svFSIplus-package/svFSIplus/build-petsc/svFSI-prefix/src/svFSI-stamp${cfgdir}") # cfgdir has leading slash
endif()
