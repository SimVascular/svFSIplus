
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(${SV_SOURCE_DIR}/ThirdParty/parmetis_svfsi/simvascular_parmetis_svfsi/ParMETISLib)
include_directories(${SV_SOURCE_DIR}/ThirdParty/tetgen/simvascular_tetgen)
include_directories(${SV_SOURCE_DIR}/ThirdParty/tinyxml/simvascular_tinyxml)
include_directories(${MPI_C_INCLUDE_PATH})

# Find Trilinos package if requested
if(SV_USE_TRILINOS)

  # source: https://github.com/trilinos/Trilinos_tutorial/wiki/CMakeFindPackageTrilinosExample

  # It may be required to add a path for Trilinos
  # This can be accomplished by adding the Trilinos
  # -DCMAKE_PREFIX_PATH=~/sfw/Trilinos_build/lib/cmake/Trilinos
  FIND_PACKAGE(Trilinos)

  # If FIND_PACKAGE successfully found your Trilinos install, it will
  # set the Boolean flag Trilinos_FOUND.  The following IF statement
  # fails with a FATAL_ERROR if Trilinos was not found.  If it _was_
  # found, it prints out the values of some Trilinos configuration
  # details.  You may find them useful for building your application
  # that uses Trilinos.
  IF(Trilinos_FOUND)
    MESSAGE("\nFound Trilinos!  Here are the details: ")
    MESSAGE("   Trilinos_DIR = ${Trilinos_DIR}")
    MESSAGE("   Trilinos_VERSION = ${Trilinos_VERSION}")
    MESSAGE("   Trilinos_PACKAGE_LIST = ${Trilinos_PACKAGE_LIST}")
    MESSAGE("   Trilinos_LIBRARIES = ${Trilinos_LIBRARIES}")
    MESSAGE("   Trilinos_INCLUDE_DIRS = ${Trilinos_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_TPL_LIST = ${Trilinos_TPL_LIST}")
    MESSAGE("   Trilinos_TPL_INCLUDE_DIRS = ${Trilinos_TPL_INCLUDE_DIRS}")
    MESSAGE("   Trilinos_TPL_LIBRARIES = ${Trilinos_TPL_LIBRARIES}")
    MESSAGE("   Trilinos_BUILD_SHARED_LIBS = ${Trilinos_BUILD_SHARED_LIBS}")
    MESSAGE("   Trilinos_CXX_COMPILER = ${Trilinos_CXX_COMPILER}")
    MESSAGE("   Trilinos_C_COMPILER = ${Trilinos_C_COMPILER}")
    MESSAGE("   Trilinos_Fortran_COMPILER = ${Trilinos_Fortran_COMPILER}")
    MESSAGE("   Trilinos_CXX_COMPILER_FLAGS = ${Trilinos_CXX_COMPILER_FLAGS}")
    MESSAGE("   Trilinos_C_COMPILER_FLAGS = ${Trilinos_C_COMPILER_FLAGS}")
    MESSAGE("   Trilinos_Fortran_COMPILER_FLAGS = ${Trilinos_Fortran_COMPILER_FLAGS}")
    MESSAGE("   Trilinos_LINKER = ${Trilinos_LINKER}")
    MESSAGE("   Trilinos_EXTRA_LD_FLAGS = ${Trilinos_EXTRA_LD_FLAGS}")
    MESSAGE("   Trilinos_AR = ${Trilinos_AR}")
    MESSAGE("End of Trilinos details\n")

    if(NOT Trilinos_LIBRARIES)
      MESSAGE(WARNING "Trilinos_FOUND is true but Trilinos_LIBRARIES is empty. This suggests that FIND_PACKAGE(Trilinos) failed quietly. Proceed with caution.")
    endif()

    set(USE_TRILINOS 1)
  ELSE()
    MESSAGE(WARNING "Could not find Trilinos. Compiling svFSI without Trilinos.")
  ENDIF()

endif()

# add trilinos flags and defines
if(USE_TRILINOS)
  ADD_DEFINITIONS(-DWITH_TRILINOS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif()

# Build with the PETSc linear algebra package.
#
if(NOT "${SV_PETSC_DIR}" STREQUAL "")

   if(NOT EXISTS ${SV_PETSC_DIR})
      message(ERROR "The PETSc directory ${SV_PETSC_DIR} could not be found.")
   endif()

   if(NOT EXISTS "${SV_PETSC_DIR}/include")
      message(ERROR "The PETSc include directory ${SV_PETSC_DIR}/include could not be found.")
   endif()

   if(NOT EXISTS "${SV_PETSC_DIR}/lib")
      message(ERROR "The PETSc library directory ${SV_PETSC_DIR}/lib could not be found.")
   endif()

  set(PETSC_INCLUDE_DIRS "${SV_PETSC_DIR}/include;${SV_PETSC_DIR}/../include")
  set(PETSC_LIBRARY_DIRS "-L${SV_PETSC_DIR}/lib -lpetsc")
  message("\nBuilding with the PETSc package; the include and libraries directories are: ")
  message("  PETSC_LIBRARY_DIRS = ${PETSC_LIBRARY_DIRS}")
  message("  PETSC_INCLUDE_DIRS = ${PETSC_INCLUDE_DIRS}")

  # Set PETSc include directory.
  include_directories(${PETSC_INCLUDE_DIRS})

  # Set C++ directive to use PETSc in svFSIplus code.
  ADD_DEFINITIONS(-DWITH_PETSC)
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pedantic -std=c99")

  set(USE_PETSC 1)

else()

  set(USE_PETSC 0)

endif()

if(ENABLE_ARRAY_INDEX_CHECKING)
  ADD_DEFINITIONS(-DENABLE_ARRAY_INDEX_CHECKING)
endif()

# svFSI requires LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Including VTK produces warnings when compiling Fortran code.
#
#   -iframe is valid for C/C++/ObjC/ObjC++ but not for Fortran
#
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(lib ${SV_LIB_SVFSI_NAME})

set(CSRCS 
  Array3.h Array3.cpp 
  Array.h Array.cpp
  Tensor4.h Tensor4.cpp
  Vector.h Vector.cpp 

  lapack_defs.h

  DebugMsg.h 
  Parameters.h Parameters.cpp
  Simulation.h Simulation.cpp
  SimulationLogger.h
  VtkData.h VtkData.cpp

  all_fun.h all_fun.cpp
  baf_ini.h baf_ini.cpp
  bf.h bf.cpp
  cep.h cep.cpp
  cep_ion.h cep_ion.cpp
  cmm.h cmm.cpp
  consts.h consts.cpp
  contact.h contact.cpp
  distribute.h distribute.cpp
  eq_assem.h eq_assem.cpp
  fluid.h fluid.cpp
  fsi.h fsi.cpp
  fs.h fs.cpp
  fft.h fft.cpp
  heatf.h heatf.cpp
  heats.h heats.cpp
  initialize.h initialize.cpp
  l_elas.h l_elas.cpp
  lhsa.h lhsa.cpp
  ls.h ls.cpp
  main.cpp
  mat_fun.h mat_fun.cpp
  mat_fun_carray.h mat_fun_carray.cpp
  mat_models.h mat_models.cpp
  mesh.h mesh.cpp
  nn.h nn.cpp
  output.h output.cpp
  load_msh.h load_msh.cpp
  pic.h pic.cpp
  post.h post.cpp
  read_files.h read_files.cpp
  read_msh.h read_msh.cpp
  remesh.h remesh.cpp
  remeshTet.cpp
  set_bc.h set_bc.cpp
  shells.h shells.cpp
  stokes.h stokes.cpp
  sv_struct.h sv_struct.cpp
  txt.h txt.cpp
  utils.h utils.cpp
  ustruct.h ustruct.cpp
  vtk_xml.h vtk_xml.cpp
  vtk_xml_parser.h vtk_xml_parser.cpp

  CepMod.h CepMod.cpp
  CepModAp.h CepModAp.cpp
  CepModBo.h CepModBo.cpp
  CepModFn.h CepModFn.cpp
  CepModTtp.h CepModTtp.cpp
  ChnlMod.h ChnlMod.cpp
  CmMod.h CmMod.cpp
  ComMod.h ComMod.cpp
  Timer.h

  SPLIT.c
)

  # Set PETSc interace code.
  if(USE_PETSC)
     #set(CSRCS ${CSRCS} petsc_linear_solver.c)
  endif()

if(USE_TRILINOS)
  set(CSRCS ${CSRCS} trilinos_linear_solver.cpp)

  # trilinos directories and libraries
  include_directories(${Trilinos_TPL_INCLUDE_DIRS})
  include_directories(${Trilinos_INCLUDE_DIRS})
  link_directories(${Trilinos_LIBRARY_DIRS})
  link_directories(${Trilinos_TPL_LIBRARY_DIRS})
endif()

simvascular_add_executable(${SV_SVFSI_EXE}
  SRCS ${CSRCS}
  DEV_SCRIPT_NAME "mysvfsi"
  INSTALL_SCRIPT_NAME "svfsi"
  INSTALL_COMP CoreExecutables
  INSTALL_DESTINATION ${SV_INSTALL_RUNTIME_DIR})

message("##### TINY  ${TINYXML_LIBRARY_NAME}")
message("##### TET ${TETGEN_LIBRARY_NAME}")

target_link_libraries(${SV_SVFSI_EXE}
  ${GLOBAL_LIBRARIES}
  ${INTELRUNTIME_LIBRARIES}
  ${ZLIB_LIBRARY}
  ${BLAS_LIBRARIES}
  ${LAPACK_LIBRARIES}
  ${METIS_SVFSI_LIBRARY_NAME}
  ${PARMETIS_SVFSI_LIBRARY_NAME}
  ${TETGEN_LIBRARY_NAME}
  ${TINYXML_LIBRARY_NAME}
  ${SV_LIB_SVFSILS_NAME}${SV_MPI_NAME_EXT}
  ${VTK_LIBRARIES}
  )

# extra MPI libraries only if there are not set to NOT_FOUND or other null
if(SV_MPI_EXTRA_LIBRARY)
  target_link_libraries(${SV_SVFSI_EXE} ${SV_MPI_EXTRA_LIBRARY})
endif()

if(USE_TRILINOS)
  target_link_libraries(${SV_SVFSI_EXE} ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES})
endif()

if(USE_PETSC)
  target_link_libraries(${SV_SVFSI_EXE} ${PETSC_LIBRARY_DIRS})
endif()

# coverage
if(ENABLE_COVERAGE)
  # set compiler flags
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -coverage")

  # find required tools
  find_program(LCOV lcov REQUIRED)
  find_program(GENHTML genhtml REQUIRED)

  # add coverage target
  add_custom_target(coverage
    # gather data
    COMMAND ${LCOV} --directory . --capture --output-file coverage.info
    # exclude externals
    COMMAND ${LCOV} --remove coverage.info -o coverage.info '/usr/*' '/opt/*' '${PROJECT_SOURCE_DIR}/ThirdParty/*' '/Library/*' 'v1/*'
    # generate report
    COMMAND ${GENHTML} --demangle-cpp -o coverage coverage.info
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
endif()

# unit tests and Google Test
if(ENABLE_UNIT_TEST)
  # install Google Test
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
  )
  FetchContent_MakeAvailable(googletest)
  enable_testing()
  include(GoogleTest)

  # add hello_test, a sample test file to check the installation of Google Test
  add_executable(
    hello_test
    ../../../tests/unitTests/hello_test.cpp
  )
  target_link_libraries(
    hello_test
    GTest::gtest_main
  )
  gtest_discover_tests(hello_test)


  # add test.cpp for unit test

  # remove the main.cpp and add test.cpp
  list(APPEND CSRCS "../../../tests/unitTests/test.cpp")
  list(REMOVE_ITEM CSRCS "main.cpp")
  
  # include source files (same as what svFSI does except for main.cpp)
  add_executable(runUnitTest ${CSRCS})
  
  # libraries
  target_link_libraries(runUnitTest
    ${METIS_SVFSI_LIBRARY_NAME}
    ${PARMETIS_SVFSI_LIBRARY_NAME}
    ${TETGEN_LIBRARY_NAME}
    ${TINYXML_LIBRARY_NAME}
    ${SV_LIB_SVFSILS_NAME}${SV_MPI_NAME_EXT}
    ${VTK_LIBRARIES}
  )

  # link Google Test
  target_link_libraries(
    runUnitTest
    GTest::gtest_main
  )
  gtest_discover_tests(runUnitTest)

endif()

