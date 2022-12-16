#-----------------------------------------------------------------------------
# Setup SimVascular Install options and directories
#-----------------------------------------------------------------------------
if(NOT SV_INSTALL_HOME_DIR)
  set(SV_INSTALL_HOME_DIR ".")
endif()

if(NOT SV_INSTALL_TIMESTAMP_DIR)
  set(SV_INSTALL_TIMESTAMP_DIR "timestamp")
endif()

if(NOT SV_INSTALL_SCRIPT_DIR)
  set(SV_INSTALL_SCRIPT_DIR ".")
endif()

if(NOT SV_INSTALL_RUNTIME_DIR)
  set(SV_INSTALL_RUNTIME_DIR bin)
endif()

if(NOT SV_INSTALL_LIBRARY_DIR)
  set(SV_INSTALL_LIBRARY_DIR lib)
endif()

if(NOT SV_INSTALL_ARCHIVE_DIR)
  set(SV_INSTALL_ARCHIVE_DIR lib)
endif()

if(NOT SV_INSTALL_INCLUDE_DIR)
  set(SV_INSTALL_INCLUDE_DIR include)
endif()

if(NOT SV_INSTALL_DATA_DIR)
  set(SV_INSTALL_DATA_DIR data)
endif()

if(NOT SV_INSTALL_DOC_DIR)
 set(SV_INSTALL_DOC_DIR doc)
endif()

if(NOT SV_INSTALL_DOXYGEN_DIR)
  set(SV_INSTALL_DOXYGEN_DIR ${SV_INSTALL_DOC_DIR}/doxygen)
endif()

#-----------------------------------------------------------------------------
# Third Party install locations
#-----------------------------------------------------------------------------
if(NOT SV_INSTALL_MPI_RUNTIME_DIR)
  set(SV_INSTALL_MPI_RUNTIME_DIR lib)
endif()

if(NOT SV_INSTALL_MPI_LIBRARY_DIR)
  set(SV_INSTALL_MPI_LIBRARY_DIR lib)
endif()

if(NOT SV_INSTALL_MPI_EXE_DIR)
  set(SV_INSTALL_MPI_EXE_DIR ".")
endif()

if(NOT SV_INSTALL_EXTERNAL_EXE_DIR)
  set(SV_INSTALL_EXTERNAL_EXE_DIR ".")
endif()

getListOfVarsPrefix("SV_INSTALL" _VARLIST)
foreach(_var ${_VARLIST})
  string(REPLACE "SV_INSTALL" "SV_INSTALL_NATIVE" _var_native ${_var})
  if(${_var})
    file(TO_NATIVE_PATH ${${_var}} ${_var_native})
    dev_message("${_var_native} ${${_var_native}}")
  else()
    dev_message("NO ${_var}")
  endif()
endforeach()

if(SV_DEVELOPER_OUTPUT)
  getListOfVarsPrefix("SV_INSTALL" _VARLIST)
  list(INSERT _VARLIST 0 CMAKE_INSTALL_PREFIX)
  print_vars(_VARLIST)
endif()


#-----------------------------------------------------------------------------
# Setup Output directories for compiling
#-----------------------------------------------------------------------------
if(NOT DEFINED OUTBIN_DIR OR NOT DEFINED OUTLIB_DIR)
  set(OUTBIN_DIR "${SV_BINARY_DIR}/bin")
  set(OUTLIB_DIR ${SV_BINARY_DIR}/lib)
endif()

if(NOT DEFINED SV_DEVELOPER_SCRIPT_DIR)
  set(SV_DEVELOPER_SCRIPT_DIR "${SV_BINARY_DIR}")
endif()

if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${OUTBIN_DIR}")
endif()

if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${OUTLIB_DIR}")
endif()

if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${OUTLIB_DIR}")
endif()

mark_as_advanced(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  CMAKE_LIBRARY_OUTPUT_DIRECTORY
  CMAKE_ARCHIVE_OUTPUT_DIRECTORY)

if(SV_DEVELOPER_OUTPUT)
  set(_VARLIST OUTBIN_DIR
    OUTLIB_DIR
    SCRIPT_DIR)
  print_vars(_VARLIST)
endif()

