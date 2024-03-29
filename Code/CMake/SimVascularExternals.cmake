#-----------------------------------------------------------------------------
# Setting up default directories for SimVascular Externals
# Note: *These can be changes by the user if they have moved the sv_extern directories*
# This is part of the legacy build system
find_package(Doxygen)
if(DOXYGEN_FOUND)
  file(TO_NATIVE_PATH "${SV_BINARY_DIR}/Doxygen/" SV_DOCS_DIR_WORK)
  set(SV_DOCS_DIR ${SV_DOCS_DIR_WORK} CACHE PATH "Location to place docs")
  configure_file(${SV_SOURCE_DIR}/../Documentation/Doxyfile
    ${SV_BINARY_DIR}/Doxyfile @ONLY)
  add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${SV_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif(DOXYGEN_FOUND)

#-----------------------------------------------------------------------------
# Enable Intel Runtime libs if we need or want them
if(SV_USE_INTEL)
  simvascular_external(INTELRUNTIME SYSTEM_DEFAULT SHARED_LIB)
else()
    unset(INTELRUNTIME_LIBRARIES CACHE)
endif()


#-----------------------------------------------------------------------------
# tkcximage (Legacy)
if(WIN32)
  if(SV_USE_TKCXIMAGE)
    find_library(TKCXIMAGE_DLL tkcximage)
    if(TKCXIMAGE_DLL)
      set(TKCXIMAGE_DLL_LIBRARY ${TKCXIMAGE_DLL})
      get_filename_component(TKCXIMAGE_DLL_PATH ${TKCXIMAGE_DLL} DIRECTORY CACHE)
      set(SV_EXTERNAL_SHARED_LIBS ${SV_EXTERNAL_SHARED_LIBS} "TKCXIMAGE")
    endif()
  endif()
endif()


#-----------------------------------------------------------------------------
# Enable MPI
if(SV_USE_MPI)
  if (NOT SV_USE_DUMMY_MPI)
    simvascular_external(MPI SYSTEM_DEFAULT SHARED_LIB)
    if(MPI_FOUND)
      get_filename_component(MPI_LIBRARY_DIR ${MPI_LIBRARY} PATH)
    endif()
    if(WIN32)
      find_library(MPI_fmpich2_LIBRARY NAMES fmpich2 HINTS ${MPI_LIBRARY_DIR})
      set(MPI_EXTRA_LIBRARY ${MPI_EXTRA_LIBRARY} ${MPI_fmpich2_LIBRARY} ${MPI_CXX_LIBRARIES})
      #message("${MPI_EXTRA_LIBRARY}")
    endif()

    if(SV_DEVELOPER_OUTPUT)
      #getListOfVarsPrefix("MPI" _VARLIST)
      #print_vars(_VARLIST)
    endif()
    if(SV_USE_MSMPI)
      # TODO(jmerkow): Change this.
      set(SV_MPI_DIR "${CMAKE_CURRENT_SOURCE_DIR}/ThirdParty/msmpi/")
      set(SV_MPI_LIB_DIR  "${SV_MPI_DIR}/Lib/x64/")
      set(SV_MPI_INCLUDE_PATH "${SV_MPI_DIR}/Include/;${SV_MPI_DIR}/Include/x64/")
      set(SV_MPI_EXTRA_LIBRARY "")
    else()
      set(SV_MPI_EXTRA_LIBRARY ${MPI_EXTRA_LIBRARY})
    endif()
    include_directories(${SV_MPI_INCLUDE_PATH})
  else()
    set(SV_MPI_EXTRA_LIBRARY lib_extra_simvascular_dummympi)
  endif()
endif()
