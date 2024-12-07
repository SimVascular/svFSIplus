# Initial SV Options
#-----------------------------------------------------------------------------
# Developer flag (Output extra info during configure)
option(SV_DEVELOPER_OUTPUT "This is a developer mode to print extra messages during configure" OFF)

set(SV_BUILD_TYPE "CMAKE" CACHE STRING "Designate CMAKE build" FORCE)
set_property(CACHE SV_BUILD_TYPE PROPERTY STRINGS CMAKE)
mark_as_advanced(SV_BUILD_TYPE)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Distribution
option(SV_ENABLE_DISTRIBUTION "Distribute" OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Enable Testing
option(BUILD_TESTING "Build ${PROJECT_NAME} testing" OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Libs - SHARED or STATIC
option(BUILD_SHARED_LIBS "Build ${PROJECT_NAME} as shared libraries." OFF)

set(SV_LIBRARY_TYPE "STATIC" CACHE STRING "Options are STATIC or SHARED" FORCE)
set_property(CACHE SV_LIBRARY_TYPE PROPERTY STRINGS STATIC SHARED)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# SimVascular Build options
option(SV_SUPPRESS_WARNINGS "Option to suppress all compiler warnings while compiling" ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# General Options
option(SV_USE_MPI "Use MSMPI" ON)
option(SV_USE_MSMPI "Use MSMPI" OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ThirdParty
option(SV_USE_ZLIB "Use ZLib" ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Remaining optional dependencies
#-----------------------------------------------------------------------------
# Enable Intel Runtime libs if we need or want them
option(SV_USE_INTEL "Add Intel Runtime Libraries (these may be needed by some libraries)" OFF)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# All OS
option(SV_USE_NOTIMER "Use notimer" ON)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Solver Build Options (Modules)
option(SV_USE_INTERNAL_EIGEN "Use Eigen headers" ON)
option(SV_USE_INTERNAL_METIS "Use the internal metis library" ON)
option(SV_USE_INTERNAL_GKLIB "Use the internal GKlib library" ON)
option(SV_USE_INTERNAL_PARMETIS "Use the internal parmetis library" ON)
option(SV_USE_TETGEN "Use tetgen library" ON)

#-----------------------------------------------------------------------------
# Externals
set(SV_EXTERNALS_TOPLEVEL_DIR "${CMAKE_BINARY_DIR}/sv_externals" CACHE PATH "Externals toplevel directory")
set(SV_EXTERNALS_SRC_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/src" CACHE PATH "Externals toplevel src dir")
set(SV_EXTERNALS_BLD_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/build" CACHE PATH "Externals toplevel build dir")
set(SV_EXTERNALS_PFX_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/prefix" CACHE PATH "Externals toplevel prefix dir")
set(SV_EXTERNALS_BIN_DIR "${SV_EXTERNALS_TOPLEVEL_DIR}/bin/${SV_COMPILER_DIR}/${SV_COMPILER_VERSION_DIR}/${SV_ARCH_DIR}" CACHE PATH "Externals toplevel bin dir")

set(SV_EXTERNALS_INSTALL_PREFIX "sv_externals" CACHE PATH "Externals toplevel directory")

option(SV_EXTERNALS_USE_TOPLEVEL_DIR "If ON, SV_EXTERNALS_TOPLEVEL_DIR will be used as location for external packages" OFF)


#-----------------------------------------------------------------------------
# The internal linear solver is always on 
#-----------------------------------------------------------------------------
set(USE_LINEAR_SOLVER 1)
set(LINEAR_SOLVER_BUILD_TYPE "Source")

#-----------------------------------------------------------------------------
# WIN32
option(SV_USE_WIN32_REGISTRY "Use Windows registry to obtain certain settings (install mode)" OFF)
mark_as_advanced(SV_USE_WIN32_REGISTRY)
#-----------------------------------------------------------------------------

