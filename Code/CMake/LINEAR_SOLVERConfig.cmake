

set(LINEAR_SOLVER_DEFINITIONS "")
set(LINEAR_SOLVER_NEEDED_LIBS LINEAR_SOLVER)

if(IS64 AND (WIN32 OR CYGWIN))
unset(LINEAR_SOLVER_INCLUDE_DIR)
	set(LINEAR_SOLVER_PATH_PREFIX "LINEAR_SOLVER-2013.08.10/win/x64")
	set(LINEAR_SOLVER_FULL_PATH "${LicensedLibs_Bin_Directory}${LINEAR_SOLVER_PATH_PREFIX}")
	set(LINEAR_SOLVER_LIB_DIR "${LINEAR_SOLVER_FULL_PATH}" CACHE TYPE PATH)
	set(LINEAR_SOLVER_POSSIBLE_INCLUDE_DIR "${LINEAR_SOLVER_FULL_PATH}" CACHE TYPE PATH)
endif()

find_path(LINEAR_SOLVER_INCLUDE_DIR LINEAR_SOLVER.h HINTS ${LINEAR_SOLVER_POSSIBLE_INCLUDE_DIR})

GENLIBS(LINEAR_SOLVER_LIBRARY "${LINEAR_SOLVER_NEEDED_LIBS}" "LINEAR_SOLVER" "${LINEAR_SOLVER_LIB_DIR}")
set(LINEAR_SOLVER_LIBRARY ${LINEAR_SOLVER_LINEAR_SOLVER_LIBRARY})

include_directories(${LINEAR_SOLVER_INCLUDE_DIR})
link_directories(${LINEAR_SOLVER_LIB_DIR})
add_definitions(${LINEAR_SOLVER_DEFINITIONS})


