
# These variables are set in SimVascularOptions.cmake.
#
if(SV_USE_INTERNAL_GKLIB)
  set(USE_INTERNAL_GKLIB ON)
  simvascular_third_party(gklib_internal)
  # require to be built here
  set(GKLIB__LIBRARY ${SV_LIB_THIRDPARTY_GKLIB_INTERNAL_NAME})
endif()

if(SV_USE_INTERNAL_METIS)
  set(USE_INTERNAL_METIS ON)
  simvascular_third_party(metis_internal)
  # require to be built here
  set(METIS_INTERNAL_LIBRARY ${SV_LIB_THIRDPARTY_METIS_INTERNAL_NAME})
endif()

if(SV_USE_INTERNAL_PARMETIS)
  set(USE_INTERNAL_PARMETIS ON)
  simvascular_third_party(parmetis_internal)
  # require to be built here 
  set(PARMETIS_INTERNAL_LIBRARY ${SV_LIB_THIRDPARTY_PARMETIS_INTERNAL_NAME})
endif()

# TETGEN
if(SV_USE_TETGEN)
  set(USE_TETGEN ON)
  simvascular_third_party(tetgen)
  # require to be built here 
  set(TETGEN_LIBRARY ${SV_LIB_THIRDPARTY_TETGEN_NAME})
endif()

#-----------------------------------------------------------------------------
# EIGEN
if(SV_USE_EIGEN)
  set(USE_EIGEN ON)
  simvascular_third_party(eigen)
  #find_package(Eigen)
  # require to be built here 
  #set(TETGEN_LIBRARY ${SV_LIB_THIRDPARTY_TETGEN_NAME})
endif()

# TINYXML
simvascular_third_party(tinyxml)
set(TINYXML_LIBRARY ${SV_LIB_THIRDPARTY_TINYXML_NAME})



