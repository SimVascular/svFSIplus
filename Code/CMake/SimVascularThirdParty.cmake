#-----------------------------------------------------------------------------
# GKLIB_SVFSI
if(SV_USE_GKLIB_SVFSI)
  set(USE_GKLIB_SVFSI ON)
  simvascular_third_party(gklib_svfsi)
  # require to be built here
  set(GKLIB__SVFSI_LIBRARY ${SV_LIB_THIRDPARTY_GKLIB_SVFSI_NAME})
endif()

# METIS_SVFSI
if(SV_USE_METIS_SVFSI)
  set(USE_METIS_SVFSI ON)
  simvascular_third_party(metis_svfsi)
  # require to be built here
  set(METIS_SVFSI_LIBRARY ${SV_LIB_THIRDPARTY_METIS_SVFSI_NAME})
endif()

#-----------------------------------------------------------------------------
# PARMETIS_SVFSI
if(SV_USE_PARMETIS_SVFSI)
  set(USE_PARMETIS_SVFSI ON)
  simvascular_third_party(parmetis_svfsi)
  # require to be built here 
  set(PARMETIS_SVFSI_LIBRARY ${SV_LIB_THIRDPARTY_PARMETIS_SVFSI_NAME})
endif()

#-----------------------------------------------------------------------------
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

# EXPRTK
simvascular_third_party(exprtk)
set(EXPRTK_LIBRARY ${SV_LIB_THIRDPARTY_EXPRTK_NAME})
