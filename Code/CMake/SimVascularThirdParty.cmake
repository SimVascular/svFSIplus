
# The Boolean SV_USE variables are set in SimVascularOptions.cmake.
#
# The SV_LIB_THIRDPARTY_*_INTERNAL_NAME is set in SimVascularInternals.cmake.
#
# The *_internal name must match the directory name in Code/ThirdParty
#
if(SV_USE_INTERNAL_GKLIB)
  set(USE_INTERNAL_GKLIB ON)
  simvascular_third_party(gklib_internal)
endif()

if(SV_USE_INTERNAL_METIS)
  set(USE_INTERNAL_METIS ON)
  simvascular_third_party(metis_internal)
endif()

if(SV_USE_INTERNAL_PARMETIS)
  set(USE_INTERNAL_PARMETIS ON)
  simvascular_third_party(parmetis_internal)
endif()

# TETGEN
if(SV_USE_TETGEN)
  set(USE_TETGEN ON)
  simvascular_third_party(tetgen)
endif()

#-----------------------------------------------------------------------------
# EIGEN
if(SV_USE_EIGEN)
  set(USE_EIGEN ON)
  simvascular_third_party(eigen)
endif()

# TINYXML
simvascular_third_party(tinyxml)



