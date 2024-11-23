
# These names are used to create variables used to store the package library name.
#
#   THIRDPARTY_METIS_INTERNAL is used to create the variable 
#   SV_LIB_THIRDPARTY_METIS_INTERNAL_NAME 
#
#   SV_LIB_THIRDPARTY_METIS_INTERNAL_NAME is later used to set the value 
#   for the METIS_INTERNAL_LIBRARY_NAME variable in SimVascularThirdParty.cmake.
# 
set(SV_LIBS
  THIRDPARTY_METIS_INTERNAL
  THIRDPARTY_PARMETIS_INTERNAL
  THIRDPARTY_GKLIB_INTERNAL
  THIRDPARTY_TETGEN
  THIRDPARTY_TINYXML
  THIRDPARTY_ZLIB
  SVFSILS
  SVFSI_CINTERFACE
  SVFSILS_CINTERFACE
)

foreach(lib ${SV_LIBS})
  string(TOLOWER "_SIMVASCULAR_${lib}" SV_LIB_${lib}_NAME)
endforeach()
