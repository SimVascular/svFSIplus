
#include_directories(./)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../gklib_svfsi/simvascular_gklib_svfsi)
#include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../metis_svfsi/simvascular_metis_svfsi/METISLib)
#include_directories(${MPI_C_INCLUDE_PATH})

set(CSRCS akwayfm.c ametis.c balancemylink.c comm.c
  csrmatch.c ctrl.c debug.c diffutil.c frename.c gklib.c
  gkmetis.c gkmpi.c graph.c initbalance.c initmsection.c initpart.c
  kmetis.c kwayrefine.c match.c mdiffusion.c mesh.c mmetis.c
  move.c msetup.c node_refine.c ometis.c pspases.c redomylink.c
  remap.c renumber.c rmetis.c selectq.c serial.c stat.c
  timer.c util.c wave.c weird.c wspace.c xyzpart.c
)

if(SV_USE_THIRDPARTY_SHARED_LIBRARIES)
  add_library(${PARMETIS_SVFSI_LIBRARY_NAME} SHARED ${CSRCS})
else()
  add_library(${PARMETIS_SVFSI_LIBRARY_NAME} STATIC ${CSRCS})
endif()

set(third_party_dir ${CMAKE_CURRENT_SOURCE_DIR}/../../..)
target_include_directories(${PARMETIS_SVFSI_LIBRARY_NAME} PUBLIC ./)
target_include_directories(${PARMETIS_SVFSI_LIBRARY_NAME} PUBLIC ${third_party_dir}/gklib_svfsi/simvascular_gklib_svfsi)
target_include_directories(${PARMETIS_SVFSI_LIBRARY_NAME} PUBLIC ${third_party_dir}/metis_svfsi/simvascular_metis_svfsi/METISLib)
target_include_directories(${PARMETIS_SVFSI_LIBRARY_NAME} PUBLIC ${MPI_C_INCLUDE_PATH}) 


target_link_libraries(${PARMETIS_SVFSI_LIBRARY_NAME} ${METIS_SVFSI_LIBRARY_NAME}  ${GKLIB_SVFSI_LIBRARY_NAME})

if(SV_INSTALL_LIBS)
  install(TARGETS ${PARMETIS_SVFSI_LIBRARY_NAME}
  RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT ThirdPartyExecutables
  LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT ThirdPartyLibraries
  ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT ThirdPartyLibraries)
endif()

if(SV_INSTALL_HEADERS)

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/temp)
  set(TEMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/temp)
  set(HDRS)
  #copy and rename files
  file(READ "parmetislib.h" METIS_H_TEMP)

  set(filename "macros.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "defs.h")
  file(COPY ${filename} DESTINATION ${TEMP_DIR})
  file(RENAME ${TEMP_DIR}/${filename}
       ${TEMP_DIR}/metis_${filename} )
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "proto.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "rename.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "stdheaders.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "struct.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "parmetis.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi RENAME parmetis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} parmetis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename ${TEMP_DIR}/parmetislib.h)
  file(WRITE ${filename} ${METIS_H_TEMP})
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/parmetis_svfsi
    COMPONENT ThirdPartyHeaders)

endif()
