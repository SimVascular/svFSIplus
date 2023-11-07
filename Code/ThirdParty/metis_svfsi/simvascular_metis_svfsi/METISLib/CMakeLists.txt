include_directories(./)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../gklib_svfsi/simvascular_gklib_svfsi)
include_directories(${MPI_C_INCLUDE_PATH})

set(CSRCS auxapi.c balance.c bucketsort.c checkgraph.c coarsen.c
  compress.c contig.c debug.c fm.c fortran.c frename.c gklib.c
  graph.c initpart.c kmetis.c kwayfm.c kwayrefine.c mcutil.c
  mesh.c meshpart.c minconn.c mincover.c mmd.c ometis.c
  options.c parmetis.c pmetis.c refine.c separator.c sfm.c
  srefine.c stat.c timing.c util.c wspace.c
)

if(SV_USE_THIRDPARTY_SHARED_LIBRARIES)
  add_library(${METIS_SVFSI_LIBRARY_NAME} SHARED ${CSRCS})
else()
  add_library(${METIS_SVFSI_LIBRARY_NAME} STATIC ${CSRCS})
endif()

if(SV_INSTALL_LIBS)
  install(TARGETS ${METIS_SVFSI_LIBRARY_NAME}
  RUNTIME DESTINATION ${SV_INSTALL_RUNTIME_DIR} COMPONENT ThirdPartyExecutables
  LIBRARY DESTINATION ${SV_INSTALL_LIBRARY_DIR} COMPONENT ThirdPartyLibraries
  ARCHIVE DESTINATION ${SV_INSTALL_ARCHIVE_DIR} COMPONENT ThirdPartyLibraries)
endif()

if(SV_INSTALL_HEADERS)

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/temp)
  set(TEMP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/temp)
  set(HDRS)
  #copy and rename files
  file(READ "metis.h" METIS_H_TEMP)

  set(filename "macros.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "defs.h")
  file(COPY ${filename} DESTINATION ${TEMP_DIR})
  file(RENAME ${TEMP_DIR}/${filename}
       ${TEMP_DIR}/metis_${filename} )
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "proto.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "rename.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "stdheaders.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "struct.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename "parmetis.h")
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi RENAME metis_svfsi_${filename}
    COMPONENT ThirdPartyHeaders)
  string(REGEX REPLACE ${filename} metis_svfsi_${filename} METIS_H_TEMP ${METIS_H_TEMP})

  set(filename ${TEMP_DIR}/metis.h)
  file(WRITE ${filename} ${METIS_H_TEMP})
  install(FILES ${filename} DESTINATION ${SV_INSTALL_INCLUDE_DIR}/thirdparty/metis_svfsi
    COMPONENT ThirdPartyHeaders)
endif()


