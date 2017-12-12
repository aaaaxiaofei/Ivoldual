PROJECT(IVOLDUAL)

#---------------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

SET(NRRD_LIBDIR "~/programs/ijk/lib")

#---------------------------------------------------------

IF (NOT CMAKE_BUILD_TYPE)
  SET (CMAKE_BUILD_TYPE Release CACHE STRING 
       "Default build type: Release" FORCE)
ENDIF (NOT CMAKE_BUILD_TYPE)

LINK_DIRECTORIES("${NRRD_LIBDIR}")
LINK_LIBRARIES(NrrdIO z)

# Include random library.
SET(CMAKE_CXX_FLAGS "-std=c++11")

ADD_EXECUTABLE(ivoldual ivoldual_main.cxx ivoldualIO.cxx isodual.cxx
                        ivoldual.cxx ijkdual_datastruct.cxx
                        ivoldual_datastruct.cxx ivoldual_triangulate.cxx 
                        ivoldualtable.cxx 
                        ivoldual_compute.cxx ivoldual_query.cxx 
                        ivoldual_move.cxx ivoldual_reposition.cxx)


ADD_CUSTOM_TARGET(tar WORKING_DIRECTORY . COMMAND tar cvfh ivoldual.tar *.cxx *.h *.txx CMakeLists.txt ivoldual_doxygen.config)

ADD_CUSTOM_TARGET(doc COMMAND doxygen ivoldual_doxygen.config)
