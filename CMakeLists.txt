cmake_minimum_required(VERSION 2.8.3)
project (uclales-cloud-tracking Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set (NETCDF_F03 "YES")
set (NETCDF_F90 "YES")
set (NETCDF_NC4 "YES")
FIND_PACKAGE(NetCDF)

add_executable(main
        src/constants.f90
        src/modnetcdf.f90
        src/constants.f90
        src/tracking_common.f90
        src/modarray.f90
        src/modstatistics.f90
        src/modtrack_cell_splitting.f90
        src/modtrack.f90
        src/tracking_data.f90
        src/field_loader.f90
        src/main.f90)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -g -O0")
set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS_RELEASE}")

if("$ENV{FC}" MATCHES "ifort")
  set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -check all -traceback -debug all")
endif()

target_link_libraries(main ${NETCDF_LIBRARIES})
include_directories(${NETCDF_INCLUDES})

# using pytest for tests for so we add a custom command instead of using
# cmake's inbuilt test system. pytest looks for the `test/` directory so needs
# to be working in the source root
if(${CMAKE_VERSION} VERSION_LESS "3.2.0")
  # the python test script will have to try and work out for itself where the
  # binary is stored
  add_custom_target(check
    COMMAND python -m pytest
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS main
  )
else()
  add_custom_target(check
    -E env "CMAKE_BINARY_DIR=${CMAKE_BINARY_DIR}"
    COMMAND python -m pytest
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    DEPENDS main
  )
endif()


#add_executable(orig_tracking_time
        #orig_src/tracking_time.f90)
#target_link_libraries(orig_tracking_time ${NETCDF_LIBRARIES_F90})


# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" OFF)

# check if Doxygen is installed
find_package(Doxygen)

if (BUILD_DOC)
  if (DOXYGEN_FOUND)
      # set input and output files
      set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in)
      set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
  
      # request to configure the file
      configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
      message("Doxygen build started")
  
      # note the option ALL which allows to build the docs together with the application
      add_custom_target( doc_doxygen ALL
          COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
          COMMENT "Generating documentation with Doxygen"
          VERBATIM )
  else (DOXYGEN_FOUND)
    message("Doxygen need to be installed to generate the doxygen documentation")
  endif (DOXYGEN_FOUND)
endif()
