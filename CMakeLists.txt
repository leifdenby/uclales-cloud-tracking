cmake_minimum_required(VERSION 2.8.3)
project (uclales-cloud-tracking Fortran)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

set (NETCDF_NC4 "YES")
FIND_PACKAGE(NetCDF)

add_executable(main
        src/modnetcdf.f90
        src/main.f90)

target_link_libraries(main ${NETCDF_LIBRARIES})
