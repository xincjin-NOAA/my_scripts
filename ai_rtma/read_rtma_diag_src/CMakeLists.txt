cmake_minimum_required(VERSION 3.19)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

project(read_rtma_diag VERSION 1.0.0 LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
endif()

# RPATH settings to ensure that the executable can find shared libraries at runtime.
# This is important on HPC systems where libraries are in non-standard locations.

# Use RPATH in the build tree
set(CMAKE_SKIP_BUILD_RPATH FALSE)
# Add the paths to all found dependencies (like NetCDF) to the RPATH.
# This allows the binary to find its libraries at runtime after installation.
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# Try to find NetCDF. spack-stack often provides 'netCDFConfig.cmake' (lowercase n).
find_package(netCDF COMPONENTS Fortran)
if(NOT netCDF_FOUND)
    find_package(NetCDF REQUIRED COMPONENTS Fortran)
endif()

# Standardize target name: spack-stack typically exports NetCDF::NetCDF_Fortran
if(NOT TARGET netCDF::netcdf_fortran AND TARGET NetCDF::NetCDF_Fortran)
    add_library(netCDF::netcdf_fortran ALIAS NetCDF::NetCDF_Fortran)
endif()

# The ncdiag package internally may look for 'NetCDF' (case sensitive).
if(netCDF_FOUND AND NOT NetCDF_FOUND)
    set(NetCDF_DIR ${netCDF_DIR})
endif()

find_package(ncdiag REQUIRED)

# Big-endian unformatted binary read (matches -convert big_endian in the makefile)
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  add_compile_options(-convert big_endian)
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  add_compile_options(-fconvert=big-endian)
endif()

add_executable(read_diag_conv_rtma.exe
  kinds.F90
  read_diag_conv_rtma.f90)

target_link_libraries(read_diag_conv_rtma.exe PRIVATE 
  ncdiag::ncdiag
  netCDF::netcdf_fortran)

install(TARGETS read_diag_conv_rtma.exe RUNTIME DESTINATION bin)
