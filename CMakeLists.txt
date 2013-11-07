cmake_minimum_required(VERSION 2.8)

# Find the BEM++ library installed on your disk
find_package(Bempp REQUIRED CONFIG)

find_library(BEMPP_LIBRARY bempp PATH ${BEMPP_LIBRARY_DIR})
find_library(BEMPP_TEUCHOSCORE_LIBRARY teuchoscore PATH ${BEMPP_LIBRARY_DIR})
find_library(BEMPP_TEUCHOSPARAMETERLIST_LIBRARY teuchosparameterlist PATH ${BEMPP_LIBRARY_DIR})

# Add the BEM++ installation directory to CMake search path to ensure that
# e.g. the version of SWIG coming with BEM++ is used
set(CMAKE_PREFIX_PATH "${BEMPP_PREFIX}/bempp;${CMAKE_PREFIX_PATH}"
    CACHE STRING "")
# Make the new module be installed in the main BEM++ directory, by default
set(CMAKE_INSTALL_PREFIX "${BEMPP_PREFIX}/bempp" CACHE PATH "" FORCE)

# By default, compile the new module with optimizations turned on
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release" CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
       FORCE)
endif ()

# set(CMAKE_CXX_FLAGS "-Wall -Wnon-virtual-dtor -Wno-sign-compare")
# SET(GCC_COMPILE_FLAGS " -DNDEBUG -DBOOST_UBLAS_NDEBUG -O3")
# SET(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} ${GCC_COMPILE_FLAGS}")

# Main library
add_subdirectory(lib)

# Python Support
add_subdirectory(python)

