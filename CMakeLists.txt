# CMakeLists files in this project can
# refer to the root source directory of the project as ${PSCF_SOURCE_DIR} and
# to the root binary directory of the project as ${PSCF_BINARY_DIR}.
cmake_minimum_required (VERSION 2.8.11)
project (PSCF)

#### Based on https://cmake.org/Wiki/CMakeForFortranExamplea ####
enable_language (Fortran)

# Add our extra Find*.cmake modules
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/tools/cmake/")

# make sure that the default is a RELEASE
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
      "Choose the type of build, options are: None Debug Release."
      FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# default installation
#get_filename_component (default_prefix ".." ABSOLUTE)
#set (CMAKE_INSTALL_PREFIX ${default_prefix} CACHE STRING
#      "Choose the installation directory; by default it installs in the NORMA directory."
#      FORCE)

# FFLAGS depend on the compiler
get_filename_component (Fortran_COMPILER_NAME ${CMAKE_Fortran_COMPILER} NAME)

if (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  # gfortran
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -ffree-form")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -ffree-form -g")
elseif (Fortran_COMPILER_NAME MATCHES "ifort.*")
  # ifort (untested)
  set (CMAKE_Fortran_FLAGS_RELEASE "-f77rtl -O3")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-f77rtl -O0 -g")
elseif (Fortran_COMPILER_NAME MATCHES "g77")
  # g77
  set (CMAKE_Fortran_FLAGS_RELEASE "-funroll-all-loops -fno-f2c -O3 -m32")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-fno-f2c -O0 -g -m32")
else (Fortran_COMPILER_NAME MATCHES "gfortran.*")
  message ("CMAKE_Fortran_COMPILER full path: " ${CMAKE_Fortran_COMPILER})
  message ("Fortran compiler: " ${Fortran_COMPILER_NAME})
  message ("No optimized Fortran compiler flags are known, we just try -O2...")
  set (CMAKE_Fortran_FLAGS_RELEASE "-O2")
  set (CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
endif (Fortran_COMPILER_NAME MATCHES "gfortran.*")
#### End based on https://cmake.org/Wiki/CMakeForFortranExamplea ####


# Recurse into the "PSCF" subdirectory. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)


# Setup the packaging
set(CPACK_PACKAGE_NAME "pscf")
set(CPACK_PACKAGE_VENDOR "MSI - UMN")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "pscf - Polymer self-consistent field theory http://research.cems.umn.edu/morse/code/pscf/home.php")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "0")
# TODO: build out components "applications libraries headers" (where headers are the .mod files)
set(CPACK_COMPONENTS_ALL Unspecified)

# TODO: finish Bundle 
# TODO: finish DEB and RPM generators (http://www.vtk.org/Wiki/CMake:Component_Install_With_CPack)
# TODO: finish NSIS generator (http://www.vtk.org/Wiki/CMake:Component_Install_With_CPack)
# TODO: test compile and package on Unix and Windows
if (APPLE)
	set(CPACK_GENERATOR 
		Bundle 
		TGZ
		ZIP
	)
else (APPLE)
	if (UNIX)
		set(CPACK_GENERATOR 
			DEB
			RPM
			TGZ
			ZIP
		)
	else (UNIX)
		set(CPACK_GENERATOR 
			NSIS	
			ZIP
		)
	endif(UNIX)
endif (APPLE)

set(CPACK_BUNDLE_NAME "pscf")
# TODO: make valid Info.plist and valid Icon
set(CPACK_BUNDLE_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/tools/bundle/Info.plist)
set(CPACK_BUNDLE_ICON ${CMAKE_CURRENT_SOURCE_DIR}/tools/bundle/Icon.png)

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Evan Bollig") #required for DEB
 
include(CPack)
