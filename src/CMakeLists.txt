
# Preprocess scf_mod module
add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f90
	COMMAND	${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/scf_mod.fp.f90 > ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f90
	DEPENDS scf_mod.fp.f90
)

# Preprocess iterate_mod module
add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f90
   COMMAND ${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/iterate_mod.fp.f90 > ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f90
   DEPENDS iterate_mod.fp.f90
)

# Choose version of wrapper for FFTW library
if (USE_FFT3)
    set (FFT_FILE fft3_mod.f90)    
    find_package(FFTW)
else (USE_FFT3)
    set (FFT_FILE fft2_mod.f90)    
    find_package(FFTW)
endif (USE_FFT3)

# Create library called "pscf_lib" that contains all module files.
add_library(pscf_lib
    const_mod.f90
    io_mod.f90
    string_mod.f90
    version_mod.f90
    field_io_mod.f90
    ${FFT_FILE}
    grid_mod.f90
    basis_mod.f90
    grid_basis_mod.f90
    group_mod.f90
    unit_cell_mod.f90
    space_groups_mod.f90
    chain_mod.f90
    chemistry_mod.f90
    step_mod.f90
    ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f90
    ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f90 
    response_pd_mod.f90
    sweep_mod.f90
    extrapolate_mod.f90
    response_mod.f90
    response_step_mod.f90
    spinodal_mod.f90
)

# Add source and build directories to list of include directories
# used when building pscf_lib.  We add the build directory to ensure 
# that the compiler/archiver can find fortran module (*.mod) files 
# and preprocessor generated source (*.f90) files, all of which are 
# placed in the build directory.

# Add the source directory.
target_include_directories (pscf_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Add the out-of-source build dir
target_include_directories (pscf_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# Link pscf_lib to FFTW and Lapack libraries
target_link_libraries(pscf_lib LINK_PUBLIC ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES})

# Preprocess the main pscf program
add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pscf.f90 
   COMMAND  ${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/pscf.fp.f90 > ${CMAKE_CURRENT_BINARY_DIR}/pscf.f90 
   DEPENDS pscf.fp.f90
)

# Build executable pscf from pscf.f90
add_executable(pscf 
   ${BUILD_MACOSX_BUNDLE} WIN32
   ${CMAKE_CURRENT_BINARY_DIR}/pscf.f90 
)

target_link_libraries(pscf LINK_PUBLIC pscf_lib)

# Install target
install(TARGETS pscf
        BUNDLE DESTINATION . COMPONENT Runtime
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib/static COMPONENT Runtime
)

set(APPS ${CMAKE_CURRENT_BINARY_DIR}/pscf${MACOSX_BUNDLE_EXT})
set(DIRS "")

if (BUILD_DMG)
# Copy libraries into the DMG
    INSTALL(CODE "
       include(BundleUtilities)
       set (BU_CHMOD_BUNDLE_ITEMS ON)
       fixup_bundle(\"${APPS}\"   \"\"   \"${DIRS}\")
    " COMPONENT Runtime)
endif (BUILD_DMG)

