# TODO: move these options to the top file

if (USE_FFT3)
    set (FFT_FILE fft3_mod.f90)    
    find_package(FFTW)
else (USE_FFT3)
    set (FFT_FILE fft2_mod.f90)    
    find_package(FFTW)
endif (USE_FFT3)

if (USE_DEVEL)
    set (DEVEL "-DDEVEL=1")
endif (USE_DEVEL)

if (BUILD_DMG)
    set (BUILD_MACOSX_BUNDLE "MACOSX_BUNDLE")
    set (MACOSX_BUNDLE_EXT ".app")
endif (BUILD_DMG)

set(FORPEDO ${CMAKE_SOURCE_DIR}/tools/python/preprocess-0.6.1/preprocess.py)

find_package(LAPACK)

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f90
	COMMAND	${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/scf_mod.fp.f90 > ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f90
	DEPENDS scf_mod.fp.f90
)

add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f90
   COMMAND ${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/iterate_mod.fp.f90 > ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f90
   DEPENDS iterate_mod.fp.f90
)

# Create a library called "pscf_lib", which includes module source files.
# The extension is already found. Any number of sources could be listed here.
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

# Ensure that linker can find any headers or fortram *.mod module files,
# including preprocessor geneerated sources placed in the build directory.
# This is the in-core source dir
target_include_directories (pscf_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# This is the out-of-source build dir, with generated sources
target_include_directories (pscf_lib PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

# Link to FFTW and Lapack libraries
target_link_libraries(pscf_lib LINK_PUBLIC ${FFTW_LIBRARIES} ${LAPACK_LIBRARIES})


### PREPROCESS THE MAIN FILE ###
add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pscf.f90 
   COMMAND  ${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/pscf.fp.f90 > ${CMAKE_CURRENT_BINARY_DIR}/pscf.f90 
   DEPENDS pscf.fp.f90
)

add_executable(pscf 
   ${BUILD_MACOSX_BUNDLE} WIN32
   ${CMAKE_CURRENT_BINARY_DIR}/pscf.f90 
)

set(APPS ${CMAKE_CURRENT_BINARY_DIR}/pscf${MACOSX_BUNDLE_EXT})
set(DIRS "")

# TODO: clean up these links
target_link_libraries(pscf LINK_PUBLIC pscf_lib)

install(TARGETS pscf
        BUNDLE DESTINATION . COMPONENT Runtime
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib/static COMPONENT Runtime
)

if (BUILD_DMG)
# Copy libraries into the DMG
    INSTALL(CODE "
       include(BundleUtilities)
       set (BU_CHMOD_BUNDLE_ITEMS ON)
       fixup_bundle(\"${APPS}\"   \"\"   \"${DIRS}\")
    " COMPONENT Runtime)
endif (BUILD_DMG)

#install(TARGETS pscf
#        RUNTIME DESTINATION bin
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib/static
#)
