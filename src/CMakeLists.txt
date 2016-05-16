# TODO: move these options to the top file

if (USE_FFT3)
    set (FFT_FILE fft3_mod.f)    
    find_package(FFTW)
else (USE_FFT3)
    set (FFT_FILE fft2_mod.f)    
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
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f
	COMMAND	${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/scf_mod.fp.f > ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f
	DEPENDS scf_mod.fp.f
)

add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f 
   COMMAND ${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/iterate_mod.fp.f > ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f 
   DEPENDS iterate_mod.fp.f
)

# Create a library called "pscf_lib", which includes module source files.
# The extension is already found. Any number of sources could be listed here.
add_library(pscf_lib
    const_mod.f
    io_mod.f
    string_mod.f
    version_mod.f
    field_io_mod.f
    ${FFT_FILE}
    grid_mod.f
    basis_mod.f
    grid_basis_mod.f
    group_mod.f
    unit_cell_mod.f
    space_groups_mod.f
    chain_mod.f
    chemistry_mod.f
    step_mod.f
    ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f
    ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f 
    response_pd_mod.f
    sweep_mod.f
    extrapolate_mod.f
    response_mod.f
    response_step_mod.f
    rpa_mod.f
    spinodal_mod.f
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
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pscf.f 
   COMMAND  ${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/pscf.fp.f > ${CMAKE_CURRENT_BINARY_DIR}/pscf.f 
   DEPENDS pscf.fp.f
)

add_executable(pscf 
   ${BUILD_MACOSX_BUNDLE} WIN32
   ${CMAKE_CURRENT_BINARY_DIR}/pscf.f 
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
