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

# Create a library called "common", which includes the source file "const_mod.f".
# The extension is already found. Any number of sources could be listed here.
add_library(common 
    const_mod.f
)
# Ensure that anything that depends on "common" will be able to find any headers
# or fortram *.mod module files for linking
# This is the in-core source dir
target_include_directories (common PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
# b) this is the out-of-core source dir (i.e., generated sources)
target_include_directories (common PUBLIC ${CMAKE_CURRENT_BINARY_DIR})


##### Build subdirs ######
add_subdirectory(io)
add_subdirectory(grid)
add_subdirectory(crystal)
add_subdirectory(scf)
add_subdirectory(rpa)
add_subdirectory(iterate)
add_subdirectory(response)

### PREPROCESS THE MAIN FILE ###
add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/pscf_pd.f 
   COMMAND  ${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/pscf_pd.fp.f > ${CMAKE_CURRENT_BINARY_DIR}/pscf_pd.f 
   DEPENDS pscf_pd.fp.f
)

add_executable(pscf 
   ${BUILD_MACOSX_BUNDLE} WIN32
   ${CMAKE_CURRENT_BINARY_DIR}/pscf_pd.f 
)

set(APPS "${CMAKE_CURRENT_BINARY_DIR}/pscf${MACOSX_BUNDLE_EXT}")
set(DIRS "")

# TODO: clean up these links
target_link_libraries(pscf LINK_PUBLIC common io io_field scf scf_aux response rpa iterate)

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
