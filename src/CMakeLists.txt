# TODO: move these options to the top file
option(UseFFT3 "Use the fft3_mod file (and FFT3.0). Otherwise, use fft2_mod (and FFTW2)" ON)
option(UseDEVEL "Use the -D DEVEL option with the Fortran preprocessor" ON)

if (UseFFT3)
    set (FFT_FILE fft3_mod.f)    
    find_package(FFTW)
else (UseFFT3)
    set (FFT_FILE fft2_mod.f)    
    find_package(FFTW)
endif (UseFFT3)

if (UseDEVEL)
    set (DEVEL "-DDEVEL=1")
endif (UseDEVEL)

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
   MACOSX_BUNDLE WIN32
   ${CMAKE_CURRENT_BINARY_DIR}/pscf_pd.f 
)

if (APPLE)
  set(APPS ${CMAKE_CURRENT_BINARY_DIR}/pscf.app)
else (APPLE)
  set(APPS ${CMAKE_CURRENT_BINARY_DIR}/pscf)
endif (APPLE) 
set(DIRS "")

# TODO: clean up these links
target_link_libraries(pscf LINK_PUBLIC common io io_field scf scf_aux response rpa iterate)

install(TARGETS pscf
        BUNDLE DESTINATION . COMPONENT Runtime
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib/static COMPONENT Runtime
)

# Copy libraries into the DMG
INSTALL(CODE "
   include(BundleUtilities)
   set (BU_CHMOD_BUNDLE_ITEMS ON)
   fixup_bundle(\"${APPS}\"   \"\"   \"${DIRS}\")
" COMPONENT Runtime)

#install(TARGETS pscf
#        RUNTIME DESTINATION bin
#        LIBRARY DESTINATION lib
#        ARCHIVE DESTINATION lib/static
#)
