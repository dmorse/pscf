add_library(grid
        ${FFT_FILE}
        grid_mod.f90
)

target_include_directories (grid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (grid PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(grid LINK_PUBLIC common io crystal_aux ${FFTW_LIBRARIES})
