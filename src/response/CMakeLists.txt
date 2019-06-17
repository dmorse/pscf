add_library(response
		extrapolate_mod.f90
		response_mod.f90
		response_step_mod.f90
)

target_include_directories (response PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (response PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(response LINK_PUBLIC common io scf io_field)
