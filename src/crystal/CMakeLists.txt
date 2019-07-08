# TODO: deform_mod

add_library(crystal
	# NOTE: deform_mod is never actually linked in the original Makefile; perhaps
	# it was never meant to be seen as completed
	# deform_mod.f90
		basis_mod.f90
		grid_basis_mod.f90
)

target_include_directories (crystal PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (crystal PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(crystal LINK_PUBLIC common io grid)

add_library(crystal_aux
		group_mod.f90
		unit_cell_mod.f90
		space_groups_mod.f90
)

target_include_directories (crystal_aux PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (crystal_aux PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(crystal_aux LINK_PUBLIC common io)

