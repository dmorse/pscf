# We make two libraries because of some convoluted dependency graphs
add_library(io
		io_mod.f90
		string_mod.f90
		version_mod.f90
)

target_include_directories (io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (io PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(io LINK_PUBLIC common)

add_library(io_field
		field_io_mod.f90
)

target_include_directories (io_field PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (io_field PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(io_field LINK_PUBLIC common crystal scf)
