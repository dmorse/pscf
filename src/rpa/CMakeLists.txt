add_library(rpa
	rpa_mod.f90
	spinodal_mod.f90
)

target_include_directories (rpa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (rpa PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(rpa LINK_PUBLIC common io scf iterate ${LAPACK_LIBRARIES})
