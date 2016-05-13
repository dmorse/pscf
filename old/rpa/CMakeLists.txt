add_library(rpa
	rpa_mod.f
	spinodal_mod.f
)

target_include_directories (rpa PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (rpa PUBLIC ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(rpa LINK_PUBLIC common io scf iterate ${LAPACK_LIBRARIES})
