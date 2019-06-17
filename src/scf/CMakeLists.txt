# TODO: Cleanup libraries (can these be combined in any way? Can I rename them to clarify what is going on?)

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f90
	COMMAND	${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/scf_mod.fp.f90 > ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f90
	DEPENDS scf_mod.fp.f90
)

add_library(scf
            chain_mod.f90
            chemistry_mod.f90
            step_mod.f90
)
target_include_directories (scf PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (scf PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(scf LINK_PUBLIC common io grid)

add_library(scf_aux
            ${CMAKE_CURRENT_BINARY_DIR}/scf_mod.f90
)
target_include_directories (scf_aux PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (scf_aux PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(scf_aux LINK_PUBLIC common io grid crystal)
