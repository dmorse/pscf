add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f 
   COMMAND ${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/iterate_mod.fp.f > ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f 
   DEPENDS iterate_mod.fp.f
)

add_library(iterate
   ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f 
   response_pd_mod.f
   sweep_mod.f
)
target_include_directories (iterate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (iterate PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(iterate LINK_PUBLIC common io scf scf_aux ${LAPACK_LIBRARIES})
