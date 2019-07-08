add_custom_command(
   OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f90 
   COMMAND ${FORPEDO} ${DEVEL} ${CMAKE_CURRENT_SOURCE_DIR}/iterate_mod.fp.f90 > ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f90 
   DEPENDS iterate_mod.fp.f90
)

add_library(iterate
   ${CMAKE_CURRENT_BINARY_DIR}/iterate_mod.f90 
   response_pd_mod.f90
   sweep_mod.f90
)
target_include_directories (iterate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories (iterate PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
target_link_libraries(iterate LINK_PUBLIC common io scf scf_aux ${LAPACK_LIBRARIES})
