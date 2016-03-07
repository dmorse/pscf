add_executable(new_group
	new_group.f
)
target_link_libraries(new_group LINK_PUBLIC common io scf rpa )

install(TARGETS new_group
		RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
)
