add_executable(new_2dgroups
	new_2dgroups.f
)
target_link_libraries(new_2dgroups LINK_PUBLIC common io scf rpa)

install(TARGETS new_2dgroups
		RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
)
