add_executable(rpa_test
	test.f
)
target_link_libraries(rpa_test LINK_PUBLIC common io scf rpa)

install(TARGETS rpa_test
		RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static

)
