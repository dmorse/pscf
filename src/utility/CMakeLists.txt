# TODO: debug for gfortran


add_executable(nto2
		nto2.f
)
target_link_libraries(nto2 LINK_PUBLIC common )


add_executable(nto3
		nto3.f
)
target_link_libraries(nto3 LINK_PUBLIC common )


add_executable(switch
		switch.f
)
target_link_libraries(switch LINK_PUBLIC common )


install(TARGETS nto2 nto3 switch
		RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static
)
