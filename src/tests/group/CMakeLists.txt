add_executable(new_group
    MACOSX_BUNDLE WIN32
	new_group.f
)
target_link_libraries(new_group LINK_PUBLIC common io scf rpa )


install(TARGETS new_group
        BUNDLE DESTINATION . COMPONENT Runtime
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib/static COMPONENT Runtime
)

set(APPS ${CMAKE_CURRENT_BINARY_DIR}/new_group.app)
set(DIRS "")
# Copy libraries into the DMG
INSTALL(CODE "
   include(BundleUtilities)
   set (BU_CHMOD_BUNDLE_ITEMS ON)
   fixup_bundle(\"${APPS}\"   \"\"   \"${DIRS}\")
" COMPONENT Runtime)


