add_executable(rpa_test
    ${BUILD_MACOSX_BUNDLE} WIN32
	test.f
)
target_link_libraries(rpa_test LINK_PUBLIC common io scf rpa)

install(TARGETS rpa_test
        BUNDLE DESTINATION . COMPONENT Runtime
        RUNTIME DESTINATION bin COMPONENT Runtime
        LIBRARY DESTINATION lib COMPONENT Runtime
        ARCHIVE DESTINATION lib/static COMPONENT Runtime
)

set(APPS ${CMAKE_CURRENT_BINARY_DIR}/rpa_test${MACOSX_BUNDLE_EXT})

set(DIRS "")

if (BUILD_DMG)
# Copy libraries into the DMG
INSTALL(CODE "
   include(BundleUtilities)
   set (BU_CHMOD_BUNDLE_ITEMS ON)
   fixup_bundle(\"${APPS}\"   \"\"   \"${DIRS}\")
" COMPONENT Runtime)
endif (BUILD_DMG)

