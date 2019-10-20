file( GLOB_RECURSE DarkPhysicsNA64_tests_SOURCES *.c* )

if( BUILDCONF_SUFFIX STREQUAL "" )
    set( DarkPhysicsNA64_UNITTESTS dph-na64-ut-tests-${BUILDCONF_SUFFIX} )
else()
    set( DarkPhysicsNA64_UNITTESTS dph-na64-ut-tests )
endif()

add_executable( ${DarkPhysicsNA64_UNITTESTS} ${DarkPhysicsNA64_tests_SOURCES} )
target_include_directories( ${DarkPhysicsNA64_UNITTESTS} PUBLIC include SYSTEM ${GTEST_INCLUDE_DIRS} )
target_link_libraries( ${DarkPhysicsNA64_UNITTESTS} PUBLIC ${DarkPhysicsNA64Test_LIBRARY} ${GTEST_BOTH_LIBRARIES} )

