cmake_minimum_required( VERSION 3.2 )
project( DarkPhysicsNA64 )

find_package( GTest )

if( BUILDCONF_SUFFIX STREQUAL "" )
    set( DarkPhysicsNA64Test_LIBRARY dph-test-common-${BUILDCONF_SUFFIX} )
else()
    set( DarkPhysicsNA64Test_LIBRARY dph-test-common )
endif()
file(GLOB_RECURSE DarkPhysicsNA64Test_SOURCES common/src/*.c*)
add_library( ${DarkPhysicsNA64Test_LIBRARY} SHARED ${DarkPhysicsNA64Test_SOURCES} )
target_include_directories(${DarkPhysicsNA64Test_LIBRARY}
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/common/include>"
    PUBLIC "$<INSTALL_INTERFACE:dphmc-na64/test>")
target_link_libraries( ${DarkPhysicsNA64Test_LIBRARY}
    PUBLIC ${DarkPhysicsNA64_LIBRARY})
set_target_properties( ${DarkPhysicsNA64Test_LIBRARY} PROPERTIES
    VERSION ${DarkPhysicsNA64_VERSION}
    SOVERSION ${DarkPhysicsNA64_VERSION} )

add_subdirectory( utils )

if( GTEST_FOUND )
    add_subdirectory( units )
endif( GTEST_FOUND )

