cmake_minimum_required( VERSION 3.2 )
project( DarkPhysicsNA64 )

set( DarkPhysicsNA64_VERSION 0.0.1 )

option( TESTS "Enables testing utils build" ON )
option( BUILD_DOC "Whether to generate the documentation (Doxygen used)" ON )
#option( BUILD_EXAMPLES "Enables build of example applications" OFF )

find_package( Geant4 REQUIRED )
find_package( GSL REQUIRED )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel. By default Debug is assumed."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if( ${CMAKE_BUILD_TYPE} STREQUAL "Debug" )
    set( BUILDCONF_SUFFIX dbg )
elseif( ${CMAKE_BUILD_TYPE} STREQUAL "MinSizeRel" )
    set( BUILDCONF_SUFFIX minsize )
endif()

if( BUILDCONF_SUFFIX STREQUAL "" )
    set( DarkPhysicsNA64_LIBRARY dph-${BUILDCONF_SUFFIX} )
else()
    set( DarkPhysicsNA64_LIBRARY dph )
endif()

configure_file( include/dphmc-config.h.in include/dphmc-config.h @ONLY )

file(GLOB_RECURSE DarkPhysicsNA64_SOURCES src/*.c*)
add_library( ${DarkPhysicsNA64_LIBRARY} SHARED ${DarkPhysicsNA64_SOURCES} )
target_compile_options( ${DarkPhysicsNA64_LIBRARY} PRIVATE -Wall -Wextra )
target_include_directories( ${DarkPhysicsNA64_LIBRARY}
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>"
    PUBLIC "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>"
    PUBLIC "$<INSTALL_INTERFACE:dphmc-na64>" )
target_include_directories( ${DarkPhysicsNA64_LIBRARY}
    SYSTEM PUBLIC ${Geant4_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS} )
target_link_libraries( ${DarkPhysicsNA64_LIBRARY}
    PUBLIC ${Geant4_LIBRARIES} ${GSL_LIBRARIES} )
set_target_properties( ${DarkPhysicsNA64_LIBRARY} PROPERTIES
    VERSION ${DarkPhysicsNA64_VERSION}
    SOVERSION ${DarkPhysicsNA64_VERSION} )

#get_target_property(debug ${DarkPhysicsNA64_LIBRARY} INTERFACE_INCLUDE_DIRECTORIES)
#message(STATUS DarkPhysicsNA64_LIBRARY "${debug}")

if( BUILD_TESTS )
    add_subdirectory( tests )
endif( BUILD_TESTS )

#
# Documentation

if( BUILD_DOC )
    find_package( Doxygen QUIET )
    # check if Doxygen is installed
    if( DOXYGEN_FOUND )
        #set(GEANT4_CPP_TAGFILE ${CMAKE_BINARY_DIR}/cppreference-doxygen-web.tag.xml)
        #^^^ consider e.g.: http://www-geant4.kek.jp/Reference/10.01/geant4.tag
        set(STL_CPP_TAGFILE ${CMAKE_BINARY_DIR}/cppreference-doxygen-web.tag.xml)
        # Try to fetch the web tag file for STL references.
        if( NOT EXISTS ${STL_CPP_TAGFILE} )
            file( DOWNLOAD http://upload.cppreference.com/mwiki/images/f/f8/cppreference-doxygen-web.tag.xml
                           ${STL_CPP_TAGFILE}
                           TIMEOUT 20
                           EXPECTED_HASH SHA1=1df9e91490837590c5bb1f40c8cbad146370d00c
                           TLS_VERIFY ON )
        endif()
        configure_file( Doxyfile.in Doxyfile @ONLY )
        if( EXISTS ${STL_CPP_TAGFILE} )
            message(STATUS "STL reference links will be included in documentation.")
            # TODO: this below does not work. Figure out what is going on and
            # remove TAGFILES=@STL... crutch from Doxyfile.in.
            #file(APPEND "Doxyfile" "TAGFILES+=\"${STL_CPP_TAGFILE}=http://en.cppreference.com/w/\"")
        endif()
        # Note the option ALL which allows to build the docs together with the
        # rest of project
        add_custom_target( doc ALL
            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            COMMENT "Generating library docs with Doxygen"
            VERBATIM )
    else (DOXYGEN_FOUND)
        message(WARNING "Doxygen is not found, documentation won't be generated.")
    endif (DOXYGEN_FOUND)
endif( BUILD_DOC )

