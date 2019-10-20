cmake_minimum_required( VERSION 3.2 )
project( DarkPhysicsNA64 )

set( DarkPhysicsNA64_VERSION 0.0.1 )

option( TESTS "Enables testing utils build" OFF )
option( HTML_DOC "Whether to generate the documentation (Doxygen used)" OFF )
option( CODECOV "If in Debug, adds compiler flags for code coverage (GNU only)" OFF )
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
    # Code coverage analysis with gcov is enabled only for debug build as
    # for optimized its output may be quite misleading and generally has no
    # sense
    if( CODECOV )
        # Get the file from
        #   https://raw.githubusercontent.com/bilke/cmake-modules/master/CodeCoverage.cmake
        # and put it to the sources directory. We do not include it in repo by
        # default since code coverage tests is not something average user would
        # do.
        include(${CMAKE_SOURCE_DIR}/CodeCoverage.cmake
                RESULT_VARIABLE CODECOV_ENABLED)
        APPEND_COVERAGE_COMPILER_FLAGS()
        # TODO: line below takes no effect
        set(COVERAGE_LCOV_EXCLUDES 'tests/utils/*' 'tests/common/*')
    endif( CODECOV )
elseif(${CMAKE_BUILD_TYPE} STREQUAL "MinSizeRel")
    set( BUILDCONF_SUFFIX minsize )
else()
    set( BUILDCONF_SUFFIX "" )
endif()

if( BUILDCONF_SUFFIX STREQUAL "" )
    set( DarkPhysicsNA64_LIBRARY dphmc )
else()
    set( DarkPhysicsNA64_LIBRARY dphmc-${BUILDCONF_SUFFIX} )
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

if( TESTS )
    add_subdirectory( tests )
else()
    message( STATUS "Testing utils disabled." )
endif( TESTS )

#
# Documentation

if( HTML_DOC )
    find_package( Doxygen QUIET )
    # check if Doxygen is installed
    if( DOXYGEN_FOUND )
        set(STL_CPP_TAGFILE ${CMAKE_BINARY_DIR}/cppreference-doxygen-web.tag.xml)
        # + Try to fetch the web tag file for STL references
        if( NOT EXISTS ${STL_CPP_TAGFILE} )
            file( DOWNLOAD http://upload.cppreference.com/mwiki/images/f/f8/cppreference-doxygen-web.tag.xml
                           ${STL_CPP_TAGFILE}
                           TIMEOUT 20
                           EXPECTED_HASH SHA1=1df9e91490837590c5bb1f40c8cbad146370d00c
                           TLS_VERIFY ON )
        endif()
        # + Try to fetch the web tag file for Geant4 references
        # NOTE: Geant4 at some point screwed up their version numbering, so to
        # get a correct path, one need to apply some CMake voodoo for the
        # versions > 4.1 && 4.10
        # - prefix minor version number with zero:
        string( REGEX REPLACE
                "\\.([0-9])\\." ".0\\1."
                Geant4_CANONICAL_VERSION
                ${Geant4_VERSION} )
        # - get rid of trailing patch version number
        string( REGEX REPLACE
                "\\.[0-9]+$" ""
                Geant4_CANONICAL_VERSION
                ${Geant4_CANONICAL_VERSION} )
        set(GEANT4_CPP_TAGFILE ${CMAKE_BINARY_DIR}/geant4.${Geant4_CANONICAL_VERSION}.tag)
        if( NOT EXISTS ${GEANT4_CPP_TAGFILE} )
            # NOTE: this file may have significant size (few tens of Mbs)
            file( DOWNLOAD http://www-geant4.kek.jp/Reference/${Geant4_CANONICAL_VERSION}/geant4.tag
                           ${GEANT4_CPP_TAGFILE}
                           TIMEOUT 300 )
        endif()
        #^^^ consider e.g.: http://www-geant4.kek.jp/Reference/10.01/geant4.tag
        #                   http://www-geant4.kek.jp/Reference/"version number"/geant4.tag
        configure_file( Doxyfile.in Doxyfile @ONLY )
        if( EXISTS ${STL_CPP_TAGFILE} )
            message(STATUS "STL reference links will be included in documentation.")
            # TODO: this below does not work. Figure out what is going on and
            # remove TAGFILES=@STL... crutch from Doxyfile.in.
            #file(APPEND "Doxyfile" "TAGFILES+=\"${STL_CPP_TAGFILE}=http://en.cppreference.com/w/\"")
        else()
            message(WARNING "STL reference links will be included in documentation.")
        endif()
        if( EXISTS ${GEANT4_CPP_TAGFILE} )
            message(STATUS "Geant4 reference links (${Geant4_VERSION} -> ${Geant4_CANONICAL_VERSION}) will be included in documentation.")
        else()
            message(WARNING "Geant4 reference links will be included in documentation.")
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
endif( HTML_DOC )

