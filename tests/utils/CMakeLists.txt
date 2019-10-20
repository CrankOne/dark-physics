cmake_minimum_required( VERSION 3.2 )

file( GLOB_RECURSE DarkPhysicsNA64_test_utils_SOURCES *.c* )

file( GLOB_RECURSE DarkPhysicsNA64_test_utils_SOURCES
      "*.c*" )

# Note: produced executables are not sensitive for build suffix.
foreach(file_path ${DarkPhysicsNA64_test_utils_SOURCES})
    string( REPLACE ".cpp" "" new_name ${file_path} )
    get_filename_component(filename ${new_name} NAME)
    add_executable( ${filename} ${file_path} )
    target_include_directories(${filename}
        PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")
    target_link_libraries( ${filename}
        PUBLIC ${DarkPhysicsNA64Test_LIBRARY})
    target_include_directories(${filename}
        SYSTEM PRIVATE "${GTEST_INCLUDE_DIRS}")
endforeach()

