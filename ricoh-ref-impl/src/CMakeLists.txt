cmake_minimum_required(VERSION 2.8)

set(VERBOSE_FLAG "")
include(CMakeForceCompiler)
#option(VERBOSE "verbose debug output" OFF)

#if(VERBOSE)
#    set(VERBOSE_FLAG -DVERBOSE)
#else(VERBOSE)
#    set(VERBOSE_FLAG "")
#endif(VERBOSE)

#set(VERBOSE_FLAG -DVERBOSE)

set(CMAKE_CXX_FLAGS "-O2 -march=native")
cmake_force_cxx_compiler(/usr/local/bin/g++-4.7 g++)

include_directories(../include)

add_executable(lp main.cpp simplex.cpp time_fnc_simplex.cpp read_lp.cpp)
add_executable(lp1 main1.cpp simplex.cpp)
add_executable(lp2 main2.cpp simplex.cpp)

#---------- generated sources ----------
add_custom_command(
    OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/time_fnc_simplex.cpp ${CMAKE_SOURCE_DIR}/include/time_fnc_simplex.hpp
    COMMAND python gen-timer.py fnc_simplex ${CMAKE_SOURCE_DIR}/include/simplex.hpp ${CMAKE_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/gen
    COMMENT "${CMAKE_BINARY_DIR}/gen/gen-timer.py simplex simplex.hpp"
    DEPENDS gen-timer.py ${CMAKE_SOURCE_DIR}/include/simplex.hpp simplex.cpp
)
