cmake_minimum_required(VERSION 2.6)

project(OmmTesting)

#include files
include_directories (/opt/google-test/include)
include_directories (/home/ksb12101/openmm/OpenMM/exec/include)
INCLUDE_DIRECTORIES (/home/ksb12101/local/rapidjson/include)


#include libraries to be linked
LINK_DIRECTORIES (/usr/lib)
LINK_DIRECTORIES (/home/ksb12101/openmm/OpenMM/exec/lib)

ENABLE_TESTING()
INCLUDE(CTest)

#automatically compile and add executables by looping over files
FILE(GLOB TEST_PROGS "Test*.cpp")
FOREACH(TEST_PROG ${TEST_PROGS})
	GET_FILENAME_COMPONENT(TEST_ROOT ${TEST_PROG} NAME_WE)
	add_executable(${TEST_ROOT} ${TEST_ROOT})
	TARGET_LINK_LIBRARIES(${TEST_ROOT} gtest OpenMM)
#start adding tests from here
add_test(${TEST_ROOT} ${TEST_ROOT})
#End of foreach for test programs
ENDFOREACH(TEST_PROG ${TEST_PROGS})
