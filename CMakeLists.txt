cmake_minimum_required(VERSION 2.6)

project(OmmTesting)

#include files
include_directories (/opt/google-test/include)
include_directories (/home/ksb12101/openmm/openmm4.1/include)
INCLUDE_DIRECTORIES (/home/ksb12101/local/rapidjson/include)

#include libraries to be linked
LINK_DIRECTORIES (/usr/lib)
LINK_DIRECTORIES (/home/ksb12101/openmm/openmm4.1/lib)
ENABLE_TESTING()
INCLUDE(CTest)

add_executable(TestHello TestHello.cpp)
add_executable(TestGtest TestSimpleGtest.cpp)
TARGET_LINK_LIBRARIES(TestGtest gtest)
add_executable(TestOpenMM TestOpenMM.cpp)
TARGET_LINK_LIBRARIES(TestOpenMM OpenMM)
add_executable(TestBinProperties TestBinProperties.cpp)
TARGET_LINK_LIBRARIES(TestBinProperties gtest OpenMM)
add_executable(TestControlBerendsenBins TestControlBerendsenBins.cpp)
TARGET_LINK_LIBRARIES(TestControlBerendsenBins gtest OpenMM)
#start adding tests from here
add_test(HelloTest TestHello)
add_test(GtestTest TestGtest)
add_test(OpenMMTest TestOpenMM)
add_test(BinPropertiesTest TestBinProperties)
add_test(ControlBerendsenBinsTest TestControlBerendsenBins)

