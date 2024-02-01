cmake_minimum_required(VERSION 3.10)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Set the project name
project(UnitTests)

# Set the relative path to svFSIplus include directory
set(SVFSI_INCLUDE_DIR ../../Code/Source/svFSI)
# set(SVFSI_THIRDPARTY ../../Code/ThirdParty/msmpi/Include)

# Include the directory
include_directories(${SVFSI_INCLUDE_DIR})
# include_directories(${SVFSI_THIRDPARTY})

# Define test executable
add_executable(runUnitTests test.cpp)

# Link against Google Test and svFSIplus libraries
# target_link_libraries(runUnitTests gtest gtest_main)

# Add the test so it can be run with `make test`
# add_test(NAME that_test COMMAND runUnitTests)

