cmake_minimum_required(VERSION 3.14)
project(hyper_catalan_solver VERSION 1.0)

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find external packages
find_package(Boost REQUIRED)
find_package(Eigen3 REQUIRED)

# Include directories
include_directories(${BOOST_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIRS})

# Main executable
add_executable(hyper_catalan_solver src/main.cpp src/hyper_catalan_solver.cpp)

# Enable testing
enable_testing()

# Add Google Test
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# Configure Google Test to compile standalone
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Unit tests
add_executable(hyper_catalan_tests test/test_hyper_catalan.cpp src/hyper_catalan_solver.cpp)
target_link_libraries(hyper_catalan_tests gtest_main)

# Example tests
add_executable(example_tests test/example_tests.cpp src/hyper_catalan_solver.cpp)
target_link_libraries(example_tests gtest_main)

include(GoogleTest)
gtest_discover_tests(hyper_catalan_tests)
gtest_discover_tests(example_tests)

# Installation
install(TARGETS hyper_catalan_solver
    RUNTIME DESTINATION bin
)

# Documentation
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Hyper-Catalan Series Polynomial Solver")
set(CPACK_PACKAGE_VENDOR "Hyper-Catalan Team")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${PROJECT_VERSION_MINOR})

include(CPack) 