#!/bin/bash

# Exit on error
set -e

# Create build directory if it doesn't exist
mkdir -p build
cd build

# Configure the project with CMake
echo "Configuring project with CMake..."
cmake ..

# Build the project
echo "Building project..."
make -j$(nproc)

echo "Build completed successfully!"
echo "You can run the solver with: ./hyper_catalan_solver"
echo "You can run the tests with: ./hyper_catalan_tests"
echo "You can run the example tests with: ./example_tests" 