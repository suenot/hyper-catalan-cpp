#include <gtest/gtest.h>
#include "../src/hyper_catalan_solver.h"
#include <iostream>

// This test demonstrates the computation of Hyper-Catalan numbers
// and provides examples from the Wildberger and Rubin paper
TEST(ExampleTest, HyperCatalanNumbersExample) {
    HyperCatalanCalculator calculator;
    
    std::cout << "Hyper-Catalan Numbers Examples:" << std::endl;
    
    // Simple subdigon types and their expected Hyper-Catalan numbers
    std::vector<std::pair<std::vector<int>, HighPrecisionFloat>> examples = {
        {{1, 0, 0}, 1},     // C_{1,0,0} = 1
        {{2, 0, 0}, 2},     // C_{2,0,0} = 2
        {{1, 1, 0}, 3},     // C_{1,1,0} = 3
        {{3, 0, 0}, 5},     // C_{3,0,0} = 5
        {{0, 2, 0}, 12},    // C_{0,2,0} = 12
        {{2, 1, 0}, 12}     // C_{2,1,0} = 12
    };
    
    for (const auto& example : examples) {
        SubdigonType type(example.first);
        HighPrecisionFloat expected = example.second;
        HighPrecisionFloat calculated = calculator.calculate(type);
        
        std::cout << "C_{";
        for (size_t i = 0; i < example.first.size(); ++i) {
            std::cout << example.first[i];
            if (i < example.first.size() - 1) {
                std::cout << ",";
            }
        }
        std::cout << "} = " << calculated << std::endl;
        
        EXPECT_EQ(calculated, expected);
    }
}

// This test demonstrates the solution of a quadratic equation
// using the Hyper-Catalan series method
TEST(ExampleTest, QuadraticEquationExample) {
    // Create a solver for quadratic equations with 15 terms
    HyperCatalanPolynomialSolver solver(2, 15);
    
    std::cout << "\nQuadratic Equation Example:" << std::endl;
    std::cout << "x² - 2x - 3 = 0" << std::endl;
    
    // Define the quadratic: x² - 2x - 3 = 0
    // Coefficient form: -3 - 2x + x²
    std::vector<HighPrecisionFloat> coefficients = {
        -3,   // x⁰
        -2,   // x¹
        1     // x²
    };
    
    // Expected roots: -1 and 3
    
    // Solve using the Hyper-Catalan series
    HighPrecisionFloat series_root;
    try {
        series_root = solver.solvePolynomial(coefficients);
        std::cout << "Root approximation from Hyper-Catalan series: " << series_root << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        FAIL() << "Failed to solve with Hyper-Catalan series";
    }
    
    // Refine using bootstrap method
    HighPrecisionFloat bootstrap_root = solver.bootstrapRoot(coefficients, series_root, 10);
    std::cout << "Refined root using bootstrap method: " << bootstrap_root << std::endl;
    
    // Check which root we found (should be either -1 or 3)
    bool is_valid_root = (abs(bootstrap_root - 3) < 0.1) || (abs(bootstrap_root + 1) < 0.1);
    EXPECT_TRUE(is_valid_root) << "Expected root to be either -1 or 3, but got " << bootstrap_root;
    
    // Verify the polynomial evaluates to approximately zero at the root
    HighPrecisionFloat polynomial_value = 0;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        polynomial_value += coefficients[i] * pow(bootstrap_root, i);
    }
    
    std::cout << "Polynomial evaluated at root: " << polynomial_value << std::endl;
    EXPECT_LT(abs(polynomial_value), 1e-10);
}

// This test demonstrates the cubic equation example from the paper
TEST(ExampleTest, CubicEquationExample) {
    // Create a solver for cubic equations with 15 terms
    HyperCatalanPolynomialSolver solver(3, 15);
    
    std::cout << "\nCubic Equation Example:" << std::endl;
    std::cout << "x³ - 6x² + 11x - 6 = 0" << std::endl;
    
    // Define the cubic: x³ - 6x² + 11x - 6 = 0 (factored as (x-1)(x-2)(x-3) = 0)
    std::vector<HighPrecisionFloat> coefficients = {
        -6,   // x⁰
        11,   // x¹
        -6,   // x²
        1     // x³
    };
    
    // Expected roots: 1, 2, and 3
    
    // Solve using the Hyper-Catalan series
    HighPrecisionFloat series_root;
    try {
        series_root = solver.solvePolynomial(coefficients);
        std::cout << "Root approximation from Hyper-Catalan series: " << series_root << std::endl;
    } catch (const std::exception& e) {
        std::cout << "Error: " << e.what() << std::endl;
        FAIL() << "Failed to solve with Hyper-Catalan series";
    }
    
    // Refine using bootstrap method
    HighPrecisionFloat bootstrap_root = solver.bootstrapRoot(coefficients, series_root, 10);
    std::cout << "Refined root using bootstrap method: " << bootstrap_root << std::endl;
    
    // Check which root we found (should be one of 1, 2, or 3)
    bool is_valid_root = (abs(bootstrap_root - 1) < 0.1) || 
                         (abs(bootstrap_root - 2) < 0.1) || 
                         (abs(bootstrap_root - 3) < 0.1);
    
    EXPECT_TRUE(is_valid_root) << "Expected root to be 1, 2, or 3, but got " << bootstrap_root;
    
    // Verify the polynomial evaluates to approximately zero at the root
    HighPrecisionFloat polynomial_value = 0;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        polynomial_value += coefficients[i] * pow(bootstrap_root, i);
    }
    
    std::cout << "Polynomial evaluated at root: " << polynomial_value << std::endl;
    EXPECT_LT(abs(polynomial_value), 1e-10);
}

// Example of finding multiple roots by using different initial guesses
TEST(ExampleTest, MultipleRootsExample) {
    HyperCatalanPolynomialSolver solver(3, 15);
    
    std::cout << "\nMultiple Roots Example (using different initial guesses):" << std::endl;
    std::cout << "x³ - 6x² + 11x - 6 = 0 with roots 1, 2, and 3" << std::endl;
    
    // Define the cubic: x³ - 6x² + 11x - 6 = 0 (factored as (x-1)(x-2)(x-3) = 0)
    std::vector<HighPrecisionFloat> coefficients = {
        -6,   // x⁰
        11,   // x¹
        -6,   // x²
        1     // x³
    };
    
    // Try different initial guesses to find all three roots
    std::vector<HighPrecisionFloat> initial_guesses = {0.5, 1.5, 2.5};
    std::vector<HighPrecisionFloat> expected_roots = {1.0, 2.0, 3.0};
    
    for (size_t i = 0; i < initial_guesses.size(); ++i) {
        HighPrecisionFloat guess = initial_guesses[i];
        HighPrecisionFloat expected = expected_roots[i];
        
        HighPrecisionFloat root = solver.bootstrapRoot(coefficients, guess, 10);
        
        std::cout << "Initial guess: " << guess << ", Found root: " << root 
                  << ", Expected: " << expected << std::endl;
        
        EXPECT_NEAR(static_cast<double>(root), static_cast<double>(expected), 0.1);
    }
} 