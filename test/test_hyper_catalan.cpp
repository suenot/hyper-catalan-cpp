#include <gtest/gtest.h>
#include "../src/hyper_catalan_solver.h"
#include <cmath>

// Test SubdigonType calculations
TEST(SubdigonTypeTest, Calculations) {
    // Create a subdigon type (2, 1, 0) - meaning 2 digons, 1 trigon, 0 tetragons
    std::vector<int> values = {2, 1, 0};
    SubdigonType type(values);
    
    // Check faces, edges, and vertices calculations
    EXPECT_EQ(type.faces(), 3); // 2 + 1 + 0 = 3 faces
    EXPECT_EQ(type.edges(), 4); // (2*2 + 3*1 + 4*0)/2 = 7/2 = 3.5 -> 4 edges (integer rounding up for SubdigonType calculations)
    EXPECT_EQ(type.vertices(), 3); // edges - faces + 2 = 4 - 3 + 2 = 3 vertices
}

// Test HyperCatalanCalculator
TEST(HyperCatalanCalculatorTest, BasicCalculations) {
    HyperCatalanCalculator calculator;
    
    // Test with simple subdigon type (1, 0, 0) - 1 digon, 0 trigons, 0 tetragons
    SubdigonType type1({1, 0, 0});
    HighPrecisionFloat result1 = calculator.calculate(type1);
    EXPECT_EQ(result1, 1); // The Hyper-Catalan number for (1,0,0) should be 1
    
    // Test with subdigon type (2, 0, 0) - 2 digons, 0 trigons, 0 tetragons
    SubdigonType type2({2, 0, 0});
    HighPrecisionFloat result2 = calculator.calculate(type2);
    EXPECT_EQ(result2, 2); // The Hyper-Catalan number for (2,0,0) should be 2
    
    // Test caching - calling calculate on the same type again should use cached result
    HighPrecisionFloat cached_result = calculator.calculate(type1);
    EXPECT_EQ(cached_result, result1);
}

// Test polynomial solver with known polynomial
TEST(HyperCatalanPolynomialSolverTest, SolveKnownPolynomial) {
    // Create a solver for degree 3 polynomial with 20 terms
    HyperCatalanPolynomialSolver solver(3, 20);
    
    // Define a polynomial with known roots: (x-1)(x-2)(x-3) = x³-6x²+11x-6
    std::vector<HighPrecisionFloat> coefficients = {
        -6,   // x⁰
        11,   // x¹
        -6,   // x²
        1     // x³
    };
    
    // Solve the polynomial
    HighPrecisionFloat root = solver.solvePolynomial(coefficients);
    
    // Refine using bootstrap with 10 iterations
    HighPrecisionFloat refined_root = solver.bootstrapRoot(coefficients, root, 10);
    
    // Verify the root is approximately one of the expected values (1, 2, or 3)
    bool is_valid_root = false;
    std::vector<HighPrecisionFloat> expected_roots = {1, 2, 3};
    
    for (const auto& expected : expected_roots) {
        if (abs(refined_root - expected) < 0.01) {
            is_valid_root = true;
            break;
        }
    }
    
    EXPECT_TRUE(is_valid_root) << "Expected root to be one of {1, 2, 3}, but got " << refined_root;
    
    // Verify the polynomial evaluates to approximately zero at the root
    HighPrecisionFloat polynomial_value = 0;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        polynomial_value += coefficients[i] * pow(refined_root, i);
    }
    
    EXPECT_LT(abs(polynomial_value), 1e-10) << "Polynomial should evaluate close to zero at root";
}

// Test bootstrap method with different initial guesses
TEST(HyperCatalanPolynomialSolverTest, BootstrapMethod) {
    HyperCatalanPolynomialSolver solver(3, 20);
    
    // Define a simple quadratic: x² - 4 = 0, with roots +2 and -2
    std::vector<HighPrecisionFloat> coefficients = {
        -4,   // x⁰
        0,    // x¹
        1     // x²
    };
    
    // Test with positive initial guess
    HighPrecisionFloat positive_guess = solver.bootstrapRoot(coefficients, 1.5, 10);
    EXPECT_NEAR(static_cast<double>(positive_guess), 2.0, 0.0001);
    
    // Test with negative initial guess
    HighPrecisionFloat negative_guess = solver.bootstrapRoot(coefficients, -1.5, 10);
    EXPECT_NEAR(static_cast<double>(negative_guess), -2.0, 0.0001);
}

// Test the solver with polynomials of different degrees
TEST(HyperCatalanPolynomialSolverTest, DifferentDegrees) {
    // Test linear equation: 2x - 4 = 0, root = 2
    {
        HyperCatalanPolynomialSolver solver(1, 20);
        std::vector<HighPrecisionFloat> coefficients = {-4, 2};
        HighPrecisionFloat root = solver.bootstrapRoot(coefficients, 1.0, 10);
        EXPECT_NEAR(static_cast<double>(root), 2.0, 0.0001);
    }
    
    // Test cubic equation: x³ - x = 0, roots = {-1, 0, 1}
    {
        HyperCatalanPolynomialSolver solver(3, 20);
        std::vector<HighPrecisionFloat> coefficients = {0, -1, 0, 1};
        
        // Try to find the root at 1
        HighPrecisionFloat root1 = solver.bootstrapRoot(coefficients, 0.5, 10);
        EXPECT_NEAR(static_cast<double>(root1), 1.0, 0.0001);
        
        // Try to find the root at -1
        HighPrecisionFloat root2 = solver.bootstrapRoot(coefficients, -0.5, 10);
        EXPECT_NEAR(static_cast<double>(root2), -1.0, 0.0001);
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
} 