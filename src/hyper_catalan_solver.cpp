#include "hyper_catalan_solver.h"

// Implementation of functions (if needed)
// Most of the implementation is already in the header file

// Example usage function
HighPrecisionFloat solveExamplePolynomial() {
    // Create a solver for degree 5 polynomials with 20 max terms
    HyperCatalanPolynomialSolver solver(5, 20);
    
    // Define coefficients for a sample polynomial: x³ - 6x² + 11x - 6 = 0
    std::vector<HighPrecisionFloat> coefficients = {
        -6,   // x⁰
        11,   // x¹
        -6,   // x²
        1     // x³
    };
    
    // Solve the polynomial
    try {
        HighPrecisionFloat root = solver.solvePolynomial(coefficients);
        
        // Refine using Newton's method
        HighPrecisionFloat refined_root = solver.bootstrapRoot(coefficients, root, 10);
        
        return refined_root;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 0;
    }
} 