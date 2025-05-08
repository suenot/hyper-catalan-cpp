#include "hyper_catalan_solver.h"
#include <iostream>
#include <vector>
#include <string>

int main(int argc, char* argv[]) {
    std::cout << "Hyper-Catalan Series Polynomial Solver" << std::endl;
    std::cout << "Based on 'A Hyper-Catalan Series Solution to Polynomial Equations, and the Geode'" << std::endl;
    std::cout << "by N.J. Wildberger and K.W. Rubin" << std::endl;
    std::cout << "------------------------------------------------" << std::endl;
    
    try {
        // Get the degree of the polynomial
        int degree;
        std::cout << "Enter the degree of polynomial: ";
        std::cin >> degree;
        
        std::vector<HighPrecisionFloat> coefficients(degree + 1);
        std::cout << "Enter coefficients from c₀ to c" << degree << " (constant term first): " << std::endl;
        
        for (int i = 0; i <= degree; ++i) {
            double coeff;
            std::cout << "c" << i << ": ";
            std::cin >> coeff;
            coefficients[i] = HighPrecisionFloat(coeff);
        }
        
        HyperCatalanPolynomialSolver solver(degree, 20);
        
        // Solve using Hyper-Catalan series
        try {
            HighPrecisionFloat series_root = solver.solvePolynomial(coefficients);
            std::cout << "Root by Hyper-Catalan series: " << series_root << std::endl;
            
            // Get initial guess and iterations for bootstrap method
            HighPrecisionFloat initial_guess = series_root;
            std::cout << "Enter initial guess for bootstrap method (default: " << series_root << "): ";
            std::string input;
            std::cin.ignore();
            std::getline(std::cin, input);
            if (!input.empty()) {
                initial_guess = HighPrecisionFloat(input);
            }
            
            int iterations;
            std::cout << "Enter number of iterations for bootstrap method: ";
            std::cin >> iterations;
            
            HighPrecisionFloat bootstrap_root = solver.bootstrapRoot(coefficients, initial_guess, iterations);
            std::cout << "Root by bootstrap method: " << bootstrap_root << std::endl;
            
            // Check accuracy of solution
            HighPrecisionFloat bootstrap_error = 0;
            for (int i = 0; i <= degree; ++i) {
                bootstrap_error += coefficients[i] * pow(bootstrap_root, i);
            }
            std::cout << "Error: " << abs(bootstrap_error) << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error solving polynomial: " << e.what() << std::endl;
            return 1;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
} 