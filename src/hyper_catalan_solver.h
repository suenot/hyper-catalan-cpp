#pragma once

#include <vector>
#include <unordered_map>
#include <cmath>
#include <iostream>
#include <functional>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <Eigen/Dense>

// High precision floating point type
using HighPrecisionFloat = boost::multiprecision::cpp_dec_float_50;

// Structure for storing subdigon type (m2, m3, m4, ...)
struct SubdigonType {
    std::vector<int> m;
    
    SubdigonType(const std::vector<int>& values) : m(values) {}
    
    bool operator==(const SubdigonType& other) const {
        return m == other.m;
    }
    
    // Calculate number of faces
    int faces() const {
        int sum = 0;
        for (size_t i = 0; i < m.size(); ++i) {
            sum += m[i];
        }
        return sum;
    }
    
    // Calculate number of edges
    int edges() const {
        int sum = 0;
        for (size_t i = 0; i < m.size(); ++i) {
            sum += (i + 2) * m[i];
        }
        return sum / 2;
    }
    
    // Calculate number of vertices
    int vertices() const {
        return edges() - faces() + 2;
    }
};

// Hash function for SubdigonType
namespace std {
    template<>
    struct hash<SubdigonType> {
        size_t operator()(const SubdigonType& type) const {
            size_t hash_val = 0;
            for (size_t i = 0; i < type.m.size(); ++i) {
                hash_val = hash_val * 31 + std::hash<int>()(type.m[i]);
            }
            return hash_val;
        }
    };
}

// Class for calculating Hyper-Catalan numbers
class HyperCatalanCalculator {
private:
    std::unordered_map<SubdigonType, HighPrecisionFloat> cache;
    
    // Calculate factorial
    HighPrecisionFloat factorial(int n) {
        HighPrecisionFloat result = 1;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }
    
public:
    HyperCatalanCalculator() {}
    
    // Calculate Hyper-Catalan number for given subdigon type
    HighPrecisionFloat calculate(const SubdigonType& type) {
        // Check cache first
        if (cache.find(type) != cache.end()) {
            return cache[type];
        }
        
        // Calculate edges count: 2*m₂ + 3*m₃ + 4*m₄ + ...
        int e = 0;
        for (size_t i = 0; i < type.m.size(); ++i) {
            e += (i + 2) * type.m[i];
        }
        
        // Calculate vertices count: 1 + m₂ + 2*m₃ + 3*m₄ + ...
        int v = 1;
        for (size_t i = 0; i < type.m.size(); ++i) {
            v += i * type.m[i];
        }
        
        // Calculate Hyper-Catalan number using formula from Theorem 5
        HighPrecisionFloat numerator = factorial(e);
        HighPrecisionFloat denominator = factorial(v);
        
        for (size_t i = 0; i < type.m.size(); ++i) {
            denominator *= factorial(type.m[i]);
        }
        
        HighPrecisionFloat result = numerator / denominator;
        
        // Store in cache
        cache[type] = result;
        
        return result;
    }
};

// Class for solving polynomial equations using Hyper-Catalan series
class HyperCatalanPolynomialSolver {
private:
    int max_degree;
    int max_terms;
    HyperCatalanCalculator calculator;
    
    // Generate all possible subdigon types with given total faces and max polygon size
    std::vector<std::vector<int>> generateTypes(int total_faces, int max_polygon_size) {
        std::vector<std::vector<int>> results;
        std::vector<int> current(max_polygon_size, 0);
        generateTypesRecursive(results, current, total_faces, 0, max_polygon_size);
        return results;
    }
    
    // Recursive helper for generating subdigon types
    void generateTypesRecursive(std::vector<std::vector<int>>& results, 
                               std::vector<int>& current, 
                               int remaining_faces, 
                               int index, 
                               int max_polygon_size) {
        if (index == max_polygon_size) {
            if (remaining_faces == 0) {
                results.push_back(current);
            }
            return;
        }
        
        for (int i = 0; i <= remaining_faces; ++i) {
            current[index] = i;
            generateTypesRecursive(results, current, remaining_faces - i, index + 1, max_polygon_size);
        }
    }
    
    // Solve polynomial in geometric form: 1 - a + t₂a² + t₃a³ + ... = 0
    HighPrecisionFloat solveGeometricForm(const std::vector<HighPrecisionFloat>& t_coefficients) {
        HighPrecisionFloat result = 0;
        
        // Iterate through all possible types of subdigons up to max_terms
        for (int total_faces = 0; total_faces < max_terms; ++total_faces) {
            std::vector<std::vector<int>> types = generateTypes(total_faces, max_degree - 1);
            
            for (const auto& type_vec : types) {
                SubdigonType type(type_vec);
                
                // Calculate Hyper-Catalan number
                HighPrecisionFloat C_m = calculator.calculate(type);
                
                // Calculate product t₂^m₂ · t₃^m₃ · t₄^m₄ · ...
                HighPrecisionFloat term_product = 1;
                for (size_t i = 0; i < type.m.size(); ++i) {
                    if (type.m[i] > 0 && i + 2 < t_coefficients.size()) {
                        term_product *= pow(t_coefficients[i + 2], type.m[i]);
                    }
                }
                
                result += C_m * term_product;
            }
        }
        
        return result;
    }
    
public:
    HyperCatalanPolynomialSolver(int degree, int terms) 
        : max_degree(degree), max_terms(terms), calculator() {}
    
    // Solve general polynomial equation: c₀ + c₁x + c₂x² + ... = 0
    HighPrecisionFloat solvePolynomial(const std::vector<HighPrecisionFloat>& coefficients) {
        if (coefficients.size() < 2) {
            throw std::invalid_argument("Polynomial must be at least of degree 1");
        }
        
        // Convert to geometric form: 1 - a + t₂a² + t₃a³ + ... = 0
        std::vector<HighPrecisionFloat> geometric_coeffs;
        geometric_coeffs.push_back(1);  // Constant 1
        geometric_coeffs.push_back(-1); // Coefficient for a¹
        
        for (size_t i = 2; i < coefficients.size(); ++i) {
            geometric_coeffs.push_back(coefficients[i] / coefficients[1]);
        }
        
        // Solve using Hyper-Catalan series
        HighPrecisionFloat root = solveGeometricForm(geometric_coeffs);
        
        // Convert back to original polynomial root
        return -coefficients[0] / (coefficients[1] * root);
    }
    
    // Bootstrap root approximation using Newton's method
    HighPrecisionFloat bootstrapRoot(const std::vector<HighPrecisionFloat>& coefficients, 
                                    HighPrecisionFloat initial_guess,
                                    int iterations,
                                    HighPrecisionFloat epsilon = 1e-15) {
        // Create polynomial function
        auto polynomial_function = [&coefficients](HighPrecisionFloat x) {
            HighPrecisionFloat result = 0;
            for (size_t i = 0; i < coefficients.size(); ++i) {
                result += coefficients[i] * pow(x, i);
            }
            return result;
        };
        
        // Create derivative function
        auto derivative_function = [&coefficients](HighPrecisionFloat x) {
            HighPrecisionFloat result = 0;
            for (size_t i = 1; i < coefficients.size(); ++i) {
                result += i * coefficients[i] * pow(x, i - 1);
            }
            return result;
        };
        
        // Apply Newton's method
        HighPrecisionFloat x = initial_guess;
        for (int i = 0; i < iterations; ++i) {
            HighPrecisionFloat f_x = polynomial_function(x);
            if (abs(f_x) < epsilon) {
                break;
            }
            
            HighPrecisionFloat df_x = derivative_function(x);
            if (abs(df_x) < epsilon) {
                break;
            }
            
            x = x - f_x / df_x;
        }
        
        return x;
    }
}; 