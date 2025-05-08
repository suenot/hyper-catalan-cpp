#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

// Простая реализация для тестирования погрешности
double solvePolynomial(const std::vector<double>& coefficients, double initial_guess, int iterations) {
    // Функция для вычисления значения полинома
    auto polynomial_function = [&coefficients](double x) {
        double result = 0.0;
        for (size_t i = 0; i < coefficients.size(); ++i) {
            result += coefficients[i] * std::pow(x, i);
        }
        return result;
    };
    
    // Функция для вычисления производной полинома
    auto derivative_function = [&coefficients](double x) {
        double result = 0.0;
        for (size_t i = 1; i < coefficients.size(); ++i) {
            result += i * coefficients[i] * std::pow(x, i - 1);
        }
        return result;
    };
    
    // Метод Ньютона для уточнения корня
    double x = initial_guess;
    const double epsilon = 1e-15;
    
    for (int i = 0; i < iterations; ++i) {
        double f_x = polynomial_function(x);
        if (std::fabs(f_x) < epsilon) {
            break;
        }
        
        double df_x = derivative_function(x);
        if (std::fabs(df_x) < epsilon) {
            break;
        }
        
        double delta = f_x / df_x;
        x = x - delta;
        
        // Добавляем вывод для отладки
        // std::cout << "Итерация " << i << ": x = " << x << ", f(x) = " << f_x << ", delta = " << delta << std::endl;
        
        if (std::fabs(delta) < epsilon) {
            break;
        }
    }
    
    return x;
}

// Функция для вычисления значения полинома
double evaluatePolynomial(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * std::pow(x, i);
    }
    return result;
}

int main() {
    std::cout << std::setprecision(15);
    
    std::cout << "Проверка погрешности для известных многочленов" << std::endl;
    std::cout << "==============================================" << std::endl;
    
    // Тест 1: x^2 - 4 = 0 (корни: -2, 2)
    {
        std::vector<double> coefficients = {-4, 0, 1};
        
        double root_positive = solvePolynomial(coefficients, 1.0, 10);
        double root_negative = solvePolynomial(coefficients, -1.0, 10);
        
        double error_positive = std::fabs(evaluatePolynomial(coefficients, root_positive));
        double error_negative = std::fabs(evaluatePolynomial(coefficients, root_negative));
        
        std::cout << "\nТест 1: x^2 - 4 = 0" << std::endl;
        std::cout << "Положительный корень: " << root_positive << " (ожидается 2)" << std::endl;
        std::cout << "Погрешность: " << error_positive << std::endl;
        
        std::cout << "Отрицательный корень: " << root_negative << " (ожидается -2)" << std::endl;
        std::cout << "Погрешность: " << error_negative << std::endl;
    }
    
    // Тест 2: x^3 - 6x^2 + 11x - 6 = 0 (корни: 1, 2, 3)
    {
        std::vector<double> coefficients = {-6, 11, -6, 1};
        
        // Используем более точные начальные приближения
        double root1 = solvePolynomial(coefficients, 0.8, 10);
        double root2 = solvePolynomial(coefficients, 2.2, 10);
        double root3 = solvePolynomial(coefficients, 3.2, 10);
        
        double error1 = std::fabs(evaluatePolynomial(coefficients, root1));
        double error2 = std::fabs(evaluatePolynomial(coefficients, root2));
        double error3 = std::fabs(evaluatePolynomial(coefficients, root3));
        
        std::cout << "\nТест 2: x^3 - 6x^2 + 11x - 6 = 0" << std::endl;
        std::cout << "Корень 1: " << root1 << " (ожидается 1)" << std::endl;
        std::cout << "Погрешность: " << error1 << std::endl;
        
        std::cout << "Корень 2: " << root2 << " (ожидается 2)" << std::endl;
        std::cout << "Погрешность: " << error2 << std::endl;
        
        std::cout << "Корень 3: " << root3 << " (ожидается 3)" << std::endl;
        std::cout << "Погрешность: " << error3 << std::endl;
    }
    
    // Тест 3: x^4 - 10x^2 + 9 = 0 (корни: -3, -1, 1, 3)
    {
        std::vector<double> coefficients = {9, 0, -10, 0, 1};
        
        double root1 = solvePolynomial(coefficients, -3.5, 10);
        double root2 = solvePolynomial(coefficients, -0.5, 10);
        double root3 = solvePolynomial(coefficients, 0.5, 10);
        double root4 = solvePolynomial(coefficients, 3.5, 10);
        
        double error1 = std::fabs(evaluatePolynomial(coefficients, root1));
        double error2 = std::fabs(evaluatePolynomial(coefficients, root2));
        double error3 = std::fabs(evaluatePolynomial(coefficients, root3));
        double error4 = std::fabs(evaluatePolynomial(coefficients, root4));
        
        std::cout << "\nТест 3: x^4 - 10x^2 + 9 = 0" << std::endl;
        std::cout << "Корень 1: " << root1 << " (ожидается -3)" << std::endl;
        std::cout << "Погрешность: " << error1 << std::endl;
        
        std::cout << "Корень 2: " << root2 << " (ожидается -1)" << std::endl;
        std::cout << "Погрешность: " << error2 << std::endl;
        
        std::cout << "Корень 3: " << root3 << " (ожидается 1)" << std::endl;
        std::cout << "Погрешность: " << error3 << std::endl;
        
        std::cout << "Корень 4: " << root4 << " (ожидается 3)" << std::endl;
        std::cout << "Погрешность: " << error4 << std::endl;
    }
    
    // Тест 4: x^5 - x - 1 = 0 (имеет единственный действительный корень около 1.16)
    {
        std::vector<double> coefficients = {-1, -1, 0, 0, 0, 1};
        
        double root = solvePolynomial(coefficients, 1.0, 15);
        
        double error = std::fabs(evaluatePolynomial(coefficients, root));
        
        std::cout << "\nТест 4: x^5 - x - 1 = 0" << std::endl;
        std::cout << "Корень: " << root << " (приблизительно 1.16)" << std::endl;
        std::cout << "Погрешность: " << error << std::endl;
    }
    
    // Тест 5: x^10 - 1 = 0 (корни на единичной окружности в комплексной плоскости)
    {
        std::vector<double> coefficients(11, 0);
        coefficients[0] = -1;
        coefficients[10] = 1;
        
        double root_pos = solvePolynomial(coefficients, 0.9, 20);
        double root_neg = solvePolynomial(coefficients, -0.9, 20);
        
        double error_pos = std::fabs(evaluatePolynomial(coefficients, root_pos));
        double error_neg = std::fabs(evaluatePolynomial(coefficients, root_neg));
        
        std::cout << "\nТест 5: x^10 - 1 = 0" << std::endl;
        std::cout << "Положительный корень: " << root_pos << " (ожидается 1)" << std::endl;
        std::cout << "Погрешность: " << error_pos << std::endl;
        
        std::cout << "Отрицательный корень: " << root_neg << " (ожидается -1)" << std::endl;
        std::cout << "Погрешность: " << error_neg << std::endl;
    }
    
    return 0;
} 