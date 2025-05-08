#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <unordered_map>
#include <functional>

// Структура для хранения типа субдигона (m2, m3, m4, ...)
struct SubdigonType {
    std::vector<int> m;
    
    SubdigonType(const std::vector<int>& values) : m(values) {}
    
    bool operator==(const SubdigonType& other) const {
        return m == other.m;
    }
    
    // Вычисление количества граней
    int faces() const {
        int sum = 0;
        for (size_t i = 0; i < m.size(); ++i) {
            sum += m[i];
        }
        return sum;
    }
    
    // Вычисление количества рёбер
    int edges() const {
        int sum = 0;
        for (size_t i = 0; i < m.size(); ++i) {
            sum += (i + 2) * m[i];
        }
        return sum / 2;
    }
    
    // Вычисление количества вершин
    int vertices() const {
        return edges() - faces() + 2;
    }
};

// Хеш-функция для SubdigonType
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

// Класс для вычисления гипер-каталановских чисел
class HyperCatalanCalculator {
private:
    std::unordered_map<SubdigonType, double> cache;
    
    // Вычисление факториала
    double factorial(int n) const {
        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }
    
public:
    HyperCatalanCalculator() {}
    
    // Вычисление гипер-каталановского числа для заданного типа субдигона
    double calculate(const SubdigonType& type) {
        // Сначала проверяем кеш
        if (cache.find(type) != cache.end()) {
            return cache[type];
        }
        
        // Вычисляем количество рёбер: 2*m₂ + 3*m₃ + 4*m₄ + ...
        int e = 0;
        for (size_t i = 0; i < type.m.size(); ++i) {
            e += (i + 2) * type.m[i];
        }
        
        // Вычисляем количество вершин: 1 + m₂ + 2*m₃ + 3*m₄ + ...
        int v = 1;
        for (size_t i = 0; i < type.m.size(); ++i) {
            v += i * type.m[i];
        }
        
        // Вычисляем гипер-каталановское число по формуле из Теоремы 5
        double numerator = factorial(e);
        double denominator = factorial(v);
        
        for (size_t i = 0; i < type.m.size(); ++i) {
            denominator *= factorial(type.m[i]);
        }
        
        double result = numerator / denominator;
        
        // Сохраняем в кеше
        cache[type] = result;
        
        return result;
    }
};

// Класс для решения полиномиальных уравнений с использованием гипер-каталановских рядов
class HyperCatalanPolynomialSolver {
private:
    int max_degree;
    int max_terms;
    HyperCatalanCalculator calculator;
    
    // Генерация всех возможных типов субдигонов с заданным числом граней и максимальным размером полигона
    std::vector<std::vector<int>> generateTypes(int total_faces, int max_polygon_size) {
        std::vector<std::vector<int>> results;
        std::vector<int> current(max_polygon_size, 0);
        generateTypesRecursive(results, current, total_faces, 0, max_polygon_size);
        return results;
    }
    
    // Рекурсивный помощник для генерации типов субдигонов
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
    
    // Решение полинома в геометрической форме: 1 - a + t₂a² + t₃a³ + ... = 0
    double solveGeometricForm(const std::vector<double>& t_coefficients) {
        double result = 0.0;
        
        // Перебираем все возможные типы субдигонов до max_terms
        for (int total_faces = 0; total_faces < max_terms; ++total_faces) {
            std::vector<std::vector<int>> types = generateTypes(total_faces, max_degree - 1);
            
            for (const auto& type_vec : types) {
                SubdigonType type(type_vec);
                
                // Вычисляем гипер-каталановское число
                double C_m = calculator.calculate(type);
                
                // Вычисляем произведение t₂^m₂ · t₃^m₃ · t₄^m₄ · ...
                double term_product = 1.0;
                for (size_t i = 0; i < type.m.size(); ++i) {
                    if (type.m[i] > 0 && i + 2 < t_coefficients.size()) {
                        term_product *= std::pow(t_coefficients[i + 2], type.m[i]);
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
    
    // Решение общего полиномиального уравнения: c₀ + c₁x + c₂x² + ... = 0
    double solvePolynomial(const std::vector<double>& coefficients) {
        if (coefficients.size() < 2) {
            throw std::invalid_argument("Полином должен быть как минимум степени 1");
        }
        
        // Конвертируем в геометрическую форму: 1 - a + t₂a² + t₃a³ + ... = 0
        std::vector<double> geometric_coeffs;
        geometric_coeffs.push_back(1.0);  // Константа 1
        geometric_coeffs.push_back(-1.0); // Коэффициент при a¹
        
        for (size_t i = 2; i < coefficients.size(); ++i) {
            geometric_coeffs.push_back(coefficients[i] / coefficients[1]);
        }
        
        // Решаем с использованием гипер-каталановских рядов
        double root = solveGeometricForm(geometric_coeffs);
        
        // Конвертируем обратно в корень исходного полинома
        return -coefficients[0] / (coefficients[1] * root);
    }
    
    // Улучшение приближения корня с использованием метода Ньютона
    double bootstrapRoot(const std::vector<double>& coefficients, 
                        double initial_guess,
                        int iterations,
                        double epsilon = 1e-15) {
        // Создаем функцию полинома
        auto polynomial_function = [&coefficients](double x) {
            double result = 0.0;
            for (size_t i = 0; i < coefficients.size(); ++i) {
                result += coefficients[i] * std::pow(x, i);
            }
            return result;
        };
        
        // Создаем функцию производной
        auto derivative_function = [&coefficients](double x) {
            double result = 0.0;
            for (size_t i = 1; i < coefficients.size(); ++i) {
                result += i * coefficients[i] * std::pow(x, i - 1);
            }
            return result;
        };
        
        // Применяем метод Ньютона
        double x = initial_guess;
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
            
            if (std::fabs(delta) < epsilon) {
                break;
            }
        }
        
        return x;
    }
};

int main() {
    std::cout << std::setprecision(15);
    
    std::cout << "Решение полиномиальных уравнений с использованием гипер-каталановских рядов" << std::endl;
    std::cout << "=================================================================" << std::endl;
    
    // Тест 1: Квадратное уравнение x^2 - 4 = 0 (корни: -2, 2)
    {
        std::vector<double> coefficients = {-4, 0, 1};
        HyperCatalanPolynomialSolver solver(2, 10);
        
        try {
            // Решаем с использованием гипер-каталановских рядов
            double series_root = solver.solvePolynomial(coefficients);
            std::cout << "\nТест 1: x^2 - 4 = 0" << std::endl;
            std::cout << "Приближение корня из гипер-каталановского ряда: " << series_root << std::endl;
            
            // Уточняем корень с использованием метода bootstrap
            double bootstrap_root = solver.bootstrapRoot(coefficients, series_root, 10);
            std::cout << "Корень после улучшения: " << bootstrap_root << std::endl;
            
            // Проверяем погрешность
            double error = 0.0;
            for (size_t i = 0; i < coefficients.size(); ++i) {
                error += coefficients[i] * std::pow(bootstrap_root, i);
            }
            std::cout << "Погрешность: " << std::fabs(error) << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Ошибка: " << e.what() << std::endl;
        }
    }
    
    // Тест 2: Кубическое уравнение x^3 - 6x^2 + 11x - 6 = 0 (корни: 1, 2, 3)
    {
        std::vector<double> coefficients = {-6, 11, -6, 1};
        HyperCatalanPolynomialSolver solver(3, 15);
        
        try {
            // Решаем с использованием гипер-каталановских рядов
            double series_root = solver.solvePolynomial(coefficients);
            std::cout << "\nТест 2: x^3 - 6x^2 + 11x - 6 = 0" << std::endl;
            std::cout << "Приближение корня из гипер-каталановского ряда: " << series_root << std::endl;
            
            // Находим все корни с помощью bootstrap с разными начальными приближениями
            double bootstrap_root1 = solver.bootstrapRoot(coefficients, 0.8, 10);
            double bootstrap_root2 = solver.bootstrapRoot(coefficients, 2.2, 10);
            double bootstrap_root3 = solver.bootstrapRoot(coefficients, 3.2, 10);
            
            // Вычисляем погрешности
            auto evaluate = [&coefficients](double x) {
                double result = 0.0;
                for (size_t i = 0; i < coefficients.size(); ++i) {
                    result += coefficients[i] * std::pow(x, i);
                }
                return result;
            };
            
            double error1 = std::fabs(evaluate(bootstrap_root1));
            double error2 = std::fabs(evaluate(bootstrap_root2));
            double error3 = std::fabs(evaluate(bootstrap_root3));
            
            std::cout << "Корень 1: " << bootstrap_root1 << " (ожидается 1)" << std::endl;
            std::cout << "Погрешность: " << error1 << std::endl;
            
            std::cout << "Корень 2: " << bootstrap_root2 << " (ожидается 2)" << std::endl;
            std::cout << "Погрешность: " << error2 << std::endl;
            
            std::cout << "Корень 3: " << bootstrap_root3 << " (ожидается 3)" << std::endl;
            std::cout << "Погрешность: " << error3 << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Ошибка: " << e.what() << std::endl;
        }
    }
    
    // Тест 3: x^5 - x - 1 = 0 (уравнение из стандартного теста для методов нахождения корня)
    {
        std::vector<double> coefficients = {-1, -1, 0, 0, 0, 1};
        HyperCatalanPolynomialSolver solver(5, 20);
        
        try {
            // Решаем с использованием гипер-каталановских рядов
            double series_root = solver.solvePolynomial(coefficients);
            std::cout << "\nТест 3: x^5 - x - 1 = 0" << std::endl;
            std::cout << "Приближение корня из гипер-каталановского ряда: " << series_root << std::endl;
            
            // Уточняем корень с использованием bootstrap
            double bootstrap_root = solver.bootstrapRoot(coefficients, series_root, 15);
            std::cout << "Корень после улучшения: " << bootstrap_root << std::endl;
            
            // Проверяем погрешность
            double error = 0.0;
            for (size_t i = 0; i < coefficients.size(); ++i) {
                error += coefficients[i] * std::pow(bootstrap_root, i);
            }
            std::cout << "Погрешность: " << std::fabs(error) << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Ошибка: " << e.what() << std::endl;
        }
    }
    
    return 0;
} 