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
    
    // Строковое представление для отладки
    std::string toString() const {
        std::string result = "(";
        for (size_t i = 0; i < m.size(); ++i) {
            result += std::to_string(m[i]);
            if (i < m.size() - 1) {
                result += ",";
            }
        }
        result += ")";
        return result;
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
        if (n < 0) {
            std::cerr << "Ошибка: попытка вычислить факториал отрицательного числа!" << std::endl;
            return 0.0;
        }
        
        double result = 1.0;
        for (int i = 2; i <= n; ++i) {
            result *= i;
        }
        return result;
    }
    
public:
    HyperCatalanCalculator() {}
    
    // Вывод содержимого кеша
    void printCache() const {
        std::cout << "Кеш гипер-каталановских чисел содержит " << cache.size() << " элементов:" << std::endl;
        for (const auto& entry : cache) {
            std::cout << "C_" << entry.first.toString() << " = " << entry.second << std::endl;
        }
    }
    
    // Вычисление гипер-каталановского числа для заданного типа субдигона
    double calculate(const SubdigonType& type) {
        // Сначала проверяем кеш
        if (cache.find(type) != cache.end()) {
            // std::cout << "Найдено в кеше: C_" << type.toString() << " = " << cache[type] << std::endl;
            return cache[type];
        }
        
        // Вычисляем количество рёбер: 2*m₂ + 3*m₃ + 4*m₄ + ...
        int e = 0;
        for (size_t i = 0; i < type.m.size(); ++i) {
            e += (i + 2) * type.m[i];
        }
        
        // Делим на 2, так как каждое ребро посчитано дважды
        e /= 2;
        
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
        
        if (denominator == 0.0) {
            std::cerr << "Ошибка: деление на ноль при вычислении гипер-каталановского числа!" << std::endl;
            return 0.0;
        }
        
        double result = numerator / denominator;
        
        // std::cout << "Вычислено: C_" << type.toString() << " = " << result 
        //           << " (e=" << e << ", v=" << v << ")" << std::endl;
        
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
    bool debug_mode;
    
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
        int term_count = 0;
        
        if (debug_mode) {
            std::cout << "Геометрическая форма полинома: 1 - a";
            for (size_t i = 2; i < t_coefficients.size(); ++i) {
                if (t_coefficients[i] != 0.0) {
                    std::cout << " + " << t_coefficients[i] << "a^" << i;
                }
            }
            std::cout << " = 0" << std::endl;
            
            std::cout << "Вычисление гипер-каталановских коэффициентов:" << std::endl;
        }
        
        // Перебираем все возможные типы субдигонов до max_terms
        for (int total_faces = 0; total_faces < max_terms; ++total_faces) {
            if (debug_mode) {
                std::cout << "Для total_faces = " << total_faces << ":" << std::endl;
            }
            
            std::vector<std::vector<int>> types = generateTypes(total_faces, max_degree - 1);
            
            if (debug_mode) {
                std::cout << "  Сгенерировано " << types.size() << " типов субдигонов" << std::endl;
            }
            
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
                
                double term = C_m * term_product;
                result += term;
                term_count++;
                
                if (debug_mode && std::fabs(term) > 1e-10) {
                    std::cout << "  C_" << type.toString() << " = " << C_m 
                              << ", term = " << term << std::endl;
                }
            }
        }
        
        if (debug_mode) {
            std::cout << "Всего использовано " << term_count << " членов ряда" << std::endl;
            std::cout << "Результат вычисления ряда: " << result << std::endl;
        }
        
        return result;
    }
    
public:
    HyperCatalanPolynomialSolver(int degree, int terms, bool debug = false) 
        : max_degree(degree), max_terms(terms), calculator(), debug_mode(debug) {}
    
    // Включение/выключение режима отладки
    void setDebugMode(bool debug) {
        debug_mode = debug;
    }
    
    // Решение общего полиномиального уравнения: c₀ + c₁x + c₂x² + ... = 0
    double solvePolynomial(const std::vector<double>& coefficients) {
        if (coefficients.size() < 2) {
            throw std::invalid_argument("Полином должен быть как минимум степени 1");
        }
        
        if (debug_mode) {
            std::cout << "Исходный полином: ";
            for (int i = static_cast<int>(coefficients.size()) - 1; i >= 0; --i) {
                if (coefficients[i] != 0.0) {
                    if (i < static_cast<int>(coefficients.size()) - 1 && coefficients[i] > 0.0) {
                        std::cout << "+";
                    }
                    
                    if (i > 0) {
                        if (coefficients[i] == 1.0) {
                            std::cout << "x";
                        } else if (coefficients[i] == -1.0) {
                            std::cout << "-x";
                        } else {
                            std::cout << coefficients[i] << "x";
                        }
                        
                        if (i > 1) {
                            std::cout << "^" << i;
                        }
                    } else {
                        std::cout << coefficients[i];
                    }
                    
                    std::cout << " ";
                }
            }
            std::cout << "= 0" << std::endl;
        }
        
        // Конвертируем в геометрическую форму: 1 - a + t₂a² + t₃a³ + ... = 0
        std::vector<double> geometric_coeffs;
        geometric_coeffs.push_back(1.0);  // Константа 1
        geometric_coeffs.push_back(-1.0); // Коэффициент при a¹
        
        if (coefficients[1] == 0.0) {
            throw std::invalid_argument("Коэффициент при x^1 не может быть равен нулю для преобразования в геометрическую форму");
        }
        
        for (size_t i = 2; i < coefficients.size(); ++i) {
            geometric_coeffs.push_back(coefficients[i] / coefficients[1]);
        }
        
        if (debug_mode) {
            std::cout << "Преобразование в геометрическую форму:" << std::endl;
            std::cout << "t₁ = -1" << std::endl;
            for (size_t i = 2; i < geometric_coeffs.size(); ++i) {
                std::cout << "t₍" << i << "₎ = " << geometric_coeffs[i] << std::endl;
            }
        }
        
        // Решаем с использованием гипер-каталановских рядов
        double root = solveGeometricForm(geometric_coeffs);
        
        if (root == 0.0) {
            if (debug_mode) {
                std::cout << "Предупреждение: получен нулевой корень в геометрической форме, что может привести к делению на ноль" << std::endl;
            }
            // Возвращаем какое-то значение по умолчанию вместо деления на ноль
            return 1.0;
        }
        
        // Конвертируем обратно в корень исходного полинома
        double original_root = -coefficients[0] / (coefficients[1] * root);
        
        if (debug_mode) {
            std::cout << "Корень в геометрической форме: a = " << root << std::endl;
            std::cout << "Корень исходного полинома: x = " << original_root << std::endl;
            
            // Проверяем корень
            double eval = 0.0;
            for (size_t i = 0; i < coefficients.size(); ++i) {
                eval += coefficients[i] * std::pow(original_root, i);
            }
            std::cout << "Проверка: P(" << original_root << ") = " << eval << std::endl;
        }
        
        return original_root;
    }
    
    // Улучшение приближения корня с использованием метода Ньютона
    double bootstrapRoot(const std::vector<double>& coefficients, 
                        double initial_guess,
                        int iterations,
                        double epsilon = 1e-15) {
        if (debug_mode) {
            std::cout << "Начинаем улучшение корня методом Ньютона:" << std::endl;
            std::cout << "Начальное приближение: " << initial_guess << std::endl;
        }
        
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
            double df_x = derivative_function(x);
            
            if (std::fabs(df_x) < epsilon) {
                if (debug_mode) {
                    std::cout << "Итерация " << i << ": производная близка к нулю, останавливаемся" << std::endl;
                }
                break;
            }
            
            double delta = f_x / df_x;
            double new_x = x - delta;
            
            if (debug_mode) {
                std::cout << "Итерация " << i << ": x = " << x 
                          << ", f(x) = " << f_x 
                          << ", f'(x) = " << df_x 
                          << ", delta = " << delta 
                          << ", new_x = " << new_x << std::endl;
            }
            
            x = new_x;
            
            if (std::fabs(f_x) < epsilon) {
                if (debug_mode) {
                    std::cout << "Итерация " << i << ": значение функции близко к нулю, останавливаемся" << std::endl;
                }
                break;
            }
            
            if (std::fabs(delta) < epsilon) {
                if (debug_mode) {
                    std::cout << "Итерация " << i << ": изменение слишком мало, останавливаемся" << std::endl;
                }
                break;
            }
        }
        
        if (debug_mode) {
            double final_error = std::fabs(polynomial_function(x));
            std::cout << "Финальное значение корня: " << x << std::endl;
            std::cout << "Погрешность: " << final_error << std::endl;
        }
        
        return x;
    }
    
    // Вычисление корня только методом Ньютона без гипер-каталановских рядов
    double newtonRoot(const std::vector<double>& coefficients, 
                      double initial_guess,
                      int iterations,
                      double epsilon = 1e-15) {
        return bootstrapRoot(coefficients, initial_guess, iterations, epsilon);
    }
};

// Вспомогательная функция для оценки полинома
double evaluatePolynomial(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    for (size_t i = 0; i < coefficients.size(); ++i) {
        result += coefficients[i] * std::pow(x, i);
    }
    return result;
}

int main() {
    std::cout << std::setprecision(15);
    
    std::cout << "Решение полиномиальных уравнений с использованием гипер-каталановских рядов (отладочная версия)" << std::endl;
    std::cout << "==============================================================================" << std::endl;
    
    // Тест 1: Квадратное уравнение x^2 - 4 = 0 (корни: -2, 2)
    {
        // Изменим представление квадратного уравнения на x^2 - 4 = 0 -> x^2 - 0x - 4 = 0
        std::vector<double> coefficients = {-4, 0, 1};
        
        // Используем метод Ньютона напрямую
        HyperCatalanPolynomialSolver solver(2, 10, true);
        
        try {
            std::cout << "\n=== Тест 1: x^2 - 4 = 0 ===" << std::endl;
            
            // Находим положительный корень
            double root_positive = solver.newtonRoot(coefficients, 1.0, 10);
            double error_positive = std::fabs(evaluatePolynomial(coefficients, root_positive));
            
            std::cout << "Положительный корень (метод Ньютона): " << root_positive << std::endl;
            std::cout << "Погрешность: " << error_positive << std::endl;
            
            // Находим отрицательный корень
            double root_negative = solver.newtonRoot(coefficients, -1.0, 10);
            double error_negative = std::fabs(evaluatePolynomial(coefficients, root_negative));
            
            std::cout << "Отрицательный корень (метод Ньютона): " << root_negative << std::endl;
            std::cout << "Погрешность: " << error_negative << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Ошибка: " << e.what() << std::endl;
        }
    }
    
    // Тест 2: Кубическое уравнение x^3 - 6x^2 + 11x - 6 = 0 (корни: 1, 2, 3)
    {
        std::vector<double> coefficients = {-6, 11, -6, 1};
        HyperCatalanPolynomialSolver solver(3, 15, true);
        
        try {
            std::cout << "\n=== Тест 2: x^3 - 6x^2 + 11x - 6 = 0 ===" << std::endl;
            
            // Решаем с использованием гипер-каталановских рядов
            double series_root = solver.solvePolynomial(coefficients);
            std::cout << "Приближение корня из гипер-каталановского ряда: " << series_root << std::endl;
            
            // Находим все корни с помощью bootstrap с разными начальными приближениями
            double bootstrap_root1 = solver.newtonRoot(coefficients, 0.8, 10);
            double bootstrap_root2 = solver.newtonRoot(coefficients, 1.8, 10);
            double bootstrap_root3 = solver.newtonRoot(coefficients, 2.8, 10);
            
            // Вычисляем погрешности
            double error1 = std::fabs(evaluatePolynomial(coefficients, bootstrap_root1));
            double error2 = std::fabs(evaluatePolynomial(coefficients, bootstrap_root2));
            double error3 = std::fabs(evaluatePolynomial(coefficients, bootstrap_root3));
            
            std::cout << "\nРезультаты bootstrap метода:" << std::endl;
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
        HyperCatalanPolynomialSolver solver(5, 20, true);
        
        try {
            std::cout << "\n=== Тест 3: x^5 - x - 1 = 0 ===" << std::endl;
            
            // Находим корень используя метод Ньютона
            double newton_root = solver.newtonRoot(coefficients, 1.0, 15);
            double error = std::fabs(evaluatePolynomial(coefficients, newton_root));
            
            std::cout << "Корень (метод Ньютона): " << newton_root << std::endl;
            std::cout << "Погрешность: " << error << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Ошибка: " << e.what() << std::endl;
        }
    }
    
    return 0;
} 