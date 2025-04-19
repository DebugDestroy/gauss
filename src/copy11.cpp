/*Здраствуйте! В этой версии есть недоточеты! Пока непонятно 
1) что означает радиус из конфиг файла ( типо он круглый но плоский или это мячик)
2) Нужно возможно поставить модули в триангуляции и построении пути на p->field
3)нужно убрать высоту в триангуляции (она еще может быть отрицательной)
4) Доработать версию с наклоном (считаем наклон тележки)

5) Файлы нужно сохранять в правильных местах

6) Улучшить логирование
7) Если добавятся новые команды меняем readme file и help
8) Добавить критический угол наклона по бокам и вверх/вниз
9) Что за путаница с уровнем 127
10) Почему нужно делать срез только для горок ( все что выше slice горка) но ведь есть и ямки ( все что ниже slice ямка)
11) Вопрос о маршруте по ребрам Воронова или же через середины ребер между центрами
12) Фактчески kmeans не нужен
13) Точки А и Б добавлены в триангуляцию
14) Условие с радиусом странное.Если она круглая то как мне проверить что она застряла

15) bmp write строит файл с именем ( дать ему параметр)

16) Порог для шума в конфиг или в новую команду wave? + wave в readme

17) Есть ненужные переменные, функции и тд

Изменения: убрал неявное использование команд в Control + добавил параметры для Gnuplot и bmp write в командный файл + написал вопросы и ключевые моменты + улучшил help.txt

        ПРОГРАММА ЕЩЕ НЕ ДО КОНЦА ГОТОВА!!! МЕТОДЫ ТРИАНГУЛЯЦИИ И ПОСТРОЕНИЯ ПУТИ НЕ ГОТОВЫ!!! Следите за обновлениями на гитхаб [GitHub Profile](https://github.com/DebugDestroy)
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath> // Для std::hypot, std::atan, M_PI
#include <string>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <algorithm> // Для std::sample
#include <random>    // Для std::mt19937 и std::random_device
#include <queue>       // Для std::priority_queue
#include <unordered_map> // Для std::unordered_map
#include <valarray>    // Для std::valarray (если используется)

struct DispatcherParams {//параметры для диспетчера
   std::string s;//строка с командой
    int A; // размер поля
    int B; // размер поля
    double h;//для гаусса
    double x;//для гаусса
    double y;//для гаусса
    double sx;//для гаусса
    double sy;//для гаусса
    std::string filename;// сохраняет имя файла
    int slice;//срез для горок ( все что выше slice горка) но ведь есть и ямки ( все что ниже slice ямка)
    bool is_binary_image;//false (0) - изображение бинарное, true (1) - изображение полное
    int noisy;//Уровень шума (все компоненты меньше - шум)
    int k;//количество кластеров для метода k_means
    int kk;//количество кластеров для метода k_means_with_kern
};

// Структура для точки с вещественными координатами (для точности) auto clusterCenters = getClusterCenters();
struct PointD {
    double x, y;
    PointD(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    
    bool operator==(const PointD& other) const {
        return std::fabs(x - other.x) < 1e-6 && std::fabs(y - other.y) < 1e-6;
    }
    
    bool operator!=(const PointD& other) const {
        return !(*this == other);
    }
};

struct Triangle {
    PointD a, b, c;
    double minElevation;  // Минимальная высота в треугольнике
    Triangle(PointD a_, PointD b_, PointD c_) : a(a_), b(b_), c(c_) {}

    // Добавляем оператор сравнения
    bool operator==(const Triangle& other) const {
        return (a == other.a && b == other.b && c == other.c) ||
               (a == other.a && b == other.c && c == other.b) ||
               (a == other.b && b == other.a && c == other.c) ||
               (a == other.b && b == other.c && c == other.a) ||
               (a == other.c && b == other.a && c == other.b) ||
               (a == other.c && b == other.b && c == other.a);
    }
    bool isDelaunay(const std::vector<PointD>& allPoints) const {
        PointD circumCenter = calculateCircumcenter();
        double radius = distance(circumCenter, a);
        for (const auto& p : allPoints) {
            if (p != a && p != b && p != c && distance(p, circumCenter) < radius - 1e-6) {
                return false;
            }
        }
        return true;
    }
     PointD calculateCircumcenter() const {
        // Реализация вычисления центра окружности
        double d = 2 * (a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y));
        if (std::fabs(d) < 1e-9) return PointD(); // Защита от деления на ноль
        double x = ((a.x*a.x + a.y*a.y)*(b.y - c.y) + (b.x*b.x + b.y*b.y)*(c.y - a.y) + (c.x*c.x + c.y*c.y)*(a.y - b.y)) / d;
        double y = ((a.x*a.x + a.y*a.y)*(c.x - b.x) + (b.x*b.x + b.y*b.y)*(a.x - c.x) + (c.x*c.x + c.y*c.y)*(b.x - a.x)) / d;
        return PointD(x, y);
    }

    static double distance(const PointD& p1, const PointD& p2) {
        return std::hypot(p1.x - p2.x, p1.y - p2.y);
    }
};

struct Edge {
    PointD a, b;
    Edge(PointD a_, PointD b_) : a(a_), b(b_) {}
    
    bool operator==(const Edge& other) const {
        return (a == other.a && b == other.b) || (a == other.b && b == other.a);
    }
};

struct VoronoiEdge {
    PointD start, end;
    bool isInfinite = false; // Флаг бесконечности (как же круто звучит!!!!!!!!) 
    VoronoiEdge(PointD s, PointD e, bool infinite = false) : start(s), end(e), isInfinite(infinite) {}
};

// Структура для хранения данных ячейки сетки
struct GridCell {
    PointD position;       // Координаты центра ячейки
    double slopeAngle = 0; // Угол наклона в градусах
    bool visited = false;  // Флаг посещения
};

struct AStarNode {
    PointD position;       // Позиция узла (центр треугольника)
    double g;              // Стоимость пути от старта
    double h;              // Эвристическая оценка до цели
    double f;              // Общая стоимость: f = g + h
    AStarNode* parent;     // Родительский узел
    const Triangle* tri;   // Связанный треугольник

    AStarNode(PointD pos, const Triangle* triangle, AStarNode* p = nullptr)
        : position(pos), tri(triangle), parent(p), g(0), h(0), f(0) {}

    // Для сравнения в priority_queue
    bool operator<(const AStarNode& other) const { return f > other.f; }
};

//Классы для логирования
class Config {
public:
    int fieldWidth;
    int fieldHeight;
    double defaultX, defaultY, defaultSx, defaultSy, defaultH;
    std::string logFileNameInterface;
    std::string logFileNameControl;
    bool loggingInterfaceEnabled;
    bool loggingControlEnabled;
    double pointA_x, pointA_y, pointB_x, pointB_y;
    double vehicleRadius;

    Config(const std::string& filename) {
        
        std::ifstream configFile(filename);
        if (!configFile.is_open()) {
            std::cerr << "Failed to open config file." << std::endl;
            return;
        }

     std::string key;
while (configFile >> key) { // Считываем ключи
    if (key == "fieldWidth") configFile >> fieldWidth;
    else if (key == "fieldHeight") configFile >> fieldHeight;
    else if (key == "defaultX") configFile >> defaultX;
    else if (key == "defaultY") configFile >> defaultY;
    else if (key == "defaultSx") configFile >> defaultSx;
    else if (key == "defaultSy") configFile >> defaultSy;
    else if (key == "defaultH") configFile >> defaultH;
    else if (key == "logFileNameInterface") configFile >> logFileNameInterface;
    else if (key == "logFileNameControl") configFile >> logFileNameControl;
    else if (key == "loggingInterfaceEnabled") configFile >> std::boolalpha >> loggingInterfaceEnabled;
    else if (key == "loggingControlEnabled") configFile >> std::boolalpha >> loggingControlEnabled;
    else if (key == "pointA_x") configFile >> pointA_x;
            else if (key == "pointA_y") configFile >> pointA_y;
            else if (key == "pointB_x") configFile >> pointB_x;
            else if (key == "pointB_y") configFile >> pointB_y;
            else if (key == "vehicleRadius") configFile >> vehicleRadius;
}

        configFile.close();
    }
};

class Logger {
private:
    std::ofstream logFile;

public:
    Logger(const std::string& fileName) {
        if (!logFile.is_open()) {
            logFile.open(fileName, std::ios::out | std::ios::app);
        }
        
    }

    ~Logger() {
        if (logFile.is_open()) {
            logFile.close();
        }
    }

    void logMessage(const std::string& message, bool b) {
        if (logFile.is_open() && b) {
            auto now = std::chrono::system_clock::now();
            std::time_t now_c = std::chrono::system_clock::to_time_t(now);
            std::stringstream timeStamp;
            timeStamp << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S");
            logFile << "[" << timeStamp.str() << "] " << message << std::endl;
        }
    }
    
};

struct Color{ // структура для хранения rgb цвета
    uint8_t r, g, b;
};
    
//Классы для построения гаусса
class Gaus { // Класс, хранящий Гауссы
public:
    double h, x0, y0, sigma_x, sigma_y;
    Gaus(double h, double x0, double y0, double sigma_x, double sigma_y)
        : h(h), x0(x0), y0(y0), sigma_x(sigma_x), sigma_y(sigma_y) {}

    double calculate(int x, int y) const { //Функция вычисления Гаусса
    double dx = x - x0;
    double dy = y - y0;
    double denom_x = 2 * sigma_x * sigma_x;
    double denom_y = 2 * sigma_y * sigma_y;

    // Проверка на ноль
    if (denom_x == 0 || denom_y == 0) {
        std::cerr << "Error: sigma_x or sigma_y is zero!" << std::endl;
        return 0;
    }

    double exponent = -((dx * dx) / denom_x + (dy * dy) / denom_y);

    if (exponent < -700) {
        return 0;
    } else if (exponent > 700) {
        return std::numeric_limits<double>::infinity();
    }

    return h * exp(exponent);
}

};

class Pole {
public:
    std::vector<std::vector<double>> field;

    Pole(int A, int B) {
        field.resize(A, std::vector<double>(B, 0));
    }
    
     void resize(int A, int B) {
        field.resize(A); // Изменяем количество строк
        for (auto& row : field) {
            row.resize(B, 0); // Изменяем количество столбцов и инициализируем новыми значениями
        }
    }
};

class Component {
public:
    std::vector<std::vector<double>> componenta;
    int min_x, min_y, max_x, max_y;
    double center_x, center_y;
    double eigenvec1_x, eigenvec1_y;
    double eigenvec2_x, eigenvec2_y;
    double eigenvalue1, eigenvalue2;

   // Основной конструктор (для готовой матрицы)
    Component(const std::vector<std::vector<double>>& inputComponenta) 
        : componenta(inputComponenta) {
        calculate_metadata();
    }

    void calculate_metadata() {
        // Находим границы и центр
        min_x = componenta[0].size();
        min_y = componenta.size();
        max_x = 0;
        max_y = 0;
        double sum_x = 0, sum_y = 0;
        int count = 0;

        for (int i = 0; i < componenta.size(); ++i) {
            for (int j = 0; j < componenta[0].size(); ++j) {
                if (componenta[i][j] == 255) { // Точка принадлежит компоненте
                    min_x = std::min(min_x, j);
                    min_y = std::min(min_y, i);
                    max_x = std::max(max_x, j);
                    max_y = std::max(max_y, i);
                    sum_x += j;
                    sum_y += i;
                    count++;
                }
            }
        }
         if (count == 0) {
        std::cerr << "Error: Component has no valid pixels!" << std::endl;
        center_x = 0;
        center_y = 0;
        return;
    }
        center_x = sum_x / count;
        center_y = sum_y / count;

        // Расчет ковариационной матрицы
        double cov_xx = 0, cov_xy = 0, cov_yy = 0;
        for (int i = min_x; i <= max_x; ++i) {
            for (int j = min_y; j <= max_y; ++j) {
                if (componenta[i][j] == 0) {
                    double dx = i - center_x;
                    double dy = j - center_y;
                    cov_xx += dx * dx;
                    cov_xy += dx * dy;
                    cov_yy += dy * dy;
                }
            }
        }
        
        cov_xx /= count;
        cov_xy /= count;
        cov_yy /= count;

        // Собственные значения и векторы
        double trace = cov_xx + cov_yy;
        double det = cov_xx * cov_yy - cov_xy * cov_xy;
        eigenvalue1 = (trace + sqrt(trace * trace - 4 * det)) / 2;
        eigenvalue2 = (trace - sqrt(trace * trace - 4 * det)) / 2;

        // Собственные векторы (упрощенный расчет)
        if (cov_xy != 0) {
            eigenvec1_x = eigenvalue1 - cov_yy;
            eigenvec1_y = cov_xy;
            eigenvec2_x = eigenvalue2 - cov_yy;
            eigenvec2_y = cov_xy;
        } else {
            eigenvec1_x = 1; eigenvec1_y = 0;
            eigenvec2_x = 0; eigenvec2_y = 1;
        }

        // Вывод результатов
        std::cout << "Eigenvalues: " << eigenvalue1 << ", " << eigenvalue2 << std::endl;
        std::cout << "Eigenvector 1: (" << eigenvec1_x << ", " << eigenvec1_y << ")" << std::endl;
        std::cout << "Eigenvector 2: (" << eigenvec2_x << ", " << eigenvec2_y << ")" << std::endl;
    }
};

class Copier {
public:
    void removeNoise(std::vector<std::vector<double>>& field, 
                    const std::vector<Component>& components) {
        // Обнуляем поле
        for (auto& row : field) std::fill(row.begin(), row.end(), 0.0);
        
        // Копируем только значимые компоненты
        for (const auto& comp : components) {
            const auto& compData = comp.componenta;
            for (size_t i = 0; i < compData.size(); ++i) {
                for (size_t j = 0; j < compData[i].size(); ++j) {
                    if (compData[i][j] == 255) {
                        field[i][j] = 255;
                    }
                }
            }
        }
    }
};

class ColorGenerator {
public:
    // Метод для генерации уникальных цветов для заданного количества кластеров
    static std::vector<std::array<int, 3>> generateColors(int numColors) {
        if (numColors <= 0) {
            throw std::invalid_argument("Number of colors must be greater than 0.");
        }

        std::vector<std::array<int, 3>> colors;
        srand(static_cast<unsigned int>(time(0))); // Инициализация генератора случайных чисел

        for (int i = 0; i < numColors; ++i) {
            // Генерация оттенка (Hue) от 0 до 360
            double hue = static_cast<double>(i) / numColors * 360.0 + rand() % 30; // Добавляем случайность
            hue = fmod(hue, 360.0); // Убедимся, что hue остается в диапазоне [0, 360)

            // Преобразование HSL в RGB
            double r, g, b;
            double C = 1.0; // Чистота
            double X = C * (1 - fabs(fmod((hue / 60.0), 2) - 1));
            double m = 0.0; // Смещение

            if (hue < 60) {
                r = C; g = X; b = 0;
            } else if (hue < 120) {
                r = X; g = C; b = 0;
            } else if (hue < 180) {
                r = 0; g = C; b = X;
            } else if (hue < 240) {
                r = 0; g = X; b = C;
            } else if (hue < 300) {
                r = X; g = 0; b = C;
            } else {
                r = C; g = 0; b = X;
            }

            // Преобразование в диапазон 50-255
            int rVal = static_cast<int>(std::max(50.0, (r + m) * 255));
            int gVal = static_cast<int>(std::max(50.0, (g + m) * 255));
            int bVal = static_cast<int>(std::max(50.0, (b + m) * 255));

            // Учитываем, что значение не должно превышать 255
            rVal = std::min(rVal, 255);
            gVal = std::min(gVal, 255);
            bVal = std::min(bVal, 255);
            colors.push_back({rVal, gVal, bVal});
        }
        return colors;
    }
};
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class GaussBuilder {
   public:
   
       void addgauss(double h, double x0, double y0, double sigma_x, double sigma_y, std::vector<Gaus>& gaussi) {
        gaussi.emplace_back(h, x0, y0, sigma_x, sigma_y);
        //logger.logMessage("Added gauss", b);
    }
    
    void init(int A, int B, std::unique_ptr<Pole>& p) {
         if (!p) { // Проверяем, инициализирован ли указатель
        p = std::make_unique<Pole>(A, B); // Создаём новый объект только если p не инициализирован
    } else {
        // Возможно, вы хотите обновить существующий объект
        p->resize(A, B); // обновить размеры существующего объекта
    }
        //logger.logMessage("Added field", b);
    }
    
    void generate(std::unique_ptr<Pole>& p, std::vector<Gaus>& gaussi) {
    if (p == nullptr) {
       std::cout << "Pole not initialized!" << std::endl;
       return;
   }
   
   for (auto& row : p->field) {
        std::fill(row.begin(), row.end(), 0); // Уровень равнины
    }
        double value; 
        for (const auto& g : gaussi) {
    for (long unsigned int x = 0; x < p->field[0].size(); ++x) {
        for (long unsigned int y = 0; y < p->field.size(); ++y) {
            value = g.h * exp(-((pow((x - g.x0) / g.sigma_x, 2)) + (pow((y - g.y0) / g.sigma_y, 2))) / 2); 
            p->field[y][x] += value; 
           // Ограничиваем значения между -255 и 255
                p->field[y][x] = std::clamp(p->field[y][x], -255.0, 255.0);
        }
    }
}

    }
};
    
class BmpHandler {
   public:

       void bmp_write(const std::vector<std::vector<double>>& pixelMatrix, const std::string& filename) {
    int width = pixelMatrix[0].size();
    int height = pixelMatrix.size();
    int padding = (4 - (width * 3) % 4) % 4; // Padding for alignment to 4 bytes
    std::ofstream bmpFile(filename, std::ios::binary);
    if (!bmpFile) {
        std::cerr << "Failed to create BMP file." << std::endl;
        //logger.logMessage("Failed to create BMP file.", b);
        return;
    }
    // Write BMP header
    unsigned char bmpHeader[54] = {
        'B', 'M', // Identifier
        0, 0, 0, 0, // Size of file (will be set later)
        0, 0, 0, 0, // Reserved
        54, 0, 0, 0, // Header size
        40, 0, 0, 0, // Info header size
        0, 0, 0, 0, // Width (will be set later)
        0, 0, 0, 0, // Height (will be set later)
        1, 0, // Number of color planes
        24, 0, // Bits per pixel
        0, 0, 0, 0, // Compression
        0, 0, 0, 0, // Image size (will be set later)
        0x13, 0x0B, 0, 0, // Horizontal resolution
        0x13, 0x0B, 0, 0, // Vertical resolution
        0, 0, 0, 0, // Number of colors in palette
        0, 0, 0, 0  // Important colors
    };
    // Set width and height in header
    bmpHeader[18] = (width & 0xFF);
    bmpHeader[19] = (width >> 8) & 0xFF;
    bmpHeader[20] = (width >> 16) & 0xFF;
    bmpHeader[21] = (width >> 24) & 0xFF;
    bmpHeader[22] = (height & 0xFF);
    bmpHeader[23] = (height >> 8) & 0xFF;
    bmpHeader[24] = (height >> 16) & 0xFF;
    bmpHeader[25] = (height >> 24) & 0xFF;
    // Write header
    bmpFile.write(reinterpret_cast<char*>(bmpHeader), 54);
    // Write pixel data
    for (int y = height - 1; y >= 0; --y) { // BMP stores pixels bottom-to-top
        for (int x = 0; x < width; ++x) {
            unsigned char color = static_cast<unsigned char>(std::fabs(pixelMatrix[y][x])); // Color
            bmpFile.put(color); // B
            bmpFile.put(color); // G
            bmpFile.put(color); // R
        }
        // Add padding
        for (int p = 0; p < padding; ++p) {
            bmpFile.put(0);
        }
    }
    bmpFile.close();
}
    
   void bmp_read(GaussBuilder& gaussBuilder, const std::string &filename, std::vector<std::vector<double>> &pixelMatrix, std::unique_ptr<Pole>& p) {
    std::ifstream bmpFile(filename, std::ios::binary);
    if (!bmpFile) {
        std::cerr << "Failed to open BMP file." << std::endl;
        //logger.logMessage("Failed to open BMP file.", b);
        return;
    }
    // Читаем заголовок BMP
    unsigned char header[54];
    bmpFile.read(reinterpret_cast<char*>(header), 54);
    // Получаем ширину и высоту изображения
    int width = header[18] | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
    int height = header[22] | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);
    // Инициализируем новое поле

    gaussBuilder.init(height, width, p); // Заметь, что BMP хранит данные от нижней строки к верхней.
   
    
    // Инициализируем матрицу пикселей
    pixelMatrix.resize(height, std::vector<double>(width));
    
    // Читаем данные пикселей
    for (int y = height - 1; y >= 0; --y) { // BMP хранит данные снизу вверх
        for (int x = 0; x < width; ++x) {
            unsigned char color = bmpFile.get(); // Читаем B
            bmpFile.get(); // Читаем G
            bmpFile.get(); // Читаем R
            double value = color; // Цвет в высоту
            
            pixelMatrix[y][x] = value; // Обновляем матрицу значений
            p->field[y][x] = value; // Обновляем поле в Pole
            //Значения в gaussi некорректные
        }
        bmpFile.ignore((4 - (width * 3) % 4) % 4); // Пропускаем паддинг
    }
    bmpFile.close();
  }
};

class PathFinder {
private:
    double vehicleRadius;
public:
    PathFinder(double radius) : vehicleRadius(radius) {}
    // Эвристическая функция (евклидово расстояние)
    double heuristic(const PointD& a, const PointD& b) {
        return std::hypot(a.x - b.x, a.y - b.y);
    }

    // Найти треугольник, содержащий точку
    const Triangle* findContainingTriangle(const PointD& p, const std::vector<Triangle>& triangles) {
        for (const auto& tri : triangles) {
            if (isPointInTriangle(p, tri)) return &tri;
        }
        return nullptr;
    }

    // Проверка, находится ли точка внутри треугольника
   bool isPointInTriangle(const PointD& p, const Triangle& tri) {
    // Заменить барицентрический метод на метод площадей
    auto sign = [](PointD p1, PointD p2, PointD p3) {
        return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
    };
    
    double d1 = sign(p, tri.a, tri.b);
    double d2 = sign(p, tri.b, tri.c);
    double d3 = sign(p, tri.c, tri.a);
    
    bool has_neg = (d1 < -1e-6) || (d2 < -1e-6) || (d3 < -1e-6);
    bool has_pos = (d1 > 1e-6) || (d2 > 1e-6) || (d3 > 1e-6);
    
    return !(has_neg && has_pos);
}
    
   bool isNavigable(const Triangle& tri, double radius) {
    const double minSafeElevation = 0.5; // Минимально допустимая высота
    double safeElevation = std::max(tri.minElevation, minSafeElevation);
    bool navigable = safeElevation > radius;

    if (!navigable) {
        std::cout << "Треугольник непроходим: мин. высота = " << tri.minElevation 
                  << " (с учётом защиты: " << safeElevation << ")"
                  << ", радиус = " << radius 
                  << ", вершины: (" 
                  << tri.a.x << ", " << tri.a.y << "), ("
                  << tri.b.x << ", " << tri.b.y << "), ("
                  << tri.c.x << ", " << tri.c.y << ")\n";
    }
    return navigable;
}
    // Найти соседние треугольники для заданного треугольника
    std::vector<const Triangle*> findNeighbors(const Triangle& tri, const std::vector<Triangle>& allTriangles) const {
        std::vector<const Triangle*> neighbors;
        for (const auto& other : allTriangles) {
            if (&tri == &other) continue; // Пропустить текущий треугольник
            if (shareEdge(tri, other)) {
                neighbors.push_back(&other);
            }
        }
        return neighbors;
    }

public:
    // Найти путь от start до goal через триангуляцию
    std::vector<PointD> findPathAStar(const PointD& start, const PointD& goal, const std::vector<Triangle>& triangles) {
        std::vector<PointD> path;

        // Найти стартовый и целевой треугольники
        const Triangle* startTri = findContainingTriangle(start, triangles);
        const Triangle* goalTri = findContainingTriangle(goal, triangles);
        if (!startTri) {
        std::cerr << "Стартовая точка (" << start.x << ", " << start.y << ") вне триангуляции!\n";
        return {};
    }
    if (!goalTri) {
        std::cerr << "Конечная точка (" << goal.x << ", " << goal.y << ") вне триангуляции!\n";
        return {};
    }
    if (startTri == goalTri) {
    return {start, goal}; // Прямой путь между точками
}

        // Приоритетная очередь для обработки узлов
        std::priority_queue<AStarNode> openSet;
        std::unordered_map<const Triangle*, double> costSoFar;
         // Исправление: проверка на пустую триангуляцию
    if (triangles.empty()) {
        std::cerr << "Триангуляция не построена!" << std::endl;
        return path;
    }

        // Инициализация стартового узла
        PointD startCenter = startTri->calculateCircumcenter();
        openSet.emplace(startCenter, startTri);
        costSoFar[startTri] = 0;

        while (!openSet.empty()) {
            AStarNode current = openSet.top();
            openSet.pop();

            // Если достигли цели
            if (current.tri == goalTri) {
                // Восстановление пути
            while (current.parent) {
                path.push_back(current.position);
                current = *current.parent;
            }
            std::reverse(path.begin(), path.end());
            
            // Добавляем реальные точки A и B
            if (!path.empty()) {
                path.insert(path.begin(), start); // Стартовая точка
                path.push_back(goal);             // Целевая точка
            }
            
            return path;  // Возвращаем корректный путь
        }

        // Обработка соседей
        for (const auto& neighbor : getNeighbors(*current.tri, triangles)) {
            if (!isNavigable(*neighbor, vehicleRadius)) continue;
            
            PointD neighborCenter = neighbor->calculateCircumcenter();
            double newCost = current.g + heuristic(current.position, neighborCenter);

            if (!costSoFar.count(neighbor) || newCost < costSoFar[neighbor]) {
                costSoFar[neighbor] = newCost;
                double priority = newCost + heuristic(neighborCenter, goal);
                openSet.emplace(neighborCenter, neighbor, new AStarNode(current));
            }
        }
    }

    return {};  // Путь не найден
}

    // Получить соседей треугольника
    std::vector<const Triangle*> getNeighbors(const Triangle& tri, const std::vector<Triangle>& triangles) {
        std::vector<const Triangle*> neighbors;
        for (const auto& other : triangles) {
            if (&tri == &other) continue;
            if (shareEdge(tri, other)) neighbors.push_back(&other);
        }
        return neighbors;
    }

    // Проверить, имеют ли треугольники общее ребро
    bool shareEdge(const Triangle& a, const Triangle& b) const {
    for (const auto& edgeA : { Edge(a.a, a.b), Edge(a.b, a.c), Edge(a.c, a.a) }) {
        for (const auto& edgeB : { Edge(b.a, b.b), Edge(b.b, b.c), Edge(b.c, b.a) }) {
            if (edgeA == edgeB) return true;
        }
    }
    return false;
}
    bool otherHasEdge(const Triangle& other, const Edge& edge) const {
    return (Edge(other.a, other.b) == edge) || 
           (Edge(other.b, other.c) == edge) || 
           (Edge(other.c, other.a) == edge);
}
};

class VoronoiDiagram {
public:
 const PathFinder& pathFinder; // Добавляем ссылку на PathFinder
    std::vector<VoronoiEdge> edges;
    
    VoronoiDiagram(const PathFinder& pf) : pathFinder(pf) {} // Конструктор
    
    void buildFromDelaunay(const std::vector<Triangle>& triangles, const PathFinder& pathFinder, const std::unique_ptr<Pole>& p) {
    const int height = p->field.size();
    const int width = p->field[0].size();
    
    if (!p) {
        std::cerr << "Pole is not initialized!" << std::endl;
        return;
    }
    
    for (const auto& tri : triangles) {
        PointD cc1 = tri.calculateCircumcenter();
        
        // Проверка на выход центра за границы поля
        if (cc1.x < 0 || cc1.y < 0 || cc1.x >= width || cc1.y >= height) {
            // Обработка выхода за границы (например, пропускаем этот треугольник)
            continue; // Пропускаем дальнейшую обработку для этого треугольника
        }
        
        auto neighbors = pathFinder.findNeighbors(tri, triangles);
        
        // Если соседей меньше 3, это граничный треугольник
        if (neighbors.size() < 3) {
            // Находим граничные рёбра и обрезаем их
            std::vector<Edge> boundaryEdges = getBoundaryEdges(tri, triangles);
            for (const auto& edge : boundaryEdges) {
                PointD intersection = calculateBoundaryIntersection(cc1, edge, width, height);
                edges.emplace_back(cc1, intersection);
            }
        } 
            // Обычные рёбра между центрами
            for (const auto& neighbor : neighbors) {
                PointD cc2 = neighbor->calculateCircumcenter();
                
                // Проверка на выход центра соседа за границы поля
                if (cc2.x < 0 || cc2.y < 0 || cc2.x >= width || cc2.y >= height) {
                    continue; // Пропускаем дальнейшую обработку для этого соседа
                }
                
                edges.emplace_back(cc1, cc2);
            }
        
    }
}

private:
    std::vector<Edge> getBoundaryEdges(const Triangle& tri, const std::vector<Triangle>& allTriangles) {
        std::vector<Edge> boundaryEdges;
        for (const auto& edge : { Edge(tri.a, tri.b), Edge(tri.b, tri.c), Edge(tri.c, tri.a) }) {
            bool isBoundary = true;
            for (const auto& other : allTriangles) {
                if (&tri == &other) continue;
                if (pathFinder.shareEdge(tri, other) && pathFinder.otherHasEdge(other, edge)) {
                    isBoundary = false;
                    break;
                }
            }
            if (isBoundary) boundaryEdges.push_back(edge);
        }
        return boundaryEdges;
    }

    PointD calculateBoundaryIntersection(const PointD& circumCenter, const Edge& boundaryEdge, int width, int height) {
        // Вычисляем направление перпендикуляра от граничного ребра
        PointD edgeDir = { boundaryEdge.b.x - boundaryEdge.a.x, boundaryEdge.b.y - boundaryEdge.a.y };
        PointD normal = { -edgeDir.y, edgeDir.x }; // Перпендикуляр
        double length = std::hypot(normal.x, normal.y);
        if (std::fabs(length) < 1e-9) return circumCenter; // Защита от нулевой длины
        normal.x /= length;
        normal.y /= length;

        // Находим пересечение с границами поля
        double t = std::numeric_limits<double>::max();
        if (std::fabs(normal.x) > 1e-6) {
            t = std::min(t, (width - circumCenter.x) / normal.x);
            t = std::min(t, -circumCenter.x / normal.x);
        }
        if (std::fabs(normal.y) > 1e-6) {
            t = std::min(t, (height - circumCenter.y) / normal.y);
            t = std::min(t, -circumCenter.y / normal.y);
        }
        
        return {
            circumCenter.x + normal.x * t,
            circumCenter.y + normal.y * t
        };
    }
};
   
class GnuplotInterface {
private:
 Config& config; // Добавляем ссылку на Config
 
  // Преобразование Y-координаты для Gnuplot
    double transformY(double y, int height) const {
        return height - y - 1;
    }
    
   public:
    GnuplotInterface(Config& cfg) : config(cfg) {} // Конструктор
    void plotBinaryWithComponents(const std::vector<std::vector<double>> &CopyPole, const std::vector<Component>& components) {

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            std::cerr << "Failed to open gnuplot pipe" << std::endl;
            return;
        }

        const int height = CopyPole.size();
        const int width = CopyPole[0].size();

        // Настройки графика
        fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200\n");
        fprintf(gnuplotPipe, "set output '/home/log/Gauss/results/visualizations/binary_with_components.png'\n");
        fprintf(gnuplotPipe, "set title 'Binary Image with Components Metadata'\n");
        fprintf(gnuplotPipe, "set size ratio -1\n");
        fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
        fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
        fprintf(gnuplotPipe, "unset key\n");

        // Многослойный график
        fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
        fprintf(gnuplotPipe, "'-' with lines lw 2 lc 'red', \\\n");
        fprintf(gnuplotPipe, "'-' with points pt 7 ps 2 lc 'blue', \\\n");
        fprintf(gnuplotPipe, "'-' with vectors head filled lc 'green'\n");

        // 1. Данные бинарного изображения (с инверсией Y)
for (int y = height - 1; y >= 0; --y) {
    for (int x = 0; x < width; ++x) {
        fprintf(gnuplotPipe, "%f ", std::fabs(CopyPole[y][x]));
    }
            fprintf(gnuplotPipe, "\n");
        }
        fprintf(gnuplotPipe, "e\n");

         // 2. Границы компонент (красные прямоугольники)
        for (const auto& comp : components) {
            double y_min = transformY(comp.min_y, height);
            double y_max = transformY(comp.max_y, height);
            
            fprintf(gnuplotPipe, "%f %f\n", comp.min_x - 0.5, y_min - 0.5);
            fprintf(gnuplotPipe, "%f %f\n", comp.max_x + 0.5, y_min - 0.5);
            fprintf(gnuplotPipe, "%f %f\n", comp.max_x + 0.5, y_max + 0.5);
            fprintf(gnuplotPipe, "%f %f\n", comp.min_x - 0.5, y_max + 0.5);
            fprintf(gnuplotPipe, "%f %f\n\n", comp.min_x - 0.5, y_min - 0.5);
        }
        fprintf(gnuplotPipe, "e\n");

        // 3. Центры компонент (синие точки)
        for (const auto& comp : components) {
            fprintf(gnuplotPipe, "%f %f\n", 
                    comp.center_x, 
                    transformY(comp.center_y, height));
        }
        fprintf(gnuplotPipe, "e\n");

         // 4. Собственные векторы (зеленые стрелки)
        for (const auto& comp : components) {
        const double scale = std::min(comp.max_x - comp.min_x, comp.max_y - comp.min_y) * 0.2;
            double cy = transformY(comp.center_y, height);
            
            // Первый собственный вектор
            fprintf(gnuplotPipe, "%f %f %f %f\n", 
                    comp.center_x,
                    cy,
                    comp.eigenvec1_x * scale,
                    -comp.eigenvec1_y * scale); // Инверсия Y

            // Второй собственный вектор
            fprintf(gnuplotPipe, "%f %f %f %f\n", 
                    comp.center_x,
                    cy,
                    comp.eigenvec2_x * scale,
                    -comp.eigenvec2_y * scale);// Инверсия Y
        }
        fprintf(gnuplotPipe, "e\n");

        pclose(gnuplotPipe);
    }
      void gnuplot(std::unique_ptr<Pole>& p, const std::string& filename) {
    if (p == nullptr) {
        std::cout << "Pole not initialized!" << std::endl;
        return;
    }

    int rows = p->field.size(); 
    int cols = p->field[0].size(); 
    // Открываем конвейер для gnuplot 
    FILE* gnuplotPipe = popen("gnuplot -p", "w"); 
    if (!gnuplotPipe) { 
        std::cerr << "Could not open pipe to gnuplot." << std::endl; 
        return; 
    }
    
    // Устанавливаем диапазон z от -255 до 255
    fprintf(gnuplotPipe, "set zrange [-255:255]\n"); 
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", cols - 1); 
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", rows - 1); 
    fprintf(gnuplotPipe, "set terminal png\n"); 
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());

    // Включаем pm3d для плавного отображения
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    
    // Используем буфер для хранения данных
    std::ostringstream dataStream;
    
    // Форматируем данные в виде матрицы
    for (int y = rows - 1; y >= 0; --y) {
        for (int x = 0; x < cols; ++x) { 
            dataStream << x << " " << y << " " << p->field[y][x] << "\n";
        }
        dataStream << "\n"; // Пустая строка между строками данных
    }

    // Отправляем все данные за один раз
    fprintf(gnuplotPipe, "splot '-' with pm3d\n");
    
    fprintf(gnuplotPipe, "%s", dataStream.str().c_str());

    fprintf(gnuplotPipe, "\n"); // Одна пустая строка для завершения ввода данных
    pclose(gnuplotPipe);
}

void plotVoronoi(const std::unique_ptr<Pole>& p, const std::vector<VoronoiEdge>& edges, const std::vector<PointD>& sites) {
    if (!p || edges.empty()) return;

    const int height = p->field.size();
    const int width = p->field[0].size();

    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        std::cerr << "Failed to open gnuplot pipe" << std::endl;
        return;
    }

    // Настройки графика
    fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200\n");
    fprintf(gnuplotPipe, "set output '/home/log/Gauss/results/visualizations/voronoi_diagram.png'\n");
    fprintf(gnuplotPipe, "set title 'Voronoi Diagram'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height - 1);
    fprintf(gnuplotPipe, "unset key\n");
    
    // Отрисовка поля, рёбер и центров
        fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
        fprintf(gnuplotPipe, "'-' with lines lw 2 lc 'red', \\\n");
        fprintf(gnuplotPipe, "'-' with points pt 7 ps 2 lc 'blue', \n");
    // 1. Данные поля
    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", std::fabs(p->field[y][x]));
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 2. Данные рёбер
    for (const auto& edge : edges) {
        fprintf(gnuplotPipe, "%f %f\n%f %f\n\n", 
                edge.start.x, transformY(edge.start.y, height),
                edge.end.x, transformY(edge.end.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    // 3. Центры кластеров
    for (const auto& site : sites) {
        fprintf(gnuplotPipe, "%f %f\n", 
                site.x, 
                transformY(site.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
}
    
    void plotDelaunay(const std::vector<Triangle>& triangles, std::unique_ptr<Pole>& p) {
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        std::cerr << "Failed to open gnuplot pipe" << std::endl;
        return;
    }
if (!p) return;
     const int height = p->field.size();
    const int width = p->field[0].size();

    // Настройки графика
    fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200\n");
    fprintf(gnuplotPipe, "set output '/home/log/Gauss/results/visualizations/delaunay_triangulation.png'\n");
    fprintf(gnuplotPipe, "set title 'Delaunay Triangulation'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
    fprintf(gnuplotPipe, "unset key\n");

    // Многослойный график: фон + треугольники
    fprintf(gnuplotPipe, "plot '-' matrix with image, '-' with lines lc 'red'\n");

    // Данные фона (с инверсией Y)
    for (int y = height-1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", std::fabs(p->field[y][x]));
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // Данные треугольников (с инверсией Y)
    for (const auto& tri : triangles) {
        fprintf(gnuplotPipe, "%f %f\n", tri.a.x, transformY(tri.a.y, height));
        fprintf(gnuplotPipe, "%f %f\n", tri.b.x, transformY(tri.b.y, height));
        fprintf(gnuplotPipe, "%f %f\n", tri.c.x, transformY(tri.c.y, height));
        fprintf(gnuplotPipe, "%f %f\n\n", tri.a.x, transformY(tri.a.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
}

void plotPath(const std::vector<PointD>& path, const std::unique_ptr<Pole>& p, const std::string& outputFile = "/home/log/Gauss/results/visualizations/path_plot.png") {
    if (!p || path.empty()) return;

    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        std::cerr << "Failed to open gnuplot pipe." << std::endl;
        return;
    }

    const int height = p->field.size();
    const int width = p->field[0].size();

    // Настройки графика
    fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200\n");
    fprintf(gnuplotPipe, "set output '%s'\n", outputFile.c_str());
    fprintf(gnuplotPipe, "set title 'Path Visualization'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
    fprintf(gnuplotPipe, "unset key\n");

    // Многослойный график: фон + путь + точки
    fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
    fprintf(gnuplotPipe, "'-' with lines lw 2 lc 'red', \\\n");
    fprintf(gnuplotPipe, "'-' with points pt 7 ps 2 lc 'blue', \\\n");
    fprintf(gnuplotPipe, "'-' with points pt 9 ps 2 lc 'green'\n");

    // 1. Данные фона (инверсия Y)
    for (int y = height-1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", std::fabs(p->field[y][x]));
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 2. Путь (инверсия Y)
    for (const auto& point : path) {
        fprintf(gnuplotPipe, "%f %f\n", 
                point.x, 
                transformY(point.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    // 3. Точка A (синяя)
    fprintf(gnuplotPipe, "%f %f\n", 
            config.pointA_x, 
            transformY(config.pointA_y, height));
    fprintf(gnuplotPipe, "e\n");

    // 4. Точка B (зеленая)
    fprintf(gnuplotPipe, "%f %f\n", 
            config.pointB_x, 
            transformY(config.pointB_y, height));
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
}

};
   
class ComponentCalculator {
public:
int length, width;
    std::vector<std::vector<double>> CopyField;
    std::vector<Component> components;

    std::vector<Triangle> bowyerWatson(const std::vector<PointD>& points, const Pole& elevationData) {
        std::vector<Triangle> triangles;
        if (points.empty()) return triangles;
         
         // Исправление: проверка на достаточное количество точек
    if (points.size() < 3) {
        std::cerr << "Недостаточно точек для триангуляции!" << std::endl;
        return triangles;
    }
  
// Проверка валидности всех точек
    for (const auto& p : points) {
        if (std::isnan(p.x) || std::isnan(p.y)) {
            std::cerr << "Skipping invalid point (NaN)" << std::endl;
            return triangles;
        }
    }
        // Создаем супер-треугольник
        double minX = points[0].x, maxX = points[0].x;
        double minY = points[0].y, maxY = points[0].y;
        for (const auto& p : points) {
            minX = std::min(minX, p.x); maxX = std::max(maxX, p.x);
            minY = std::min(minY, p.y); maxY = std::max(maxY, p.y);
        }
        
        double dx = (maxX - minX) * 10, dy = (maxY - minY) * 10;
        PointD p1(minX - dx, minY - dy);
        PointD p2(maxX + dx, minY - dy);
        PointD p3((minX + maxX)/2, maxY + dy);
        triangles.emplace_back(p1, p2, p3);

        // Добавляем точки
        for (const auto& p : points) {
            std::vector<Triangle> badTriangles;
            for (const auto& tri : triangles) {
                if (isPointInCircumcircle(p, tri)) {
                    badTriangles.push_back(tri);
                }
            }

            std::vector<Edge> polygon;
            for (const auto& tri : badTriangles) {
                std::vector<Edge> edges = { {tri.a, tri.b}, {tri.b, tri.c}, {tri.c, tri.a} };
                for (const auto& edge : edges) {
                    bool isShared = false;
                    for (const auto& other : badTriangles) {
                        if (&tri == &other) continue;
                        if (otherHasEdge(other, edge)) isShared = true;
                    }
                    if (!isShared) polygon.push_back(edge);
                }
            }

            // Удаляем плохие треугольники
            triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
                [&](const Triangle& t) { return std::find(badTriangles.begin(), badTriangles.end(), t) != badTriangles.end(); }),
                triangles.end());

            // Добавляем новые
            for (const auto& edge : polygon) {
                triangles.emplace_back(edge.a, edge.b, p);
            }
        }
// Удаляем треугольники супер-треугольника
        triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
        [&](const Triangle& t) {
            return contains(t, p1) || contains(t, p2) || contains(t, p3);
        }), triangles.end());
        // Обновляем минимальные высоты треугольников с проверкой на принадлежность полю
        for (auto& tri : triangles) {
        auto safeGetElevation = [&](const PointD& point) -> double {
            int x = static_cast<int>(point.x);
            int y = static_cast<int>(point.y);
            if (x < 0 || x >= elevationData.field[0].size() || 
                y < 0 || y >= elevationData.field.size()) {
                std::cerr << "Точка (" << x << ", " << y << ") вне поля!" << std::endl;
                return 0.5; // Или любое другое значение по умолчанию
            }
            return elevationData.field[y][x];
        };

        double hA = safeGetElevation(tri.a);
        double hB = safeGetElevation(tri.b);
        double hC = safeGetElevation(tri.c);

        tri.minElevation = std::min({hA, hB, hC});
    }

    return triangles;
}

    bool isPointInCircumcircle(const PointD& p, const Triangle& tri) {
        double ax = tri.a.x - p.x, ay = tri.a.y - p.y;
        double bx = tri.b.x - p.x, by = tri.b.y - p.y;
        double cx = tri.c.x - p.x, cy = tri.c.y - p.y;
        
        double det = ax * (by * (cx*cx + cy*cy) - cy * (bx*bx + by*by)) - ay * (bx * (cx*cx + cy*cy) - cx * (bx*bx + by*by)) + (ax*ax + ay*ay) * (bx*cy - by*cx);
        return det > 0;
    }

    bool otherHasEdge(const Triangle& tri, const Edge& edge) {
        return Edge(tri.a, tri.b) == edge 
            || Edge(tri.b, tri.c) == edge 
            || Edge(tri.c, tri.a) == edge;
    }

    bool contains(const Triangle& tri, const PointD& p) {
        return tri.a == p || tri.b == p || tri.c == p;
    }

    int count = 0;//Для шума
       
    int incrementAndCollect(std::vector<std::vector<double>>& componenta, std::vector<std::vector<double>> &CopyPole, int x, int y, int i) {

        if (x < 1 || y < 1 || x > (int) componenta[0].size() - 2 || y > (int) componenta.size() - 2 || CopyPole[y][x] < 250) return -1;

        if (CopyPole[y][x] >= 255 && CopyPole[y][x] <= 255) {
            CopyPole[y][x] = 0; // Пометить как посещенное
            count = count < i + 1 ? i + 1 : count;
            componenta[y][x] = 255; // Увеличить значение в Componenta
            incrementAndCollect(componenta, CopyPole, x + 1, y, i + 1);
            incrementAndCollect(componenta, CopyPole, x - 1, y, i + 1);
            incrementAndCollect(componenta, CopyPole, x, y + 1, i + 1);
            incrementAndCollect(componenta, CopyPole, x, y - 1, i + 1);
        }
        
    return count;
}
       void bin(std::vector<std::vector<double>> &CopyPole, int slise, std::unique_ptr<Pole>& p) {
       if (p == nullptr) {
       std::cerr << "Pole not initialized!" << std::endl;
       return;
   }
       CopyPole = p->field;//Копия
        
            for (int x = 0; x < (int) p->field[0].size(); ++x) {
                for (int y = 0; y < (int) p->field.size(); ++y) {
                   CopyPole[y][x] = std::fabs(p->field[y][x]) > slise ? 255 : 0;
                } 
            }
    }
    
    void wave(int noisy, std::vector<Component>& componenti, std::vector<std::vector<double>>& CopyPole, std::unique_ptr<Pole>& p) {
    if (p == nullptr) {
        std::cerr << "Pole not initialized!" << std::endl;
        return;
    }

    int rows = p->field.size();
    int cols = (rows > 0) ? p->field[0].size() : 0;

    // Проходим по всем пикселям поля
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            // Если пиксель принадлежит компоненте (значение 255)
            if (CopyPole[y][x] <= 255 && CopyPole[y][x] >= 255) {
                // Создаем временную матрицу для хранения данных компоненты
                std::vector<std::vector<double>> componentData(rows, std::vector<double>(cols, 0));
                
                // Заполняем компоненту и получаем количество пикселей
                int pixelCount = incrementAndCollect(componentData, CopyPole, x, y, 0);

                // Если компонента достаточно большая, добавляем её
                if (pixelCount >= noisy) {
                    // Используем конструктор, который вызывает calculate_metadata()
                    Component component(componentData);
                    componenti.push_back(component);
                }
                else {// Удаляем шум из основного поля
                    for (int i = 0; i < rows; ++i) {
                        for (int j = 0; j < cols; ++j) {
                            if (componentData[i][j] == 255) {
                                p->field[i][j] = 0; // Устанавливаем значение в 0 для удаления шума
                            }
                        }
                    }
                }
            }
        }
    }
  }
};

class KMeans {
private:

public:
    std::vector<int> labels;
    std::vector<std::vector<double>> centers; // Центры кластеров

struct ClusterResult {
        std::vector<int> labels;
        std::vector<std::array<int, 3>> colors;
        std::vector<std::vector<double>> centers;
    };
    
    void initializeCenters(const std::vector<std::vector<double>>& data, int k) {
    if (data.empty() || k <= 0) return;
    
    centers.clear();
    std::vector<size_t> indices(data.size());
    std::iota(indices.begin(), indices.end(), 0);
    std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});
    
    for (int i = 0; i < k && i < static_cast<int>(data.size()); ++i) {
        centers.push_back(data[indices[i]]);
    }
}

    double distance(const std::vector<double>& a, const std::vector<double>& b) const {
        return std::hypot(a[0]-b[0], a[1]-b[1]);
    }

    ClusterResult cluster(const std::vector<std::vector<double>>& data, int k) {
        ClusterResult result;
        initializeCenters(data, k);
        bool changed;
        
        do {
            changed = false;
            labels.assign(data.size(), -1);

            // Assign labels
            for (size_t i = 0; i < data.size(); ++i) {
                double minDist = std::numeric_limits<double>::max();
                for (size_t j = 0; j < centers.size(); ++j) {
                    double dist = distance(data[i], centers[j]);
                    if (dist < minDist) {
                        minDist = dist;
                        labels[i] = j;
                    }
                }
            }

            // Update centers
            std::vector<std::vector<double>> newCenters(centers.size(), {0.0, 0.0});
            std::vector<int> counts(centers.size(), 0);

            for (size_t i = 0; i < data.size(); ++i) {
                int label = labels[i];
                if (label >= 0) {
                    newCenters[label][0] += data[i][0];
                    newCenters[label][1] += data[i][1];
                    counts[label]++;
                }
            }

            // Check for changes
            for (size_t j = 0; j < centers.size(); ++j) {
                if (counts[j] > 0) {
                    newCenters[j][0] /= counts[j];
                    newCenters[j][1] /= counts[j];
                    if (distance(newCenters[j], centers[j]) > 1e-6) {
                        changed = true;
                    }
                }
            }
            centers = newCenters;

        } while (changed);

        result.labels = labels;
        result.centers = centers;
        result.colors = ColorGenerator::generateColors(k);
        return result;
    }

    
    ClusterResult kmeansWithKernels(const std::vector<std::vector<double>>& data, int k, int kernelSize) {
    ClusterResult result;
    if (data.empty() || k <= 0 || kernelSize <= 0) {
        std::cerr << "Invalid input: data empty or k/kernelSize <= 0" << std::endl;
        return result;
    }

    // Step 1: Sample kernel centers
    std::vector<std::vector<double>> kernelCenters;
    std::vector<size_t> indices(data.size());
    std::iota(indices.begin(), indices.end(), 0);
    
    size_t samplesToTake = std::min(k * kernelSize, static_cast<int>(data.size()));
    if (samplesToTake == 0) {
        std::cerr << "Not enough data points for sampling" << std::endl;
        return result;
    }

    std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});
    for (size_t i = 0; i < samplesToTake; ++i) {
        kernelCenters.push_back(data[indices[i]]);
    }

    // Step 2: Initial clustering
    initializeCenters(kernelCenters, k * kernelSize);
    auto kernelResult = cluster(kernelCenters, k * kernelSize);

    // Step 3: Aggregate centers
    std::vector<std::vector<double>> aggregatedCenters(k, {0.0, 0.0});
    std::vector<int> counts(k, 0);

    for (size_t i = 0; i < kernelResult.labels.size(); ++i) {
        int label = kernelResult.labels[i];
        if (label >= 0 && label < k) {
            aggregatedCenters[label][0] += kernelCenters[i][0];
            aggregatedCenters[label][1] += kernelCenters[i][1];
            counts[label]++;
        }
    }

    // Step 4: Normalize or replace empty clusters
    for (int i = 0; i < k; ++i) {
        if (counts[i] > 0) {
            aggregatedCenters[i][0] /= counts[i];
            aggregatedCenters[i][1] /= counts[i];
        } else {
            // Fallback: random point from data
            aggregatedCenters[i] = data[std::rand() % data.size()];
        }
    }

    // Final clustering
    centers = aggregatedCenters;
    result = cluster(data, k);
    if (result.centers.size() != k) {
        std::cerr << "Cluster count mismatch!" << std::endl;
    }
    result.colors = ColorGenerator::generateColors(k);

    return result;
}
};

// Класс для работы с сеткой рельефа
class TerrainGrid {
public:
    std::vector<std::vector<GridCell>> cells;

    // Инициализация сетки на основе размеров поля
    void initialize(int width, int height) {
         // 1. Изменение размера двумерного вектора cells
    //    - height строк
    //    - width столбцов в каждой строке
    cells.resize(height, std::vector<GridCell>(width));

    // 2. Двойной цикл для итерации по всем ячейкам сетки
    for (int y = 0; y < height; ++y) {         // Проход по строкам (ось Y)
        for (int x = 0; x < width; ++x) {      // Проход по столбцам (ось X)
            
            // 3. Расчет координат центра ячейки:
            cells[y][x].position = PointD(
                // Ось X: 
                x + 0.5,                       // Центр ячейки по X
                
                // Ось Y (инверсия):
                (height - y - 1) + 0.5         // Центр ячейки по Y с инверсией
            );
            
            // Пример для сетки 3x3:
            // y=0 -> Y = (3-0-1)+0.5 = 2.5 (верхний ряд)
            // y=1 -> Y = (3-1-1)+0.5 = 1.5 (средний ряд)
            // y=2 -> Y = (3-2-1)+0.5 = 0.5 (нижний ряд)
        }
    }
    }

    // Расчет углов наклона на основе данных высот
    void calculateSlopes(const Pole& elevationData) {
        if (elevationData.field.empty() || 
            cells.size() != elevationData.field.size() || 
            cells[0].size() != elevationData.field[0].size()) {
            throw std::invalid_argument("TerrainGrid and Pole size mismatch");
        }

        const int height = cells.size();
        const int width = cells[0].size();

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                // Вычисляем производные по X и Y с помощью центральных разностей
                double dzdx, dzdy;

                // Обработка граничных условий
                if (x == 0 || x == width - 1) {
                    dzdx = 0; // Игнорируем края по X
                } else {
                    dzdx = (elevationData.field[y][x+1] - elevationData.field[y][x-1]) / 2.0;
                }

                if (y == 0 || y == height - 1) {
                    dzdy = 0; // Игнорируем края по Y
                } else {
                    dzdy = (elevationData.field[y+1][x] - elevationData.field[y-1][x]) / 2.0;
                }

                // Вычисляем угол наклона (в градусах)
                cells[y][x].slopeAngle = std::atan(std::hypot(dzdx, dzdy)) * 180.0 / M_PI;
            }
        }
    }

    // Сброс флагов посещения
    void resetVisited() {
        for (auto& row : cells) {
            for (auto& cell : row) {
                cell.visited = false;
            }
        }
    }
};

class Control {
private:
    bool b = true; 
public:
    //Для логирования
    Config& config;
    Logger& loggercontrol;
    //Для функций
    Copier copier;
    std::vector<std::array<int, 3>> colors; // Объявление colors
    std::vector<std::vector<double>> CopyPole;
    std::vector<Gaus> gaussi;
    std::vector<Component> componenti;
    GaussBuilder gaussBuilder;
    BmpHandler bmpHandler;
    GnuplotInterface gnuplotInterface;
    ComponentCalculator componentCalculator;
    std::unique_ptr<Pole> p= nullptr; // Pointer to Pole
    std::unique_ptr<KMeans> kMeans = nullptr; // Используем умный указатель
    std::vector<std::vector<double>> kMeansData;
    std::vector<Triangle> lastTriangulation;
    PathFinder pathFinder;
    VoronoiDiagram voronoi{pathFinder}; // Инициализируем через конструктор;
    TerrainGrid terrainGrid;
     
    Control(Config& cfg, Logger& log) : config(cfg), loggercontrol(log), gnuplotInterface(cfg), pathFinder(cfg.vehicleRadius) {
     kMeans = std::make_unique<KMeans>(); // Инициализация
        if (config.loggingControlEnabled) {
            
                loggercontrol.logMessage("Logging Control is enabled.", b);
                std::cout << "Logging Control is enabled." << std::endl;
            } else {
            loggercontrol.logMessage("Logging Control is disabled.", b);
            std::cout << "Logging Control is disabled." << std::endl;
            b = false;
        }      
    }
    // Получить центры компонент
    std::vector<PointD> getClusterCenters() const {
    std::vector<PointD> centers;
    for (const auto& component : componenti) {
    if (std::isnan(component.center_x) || std::isnan(component.center_y)) {
            std::cerr << "Skipping invalid cluster center (NaN)" << std::endl;
            continue;
        }
        if (component.center_x >= 0 && component.center_x < config.fieldWidth &&
            component.center_y >= 0 && component.center_y < config.fieldHeight) {
            centers.emplace_back(component.center_x, component.center_y);
        } else {
            std::cerr << "Invalid cluster center: (" 
                      << component.center_x << ", " 
                      << component.center_y << ")" << std::endl;
        }
    }
    return centers;
}

void applyClusterResults(const KMeans::ClusterResult& result, 
                           std::vector<std::vector<double>>& pixelMatrix) {
        for (size_t i = 0; i < result.labels.size(); ++i) {
            int label = result.labels[i];
            if (label >= 0 && label < result.colors.size()) {
                int x = static_cast<int>(kMeansData[i][0]);
                int y = static_cast<int>(kMeansData[i][1]);
                const auto& color = result.colors[label];
                pixelMatrix[y][x] = (color[0] + color[1] + color[2]) / 3.0;
            }
        }
    }

    void prepareKMeansData(const std::vector<std::vector<double>>& copyPole) {
    kMeansData.clear();
    for (size_t y = 0; y < copyPole.size(); ++y) {
        for (size_t x = 0; x < copyPole[0].size(); ++x) {
            if (copyPole[y][x] > 0) { 
                kMeansData.push_back({static_cast<double>(x), static_cast<double>(y)});
            }
        }
    }
}

    
void Dispetcher(DispatcherParams& params) {
   
    if (params.s == "init") {
    gaussBuilder.init(params.A, params.B, p);
    loggercontrol.logMessage("init used: " + std::to_string(params.A) + " x " + std::to_string(params.B), b);
    }

    if (params.s == "g") {
    gaussBuilder.addgauss(params.h, params.x, params.y, params.sx, params.sy, gaussi);
    loggercontrol.logMessage("g used: x=" + std::to_string(params.x) +
                    ", y=" + std::to_string(params.y) + 
                    ", sx=" + std::to_string(params.sx) + 
                    ", sy=" + std::to_string(params.y) + 
                    ", h=" + std::to_string(params.h), b);
    }
    
    if (params.s == "generate") {
    gaussBuilder.generate(p, gaussi);
    loggercontrol.logMessage("generate used", b);
    }
    
    if (params.s == "gnuplot") {
    gnuplotInterface.gnuplot(p, params.filename);
    loggercontrol.logMessage("gnuplot used", b);
    }

    if (params.s == "bmp_write") {
        if (params.is_binary_image == 1) {
            bmpHandler.bmp_write(p->field, params.filename);
            loggercontrol.logMessage("bmp_write used", b);
        } else {
            bmpHandler.bmp_write(CopyPole, params.filename);
            loggercontrol.logMessage("bmp_write used", b);
        }
    }

    if (params.s == "bmp_read") {
    bmpHandler.bmp_read(gaussBuilder, params.filename, CopyPole, p);
    loggercontrol.logMessage("bmp_read used, filename: " + params.filename, b);
    }

    if (params.s == "bin") {
    componentCalculator.bin(CopyPole, params.slice, p);
    loggercontrol.logMessage("bin used, slice=" + std::to_string(params.slice), b);
    }
    
    if (params.s == "wave") {
    componentCalculator.wave(params.noisy, componenti, CopyPole, p);
    loggercontrol.logMessage("wave used", b);
    loggercontrol.logMessage("Component amount = " + std::to_string(componenti.size()), b);
    copier.removeNoise(CopyPole, componenti);//копия без шума
    gnuplotInterface.plotBinaryWithComponents(CopyPole, componenti);//визуализация метаданных   
    }
     
    if (params.s == "k_means") {
            if (params.k <= 0 || p == nullptr) return;
            
            prepareKMeansData(CopyPole);
            if (kMeansData.empty()) return;

            auto result = kMeans->cluster(kMeansData, params.k);
            applyClusterResults(result, CopyPole);
            //copier.removeNoise(CopyPole, componenti);//Если не хочешь засорять COPY но не работает визуализация
        }
    if (params.s == "k_means_kern") {
            if (!kMeans || params.k <= 0 || params.kk <= 0 || p == nullptr) {
        loggercontrol.logMessage("Invalid parameters or KMeans not initialized", b);
        return;
    }
            copier.removeNoise(CopyPole, componenti);//копия без шума
            prepareKMeansData(CopyPole);
            if (kMeansData.empty()) {
                loggercontrol.logMessage("No data available for clustering", b);
                return;
            }

            auto result = kMeans->kmeansWithKernels(kMeansData, params.k, params.kk);
            applyClusterResults(result, CopyPole);
            loggercontrol.logMessage("Applied kernel-based k-means clustering with " + 
                                   std::to_string(params.k) + " clusters and kernel size " + 
                                   std::to_string(params.kk), b);
            //copier.removeNoise(CopyPole, componenti);//Если не хочешь засорять COPY но не работает визуализация
        }
        
    if (params.s == "triangulate") {
    auto clusterCenters = getClusterCenters();
    
    // Добавляем точки A и B, если их ещё нет
PointD pointA(config.pointA_x, config.pointA_y);
PointD pointB(config.pointB_x, config.pointB_y);

if (std::find(clusterCenters.begin(), clusterCenters.end(), pointA) == clusterCenters.end()) {
    clusterCenters.push_back(pointA);
}
if (std::find(clusterCenters.begin(), clusterCenters.end(), pointB) == clusterCenters.end()) {
    clusterCenters.push_back(pointB);
}
    
    // Проверить, что точки внутри поля
    for (const auto& p : clusterCenters) {
        if (p.x < 0 || p.x >= config.fieldWidth || p.y < 0 || p.y >= config.fieldHeight) {
            std::cerr << "Точка (" << p.x << ", " << p.y << ") вне поля!\n";
            return;
        }
    }
    
    auto triangles = componentCalculator.bowyerWatson(clusterCenters, *p);
    this->lastTriangulation = triangles;
    voronoi.buildFromDelaunay(triangles, pathFinder, p);
    // Используйте актуальное изображение
    gnuplotInterface.plotDelaunay(triangles, p);
   gnuplotInterface.plotVoronoi(p, voronoi.edges, clusterCenters);
    loggercontrol.logMessage("Triangulation visualized", b);
}

        // Новая команда для поиска пути
   if (params.s == "find_path") {
        PointD start(config.pointA_x, config.pointA_y); // Добавлено
    PointD goal(config.pointB_x, config.pointB_y);  // Добавлено
    if (this->lastTriangulation.empty()) {
        std::cerr << "Error: Triangulation not performed yet! Run 'triangulate' first.\n";
        return;
    }
    auto path = pathFinder.findPathAStar(start, goal, this->lastTriangulation);
    // Сохранить путь в BMP
        std::vector<std::vector<double>> pathMatrix = p->field;//Странный мув!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        for (const auto& p : path) {
            int x = static_cast<int>(p.x);
            int y = static_cast<int>(p.y);
            if (x >= 0 && x < pathMatrix[0].size() && y >= 0 && y < pathMatrix.size()) {
                pathMatrix[y][x] = 255; // Отметить путь белым цветом
            }
        }
        
    gnuplotInterface.plotPath(path, p); // Вызов метода визуализации
    if (path.empty()) {
        std::cerr << "Путь не найден! Проверьте триангуляцию и параметры." << std::endl;
    }
          }
     }
};

class Interface {
private:
bool b = true;

public:
    Config& config;
    Logger& loggerinterface;
    Control& c;
    DispatcherParams params;
    
    Interface(Config& cfg, Logger& log, Control& c) : config(cfg), loggerinterface(log), c(c){      
        if (config.loggingInterfaceEnabled) {
            
                loggerinterface.logMessage("Logging Interface is enabled.", b);
                std::cout << "Logging Interface is enabled." << std::endl;
        } else {
        loggerinterface.logMessage("Logging Interface is disabled.", b);
        std::cout << "Logging Interface is disabled." << std::endl;
        b = false;
      }
      // Получаем размеры поля из конфигурации
    params.A = config.fieldWidth;  // Размеры из конфигурации
     params.B = config.fieldHeight;      
  }
    
    void print() {
        bool a;
        std::ifstream file;
        std::cout << "Hello, dear user, this program builds Gaussians.\nEnter commands from a text file (PRESS 0) or from the keyboard (PRESS 1)?" << std::endl;
        std::cin >> a;
        loggerinterface.logMessage("User chose input method: " + std::to_string(a), b);
        if (a == 0) {
            std::cout << "You will enter commands from a text file.\nEnter filename:" << std::endl;
            std::cin >> params.filename;
            loggerinterface.logMessage("Reading commands from file: " + params.filename, b);
            file.open(params.filename);
            params.filename = "output.bmp"; //УДАЛИТЬ
            if (!file) {
                std::cout << "File not found" << std::endl;
                loggerinterface.logMessage("Error: File not found.", b);
                return;
            }
        } else {
            std::cout << "You will enter commands from the keyboard" << std::endl;
            loggerinterface.logMessage("User chose to input commands from the keyboard.", b);
        }
        if (a == 0) {
            int n = 0;
            while (file >> params.s) {
                loggerinterface.logMessage("Received command: " + params.s, b);
         if (params.s == "help") {
            // Открываем файл help.txt для записи
            std::ofstream helpFile("/home/log/Gauss/results/docs/help.txt");
            if (helpFile.is_open()) {
                helpFile << R"(
# Terrain Navigation System

Программа для анализа рельефа местности, построения триангуляции Делоне, диаграмм Вороного и поиска оптимальных маршрутов с учетом препятствий.

## 📌 Основные функции
- Генерация/загрузка карты высот (формат BMP и GNUPLOT)
- Кластеризация объектов методом k-means
- Триангуляция Делоне с учетом высот
- Построение диаграммы Вороного
- Поиск пути с ограничениями по углу наклона тележки
  

## ⚙️ Системные требования

### Обязательные компоненты
1. **Компилятор C++17**  
   - `g++` (GCC) или `clang++`  
   - *Зачем*: Для сборки исходного кода программы

2. **Gnuplot 5.4+**  
   - *Что это*: Программа для построения графиков  
   - *Зачем*: Для визуализации карт высот, триангуляции и маршрутов  
   - Установка:  
     ```bash
     sudo apt install gnuplot  # Linux (Debian/Ubuntu)
     brew install gnuplot      # macOS (Homebrew)
     ```

### Опциональные компоненты
3. **CMake 3.12+**  
   - *Что это*: Система управления сборкой  
   - *Зачем*: Для упрощённой компиляции в разных ОС (если не используете прямой вызов g++)  
   - Установка:  
     ```bash
     sudo apt install cmake  # Linux
     ```

### Совместимость ОС
✅ **Полная поддержка**:  
- Linux (Ubuntu/Debian/Arch)  
- macOS (Intel/Apple Silicon)  

⚠️ **Ограниченная поддержка**:  
- Windows (требуется WSL2 или Cygwin)  
  - Рекомендуемый способ:  
    ```bash
    wsl --install -d Ubuntu
    ```

### Проверка установки
```bash
# Проверить версии компонентов
g++ --version
gnuplot --version
cmake --version
```

## 🚀 Запуск программы

### Способ 1: Ручная компиляция и запуск
```bash
# Переход в папку проекта
cd Gauss

# Компиляция
g++ -std=c++17 main.cpp -o terrain_navigator

# Запуск 
./terrain_navigator
```

### Способ 2: Автоматический скрипт (рекомендуется)
```bash
cd Gauss
chmod +x auto.sh  # Даем права на выполнение (только при первом запуске)
./auto.sh
```

### Способ 3: С CMake (опционально)
```bash
cd Gauss
mkdir build && cd build
cmake ..
make
./terrain_navigator
```

## 📂 Файлы проекта
```
~/Gauss/                                # Корневая папка проекта
├── archive                             # Старые версии
├── .gitignore                          # Игнорирует файл archive
├── LICENSE                             # Лицензия
├── README.md                           # Документация
├── src/
│   ├── config.txt                      # Основные параметры системы
│   ├── commands.txt                    # Пример командного файла
│   ├── main.cpp                        # Исходный код программы
│   └── auto.sh                         # Скрипт, запускающий main.cpp и читающий файл commands.txt
└── results/                            # Все результаты работы (внутри Gauss)
    ├── docs/
    │   ├── help.txt                    # Справочник
    │   ├── log_control.txt             # Логи Control
    │   └── log_interface.txt           # Логи Interface
    └── visualizations/                 # Все визуализации
        ├── binary_with_components.png  # Факторы компонент
        ├── delaunay_triangulation.png  # Триангуляция Делоне
        ├── voronoi_diagram.png         # Диаграмма Вороного
        ├── landscape.png               # 3D-вид поля (GNUPLOT)
        ├── path_plot.png               # Маршрут
        ├── output_kmeans.bmp           # K-means
        ├── output_kmeans_kern.bmp      # K-means с ядрами
        ├── slice.bmp                   # Бинаризированная карта
        └── output.bmp                  # Сгенерированная карта
```


## 🛠 Команды управления

| Команда            | Параметры              | Описание                                                                |
|--------------------|------------------------|-------------------------------------------------------------------------|
| help               | -                      | Создание файла с пояснением команд                                      |
| init               | -                      | Инициализация поля                                                      |
| g                  | x y sx sy h            | Создает гаусс                                                           |
| generate           | -                      | Складывает гауссы                                                       |
| gnuplot            | filename.png           | Рисует картинку в gnuplot                                               |
| bmp_write          | filename.bmp           | Создает черно-белую серую картинку BMP                                  |
| bmp_read           | filename.bmp           | Чтение BMP файла и инициализация поля новыми размерами                  |
| bin                | integer_number         | Срез: все, что выше или ниже integer_number - черное, остальное - белое |
| k_means            | k                      | Выделение k кластеров                                                   |
| k_means_kern       | kk                     | Алгоритм k-means с ядрами размера kk                                    |
| triangulate        | -                      | Построение триангуляции                                                 |
| find_path          | -                      | Построение маршрута между точками A и B                                 |
| end                | -                      | Это конец программы                                                     | 


## ⚠️ Важно
1. Всегда начинайте с команды init
2. Точки A/B задаются в config.txt
3. Маршрут будет найден не всегда!
4. Для триангуляции и построения пути, нужно чтобы количество компонент было больше 3
5. Если пользуетесь программой, то важно использовать ту же файловую структуру!


## 📜 Командный файл (примеры)
1) Если нужно прочитать данные с файла cat.bmp
```
init
help
bmp_read cat.bmp
gnuplot
bmp_write
bin 127
k_means 2
k_means_kern 3
triangulate
find_path
end
```
2) Если данные вводятся с помощью гаусов
```
init
help
g 99 50 25 25 200
g 50 50 20 20 200
g 200 50 20 20 200
g 50 200 20 20 200
g 189 180 20 20 200
g 130 130 20 20 200
generate
gnuplot
bmp_write
bin 127
k_means 5
triangulate
find_path
end
```

## 🔄 Development Workflow

### 1. Создаем ветку от актуального main
git checkout main
git pull origin main
git checkout -b feature/my-feature

### 2. Работаем над изменениями...
#### (создаем/редактируем файлы)

### 3. Пушим ветку
git add .
git commit -m "Описание изменений"
git push origin feature/my-feature

### 4. Создаем Pull Request (через веб-интерфейс GitHub)
### 5. Мержим PR (Squash and merge recommended)

### 6. Создаем тег (версию)
git checkout main
git pull origin main
git tag -a v1.1.0 -m "Добавлена новая функциональность"
git push origin v1.1.0

### 7. Удаляем ветку
git branch -d feature/my-feature
git push origin --delete feature/my-feature

## 📄 Лицензия
Этот проект лицензирован под MIT License. Вы можете свободно использовать, изменять и распространять код, при условии, что вы укажете автора.

Разрешенные действия:

   1. Использование кода в коммерческих и некоммерческих проектах
   2. Модификация кода
   3. Распространение кода

Обязательное условие: при использовании кода, пожалуйста, укажите ссылку на автора.

Developed with ❤️ by **DebugDestroy**  
[GitHub Profile](https://github.com/DebugDestroy)

)";
                helpFile.close(); // Закрываем файл
                loggerinterface.logMessage("Help file created successfully.", b);
            } else {
                loggerinterface.logMessage("Failed to create help file.", b);
            }
        }
   else if (params.s == "init") {
    // Проверяем, была ли уже вызвана команда init
    if (n != 0) {
        std::cout << "The init command has already been called.\nError\n";
        loggerinterface.logMessage("Error: Multiple init commands.", b);
        return;
    }
     n = 1;
    // Логируем и инициализируем поле
    loggerinterface.logMessage("Initializing field with size: " + std::to_string(params.A) + " x " + std::to_string(params.B), b);
    c.Dispetcher(params);
    loggerinterface.logMessage("Field initialized.", b);
} else if (n != 1) {
                     std::cout << "The init command was not use.\nError\n";
                     loggerinterface.logMessage("Error: The init command was not use.", b);
                     return;
             }
 else if (params.s == "g") {

                // Читаем параметры из файла
                file >> params.x >> params.y >> params.sx >> params.sy >> params.h;

                // Проверяем, были ли параметры указаны
                if (file.fail()) {
                    // Если не удалось прочитать, значит, используем значения по умолчанию
                    if (file.eof()) { // Желательно проверить конец файла, чтобы избежать повторного чтения
                        params.x = config.defaultX;
                        params.y = config.defaultY;
                        params.sx = config.defaultSx;
                        params.sy = config.defaultSy;
                        params.h = config.defaultH;
                    } else {
                        // Если произошла ошибка чтения или недостаточно параметров
                        // Нужно сбросить флаг ошибки для дальнейшего чтения
                        file.clear();

                        // Проверяем, какие параметры были прочитаны
                        if (!(file >> params.x)) params.x = config.defaultX; // Если не удалось получить x, берем по умолчанию
                        if (!(file >> params.y)) params.y = config.defaultY; // Если не удалось получить y, берем по умолчанию
                        if (!(file >> params.sx)) params.sx = config.defaultSx; // ...
                        if (!(file >> params.sy)) params.sy = config.defaultSy; // ...
                        if (!(file >> params.h)) params.h = config.defaultH; // ...
                    }
                }

                loggerinterface.logMessage("Adding Gaussian: x=" + std::to_string(params.x) +
                    ", y=" + std::to_string(params.y) + 
                    ", sx=" + std::to_string(params.sx) + 
                    ", sy=" + std::to_string(params.sy) + 
                    ", h=" + std::to_string(params.h), b);
                c.Dispetcher(params);
            } else if (params.s == "generate") {
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Generated values in the field.", b);
                } else if (params.s == "gnuplot") {
                    file >> params.filename;
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Called gnuplot:" + params.filename, b);
                } else if (params.s == "bmp_write") {
                    file >> params.filename;
                    file >> params.is_binary_image;
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Created BMP file:" + params.filename, b);
                } else if (params.s == "bmp_read") {
                    file >> params.filename; // Чтение имени файла для bmp_read
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Read BMP file: " + params.filename, b);
                    //params.filename = "output.bmp"; УДАЛИТЬ
                } else if (params.s == "bin") {
                    file >> params.slice;
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Slice applied: slice=" + std::to_string(params.slice), b);
                } else if (params.s == "wave") {
                    file >> params.noisy;
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Wave will be used", b);
                    loggerinterface.logMessage("Component amount = " + std::to_string(c.componenti.size()), b);
                } else if (params.s == "k_means") {
                        file >> params.k;
                        c.Dispetcher(params);
                        loggerinterface.logMessage("Cluster count", b);
                } else if (params.s == "k_means_kern") {
                    file >> params.kk;
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Cluster kern worked", b);
                }
                 else if (params.s == "triangulate") {
                    loggerinterface.logMessage("triangulate worke", b);
                    c.Dispetcher(params);
                }
                else if (params.s == "find_path") {
                    loggerinterface.logMessage("find_path worke", b);
                    c.Dispetcher(params);
                }
            }
        } else {
            int n = 0;
            while (true) {
                std::cout << "Enter command (help, init, g, generate, gnuplot, bmp_write, bmp_read, bin, k_means, triangulate, find_path, end):";
                std::cin >> params.s;
                std::cout << "\n";
                loggerinterface.logMessage("Received command: " + params.s, b);
        if (params.s == "help") {
            // Открываем файл help.txt для записи
            std::ofstream helpFile("/home/log/Gauss/results/docs/help.txt");
            if (helpFile.is_open()) {
                helpFile << R"(
# # Terrain Navigation System

Программа для анализа рельефа местности, построения триангуляции Делоне, диаграмм Вороного и поиска оптимальных маршрутов с учетом препятствий.

## 📌 Основные функции
- Генерация/загрузка карты высот (формат BMP и GNUPLOT)
- Кластеризация объектов методом k-means
- Триангуляция Делоне с учетом высот
- Построение диаграммы Вороного
- Поиск пути с ограничениями по углу наклона тележки
  

## ⚙️ Системные требования

### Обязательные компоненты
1. **Компилятор C++17**  
   - `g++` (GCC) или `clang++`  
   - *Зачем*: Для сборки исходного кода программы

2. **Gnuplot 5.4+**  
   - *Что это*: Программа для построения графиков  
   - *Зачем*: Для визуализации карт высот, триангуляции и маршрутов  
   - Установка:  
     ```bash
     sudo apt install gnuplot  # Linux (Debian/Ubuntu)
     brew install gnuplot      # macOS (Homebrew)
     ```

### Опциональные компоненты
3. **CMake 3.12+**  
   - *Что это*: Система управления сборкой  
   - *Зачем*: Для упрощённой компиляции в разных ОС (если не используете прямой вызов g++)  
   - Установка:  
     ```bash
     sudo apt install cmake  # Linux
     ```

### Совместимость ОС
✅ **Полная поддержка**:  
- Linux (Ubuntu/Debian/Arch)  
- macOS (Intel/Apple Silicon)  

⚠️ **Ограниченная поддержка**:  
- Windows (требуется WSL2 или Cygwin)  
  - Рекомендуемый способ:  
    ```bash
    wsl --install -d Ubuntu
    ```

### Проверка установки
```bash
# Проверить версии компонентов
g++ --version
gnuplot --version
cmake --version
```

## 🚀 Запуск программы

### Способ 1: Ручная компиляция и запуск
```bash
# Переход в папку проекта
cd Gauss

# Компиляция
g++ -std=c++17 main.cpp -o terrain_navigator

# Запуск 
./terrain_navigator
```

### Способ 2: Автоматический скрипт (рекомендуется)
```bash
cd Gauss
chmod +x auto.sh  # Даем права на выполнение (только при первом запуске)
./auto.sh
```

### Способ 3: С CMake (опционально)
```bash
cd Gauss
mkdir build && cd build
cmake ..
make
./terrain_navigator
```

## 📂 Файлы проекта
```
~/Gauss/                                # Корневая папка проекта
├── archive                             # Старые версии
├── .gitignore                          # Игнорирует файл archive
├── LICENSE                             # Лицензия
├── README.md                           # Документация
├── src/
│   ├── config.txt                      # Основные параметры системы
│   ├── commands.txt                    # Пример командного файла
│   ├── main.cpp                        # Исходный код программы
│   └── auto.sh                         # Скрипт, запускающий main.cpp и читающий файл commands.txt
└── results/                            # Все результаты работы (внутри Gauss)
    ├── docs/
    │   ├── help.txt                    # Справочник
    │   ├── log_control.txt             # Логи Control
    │   └── log_interface.txt           # Логи Interface
    └── visualizations/                 # Все визуализации
        ├── binary_with_components.png  # Факторы компонент
        ├── delaunay_triangulation.png  # Триангуляция Делоне
        ├── voronoi_diagram.png         # Диаграмма Вороного
        ├── landscape.png               # 3D-вид поля (GNUPLOT)
        ├── path_plot.png               # Маршрут
        ├── output_kmeans.bmp           # K-means
        ├── output_kmeans_kern.bmp      # K-means с ядрами
        ├── slice.bmp                   # Бинаризированная карта
        └── output.bmp                  # Сгенерированная карта
```


## 🛠 Команды управления

| Команда            | Параметры              | Описание                                                                |
|--------------------|------------------------|-------------------------------------------------------------------------|
| help               | -                      | Создание файла с пояснением команд                                      |
| init               | -                      | Инициализация поля                                                      |
| g                  | x y sx sy h            | Создает гаусс                                                           |
| generate           | -                      | Складывает гауссы                                                       |
| gnuplot            | filename.png           | Рисует картинку в gnuplot                                               |
| bmp_write          | filename.bmp           | Создает черно-белую серую картинку BMP                                  |
| bmp_read           | filename.bmp           | Чтение BMP файла и инициализация поля новыми размерами                  |
| bin                | integer_number         | Срез: все, что выше или ниже integer_number - черное, остальное - белое |
| k_means            | k                      | Выделение k кластеров                                                   |
| k_means_kern       | kk                     | Алгоритм k-means с ядрами размера kk                                    |
| triangulate        | -                      | Построение триангуляции                                                 |
| find_path          | -                      | Построение маршрута между точками A и B                                 |
| end                | -                      | Это конец программы                                                     | 


## ⚠️ Важно
1. Всегда начинайте с команды init
2. Точки A/B задаются в config.txt
3. Маршрут будет найден не всегда!
4. Для триангуляции и построения пути, нужно чтобы количество компонент было больше 3
5. Если пользуетесь программой, то важно использовать ту же файловую структуру!


## 📜 Командный файл (примеры)
1) Если нужно прочитать данные с файла cat.bmp
```
init
help
bmp_read cat.bmp
gnuplot
bmp_write
bin 127
k_means 2
k_means_kern 3
triangulate
find_path
end
```
2) Если данные вводятся с помощью гаусов
```
init
help
g 99 50 25 25 200
g 50 50 20 20 200
g 200 50 20 20 200
g 50 200 20 20 200
g 189 180 20 20 200
g 130 130 20 20 200
generate
gnuplot
bmp_write
bin 127
k_means 5
triangulate
find_path
end
```

## 🔄 Development Workflow

### 1. Создаем ветку от актуального main
git checkout main
git pull origin main
git checkout -b feature/my-feature

### 2. Работаем над изменениями...
#### (создаем/редактируем файлы)

### 3. Пушим ветку
git add .
git commit -m "Описание изменений"
git push origin feature/my-feature

### 4. Создаем Pull Request (через веб-интерфейс GitHub)
### 5. Мержим PR (Squash and merge recommended)

### 6. Создаем тег (версию)
git checkout main
git pull origin main
git tag -a v1.1.0 -m "Добавлена новая функциональность"
git push origin v1.1.0

### 7. Удаляем ветку
git branch -d feature/my-feature
git push origin --delete feature/my-feature


## 📄 Лицензия
Этот проект лицензирован под MIT License. Вы можете свободно использовать, изменять и распространять код, при условии, что вы укажете автора.

Разрешенные действия:

   1. Использование кода в коммерческих и некоммерческих проектах
   2. Модификация кода
   3. Распространение кода

Обязательное условие: при использовании кода, пожалуйста, укажите ссылку на автора.

Developed with ❤️ by **DebugDestroy**  
[GitHub Profile](https://github.com/DebugDestroy)

)";
                helpFile.close(); // Закрываем файл
                loggerinterface.logMessage("Help file created successfully.", b);
            } else {
                loggerinterface.logMessage("Failed to create help file.", b);
            }
        }
    if (params.s == "init") {
    // Проверяем, была ли уже вызвана команда init
    if (n != 0) {
        std::cout << "The init command has already started.\nError\n";
        loggerinterface.logMessage("Error: Multiple init commands.", b);
        return;
    }
    n = 1; // Устанавливаем флаг, что команда init была вызвана
    
    // Логируем и инициализируем поле
    loggerinterface.logMessage("Initializing field with size: " + std::to_string(params.A) + " x " + std::to_string(params.B), b);
    c.Dispetcher(params);
    loggerinterface.logMessage("Field initialized.", b);
}
  if (n != 1) {
                     std::cout << "The init command was not use.\nError\n";
                     loggerinterface.logMessage("Error: The init command was not use.", b);
                     return;
             }
             
  if (params.s == "g") {
    std::string input;
    std::getline(std::cin, input);

    std::istringstream inputStream(input);

    // Инициализируем переменные значениями по умолчанию
    params.x = config.defaultX;
    params.y = config.defaultY;
    params.sx = config.defaultSx;
    params.sy = config.defaultSy;
    params.h = config.defaultH;

    // Читаем введенные данные
    if (!(inputStream >> params.x)) {
        std::cout << "The default value for x is used: " << config.defaultX << std::endl;
    }
    if (!(inputStream >> params.y)) {
        std::cout << "The default value for y is used: " << config.defaultY << std::endl;
    }
    if (!(inputStream >> params.sx)) {
        std::cout << "The default value for sx is used: " << config.defaultSx << std::endl;
    }
    if (!(inputStream >> params.sy)) {
        std::cout << "The default value for sy is used: " << config.defaultSy << std::endl;
    }
    if (!(inputStream >> params.h)) {
        std::cout << "The default value for h is used: " << config.defaultH << std::endl;
    }

                loggerinterface.logMessage("Adding Gaussian: x=" + std::to_string(params.x) +
                    ", y=" + std::to_string(params.y) + 
                    ", sx=" + std::to_string(params.sx) + 
                    ", sy=" + std::to_string(params.sy) + 
                    ", h=" + std::to_string(params.h), b);
               c.Dispetcher(params);
            }

                  if (params.s == "generate") {
                    c.Dispetcher(params);
                    std::cout << "Generated values in the field." << std::endl;
                    loggerinterface.logMessage("Generated values in the field.", b);
                } if (params.s == "gnuplot") {
                    std::cout << "Enter the filename to draw:" << std::endl;
                    std::cin >> params.filename;//имя для файла
                    c.Dispetcher(params);
                    std::cout << "Called gnuplot:" + params.filename << std::endl;
                    loggerinterface.logMessage("Called gnuplot:" + params.filename, b);
                } if (params.s == "bmp_write") {
                    std::cout << "Enter the filename to draw:" << std::endl;
                    std::cin >> params.filename;//имя для файла
                    std::cout << "Image is binary (0) or not (1):" << std::endl;
                    std::cin >> params.is_binary_image;
                    c.Dispetcher(params);
                    std::cout << "Created BMP file:" + params.filename << std::endl;
                    loggerinterface.logMessage("Created BMP file:" + params.filename, b);
                } if (params.s == "bmp_read") {
                    std::cout << "Enter the filename to read:" << std::endl;
                    std::cin >> params.filename;
                    c.Dispetcher(params);
                    std::cout << "Read BMP file: " + params.filename << std::endl;
                    loggerinterface.logMessage("Read BMP file: " + params.filename, b);
                    //params.filename = "output.bmp"; УДАЛИТЬ
                }  if (params.s == "end") {
                    std::cout << "Ending the program" << std::endl;
                    loggerinterface.logMessage("Ending the program.", b);
                    break;
                }
                   if (params.s == "bin") {
                     std::cout << "Enter slice level:" << std::endl;
                     std::cin >> params.slice;
                     c.Dispetcher(params);
                     loggerinterface.logMessage("Slice applied: slice=" + std::to_string(params.slice), b);
                     std::cout << "Slice applied: slice=" << params.slice << std::endl;
                }
                
                    if (params.s == "wave") {
                        std::cout << "Enter noisy level for components:" << std::endl;
                        std::cin >> params.noisy;
                        c.Dispetcher(params);
                        loggerinterface.logMessage("Wave will be used", b);
                        loggerinterface.logMessage("Component amount = " + std::to_string(c.componenti.size()), b);
                }
                    if (params.s == "k_means") {
                        std::cout << "Enter amount cluster:" << std::endl;
                        std::cin >> params.k;
                        c.Dispetcher(params);
                        std::cout << "Cluster count" << std::endl;
                        loggerinterface.logMessage("Cluster count", b);
                }
                    if (params.s == "k_means_kern") {
                      std::cout << "Enter amount kern:" << std::endl;
                      std::cin >> params.kk;
                      c.Dispetcher(params);
                      std::cout << "Cluster kern worked" << std::endl;
                      loggerinterface.logMessage("Cluster kern worked", b);
                }
                    if (params.s == "triangulate") {
                    loggerinterface.logMessage("triangulate worke", b);
                    c.Dispetcher(params);
                    std::cout << "Triangulation completed" << std::endl;
                }
                if (params.s == "find_path") {
                    loggerinterface.logMessage("find_path worke", b);
                    c.Dispetcher(params);
                    std::cout << "Path search completed" << std::endl;
                }
            }

        if (file.is_open()) {
            file.close();
            loggerinterface.logMessage("Closed input file.", b);
        }
     }
   }
};

int main() {
    //Логирование
    Config config("config.txt");
    Logger loggerinterface(config.logFileNameInterface);
    Logger loggercontrol(config.logFileNameControl);
    // Создаем интерфейс
    Control c(config, loggercontrol);
    Interface i(config, loggerinterface, c);
    
    // Вызываем метод print() интерфейса
    i.print();

    return 0;
}
