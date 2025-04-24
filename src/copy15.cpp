/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2)Возможно новые требования
3) Оптимизация

Изменения: Улучшена визуализация (теперь окружности в трехмерном рисунке зависят от радиуса тележки, а в двумерном рисуется сечение окружности(если оно есть, иначе линии)) + 
очень сильно улучшено логирование + добавлены уровни логирования + упрощен интерфейс + исправил небольшой недочет в wave + убрал дублирование кода

                           Программа готова!!! Но всю равно следите за обновлениями на гитхаб [GitHub Profile](https://github.com/DebugDestroy)
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

enum class BmpWriteMode {//для bmp_write
    Binary,  // Бинарное изображение
    Full     // Полное поле
};

// Уровни логирования от самого детального до самого важного
enum class LogLevel {
    Trace,    // Максимальная детализация (трассировка выполнения)
    Debug,    // Отладочная информация
    Info,     // Информационные сообщения
    Warning,  // Предупреждения
    Error,    // Ошибки (не в алгоритме)
    Critical, // Критические ошибки (в алгоритме)
    Off       // Логирование отключено
};


enum class ThresholdMode {// Для bin
    All,     // и ямы, и горы (по модулю)
    Peaks,   // только горы (без модуля, значение > slice)
    Valleys  // только ямы (без модуля, значение < -slice)
};

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
    BmpWriteMode bmp_mode;  // Режим записи BMP
    int slice;//срез 
    ThresholdMode bin_mode;// Для выбора режима bin
    int noisy;//Уровень шума (все компоненты меньше - шум)
    int k;//количество кластеров для метода k_means
    int kk;//количество кластеров для метода k_means_with_kern
    double pointA_x, pointA_y, pointB_x, pointB_y;// точки A и B для путя
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
    static bool areCollinear(const PointD& a, const PointD& b, const PointD& c) {
    double area = (b.y - a.y) * (c.x - b.x) - (c.y - b.y) * (b.x - a.x);
    return std::fabs(area) < 1e-6; // Проверяем, близко ли значение к нулю
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
    
    double length() const {
        return std::hypot(a.x - b.x, a.y - b.y);
    }
};

struct VoronoiEdge {
    PointD start, end;
    VoronoiEdge(PointD s, PointD e) : start(s), end(e) {}
};

// Структура для хранения данных ячейки сетки
struct GridCell {
    PointD position;       // Координаты центра ячейки
    double height = 0;     // Высота в этой точке
    double slopeX = 0;     // Наклон по оси X (в радианах)
    double slopeY = 0;     // Наклон по оси Y (в радианах)
    bool visited = false;
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
    int fieldWidth, fieldHeight;
    double defaultX, defaultY, defaultSx, defaultSy, defaultH;
    std::string FiltrationLogLevelInterface, FiltrationLogLevelControl;
    std::string logFileNameInterface;
    std::string logFileNameControl;
    std::string defaultGnuplot, defaultPlotMetedata, defaultPlotVoronoi, defaultPlotDelaunay, defaultPlotPath, defaultWrite, defaultRead, defaultWriteModeImage, defaultBinMode, defaultPlot3DPath;
    int defaultSlice, defaultNoisy, defaultKlaster, defaultKlasterKern;
    double defaultpointA_x, defaultpointA_y, defaultpointB_x, defaultpointB_y;
    double vehicleRadius, maxSideAngle, maxUpDownAngle;

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
    else if (key == "defaultGnuplot") configFile >> defaultGnuplot;  
    else if (key == "defaultPlotMetedata") configFile >> defaultPlotMetedata;
    else if (key == "defaultPlotVoronoi") configFile >> defaultPlotVoronoi;
    else if (key == "defaultPlotDelaunay") configFile >> defaultPlotDelaunay;
    else if (key == "defaultPlotPath") configFile >> defaultPlotPath;   
    else if (key == "defaultWrite") configFile >> defaultWrite;
    else if (key == "defaultWriteModeImage") configFile >> defaultWriteModeImage;
    else if (key == "defaultRead") configFile >> defaultRead;
    else if (key == "defaultSlice") configFile >> defaultSlice;
    else if (key == "defaultBinMode") configFile >> defaultBinMode;
    else if (key == "defaultNoisy") configFile >> defaultNoisy;
    else if (key == "defaultKlaster") configFile >> defaultKlaster;
    else if (key == "defaultKlasterKern") configFile >> defaultKlasterKern;  
    else if (key == "logFileNameInterface") configFile >> logFileNameInterface;
    else if (key == "logFileNameControl") configFile >> logFileNameControl;
    else if (key == "FiltrationLogLevelInterface") configFile >> FiltrationLogLevelInterface;
    else if (key == "FiltrationLogLevelControl") configFile >> FiltrationLogLevelControl;
    else if (key == "defaultpointA_x") configFile >> defaultpointA_x;
    else if (key == "defaultpointA_y") configFile >> defaultpointA_y;
    else if (key == "defaultpointB_x") configFile >> defaultpointB_x;
    else if (key == "defaultpointB_y") configFile >> defaultpointB_y;
    else if (key == "vehicleRadius") configFile >> vehicleRadius;    
    else if (key == "maxSideAngle") configFile >> maxSideAngle;
    else if (key == "maxUpDownAngle") configFile >> maxUpDownAngle;
    else if (key == "defaultPlot3DPath") configFile >> defaultPlot3DPath;
}

        configFile.close();
    }
};

class Logger {
private:
    std::ofstream logFile;
    LogLevel currentLogLevel;
    std::string moduleName;
    
    // Преобразуем уровень логирования в строку
    static const char* levelToString(LogLevel level) {
        switch(level) {
            case LogLevel::Trace:    return "TRACE";
            case LogLevel::Debug:   return "DEBUG";
            case LogLevel::Info:    return "INFO";
            case LogLevel::Warning: return "WARNING";
            case LogLevel::Error:   return "ERROR";
            case LogLevel::Critical: return "CRITICAL";
            default:                return "UNKNOWN";
        }
    }
    
    // Преобразуем строку из конфига в LogLevel
    static LogLevel stringToLevel(const std::string& levelStr) {
        if (levelStr == "TRACE")    return LogLevel::Trace;
        if (levelStr == "DEBUG")    return LogLevel::Debug;
        if (levelStr == "INFO")     return LogLevel::Info;
        if (levelStr == "WARNING")  return LogLevel::Warning;
        if (levelStr == "ERROR")    return LogLevel::Error;
        if (levelStr == "CRITICAL") return LogLevel::Critical;
        if (levelStr == "OFF")      return LogLevel::Off;
        
        // По умолчанию INFO, если неизвестный уровень
        return LogLevel::Info;
    }

public:
    Logger(const std::string& fileName, const std::string& module, const std::string& configLevel)
        : moduleName(module), currentLogLevel(stringToLevel(configLevel)) {
        if (currentLogLevel != LogLevel::Off) {
            logFile.open(fileName, std::ios::out | std::ios::app);
            if (!logFile.is_open()) {
                std::cerr << "Failed to open log file: " << fileName << std::endl;
            }
        }
    }

    ~Logger() {
        if (logFile.is_open()) {
            logFile.close();
        }
    }
    
    LogLevel getLogLevel() const {
        return currentLogLevel;
    }
    
    void logMessage(LogLevel level, const std::string& message) {
        // Проверяем, нужно ли логировать это сообщение
        if (level < currentLogLevel || !logFile.is_open()) {
            return;
        }
        
        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        
        logFile << "[" << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << "] "
                << "[" << levelToString(level) << "] "
                << "[" << moduleName << "] "
                << message << std::endl;
    }
    
    // Удобные методы для каждого уровня
    void trace(const std::string& message) {
        logMessage(LogLevel::Trace, message);
    }
    
    void debug(const std::string& message) {
        logMessage(LogLevel::Debug, message);
    }
    
    void info(const std::string& message) {
        logMessage(LogLevel::Info, message);
    }
    
    void warning(const std::string& message) {
        logMessage(LogLevel::Warning, message);
    }
    
    void error(const std::string& message) {
        logMessage(LogLevel::Error, message);
    }
    
    void critical(const std::string& message) {
        logMessage(LogLevel::Critical, message);
    }
};

struct Color{ // структура для хранения rgb цвета
    uint8_t r, g, b;
};
    
//Классы для построения гаусса
struct Gaus { // Структура, хранящая Гауссы
public:
    double h, x0, y0, sigma_x, sigma_y;
    Gaus(double h, double x0, double y0, double sigma_x, double sigma_y) : h(h), x0(x0), y0(y0), sigma_x(sigma_x), sigma_y(sigma_y) {}
};

class Pole {
private:
    Logger& logger;

    void logResizeOperation(int old_rows, int old_cols, int new_rows, int new_cols) const {
        logger.debug(std::string("[Pole] Resizing field from ") +
                    std::to_string(old_rows) + std::string("x") + std::to_string(old_cols) + 
                    std::string(" to ") + std::to_string(new_rows) + std::string("x") + std::to_string(new_cols));
    }

public:
    std::vector<std::vector<double>> field;

    Pole(int A, int B, Logger& lg) : logger(lg) {
        logger.trace(std::string("[Pole] Constructor called with dimensions ") +
                   std::to_string(A) + std::string("x") + std::to_string(B));
        
        if (A <= 0 || B <= 0) {
            logger.error(std::string("[Pole] Invalid dimensions: ") +
                       std::to_string(A) + std::string("x") + std::to_string(B));
            throw std::invalid_argument("Field dimensions must be positive");
        }

        const size_t old_rows = field.size();
        const size_t old_cols = old_rows > 0 ? field[0].size() : 0;
        
        field.resize(A, std::vector<double>(B, 0));
        
        logger.info(std::string("[Pole] Field initialized successfully. Dimensions: ") +
                   std::to_string(field.size()) + std::string("x") + 
                   std::to_string(field.empty() ? 0 : field[0].size()));
    }
    
    void resize(int A, int B) {
        logger.trace("[Pole::resize] Starting resize operation");
        
        if (A <= 0 || B <= 0) {
            logger.error(std::string("[Pole::resize] Invalid target dimensions: ") +
                       std::to_string(A) + std::string("x") + std::to_string(B));
            throw std::invalid_argument("Field dimensions must be positive");
        }

        const size_t old_rows = field.size();
        const size_t old_cols = old_rows > 0 ? field[0].size() : 0;
        
        logResizeOperation(old_rows, old_cols, A, B);
        
        // Сохраняем данные при уменьшении размера
        if (A < static_cast<int>(old_rows)) {
            logger.debug(std::string("[Pole::resize] Truncating rows from ") + 
                       std::to_string(old_rows) + std::string(" to ") + std::to_string(A));
        }
        
        field.resize(A);
        
        for (auto& row : field) {
            if (B < static_cast<int>(row.size())) {
                logger.debug(std::string("[Pole::resize] Truncating columns from ") +
                           std::to_string(row.size()) + std::string(" to ") + std::to_string(B));
            }
            row.resize(B, 0);
        }

        logger.info(std::string("[Pole::resize] Resize completed. New dimensions: ") +
                   std::to_string(field.size()) + std::string("x") + 
                   std::to_string(field.empty() ? 0 : field[0].size()));
    }

    void logFieldStats() {
        logger.debug(std::string("[Pole] Current field statistics:\n") +
             std::string("  Dimensions: ") + std::to_string(field.size()) + "x" + 
             std::to_string(field.empty() ? 0 : field[0].size()) + "\n" +
             std::string("  Memory usage: ~") + 
             std::to_string(field.size() * (field.empty() ? 0 : field[0].size()) * sizeof(double)) + 
             " bytes");
    }
};

class Component {
private:
     Logger& logger;
    
    // Полный лог всех полей компонента
    void logComponentStats() {
        std::ostringstream oss;
        oss << "Full component stats:\n"
            << "  Bounds: x[" << min_x << ".." << max_x << "], y[" << min_y << ".." << max_y << "]\n"
            << "  Center: (" << center_x << ", " << center_y << ")\n"
            << "  Pixels: " << pixelCount << "\n"
            << "  Eigenvalues: " << eigenvalue1 << ", " << eigenvalue2 << "\n"
            << "  Eigenvector1: (" << eigenvec1_x << ", " << eigenvec1_y << ")\n"
            << "  Eigenvector2: (" << eigenvec2_x << ", " << eigenvec2_y << ")";
            
        logger.info("[Component] " + oss.str());
    }
    
public:
    std::vector<std::vector<double>> componenta;
    int min_x, min_y, max_x, max_y;
    double center_x, center_y;
    double eigenvec1_x, eigenvec1_y;
    double eigenvec2_x, eigenvec2_y;
    double eigenvalue1, eigenvalue2;
    int pixelCount;  // Поле для хранения количества пикселей
    
   // Основной конструктор (для готовой матрицы)
    Component(Logger& lg, const std::vector<std::vector<double>>& inputComponenta, int count) : logger(lg), componenta(inputComponenta), pixelCount(count) {
        logger.trace("[Component] Constructor started");
        calculate_metadata();
        logger.debug(std::string("[Component] Created with pixel count: ") + std::to_string(pixelCount));
        logComponentStats(); // Теперь печатает все поля
    }

    void calculate_metadata() {
    // Находим границы и центр
    min_x = componenta[0].size();
    min_y = componenta.size();
    max_x = 0;
    max_y = 0;
    double sum_x = 0, sum_y = 0;
    int count = 0;
    
    logger.trace("[Component::calculate_metadata] Starting calculate_metadata");
    logger.debug("[Component::calculate_metadata] Scanning component pixels...");
    for (int y = 0; y < componenta.size(); ++y) {
        for (int x = 0; x < componenta[0].size(); ++x) {
            if (componenta[y][x] == 255) {
                min_x = std::min(min_x, x);
                min_y = std::min(min_y, y);
                max_x = std::max(max_x, x);
                max_y = std::max(max_y, y);
                sum_x += x;
                sum_y += y;
                count++;
            }
        }
    }

    if (count == 0) {
        logger.warning("[Component::calculate_metadata] Component has no valid pixels!");
        center_x = center_y = 0;
        return;
    }

    center_x = sum_x / count;
    center_y = sum_y / count;
    
    logger.debug(std::string("[Component::calculate_metadata] Center calculated: (") + 
                    std::to_string(center_x) + std::string(", ") + 
                    std::to_string(center_y) + std::string(")"));

    // Корректный расчет ковариационной матрицы
    double cov_xx = 0, cov_xy = 0, cov_yy = 0;
    logger.trace("[Component::calculate_metadata] Calculating covariance matrix...");
    for (int y = min_y; y <= max_y; ++y) {
        for (int x = min_x; x <= max_x; ++x) {
            if (componenta[y][x] == 255) {  // Исправлено условие
                double dx = x - center_x;
                double dy = y - center_y;
                cov_xx += dx * dx;
                cov_xy += dx * dy;
                cov_yy += dy * dy;
            }
        }
    }

    cov_xx /= count;
    cov_xy /= count;
    cov_yy /= count;
    
    logger.debug(std::string("[Component::calculate_metadata] Covariance matrix: xx=") + std::to_string(cov_xx) +
                    std::string(", xy=") + std::to_string(cov_xy) +
                    std::string(", yy=") + std::to_string(cov_yy));
                    
    // Собственные значения
    double trace = cov_xx + cov_yy;
    double det = cov_xx * cov_yy - cov_xy * cov_xy;
    eigenvalue1 = (trace + sqrt(trace * trace - 4 * det)) / 2;
    eigenvalue2 = (trace - sqrt(trace * trace - 4 * det)) / 2;
    
    logger.debug(std::string("[Component::calculate_metadata] Eigenvalues calculated: λ1=") + 
                    std::to_string(eigenvalue1) + 
                    std::string(", λ2=") + std::to_string(eigenvalue2));
                    
    // Корректный расчет собственных векторов
    if (fabs(cov_xy) > 1e-10) {
        logger.trace("[Component::calculate_metadata] Calculating eigenvectors for non-diagonal matrix");
        // Первый собственный вектор (для eigenvalue1)
        eigenvec1_x = cov_xy;
        eigenvec1_y = eigenvalue1 - cov_xx;
        
        // Второй собственный вектор (для eigenvalue2)
        eigenvec2_x = cov_xy;
        eigenvec2_y = eigenvalue2 - cov_xx;
    } else {
        logger.trace("[Component::calculate_metadata] Using identity matrix for eigenvectors (diagonal case)");
        eigenvec1_x = 1; eigenvec1_y = 0;
        eigenvec2_x = 0; eigenvec2_y = 1;
    }

    // Нормализация векторов
    double norm1 = sqrt(eigenvec1_x * eigenvec1_x + eigenvec1_y * eigenvec1_y);
    double norm2 = sqrt(eigenvec2_x * eigenvec2_x + eigenvec2_y * eigenvec2_y);
    
    eigenvec1_x /= norm1; eigenvec1_y /= norm1;
    eigenvec2_x /= norm2; eigenvec2_y /= norm2;
    
    logger.debug(std::string("[Component::calculate_metadata] Normalized eigenvectors: ") +
                    "v1=(" + std::to_string(eigenvec1_x) + "," + std::to_string(eigenvec1_y) + "), " +
                    "v2=(" + std::to_string(eigenvec2_x) + "," + std::to_string(eigenvec2_y) + ")");
        
        logger.trace("[Component::calculate_metadata] calculate_metadata completed");
}
};

class Copier {
private:
    Logger& logger;

public:
    Copier(Logger& lg) : logger(lg) {}

    void removeNoise(std::vector<std::vector<double>>& field, const std::vector<Component>& components) {
        logger.trace("[Copier::removeNoise] Starting noise removal");
        
        // Логируем исходные параметры
             logger.debug(std::string("[Copier::removeNoise] Input parameters: ") +
             std::to_string(field.size()) + "x" + 
             std::to_string(field.empty() ? 0 : field[0].size()) + " field, " +
             std::to_string(components.size()) + " components");
             
        // Обнуляем поле
        size_t total_zeroed = 0;
        for (auto& row : field) {
            total_zeroed += row.size();
            std::fill(row.begin(), row.end(), 0.0);
        }
        logger.debug(std::string("[Copier::removeNoise] Zeroed ") + std::to_string(total_zeroed) + std::string(" pixels"));

        // Копируем только значимые компоненты
        size_t total_copied = 0;
        size_t components_processed = 0;
        
        for (const auto& comp : components) {
            const auto& compData = comp.componenta;
            size_t component_pixels = 0;
            
            for (size_t i = 0; i < compData.size(); ++i) {
                for (size_t j = 0; j < compData[i].size(); ++j) {
                    if (compData[i][j] == 255) {
                        field[i][j] = 255;
                        component_pixels++;
                        total_copied++;
                    }
                }
            }
            
            components_processed++;
            logger.debug(std::string("[Copier::removeNoise] Component #") + std::to_string(components_processed) + 
                       std::string(": copied ") + std::to_string(component_pixels) + std::string(" pixels"));
        }

        // Итоговая статистика
        logger.info(std::string("[Copier::removeNoise] Completed. Stats:\n") +
                  std::string("  Total components processed: ") + std::to_string(components_processed) + std::string("\n") +
                  std::string("  Total pixels copied: ") + std::to_string(total_copied) + std::string("\n") +
                  std::string("  Zeroed pixels: ") + std::to_string(total_zeroed - total_copied));
    }
};

class ColorGenerator {
private:
    static void logColorConversion(double hue, const std::array<int, 3>& rgb, Logger& logger) {
        logger.debug(std::string("[ColorGenerator] Generated color - HSL: ") + 
                    std::to_string(hue) + "°, RGB: (" +
                    std::to_string(rgb[0]) + ", " +
                    std::to_string(rgb[1]) + ", " +
                    std::to_string(rgb[2]) + ")");
    }

public:
    static std::vector<std::array<int, 3>> generateColors(int numColors, Logger& logger) {
        logger.trace(std::string("[ColorGenerator::generateColors] Starting color generation for ") + 
                   std::to_string(numColors) + " colors");

        if (numColors <= 0) {
            logger.error(std::string("[ColorGenerator::generateColors] Invalid number of colors: ") + 
                        std::to_string(numColors));
            throw std::invalid_argument("Number of colors must be greater than 0.");
        }

        std::vector<std::array<int, 3>> colors;
        unsigned int seed = static_cast<unsigned int>(time(0));
        srand(seed);
        logger.debug(std::string("[ColorGenerator] Seed initialized: ") + std::to_string(seed));

        for (int i = 0; i < numColors; ++i) {
            double hue = static_cast<double>(i) / numColors * 360.0 + rand() % 30;
            hue = fmod(hue, 360.0);

            // HSL to RGB conversion
            double r, g, b;
            double C = 1.0;
            double X = C * (1 - fabs(fmod((hue / 60.0), 2) - 1));
            double m = 0.0;

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

            // Scale to 50-255 range
            std::array<int, 3> rgb = {
                std::min(255, static_cast<int>(std::max(50.0, (r + m) * 255))),
                std::min(255, static_cast<int>(std::max(50.0, (g + m) * 255))),
                std::min(255, static_cast<int>(std::max(50.0, (b + m) * 255)))
            };

            colors.push_back(rgb);
            logColorConversion(hue, rgb, logger);
        }

        logger.info(std::string("[ColorGenerator::generateColors] Successfully generated ") + 
                   std::to_string(colors.size()) + " distinct colors");
        return colors;
    }
};

class GaussBuilder {
private:
    Logger& logger;

    void logGaussParameters(const Gaus& g) const {
        logger.debug(std::string("[GaussBuilder] Added Gaussian - ") +
                   "h: " + std::to_string(g.h) + 
                   ", x0: " + std::to_string(g.x0) +
                   ", y0: " + std::to_string(g.y0) +
                   ", σx: " + std::to_string(g.sigma_x) +
                   ", σy: " + std::to_string(g.sigma_y));
    }

public:
    GaussBuilder(Logger& lg) : logger(lg) {}

    void addgauss(double h, double x0, double y0, double sigma_x, double sigma_y, std::vector<Gaus>& gaussi) {
        logger.trace("[GaussBuilder::addgauss] Adding new Gaussian distribution");
        
        gaussi.emplace_back(h, x0, y0, sigma_x, sigma_y);
        logGaussParameters(gaussi.back());
    }
    
    void init(int A, int B, std::unique_ptr<Pole>& p) {
        logger.trace(std::string("[GaussBuilder::init] Initializing field with size ") +
                   std::to_string(A) + "x" + std::to_string(B));
        
        if (!p) {
            logger.debug("[GaussBuilder::init] Creating new Pole object");
            p = std::make_unique<Pole>(A, B, logger);
        } else {
            logger.debug("[GaussBuilder::init] Resizing existing Pole object");
            p->resize(A, B);
        }
        
        logger.info("[GaussBuilder::init] Field initialized successfully");
    }
    
    void generate(std::unique_ptr<Pole>& p, std::vector<Gaus>& gaussi) {
        logger.trace("[GaussBuilder::generate] Starting field generation");
        
        if (p == nullptr) {
            logger.error("[GaussBuilder::generate] Pole not initialized!");
            return;
        }
        
        const size_t width = p->field[0].size();
        const size_t height = p->field.size();
        logger.debug(std::string("[GaussBuilder::generate] Generating field ") +
                   std::to_string(width) + "x" + std::to_string(height) + 
                   " with " + std::to_string(gaussi.size()) + " Gaussians");

        // Заполнение базовым значением
        for (auto& row : p->field) {
            std::fill(row.begin(), row.end(), 127);
        }
        logger.debug("[GaussBuilder::generate] Base level set to 127");

        // Применение гауссиан
        size_t total_pixels = 0;
        size_t clamped_pixels = 0;
        
        for (const auto& g : gaussi) {
            logger.trace(std::string("[GaussBuilder::generate] Applying Gaussian at (") +
                       std::to_string(g.x0) + "," + std::to_string(g.y0) + ")");
            
            for (size_t x = 0; x < width; ++x) {
                for (size_t y = 0; y < height; ++y) {
                    double value = g.h * exp(-((pow((x - g.x0)/g.sigma_x, 2) + 
                                             pow((y - g.y0)/g.sigma_y, 2))/2));
                    p->field[y][x] += value;
                    
                    // Логирование clamping
                    if (p->field[y][x] <= 0.0 || p->field[y][x] >= 255.0) {
                        clamped_pixels++;
                    }
                    total_pixels++;
                    
                    p->field[y][x] = std::clamp(p->field[y][x], 0.0, 255.0);
                }
            }
        }

        logger.info(std::string("[GaussBuilder::generate] Generation completed\n") +
                  "  Total pixels processed: " + std::to_string(total_pixels) + "\n" +
                  "  Clamped pixels: " + std::to_string(clamped_pixels) + "\n" +
                  "  Clamping ratio: " + 
                  std::to_string((double)clamped_pixels/total_pixels*100) + "%");
    }
};
    
class BmpHandler {
private:
    Logger& logger;

    void logBmpHeader(const unsigned char* header) const {
        std::ostringstream oss;
        oss << "[BmpHandler] BMP Header Info:\n"
            << "  File Size: " << *reinterpret_cast<const uint32_t*>(&header[2]) << " bytes\n"
            << "  Image Size: " << *reinterpret_cast<const uint32_t*>(&header[34]) << " bytes\n"
            << "  Dimensions: " << *reinterpret_cast<const int32_t*>(&header[18]) 
            << "x" << *reinterpret_cast<const int32_t*>(&header[22]) << "\n"
            << "  Bits Per Pixel: " << *reinterpret_cast<const uint16_t*>(&header[28]);
        logger.debug(oss.str());
    }

public:
    BmpHandler(Logger& lg) : logger(lg) {}

    void bmp_write(const std::vector<std::vector<double>>& pixelMatrix, const std::string& filename) {
        logger.trace("[BmpHandler::bmp_write] Starting BMP write operation");
        
        if (pixelMatrix.empty() || pixelMatrix[0].empty()) {
            logger.error("[BmpHandler::bmp_write] Empty pixel matrix provided");
            return;
        }

        int width = pixelMatrix[0].size();
        int height = pixelMatrix.size();
        int padding = (4 - (width * 3) % 4) % 4;

        logger.debug(std::string("[BmpHandler::bmp_write] Creating BMP file: ") + filename + 
                   " with dimensions " + std::to_string(width) + "x" + std::to_string(height) +
                   ", padding: " + std::to_string(padding) + " bytes");

        std::ofstream bmpFile(filename, std::ios::binary);
        if (!bmpFile) {
            logger.error(std::string("[BmpHandler::bmp_write] Failed to create BMP file: ") + filename);
            return;
        }

        // Prepare and write header
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
    logger.trace("[BmpHandler::bmp_write] Writing BMP header");
        bmpFile.write(reinterpret_cast<char*>(bmpHeader), 54);

        // Write pixel data
        size_t pixels_written = 0;
        size_t clamped_pixels = 0;
        logger.trace("[BmpHandler::bmp_write] Writing pixel data (bottom-to-top)");
        
        for (int y = height - 1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                double pixelValue = pixelMatrix[y][x];
                if (pixelValue < 0.0 || pixelValue > 255.0) {
                    clamped_pixels++;
                    pixelValue = std::clamp(pixelValue, 0.0, 255.0);
                }
                unsigned char color = static_cast<unsigned char>(pixelValue);
                bmpFile.put(color).put(color).put(color);
                pixels_written++;
            }
            // Write padding
            for (int p = 0; p < padding; ++p) {
                bmpFile.put(0);
            }
        }

        bmpFile.close();
        logger.info(std::string("[BmpHandler::bmp_write] BMP file successfully written\n") +
                  "  Total pixels: " + std::to_string(pixels_written) + "\n" +
                  "  Clamped pixels: " + std::to_string(clamped_pixels) + "\n" +
                  "  File size: ~" + std::to_string(54 + (width*3 + padding)*height) + " bytes");
    }
    
   void bmp_read(GaussBuilder& gaussBuilder, const std::string& filename, 
                 std::vector<std::vector<double>>& pixelMatrix, std::unique_ptr<Pole>& p) {
        logger.trace(std::string("[BmpHandler::bmp_read] Starting BMP read operation: ") + filename);
        
        std::ifstream bmpFile(filename, std::ios::binary);
        if (!bmpFile) {
            logger.error(std::string("[BmpHandler::bmp_read] Failed to open BMP file: ") + filename);
            return;
        }

        // Read header
        unsigned char header[54];
        bmpFile.read(reinterpret_cast<char*>(header), 54);
        logBmpHeader(header);

        int width = header[18] | (header[19] << 8) | (header[20] << 16) | (header[21] << 24);
        int height = header[22] | (header[23] << 8) | (header[24] << 16) | (header[25] << 24);
        int padding = (4 - (width * 3) % 4) % 4;

        logger.debug(std::string("[BmpHandler::bmp_read] Image dimensions: ") + 
                   std::to_string(width) + "x" + std::to_string(height) +
                   ", padding: " + std::to_string(padding) + " bytes");

        // Initialize field
        logger.trace("[BmpHandler::bmp_read] Initializing field through GaussBuilder");
        gaussBuilder.init(height, width, p);
        pixelMatrix.resize(height, std::vector<double>(width));

        // Read pixel data
        size_t pixels_read = 0;
        double min_val = 255.0, max_val = 0.0;
        logger.trace("[BmpHandler::bmp_read] Reading pixel data (bottom-to-top)");
        
        for (int y = height - 1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                unsigned char color = bmpFile.get();
                bmpFile.get(); // G
                bmpFile.get(); // R
                
                double value = static_cast<double>(color);
                pixelMatrix[y][x] = value;
                p->field[y][x] = value;
                
                min_val = std::min(min_val, value);
                max_val = std::max(max_val, value);
                pixels_read++;
            }
            bmpFile.ignore(padding);
        }

        bmpFile.close();
        logger.info(std::string("[BmpHandler::bmp_read] BMP file successfully loaded\n") +
                  "  Pixels read: " + std::to_string(pixels_read) + "\n" +
                  "  Value range: [" + std::to_string(min_val) + ".." + 
                  std::to_string(max_val) + "]");
    }
};

// Класс для работы с сеткой рельефа
class TerrainGrid {
private:
    Logger& logger;

    void logGridInitialization(int width, int height) const {
        std::ostringstream oss;
        oss << "[TerrainGrid] Grid initialized:\n"
            << "  Dimensions: " << width << "x" << height << "\n"
            << "  Cell spacing: 1.0 units\n"
            << "  Y-axis inverted (top-left origin)";
        logger.debug(oss.str());
    }

    void logSlopeCalculation(int x, int y, double slopeX, double slopeY) const {
        if (logger.getLogLevel() <= LogLevel::Debug) {
            std::ostringstream oss;
            oss << "[TerrainGrid] Slope at (" << x << "," << y << "):\n"
                << "  slopeX: " << slopeX << " (" << atan(slopeX) * 180/M_PI << "°)\n"
                << "  slopeY: " << slopeY << " (" << atan(slopeY) * 180/M_PI << "°)";
            logger.debug(oss.str());
        }
    }

public:
    std::vector<std::vector<GridCell>> cells;

    TerrainGrid(int width, int height, Logger& lg) : logger(lg) {
        logger.trace("[TerrainGrid] Constructor called");
        initialize(width, height);
    }

    void initialize(int width, int height) {
        logger.trace("[TerrainGrid::initialize] Initializing grid");
        
        if (width <= 0 || height <= 0) {
            logger.error(std::string("[TerrainGrid::initialize] Invalid dimensions: ") + 
                       std::to_string(width) + "x" + std::to_string(height));
            throw std::invalid_argument("Grid dimensions must be positive");
        }

        cells.resize(height, std::vector<GridCell>(width));

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                cells[y][x].position = PointD(x + 0.5, (height - y - 1) + 0.5);
                
                if (x == 0 && y == 0) {  // Логируем только первую ячейку для примера
                    logger.trace(std::string("[TerrainGrid] Sample cell [0][0] position: (") + 
                               std::to_string(cells[y][x].position.x) + ", " +
                               std::to_string(cells[y][x].position.y) + ")");
                }
            }
        }

        logGridInitialization(width, height);
        logger.info("[TerrainGrid] Grid initialization completed");
    }
    
    std::pair<double, double> getEdgeSlopes(const Edge& edge) const {
        logger.trace(std::string("[TerrainGrid::getEdgeSlopes] Calculating slopes for edge: (") +
                   std::to_string(edge.a.x) + "," + std::to_string(edge.a.y) + ") -> (" +
                   std::to_string(edge.b.x) + "," + std::to_string(edge.b.y) + ")");

        PointD dir = {edge.b.x - edge.a.x, edge.b.y - edge.a.y};
        double length = std::hypot(dir.x, dir.y);
        
        if (length < 1e-6) {
            logger.debug("[TerrainGrid::getEdgeSlopes] Edge has zero length");
            return {0, 0};
        }

        dir.x /= length;
        dir.y /= length;
        int x = static_cast<int>((edge.a.x + edge.b.x) / 2);
        int y = static_cast<int>((edge.a.y + edge.b.y) / 2);

        if (x <= 0 || y <= 0 || x >= cells[0].size()-1 || y >= cells.size()-1) {
            logger.warning(std::string("[TerrainGrid::getEdgeSlopes] Edge center out of bounds: (") +
                         std::to_string(x) + "," + std::to_string(y) + ")");
            return {0, 0};
        }

        double forwardSlope = atan(cells[y][x].slopeX * dir.x + cells[y][x].slopeY * dir.y) * 180.0 / M_PI;
        double sideSlope = atan(cells[y][x].slopeX * (-dir.y) + cells[y][x].slopeY * dir.x) * 180.0 / M_PI;

        logger.debug(std::string("[TerrainGrid::getEdgeSlopes] Calculated slopes:\n") +
                   "  Forward: " + std::to_string(forwardSlope) + "°\n" +
                   "  Side: " + std::to_string(sideSlope) + "°");
        
        return {forwardSlope, sideSlope};
    }

    void calculateSlopes(const Pole& elevationData) {
        logger.trace("[TerrainGrid::calculateSlopes] Starting slope calculations");
        
        if (elevationData.field.empty()) {
            logger.error("[TerrainGrid::calculateSlopes] Empty elevation data provided");
            throw std::invalid_argument("Pole data is empty");
        }

        const int height = cells.size();
        const int width = cells[0].size();

        if (height != elevationData.field.size() || width != elevationData.field[0].size()) {
            logger.error(std::string("[TerrainGrid::calculateSlopes] Dimension mismatch: ") +
                       std::to_string(width) + "x" + std::to_string(height) + " vs " +
                       std::to_string(elevationData.field[0].size()) + "x" + 
                       std::to_string(elevationData.field.size()));
            throw std::invalid_argument("Grid and elevation data dimensions mismatch");
        }

        size_t processed_cells = 0;
        for (int y = 1; y < height-1; ++y) {
            for (int x = 1; x < width-1; ++x) {
                cells[y][x].height = elevationData.field[y][x];
                cells[y][x].slopeX = (elevationData.field[y][x+1] - elevationData.field[y][x-1]) / 2.0;
                cells[y][x].slopeY = (elevationData.field[y+1][x] - elevationData.field[y-1][x]) / 2.0;
                processed_cells++;

                if (x == 1 && y == 1) {  // Логируем только первую ячейку для примера
                    logSlopeCalculation(x, y, cells[y][x].slopeX, cells[y][x].slopeY);
                }
            }
        }

        logger.info(std::string("[TerrainGrid::calculateSlopes] Completed\n") +
                  "  Processed cells: " + std::to_string(processed_cells) + "\n" +
                  "  Border cells skipped: " + 
                  std::to_string(width*height - processed_cells));
    }

    void resetVisited() {
        size_t reset_count = 0;
        for (auto& row : cells) {
            for (auto& cell : row) {
                if (cell.visited) {
                    cell.visited = false;
                    reset_count++;
                }
            }
        }
        logger.debug(std::string("[TerrainGrid::resetVisited] Reset ") + 
                    std::to_string(reset_count) + " visited cells");
    }
};

class PathFinder {
private:
    const Config& config;
    Logger& logger;

    void logPoint(const std::string& prefix, const PointD& p) const {
        logger.debug(prefix + " (" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")");
    }

    void logTriangle(const std::string& prefix, const Triangle* tri) const {
        if (tri) {
            std::ostringstream oss;
            oss << prefix << " ["
                << "A(" << tri->a.x << "," << tri->a.y << "), "
                << "B(" << tri->b.x << "," << tri->b.y << "), "
                << "C(" << tri->c.x << "," << tri->c.y << ")]";
            logger.debug(oss.str());
        } else {
            logger.debug(prefix + " [не найден]");
        }
    }

    void logEdge(const Edge& edge) const {
        std::ostringstream oss;
        oss << "Ребро: (" << edge.a.x << "," << edge.a.y << ")-(" 
            << edge.b.x << "," << edge.b.y << "), длина: " << edge.length();
        logger.trace(oss.str());
    }

public:
    PathFinder(const Config& cfg, Logger& lg) : config(cfg), logger(lg) {
        logger.trace("[PathFinder] Инициализация поисковика пути");
    }
    
    double heuristic(const PointD& a, const PointD& b) {
        double dist = std::hypot(a.x - b.x, a.y - b.y);
        logger.trace(std::string("[PathFinder::heuristic] Расстояние между (") + 
                   std::to_string(a.x) + "," + std::to_string(a.y) + ") и (" +
                   std::to_string(b.x) + "," + std::to_string(b.y) + "): " + 
                   std::to_string(dist));
        return dist;
    }

    const Triangle* findContainingTriangle(const PointD& p, const std::vector<Triangle>& triangles) {
        logger.trace(std::string("[PathFinder::findContainingTriangle] Поиск треугольника для точки (") + 
                   std::to_string(p.x) + "," + std::to_string(p.y) + ")");
        
        for (const auto& tri : triangles) {
            if (isPointInTriangle(p, tri)) {
                logTriangle("[PathFinder::findContainingTriangle] Найден содержащий треугольник", &tri);
                return &tri;
            }
        }
        
        logger.warning("[PathFinder::findContainingTriangle] Точка не принадлежит ни одному треугольнику");
        return nullptr;
    }

    bool isPointInTriangle(const PointD& p, const Triangle& tri) {
        auto sign = [](PointD p1, PointD p2, PointD p3) {
            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
        };
        
        double d1 = sign(p, tri.a, tri.b);
        double d2 = sign(p, tri.b, tri.c);
        double d3 = sign(p, tri.c, tri.a);
        
        bool has_neg = (d1 < -1e-6) || (d2 < -1e-6) || (d3 < -1e-6);
        bool has_pos = (d1 > 1e-6) || (d2 > 1e-6) || (d3 > 1e-6);
        
        bool result = !(has_neg && has_pos);
        logger.trace(std::string("[PathFinder::isPointInTriangle] Точка (") + 
                   std::to_string(p.x) + "," + std::to_string(p.y) + ") " + 
                   (result ? "внутри" : "вне") + " треугольника");
        return result;
    }
    
    bool isNavigable(const Edge& edge, const TerrainGrid& terrainGrid) const {
        logEdge(edge);
        auto [forwardAngle, sideAngle] = terrainGrid.getEdgeSlopes(edge);
        
        bool navigable = abs(forwardAngle) <= config.maxUpDownAngle && 
               abs(sideAngle) <= config.maxSideAngle &&
               edge.length() > config.vehicleRadius * 2;
        
        if (!navigable) {
            std::ostringstream oss;
            oss << "[PathFinder::isNavigable] Ребро НЕ проходимо:\n"
                << "  Уклон вперед: " << forwardAngle << "° (макс. " << config.maxUpDownAngle << "°)\n"
                << "  Боковой уклон: " << sideAngle << "° (макс. " << config.maxSideAngle << "°)\n"
                << "  Длина: " << edge.length() << " (мин. " << config.vehicleRadius * 2 << ")";
            logger.debug(oss.str());
        }
        
        return navigable;
    }
    
    std::vector<const Triangle*> findNeighbors(const Triangle& tri, const std::vector<Triangle>& allTriangles) const {
        logger.trace("[PathFinder::findNeighbors] Поиск соседей треугольника");
        std::vector<const Triangle*> neighbors;
        
        for (const auto& other : allTriangles) {
            if (&tri == &other) continue;
            if (shareEdge(tri, other)) {
                neighbors.push_back(&other);
                logger.trace("[PathFinder::findNeighbors] Найден соседний треугольник");
            }
        }
        
        logger.debug(std::string("[PathFinder::findNeighbors] Найдено ") + std::to_string(neighbors.size()) + " соседей");
        return neighbors;
    }

public:
    std::vector<PointD> findPathAStar(const PointD& start, const PointD& goal, 
                                      const std::vector<Triangle>& triangles, 
                                      const TerrainGrid& terrainGrid,
                                      const std::vector<std::vector<double>>& binaryMap,
                                      double slice, const std::unique_ptr<Pole>& p) {
        logger.info("[PathFinder::findPathAStar] Начало поиска пути");
        logPoint("Стартовая точка:", start);
        logPoint("Целевая точка:", goal);
        
        std::vector<PointD> path;
        const Triangle* startTri = findContainingTriangle(start, triangles);
        const Triangle* goalTri = findContainingTriangle(goal, triangles);
        
        if (!startTri || !goalTri) {
            logger.error("[PathFinder::findPathAStar] Старт или цель вне триангуляции");
            logTriangle("Стартовый треугольник:", startTri);
            logTriangle("Целевой треугольник:", goalTri);
            return {};
        }

        if (startTri == goalTri) {
            logger.info("[PathFinder::findPathAStar] Старт и цель в одном треугольнике");
            path = {start, goal};
            if (!checkPathFeasibility(path, terrainGrid, binaryMap, slice, p)) {
                logger.warning("[PathFinder::findPathAStar] Прямой путь непроходим");
                return {};
            }
            logger.info("[PathFinder::findPathAStar] Прямой путь проходим");
            return path;
        }

        logger.debug("[PathFinder::findPathAStar] Запуск алгоритма A*");
        std::priority_queue<AStarNode> openSet;
        std::unordered_map<const Triangle*, double> costSoFar;

        PointD startCenter = startTri->calculateCircumcenter();
        openSet.emplace(startCenter, startTri);
        costSoFar[startTri] = 0;
        
        logger.debug("[PathFinder::findPathAStar] Начальный узел добавлен в очередь");

        while (!openSet.empty()) {
            AStarNode current = openSet.top();
            openSet.pop();

            if (current.tri == goalTri) {
                logger.info("[PathFinder::findPathAStar] Достигнут целевой треугольник");
                while (current.parent) {
                    path.push_back(current.position);
                    current = *current.parent;
                }
                std::reverse(path.begin(), path.end());
                
                if (!path.empty()) {
                    path.insert(path.begin(), start);
                    path.push_back(goal);
                    
                    logger.debug(std::string("[PathFinder::findPathAStar] Найден путь из ") + 
                               std::to_string(path.size()) + " точек");
                    
                    if (!checkPathFeasibility(path, terrainGrid, binaryMap, slice, p)) {
                        logger.warning("[PathFinder::findPathAStar] Путь непроходим");
                        return {};
                    }
                    logger.info("[PathFinder::findPathAStar] Путь проходим");
                }
                return path;
            }

            for (const auto& neighbor : getNeighbors(*current.tri, triangles)) {
                Edge edgeToCheck(current.position, neighbor->calculateCircumcenter());
                
                if (!isNavigable(edgeToCheck, terrainGrid)) {
                    logger.trace("[PathFinder::findPathAStar] Ребро непроходимо, пропускаем");
                    continue;
                }

                PointD neighborCenter = neighbor->calculateCircumcenter();
                double newCost = current.g + heuristic(current.position, neighborCenter);

                if (!costSoFar.count(neighbor) || newCost < costSoFar[neighbor]) {
                    costSoFar[neighbor] = newCost;
                    double priority = newCost + heuristic(neighborCenter, goal);
                    openSet.emplace(neighborCenter, neighbor, new AStarNode(current));
                    
                    logger.trace(std::string("[PathFinder::findPathAStar] Добавлен новый узел в очередь: (") +
                               std::to_string(neighborCenter.x) + "," + 
                               std::to_string(neighborCenter.y) + ")");
                }
            }
        }

        logger.warning("[PathFinder::findPathAStar] Путь не найден");
        return {};
    }

    bool checkPathFeasibility(const std::vector<PointD>& path,
                              const TerrainGrid& terrainGrid,
                              const std::vector<std::vector<double>>& binaryMap,
                              double slice, const std::unique_ptr<Pole>& p) const {
        logger.info("[PathFinder::checkPathFeasibility] Проверка проходимости пути");
        logger.debug(std::string("Сегментов пути: ") + std::to_string(path.size()-1));

        for (size_t i = 0; i < path.size() - 1; ++i) {
            std::ostringstream segInfo;
            segInfo << "Сегмент " << i+1 << ": (" 
                   << path[i].x << "," << path[i].y << ") -> (" 
                   << path[i+1].x << "," << path[i+1].y << ")";
            logger.trace(segInfo.str());
            
            const auto linePixels = bresenhamLine(path[i], path[i+1]);
            
            for (const auto& pixel : linePixels) {
                auto [forwardAngle, sideAngle] = getVehicleSlopeAngles(pixel, path[i+1], terrainGrid);
                
                if (forwardAngle > config.maxUpDownAngle || sideAngle > config.maxSideAngle) {
                    std::ostringstream oss;
                    oss << "Превышены углы наклона в точке (" << pixel.x << "," << pixel.y << "):\n"
                        << "  Уклон вперед: " << forwardAngle << "° (макс. " << config.maxUpDownAngle << "°)\n"
                        << "  Боковой уклон: " << sideAngle << "° (макс. " << config.maxSideAngle << "°)";
                    logger.debug(oss.str());
                    return false;
                }

                if (!isVehicleRadiusValid(pixel, binaryMap, p, slice)) {
                    logger.debug(std::string("Столкновение в точке (") + 
                               std::to_string(pixel.x) + "," + 
                               std::to_string(pixel.y) + ")");
                    return false;
                }
            }
        }
        
        logger.info("[PathFinder::checkPathFeasibility] Путь проходим");
        return true;
    }

std::vector<PointD> bresenhamLine(const PointD& start, const PointD& end) const {
    logger.trace(std::string("[PathFinder::bresenhamLine] Построение линии от (") + 
               std::to_string(start.x) + "," + std::to_string(start.y) + ") до (" +
               std::to_string(end.x) + "," + std::to_string(end.y) + ")");
    
    std::vector<PointD> linePoints;
    int x0 = static_cast<int>(start.x);
    int y0 = static_cast<int>(start.y);
    int x1 = static_cast<int>(end.x);
    int y1 = static_cast<int>(end.y);

    logger.debug(std::string("[PathFinder::bresenhamLine] Целочисленные координаты: (") + 
               std::to_string(x0) + "," + std::to_string(y0) + ") -> (" +
               std::to_string(x1) + "," + std::to_string(y1) + ")");

    int dx = abs(x1 - x0);
    int dy = -abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1;
    int sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;

    logger.trace(std::string("[PathFinder::bresenhamLine] Параметры алгоритма: dx=") + std::to_string(dx) + 
               ", dy=" + std::to_string(dy) + ", err=" + std::to_string(err));

    size_t pointCount = 0;
    while (true) {
        linePoints.emplace_back(x0, y0);
        pointCount++;
        
        if (x0 == x1 && y0 == y1) {
            logger.debug(std::string("[PathFinder::bresenhamLine] Линия завершена, точек: ") + std::to_string(pointCount));
            break;
        }
        
        int e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
            logger.trace(std::string("[PathFinder::bresenhamLine] Шаг по X: x=") + std::to_string(x0));
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
            logger.trace(std::string("[PathFinder::bresenhamLine] Шаг по Y: y=") + std::to_string(y0));
        }
    }
    
    return linePoints;
}

std::pair<double, double> getVehicleSlopeAngles(
    const PointD& pixel, 
    const PointD& nextPixel, 
    const TerrainGrid& grid
) const {
    logger.trace(std::string("[PathFinder::getVehicleSlopeAngles] Расчет углов наклона для точки (") +
               std::to_string(pixel.x) + "," + std::to_string(pixel.y) + ")");

    if (pixel == nextPixel) {
        logger.debug("[PathFinder::getVehicleSlopeAngles] Точки совпадают, углы = 0");
        return {0, 0};
    }

    int x = static_cast<int>(pixel.x);
    int y = static_cast<int>(pixel.y);
    
    if (x <= 0 || y <= 0 || x >= grid.cells[0].size()-1 || y >= grid.cells.size()-1) {
        logger.warning(std::string("[PathFinder::getVehicleSlopeAngles] Точка на границе сетки (") +
                     std::to_string(x) + "," + std::to_string(y) + "), помечаем как непроходимую");
        return {config.maxUpDownAngle + 1, config.maxSideAngle + 1};
    }

    PointD dir = {nextPixel.x - pixel.x, nextPixel.y - pixel.y};
    double length = std::hypot(dir.x, dir.y);
    
    if (length < 1e-6) {
        logger.debug("[PathFinder::getVehicleSlopeAngles] Нулевой вектор направления");
        return {0, 0};
    }
    
    dir.x /= length;
    dir.y /= length;
    logger.trace(std::string("[PathFinder::getVehicleSlopeAngles] Нормализованный вектор направления: (") +
               std::to_string(dir.x) + "," + std::to_string(dir.y) + ")");

    double z = grid.cells[y][x].height;
    double z_right = grid.cells[y][x+1].height;
    double z_left = grid.cells[y][x-1].height;
    double z_top = grid.cells[y-1][x].height;
    double z_bottom = grid.cells[y+1][x].height;

    logger.trace(std::string("[PathFinder::getVehicleSlopeAngles] Высоты вокруг точки: центр=") + std::to_string(z) +
               ", слева=" + std::to_string(z_left) + ", справа=" + std::to_string(z_right) +
               ", сверху=" + std::to_string(z_top) + ", снизу=" + std::to_string(z_bottom));

    double dzdx = (z_right - z_left) / 2.0;
    double dzdy = (z_bottom - z_top) / 2.0;
    logger.trace(std::string("[PathFinder::getVehicleSlopeAngles] Градиенты: dzdx=") + std::to_string(dzdx) +
               ", dzdy=" + std::to_string(dzdy));

    double forwardAngle = std::atan(dzdx * dir.x + dzdy * dir.y) * 180.0 / M_PI;
    PointD perp = {-dir.y, dir.x};
    double sideAngle = std::atan(dzdx * perp.x + dzdy * perp.y) * 180.0 / M_PI;

    logger.debug(std::string("[PathFinder::getVehicleSlopeAngles] Результат: forward=") + std::to_string(forwardAngle) +
               "°, side=" + std::to_string(sideAngle) + "°");
    
    return {std::abs(forwardAngle), std::abs(sideAngle)};
}

bool isVehicleRadiusValid(const PointD& pixel, 
                         const std::vector<std::vector<double>>& binaryMap,
                         const std::unique_ptr<Pole>& elevationData,
                         double sliceLevel) const {
    logger.trace(std::string("[PathFinder::isVehicleRadiusValid] Проверка радиуса в точке (") +
               std::to_string(pixel.x) + "," + std::to_string(pixel.y) + ")");
    
    const double x = pixel.x;
    const double y = pixel.y;
    const double currentHeight = elevationData->field[static_cast<int>(y)][static_cast<int>(x)];
    
    logger.debug(std::string("[PathFinder::isVehicleRadiusValid] Высота в точке: ") + std::to_string(currentHeight) +
               ", уровень среза: " + std::to_string(sliceLevel));

    const double heightDiff = std::abs(currentHeight - sliceLevel);
    if (heightDiff >= config.vehicleRadius) {
        logger.debug(std::string("[PathFinder::isVehicleRadiusValid] Высота вне радиуса (") + 
                   std::to_string(heightDiff) + " >= " + std::to_string(config.vehicleRadius) + 
                   "), столкновений нет");
        return true;
    }

    const double effectiveRadius = std::sqrt(config.vehicleRadius * config.vehicleRadius - heightDiff * heightDiff);
    const int radiusPixels = static_cast<int>(std::ceil(effectiveRadius));
    
    logger.trace(std::string("[PathFinder::isVehicleRadiusValid] Эффективный радиус: ") + 
               std::to_string(effectiveRadius) + " пикселей (" + 
               std::to_string(radiusPixels) + " целых)");

    const double radiusSquared = effectiveRadius * effectiveRadius;
    size_t checkedPixels = 0;
    size_t collisionPixels = 0;

    for (int dy = -radiusPixels; dy <= radiusPixels; ++dy) {
        for (int dx = -radiusPixels; dx <= radiusPixels; ++dx) {
            if (dx*dx + dy*dy > radiusSquared) continue;
            
            const int nx = static_cast<int>(x) + dx;
            const int ny = static_cast<int>(y) + dy;
            checkedPixels++;
            
            if (nx >= 0 && ny >= 0 && 
                nx < static_cast<int>(binaryMap[0].size()) && 
                ny < static_cast<int>(binaryMap.size()) &&
                binaryMap[ny][nx] == 255) {
                collisionPixels++;
                logger.trace(std::string("[PathFinder::isVehicleRadiusValid] Столкновение в (") +
                           std::to_string(nx) + "," + std::to_string(ny) + ")");
            }
        }
    }

    bool result = (collisionPixels == 0);
    logger.debug(std::string("[PathFinder::isVehicleRadiusValid] Результат проверки: ") + 
               std::string(result ? "проходимо" : "столкновение") + 
               ", проверено пикселей: " + std::to_string(checkedPixels) +
               ", столкновений: " + std::to_string(collisionPixels));
    
    return result;
}

   std::vector<const Triangle*> getNeighbors(const Triangle& tri, const std::vector<Triangle>& triangles) {
    logger.trace(std::string("[PathFinder::getNeighbors] Поиск соседей треугольника [") +
               "A(" + std::to_string(tri.a.x) + "," + std::to_string(tri.a.y) + "), " +
               "B(" + std::to_string(tri.b.x) + "," + std::to_string(tri.b.y) + "), " +
               "C(" + std::to_string(tri.c.x) + "," + std::to_string(tri.c.y) + ")]");
    
    std::vector<const Triangle*> neighbors;
    size_t totalChecked = 0;
    size_t edgesMatched = 0;

    for (const auto& other : triangles) {
        if (&tri == &other) {
            logger.trace("[PathFinder::getNeighbors] Пропуск самосравнения");
            continue;
        }
        
        totalChecked++;
        if (shareEdge(tri, other)) {
            neighbors.push_back(&other);
            edgesMatched++;
            logger.debug(std::string("[PathFinder::getNeighbors] Найдено общее ребро с треугольником [") +
                       "A(" + std::to_string(other.a.x) + "," + std::to_string(other.a.y) + "), " +
                       "B(" + std::to_string(other.b.x) + "," + std::to_string(other.b.y) + "), " +
                       "C(" + std::to_string(other.c.x) + "," + std::to_string(other.c.y) + ")]");
        }
    }

    logger.info(std::string("[PathFinder::getNeighbors] Результаты поиска: ") +
              std::to_string(neighbors.size()) + " соседей из " +
              std::to_string(totalChecked) + " проверенных треугольников, " +
              std::to_string(edgesMatched) + " общих ребер");
    
    return neighbors;
}

bool shareEdge(const Triangle& a, const Triangle& b) const {
    logger.trace("[PathFinder::shareEdge] Проверка общего ребра между треугольниками");
    
    const std::array<Edge, 3> edgesA = { Edge(a.a, a.b), Edge(a.b, a.c), Edge(a.c, a.a) };
    const std::array<Edge, 3> edgesB = { Edge(b.a, b.b), Edge(b.b, b.c), Edge(b.c, b.a) };

    for (size_t i = 0; i < 3; ++i) {
        const Edge& edgeA = edgesA[i];
        for (size_t j = 0; j < 3; ++j) {
            const Edge& edgeB = edgesB[j];
            if (edgeA == edgeB) {
                logger.debug(std::string("[PathFinder::shareEdge] Найдено общее ребро: (") +
                           std::to_string(edgeA.a.x) + "," + std::to_string(edgeA.a.y) + ")-(" +
                           std::to_string(edgeA.b.x) + "," + std::to_string(edgeA.b.y) + ")");
                return true;
            }
        }
    }

    logger.trace("[PathFinder::shareEdge] Общих ребер не обнаружено");
    return false;
}

bool otherHasEdge(const Triangle& other, const Edge& edge) const {
    logger.trace("[PathFinder::otherHasEdge] Проверка наличия ребра в треугольнике");
    
    const Edge edges[3] = { Edge(other.a, other.b), Edge(other.b, other.c), Edge(other.c, other.a) };
    
    for (int i = 0; i < 3; ++i) {
        if (edges[i] == edge) {
            logger.debug(std::string("[PathFinder::otherHasEdge] Ребро найдено: (") +
                       std::to_string(edge.a.x) + "," + std::to_string(edge.a.y) + ")-(" +
                       std::to_string(edge.b.x) + "," + std::to_string(edge.b.y) + ")");
            return true;
        }
    }

    logger.trace("[PathFinder::otherHasEdge] Ребро отсутствует в треугольнике");
    return false;
}
};

class VoronoiDiagram {
private:
    Logger& logger;

    void logEdge(const Edge& edge, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Ребро: (" << edge.a.x << "," << edge.a.y << ")-(" 
            << edge.b.x << "," << edge.b.y << "), длина: " << edge.length();
        logger.trace(oss.str());
    }

    void logPoint(const PointD& p, const std::string& prefix = "") const {
        logger.trace(prefix + "Точка (" + std::to_string(p.x) + "," + std::to_string(p.y) + ")");
    }

public:
    VoronoiDiagram(Logger& lg) : logger(lg) {
        logger.trace("[VoronoiDiagram] Инициализация диаграммы Вороного");
    }

    void buildFromDelaunay(const std::vector<Triangle>& triangles, 
                          const PathFinder& pathFinder, 
                          const std::unique_ptr<Pole>& p,
                          std::vector<VoronoiEdge>& edges) {
        logger.info("[VoronoiDiagram::buildFromDelaunay] Построение диаграммы из триангуляции Делоне");
        
        const int height = p->field.size();
        const int width = p->field[0].size();
        edges.clear();

        if (!p) {
            logger.error("[VoronoiDiagram::buildFromDelaunay] Ошибка: данные высот не инициализированы!");
            return;
        }

        logger.debug(std::string("[VoronoiDiagram::buildFromDelaunay] Размер области: ") + 
                   std::to_string(width) + "x" + std::to_string(height) + 
                   ", треугольников: " + std::to_string(triangles.size()));

        for (const auto& tri : triangles) {
            PointD cc1 = tri.calculateCircumcenter();
            bool cc1_valid = isPointInsideField(cc1, width, height);
            
            logger.debug(std::string("[VoronoiDiagram::buildFromDelaunay] Обработка треугольника с центром (") +
                       std::to_string(cc1.x) + "," + std::to_string(cc1.y) + "), " +
                       (cc1_valid ? "внутри" : "снаружи") + " области");

            auto neighbors = pathFinder.findNeighbors(tri, triangles);
            logger.trace(std::string("[VoronoiDiagram::buildFromDelaunay] Найдено ") + 
                       std::to_string(neighbors.size()) + " соседей");

            // Обработка обычных ребер
            for (const auto& neighbor : neighbors) {
                PointD cc2 = neighbor->calculateCircumcenter();
                bool cc2_valid = isPointInsideField(cc2, width, height);

                logger.trace(std::string("[VoronoiDiagram::buildFromDelaunay] Соседний центр (") +
                           std::to_string(cc2.x) + "," + std::to_string(cc2.y) + "), " +
                           (cc2_valid ? "внутри" : "снаружи") + " области");

                if (cc1_valid && cc2_valid) {
                    edges.emplace_back(cc1, cc2);
                    logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлено полное ребро Вороного");
                } 
                else if (cc1_valid || cc2_valid) {
                    auto clipped = clipEdgeToField(cc1, cc2, width, height);
                    if (!clipped.empty()) {
                        edges.emplace_back(clipped[0], clipped[1]);
                        logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлено обрезанное ребро Вороного");
                    } else {
                        logger.trace("[VoronoiDiagram::buildFromDelaunay] Ребро полностью вне области после отсечения");
                    }
                }
            }

            // Обработка граничных треугольников
            if (neighbors.size() < 3) {
                logger.trace("[VoronoiDiagram::buildFromDelaunay] Обработка граничного треугольника");
                handleBoundaryTriangle(tri, cc1, width, height, edges, cc1_valid, triangles, pathFinder);
            }
        }

        logger.info(std::string("[VoronoiDiagram::buildFromDelaunay] Построение завершено, ребер: ") + 
                  std::to_string(edges.size()));
    }

private:
    std::vector<PointD> clipEdgeToField(const PointD& p1, const PointD& p2, int width, int height) {
        logger.trace(std::string("[VoronoiDiagram::clipEdgeToField] Отсечение ребра (") +
                   std::to_string(p1.x) + "," + std::to_string(p1.y) + ")-(" +
                   std::to_string(p2.x) + "," + std::to_string(p2.y) + ")");

        auto code = [&](const PointD& p) {
            int c = 0;
            if (p.x < 0) c |= 1;
            if (p.x > width) c |= 2;
            if (p.y < 0) c |= 4;
            if (p.y > height) c |= 8;
            return c;
        };

        int code1 = code(p1);
        int code2 = code(p2);
        PointD a = p1, b = p2;

        logger.trace(std::string("[VoronoiDiagram::clipEdgeToField] Коды: p1=") + std::to_string(code1) + 
                   ", p2=" + std::to_string(code2));

        while (true) {
            if (!(code1 | code2)) {
                logger.trace("[VoronoiDiagram::clipEdgeToField] Ребро полностью внутри области");
                return {a, b};
            }
            if (code1 & code2) {
                logger.trace("[VoronoiDiagram::clipEdgeToField] Ребро полностью снаружи области");
                return {};
            }

            int outcode = code1 ? code1 : code2;
            PointD p;

            if (outcode & 8) {
                p.x = a.x + (b.x - a.x) * (height - a.y) / (b.y - a.y);
                p.y = height;
                logger.trace("[VoronoiDiagram::clipEdgeToField] Пересечение с верхней границей");
            }
            else if (outcode & 4) {
                p.x = a.x + (b.x - a.x) * (-a.y) / (b.y - a.y);
                p.y = 0;
                logger.trace("[VoronoiDiagram::clipEdgeToField] Пересечение с нижней границей");
            }
            else if (outcode & 2) {
                p.y = a.y + (b.y - a.y) * (width - a.x) / (b.x - a.x);
                p.x = width;
                logger.trace("[VoronoiDiagram::clipEdgeToField] Пересечение с правой границей");
            }
            else if (outcode & 1) {
                p.y = a.y + (b.y - a.y) * (-a.x) / (b.x - a.x);
                p.x = 0;
                logger.trace("[VoronoiDiagram::clipEdgeToField] Пересечение с левой границей");
            }

            if (outcode == code1) {
                a = p;
                code1 = code(a);
            } else {
                b = p;
                code2 = code(b);
            }
        }
    }

    bool isPointInsideField(const PointD& p, int width, int height) {
        bool inside = p.x >= 0 && p.y >= 0 && p.x < width && p.y < height;
        logger.trace(std::string("[VoronoiDiagram::isPointInsideField] Точка (") +
                    std::to_string(p.x) + "," + std::to_string(p.y) + ") " +
                    (inside ? "внутри" : "снаружи") + " области");
        return inside;
    }

    void handleBoundaryTriangle(const Triangle& tri, const PointD& cc, int width, int height, 
                               std::vector<VoronoiEdge>& edges, bool cc_valid,
                               const std::vector<Triangle>& allTriangles,  
                               const PathFinder& pathFinder) {
        if (!cc_valid) {
            logger.trace("[VoronoiDiagram::handleBoundaryTriangle] Центр снаружи, пропуск");
            return;
        }

        auto boundaryEdges = getBoundaryEdges(tri, allTriangles, pathFinder);
        logger.debug(std::string("[VoronoiDiagram::handleBoundaryTriangle] Найдено ") +
                   std::to_string(boundaryEdges.size()) + " граничных ребер");

        for (const auto& edge : boundaryEdges) {
            PointD boundaryPoint = calculateBoundaryIntersection(cc, edge, width, height);
            edges.emplace_back(cc, boundaryPoint);
            logger.debug("[VoronoiDiagram::handleBoundaryTriangle] Добавлено граничное ребро Вороного");
        }
    }

    std::vector<Edge> getBoundaryEdges(const Triangle& tri, 
                                      const std::vector<Triangle>& allTriangles,  
                                      const PathFinder& pathFinder) {
        std::vector<Edge> boundaryEdges;
        logger.trace("[VoronoiDiagram::getBoundaryEdges] Поиск граничных ребер треугольника");

        for (const auto& edge : { Edge(tri.a, tri.b), Edge(tri.b, tri.c), Edge(tri.c, tri.a) }) {
            bool isBoundary = true;
            logEdge(edge, "Проверка ребра: ");

            for (const auto& other : allTriangles) {
                if (&tri == &other) continue;
                if (pathFinder.shareEdge(tri, other) && pathFinder.otherHasEdge(other, edge)) {
                    isBoundary = false;
                    break;
                }
            }

            if (isBoundary) {
                boundaryEdges.push_back(edge);
                logger.trace("[VoronoiDiagram::getBoundaryEdges] Найдено граничное ребро");
            }
        }

        return boundaryEdges;
    }

    PointD calculateBoundaryIntersection(const PointD& circumCenter, 
                                        const Edge& boundaryEdge, 
                                        int width, int height) {
        logger.trace("[VoronoiDiagram::calculateBoundaryIntersection] Вычисление пересечения с границей");

        PointD edgeDir = { boundaryEdge.b.x - boundaryEdge.a.x, boundaryEdge.b.y - boundaryEdge.a.y };
        PointD normal = { -edgeDir.y, edgeDir.x };
        double length = std::hypot(normal.x, normal.y);
        
        if (std::fabs(length) < 1e-9) {
            logger.warning("[VoronoiDiagram::calculateBoundaryIntersection] Нулевая длина нормали!");
            return circumCenter;
        }
        
        normal.x /= length;
        normal.y /= length;

        double t = std::numeric_limits<double>::max();
        if (std::fabs(normal.x) > 1e-6) {
            t = std::min(t, (width - circumCenter.x) / normal.x);
            t = std::min(t, -circumCenter.x / normal.x);
        }
        if (std::fabs(normal.y) > 1e-6) {
            t = std::min(t, (height - circumCenter.y) / normal.y);
            t = std::min(t, -circumCenter.y / normal.y);
        }
        
        PointD result = {
            circumCenter.x + normal.x * t,
            circumCenter.y + normal.y * t
        };

        logger.debug(std::string("[VoronoiDiagram::calculateBoundaryIntersection] Точка пересечения: (") +
                    std::to_string(result.x) + "," + std::to_string(result.y) + ")");
        
        return result;
    }
};
   
class GnuplotInterface {
private:
    Logger& logger;

    double transformY(double y, int height) const {
        double transformed = height - y - 1;
        logger.trace(std::string("[GnuplotInterface::transformY] Преобразование Y: ") + 
                   std::to_string(y) + " -> " + std::to_string(transformed));
        return transformed;
    }

    void logPlotStart(const std::string& plotType, const std::string& filename) const {
        logger.info(std::string("[GnuplotInterface] Начало визуализации: ") + plotType);
        logger.debug(std::string("Файл вывода: ") + filename);
    }

    void logPlotEnd(const std::string& plotType) const {
        logger.info(std::string("[GnuplotInterface] Визуализация завершена: ") + plotType);
    }

public:
    GnuplotInterface(Logger& lg) : logger(lg) {
        logger.trace("[GnuplotInterface] Инициализация интерфейса Gnuplot");
    }

    void plotBinaryWithComponents(const std::vector<std::vector<double>>& CopyPole, 
                                const std::vector<Component>& components, 
                                const std::string& filename) {
        logPlotStart("BinaryWithComponents", filename);
        
        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotBinaryWithComponents] Ошибка открытия gnuplot pipe");
            return;
        }

        const int height = CopyPole.size();
        const int width = CopyPole[0].size();
        logger.debug("Размер данных: " + std::to_string(width) + "x" + std::to_string(height) + 
                   ", компонентов: " + std::to_string(components.size()));

        // Настройки графика
        fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200\n");
        fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
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
        fprintf(gnuplotPipe, "%f ", CopyPole[y][x]);
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
    double cy = transformY(comp.center_y, height);
    
    // Масштабируем собственные вектора по собственным значениям
    double scale_factor = 1; // Общий масштаб для визуализации
    double vec1_scale = scale_factor * sqrt(comp.eigenvalue1);
    double vec2_scale = scale_factor * sqrt(comp.eigenvalue2);
    
    // Первый собственный вектор
    fprintf(gnuplotPipe, "%f %f %f %f\n", 
            comp.center_x,
            cy,
            comp.eigenvec1_x * vec1_scale,
            -comp.eigenvec1_y * vec1_scale);

    // Второй собственный вектор
    fprintf(gnuplotPipe, "%f %f %f %f\n", 
            comp.center_x,
            cy,
            comp.eigenvec2_x * vec2_scale,
            -comp.eigenvec2_y * vec2_scale);
}
fprintf(gnuplotPipe, "e\n");

        pclose(gnuplotPipe);
        logPlotEnd("BinaryWithComponents");
    }
      void gnuplot(std::unique_ptr<Pole>& p, const std::string& filename) {
        logPlotStart("HeightMap3D", filename);
        
        if (p == nullptr) {
            logger.error("[GnuplotInterface::gnuplot] Ошибка: данные высот не инициализированы");
            return;
        }

        int rows = p->field.size();
        int cols = p->field[0].size();
        logger.debug(std::string("Размер сетки: ") + std::to_string(cols) + "x" + std::to_string(rows));
        
        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::gnuplot] Ошибка открытия gnuplot pipe");
            return;
        }
    
    // Настройки 3D графика
    fprintf(gnuplotPipe, "set terminal pngcairo enhanced size 1600,1200\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "set zlabel 'Height'\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", cols - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", rows - 1);
    fprintf(gnuplotPipe, "set zrange [*:*]\n"); // Автомасштабирование по Z
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set view 60, 30, 1, 1\n"); // Угол обзора
    
    // Формат данных: x y z
    fprintf(gnuplotPipe, "splot '-' with pm3d title 'Height Map'\n");
    
    // Записываем данные в правильном порядке
    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            // Преобразуем координаты для правильной ориентации
            fprintf(gnuplotPipe, "%d %d %f\n", x, rows-1-y, p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n"); // Пустая строка между слоями Y
    }
    
    fprintf(gnuplotPipe, "e\n"); // Конец данных
    pclose(gnuplotPipe);
    logPlotEnd("HeightMap3D");
}

void plotVoronoi(const std::unique_ptr<Pole>& p, 
                    const std::vector<VoronoiEdge>& edges, 
                    const std::vector<PointD>& sites, 
                    const std::string& filename) {
        logPlotStart("VoronoiDiagram", filename);
        
        if (!p || edges.empty()) {
            logger.warning("[GnuplotInterface::plotVoronoi] Нет данных для визуализации");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("Диаграмма Вороного: ") + std::to_string(edges.size()) + " ребер, " + 
                   std::to_string(sites.size()) + " сайтов");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotVoronoi] Ошибка открытия gnuplot pipe");
            return;
        }

    // Настройки графика с контрастными цветами для красного фона
    fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set title 'Voronoi Diagram on Red Field'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width - 1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height - 1);
    fprintf(gnuplotPipe, "unset key\n");
    
    // Цветовая схема:
    fprintf(gnuplotPipe, "set style line 1 lc rgb '#00FF00' lw 2    # Ярко-зеленые ребра\n");
    fprintf(gnuplotPipe, "set style line 2 lc rgb '#FFFFFF' pt 7 ps 2 # Белые центры\n");

    // Многослойный график:
    fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
    fprintf(gnuplotPipe, "'-' with lines ls 1, \\\n");  // Ребра
    fprintf(gnuplotPipe, "'-' with points ls 2\n");     // Центры

    // 1. Данные поля (красный фон)
    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 2. Ребра Вороного (ярко-зеленые)
    for (const auto& edge : edges) {
        fprintf(gnuplotPipe, "%f %f\n%f %f\n\n", 
                edge.start.x, transformY(edge.start.y, height),
                edge.end.x, transformY(edge.end.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    // 3. Центры (белые точки)
    for (const auto& site : sites) {
        fprintf(gnuplotPipe, "%f %f\n", 
                site.x, 
                transformY(site.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("VoronoiDiagram");
}
    
    void plotDelaunay(const std::vector<Triangle>& triangles, 
                     std::unique_ptr<Pole>& p, 
                     const std::string& filename) {
        logPlotStart("DelaunayTriangulation", filename);
        
        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotDelaunay] Ошибка открытия gnuplot pipe");
            return;
        }
        
        if (!p) {
            logger.error("[GnuplotInterface::plotDelaunay] Нет данных высот");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("Триангуляция Делоне: ") + std::to_string(triangles.size()) + " треугольников");

    // Улучшенные настройки графика
    fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200 enhanced font 'Arial,12'\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set title 'Delaunay Triangulation'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
    fprintf(gnuplotPipe, "unset key\n");
    
    // Цветовая схема
    fprintf(gnuplotPipe, "set style line 1 lc rgb '#00FF00' lw 1.5\n");   // Ярко-зеленые линии
    
    // Многослойный график: фон + треугольники
    fprintf(gnuplotPipe, "plot '-' matrix with image, '-' with lines ls 1\n");

    // 1. Данные фона (с инверсией Y)
    for (int y = height-1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 2. Данные треугольников (с инверсией Y)
    for (const auto& tri : triangles) {
        fprintf(gnuplotPipe, "%f %f\n", tri.a.x, transformY(tri.a.y, height));
        fprintf(gnuplotPipe, "%f %f\n", tri.b.x, transformY(tri.b.y, height));
        fprintf(gnuplotPipe, "%f %f\n", tri.c.x, transformY(tri.c.y, height));
        fprintf(gnuplotPipe, "%f %f\n\n", tri.a.x, transformY(tri.a.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("DelaunayTriangulation");
}

void plotPath(const std::vector<PointD>& path, 
                 const std::unique_ptr<Pole>& p, 
                 const std::string& filename, 
                 const DispatcherParams& params, 
                 const PathFinder& pathFinder, 
                 const double Radius) {
        logPlotStart("PathVisualization", filename);
        
        if (!p || path.empty()) {
            logger.warning("[GnuplotInterface::plotPath] Нет данных пути для визуализации");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("Визуализация пути: ") + std::to_string(path.size()) + " точек");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotPath] Ошибка открытия gnuplot pipe");
            return;
        }

    // Настройки графика
    fprintf(gnuplotPipe, "set terminal pngcairo size 1600,1200\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set title 'Path Visualization'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
    fprintf(gnuplotPipe, "unset key\n");

    // 1. Собираем все пиксели пути
    std::vector<PointD> pathPixels;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto segment = pathFinder.bresenhamLine(path[i], path[i+1]);
        pathPixels.insert(pathPixels.end(), segment.begin(), segment.end());
    }

    // 2. Создаем команды для окружностей
    for (const auto& pixel : pathPixels) {
        const int x = static_cast<int>(pixel.x);
        const int y = static_cast<int>(pixel.y);
        const double currentHeight = p->field[y][x];
        const double heightDiff = std::abs(currentHeight - params.slice);
        
        if (heightDiff < Radius) {
            const double radius = std::sqrt(
                Radius * Radius - 
                heightDiff * heightDiff
            );
            
            fprintf(gnuplotPipe, "set object circle at %d,%d size %f fc rgb '#004400' fs transparent solid 0.5 front\n",
                   static_cast<int>(pixel.x), static_cast<int>(transformY(pixel.y, height)), radius);
        }
    }
    
    // 3. Подписи START/END
    fprintf(gnuplotPipe, "set label 'START' at %d,%d front\n", 
            static_cast<int>(params.pointA_x), static_cast<int>(transformY(params.pointA_y, height)));
    fprintf(gnuplotPipe, "set label 'END' at %d,%d front\n",
            static_cast<int>(params.pointB_x), static_cast<int>(transformY(params.pointB_y, height)));
    
    // Многослойный график: фон + путь + точки
    fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
    fprintf(gnuplotPipe, "'-' with lines lw 2 lc rgb '#00FF0080', \\\n"); // Светло-зеленый 
    fprintf(gnuplotPipe, "'-' with points pt 7 ps 2 lc 'blue', \\\n");
    fprintf(gnuplotPipe, "'-' with points pt 9 ps 2 lc 'purple'\n");

    // 4. Данные фона (инверсия Y)
    for (int y = height-1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    // 5. Путь (инверсия Y)
    for (const auto& point : path) {
        fprintf(gnuplotPipe, "%f %f\n", 
                point.x, 
                transformY(point.y, height));
    }
    fprintf(gnuplotPipe, "e\n");

    // 6. Точка A (синяя)
    fprintf(gnuplotPipe, "%f %f\n", 
            params.pointA_x, 
            transformY(params.pointA_y, height));
    fprintf(gnuplotPipe, "e\n");

    // 7. Точка B (фиолетовая)
    fprintf(gnuplotPipe, "%f %f\n", 
            params.pointB_x, 
            transformY(params.pointB_y, height));
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("PathVisualization");
}

void plotInteractive3DPath(const std::vector<PointD>& path, 
                              const std::unique_ptr<Pole>& p, 
                              const PointD& start, 
                              const PointD& end, 
                              const PathFinder& pathFinder,  
                              const double sphereRadius) {
        logPlotStart("Interactive3DPath", "интерактивное окно");
        
        if (!p || path.empty()) {
            logger.warning("[GnuplotInterface::plotInteractive3DPath] Нет данных пути");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("Интерактивный 3D путь: ") + std::to_string(path.size()) + " точек");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plotInteractive3DPath] Ошибка открытия gnuplot pipe");
            return;
        }
    
    // 1. Настройки графика (должны идти первыми!)
    fprintf(gnuplotPipe, "set terminal qt size 1600,1200 enhanced\n");
    fprintf(gnuplotPipe, "set title '3D Path Visualization (Use mouse to rotate)'\n");
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "set zlabel 'Height'\n");
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set view 60, 30\n");
    fprintf(gnuplotPipe, "set mouse\n");  // Включаем интерактивное управление
    fprintf(gnuplotPipe, "set key outside\n");
    
    // 2. Создаем проекцию пути (все пиксели между узлами)
    std::vector<std::pair<int, int>> pathProjection;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto linePixels = pathFinder.bresenhamLine(path[i], path[i+1]); // Алгоритм Брезенхема для рисования линии
        for (const auto& pixel : linePixels) {
            pathProjection.emplace_back(static_cast<int>(pixel.x), static_cast<int>(pixel.y));
        }
    }
    
    // 3. Проекция пути на поле (зеленые точки)
    for (const auto& [x, y] : pathProjection) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#00FF00' fs solid front\n",
                x, static_cast<int>(transformY(y, height)), p->field[y][x], sphereRadius);
    }

    // 4. Точки A и B
    int startX = static_cast<int>(start.x);
    int startY = static_cast<int>(start.y);
    int endX = static_cast<int>(end.x);
    int endY = static_cast<int>(end.y);
    
    if (startX >= 0 && startY >= 0 && startX < width && startY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#FF0000' fs solid front\n",
                startX, static_cast<int>(transformY(startY, height)), 
                p->field[startY][startX], sphereRadius);
    }
    
    fprintf(gnuplotPipe, "set label 'START' at %d,%d,%f front\n",
                startX, static_cast<int>(transformY(startY, height)), p->field[startY][startX] + 1);
    
    if (endX >= 0 && endY >= 0 && endX < width && endY < height) {
         fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#0000FF' fs solid front\n",
                endX, static_cast<int>(transformY(endY, height)), 
                p->field[endY][endX], sphereRadius);
    }
    
     fprintf(gnuplotPipe, "set label 'END' at %d,%d,%f front\n",
                endX, static_cast<int>(transformY(endY, height)), p->field[endY][endX] + 1);

    // 5. Рисуем поверхность
    fprintf(gnuplotPipe, "splot '-' with pm3d title 'Terrain'\n");
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%d %d %f\n", x, height-1-y, p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");
    fprintf(gnuplotPipe, "pause mouse close\n");
    pclose(gnuplotPipe);
    logPlotEnd("Interactive3DPath");
}

void plot3DPath(const std::vector<PointD>& path, 
                   const std::unique_ptr<Pole>& p, 
                   const std::string& filename, 
                   const PointD& start, 
                   const PointD& end, 
                   const PathFinder& pathFinder, 
                   const double sphereRadius) {
        logPlotStart("3DPathProjection", filename);
        
        if (!p || path.empty()) {
            logger.warning("[GnuplotInterface::plot3DPath] Нет данных пути");
            return;
        }

        const int height = p->field.size();
        const int width = p->field[0].size();
        logger.debug(std::string("3D проекция пути: ") + std::to_string(path.size()) + " точек");

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            logger.error("[GnuplotInterface::plot3DPath] Ошибка открытия gnuplot pipe");
            return;
        }

    // 1. Создаем проекцию пути (все пиксели между узлами)
    std::vector<std::pair<int, int>> pathProjection;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto linePixels = pathFinder.bresenhamLine(path[i], path[i+1]); // Алгоритм Брезенхема для рисования линии
        for (const auto& pixel : linePixels) {
            pathProjection.emplace_back(static_cast<int>(pixel.x), static_cast<int>(pixel.y));
        }
    }

    // 2. Настройки графика
    fprintf(gnuplotPipe, "set terminal pngcairo enhanced size 1600,1200\n");
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set title '3D Path Projection'\n");
    fprintf(gnuplotPipe, "set pm3d\n");
    fprintf(gnuplotPipe, "set hidden3d\n");
    fprintf(gnuplotPipe, "set view 60, 30\n");
    
    // 3. Проекция пути на поле (зеленые точки)
    for (const auto& [x, y] : pathProjection) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#00FF00' fs solid front\n",
                x, static_cast<int>(transformY(y, height)), p->field[y][x], sphereRadius);
    }

    // 4. Точки A и B
    int startX = static_cast<int>(start.x);
    int startY = static_cast<int>(start.y);
    int endX = static_cast<int>(end.x);
    int endY = static_cast<int>(end.y);
    
    if (startX >= 0 && startY >= 0 && startX < width && startY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#FF0000' fs solid front\n",
                startX, static_cast<int>(transformY(startY, height)), 
                p->field[startY][startX] , sphereRadius);
    }
    
    fprintf(gnuplotPipe, "set label 'START' at %d,%d,%f front\n",
                startX, static_cast<int>(transformY(startY, height)), p->field[startY][startX] + 1);
    
    if (endX >= 0 && endY >= 0 && endX < width && endY < height) {
         fprintf(gnuplotPipe, "set object circle at %d,%d,%f size %f fc rgb '#0000FF' fs solid front\n",
                endX, static_cast<int>(transformY(endY, height)), 
                p->field[endY][endX], sphereRadius);
    }
    
     fprintf(gnuplotPipe, "set label 'END' at %d,%d,%f front\n",
                endX, static_cast<int>(transformY(endY, height)), p->field[endY][endX] + 1);

    // 5. Рисуем поверхность
    fprintf(gnuplotPipe, "splot '-' with pm3d title 'Terrain'\n");
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%d %d %f\n", x, height-1-y, p->field[y][x]);
        }
        fprintf(gnuplotPipe, "\n");
    }
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
    logPlotEnd("3DPathProjection");
}

};
   
class ComponentCalculator {
private:
    Logger& logger;

    void logTriangle(const Triangle& tri, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Треугольник ["
            << "A(" << tri.a.x << "," << tri.a.y << "), "
            << "B(" << tri.b.x << "," << tri.b.y << "), "
            << "C(" << tri.c.x << "," << tri.c.y << ")]";
        logger.debug(oss.str());
    }

    void logPoint(const PointD& p, const std::string& prefix = "") const {
        logger.trace(prefix + "Точка (" + std::to_string(p.x) + "," + std::to_string(p.y) + ")");
    }

    void logEdge(const Edge& edge, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Ребро: (" << edge.a.x << "," << edge.a.y << ")-(" 
            << edge.b.x << "," << edge.b.y << ")";
        logger.trace(oss.str());
    }

public:
    ComponentCalculator(Logger& lg) : logger(lg) {
        logger.trace("[ComponentCalculator] Инициализация калькулятора компонент");
    }

    std::vector<Triangle> bowyerWatson(const std::vector<PointD>& points) {
        logger.info("[ComponentCalculator::bowyerWatson] Начало триангуляции Боуера-Ватсона");
        logger.debug(std::string("Количество точек: ") + std::to_string(points.size()));
        
        std::vector<Triangle> triangles;

        // Проверка на минимальное количество точек
        if (points.size() < 3) {
            logger.error("[ComponentCalculator::bowyerWatson] Недостаточно точек для триангуляции (меньше 3)");
            return triangles;
        }

        // Проверка на NaN-точки
        for (const auto& p : points) {
            if (std::isnan(p.x) || std::isnan(p.y)) {
                logger.error("[ComponentCalculator::bowyerWatson] Обнаружена некорректная точка (NaN)");
                return triangles;
            }
        }

        // Проверка на коллинеарность
        if (std::all_of(points.begin(), points.end(), 
                       [&](const PointD& p) { return Triangle::areCollinear(points[0], points[1], p); })) {
            logger.error("[ComponentCalculator::bowyerWatson] Все точки коллинеарны, триангуляция невозможна");
            return triangles;
        }

        // Создаём супер-треугольник
        auto [minX, maxX] = std::minmax_element(points.begin(), points.end(), 
                                              [](auto& a, auto& b) { return a.x < b.x; });
        auto [minY, maxY] = std::minmax_element(points.begin(), points.end(), 
                                              [](auto& a, auto& b) { return a.y < b.y; });

        double dx = (maxX->x - minX->x) * 10;
        double dy = (maxY->y - minY->y) * 10;
        PointD p1(minX->x - dx, minY->y - dy);
        PointD p2(maxX->x + dx, minY->y - dy);
        PointD p3((minX->x + maxX->x) / 2, maxY->y + dy);

        triangles.emplace_back(p1, p2, p3);
        logger.info("[ComponentCalculator::bowyerWatson] Создан супер-треугольник");
        logTriangle(triangles.back(), "Супер-треугольник: ");

        // Основной алгоритм
        size_t pointIndex = 0;
        for (const auto& point : points) {
            logger.trace(std::string("[ComponentCalculator::bowyerWatson] Обработка точки #") + 
                        std::to_string(++pointIndex) + " (" + 
                        std::to_string(point.x) + "," + std::to_string(point.y) + ")");

            std::vector<Triangle> badTriangles;
            std::copy_if(triangles.begin(), triangles.end(), std::back_inserter(badTriangles),
                         [&](const Triangle& tri) { return isPointInCircumcircle(point, tri); });

            logger.debug(std::string("[ComponentCalculator::bowyerWatson] Найдено ") + 
                       std::to_string(badTriangles.size()) + " плохих треугольников");

            std::vector<Edge> polygonEdges;
            for (const auto& tri : badTriangles) {
                for (const auto& edge : { Edge{tri.a, tri.b}, Edge{tri.b, tri.c}, Edge{tri.c, tri.a} }) {
                    bool isShared = std::any_of(badTriangles.begin(), badTriangles.end(),
                        [&](const Triangle& other) { 
                            return !(tri == other) && hasEdge(other, edge); 
                        });
                    
                    if (!isShared) {
                        polygonEdges.push_back(edge);
                        logEdge(edge, "Добавлено ребро в полигон: ");
                    }
                }
            }

            // Удаляем плохие треугольники
            triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
                [&](const Triangle& t) { 
                    return std::find(badTriangles.begin(), badTriangles.end(), t) != badTriangles.end();
                }), triangles.end());

            // Добавляем новые треугольники
            for (const auto& edge : polygonEdges) {
                triangles.emplace_back(edge.a, edge.b, point);
                logTriangle(triangles.back(), "Добавлен новый треугольник: ");
            }
        }

        // Удаляем треугольники, связанные с супер-треугольником
        size_t initialCount = triangles.size();
        triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
            [&](const Triangle& t) {
                return t.a == p1 || t.a == p2 || t.a == p3 ||
                       t.b == p1 || t.b == p2 || t.b == p3 ||
                       t.c == p1 || t.c == p2 || t.c == p3;
            }), triangles.end());

        logger.info("[ComponentCalculator::bowyerWatson] Триангуляция завершена");
        logger.debug(std::string("Удалено ") + std::to_string(initialCount - triangles.size()) + 
                    " треугольников, связанных с супер-треугольником");
        logger.debug(std::string("Итоговое количество треугольников: ") + std::to_string(triangles.size()));

        return triangles;
    }

    bool isPointInCircumcircle(const PointD& p, const Triangle& tri) {
        double ax = tri.a.x - p.x, ay = tri.a.y - p.y;
        double bx = tri.b.x - p.x, by = tri.b.y - p.y;
        double cx = tri.c.x - p.x, cy = tri.c.y - p.y;
        
        double det = ax * (by * (cx*cx + cy*cy) - cy * (bx*bx + by*by)) 
                   - ay * (bx * (cx*cx + cy*cy) - cx * (bx*bx + by*by)) 
                   + (ax*ax + ay*ay) * (bx*cy - by*cx);

        bool inside = det > 0;
        logger.trace(std::string("[ComponentCalculator::isPointInCircumcircle] Точка (") + 
                    std::to_string(p.x) + "," + std::to_string(p.y) + ") " +
                    (inside ? "внутри" : "вне") + " окружности треугольника");
        logTriangle(tri, "Проверяемый треугольник: ");

        return inside;
    }

    bool hasEdge(const Triangle& tri, const Edge& edge) {
        bool has = Edge(tri.a, tri.b) == edge 
                || Edge(tri.b, tri.c) == edge 
                || Edge(tri.c, tri.a) == edge;

        logger.trace(std::string("[ComponentCalculator::hasEdge] Треугольник ") + 
                    std::string(has ? "содержит" : "не содержит") + " ребро");
        logTriangle(tri);
        logEdge(edge);

        return has;
    }
       
    int incrementAndCollect(std::vector<std::vector<double>>& componenta, 
                          std::vector<std::vector<double>>& CopyPole, 
                          int x, int y, int i, int& pixelCount) {
        logger.trace(std::string("[ComponentCalculator::incrementAndCollect] Проверка пикселя (") + 
                    std::to_string(x) + "," + std::to_string(y) + "), глубина=" + 
                    std::to_string(i));

        if (x < 1 || y < 1 || x > (int)componenta[0].size() - 2 || 
            y > (int)componenta.size() - 2 || CopyPole[y][x] < 250) {
            logger.trace("[ComponentCalculator::incrementAndCollect] Выход за границы или неподходящее значение");
            return -1;
        }

        if (CopyPole[y][x] >= 255 && CopyPole[y][x] <= 255) {
            CopyPole[y][x] = 0;
            pixelCount++;
            componenta[y][x] = 255;
            
            logger.trace(std::string("[ComponentCalculator::incrementAndCollect] Обработан пиксель (") + 
                        std::to_string(x) + "," + std::to_string(y) + "), счетчик=" + 
                        std::to_string(pixelCount));

            // Рекурсивные вызовы
            incrementAndCollect(componenta, CopyPole, x + 1, y, i + 1, pixelCount);
            incrementAndCollect(componenta, CopyPole, x - 1, y, i + 1, pixelCount);
            incrementAndCollect(componenta, CopyPole, x, y + 1, i + 1, pixelCount);
            incrementAndCollect(componenta, CopyPole, x, y - 1, i + 1, pixelCount);
        }
        
        return pixelCount;
    }

    void bin(std::vector<std::vector<double>>& CopyPole, 
            int slice, 
            std::unique_ptr<Pole>& p, 
            ThresholdMode mode) {
        logger.info(std::string("[ComponentCalculator::bin] Бинаризация данных, slice=") + 
                  std::to_string(slice) + ", mode=" + 
                  (mode == ThresholdMode::Peaks ? "Peaks" : 
                   mode == ThresholdMode::Valleys ? "Valleys" : "All"));

        if (p == nullptr) {
            logger.error("[ComponentCalculator::bin] Ошибка: данные высот не инициализированы!");
            return;
        }
        
        CopyPole = p->field;
        int symmetric_slice = 2*127 - slice;

        for (int x = 0; x < (int)p->field[0].size(); ++x) {
            for (int y = 0; y < (int)p->field.size(); ++y) {
                double value = p->field[y][x];
                switch (mode) {
                    case ThresholdMode::All: {
                        bool is_peak = (slice >= 127) ? (value > slice) : (value > symmetric_slice);
                        bool is_valley = (slice >= 127) ? (value < symmetric_slice) : (value < slice);
                        CopyPole[y][x] = (is_peak || is_valley) ? 255 : 0;
                        break;
                    }
                    case ThresholdMode::Peaks:
                        CopyPole[y][x] = (slice >= 127) ? (value > slice) : (value > symmetric_slice);
                        break;
                    case ThresholdMode::Valleys:
                        CopyPole[y][x] = (slice >= 127) ? (value < symmetric_slice) : (value < slice);
                        break;
                }
            } 
        }

        logger.info("[ComponentCalculator::bin] Бинаризация завершена");
    }
    
    void wave(int noisy, 
             std::vector<Component>& componenti, 
             std::vector<std::vector<double>>& CopyPole, 
             std::unique_ptr<Pole>& p) {
        logger.info(std::string("[ComponentCalculator::wave] Начало волнового алгоритма, порог=") + 
                  std::to_string(noisy));

        if (p == nullptr) {
            logger.error("[ComponentCalculator::wave] Ошибка: данные высот не инициализированы!");
            return;
        }

        int rows = p->field.size();
        int cols = (rows > 0) ? p->field[0].size() : 0;
        std::vector<Component> noiseComponents;

        logger.debug(std::string("Размер данных: ") + std::to_string(cols) + "x" + std::to_string(rows));

        for (int y = 2; y < rows - 2; ++y) {
            for (int x = 2; x < cols - 2; ++x) {
                if (CopyPole[y][x] <= 255 && CopyPole[y][x] >= 255) {
                    std::vector<std::vector<double>> componentData(rows, std::vector<double>(cols, 0));
                    int pixelCount = 0;
                    incrementAndCollect(componentData, CopyPole, x, y, 0, pixelCount);

                    if (pixelCount >= noisy) {
                        Component component(logger, componentData, pixelCount);
                        componenti.push_back(component);
                        
                        logger.debug(std::string("[ComponentCalculator::wave] Значимая компонента: ") +
                                   "pixels=" + std::to_string(pixelCount) +
                                   ", center=(" + std::to_string(component.center_x) + "," + 
                                   std::to_string(component.center_y) + ")" +
                                   ", size=(" + std::to_string(component.max_x - component.min_x) + 
                                   "x" + std::to_string(component.max_y - component.min_y) + ")");
                    } else {
                        Component noiseComponent(logger, componentData, pixelCount);
                        noiseComponents.push_back(noiseComponent);
                        
                        logger.trace(std::string("[ComponentCalculator::wave] Шумовая компонента: ") +
                                    "pixels=" + std::to_string(pixelCount) +
                                    ", center=(" + std::to_string(noiseComponent.center_x) + "," + 
                                    std::to_string(noiseComponent.center_y) + ")");

                        // Удаляем шум из основного поля
                        for (int i = 0; i < rows; ++i) {
                            for (int j = 0; j < cols; ++j) {
                                if (componentData[i][j] == 255) {
                                    p->field[i][j] = 127;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        logger.info("[ComponentCalculator::wave] Обработка завершена");
        logger.debug(std::string("Найдено значимых компонент: ") + std::to_string(componenti.size()) +
                   ", шумовых компонент: " + std::to_string(noiseComponents.size()));
    }
};

class KMeans {
private:
    Logger& logger;

    void logCenters(const std::vector<std::vector<double>>& centers, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Центры кластеров (" << centers.size() << "):";
        for (size_t i = 0; i < centers.size(); ++i) {
            oss << "\n  Кластер " << i << ": (" << centers[i][0] << ", " << centers[i][1] << ")";
        }
        logger.debug(oss.str());
    }

    void logPoint(const std::vector<double>& point, const std::string& prefix = "") const {
        logger.trace(prefix + "Точка (" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ")");
    }

public:
    std::vector<int> labels;
    std::vector<std::vector<double>> centers; // Центры кластеров

    struct ClusterResult {
        std::vector<int> labels;
        std::vector<std::array<int, 3>> colors;
        std::vector<std::vector<double>> centers;
    };

    KMeans(Logger& lg) : logger(lg) {
        logger.trace("[KMeans] Инициализация алгоритма K-средних");
    }

    void initializeCenters(const std::vector<std::vector<double>>& data, int k) {
        logger.info("[KMeans::initializeCenters] Инициализация центров кластеров");
        logger.debug(std::string("Количество кластеров: ") + std::to_string(k) + 
                   ", количество точек: " + std::to_string(data.size()));

        if (data.empty() || k <= 0) {
            logger.error("[KMeans::initializeCenters] Ошибка: пустые данные или k <= 0");
            return;
        }
        
        centers.clear();
        std::vector<size_t> indices(data.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});
        
        for (int i = 0; i < k && i < static_cast<int>(data.size()); ++i) {
            centers.push_back(data[indices[i]]);
            logPoint(data[indices[i]], "Выбран начальный центр: ");
        }

        logCenters(centers, "Инициализированы центры: ");
    }

    double distance(const std::vector<double>& a, const std::vector<double>& b) const {
        double dist = std::hypot(a[0]-b[0], a[1]-b[1]);
        logger.trace(std::string("[KMeans::distance] Расстояние между (") + 
                   std::to_string(a[0]) + "," + std::to_string(a[1]) + ") и (" +
                   std::to_string(b[0]) + "," + std::to_string(b[1]) + ") = " +
                   std::to_string(dist));
        return dist;
    }

    ClusterResult cluster(const std::vector<std::vector<double>>& data, int k) {
        logger.info("[KMeans::cluster] Начало кластеризации");
        logger.debug(std::string("Количество кластеров: ") + std::to_string(k) + 
                   ", точек: " + std::to_string(data.size()));

        ClusterResult result;
        initializeCenters(data, k);
        bool changed;
        int iteration = 0;
        
        do {
            iteration++;
            logger.debug(std::string("[KMeans::cluster] Итерация ") + std::to_string(iteration));
            changed = false;
            labels.assign(data.size(), -1);

            // Назначение меток
            for (size_t i = 0; i < data.size(); ++i) {
                double minDist = std::numeric_limits<double>::max();
                for (size_t j = 0; j < centers.size(); ++j) {
                    double dist = distance(data[i], centers[j]);
                    if (dist < minDist) {
                        minDist = dist;
                        labels[i] = j;
                    }
                }
                logger.trace(std::string("[KMeans::cluster] Точка ") + std::to_string(i) + 
                            " отнесена к кластеру " + std::to_string(labels[i]));
            }

            // Обновление центров
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

            // Проверка изменений
            for (size_t j = 0; j < centers.size(); ++j) {
                if (counts[j] > 0) {
                    newCenters[j][0] /= counts[j];
                    newCenters[j][1] /= counts[j];
                    double centerDist = distance(newCenters[j], centers[j]);
                    if (centerDist > 1e-6) {
                        changed = true;
                    }
                    logger.debug(std::string("[KMeans::cluster] Кластер ") + std::to_string(j) + 
                               ": точек=" + std::to_string(counts[j]) + 
                               ", смещение центра=" + std::to_string(centerDist));
                } else {
                    logger.warning(std::string("[KMeans::cluster] Кластер ") + std::to_string(j) + 
                                 " остался без точек!");
                }
            }
            
            centers = newCenters;
            logCenters(centers, "Обновленные центры: ");

        } while (changed);

        logger.info(std::string("[KMeans::cluster] Кластеризация завершена за ") + 
                   std::to_string(iteration) + " итераций");

        result.labels = labels;
        result.centers = centers;
        result.colors = ColorGenerator::generateColors(k, logger);
        
        logger.debug(std::string("[KMeans::cluster] Сгенерировано ") + 
                    std::to_string(result.colors.size()) + " цветов для кластеров");
        logCenters(result.centers, "Финальные центры кластеров: ");
        
        return result;
    }

    ClusterResult kmeansWithKernels(const std::vector<std::vector<double>>& data, int k, int kernelSize) {
        logger.info("[KMeans::kmeansWithKernels] Начало кластеризации с ядрами");
        logger.debug(std::string("Параметры: k=") + std::to_string(k) + 
                   ", kernelSize=" + std::to_string(kernelSize) + 
                   ", точек=" + std::to_string(data.size()));

        ClusterResult result;
        if (data.empty() || k <= 0 || kernelSize <= 0) {
            logger.error("[KMeans::kmeansWithKernels] Ошибка: пустые данные или неверные параметры");
            return result;
        }

        // Шаг 1: Выборка ядерных центров
        std::vector<std::vector<double>> kernelCenters;
        std::vector<size_t> indices(data.size());
        std::iota(indices.begin(), indices.end(), 0);
        
        size_t samplesToTake = std::min(k * kernelSize, static_cast<int>(data.size()));
        if (samplesToTake == 0) {
            logger.error("[KMeans::kmeansWithKernels] Ошибка: недостаточно точек для выборки");
            return result;
        }

        std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});
        for (size_t i = 0; i < samplesToTake; ++i) {
            kernelCenters.push_back(data[indices[i]]);
        }

        logger.debug(std::string("[KMeans::kmeansWithKernels] Выбрано ") + 
                    std::to_string(kernelCenters.size()) + " ядерных точек");

        // Шаг 2: Первичная кластеризация
        initializeCenters(kernelCenters, k * kernelSize);
        auto kernelResult = cluster(kernelCenters, k * kernelSize);

        // Шаг 3: Агрегация центров
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

        // Шаг 4: Нормализация или замена пустых кластеров
        for (int i = 0; i < k; ++i) {
            if (counts[i] > 0) {
                aggregatedCenters[i][0] /= counts[i];
                aggregatedCenters[i][1] /= counts[i];
                logger.debug(std::string("[KMeans::kmeansWithKernels] Агрегированный центр ") + 
                            std::to_string(i) + " на основе " + 
                            std::to_string(counts[i]) + " точек");
            } else {
                aggregatedCenters[i] = data[std::rand() % data.size()];
                logger.warning(std::string("[KMeans::kmeansWithKernels] Кластер ") + 
                             std::to_string(i) + " был пуст, заменен случайной точкой");
            }
        }

        // Финальная кластеризация
        centers = aggregatedCenters;
        result = cluster(data, k);
        if (result.centers.size() != k) {
            logger.error("[KMeans::kmeansWithKernels] Несоответствие количества кластеров!");
        }
        result.colors = ColorGenerator::generateColors(k, logger);

        logger.info("[KMeans::kmeansWithKernels] Кластеризация с ядрами завершена");
        logCenters(result.centers, "Финальные центры: ");
        
        return result;
    }
};

class Control {
private:
    void logOperation(LogLevel level, const std::string& operation, const std::string& details = "") {
        std::string message = std::string("Operation: ") + operation;
        if (!details.empty()) {
            message += std::string(", ") + details;
        }
        logger.logMessage(level, message);
    }

public:
    Config& config;
    Logger& logger;
    
    // Компоненты системы
    Copier copier;
    std::vector<std::array<int, 3>> colors;
    std::vector<std::vector<double>> CopyPole;
    std::vector<Gaus> gaussi;
    std::vector<Component> componenti;
    GaussBuilder gaussBuilder;
    BmpHandler bmpHandler;
    GnuplotInterface gnuplotInterface;
    ComponentCalculator componentCalculator;
    TerrainGrid terrainGrid;
    std::unique_ptr<Pole> p = nullptr;
    std::unique_ptr<KMeans> kMeans = nullptr;
    std::vector<std::vector<double>> kMeansData;
    std::vector<Triangle> lastTriangulation;
    std::vector<VoronoiEdge> voronoiEdges;
    PathFinder pathFinder;
    VoronoiDiagram voronoi;
    std::vector<std::vector<PointD>> paths;
    std::vector<PointD> clusterCenters;
    std::vector<PointD> path;

    Control(Config& cfg, Logger& log) 
        : config(cfg), 
          logger(log),
          copier(log),
          gaussBuilder(log),
          bmpHandler(log),
          gnuplotInterface(log),
          componentCalculator(log),
          terrainGrid(cfg.fieldWidth, cfg.fieldHeight, log),
          pathFinder(cfg, log),
          voronoi(log) {
        
        kMeans = std::make_unique<KMeans>(log);
        
        logger.info("Control system initialized");
        logger.debug(std::string("Field dimensions: ") + std::to_string(cfg.fieldWidth) + "x" + 
                    std::to_string(cfg.fieldHeight));
    }

    std::vector<PointD> getClusterCenters() const {
        std::vector<PointD> centers;
        for (const auto& component : componenti) {
            if (std::isnan(component.center_x) || std::isnan(component.center_y)) {
                logger.logMessage(LogLevel::Warning, "Skipping invalid cluster center (NaN)");
                continue;
            }
            if (component.center_x >= 0 && component.center_x < config.fieldWidth &&
                component.center_y >= 0 && component.center_y < config.fieldHeight) {
                centers.emplace_back(component.center_x, component.center_y);
            } else {
                logger.logMessage(LogLevel::Error, 
                    std::string("Invalid cluster center: (") + std::to_string(component.center_x) + 
                    ", " + std::to_string(component.center_y) + ")");
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
        size_t validPoints = 0;
        
        for (size_t y = 0; y < copyPole.size(); ++y) {
            for (size_t x = 0; x < copyPole[0].size(); ++x) {
                if (copyPole[y][x] > 0) { 
                    kMeansData.push_back({static_cast<double>(x), static_cast<double>(y)});
                    validPoints++;
                }
            }
        }
            logger.logMessage(LogLevel::Debug, 
                std::string("Prepared ") + std::to_string(validPoints) + " points for K-means");
    }
    
    void Dispetcher(DispatcherParams& params) {     
            logger.logMessage(LogLevel::Info, std::string("Processing command: ") + params.s);
               
        if (params.s == "init") {
            gaussBuilder.init(params.A, params.B, p);
            logOperation(LogLevel::Info, std::string("init"), std::string("size: ") + std::to_string(params.A) + "x" + std::to_string(params.B));
        }

        if (params.s == "g") {
            gaussBuilder.addgauss(params.h, params.x, params.y, params.sx, params.sy, gaussi);
            logOperation(LogLevel::Info, std::string("addgauss"), 
                std::string("x=") + std::to_string(params.x) + 
                ", y=" + std::to_string(params.y) + 
                ", h=" + std::to_string(params.h));
        }
        
        if (params.s == "generate") {
            gaussBuilder.generate(p, gaussi);
            logOperation(LogLevel::Info, std::string("generate"));
        }
        
        if (params.s == "gnuplot") {
            gnuplotInterface.gnuplot(p, params.filename);
            logOperation(LogLevel::Info, std::string("gnuplot"), std::string("file: ") + params.filename);
        }

        if (params.s == "PlotMetedata") {
            gnuplotInterface.plotBinaryWithComponents(CopyPole, componenti, params.filename);
            logOperation(LogLevel::Info, std::string("PlotMetedata"), std::string("file: ") + params.filename);
        }
        
        if (params.s == "PlotVoronoi") {
            gnuplotInterface.plotVoronoi(p, voronoiEdges, clusterCenters, params.filename);
            logOperation(LogLevel::Info, std::string("PlotVoronoi"), std::string("file: ") + params.filename);
        }
        
        if (params.s == "PlotDelaunay") {
            gnuplotInterface.plotDelaunay(lastTriangulation, p, params.filename);
            logOperation(LogLevel::Info, std::string("PlotDelaunay"), std::string("file: ") + params.filename);
        }
        
        if (params.s == "PlotPath") {
            gnuplotInterface.plotPath(path, p, params.filename, params, pathFinder, config.vehicleRadius);
            logOperation(LogLevel::Info, std::string("PlotPath"), std::string("file: ") + params.filename);
        }

        if (params.s == "bmp_write") {
            if (params.bmp_mode == BmpWriteMode::Full) {
                bmpHandler.bmp_write(p->field, params.filename);
            } else {
                bmpHandler.bmp_write(CopyPole, params.filename);
            }
            logOperation(LogLevel::Info, std::string("bmp_write"), 
                std::string("mode: ") + std::string(params.bmp_mode == BmpWriteMode::Full ? "full" : "binary") +
                ", file: " + params.filename);
        }

        if (params.s == "bmp_read") {
            bmpHandler.bmp_read(gaussBuilder, params.filename, CopyPole, p);
            logOperation(LogLevel::Info, std::string("bmp_read"), std::string("file: ") + params.filename);
        }

        if (params.s == "bin") {
            componentCalculator.bin(CopyPole, params.slice, p, params.bin_mode);
            logOperation(LogLevel::Info, std::string("bin"), 
                std::string("slice=") + std::to_string(params.slice) + 
                ", mode=" + (params.bin_mode == ThresholdMode::Peaks ? "peaks" : 
                            params.bin_mode == ThresholdMode::Valleys ? "valleys" : "all"));
        }
        
        if (params.s == "wave") {
            componentCalculator.wave(params.noisy, componenti, CopyPole, p);
            logOperation(LogLevel::Info, std::string("wave"), 
                std::string("noisy=") + std::to_string(params.noisy) + 
                ", components=" + std::to_string(componenti.size()));
            copier.removeNoise(CopyPole, componenti);
        }
         
        if (params.s == "k_means") {
            if (params.k <= 0 || p == nullptr) {
                logger.logMessage(LogLevel::Error, "Invalid parameters for k_means");
                return;
            }
            
            copier.removeNoise(CopyPole, componenti);
            prepareKMeansData(CopyPole);
            if (kMeansData.empty()) {
                logger.logMessage(LogLevel::Warning, "No data for k_means clustering");
                return;
            }

            auto result = kMeans->cluster(kMeansData, params.k);
            applyClusterResults(result, CopyPole);
            logOperation(LogLevel::Info, std::string("k_means"), std::string("k=") + std::to_string(params.k));
        }
        
        if (params.s == "k_means_kern") {
            if (!kMeans || params.k <= 0 || params.kk <= 0 || p == nullptr) {
                logger.logMessage(LogLevel::Error, "Invalid parameters for k_means_kern");
                return;
            }
            
            copier.removeNoise(CopyPole, componenti);
            prepareKMeansData(CopyPole);
            if (kMeansData.empty()) {
                logger.logMessage(LogLevel::Warning, "No data for k_means_kern clustering");
                return;
            }

            auto result = kMeans->kmeansWithKernels(kMeansData, params.k, params.kk);
            applyClusterResults(result, CopyPole);
            logOperation(LogLevel::Info, std::string("k_means_kern"), 
                std::string("k=") + std::to_string(params.k) + 
                ", kk=" + std::to_string(params.kk));
        }
        
        if (params.s == "triangulate") {
            clusterCenters = getClusterCenters();     
            lastTriangulation = componentCalculator.bowyerWatson(clusterCenters);
            voronoi.buildFromDelaunay(lastTriangulation, pathFinder, p, voronoiEdges);
            logOperation(LogLevel::Info, std::string("triangulate"), 
                std::string("clusters=") + std::to_string(clusterCenters.size()) + 
                ", triangles=" + std::to_string(lastTriangulation.size()));
        }

        if (params.s == "find_path") {
            if (!p) {
                logger.logMessage(LogLevel::Error, "Pole not initialized for find_path");
                return;
            }
            
            terrainGrid.calculateSlopes(*p);
            PointD start(params.pointA_x, params.pointA_y);
            PointD goal(params.pointB_x, params.pointB_y);
            
            if (lastTriangulation.empty()) {
                logger.logMessage(LogLevel::Error, "No triangulation for find_path");
                return;
            }
            
            copier.removeNoise(CopyPole, componenti);
            path = pathFinder.findPathAStar(start, goal, lastTriangulation, terrainGrid, CopyPole, params.slice, p);
            
            if (path.empty()) {
                logger.logMessage(LogLevel::Warning, "Path not found");
            } else {
                paths.push_back(path);
                logOperation(LogLevel::Info, std::string("find_path"), 
                    std::string("from (") + std::to_string(start.x) + "," + std::to_string(start.y) + ")" +
                    " to (" + std::to_string(goal.x) + "," + std::to_string(goal.y) + ")" +
                    ", points=" + std::to_string(path.size()));
            }
        }
        
        if (params.s == "Plot3DPath") {
            PointD start(params.pointA_x, params.pointA_y);
            PointD end(params.pointB_x, params.pointB_y);
            gnuplotInterface.plot3DPath(path, p, params.filename, start, end, pathFinder, config.vehicleRadius);
            logOperation(LogLevel::Info, std::string("Plot3DPath"), std::string("file: ") + params.filename);
        }
        
        if (params.s == "plotInteractive3DPath") {
            PointD start(params.pointA_x, params.pointA_y);
            PointD end(params.pointB_x, params.pointB_y);
            gnuplotInterface.plotInteractive3DPath(path, p, start, end, pathFinder, config.vehicleRadius);
            logOperation(LogLevel::Info, std::string("plotInteractive3DPath"));
        }
    }
};

class Interface {
private:
    Config& config;
    Logger& logger;
    Control& control;
    DispatcherParams params;
    int n = 0; // Флаг для команды init

    // Приватная функция для вывода справки
    void showHelp() {
    logger.info("Showing help information");
    
    // Содержимое справки
    const std::string helpContent = R"(
# Terrain Navigation System

Программа для анализа рельефа местности, построения триангуляции Делоне, диаграмм Вороного и поиска оптимальных маршрутов с учетом препятствий.

## 📌 Основные функции
- Генерация/загрузка карты высот (формат BMP и GNUPLOT)
- Кластеризация объектов методом k-means
- Триангуляция Делоне с учетом высот
- Построение диаграммы Вороного
- Поиск пути с ограничениями по углу наклона тележки и ее радиуса
  

## 💻️ Системные требования

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
│   ├── commandsGauss.txt               # Пример командного файла, генерирующий поле гаусами
│   ├── commandsRead.txt                # Пример командного файла, генерирующий поле картинкой Read.bmp
│   ├── main.cpp                        # Исходный код программы
│   └── auto.sh                         # Скрипт, запускающий main.cpp и читающий файл commandsGauss.txt
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
        ├── path_plot2D.png             # Маршрут 2D
        ├── path_plot3D.png             # Маршрут 3D
        ├── output_kmeans.bmp           # K-means
        ├── output_kmeans_kern.bmp      # K-means с ядрами
        ├── slice.bmp                   # Бинаризированная карта
        ├── Read.bmp                    # Поле для чтения
        └── output.bmp                  # Сгенерированная карта
```


## 🛠 Команды управления (для командного файла command.txt)

| Команда              | Параметры                      | Описание                                                                 |
|----------------------|--------------------------------|--------------------------------------------------------------------------|
| help                 | -                              | Создание файла с пояснением команд                                       |
| init                 | -                              | Инициализация поля                                                       |
| g                    | x y sx sy h                    | Создает гаусс с центром (x,y), размерами (sx,sy) и высотой h             |
| generate             | -                              | Складывает все добавленные гауссы в итоговое поле                        |
| gnuplot              | filename.png                   | Сохраняет 3D-визуализацию поля в PNG файл                                |
| PlotMetedata         | filename.png                   | Визуализирует метаданные компонент с границами и центрами                |
| PlotVoronoi          | filename.png                   | Строит диаграмму Вороного по текущей триангуляции                        |
| PlotDelaunay         | filename.png                   | Визуализирует триангуляцию Делоне                                        |
| PlotPath             | filename.png                   | Отображает найденный путь между точками A и B                            |
| bmp_write            | filename.bmp [Full/Binary]     | Сохраняет поле в BMP: Full - полное, Binary - бинаризованное             |
| bmp_read             | filename.bmp                   | Загружает поле из BMP файла                                              |
| bin                  | slice [Peaks/Valleys/All]      | Бинаризация: Peaks - только пики, Valleys - впадины, All - по модулю     |
| wave                 | noisy                          | Удаляет компоненты размером ≤ noisy как шум                              |
| k_means              | k                              | Кластеризует данные в k кластеров                                        |
| k_means_kern         | kk                             | Кластеризация с ядрами размера kk                                        |
| triangulate          | -                              | Строит триангуляцию Делоне по центрам компонент                          |
| find_path            | Ax Ay Bx By                    | Ищет путь между точками A и B через триангуляцию                         |
| Plot3DPath           | filename.png                   | Сохраняет 3D-визуализацию путя в PNG файл                                |
| plotInteractive3DPath| -                              | Интерактвный 3D режим с путем                                            |
| end                  | -                              | Завершает работу программы                                               |


## ⚙️ Параметры конфигурационного файла (config.txt)

| Параметр                     | Значение                                       | Описание                                                                         |
|------------------------------|------------------------------------------------|----------------------------------------------------------------------------------|
| fieldWidth                   | `fieldWidth`                                   | Ширина рабочего поля в пикселях                                                  |  
| fieldHeight                  | `fieldHeight`                                  | Высота рабочего поля в пикселях                                                  |
| defaultX                     | `defaultX`                                     | Стандартная X-координата центра гауссова распределения по умолчанию              |
| defaultY                     | `defaultY`                                     | Стандартная Y-координата центра гауссова распределения по умолчанию              |
| defaultSx                    | `defaultSx`                                    | Стандартное отклонение по оси X по умолчанию                                     |
| defaultSy                    | `defaultSy`                                    | Стандартное отклонение по оси Y по умолчанию                                     |
| defaultH                     | `defaultH`                                     | Стандартная высота гауссова распределения по умолчанию                           |
| defaultGnuplot               | `filename_gnuplot.png`                         | Путь к файлу для сохранения 3D-визуализации по умолчанию                         |
| defaultPlotMetedata          | `filename_metadata.png`                        | Путь к файлу для визуализации метаданных компонент по умолчанию                  |
| defaultPlotVoronoi           | `filename_voronoi.png`                         | Путь к файлу для диаграммы Вороного по умолчанию                                 |
| defaultPlotDelaunay          | `filename_delaunay.png`                        | Путь к файлу для триангуляции Делоне по умолчанию                                |
| defaultPlotPath              | `filename_path.png`                            | Путь к файлу для визуализации маршрута по умолчанию                              |
| defaultWrite                 | `filename_write.bmp`                           | Путь к файлу для сохранения BMP-изображения по умолчанию                         |
| defaultWriteModeImage        | `writeMode`                                    | Режим сохранения BMP (Full/Binary) по умолчанию                                  |
| defaultRead                  | `filename_read.bmp`                            | Путь к файлу для загрузки BMP-изображения по умолчанию                           |
| defaultSlice                 | `defaultSlice`                                 | Порог бинаризации по умолчанию                                                   |
| defaultBinMode               | `binMode`                                      | Режим бинаризации (Peaks/Valleys/All) по умолчанию                               |
| defaultNoisy                 | `defaultNoisy`                                 | Порог для удаления шумовых компонент по умолчанию                                |
| defaultKlaster               | `defaultKlaster`                               | Количество кластеров для k-mean по умолчанию                                     |
| defaultKlasterKern           | `defaultKlasterKern`                           | Размер ядра для кластеризации по умолчанию                                       |
| defaultpointA_x              | `pointA_x`                                     | X-координата точки A для поиска пути по умолчанию                                |
| defaultpointA_y              | `pointA_y`                                     | Y-координата точки A для поиска пути по умолчанию                                |
| defaultpointB_x              | `pointB_x`                                     | X-координата точки B для поиска пути по умолчанию                                |
| defaultpointB_y              | `pointB_y`                                     | Y-координата точки B для поиска пути по умолчанию                                |
| defaultPlot3DPath            | `filename_plot3dpath.png`                      | Путь к файлу для 3D-визуализации маршрута по умолчанию                           |
| vehicleRadius                | `vehicleRadius`                                | Радиус транспортного средства                                                    |
| maxSideAngle                 | `maxSideAngle`                                 | Максимальный угол поворота вбок (градусы)                                        |
| maxUpDownAngle               | `maxUpDownAngle`                               | Максимальный угол наклона вверх/вниз (градусы)                                   |
| logFileNameInterface         | `filename_log_interface.txt`                   | Путь к лог-файлу интерфейса                                                      |
| logFileNameControl           | `filename_log_control.txt`                     | Путь к лог-файлу управления                                                      |
| FiltrationLogLevelInterface  | `logLevelInterface`                            | Уровень логирования интерфейса (TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF)     |
| FiltrationLogLevelControl    | `logLevelControl`                              | Уровень логирования управления (TRACE/DEBUG/INFO/WARNING/ERROR/CRITICAL/OFF)     |


## ⚠️ Важно
1. Всегда начинайте с команды init
2. Точки A/B задаются в config.txt
3. Маршрут будет найден не всегда!
4. Для триангуляции и построения пути, нужно чтобы количество компонент было больше 3
5. Если пользуетесь программой, то важно использовать ту же файловую структуру!
6. Путь к файлу пишем полностью
7. Важен порядок команд, не забывайте делать картинки после команд
8. Для команды bmp_write не полагайтесь на значения по умолчанию
9. Точки A и B должны попадать в триангуляцию
10. Уровень "равнины" = 127, чтобы метод записи поля по гаусам согласовался с записью по картике
11. Условия проходимости: 1)Если угол наклона в любом пикселе путя превосходит допустимый (по направлению или вбок) 2) Расстояние до препятсвия на срезе меньше радиуса

## 📜 Командный файл (примеры)
1) Если нужно прочитать данные с файла Read.bmp
```
init
help
bmp_read /home/log/Gauss/results/visualizations/Read.bmp
gnuplot /home/log/Gauss/results/visualizations/gnuplot.png
bmp_write /home/log/Gauss/results/visualizations/Polew.bmp Full
bin 20 All
bmp_write /home/log/Gauss/results/visualizations/Slice.bmp Binary
wave 10
PlotMetedata /home/log/Gauss/results/visualizations/Metadata.png
k_means 5
bmp_write /home/log/Gauss/results/visualizations/kmeans.bmp Binary
triangulate
PlotVoronoi /home/log/Gauss/results/visualizations/Diagramma_Voronova.png
PlotDelaunay /home/log/Gauss/results/visualizations/Triangulation_Delone.png
find_path
PlotPath /home/log/Gauss/results/visualizations/Path.png
end
```
2) Если данные вводятся с помощью гаусов
```
init
help
g 99 50 25 25 -50
g 50 50 20 20 50
g 200 50 20 20 -50
g 50 200 20 20 50
g 189 180 20 20 -50
g 130 130 20 20 50
generate
gnuplot /home/log/Gauss/results/visualizations/gnuplot.png
bmp_write /home/log/Gauss/results/visualizations/Pole.bmp Full
bin 20 All
bmp_write /home/log/Gauss/results/visualizations/Slice.bmp Binary
wave 10
PlotMetedata /home/log/Gauss/results/visualizations/Metadata.png
k_means 5
bmp_write /home/log/Gauss/results/visualizations/kmeans.bmp Binary
triangulate
PlotVoronoi /home/log/Gauss/results/visualizations/Diagramma_Voronova.png
PlotDelaunay /home/log/Gauss/results/visualizations/Triangulation_Delone.png
find_path
PlotPath /home/log/Gauss/results/visualizations/Path.png
end
```

## 📃️ Конфигурационный файл (пример)
```
fieldWidth 250
fieldHeight 250
defaultX 50.0
defaultY 50.0
defaultSx 20.0
defaultSy 20.0
defaultH 200.0
defaultGnuplot /home/log/Gauss/results/visualizations/Gnuplot.png
defaultPlotMetedata /home/log/Gauss/results/visualizations/Metadata.png
defaultPlotVoronoi /home/log/Gauss/results/visualizations/Voronoi.png
defaultPlotDelaunay /home/log/Gauss/results/visualizations/Delaunay.png
defaultPlotPath /home/log/Gauss/results/visualizations/Path.png
defaultWrite /home/log/Gauss/results/visualizations/Write.bmp
defaultWriteModeImage Full
defaultRead /home/log/Gauss/results/visualizations/Read.bmp
defaultSlice 127
defaultBinMode All
defaultNoisy 10
defaultKlaster 5
defaultKlasterKern 5
defaultpointA_x 150.0
defaultpointA_y 150.0
defaultpointB_x 160.0
defaultpointB_y 160.0
defaultPlot3DPath /home/log/Gauss/results/visualizations/Plot3DPath.png
vehicleRadius 5
maxSideAngle 90.0
maxUpDownAngle 90.0
logFileNameInterface /home/log/Gauss/results/docs/loginterface.txt
logFileNameControl /home/log/Gauss/results/docs/logcontrol.txt
FiltrationLogLevelInterface INFO
FiltrationLogLevelControl INFO
```

##### ======================================
## 🔄 SAFE WORKFLOW v1.0
##### ======================================
#### Философия: "Мои изменения — священны 😤️, main обновляется без боли 😇️"

### 1️⃣ Начало работы (без опасного pull!)
```
git checkout main
git checkout -b feature/improved-structure  # Создаем ветку ОТ локального main
```

### 2️⃣ Работа с файлами
#### Редактируете файлы -> сохраняете -> проверяете:
```
git status
```

### 3️⃣ Фиксация изменений в ветке
```
git add .
git commit -m "Мои изменения ..."
git push origin feature/improved-structure
```

### 4️⃣ Лучший способ обновить main:

#### --- Способ 1: Идеальный мир (без конфликтов) ---
```
git checkout main
git merge --squash feature/improved-structure #Лень решать конфликты? Пишем "git reset --merge" после и идем ко 2 способу!
git commit -m "РЕЛИЗ: Новая структура проекта"
git push origin main
```

#### --- Способ 2: Жёсткий reset (когда конфликты лень решать) ---
```
git checkout main
git reset --hard feature/improved-structure
git push origin main --force-with-lease
```

#### --- Способ 3: Аккуратный мерж (если main меняли другие) ---
```
git fetch origin
git merge origin/main --no-commit
[ручное разрешение конфликтов]
git commit -m "Мерж улучшенной структуры"
git push origin main
```

### 5️⃣ Создание тега
```
git tag -a v2.0.0 -m "Обновленная структура проекта"
git push origin v2.0.0
```

### 6️⃣ Очистка
```
git branch -D feature/improved-structure
git push origin --delete feature/improved-structure
```


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
    // 1. Вывод в консоль
    std::cout << helpContent << std::endl;

    // 2. Сохранение в файл
    std::ofstream helpFile("/home/log/Gauss/results/docs/help.txt");
    
    if (helpFile) {
        helpFile << helpContent;
        logger.info("Help file created successfully at /home/log/Gauss/results/docs/help.txt");
    } else {
        logger.error("Failed to create help file at /home/log/Gauss/results/docs/help.txt");
        std::cerr << "Error: Could not create help file at /home/log/Gauss/results/docs/help.txt" << std::endl;
        
        // Попытка создать директории, если они не существуют
        system("mkdir -p /home/log/Gauss/results/docs/");
        
        // Повторная попытка
        helpFile.open("/home/log/Gauss/results/docs/help.txt");
        if (helpFile) {
            helpFile << helpContent;
            logger.info("Help file created successfully after creating directories");
        } else {
            logger.error("Failed to create help file even after creating directories");
        }
    }
}

    // Обработка команд из файла
    void processFileCommands(std::ifstream& file) {
        int commandCount = 0;
        
        while (file >> params.s) {
            commandCount++;
            logger.info(std::string("Processing command #") + std::to_string(commandCount) + ": " + params.s);
            
            if (params.s == "help") {
                showHelp();
                continue;
            }
            
            if (params.s == "end") {
                logger.info("Ending the program.");
                break;
            }
            
            if (!processCommand(file, true)) {
                break;
            }
        }
    }

    // Обработка команд с клавиатуры
    void processKeyboardCommands() {
        
        while (true) {
            std::cout << "Enter command (help, init, g, generate, gnuplot, bmp_write, bmp_read, bin, k_means, triangulate, find_path, end): ";
            std::cin >> params.s;
            std::cout << "\n";
            logger.info(std::string("Received command: ") + params.s);
            
            if (params.s == "help") {
                showHelp();
                continue;
            }
            
            if (params.s == "end") {
                std::cout << "Ending the program" << std::endl;
                logger.info("Ending the program.");
                break;
            }
            
            if (!processCommand(std::cin, false)) {
                break;
            }
        }
    }

    // Общая обработка команды (для файла и клавиатуры)
    bool processCommand(std::istream& input, bool fromFile) {
    std::string line;
    std::string modeWrite, modeBin;
    
    if (params.s == "init") {
        if (n != 0) {
            std::cout << "The init command has already been called.\nError\n";
            logger.error("Error: Multiple init commands.");
            return false;
        }
        n = 1;
        logger.info(std::string("Initializing field with size: ") + std::to_string(params.A) + " x " + std::to_string(params.B));
        control.Dispetcher(params);
        logger.info("Field initialized.");
    }
    else if (n != 1) {
        std::cout << "The init command was not used.\nError\n";
        logger.error("Error: The init command was not used.");
        return false;
    }
    else if (params.s == "g") {
        params.x = config.defaultX;
        params.y = config.defaultY;
        params.sx = config.defaultSx;
        params.sy = config.defaultSy;
        params.h = config.defaultH;

        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.x >> params.y >> params.sx >> params.sy >> params.h;
        } else {
            std::cout << "Enter parameters (x y sx sy h): ";
            std::cin >> params.x >> params.y >> params.sx >> params.sy >> params.h;
        }

        logger.info(std::string("Adding Gaussian: x=") + std::to_string(params.x) + 
                  ", y=" + std::to_string(params.y) + 
                  ", sx=" + std::to_string(params.sx) + 
                  ", sy=" + std::to_string(params.sy) + 
                  ", h=" + std::to_string(params.h));
        control.Dispetcher(params);
    }
    else if (params.s == "generate") {
        logger.info("Generating field by summing all Gaussians");
        control.Dispetcher(params);
        logger.info("Field generation completed");
    }
    else if (params.s == "gnuplot") {
        params.filename = config.defaultGnuplot;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
        } else {
            std::cout << "Enter output filename: ";
            std::cin >> params.filename;
        }

        logger.info(std::string("Calling gnuplot with filename: ") + params.filename);
        control.Dispetcher(params);
        logger.info("Gnuplot visualization completed");
    }
    else if (params.s == "PlotMetedata") {
        params.filename = config.defaultPlotMetedata;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
        } else {
            std::cout << "Enter output filename: ";
            std::cin >> params.filename;
        }

        logger.info(std::string("Plotting metadata to: ") + params.filename);
        control.Dispetcher(params);
        logger.info("Metadata plotting completed");
    }
    else if (params.s == "PlotVoronoi") {
        params.filename = config.defaultPlotVoronoi;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
        } else {
            std::cout << "Enter output filename: ";
            std::cin >> params.filename;
        }

        logger.info(std::string("Plotting Voronoi diagram to: ") + params.filename);
        control.Dispetcher(params);
        logger.info("Voronoi diagram plotting completed");
    }
    else if (params.s == "PlotDelaunay") {
        params.filename = config.defaultPlotDelaunay;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
        } else {
            std::cout << "Enter output filename: ";
            std::cin >> params.filename;
        }

        logger.info(std::string("Plotting Delaunay triangulation to: ") + params.filename);
        control.Dispetcher(params);
        logger.info("Delaunay triangulation plotting completed");
    }
    else if (params.s == "PlotPath") {
        params.filename = config.defaultPlotPath;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
        } else {
            std::cout << "Enter output filename: ";
            std::cin >> params.filename;
        }

        logger.info(std::string("Plotting path to: ") + params.filename);
        control.Dispetcher(params);
        logger.info("Path plotting completed");
    }
    else if (params.s == "bmp_write") {
        params.filename = config.defaultWrite;
        modeWrite = config.defaultWriteModeImage;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename >> modeWrite;
        } else {
            std::cout << "Enter output filename and mode (Full/Binary): ";
            std::cin >> params.filename >> modeWrite;
        }

        if (modeWrite == "Full") {
            params.bmp_mode = BmpWriteMode::Full;
        } else {
            params.bmp_mode = BmpWriteMode::Binary;
        }

        logger.info(std::string("Writing BMP file: ") + params.filename + " with mode: " + modeWrite);
        control.Dispetcher(params);
        logger.info("BMP writing completed");
    }
    else if (params.s == "bmp_read") {
        params.filename = config.defaultRead;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
        } else {
            std::cout << "Enter input filename: ";
            std::cin >> params.filename;
        }

        logger.info(std::string("Reading BMP file: ") + params.filename);
        control.Dispetcher(params);
        logger.info("BMP reading completed");
    }
    else if (params.s == "bin") {
        params.slice = config.defaultSlice;
        modeBin = config.defaultBinMode;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.slice >> modeBin;
        } else {
            std::cout << "Enter slice level and mode (Peaks/Valleys/All): ";
            std::cin >> params.slice >> modeBin;
        }

        if (modeBin == "Peaks") {
            params.bin_mode = ThresholdMode::Peaks;
        } else if (modeBin == "Valleys") {
            params.bin_mode = ThresholdMode::Valleys;
        } else {
            params.bin_mode = ThresholdMode::All;
        }

        logger.info(std::string("Applying binary filter with slice: ") + std::to_string(params.slice) + 
                   " and mode: " + modeBin);
        control.Dispetcher(params);
        logger.info("Binary filtering completed");
    }
    else if (params.s == "wave") {
        params.noisy = config.defaultNoisy;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.noisy;
        } else {
            std::cout << "Enter noisy level: ";
            std::cin >> params.noisy;
        }

        logger.info(std::string("Applying wave filter with noisy level: ") + std::to_string(params.noisy));
        control.Dispetcher(params);
        logger.info(std::string("Wave filtering completed. Components count: ") + 
                   std::to_string(control.componenti.size()));
    }
    else if (params.s == "k_means") {
        params.k = config.defaultKlaster;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.k;
        } else {
            std::cout << "Enter cluster count: ";
            std::cin >> params.k;
        }

        logger.info(std::string("Running k-means with k: ") + std::to_string(params.k));
        control.Dispetcher(params);
        logger.info("k-means clustering completed");
    }
    else if (params.s == "k_means_kern") {
        params.kk = config.defaultKlasterKern;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.kk;
        } else {
            std::cout << "Enter kernel size: ";
            std::cin >> params.kk;
        }

        logger.info(std::string("Running k-means with kernel size: ") + std::to_string(params.kk));
        control.Dispetcher(params);
        logger.info("k-means with kernels completed");
    }
    else if (params.s == "triangulate") {
        logger.info("Starting Delaunay triangulation");
        control.Dispetcher(params);
        logger.info("Delaunay triangulation completed");
    }
    else if (params.s == "find_path") {
        params.pointA_x = config.defaultpointA_x;
        params.pointA_y = config.defaultpointA_y;
        params.pointB_x = config.defaultpointB_x;
        params.pointB_y = config.defaultpointB_y;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.pointA_x >> params.pointA_y >> params.pointB_x >> params.pointB_y;
        } else {
            std::cout << "Enter points A(x y) and B(x y): ";
            std::cin >> params.pointA_x >> params.pointA_y >> params.pointB_x >> params.pointB_y;
        }

        logger.info(std::string("Finding path from (") + std::to_string(params.pointA_x) + "," + 
                   std::to_string(params.pointA_y) + ") to (" + 
                   std::to_string(params.pointB_x) + "," + 
                   std::to_string(params.pointB_y) + ")");
        control.Dispetcher(params);
        logger.info("Path finding completed");
    }
    else if (params.s == "Plot3DPath") {
        params.filename = config.defaultPlot3DPath;
        
        if (fromFile) {
            std::getline(input, line);
            std::istringstream iss(line);
            iss >> params.filename;
        } else {
            std::cout << "Enter output filename: ";
            std::cin >> params.filename;
        }

        logger.info(std::string("Plotting 3D path to: ") + params.filename);
        control.Dispetcher(params);
        logger.info("3D path plotting completed");
    }
    else if (params.s == "plotInteractive3DPath") {
        logger.info("Starting interactive 3D path visualization");
        control.Dispetcher(params);
        logger.info("Interactive 3D visualization completed");
    }
    else {
        std::cout << "Unknown command: " << params.s << std::endl;
        logger.warning(std::string("Unknown command received: ") + params.s);
        return false;
    }
    
    return true;
}

public:
    Interface(Config& cfg, Logger& log, Control& c) 
        : config(cfg), logger(log), control(c) {
        logger.info("Interface initialized");
        logger.debug(std::string("Field dimensions: ") + std::to_string(config.fieldWidth) + 
                   "x" + std::to_string(config.fieldHeight));
        
        params.A = config.fieldWidth;
        params.B = config.fieldHeight;
    }
    
    void print() {
        logger.info("Starting user interface");
        
        bool useFileInput;
        std::string commandfilename;
        
        std::cout << "Hello, dear user, this program builds Gaussians.\n"
                  << "Enter commands from a text file (PRESS 1) or from the keyboard (PRESS 0)?" 
                  << std::endl;
        std::cin >> useFileInput;
        
        logger.info(std::string("User chose input method: ") + 
                  std::string(useFileInput ? "file" : "keyboard"));
        
        if (useFileInput) {
            std::cout << "You will enter commands from a text file.\nEnter filename:" << std::endl;
            std::cin >> commandfilename;
            
            logger.info(std::string("Attempting to open command file: ") + commandfilename);
            std::ifstream file(commandfilename);
            
            if (!file) {
                logger.error(std::string("Failed to open command file: ") + commandfilename);
                std::cout << "File not found" << std::endl;
                return;
            }
            
            processFileCommands(file);
            file.close();
            logger.info("Closed input file.");
        } else {
            logger.info("User selected keyboard input mode");
            std::cout << "You will enter commands from the keyboard" << std::endl;
            processKeyboardCommands();
        }
    }
};

int main() {
    //Логирование
    Config config("config.txt");
    Logger loggerinterface(config.logFileNameInterface, "Interface", config.FiltrationLogLevelInterface);
    Logger loggercontrol(config.logFileNameControl, "Control", config.FiltrationLogLevelControl);
    // Создаем интерфейс
    Control c(config, loggercontrol);
    Interface i(config, loggerinterface, c);
    
    // Вызываем метод print() интерфейса
    i.print();

    return 0;
}
