/*Здраствуйте! В этой версии есть недоточеты! Осталось
1) Если добавятся новые команды меняем readme file и help 
2)Улучшить логирование
3) Возможно новые требования
4)Оптимизация

Изменения:  храним список путей (чтобы в начале применили команду path сколько душе угодно, а потом их распечатали разом, но эта распечатка не реализована), исправлена диаграмма воронова + 
улучшена визуализация + теперь условия проходимости коррекны ( условие с радиусом и с наклонами в стороны) + добавлены параметры в команду path + улучшено логирование (но не до конца) + 
 правильная визуализация метаданных компонент + убрал точки A и B из триангулизации + Визуализация маршрута 3D новая команда + углы правильно считаются и логируются

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
    std::string logFileNameInterface;
    std::string logFileNameControl;
    std::string defaultGnuplot, defaultPlotMetedata, defaultPlotVoronoi, defaultPlotDelaunay, defaultPlotPath, defaultWrite, defaultRead, defaultWriteModeImage, defaultBinMode, defaultPlot3DPath;
    int defaultSlice, defaultNoisy, defaultKlaster, defaultKlasterKern;
    bool loggingInterfaceEnabled;
    bool loggingControlEnabled;
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
    else if (key == "loggingInterfaceEnabled") configFile >> std::boolalpha >> loggingInterfaceEnabled;
    else if (key == "loggingControlEnabled") configFile >> std::boolalpha >> loggingControlEnabled;
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
    Gaus(double h, double x0, double y0, double sigma_x, double sigma_y) : h(h), x0(x0), y0(y0), sigma_x(sigma_x), sigma_y(sigma_y) {}
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
    int pixelCount;  // Поле для хранения количества пикселей
    
   // Основной конструктор (для готовой матрицы)
    Component(const std::vector<std::vector<double>>& inputComponenta, int count) : componenta(inputComponenta), pixelCount(count) {
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
        std::cerr << "Error: Component has no valid pixels!" << std::endl;
        center_x = center_y = 0;
        return;
    }

    center_x = sum_x / count;
    center_y = sum_y / count;

    // Корректный расчет ковариационной матрицы
    double cov_xx = 0, cov_xy = 0, cov_yy = 0;
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

    // Собственные значения
    double trace = cov_xx + cov_yy;
    double det = cov_xx * cov_yy - cov_xy * cov_xy;
    eigenvalue1 = (trace + sqrt(trace * trace - 4 * det)) / 2;
    eigenvalue2 = (trace - sqrt(trace * trace - 4 * det)) / 2;

    // Корректный расчет собственных векторов
    if (fabs(cov_xy) > 1e-10) {
        // Первый собственный вектор (для eigenvalue1)
        eigenvec1_x = cov_xy;
        eigenvec1_y = eigenvalue1 - cov_xx;
        
        // Второй собственный вектор (для eigenvalue2)
        eigenvec2_x = cov_xy;
        eigenvec2_y = eigenvalue2 - cov_xx;
    } else {
        eigenvec1_x = 1; eigenvec1_y = 0;
        eigenvec2_x = 0; eigenvec2_y = 1;
    }

    // Нормализация векторов
    double norm1 = sqrt(eigenvec1_x * eigenvec1_x + eigenvec1_y * eigenvec1_y);
    double norm2 = sqrt(eigenvec2_x * eigenvec2_x + eigenvec2_y * eigenvec2_y);
    
    eigenvec1_x /= norm1; eigenvec1_y /= norm1;
    eigenvec2_x /= norm2; eigenvec2_y /= norm2;
}
};

class Copier {
public:
    void removeNoise(std::vector<std::vector<double>>& field, const std::vector<Component>& components) {
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
        std::fill(row.begin(), row.end(), 127); // Уровень равнины
    }
        double value; 
        for (const auto& g : gaussi) {
    for (long unsigned int x = 0; x < p->field[0].size(); ++x) {
        for (long unsigned int y = 0; y < p->field.size(); ++y) {
            value = g.h * exp(-((pow((x - g.x0) / g.sigma_x, 2)) + (pow((y - g.y0) / g.sigma_y, 2))) / 2); 
            p->field[y][x] += value;
           // Ограничиваем значения между 0 и 255
                p->field[y][x] = std::clamp(p->field[y][x], 0.0, 255.0);
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
            unsigned char color = static_cast<unsigned char>(pixelMatrix[y][x]); // Color
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

// Класс для работы с сеткой рельефа
class TerrainGrid {
public:
TerrainGrid(int width, int height) {
        initialize(width, height);
    }
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
    
    std::pair<double, double> getEdgeSlopes(const Edge& edge) const {
    // Вычисляем направление ребра
    PointD dir = {
        edge.b.x - edge.a.x,
        edge.b.y - edge.a.y
    };
    double length = std::hypot(dir.x, dir.y);
    if (length < 1e-6) return {0, 0};
    
    // Нормализуем направление
    dir.x /= length;
    dir.y /= length;
    
    // Средняя точка ребра
    int x = static_cast<int>((edge.a.x + edge.b.x) / 2);
    int y = static_cast<int>((edge.a.y + edge.b.y) / 2);
    
    // Проверяем границы
    if (x <= 0 || y <= 0 || x >= cells[0].size()-1 || y >= cells.size()-1) {
        return {0, 0};
    }

    // Вычисляем углы наклона
    double forwardSlope = atan(cells[y][x].slopeX * dir.x + cells[y][x].slopeY * dir.y);
    double sideSlope = atan(cells[y][x].slopeX * (-dir.y) + cells[y][x].slopeY * dir.x);
    
    return {
        forwardSlope * 180.0 / M_PI,  // Преобразуем в градусы
        sideSlope * 180.0 / M_PI
    };
}

    // Расчет углов наклона на основе данных высот
    void calculateSlopes(const Pole& elevationData) {
    if (elevationData.field.empty()) {
        throw std::invalid_argument("Pole data is empty");
    }

    const int height = cells.size();
    const int width = cells[0].size();

    for (int y = 1; y < height-1; ++y) {
        for (int x = 1; x < width-1; ++x) {
            // Сохраняем высоту
            cells[y][x].height = elevationData.field[y][x];
            
            // Вычисляем производные (центральные разности)
            cells[y][x].slopeX = (elevationData.field[y][x+1] - elevationData.field[y][x-1]) / 2.0;
            cells[y][x].slopeY = (elevationData.field[y+1][x] - elevationData.field[y-1][x]) / 2.0;
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

class PathFinder {
private:
    const Config& config;
public:
    PathFinder(const Config& cfg) : config(cfg) {}
    
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
    
  bool isNavigable(const Edge& edge, const TerrainGrid& terrainGrid) const {
    auto [forwardAngle, sideAngle] = terrainGrid.getEdgeSlopes(edge);
    
    return abs(forwardAngle) <= config.maxUpDownAngle && 
           abs(sideAngle) <= config.maxSideAngle &&
           edge.length() > config.vehicleRadius * 2;
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
    std::vector<PointD> findPathAStar(const PointD& start, const PointD& goal, 
                                            const std::vector<Triangle>& triangles, 
                                            const TerrainGrid& terrainGrid,
                                            const std::vector<std::vector<double>>& binaryMap,
                                            Logger& logger, bool loggingEnabled) {  // Добавляем параметр логгера
    std::vector<PointD> path;

    logger.logMessage("=== НАЧАЛО ПОИСКА ПУТИ ===", loggingEnabled);
    logger.logMessage("Старт: (" + std::to_string(start.x) + ", " + std::to_string(start.y) + ")", loggingEnabled);
    logger.logMessage("Цель: (" + std::to_string(goal.x) + ", " + std::to_string(goal.y) + ")", loggingEnabled);

    // Найти стартовый и целевой треугольники
    const Triangle* startTri = findContainingTriangle(start, triangles);
    const Triangle* goalTri = findContainingTriangle(goal, triangles);
    
    if (!startTri) {
        std::string msg = "Стартовая точка (" + std::to_string(start.x) + ", " + std::to_string(start.y) + ") вне триангуляции!";
        logger.logMessage(msg, loggingEnabled);
        return {};
    }
    if (!goalTri) {
        std::string msg = "Конечная точка (" + std::to_string(goal.x) + ", " + std::to_string(goal.y) + ") вне триангуляции!";
        logger.logMessage(msg, loggingEnabled);
        return {};
    }

    // Проверка для случая, когда start и goal в одном треугольнике
    if (startTri == goalTri) {
        logger.logMessage("Старт и цель находятся в одном треугольнике", loggingEnabled);
        path = {start, goal};
        if (!checkPathFeasibility(path, terrainGrid, binaryMap, logger, loggingEnabled)) {  // Передаем логгер
            logger.logMessage("Прямой путь НЕ проходим!", loggingEnabled);
            return {};
        }
        logger.logMessage("Прямой путь проходим", loggingEnabled);
        return path;
    }

    // Основной алгоритм A*
    logger.logMessage("Запуск алгоритма A*...", loggingEnabled);
    std::priority_queue<AStarNode> openSet;
    std::unordered_map<const Triangle*, double> costSoFar;

    if (triangles.empty()) {
        logger.logMessage("Ошибка: триангуляция не построена!", loggingEnabled);
        return {};
    }

    PointD startCenter = startTri->calculateCircumcenter();
    openSet.emplace(startCenter, startTri);
    costSoFar[startTri] = 0;
    logger.logMessage("Начальный узел добавлен в очередь", loggingEnabled);

    while (!openSet.empty()) {
        AStarNode current = openSet.top();
        openSet.pop();

        if (current.tri == goalTri) {
            logger.logMessage("Целевой треугольник достигнут!", loggingEnabled);
            // Восстановление пути
            while (current.parent) {
                path.push_back(current.position);
                current = *current.parent;
            }
            std::reverse(path.begin(), path.end());
            
            // Добавляем старт и цель
            if (!path.empty()) {
                path.insert(path.begin(), start);
                path.push_back(goal);
                
                logger.logMessage("Проверка проходимости пути...", loggingEnabled);
                if (!checkPathFeasibility(path, terrainGrid, binaryMap, logger, loggingEnabled)) {
                    logger.logMessage("Путь НЕ проходим!", loggingEnabled);
                    return {};
                }
                logger.logMessage("Путь успешно проверен и проходим!", loggingEnabled);
            }
            return path;
        }

        // Обработка соседей
        for (const auto& neighbor : getNeighbors(*current.tri, triangles)) {
            Edge edgeToCheck(current.position, neighbor->calculateCircumcenter());
            if (!isNavigable(edgeToCheck, terrainGrid)) {
                logger.logMessage("Ребро (" + std::to_string(edgeToCheck.a.x) + "," + std::to_string(edgeToCheck.a.y) + 
                                ")-(" + std::to_string(edgeToCheck.b.x) + "," + std::to_string(edgeToCheck.b.y) + 
                                ") непроходимо", loggingEnabled);
                continue;
            }

            PointD neighborCenter = neighbor->calculateCircumcenter();
            double newCost = current.g + heuristic(current.position, neighborCenter);

            if (!costSoFar.count(neighbor) || newCost < costSoFar[neighbor]) {
                costSoFar[neighbor] = newCost;
                double priority = newCost + heuristic(neighborCenter, goal);
                openSet.emplace(neighborCenter, neighbor, new AStarNode(current));
                logger.logMessage("Добавлен новый узел в очередь: (" + 
                                std::to_string(neighborCenter.x) + "," + 
                                std::to_string(neighborCenter.y) + ")", loggingEnabled);
            }
        }
    }

    logger.logMessage("Путь не найден!", loggingEnabled);
    return {};
}

bool checkPathFeasibility(const std::vector<PointD>& path,
                                    const TerrainGrid& terrainGrid,
                                    const std::vector<std::vector<double>>& binaryMap,
                                    Logger& logger, bool loggingEnabled) const {
    logger.logMessage("=== ПРОВЕРКА ПРОХОДИМОСТИ ПУТИ ===", loggingEnabled);
    logger.logMessage("Всего сегментов: " + std::to_string(path.size()-1), loggingEnabled);

    for (size_t i = 0; i < path.size() - 1; ++i) {
        logger.logMessage("Анализ сегмента " + std::to_string(i+1) + ": (" + 
                        std::to_string(path[i].x) + "," + std::to_string(path[i].y) + ") -> (" +
                        std::to_string(path[i+1].x) + "," + std::to_string(path[i+1].y) + ")", loggingEnabled);
        
        const auto linePixels = bresenhamLine(path[i], path[i+1]);
        logger.logMessage("Пикселей в сегменте: " + std::to_string(linePixels.size()), loggingEnabled);
        
        for (const auto& pixel : linePixels) {
            // 1. Проверка углов наклона
            auto [forwardAngle, sideAngle] = getVehicleSlopeAngles(pixel, path[i+1], terrainGrid);
            
            std::string angleInfo = "Пиксель (" + std::to_string(pixel.x) + "," + std::to_string(pixel.y) + 
                                  "): forward=" + std::to_string(forwardAngle) + 
                                  "°, side=" + std::to_string(sideAngle) + "°";
            logger.logMessage(angleInfo, loggingEnabled);

            if (forwardAngle > config.maxUpDownAngle || sideAngle > config.maxSideAngle) {
                std::string msg = "ПРЕВЫШЕН УГОЛ НАКЛОНА в пикселе (" + 
                                std::to_string(pixel.x) + "," + std::to_string(pixel.y) + 
                                "): forward=" + std::to_string(forwardAngle) + 
                                " (max " + std::to_string(config.maxUpDownAngle) + 
                                "), side=" + std::to_string(sideAngle) + 
                                " (max " + std::to_string(config.maxSideAngle) + ")";
                logger.logMessage(msg, loggingEnabled);
                return false;
            }

            // 2. Проверка радиуса тележки
            if (!isVehicleRadiusValid(pixel, binaryMap)) {
                std::string msg = "СТОЛКНОВЕНИЕ в пикселе (" + 
                                std::to_string(pixel.x) + "," + std::to_string(pixel.y) + 
                                "), радиус тележки: " + std::to_string(config.vehicleRadius);
                logger.logMessage(msg, loggingEnabled);
                return false;
            }
        }
    }
    
    logger.logMessage("Все проверки пройдены успешно!", loggingEnabled);
    return true;
}

std::vector<PointD> bresenhamLine(const PointD& start, const PointD& end) const {
    std::vector<PointD> linePoints;
    
    int x0 = static_cast<int>(start.x);
    int y0 = static_cast<int>(start.y);
    int x1 = static_cast<int>(end.x);
    int y1 = static_cast<int>(end.y);

    int dx = abs(x1 - x0);
    int dy = -abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1;
    int sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;

    while (true) {
        linePoints.emplace_back(x0, y0);
        if (x0 == x1 && y0 == y1) break;
        
        int e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
        }
    }
    return linePoints;
}

std::pair<double, double> getVehicleSlopeAngles(
    const PointD& pixel, 
    const PointD& nextPixel, 
    const TerrainGrid& grid
) const {
    if (pixel == nextPixel) return {0, 0};

    // Проверяем границы сетки
    int x = static_cast<int>(pixel.x);
    int y = static_cast<int>(pixel.y);
    if (x <= 0 || y <= 0 || x >= grid.cells[0].size()-1 || y >= grid.cells.size()-1) {
        return {90, 90}; // Края сетки считаем непроходимыми
    }

    // Вектор направления движения (нормализованный)
    PointD dir = {
        nextPixel.x - pixel.x,
        nextPixel.y - pixel.y
    };
    double length = std::hypot(dir.x, dir.y);
    if (length < 1e-6) return {0, 0};
    dir.x /= length;
    dir.y /= length;

    // Получаем высоты в текущей точке и соседних точках
    double z = grid.cells[y][x].height; 
    double z_right = grid.cells[y][x+1].height;
    double z_left = grid.cells[y][x-1].height;
    double z_top = grid.cells[y-1][x].height;
    double z_bottom = grid.cells[y+1][x].height;

    // Вычисляем производные по x и y
    double dzdx = (z_right - z_left) / 2.0;
    double dzdy = (z_bottom - z_top) / 2.0;

    // Угол наклона по направлению движения (forward)
    double forwardAngle = std::atan(dzdx * dir.x + dzdy * dir.y) * 180.0 / M_PI;

    // Перпендикулярное направление (для side angle)
    PointD perp = {-dir.y, dir.x};
    double sideAngle = std::atan(dzdx * perp.x + dzdy * perp.y) * 180.0 / M_PI;

    return {
        std::abs(forwardAngle),  // Угол вдоль направления движения
        std::abs(sideAngle)      // Угол поперек направления движения
    };
}

bool isVehicleRadiusValid(const PointD& pixel, const std::vector<std::vector<double>>& binaryMap) const {
    int x = static_cast<int>(pixel.x);
    int y = static_cast<int>(pixel.y);
    int radius = static_cast<int>(config.vehicleRadius);

    // Проверяем все пиксели в радиусе
    for (int dy = -radius; dy <= radius; ++dy) {
        for (int dx = -radius; dx <= radius; ++dx) {
            int nx = x + dx;
            int ny = y + dy;
            if (nx >= 0 && ny >= 0 && nx < binaryMap[0].size() && ny < binaryMap.size()) {
                if (binaryMap[ny][nx] == 255) { // Препятствие
                    double dist = std::hypot(dx, dy);
                    if (dist < config.vehicleRadius) {
                        return false;
                    }
                }
            }
        }
    }
    return true;
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
    
    void buildFromDelaunay(const std::vector<Triangle>& triangles, const PathFinder& pathFinder, const std::unique_ptr<Pole>& p, std::vector<VoronoiEdge>& edges) {
    const int height = p->field.size();
    const int width = p->field[0].size();
    edges.clear();

    if (!p) {
        std::cerr << "Pole is not initialized!" << std::endl;
        return;
    }

    for (const auto& tri : triangles) {
            PointD cc1 = tri.calculateCircumcenter();
            bool cc1_valid = isPointInsideField(cc1, width, height);
            
            auto neighbors = pathFinder.findNeighbors(tri, triangles);

            // Обработка обычных ребер
            for (const auto& neighbor : neighbors) {
                PointD cc2 = neighbor->calculateCircumcenter();
                bool cc2_valid = isPointInsideField(cc2, width, height);

                if (cc1_valid && cc2_valid) {
                    // Оба центра внутри - просто добавляем ребро
                    edges.emplace_back(cc1, cc2);
                } 
                else if (cc1_valid || cc2_valid) {
                    // Только один центр внутри - обрезаем ребро
                    auto clipped = clipEdgeToField(cc1, cc2, width, height);
                    if (!clipped.empty()) {
                        edges.emplace_back(clipped[0], clipped[1]);
                    }
                }
                // Оба центра снаружи - игнорируем
            }

            // Обработка граничных треугольников
            if (neighbors.size() < 3) {
                handleBoundaryTriangle(tri, cc1, width, height, edges, cc1_valid, triangles, pathFinder);
            }
        }
    }

private:
    std::vector<PointD> clipEdgeToField(const PointD& p1, const PointD& p2, int width, int height) {
    // Реализация алгоритма Коэна-Сазерленда для отсечения отрезка
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

    while (true) {
        if (!(code1 | code2)) return {a, b}; // Полностью внутри
        if (code1 & code2) return {}; // Полностью снаружи

        int outcode = code1 ? code1 : code2;
        PointD p;

        // Находим точку пересечения
        if (outcode & 8) { // Верхняя граница
            p.x = a.x + (b.x - a.x) * (height - a.y) / (b.y - a.y);
            p.y = height;
        }
        else if (outcode & 4) { // Нижняя граница
            p.x = a.x + (b.x - a.x) * (-a.y) / (b.y - a.y);
            p.y = 0;
        }
        else if (outcode & 2) { // Правая граница
            p.y = a.y + (b.y - a.y) * (width - a.x) / (b.x - a.x);
            p.x = width;
        }
        else if (outcode & 1) { // Левая граница
            p.y = a.y + (b.y - a.y) * (-a.x) / (b.x - a.x);
            p.x = 0;
        }

        // Обновляем внешнюю точку
        if (outcode == code1) {
            a = p;
            code1 = code(a);
        }
        else {
            b = p;
            code2 = code(b);
        }
    }
}

bool isPointInsideField(const PointD& p, int width, int height) {
        return p.x >= 0 && p.y >= 0 && p.x < width && p.y < height;
    }

    void handleBoundaryTriangle(const Triangle& tri, const PointD& cc, int width, int height, std::vector<VoronoiEdge>& edges, bool cc_valid,
                                const std::vector<Triangle>& allTriangles,  const PathFinder& pathFinder) {
        if (!cc_valid) return;

        for (const auto& edge : getBoundaryEdges(tri, allTriangles, pathFinder)) {
            PointD boundaryPoint = calculateBoundaryIntersection(cc, edge, width, height);
            edges.emplace_back(cc, boundaryPoint);
        }
    }

    std::vector<Edge> getBoundaryEdges(const Triangle& tri, const std::vector<Triangle>& allTriangles,  const PathFinder& pathFinder) {
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
 
  // Преобразование Y-координаты для Gnuplot
    double transformY(double y, int height) const {
        return height - y - 1;
    }
    
   public:
    void plotBinaryWithComponents(const std::vector<std::vector<double>> &CopyPole, const std::vector<Component>& components, const std::string& filename) {

        FILE* gnuplotPipe = popen("gnuplot -persist", "w");
        if (!gnuplotPipe) {
            std::cerr << "Failed to open gnuplot pipe" << std::endl;
            return;
        }

        const int height = CopyPole.size();
        const int width = CopyPole[0].size();

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
    }
      void gnuplot(std::unique_ptr<Pole>& p, const std::string& filename) {
    if (p == nullptr) {
        std::cout << "Pole not initialized!" << std::endl;
        return;
    }

    int rows = p->field.size();
    int cols = p->field[0].size();
    
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        std::cerr << "Could not open pipe to gnuplot." << std::endl;
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
}

void plotVoronoi(const std::unique_ptr<Pole>& p, const std::vector<VoronoiEdge>& edges, const std::vector<PointD>& sites, const std::string& filename) {
    if (!p || edges.empty()) return;

    const int height = p->field.size();
    const int width = p->field[0].size();

    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        std::cerr << "Failed to open gnuplot pipe" << std::endl;
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
}
    
    void plotDelaunay(const std::vector<Triangle>& triangles, std::unique_ptr<Pole>& p, const std::string& filename) {
    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        std::cerr << "Failed to open gnuplot pipe" << std::endl;
        return;
    }
    
    if (!p) return;
    const int height = p->field.size();
    const int width = p->field[0].size();

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
}

void plotPath(const std::vector<PointD>& path, const std::unique_ptr<Pole>& p, const std::string& filename, DispatcherParams& params) {
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
    fprintf(gnuplotPipe, "set output '%s'\n", filename.c_str());
    fprintf(gnuplotPipe, "set title 'Path Visualization'\n");
    fprintf(gnuplotPipe, "set size ratio -1\n");
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", width-1);
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", height-1);
    fprintf(gnuplotPipe, "unset key\n");

    // Многослойный график: фон + путь + точки
    fprintf(gnuplotPipe, "plot '-' matrix with image, \\\n");
    fprintf(gnuplotPipe, "'-' with lines lw 2 lc 'green', \\\n");
    fprintf(gnuplotPipe, "'-' with points pt 7 ps 2 lc 'blue', \\\n");
    fprintf(gnuplotPipe, "'-' with points pt 9 ps 2 lc 'purple'\n");

    // 1. Данные фона (инверсия Y)
    for (int y = height-1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            fprintf(gnuplotPipe, "%f ", p->field[y][x]);
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
            params.pointA_x, 
            transformY(params.pointA_y, height));
    fprintf(gnuplotPipe, "e\n");

    // 4. Точка B (зеленая)
    fprintf(gnuplotPipe, "%f %f\n", 
            params.pointB_x, 
            transformY(params.pointB_y, height));
    fprintf(gnuplotPipe, "e\n");

    pclose(gnuplotPipe);
}

void plotInteractive3DPath(const std::vector<PointD>& path, const std::unique_ptr<Pole>& p, const PointD& start, const PointD& end) {
    if (!p || path.empty()) return;

    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        std::cerr << "Failed to open gnuplot pipe." << std::endl;
        return;
    }

    const int height = p->field.size();
    const int width = p->field[0].size();
    
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
        PointD from = path[i];
        PointD to = path[i+1];
        
        // Алгоритм Брезенхема для рисования линии
        int x0 = static_cast<int>(from.x);
        int y0 = static_cast<int>(from.y);
        int x1 = static_cast<int>(to.x);
        int y1 = static_cast<int>(to.y);
        
        int dx = abs(x1 - x0);
        int dy = abs(y1 - y0);
        int sx = (x0 < x1) ? 1 : -1;
        int sy = (y0 < y1) ? 1 : -1;
        int err = dx - dy;
        
        while (true) {
            if (x0 >= 0 && y0 >= 0 && x0 < width && y0 < height) {
                pathProjection.emplace_back(x0, y0);
            }
            
            if (x0 == x1 && y0 == y1) break;
            
            int e2 = 2 * err;
            if (e2 > -dy) {
                err -= dy;
                x0 += sx;
            }
            if (e2 < dx) {
                err += dx;
                y0 += sy;
            }
        }
    }
    
    // 3. Проекция пути на поле (зеленые точки)
    for (const auto& [x, y] : pathProjection) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size screen 0.003 fc rgb '#00FF00' fs solid front\n",
                x, static_cast<int>(transformY(y, height)), p->field[y][x] + 0.1);
    }

    // 4. Точки A и B
    int startX = static_cast<int>(start.x);
    int startY = static_cast<int>(start.y);
    int endX = static_cast<int>(end.x);
    int endY = static_cast<int>(end.y);
    
    if (startX >= 0 && startY >= 0 && startX < width && startY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size screen 0.008 fc rgb '#FF0000' fs solid front\n",
                startX, static_cast<int>(transformY(startY, height)), 
                p->field[startY][startX] + 0.2);
    }
    
    fprintf(gnuplotPipe, "set label 'START' at %d,%d,%f front\n",
                startX, static_cast<int>(transformY(startY, height)), p->field[startY][startX] + 1);
    
    if (endX >= 0 && endY >= 0 && endX < width && endY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size screen 0.008 fc rgb '#0000FF' fs solid front\n",
                endX, static_cast<int>(transformY(endY, height)), 
                p->field[endY][endX] + 0.2);
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
}

void plot3DPath(const std::vector<PointD>& path, const std::unique_ptr<Pole>& p, const std::string& filename, const PointD& start, const PointD& end) {
    if (!p || path.empty()) return;

    FILE* gnuplotPipe = popen("gnuplot -persist", "w");
    if (!gnuplotPipe) {
        std::cerr << "Failed to open gnuplot pipe." << std::endl;
        return;
    }

    const int height = p->field.size();
    const int width = p->field[0].size();

    // 1. Создаем проекцию пути (все пиксели между узлами)
    std::vector<std::pair<int, int>> pathProjection;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        PointD from = path[i];
        PointD to = path[i+1];
        
        // Алгоритм Брезенхема для рисования линии
        int x0 = static_cast<int>(from.x);
        int y0 = static_cast<int>(from.y);
        int x1 = static_cast<int>(to.x);
        int y1 = static_cast<int>(to.y);
        
        int dx = abs(x1 - x0);
        int dy = abs(y1 - y0);
        int sx = (x0 < x1) ? 1 : -1;
        int sy = (y0 < y1) ? 1 : -1;
        int err = dx - dy;
        
        while (true) {
            if (x0 >= 0 && y0 >= 0 && x0 < width && y0 < height) {
                pathProjection.emplace_back(x0, y0);
            }
            
            if (x0 == x1 && y0 == y1) break;
            
            int e2 = 2 * err;
            if (e2 > -dy) {
                err -= dy;
                x0 += sx;
            }
            if (e2 < dx) {
                err += dx;
                y0 += sy;
            }
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
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size screen 0.003 fc rgb '#00FF00' fs solid front\n",
                x, static_cast<int>(transformY(y, height)), p->field[y][x] + 0.1);
    }

    // 4. Точки A и B
    int startX = static_cast<int>(start.x);
    int startY = static_cast<int>(start.y);
    int endX = static_cast<int>(end.x);
    int endY = static_cast<int>(end.y);
    
    if (startX >= 0 && startY >= 0 && startX < width && startY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size screen 0.008 fc rgb '#FF0000' fs solid front\n",
                startX, static_cast<int>(transformY(startY, height)), 
                p->field[startY][startX] + 0.2);
    }
    
    fprintf(gnuplotPipe, "set label 'START' at %d,%d,%f front\n",
                startX, static_cast<int>(transformY(startY, height)), p->field[startY][startX] + 1);
    
    if (endX >= 0 && endY >= 0 && endX < width && endY < height) {
        fprintf(gnuplotPipe, "set object circle at %d,%d,%f size screen 0.008 fc rgb '#0000FF' fs solid front\n",
                endX, static_cast<int>(transformY(endY, height)), 
                p->field[endY][endX] + 0.2);
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
}

};
   
class ComponentCalculator {
public:
    std::vector<Triangle> bowyerWatson(const std::vector<PointD>& points) {
    std::vector<Triangle> triangles;

    // Проверка на минимальное количество точек
    if (points.size() < 3) {
        std::cerr << "Недостаточно точек для триангуляции!" << std::endl;
        return triangles;
    }

    // Проверка на NaN-точки и коллинеарность
    for (const auto& p : points) {
        if (std::isnan(p.x) || std::isnan(p.y)) {
            std::cerr << "Обнаружена некорректная точка (NaN)!" << std::endl;
            return triangles;
        }
    }

    // Проверка на коллинеарность всех точек
    if (std::all_of(points.begin(), points.end(), 
                   [&](const PointD& p) { return Triangle::areCollinear(points[0], points[1], p); })) {
        std::cerr << "Все точки коллинеарны, триангуляция невозможна!" << std::endl;
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

    // Основной алгоритм
    for (const auto& point : points) {
        std::vector<Triangle> badTriangles;
        std::copy_if(triangles.begin(), triangles.end(), std::back_inserter(badTriangles),
                     [&](const Triangle& tri) { return isPointInCircumcircle(point, tri); });

        std::vector<Edge> polygonEdges;
        for (const auto& tri : badTriangles) {
            for (const auto& edge : { Edge{tri.a, tri.b}, Edge{tri.b, tri.c}, Edge{tri.c, tri.a} }) {
                bool isShared = std::any_of(badTriangles.begin(), badTriangles.end(),
                    [&](const Triangle& other) { 
                        return !(tri == other) && hasEdge(other, edge); 
                    });
                
                if (!isShared) {
                    polygonEdges.push_back(edge);
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
        }
    }

    // Удаляем треугольники, связанные с супер-треугольником
    triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
        [&](const Triangle& t) {
            return t.a == p1 || t.a == p2 || t.a == p3 ||
                   t.b == p1 || t.b == p2 || t.b == p3 ||
                   t.c == p1 || t.c == p2 || t.c == p3;
        }), triangles.end());

    return triangles;
}

bool isPointInCircumcircle(const PointD& p, const Triangle& tri) {
    double ax = tri.a.x - p.x, ay = tri.a.y - p.y;
    double bx = tri.b.x - p.x, by = tri.b.y - p.y;
    double cx = tri.c.x - p.x, cy = tri.c.y - p.y;
    
    double det = ax * (by * (cx*cx + cy*cy) - cy * (bx*bx + by*by)) 
               - ay * (bx * (cx*cx + cy*cy) - cx * (bx*bx + by*by)) 
               + (ax*ax + ay*ay) * (bx*cy - by*cx);
    return det > 0;
}

bool hasEdge(const Triangle& tri, const Edge& edge) {
    return Edge(tri.a, tri.b) == edge 
        || Edge(tri.b, tri.c) == edge 
        || Edge(tri.c, tri.a) == edge;
}
       
    int incrementAndCollect(std::vector<std::vector<double>>& componenta, std::vector<std::vector<double>>& CopyPole, int x, int y, int i, int& pixelCount) {
    if (x < 1 || y < 1 || x > (int)componenta[0].size() - 2 || 
        y > (int)componenta.size() - 2 || CopyPole[y][x] < 250) {
        return -1;
    }

    if (CopyPole[y][x] >= 255 && CopyPole[y][x] <= 255) {
        CopyPole[y][x] = 0; // Пометить как посещенное
        pixelCount++; // Увеличиваем счетчик пикселей
        componenta[y][x] = 255;
        
        // Рекурсивный вызов для соседних пикселей
        incrementAndCollect(componenta, CopyPole, x + 1, y, i + 1, pixelCount);
        incrementAndCollect(componenta, CopyPole, x - 1, y, i + 1, pixelCount);
        incrementAndCollect(componenta, CopyPole, x, y + 1, i + 1, pixelCount);
        incrementAndCollect(componenta, CopyPole, x, y - 1, i + 1, pixelCount);
    }
    
    return pixelCount;
}
       void bin(std::vector<std::vector<double>>& CopyPole, int slice, std::unique_ptr<Pole>& p, ThresholdMode mode) {
    if (p == nullptr) {
        std::cerr << "Pole not initialized!" << std::endl;
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
}
    
    void wave(int noisy, std::vector<Component>& componenti, std::vector<std::vector<double>>& CopyPole, std::unique_ptr<Pole>& p, Logger& logger, bool loggingEnabled) {
    if (p == nullptr) {
        std::cerr << "Pole not initialized!" << std::endl;
        return;
    }

    int rows = p->field.size();
    int cols = (rows > 0) ? p->field[0].size() : 0;
    std::vector<Component> noiseComponents; // Для хранения шумовых компонент

    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {
            if (CopyPole[y][x] <= 255 && CopyPole[y][x] >= 255) {
                std::vector<std::vector<double>> componentData(rows, std::vector<double>(cols, 0));
                int pixelCount = 0;
                incrementAndCollect(componentData, CopyPole, x, y, 0, pixelCount);

                if (pixelCount >= noisy) {
                    Component component(componentData, pixelCount);
                    componenti.push_back(component);
                    
                    // Логируем метаданные значимой компоненты
               logger.logMessage( "Significant component: pixels=" + std::to_string(pixelCount) +
                                  ", center=(" + std::to_string(component.center_x) + "," + std::to_string(component.center_y) +
                                  "), size=(" + std::to_string(component.max_x - component.min_x) + 
                                  "x" + std::to_string(component.max_y - component.min_y) + ")",
                                   loggingEnabled);
                } else {
                    Component noiseComponent(componentData, pixelCount);
                    noiseComponents.push_back(noiseComponent);
                    
                    // Логируем метаданные шумовой компоненты  
                     logger.logMessage( "Noise component: pixels=" + std::to_string(pixelCount) +
                                        ", center=(" + std::to_string(noiseComponent.center_x) + "," + std::to_string(noiseComponent.center_y) +
                                        "), size=(" + std::to_string(noiseComponent.max_x - noiseComponent.min_x) + 
                                        "x" + std::to_string(noiseComponent.max_y - noiseComponent.min_y) + ")",
                                        loggingEnabled);
                    
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
    
        logger.logMessage( "Wave processing complete. Significant components: " + std::to_string(componenti.size()) +
                           ", Noise components: " + std::to_string(noiseComponents.size()),
                           loggingEnabled);
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
    TerrainGrid terrainGrid;
    std::unique_ptr<Pole> p= nullptr; // Pointer to Pole
    std::unique_ptr<KMeans> kMeans = nullptr; // Используем умный указатель
    std::vector<std::vector<double>> kMeansData;
    std::vector<Triangle> lastTriangulation; // Результат триангуляции
    std::vector<VoronoiEdge> voronoiEdges;  // Храним ребра здесь
    PathFinder pathFinder{config};
    VoronoiDiagram voronoi;
    std::vector<std::vector<PointD>> paths;
    
    std::vector<PointD> clusterCenters;      // Центры кластеров
    std::vector<PointD> path;                // Найденный путь 
    
    Control(Config& cfg, Logger& log) : config(cfg), loggercontrol(log), terrainGrid(cfg.fieldWidth, cfg.fieldHeight),  // Инициализация с параметрами
                    pathFinder(cfg) {
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

    if (params.s == "PlotMetedata") {
    gnuplotInterface.plotBinaryWithComponents(CopyPole, componenti, params.filename);//визуализация метаданных 
    loggercontrol.logMessage("PlotMetedata used", b);
    }
    
    if (params.s == "PlotVoronoi") {
    gnuplotInterface.plotVoronoi(p, voronoiEdges, clusterCenters, params.filename);
    loggercontrol.logMessage("PlotVoronoi used", b);
    }
    
    if (params.s == "PlotDelaunay") {
    gnuplotInterface.plotDelaunay(lastTriangulation, p, params.filename);
    loggercontrol.logMessage("PlotDelaunay used", b);
    }
    
    if (params.s == "PlotPath") {
    gnuplotInterface.plotPath(path, p, params.filename, params); // Вызов метода визуализации
    loggercontrol.logMessage("PlotPath used", b);
    }

    if (params.s == "bmp_write") {
    if (params.bmp_mode == BmpWriteMode::Full) {
        bmpHandler.bmp_write(p->field, params.filename);
    } else {
        bmpHandler.bmp_write(CopyPole, params.filename);
    }
    loggercontrol.logMessage("bmp_write used, mode: " + std::string((params.bmp_mode == BmpWriteMode::Full ? "full" : "binary")), b);
}

    if (params.s == "bmp_read") {
    bmpHandler.bmp_read(gaussBuilder, params.filename, CopyPole, p);
    loggercontrol.logMessage("bmp_read used, filename: " + params.filename, b);
    }

    if (params.s == "bin") {
    componentCalculator.bin(CopyPole, params.slice, p, params.bin_mode);
    loggercontrol.logMessage("bin used, slice=" + std::to_string(params.slice), b);
    }
    
    if (params.s == "wave") {
    componentCalculator.wave(params.noisy, componenti, CopyPole, p, loggercontrol, b);
    loggercontrol.logMessage("Wave used with noisy threshold: " + std::to_string(params.noisy), b);
    loggercontrol.logMessage("Component amount = " + std::to_string(componenti.size()), b);
    copier.removeNoise(CopyPole, componenti); //копия без шума 
    }
     
    if (params.s == "k_means") {
            if (params.k <= 0 || p == nullptr) return;
            
            copier.removeNoise(CopyPole, componenti);//копия без шума
            prepareKMeansData(CopyPole);
            if (kMeansData.empty()) return;

            auto result = kMeans->cluster(kMeansData, params.k);
            applyClusterResults(result, CopyPole);
            
    // Логирование параметров и центров кластеров
    std::string logMessage = "k_means clustering: k=" + std::to_string(params.k) + 
                            ", centers=[";
    for (const auto& center : result.centers) {
        logMessage += "(" + std::to_string(center[0]) + "," + std::to_string(center[1]) + "),";
    }
    if (!result.centers.empty()) logMessage.pop_back(); // Удаляем последнюю запятую
    logMessage += "]";
    loggercontrol.logMessage(logMessage, b);
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
            // Логирование параметров и центров кластеров одной строкой
    std::string logMessage = "k_means_kern clustering: k=" + std::to_string(params.k) + 
                            ", kk=" + std::to_string(params.kk) + 
                            ", centers=[";
    for (const auto& center : result.centers) {
        logMessage += "(" + std::to_string(center[0]) + "," + std::to_string(center[1]) + "),";
    }
    if (!result.centers.empty()) logMessage.pop_back(); // Удаляем последнюю запятую
    logMessage += "]";
    loggercontrol.logMessage(logMessage, b);
    }
        
    if (params.s == "triangulate") {
        clusterCenters = getClusterCenters();     
    lastTriangulation = componentCalculator.bowyerWatson(clusterCenters);
    voronoi.buildFromDelaunay(lastTriangulation, pathFinder, p, voronoiEdges);
}

        // Новая команда для поиска пути
   if (params.s == "find_path") {
       if (p) {
    terrainGrid.calculateSlopes(*p);
} else {
    loggercontrol.logMessage("Error: Pole is not initialized!", b);
}
        PointD start(params.pointA_x, params.pointA_y); // Добавлено
    PointD goal(params.pointB_x, params.pointB_y);  // Добавлено
    if (lastTriangulation.empty()) {
        std::cerr << "Error: Triangulation not performed yet! Run 'triangulate' first.\n";
        return;
    }
    copier.removeNoise(CopyPole, componenti);//копия без шума
    
    path = pathFinder.findPathAStar(start, goal, lastTriangulation, terrainGrid, CopyPole, loggercontrol, b);
    
    if (path.empty()) {
        std::cerr << "Путь не найден! Проверьте триангуляцию и параметры." << std::endl;
    } else paths.push_back(path);
    
    }
    
    if (params.s == "Plot3DPath") {
    PointD start(params.pointA_x, params.pointA_y);
    PointD end(params.pointB_x, params.pointB_y);
    gnuplotInterface.plot3DPath(path, p, params.filename, start, end);
    loggercontrol.logMessage("Plot3DPath used", b);
    }
    
    if (params.s == "plotInteractive3DPath") {
    PointD start(params.pointA_x, params.pointA_y);
    PointD end(params.pointB_x, params.pointB_y);
    gnuplotInterface.plotInteractive3DPath(path, p, start, end);
    loggercontrol.logMessage("plotInteractive3DPath used", b);
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
        std::string commandfilename;
        std::ifstream file;
        std::string line;
        std::cout << "Hello, dear user, this program builds Gaussians.\nEnter commands from a text file (PRESS 0) or from the keyboard (PRESS 1)?" << std::endl;
        std::cin >> a;
        loggerinterface.logMessage("User chose input method: " + std::to_string(a), b);
        if (a == 0) {
            std::cout << "You will enter commands from a text file.\nEnter filename:" << std::endl;
            std::cin >> commandfilename;
            loggerinterface.logMessage("Reading commands from file: " + commandfilename, b);
            file.open(commandfilename);
            
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
- Поиск пути с ограничениями по углу наклона тележки и ее радиуса
  

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


## 🛠 Команды управления

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
    // Установка значений по умолчанию
    params.x = config.defaultX;
    params.y = config.defaultY;
    params.sx = config.defaultSx;
    params.sy = config.defaultSy;
    params.h = config.defaultH;

    // Чтение строки из файла
    if (std::getline(file, line)) {
        std::istringstream iss(line);
        
        // Чтение параметров команды
        iss >> params.x;  // Если не прочитается, останется defaultX
        iss >> params.y;  // Аналогично для остальных параметров
        iss >> params.sx; //...
        iss >> params.sy; //...
        iss >> params.h; //...
    }

    // Логирование и вызов обработчика
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
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.filename = config.defaultGnuplot;
                    
                    // Чтение параметров команды
                    iss >> params.filename;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Called gnuplot:" + params.filename, b);
                } else if (params.s == "PlotMetedata") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.filename = config.defaultPlotMetedata;
                    
                    // Чтение параметров команды
                    iss >> params.filename;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Called PlotMetedata:" + params.filename, b);
                } else if (params.s == "PlotVoronoi") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.filename = config.defaultPlotVoronoi;
                    
                    // Чтение параметров команды
                    iss >> params.filename;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Called PlotVoronoi:" + params.filename, b);
                } else if (params.s == "PlotDelaunay") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.filename = config.defaultPlotDelaunay;
                    
                    // Чтение параметров команды
                    iss >> params.filename;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Called PlotDelaunay:" + params.filename, b);
                } else if (params.s == "PlotPath") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.filename = config.defaultPlotPath;
                    
                    // Чтение параметров команды
                    iss >> params.filename;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Called PlotPath:" + params.filename, b);
                } else if (params.s == "bmp_write") {
                    std::string modeWrite;
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.filename = config.defaultWrite;
                    modeWrite = config.defaultWriteModeImage;
                    
                    // Чтение параметров команды
                    iss >> params.filename;
                    iss >> modeWrite;
                    
                    if (modeWrite == "Full") {
                params.bmp_mode = BmpWriteMode::Full;
            } else params.bmp_mode = BmpWriteMode::Binary;
                
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Created BMP file:" + params.filename + ", Mode image:" + modeWrite, b);
                } else if (params.s == "bmp_read") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.filename = config.defaultRead;
                    
                    // Чтение параметров команды
                    iss >> params.filename;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Read BMP file: " + params.filename, b);
                } else if (params.s == "bin") {
                    std::string modeBin;
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.slice = config.defaultSlice;
                    modeBin = config.defaultBinMode;
                    
                    // Чтение параметров команды
                    iss >> params.slice;
                    iss >> modeBin;
                    
            if (modeBin == "Peaks") {
                params.bin_mode = ThresholdMode::Peaks;
            } else if (modeBin == "Valleys") {
                params.bin_mode = ThresholdMode::Valleys;
            } else {
                params.bin_mode = ThresholdMode::All;
            }
        
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Slice applied: slice=" + std::to_string(params.slice) + ", Bin Mode:" + modeBin, b);
                } else if (params.s == "wave") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.noisy = config.defaultNoisy;
                    
                    // Чтение параметров команды
                    iss >> params.noisy;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Wave will be used", b);
                    loggerinterface.logMessage("Component amount = " + std::to_string(c.componenti.size()), b);
                } else if (params.s == "k_means") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.k = config.defaultKlaster;
                    
                    // Чтение параметров команды
                    iss >> params.k;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Cluster count", b);
                } else if (params.s == "k_means_kern") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.kk = config.defaultKlasterKern;
                    
                    // Чтение параметров команды
                    iss >> params.kk;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Cluster kern worked", b);
                } else if (params.s == "triangulate") {
                     loggerinterface.logMessage("triangulate worke", b);
                     c.Dispetcher(params);
                 
                } else if (params.s == "find_path") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                       
                    // Установка значений по умолчанию
                    params.pointA_x = config.defaultpointA_x;
                    params.pointA_y = config.defaultpointA_y;
                    params.pointB_x = config.defaultpointB_x;
                    params.pointB_y = config.defaultpointB_y;
                    // Чтение параметров команды
                    iss >> params.pointA_x;  // Если не прочитается, останется defaultpointA_x
                    iss >> params.pointA_y;  // Аналогично для остальных параметров
                    iss >> params.pointB_x; //...
                    iss >> params.pointB_y; 
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("find_path worked", b);
                } else if (params.s == "Plot3DPath") {
                    std::getline(file, line);
                    std::istringstream iss(line);
                    
                    // Установка значений по умолчанию
                    params.filename = config.defaultPlot3DPath;
                    
                    // Чтение параметров команды
                    iss >> params.filename;
                    
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Called Plot3DPath:" + params.filename, b);
                } else if (params.s == "plotInteractive3DPath") {
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Worked plotInteractive3DPath", b);
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
# Terrain Navigation System

Программа для анализа рельефа местности, построения триангуляции Делоне, диаграмм Вороного и поиска оптимальных маршрутов с учетом препятствий.

## 📌 Основные функции
- Генерация/загрузка карты высот (формат BMP и GNUPLOT)
- Кластеризация объектов методом k-means
- Триангуляция Делоне с учетом высот
- Построение диаграммы Вороного
- Поиск пути с ограничениями по углу наклона тележки и ее радиуса
  

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


## 🛠 Команды управления

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
    std::getline(std::cin, line);
    std::istringstream iss(line);
    
    // Инициализируем переменные значениями по умолчанию
    params.x = config.defaultX;
    params.y = config.defaultY;
    params.sx = config.defaultSx;
    params.sy = config.defaultSy;
    params.h = config.defaultH;
    
        
    // Читаем введенные данные
    if (!(iss >> params.x)) {
        std::cout << "The default value for x is used: " << config.defaultX << std::endl;
    }
    if (!(iss >> params.y)) {
        std::cout << "The default value for y is used: " << config.defaultY << std::endl;
    }
    if (!(iss >> params.sx)) {
        std::cout << "The default value for sx is used: " << config.defaultSx << std::endl;
    }
    if (!(iss >> params.sy)) {
        std::cout << "The default value for sy is used: " << config.defaultSy << std::endl;
    }
    if (!(iss >> params.h)) {
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
                    std::cout << "Enter the filename(full path) to draw:" << std::endl;
                    std::cin >> params.filename;//имя для файла
                    c.Dispetcher(params);
                    std::cout << "Called gnuplot:" + params.filename << std::endl;
                    loggerinterface.logMessage("Called gnuplot:" + params.filename, b);
                } if (params.s == "PlotMetedata") {
                    std::cout << "Enter the filename(full path) to draw:" << std::endl;
                    std::cin >> params.filename;//имя для файла
                    c.Dispetcher(params);
                    std::cout << "Called PlotMetedata:" + params.filename << std::endl;
                    loggerinterface.logMessage("Called PlotMetedata:" + params.filename, b);
                } if (params.s == "PlotVoronoi") {
                    std::cout << "Enter the filename(full path) to draw:" << std::endl;
                    std::cin >> params.filename;//имя для файла
                    c.Dispetcher(params);
                    std::cout << "Called PlotVoronoi:" + params.filename << std::endl;
                    loggerinterface.logMessage("Called PlotVoronoi:" + params.filename, b);
                } if (params.s == "PlotDelaunay") {
                    std::cout << "Enter the filename(full path) to draw:" << std::endl;
                    std::cin >> params.filename;//имя для файла
                    c.Dispetcher(params);
                    std::cout << "Called PlotDelaunay:" + params.filename << std::endl;
                    loggerinterface.logMessage("Called PlotDelaunay:" + params.filename, b);
                } if (params.s == "PlotPath") {
                    std::cout << "Enter the filename(full path) to draw:" << std::endl;
                    std::cin >> params.filename;//имя для файла
                    c.Dispetcher(params);
                    std::cout << "Called PlotPath:" + params.filename << std::endl;
                    loggerinterface.logMessage("Called PlotPath:" + params.filename, b);
                } if (params.s == "bmp_write") {
                    std::string modeWrite;
                    
                    std::cout << "Enter the filename(full path) to draw:" << std::endl;
                    std::cin >> params.filename;//имя для файла
                    std::cout << "Image is binary (Binary) or not (Full):" << std::endl;
                    std::cin >> modeWrite;
                    
                    if (modeWrite == "Full") {
                params.bmp_mode = BmpWriteMode::Full;
            } else params.bmp_mode = BmpWriteMode::Binary;
            
                    c.Dispetcher(params);
                    std::cout << "Created BMP file:" + params.filename + ", bmp_mode:" + modeWrite << std::endl;
                    loggerinterface.logMessage("Created BMP file:" + params.filename, b);
                } if (params.s == "bmp_read") {
                    std::cout << "Enter the filename(full path) to read:" << std::endl;
                    std::cin >> params.filename;
                    c.Dispetcher(params);
                    std::cout << "Read BMP file: " + params.filename << std::endl;
                    loggerinterface.logMessage("Read BMP file: " + params.filename, b);
                }  if (params.s == "end") {
                    std::cout << "Ending the program" << std::endl;
                    loggerinterface.logMessage("Ending the program.", b);
                    break;
                }
                   if (params.s == "bin") {
                     std::string modeBin;
                     
                     std::cout << "Enter slice level:" << std::endl;
                     std::cin >> params.slice;
                     std::cout << "Enter Bin Mode (Peaks, Valleys, All): " << std::endl;
                     std::cin >> modeBin;
                     
            if (modeBin == "Peaks") {
                params.bin_mode = ThresholdMode::Peaks;
            } else if (modeBin == "Valleys") {
                params.bin_mode = ThresholdMode::Valleys;
            } else {
                params.bin_mode = ThresholdMode::All;
            }
            
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
                        std::cout << "Enter amount pointA_x:" << std::endl;
                        std::cin >> params.pointA_x;
                        std::cout << "Enter amount pointA_y:" << std::endl;
                        std::cin >> params.pointA_y;
                        std::cout << "Enter amount pointB_x:" << std::endl;
                        std::cin >> params.pointB_x;
                        std::cout << "Enter amount pointB_y:" << std::endl;
                        std::cin >> params.pointB_y;
                        
                        c.Dispetcher(params);
                        std::cout << "Path search completed" << std::endl;
                        loggerinterface.logMessage("find_path worked", b);
                }
                    if (params.s == "Plot3DPath") {
                        std::cout << "Enter the filename(full path) to draw:" << std::endl;
                        std::cin >> params.filename;//имя для файла
                        c.Dispetcher(params);
                        std::cout << "Called Plot3DPath:" + params.filename << std::endl;
                        loggerinterface.logMessage("Called Plot3DPath:" + params.filename, b);
                }
                    if (params.s == "plotInteractive3DPath") {
                    c.Dispetcher(params);
                    loggerinterface.logMessage("Worked plotInteractive3DPath", b);
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
