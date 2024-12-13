
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <algorithm> // Для std::sample
#include <random>    // Для std::mt19937 и std::random_device

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
    
//Классы для построения гаусса
class Gaus {
public:
    double h, x0, y0, sx, sy;

    Gaus(double h, double x0, double y0, double sx, double sy)
        : h(h), x0(x0), y0(y0), sx(sx), sy(sy) {}
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
    
    Component(const std::vector<std::vector<double>>& inputComponenta) : componenta(inputComponenta) {}
    
    Component(int A, int B) {
        componenta.resize(A, std::vector<double>(B, 0));
    }
};

class Copier {
public:
    // Метод для копирования не нулевых значений из списка компонентов в Copypole
    void copyNonZeroValues(std::vector<std::vector<double>>& copypole, const std::vector<Component>& components) {
        // Проверяем, что размеры Copypole и компонентов совпадают
        if (components.empty() || copypole.empty() || 
            copypole.size() != components[0].componenta.size() || 
            copypole[0].size() != components[0].componenta[0].size()) {
            std::cerr << "Размеры Copypole и компонентов не совпадают!" << std::endl;
            return;
        }

        // Проходим по каждому компоненту и заполняем copypole
        for (const auto& component : components) {
            for (size_t i = 0; i < component.componenta.size(); ++i) {
                for (size_t j = 0; j < component.componenta[i].size(); ++j) {
                    if (component.componenta[i][j] != 0.0) { // Проверяем на не нулевое значение
                        copypole[i][j] = component.componenta[i][j]; // Записываем значение в copypole
                    }
                    if (copypole[i][j] < 255) {
                        copypole[i][j] = 0;
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
   
        double value; 
        for (const auto& g : gaussi) {
            for (long unsigned int x = 0; x < p->field[0].size(); ++x) {
                for (long unsigned int y = 0; y < p->field.size(); ++y) {
                    value = g.h * exp(-((pow((x - g.x0) / g.sx, 2)) + (pow((y - g.y0) / g.sy, 2))) / 2); 
                    p->field[y][x] += value; 
                    if (p->field[y][x] > 255) {
                        p->field[y][x] = 255;
                    }
                    if (p->field[y][x] < 0) {
                        p->field[y][x] = 0;
                    }
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
            unsigned char color = 255 - static_cast<unsigned char>(pixelMatrix[y][x]); // Color
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
            double value = 255 - color; // Цвет в высоту
            pixelMatrix[y][x] = value; // Обновляем матрицу значений
            p->field[y][x] = value; // Обновляем поле в Pole
            //Значения в gaussi некорректные
        }
        bmpFile.ignore((4 - (width * 3) % 4) % 4); // Пропускаем паддинг
    }
    bmpFile.close();
  }
};
   
class GnuplotInterface {
   public:
   
      void gnuplot(std::unique_ptr<Pole>& p) {
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
    
    // Устанавливаем диапазон z от 0 до 255
    fprintf(gnuplotPipe, "set zrange [0:255]\n"); 
    fprintf(gnuplotPipe, "set xrange [0:%d]\n", cols - 1); 
    fprintf(gnuplotPipe, "set yrange [0:%d]\n", rows - 1); 
    fprintf(gnuplotPipe, "set terminal png\n"); 
    fprintf(gnuplotPipe, "set output 'landscape.png'\n"); 
    fprintf(gnuplotPipe, "splot '-' with lines\n"); 

    // Используем буфер для хранения данных
    std::ostringstream dataStream;
    for (int y = 0; y < rows; ++y) { 
        for (int x = 0; x < cols; ++x) { 
            dataStream << x << " " << y << " " << p->field[y][x] << "\n";
        }
    }

    // Отправляем все данные за один раз
    fprintf(gnuplotPipe, "%s", dataStream.str().c_str());

    fprintf(gnuplotPipe, "\n"); // Одна пустая строка для завершения ввода данных
    pclose(gnuplotPipe);
}



};
   
class ComponentCalculator {
   public:
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
       void bin(BmpHandler& bmpHandler, std::vector<std::vector<double>> &CopyPole, int slise, std::unique_ptr<Pole>& p) {
       if (p == nullptr) {
       std::cerr << "Pole not initialized!" << std::endl;
       return;
   }
       CopyPole = p->field;//Копия
        
            for (int x = 0; x < (int) p->field[0].size(); ++x) {
                for (int y = 0; y < (int) p->field.size(); ++y) {
                   CopyPole[y][x] = p->field[y][x] > slise ? 255 : 0;
                } 
            }
        
        bmpHandler.bmp_write(CopyPole, "slise.bmp");
        //loggerinterface.logMessage("Created BMP file.", b);
    }
    
    void wave(std::vector<Component>& componenti, std::vector<std::vector<double>> &CopyPole, std::unique_ptr<Pole>& p) {
    if (p == nullptr) {
       std::cerr << "Pole not initialized!" << std::endl;
       return;
   }
    
            for (int y = 0; y < (int) p->field.size(); ++y) {
                for (int x = 0; x < (int) p->field[y].size(); ++x) {
                    if (CopyPole[y][x] <= 255 && CopyPole[y][x] >= 255) {
                        count = 0;
                        Component Componenta(p->field.size(), p->field[0].size());
                        if (incrementAndCollect(Componenta.componenta, CopyPole, x, y, 1) >= 10) {
                         componenti.emplace_back(Componenta);
                       }
                   }
               }
           }

        //loggerinterface.logMessage("Wave used, amount component = " + std::to_string(componenti.size()), b);
   }
};

class KMeans {
private:

public:
    std::vector<int> labels;
    std::vector<std::vector<double>> centers; // Центры кластеров

    // Метод для инициализации центров случайными пикселями
    void initializeCenters(const std::vector<std::vector<double>>& data, int k) {
        std::vector<int> indices(data.size());
        std::iota(indices.begin(), indices.end(), 0); // Заполняем индексы от 0 до data.size() - 1

        std::shuffle(indices.begin(), indices.end(), std::mt19937{ std::random_device{}() }); // Перемешиваем индексы

        centers.clear();
        for (int i = 0; i < k; ++i) {
            if (i < indices.size()) {
                centers.push_back({data[indices[i]][0], data[indices[i]][1]}); // Добавляем только первые два элемента
            }
        }
    }

    // Метод для вычисления расстояния между двумя точками
    double distance(const std::vector<double>& a, const std::vector<double>& b) {
        double sum = 0.0;
        for (size_t i = 0; i < 2; ++i) {
            sum += (a[i] - b[i]) * (a[i] - b[i]);
        }
        return std::sqrt(sum);
    }

    // Метод для выполнения кластеризации
    void cluster(const std::vector<std::vector<double>>& data, std::vector<std::vector<double>>& center, int k) {
        labels.assign(data.size(), -1); // Массив меток для данных
        bool changed;

        do {
            changed = false;

            // Присваиваем метки кластерам
            for (size_t i = 0; i < data.size(); ++i) {
                 if (data[i].size() > 2 && data[i][2] <= 255 && data[i][2] >= 255) {
                int nearestCenter = -1;
                double minDistance = std::numeric_limits<double>::max();

                for (size_t j = 0; j < center.size(); ++j) {
                    double dist = distance(data[i], center[j]);
                    if (dist < minDistance) {
                        minDistance = dist;
                        nearestCenter = j;
                    }
                }

                if (labels[i] != nearestCenter) {
                    labels[i] = nearestCenter;
                    changed = true; // Если метка изменилась, помечаем, что произошло изменение
                }
            }
        }

            // Перерасчет центров
            std::vector<std::vector<double>> newCenters(k, std::vector<double>(2, 0.0));
            std::vector<int> counts(k, 0);

            for (size_t i = 0; i < data.size(); ++i) {
                if (labels[i] != -1) {
                    for (size_t j = 0; j < 2; ++j) {
                        newCenters[labels[i]][j] += data[i][j];
                    }
                    counts[labels[i]]++;
                }
            }

            for (size_t j = 0; j < k; ++j) {
                if (counts[j] > 0) {
                    for (size_t m = 0; m < newCenters[j].size(); ++m) {
                        newCenters[j][m] /= counts[j]; // Среднее значение
                    }
                }
            }

            center = newCenters; // Обновляем центры

        } while (changed); // Повторяем, пока изменения происходят
    }
    
    // Новый метод с ядрами
    void kmeansWithKernels(const std::vector<std::vector<double>>& data, int k, int p) {
    // Создаем вектор cent с размером как у data
std::vector<std::vector<double>> cent(data.size(), std::vector<double>(3, 0)); // Инициализируем координаты и значение

    initializeCenters(data, k*p);
    // Предполагаем, что centers - это вектор, где каждая строка содержит координаты пикселей (например, {i, j})
for (const auto& center : centers) {
    // Определяем координаты пикселя
    size_t i = static_cast<size_t>(center[0]); // x координаты
    size_t j = static_cast<size_t>(center[1]); // y координаты

    // Проверяем, находятся ли координаты в пределах допустимого диапазона
    if (i < data.size() && j < data[0].size()) {
        // Присваиваем значение 255 в cent для пикселей из centers
        cent[i][2] = 255; // значение для пикселей из centers
    }
}

for (size_t i = 0; i < data.size(); ++i) {
    for (size_t j = 0; j < data[i].size(); ++j) {
            cent[i][0] = static_cast<double>(i); // x координаты
            cent[i][1] = static_cast<double>(j); // y координаты
        }
    }
    
    initializeCenters(cent, k);
    
    cluster(cent, centers, k);

    // После обновления центров можно снова применить обычный K-Means
    cluster(data, centers, k);
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
    Copier copy;
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
     
    Control(Config& cfg, Logger& log) : config(cfg), loggercontrol(log) {
        if (config.loggingControlEnabled) {
            
                loggercontrol.logMessage("Logging Control is enabled.", b);
                std::cout << "Logging Control is enabled." << std::endl;
            } else {
            loggercontrol.logMessage("Logging Control is disabled.", b);
            std::cout << "Logging Control is disabled." << std::endl;
            b = false;
        }      
    }
    
void Dispetcher(const std::string& command, int A, int B, double h, double x0, double y0, double sigma_x, double sigma_y, 
std::vector<std::vector<double>>& pixelMatrix, const std::string& filename, int slice, int klaster, int pixelsInKernel) {
   
    if (command == "init") {
    gaussBuilder.init(A, B, p);
    loggercontrol.logMessage("init used: " + std::to_string(A) + " x " + std::to_string(B), b);
    }

    if (command == "g") {
    gaussBuilder.addgauss(h, x0, y0, sigma_x, sigma_y, gaussi);
    loggercontrol.logMessage("g used: x=" + std::to_string(x0) +
                    ", y=" + std::to_string(y0) + 
                    ", sx=" + std::to_string(sigma_x) + 
                    ", sy=" + std::to_string(sigma_y) + 
                    ", h=" + std::to_string(h), b);
    }
    
    if (command == "generate") {
    gaussBuilder.generate(p, gaussi);
    loggercontrol.logMessage("generate used", b);
    }
    
    if (command == "gnuplot") {
    gnuplotInterface.gnuplot(p);
    loggercontrol.logMessage("gnuplot used", b);
    }

    if (command == "bmp_write") {
    bmpHandler.bmp_write(pixelMatrix, filename);
    loggercontrol.logMessage("bmp_write used", b);
    }

    if (command == "bmp_read") {
    bmpHandler.bmp_read(gaussBuilder, filename, pixelMatrix, p);
    loggercontrol.logMessage("bmp_read used, filename: " + filename, b);
    }

    if (command == "bin") {
    componentCalculator.bin(bmpHandler, pixelMatrix, slice, p);
    loggercontrol.logMessage("bin used, slice=" + std::to_string(slice), b);
    componentCalculator.wave(componenti, pixelMatrix, p);
    loggercontrol.logMessage("wave used", b);
    loggercontrol.logMessage("Component amount = " + std::to_string(componenti.size()), b);
    copy.copyNonZeroValues(pixelMatrix, componenti);//копия без шума
    copy.copyNonZeroValues(p->field, componenti);//поле без шума
    }
    
    if (command == "k_means") {
            if (klaster <= 0) {
                loggercontrol.logMessage("Error: klaster must be greater than 0", b);
                return; // Проверка на корректность klaster
            }
            
            if (pixelMatrix.empty()) {
        std::cerr << "Error: pixelMatrix is not initialized!" << std::endl;
        return;
    }


            // Генерация уникальных цветов
            colors = ColorGenerator::generateColors(klaster);

             // Преобразуем pixelMatrix в формат, подходящий для K-средних
        for (size_t i = 0; i < pixelMatrix.size(); ++i) {
            for (size_t j = 0; j < pixelMatrix[i].size(); ++j) {
                    kMeansData.push_back({static_cast<double>(i), static_cast<double>(j), pixelMatrix[i][j]});
            }
        }

        // Если нет пикселей для кластеризации, выходим
        if (kMeansData.empty()) {
            loggercontrol.logMessage("Error: No colored pixels found for clustering", b);
            return;
        }

        kMeans = std::make_unique<KMeans>(); // Используем умный указатель
        loggercontrol.logMessage("k_means used", b);
        kMeans->initializeCenters(kMeansData, klaster); // Инициализируем центры
        kMeans->cluster(kMeansData, kMeans->centers, klaster); // Выполняем кластеризацию

        // Обновляем pixelMatrix, используя метки кластеров
for (size_t i = 0; i < kMeansData.size(); ++i) {
    int clusterLabel = kMeans->labels[i]; // Получаем метку кластера из массива labels
    if (clusterLabel != -1) { // Если метка действительна
        auto color = colors[clusterLabel]; // Получаем цвет для кластера
        pixelMatrix[static_cast<int>(kMeansData[i][0])][static_cast<int>(kMeansData[i][1])] = (color[0] + color[1] + color[2]) / 3.0; // Сохраняем среднее значение в pixelMatrix
    } else {
        pixelMatrix[static_cast<int>(kMeansData[i][0])][static_cast<int>(kMeansData[i][1])] = 0;
    }
}

        // Создаем черно-белую картинку (BMP)
        bmpHandler.bmp_write(pixelMatrix, "output_kmeans.bmp");
        }
        if (command == "k_means_kern") {
         // Генерация уникальных цветов
            colors = ColorGenerator::generateColors(klaster);
             // Преобразуем pixelMatrix в формат, подходящий для K-средних
        for (size_t i = 0; i < pixelMatrix.size(); ++i) {
            for (size_t j = 0; j < pixelMatrix[i].size(); ++j) {
                    kMeansData.push_back({static_cast<double>(i), static_cast<double>(j), pixelMatrix[i][j]});
            }
        }

        // Если нет пикселей для кластеризации, выходим
        if (kMeansData.empty()) {
            loggercontrol.logMessage("Error: No colored pixels found for clustering", b);
            return;
        }
            kMeans->kmeansWithKernels(kMeansData, klaster, pixelsInKernel);
             // Обновляем pixelMatrix, используя метки кластеров
for (size_t i = 0; i < kMeansData.size(); ++i) {
    int clusterLabel = kMeans->labels[i]; // Получаем метку кластера из массива labels
    if (clusterLabel != -1) { // Если метка действительна
        auto color = colors[clusterLabel]; // Получаем цвет для кластера
        pixelMatrix[static_cast<int>(kMeansData[i][0])][static_cast<int>(kMeansData[i][1])] = (color[0] + color[1] + color[2]) / 3.0; // Сохраняем среднее значение в pixelMatrix
    } else {
        pixelMatrix[static_cast<int>(kMeansData[i][0])][static_cast<int>(kMeansData[i][1])] = 0;
    }
}

        // Создаем черно-белую картинку (BMP)
        bmpHandler.bmp_write(pixelMatrix, "output_kmeansWithKernels.bmp");
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
    
    Interface(Config& cfg, Logger& log, Control& c) : config(cfg), loggerinterface(log), c(c){       
        if (config.loggingInterfaceEnabled) {
            
                loggerinterface.logMessage("Logging Interface is enabled.", b);
                std::cout << "Logging Interface is enabled." << std::endl;
        } else {
        loggerinterface.logMessage("Logging Interface is disabled.", b);
        std::cout << "Logging Interface is disabled." << std::endl;
        b = false;
      }        
  }
  
  // Получаем размеры поля из конфигурации
    int A = config.fieldWidth;  // Размеры из конфигурации
    int B = config.fieldHeight;
    int slice = 0;
    int k = 0;//cluster
    int p = 0;//cluster kern
    
    void print() {
        double x = 0, y = 0, sx = 0, sy = 0, h = 0;
        std::string s;
        bool a;
        std::string filename;
        std::ifstream file;
        std::cout << "Hello, dear user, this program builds Gaussians.\nEnter commands from a text file (PRESS 0) or from the keyboard (PRESS 1)?" << std::endl;
        std::cin >> a;
        loggerinterface.logMessage("User chose input method: " + std::to_string(a), b);
        if (a == 0) {
            std::cout << "You will enter commands from a text file.\nEnter filename:" << std::endl;
            std::cin >> filename;
            loggerinterface.logMessage("Reading commands from file: " + filename, b);
            file.open(filename);
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
            while (file >> s) {
                loggerinterface.logMessage("Received command: " + s, b);
         if (s == "help") {
            // Открываем файл help.txt для записи
            std::ofstream helpFile("help.txt");
            if (helpFile.is_open()) {
                helpFile << R"(
Компиляция программы:
Сначала: g++ -std=c++17 Fail.cpp
Потом: ./a.out

Доступные команды:
help - создание файла с пояснением команд
init - инициализация поля. Первая команда!!!
g (x,y,sx,sy,h) - создает гаусс с 5 параметрами
generate - складывает гауссы
gnuplot - рисует картинку в gnuplot
bmp_write - создает черно-белую серую картинку bmp
bmp_read (FILENAME) - чтение bmp файла и инициализация поля новыми размерами
bin (INTEGER NUMBER) - срез: все, что выше черное, все что ниже белое, считается число компонентов, каждая компонента запомнена
k_means (k) - выделение k кластеров
k_means_kern (p) - алгоритм kmeans с ядрами размера p
end - это конец программы

Примеры:
1) Если нужно прочитать данные с файла cat.bmp

init
help
bmp_read cat.bmp
gnuplot
bmp_write
bin 127
k_means 2
k_means_kern 3
end

2) Если данные вводятся с помощью гаусов

init
help
g 50 50 20 20 200
generate
gnuplot
bmp_write
bin 127
k_means 2
k_means_kern 3
end

Да? Нет? Понятно?
)";
                helpFile.close(); // Закрываем файл
                loggerinterface.logMessage("Help file created successfully.", b);
            } else {
                loggerinterface.logMessage("Failed to create help file.", b);
            }
        }
   else if (s == "init") {
    // Проверяем, была ли уже вызвана команда init
    if (n != 0) {
        std::cout << "The init command has already been called.\nError\n";
        loggerinterface.logMessage("Error: Multiple init commands.", b);
        return;
    }
     n = 1;
    // Логируем и инициализируем поле
    loggerinterface.logMessage("Initializing field with size: " + std::to_string(A) + " x " + std::to_string(B), b);
    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
    loggerinterface.logMessage("Field initialized.", b);
} else if (n != 1) {
                     std::cout << "The init command was not use.\nError\n";
                     loggerinterface.logMessage("Error: The init command was not use.", b);
                     return;
             }
 else if (s == "g") {

                // Читаем параметры из файла
                file >> x >> y >> sx >> sy >> h;

                // Проверяем, были ли параметры указаны
                if (file.fail()) {
                    // Если не удалось прочитать, значит, используем значения по умолчанию
                    if (file.eof()) { // Желательно проверить конец файла, чтобы избежать повторного чтения
                        x = config.defaultX;
                        y = config.defaultY;
                        sx = config.defaultSx;
                        sy = config.defaultSy;
                        h = config.defaultH;
                    } else {
                        // Если произошла ошибка чтения или недостаточно параметров
                        // Нужно сбросить флаг ошибки для дальнейшего чтения
                        file.clear();

                        // Проверяем, какие параметры были прочитаны
                        if (!(file >> x)) x = config.defaultX; // Если не удалось получить x, берем по умолчанию
                        if (!(file >> y)) y = config.defaultY; // Если не удалось получить y, берем по умолчанию
                        if (!(file >> sx)) sx = config.defaultSx; // ...
                        if (!(file >> sy)) sy = config.defaultSy; // ...
                        if (!(file >> h)) h = config.defaultH; // ...
                    }
                }

                loggerinterface.logMessage("Adding Gaussian: x=" + std::to_string(x) +
                    ", y=" + std::to_string(y) + 
                    ", sx=" + std::to_string(sx) + 
                    ", sy=" + std::to_string(sy) + 
                    ", h=" + std::to_string(h), b);
                c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
            } else if (s == "generate") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                    loggerinterface.logMessage("Generated values in the field.", b);
                } else if (s == "gnuplot") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                    loggerinterface.logMessage("Called gnuplot: landscape.png", b);
                } else if (s == "bmp_write") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                    loggerinterface.logMessage("Created BMP file: output.bmp", b);
                } else if (s == "bmp_read") {
                    file >> filename; // Чтение имени файла для bmp_read
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, p);
                    loggerinterface.logMessage("Read BMP file: " + filename, b);
                } else if (s == "bin") {
                    file >> slice;
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                    loggerinterface.logMessage("Slice applied: slice=" + std::to_string(slice), b);
                    loggerinterface.logMessage("Wave will be used", b);
                    loggerinterface.logMessage("Component amount = " + std::to_string(c.componenti.size()), b);
                } else if (s == "k_means") {
                        file >> k;
                        c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                        loggerinterface.logMessage("Cluster count", b);
                } else if (s == "k_means_kern") {
                    file >> p;
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                    loggerinterface.logMessage("Cluster kern worked", b);
                }
            }
        } else {
            int n = 0;
            while (true) {
                std::cout << "Enter command (help, init, g, generate, gnuplot, bmp_write, bmp_read, bin, k_means, end):";
                std::cin >> s;
                std::cout << "\n";
                loggerinterface.logMessage("Received command: " + s, b);
        if (s == "help") {
            // Открываем файл help.txt для записи
            std::ofstream helpFile("help.txt");
            if (helpFile.is_open()) {
                helpFile << R"(
Компиляция программы:
Сначала: g++ -std=c++17 Fail.cpp
Потом: ./a.out

Доступные команды:
help - создание файла с пояснением команд
init - инициализация поля. Первая команда!!!
g (x,y,sx,sy,h) - создает гаусс с 5 параметрами
generate - складывает гауссы
gnuplot - рисует картинку в gnuplot
bmp_write - создает черно-белую серую картинку bmp
bmp_read (FILENAME) - чтение bmp файла и инициализация поля новыми размерами
bin (INTEGER NUMBER) - срез: все, что выше черное, все что ниже белое, считается число компонентов, каждая компонента запомнена
k_means (k) - выделение k кластеров
k_means_kern (p) - алгоритм kmeans с ядрами размера p
end - это конец программы

Примеры:
1) Если нужно прочитать данные с файла cat.bmp

init
help
bmp_read cat.bmp
gnuplot
bmp_write
bin 127
k_means 2
k_means_kern 3
end

2) Если данные вводятся с помощью гаусов

init
help
g 50 50 20 20 200
generate
gnuplot
bmp_write
bin 127
k_means 2
k_means_kern 3
end

Да? Нет? Понятно?
)";
                helpFile.close(); // Закрываем файл
                loggerinterface.logMessage("Help file created successfully.", b);
            } else {
                loggerinterface.logMessage("Failed to create help file.", b);
            }
        }
    if (s == "init") {
    // Проверяем, была ли уже вызвана команда init
    if (n != 0) {
        std::cout << "The init command has already started.\nError\n";
        loggerinterface.logMessage("Error: Multiple init commands.", b);
        return;
    }
    n = 1; // Устанавливаем флаг, что команда init была вызвана
    
    // Логируем и инициализируем поле
    loggerinterface.logMessage("Initializing field with size: " + std::to_string(A) + " x " + std::to_string(B), b);
    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
    loggerinterface.logMessage("Field initialized.", b);
}
  if (n != 1) {
                     std::cout << "The init command was not use.\nError\n";
                     loggerinterface.logMessage("Error: The init command was not use.", b);
                     return;
             }
             
  if (s == "g") {
    std::string input;
    std::getline(std::cin, input);

    std::istringstream inputStream(input);

    // Инициализируем переменные значениями по умолчанию
    x = config.defaultX;
    y = config.defaultY;
    sx = config.defaultSx;
    sy = config.defaultSy;
    h = config.defaultH;

    // Читаем введенные данные
    if (!(inputStream >> x)) {
        std::cout << "The default value for x is used: " << config.defaultX << std::endl;
    }
    if (!(inputStream >> y)) {
        std::cout << "The default value for y is used: " << config.defaultY << std::endl;
    }
    if (!(inputStream >> sx)) {
        std::cout << "The default value for sx is used: " << config.defaultSx << std::endl;
    }
    if (!(inputStream >> sy)) {
        std::cout << "The default value for sy is used: " << config.defaultSy << std::endl;
    }
    if (!(inputStream >> h)) {
        std::cout << "The default value for h is used: " << config.defaultH << std::endl;
    }

                loggerinterface.logMessage("Adding Gaussian: x=" + std::to_string(x) +
                    ", y=" + std::to_string(y) + 
                    ", sx=" + std::to_string(sx) + 
                    ", sy=" + std::to_string(sy) + 
                    ", h=" + std::to_string(h), b);
               c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
            }

                  if (s == "generate") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                    std::cout << "Generated values in the field." << std::endl;
                    loggerinterface.logMessage("Generated values in the field.", b);
                } if (s == "gnuplot") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                    std::cout << "Called gnuplot: landscape.png" << std::endl;
                    loggerinterface.logMessage("Called gnuplot: landscape.png", b);
                } if (s == "bmp_write") {
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                    std::cout << "Created BMP file: output.bmp" << std::endl;
                    loggerinterface.logMessage("Created BMP file: output.bmp", b);
                } if (s == "bmp_read") {
                    std::cout << "Enter the filename to read:" << std::endl;
                    std::cin >> filename;
                    c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, filename, slice, k, p);
                    std::cout << "Read BMP file: " + filename << std::endl;
                    loggerinterface.logMessage("Read BMP file: " + filename, b);
                }  if (s == "end") {
                    std::cout << "Ending the program" << std::endl;
                    loggerinterface.logMessage("Ending the program.", b);
                    break;
                }
                   if (s == "bin") {
                     std::cout << "Enter slice level:" << std::endl;
                     std::cin >> slice;
                     c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                     loggerinterface.logMessage("Slice applied: slice=" + std::to_string(slice), b);
                     std::cout << "Slice applied: slice=" << slice << std::endl;
                     loggerinterface.logMessage("Wave will be used", b);
                     loggerinterface.logMessage("Component amount = " + std::to_string(c.componenti.size()), b);
                }
                    if (s == "k_means") {
                        std::cout << "Enter amount cluster:" << std::endl;
                        std::cin >> k;
                        c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                        std::cout << "Cluster count" << std::endl;
                        loggerinterface.logMessage("Cluster count", b);
                }
                    if (s == "k_means_kern") {
                      std::cout << "Enter amount kern:" << std::endl;
                      std::cin >> p;
                      c.Dispetcher(s, A, B, h, x, y, sx, sy, c.p->field, "output.bmp", slice, k, p);
                      std::cout << "Cluster kern worked" << std::endl;
                      loggerinterface.logMessage("Cluster kern worked", b);
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
