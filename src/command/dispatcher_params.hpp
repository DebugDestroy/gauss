#pragma once
#include <string>
#include "io/bmp_handler.hpp"  // Для BmpWriteMode
#include "algorithms/components/binary.hpp"  // Для ThresholdMode

namespace command {
struct DispatcherParams {
    std::string command;                                  // Команда для выполнения (help, init, g, generate и т.д.)
    
    int fieldWidth, fieldHeight;                          // Размеры поля (ширина × высота в пикселях)
    
    // ----- addgauss -----
    double height;                                        // Высота гауссова распределения (в условных единицах)
    double centerX, centerY;                              // Координаты центра (x, y)
    double sigmaX, sigmaY;                                // Разброс по осям (σx, σy)
    
    // ----- Параметры для g_auto -----
    double xmin, xmax;                                    // Диапазон X координаты центра
    double ymin, ymax;                                    // Диапазон Y координаты центра
    double sx_min, sx_max;                                // Диапазон σx
    double sy_min, sy_max;                                // Диапазон σy
    double h_min, h_max;                                  // Диапазон высоты
    int count_min, count_max;                             // Диапазон количества генерируемых гауссов
    
    // ----- Ввод / вывод -----
    std::string filename;                                 // Имя файла для операций ввода/вывода
    io::BmpWriteMode bmpWriteMode;                        // Режим записи BMP (Full/Binary)
    
    // ----- binary -----
    int threshold;                                        // Порог бинаризации (0-255)    
    algorithms::components::ThresholdMode thresholdMode;  // Режим бинаризации (Peaks/Valleys/All)
    
    // ----- wave -----
    int noiseLevel;                                       // Максимальный размер шумовых компонент (в пикселях)
    
    // ----- kmeans -----
    int clusterCount;                                     // Количество кластеров для k-means
    int kernelSize;                                       // Размер ядра для кластеризации
    
    // ----- Путь -----
    double startPointX, startPointY;                      // Координаты начальной точки маршрута (Ax, Ay)
    double endPointX, endPointY;                          // Координаты конечной точки маршрута (Bx, By)
};
}
