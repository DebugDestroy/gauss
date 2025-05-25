#pragma once
#include <string>
#include "io/bmp_handler.hpp"  // Для BmpWriteMode
#include "algorithms/components/binary.hpp"  // Для ThresholdMode

namespace command {
struct DispatcherParams {
    std::string command;                                  // Команда для выполнения (help, init, g, generate и т.д.)
    int fieldWidth, fieldHeight;                          // Размеры поля (ширина × высота в пикселях)
    double height;                                        // Высота гауссова распределения (в условных единицах)
    double centerX, centerY;                              // Координаты центра (x, y)
    double sigmaX, sigmaY;                                // Разброс по осям (σx, σy)
    std::string filename;                                 // Имя файла для операций ввода/вывода
    io::BmpWriteMode bmpWriteMode;                        // Режим записи BMP (Full/Binary)
    int threshold;                                        // Порог бинаризации (0-255)
    algorithms::components::ThresholdMode thresholdMode;  // Режим бинаризации (Peaks/Valleys/All)
    int noiseLevel;                                       // Максимальный размер шумовых компонент (в пикселях)
    int clusterCount;                                     // Количество кластеров для k-means
    int kernelSize;                                       // Размер ядра для кластеризации
    double startPointX, startPointY;                      // Координаты начальной точки маршрута (Ax, Ay)
    double endPointX, endPointY;                          // Координаты конечной точки маршрута (Bx, By)
};
}