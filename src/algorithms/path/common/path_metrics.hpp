#pragma once
#include <cstddef>
#include <string>

namespace algorithms::path {

struct PathMetrics {
    std::string environment;                        // Рабочая среда
    std::string algorithmName;                      // Имя алгортма
    double executionTimeMs = 0.0;                   // Время работы алгоритма
    size_t pathNodes = 0;                           // Количество узлов в пути
    size_t expandedNodes = 0;                       // Количество обработанных узлов
    double euclideanLength = 0.0;                   // Евклидова длина пути
    size_t pixelLength = 0;                         // Длина пути в пикселях (Bresenham)
    bool pathFound = false;                         // Найден ли путь
    double minObstacleDistance = 0;                 // Минимальное Евклидово расстояние до препятсвия
    int minObstacleDistancePixel = 0;               // Минимальное пиксельное расстояние до препятсвия
    double maxSideAngle = 0;                        // Максимальный наклон вбок
    double maxUpDownAngle = 0;                      // Максимальный наклон вперед/назад
};

} // namespace algorithms::path
