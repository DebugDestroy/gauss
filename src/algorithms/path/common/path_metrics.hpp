#pragma once
#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/geometry/bresenham_line.hpp"

#include <vector>
#include <cmath>
#include <string>
#include <chrono>

namespace algorithms::path {

struct PathMetrics {
    double timeMs = 0.0;                           // Время работы алгоритма
    size_t pathNodes = 0;                           // Количество узлов в пути
    size_t expandedNodes = 0;                       // Количество обработанных узлов
    double euclideanLength = 0.0;                  // Евклидова длина пути
    size_t pixelLength = 0;                         // Длина пути в пикселях (Bresenham)

    std::chrono::high_resolution_clock::time_point startTime;

    // Начало измерения времени
    void startTimer() {
        startTime = std::chrono::high_resolution_clock::now();
    }

    // Вычислить метрики по пути
    void computeFromPath(const std::vector<geometry::PointD>& path) {
        pathNodes = path.size();
        euclideanLength = 0.0;
        pixelLength = 0;

        for (size_t i = 1; i < path.size(); ++i) {
            const auto& p0 = path[i - 1];
            const auto& p1 = path[i];

            euclideanLength += std::hypot(p1.x - p0.x, p1.y - p0.y);
            pixelLength += geometry::bresenhamLine(p0, p1).size();
        }
    }

    // Завершить таймер и залогировать
    void finishAndLog(core::Logger& logger, const std::string& label) {
        auto endTime = std::chrono::high_resolution_clock::now();
        timeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();

        logger.info("[PathMetrics][" + label + "] "
            "time=" + std::to_string(timeMs) + " ms, "
            "path nodes=" + std::to_string(pathNodes) + ", "
            "expanded nodes=" + std::to_string(expandedNodes) + ", "
            "euclidean length=" + std::to_string(euclideanLength) + ", "
            "pixel length=" + std::to_string(pixelLength));
            
         reset();
    }

    // Обнулить все метрики
    void reset() {
        timeMs = 0.0;
        pathNodes = 0;
        expandedNodes = 0;
        euclideanLength = 0.0;
        pixelLength = 0;
    }
};

} // namespace algorithms::path
