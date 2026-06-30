#include "statistics/statistics_manager.hpp"
#include "algorithms/geometry/bresenham_line.hpp" 

#include <vector>
#include <cmath>
#include <string>
#include <chrono>

namespace statistics {
    // Начало измерения времени
    void StatisticsManager::startTimer() {
        startTime = std::chrono::high_resolution_clock::now();
    }

    // Вычислить метрики по пути
    void StatisticsManager::computePathMetrics(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::Pixel>& path)
    {
        metrics.pathNodes = path.size();
        metrics.euclideanLength = 0.0;
        metrics.pixelLength = 0;

        for (size_t i = 1; i < path.size(); ++i) {
            const auto& p0 = path[i - 1];
            const auto& p1 = path[i];

            metrics.euclideanLength += std::hypot(p1.x - p0.x, p1.y - p0.y);
            metrics.pixelLength += algorithms::geometry::bresenhamLine(p0, p1).size();
        }
    }

    // Завершить таймер и посчитать время
    void  StatisticsManager::finishTimer(algorithms::path::PathMetrics& metrics) {
        auto endTime = std::chrono::high_resolution_clock::now();
        metrics.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    }
    
    // Сбросить метрики
    void StatisticsManager::reset(algorithms::path::PathMetrics& metrics) {
        metrics.algorithmName.clear();
        metrics.executionTimeMs = 0.0;
        metrics.pathNodes = 0;
        metrics.expandedNodes = 0;
        metrics.euclideanLength = 0.0;
        metrics.pixelLength = 0;
        metrics.pathFound = false;
    }
}
