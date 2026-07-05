#pragma once
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/path/common/path_metrics.hpp"
#include "algorithms/path/common/conditions.hpp"

#include <vector>
#include <chrono>

namespace statistics {
class StatisticsManager {
public:
    // Начать отсчет
    void startTimer();
    
    // Завершить таймер и посчитать время
    void finishTimer(algorithms::path::PathMetrics& metrics);

    // Вычислить длину пути
    void computePathLength(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::Pixel>& path);
    
    // Вычислить максимальные углы вбок и вперед/назад
    void computeMaxTerrainAngles(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::Pixel>& path,
        const std::vector<std::vector<double>>& field,
        int vehicleRadius);
        
    // Вычислить минимальное расстояние до препятствия в пикселях и евклидово
    void computeMinObstacleDistance(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::Pixel>& path,
        const std::vector<std::vector<double>>& binaryMap,
        const algorithms::path::common::Conditions& conds);
          
    // Сбросить метрики
    void reset(algorithms::path::PathMetrics& metrics);
            
private:
    std::chrono::high_resolution_clock::time_point startTime;
};

} // namespace statistics
