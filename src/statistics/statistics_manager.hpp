#pragma once
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/path/common/path_metrics.hpp"
#include "algorithms/path/common/collision.hpp"
#include "algorithms/geometry/bresenham_line.hpp" 
#include "algorithms/geometry/math.hpp"

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
    template<class Point>
    void computePathLength(algorithms::path::PathMetrics& metrics, const std::vector<Point>& path)
    {
        metrics.pathNodes = path.size();
        metrics.euclideanLength = 0.0;
        metrics.pixelLength = 0;

        for (size_t i = 1; i < path.size(); ++i) {
            const auto& p0 = path[i - 1];
            const auto& p1 = path[i];

            metrics.euclideanLength += std::hypot(p1.x - p0.x, p1.y - p0.y);
            metrics.pixelLength +=
                algorithms::geometry::bresenhamLine(
                    algorithms::geometry::toPixel(p0), 
                    algorithms::geometry::toPixel(p1)
                ).size() - 1;
        }
        if (!path.empty())
            metrics.pixelLength += 1;
    }
    
    // Вычислить максимальные углы вбок и вперед/назад для пиксельного пути
    void computeMaxTerrainAngles(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::Pixel>& path,
        const std::vector<std::vector<double>>& field,
        int vehicleRadius);
    
    // Вычислить минимальное расстояние до препятствия в пикселях и евклидово для пиксельного пути
    void computeMinObstacleDistance(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::Pixel>& path,
        const std::vector<std::vector<double>>& binaryMap);
        
    // Вычислить максимальные углы вбок и вперед/назад для непрерывного пути
    void computeMaxTerrainAngles(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::PointD>& path,
        const std::vector<algorithms::gauss::Gaus>& gaussi,
        double vehicleRadius,
        double interpEdge);    
        
    // Вычислить минимальное расстояние до препятствия в пикселях и евклидово для непрерывного пути
    void computeMinObstacleDistance(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::PointD>& path,
        const std::vector<algorithms::gauss::Gaus>& gaussi,
        int fieldWidth,
        int fieldHeight,
        double heightThreshold,
        double interpEdge,
        double interpolationCollision,
        double interpAngle);
          
    // Сбросить метрики
    void reset(algorithms::path::PathMetrics& metrics);
            
private:
    std::chrono::high_resolution_clock::time_point startTime;
};

} // namespace statistics
