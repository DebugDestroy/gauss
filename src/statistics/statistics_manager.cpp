#include "statistics/statistics_manager.hpp"
#include "algorithms/geometry/bresenham_line.hpp" 
#include "algorithms/kinematics/incline_angle.hpp"

#include <vector>
#include <cmath>
#include <string>
#include <chrono>
#include <limits>

namespace statistics {
    // Начало измерения времени
    void StatisticsManager::startTimer() {
        startTime = std::chrono::high_resolution_clock::now();
    }

    // Вычислить длину пути
    void StatisticsManager::computePathLength(
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
    
// Вычислить максимальные углы вбок и вперед/назад
void StatisticsManager::computeMaxTerrainAngles(
    algorithms::path::PathMetrics& metrics,
    const std::vector<algorithms::geometry::Pixel>& path,
    const std::vector<std::vector<double>>& field,
    int vehicleRadius)
{
    if (path.size() < 2)
        return;

    metrics.maxSideAngle = 0.0;
    metrics.maxUpDownAngle = 0.0;

    for (size_t i = 0; i + 1 < path.size(); ++i)
    {
        algorithms::geometry::PointD dir{
            static_cast<double>(path[i + 1].x - path[i].x),
            static_cast<double>(path[i + 1].y - path[i].y)
        };

        auto angles = algorithms::kinematics::calculateVehicleAngles(
            field,
            path[i],
            dir,
            vehicleRadius);

        metrics.maxSideAngle =
            std::max(metrics.maxSideAngle,
                     std::fabs(angles.sideAngle));

        metrics.maxUpDownAngle =
            std::max(metrics.maxUpDownAngle,
                     std::fabs(angles.upDownAngle));
    }
}
        
void StatisticsManager::computeMinObstacleDistance(
    algorithms::path::PathMetrics& metrics,
    const std::vector<algorithms::geometry::Pixel>& path,
    const std::vector<std::vector<double>>& binaryMap,
    const algorithms::path::common::Conditions& conds)
{
    if (path.empty())
        return;

    double minEuclid = std::numeric_limits<double>::infinity();
    int minPixel = std::numeric_limits<int>::max();

    for (const auto& p : path)
    {
        double d = conds.minObstacleDistance(p, binaryMap);
        int dp   = conds.minObstacleDistancePixel(p, binaryMap);

        minEuclid = std::min(minEuclid, d);
        minPixel  = std::min(minPixel, dp);
    }

    metrics.minObstacleDistance = minEuclid;
    metrics.minObstacleDistancePixel = minPixel;
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
        metrics.minObstacleDistance = 0;
        metrics.minObstacleDistancePixel = 0; 
        metrics.maxSideAngle = 0;
        metrics.maxUpDownAngle = 0;
    }
}
