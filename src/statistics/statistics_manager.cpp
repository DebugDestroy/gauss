#include "statistics/statistics_manager.hpp"
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
        
        const auto segmentPixels = algorithms::geometry::bresenhamLine(path[i], path[i + 1]);
        
        if (segmentPixels.size() < 2)
            continue;
            
        for (const auto& center : segmentPixels)
        {     
            auto angles = algorithms::kinematics::calculateVehicleAngles(
                field,
                center,
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
}

// Вычислить минимальное расстояние до препятствия в пикселях и евклидово для пиксельного пути        
void StatisticsManager::computeMinObstacleDistance(
    algorithms::path::PathMetrics& metrics,
    const std::vector<algorithms::geometry::Pixel>& path,
    const std::vector<std::vector<double>>& binaryMap)
{
    if (path.size() < 2)
        return;

    double minEuclid = std::numeric_limits<double>::infinity();
    int minPixel = std::numeric_limits<int>::max();

    for (size_t i = 0; i + 1 < path.size(); ++i)
    {
        const auto segmentPixels = algorithms::geometry::bresenhamLine(path[i], path[i + 1]);
            
        for (const auto& center : segmentPixels)
        {     
            auto result = algorithms::path::common::minObstacleDistance(center, binaryMap);

            minEuclid = std::min(minEuclid, result.euclidean);
            minPixel  = std::min(minPixel, result.pixel);
        }
    }

    metrics.minObstacleDistance = minEuclid;
    metrics.minObstacleDistancePixel = minPixel;
}

    // Завершить таймер и посчитать время
    void  StatisticsManager::finishTimer(algorithms::path::PathMetrics& metrics) {
        auto endTime = std::chrono::high_resolution_clock::now();
        metrics.executionTimeMs = std::chrono::duration<double, std::milli>(endTime - startTime).count();
    }
    
    // Вычислить максимальные углы вбок и вперед/назад для непрерывного пути
    void StatisticsManager::computeMaxTerrainAngles(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::PointD>& path,
        const algorithms::gauss::GaussBuilder& gaussBuilder,
        const std::vector<algorithms::gauss::Gaus>& gaussi,
        double vehicleRadius,
        double interpEdge)
    {
        if (path.size() < 2)
            return;

        metrics.maxSideAngle = 0.0;
        metrics.maxUpDownAngle = 0.0;
        double length;
        algorithms::geometry::PointD dir;
        algorithms::geometry::PointD direction; // Нормализованное направление движения
        algorithms::geometry::PointD center; 
        int steps;
        double t;
        
        for (size_t segment = 0; segment + 1 < path.size(); ++segment)
        {
            dir = {
                path[segment + 1].x - path[segment].x,
                path[segment + 1].y - path[segment].y
            };
            
            length = std::hypot(dir.x, dir.y);

            // Ребро нулевой длины
            if (length < core::EPSILON)
                continue;

            // Нормализованное направление движения
            direction = {
                dir.x / length,
                dir.y / length
            };
    
            steps = std::max(
                1,
                static_cast<int>(std::ceil(length / interpEdge))
            );
    
            for (int step = 0; step <= steps; ++step)
            {
                t = static_cast<double>(step) / steps;

                center = {
                    path[segment].x + t * dir.x,
                    path[segment].y + t * dir.y
                }; 

                auto angles =
                    algorithms::kinematics::calculateVehicleAnglesContinuous(
                        gaussBuilder,
                        gaussi,
                        center,
                        direction,
                        vehicleRadius);

                metrics.maxSideAngle =
                    std::max(metrics.maxSideAngle,
                        std::fabs(angles.sideAngle));

                metrics.maxUpDownAngle =
                    std::max(metrics.maxUpDownAngle,
                        std::fabs(angles.upDownAngle));
            }
        }
    }
    
    // Вычислить минимальное расстояние до препятствия в пикселях и евклидово для непрерывного пути
    void StatisticsManager::computeMinObstacleDistance(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::PointD>& path,
        const algorithms::gauss::GaussBuilder& gaussBuilder,
        const std::vector<algorithms::gauss::Gaus>& gaussi,
        int fieldWidth,
        int fieldHeight,
        double heightThreshold,
        double interpEdge,
        double interpolationCollision,
        double interpAngle)
    {
        if (path.size() < 2)
            return;

        double minEuclid = std::numeric_limits<double>::infinity();
        int minPixel = std::numeric_limits<int>::max();
        
        double length;
        algorithms::geometry::PointD dir;
        algorithms::geometry::PointD center; 
        int steps;
        double t;
        
        for (size_t segment = 0; segment + 1 < path.size(); ++segment)
        {
            length = algorithms::geometry::distance(path[segment], path[segment + 1]);
            
            dir = {
                path[segment + 1].x - path[segment].x,
                path[segment + 1].y - path[segment].y
            };
            
            // Ребро нулевой длины
            if (length < core::EPSILON)
            {
                auto result = algorithms::path::common::minObstacleDistance(path[segment], gaussBuilder, gaussi, fieldWidth, fieldHeight, heightThreshold, interpolationCollision, interpAngle);
                minEuclid = std::min(minEuclid, result.euclidean);
                minPixel  = std::min(minPixel, result.pixel);
                continue;
            }
            
            steps = std::max(
                1,
                static_cast<int>(std::ceil(length / interpEdge))
            );
    
            for (int step = 0; step <= steps; ++step)
            {
                t = static_cast<double>(step) / steps;

                center = {
                    path[segment].x + t * dir.x,  
                    path[segment].y + t * dir.y
                };
                
                auto result = algorithms::path::common::minObstacleDistance(center, gaussBuilder, gaussi, fieldWidth, fieldHeight, heightThreshold, interpolationCollision, interpAngle);

                minEuclid = std::min(minEuclid, result.euclidean);
                minPixel  = std::min(minPixel, result.pixel);
            }
        }
        metrics.minObstacleDistance = minEuclid;
        metrics.minObstacleDistancePixel = minPixel;
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
