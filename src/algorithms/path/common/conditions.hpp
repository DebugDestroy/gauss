#pragma once

#include <vector>
#include <memory>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/path/common/grid.hpp"

namespace algorithms::path::common {

class Conditions {
private:
    core::Logger& logger;
    
public:
    Conditions(core::Logger& lg);
    
    // Быстрая проверка есть ли столкновение в точке
    bool isVehicleRadiusValid(const algorithms::geometry::Pixel& pixel, 
                         const std::vector<std::vector<double>>& binaryMap,
                         int vehicleRadius) const;
                         
    // Расстояние до препятсвия от центра Евклидово
    double minObstacleDistance(
        const algorithms::geometry::Pixel& center,
        const std::vector<std::vector<double>>& binaryMap) const;
        
    // Расстояние до препятсвия от центра пиксельное
    int minObstacleDistancePixel(
        const algorithms::geometry::Pixel& center,
        const std::vector<std::vector<double>>& binaryMap) const;
                              
    // Проверка, проходимо ли ребро (уклоны, крены, коллизии)
    bool isEdgeNavigable(const algorithms::geometry::PixelEdge& edge, 
                         const std::vector<std::vector<double>>& field,
                         const std::vector<std::vector<double>>& binaryMap,
                         int vehicleRadius,
                         double maxSideAngle,
                         double maxUpDownAngle) const;
                         
    // Проверка, проходима ли ячейка
    bool isCellFree(const algorithms::path::common::GridCell& cell, 
                    const std::vector<std::vector<double>>& binaryMap,
                    int cellSize,
                    int noisy) const;
    };
}
