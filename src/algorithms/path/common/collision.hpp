#pragma once

#include <vector>

#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/gauss/gauss_builder.hpp"

namespace algorithms::path::common {

struct ObstacleDistance
{
    double euclidean;
    int pixel;
};

// Быстрая проверка есть ли столкновение в точке
bool isVehicleRadiusValid(const algorithms::geometry::Pixel& pixel, 
                         const std::vector<std::vector<double>>& binaryMap,
                         int vehicleRadius);
                         
// Расстояние до препятсвия от центра для дискретных
ObstacleDistance minObstacleDistance(
        const algorithms::geometry::Pixel& center,
        const std::vector<std::vector<double>>& binaryMap);
        
// Расстояние до препятсвия от центра для непрерывных
ObstacleDistance minObstacleDistance(
        const algorithms::geometry::PointD& center,
        const algorithms::gauss::GaussBuilder& gaussBuilder,
        const std::vector<algorithms::gauss::Gaus>& gaussi,
        int fieldWidth,
        int fieldHeight,
        double heightThreshold,
        double interpolationCollision,
        double interpAngle);

// Проверка на столкновения в непр случае
bool checkPointContinuous(
    const algorithms::gauss::GaussBuilder& gaussBuilder,
    const std::vector<algorithms::gauss::Gaus>& gaussi,
    const algorithms::geometry::PointD& center,
    int fieldWidth,
    int fieldHeight,
    double heightThreshold,
    double vehicleRadius,
    double interpCollision,
    double interpAngle);
                    
}
