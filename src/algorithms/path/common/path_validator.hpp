#pragma once

#include <vector>
#include <memory>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/path/common/grid.hpp"
#include "algorithms/gauss/gauss_builder.hpp"

namespace algorithms::path::common {

class PathValidator {
private:
    core::Logger& logger;
    
public:
    PathValidator(core::Logger& lg);
                              
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
                    std::size_t noisy) const;
    
    // Проверка, проходимо ли ребро (уклоны, крены, коллизии) в непр случае
    bool isEdgeValidContinuous(
    const std::vector<algorithms::gauss::Gaus>& gaussi,
    const algorithms::geometry::PointD& A,
    const algorithms::geometry::PointD& B,
    int fieldWidth,
    int fieldHeight,
    double heightThreshold,
    double vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle,
    double interpolationEdge,
    double interpolationCollision,
    double interpolationAngle) const;
};

}
