#pragma once

#include <vector>
#include <memory>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::path::common {

class Conditions {
private:
    core::Logger& logger;

public:
    Conditions(core::Logger& lg);

    // Проверка на коллизию по радиусу тележки
    bool isVehicleRadiusValid(const algorithms::geometry::Pixel& pixel, 
                              const std::vector<std::vector<double>>& binaryMap,
                              int vehicleRadius) const;

    // Проверка, проходимо ли ребро (уклоны, крены, коллизии)
    bool isEdgeNavigable(const algorithms::geometry::PixelEdge& edge, 
                         const std::vector<std::vector<double>>& field,
                         const std::vector<std::vector<double>>& binaryMap,
                         int vehicleRadius,
                         double maxSideAngle,
                         double maxUpDownAngle) const;
    };
}
