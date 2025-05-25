#pragma once

#include <vector>
#include <memory>

#include "core/config.hpp"
#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/gauss/pole.hpp"

namespace algorithms::path::a_star {

class Conditions {
private:
    const core::Config& config;
    core::Logger& logger;

public:
    Conditions(const core::Config& cfg, core::Logger& lg);

    // Проверка на коллизию по радиусу тележки
    bool isVehicleRadiusValid(const algorithms::geometry::PointD& pixel, 
                              const std::vector<std::vector<double>>& binaryMap) const;

    // Проверка, проходимо ли ребро (уклоны, крены, коллизии)
    bool isEdgeNavigable(const algorithms::geometry::Edge& edge, 
                         const std::unique_ptr<algorithms::gauss::Pole>& p,
                         const std::vector<std::vector<double>>& binaryMap) const;
};

}
