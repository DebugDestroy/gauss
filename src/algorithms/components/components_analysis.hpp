#pragma once

#include <vector>
#include "core/logger.hpp"
#include "core/constants.hpp"
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::components {
class Component {
private:
    core::Logger& logger;

    void logComponentStats();  // Логгирование параметров компоненты

public:
    std::vector<algorithms::geometry::Pixel> componenta;
    int min_x, min_y, max_x, max_y;
    double center_x, center_y;
    double eigenvec1_x, eigenvec1_y;
    double eigenvec2_x, eigenvec2_y;
    double eigenvalue1, eigenvalue2;
    int pixelCount;

    Component(core::Logger& lg, const std::vector<algorithms::geometry::Pixel>& inputComponenta);
    void calculate_metadata();
};
}
