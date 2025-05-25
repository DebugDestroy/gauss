#pragma once

#include <vector>
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::geometry {
    std::vector<PointD> bresenhamLine(const PointD& start, const PointD& end);
}
