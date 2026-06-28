#pragma once

#include <vector>
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::geometry {
    std::vector<Pixel> bresenhamLine(const Pixel& start, const Pixel& end);
}
