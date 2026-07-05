#pragma once

#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::path::common {

// Эвристическая евклидова функция расстояния
    inline double euclidean(const algorithms::geometry::Pixel& a, const algorithms::geometry::Pixel& b) {
return std::hypot(a.x - b.x, a.y - b.y);
    }
    
    inline double euclidean(const algorithms::path::common::GridCell& a, const algorithms::path::common::GridCell& b) {
return std::hypot(a.row - b.row, a.col - b.col);
    }

}
