#include "algorithms/path/a_star/heuristics.hpp"
#include "core/logger.hpp"

#include <cmath>
#include <string>

namespace algorithms::path::a_star {
   
    double euclidean(const algorithms::geometry::PointD& a, const algorithms::geometry::PointD& b) {
        double dist = std::hypot(a.x - b.x, a.y - b.y);
        return dist;
    }
}
