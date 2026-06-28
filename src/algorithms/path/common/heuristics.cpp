#include "algorithms/path/common/heuristics.hpp"
#include "core/logger.hpp"

#include <cmath>
#include <string>

namespace algorithms::path::common {
   
    double euclidean(const algorithms::geometry::Pixel& a, const algorithms::geometry::Pixel& b) {
        double dist = std::hypot(a.x - b.x, a.y - b.y);
        return dist;
    }
}
