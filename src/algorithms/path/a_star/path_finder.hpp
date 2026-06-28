#pragma once

#include <vector>
#include <unordered_map>
#include <memory>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::path::a_star {

struct AStarNode {
    algorithms::geometry::Pixel position;
    double gScore;
    double fScore;

    bool operator>(const AStarNode& other) const {
        return fScore > other.fScore;
    }
};

class PathFinder {
private:
    core::Logger& logger;

public:
    explicit PathFinder(core::Logger& lg);

    std::vector<algorithms::geometry::Pixel> findPathAStar(
        const algorithms::geometry::Pixel& start,
        const algorithms::geometry::Pixel& goal,
        const std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph);
};

}
