#pragma once

#include <vector>
#include <unordered_map>
#include <memory>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::path::dekstra {

struct DijkstraNode {
    algorithms::geometry::Pixel position;
    double gScore;

    bool operator>(const DijkstraNode& other) const {
        return gScore > other.gScore;
    }
};

class PathFinder {
private:
    core::Logger& logger;

public:
    explicit PathFinder(core::Logger& lg);

    std::vector<algorithms::geometry::Pixel> findPathDijkstra(
        const algorithms::geometry::Pixel& start,
        const algorithms::geometry::Pixel& goal,
        const std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph);
};

}
