#pragma once

#include <vector>
#include <unordered_map>
#include <memory>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/path/common/path_metrics.hpp"
#include "algorithms/path/common/grid.hpp"

namespace algorithms::path::a_star {

struct NodeGraph {
    algorithms::geometry::Pixel position;
    double gScore;
    double fScore;

    bool operator>(const NodeGraph& other) const {
        return fScore > other.fScore;
    }
};

struct NodeGrid {
    std::size_t idx;
    double gScore;
    double fScore;

    bool operator>(const NodeGrid& other) const {
        return fScore > other.fScore;
    }
};

class PathFinder {
private:
    core::Logger& logger;

public:
    explicit PathFinder(core::Logger& lg);

    std::vector<algorithms::geometry::Pixel> findPathAStarGraph(
        const algorithms::geometry::Pixel& start,
        const algorithms::geometry::Pixel& goal,
        const std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph,
        algorithms::path::PathMetrics& metrics);
        
    std::vector<algorithms::path::common::GridCell> findPathAStarGrid(
        const algorithms::path::common::GridCell& startCell,
        const algorithms::path::common::GridCell& endCell,
        const algorithms::path::common::Grid& grid,
        algorithms::path::PathMetrics& metrics);
};

}
