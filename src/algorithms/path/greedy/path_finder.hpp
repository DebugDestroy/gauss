#pragma once

#include <vector>
#include <memory>
#include <unordered_map>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/path/common/path_metrics.hpp"
#include "algorithms/path/common/grid.hpp"

namespace algorithms::path::greedy {

class PathFinder {
private:
    core::Logger& logger;

public:
    explicit PathFinder(core::Logger& lg);

    std::vector<algorithms::geometry::Pixel> findPathGreedyGraph(
        const algorithms::geometry::Pixel& start,
        const algorithms::geometry::Pixel& goal,
        const std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph,
        algorithms::path::PathMetrics& metrics);
        
    std::vector<algorithms::path::common::GridCell> findPathGreedyGrid(
        const algorithms::path::common::GridCell& startCell,
        const algorithms::path::common::GridCell& endCell,
        const algorithms::path::common::Grid& grid,
        algorithms::path::PathMetrics& metrics);
};

}
