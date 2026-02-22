#pragma once

#include <vector>
#include <memory>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/gauss/pole.hpp"
#include "algorithms/path/common/graph.hpp"
#include "algorithms/path/common/conditions.hpp"

namespace algorithms::path::greedy {

class PathFinder {
private:
    core::Logger& logger;

public:
    explicit PathFinder(core::Logger& lg);

    std::vector<algorithms::geometry::PointD> findPathGreedy(
        const algorithms::geometry::PointD& start,
        const algorithms::geometry::PointD& goal,
        std::vector<algorithms::geometry::Edge> voronoiEdges,
        algorithms::path::common::Graph& graph,
        const algorithms::path::common::Conditions& conds,
        const std::vector<std::vector<double>>& binaryMap,
        const std::unique_ptr<algorithms::gauss::Pole>& elevationData
    );
};

}
