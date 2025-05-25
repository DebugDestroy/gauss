#pragma once

#include <vector>
#include <unordered_map>
#include <memory>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/gauss/pole.hpp"
#include "algorithms/path/a_star/graph.hpp"
#include "algorithms/path/a_star/conditions.hpp"

namespace algorithms::path::a_star {

struct AStarNode {
    algorithms::geometry::PointD position;
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

    std::vector<algorithms::geometry::PointD> findPathAStar(
        const algorithms::geometry::PointD& start,
        const algorithms::geometry::PointD& goal,
        std::vector<algorithms::geometry::Edge> voronoiEdges,
        algorithms::path::a_star::Graph& graph,
        const algorithms::path::a_star::Conditions& conds,
        const std::vector<std::vector<double>>& binaryMap,
        const std::unique_ptr<algorithms::gauss::Pole>& elevationData);
};

}
