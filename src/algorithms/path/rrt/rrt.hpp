#pragma once

#include <vector>
#include <random>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/path/common/path_metrics.hpp"
#include "algorithms/path/common/path_validator.hpp"
#include "algorithms/gauss/gauss_builder.hpp"

namespace algorithms::path::rrt {

struct RRTNode
{
    algorithms::geometry::PointD point;
    int parent;          // индекс родителя
};

struct RRTResult
{
    std::vector<algorithms::geometry::PointD> path;
    std::vector<RRTNode> tree;
};

class PathFinder {
private:
    core::Logger& logger;
    std::mt19937 &gen;
    std::uniform_real_distribution<double> xDist;
    std::uniform_real_distribution<double> yDist;
    std::uniform_real_distribution<double> probDist{0.0, 1.0};
    
    algorithms::geometry::PointD steer(
    const algorithms::geometry::PointD& from,
    const algorithms::geometry::PointD& to,
    double step);
    
    algorithms::geometry::PointD randomPoint(const algorithms::geometry::PointD& goal,
                   double goalBias);
    
    int nearestNode(
        const std::vector<RRTNode>& tree,
        const algorithms::geometry::PointD& point) const;
    
    std::vector<algorithms::geometry::PointD> restorePath(
        const std::vector<RRTNode>& tree,
        int goalIndex) const;
    
public:
    explicit PathFinder(core::Logger& lg, std::mt19937& generator);

    RRTResult findPathRRT(
        const algorithms::geometry::PointD& start,
        const algorithms::geometry::PointD& goal,
        const std::vector<algorithms::gauss::Gaus>& gaussi,
        int fieldWidth,
        int fieldHeight,
        double heightThreshold,
        double vehicleRadius,
        double maxSideAngle,
        double maxUpDownAngle,
        double interpolationEdge,
        double interpolationCollision,
        double interpolationAngle,
        std::size_t maxIterations,
        double step,
        double goalRadius,
        double goalBias,
        const algorithms::path::common::PathValidator& conds,
        algorithms::path::PathMetrics& metrics);
};

}
