#pragma once

#include <vector>
#include <memory>
#include <string>
#include "core/logger.hpp"
#include "core/constants.hpp"
#include "algorithms/gauss/pole.hpp"
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::geometry {

class VoronoiDiagram {
public:
    VoronoiDiagram(core::Logger& logger);

    void buildFromDelaunay(const std::vector<Triangle>& triangles,
                           const std::unique_ptr<algorithms::gauss::Pole>& p,
                           std::vector<Edge>& edges);

private:
    core::Logger& logger;

    void logEdge(const Edge& edge, const std::string& prefix = "") const;

    PointD findBoundaryIntersection(const PointD& start, const PointD& dir, int width, int height) const;

    std::vector<Edge> getConvexHullEdges(const Triangle& tri, const std::vector<Triangle>& allTriangles) const;
};

}
