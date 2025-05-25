#pragma once

#include <vector>
#include <string>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::geometry {

class Triangulator {
private:
    core::Logger& logger;

    void logTriangle(const Triangle& tri, const std::string& prefix = "") const;
    void logEdge(const Edge& edge, const std::string& prefix = "") const;

public:
    explicit Triangulator(core::Logger& lg);

    std::vector<Triangle> bowyerWatson(
        const std::vector<PointD>& points,
        int width,
        int height
    );
};

}
