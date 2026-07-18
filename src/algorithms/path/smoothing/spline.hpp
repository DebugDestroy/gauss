#pragma once

#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/gauss/gauss_builder.hpp"
#include "algorithms/path/common/path_validator.hpp"

#include <vector>

namespace algorithms::path::smoothing
{

class Spline
{
public:

    static void
    spline(
        std::vector<algorithms::geometry::PointD>& path,
        const algorithms::gauss::GaussBuilder& gaussBuilder,
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
        const algorithms::path::common::PathValidator& validator,
        std::size_t samplesPerSegment);


private:

    static algorithms::geometry::PointD catmullRom(
        const algorithms::geometry::PointD& p0,
        const algorithms::geometry::PointD& p1,
        const algorithms::geometry::PointD& p2,
        const algorithms::geometry::PointD& p3,
        double t);

};

}
