#pragma once

#include <vector>

#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/gauss/gauss_builder.hpp"
#include "algorithms/path/common/path_validator.hpp"

namespace algorithms::path::smoothing
{

class Shortcut
{
public:

    static void
    shortcutDiscrete(
        std::vector<algorithms::geometry::Pixel>& path,
        const std::vector<std::vector<double>>& field,
        const std::vector<std::vector<double>>& binaryMap,
        const algorithms::path::common::PathValidator& validator,
        int vehicleRadius,
        double maxSideAngle,
        double maxUpDownAngle);

    static void
    shortcutContinuous(
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
        const algorithms::path::common::PathValidator& validator);

};

}
