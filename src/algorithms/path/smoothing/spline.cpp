#include "algorithms/path/smoothing/spline.hpp"

#include <cmath>

namespace algorithms::path::smoothing
{

algorithms::geometry::PointD Spline::catmullRom(
    const algorithms::geometry::PointD& p0,
    const algorithms::geometry::PointD& p1,
    const algorithms::geometry::PointD& p2,
    const algorithms::geometry::PointD& p3,
    double t)
{
    const double t2 = t * t;
    const double t3 = t2 * t;

    return {
        0.5 * (
            2.0 * p1.x +
            (-p0.x + p2.x) * t +
            (2.0 * p0.x - 5.0 * p1.x + 4.0 * p2.x - p3.x) * t2 +
            (-p0.x + 3.0 * p1.x - 3.0 * p2.x + p3.x) * t3
        ),

        0.5 * (
            2.0 * p1.y +
            (-p0.y + p2.y) * t +
            (2.0 * p0.y - 5.0 * p1.y + 4.0 * p2.y - p3.y) * t2 +
            (-p0.y + 3.0 * p1.y - 3.0 * p2.y + p3.y) * t3
        )
    };
}


void
Spline::spline(
    std::vector<algorithms::geometry::PointD>& path,
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
    std::size_t samplesPerSegment)
{
    if (path.size() < 3 || samplesPerSegment == 0)
        return;


std::vector<algorithms::geometry::PointD> smooth;

smooth.reserve(
    (path.size() - 1) * samplesPerSegment + 1
);

smooth.push_back(path.front());


for (std::size_t i = 0; i + 1 < path.size(); ++i)
{
    algorithms::geometry::PointD p0 = (i == 0)
        ? path[i]
        : path[i-1];

    algorithms::geometry::PointD p1 = path[i];

    algorithms::geometry::PointD p2 = path[i+1];

    algorithms::geometry::PointD p3 = (i + 2 < path.size())
        ? path[i+2]
        : path[i+1];


    std::vector<algorithms::geometry::PointD> segment;

    segment.reserve(samplesPerSegment);


    for (std::size_t j = 1; j <= samplesPerSegment; ++j)
    {
        double t =
            static_cast<double>(j) /
            static_cast<double>(samplesPerSegment);

        segment.push_back(
            catmullRom(
                p0,
                p1,
                p2,
                p3,
                t)
        );
    }


    bool valid = true;

    algorithms::geometry::PointD previous = smooth.back();


    for (const auto& point : segment)
    {
        if (!validator.isEdgeValidContinuous(
                gaussi,
                previous,
                point,
                fieldWidth,
                fieldHeight,
                heightThreshold,
                vehicleRadius,
                maxSideAngle,
                maxUpDownAngle,
                interpolationEdge,
                interpolationCollision,
                interpolationAngle))
        {
            valid = false;
            break;
        }

        previous = point;
    }


    if (valid)
    {
        smooth.insert(
            smooth.end(),
            segment.begin(),
            segment.end()
        );
    }
    else
    {
        smooth.push_back(p2);
    }
}

path = std::move(smooth);
}

} // namespace algorithms::path::smoothing
