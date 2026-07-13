#include "algorithms/path/smoothing/shortcut.hpp"

namespace algorithms::path::smoothing
{

void
Shortcut::shortcutDiscrete(
    std::vector<algorithms::geometry::Pixel>& path,
    const std::vector<std::vector<double>>& field,
    const std::vector<std::vector<double>>& binaryMap,
    const algorithms::path::common::PathValidator& validator,
    int vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle)
{
    if (path.size() <= 2)
        return;


    std::vector<algorithms::geometry::Pixel> result;
    result.reserve(path.size());


    std::size_t i = 0;

    result.push_back(path.front());


    while (i < path.size() - 1)
    {
        bool connected = false;


        for (std::size_t j = path.size() - 1; j > i; --j)
        {
            algorithms::geometry::PixelEdge edge(path[i], path[j]);


            if (!validator.isEdgeNavigable(
                    edge,
                    field,
                    binaryMap,
                    vehicleRadius,
                    maxSideAngle,
                    maxUpDownAngle))
            {
                continue;
            }


            result.push_back(path[j]);

            i = j;

            connected = true;

            break;
        }


        if (!connected)
        {
            // Оставляем следующую вершину
            result.push_back(path[i + 1]);

            ++i;
        }
    }


    path = std::move(result);
}

void
Shortcut::shortcutContinuous(
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
    const algorithms::path::common::PathValidator& validator)
{
    if (path.size() <= 2)
        return;

    std::vector<algorithms::geometry::PointD> result;
    result.reserve(path.size());
    
    std::size_t i = 0;

    while (i < path.size() - 1)
    {
        result.push_back(path[i]);

        for (std::size_t j = path.size() - 1; j > i; --j)
        {
            if (!validator.isEdgeValidContinuous(
                    gaussi,
                    path[i],
                    path[j],
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
                continue;
            }

            i = j;
            break;
        }
    }

    result.push_back(path.back());
    path = std::move(result);
}

} // namespace algorithms::path::smoothing
