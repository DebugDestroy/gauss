#include "algorithms/path/common/collision.hpp"
#include "algorithms/geometry/geometry_structures.hpp"
#include "core/constants.hpp"

#include <vector>
#include <cmath>
#include <limits>
#include <numbers>

namespace algorithms::path::common {

 // Быстрая проверка есть ли столкновение в точке
bool isVehicleRadiusValid(const algorithms::geometry::Pixel& pixel, 
                         const std::vector<std::vector<double>>& binaryMap,
                         int vehicleRadius)
{
    const int x = pixel.x;
    const int y = pixel.y;

    const int radiusPixels = vehicleRadius;
    const int radiusSquared = vehicleRadius * vehicleRadius;

    for (int dy = -radiusPixels; dy <= radiusPixels; ++dy) {
        for (int dx = -radiusPixels; dx <= radiusPixels; ++dx) {
            if (dx * dx + dy * dy > radiusSquared) continue;

            const int nx = x + dx;
            const int ny = y + dy;
            
            if (nx < 0 || ny < 0 ||
                ny >= static_cast<int>(binaryMap.size()) ||
                nx >= static_cast<int>(binaryMap[0].size()))
            {
                return false;
            }
            
            if (std::fabs(binaryMap[static_cast<std::size_t>(ny)][static_cast<std::size_t>(nx)] - core::WHITE) < core::EPSILON) {
                return false;     
            }
        }
    }

    return true;
}
                         
// Расстояние до препятсвия от центра для дискретных
ObstacleDistance minObstacleDistance(
    const algorithms::geometry::Pixel& center,
    const std::vector<std::vector<double>>& binaryMap)
{
    const int height = static_cast<int>(binaryMap.size());
    const int width  = static_cast<int>(binaryMap[0].size());

    const int maxRadius = std::max(width, height);

    for (int r = 0; r <= maxRadius; ++r)
    {
        double minEuclid = std::numeric_limits<double>::infinity();
        bool found = false;

        for (int dy = -r; dy <= r; ++dy)
        {
            for (int dx = -r; dx <= r; ++dx)
            {
                if (r > 0 && dx * dx + dy * dy <= (r - 1) * (r - 1))
                    continue;

                if (dx * dx + dy * dy > r * r)
                    continue;

                const int x = center.x + dx;
                const int y = center.y + dy;

                bool obstacle =
                    (x < 0 || x >= width || y < 0 || y >= height);

                if (!obstacle)
                {
                    obstacle =
                        std::fabs(binaryMap[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)] - core::WHITE) < core::EPSILON;
                }

                if (obstacle)
                {
                    found = true;
                    minEuclid = std::min(
                        minEuclid,
                        std::hypot(static_cast<double>(dx),
                                   static_cast<double>(dy)));
                }
            }
        }

        if (found)
        {
            return {minEuclid, r};
        }
    }

    return {
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<int>::max()
    };
}

// Расстояние до препятсвия от центра для непрерывных
ObstacleDistance minObstacleDistance(
    const algorithms::geometry::PointD& center,
    const std::vector<algorithms::gauss::Gaus>& gaussi,
    int fieldWidth,
    int fieldHeight,
    double heightThreshold,
    double interpolationCollision,
    double interpAngle)
{
    const double maxRadius =
        static_cast<double>(std::max(fieldWidth, fieldHeight));

    auto isPointValid = [&](const algorithms::geometry::PointD& p)
    {
        if (p.x < 0.0 || p.y < 0.0 ||
            p.x >= fieldWidth ||
            p.y >= fieldHeight)
        {
            return false;
        }

        const double h =
            algorithms::gauss::GaussBuilder::heightAt(
                p,
                gaussi);

        return std::abs(h - core::MID_GRAY) < heightThreshold;
    };


    if (!isPointValid(center))
        return {0.0, 0};


    for (double radius = interpolationCollision;
         radius <= maxRadius;
         radius += interpolationCollision)
    {
        const double circumference =
            2.0 * std::numbers::pi * radius;

        const int samples = std::max(
            1,
            static_cast<int>(
                std::ceil(circumference / interpAngle))
        );


        for (int i = 0; i < samples; ++i)
        {
            const double angle =
                2.0 * std::numbers::pi * i / samples;

            algorithms::geometry::PointD p{
                center.x + radius * std::cos(angle),
                center.y + radius * std::sin(angle)
            };


            if (!isPointValid(p))
            {
                return {
                    radius,
                    static_cast<int>(std::floor(radius))
                };
            }
        }
    }


    return {
        std::numeric_limits<double>::infinity(),
        std::numeric_limits<int>::max()
    };
}

// Проверка на столкновения в непрерывном случае
bool checkPointContinuous(
    const std::vector<algorithms::gauss::Gaus>& gaussi,
    const algorithms::geometry::PointD& center,
    int fieldWidth,
    int fieldHeight,
    double heightThreshold,
    double vehicleRadius,
    double interpCollision,
    double interpAngle)
{
    auto isPointValid = [&](const algorithms::geometry::PointD& p)
    {
        if (p.x < 0.0 || p.y < 0.0 ||
            p.x >= fieldWidth ||
            p.y >= fieldHeight)
        {
            return false;
        }

        const double h =
            algorithms::gauss::GaussBuilder::heightAt(
                p,
                gaussi);

        return std::abs(h - core::MID_GRAY) < heightThreshold;
    };


    // Проверяем сам центр
    if (!isPointValid(center))
        return false;


    // Проверяем окружности радиуса interpCollision, 2*interpCollision...
    for (double radius = interpCollision;
         radius <= vehicleRadius;
         radius += interpCollision)
    {
        const double circumference =
            2.0 * std::numbers::pi * radius;

        const int samples = std::max(
            1,
            static_cast<int>(std::ceil(circumference / interpAngle))
        );

        for (int i = 0; i < samples; ++i)
        {
            const double angle =
                2.0 * std::numbers::pi * i / samples;

            algorithms::geometry::PointD p{
                center.x + radius * std::cos(angle),
                center.y + radius * std::sin(angle)
            };

            if (!isPointValid(p))
                return false;
        }
    }

    return true;
}
                    
}
