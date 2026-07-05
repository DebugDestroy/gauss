#include "algorithms/path/common/conditions.hpp"
#include "algorithms/geometry/bresenham_line.hpp"
#include "algorithms/kinematics/incline_angle.hpp"

#include <cmath>
#include <sstream>
#include <limits>

namespace algorithms::path::common {

    Conditions::Conditions(core::Logger& lg) : logger(lg) {
        logger.trace("[Conditions] Инициализация проверок проходимости пути");
    }
// -----------------------------------------------------Граф-----------------------------------------------------    
bool Conditions::isVehicleRadiusValid(const algorithms::geometry::Pixel& pixel, 
                         const std::vector<std::vector<double>>& binaryMap,
                         int vehicleRadius) const
{
    logger.trace("[PathFinder::isVehicleRadiusValid] Проверка радиуса в точке (" +
                 std::to_string(pixel.x) + "," + std::to_string(pixel.y) + ")");

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
            
            if (std::fabs(binaryMap[ny][nx] - core::WHITE) < core::EPSILON) {
                return false;     
            }
        }
    }

    return true;
}

double Conditions::minObstacleDistance(
    const algorithms::geometry::Pixel& center,
    const std::vector<std::vector<double>>& binaryMap) const
{
    const int height = binaryMap.size();
    const int width  = binaryMap[0].size();

    const int maxRadius = std::max(width, height);

    double obstacleDistance = std::numeric_limits<double>::infinity();

    logger.trace("[minObstacleDistance] start center=(" +
                std::to_string(center.x) + "," +
                std::to_string(center.y) + ")");

    for (int r = 0; r <= maxRadius; ++r)
    {
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

                if (x < 0 || x >= width || y < 0 || y >= height)
                {
                    double d = std::hypot(dx, dy);

                    logger.trace("[minObstacleDistance] boundary hit at r=" +
                                std::to_string(r) +
                                " dist=" + std::to_string(d));

                    obstacleDistance = std::min(d, obstacleDistance);
                    continue;
                }

                if (std::fabs(binaryMap[y][x] - core::WHITE) < core::EPSILON)
                {
                    double d = std::hypot(dx, dy);

                    logger.trace("[minObstacleDistance] obstacle hit at (" +
                                std::to_string(x) + "," +
                                std::to_string(y) +
                                ") r=" + std::to_string(r) +
                                " dist=" + std::to_string(d));

                    obstacleDistance = std::min(d, obstacleDistance);
                }
            }
        }

        // если уже нашли ближайшее возможное расстояние — можно выйти
        if (std::isfinite(obstacleDistance))
        {
            logger.trace("[minObstacleDistance] return early r=" +
                        std::to_string(r) +
                        " result=" + std::to_string(obstacleDistance));

            return obstacleDistance;
        }
    }

    logger.error("[minObstacleDistance] no obstacle found -> inf");

    return std::numeric_limits<double>::infinity();
}

int Conditions::minObstacleDistancePixel(
    const algorithms::geometry::Pixel& center,
    const std::vector<std::vector<double>>& binaryMap) const
{
    const int height = binaryMap.size();
    const int width  = binaryMap[0].size();

    const int maxRadius = std::max(width, height);

    logger.trace("[minObstacleDistancePixel] start center=(" +
                std::to_string(center.x) + "," +
                std::to_string(center.y) + ")");

    for (int r = 0; r <= maxRadius; ++r)
    {
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

                if (x < 0 || x >= width || y < 0 || y >= height)
                {
                    logger.trace("[minObstacleDistancePixel] boundary hit -> return r=" +
                                std::to_string(r));

                    return r;
                }

                if (std::fabs(binaryMap[y][x] - core::WHITE) < core::EPSILON)
                {
                    logger.trace("[minObstacleDistancePixel] obstacle hit at (" +
                                std::to_string(x) + "," +
                                std::to_string(y) +
                                ") -> return r=" +
                                std::to_string(r));

                    return r;
                }
            }
        }
    }

    logger.error("[minObstacleDistancePixel] no obstacle found -> INT_MAX");

    return std::numeric_limits<int>::max();
}

bool Conditions::isEdgeNavigable(
    const algorithms::geometry::PixelEdge& edge,
    const std::vector<std::vector<double>>& field,
    const std::vector<std::vector<double>>& binaryMap,
    int vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle) const
{
    const auto line = algorithms::geometry::bresenhamLine(edge.a, edge.b);

    if (line.size() < 2) {
        logger.debug("Слишком короткое ребро для проверки");
        return false;
    }

    logger.trace("Длина ребра: " + std::to_string(line.size()) + " точек");

    algorithms::geometry::PointD dir{
        static_cast<double>(line.back().x - line.front().x),
        static_cast<double>(line.back().y - line.front().y)
    };

    logger.debug("--- Начало проверки ребра ---");
    logger.debug("Параметры тележки: радиус=" +
                 std::to_string(vehicleRadius) +
                 ", maxUpDownAngle=" +
                 std::to_string(maxUpDownAngle) +
                 "°, maxSideAngle=" +
                 std::to_string(maxSideAngle) + "°");

    for (const auto& center : line)
    {
        // Сначала проверяем столкновения
        if (!isVehicleRadiusValid(center, binaryMap, vehicleRadius))
        {
            logger.debug("Столкновение в (" +
                         std::to_string(center.x) + ", " +
                         std::to_string(center.y) + ")");
            return false;
        }

        // Затем вычисляем углы
        const auto angles =
            algorithms::kinematics::calculateVehicleAngles(
                field,
                center,
                dir,
                vehicleRadius);

        logger.trace("Углы в (" +
                     std::to_string(center.x) + ", " +
                     std::to_string(center.y) +
                     "): side=" +
                     std::to_string(angles.sideAngle) +
                     "°, upDown=" +
                     std::to_string(angles.upDownAngle) + "°");

        if (std::fabs(angles.upDownAngle) > maxUpDownAngle)
        {
            logger.debug("Превышен продольный наклон: " +
                         std::to_string(angles.upDownAngle) + "°");
            return false;
        }

        if (std::fabs(angles.sideAngle) > maxSideAngle)
        {
            logger.debug("Превышен боковой наклон: " +
                         std::to_string(angles.sideAngle) + "°");
            return false;
        }
    }

    logger.debug("--- Ребро проходимо ---");
    return true;
}

// -----------------------------------------------------Сетка-----------------------------------------------------
    bool Conditions::isCellFree(const algorithms::path::common::GridCell& cell, 
                    const std::vector<std::vector<double>>& binaryMap,
                    int cellSize,
                    int noisy) const
{
    int dangerousPixels = 0;

    const int x0 = cell.col * cellSize;
    const int y0 = cell.row * cellSize;

    for (int y = y0; y < y0 + cellSize; ++y) {
        for (int x = x0; x < x0 + cellSize; ++x) {

            if (std::fabs(binaryMap[y][x] - core::WHITE) < core::EPSILON) {
                ++dangerousPixels;

                if (dangerousPixels >= noisy)
                    return false;
            }
        }
    }

    return true;
}

}
