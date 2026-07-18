#include "algorithms/path/common/path_validator.hpp"
#include "algorithms/geometry/bresenham_line.hpp"
#include "algorithms/kinematics/incline_angle.hpp"
#include "algorithms/path/common/collision.hpp"
#include "core/constants.hpp"

#include <cmath>
#include <sstream>
#include <limits>

namespace algorithms::path::common {

    PathValidator::PathValidator(core::Logger& lg) : logger(lg) {
        logger.trace("[PathValidator] Инициализация проверок проходимости пути");
    }
// -----------------------------------------------------Граф-----------------------------------------------------    
bool PathValidator::isEdgeNavigable(
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
        if (!algorithms::path::common::isVehicleRadiusValid(center, binaryMap, vehicleRadius))
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
    bool PathValidator::isCellFree(const algorithms::path::common::GridCell& cell, 
                    const std::vector<std::vector<double>>& binaryMap,
                    int cellSize,
                    std::size_t noisy) const
{
    std::size_t dangerousPixels = 0;

    const int x0 = cell.col * cellSize;
    const int y0 = cell.row * cellSize;

    for (int y = y0; y < y0 + cellSize; ++y) {
        for (int x = x0; x < x0 + cellSize; ++x) {

            if (std::fabs(binaryMap[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)] - core::WHITE) < core::EPSILON) {
                ++dangerousPixels;

                if (dangerousPixels >= noisy)
                    return false;
            }
        }
    }

    return true;
}

// -----------------------------------------------------Непрерывно-----------------------------------------------------
bool PathValidator::isEdgeValidContinuous(
    const algorithms::gauss::GaussBuilder& gaussBuilder,
    const std::vector<algorithms::gauss::Gaus>& gaussi,
    const algorithms::geometry::PointD& A,
    const algorithms::geometry::PointD& B,
    int fieldWidth,
    int fieldHeight,
    double heightThreshold,
    double vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle,
    double interpolationEdge,
    double interpolationCollision,
    double interpolationAngle) const
{
    using algorithms::geometry::PointD;

    PointD dir{
        B.x - A.x,
        B.y - A.y
    };

    double length = std::hypot(dir.x, dir.y);

    // Ребро нулевой длины
    if (length < core::EPSILON) {
        return algorithms::path::common::checkPointContinuous(
            gaussBuilder,
            gaussi,
            A,
            fieldWidth,
            fieldHeight,
            heightThreshold,
            vehicleRadius,
            interpolationCollision,
            interpolationAngle);
    }

    // Нормализованное направление движения
    PointD direction{
        dir.x / length,
        dir.y / length
    };

    int steps = std::max(
        1,
        static_cast<int>(std::ceil(length / interpolationEdge))
    );

    for (int i = 0; i <= steps; ++i)
    {

        double t = static_cast<double>(i) / steps;

        PointD center{
            A.x + t * dir.x,
            A.y + t * dir.y
        };

        // Проверяем столкновение
        if (!algorithms::path::common::checkPointContinuous(
                gaussBuilder,
                gaussi,
                center,
                fieldWidth,
                fieldHeight,
                heightThreshold,
                vehicleRadius,
                interpolationCollision,
                interpolationAngle))
        {
            return false;
        }

        // Проверяем наклон тележки
        auto angles =
            algorithms::kinematics::calculateVehicleAnglesContinuous(
                gaussBuilder,
                gaussi,
                center,
                direction,
                vehicleRadius);

        if (std::abs(angles.sideAngle) > maxSideAngle)
            return false;

        if (std::abs(angles.upDownAngle) > maxUpDownAngle)
            return false;
    }

    return true;
}

}
