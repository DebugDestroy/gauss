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

    size_t checkedPixels = 0;
    size_t collisionPixels = 0;

    for (int dy = -radiusPixels; dy <= radiusPixels; ++dy) {
        for (int dx = -radiusPixels; dx <= radiusPixels; ++dx) {
            if (dx * dx + dy * dy > radiusSquared) continue;

            const int nx = x + dx;
            const int ny = y + dy;
            ++checkedPixels;

            if (nx >= 0 && ny >= 0 &&
                ny < static_cast<int>(binaryMap.size()) &&
                nx < static_cast<int>(binaryMap[0].size()) &&
                std::fabs(binaryMap[ny][nx] - core::WHITE) < core::EPSILON) {

                ++collisionPixels;
                logger.trace("[PathFinder::isVehicleRadiusValid] Столкновение в (" +
                             std::to_string(nx) + "," + std::to_string(ny) + ")");
            }
        }
    }

    bool result = (collisionPixels == 0);
    logger.debug("[PathFinder::isVehicleRadiusValid] Результат: " + 
                 std::string(result ? "проходимо" : "столкновение") +
                 ", проверено пикселей: " + std::to_string(checkedPixels) + 
                 ", столкновений: " + std::to_string(collisionPixels));

    return result;
}

    bool Conditions::isEdgeNavigable(const algorithms::geometry::PixelEdge& edge, 
                    const std::vector<std::vector<double>>& field,
                    const std::vector<std::vector<double>>& binaryMap,
                    int vehicleRadius,
                    double maxSideAngle,
                    double maxUpDownAngle) const 
{
                    
    const auto line = algorithms::geometry::bresenhamLine(edge.a, edge.b);
    
     // Проверка на минимальную длину ребра
        if (line.size() < 2) {
            logger.debug("Слишком короткое ребро для проверки");
            return false;
        }
        logger.trace("Длина ребра: " + std::to_string(line.size()) + " точек");
        
    // Вектор пути
    algorithms::geometry::PointD dir = 
        {static_cast<double>(line.back().x - line.front().x), 
            static_cast<double>(line.back().y - line.front().y)};
    double dirLen = std::hypot(dir.x, dir.y);

    // Защита от деления на ноль
    if (dirLen < core::EPSILON) {
        logger.debug("Нулевая длина направления");
        return false;
    }

    dir.x /= dirLen;
    dir.y /= dirLen;
    algorithms::geometry::PointD perp = {-dir.y, dir.x}; // Единичный перпендикуляр
    const double L = vehicleRadius / std::sqrt(2.0);
    const double W = vehicleRadius / std::sqrt(2.0);

auto isInside = [&](const algorithms::geometry::Pixel& p) {
    return p.x >= 0 &&
           p.y >= 0 &&
           p.x < field[0].size() &&
           p.y < field.size();
};

    logger.debug("--- Начало проверки ребра ---");
    logger.debug("Параметры тележки: радиус=" + std::to_string(vehicleRadius) + 
                ", max_уклон=" + std::to_string(maxUpDownAngle) + 
                "°, max_крен=" + std::to_string(maxSideAngle) + "°");
    
    for (size_t i = 0; i < line.size(); ++i) {
        const algorithms::geometry::Pixel& center = line[i];
            algorithms::geometry::Pixel frontLeft = {
    static_cast<int>(std::round(center.x + L * dir.x + W * perp.x)),
    static_cast<int>(std::round(center.y + L * dir.y + W * perp.y))
};

algorithms::geometry::Pixel frontRight = {
    static_cast<int>(std::round(center.x + L * dir.x - W * perp.x)),
    static_cast<int>(std::round(center.y + L * dir.y - W * perp.y))
};

algorithms::geometry::Pixel rearLeft = {
    static_cast<int>(std::round(center.x - L * dir.x + W * perp.x)),
    static_cast<int>(std::round(center.y - L * dir.y + W * perp.y))
};

algorithms::geometry::Pixel rearRight = {
    static_cast<int>(std::round(center.x - L * dir.x - W * perp.x)),
    static_cast<int>(std::round(center.y - L * dir.y - W * perp.y))
};
        // Проверка границ массива
            if (!isInside(frontLeft) ||
    !isInside(frontRight) ||
    !isInside(rearLeft) ||
    !isInside(rearRight))
{
    logger.debug("Колесо вне поля");
    return false;
}

double hFL = field[frontLeft.y][frontLeft.x];
double hFR = field[frontRight.y][frontRight.x];
double hRL = field[rearLeft.y][rearLeft.x];
double hRR = field[rearRight.y][rearRight.x];
double hFront = (hFL + hFR) / 2.0;
double hRear  = (hRL + hRR) / 2.0;

double pitchAngle = algorithms::kinematics::calculateAngle(hFront, hRear, 2.0 * L);
    if (std::fabs(pitchAngle) > maxUpDownAngle)
    return false;
    
    double hLeft  = (hFL + hRL) / 2.0;
double hRight = (hFR + hRR) / 2.0;

double rollAngle  = algorithms::kinematics::calculateAngle(hLeft, hRight, 2.0 * W);
    if (std::fabs(rollAngle) > maxSideAngle)
    return false;
            

        // Проверка коллизий
        bool collisionCheck = isVehicleRadiusValid(center, binaryMap, vehicleRadius);
        logger.trace("  Коллизия: " + std::string(collisionCheck ? "нет" : "есть"));
        
        if (!collisionCheck) {
            logger.debug("  ! СТОЛКНОВЕНИЕ в точке (" + 
                       std::to_string(center.x) + "," + std::to_string(center.y) + ")");
            return false;
        }
    }
    
    logger.debug("--- Ребро проходимо ---");
    return true;
}
}
