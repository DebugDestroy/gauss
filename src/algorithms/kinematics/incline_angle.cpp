#include "algorithms/kinematics/incline_angle.hpp"
#include <cmath>

namespace algorithms::kinematics {
double calculateAngle(double h1,
                      double h2,
                      double distance)
{
    return std::atan2(h1 - h2, distance)
           * 180.0 / M_PI;
}

VehicleAngles calculateVehicleAngles(
    const std::vector<std::vector<double>>& field,
    const algorithms::geometry::Pixel& center,
    const algorithms::geometry::PointD& dir,
    int vehicleRadius)
{
    algorithms::geometry::PointD direction = dir;
    double dirLen = std::hypot(direction.x, direction.y);
    
    // Защита от деления на ноль
    if (dirLen < core::EPSILON) {
        return {0.0, 0.0}; // Угол нулевой так как точки очень близко
    }
    
    direction.x /= dirLen;
    direction.y /= dirLen;
    algorithms::geometry::PointD perp = {-direction.y, direction.x}; // Единичный перпендикуляр
    const double L = vehicleRadius / std::sqrt(2.0);
    const double W = vehicleRadius / std::sqrt(2.0);
    
    algorithms::geometry::Pixel frontLeft = {
    static_cast<int>(std::round(center.x + L * direction.x + W * perp.x)),
    static_cast<int>(std::round(center.y + L * direction.y + W * perp.y))
};

algorithms::geometry::Pixel frontRight = {
    static_cast<int>(std::round(center.x + L * direction.x - W * perp.x)),
    static_cast<int>(std::round(center.y + L * direction.y - W * perp.y))
};

algorithms::geometry::Pixel rearLeft = {
    static_cast<int>(std::round(center.x - L * direction.x + W * perp.x)),
    static_cast<int>(std::round(center.y - L * direction.y + W * perp.y))
};

algorithms::geometry::Pixel rearRight = {
    static_cast<int>(std::round(center.x - L * direction.x - W * perp.x)),
    static_cast<int>(std::round(center.y - L * direction.y - W * perp.y))
};

auto isInside = [&](const algorithms::geometry::Pixel& p) {
    return p.x >= 0 &&
           p.y >= 0 &&
           p.x < field[0].size() &&
           p.y < field.size();
};
        
// Проверка границ массива
            if (!isInside(frontLeft) ||
    !isInside(frontRight) ||
    !isInside(rearLeft) ||
    !isInside(rearRight))
{
    return {0.0, 0.0}; // Вышли за границы значит угол нулевой
}

double hFL = field[frontLeft.y][frontLeft.x];
double hFR = field[frontRight.y][frontRight.x];
double hRL = field[rearLeft.y][rearLeft.x];
double hRR = field[rearRight.y][rearRight.x];

double hFront = (hFL + hFR) / 2.0;
double hRear  = (hRL + hRR) / 2.0;

double upDownAngle = calculateAngle(hFront, hRear, 2.0 * L);

double hLeft  = (hFL + hRL) / 2.0;
double hRight = (hFR + hRR) / 2.0;

double sideAngle  = calculateAngle(hLeft, hRight, 2.0 * W);

return {sideAngle, upDownAngle};
}

}
