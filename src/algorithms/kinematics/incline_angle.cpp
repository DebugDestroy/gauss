#include "algorithms/kinematics/incline_angle.hpp"
#include "algorithms/geometry/math.hpp"

#include <cmath>
#include <numbers>

namespace algorithms::kinematics {
double calculateAngle(double h1,
                      double h2,
                      double distance)
{
    return std::atan2(h1 - h2, distance)
           * 180.0 / std::numbers::pi;
}

VehicleAngles calculateVehicleAngles(
    const std::vector<std::vector<double>>& field,
    const algorithms::geometry::Pixel& center,
    const algorithms::geometry::PointD& dir,
    int vehicleRadius)
{
    const int width = static_cast<int>(field[0].size());
    const int height = static_cast<int>(field.size());
    
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
        
// Проверка границ массива
            if (!algorithms::geometry::isPointInsideField(frontLeft, width, height) ||
    !algorithms::geometry::isPointInsideField(frontRight, width, height) ||
    !algorithms::geometry::isPointInsideField(rearLeft, width, height) ||
    !algorithms::geometry::isPointInsideField(rearRight, width, height))
{
    return {0.0, 0.0}; // Вышли за границы значит угол нулевой
}

double hFL = field[static_cast<std::size_t>(frontLeft.y)]
                  [static_cast<std::size_t>(frontLeft.x)];

double hFR = field[static_cast<std::size_t>(frontRight.y)]
                  [static_cast<std::size_t>(frontRight.x)];

double hRL = field[static_cast<std::size_t>(rearLeft.y)]
                  [static_cast<std::size_t>(rearLeft.x)];

double hRR = field[static_cast<std::size_t>(rearRight.y)]
                  [static_cast<std::size_t>(rearRight.x)];

double hFront = (hFL + hFR) / 2.0;
double hRear  = (hRL + hRR) / 2.0;

double upDownAngle = calculateAngle(hFront, hRear, 2.0 * L);

double hLeft  = (hFL + hRL) / 2.0;
double hRight = (hFR + hRR) / 2.0;

double sideAngle  = calculateAngle(hLeft, hRight, 2.0 * W);

return {sideAngle, upDownAngle};
}

VehicleAngles calculateVehicleAnglesContinuous(
    const std::vector<algorithms::gauss::Gaus>& gaussi,
    const algorithms::geometry::PointD& center,
    const algorithms::geometry::PointD& direction, // единичное направление
    double vehicleRadius)
{

    algorithms::geometry::PointD perp = {-direction.y, direction.x};

    const double L = vehicleRadius / std::sqrt(2.0);
    const double W = vehicleRadius / std::sqrt(2.0);

    // точки корпуса (в непрерывной системе — тоже PointD)
    algorithms::geometry::PointD frontLeft  = {
        center.x + L * direction.x + W * perp.x,
        center.y + L * direction.y + W * perp.y
    };

    algorithms::geometry::PointD frontRight = {
        center.x + L * direction.x - W * perp.x,
        center.y + L * direction.y - W * perp.y
    };

    algorithms::geometry::PointD rearLeft   = {
        center.x - L * direction.x + W * perp.x,
        center.y - L * direction.y + W * perp.y
    };

    algorithms::geometry::PointD rearRight  = {
        center.x - L * direction.x - W * perp.x,
        center.y - L * direction.y - W * perp.y
    };
    
    // высоты через гауссы
    auto hFL = algorithms::gauss::GaussBuilder::heightAt(frontLeft, gaussi);
    auto hFR = algorithms::gauss::GaussBuilder::heightAt(frontRight, gaussi);
    auto hRL = algorithms::gauss::GaussBuilder::heightAt(rearLeft, gaussi);
    auto hRR = algorithms::gauss::GaussBuilder::heightAt(rearRight, gaussi);

    double hFront = (hFL + hFR) / 2.0;
    double hRear  = (hRL + hRR) / 2.0;

    double upDownAngle = calculateAngle(hFront, hRear, 2.0 * L);

    double hLeft  = (hFL + hRL) / 2.0;
    double hRight = (hFR + hRR) / 2.0;

    double sideAngle = calculateAngle(hLeft, hRight, 2.0 * W);

    return {sideAngle, upDownAngle};
}

}
