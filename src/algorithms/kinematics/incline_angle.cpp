#include "algorithms/kinematics/incline_angle.hpp"
#include <cmath>
#include <limits>

namespace algorithms::kinematics {
double calculateWheelAngle(const algorithms::geometry::PointD& center, 
                          const algorithms::geometry::PointD& wheel,
                          const std::unique_ptr<algorithms::gauss::Pole>& p) {
    // Проверка границ                      
    if (wheel.x < 0 || wheel.y < 0 || 
    wheel.x >= p->field[0].size() || 
    wheel.y >= p->field.size()) {
    return std::numeric_limits<double>::max(); // Помечаем как непроходимое
}

    // Получаем высоты в точках
    double h_center = p->field[static_cast<int>(center.y)][static_cast<int>(center.x)];
    double h_wheel = p->field[static_cast<int>(wheel.y)][static_cast<int>(wheel.x)];
   
    // Расчет угла в градусах
    return std::atan2(h_wheel - h_center, 
                     std::hypot(wheel.x - center.x, wheel.y - center.y)) * 180.0 / M_PI;
}
}
