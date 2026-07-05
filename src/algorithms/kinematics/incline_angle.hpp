#pragma once
#include "algorithms/geometry/geometry_structures.hpp"

#include <vector>

namespace algorithms::kinematics {

struct VehicleAngles {
    double sideAngle;   // наклон вбок
    double upDownAngle;    // наклон вперед/назад
};

// Вычисляет угол наклона между точками по высоте
double calculateAngle(double h1,
                      double h2,
                      double distance);

// Вычисляет VehicleAngles                       
VehicleAngles calculateVehicleAngles(
    const std::vector<std::vector<double>>& field,
    const algorithms::geometry::Pixel& center,
    const algorithms::geometry::PointD& dir,
    int vehicleRadius);

}
