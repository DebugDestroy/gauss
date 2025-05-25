#pragma once

#include <memory>
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/gauss/pole.hpp"

namespace algorithms::kinematics {

// Вычисляет угол наклона между центром и колесом по данным высот
double calculateWheelAngle(const algorithms::geometry::PointD& center, 
                           const algorithms::geometry::PointD& wheel,
                           const std::unique_ptr<algorithms::gauss::Pole>& p);

}