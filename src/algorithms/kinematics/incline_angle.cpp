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
}
