#pragma once

#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::path::a_star {

// Эвристическая евклидова функция расстояния
double euclidean(const algorithms::geometry::PointD& a, const algorithms::geometry::PointD& b);

}
