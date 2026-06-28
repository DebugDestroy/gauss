#pragma once

#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::path::common {

// Эвристическая евклидова функция расстояния
double euclidean(const algorithms::geometry::Pixel& a, const algorithms::geometry::Pixel& b);

}
