#pragma once

#include <unordered_map>
#include <vector>
#include <memory>

#include "algorithms/geometry/geometry_structures.hpp" // для PointD, Triangle, Edge
#include "algorithms/path/a_star/conditions.hpp"              // Conditions
#include "algorithms/gauss/pole.hpp"

namespace algorithms::path::a_star {

class Graph {
private:
    core::Logger& logger;

public:
    explicit Graph(core::Logger& lg);
        
// Возвращает граф достижимости на основе навигабельных рёбер
std::unordered_map<algorithms::geometry::PointD, std::vector<algorithms::geometry::PointD>> buildGraphFromEdges(
    const std::vector<algorithms::geometry::Edge>& edges,
    const std::vector<std::vector<double>>& binaryMap,
    const std::unique_ptr<algorithms::gauss::Pole>& elevationData,
    const Conditions& conds);

// Ищет ближайшую точку Вороного (узел графа) к произвольной точке
algorithms::geometry::PointD findClosestVoronoiNode(
    const algorithms::geometry::PointD& point,
    const std::unordered_map<algorithms::geometry::PointD, std::vector<algorithms::geometry::PointD>>& graph) const;

};
}
