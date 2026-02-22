#pragma once

#include <unordered_map>
#include <vector>
#include <memory>

#include "algorithms/geometry/geometry_structures.hpp" // для PointD, Triangle, Edge
#include "algorithms/path/common/conditions.hpp"              // Conditions
#include "algorithms/gauss/pole.hpp"

namespace algorithms::path::common {

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

// Добавляет точку в граф и добавляет проходимые ребра от этой точки к другим вершинам
bool connectPointToGraph(
    std::unordered_map<algorithms::geometry::PointD,
                       std::vector<algorithms::geometry::PointD>>& graph,
    const algorithms::geometry::PointD& point,
    const std::vector<std::vector<double>>& binaryMap,
    const std::unique_ptr<algorithms::gauss::Pole>& elevationData,
    const Conditions& conds);
    
// Ищет ближайшую точку Вороного (узел графа) к произвольной точке
algorithms::geometry::PointD findClosestVoronoiNode(
    const algorithms::geometry::PointD& point,
    const std::unordered_map<algorithms::geometry::PointD, std::vector<algorithms::geometry::PointD>>& graph) const;

};
}
