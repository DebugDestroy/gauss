#pragma once

#include <unordered_map>
#include <vector>
#include <memory>

#include "algorithms/geometry/geometry_structures.hpp" // для PointD, Triangle, Edge...
#include "algorithms/path/common/conditions.hpp"              // Conditions

namespace algorithms::path::common {

enum class ConnectMode {
    Nearest,     // подключить к ближайшей вершине
    NearestK,    // подключить к k ближайшим вершинам
    All          // подключить ко всем достижимым вершинам
};

class Graph {
private:
    core::Logger& logger;

public:
    explicit Graph(core::Logger& lg);
        
// Возвращает граф достижимости на основе навигабельных рёбер
std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>> buildGraphFromEdges(
    const std::vector<algorithms::geometry::Edge>& edges,
    const std::vector<std::vector<double>>& binaryMap,
    const std::vector<std::vector<double>>& field,
    const Conditions& conds,
    int vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle);

// Добавляет точку в граф и добавляет проходимые ребра от этой точки к другим вершинам
void connectPointToGraph(
    std::unordered_map<algorithms::geometry::Pixel,
                       std::vector<algorithms::geometry::Pixel>>& graph,
    const algorithms::geometry::Pixel& point,
    const std::vector<std::vector<double>>& binaryMap,
    const std::vector<std::vector<double>>& field,
    const Conditions& conds,
    int vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle,
    ConnectMode mode,
    int nearestVerticesCount = 1);
};
}
