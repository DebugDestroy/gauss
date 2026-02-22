#include "graph.hpp"
#include "utils/hash.hpp"

#include <cmath>
#include <limits>
#include <string>

namespace algorithms::path::common {

Graph::Graph(core::Logger& lg) : logger(lg) {}

// Перегрузка: принимает вектор Edge
std::unordered_map<algorithms::geometry::PointD, std::vector<algorithms::geometry::PointD>> Graph::buildGraphFromEdges(
    const std::vector<algorithms::geometry::Edge>& edges,
    const std::vector<std::vector<double>>& binaryMap,
    const std::unique_ptr<algorithms::gauss::Pole>& elevationData,
    const Conditions& conds)
{
    std::unordered_map<algorithms::geometry::PointD, std::vector<algorithms::geometry::PointD>> graph;

    for (const auto& edge : edges) {
        if (!conds.isEdgeNavigable(edge, elevationData, binaryMap)) {
            logger.trace("[PathFinder::buildGraphFromEdges] Ребро непроходимо, пропуск: (" +
                         std::to_string(edge.a.x) + ", " + std::to_string(edge.a.y) + ") -> (" +
                         std::to_string(edge.b.x) + ", " + std::to_string(edge.b.y) + ")");
            continue;
        }

        if (!conds.isVehicleRadiusValid(edge.a, binaryMap)) {
            logger.trace("[PathFinder::buildGraphFromEdges] Точка edge.a непригодна для радиуса машины: (" +
                         std::to_string(edge.a.x) + ", " + std::to_string(edge.a.y) + ")");
            continue;
        }

        if (!conds.isVehicleRadiusValid(edge.b, binaryMap)) {
            logger.trace("[PathFinder::buildGraphFromEdges] Точка edge.b непригодна для радиуса машины: (" +
                         std::to_string(edge.b.x) + ", " + std::to_string(edge.b.y) + ")");
            continue;
        }

        graph[edge.a].push_back(edge.b);
        graph[edge.b].push_back(edge.a);
    }

    return graph;
}

bool Graph::connectPointToGraph(
    std::unordered_map<algorithms::geometry::PointD, std::vector<algorithms::geometry::PointD>>& graph,
    const algorithms::geometry::PointD& p,
    const std::vector<std::vector<double>>& binaryMap,
    const std::unique_ptr<gauss::Pole>& elevationData,
    const Conditions& conds)
{
    if (!conds.isVehicleRadiusValid(p, binaryMap)) {
        logger.warning("[Graph::connectPointToGraph] Точка непригодна");
        return false;
    }

    // Добавляем точку заранее, чтобы не было вставки во время итерации
    graph.emplace(p, std::vector<algorithms::geometry::PointD>{});

    bool connected = false;

    for (const auto& [node, _] : graph)
    {
        if (node == p)
            continue;

        geometry::Edge e(p, node);

        if (!conds.isEdgeNavigable(e, elevationData, binaryMap))
            continue;

        graph[p].push_back(node);
        graph[node].push_back(p);

        connected = true;
    }

    if (!connected) {
        logger.warning("[Graph::connectPointToGraph] Нет допустимых соединений");
    }

    return connected;
}

 algorithms::geometry::PointD Graph::findClosestVoronoiNode(const algorithms::geometry::PointD& point, const std::unordered_map<algorithms::geometry::PointD, std::vector<algorithms::geometry::PointD>>& graph) const 
 {
    double bestDist = std::numeric_limits<double>::max();
    algorithms::geometry::PointD closest = point;

    for (const auto& [node, _] : graph) {
        double dist = std::hypot(point.x - node.x, point.y - node.y);
        if (dist < bestDist) {
            bestDist = dist;
            closest = node;
        }
    }
    return closest;
 }
}
