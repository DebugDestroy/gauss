#include "graph.hpp"
#include "utils/hash.hpp"
#include "algorithms/geometry/math.hpp"
#include "algorithms/path/common/collision.hpp"

#include <cmath>
#include <limits>
#include <string>
#include <algorithm> // для sort

namespace algorithms::path::common {

Graph::Graph(core::Logger& lg) : logger(lg) {}

// Перегрузка: принимает вектор Edge
std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>> Graph::buildGraphFromEdges(
    const std::vector<algorithms::geometry::Edge>& edges,
    const std::vector<std::vector<double>>& binaryMap,
    const std::vector<std::vector<double>>& field,
    const PathValidator& conds,
    int vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle)
{
    std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>> graph;

    for (const auto& edge : edges) {
        logger.trace(
        "[Graph] Continuous edge: (" +
        std::to_string(edge.a.x) + ", " +
        std::to_string(edge.a.y) + ") -> (" +
        std::to_string(edge.b.x) + ", " +
        std::to_string(edge.b.y) + ")");
        
        algorithms::geometry::Pixel a = algorithms::geometry::toPixel(edge.a);
        algorithms::geometry::Pixel b = algorithms::geometry::toPixel(edge.b);

        algorithms::geometry::PixelEdge pe(a,b);
        
        logger.trace(
        "[Graph] Discretized edge: (" +
        std::to_string(pe.a.x) + ", " +
        std::to_string(pe.a.y) + ") -> (" +
        std::to_string(pe.b.x) + ", " +
        std::to_string(pe.b.y) + ")");
        
        if (!conds.isEdgeNavigable(pe, field, binaryMap, vehicleRadius, maxSideAngle, maxUpDownAngle)) {
            logger.trace("[PathFinder::buildGraphFromEdges] Ребро непроходимо, пропуск: (" +
                         std::to_string(pe.a.x) + ", " + std::to_string(pe.a.y) + ") -> (" +
                         std::to_string(pe.b.x) + ", " + std::to_string(pe.b.y) + ")");
            continue;
        }

        graph[pe.a].push_back(pe.b);
        graph[pe.b].push_back(pe.a);
    }

    return graph;
}

void Graph::connectPointToGraph(
    std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph,
    const algorithms::geometry::Pixel& p,
    const std::vector<std::vector<double>>& binaryMap,
    const std::vector<std::vector<double>>& field,
    const PathValidator& conds,
    int vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle,
    ConnectMode mode,
    std::size_t nearestVerticesCount)
{
    if (!algorithms::path::common::isVehicleRadiusValid(p, binaryMap, vehicleRadius)) {
        logger.warning("[Graph::connectPointToGraph] Точка непригодна");
        return;
    }
    
    if (graph.contains(p)) {
        logger.warning("[Graph::connectPointToGraph] Точка уже присутствует в графе");
        return;
    }
    
    if (graph.empty()) {
    logger.warning("[Graph::connectPointToGraph] Граф пустой");
    return;
}
    
    std::vector<std::pair<double, algorithms::geometry::Pixel>> candidates;

    for (const auto& [node, _] : graph)
    {
        algorithms::geometry::PixelEdge e(p, node);

        if (!conds.isEdgeNavigable(e, field, binaryMap, vehicleRadius, maxSideAngle, maxUpDownAngle))
            continue;

        double dist = std::hypot(
        p.x - node.x,
        p.y - node.y);
    logger.debug(
        "[Graph::connectPointToGraph] Ребро принято, расстояние = " +
        std::to_string(dist));
    candidates.emplace_back(dist, node);
    }
    logger.debug(
    "[Graph::connectPointToGraph] Найдено кандидатов: " +
    std::to_string(candidates.size()));
     if (candidates.empty()) {
        logger.warning("[Graph::connectPointToGraph] Нет допустимых соединений");
        return;
    }
    
    std::sort(candidates.begin(),
          candidates.end(),
          [](const auto& a, const auto& b)
          {
              return a.first < b.first;
          });
          
    // Добавляем вершину
    graph.emplace(p, std::vector<algorithms::geometry::Pixel>{});

    size_t connectionsCount = 0;
    
    switch (mode)
    {
        case ConnectMode::All:
        {
            connectionsCount = candidates.size();
            break;
        }

        case ConnectMode::Nearest:
        {
            connectionsCount = 1;
            break;
        }

        case ConnectMode::NearestK:
        {
 if (nearestVerticesCount <= 0)
    {
        logger.warning(
            "[Graph::connectPointToGraph] Некорректное число ближайших вершин: " +
            std::to_string(nearestVerticesCount) +
            ". Используется 1");

        connectionsCount = std::min<size_t>(1, candidates.size());
    }
    else
    {
        connectionsCount =
            std::min(nearestVerticesCount,
                     candidates.size());

        if (connectionsCount < nearestVerticesCount)
        {
            logger.warning(
                "[Graph::connectPointToGraph] Запрошено " +
                std::to_string(nearestVerticesCount) +
                " вершин, но доступно только " +
                std::to_string(connectionsCount)+
            ". Используется " +
                std::to_string(connectionsCount));
        }
    }
    break;
    }
    default:
    logger.warning("[Graph::connectPointToGraph] Неизвестный режим");
    return;
    }
    for (size_t i = 0; i < connectionsCount; ++i)
    {
        const auto& node = candidates[i].second;
        graph[p].push_back(node);
        graph[node].push_back(p);
    }

    logger.debug(
    "[Graph::connectPointToGraph] Точка (" +
    std::to_string(p.x) + ", " +
    std::to_string(p.y) +
    ") подключена к " +
    std::to_string(connectionsCount) +
    " вершинам");

    return;
}
}
