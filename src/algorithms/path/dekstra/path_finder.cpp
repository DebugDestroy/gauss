#include "algorithms/path/dekstra/path_finder.hpp"
#include "utils/hash.hpp"
#include "algorithms/path/common/path_metrics.hpp"

#include <queue>
#include <cmath>
#include <algorithm>
#include <string>
#include <unordered_set>

namespace algorithms::path::dekstra {

PathFinder::PathFinder(core::Logger& lg)
    : logger(lg)
{
    logger.trace("[PathFinder] Инициализация поисковика пути Dijkstra");
}

std::vector<algorithms::geometry::PointD> PathFinder::findPathDijkstra(
    const algorithms::geometry::PointD& start,
    const algorithms::geometry::PointD& goal,
    std::vector<algorithms::geometry::Edge> voronoiEdges,
    algorithms::path::common::Graph& graph,
    const algorithms::path::common::Conditions& conds,
    const std::vector<std::vector<double>>& binaryMap,
    const std::unique_ptr<algorithms::gauss::Pole>& elevationData)
{
    logger.info("[PathFinder::findPathDijkstra] Начало поиска пути...");
    
    PathMetrics metrics;
    metrics.startTimer();
    
    // Строим локальный граф
auto localGraph = graph.buildGraphFromEdges(
    voronoiEdges, binaryMap, elevationData, conds);

logger.debug("[PathFinder::findPathDijkstra] Подключаем старт и цель ко всем видимым вершинам");

if (!graph.connectPointToGraph(localGraph, start, binaryMap, elevationData, conds)) { // Добавили старт
    logger.warning("[PathFinder::findPathDijkstra] Старт не подключен в граф!");
    
metrics.finishAndLog(logger, "Dijkstra (Failed)");

    return {};
}

if (!graph.connectPointToGraph(localGraph, goal, binaryMap, elevationData, conds)) { // Добавили финиш
    logger.warning("[PathFinder::findPathDijkstra] Финиш не подключен в граф!");
    
metrics.finishAndLog(logger, "Dijkstra (Failed)");

    return {};
}

    // =============================
    //         Dijkstra
    // =============================

    std::priority_queue<
        DijkstraNode,
        std::vector<DijkstraNode>,
        std::greater<DijkstraNode>
    > openSet;

    std::unordered_set<algorithms::geometry::PointD> closedSet;
    std::unordered_map<algorithms::geometry::PointD,
                       algorithms::geometry::PointD> cameFrom;
    std::unordered_map<algorithms::geometry::PointD,
                       double> gScore;

    gScore[start] = 0.0;
    openSet.push({start, 0.0});

    while (!openSet.empty()) {

        DijkstraNode current = openSet.top();
        openSet.pop();
        
        logger.debug("[DEKSTRA] Current: (" +
             std::to_string(current.position.x) + ", " +
             std::to_string(current.position.y) + ")");
             
        if (closedSet.count(current.position))
            continue;
            
        metrics.expandedNodes++; 

        if (current.position == goal) {

            logger.info(
                "[PathFinder::findPathDijkstra] Цель достигнута!");

            std::vector<algorithms::geometry::PointD> path;

            for (algorithms::geometry::PointD node = current.position;
                 cameFrom.count(node);
                 node = cameFrom[node])
            {
                path.push_back(node);
            }

            path.push_back(start);
            std::reverse(path.begin(), path.end());
            
metrics.computeFromPath(path);
metrics.finishAndLog(logger, "Dijkstra");
    
            return path;
        }

        closedSet.insert(current.position);

        for (const auto& neighbor :
             localGraph[current.position])
        {
            if (closedSet.count(neighbor))
                continue;

            double tentativeG =
                gScore[current.position] + std::hypot(current.position.x - neighbor.x, current.position.y - neighbor.y); 

            if (!gScore.count(neighbor) ||
                tentativeG < gScore[neighbor])
            {
                cameFrom[neighbor] = current.position;
                gScore[neighbor] = tentativeG;
                openSet.push({neighbor, gScore[neighbor]});
            }
        }
    }

    logger.warning(
        "[PathFinder::findPathDijkstra] Путь не найден!");

metrics.finishAndLog(logger, "Dijkstra (Failed)");

    return {};
}

}
