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

std::vector<algorithms::geometry::Pixel> PathFinder::findPathDijkstra(
    const algorithms::geometry::Pixel& start,
        const algorithms::geometry::Pixel& goal,
        const std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph)
{
    logger.info("[PathFinder::findPathDijkstra] Начало поиска пути...");
    
    PathMetrics metrics;
    metrics.startTimer();

    // =============================
    //         Dijkstra
    // =============================
if (graph.empty()) {
    logger.warning(
        "[PathFinder::findPathAStar] Граф пуст");
    return {};
}

if (!graph.contains(start)) {
    logger.warning(
        "[PathFinder::findPathAStar] Старт отсутствует в графе");
    return {};
}

if (!graph.contains(goal)) {
    logger.warning(
        "[PathFinder::findPathAStar] Финиш отсутствует в графе");
    return {};
}
    std::priority_queue<
        DijkstraNode,
        std::vector<DijkstraNode>,
        std::greater<DijkstraNode>
    > openSet;

    std::unordered_set<algorithms::geometry::Pixel> closedSet;
    std::unordered_map<algorithms::geometry::Pixel,
                       algorithms::geometry::Pixel> cameFrom;
    std::unordered_map<algorithms::geometry::Pixel,
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

            std::vector<algorithms::geometry::Pixel> path;

            for (algorithms::geometry::Pixel node = current.position;
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
             graph.at(current.position))
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
