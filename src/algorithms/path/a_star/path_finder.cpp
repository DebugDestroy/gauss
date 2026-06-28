#include "algorithms/path/a_star/path_finder.hpp"
#include "algorithms/path/common/heuristics.hpp"
#include "utils/hash.hpp"
#include "algorithms/path/common/path_metrics.hpp"

#include <queue>
#include <cmath>
#include <algorithm>
#include <string>
#include <unordered_set>

namespace algorithms::path::a_star {

PathFinder::PathFinder(core::Logger& lg) : logger(lg) {
    logger.trace("[PathFinder] Инициализация поисковика пути A*");
}

std::vector<algorithms::geometry::Pixel> PathFinder::findPathAStar(
    const algorithms::geometry::Pixel& start,
    const algorithms::geometry::Pixel& goal,
    const  std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph)
{
    logger.info("[PathFinder::findPathAStar] Начало поиска пути...");
    
    PathMetrics metrics;
    metrics.startTimer();
    
    // =============================
    //              A*
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

    //ОБЪВЛЕНИЕ ПЕРЕМЕННЫХ
    std::priority_queue<AStarNode, std::vector<AStarNode>, std::greater<AStarNode>> openSet;
    std::unordered_set<algorithms::geometry::Pixel> closedSet;
    std::unordered_map<algorithms::geometry::Pixel, algorithms::geometry::Pixel> cameFrom;
    std::unordered_map<algorithms::geometry::Pixel, double> gScore;
    std::unordered_map<algorithms::geometry::Pixel, double> fScore;

    gScore[start] = 0.0;
    fScore[start] = algorithms::path::common::euclidean(start, goal);
    openSet.push({start, 0.0, fScore[start]});

    while (!openSet.empty()) {
        AStarNode current = openSet.top();
        openSet.pop();
        
        logger.debug("[ASTAR] Current: (" +
             std::to_string(current.position.x) + ", " +
             std::to_string(current.position.y) + ")");
             
        if (closedSet.count(current.position))
           continue;
           
        metrics.expandedNodes++; 
        
        if (current.position == goal) {
            logger.info("[PathFinder::findPathAStar] Цель достигнута!");

            std::vector<algorithms::geometry::Pixel> path;
            for (algorithms::geometry::Pixel node = current.position; cameFrom.count(node); node = cameFrom[node])
                path.push_back(node);
            path.push_back(start);
            std::reverse(path.begin(), path.end());
            
metrics.computeFromPath(path);
metrics.finishAndLog(logger, "A*");

            return path;
        }
        
        closedSet.insert(current.position);

        for (const algorithms::geometry::Pixel& neighbor : graph.at(current.position)) {
            double tentativeG = gScore[current.position] + std::hypot(current.position.x - neighbor.x, current.position.y - neighbor.y);

            if (closedSet.count(neighbor))
            continue;
            
            if (!gScore.count(neighbor) || tentativeG < gScore[neighbor]) {
                cameFrom[neighbor] = current.position;
                gScore[neighbor] = tentativeG;
                fScore[neighbor] = tentativeG + algorithms::path::common::euclidean(neighbor, goal);
                openSet.push({neighbor, gScore[neighbor], fScore[neighbor]});
            }
        }
    }

    logger.warning("[PathFinder::findPathAStar] Путь не найден!");
    
metrics.finishAndLog(logger, "A* (Failed)");

    return {};
}

}
