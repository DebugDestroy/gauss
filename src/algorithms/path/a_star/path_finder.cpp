#include "algorithms/path/a_star/path_finder.hpp"
#include "algorithms/path/common/heuristics.hpp"
#include "utils/hash.hpp"

#include <queue>
#include <cmath>
#include <algorithm>
#include <string>
#include <unordered_set>

namespace algorithms::path::a_star {

PathFinder::PathFinder(core::Logger& lg) : logger(lg) {
    logger.trace("[PathFinder] Инициализация поисковика пути A*");
}

std::vector<algorithms::geometry::Pixel> PathFinder::findPathAStarGraph(
    const algorithms::geometry::Pixel& start,
    const algorithms::geometry::Pixel& goal,
    const  std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph,
    algorithms::path::PathMetrics& metrics)
{
    logger.info("[PathFinder::findPathAStar] Начало поиска пути...");
    
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
    std::priority_queue<NodeGraph, std::vector<NodeGraph>, std::greater<NodeGraph>> openSet;
    std::unordered_set<algorithms::geometry::Pixel> closedSet;
    std::unordered_map<algorithms::geometry::Pixel, algorithms::geometry::Pixel> cameFrom;
    std::unordered_map<algorithms::geometry::Pixel, double> gScore;
    std::unordered_map<algorithms::geometry::Pixel, double> fScore;

    gScore[start] = 0.0;
    fScore[start] = algorithms::path::common::euclidean(start, goal);
    openSet.push({start, 0.0, fScore[start]});

    while (!openSet.empty()) {
        NodeGraph current = openSet.top();
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
            metrics.pathFound = true;
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

    return {};
}

std::vector<algorithms::path::common::GridCell>
PathFinder::findPathAStarGrid(
        const algorithms::path::common::GridCell& startCell,
        const algorithms::path::common::GridCell& endCell,
        const algorithms::path::common::Grid& grid,
        algorithms::path::PathMetrics& metrics)
{
    logger.debug("[PathFinder::findPathAStarGrid] Start A* on grid");

    if (grid.cells.empty()) {
        logger.warning("[A* Grid] Grid is empty");
        return {};
    }
    
    auto index = [&](int row, int col) {
        return static_cast<std::size_t>(row * grid.cols + col);
    };
    
    std::size_t startIdx = index(startCell.row, startCell.col);
    std::size_t goalIdx  = index(endCell.row, endCell.col);
    
    if (!grid.cells[startIdx].traversable ||
    !grid.cells[goalIdx].traversable)
    {
        logger.warning("[A* Grid] Start or goal cell is blocked");
        return {};
    }

    std::priority_queue<
        NodeGrid,
        std::vector<NodeGrid>,
        std::greater<NodeGrid>
    > openSet;
    
    const std::size_t nodeCount = grid.cells.size();
    std::vector<bool> closedSet(nodeCount, false);
    std::vector<int> cameFrom(nodeCount, -1);
    std::vector<double> gScore(nodeCount, std::numeric_limits<double>::infinity());
    std::vector<double> fScore(nodeCount, std::numeric_limits<double>::infinity());

    gScore[startIdx] = 0.0;
    fScore[startIdx] = algorithms::path::common::euclidean(startCell, endCell);

    openSet.push({startIdx, 0.0, fScore[startIdx]});

    while (!openSet.empty()) {

        NodeGrid current = openSet.top();
        openSet.pop();

        std::size_t curIdx = current.idx;
        const algorithms::path::common::GridCell& curCell = grid.cells[curIdx];

        if (closedSet[curIdx])
            continue;

        metrics.expandedNodes++;

        if (curIdx == goalIdx) {

            logger.debug("[A* Grid] Goal reached");

            std::vector<algorithms::path::common::GridCell> path;

            for (int at = static_cast<int>(curIdx);
                 at != -1;
                 at = cameFrom[static_cast<std::size_t>(at)])
            {
                  path.push_back(grid.cells[static_cast<std::size_t>(at)]);
            }

            std::reverse(path.begin(), path.end());

            metrics.pathFound = true;
            return path;
        }

        closedSet[curIdx] = true;

        auto neighbours = getNeighbours(grid, curCell);

        for (const algorithms::path::common::GridCell& neighbor : neighbours)
        {
            if (!neighbor.traversable)
                continue;
        
            std::size_t nIdx = index(neighbor.row, neighbor.col);

            if (closedSet[nIdx])
                continue;

            double stepCost =
                (neighbor.row != curCell.row && neighbor.col != curCell.col)
                    ? std::sqrt(2.0)
                    : 1.0;

            double tentativeG = gScore[curIdx] + stepCost;

            if (tentativeG < gScore[nIdx]) {

                cameFrom[nIdx] = static_cast<int>(curIdx);
                gScore[nIdx] = tentativeG;
                fScore[nIdx] = tentativeG +
                    algorithms::path::common::euclidean(neighbor, endCell);

                openSet.push({
                    nIdx,
                    tentativeG,
                    fScore[nIdx]
                });
            }
        }

    }

    logger.warning("[A* Grid] Path not found");
    return {};
}

}
