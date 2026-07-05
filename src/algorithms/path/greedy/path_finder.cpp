#include "algorithms/path/greedy/path_finder.hpp"
#include "utils/hash.hpp"
#include "algorithms/path/common/path_metrics.hpp"
#include "algorithms/path/common/heuristics.hpp"

#include <algorithm>
#include <unordered_set>

namespace algorithms::path::greedy {

PathFinder::PathFinder(core::Logger& lg)
    : logger(lg)
{
    logger.trace("[PathFinder] Инициализация Greedy поисковика пути");
}

std::vector<algorithms::geometry::Pixel> PathFinder::findPathGreedyGraph(
    const algorithms::geometry::Pixel& start,
    const algorithms::geometry::Pixel& goal,
    const std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph,
    algorithms::path::PathMetrics& metrics)
{
    logger.info("[PathFinder::findPathGreedy] Начало поиска пути...");
    
    // =============================
    //            Greedy
    // =============================
    if (graph.empty()) {
        logger.warning(
            "[PathFinder::findPathGreedy] Граф пуст");
        return {};
    }

    if (!graph.contains(start)) {
        logger.warning(
            "[PathFinder::findPathGreedy] Старт отсутствует в графе");
        return {};
    }

    if (!graph.contains(goal)) {
        logger.warning(
            "[PathFinder::findPathGreedy] Финиш отсутствует в графе");
        return {};
    }
    
    std::vector<algorithms::geometry::Pixel> path;
    algorithms::geometry::Pixel current = start;
    path.push_back(current);

    std::unordered_set<algorithms::geometry::Pixel> visited;
    visited.insert(current);

    while (!(current == goal)) {
        const auto& neighbors = graph.at(current);
        std::vector<algorithms::geometry::Pixel> candidates;
        metrics.expandedNodes++; 

for (const auto& n : neighbors) {
    if (!visited.count(n)) {
        candidates.push_back(n);
    }
}

logger.debug("[Greedy] Current: (" +
             std::to_string(current.x) + ", " +
             std::to_string(current.y) + ")");

if (candidates.empty()) {
    logger.warning("[Greedy] Все соседи посещены, тупик");
    
    return path;
}

auto nextIt = std::min_element(
    candidates.begin(),
    candidates.end(),
    [&](const auto& a, const auto& b) {
        return algorithms::path::common::euclidean(a, goal) <
               algorithms::path::common::euclidean(b, goal);
    }
);

        current = *nextIt;
        path.push_back(current);
        visited.insert(current);
    }

    logger.info("[PathFinder::findPathGreedy] Цель достигнута!");
    metrics.pathFound = true;
    return path;
}

std::vector<algorithms::path::common::GridCell>
PathFinder::findPathGreedyGrid(
    const algorithms::path::common::GridCell& startCell,
    const algorithms::path::common::GridCell& endCell,
    const algorithms::path::common::Grid& grid,
    algorithms::path::PathMetrics& metrics)
{
    logger.info("[PathFinder::findPathGreedyGrid] Start Greedy on grid");

    if (grid.cells.empty()) {
        logger.warning("[Greedy Grid] Grid is empty");
        return {};
    }

    auto index = [&](int row, int col) {
        return row * grid.cols + col;
    };

    int startIdx = index(startCell.row, startCell.col);
    int goalIdx  = index(endCell.row, endCell.col);

    if (!grid.cells[startIdx].traversable ||
        !grid.cells[goalIdx].traversable)
    {
        logger.warning("[Greedy Grid] Start or goal cell is blocked");
        return {};
    }

    std::vector<algorithms::path::common::GridCell> path;

    algorithms::path::common::GridCell current = startCell;
    path.push_back(current);

    std::vector<bool> visited(grid.cells.size(), false);
    visited[startIdx] = true;

    while (true)
    {
        int curIdx = index(current.row, current.col);

        if (curIdx == goalIdx)
            break;

        metrics.expandedNodes++;

        auto neighbours = getNeighbours(grid, current);

        std::vector<algorithms::path::common::GridCell> candidates;

        for (const algorithms::path::common::GridCell& n : neighbours)
        {
            if (!n.traversable)
                continue;
                
            int nIdx = index(n.row, n.col);

            if (visited[nIdx])
                continue;

            candidates.push_back(n);
        }

        logger.debug("[Greedy] Current: (" +
                     std::to_string(current.row) + ", " +
                     std::to_string(current.col) + ")");

        if (candidates.empty()) {
            logger.warning("[Greedy Grid] Dead end");
            return path;
        }

        auto nextIt = std::min_element(
            candidates.begin(),
            candidates.end(),
            [&](const algorithms::path::common::GridCell& a, const algorithms::path::common::GridCell& b)
            {
                return algorithms::path::common::euclidean(a, endCell) <
                       algorithms::path::common::euclidean(b, endCell);
            });

        current = *nextIt;
        visited[index(current.row, current.col)] = true;
        path.push_back(current);
    }

    logger.info("[Greedy Grid] Goal reached");
    metrics.pathFound = true;

    return path;
}

}
