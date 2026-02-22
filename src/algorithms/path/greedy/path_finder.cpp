#include "algorithms/path/greedy/path_finder.hpp"
#include "utils/hash.hpp"
#include "algorithms/path/common/path_metrics.hpp"

#include <algorithm>
#include <unordered_set>

namespace algorithms::path::greedy {

PathFinder::PathFinder(core::Logger& lg)
    : logger(lg)
{
    logger.trace("[PathFinder] Инициализация Greedy поисковика пути");
}

std::vector<algorithms::geometry::PointD> PathFinder::findPathGreedy(
    const algorithms::geometry::PointD& start,
    const algorithms::geometry::PointD& goal,
    std::vector<algorithms::geometry::Edge> voronoiEdges,
    algorithms::path::common::Graph& graph,
    const algorithms::path::common::Conditions& conds,
    const std::vector<std::vector<double>>& binaryMap,
    const std::unique_ptr<algorithms::gauss::Pole>& elevationData)
{
    logger.info("[PathFinder::findPathGreedy] Начало поиска пути...");
    
    PathMetrics metrics;
    metrics.startTimer();
    
        // Строим локальный граф
auto localGraph = graph.buildGraphFromEdges(
    voronoiEdges, binaryMap, elevationData, conds);

logger.debug("[PathFinder::findPathGreedy] Подключаем старт и цель ко всем видимым вершинам");

if (!graph.connectPointToGraph(localGraph, start, binaryMap, elevationData, conds)) { // Добавили старт
    logger.warning("[PathFinder::findPathGreedy] Старт не подключен в граф!");
    
metrics.finishAndLog(logger, "Greedy (Failed)");

    return {};
}

if (!graph.connectPointToGraph(localGraph, goal, binaryMap, elevationData, conds)) { // Добавили финиш
    logger.warning("[PathFinder::findPathGreedy] Финиш не подключен в граф!");

metrics.finishAndLog(logger, "Greedy (Failed)");

    return {};
}
    
    // =============================
    //            Greedy
    // =============================

    std::vector<algorithms::geometry::PointD> path;
    algorithms::geometry::PointD current = start;
    path.push_back(current);

    std::unordered_set<algorithms::geometry::PointD> visited;
    visited.insert(current);

    while (!(current == goal)) {
        const auto& neighbors = localGraph[current];
        std::vector<algorithms::geometry::PointD> candidates;
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
    
metrics.computeFromPath(path);
metrics.finishAndLog(logger, "Greedy (Failed)");
    
    return path;
}

auto nextIt = std::min_element(
    candidates.begin(),
    candidates.end(),
    [&](const auto& a, const auto& b) {
        return std::hypot(a.x - goal.x, a.y - goal.y) <
               std::hypot(b.x - goal.x, b.y - goal.y);
    }
);

        current = *nextIt;
        path.push_back(current);
        visited.insert(current);
    }

    logger.info("[PathFinder::findPathGreedy] Цель достигнута!");
    
metrics.computeFromPath(path);
metrics.finishAndLog(logger, "Greedy");
    
    return path;
}

}
