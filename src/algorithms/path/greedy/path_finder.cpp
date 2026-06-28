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

std::vector<algorithms::geometry::Pixel> PathFinder::findPathGreedy(
    const algorithms::geometry::Pixel& start,
        const algorithms::geometry::Pixel& goal,
        const std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph)
{
    logger.info("[PathFinder::findPathGreedy] Начало поиска пути...");
    
    PathMetrics metrics;
    metrics.startTimer();
    
    // =============================
    //            Greedy
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
