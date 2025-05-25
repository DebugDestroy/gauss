#include "algorithms/path/a_star/path_finder.hpp"
#include "algorithms/path/a_star/heuristics.hpp"
#include "utils/hash.hpp"

#include <queue>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <string>

namespace algorithms::path::a_star {

PathFinder::PathFinder(core::Logger& lg) : logger(lg) {
    logger.trace("[PathFinder] Инициализация поисковика пути");
}

std::vector<algorithms::geometry::PointD> PathFinder::findPathAStar(
    const algorithms::geometry::PointD& start,
    const algorithms::geometry::PointD& goal,
    std::vector<algorithms::geometry::Edge> voronoiEdges,
    Graph& graph,
    const Conditions& conds,
    const std::vector<std::vector<double>>& binaryMap,
    const std::unique_ptr<algorithms::gauss::Pole>& elevationData)
{
    logger.info("[PathFinder::findPathAStar] Начало поиска пути...");
    
    auto localGraph = graph.buildGraphFromEdges(voronoiEdges, binaryMap, elevationData, conds);
    
    algorithms::geometry::PointD startNode = graph.findClosestVoronoiNode(start, localGraph);
    algorithms::geometry::PointD goalNode = graph.findClosestVoronoiNode(goal, localGraph);

    logger.debug("[PathFinder::findPathAStar] Ближайший к старту: (" + std::to_string(startNode.x) + ", " + std::to_string(startNode.y) + ")");
    logger.debug("[PathFinder::findPathAStar] Ближайший к цели: (" + std::to_string(goalNode.x) + ", " + std::to_string(goalNode.y) + ")");

    auto addTemporaryEdge = [&](const algorithms::geometry::PointD& a, const algorithms::geometry::PointD& b) -> bool {
        if (!conds.isEdgeNavigable(algorithms::geometry::Edge(a, b), elevationData, binaryMap) ||
            !conds.isVehicleRadiusValid(a, binaryMap) ||
            !conds.isVehicleRadiusValid(b, binaryMap)) {
            logger.warning("[PathFinder::findPathAStar] Временное ребро непригодно: (" +
                           std::to_string(a.x) + ", " + std::to_string(a.y) + ") -> (" +
                           std::to_string(b.x) + ", " + std::to_string(b.y) + ")");
            return false;
        }
        localGraph[a].push_back(b);
        localGraph[b].push_back(a);
        return true;
    };

    if (!addTemporaryEdge(start, startNode) || !addTemporaryEdge(goal, goalNode)) {
        logger.warning("[PathFinder::findPathAStar] Не удалось добавить стартовое или конечное ребро");
        return {};
    }

    std::priority_queue<AStarNode, std::vector<AStarNode>, std::greater<AStarNode>> openSet;
    std::unordered_map<algorithms::geometry::PointD, algorithms::geometry::PointD> cameFrom;
    std::unordered_map<algorithms::geometry::PointD, double> gScore;
    std::unordered_map<algorithms::geometry::PointD, double> fScore;

    gScore[start] = 0.0;
    fScore[start] = euclidean(start, goal);
    openSet.push({start, 0.0, fScore[start]});

    while (!openSet.empty()) {
        AStarNode current = openSet.top();
        openSet.pop();

        if (current.position == goal) {
            logger.info("[PathFinder::findPathAStar] Цель достигнута!");

            std::vector<algorithms::geometry::PointD> path;
            for (algorithms::geometry::PointD node = current.position; cameFrom.count(node); node = cameFrom[node])
                path.push_back(node);
            path.push_back(start);
            std::reverse(path.begin(), path.end());
            return path;
        }

        for (const algorithms::geometry::PointD& neighbor : localGraph[current.position]) {
            double tentativeG = gScore[current.position] + euclidean(current.position, neighbor);

            if (!gScore.count(neighbor) || tentativeG < gScore[neighbor]) {
                cameFrom[neighbor] = current.position;
                gScore[neighbor] = tentativeG;
                fScore[neighbor] = tentativeG + euclidean(neighbor, goal);
                openSet.push({neighbor, gScore[neighbor], fScore[neighbor]});
            }
        }
    }

    logger.warning("[PathFinder::findPathAStar] Путь не найден!");
    return {};
}

}
