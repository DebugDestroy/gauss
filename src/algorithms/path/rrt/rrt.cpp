#include "algorithms/path/rrt/rrt.hpp"
#include "algorithms/geometry/math.hpp"
#include "algorithms/path/common/collision.hpp"

#include <algorithm>
#include <utility>

namespace algorithms::path::rrt {

PathFinder::PathFinder(core::Logger& lg,
           std::mt19937& generator)
    : logger(lg),
      gen(generator)
{}

algorithms::geometry::PointD PathFinder::steer(
    const algorithms::geometry::PointD& from,
    const algorithms::geometry::PointD& to,
    double step)
{
    algorithms::geometry::PointD dir = to - from;

    double len = std::hypot(dir.x, dir.y);

    if (len < step)
        return to;

    dir /= len;

    return {
        from.x + dir.x * step,
        from.y + dir.y * step
    };
}

algorithms::geometry::PointD PathFinder::randomPoint(const algorithms::geometry::PointD& goal,
                               double goalBias)
{
    if (probDist(gen) < goalBias)
        return goal;

    return {xDist(gen), yDist(gen)};
}

int PathFinder::nearestNode(
    const std::vector<RRTNode>& tree,
    const algorithms::geometry::PointD& point) const
{
    int bestIndex = 0;
    double bestDist2 = std::numeric_limits<double>::max();

    for (std::size_t i = 0; i < tree.size(); ++i)
    {
        const double dx = tree[i].point.x - point.x;
        const double dy = tree[i].point.y - point.y;

        const double dist2 = dx * dx + dy * dy;

        if (dist2 < bestDist2)
        {
            bestDist2 = dist2;
            bestIndex = static_cast<int>(i);
        }
    }

    return bestIndex;
}

std::vector<algorithms::geometry::PointD>
PathFinder::restorePath(
    const std::vector<RRTNode>& tree,
    int goalIndex) const
{
    std::vector<algorithms::geometry::PointD> path;

    int current = goalIndex;

    while (current != -1)
    {
        const auto& p = tree[static_cast<std::size_t>(current)].point;

        path.push_back(p);

        current = tree[static_cast<std::size_t>(current)].parent;
    }

    std::reverse(path.begin(), path.end());

    return path;
}

RRTResult PathFinder::findPathRRT(
    const algorithms::geometry::PointD& start,
    const algorithms::geometry::PointD& goal,
    const std::vector<algorithms::gauss::Gaus>& gaussi,
    int fieldWidth,
    int fieldHeight,
    double heightThreshold,
    double vehicleRadius,
    double maxSideAngle,
    double maxUpDownAngle,
    double interpolationEdge,
    double interpolationCollision,
    double interpolationAngle,
    std::size_t maxIterations,
    double step,
    double goalRadius,
    double goalBias,
    const algorithms::path::common::PathValidator& conds,
    PathMetrics& metrics)
{
    std::vector<RRTNode> tree;
    RRTResult result;
    
    logger.debug("[RRT] Starting search");
    logger.debug("Start = (" + std::to_string(start.x) + ", " + std::to_string(start.y) + ")");
    logger.debug("Goal  = (" + std::to_string(goal.x) + ", " + std::to_string(goal.y) + ")");
    logger.debug("Iterations = " + std::to_string(maxIterations));
    logger.debug("Step = " + std::to_string(step));
    logger.debug("Goal radius = " + std::to_string(goalRadius));
    
    if (!algorithms::path::common::checkPointContinuous(
            gaussi,
            start,
            fieldWidth,
            fieldHeight,
            heightThreshold,
            vehicleRadius,
            interpolationCollision,
            interpolationAngle))
    {
        logger.error("[RRT] Invalid start position");
        metrics.pathFound = false;
        return result;
    }

    if (!algorithms::path::common::checkPointContinuous(
            gaussi,
            goal,
            fieldWidth,
            fieldHeight,
            heightThreshold,
            vehicleRadius,
            interpolationCollision,
            interpolationAngle))
    {
        logger.error("[RRT] Invalid goal position");
        metrics.pathFound = false;
        return result;
    }
    
    xDist.param(
        std::uniform_real_distribution<double>::param_type(
            0.0,
            static_cast<double>(fieldWidth)
        )
    );

    yDist.param(
        std::uniform_real_distribution<double>::param_type(
            0.0,
            static_cast<double>(fieldHeight)
        )
    );
    
    tree.push_back({start, -1});

    for (size_t iter = 0; iter < maxIterations; ++iter)
    {
        // 1. Генерируем случайную цель
        algorithms::geometry::PointD qRand = randomPoint(goal, goalBias);

        // 2. Ищем ближайшую вершину дерева
        int nearest = nearestNode(tree, qRand);

        // 3. Делаем шаг в сторону случайной точки
        algorithms::geometry::PointD qNew =
            steer(tree[static_cast<std::size_t>(nearest)].point,
                  qRand,
                  step);
                  
        // 4. Проверяем новое ребро
        if (!conds.isEdgeValidContinuous(
                gaussi,
                tree[static_cast<std::size_t>(nearest)].point,
                qNew,
                fieldWidth,
                fieldHeight,
                heightThreshold,
                vehicleRadius,
                maxSideAngle,
                maxUpDownAngle,
                interpolationEdge,
                interpolationCollision,
                interpolationAngle))
        {
            continue;
        }

        // 5. Добавляем вершину
        tree.push_back({qNew, nearest});
        metrics.expandedNodes++; 
        
        int newIndex = static_cast<int>(tree.size()) - 1;

        // 6. Проверяем достижение цели
        if (distance(qNew, goal) <= goalRadius)
        {
            if (conds.isEdgeValidContinuous(
                gaussi,
                qNew,
                goal,
                fieldWidth,
                fieldHeight,
                heightThreshold,
                vehicleRadius,
                maxSideAngle,
                maxUpDownAngle,
                interpolationEdge,
                interpolationCollision,
                interpolationAngle))
            {
            tree.push_back({goal, newIndex});
            metrics.pathFound = true;
            result.path = restorePath(tree, static_cast<int>(tree.size()) - 1);
            result.tree = std::move(tree);
            return result;
            }
        }
    }
    
    result.tree = std::move(tree);
    metrics.pathFound = false;
    return result;
}
}
