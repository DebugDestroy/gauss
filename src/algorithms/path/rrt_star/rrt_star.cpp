#include "algorithms/path/rrt_star/rrt_star.hpp"
#include "algorithms/geometry/math.hpp"
#include "algorithms/path/common/collision.hpp"

#include <algorithm>
#include <utility>

namespace algorithms::path::rrt_star {

PathFinder::PathFinder(core::Logger& lg,
           std::mt19937& generator)
    : logger(lg),
      gen(generator)
{}

double PathFinder::calculateRadius(std::size_t nodes,
                       double gammaConstant,
                       double maxRadius)
{
    if (nodes < 2)
        return maxRadius;
        
    double r =
        gammaConstant * std::sqrt(
            std::log(static_cast<double>(nodes)) /
            static_cast<double>(nodes)
        );

    return std::min(r,maxRadius);
}

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
    const std::vector<RRTStarNode>& tree,
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

std::vector<int> PathFinder::nearNodes(
    const std::vector<RRTStarNode>& tree,
    const algorithms::geometry::PointD& point,
    double radius) const
{
    std::vector<int> result;
    int nearest = nearestNode(tree, point);
    result.push_back(nearest);
    
    const double radiusSquared = radius * radius;

    for (std::size_t i = 0; i < tree.size(); ++i)
    {
        if (static_cast<int>(i) == nearest)
            continue;
            
        const double dx = tree[i].point.x - point.x;
        const double dy = tree[i].point.y - point.y;

        const double distSquared = dx * dx + dy * dy;

        if (distSquared <= radiusSquared)
        {
            result.push_back(static_cast<int>(i));
        }
    }

    return result;
}

int PathFinder::chooseParent(
    const std::vector<RRTStarNode>& tree,
    const std::vector<int>& neighbors,
    const algorithms::geometry::PointD& newPoint,
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
    const algorithms::path::common::PathValidator& conds)
{
    int bestParent = -1;

    double bestCost =
        std::numeric_limits<double>::max();


    for (const int index : neighbors)
    {
        const auto& node =
            tree[static_cast<std::size_t>(index)];


        if (!conds.isEdgeValidContinuous(
                gaussi,
                node.point,
                newPoint,
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


        const double newCost = node.cost + algorithms::geometry::distance(node.point, newPoint);


        if (newCost < bestCost)
        {
            bestCost = newCost;
            bestParent = index;
        }
    }


    return bestParent;
}

void PathFinder::updateChildrenCost(
    std::vector<RRTStarNode>& tree,
    int parentIndex)
{
    const auto& parent =
        tree[static_cast<std::size_t>(parentIndex)];


    for (std::size_t i = 0; i < tree.size(); ++i)
    {
        if (static_cast<int>(i) == parentIndex)
            continue;


        auto& child = tree[i];


        if (child.parent != parentIndex)
            continue;


        const double distance = algorithms::geometry::distance(child.point, parent.point);


        child.cost =
            parent.cost + distance;


        updateChildrenCost(
            tree,
            static_cast<int>(i)
        );
    }
}

void PathFinder::rewire(
    std::vector<RRTStarNode>& tree,
    const std::vector<int>& neighbors,
    int newIndex,
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
    const algorithms::path::common::PathValidator& conds)
{
    const auto newNode =
        tree[static_cast<std::size_t>(newIndex)];


    for (const int index : neighbors)
    {
        // Нельзя переподключать самого себя
        if (index == newIndex)
            continue;


        auto& node =
            tree[static_cast<std::size_t>(index)];


        const double distance = algorithms::geometry::distance(newNode.point, node.point);

        const double newCost =
            newNode.cost + distance;


        if (newCost >= node.cost)
            continue;


        // Проверяем новое ребро
        if (!conds.isEdgeValidContinuous(
                gaussi,
                newNode.point,
                node.point,
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


        // Улучшили путь
        node.parent = newIndex;
        node.cost = newCost;
        
        // Обновили стоимость детишек
        updateChildrenCost(
            tree,
            index
        );
    }
}

std::vector<algorithms::geometry::PointD>
PathFinder::restorePath(
    const std::vector<RRTStarNode>& tree,
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

RRTStarResult PathFinder::findPathRRTStar(
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
    double maxFindRadius,
    double gammaConstant,
    double goalRadius,
    double goalBias,
    const algorithms::path::common::PathValidator& conds,
    PathMetrics& metrics)
{
    std::vector<RRTStarNode> tree;
    RRTStarResult result;
    int bestGoalParent = -1;
    double bestGoalCost = std::numeric_limits<double>::infinity();

    logger.debug("[RRTStar] Starting search");
    logger.debug("Start = (" + std::to_string(start.x) + ", " + std::to_string(start.y) + ")");
    logger.debug("Goal  = (" + std::to_string(goal.x) + ", " + std::to_string(goal.y) + ")");
    logger.debug("Iterations = " + std::to_string(maxIterations));
    logger.debug("Step = " + std::to_string(step));
    logger.debug("Find Radius = " + std::to_string(maxFindRadius));
    logger.debug("Gamma constant = " + std::to_string(gammaConstant));
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
        logger.error("[RRTStar] Invalid start position");
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
        logger.error("[RRTStar] Invalid goal position");
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
    
    tree.push_back({start, -1, 0.0});

    for (size_t iter = 0; iter < maxIterations; ++iter)
    {
        if (iter % 10000 == 0)
            logger.debug("iter = " + std::to_string(iter));
            
        // 1. Генерируем случайную цель
        algorithms::geometry::PointD qRand = randomPoint(goal, goalBias);

        // 2. Ищем ближайшую вершину дерева
        int nearest = nearestNode(tree, qRand);

        // 3. Делаем шаг в сторону случайной точки
        algorithms::geometry::PointD qNew =
            steer(tree[static_cast<std::size_t>(nearest)].point,
                  qRand,
                  step);
        
        // Если точка недостижима, то ищем другую
        if (!algorithms::path::common::checkPointContinuous(
                gaussi,
                qNew,
                fieldWidth,
                fieldHeight,
                heightThreshold,
                vehicleRadius,
                interpolationCollision,
                interpolationAngle))
            continue;
    
        // 4. Пересчитываем find radius
        double radius = calculateRadius(tree.size(), gammaConstant, maxFindRadius);
        
        // 5. Ищем соседей в круге радиуса radius (и всегда добавляем ближнего)
        auto neighbors = nearNodes(
            tree,
            qNew,           
            radius
        );
                
        // 6. Ищем лучшего соседа
       int parent = chooseParent(
                        tree,
                        neighbors,
                        qNew,
                        gaussi,
                        fieldWidth,
                        fieldHeight,
                        heightThreshold,
                        vehicleRadius,
                        maxSideAngle,
                        maxUpDownAngle,
                        interpolationEdge,
                        interpolationCollision,
                        interpolationAngle,
                        conds);
                     
         // Если нет соседей то ищем новую точку
         if (parent == -1)
             continue;
        
        // 7. Добавляем вершину
        // Стоимость новой вершины
        double cost =
            tree[static_cast<std::size_t>(parent)].cost +
            algorithms::geometry::distance(
                tree[static_cast<std::size_t>(parent)].point,
                qNew
            );
            
        tree.push_back({qNew, parent, cost});
        metrics.expandedNodes++; 
        
        int newIndex = static_cast<int>(tree.size()) - 1;
        
        // 8. Обновляем дерево
        rewire(
            tree,
            neighbors,
            newIndex,
            gaussi,
            fieldWidth,
            fieldHeight,
            heightThreshold,
            vehicleRadius,
            maxSideAngle,
            maxUpDownAngle,
            interpolationEdge,
            interpolationCollision,
            interpolationAngle,
            conds);
        
        double distGoal = algorithms::geometry::distance(qNew, goal);
        
        // 9. Проверяем достижение цели
        if (distGoal <= goalRadius)
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
                double goalCost = tree[static_cast<std::size_t>(newIndex)].cost + distGoal;
                    
                if(goalCost < bestGoalCost)
                {
                    bestGoalParent = newIndex;
                    bestGoalCost = goalCost;
                    
                    logger.debug(
                        "[RRT*] Better path found. "
                        "iter=" + std::to_string(iter) +
                        ", cost=" + std::to_string(bestGoalCost)
                    );
                }
            }
        }
    }
    
    if(bestGoalParent != -1)
    {
        tree.push_back({goal, bestGoalParent, bestGoalCost});
        int goalIndex = static_cast<int>(tree.size()) - 1;
        result.path = restorePath(tree, goalIndex);
        metrics.pathFound = true;
    }
    else
    {
        metrics.pathFound = false;
    }
    
    result.tree = std::move(tree);
    return result;
}
}
