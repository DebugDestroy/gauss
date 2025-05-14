#pragma once
#include <vector>
#include <queue> // для std::priority_queue
#include <unordered_map> // для std::unordered_map
#include <cmath> // для std::hypot, std::atan, M_PI
#include <algorithm> // для std::find_if
#include <string>
#include <sstream>
#include <array>
#include <map>

// Локальные заголовки
#include "core/Constants.hpp"  // Подключаем константы
#include "services/Geometry.hpp" // для PointD, Triangle, Edge

struct AStarNode {
    PointD position;
    double gScore;
    double fScore;

    bool operator>(const AStarNode& other) const {
        return fScore > other.fScore;
    }
}; 

class PathFinder {
private:
    const Config& config;
    Logger& logger;

    double heuristic(const PointD& a, const PointD& b) {
        double dist = std::hypot(a.x - b.x, a.y - b.y);
        logger.trace(std::string("[PathFinder::heuristic] Расстояние между (") + 
                   std::to_string(a.x) + "," + std::to_string(a.y) + ") и (" +
                   std::to_string(b.x) + "," + std::to_string(b.y) + "): " + 
                   std::to_string(dist));
        return dist;
    }
    
    bool isVehicleRadiusValid(const PointD& pixel, 
                         const std::vector<std::vector<double>>& binaryMap) const
{
    logger.trace("[PathFinder::isVehicleRadiusValid] Проверка радиуса в точке (" +
                 std::to_string(pixel.x) + "," + std::to_string(pixel.y) + ")");

    const double x = pixel.x;
    const double y = pixel.y;

    const double effectiveRadius = config.vehicleRadius;
    const int radiusPixels = static_cast<int>(std::ceil(effectiveRadius));
    const double radiusSquared = effectiveRadius * effectiveRadius;

    size_t checkedPixels = 0;
    size_t collisionPixels = 0;

    for (int dy = -radiusPixels; dy <= radiusPixels; ++dy) {
        for (int dx = -radiusPixels; dx <= radiusPixels; ++dx) {
            if (dx * dx + dy * dy > radiusSquared) continue;

            const int nx = static_cast<int>(x) + dx;
            const int ny = static_cast<int>(y) + dy;
            ++checkedPixels;

            if (nx >= 0 && ny >= 0 &&
                ny < static_cast<int>(binaryMap.size()) &&
                nx < static_cast<int>(binaryMap[0].size()) &&
                std::fabs(binaryMap[ny][nx] - Constants::WHITE) < Constants::EPSILON) {

                ++collisionPixels;
                logger.trace("[PathFinder::isVehicleRadiusValid] Столкновение в (" +
                             std::to_string(nx) + "," + std::to_string(ny) + ")");
            }
        }
    }

    bool result = (collisionPixels == 0);
    logger.debug("[PathFinder::isVehicleRadiusValid] Результат: " + 
                 std::string(result ? "проходимо" : "столкновение") +
                 ", проверено пикселей: " + std::to_string(checkedPixels) + 
                 ", столкновений: " + std::to_string(collisionPixels));

    return result;
}

bool isEdgeNavigable(const Edge& edge, 
                    const std::unique_ptr<Pole>& p,
                    const std::vector<std::vector<double>>& binaryMap) const {
    const auto line = edge.bresenhamLine(edge.a, edge.b);
    const double r = config.vehicleRadius;
    
    // Вектор пути и ортогональный ему
    const PointD dir = {line.back().x - line.front().x, line.back().y - line.front().y};
    PointD perp = {-dir.y, dir.x};
    double perpLen = std::hypot(perp.x, perp.y);
    
    // Защита от деления на ноль
        if (perpLen < Constants::EPSILON) {
            logger.error("Нулевая длина перпендикуляра");
            return false;
        }
        
    perp.x = perp.x / perpLen * r;
    perp.y = perp.y / perpLen * r;

    logger.debug("--- Начало проверки ребра ---");
    logger.debug("Параметры тележки: радиус=" + std::to_string(r) + 
                ", max_уклон=" + std::to_string(config.maxUpDownAngle) + 
                "°, max_крен=" + std::to_string(config.maxSideAngle) + "°");
    // Проверка на минимальную длину ребра
        if (line.size() < 2) {
            logger.error("Слишком короткое ребро для проверки");
            return false;
        }
        logger.trace("Длина ребра: " + std::to_string(line.size()) + " точек");
    
    for (size_t i = 0; i < line.size(); ++i) {
        const PointD& center = line[i];
        double h_center; // Высота центра тележки
        // Проверка границ массива
            if (center.y < 0 || center.x < 0 || 
                center.y >= p->field.size() || 
                center.x >= p->field[0].size()) {
                logger.error("Точка за границами поля: (" + 
                           std::to_string(center.x) + "," + 
                           std::to_string(center.y) + ")");
                return false;
            }
            
        h_center = p->field[static_cast<int>(center.y)][static_cast<int>(center.x)];
        
        std::ostringstream pointHeader;
        pointHeader << "Точка [" << i << "/" << line.size()-1 << "] (" 
                   << center.x << "," << center.y << "): "
                   << "высота=" << h_center;
        logger.trace(pointHeader.str());

        // Проверка переднего колеса
        if (i + r < line.size()) {
            const PointD& frontWheel = line[i + static_cast<size_t>(r)];
            double h_front = p->field[static_cast<int>(frontWheel.y)][static_cast<int>(frontWheel.x)];
            double frontAngle = calculateWheelAngle(center, frontWheel, p);
            
            logger.trace("  Переднее колесо (" + 
                        std::to_string(frontWheel.x) + "," + std::to_string(frontWheel.y) + 
                        "): угол=" + std::to_string(frontAngle) + 
                        "°, высота=" + std::to_string(h_front));
            
            if (std::fabs(frontAngle) > config.maxUpDownAngle) {
                logger.debug("  ! ПРЕВЫШЕНИЕ: передний угол " + std::to_string(frontAngle) + 
                           "° > допустимого " + std::to_string(config.maxUpDownAngle) + "°");
                return false;
            }
        }

        // Боковые колеса
        PointD leftWheel = {std::round(center.x + perp.x), std::round(center.y + perp.y)};
        PointD rightWheel = {std::round(center.x - perp.x), std::round(center.y - perp.y)};

        double leftAngle = 0, rightAngle = 0;
        if (leftWheel.x >= 0 && leftWheel.y >= 0 && 
            leftWheel.x < p->field[0].size() && leftWheel.y < p->field.size()) {
            double h_left = p->field[static_cast<int>(leftWheel.y)][static_cast<int>(leftWheel.x)];
            leftAngle = calculateWheelAngle(center, leftWheel, p);
            logger.trace("  Левое колесо (" + 
                        std::to_string(leftWheel.x) + "," + std::to_string(leftWheel.y) + 
                        "): угол=" + std::to_string(leftAngle) + 
                        "°, высота=" + std::to_string(h_left));
        }

        if (rightWheel.x >= 0 && rightWheel.y >= 0 && 
            rightWheel.x < p->field[0].size() && rightWheel.y < p->field.size()) {
            double h_right = p->field[static_cast<int>(rightWheel.y)][static_cast<int>(rightWheel.x)];
            rightAngle = calculateWheelAngle(center, rightWheel, p);
            logger.trace("  Правое колесо (" + 
                        std::to_string(rightWheel.x) + "," + std::to_string(rightWheel.y) + 
                        "): угол=" + std::to_string(rightAngle) + 
                        "°, высота=" + std::to_string(h_right));
        }

        if (std::fabs(leftAngle) > config.maxSideAngle || 
            std::fabs(rightAngle) > config.maxSideAngle) {
            logger.debug("  ! ПРЕВЫШЕНИЕ: боковые углы L=" + std::to_string(leftAngle) + 
                       "° R=" + std::to_string(rightAngle) + 
                       " > допустимого " + std::to_string(config.maxSideAngle) + "°");
            return false;
        }

        // Проверка коллизий
        bool collisionCheck = isVehicleRadiusValid(center, binaryMap);
        logger.trace("  Коллизия: " + std::string(collisionCheck ? "нет" : "есть"));
        
        if (!collisionCheck) {
            logger.debug("  ! СТОЛКНОВЕНИЕ в точке (" + 
                       std::to_string(center.x) + "," + std::to_string(center.y) + ")");
            return false;
        }
    }
    
    logger.debug("--- Ребро проходимо ---");
    return true;
}

double calculateWheelAngle(const PointD& center, 
                          const PointD& wheel,
                          const std::unique_ptr<Pole>& p) const {
    // Получаем высоты в точках
    double h_center = p->field[static_cast<int>(center.y)][static_cast<int>(center.x)];
    double h_wheel = p->field[static_cast<int>(wheel.y)][static_cast<int>(wheel.x)];
    
    if (wheel.x < 0 || wheel.y < 0 || 
    wheel.x >= p->field[0].size() || 
    wheel.y >= p->field.size()) {
    return std::numeric_limits<double>::max(); // Помечаем как непроходимое
}
    // Расчет угла в градусах
    return std::atan2(h_wheel - h_center, 
                     std::hypot(wheel.x - center.x, wheel.y - center.y)) * 180.0 / M_PI;
}

PointD findClosestVoronoiNode(const PointD& point, const std::unordered_map<PointD, std::vector<PointD>>& graph) {
    double bestDist = std::numeric_limits<double>::max();
    PointD closest = point;

    for (const auto& [node, _] : graph) {
        double dist = std::hypot(point.x - node.x, point.y - node.y);
        if (dist < bestDist) {
            bestDist = dist;
            closest = node;
        }
    }
    return closest;
}

// Перегрузка: принимает вектор Edge
std::unordered_map<PointD, std::vector<PointD>> buildGraphFromEdges(
    const std::vector<Edge>& edges,
    const std::vector<std::vector<double>>& binaryMap,
    const std::unique_ptr<Pole>& elevationData)
{
    std::unordered_map<PointD, std::vector<PointD>> graph;

    for (const auto& edge : edges) {
        if (!isEdgeNavigable(edge, elevationData, binaryMap)) {
            logger.trace("[PathFinder::buildGraphFromEdges] Ребро непроходимо, пропуск: (" +
                         std::to_string(edge.a.x) + ", " + std::to_string(edge.a.y) + ") -> (" +
                         std::to_string(edge.b.x) + ", " + std::to_string(edge.b.y) + ")");
            continue;
        }

        if (!isVehicleRadiusValid(edge.a, binaryMap)) {
            logger.trace("[PathFinder::buildGraphFromEdges] Точка edge.a непригодна для радиуса машины: (" +
                         std::to_string(edge.a.x) + ", " + std::to_string(edge.a.y) + ")");
            continue;
        }

        if (!isVehicleRadiusValid(edge.b, binaryMap)) {
            logger.trace("[PathFinder::buildGraphFromEdges] Точка edge.b непригодна для радиуса машины: (" +
                         std::to_string(edge.b.x) + ", " + std::to_string(edge.b.y) + ")");
            continue;
        }

        graph[edge.a].push_back(edge.b);
        graph[edge.b].push_back(edge.a);
    }

    return graph;
}


public:
    PathFinder(const Config& cfg, Logger& lg) : config(cfg), logger(lg) {
        logger.trace("[PathFinder] Инициализация поисковика пути");
    }

    std::vector<PointD> findPathAStar(
    const PointD& start,
    const PointD& goal,
    std::vector<Edge>& voronoiEdges,
    const std::vector<std::vector<double>>& binaryMap,
    const std::unique_ptr<Pole>& elevationData)
{
    logger.info("[PathFinder::findPathAStar] Начало поиска пути...");

    auto graph = buildGraphFromEdges(voronoiEdges, binaryMap, elevationData);
    logger.debug("[PathFinder::findPathAStar] Граф построен, количество узлов: " + std::to_string(graph.size()));

    PointD startNode = findClosestVoronoiNode(start, graph);
    PointD goalNode = findClosestVoronoiNode(goal, graph);
    logger.debug("[PathFinder::findPathAStar] Ближайший к старту: (" + std::to_string(startNode.x) + ", " + std::to_string(startNode.y) + ")");
    logger.debug("[PathFinder::findPathAStar] Ближайший к цели: (" + std::to_string(goalNode.x) + ", " + std::to_string(goalNode.y) + ")");

    // Добавим временные рёбра с проверками проходимости
auto addTemporaryEdge = [&](const PointD& a, const PointD& b) -> bool {
    if (!isEdgeNavigable(Edge(a, b), elevationData, binaryMap)) {
        logger.warning("[PathFinder::findPathAStar] Временное ребро непроходимо: (" +
                       std::to_string(a.x) + ", " + std::to_string(a.y) + ") -> (" +
                       std::to_string(b.x) + ", " + std::to_string(b.y) + ")");
        return false;
    }
    if (!isVehicleRadiusValid(a, binaryMap) || !isVehicleRadiusValid(b, binaryMap)) {
        logger.warning("[PathFinder::findPathAStar] Временное ребро вне допустимого радиуса: (" +
                       std::to_string(a.x) + ", " + std::to_string(a.y) + ") -> (" +
                       std::to_string(b.x) + ", " + std::to_string(b.y) + ")");
        return false;
    }
    graph[a].push_back(b);
    graph[b].push_back(a);
    return true;
};

if (!addTemporaryEdge(start, startNode) || !addTemporaryEdge(goal, goalNode)) {
    logger.warning("[PathFinder::findPathAStar] Не удалось добавить стартовое или конечное ребро");
    return {};
}

    std::priority_queue<AStarNode, std::vector<AStarNode>, std::greater<AStarNode>> openSet;
    std::unordered_map<PointD, PointD> cameFrom;
    std::unordered_map<PointD, double> gScore;
    std::unordered_map<PointD, double> fScore;

    gScore[start] = 0.0;
    fScore[start] = heuristic(start, goal);
    openSet.push({start, 0.0, fScore[start]});

    while (!openSet.empty()) {
        AStarNode current = openSet.top();
        openSet.pop();

        logger.trace("[PathFinder::findPathAStar] Обработка узла: (" +
                     std::to_string(current.position.x) + ", " + std::to_string(current.position.y) + ")");

        if (current.position == goal) {
            logger.info("[PathFinder::findPathAStar] Цель достигнута!");

            std::vector<PointD> path;
            PointD node = current.position;
            while (cameFrom.count(node)) {
                path.push_back(node);
                node = cameFrom[node];
            }
            path.push_back(start);
            std::reverse(path.begin(), path.end());

            return path;
        }

        for (const auto& neighbor : graph[current.position]) {
            Edge edge(current.position, neighbor);

            double tentativeG = gScore[current.position] + heuristic(current.position, neighbor);

            if (!gScore.count(neighbor) || tentativeG < gScore[neighbor]) {
                cameFrom[neighbor] = current.position;
                gScore[neighbor] = tentativeG;
                fScore[neighbor] = tentativeG + heuristic(neighbor, goal);
                openSet.push({neighbor, gScore[neighbor], fScore[neighbor]});
                logger.trace("[PathFinder::findPathAStar] Добавление в очередь: (" +
                             std::to_string(neighbor.x) + ", " + std::to_string(neighbor.y) + "), f = " +
                             std::to_string(fScore[neighbor]));
            }
        }
    }

    logger.warning("[PathFinder::findPathAStar] Путь не найден!");
    return {};
}
};
