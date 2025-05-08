#pragma once
#include <vector>
#include <queue> // для std::priority_queue
#include <unordered_map> // для std::unordered_map
#include <cmath> // для std::hypot, std::atan, M_PI
#include <algorithm> // для std::find_if
#include <string>
#include <sstream>
#include <array>

// Локальные заголовки
#include "core/Constants.hpp"  // Подключаем константы
#include "services/Geometry.hpp" // для PointD, Triangle, Edge

struct AStarNode {
    PointD position;       // Позиция узла (центр треугольника)
    const Triangle* tri;   // Связанный треугольник
    AStarNode* parent;     // Родительский узел
    double g;              // Стоимость пути от старта
    double h;              // Эвристическая оценка до цели
    double f;              // Общая стоимость: f = g + h

    AStarNode(PointD pos, const Triangle* triangle, AStarNode* p = nullptr, double cost = 0, double heuristic_val = 0)
    : position(pos), tri(triangle), parent(p), 
      g(cost), h(heuristic_val), f(g + h) {}

    // Для сравнения в priority_queue
    bool operator<(const AStarNode& other) const { return f > other.f; }
};
    

class PathFinder {
private:
    const Config& config;
    Logger& logger;

std::vector<PointD> reconstructPath(const AStarNode& endNode, 
                                   const PointD& start, 
                                   const PointD& goal) const {
    std::vector<PointD> path;
    const AStarNode* current = &endNode;
    while (current) {
        path.push_back(current->position);
        current = current->parent;
    }
    std::reverse(path.begin(), path.end());
    path.insert(path.begin(), start);
    path.push_back(goal);
    return path;
}

    double heuristic(const PointD& a, const PointD& b) {
        double dist = std::hypot(a.x - b.x, a.y - b.y);
        logger.trace(std::string("[PathFinder::heuristic] Расстояние между (") + 
                   std::to_string(a.x) + "," + std::to_string(a.y) + ") и (" +
                   std::to_string(b.x) + "," + std::to_string(b.y) + "): " + 
                   std::to_string(dist));
        return dist;
    }
    
    bool isVehicleRadiusValid(const PointD& pixel, 
                         const std::vector<std::vector<double>>& binaryMap,
                         const std::unique_ptr<Pole>& elevationData,
                         double sliceLevel) const {
    logger.trace(std::string("[PathFinder::isVehicleRadiusValid] Проверка радиуса в точке (") +
               std::to_string(pixel.x) + "," + std::to_string(pixel.y) + ")");
    
    const double x = pixel.x;
    const double y = pixel.y;
    const double currentHeight = elevationData->field[static_cast<int>(y)][static_cast<int>(x)];
    
    logger.debug(std::string("[PathFinder::isVehicleRadiusValid] Высота в точке: ") + std::to_string(currentHeight) +
               ", уровень среза: " + std::to_string(sliceLevel));

    const double heightDiff = std::fabs(currentHeight - sliceLevel);
    if (heightDiff >= config.vehicleRadius) {
        logger.debug(std::string("[PathFinder::isVehicleRadiusValid] Тележка ниже среза (") + 
                   std::to_string(heightDiff) + " >= " + std::to_string(config.vehicleRadius) + 
                   "), столкновений нет");
        return true;
    }

    const double effectiveRadius = std::sqrt(config.vehicleRadius * config.vehicleRadius - heightDiff * heightDiff);
    const int radiusPixels = static_cast<int>(std::ceil(effectiveRadius));
    
    logger.trace(std::string("[PathFinder::isVehicleRadiusValid] Эффективный радиус: ") + 
               std::to_string(effectiveRadius) + " пикселей (" + 
               std::to_string(radiusPixels) + " целых)");

    const double radiusSquared = effectiveRadius * effectiveRadius;
    size_t checkedPixels = 0;
    size_t collisionPixels = 0;

    for (int dy = -radiusPixels; dy <= radiusPixels; ++dy) {
        for (int dx = -radiusPixels; dx <= radiusPixels; ++dx) {
            if (dx*dx + dy*dy > radiusSquared) continue;
            
            const int nx = static_cast<int>(x) + dx;
            const int ny = static_cast<int>(y) + dy;
            checkedPixels++;
            
            if (nx >= 0 && ny >= 0 && 
                nx < static_cast<int>(binaryMap[0].size()) && 
                ny < static_cast<int>(binaryMap.size()) &&
                std::fabs(binaryMap[ny][nx] - Constants::WHITE) < Constants::EPSILON) {
                collisionPixels++;
                logger.trace(std::string("[PathFinder::isVehicleRadiusValid] Столкновение в (") +
                           std::to_string(nx) + "," + std::to_string(ny) + ")");
            }
        }
    }

    bool result = (collisionPixels == 0);
    logger.debug(std::string("[PathFinder::isVehicleRadiusValid] Результат проверки: ") + 
               std::string(result ? "проходимо" : "столкновение") + 
               ", проверено пикселей: " + std::to_string(checkedPixels) +
               ", столкновений: " + std::to_string(collisionPixels));
    
    return result;
}

bool isEdgeNavigable(const Edge& edge, 
                    const std::unique_ptr<Pole>& p,
                    const std::vector<std::vector<double>>& binaryMap,
                    double sliceLevel) const {
    const auto line = edge.bresenhamLine(edge.a, edge.b);
    const double r = config.vehicleRadius;
    
    // Вектор пути и ортогональный ему
    const PointD dir = {line.back().x - line.front().x, line.back().y - line.front().y};
    PointD perp = {-dir.y, dir.x};
    double perpLen = std::hypot(perp.x, perp.y);
    perp.x = perp.x / perpLen * r;
    perp.y = perp.y / perpLen * r;

    logger.debug("--- Начало проверки ребра ---");
    logger.debug("Параметры тележки: радиус=" + std::to_string(r) + 
                ", max_уклон=" + std::to_string(config.maxUpDownAngle) + 
                "°, max_крен=" + std::to_string(config.maxSideAngle) + "°");
    
    for (size_t i = 0; i < line.size(); ++i) {
        const PointD& center = line[i];
        double h_center = p->field[static_cast<int>(center.y)][static_cast<int>(center.x)];
        
        std::ostringstream pointHeader;
        pointHeader << "Точка [" << i << "/" << line.size()-1 << "] (" 
                   << center.x << "," << center.y << "): "
                   << "высота=" << h_center << ", slice=" << sliceLevel;
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
        bool collisionCheck = isVehicleRadiusValid(center, binaryMap, p, sliceLevel);
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

public:
    PathFinder(const Config& cfg, Logger& lg) : config(cfg), logger(lg) {
        logger.trace("[PathFinder] Инициализация поисковика пути");
    }

    std::vector<PointD> findPathAStar(const PointD& start, const PointD& goal, 
                                 const std::vector<Triangle>& triangles,
                                 const std::vector<std::vector<double>>& binaryMap,
                                 double slice, const std::unique_ptr<Pole>& p) {
    logger.info("[PathFinder::findPathAStar] Начало поиска пути");
    
    const Triangle* startTri = Triangle::findContainingTriangle(start, triangles);
    const Triangle* goalTri = Triangle::findContainingTriangle(goal, triangles);
    
    if (!startTri || !goalTri) {
        logger.error("Старт или цель вне триангуляции");
        return {};
    }

    // Проверка для точек в одном треугольнике
    if (startTri == goalTri) {
        Edge directEdge(start, goal);
        if (!isEdgeNavigable(directEdge, p, binaryMap, slice)) {
            logger.warning("Прямой путь в треугольнике непроходим");
            return {};
        }
        return {start, goal};
    }

    // Основной алгоритм A*
    std::priority_queue<AStarNode> openSet;
    std::unordered_map<const Triangle*, double> costSoFar;
    openSet.emplace(startTri->calculateCircumcenter(), startTri);
    costSoFar[startTri] = 0;

    while (!openSet.empty()) {
        AStarNode current = openSet.top();
        openSet.pop();

        if (current.tri == goalTri) {
            std::vector<PointD> path = reconstructPath(current, start, goal);
            logger.info("Путь найден. Длина: " + std::to_string(path.size()));
            return path;
        }

        for (const auto& neighbor : Triangle::getNeighbors(*current.tri, triangles)) {
            Edge edge(current.position, neighbor->calculateCircumcenter());
            
            if (!isEdgeNavigable(edge, p, binaryMap, slice)) {
                continue; // Непроходимые рёбра отсекаются
            }

            double newCost = current.g + edge.length();
            if (!costSoFar.count(neighbor) || newCost < costSoFar[neighbor]) {
                costSoFar[neighbor] = newCost;
                double priority = newCost + heuristic(neighbor->calculateCircumcenter(), goal);
                openSet.emplace(neighbor->calculateCircumcenter(), neighbor, new AStarNode(current), priority);
            }
        }
    }

    logger.warning("Путь не найден");
    return {};
}

};
