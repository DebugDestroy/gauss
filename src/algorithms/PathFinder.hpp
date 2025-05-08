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

    void logPoint(const std::string& prefix, const PointD& p) const {
        logger.debug(prefix + " (" + std::to_string(p.x) + ", " + std::to_string(p.y) + ")");
    }

    void logTriangle(const std::string& prefix, const Triangle* tri) const {
        if (tri) {
            std::ostringstream oss;
            oss << prefix << " ["
                << "A(" << tri->a.x << "," << tri->a.y << "), "
                << "B(" << tri->b.x << "," << tri->b.y << "), "
                << "C(" << tri->c.x << "," << tri->c.y << ")]";
            logger.debug(oss.str());
        } else {
            logger.debug(prefix + " [не найден]");
        }
    }

    void logEdge(const Edge& edge) const {
        std::ostringstream oss;
        oss << "Ребро: (" << edge.a.x << "," << edge.a.y << ")-(" 
            << edge.b.x << "," << edge.b.y << "), длина: " << edge.length();
        logger.trace(oss.str());
    }

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
    
    const Triangle* findContainingTriangle(const PointD& p, const std::vector<Triangle>& triangles) {
        logger.trace(std::string("[PathFinder::findContainingTriangle] Поиск треугольника для точки (") + 
                   std::to_string(p.x) + "," + std::to_string(p.y) + ")");
        
        for (const auto& tri : triangles) {
            if (isPointInTriangle(p, tri)) {
                logTriangle("[PathFinder::findContainingTriangle] Найден содержащий треугольник", &tri);
                return &tri;
            }
        }
        
        logger.warning("[PathFinder::findContainingTriangle] Точка не принадлежит ни одному треугольнику");
        return nullptr;
    }

    bool isPointInTriangle(const PointD& p, const Triangle& tri) {
        auto sign = [](PointD p1, PointD p2, PointD p3) {
            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
        };
        
        double d1 = sign(p, tri.a, tri.b);
        double d2 = sign(p, tri.b, tri.c);
        double d3 = sign(p, tri.c, tri.a);
        
        bool has_neg = (d1 < -Constants::EPSILON) || (d2 < -Constants::EPSILON) || (d3 < -Constants::EPSILON);
        bool has_pos = (d1 > Constants::EPSILON) || (d2 > Constants::EPSILON) || (d3 > Constants::EPSILON);
        
        bool result = !(has_neg && has_pos);
        logger.trace(std::string("[PathFinder::isPointInTriangle] Точка (") + 
                   std::to_string(p.x) + "," + std::to_string(p.y) + ") " + 
                   (result ? "внутри" : "вне") + " треугольника");
        return result;
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

   std::vector<const Triangle*> getNeighbors(const Triangle& tri, const std::vector<Triangle>& triangles) {
    logger.trace(std::string("[PathFinder::getNeighbors] Поиск соседей треугольника [") +
               "A(" + std::to_string(tri.a.x) + "," + std::to_string(tri.a.y) + "), " +
               "B(" + std::to_string(tri.b.x) + "," + std::to_string(tri.b.y) + "), " +
               "C(" + std::to_string(tri.c.x) + "," + std::to_string(tri.c.y) + ")]");
    
    std::vector<const Triangle*> neighbors;
    size_t totalChecked = 0;
    size_t edgesMatched = 0;

    for (const auto& other : triangles) {
        if (&tri == &other) {
            logger.trace("[PathFinder::getNeighbors] Пропуск самосравнения");
            continue;
        }
        
        totalChecked++;
        if (shareEdge(tri, other)) {
            neighbors.push_back(&other);
            edgesMatched++;
            logger.debug(std::string("[PathFinder::getNeighbors] Найдено общее ребро с треугольником [") +
                       "A(" + std::to_string(other.a.x) + "," + std::to_string(other.a.y) + "), " +
                       "B(" + std::to_string(other.b.x) + "," + std::to_string(other.b.y) + "), " +
                       "C(" + std::to_string(other.c.x) + "," + std::to_string(other.c.y) + ")]");
        }
    }

    logger.info(std::string("[PathFinder::getNeighbors] Результаты поиска: ") +
              std::to_string(neighbors.size()) + " соседей из " +
              std::to_string(totalChecked) + " проверенных треугольников, " +
              std::to_string(edgesMatched) + " общих ребер");
    
    return neighbors;
}

bool isEdgeNavigable(const Edge& edge, 
                    const std::unique_ptr<Pole>& p,
                    const std::vector<std::vector<double>>& binaryMap,
                    double sliceLevel) const {
    const auto line = bresenhamLine(edge.a, edge.b);
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
    
    bool shareEdge(const Triangle& a, const Triangle& b) const {
    logger.trace("[PathFinder::shareEdge] Проверка общего ребра между треугольниками");
    
    const std::array<Edge, 3> edgesA = { Edge(a.a, a.b), Edge(a.b, a.c), Edge(a.c, a.a) };
    const std::array<Edge, 3> edgesB = { Edge(b.a, b.b), Edge(b.b, b.c), Edge(b.c, b.a) };

    for (size_t i = 0; i < 3; ++i) {
        const Edge& edgeA = edgesA[i];
        for (size_t j = 0; j < 3; ++j) {
            const Edge& edgeB = edgesB[j];
            if (edgeA == edgeB) {
                logger.debug(std::string("[PathFinder::shareEdge] Найдено общее ребро: (") +
                           std::to_string(edgeA.a.x) + "," + std::to_string(edgeA.a.y) + ")-(" +
                           std::to_string(edgeA.b.x) + "," + std::to_string(edgeA.b.y) + ")");
                return true;
            }
        }
    }

    logger.trace("[PathFinder::shareEdge] Общих ребер не обнаружено");
    return false;
}

bool otherHasEdge(const Triangle& other, const Edge& edge) const {
    logger.trace("[PathFinder::otherHasEdge] Проверка наличия ребра в треугольнике");
    
    const Edge edges[3] = { Edge(other.a, other.b), Edge(other.b, other.c), Edge(other.c, other.a) };
    
    for (int i = 0; i < 3; ++i) {
        if (edges[i] == edge) {
            logger.debug(std::string("[PathFinder::otherHasEdge] Ребро найдено: (") +
                       std::to_string(edge.a.x) + "," + std::to_string(edge.a.y) + ")-(" +
                       std::to_string(edge.b.x) + "," + std::to_string(edge.b.y) + ")");
            return true;
        }
    }

    logger.trace("[PathFinder::otherHasEdge] Ребро отсутствует в треугольнике");
    return false;
}

    std::vector<const Triangle*> findNeighbors(const Triangle& tri, const std::vector<Triangle>& allTriangles) const {
        logger.trace("[PathFinder::findNeighbors] Поиск соседей треугольника");
        std::vector<const Triangle*> neighbors;
        
        for (const auto& other : allTriangles) {
            if (&tri == &other) continue;
            if (shareEdge(tri, other)) {
                neighbors.push_back(&other);
                logger.trace("[PathFinder::findNeighbors] Найден соседний треугольник");
            }
        }
        
        logger.debug(std::string("[PathFinder::findNeighbors] Найдено ") + std::to_string(neighbors.size()) + " соседей");
        return neighbors;
    }

    std::vector<PointD> findPathAStar(const PointD& start, const PointD& goal, 
                                 const std::vector<Triangle>& triangles,
                                 const std::vector<std::vector<double>>& binaryMap,
                                 double slice, const std::unique_ptr<Pole>& p) {
    logger.info("[PathFinder::findPathAStar] Начало поиска пути");
    
    const Triangle* startTri = findContainingTriangle(start, triangles);
    const Triangle* goalTri = findContainingTriangle(goal, triangles);
    
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

        for (const auto& neighbor : getNeighbors(*current.tri, triangles)) {
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

std::vector<PointD> bresenhamLine(const PointD& start, const PointD& end) const {
    logger.trace(std::string("[PathFinder::bresenhamLine] Построение линии от (") + 
               std::to_string(start.x) + "," + std::to_string(start.y) + ") до (" +
               std::to_string(end.x) + "," + std::to_string(end.y) + ")");
    
    std::vector<PointD> linePoints;
    int x0 = static_cast<int>(start.x);
    int y0 = static_cast<int>(start.y);
    int x1 = static_cast<int>(end.x);
    int y1 = static_cast<int>(end.y);

    logger.debug(std::string("[PathFinder::bresenhamLine] Целочисленные координаты: (") + 
               std::to_string(x0) + "," + std::to_string(y0) + ") -> (" +
               std::to_string(x1) + "," + std::to_string(y1) + ")");

    int dx = std::abs(x1 - x0);
    int dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1;
    int sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;

    logger.trace(std::string("[PathFinder::bresenhamLine] Параметры алгоритма: dx=") + std::to_string(dx) + 
               ", dy=" + std::to_string(dy) + ", err=" + std::to_string(err));

    size_t pointCount = 0;
    while (true) {
        linePoints.emplace_back(x0, y0);
        pointCount++;
        
        if (x0 == x1 && y0 == y1) {
            logger.debug(std::string("[PathFinder::bresenhamLine] Линия завершена, точек: ") + std::to_string(pointCount));
            break;
        }
        
        int e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
            logger.trace(std::string("[PathFinder::bresenhamLine] Шаг по X: x=") + std::to_string(x0));
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
            logger.trace(std::string("[PathFinder::bresenhamLine] Шаг по Y: y=") + std::to_string(y0));
        }
    }
    
    return linePoints;
}

};
