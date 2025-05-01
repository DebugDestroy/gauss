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
#include "core/Config.hpp"
#include "core/Logger.hpp"
#include "core/Geometry.hpp" // для PointD, Triangle, Edge
#include "algorithms/TerrainGrid.hpp" // для TerrainGrid
#include "core/Pole.hpp" // для std::unique_ptr<Pole>

struct AStarNode {
    PointD position;       // Позиция узла (центр треугольника)
    const Triangle* tri;   // Связанный треугольник
    AStarNode* parent;     // Родительский узел
    double g;              // Стоимость пути от старта
    double h;              // Эвристическая оценка до цели
    double f;              // Общая стоимость: f = g + h

    AStarNode(PointD pos, const Triangle* triangle, AStarNode* p = nullptr)
        : position(pos), tri(triangle), parent(p), g(0), h(0), f(0) {}

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

public:
    PathFinder(const Config& cfg, Logger& lg) : config(cfg), logger(lg) {
        logger.trace("[PathFinder] Инициализация поисковика пути");
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
        
        bool has_neg = (d1 < -1e-6) || (d2 < -1e-6) || (d3 < -1e-6);
        bool has_pos = (d1 > 1e-6) || (d2 > 1e-6) || (d3 > 1e-6);
        
        bool result = !(has_neg && has_pos);
        logger.trace(std::string("[PathFinder::isPointInTriangle] Точка (") + 
                   std::to_string(p.x) + "," + std::to_string(p.y) + ") " + 
                   (result ? "внутри" : "вне") + " треугольника");
        return result;
    }
    
    bool isNavigable(const Edge& edge, const TerrainGrid& terrainGrid) const {
        logEdge(edge);
        auto [forwardAngle, sideAngle] = terrainGrid.getEdgeSlopes(edge);
        
        bool navigable = abs(forwardAngle) <= config.maxUpDownAngle && 
               abs(sideAngle) <= config.maxSideAngle &&
               edge.length() > config.vehicleRadius * 2;
        
        if (!navigable) {
            std::ostringstream oss;
            oss << "[PathFinder::isNavigable] Ребро НЕ проходимо:\n"
                << "  Уклон вперед: " << forwardAngle << "° (макс. " << config.maxUpDownAngle << "°)\n"
                << "  Боковой уклон: " << sideAngle << "° (макс. " << config.maxSideAngle << "°)\n"
                << "  Длина: " << edge.length() << " (мин. " << config.vehicleRadius * 2 << ")";
            logger.debug(oss.str());
        }
        
        return navigable;
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
                                      const TerrainGrid& terrainGrid,
                                      const std::vector<std::vector<double>>& binaryMap,
                                      double slice, const std::unique_ptr<Pole>& p) {
        logger.info("[PathFinder::findPathAStar] Начало поиска пути");
        logPoint("Стартовая точка:", start);
        logPoint("Целевая точка:", goal);
        
        std::vector<PointD> path;
        const Triangle* startTri = findContainingTriangle(start, triangles);
        const Triangle* goalTri = findContainingTriangle(goal, triangles);
        
        if (!startTri || !goalTri) {
            logger.error("[PathFinder::findPathAStar] Старт или цель вне триангуляции");
            logTriangle("Стартовый треугольник:", startTri);
            logTriangle("Целевой треугольник:", goalTri);
            return {};
        }

        if (startTri == goalTri) {
            logger.info("[PathFinder::findPathAStar] Старт и цель в одном треугольнике");
            path = {start, goal};
            if (!checkPathFeasibility(path, terrainGrid, binaryMap, slice, p)) {
                logger.warning("[PathFinder::findPathAStar] Прямой путь непроходим");
                return {};
            }
            logger.info("[PathFinder::findPathAStar] Прямой путь проходим");
            return path;
        }

        logger.debug("[PathFinder::findPathAStar] Запуск алгоритма A*");
        std::priority_queue<AStarNode> openSet;
        std::unordered_map<const Triangle*, double> costSoFar;

        PointD startCenter = startTri->calculateCircumcenter();
        openSet.emplace(startCenter, startTri);
        costSoFar[startTri] = 0;
        
        logger.debug("[PathFinder::findPathAStar] Начальный узел добавлен в очередь");

        while (!openSet.empty()) {
            AStarNode current = openSet.top();
            openSet.pop();

            if (current.tri == goalTri) {
                logger.info("[PathFinder::findPathAStar] Достигнут целевой треугольник");
                while (current.parent) {
                    path.push_back(current.position);
                    current = *current.parent;
                }
                std::reverse(path.begin(), path.end());
                
                if (!path.empty()) {
                    path.insert(path.begin(), start);
                    path.push_back(goal);
                    
                    logger.debug(std::string("[PathFinder::findPathAStar] Найден путь из ") + 
                               std::to_string(path.size()) + " точек");
                    
                    if (!checkPathFeasibility(path, terrainGrid, binaryMap, slice, p)) {
                        logger.warning("[PathFinder::findPathAStar] Путь непроходим");
                        return {};
                    }
                    logger.info("[PathFinder::findPathAStar] Путь проходим");
                }
                return path;
            }

            for (const auto& neighbor : getNeighbors(*current.tri, triangles)) {
                Edge edgeToCheck(current.position, neighbor->calculateCircumcenter());
                
                if (!isNavigable(edgeToCheck, terrainGrid)) {
                    logger.trace("[PathFinder::findPathAStar] Ребро непроходимо, пропускаем");
                    continue;
                }

                PointD neighborCenter = neighbor->calculateCircumcenter();
                double newCost = current.g + heuristic(current.position, neighborCenter);

                if (!costSoFar.count(neighbor) || newCost < costSoFar[neighbor]) {
                    costSoFar[neighbor] = newCost;
                    //double priority = newCost + heuristic(neighborCenter, goal); !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
                    openSet.emplace(neighborCenter, neighbor, new AStarNode(current));
                    
                    logger.trace(std::string("[PathFinder::findPathAStar] Добавлен новый узел в очередь: (") +
                               std::to_string(neighborCenter.x) + "," + 
                               std::to_string(neighborCenter.y) + ")");
                }
            }
        }

        logger.warning("[PathFinder::findPathAStar] Путь не найден");
        return {};
    }

    bool checkPathFeasibility(const std::vector<PointD>& path,
                              const TerrainGrid& terrainGrid,
                              const std::vector<std::vector<double>>& binaryMap,
                              double slice, const std::unique_ptr<Pole>& p) const {
        logger.info("[PathFinder::checkPathFeasibility] Проверка проходимости пути");
        logger.debug(std::string("Сегментов пути: ") + std::to_string(path.size()-1));

        for (size_t i = 0; i < path.size() - 1; ++i) {
            std::ostringstream segInfo;
            segInfo << "Сегмент " << i+1 << ": (" 
                   << path[i].x << "," << path[i].y << ") -> (" 
                   << path[i+1].x << "," << path[i+1].y << ")";
            logger.trace(segInfo.str());
            
            const auto linePixels = bresenhamLine(path[i], path[i+1]);
            
            for (const auto& pixel : linePixels) {
                auto [forwardAngle, sideAngle] = getVehicleSlopeAngles(pixel, path[i+1], terrainGrid);
                
                if (forwardAngle > config.maxUpDownAngle || sideAngle > config.maxSideAngle) {
                    std::ostringstream oss;
                    oss << "Превышены углы наклона в точке (" << pixel.x << "," << pixel.y << "):\n"
                        << "  Уклон вперед: " << forwardAngle << "° (макс. " << config.maxUpDownAngle << "°)\n"
                        << "  Боковой уклон: " << sideAngle << "° (макс. " << config.maxSideAngle << "°)";
                    logger.debug(oss.str());
                    return false;
                }

                if (!isVehicleRadiusValid(pixel, binaryMap, p, slice)) {
                    logger.debug(std::string("Столкновение в точке (") + 
                               std::to_string(pixel.x) + "," + 
                               std::to_string(pixel.y) + ")");
                    return false;
                }
            }
        }
        
        logger.info("[PathFinder::checkPathFeasibility] Путь проходим");
        return true;
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

    int dx = abs(x1 - x0);
    int dy = -abs(y1 - y0);
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

std::pair<double, double> getVehicleSlopeAngles(
    const PointD& pixel, 
    const PointD& nextPixel, 
    const TerrainGrid& grid
) const {
    logger.trace(std::string("[PathFinder::getVehicleSlopeAngles] Расчет углов наклона для точки (") +
               std::to_string(pixel.x) + "," + std::to_string(pixel.y) + ")");

    if (pixel == nextPixel) {
        logger.debug("[PathFinder::getVehicleSlopeAngles] Точки совпадают, углы = 0");
        return {0, 0};
    }

    int x = static_cast<int>(pixel.x);
    int y = static_cast<int>(pixel.y);
    
    if (x <= 0 || y <= 0 || x >= static_cast<int>(grid.cells[0].size()-1) || y >= static_cast<int>(grid.cells.size()-1)) {
        logger.warning(std::string("[PathFinder::getVehicleSlopeAngles] Точка на границе сетки (") +
                     std::to_string(x) + "," + std::to_string(y) + "), помечаем как непроходимую");
        return {config.maxUpDownAngle + 1, config.maxSideAngle + 1};
    }

    PointD dir = {nextPixel.x - pixel.x, nextPixel.y - pixel.y};
    double length = std::hypot(dir.x, dir.y);
    
    if (length < 1e-6) {
        logger.debug("[PathFinder::getVehicleSlopeAngles] Нулевой вектор направления");
        return {0, 0};
    }
    
    dir.x /= length;
    dir.y /= length;
    logger.trace(std::string("[PathFinder::getVehicleSlopeAngles] Нормализованный вектор направления: (") +
               std::to_string(dir.x) + "," + std::to_string(dir.y) + ")");

    double z = grid.cells[y][x].height;
    double z_right = grid.cells[y][x+1].height;
    double z_left = grid.cells[y][x-1].height;
    double z_top = grid.cells[y-1][x].height;
    double z_bottom = grid.cells[y+1][x].height;

    logger.trace(std::string("[PathFinder::getVehicleSlopeAngles] Высоты вокруг точки: центр=") + std::to_string(z) +
               ", слева=" + std::to_string(z_left) + ", справа=" + std::to_string(z_right) +
               ", сверху=" + std::to_string(z_top) + ", снизу=" + std::to_string(z_bottom));

    double dzdx = (z_right - z_left) / 2.0;
    double dzdy = (z_bottom - z_top) / 2.0;
    logger.trace(std::string("[PathFinder::getVehicleSlopeAngles] Градиенты: dzdx=") + std::to_string(dzdx) +
               ", dzdy=" + std::to_string(dzdy));

    double forwardAngle = std::atan(dzdx * dir.x + dzdy * dir.y) * 180.0 / M_PI;
    PointD perp = {-dir.y, dir.x};
    double sideAngle = std::atan(dzdx * perp.x + dzdy * perp.y) * 180.0 / M_PI;

    logger.debug(std::string("[PathFinder::getVehicleSlopeAngles] Результат: forward=") + std::to_string(forwardAngle) +
               "°, side=" + std::to_string(sideAngle) + "°");
    
    return {std::abs(forwardAngle), std::abs(sideAngle)};
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

    const double heightDiff = std::abs(currentHeight - sliceLevel);
    if (heightDiff >= config.vehicleRadius) {
        logger.debug(std::string("[PathFinder::isVehicleRadiusValid] Высота вне радиуса (") + 
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
                binaryMap[ny][nx] <= 255 && binaryMap[ny][nx] >= 255) {
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
};
