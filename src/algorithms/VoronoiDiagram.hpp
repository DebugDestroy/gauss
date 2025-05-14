#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <cmath> // для std::hypot, std::fabs

// Локальные заголовки
#include "core/Constants.hpp"  // Подключаем константы
#include "services/Geometry.hpp" // для PointD, Edge, Triangle

class VoronoiDiagram {
private:
    Logger& logger;

    void logEdge(const Edge& edge, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Ребро: (" << edge.a.x << "," << edge.a.y << ")-(" 
            << edge.b.x << "," << edge.b.y << "), длина: " << edge.length();
        logger.trace(oss.str());
    }
    
    bool isPointInsideField(const PointD& p, int width, int height) {
        bool inside = p.x >= 0 && p.y >= 0 && p.x < width && p.y < height;
        logger.trace(std::string("[VoronoiDiagram::isPointInsideField] Точка (") +
                    std::to_string(p.x) + "," + std::to_string(p.y) + ") " +
                    (inside ? "внутри" : "снаружи") + " области");
        return inside;
    }
    
    // Нахождение пересечения с границей (простая версия в реализации)
    PointD findBoundaryIntersection(const PointD& start, const PointD& dir, int width, int height) const {
    // Нормализуем направление
    double len = std::hypot(dir.x, dir.y);
    if (len < Constants::EPSILON) {
        logger.warning("[VoronoiDiagram::findBoundaryIntersection] Нулевое направление, возвращается стартовая точка");
        return start;
    }

    PointD unitDir = { dir.x / len, dir.y / len };

    double t = std::numeric_limits<double>::max();

    // Проверяем пересечения с 4 границами
    if (unitDir.x != 0) {
        if (unitDir.x > 0) t = std::min(t, (width - 1 - start.x) / unitDir.x);  // Правая граница
        else               t = std::min(t, -start.x / unitDir.x);              // Левая граница
    }

    if (unitDir.y != 0) {
        if (unitDir.y > 0) t = std::min(t, (height - 1 - start.y) / unitDir.y); // Верхняя граница
        else                t = std::min(t, -start.y / unitDir.y);              // Нижняя граница
    }

    // Проверка t ещё раз не нужна — если длина была > 0, t гарантированно установится
    PointD result = { start.x + unitDir.x * t, start.y + unitDir.y * t };

    logger.trace(std::string("[VoronoiDiagram::findBoundaryIntersection] Пересечение с границей: (") +
                 std::to_string(result.x) + "," + std::to_string(result.y) + ")");
    return result;
}
    
    // Получение рёбер выпуклой оболочки
    std::vector<Edge> getConvexHullEdges(const Triangle& tri, const std::vector<Triangle>& allTriangles) const {
        std::vector<Edge> hullEdges;
        
        for (const Edge& edge : { Edge(tri.a, tri.b), Edge(tri.b, tri.c), Edge(tri.c, tri.a) }) {
            bool isShared = false;
            for (const auto& other : allTriangles) {
                if (&other == &tri) continue;
                if (Triangle::otherHasEdge(other, edge)) {
                    isShared = true;
                    break;
                }
            }
            if (!isShared) hullEdges.push_back(edge);
        }
        
        return hullEdges;
    }

public:
    VoronoiDiagram(Logger& lg) : logger(lg) {
        logger.trace("[VoronoiDiagram] Инициализация диаграммы Вороного");
    }

    void buildFromDelaunay(const std::vector<Triangle>& triangles,
                      const std::unique_ptr<Pole>& p,
                      std::vector<Edge>& edges) {
    logger.info("[VoronoiDiagram::buildFromDelaunay] Начало построения диаграммы Вороного из триангуляции Делоне");
    if (!p) {
        logger.error("[VoronoiDiagram::buildFromDelaunay] Ошибка: данные высот не инициализированы!");
        return;
    }
    
    const int width = p->field[0].size();
    const int height = p->field.size();
    edges.clear();

    logger.debug(std::string("[VoronoiDiagram::buildFromDelaunay] Размер области: ") + 
               std::to_string(width) + "x" + std::to_string(height) + 
               ", треугольников: " + std::to_string(triangles.size()));

for (const auto& tri : triangles) {
    PointD cc = tri.calculateCircumcenter();
    bool ccValid = isPointInsideField(cc, width, height);
    
    logger.debug(std::string("[VoronoiDiagram::buildFromDelaunay] Обработка треугольника с центром (") +
               std::to_string(cc.x) + "," + std::to_string(cc.y) + "), " +
               (ccValid ? "внутри" : "снаружи") + " области");

    // Обработка внутренних рёбер
    auto neighbors = Triangle::findNeighbors(tri, triangles);
    logger.trace(std::string("[VoronoiDiagram::buildFromDelaunay] Найдено ") + 
               std::to_string(neighbors.size()) + " соседей");

   for (const auto& neighbor : neighbors) {
    auto neighborCC = neighbor->calculateCircumcenter();
    bool neighborCCValid = isPointInsideField(neighborCC, width, height);

    if (ccValid && neighborCCValid) {
        if (!isPointInsideField(cc, width, height) || !isPointInsideField(neighborCC, width, height)) {
logger.warning("!!! ВНИМАНИЕ: добавляется ребро снаружи: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(neighborCC.x) + ", " + std::to_string(neighborCC.y) + ")");
}
logger.trace("!!! ВНИМАНИЕ: добавляется ребро: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(neighborCC.x) + ", " + std::to_string(neighborCC.y) + ")");
        edges.emplace_back(cc, neighborCC);
        logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлено внутреннее ребро Вороного");
        continue;
    }

    if (ccValid && !neighborCCValid) {
        PointD direction = { neighborCC.x - cc.x, neighborCC.y - cc.y };
        PointD clipped = findBoundaryIntersection(cc, direction, width, height);
        if (!isPointInsideField(cc, width, height) || !isPointInsideField(clipped, width, height)) {
logger.warning("!!! ВНИМАНИЕ: добавляется ребро снаружи: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(clipped.x) + ", " + std::to_string(clipped.y) + ")");
}
logger.trace("!!! ВНИМАНИЕ: добавляется ребро: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(clipped.x) + ", " + std::to_string(clipped.y) + ")");
        edges.emplace_back(cc, clipped);
        logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлен обрезанный луч Вороного");
        continue;
    }
}


    // Обработка граничных рёбер
    if (ccValid && neighbors.size() < 3) {
        logger.trace("[VoronoiDiagram::buildFromDelaunay] Обработка граничного треугольника");
        
        auto hullEdges = getConvexHullEdges(tri, triangles);
        logger.trace(std::string("[VoronoiDiagram::buildFromDelaunay] Найдено ") + 
                   std::to_string(hullEdges.size()) + " рёбер выпуклой оболочки");

        for (const auto& hullEdge : hullEdges) {
            PointD mid = hullEdge.midPoint();
            PointD outwardDir = hullEdge.perpendicular();

            // Нормализуем
            double len = std::hypot(outwardDir.x, outwardDir.y);
            if (len > Constants::EPSILON) {
                outwardDir.x /= len;
                outwardDir.y /= len;
            }

            PointD boundaryPt = findBoundaryIntersection(mid, outwardDir, width, height);

            if (isPointInsideField(boundaryPt, width, height)) {
            if (!isPointInsideField(cc, width, height) || !isPointInsideField(boundaryPt, width, height)) {
logger.warning("!!! ВНИМАНИЕ: добавляется ребро снаружи: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(boundaryPt.x) + ", " + std::to_string(boundaryPt.y) + ")");
}
logger.trace("!!! ВНИМАНИЕ: добавляется ребро: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(boundaryPt.x) + ", " + std::to_string(boundaryPt.y) + ")");
                edges.emplace_back(cc, boundaryPt);
                logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлено ребро Вороного к границе области");
            } else {
                logger.warning("[VoronoiDiagram::buildFromDelaunay] Не удалось найти валидную граничную точку");
            }
        }
    }
}

logger.info(std::string("[VoronoiDiagram::buildFromDelaunay] Построение завершено, всего ребер: ") + 
          std::to_string(edges.size()));
          for (const auto& correctedge : edges) {
              logEdge(correctedge, "!!! ВНИМАНИЕ: добавлено ребро:");
              }          
          
}
};
