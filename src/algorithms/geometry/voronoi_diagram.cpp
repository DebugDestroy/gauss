#include "algorithms/geometry/voronoi_diagram.hpp"
#include "algorithms/geometry/math.hpp"

#include <sstream>
#include <cmath>
#include <limits>

namespace algorithms::geometry {

    void VoronoiDiagram::logEdge(const Edge& edge, const std::string& prefix) const {
        std::ostringstream oss;
        oss << prefix << "Ребро: (" << edge.a.x << "," << edge.a.y << ")-(" 
            << edge.b.x << "," << edge.b.y << "), длина: " << length(edge);
        logger.trace(oss.str());
    }
    
    // Нахождение пересечения с границей (простая версия в реализации)
    PointD VoronoiDiagram::findBoundaryIntersection(const PointD& start, const PointD& dir, int width, int height) const {
    // Нормализуем направление
    double len = std::hypot(dir.x, dir.y);
    if (len < core::EPSILON) {
        logger.warning("[VoronoiDiagram::findBoundaryIntersection] Нулевое направление, возвращается стартовая точка");
        return start;
    }

    PointD unitDir = { dir.x / len, dir.y / len };

    double t = std::numeric_limits<double>::max();

    // Проверяем пересечения с 4 границами
    if (std::fabs(unitDir.x) > core::EPSILON) {
        if (unitDir.x > 0) t = std::min(t, (width - 1 - start.x) / unitDir.x);  // Правая граница
        else               t = std::min(t, -start.x / unitDir.x);              // Левая граница
    }

    if (std::fabs(unitDir.y) > core::EPSILON) {
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
    std::vector<Edge> VoronoiDiagram::getConvexHullEdges(const Triangle& tri, const std::vector<Triangle>& allTriangles) const {
        std::vector<Edge> hullEdges;
        
        for (const Edge& edge : { Edge(tri.a, tri.b), Edge(tri.b, tri.c), Edge(tri.c, tri.a) }) {
            bool isShared = false;
            for (const auto& other : allTriangles) {
                if (&other == &tri) continue;
                if (otherHasEdge(other, edge)) {
                    isShared = true;
                    break;
                }
            }
            if (!isShared) hullEdges.push_back(edge);
        }
        
        return hullEdges;
    }

    VoronoiDiagram::VoronoiDiagram(core::Logger& lg) : logger(lg) {
        logger.trace("[VoronoiDiagram] Инициализация диаграммы Вороного");
    }

    void VoronoiDiagram::buildFromDelaunay(const std::vector<Triangle>& triangles,
                      const std::unique_ptr<algorithms::gauss::Pole>& p,
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
    PointD cc = calculateCircumcenter(tri.a, tri.b, tri.c);
    bool ccValid = isPointInsideField(cc, width, height);
    
    logger.debug(std::string("[VoronoiDiagram::buildFromDelaunay] Обработка треугольника с центром (") +
               std::to_string(cc.x) + "," + std::to_string(cc.y) + "), " +
               (ccValid ? "внутри" : "снаружи") + " области");

    // Обработка внутренних рёбер
    auto neighbors = findNeighbors(tri, triangles);
    logger.trace(std::string("[VoronoiDiagram::buildFromDelaunay] Найдено ") + 
               std::to_string(neighbors.size()) + " соседей");

   for (const auto& neighbor : neighbors) {
    auto neighborCC = calculateCircumcenter(neighbor->a, neighbor->b, neighbor->c);
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
            PointD mid = midPoint(hullEdge);
            PointD outwardDir = perpendicular(hullEdge);

            // Нормализуем
            double len = std::hypot(outwardDir.x, outwardDir.y);
            if (len > core::EPSILON) {
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
}
