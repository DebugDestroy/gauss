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
        logger.debug(oss.str());
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
                      const std::vector<std::vector<double>>& field,
                      std::vector<Edge>& edges) {
    logger.info("[VoronoiDiagram::buildFromDelaunay] Начало построения диаграммы Вороного из триангуляции Делоне");
    if (field.empty()) {
        logger.error("[VoronoiDiagram::buildFromDelaunay] Ошибка: данные высот не инициализированы!");
        return;
    }
    
    if (triangles.empty()) {
        logger.error("[VoronoiDiagram::buildFromDelaunay] No triangulation for VoronoiDiagram");
        return;
    }
    
    const int width = static_cast<int>(field[0].size());
    const int height = static_cast<int>(field.size());
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
Edge newEdge(cc, neighborCC);

if (!containsEdge(edges, newEdge)) {
    logger.trace("Добавляется внутреннее ребро: (" +
                 std::to_string(cc.x) + ", " +
                 std::to_string(cc.y) + ") (" +
                 std::to_string(neighborCC.x) + ", " +
                 std::to_string(neighborCC.y) + ")");

    edges.push_back(newEdge);

    logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлено внутреннее ребро Вороного");
}
else {
    logger.trace("[VoronoiDiagram::buildFromDelaunay] Пропущено дублирующееся ребро");
}
        continue;
    }

    if (ccValid && !neighborCCValid) {
        PointD direction = { neighborCC.x - cc.x, neighborCC.y - cc.y };
        PointD clipped = findBoundaryIntersection(cc, direction, width, height);
        if (!isPointInsideField(cc, width, height) || !isPointInsideField(clipped, width, height)) {
logger.warning("!!! ВНИМАНИЕ: добавляется ребро снаружи: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(clipped.x) + ", " + std::to_string(clipped.y) + ")");
}
logger.trace("!!! ВНИМАНИЕ: добавляется ребро: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(clipped.x) + ", " + std::to_string(clipped.y) + ")");
Edge newEdge(cc, clipped);

if (!containsEdge(edges, newEdge)) {
    edges.push_back(newEdge);

    logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлен обрезанный луч Вороного");
}
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
            // 1. считаем инцентр треугольника
            PointD incenter = calculateIncenter(tri);

            // 2. направление "наружу" (локально)
            PointD referenceDir = {
                mid.x - incenter.x,
                mid.y - incenter.y
            };

            // 3. нормаль к ребру (две возможные)
            PointD normal = perpendicular(hullEdge);

            // 4. выбираем правильную ориентацию нормали
            PointD outwardNormal = makeCodirected(normal, referenceDir);
            
            // Нормализуем
            double len = std::hypot(outwardNormal.x, outwardNormal.y);
            if (len > core::EPSILON) {
                outwardNormal.x /= len;
                outwardNormal.y /= len;
            }

            PointD boundaryPt = findBoundaryIntersection(mid, outwardNormal, width, height);

            if (isPointInsideField(boundaryPt, width, height)) {
            if (!isPointInsideField(cc, width, height) || !isPointInsideField(boundaryPt, width, height)) {
logger.warning("!!! ВНИМАНИЕ: добавляется ребро снаружи: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(boundaryPt.x) + ", " + std::to_string(boundaryPt.y) + ")");
}
logger.trace("!!! ВНИМАНИЕ: добавляется ребро: (" + std::to_string(cc.x) + ", " + std::to_string(cc.y) + ") (" + std::to_string(boundaryPt.x) + ", " + std::to_string(boundaryPt.y) + ")");
Edge newEdge(cc, boundaryPt);

if (!containsEdge(edges, newEdge)) {
    edges.push_back(newEdge);

    logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлено ребро Вороного к границе области");
}
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
          logger.info("[VoronoiDiagram] Итоговый граф Вороного:");


}
}
