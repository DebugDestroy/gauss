#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <cmath> // для std::hypot, std::fabs

// Локальные заголовки
#include "core/Constants.hpp"  // Подключаем константы
#include "core/Logger.hpp"
#include "services/Geometry.hpp" // для PointD, Edge, Triangle
#include "services/Pole.hpp" // для std::unique_ptr<Pole>
#include "algorithms/PathFinder.hpp" // для PathFinder

class VoronoiDiagram {
private:
    Logger& logger;

    void logEdge(const Edge& edge, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Ребро: (" << edge.a.x << "," << edge.a.y << ")-(" 
            << edge.b.x << "," << edge.b.y << "), длина: " << edge.length();
        logger.trace(oss.str());
    }

    void logPoint(const PointD& p, const std::string& prefix = "") const {
        logger.trace(prefix + "Точка (" + std::to_string(p.x) + "," + std::to_string(p.y) + ")");
    }

public:
    VoronoiDiagram(Logger& lg) : logger(lg) {
        logger.trace("[VoronoiDiagram] Инициализация диаграммы Вороного");
    }

    void buildFromDelaunay(const std::vector<Triangle>& triangles, 
                          const PathFinder& pathFinder, 
                          const std::unique_ptr<Pole>& p,
                          std::vector<VoronoiEdge>& edges) {
        logger.info("[VoronoiDiagram::buildFromDelaunay] Построение диаграммы из триангуляции Делоне");
        
        const int height = p->field.size();
        const int width = p->field[0].size();
        edges.clear();

        if (!p) {
            logger.error("[VoronoiDiagram::buildFromDelaunay] Ошибка: данные высот не инициализированы!");
            return;
        }

        logger.debug(std::string("[VoronoiDiagram::buildFromDelaunay] Размер области: ") + 
                   std::to_string(width) + "x" + std::to_string(height) + 
                   ", треугольников: " + std::to_string(triangles.size()));

        for (const auto& tri : triangles) {
            PointD cc1 = tri.calculateCircumcenter();
            bool cc1_valid = isPointInsideField(cc1, width, height);
            
            logger.debug(std::string("[VoronoiDiagram::buildFromDelaunay] Обработка треугольника с центром (") +
                       std::to_string(cc1.x) + "," + std::to_string(cc1.y) + "), " +
                       (cc1_valid ? "внутри" : "снаружи") + " области");

            auto neighbors = pathFinder.findNeighbors(tri, triangles);
            logger.trace(std::string("[VoronoiDiagram::buildFromDelaunay] Найдено ") + 
                       std::to_string(neighbors.size()) + " соседей");

            // Обработка обычных ребер
            for (const auto& neighbor : neighbors) {
                PointD cc2 = neighbor->calculateCircumcenter();
                bool cc2_valid = isPointInsideField(cc2, width, height);

                logger.trace(std::string("[VoronoiDiagram::buildFromDelaunay] Соседний центр (") +
                           std::to_string(cc2.x) + "," + std::to_string(cc2.y) + "), " +
                           (cc2_valid ? "внутри" : "снаружи") + " области");

                if (cc1_valid && cc2_valid) {
                    edges.emplace_back(cc1, cc2);
                    logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлено полное ребро Вороного");
                } 
                else if (cc1_valid || cc2_valid) {
                    auto clipped = clipEdgeToField(cc1, cc2, width, height);
                    if (!clipped.empty()) {
                        edges.emplace_back(clipped[0], clipped[1]);
                        logger.debug("[VoronoiDiagram::buildFromDelaunay] Добавлено обрезанное ребро Вороного");
                    } else {
                        logger.trace("[VoronoiDiagram::buildFromDelaunay] Ребро полностью вне области после отсечения");
                    }
                }
            }

            // Обработка граничных треугольников
            if (neighbors.size() < 3) {
                logger.trace("[VoronoiDiagram::buildFromDelaunay] Обработка граничного треугольника");
                handleBoundaryTriangle(tri, cc1, width, height, edges, cc1_valid, triangles, pathFinder);
            }
        }

        logger.info(std::string("[VoronoiDiagram::buildFromDelaunay] Построение завершено, ребер: ") + 
                  std::to_string(edges.size()));
    }

private:
    std::vector<PointD> clipEdgeToField(const PointD& p1, const PointD& p2, int width, int height) {
        logger.trace(std::string("[VoronoiDiagram::clipEdgeToField] Отсечение ребра (") +
                   std::to_string(p1.x) + "," + std::to_string(p1.y) + ")-(" +
                   std::to_string(p2.x) + "," + std::to_string(p2.y) + ")");

        auto code = [&](const PointD& p) {
            int c = 0;
            if (p.x < 0) c |= 1;
            if (p.x > width) c |= 2;
            if (p.y < 0) c |= 4;
            if (p.y > height) c |= 8;
            return c;
        };

        int code1 = code(p1);
        int code2 = code(p2);
        PointD a = p1, b = p2;

        logger.trace(std::string("[VoronoiDiagram::clipEdgeToField] Коды: p1=") + std::to_string(code1) + 
                   ", p2=" + std::to_string(code2));

        while (true) {
            if (!(code1 | code2)) {
                logger.trace("[VoronoiDiagram::clipEdgeToField] Ребро полностью внутри области");
                return {a, b};
            }
            if (code1 & code2) {
                logger.trace("[VoronoiDiagram::clipEdgeToField] Ребро полностью снаружи области");
                return {};
            }

            int outcode = code1 ? code1 : code2;
            PointD p;

            if (outcode & 8) {
                p.x = a.x + (b.x - a.x) * (height - a.y) / (b.y - a.y);
                p.y = height;
                logger.trace("[VoronoiDiagram::clipEdgeToField] Пересечение с верхней границей");
            }
            else if (outcode & 4) {
                p.x = a.x + (b.x - a.x) * (-a.y) / (b.y - a.y);
                p.y = 0;
                logger.trace("[VoronoiDiagram::clipEdgeToField] Пересечение с нижней границей");
            }
            else if (outcode & 2) {
                p.y = a.y + (b.y - a.y) * (width - a.x) / (b.x - a.x);
                p.x = width;
                logger.trace("[VoronoiDiagram::clipEdgeToField] Пересечение с правой границей");
            }
            else if (outcode & 1) {
                p.y = a.y + (b.y - a.y) * (-a.x) / (b.x - a.x);
                p.x = 0;
                logger.trace("[VoronoiDiagram::clipEdgeToField] Пересечение с левой границей");
            }

            if (outcode == code1) {
                a = p;
                code1 = code(a);
            } else {
                b = p;
                code2 = code(b);
            }
        }
    }

    bool isPointInsideField(const PointD& p, int width, int height) {
        bool inside = p.x >= 0 && p.y >= 0 && p.x < width && p.y < height;
        logger.trace(std::string("[VoronoiDiagram::isPointInsideField] Точка (") +
                    std::to_string(p.x) + "," + std::to_string(p.y) + ") " +
                    (inside ? "внутри" : "снаружи") + " области");
        return inside;
    }

    void handleBoundaryTriangle(const Triangle& tri, const PointD& cc, int width, int height, 
                               std::vector<VoronoiEdge>& edges, bool cc_valid,
                               const std::vector<Triangle>& allTriangles,  
                               const PathFinder& pathFinder) {
        if (!cc_valid) {
            logger.trace("[VoronoiDiagram::handleBoundaryTriangle] Центр снаружи, пропуск");
            return;
        }

        auto boundaryEdges = getBoundaryEdges(tri, allTriangles, pathFinder);
        logger.debug(std::string("[VoronoiDiagram::handleBoundaryTriangle] Найдено ") +
                   std::to_string(boundaryEdges.size()) + " граничных ребер");

        for (const auto& edge : boundaryEdges) {
            PointD boundaryPoint = calculateBoundaryIntersection(cc, edge, width, height);
            edges.emplace_back(cc, boundaryPoint);
            logger.debug("[VoronoiDiagram::handleBoundaryTriangle] Добавлено граничное ребро Вороного");
        }
    }

    std::vector<Edge> getBoundaryEdges(const Triangle& tri, 
                                      const std::vector<Triangle>& allTriangles,  
                                      const PathFinder& pathFinder) {
        std::vector<Edge> boundaryEdges;
        logger.trace("[VoronoiDiagram::getBoundaryEdges] Поиск граничных ребер треугольника");

        for (const auto& edge : { Edge(tri.a, tri.b), Edge(tri.b, tri.c), Edge(tri.c, tri.a) }) {
            bool isBoundary = true;
            logEdge(edge, "Проверка ребра: ");

            for (const auto& other : allTriangles) {
                if (&tri == &other) continue;
                if (pathFinder.shareEdge(tri, other) && pathFinder.otherHasEdge(other, edge)) {
                    isBoundary = false;
                    break;
                }
            }

            if (isBoundary) {
                boundaryEdges.push_back(edge);
                logger.trace("[VoronoiDiagram::getBoundaryEdges] Найдено граничное ребро");
            }
        }

        return boundaryEdges;
    }

    PointD calculateBoundaryIntersection(const PointD& circumCenter, 
                                        const Edge& boundaryEdge, 
                                        int width, int height) {
        logger.trace("[VoronoiDiagram::calculateBoundaryIntersection] Вычисление пересечения с границей");

        PointD edgeDir = { boundaryEdge.b.x - boundaryEdge.a.x, boundaryEdge.b.y - boundaryEdge.a.y };
        PointD normal = { -edgeDir.y, edgeDir.x };
        double length = std::hypot(normal.x, normal.y);
        
        if (std::fabs(length) < Constants::EPSILON) {
            logger.warning("[VoronoiDiagram::calculateBoundaryIntersection] Нулевая длина нормали!");
            return circumCenter;
        }
        
        normal.x /= length;
        normal.y /= length;

        double t = std::numeric_limits<double>::max();
        if (std::fabs(normal.x) > Constants::EPSILON) {
            t = std::min(t, (width - circumCenter.x) / normal.x);
            t = std::min(t, -circumCenter.x / normal.x);
        }
        if (std::fabs(normal.y) > Constants::EPSILON) {
            t = std::min(t, (height - circumCenter.y) / normal.y);
            t = std::min(t, -circumCenter.y / normal.y);
        }
        
        PointD result = {
            circumCenter.x + normal.x * t,
            circumCenter.y + normal.y * t
        };

        logger.debug(std::string("[VoronoiDiagram::calculateBoundaryIntersection] Точка пересечения: (") +
                    std::to_string(result.x) + "," + std::to_string(result.y) + ")");
        
        return result;
    }
};
