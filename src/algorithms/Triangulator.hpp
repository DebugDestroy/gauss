#pragma once
#include <vector>
#include <algorithm> // для std::remove_if, std::find, std::all_of
#include <cmath> // для std::isnan, std::sqrt
#include <string>
#include <sstream>
#include <limits> // для std::numeric_limits

// Локальные заголовки
#include "core/Logger.hpp"
#include "services/Geometry.hpp" // для PointD, Triangle, Edge

class Triangulator {
private:
    Logger& logger;

    void logTriangle(const Triangle& tri, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Треугольник ["
            << "A(" << tri.a.x << "," << tri.a.y << "), "
            << "B(" << tri.b.x << "," << tri.b.y << "), "
            << "C(" << tri.c.x << "," << tri.c.y << ")]";
        logger.debug(oss.str());
    }

    void logEdge(const Edge& edge, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Ребро: (" << edge.a.x << "," << edge.a.y << ")-(" 
            << edge.b.x << "," << edge.b.y << ")";
        logger.trace(oss.str());
    }
    
    bool isPointInCircumcircle(const PointD& p, const Triangle& tri) {
        double ax = tri.a.x - p.x, ay = tri.a.y - p.y;
        double bx = tri.b.x - p.x, by = tri.b.y - p.y;
        double cx = tri.c.x - p.x, cy = tri.c.y - p.y;
        
        double det = ax * (by * (cx*cx + cy*cy) - cy * (bx*bx + by*by)) 
                   - ay * (bx * (cx*cx + cy*cy) - cx * (bx*bx + by*by)) 
                   + (ax*ax + ay*ay) * (bx*cy - by*cx);

        bool inside = det > 0;
        logger.trace(std::string("[Triangulator::isPointInCircumcircle] Точка (") + 
                    std::to_string(p.x) + "," + std::to_string(p.y) + ") " +
                    (inside ? "внутри" : "вне") + " окружности треугольника");
        logTriangle(tri, "Проверяемый треугольник: ");

        return inside;
    }

    bool hasEdge(const Triangle& tri, const Edge& edge) {
        bool has = Edge(tri.a, tri.b) == edge 
                || Edge(tri.b, tri.c) == edge 
                || Edge(tri.c, tri.a) == edge;

        logger.trace(std::string("[Triangulator::hasEdge] Треугольник ") + 
                    std::string(has ? "содержит" : "не содержит") + " ребро");
        logTriangle(tri);
        logEdge(edge);

        return has;
    }
public:
    Triangulator(Logger& lg) : logger(lg) {
        logger.trace("[Triangulator] Инициализация калькулятора компонент");
    }

    std::vector<Triangle> bowyerWatson(const std::vector<PointD>& points) {
        logger.info("[Triangulator::bowyerWatson] Начало триангуляции Боуера-Ватсона");
        logger.debug(std::string("Количество точек: ") + std::to_string(points.size()));
        
        std::vector<Triangle> triangles;

        // Проверка на минимальное количество точек
        if (points.size() < 3) {
            logger.error("[Triangulator::bowyerWatson] Недостаточно точек для триангуляции (меньше 3)");
            return triangles;
        }

        // Проверка на NaN-точки
        for (const auto& p : points) {
            if (std::isnan(p.x) || std::isnan(p.y)) {
                logger.error("[Triangulator::bowyerWatson] Обнаружена некорректная точка (NaN)");
                return triangles;
            }
        }

        // Проверка на коллинеарность
        if (std::all_of(points.begin(), points.end(), 
                       [&](const PointD& p) { return Triangle::areCollinear(points[0], points[1], p); })) {
            logger.error("[Triangulator::bowyerWatson] Все точки коллинеарны, триангуляция невозможна");
            return triangles;
        }

        // Создаём супер-треугольник
        auto [minX, maxX] = std::minmax_element(points.begin(), points.end(), 
                                              [](auto& a, auto& b) { return a.x < b.x; });
        auto [minY, maxY] = std::minmax_element(points.begin(), points.end(), 
                                              [](auto& a, auto& b) { return a.y < b.y; });

        double dx = (maxX->x - minX->x) * 10;
        double dy = (maxY->y - minY->y) * 10;
        PointD p1(minX->x - dx, minY->y - dy);
        PointD p2(maxX->x + dx, minY->y - dy);
        PointD p3((minX->x + maxX->x) / 2, maxY->y + dy);

        triangles.emplace_back(p1, p2, p3);
        logger.info("[Triangulator::bowyerWatson] Создан супер-треугольник");
        logTriangle(triangles.back(), "Супер-треугольник: ");

        // Основной алгоритм
        size_t pointIndex = 0;
        for (const auto& point : points) {
            logger.trace(std::string("[Triangulator::bowyerWatson] Обработка точки #") + 
                        std::to_string(++pointIndex) + " (" + 
                        std::to_string(point.x) + "," + std::to_string(point.y) + ")");

            std::vector<Triangle> badTriangles;
            std::copy_if(triangles.begin(), triangles.end(), std::back_inserter(badTriangles),
                         [&](const Triangle& tri) { return isPointInCircumcircle(point, tri); });

            logger.debug(std::string("[Triangulator::bowyerWatson] Найдено ") + 
                       std::to_string(badTriangles.size()) + " плохих треугольников");

            std::vector<Edge> polygonEdges;
            for (const auto& tri : badTriangles) {
                for (const auto& edge : { Edge{tri.a, tri.b}, Edge{tri.b, tri.c}, Edge{tri.c, tri.a} }) {
                    bool isShared = std::any_of(badTriangles.begin(), badTriangles.end(),
                        [&](const Triangle& other) { 
                            return !(tri == other) && hasEdge(other, edge); 
                        });
                    
                    if (!isShared) {
                        polygonEdges.push_back(edge);
                        logEdge(edge, "Добавлено ребро в полигон: ");
                    }
                }
            }

            // Удаляем плохие треугольники
            triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
                [&](const Triangle& t) { 
                    return std::find(badTriangles.begin(), badTriangles.end(), t) != badTriangles.end();
                }), triangles.end());

            // Добавляем новые треугольники
            for (const auto& edge : polygonEdges) {
                triangles.emplace_back(edge.a, edge.b, point);
                logTriangle(triangles.back(), "Добавлен новый треугольник: ");
            }
        }

        // Удаляем треугольники, связанные с супер-треугольником
        size_t initialCount = triangles.size();
        triangles.erase(std::remove_if(triangles.begin(), triangles.end(),
            [&](const Triangle& t) {
                return t.a == p1 || t.a == p2 || t.a == p3 ||
                       t.b == p1 || t.b == p2 || t.b == p3 ||
                       t.c == p1 || t.c == p2 || t.c == p3;
            }), triangles.end());

        logger.info("[Triangulator::bowyerWatson] Триангуляция завершена");
        logger.debug(std::string("Удалено ") + std::to_string(initialCount - triangles.size()) + 
                    " треугольников, связанных с супер-треугольником");
        logger.debug(std::string("Итоговое количество треугольников: ") + std::to_string(triangles.size()));

        return triangles;
    }
};
