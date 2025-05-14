#pragma once
#include <vector>
#include <algorithm> // для std::remove_if, std::find, std::all_of
#include <cmath> // для std::isnan, std::sqrt
#include <string>
#include <sstream>
#include <limits> // для std::numeric_limits

// Локальные заголовки
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

public:
    Triangulator(Logger& lg) : logger(lg) {
        logger.trace("[Triangulator] Инициализация калькулятора компонент");
    }

    std::vector<Triangle> bowyerWatson(const std::vector<PointD>& points, const int width, const int height) {
         std::vector<Triangle> triangles;// Треугольники

        // Вычисляем размер супер-треугольника (2 * максимальная сторона)
        const double super_size = 2 * std::max(width, height);

         // Создаём вершины треугольника (в целочисленных координатах для точности)
         PointD p1(0, 0);                       // Левый нижний угол
         PointD p2(0, super_size);               // Левый верхний угол
         PointD p3(super_size, 0);               // Правый нижний угол
         
        logger.info("[Triangulator::bowyerWatson] Начало триангуляции Боуера-Ватсона");
        logger.debug(std::string("Количество точек: ") + std::to_string(points.size()));

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

// Создаём и логируем супер-треугольник
triangles.emplace_back(p1, p2, p3);
logger.info("[Triangulator] Создан супер-треугольник с вершинами:");
logger.info("  A(0, 0)");
logger.info("  B(0, " + std::to_string(super_size) + ")");
logger.info("  C(" + std::to_string(super_size) + ", 0)");

        // Основной алгоритм
        size_t pointIndex = 0;
        for (const auto& point : points) {
            
            logger.trace(std::string("[Triangulator::bowyerWatson] Обработка точки #") + 
                        std::to_string(++pointIndex) + " (" + 
                        std::to_string(point.x) + "," + std::to_string(point.y) + ")");

            std::vector<Triangle> badTriangles;
            std::copy_if(triangles.begin(), triangles.end(), std::back_inserter(badTriangles),
                         [&](const Triangle& tri) { return Triangle::isPointInCircumcircle(point, tri); });

            logger.debug(std::string("[Triangulator::bowyerWatson] Найдено ") + 
                       std::to_string(badTriangles.size()) + " плохих треугольников");

            std::vector<Edge> polygonEdges;
            for (const auto& tri : badTriangles) {
                for (const auto& edge : { Edge{tri.a, tri.b}, Edge{tri.b, tri.c}, Edge{tri.c, tri.a} }) {
                    bool isShared = std::any_of(badTriangles.begin(), badTriangles.end(),
                        [&](const Triangle& other) { 
                            return !(tri == other) && Triangle::otherHasEdge(other, edge); 
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
