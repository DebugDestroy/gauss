#pragma once
#include <cmath>    // Для std::fabs, std::hypot
#include <vector>   // Для использования в других классах
#include <string>   // Для возможного использования строк

// Локальные заголовки
#include "core/Constants.hpp"  // Подключаем константы

// Структура для точки с вещественными координатами (для точности) auto clusterCenters = getClusterCenters();
struct PointD {
    double x, y;
    PointD(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    
    bool operator==(const PointD& other) const {
        return std::fabs(x - other.x) < Constants::EPSILON && std::fabs(y - other.y) < Constants::EPSILON;
    }
    
    bool operator!=(const PointD& other) const {
        return !(*this == other);
    }
};

struct Edge {
    PointD a, b;
    Edge(PointD a_, PointD b_) : a(a_), b(b_) {}
    
    bool operator==(const Edge& other) const {
        return (a == other.a && b == other.b) || (a == other.b && b == other.a);
    }
    
    double length() const {
        return std::hypot(a.x - b.x, a.y - b.y);
    }
    
   static std::vector<PointD> bresenhamLine(const PointD& start, const PointD& end) {
    
    std::vector<PointD> linePoints;
    int x0 = static_cast<int>(start.x);
    int y0 = static_cast<int>(start.y);
    int x1 = static_cast<int>(end.x);
    int y1 = static_cast<int>(end.y);

    int dx = std::abs(x1 - x0);
    int dy = -std::abs(y1 - y0);
    int sx = x0 < x1 ? 1 : -1;
    int sy = y0 < y1 ? 1 : -1;
    int err = dx + dy;

    size_t pointCount = 0;
    while (true) {
        linePoints.emplace_back(x0, y0);
        pointCount++;
        
        if (x0 == x1 && y0 == y1) {
            break;
        }
        
        int e2 = 2 * err;
        if (e2 >= dy) {
            err += dy;
            x0 += sx;
        }
        if (e2 <= dx) {
            err += dx;
            y0 += sy;
        }
    }
    
    return linePoints;
}
};

struct VoronoiEdge {
    PointD start, end;
    VoronoiEdge(PointD s, PointD e) : start(s), end(e) {}
};

struct Triangle {
    PointD a, b, c;
    Triangle(PointD a_, PointD b_, PointD c_) : a(a_), b(b_), c(c_) {}

    // Добавляем оператор сравнения
    bool operator==(const Triangle& other) const {
        return (a == other.a && b == other.b && c == other.c) ||
               (a == other.a && b == other.c && c == other.b) ||
               (a == other.b && b == other.a && c == other.c) ||
               (a == other.b && b == other.c && c == other.a) ||
               (a == other.c && b == other.a && c == other.b) ||
               (a == other.c && b == other.b && c == other.a);
    }
    static bool areCollinear(const PointD& a, const PointD& b, const PointD& c) {
    double area = (b.y - a.y) * (c.x - b.x) - (c.y - b.y) * (b.x - a.x);
    return std::fabs(area) < Constants::EPSILON; // Проверяем, близко ли значение к нулю
}

     PointD calculateCircumcenter() const {
        // Реализация вычисления центра окружности
        double d = 2 * (a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y));
        if (std::fabs(d) < Constants::EPSILON) return PointD(); // Защита от деления на ноль
        double x = ((a.x*a.x + a.y*a.y)*(b.y - c.y) + (b.x*b.x + b.y*b.y)*(c.y - a.y) + (c.x*c.x + c.y*c.y)*(a.y - b.y)) / d;
        double y = ((a.x*a.x + a.y*a.y)*(c.x - b.x) + (b.x*b.x + b.y*b.y)*(a.x - c.x) + (c.x*c.x + c.y*c.y)*(b.x - a.x)) / d;
        return PointD(x, y);
    }

    static double distance(const PointD& p1, const PointD& p2) {
        return std::hypot(p1.x - p2.x, p1.y - p2.y);
    }

   static const Triangle* findContainingTriangle(const PointD& p, const std::vector<Triangle>& triangles) {
        for (const auto& tri : triangles) {
            if (isPointInTriangle(p, tri)) {
                return &tri;
            }
        }
        return nullptr;
    }
    
   static bool isPointInTriangle(const PointD& p, const Triangle& tri) {
        auto sign = [](PointD p1, PointD p2, PointD p3) {
            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
        };
        
        double d1 = sign(p, tri.a, tri.b);
        double d2 = sign(p, tri.b, tri.c);
        double d3 = sign(p, tri.c, tri.a);
        
        bool has_neg = (d1 < -Constants::EPSILON) || (d2 < -Constants::EPSILON) || (d3 < -Constants::EPSILON);
        bool has_pos = (d1 > Constants::EPSILON) || (d2 > Constants::EPSILON) || (d3 > Constants::EPSILON);
        
        bool result = !(has_neg && has_pos);
        return result;
    }
    
   static std::vector<const Triangle*> getNeighbors(const Triangle& tri, const std::vector<Triangle>& triangles) {
    
    std::vector<const Triangle*> neighbors;

    for (const auto& other : triangles) {
        if (&tri == &other) {
            continue;
        }
        
        if (shareEdge(tri, other)) {
            neighbors.push_back(&other);
        }
    }
    
    return neighbors;
}

    static std::vector<const Triangle*> findNeighbors(const Triangle& tri, const std::vector<Triangle>& allTriangles) {
        std::vector<const Triangle*> neighbors;
        
        for (const auto& other : allTriangles) {
            if (&tri == &other) continue;
            if (shareEdge(tri, other)) {
                neighbors.push_back(&other);
            }
        }
        
        return neighbors;
    }
    
static bool shareEdge(const Triangle& a, const Triangle& b) {
    
    const std::array<Edge, 3> edgesA = { Edge(a.a, a.b), Edge(a.b, a.c), Edge(a.c, a.a) };
    const std::array<Edge, 3> edgesB = { Edge(b.a, b.b), Edge(b.b, b.c), Edge(b.c, b.a) };

    for (size_t i = 0; i < 3; ++i) {
        const Edge& edgeA = edgesA[i];
        for (size_t j = 0; j < 3; ++j) {
            const Edge& edgeB = edgesB[j];
            if (edgeA == edgeB) {
                return true;
            }
        }
    }
    return false;
}

static bool otherHasEdge(const Triangle& other, const Edge& edge) {
    const Edge edges[3] = { Edge(other.a, other.b), Edge(other.b, other.c), Edge(other.c, other.a) };
    
    for (int i = 0; i < 3; ++i) {
        if (edges[i] == edge) {
            return true;
        }
    }
    return false;
}
};
