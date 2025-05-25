#include "algorithms/geometry/math.hpp"
#include "core/constants.hpp"
#include <cmath>
#include <array>

namespace algorithms::geometry {
    double distance(const PointD& p1, const PointD& p2) {
        return std::hypot(p1.x - p2.x, p1.y - p2.y);
}
     double length(const Edge& e) {
    return distance(e.a, e.b);
}

 PointD midPoint(const Edge& e) {
    return PointD((e.a.x + e.b.x) / 2.0, (e.a.y + e.b.y) / 2.0);
}

 PointD perpendicular(const Edge& e) {
    return PointD(-(e.b.y - e.a.y), e.b.x - e.a.x);
}

    bool areCollinear(const PointD& a, const PointD& b, const PointD& c) {
    double area = (b.y - a.y) * (c.x - b.x) - (c.y - b.y) * (b.x - a.x);
    return std::fabs(area) < core::EPSILON; // Проверяем, близко ли значение к нулю
}

    bool isPointInsideField(const PointD& p, int width, int height) {
        bool inside = p.x >= 0 && p.y >= 0 && p.x < width && p.y < height;
        return inside;
    }

     bool shareEdge(const Triangle& a, const Triangle& b) {
    
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

    bool otherHasEdge(const Triangle& other, const Edge& edge) {
    const Edge edges[3] = { Edge(other.a, other.b), Edge(other.b, other.c), Edge(other.c, other.a) };
    
    for (int i = 0; i < 3; ++i) {
        if (edges[i] == edge) {
            return true;
        }
    }
    return false;
}
   
    bool isPointInTriangle(const PointD& p, const Triangle& tri) {
        auto sign = [](PointD p1, PointD p2, PointD p3) {
            return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y);
        };
        
        double d1 = sign(p, tri.a, tri.b);
        double d2 = sign(p, tri.b, tri.c);
        double d3 = sign(p, tri.c, tri.a);
        
        bool has_neg = (d1 < -core::EPSILON) || (d2 < -core::EPSILON) || (d3 < -core::EPSILON);
        bool has_pos = (d1 > core::EPSILON) || (d2 > core::EPSILON) || (d3 > core::EPSILON);
        
        bool result = !(has_neg && has_pos);
        return result;
    }
    
    const Triangle* findContainingTriangle(const PointD& p, const std::vector<Triangle>& triangles) {
        for (const auto& tri : triangles) {
            if (isPointInTriangle(p, tri)) {
                return &tri;
            }
        }
        return nullptr;
    }

    std::vector<const Triangle*> findNeighbors(const Triangle& tri, const std::vector<Triangle>& allTriangles) {
        std::vector<const Triangle*> neighbors;
        
        for (const auto& other : allTriangles) {
            if (&tri == &other) continue;
            if (shareEdge(tri, other)) {
                neighbors.push_back(&other);
            }
        }
        
        return neighbors;
    }

    bool isPointInCircumcircle(const PointD& p, const Triangle& tri) {
    const PointD& a = tri.a;
    const PointD& b = tri.b;
    const PointD& c = tri.c;

    // 1. Проверка на вырожденный треугольник
    if (areCollinear(a, b, c)) {
        // Вырожденный треугольник - проверка окружности невозможна
        return false;
    }

    // 2. Вычисление центра окружности
    const double D = 2 * (a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y));
    if (std::fabs(D) < core::EPSILON) {
        //Ошибка вычисления центра окружности
        return false;
    }

    const PointD center = {
        ((a.x*a.x + a.y*a.y) * (b.y - c.y) + 
         (b.x*b.x + b.y*b.y) * (c.y - a.y) + 
         (c.x*c.x + c.y*c.y) * (a.y - b.y)) / D,
         
        ((a.x*a.x + a.y*a.y) * (c.x - b.x) + 
         (b.x*b.x + b.y*b.y) * (a.x - c.x) + 
         (c.x*c.x + c.y*c.y) * (b.x - a.x)) / D
    };

    // 3. Вычисление квадрата радиуса
    const double radius_sq = (a.x - center.x)*(a.x - center.x) + 
                           (a.y - center.y)*(a.y - center.y);

    // 4. Вычисление квадрата расстояния от точки до центра
    const double dist_sq = (p.x - center.x)*(p.x - center.x) + 
                         (p.y - center.y)*(p.y - center.y);

    // 5. Epsilon-сравнение
    return dist_sq < radius_sq - core::EPSILON;
}

 PointD calculateCircumcenter(const PointD& a, const PointD& b, const PointD& c) {
    double d = 2 * (a.x * (b.y - c.y) +
                    b.x * (c.y - a.y) +
                    c.x * (a.y - b.y));
    
    if (std::fabs(d) < core::EPSILON)
        return PointD(); // Центр не определён — вершины почти на одной прямой

    double x = ((a.x * a.x + a.y * a.y) * (b.y - c.y) +
                (b.x * b.x + b.y * b.y) * (c.y - a.y) +
                (c.x * c.x + c.y * c.y) * (a.y - b.y)) / d;

    double y = ((a.x * a.x + a.y * a.y) * (c.x - b.x) +
                (b.x * b.x + b.y * b.y) * (a.x - c.x) +
                (c.x * c.x + c.y * c.y) * (b.x - a.x)) / d;

    return PointD(x, y);
}
}
