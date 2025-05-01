#pragma once
#include <cmath>    // Для std::fabs, std::hypot
#include <vector>   // Для использования в других классах
#include <string>   // Для возможного использования строк

// Структура для точки с вещественными координатами (для точности) auto clusterCenters = getClusterCenters();
struct PointD {
    double x, y;
    PointD(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    
    bool operator==(const PointD& other) const {
        return std::fabs(x - other.x) < 1e-6 && std::fabs(y - other.y) < 1e-6;
    }
    
    bool operator!=(const PointD& other) const {
        return !(*this == other);
    }
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
    return std::fabs(area) < 1e-6; // Проверяем, близко ли значение к нулю
}

     PointD calculateCircumcenter() const {
        // Реализация вычисления центра окружности
        double d = 2 * (a.x*(b.y - c.y) + b.x*(c.y - a.y) + c.x*(a.y - b.y));
        if (std::fabs(d) < 1e-9) return PointD(); // Защита от деления на ноль
        double x = ((a.x*a.x + a.y*a.y)*(b.y - c.y) + (b.x*b.x + b.y*b.y)*(c.y - a.y) + (c.x*c.x + c.y*c.y)*(a.y - b.y)) / d;
        double y = ((a.x*a.x + a.y*a.y)*(c.x - b.x) + (b.x*b.x + b.y*b.y)*(a.x - c.x) + (c.x*c.x + c.y*c.y)*(b.x - a.x)) / d;
        return PointD(x, y);
    }

    static double distance(const PointD& p1, const PointD& p2) {
        return std::hypot(p1.x - p2.x, p1.y - p2.y);
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
};

struct VoronoiEdge {
    PointD start, end;
    VoronoiEdge(PointD s, PointD e) : start(s), end(e) {}
};
