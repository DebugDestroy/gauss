#pragma once
#include <optional>  // Для std::optional
#include <cmath>    // Для std::fabs, std::hypot
#include <vector>   // Для использования в других классах
#include <string>   // Для возможного использования строк

// Локальные заголовки
#include "core/constants.hpp"  // Подключаем константы

namespace algorithms::geometry {
// Структура для точки с вещественными координатами (для точности) auto clusterCenters = getClusterCenters();
struct PointD {
    double x, y;
    PointD(double x_ = 0, double y_ = 0) : x(x_), y(y_) {}
    
    bool operator==(const PointD& other) const
{
    return std::abs(x - other.x) < core::EPSILON &&
           std::abs(y - other.y) < core::EPSILON;
}
    
    bool operator!=(const PointD& other) const {
        return !(*this == other);
    }
    
    PointD operator-(const PointD& other) const {
        return {x - other.x, y - other.y};
    }
    
    PointD& operator/=(double k) {
        x /= k;
        y /= k;
        return *this;
    }
};

struct Edge {
    PointD a, b;
    Edge(PointD a_, PointD b_) : a(a_), b(b_) {}
    
    bool operator==(const Edge& other) const {
        return (a == other.a && b == other.b) || (a == other.b && b == other.a);
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
};

struct Pixel {
    int x, y;

    Pixel(int x_ = 0, int y_ = 0)
        : x(x_), y(y_) {}

    bool operator==(const Pixel& other) const {
        return x == other.x && y == other.y;
    }

    bool operator!=(const Pixel& other) const {
        return !(*this == other);
    }
    
    Pixel operator-(const Pixel& other) const {
        return {x - other.x, y - other.y};
    }
};

struct PixelEdge
{
    Pixel a, b;

    PixelEdge(Pixel a_ = {}, Pixel b_ = {})
        : a(a_), b(b_) {}

    bool operator==(const PixelEdge& other) const
    {
        return (a == other.a && b == other.b) ||
               (a == other.b && b == other.a);
    }
};


}
