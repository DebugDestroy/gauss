#pragma once
#include "algorithms/geometry/geometry_structures.hpp"
#include <vector>

namespace algorithms::geometry {
double distance(const PointD& p1, const PointD& p2); // Возвращает расстояние между двумя точками

double length(const Edge& e); // Длина ребра (расстояние между концами)

PointD midPoint(const Edge& e); // Середина ребра — типа точка между двумя концами

PointD perpendicular(const Edge& e); // Перпендикулярный вектор к ребру (не нормированный)

bool areCollinear(const PointD& a, const PointD& b, const PointD& c); // Проверяет, лежат ли 3 точки на одной прямой

bool isPointInsideField(const PointD& p, int width, int height); // Проверяет, внутри ли поле p (размеры задаются)

bool shareEdge(const Triangle& a, const Triangle& b); // Проверяет, есть ли общее ребро у двух треугольников

bool otherHasEdge(const Triangle& other, const Edge& edge); // Есть ли такое ребро в другом треугольнике

bool isPointInTriangle(const PointD& p, const Triangle& tri); // Внутри ли точка p заданного треугольника

const Triangle* findContainingTriangle(const PointD& p, const std::vector<Triangle>& triangles); // Находит треугольник, внутри которого точка

std::vector<const Triangle*> findNeighbors(const Triangle& tri, const std::vector<Triangle>& allTriangles); // Возвращает соседей по общему ребру

bool isPointInCircumcircle(const PointD& p, const Triangle& tri); // Проверяет, лежит ли точка внутри описанной окружности треугольника

PointD calculateCircumcenter(const PointD& a, const PointD& b, const PointD& c); // Центр описанной окружности треугольника ABC
}
