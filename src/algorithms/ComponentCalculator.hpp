#pragma once
#include <vector>
#include <algorithm> // для std::remove_if, std::find, std::all_of
#include <cmath> // для std::isnan, std::sqrt
#include <numeric> // для std::iota
#include <string>
#include <sstream>
#include <limits> // для std::numeric_limits

// Локальные заголовки
#include "core/Logger.hpp"
#include "core/Geometry.hpp" // для PointD, Triangle, Edge
#include "core/Pole.hpp" // для std::unique_ptr<Pole>
#include "algorithms/Component.hpp" // для Component, ThresholdMode

class ComponentCalculator {
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

    void logPoint(const PointD& p, const std::string& prefix = "") const {
        logger.trace(prefix + "Точка (" + std::to_string(p.x) + "," + std::to_string(p.y) + ")");
    }

    void logEdge(const Edge& edge, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Ребро: (" << edge.a.x << "," << edge.a.y << ")-(" 
            << edge.b.x << "," << edge.b.y << ")";
        logger.trace(oss.str());
    }

public:
    ComponentCalculator(Logger& lg) : logger(lg) {
        logger.trace("[ComponentCalculator] Инициализация калькулятора компонент");
    }

    std::vector<Triangle> bowyerWatson(const std::vector<PointD>& points) {
        logger.info("[ComponentCalculator::bowyerWatson] Начало триангуляции Боуера-Ватсона");
        logger.debug(std::string("Количество точек: ") + std::to_string(points.size()));
        
        std::vector<Triangle> triangles;

        // Проверка на минимальное количество точек
        if (points.size() < 3) {
            logger.error("[ComponentCalculator::bowyerWatson] Недостаточно точек для триангуляции (меньше 3)");
            return triangles;
        }

        // Проверка на NaN-точки
        for (const auto& p : points) {
            if (std::isnan(p.x) || std::isnan(p.y)) {
                logger.error("[ComponentCalculator::bowyerWatson] Обнаружена некорректная точка (NaN)");
                return triangles;
            }
        }

        // Проверка на коллинеарность
        if (std::all_of(points.begin(), points.end(), 
                       [&](const PointD& p) { return Triangle::areCollinear(points[0], points[1], p); })) {
            logger.error("[ComponentCalculator::bowyerWatson] Все точки коллинеарны, триангуляция невозможна");
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
        logger.info("[ComponentCalculator::bowyerWatson] Создан супер-треугольник");
        logTriangle(triangles.back(), "Супер-треугольник: ");

        // Основной алгоритм
        size_t pointIndex = 0;
        for (const auto& point : points) {
            logger.trace(std::string("[ComponentCalculator::bowyerWatson] Обработка точки #") + 
                        std::to_string(++pointIndex) + " (" + 
                        std::to_string(point.x) + "," + std::to_string(point.y) + ")");

            std::vector<Triangle> badTriangles;
            std::copy_if(triangles.begin(), triangles.end(), std::back_inserter(badTriangles),
                         [&](const Triangle& tri) { return isPointInCircumcircle(point, tri); });

            logger.debug(std::string("[ComponentCalculator::bowyerWatson] Найдено ") + 
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

        logger.info("[ComponentCalculator::bowyerWatson] Триангуляция завершена");
        logger.debug(std::string("Удалено ") + std::to_string(initialCount - triangles.size()) + 
                    " треугольников, связанных с супер-треугольником");
        logger.debug(std::string("Итоговое количество треугольников: ") + std::to_string(triangles.size()));

        return triangles;
    }

    bool isPointInCircumcircle(const PointD& p, const Triangle& tri) {
        double ax = tri.a.x - p.x, ay = tri.a.y - p.y;
        double bx = tri.b.x - p.x, by = tri.b.y - p.y;
        double cx = tri.c.x - p.x, cy = tri.c.y - p.y;
        
        double det = ax * (by * (cx*cx + cy*cy) - cy * (bx*bx + by*by)) 
                   - ay * (bx * (cx*cx + cy*cy) - cx * (bx*bx + by*by)) 
                   + (ax*ax + ay*ay) * (bx*cy - by*cx);

        bool inside = det > 0;
        logger.trace(std::string("[ComponentCalculator::isPointInCircumcircle] Точка (") + 
                    std::to_string(p.x) + "," + std::to_string(p.y) + ") " +
                    (inside ? "внутри" : "вне") + " окружности треугольника");
        logTriangle(tri, "Проверяемый треугольник: ");

        return inside;
    }

    bool hasEdge(const Triangle& tri, const Edge& edge) {
        bool has = Edge(tri.a, tri.b) == edge 
                || Edge(tri.b, tri.c) == edge 
                || Edge(tri.c, tri.a) == edge;

        logger.trace(std::string("[ComponentCalculator::hasEdge] Треугольник ") + 
                    std::string(has ? "содержит" : "не содержит") + " ребро");
        logTriangle(tri);
        logEdge(edge);

        return has;
    }
       
    int incrementAndCollect(std::vector<std::vector<double>>& componenta, 
                          std::vector<std::vector<double>>& CopyPole, 
                          int x, int y, int i, int& pixelCount) {
        logger.trace(std::string("[ComponentCalculator::incrementAndCollect] Проверка пикселя (") + 
                    std::to_string(x) + "," + std::to_string(y) + "), глубина=" + 
                    std::to_string(i));

        if (x < 1 || y < 1 || x > (int)componenta[0].size() - 2 || 
            y > (int)componenta.size() - 2 || CopyPole[y][x] < 250) {
            logger.trace("[ComponentCalculator::incrementAndCollect] Выход за границы или неподходящее значение");
            return -1;
        }

        if (CopyPole[y][x] >= 255 && CopyPole[y][x] <= 255) {
            CopyPole[y][x] = 0;
            pixelCount++;
            componenta[y][x] = 255;
            
            logger.trace(std::string("[ComponentCalculator::incrementAndCollect] Обработан пиксель (") + 
                        std::to_string(x) + "," + std::to_string(y) + "), счетчик=" + 
                        std::to_string(pixelCount));

            // Рекурсивные вызовы
            incrementAndCollect(componenta, CopyPole, x + 1, y, i + 1, pixelCount);
            incrementAndCollect(componenta, CopyPole, x - 1, y, i + 1, pixelCount);
            incrementAndCollect(componenta, CopyPole, x, y + 1, i + 1, pixelCount);
            incrementAndCollect(componenta, CopyPole, x, y - 1, i + 1, pixelCount);
        }
        
        return pixelCount;
    }

    void bin(std::vector<std::vector<double>>& CopyPole, 
            int slice, 
            std::unique_ptr<Pole>& p, 
            ThresholdMode mode) {
        logger.info(std::string("[ComponentCalculator::bin] Бинаризация данных, slice=") + 
                  std::to_string(slice) + ", mode=" + 
                  (mode == ThresholdMode::Peaks ? "Peaks" : 
                   mode == ThresholdMode::Valleys ? "Valleys" : "All"));

        if (p == nullptr) {
            logger.error("[ComponentCalculator::bin] Ошибка: данные высот не инициализированы!");
            return;
        }
        
        CopyPole = p->field;
        int symmetric_slice = 2*127 - slice;

        for (int x = 0; x < (int)p->field[0].size(); ++x) {
            for (int y = 0; y < (int)p->field.size(); ++y) {
                double value = p->field[y][x];
                switch (mode) {
                    case ThresholdMode::All: {
                        bool is_peak = (slice >= 127) ? (value > slice) : (value > symmetric_slice);
                        bool is_valley = (slice >= 127) ? (value < symmetric_slice) : (value < slice);
                        CopyPole[y][x] = (is_peak || is_valley) ? 255 : 0;
                        break;
                    }
                    case ThresholdMode::Peaks:
                        CopyPole[y][x] = (slice >= 127) ? (value > slice) : (value > symmetric_slice);
                        break;
                    case ThresholdMode::Valleys:
                        CopyPole[y][x] = (slice >= 127) ? (value < symmetric_slice) : (value < slice);
                        break;
                }
            } 
        }

        logger.info("[ComponentCalculator::bin] Бинаризация завершена");
    }
    
    void wave(int noisy, 
             std::vector<Component>& componenti, 
             std::vector<std::vector<double>>& CopyPole, 
             std::unique_ptr<Pole>& p) {
        logger.info(std::string("[ComponentCalculator::wave] Начало волнового алгоритма, порог=") + 
                  std::to_string(noisy));

        if (p == nullptr) {
            logger.error("[ComponentCalculator::wave] Ошибка: данные высот не инициализированы!");
            return;
        }

        int rows = p->field.size();
        int cols = (rows > 0) ? p->field[0].size() : 0;
        std::vector<Component> noiseComponents;

        logger.debug(std::string("Размер данных: ") + std::to_string(cols) + "x" + std::to_string(rows));

        for (int y = 2; y < rows - 2; ++y) {
            for (int x = 2; x < cols - 2; ++x) {
                if (CopyPole[y][x] <= 255 && CopyPole[y][x] >= 255) {
                    std::vector<std::vector<double>> componentData(rows, std::vector<double>(cols, 0));
                    int pixelCount = 0;
                    incrementAndCollect(componentData, CopyPole, x, y, 0, pixelCount);

                    if (pixelCount >= noisy) {
                        Component component(logger, componentData, pixelCount);
                        componenti.push_back(component);
                        
                        logger.debug(std::string("[ComponentCalculator::wave] Значимая компонента: ") +
                                   "pixels=" + std::to_string(pixelCount) +
                                   ", center=(" + std::to_string(component.center_x) + "," + 
                                   std::to_string(component.center_y) + ")" +
                                   ", size=(" + std::to_string(component.max_x - component.min_x) + 
                                   "x" + std::to_string(component.max_y - component.min_y) + ")");
                    } else {
                        Component noiseComponent(logger, componentData, pixelCount);
                        noiseComponents.push_back(noiseComponent);
                        
                        logger.trace(std::string("[ComponentCalculator::wave] Шумовая компонента: ") +
                                    "pixels=" + std::to_string(pixelCount) +
                                    ", center=(" + std::to_string(noiseComponent.center_x) + "," + 
                                    std::to_string(noiseComponent.center_y) + ")");

                        // Удаляем шум из основного поля
                        for (int i = 0; i < rows; ++i) {
                            for (int j = 0; j < cols; ++j) {
                                if (componentData[i][j] <= 255 && componentData[i][j] >= 255) {
                                    p->field[i][j] = 127;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        logger.info("[ComponentCalculator::wave] Обработка завершена");
        logger.debug(std::string("Найдено значимых компонент: ") + std::to_string(componenti.size()) +
                   ", шумовых компонент: " + std::to_string(noiseComponents.size()));
    }
};
