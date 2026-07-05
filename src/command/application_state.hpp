#pragma once

#include <vector>
#include <array>
#include <memory>
#include <optional>

// algorithms::gauss
#include "algorithms/gauss/gauss_builder.hpp"

// algorithms::components
#include "algorithms/components/find_components.hpp"

// algorithms::geometry
#include "algorithms/geometry/geometry_structures.hpp"

// algorithms::path::common
#include "algorithms/path/common/path_metrics.hpp" // Для метрик
#include "algorithms/path/common/grid.hpp"

// utils
#include "utils/hash.hpp" // Для навигационного графа

namespace command {

struct ApplicationState {

    // Визуализация
    std::vector<std::array<int, 3>> colors;

    // Бинарное представление поля
    std::vector<std::vector<double>> binaryMap;

    // Гауссово поле
    std::vector<std::vector<double>> field;
    std::vector<algorithms::gauss::Gaus> gaussi;

    // Компоненты связности
    std::vector<algorithms::components::Component> components;

    // Кластеризация
    std::vector<algorithms::geometry::PointD> clusterCenters;
    std::vector<algorithms::geometry::PointD> kmeansCenters;
    std::vector<std::vector<double>> kmeansField;

    // Геометрия
    std::vector<algorithms::geometry::Triangle> lastTriangulation;
    std::vector<algorithms::geometry::Edge> voronoiEdges;
    
    // Навигационный граф
    std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>> navigationGraph;
    
    // Сетка
    algorithms::path::common::Grid grid;
    
    // Поиск пути
    std::optional<algorithms::geometry::Pixel> start;
    std::optional<algorithms::geometry::Pixel> end;
    algorithms::path::PathMetrics PathMetrics;
    std::vector<algorithms::geometry::Pixel> path;
    
    // Поиск пути на сетке
    std::optional<algorithms::path::common::GridCell> startCell;
    std::optional<algorithms::path::common::GridCell> endCell;
    std::vector<algorithms::path::common::GridCell> gridPath;
};

} // namespace command
