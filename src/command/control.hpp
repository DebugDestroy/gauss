#pragma once

#include <vector>
#include <array>
#include <memory>
#include <string>

// Локальные заголовки
// command
#include "command/dispatcher_params.hpp"

// core
#include "core/config.hpp"
#include "core/logger.hpp"

// io
#include "io/bmp_handler.hpp"

// visualization
#include "visualization/gnuplot.hpp"
#include "visualization/color.hpp"

// algorithms::gauss
#include "algorithms/gauss/pole.hpp"
#include "algorithms/gauss/gauss_builder.hpp"

// algorithms::components
#include "algorithms/components/copier_service.hpp"
#include "algorithms/components/binary.hpp"
#include "algorithms/components/find_components.hpp"
#include "algorithms/components/cluster_service.hpp"
#include "algorithms/components/kmeans.hpp"

// algorithms::geometry
#include "algorithms/geometry/triangulator.hpp"
#include "algorithms/geometry/voronoi_diagram.hpp"

// algorithms::path::common
#include "algorithms/path/common/conditions.hpp"
#include "algorithms/path/common/graph.hpp"

// algorithms::path::a_star
#include "algorithms/path/a_star/path_finder.hpp"

// algorithms::path::dekstra
#include "algorithms/path/dekstra/path_finder.hpp"

// algorithms::path::greedy
#include "algorithms/path/greedy/path_finder.hpp"

namespace command {

class Control {
private:

    void logOperation(core::LogLevel level, const std::string& operation, const std::string& details = "");

public:
    // Переменные
    std::vector<std::array<int, 3>> colors;  
    std::vector<std::vector<double>> CopyPole;
    
    // core
    core::Config& config;
    core::Logger& logger;
    
    // io
    io::BmpHandler bmpHandler;
    
    // visualization
    visualization::GnuplotInterface gnuplotInterface;
    visualization::ColorGenerator colorGenerator;
    
    // algorithms::gauss
    std::unique_ptr<algorithms::gauss::Pole> p = nullptr;
    algorithms::gauss::GaussBuilder gaussBuilder;
    std::vector<algorithms::gauss::Gaus> gaussi;
    
    // algorithms::components
    algorithms::components::Copier copier;
    algorithms::components::Binarizer binarizer;
    algorithms::components::ComponentCalculator componentCalculator;
    std::vector<algorithms::components::Component> componenti;
    algorithms::components::ClusterService clusterService;
    algorithms::components::KMeans kMeans;
    
    // algorithms::geometry
    std::vector<algorithms::geometry::PointD> clusterCenters;
    algorithms::geometry::Triangulator triangulator;
    std::vector<algorithms::geometry::Triangle> lastTriangulation;
    algorithms::geometry::VoronoiDiagram voronoi;
    std::vector<algorithms::geometry::Edge> voronoiEdges;
    algorithms::geometry::PointD start;
    algorithms::geometry::PointD end;
    std::vector<algorithms::geometry::PointD> path;
    
    // algorithms::path::common
    algorithms::path::common::Conditions conditions;
    algorithms::path::common::Graph graph;
    
    // algorithms::path::a_star
    algorithms::path::a_star::PathFinder astarFinder;
    
    // algorithms::path::dekstra
    algorithms::path::dekstra::PathFinder dekstraFinder;
    
    // algorithms::path::greedy
    algorithms::path::greedy::PathFinder greedyFinder;
    
    // Конструктор
    Control(core::Config& cfg, core::Logger& log);

    // Основной диспетчер команд
    void Dispetcher(command::DispatcherParams& params);
};

}
