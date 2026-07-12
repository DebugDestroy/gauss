#pragma once

#include <vector>
#include <array>
#include <memory>
#include <string>
#include <random>

// Локальные заголовки
// command
#include "command/dispatcher_params.hpp"
#include "command/application_state.hpp"

// core
#include "core/logger.hpp"

// io
#include "io/bmp_handler.hpp"

// visualization
#include "visualization/gnuplot.hpp"
#include "visualization/color.hpp"

// statistics
#include "statistics/statistics_manager.hpp"
#include "statistics/csv_writer.hpp"

// algorithms::gauss
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
#include "algorithms/geometry/math.hpp"

// algorithms::path::common
#include "algorithms/path/common/path_validator.hpp"
#include "algorithms/path/common/graph.hpp"
#include "algorithms/path/common/grid.hpp"
#include "algorithms/path/common/path_metrics.hpp"

// algorithms::path::a_star
#include "algorithms/path/a_star/path_finder.hpp"

// algorithms::path::dekstra
#include "algorithms/path/dekstra/path_finder.hpp"

// algorithms::path::greedy
#include "algorithms/path/greedy/path_finder.hpp"

// algorithms::path::rrt
#include "algorithms/path/rrt/rrt.hpp"

// algorithms::path::rrt_star
#include "algorithms/path/rrt_star/rrt_star.hpp"

namespace command {

class Control {
private:

    void logOperation(core::LogLevel level, const std::string& operation, const std::string& details = "");
    
    template<typename Point>
    static std::string formatPathMetricsLog(
        const algorithms::path::PathMetrics& m,
        const Point& start,
        const Point& end)
    {
        return std::string("from (") +
               std::to_string(start.x) + "," +
               std::to_string(start.y) + ")" +

               " to (" +
               std::to_string(end.x) + "," +
               std::to_string(end.y) + ")" +

               " | env=" + m.environment +
               " | algo=" + m.algorithmName +
     
               " | path_found=" + (m.pathFound ? "true" : "false") +
               " | time=" + std::to_string(m.executionTimeMs) + " ms" +

               " | nodes=" + std::to_string(m.pathNodes) +
               " | expanded=" + std::to_string(m.expandedNodes) +

               " | euclid_len=" + std::to_string(m.euclideanLength) +
               " | pixel_len=" + std::to_string(m.pixelLength) +

               " | min_obs_euc=" + std::to_string(m.minObstacleDistance) +
               " | min_obs_px=" + std::to_string(m.minObstacleDistancePixel) +

               " | max_side_angle=" + std::to_string(m.maxSideAngle) +
               " | max_updown_angle=" + std::to_string(m.maxUpDownAngle);
    }
public:
    // Генератор случайных чисел
    std::mt19937 randomGenerator;
    
    // Состояние приложения
    ApplicationState state;
    
    // core
    core::Logger& logger;
    
    // io
    io::BmpHandler bmpHandler;
    
    // visualization
    visualization::GnuplotInterface gnuplotInterface;
    visualization::ColorGenerator colorGenerator;
    
    // statistics
    statistics::StatisticsManager statisticsManager;
    statistics::CsvWriter csvWriter;
    
    // algorithms::gauss
    algorithms::gauss::GaussBuilder gaussBuilder;
    
    // algorithms::components
    algorithms::components::Copier copier;
    algorithms::components::Binarizer binarizer;
    algorithms::components::ComponentCalculator componentCalculator;
    algorithms::components::ClusterService clusterService;
    algorithms::components::KMeans kMeans;
    
    // algorithms::geometry
    algorithms::geometry::Triangulator triangulator;
    algorithms::geometry::VoronoiDiagram voronoi;
    
    // algorithms::path::common
    algorithms::path::common::PathValidator conditions;
    algorithms::path::common::Graph graph;
    algorithms::path::common::GridBuilder grid;
    
    // algorithms::path::a_star
    algorithms::path::a_star::PathFinder astarFinder;
    
    // algorithms::path::dekstra
    algorithms::path::dekstra::PathFinder dekstraFinder;
    
    // algorithms::path::greedy
    algorithms::path::greedy::PathFinder greedyFinder;
    
     // algorithms::path::rrt
    algorithms::path::rrt::PathFinder rrt;
    
    // algorithms::path::rrt_star
    algorithms::path::rrt_star::PathFinder rrtStar;
    
    // Конструктор
    Control(core::Logger& log, const std::string seedMode, uint32_t seed);

    // Основной диспетчер команд
    void Dispetcher(command::DispatcherParams& params);
};

}
