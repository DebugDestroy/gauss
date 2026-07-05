#pragma once

#include <vector>
#include <array>
#include <memory>
#include <string>

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
#include "algorithms/path/common/conditions.hpp"
#include "algorithms/path/common/graph.hpp"
#include "algorithms/path/common/grid.hpp"
#include "algorithms/path/common/path_metrics.hpp"

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
    
    std::string formatPathMetricsLog(
        const algorithms::path::PathMetrics& m,
        const algorithms::geometry::Pixel& start,
        const algorithms::geometry::Pixel& end);
public:
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
    algorithms::path::common::Conditions conditions;
    algorithms::path::common::Graph graph;
    algorithms::path::common::GridBuilder grid;
    
    // algorithms::path::a_star
    algorithms::path::a_star::PathFinder astarFinder;
    
    // algorithms::path::dekstra
    algorithms::path::dekstra::PathFinder dekstraFinder;
    
    // algorithms::path::greedy
    algorithms::path::greedy::PathFinder greedyFinder;
    
    // Конструктор
    Control(core::Logger& log);

    // Основной диспетчер команд
    void Dispetcher(command::DispatcherParams& params);
};

}
