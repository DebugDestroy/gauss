#pragma once

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>
#include <optional>

#include "command/dispatcher_params.hpp"
#include "core/logger.hpp"
#include "algorithms/components/components_analysis.hpp" // Для Component
#include "algorithms/geometry/geometry_structures.hpp"   // Для PointD, Edge
#include "algorithms/geometry/bresenham_line.hpp"        // Для bresenhamLine
#include "algorithms/path/common/grid.hpp"               // Для сетки
#include "algorithms/path/rrt/rrt.hpp"                   // Для rrt

namespace visualization {

class GnuplotInterface {
private:
    core::Logger& logger;

    void applyNiceStyle(FILE* pipe, const std::string& title, bool is3D = false, const std::string& terminal = "pngcairo");
    void logPlotStart(const std::string& plotType, const std::string& filename) const;
    void logPlotEnd(const std::string& plotType) const;

public:
    explicit GnuplotInterface(core::Logger& lg);

    void plotBinaryWithComponents(const std::vector<std::vector<double>>& binaryMap,
                                   const std::vector<algorithms::components::Component>& components,
                                   const std::string& filename);

    void gnuplot(const std::vector<std::vector<double>>& field, const std::string& filename);
    
    void plotKmeans(const std::vector<std::vector<double>>& binaryMap, const std::vector<algorithms::geometry::PointD>& centers, const std::string& filename);

    void plotVoronoi(const std::vector<std::vector<double>>& binaryMap,
                     const std::vector<algorithms::geometry::Edge>& edges,
                     const std::vector<algorithms::geometry::PointD>& sites,
                     const std::string& filename);

    void plotGraph(const std::vector<std::vector<double>>& binaryMap,
                   const std::unordered_map<algorithms::geometry::Pixel, std::vector<algorithms::geometry::Pixel>>& graph,
                   const std::string& filename,
                   std::optional<algorithms::geometry::Pixel> start = std::nullopt,
                   std::optional<algorithms::geometry::Pixel> goal = std::nullopt);
    
    void plotGrid(
    const algorithms::path::common::Grid& grid,
    const std::vector<std::vector<double>>& binaryMap,
    const std::string& filename);
    
    void plotNavGrid(
    const algorithms::path::common::Grid& grid,
    const std::vector<std::vector<double>>& binaryMap,
    const std::string& filename,
    const std::optional<algorithms::path::common::GridCell> start = std::nullopt,
    const std::optional<algorithms::path::common::GridCell> end = std::nullopt);
    
    void plotGridPath(
    const algorithms::path::common::Grid& grid,
    std::vector<algorithms::path::common::GridCell>& gridPath,
    const algorithms::path::common::GridCell& start,
    const algorithms::path::common::GridCell& end,
    const std::vector<std::vector<double>>& binaryMap,
    const std::string& filename);
                     
    void plotDelaunay(const std::vector<algorithms::geometry::Triangle>& triangles, 
                     const std::vector<std::vector<double>>& binaryMap, 
                     const std::string& filename);

    void plotPath(const std::vector<algorithms::geometry::Pixel>& path, 
                 const std::vector<std::vector<double>>& field,
                 const std::vector<std::vector<double>>& binaryMap,
                 const std::string& filename, 
                 const command::DispatcherParams& params,
                 const int Radius);
    
    void plotRRT(const std::vector<algorithms::geometry::PointD>& path, 
                 const std::vector<algorithms::path::rrt::RRTNode>& treeRRT,
                 const algorithms::geometry::PointD& start, 
                 const algorithms::geometry::PointD& end,
                 const std::vector<std::vector<double>>& binaryMap,
                 const std::string& filename);
                 
    void plotInteractive3DPath(const std::vector<algorithms::geometry::Pixel>& path, 
                              const std::vector<std::vector<double>>& field, 
                              const algorithms::geometry::Pixel& start, 
                              const algorithms::geometry::Pixel& end,
                              const int sphereRadius);

    void plot3DPath(const std::vector<algorithms::geometry::Pixel>& path, 
                   const std::vector<std::vector<double>>& field, 
                   const std::string& filename, 
                   const algorithms::geometry::Pixel& start, 
                   const algorithms::geometry::Pixel& end,
                   const int sphereRadius);
};

}
