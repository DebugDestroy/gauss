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

namespace visualization {

class GnuplotInterface {
private:
    core::Logger& logger;

    void applyNiceStyle(FILE* pipe, const std::string& title, bool is3D = false, const std::string& terminal = "pngcairo");
    double transformY(double y, int height) const;
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
                     
    void plotDelaunay(const std::vector<algorithms::geometry::Triangle>& triangles, 
                     const std::vector<std::vector<double>>& binaryMap, 
                     const std::string& filename);

    void plotPath(const std::vector<algorithms::geometry::Pixel>& path, 
                 const std::vector<std::vector<double>>& field,
                 const std::vector<std::vector<double>>& binaryMap,
                 const std::string& filename, 
                 const command::DispatcherParams& params,
                 const int Radius);

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
