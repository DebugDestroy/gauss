#pragma once

#include <vector>
#include <string>
#include <memory>

#include "command/dispatcher_params.hpp"
#include "core/logger.hpp"
#include "algorithms/components/components_analysis.hpp" // Для Component
#include "algorithms/gauss/pole.hpp"                     // Для Pole
#include "algorithms/geometry/geometry_structures.hpp"   // Для PointD, Edge
#include "algorithms/geometry/bresenham_line.hpp"        // Для bresenhamLine

namespace visualization {

class GnuplotInterface {
private:
    core::Logger& logger;

    double transformY(double y, int height) const;
    void logPlotStart(const std::string& plotType, const std::string& filename) const;
    void logPlotEnd(const std::string& plotType) const;

public:
    explicit GnuplotInterface(core::Logger& lg);

    void plotBinaryWithComponents(const std::vector<std::vector<double>>& CopyPole,
                                   const std::vector<algorithms::components::Component>& components,
                                   const std::string& filename);

    void gnuplot(std::unique_ptr<algorithms::gauss::Pole>& p, const std::string& filename);

    void plotVoronoi(const std::unique_ptr<algorithms::gauss::Pole>& p,
                     const std::vector<algorithms::geometry::Edge>& edges,
                     const std::vector<algorithms::geometry::PointD>& sites,
                     const std::string& filename);

    void plotDelaunay(const std::vector<algorithms::geometry::Triangle>& triangles, 
                     std::unique_ptr<algorithms::gauss::Pole>& p, 
                     const std::string& filename);

    void plotPath(const std::vector<algorithms::geometry::PointD>& path, 
                 const std::unique_ptr<algorithms::gauss::Pole>& p, 
                 const std::string& filename, 
                 const command::DispatcherParams& params,
                 const double Radius);

    void plotInteractive3DPath(const std::vector<algorithms::geometry::PointD>& path, 
                              const std::unique_ptr<algorithms::gauss::Pole>& p, 
                              const algorithms::geometry::PointD& start, 
                              const algorithms::geometry::PointD& end,
                              const double sphereRadius);

    void plot3DPath(const std::vector<algorithms::geometry::PointD>& path, 
                   const std::unique_ptr<algorithms::gauss::Pole>& p, 
                   const std::string& filename, 
                   const algorithms::geometry::PointD& start, 
                   const algorithms::geometry::PointD& end,
                   const double sphereRadius);
};

}
