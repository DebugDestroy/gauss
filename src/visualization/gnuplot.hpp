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

    void plotPathDiscrete(const std::vector<algorithms::geometry::Pixel>& path, 
                 const std::vector<std::vector<double>>& field,
                 const std::vector<std::vector<double>>& binaryMap,
                 const std::string& filename, 
                 const command::DispatcherParams& params,
                 const int Radius);
    
    void plotPathContinuous(
        const std::vector<algorithms::geometry::PointD>& path,
        const std::vector<algorithms::gauss::Gaus>& gaussi,
        int fieldWidth,
        int fieldHeight,
        double heightThreshold,
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

template<typename NodeType>
    void plotRRT(const std::vector<algorithms::geometry::PointD>& path, 
                 const std::vector<NodeType>& treeRRT,
                 const algorithms::geometry::PointD& start, 
                 const algorithms::geometry::PointD& end,
                 const std::vector<std::vector<double>>& binaryMap,
                 const std::string& filename)
{
    logPlotStart("RRT Animation", filename);

    if (binaryMap.empty() || treeRRT.empty())
    {
        logger.warning("[plotRRT] Empty data");
        return;
    }

    FILE* pipe = popen("gnuplot", "w");

    if (!pipe)
    {
        logger.error("[plotRRT] Cannot start gnuplot");
        return;
    }

    const int height = static_cast<int>(binaryMap.size());
    const int width  = static_cast<int>(binaryMap[0].size());

    applyNiceStyle(pipe,
                   "Rapidly Exploring Random Tree",
                   false,
                   "gif");

    fprintf(pipe,
        "set output '%s'\n",
        filename.c_str());

    fprintf(pipe, "unset key\n");
    fprintf(pipe, "set size ratio -1\n");

    fprintf(pipe,
        "set xrange [0:%d]\n"
        "set yrange [0:%d]\n",
        width - 1,
        height - 1);

    fprintf(pipe,
        "set label 1 'Start' at %f,%f tc rgb 'blue' front\n",
        start.x,
        start.y);

    fprintf(pipe,
        "set label 2 'Goal' at %f,%f tc rgb 'purple' front\n",
        end.x,
        end.y);

    // ============================
    // Постепенное построение дерева
    // ============================
    
    std::size_t step = std::max<std::size_t>(1, treeRRT.size() / 100);
        
    for (std::size_t frame = 1; frame <= treeRRT.size(); frame += step)
    {
        fprintf(pipe,
            "plot '-' matrix with image,\\\n"
            "'-' with lines lw 2 lc rgb '#00AA00',\\\n"
            "'-' with points pt 7 ps 2 lc rgb 'blue',\\\n"
            "'-' with points pt 9 ps 2 lc rgb 'purple'\n");

        //-------------------------
        // Карта
        //-------------------------

        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
            {
                fprintf(pipe,
                        "%f ",
                        binaryMap[static_cast<std::size_t>(y)]
                                 [static_cast<std::size_t>(x)]);
            }
            fprintf(pipe, "\n");
        }
        fprintf(pipe, "e\n");

        //-------------------------
        // Дерево
        //-------------------------

        for (std::size_t i = 1; i < frame; ++i)
        {
            const auto& node = treeRRT[i];

            if (node.parent < 0)
                continue;

            const auto& parent =
                treeRRT[static_cast<std::size_t>(node.parent)];

            fprintf(pipe,
                "%f %f\n",
                parent.point.x,
                parent.point.y);

            fprintf(pipe,
                "%f %f\n\n",
                node.point.x,
                node.point.y);
        }
        
        fprintf(pipe, "e\n");

        //-------------------------
        // Старт
        //-------------------------

        fprintf(pipe,
            "%f %f\n"
            "e\n",
            start.x,
            start.y);

        //-------------------------
        // Финиш
        //-------------------------

        fprintf(pipe,
            "%f %f\n"
            "e\n",
            end.x,
            end.y);
    }

    // =====================================
    // Последний кадр — рисуем путь поверх
    // =====================================
    for (int pause = 0; pause < 50; ++pause) 
    {
        fprintf(pipe,
            "plot '-' matrix with image,\\\n"
            "'-' with lines lw 2 lc rgb '#00AA00',\\\n"
            "'-' with lines lw 4 lc rgb 'red',\\\n"
            "'-' with points pt 7 ps 2 lc rgb 'blue',\\\n"
            "'-' with points pt 9 ps 2 lc rgb 'purple'\n");

        // карта
        for (int y = 0; y < height; ++y)
        {
            for (int x = 0; x < width; ++x)
                fprintf(pipe,
                        "%f ",
                        binaryMap[static_cast<std::size_t>(y)]
                                 [static_cast<std::size_t>(x)]);
            fprintf(pipe, "\n");
        }
        fprintf(pipe, "e\n");

        // дерево
        for (std::size_t i = 1; i < treeRRT.size(); ++i)
        {
            if (treeRRT[i].parent < 0)
                continue;

            const auto& p =
                treeRRT[static_cast<std::size_t>(treeRRT[i].parent)];

            fprintf(pipe,
                "%f %f\n"
                "%f %f\n\n",
                p.point.x,
                p.point.y,
                treeRRT[i].point.x,
                treeRRT[i].point.y);
        }
        fprintf(pipe, "e\n");

        // путь
        for (const auto& p : path)
            fprintf(pipe,
                    "%f %f\n",
                    p.x,
                    p.y);
        fprintf(pipe, "e\n");

        // старт
        fprintf(pipe,
            "%f %f\n"
            "e\n",
            start.x,
            start.y);

        // финиш
        fprintf(pipe,
            "%f %f\n"
            "e\n",
            end.x,
            end.y);
    }
    pclose(pipe);

    logPlotEnd("RRT Animation");
}

};

}
