#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include <stdexcept>
#include <sstream>
#include <ctime>
#include <chrono>
#include <iomanip>
#include <algorithm>
#include <random>
#include <queue>
#include <unordered_map>
#include <memory> // Для std::unique_ptr

// Локальные заголовки
#include "core/Config.hpp"
#include "core/Logger.hpp"
#include "command/DispatcherParams.hpp"  // Добавляем include
#include "algorithms/Copier.hpp"
#include "algorithms/Component.hpp"
#include "algorithms/GaussBuilder.hpp"
#include "algorithms/KMeans.hpp"
#include "algorithms/PathFinder.hpp"
#include "algorithms/TerrainGrid.hpp"
#include "algorithms/VoronoiDiagram.hpp"
#include "core/Pole.hpp"
#include "visualization/GnuplotInterface.hpp"
#include "io/BmpHandler.hpp"
#include "algorithms/ComponentCalculator.hpp"

class Control {
private:
    void logOperation(LogLevel level, const std::string& operation, const std::string& details = "") {
        std::string message = std::string("Operation: ") + operation;
        if (!details.empty()) {
            message += std::string(", ") + details;
        }
        logger.logMessage(level, message);
    }

public:
    Config& config;
    Logger& logger;
    
    // Компоненты системы
    Copier copier;
    std::vector<std::array<int, 3>> colors;
    std::vector<std::vector<double>> CopyPole;
    std::vector<Gaus> gaussi;
    std::vector<Component> componenti;
    GaussBuilder gaussBuilder;
    BmpHandler bmpHandler;
    GnuplotInterface gnuplotInterface;
    ComponentCalculator componentCalculator;
    TerrainGrid terrainGrid;
    std::unique_ptr<Pole> p = nullptr;
    std::unique_ptr<KMeans> kMeans = nullptr;
    std::vector<std::vector<double>> kMeansData;
    std::vector<Triangle> lastTriangulation;
    std::vector<VoronoiEdge> voronoiEdges;
    PathFinder pathFinder;
    VoronoiDiagram voronoi;
    std::vector<PointD> clusterCenters;
    std::vector<PointD> path;

    Control(Config& cfg, Logger& log) 
        : config(cfg), 
          logger(log),
          copier(log),
          gaussBuilder(log),
          bmpHandler(log),
          gnuplotInterface(log),
          componentCalculator(log),
          terrainGrid(cfg.fieldWidth, cfg.fieldHeight, log),
          pathFinder(cfg, log),
          voronoi(log) {
        
        kMeans = std::make_unique<KMeans>(log);
        
        logger.info("Control system initialized");
        logger.debug(std::string("Field dimensions: ") + std::to_string(cfg.fieldWidth) + "x" + 
                    std::to_string(cfg.fieldHeight));
    }

    std::vector<PointD> getClusterCenters() const {
        std::vector<PointD> centers;
        for (const auto& component : componenti) {
            if (std::isnan(component.center_x) || std::isnan(component.center_y)) {
                logger.logMessage(LogLevel::Warning, "Skipping invalid cluster center (NaN)");
                continue;
            }
            if (component.center_x >= 0 && component.center_x < config.fieldWidth &&
                component.center_y >= 0 && component.center_y < config.fieldHeight) {
                centers.emplace_back(component.center_x, component.center_y);
            } else {
                logger.logMessage(LogLevel::Error, 
                    std::string("Invalid cluster center: (") + std::to_string(component.center_x) + 
                    ", " + std::to_string(component.center_y) + ")");
            }
        }
        return centers;
    }

    void applyClusterResults(const KMeans::ClusterResult& result, 
                           std::vector<std::vector<double>>& pixelMatrix) {
        for (size_t i = 0; i < result.labels.size(); ++i) {
            int label = result.labels[i];
            if (label >= 0 && label < static_cast<int>(result.colors.size())) {
                int x = static_cast<int>(kMeansData[i][0]);
                int y = static_cast<int>(kMeansData[i][1]);
                const auto& color = result.colors[label];
                pixelMatrix[y][x] = (color[0] + color[1] + color[2]) / 3.0;
            }
        }
    }

    void prepareKMeansData(const std::vector<std::vector<double>>& copyPole) {
        kMeansData.clear();
        size_t validPoints = 0;
        
        for (size_t y = 0; y < copyPole.size(); ++y) {
            for (size_t x = 0; x < copyPole[0].size(); ++x) {
                if (copyPole[y][x] > 0) { 
                    kMeansData.push_back({static_cast<double>(x), static_cast<double>(y)});
                    validPoints++;
                }
            }
        }
            logger.logMessage(LogLevel::Debug, 
                std::string("Prepared ") + std::to_string(validPoints) + " points for K-means");
    }
    
    void Dispetcher(DispatcherParams& params) {     
            logger.logMessage(LogLevel::Info, std::string("Processing command: ") + params.s);
               
        if (params.s == "init") {
            gaussBuilder.init(params.A, params.B, p);
            logOperation(LogLevel::Info, std::string("init"), std::string("size: ") + std::to_string(params.A) + "x" + std::to_string(params.B));
        }

        if (params.s == "g") {
            gaussBuilder.addgauss(params.h, params.x, params.y, params.sx, params.sy, gaussi);
            logOperation(LogLevel::Info, std::string("addgauss"), 
                std::string("x=") + std::to_string(params.x) + 
                ", y=" + std::to_string(params.y) + 
                ", h=" + std::to_string(params.h));
        }
        
        if (params.s == "generate") {
            gaussBuilder.generate(p, gaussi);
            logOperation(LogLevel::Info, std::string("generate"));
        }
        
        if (params.s == "gnuplot") {
            gnuplotInterface.gnuplot(p, params.filename);
            logOperation(LogLevel::Info, std::string("gnuplot"), std::string("file: ") + params.filename);
        }

        if (params.s == "PlotMetedata") {
            gnuplotInterface.plotBinaryWithComponents(CopyPole, componenti, params.filename);
            logOperation(LogLevel::Info, std::string("PlotMetedata"), std::string("file: ") + params.filename);
        }
        
        if (params.s == "PlotVoronoi") {
            gnuplotInterface.plotVoronoi(p, voronoiEdges, clusterCenters, params.filename);
            logOperation(LogLevel::Info, std::string("PlotVoronoi"), std::string("file: ") + params.filename);
        }
        
        if (params.s == "PlotDelaunay") {
            gnuplotInterface.plotDelaunay(lastTriangulation, p, params.filename);
            logOperation(LogLevel::Info, std::string("PlotDelaunay"), std::string("file: ") + params.filename);
        }
        
        if (params.s == "PlotPath") {
            gnuplotInterface.plotPath(path, p, params.filename, params, pathFinder, config.vehicleRadius);
            logOperation(LogLevel::Info, std::string("PlotPath"), std::string("file: ") + params.filename);
        }

        if (params.s == "bmp_write") {
            if (params.bmp_mode == BmpWriteMode::Full) {
                bmpHandler.bmp_write(p->field, params.filename);
            } else {
                bmpHandler.bmp_write(CopyPole, params.filename);
            }
            logOperation(LogLevel::Info, std::string("bmp_write"), 
                std::string("mode: ") + std::string(params.bmp_mode == BmpWriteMode::Full ? "full" : "binary") +
                ", file: " + params.filename);
        }

        if (params.s == "bmp_read") {
            bmpHandler.bmp_read(gaussBuilder, params.filename, CopyPole, p);
            logOperation(LogLevel::Info, std::string("bmp_read"), std::string("file: ") + params.filename);
        }

        if (params.s == "bin") {
            componentCalculator.bin(CopyPole, params.slice, p, params.bin_mode);
            logOperation(LogLevel::Info, std::string("bin"), 
                std::string("slice=") + std::to_string(params.slice) + 
                ", mode=" + (params.bin_mode == ThresholdMode::Peaks ? "peaks" : 
                            params.bin_mode == ThresholdMode::Valleys ? "valleys" : "all"));
        }
        
        if (params.s == "wave") {
            componentCalculator.wave(params.noisy, componenti, CopyPole, p);
            logOperation(LogLevel::Info, std::string("wave"), 
                std::string("noisy=") + std::to_string(params.noisy) + 
                ", components=" + std::to_string(componenti.size()));
            copier.removeNoise(CopyPole, componenti);
        }
         
        if (params.s == "k_means") {
            if (params.k <= 0 || p == nullptr) {
                logger.logMessage(LogLevel::Error, "Invalid parameters for k_means");
                return;
            }
            
            copier.removeNoise(CopyPole, componenti);
            prepareKMeansData(CopyPole);
            if (kMeansData.empty()) {
                logger.logMessage(LogLevel::Warning, "No data for k_means clustering");
                return;
            }

            auto result = kMeans->cluster(kMeansData, params.k);
            applyClusterResults(result, CopyPole);
            logOperation(LogLevel::Info, std::string("k_means"), std::string("k=") + std::to_string(params.k));
        }
        
        if (params.s == "k_means_kern") {
            if (!kMeans || params.k <= 0 || params.kk <= 0 || p == nullptr) {
                logger.logMessage(LogLevel::Error, "Invalid parameters for k_means_kern");
                return;
            }
            
            copier.removeNoise(CopyPole, componenti);
            prepareKMeansData(CopyPole);
            if (kMeansData.empty()) {
                logger.logMessage(LogLevel::Warning, "No data for k_means_kern clustering");
                return;
            }

            auto result = kMeans->kmeansWithKernels(kMeansData, params.k, params.kk);
            applyClusterResults(result, CopyPole);
            logOperation(LogLevel::Info, std::string("k_means_kern"), 
                std::string("k=") + std::to_string(params.k) + 
                ", kk=" + std::to_string(params.kk));
        }
        
        if (params.s == "triangulate") {
            clusterCenters = getClusterCenters();     
            lastTriangulation = componentCalculator.bowyerWatson(clusterCenters);
            voronoi.buildFromDelaunay(lastTriangulation, pathFinder, p, voronoiEdges);
            logOperation(LogLevel::Info, std::string("triangulate"), 
                std::string("clusters=") + std::to_string(clusterCenters.size()) + 
                ", triangles=" + std::to_string(lastTriangulation.size()));
        }

        if (params.s == "find_path") {
            if (!p) {
                logger.logMessage(LogLevel::Error, "Pole not initialized for find_path");
                return;
            }
            
            terrainGrid.calculateSlopes(*p);
            PointD start(params.pointA_x, params.pointA_y);
            PointD goal(params.pointB_x, params.pointB_y);
            
            if (lastTriangulation.empty()) {
                logger.logMessage(LogLevel::Error, "No triangulation for find_path");
                return;
            }
            
            copier.removeNoise(CopyPole, componenti);
            path = pathFinder.findPathAStar(start, goal, lastTriangulation, terrainGrid, CopyPole, params.slice, p);
            
            if (path.empty()) {
                logger.logMessage(LogLevel::Warning, "Path not found");
            } else {
                logOperation(LogLevel::Info, std::string("find_path"), 
                    std::string("from (") + std::to_string(start.x) + "," + std::to_string(start.y) + ")" +
                    " to (" + std::to_string(goal.x) + "," + std::to_string(goal.y) + ")" +
                    ", points=" + std::to_string(path.size()));
            }
        }
        
        if (params.s == "Plot3DPath") {
            PointD start(params.pointA_x, params.pointA_y);
            PointD end(params.pointB_x, params.pointB_y);
            gnuplotInterface.plot3DPath(path, p, params.filename, start, end, pathFinder, config.vehicleRadius);
            logOperation(LogLevel::Info, std::string("Plot3DPath"), std::string("file: ") + params.filename);
        }
        
        if (params.s == "plotInteractive3DPath") {
            PointD start(params.pointA_x, params.pointA_y);
            PointD end(params.pointB_x, params.pointB_y);
            gnuplotInterface.plotInteractive3DPath(path, p, start, end, pathFinder, config.vehicleRadius);
            logOperation(LogLevel::Info, std::string("plotInteractive3DPath"));
        }
    }
};
