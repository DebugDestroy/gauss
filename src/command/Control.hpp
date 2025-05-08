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
#include "services/ClusterService.hpp"
#include "services/ColorGenerator.hpp"
#include "algorithms/Copier.hpp"
#include "algorithms/Component.hpp"
#include "algorithms/GaussBuilder.hpp"
#include "algorithms/KMeans.hpp"
#include "algorithms/PathFinder.hpp"
#include "algorithms/VoronoiDiagram.hpp"
#include "services/Pole.hpp"
#include "visualization/GnuplotInterface.hpp"
#include "io/BmpHandler.hpp"
#include "algorithms/ComponentCalculator.hpp"
#include "algorithms/Triangulator.hpp"

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
    KMeans kMeans;
    ColorGenerator colorGenerator;
    Triangulator triangulator;
    std::unique_ptr<Pole> p = nullptr;
    ClusterService clusterService;
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
          kMeans(log),           // Инициализируем KMeans
          colorGenerator(),       // Инициализируем ColorGenerator
          triangulator(log),
          clusterService(log, kMeans, colorGenerator),
          pathFinder(cfg, log),
          voronoi(log) {
        
        logger.info("Control system initialized");
        logger.debug(std::string("Field dimensions: ") + std::to_string(cfg.fieldWidth) + "x" + 
                    std::to_string(cfg.fieldHeight));
    }
    
    void Dispetcher(DispatcherParams& params) {     
            logger.logMessage(LogLevel::Info, std::string("Processing command: ") + params.command);
               
        if (params.command == "init") {
            gaussBuilder.init(params.fieldWidth, params.fieldHeight, p);
            logOperation(LogLevel::Info, std::string("init"), std::string("size: ") + std::to_string(params.fieldWidth) + "x" + std::to_string(params.fieldHeight));
        }

        if (params.command == "g") {
            gaussBuilder.addgauss(params.height, params.centerX, params.centerY, params.sigmaX, params.sigmaY, gaussi);
            logOperation(LogLevel::Info, std::string("addgauss"), 
                std::string("x=") + std::to_string(params.centerX) + 
                ", y=" + std::to_string(params.centerY) + 
                ", h=" + std::to_string(params.height));
        }
        
        if (params.command == "generate") {
            gaussBuilder.generate(p, gaussi);
            logOperation(LogLevel::Info, std::string("generate"));
        }
        
        if (params.command == "gnuplot") {
            gnuplotInterface.gnuplot(p, params.filename);
            logOperation(LogLevel::Info, std::string("gnuplot"), std::string("file: ") + params.filename);
        }

        if (params.command == "PlotMetedata") {
            gnuplotInterface.plotBinaryWithComponents(CopyPole, componenti, params.filename);
            logOperation(LogLevel::Info, std::string("PlotMetedata"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotVoronoi") {
            gnuplotInterface.plotVoronoi(p, voronoiEdges, clusterCenters, params.filename);
            logOperation(LogLevel::Info, std::string("PlotVoronoi"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotDelaunay") {
            gnuplotInterface.plotDelaunay(lastTriangulation, p, params.filename);
            logOperation(LogLevel::Info, std::string("PlotDelaunay"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotPath") {
            gnuplotInterface.plotPath(path, p, params.filename, params, pathFinder, config.vehicleRadius);
            logOperation(LogLevel::Info, std::string("PlotPath"), std::string("file: ") + params.filename);
        }

        if (params.command == "bmp_write") {
            if (params.bmpWriteMode == BmpWriteMode::Full) {
                bmpHandler.bmp_write(p->field, params.filename);
            } else {
                bmpHandler.bmp_write(CopyPole, params.filename);
            }
            logOperation(LogLevel::Info, std::string("bmp_write"), 
                std::string("mode: ") + std::string(params.bmpWriteMode == BmpWriteMode::Full ? "full" : "binary") +
                ", file: " + params.filename);
        }

        if (params.command == "bmp_read") {
            bmpHandler.bmp_read(gaussBuilder, params.filename, CopyPole, p);
            logOperation(LogLevel::Info, std::string("bmp_read"), std::string("file: ") + params.filename);
        }

        if (params.command == "bin") {
            componentCalculator.bin(CopyPole, params.threshold, p, params.thresholdMode);
            logOperation(LogLevel::Info, std::string("bin"), 
                std::string("threshold=") + std::to_string(params.threshold) + 
                ", mode=" + (params.thresholdMode == ThresholdMode::Peaks ? "peaks" : 
                            params.thresholdMode == ThresholdMode::Valleys ? "valleys" : "all"));
        }
        
        if (params.command == "wave") {
            componentCalculator.wave(params.noiseLevel, componenti, CopyPole, p);
            logOperation(LogLevel::Info, std::string("wave"), 
                std::string("noiseLevel=") + std::to_string(params.noiseLevel) + 
                ", components=" + std::to_string(componenti.size()));
            copier.removeNoise(CopyPole, componenti);
        }
         
        if (params.command == "k_means") {
            if (params.clusterCount <= 0 || p == nullptr) {
                logger.logMessage(LogLevel::Error, "Invalid parameters for k_means");
                return;
            }
            
            copier.removeNoise(CopyPole, componenti);
            clusterService.prepareKMeansData(CopyPole);
            if (clusterService.getKMeansData().empty()) {
                logger.logMessage(LogLevel::Warning, "No data for k_means clustering");
                return;
            }
            auto result = kMeans.cluster(clusterService.getKMeansData(), params.clusterCount);
            clusterService.applyClusterResults(result, CopyPole);
            logOperation(LogLevel::Info, std::string("k_means"), std::string("clusterCount=") + std::to_string(params.clusterCount));
        }
        
       if (params.command == "k_means_kern") {
           if (params.clusterCount <= 0 || params.kernelSize <= 0 || p == nullptr) {
               logger.logMessage(LogLevel::Error, "Invalid parameters for k_means_kern");
               return;
           }
    
             copier.removeNoise(CopyPole, componenti);
             clusterService.prepareKMeansData(CopyPole);
             if (clusterService.getKMeansData().empty()) {
                 logger.logMessage(LogLevel::Warning, "No data for k_means_kern clustering");
                 return;
             }

    auto result = kMeans.kmeansWithKernels(
        clusterService.getKMeansData(), 
        params.clusterCount, 
        params.kernelSize
    );
    
            clusterService.applyClusterResults(result, CopyPole);
            logOperation(LogLevel::Info, std::string("k_means_kern"), 
            std::string("clusterCount=") + std::to_string(params.clusterCount) + 
            ", kernelSize=" + std::to_string(params.kernelSize));
        }
        
        if (params.command == "triangulate") {
            clusterCenters = clusterService.getClusterCenters(componenti, config);  
            lastTriangulation = triangulator.bowyerWatson(clusterCenters);
            voronoi.buildFromDelaunay(lastTriangulation, pathFinder, p, voronoiEdges);
            logOperation(LogLevel::Info, std::string("triangulate"), 
                std::string("clusters=") + std::to_string(clusterCenters.size()) + 
                ", triangles=" + std::to_string(lastTriangulation.size()));
        }

        if (params.command == "find_path") {
            if (!p) {
                logger.logMessage(LogLevel::Error, "Pole not initialized for find_path");
                return;
            }
            
            PointD start(params.startPointX, params.startPointY);
            PointD goal(params.endPointX, params.endPointY);
            
            if (lastTriangulation.empty()) {
                logger.logMessage(LogLevel::Error, "No triangulation for find_path");
                return;
            }
            
            copier.removeNoise(CopyPole, componenti);
            path = pathFinder.findPathAStar(start, goal, lastTriangulation, CopyPole, params.threshold, p);
            
            if (path.empty()) {
                logger.logMessage(LogLevel::Warning, "Path not found");
            } else {
                logOperation(LogLevel::Info, std::string("find_path"), 
                    std::string("from (") + std::to_string(start.x) + "," + std::to_string(start.y) + ")" +
                    " to (" + std::to_string(goal.x) + "," + std::to_string(goal.y) + ")" +
                    ", points=" + std::to_string(path.size()));
            }
        }
        
        if (params.command == "Plot3DPath") {
            PointD start(params.startPointX, params.startPointY);
            PointD end(params.endPointX, params.endPointY);
            gnuplotInterface.plot3DPath(path, p, params.filename, start, end, pathFinder, config.vehicleRadius);
            logOperation(LogLevel::Info, std::string("Plot3DPath"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "plotInteractive3DPath") {
            PointD start(params.startPointX, params.startPointY);
            PointD end(params.endPointX, params.endPointY);
            gnuplotInterface.plotInteractive3DPath(path, p, start, end, pathFinder, config.vehicleRadius);
            logOperation(LogLevel::Info, std::string("plotInteractive3DPath"));
        }
    }
};
