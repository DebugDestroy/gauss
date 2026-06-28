#include "command/control.hpp"

namespace command {

    // Конструктор
    Control::Control(core::Config& cfg, core::Logger& log)
      :
      // core
      config(cfg),
      logger(log),

      // io
      bmpHandler(log),

      // visualization
      gnuplotInterface(log),
      colorGenerator(),

      // algorithms::gauss
      gaussBuilder(log),

      // algorithms::components
      copier(log),
      binarizer(log),
      componentCalculator(log),
      clusterService(log, kMeans, colorGenerator),
      kMeans(log),

      // algorithms::geometry
      triangulator(log),
      voronoi(log),

      // algorithms::path::common
      conditions(log),
      graph(log),
      
      // algorithms::path::a_star
      astarFinder(log),
      
      // algorithms::path::dekstra
      dekstraFinder(log),
      
      // algorithms::path::greedy
      greedyFinder(log)    
{
    logger.info("Control system initialized");
    logger.debug(std::string("Field dimensions: ") +
                 std::to_string(cfg.fieldWidth) + "x" +
                 std::to_string(cfg.fieldHeight));
}
    void Control::logOperation(core::LogLevel level, const std::string& operation, const std::string& details) {
        std::string message = std::string("Operation: ") + operation;
        if (!details.empty()) {
            message += std::string(", ") + details;
        }
        logger.logMessage(level, message);
    }
    
    void Control::Dispetcher(command::DispatcherParams& params) {
            logger.logMessage(core::LogLevel::Info, std::string("Processing command: ") + params.command);
               
        if (params.command == "init") {
            gaussBuilder.init(params.fieldWidth, params.fieldHeight, state.field);
            logOperation(core::LogLevel::Info, std::string("init"), std::string("size: ") + std::to_string(params.fieldWidth) + "x" + std::to_string(params.fieldHeight));
        }

        if (params.command == "g") {
            gaussBuilder.addgauss(params.height, params.centerX, params.centerY, params.sigmaX, params.sigmaY, state.gaussi);
            logOperation(core::LogLevel::Info, std::string("addgauss"), 
                std::string("x=") + std::to_string(params.centerX) + 
                ", y=" + std::to_string(params.centerY) + 
                ", h=" + std::to_string(params.height));
        }
        
        if (params.command == "g_auto") {
            gaussBuilder.addgaussRandom(
            params.xmin, params.xmax, 
            params.ymin, params.ymax,
            params.sx_min, params.sx_max,
            params.sy_min, params.sy_max,
            params.h_min, params.h_max,
            params.count_min, params.count_max,
            params.gAutoMode, params.seedGAuto,
            state.gaussi);
                
            logOperation(core::LogLevel::Info, std::string("g_auto"), std::string(
               "x=[" + std::to_string(params.xmin) + ", " + std::to_string(params.xmax) + "], " +
               "y=[" + std::to_string(params.ymin) + ", " + std::to_string(params.ymax) + "], " +
               "sx=[" + std::to_string(params.sx_min) + ", " + std::to_string(params.sx_max) + "], " +
               "sy=[" + std::to_string(params.sy_min) + ", " + std::to_string(params.sy_max) + "], " +
               "h=[" + std::to_string(params.h_min) + ", " + std::to_string(params.h_max) + "], " +
               "count=[" + std::to_string(params.count_min) + ", " + std::to_string(params.count_max) + "]"));
        }
        
        if (params.command == "generate") {
            gaussBuilder.generate(state.field, state.gaussi);
            logOperation(core::LogLevel::Info, std::string("generate"));
        }
        
        if (params.command == "save_g") {
            gaussBuilder.saveGaussiansToFile(params.filename, state.gaussi);
            logOperation(core::LogLevel::Info, std::string("save_g"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "gnuplot") {
            gnuplotInterface.gnuplot(state.field, params.filename);
            logOperation(core::LogLevel::Info, std::string("gnuplot"), std::string("file: ") + params.filename);
        }

        if (params.command == "PlotMetedata") {
            gnuplotInterface.plotBinaryWithComponents(state.binaryMap, state.components, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotMetedata"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotKmeans") {
            gnuplotInterface.plotKmeans(state.kmeansField, state.kmeansCenters, params.filename);
            logOperation(core::LogLevel::Info, std::string("gnuplot"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotGraph") {
            gnuplotInterface.plotGraph(state.binaryMap, 
                                       state.navigationGraph, 
                                       params.filename, 
                                       state.start,
                                       state.end);
            logOperation(core::LogLevel::Info, std::string("PlotGraph"), std::string("file: ") + params.filename);
        }        
        
        if (params.command == "PlotVoronoi") {
            gnuplotInterface.plotVoronoi(state.binaryMap, state.voronoiEdges, state.clusterCenters, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotVoronoi"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotDelaunay") {
            gnuplotInterface.plotDelaunay(state.lastTriangulation, state.binaryMap, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotDelaunay"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotPath") {
            gnuplotInterface.plotPath(state.path, state.field, state.binaryMap, params.filename, params, params.vehicleRadius);
            logOperation(core::LogLevel::Info, std::string("PlotPath"), std::string("file: ") + params.filename);
        }

        if (params.command == "bmp_write") {
            if (params.bmpWriteMode == io::BmpWriteMode::Full) {
                bmpHandler.bmp_write(state.field, params.filename);
            } else {
                bmpHandler.bmp_write(state.binaryMap, params.filename);
            }
            logOperation(core::LogLevel::Info, std::string("bmp_write"), 
                std::string("mode: ") + std::string(params.bmpWriteMode == io::BmpWriteMode::Full ? "full" : "binary") +
                ", file: " + params.filename);
        }

        if (params.command == "bmp_read") {
            bmpHandler.bmp_read(gaussBuilder, params.filename, state.field);
            logOperation(core::LogLevel::Info, std::string("bmp_read"), std::string("file: ") + params.filename);
        }

        if (params.command == "bin") {
            binarizer.bin(state.binaryMap, params.threshold, state.field, params.thresholdMode);
            logOperation(core::LogLevel::Info, std::string("bin"), 
                std::string("threshold=") + std::to_string(params.threshold) + 
                ", mode=" + (params.thresholdMode == algorithms::components::ThresholdMode::Peaks ? "peaks" : 
                            params.thresholdMode == algorithms::components::ThresholdMode::Valleys ? "valleys" : "all"));
        }
        
        if (params.command == "wave") {
            componentCalculator.wave(params.noiseLevel, state.components, state.binaryMap, state.field);
            logOperation(core::LogLevel::Info, std::string("wave"), 
                std::string("noiseLevel=") + std::to_string(params.noiseLevel) + 
                ", components=" + std::to_string(state.components.size()));
        }
         
        if (params.command == "k_means") {
            if (params.clusterCount <= 0 || state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Invalid parameters for k_means");
                return;
            }
            
            state.kmeansField = state.binaryMap;
            clusterService.prepareKMeansData(state.kmeansField);
            if (clusterService.getKMeansData().empty()) {
                logger.logMessage(core::LogLevel::Warning, "No data for k_means clustering");
                return;
            }
            auto result = kMeans.cluster(clusterService.getKMeansData(), params.clusterCount);
            state.kmeansCenters.clear();
            for (const auto& c : result.centers)
                {
                    state.kmeansCenters.emplace_back(c[0], c[1]);
                }
            clusterService.applyClusterResults(result, state.kmeansField);
            logOperation(core::LogLevel::Info, std::string("k_means"), std::string("clusterCount=") + std::to_string(params.clusterCount));
        }
        
       if (params.command == "k_means_kern") {
           if (params.clusterCount <= 0 || params.kernelSize <= 0 || state.field.empty()) {
               logger.logMessage(core::LogLevel::Error, "Invalid parameters for k_means_kern");
               return;
           }
    
             state.kmeansField = state.binaryMap;
             clusterService.prepareKMeansData(state.kmeansField);
             if (clusterService.getKMeansData().empty()) {
                 logger.logMessage(core::LogLevel::Warning, "No data for k_means_kern clustering");
                 return;
             }

    auto result = kMeans.kmeansWithKernels(
        clusterService.getKMeansData(), 
        params.clusterCount, 
        params.kernelSize
    );
        state.kmeansCenters.clear();
            for (const auto& c : result.centers)
                {
                    state.kmeansCenters.emplace_back(c[0], c[1]);
                }
            clusterService.applyClusterResults(result, state.kmeansField);
            logOperation(core::LogLevel::Info, std::string("k_means_kern"), 
            std::string("clusterCount=") + std::to_string(params.clusterCount) + 
            ", kernelSize=" + std::to_string(params.kernelSize));
        }
        
        if (params.command == "triangulate") {
            state.clusterCenters = clusterService.getClusterCenters(state.components, params.fieldWidth, params.fieldHeight);  
            state.lastTriangulation = triangulator.bowyerWatson(state.clusterCenters, params.fieldWidth, params.fieldHeight);
            logOperation(core::LogLevel::Info, std::string("triangulate"), 
                std::string("clusters=") + std::to_string(state.clusterCenters.size()) + 
                ", triangles=" + std::to_string(state.lastTriangulation.size()));
        }
        
        if (params.command == "voronoi") {
            voronoi.buildFromDelaunay(state.lastTriangulation, state.field, state.voronoiEdges);
            logOperation(core::LogLevel::Info, std::string("voronoi"));
        }
        if (params.command == "build_nav_graph") {
            state.navigationGraph = 
                graph.buildGraphFromEdges(
                    state.voronoiEdges, 
                    state.binaryMap, 
                    state.field, 
                    conditions, 
                    params.vehicleRadius,
                    params.maxSideAngle,
                    params.maxUpDownAngle);                                                                    
            logOperation(core::LogLevel::Info, std::string("build_nav_graph"));
        }
        if (params.command == "connect_to_graph") {
            state.start = algorithms::geometry::Pixel(params.startPointX, params.startPointY);
            state.end = algorithms::geometry::Pixel(params.endPointX, params.endPointY);
                
                 graph.connectPointToGraph(
                     state.navigationGraph,
                     *state.start,
                     state.binaryMap,
                     state.field,
                     conditions,
                     params.vehicleRadius,
                     params.maxSideAngle,
                     params.maxUpDownAngle,
                     params.connectMode,
                     params.nearestVerticesCount);

               graph.connectPointToGraph(
                     state.navigationGraph,
                     *state.end,
                     state.binaryMap,
                     state.field,
                     conditions,
                     params.vehicleRadius,
                     params.maxSideAngle,
                     params.maxUpDownAngle,
                     params.connectMode,
                     params.nearestVerticesCount);

           logOperation(core::LogLevel::Info, "connect_to_graph");
        }
        if (params.command == "find_path_astar") {
            if (state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for A*");
                return;
            }
            
            state.path = astarFinder.findPathAStar(*state.start, *state.end, state.navigationGraph);
            
            if (state.path.empty()) {
                logger.logMessage(core::LogLevel::Warning, "Path A* not found");
            } else {
                logOperation(core::LogLevel::Info, std::string("find_path_astar"), 
                    std::string("from (") + std::to_string(state.start->x) + "," + std::to_string(state.start->y) + ")" +
                    " to (" + std::to_string(state.end->x) + "," + std::to_string(state.end->y) + ")" +
                    ", points=" + std::to_string(state.path.size()));
            }
        }

        if (params.command == "find_path_dekstra") {
            if (state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for dekstra");
                return;
            }

            state.path = dekstraFinder.findPathDijkstra(*state.start, *state.end, state.navigationGraph);
            
            if (state.path.empty()) {
                logger.logMessage(core::LogLevel::Warning, "Path dekstra not found");
            } else {
                logOperation(core::LogLevel::Info, std::string("find_path_dekstra"), 
                    std::string("from (") + std::to_string(state.start->x) + "," + std::to_string(state.start->y) + ")" +
                    " to (" + std::to_string(state.end->x) + "," + std::to_string(state.end->y) + ")" +
                    ", points=" + std::to_string(state.path.size()));
            }
        }

        if (params.command == "find_path_greedy") {
            if (state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for greedy");
                return;
            }

            state.path = greedyFinder.findPathGreedy(*state.start, *state.end, state.navigationGraph);
            
            if (state.path.empty()) {
                logger.logMessage(core::LogLevel::Warning, "Path greedy not found");
            } else {
                logOperation(core::LogLevel::Info, std::string("find_path_greedy"), 
                    std::string("from (") + std::to_string(state.start->x) + "," + std::to_string(state.start->y) + ")" +
                    " to (" + std::to_string(state.end->x) + "," + std::to_string(state.end->y) + ")" +
                    ", points=" + std::to_string(state.path.size()));
            }
        }
        
        if (params.command == "Plot3DPath") {
            gnuplotInterface.plot3DPath(state.path, state.field, params.filename, *state.start, *state.end, params.vehicleRadius);
            logOperation(core::LogLevel::Info, std::string("Plot3DPath"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "plotInteractive3DPath") {
            gnuplotInterface.plotInteractive3DPath(state.path, state.field, *state.start, *state.end, params.vehicleRadius);
            logOperation(core::LogLevel::Info, std::string("plotInteractive3DPath"));
        }
    }
}
