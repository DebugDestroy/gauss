#include "command/control.hpp"

namespace command {

    // Конструктор
    Control::Control(core::Logger& log, const std::string seedMode, uint32_t seed)
      :
      // core
      logger(log),

      // io
      bmpHandler(log),

      // visualization
      gnuplotInterface(log),
      colorGenerator(log, randomGenerator),
      
      // statistics
      statisticsManager(),
      csvWriter(log),
      
      // algorithms::gauss
      gaussBuilder(log, randomGenerator),

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
      grid(log),
      
      // algorithms::path::a_star
      astarFinder(log),
      
      // algorithms::path::dekstra
      dekstraFinder(log),
      
      // algorithms::path::greedy
      greedyFinder(log),
      
      // algorithms::path::rrt
      rrt(log, randomGenerator),
      
      // algorithms::path::rrt_star
      rrtStar(log, randomGenerator)
{
    logger.info("Control system initialized");
    
    if (seedMode == "Random")
    {
        std::random_device rd;
        uint32_t generatedSeed = rd();

        randomGenerator.seed(generatedSeed);

        logger.info("Random generator seed = " + std::to_string(generatedSeed));
    }
    else
    {
        randomGenerator.seed(seed);

        logger.info("Random generator seed = " + std::to_string(seed));
    }
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
                                       state.startPixel,
                                       state.goalPixel);
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
        
        if (params.command == "PlotGrid") {
            gnuplotInterface.plotGrid(state.grid, state.binaryMap, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotGrid"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotNavGrid") {
            gnuplotInterface.plotNavGrid(state.grid, state.binaryMap, params.filename, state.startCell, state.goalCell);
            logOperation(core::LogLevel::Info, std::string("PlotNavGrid"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotGridPath") {
            gnuplotInterface.plotGridPath(state.grid, state.gridPath, *state.startCell, *state.goalCell, state.binaryMap, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotGridPath"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotPathDiscrete") {
            gnuplotInterface.plotPathDiscrete(state.pathPixel, state.field, state.binaryMap, params.filename, params, params.vehicleRadiusPixel);
            logOperation(core::LogLevel::Info, std::string("PlotPathDiscrete"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotPathContinuous") {
            gnuplotInterface.plotPathContinuous(state.pathWorld,
                                                state.gaussi,
                                                params.fieldWidth, params.fieldHeight,
                                                params.heightThresholdWorld,
                                                params.filename);
                                                
            logOperation(core::LogLevel::Info, std::string("PlotPathContinuous"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotRRT") {
            gnuplotInterface.plotRRT(state.pathWorld, state.treeRRT, state.startWorld, state.goalWorld, state.binaryMap, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotRRT"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotRRTStar") {
            gnuplotInterface.plotRRT(state.pathWorld, state.treeRRTStar, state.startWorld, state.goalWorld, state.binaryMap, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotRRTStar"), std::string("file: ") + params.filename);
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
            binarizer.bin(state.binaryMap, params.heightThresholdPixel, state.field);
            logOperation(core::LogLevel::Info, std::string("bin"), 
                std::string("heightThresholdPixel=") + std::to_string(params.heightThresholdPixel));
        }
        
        if (params.command == "wave") {
            componentCalculator.wave(params.waveNoisy, state.components, state.binaryMap, state.field);
            logOperation(core::LogLevel::Info, std::string("wave"), 
                std::string("waveNoisy=") + std::to_string(params.waveNoisy) + 
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
                    params.vehicleRadiusPixel,
                    params.maxSideAngle,
                    params.maxUpDownAngle);                                                                    
            logOperation(core::LogLevel::Info, std::string("build_nav_graph"));
        }
        
        if (params.command == "grid") {
            grid.buildGrid(state.grid, params.fieldWidth, params.fieldHeight, params.gridWidth);
            logOperation(core::LogLevel::Info, std::string("grid"));
        }
        
        if (params.command == "build_nav_grid") {
                grid.buildNavGrid(
                    state.grid,
                    state.binaryMap,
                    conditions,
                    params.gridNoisy);                                                                    
            logOperation(core::LogLevel::Info, std::string("build_nav_grid"));
        }
        
        if (params.command == "connect_to_grid") {
            state.startPixel = algorithms::geometry::Pixel(params.startPixelX, params.startPixelY);
            state.goalPixel = algorithms::geometry::Pixel(params.goalPixelX, params.goalPixelY);
                
                 state.startCell = grid.connectPointToGrid(state.grid, *state.startPixel);
                 state.goalCell = grid.connectPointToGrid(state.grid, *state.goalPixel);

           logOperation(core::LogLevel::Info, "connect_to_grid");
        }
        
        if (params.command == "connect_to_graph") {
            state.startPixel = algorithms::geometry::Pixel(params.startPixelX, params.startPixelY);
            state.goalPixel = algorithms::geometry::Pixel(params.goalPixelX, params.goalPixelY);
                
                 graph.connectPointToGraph(
                     state.navigationGraph,
                     *state.startPixel,
                     state.binaryMap,
                     state.field,
                     conditions,
                     params.vehicleRadiusPixel,
                     params.maxSideAngle,
                     params.maxUpDownAngle,
                     params.connectMode,
                     params.nearestVerticesCount);

               graph.connectPointToGraph(
                     state.navigationGraph,
                     *state.goalPixel,
                     state.binaryMap,
                     state.field,
                     conditions,
                     params.vehicleRadiusPixel,
                     params.maxSideAngle,
                     params.maxUpDownAngle,
                     params.connectMode,
                     params.nearestVerticesCount);

           logOperation(core::LogLevel::Info, "connect_to_graph");
        }
        
        if (params.command == "astar_graph") {
            if (state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for A*");
                return;
            }
            
            state.PathMetrics.environment = "graph";
            statisticsManager.reset(state.PathMetrics);
            state.PathMetrics.algorithmName = "A*";
            statisticsManager.startTimer();       
            
            state.pathPixel = astarFinder.findPathAStarGraph(*state.startPixel, *state.goalPixel, state.navigationGraph, state.PathMetrics);
            
            statisticsManager.finishTimer(state.PathMetrics);
            
            if (state.PathMetrics.pathFound) {
                statisticsManager.computePathLength(state.PathMetrics, state.pathPixel);
                statisticsManager.computeMaxTerrainAngles(state.PathMetrics, state.pathPixel, state.field, params.vehicleRadiusPixel);
                statisticsManager.computeMinObstacleDistance(state.PathMetrics, state.pathPixel, state.binaryMap);
            
                logOperation(
                    core::LogLevel::Info,
                    "astar_graph",
                    Control::formatPathMetricsLog(state.PathMetrics, *state.startPixel, *state.goalPixel)
                );
            } else {
                logger.logMessage(core::LogLevel::Warning,
                                  "[astar_graph] Path not found");
           }
        }

        if (params.command == "dekstra_graph") {
            if (state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for dekstra");
                return;
            }
            
            state.PathMetrics.environment = "graph";            
            statisticsManager.reset(state.PathMetrics);
            state.PathMetrics.algorithmName = "Dijkstra";
            statisticsManager.startTimer();       
            
            state.pathPixel = dekstraFinder.findPathDijkstraGraph(*state.startPixel, *state.goalPixel, state.navigationGraph, state.PathMetrics);
             
            statisticsManager.finishTimer(state.PathMetrics);
            
            if (state.PathMetrics.pathFound) {
                statisticsManager.computePathLength(state.PathMetrics, state.pathPixel);
                statisticsManager.computeMaxTerrainAngles(state.PathMetrics, state.pathPixel, state.field, params.vehicleRadiusPixel);
                statisticsManager.computeMinObstacleDistance(state.PathMetrics, state.pathPixel, state.binaryMap);
            
                logOperation(
                    core::LogLevel::Info,
                    "dekstra_graph",
                    Control::formatPathMetricsLog(state.PathMetrics, *state.startPixel, *state.goalPixel)
                );
            } else {
                logger.logMessage(core::LogLevel::Warning,
                                  "[dekstra_graph] Path not found");
           }
        }

        if (params.command == "greedy_graph") {
            if (state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for greedy");
                return;
            }
            
            state.PathMetrics.environment = "graph";
            statisticsManager.reset(state.PathMetrics);
            state.PathMetrics.algorithmName = "Greedy";
            statisticsManager.startTimer();       
            
            state.pathPixel = greedyFinder.findPathGreedyGraph(*state.startPixel, *state.goalPixel, state.navigationGraph, state.PathMetrics);
                         
            statisticsManager.finishTimer(state.PathMetrics);
            
            if (state.PathMetrics.pathFound) {
                statisticsManager.computePathLength(state.PathMetrics, state.pathPixel);
                statisticsManager.computeMaxTerrainAngles(state.PathMetrics, state.pathPixel, state.field, params.vehicleRadiusPixel);
                statisticsManager.computeMinObstacleDistance(state.PathMetrics, state.pathPixel, state.binaryMap);
            
                logOperation(
                    core::LogLevel::Info,
                    "greedy_graph",
                    Control::formatPathMetricsLog(state.PathMetrics, *state.startPixel, *state.goalPixel)
                );
            } else {
                logger.logMessage(core::LogLevel::Warning,
                                  "[greedy_graph] Path not found");
           }
        }
        
        if (params.command == "astar_grid") {
            if (state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for A*");
                return;
            }
            
            state.PathMetrics.environment = "grid";
            statisticsManager.reset(state.PathMetrics);
            state.PathMetrics.algorithmName = "A*";
            statisticsManager.startTimer();       
            
            state.gridPath = astarFinder.findPathAStarGrid(*state.startCell, *state.goalCell, state.grid, state.PathMetrics);
            
            statisticsManager.finishTimer(state.PathMetrics);
            
            state.pathPixel = algorithms::geometry::toPixelPath(state.gridPath, params.gridWidth);
             
            if (state.PathMetrics.pathFound) {
                statisticsManager.computePathLength(state.PathMetrics, state.pathPixel);
                statisticsManager.computeMaxTerrainAngles(state.PathMetrics, state.pathPixel, state.field, params.vehicleRadiusPixel);
                statisticsManager.computeMinObstacleDistance(state.PathMetrics, state.pathPixel, state.binaryMap);
            
                logOperation(
                    core::LogLevel::Info,
                    "astar_grid",
                    Control::formatPathMetricsLog(state.PathMetrics, *state.startPixel, *state.goalPixel)
                );
            } else {
                logger.logMessage(core::LogLevel::Warning,
                                  "[astar_grid] Path not found");
           }
        }

        if (params.command == "dekstra_grid") {
            if (state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for dekstra");
                return;
            }
            
            state.PathMetrics.environment = "grid";            
            statisticsManager.reset(state.PathMetrics);
            state.PathMetrics.algorithmName = "Dijkstra";
            statisticsManager.startTimer();       
            
            state.gridPath = dekstraFinder.findPathDijkstraGrid(*state.startCell, *state.goalCell, state.grid, state.PathMetrics);
            
            statisticsManager.finishTimer(state.PathMetrics);
            
            state.pathPixel = algorithms::geometry::toPixelPath(state.gridPath, params.gridWidth);
             
            if (state.PathMetrics.pathFound) {
                statisticsManager.computePathLength(state.PathMetrics, state.pathPixel);
                statisticsManager.computeMaxTerrainAngles(state.PathMetrics, state.pathPixel, state.field, params.vehicleRadiusPixel);
                statisticsManager.computeMinObstacleDistance(state.PathMetrics, state.pathPixel, state.binaryMap);
            
                logOperation(
                    core::LogLevel::Info,
                    "dekstra_grid",
                    Control::formatPathMetricsLog(state.PathMetrics, *state.startPixel, *state.goalPixel)
                );
            } else {
                logger.logMessage(core::LogLevel::Warning,
                                  "[dekstra_grid] Path not found");
           }
        }

        if (params.command == "greedy_grid") {
            if (state.field.empty()) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for greedy");
                return;
            }
            
            state.PathMetrics.environment = "grid";
            statisticsManager.reset(state.PathMetrics);
            state.PathMetrics.algorithmName = "Greedy";
            statisticsManager.startTimer();       
            
            state.gridPath = greedyFinder.findPathGreedyGrid(*state.startCell, *state.goalCell, state.grid, state.PathMetrics);
            
            statisticsManager.finishTimer(state.PathMetrics);
            
            state.pathPixel = algorithms::geometry::toPixelPath(state.gridPath, params.gridWidth);
             
            if (state.PathMetrics.pathFound) {
                statisticsManager.computePathLength(state.PathMetrics, state.pathPixel);
                statisticsManager.computeMaxTerrainAngles(state.PathMetrics, state.pathPixel, state.field, params.vehicleRadiusPixel);
                statisticsManager.computeMinObstacleDistance(state.PathMetrics, state.pathPixel, state.binaryMap);
            
                logOperation(
                    core::LogLevel::Info,
                    "greedy_grid",
                    Control::formatPathMetricsLog(state.PathMetrics, *state.startPixel, *state.goalPixel)
                );
            } else {
                logger.logMessage(core::LogLevel::Warning,
                                  "[greedy_grid] Path not found");
           }
        }
        
        if (params.command == "rrt") {
            if (state.gaussi.empty()) {
                logger.logMessage(core::LogLevel::Error, "Gaussi not initialized for rrt_star");
                return;
            }
            
            state.startWorld = algorithms::geometry::PointD(params.startWorldX, params.startWorldY);
            state.goalWorld = algorithms::geometry::PointD(params.goalWorldX, params.goalWorldY);
            
            state.PathMetrics.environment = "continuous";
            statisticsManager.reset(state.PathMetrics);
            state.PathMetrics.algorithmName = "RRT";
            statisticsManager.startTimer();       
            
            auto result = rrt.findPathRRT(state.startWorld, state.goalWorld, 
                                          state.gaussi,
                                          params.fieldWidth, params.fieldHeight,
                                          params.heightThresholdWorld,
                                          params.vehicleRadiusWorld,
                                          params.maxSideAngle, params.maxUpDownAngle,
                                          params.interpEdge, params.interpCollision, params.interpAngle,
                                          params.maxIterations,
                                          params.step,
                                          params.goalRadius,
                                          params.goalBias,
                                          conditions,
                                          state.PathMetrics);
            
            statisticsManager.finishTimer(state.PathMetrics);
            
            state.pathWorld = std::move(result.path);
            state.treeRRT   = std::move(result.tree);
             
            if (state.PathMetrics.pathFound) {
                statisticsManager.computePathLength(state.PathMetrics, 
                                                    state.pathWorld);
                                                    
                statisticsManager.computeMaxTerrainAngles(state.PathMetrics, 
                                                          state.pathWorld, 
                                                          state.gaussi, 
                                                          params.vehicleRadiusWorld, 
                                                          params.interpEdge);
                                                          
                statisticsManager.computeMinObstacleDistance(state.PathMetrics, 
                                                             state.pathWorld, 
                                                             state.gaussi, 
                                                             params.fieldWidth, 
                                                             params.fieldHeight, 
                                                             params.heightThresholdWorld,
                                                             params.interpEdge,
                                                             params.interpCollision,
                                                             params.interpAngle);
            
                logOperation(
                    core::LogLevel::Info,
                    "rrt",
                    Control::formatPathMetricsLog(state.PathMetrics, state.startWorld, state.goalWorld)
                );
            } else {
                logger.logMessage(core::LogLevel::Warning,
                                  "[rrt] Path not found");
           }
        }
        
        if (params.command == "rrt_star") {
            if (state.gaussi.empty()) {
                logger.logMessage(core::LogLevel::Error, "Gaussi not initialized for rrt_star");
                return;
            }
            
            state.startWorld = algorithms::geometry::PointD(params.startWorldX, params.startWorldY);
            state.goalWorld = algorithms::geometry::PointD(params.goalWorldX, params.goalWorldY);
            
            state.PathMetrics.environment = "continuous";
            statisticsManager.reset(state.PathMetrics);
            state.PathMetrics.algorithmName = "RRT*";
            statisticsManager.startTimer();
            
            auto result = rrtStar.findPathRRTStar(state.startWorld, state.goalWorld, 
                                                  state.gaussi,
                                                  params.fieldWidth, params.fieldHeight,
                                                  params.heightThresholdWorld,
                                                  params.vehicleRadiusWorld,
                                                  params.maxSideAngle, params.maxUpDownAngle,
                                                  params.interpEdge, params.interpCollision, params.interpAngle,
                                                  params.maxIterations,
                                                  params.step,
                                                  params.maxFindRadius,
                                                  params.gammaConstant,
                                                  params.goalRadius,
                                                  params.goalBias,
                                                  conditions,
                                                  state.PathMetrics);
            
            statisticsManager.finishTimer(state.PathMetrics);
            
            state.pathWorld = std::move(result.path);
            state.treeRRTStar   = std::move(result.tree);
             
            if (state.PathMetrics.pathFound) {
                statisticsManager.computePathLength(state.PathMetrics, 
                                                    state.pathWorld);
                                                    
                statisticsManager.computeMaxTerrainAngles(state.PathMetrics, 
                                                          state.pathWorld, 
                                                          state.gaussi, 
                                                          params.vehicleRadiusWorld, 
                                                          params.interpEdge);
                                                          
                statisticsManager.computeMinObstacleDistance(state.PathMetrics, 
                                                             state.pathWorld, 
                                                             state.gaussi, 
                                                             params.fieldWidth, 
                                                             params.fieldHeight, 
                                                             params.heightThresholdWorld,
                                                             params.interpEdge,
                                                             params.interpCollision,
                                                             params.interpAngle);
            
                logOperation(
                    core::LogLevel::Info,
                    "rrt_star",
                    Control::formatPathMetricsLog(state.PathMetrics, state.startWorld, state.goalWorld)
                );
            } else {
                logger.logMessage(core::LogLevel::Warning,
                                  "[rrt_star] Path not found");
           }
        }
        
        if (params.command == "shortcut_discrete") {
            algorithms::path::smoothing::Shortcut::shortcutDiscrete(
                    state.pathPixel,
                    state.field,
                    state.binaryMap,
                    conditions,
                    params.vehicleRadiusPixel,
                    params.maxSideAngle,
                    params.maxUpDownAngle);     
                                                                                   
            logOperation(core::LogLevel::Info, std::string("shortcut_discrete"));
        }
        
        if (params.command == "shortcut_continuous") {
            algorithms::path::smoothing::Shortcut::shortcutContinuous(
                    state.pathWorld,
                    state.gaussi,
                    params.fieldWidth, params.fieldHeight,
                    params.heightThresholdWorld,
                    params.vehicleRadiusWorld,
                    params.maxSideAngle, params.maxUpDownAngle,
                    params.interpEdge, params.interpCollision, params.interpAngle,
                    conditions);       
                                                                                 
            logOperation(core::LogLevel::Info, std::string("shortcut_continuous"));
        }
        
        if (params.command == "spline_discrete") {
            state.pathWorld = algorithms::geometry::toPointDPath(state.pathPixel);
            algorithms::path::smoothing::Spline::spline(
                    state.pathWorld,
                    state.gaussi,
                    params.fieldWidth, params.fieldHeight,
                    params.heightThresholdWorld,
                    params.vehicleRadiusWorld,
                    params.maxSideAngle, params.maxUpDownAngle,
                    params.interpEdge, params.interpCollision, params.interpAngle,
                    conditions,
                    params.samplesPerSegment);                                                                       
            logOperation(core::LogLevel::Info, std::string("spline_discrete"));
        }
        
        if (params.command == "spline_continuous") {
            algorithms::path::smoothing::Spline::spline(
                    state.pathWorld,
                    state.gaussi,
                    params.fieldWidth, params.fieldHeight,
                    params.heightThresholdWorld,
                    params.vehicleRadiusWorld,
                    params.maxSideAngle, params.maxUpDownAngle,
                    params.interpEdge, params.interpCollision, params.interpAngle,
                    conditions,
                    params.samplesPerSegment);                                                                       
            logOperation(core::LogLevel::Info, std::string("spline_continuous"));
        }
        
        if (params.command == "save_metrics") {
            csvWriter.writeMetrics(state.PathMetrics, params.filename);
            logOperation(core::LogLevel::Info, std::string("save_metrics"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "Plot3DPath") {
            gnuplotInterface.plot3DPath(state.pathPixel, state.field, params.filename, *state.startPixel, *state.goalPixel, params.vehicleRadiusPixel);
            logOperation(core::LogLevel::Info, std::string("Plot3DPath"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "plotInteractive3DPath") {
            gnuplotInterface.plotInteractive3DPath(state.pathPixel, state.field, *state.startPixel, *state.goalPixel, params.vehicleRadiusPixel);
            logOperation(core::LogLevel::Info, std::string("plotInteractive3DPath"));
        }
    }
}
