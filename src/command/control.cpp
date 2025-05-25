#include "command/control.hpp"

namespace command {

    // Конструктор
    Control::Control(core::Config& cfg, core::Logger& log)
      :
      // Переменные
      colors(),
      CopyPole(),
      
      // core
      config(cfg),
      logger(log),

      // io
      bmpHandler(log),

      // visualization
      gnuplotInterface(log),
      colorGenerator(),

      // algorithms::gauss
      p(nullptr),
      gaussBuilder(log),
      gaussi(),

      // algorithms::components
      copier(log),
      binarizer(log),
      componentCalculator(log),
      componenti(),
      clusterService(log, kMeans, colorGenerator),
      kMeans(log),

      // algorithms::geometry
      clusterCenters(),
      triangulator(log),
      lastTriangulation(),
      voronoi(log),
      voronoiEdges(),
      start(),
      end(),
      path(),

      // algorithms::path::a_star
      conditions(cfg, log),
      graph(log),
      pathFinder(log)
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
            gaussBuilder.init(params.fieldWidth, params.fieldHeight, p);
            logOperation(core::LogLevel::Info, std::string("init"), std::string("size: ") + std::to_string(params.fieldWidth) + "x" + std::to_string(params.fieldHeight));
        }

        if (params.command == "g") {
            gaussBuilder.addgauss(params.height, params.centerX, params.centerY, params.sigmaX, params.sigmaY, gaussi);
            logOperation(core::LogLevel::Info, std::string("addgauss"), 
                std::string("x=") + std::to_string(params.centerX) + 
                ", y=" + std::to_string(params.centerY) + 
                ", h=" + std::to_string(params.height));
        }
        
        if (params.command == "generate") {
            gaussBuilder.generate(p, gaussi);
            logOperation(core::LogLevel::Info, std::string("generate"));
        }
        
        if (params.command == "gnuplot") {
            gnuplotInterface.gnuplot(p, params.filename);
            logOperation(core::LogLevel::Info, std::string("gnuplot"), std::string("file: ") + params.filename);
        }

        if (params.command == "PlotMetedata") {
            gnuplotInterface.plotBinaryWithComponents(CopyPole, componenti, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotMetedata"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotVoronoi") {
            gnuplotInterface.plotVoronoi(p, voronoiEdges, clusterCenters, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotVoronoi"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotDelaunay") {
            gnuplotInterface.plotDelaunay(lastTriangulation, p, params.filename);
            logOperation(core::LogLevel::Info, std::string("PlotDelaunay"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "PlotPath") {
            gnuplotInterface.plotPath(path, p, params.filename, params, config.vehicleRadius);
            logOperation(core::LogLevel::Info, std::string("PlotPath"), std::string("file: ") + params.filename);
        }

        if (params.command == "bmp_write") {
            if (params.bmpWriteMode == io::BmpWriteMode::Full) {
                bmpHandler.bmp_write(p->field, params.filename);
            } else {
                bmpHandler.bmp_write(CopyPole, params.filename);
            }
            logOperation(core::LogLevel::Info, std::string("bmp_write"), 
                std::string("mode: ") + std::string(params.bmpWriteMode == io::BmpWriteMode::Full ? "full" : "binary") +
                ", file: " + params.filename);
        }

        if (params.command == "bmp_read") {
            bmpHandler.bmp_read(gaussBuilder, params.filename, CopyPole, p);
            logOperation(core::LogLevel::Info, std::string("bmp_read"), std::string("file: ") + params.filename);
        }

        if (params.command == "bin") {
            binarizer.bin(CopyPole, params.threshold, p, params.thresholdMode);
            logOperation(core::LogLevel::Info, std::string("bin"), 
                std::string("threshold=") + std::to_string(params.threshold) + 
                ", mode=" + (params.thresholdMode == algorithms::components::ThresholdMode::Peaks ? "peaks" : 
                            params.thresholdMode == algorithms::components::ThresholdMode::Valleys ? "valleys" : "all"));
        }
        
        if (params.command == "wave") {
            componentCalculator.wave(params.noiseLevel, componenti, CopyPole, p);
            logOperation(core::LogLevel::Info, std::string("wave"), 
                std::string("noiseLevel=") + std::to_string(params.noiseLevel) + 
                ", components=" + std::to_string(componenti.size()));
            copier.removeNoise(CopyPole, componenti);
        }
         
        if (params.command == "k_means") {
            if (params.clusterCount <= 0 || p == nullptr) {
                logger.logMessage(core::LogLevel::Error, "Invalid parameters for k_means");
                return;
            }
            
            copier.removeNoise(CopyPole, componenti);
            clusterService.prepareKMeansData(CopyPole);
            if (clusterService.getKMeansData().empty()) {
                logger.logMessage(core::LogLevel::Warning, "No data for k_means clustering");
                return;
            }
            auto result = kMeans.cluster(clusterService.getKMeansData(), params.clusterCount);
            clusterService.applyClusterResults(result, CopyPole);
            logOperation(core::LogLevel::Info, std::string("k_means"), std::string("clusterCount=") + std::to_string(params.clusterCount));
        }
        
       if (params.command == "k_means_kern") {
           if (params.clusterCount <= 0 || params.kernelSize <= 0 || p == nullptr) {
               logger.logMessage(core::LogLevel::Error, "Invalid parameters for k_means_kern");
               return;
           }
    
             copier.removeNoise(CopyPole, componenti);
             clusterService.prepareKMeansData(CopyPole);
             if (clusterService.getKMeansData().empty()) {
                 logger.logMessage(core::LogLevel::Warning, "No data for k_means_kern clustering");
                 return;
             }

    auto result = kMeans.kmeansWithKernels(
        clusterService.getKMeansData(), 
        params.clusterCount, 
        params.kernelSize
    );
    
            clusterService.applyClusterResults(result, CopyPole);
            logOperation(core::LogLevel::Info, std::string("k_means_kern"), 
            std::string("clusterCount=") + std::to_string(params.clusterCount) + 
            ", kernelSize=" + std::to_string(params.kernelSize));
        }
        
        if (params.command == "triangulate") {
            clusterCenters = clusterService.getClusterCenters(componenti, config);  
            lastTriangulation = triangulator.bowyerWatson(clusterCenters, params.fieldWidth, params.fieldHeight);
            voronoi.buildFromDelaunay(lastTriangulation, p, voronoiEdges);
            logOperation(core::LogLevel::Info, std::string("triangulate"), 
                std::string("clusters=") + std::to_string(clusterCenters.size()) + 
                ", triangles=" + std::to_string(lastTriangulation.size()));
        }

        if (params.command == "find_path") {
            if (!p) {
                logger.logMessage(core::LogLevel::Error, "Pole not initialized for find_path");
                return;
            }
            
            if (lastTriangulation.empty()) {
                logger.logMessage(core::LogLevel::Error, "No triangulation for find_path");
                return;
            }
            start = algorithms::geometry::PointD(params.startPointX, params.startPointY);
              end = algorithms::geometry::PointD(params.endPointX, params.endPointY);

            
            copier.removeNoise(CopyPole, componenti);
            path = pathFinder.findPathAStar(start, end, voronoiEdges, graph, conditions, CopyPole, p);
            
            if (path.empty()) {
                logger.logMessage(core::LogLevel::Warning, "Path not found");
            } else {
                logOperation(core::LogLevel::Info, std::string("find_path"), 
                    std::string("from (") + std::to_string(start.x) + "," + std::to_string(start.y) + ")" +
                    " to (" + std::to_string(end.x) + "," + std::to_string(end.y) + ")" +
                    ", points=" + std::to_string(path.size()));
            }
        }
        
        if (params.command == "Plot3DPath") {
            gnuplotInterface.plot3DPath(path, p, params.filename, start, end, config.vehicleRadius);
            logOperation(core::LogLevel::Info, std::string("Plot3DPath"), std::string("file: ") + params.filename);
        }
        
        if (params.command == "plotInteractive3DPath") {
            gnuplotInterface.plotInteractive3DPath(path, p, start, end, config.vehicleRadius);
            logOperation(core::LogLevel::Info, std::string("plotInteractive3DPath"));
        }
    }
}
