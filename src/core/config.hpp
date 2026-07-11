#pragma once
#include <string>
#include <fstream>
#include <stdexcept>
#include <unordered_set>

namespace core {

class Config {
private:

std::unordered_set<std::string> readKeys;
static const std::unordered_set<std::string> required;

template<typename T>
void readParameter(std::ifstream& file,
                   T& value,
                   const std::string& parameterName)
{
    if (!(file >> value))
    {
        throw std::runtime_error(
            "Invalid value for parameter '" +
            parameterName + "'");
    }
}
public:
    // HELP
    std::string defaultHelp;
    
    // FIELD PARAMETERS 
    int defaultfieldWidth = 0, defaultfieldHeight = 0;
    
    // DEFAULT GAUSSIAN PARAMETERS
    double defaultCenterX = 0, defaultCenterY = 0, defaultSigmaX = 0, defaultSigmaY = 0, defaultHeight = 0;
    
    // DEFAULT G_AUTO PARAMETERS
    std::size_t count_min = 0, count_max = 0;                             // Диапазон количества генерируемых гауссов
    double xmin = 0, xmax = 0;                                    // Диапазон X координаты центра
    double ymin = 0, ymax = 0;                                    // Диапазон Y координаты центра
    double sx_min = 0, sx_max = 0;                                // Диапазон σx
    double sy_min = 0, sy_max = 0;                                // Диапазон σy
    double h_min = 0, h_max = 0;                                  // Диапазон высоты
    
    // OUTPUT FILES
    std::string defaultGnuplot, defaultPlotMetedata, defaultPlotKmeans, defaultPlotGraph, defaultPlotVoronoi, defaultPlotDelaunay, defaultPlotGrid, defaultPlotNavGrid, defaultPlotGridPath, 
    defaultPlotPath, PlotRRT, defaultPlot3DPath;
    
    std::string defaultWrite, defaultRead, defaultWriteModeImage,
    save_g,
    defaultsave_metrics;
    
    // BINARY
    int heightThresholdPixel = 0;
    
    // WAVE
    std::size_t defaultWaveNoisy = 0;
    
    // KMEANS
    std::size_t defaultKlaster = 0, defaultKlasterKern = 0;
    
    // GRID
    int defaultgridWidth = 0;
    std::size_t defaultgridNoisy = 0;
    
    // PATHS
    int startPixelX = 0, startPixelY = 0, goalPixelX = 0, goalPixelY = 0;
    std::string defaultConnectMode;
    std::size_t defaultNearestVerticesCount = 0;
    
    int vehicleRadiusPixel = 0; 
    double maxSideAngle = 0, maxUpDownAngle = 0;
    
    // RRT
    size_t maxIterations = 0;
    double startWorldX = 0;
    double startWorldY = 0;
    double goalWorldX = 0;
    double goalWorldY = 0;
    double vehicleRadiusWorld = 0;
    double heightThresholdWorld = 0;
    double interpEdge = 0;
    double interpCollision = 0;
    double interpAngle = 0;
    double step = 0;
    double goalRadius = 0;
    double goalBias = 0;

    // LOGGER
    std::string logFileNameInterface, logFileNameControl;
    std::string FiltrationLogLevelInterface, FiltrationLogLevelControl;
    
    // SEED
    std::string seedMode;
    std::uint32_t seed = 0;
    
    // Конструктор
    explicit Config(const std::string& filename);
};

}
