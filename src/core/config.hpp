#pragma once
#include <string>
#include <fstream>
#include <stdexcept>

namespace core {

class Config {
private:

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
    double xmin = 0, xmax = 0;                                    // Диапазон X координаты центра
    double ymin = 0, ymax = 0;                                    // Диапазон Y координаты центра
    double sx_min = 0, sx_max = 0;                                // Диапазон σx
    double sy_min = 0, sy_max = 0;                                // Диапазон σy
    double h_min = 0, h_max = 0;                                  // Диапазон высоты
    int count_min = 0, count_max = 0;                             // Диапазон количества генерируемых гауссов
    std::string defaultGAutoMode;
    std::uint32_t defaultSeedGAuto = 0;
    
    // OUTPUT FILES
    std::string defaultGnuplot, defaultPlotMetedata, defaultPlotKmeans, defaultPlotGraph, defaultPlotVoronoi, defaultPlotDelaunay, defaultPlotGrid, defaultPlotNavGrid, defaultPlotGridPath, 
    defaultPlotPath, defaultPlot3DPath;
    
    std::string defaultWrite, defaultRead, defaultWriteModeImage,
    save_g,
    defaultsave_metrics;
    
    // BINARY
    std::string defaultBinMode;
    int defaultThreshold = 0;
    
    // WAVE
    int defaultWaveNoisy = 0;
    
    // KMEANS
    int defaultKlaster = 0, defaultKlasterKern = 0;
    
    // GRID
    int defaultgridWidth = 0;
    int defaultgridNoisy = 0;
    
    // PATHS
    int defaultstartPointX = 0, defaultstartPointY = 0, defaultendPointX = 0, defaultendPointY = 0;
    std::string defaultConnectMode;
    int defaultNearestVerticesCount = 0;
    
    int vehicleRadius = 0; 
    double maxSideAngle = 0, maxUpDownAngle = 0;
    
    // LOGGER
    std::string logFileNameInterface, logFileNameControl;
    std::string FiltrationLogLevelInterface, FiltrationLogLevelControl;
    
    // Конструктор
    explicit Config(const std::string& filename);
};

}
