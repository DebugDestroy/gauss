#pragma once
#include <string>

namespace core {

class Config {
public:
    // HELP
    std::string defaultHelp;
    
    // FIELD PARAMETERS 
    int fieldWidth, fieldHeight;
    
    // DEFAULT GAUSSIAN PARAMETERS
    double defaultCenterX, defaultCenterY, defaultSigmaX, defaultSigmaY, defaultHeight;
    
    // DEFAULT G_AUTO PARAMETERS
    double xmin, xmax;                                    // Диапазон X координаты центра
    double ymin, ymax;                                    // Диапазон Y координаты центра
    double sx_min, sx_max;                                // Диапазон σx
    double sy_min, sy_max;                                // Диапазон σy
    double h_min, h_max;                                  // Диапазон высоты
    int count_min, count_max;                             // Диапазон количества генерируемых гауссов
    
    // OUTPUT FILES
    std::string defaultGnuplot, defaultPlotMetedata, defaultPlotVoronoi, defaultPlotDelaunay, defaultPlotPath, defaultPlot3DPath,
    defaultWrite, defaultRead, defaultWriteModeImage,
    save_g;
    
    // BINARY
    std::string defaultBinMode;
    int defaultThreshold;
    
    // WAVE
    int defaultNoisy;
    
    // KMEANS
    int defaultKlaster, defaultKlasterKern;
    
    // PATHS
    double defaultstartPointX, defaultstartPointY, defaultendPointX, defaultendPointY;
    double vehicleRadius, maxSideAngle, maxUpDownAngle;
    
    // LOGGER
    std::string logFileNameInterface, logFileNameControl;
    std::string FiltrationLogLevelInterface, FiltrationLogLevelControl;
    
    // Конструктор
    explicit Config(const std::string& filename);
};

}
