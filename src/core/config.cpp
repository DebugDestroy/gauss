#include "core/config.hpp"
#include <fstream>
#include <iostream>

namespace core {

Config::Config(const std::string& filename) {
    std::ifstream configFile(filename);
    if (!configFile.is_open()) {
        std::cerr << "Failed to open config file." << std::endl;
        return;
    }

    std::string key;
    while (configFile >> key) {
        // HELP
        if (key == "defaultHelp") configFile >> defaultHelp;
        
        // FIELD PARAMETERS
        else if (key == "fieldWidth") configFile >> fieldWidth;
        else if (key == "fieldHeight") configFile >> fieldHeight;
        
        // DEFAULT GAUSSIAN PARAMETERS
        else if (key == "defaultCenterX") configFile >> defaultCenterX;
        else if (key == "defaultCenterY") configFile >> defaultCenterY;
        else if (key == "defaultSigmaX") configFile >> defaultSigmaX;
        else if (key == "defaultSigmaY") configFile >> defaultSigmaY;
        else if (key == "defaultHeight") configFile >> defaultHeight;
        
        // DEFAULT G_AUTO PARAMETERS
        else if (key == "xmin") configFile >> xmin;
        else if (key == "xmax") configFile >> xmax;
        else if (key == "ymin") configFile >> ymin;
        else if (key == "ymax") configFile >> ymax;
        else if (key == "sx_min") configFile >> sx_min;
        else if (key == "sx_max") configFile >> sx_max;
        else if (key == "sy_min") configFile >> sy_min;
        else if (key == "sy_max") configFile >> sy_max;
        else if (key == "h_min") configFile >> h_min;
        else if (key == "h_max") configFile >> h_max;
        else if (key == "count_min") configFile >> count_min;
        else if (key == "count_max") configFile >> count_max;
        
        // OUTPUT FILES
        else if (key == "defaultGnuplot") configFile >> defaultGnuplot;
        else if (key == "defaultPlotMetedata") configFile >> defaultPlotMetedata;
        else if (key == "defaultPlotVoronoi") configFile >> defaultPlotVoronoi;
        else if (key == "defaultPlotDelaunay") configFile >> defaultPlotDelaunay;
        else if (key == "defaultPlotPath") configFile >> defaultPlotPath;
        else if (key == "defaultPlot3DPath") configFile >> defaultPlot3DPath;
        
        else if (key == "defaultWrite") configFile >> defaultWrite;
        else if (key == "defaultWriteModeImage") configFile >> defaultWriteModeImage;
        else if (key == "defaultRead") configFile >> defaultRead;
        
        else if (key == "save_g") configFile >> save_g;
        
        // BINARY
        else if (key == "defaultThreshold") configFile >> defaultThreshold;
        else if (key == "defaultBinMode") configFile >> defaultBinMode;
        
        // WAVE
        else if (key == "defaultNoisy") configFile >> defaultNoisy;
        
        // KMEANS
        else if (key == "defaultKlaster") configFile >> defaultKlaster;
        else if (key == "defaultKlasterKern") configFile >> defaultKlasterKern;
        
        // PATHS
        else if (key == "defaultstartPointX") configFile >> defaultstartPointX;
        else if (key == "defaultstartPointY") configFile >> defaultstartPointY;
        else if (key == "defaultendPointX") configFile >> defaultendPointX;
        else if (key == "defaultendPointY") configFile >> defaultendPointY;
        
        else if (key == "vehicleRadius") configFile >> vehicleRadius;
        else if (key == "maxSideAngle") configFile >> maxSideAngle;
        else if (key == "maxUpDownAngle") configFile >> maxUpDownAngle;
        
        // LOGGER
        else if (key == "logFileNameInterface") configFile >> logFileNameInterface;
        else if (key == "logFileNameControl") configFile >> logFileNameControl;
        
        else if (key == "FiltrationLogLevelInterface") configFile >> FiltrationLogLevelInterface;
        else if (key == "FiltrationLogLevelControl") configFile >> FiltrationLogLevelControl;        
    }

    configFile.close();
}

}
