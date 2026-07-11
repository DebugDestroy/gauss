#include "core/config.hpp"
#include <iostream>

namespace core {

const std::unordered_set<std::string> Config::required = {
    "defaultHelp",

    "defaultfieldWidth",
    "defaultfieldHeight",

    "defaultCenterX",
    "defaultCenterY",
    "defaultSigmaX",
    "defaultSigmaY",
    "defaultHeight",

    "count_min",
    "count_max",
    "xmin",
    "xmax",
    "ymin",
    "ymax",
    "sx_min",
    "sx_max",
    "sy_min",
    "sy_max",
    "h_min",
    "h_max",

    "defaultGnuplot",
    "defaultPlotMetedata",
    "defaultPlotKmeans",
    "defaultPlotGraph",
    "defaultPlotVoronoi",
    "defaultPlotDelaunay",
    "defaultPlotGrid",
    "defaultPlotNavGrid",
    "defaultPlotGridPath",
    "defaultPlotPath",
    "PlotRRT",
    "defaultPlot3DPath",

    "defaultWrite",
    "defaultWriteModeImage",
    "defaultRead",

    "save_g",
    "defaultsave_metrics",

    "heightThresholdPixel",

    "defaultWaveNoisy",

    "defaultKlaster",
    "defaultKlasterKern",

    "defaultgridWidth",
    "defaultgridNoisy",

    "startPixelX",
    "startPixelY",
    "goalPixelX",
    "goalPixelY",
    "defaultConnectMode",

    "vehicleRadiusPixel",
    "maxSideAngle",
    "maxUpDownAngle",

    "logFileNameInterface",
    "logFileNameControl",

    "FiltrationLogLevelInterface",
    "FiltrationLogLevelControl",

    "seedMode",

    "maxIterations",
    "startWorldX",
    "startWorldY",
    "goalWorldX",
    "goalWorldY",
    "vehicleRadiusWorld",
    "heightThresholdWorld",
    "interpEdge",
    "interpCollision",
    "interpAngle",
    "step",
    "goalRadius",
    "goalBias"
};

Config::Config(const std::string& filename) {
    std::ifstream configFile(filename);
    if (!configFile.is_open()) {
        throw std::runtime_error(
            "Failed to open config file: " + filename);
    }

    std::string key;
    while (configFile >> key) {
         // Запись команды в список и проверка на дубликат
         if (!readKeys.insert(key).second)
         {
             throw std::runtime_error(
                 "Duplicate config parameter '" +
                 key +
                 "'");
         }
         
        // HELP
        if (key == "defaultHelp")
        readParameter(configFile, defaultHelp, "defaultHelp");
        
        // FIELD PARAMETERS
        else if (key == "defaultfieldWidth")
        readParameter(configFile, defaultfieldWidth, "defaultfieldWidth");
        else if (key == "defaultfieldHeight")
        readParameter(configFile, defaultfieldHeight, "defaultfieldHeight");
        
        // DEFAULT GAUSSIAN PARAMETERS
        else if (key == "defaultCenterX")
        readParameter(configFile, defaultCenterX, "defaultCenterX");
        else if (key == "defaultCenterY")
        readParameter(configFile, defaultCenterY, "defaultCenterY");
        else if (key == "defaultSigmaX") 
        readParameter(configFile, defaultSigmaX, "defaultSigmaX");
        else if (key == "defaultSigmaY")
        readParameter(configFile, defaultSigmaY, "defaultSigmaY");
        else if (key == "defaultHeight")
        readParameter(configFile, defaultHeight, "defaultHeight");
        
        // DEFAULT G_AUTO PARAMETERS
        else if (key == "count_min")
        readParameter(configFile, count_min, "count_min");
        else if (key == "count_max")
        readParameter(configFile, count_max, "count_max");
        else if (key == "xmin")
        readParameter(configFile, xmin, "xmin");
        else if (key == "xmax")
        readParameter(configFile, xmax, "xmax");
        else if (key == "ymin")
        readParameter(configFile, ymin, "ymin");
        else if (key == "ymax")
        readParameter(configFile, ymax, "ymax");
        else if (key == "sx_min")
        readParameter(configFile, sx_min, "sx_min");
        else if (key == "sx_max") 
        readParameter(configFile, sx_max, "sx_max");
        else if (key == "sy_min")
        readParameter(configFile, sy_min, "sy_min");
        else if (key == "sy_max")
        readParameter(configFile, sy_max, "sy_max");
        else if (key == "h_min")
        readParameter(configFile, h_min, "h_min");
        else if (key == "h_max")
        readParameter(configFile, h_max, "h_max");
        
        // OUTPUT FILES
        else if (key == "defaultGnuplot")
        readParameter(configFile, defaultGnuplot, "defaultGnuplot");
        else if (key == "defaultPlotMetedata")
        readParameter(configFile, defaultPlotMetedata, "defaultPlotMetedata");
        else if (key == "defaultPlotKmeans")
        readParameter(configFile, defaultPlotKmeans, "defaultPlotKmeans");
        else if (key == "defaultPlotGraph")
        readParameter(configFile, defaultPlotGraph, "defaultPlotGraph");
        else if (key == "defaultPlotVoronoi")
        readParameter(configFile, defaultPlotVoronoi, "defaultPlotVoronoi");
        else if (key == "defaultPlotDelaunay")
        readParameter(configFile, defaultPlotDelaunay, "defaultPlotDelaunay");
        else if (key == "defaultPlotGrid")
        readParameter(configFile, defaultPlotGrid, "defaultPlotGrid");
        else if (key == "defaultPlotNavGrid")
        readParameter(configFile, defaultPlotNavGrid, "defaultPlotNavGrid");
        else if (key == "defaultPlotGridPath")
        readParameter(configFile, defaultPlotGridPath, "defaultPlotGridPath");
        else if (key == "defaultPlotPath")
        readParameter(configFile, defaultPlotPath, "defaultPlotPath");
        else if (key == "PlotRRT")
        readParameter(configFile, PlotRRT, "PlotRRT");
        else if (key == "defaultPlot3DPath")
        readParameter(configFile, defaultPlot3DPath, "defaultPlot3DPath");
        
        else if (key == "defaultWrite")
        readParameter(configFile, defaultWrite, "defaultWrite");
        else if (key == "defaultWriteModeImage")
        readParameter(configFile, defaultWriteModeImage, "defaultWriteModeImage");
        else if (key == "defaultRead")
        readParameter(configFile, defaultRead, "defaultRead");
        
        else if (key == "save_g")
        readParameter(configFile, save_g, "save_g");
        
        else if (key == "defaultsave_metrics")
        readParameter(configFile, defaultsave_metrics, "defaultsave_metrics");
        
        // BINARY
        else if (key == "heightThresholdPixel")
        readParameter(configFile, heightThresholdPixel, "heightThresholdPixel");
        
        // WAVE
        else if (key == "defaultWaveNoisy")
        readParameter(configFile, defaultWaveNoisy, "defaultWaveNoisy");
        
        // KMEANS
        else if (key == "defaultKlaster")
        readParameter(configFile, defaultKlaster, "defaultKlaster");
        else if (key == "defaultKlasterKern")
        readParameter(configFile, defaultKlasterKern, "defaultKlasterKern");
        
        // GRID
        else if (key == "defaultgridWidth")
        readParameter(configFile, defaultgridWidth, "defaultgridWidth");
        else if (key == "defaultgridNoisy")
        readParameter(configFile, defaultgridNoisy, "defaultgridNoisy");
        
        // PATHS
        else if (key == "startPixelX")
        readParameter(configFile, startPixelX, "startPixelX");
        else if (key == "startPixelY")
        readParameter(configFile, startPixelY, "startPixelY");
        else if (key == "goalPixelX")
        readParameter(configFile, goalPixelX, "goalPixelX");
        else if (key == "goalPixelY")
        readParameter(configFile, goalPixelY, "goalPixelY");
        else if (key == "defaultConnectMode") {
            readParameter(configFile, defaultConnectMode, "defaultConnectMode");

            if (defaultConnectMode == "NearestK") {
                readParameter(configFile, defaultNearestVerticesCount, "defaultNearestVerticesCount");
            }
        }
        
        else if (key == "vehicleRadiusPixel")
        readParameter(configFile, vehicleRadiusPixel, "vehicleRadiusPixel");
        else if (key == "maxSideAngle")
        readParameter(configFile, maxSideAngle, "maxSideAngle");
        else if (key == "maxUpDownAngle")
        readParameter(configFile, maxUpDownAngle, "maxUpDownAngle");
        
        // LOGGER
        else if (key == "logFileNameInterface")
        readParameter(configFile, logFileNameInterface, "logFileNameInterface");
        else if (key == "logFileNameControl")
        readParameter(configFile, logFileNameControl, "logFileNameControl");
        
        else if (key == "FiltrationLogLevelInterface")
        readParameter(configFile, FiltrationLogLevelInterface, "FiltrationLogLevelInterface");
        else if (key == "FiltrationLogLevelControl")
        readParameter(configFile, FiltrationLogLevelControl, "FiltrationLogLevelControl");  
        
        // SEED
        else if (key == "seedMode") {
            readParameter(configFile, seedMode, "seedMode");

            if (seedMode == "Fixed") {
                readParameter(configFile, seed, "seed");
            }
        }
        
        // RRT
        else if (key == "maxIterations")
        readParameter(configFile, maxIterations, "maxIterations");
        else if (key == "startWorldX")
        readParameter(configFile, startWorldX, "startWorldX");
        else if (key == "startWorldY")
        readParameter(configFile, startWorldY, "startWorldY");
        else if (key == "goalWorldX")
        readParameter(configFile, goalWorldX, "goalWorldX");
        else if (key == "goalWorldY")
        readParameter(configFile, goalWorldY, "goalWorldY");
        else if (key == "vehicleRadiusWorld")
        readParameter(configFile, vehicleRadiusWorld, "vehicleRadiusWorld");
        else if (key == "heightThresholdWorld")
        readParameter(configFile, heightThresholdWorld, "heightThresholdWorld");
        else if (key == "interpEdge")
        readParameter(configFile, interpEdge, "interpEdge");
        else if (key == "interpCollision")
        readParameter(configFile, interpCollision, "interpCollision");
        else if (key == "interpAngle")
        readParameter(configFile, interpAngle, "interpAngle");
        else if (key == "step")
        readParameter(configFile, step, "step");
        else if (key == "goalRadius")
        readParameter(configFile, goalRadius, "goalRadius");
        else if (key == "goalBias")
        readParameter(configFile, goalBias, "goalBias");
        
        // ELSE
        else throw std::runtime_error("Unknown configuration parameter: '" + key + "'");     
    }

    configFile.close();
    
    for (const auto& name : required)
{
    if (!readKeys.contains(name))
    {
        throw std::runtime_error(
            "Required config parameter '" +
            name +
            "' is missing");
    }
}

}

}
