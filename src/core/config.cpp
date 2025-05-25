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
        if (key == "fieldWidth") configFile >> fieldWidth;
        else if (key == "fieldHeight") configFile >> fieldHeight;
        else if (key == "defaultCenterX") configFile >> defaultCenterX;
        else if (key == "defaultCenterY") configFile >> defaultCenterY;
        else if (key == "defaultSigmaX") configFile >> defaultSigmaX;
        else if (key == "defaultSigmaY") configFile >> defaultSigmaY;
        else if (key == "defaultHeight") configFile >> defaultHeight;
        else if (key == "defaultGnuplot") configFile >> defaultGnuplot;
        else if (key == "defaultPlotMetedata") configFile >> defaultPlotMetedata;
        else if (key == "defaultPlotVoronoi") configFile >> defaultPlotVoronoi;
        else if (key == "defaultPlotDelaunay") configFile >> defaultPlotDelaunay;
        else if (key == "defaultPlotPath") configFile >> defaultPlotPath;
        else if (key == "defaultWrite") configFile >> defaultWrite;
        else if (key == "defaultWriteModeImage") configFile >> defaultWriteModeImage;
        else if (key == "defaultRead") configFile >> defaultRead;
        else if (key == "defaultThreshold") configFile >> defaultThreshold;
        else if (key == "defaultBinMode") configFile >> defaultBinMode;
        else if (key == "defaultNoisy") configFile >> defaultNoisy;
        else if (key == "defaultKlaster") configFile >> defaultKlaster;
        else if (key == "defaultKlasterKern") configFile >> defaultKlasterKern;
        else if (key == "defaultHelp") configFile >> defaultHelp;
        else if (key == "logFileNameInterface") configFile >> logFileNameInterface;
        else if (key == "logFileNameControl") configFile >> logFileNameControl;
        else if (key == "FiltrationLogLevelInterface") configFile >> FiltrationLogLevelInterface;
        else if (key == "FiltrationLogLevelControl") configFile >> FiltrationLogLevelControl;
        else if (key == "defaultstartPointX") configFile >> defaultstartPointX;
        else if (key == "defaultstartPointY") configFile >> defaultstartPointY;
        else if (key == "defaultendPointX") configFile >> defaultendPointX;
        else if (key == "defaultendPointY") configFile >> defaultendPointY;
        else if (key == "vehicleRadius") configFile >> vehicleRadius;
        else if (key == "maxSideAngle") configFile >> maxSideAngle;
        else if (key == "maxUpDownAngle") configFile >> maxUpDownAngle;
        else if (key == "defaultPlot3DPath") configFile >> defaultPlot3DPath;
    }

    configFile.close();
}

}
