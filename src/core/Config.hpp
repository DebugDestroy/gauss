#pragma once
#include <fstream>   // Для std::ifstream
#include <string>    // Для std::string
#include <iostream>  // Для std::cerr

// Чтение конфигурационного файла
class Config {
public:
    int fieldWidth, fieldHeight;
    double defaultX, defaultY, defaultSx, defaultSy, defaultH;
    std::string FiltrationLogLevelInterface, FiltrationLogLevelControl;
    std::string logFileNameInterface;
    std::string logFileNameControl;
    std::string defaultHelp;
    std::string defaultGnuplot, defaultPlotMetedata, defaultPlotVoronoi, defaultPlotDelaunay, defaultPlotPath, defaultWrite, defaultRead, defaultWriteModeImage, defaultBinMode, defaultPlot3DPath;
    int defaultSlice, defaultNoisy, defaultKlaster, defaultKlasterKern;
    double defaultpointA_x, defaultpointA_y, defaultpointB_x, defaultpointB_y;
    double vehicleRadius, maxSideAngle, maxUpDownAngle;

    Config(const std::string& filename) {
        
        std::ifstream configFile(filename);
        if (!configFile.is_open()) {
            std::cerr << "Failed to open config file." << std::endl;
            return;
        }

     std::string key;
while (configFile >> key) { // Считываем ключи
    if (key == "fieldWidth") configFile >> fieldWidth;
    else if (key == "fieldHeight") configFile >> fieldHeight;
    else if (key == "defaultX") configFile >> defaultX;
    else if (key == "defaultY") configFile >> defaultY;
    else if (key == "defaultSx") configFile >> defaultSx;
    else if (key == "defaultSy") configFile >> defaultSy;
    else if (key == "defaultH") configFile >> defaultH;  
    else if (key == "defaultGnuplot") configFile >> defaultGnuplot;  
    else if (key == "defaultPlotMetedata") configFile >> defaultPlotMetedata;
    else if (key == "defaultPlotVoronoi") configFile >> defaultPlotVoronoi;
    else if (key == "defaultPlotDelaunay") configFile >> defaultPlotDelaunay;
    else if (key == "defaultPlotPath") configFile >> defaultPlotPath;   
    else if (key == "defaultWrite") configFile >> defaultWrite;
    else if (key == "defaultWriteModeImage") configFile >> defaultWriteModeImage;
    else if (key == "defaultRead") configFile >> defaultRead;
    else if (key == "defaultSlice") configFile >> defaultSlice;
    else if (key == "defaultBinMode") configFile >> defaultBinMode;
    else if (key == "defaultNoisy") configFile >> defaultNoisy;
    else if (key == "defaultKlaster") configFile >> defaultKlaster;
    else if (key == "defaultKlasterKern") configFile >> defaultKlasterKern;   
    else if (key == "defaultHelp") configFile >> defaultHelp;  
    else if (key == "logFileNameInterface") configFile >> logFileNameInterface;
    else if (key == "logFileNameControl") configFile >> logFileNameControl;
    else if (key == "FiltrationLogLevelInterface") configFile >> FiltrationLogLevelInterface;
    else if (key == "FiltrationLogLevelControl") configFile >> FiltrationLogLevelControl;
    else if (key == "defaultpointA_x") configFile >> defaultpointA_x;
    else if (key == "defaultpointA_y") configFile >> defaultpointA_y;
    else if (key == "defaultpointB_x") configFile >> defaultpointB_x;
    else if (key == "defaultpointB_y") configFile >> defaultpointB_y;
    else if (key == "vehicleRadius") configFile >> vehicleRadius;    
    else if (key == "maxSideAngle") configFile >> maxSideAngle;
    else if (key == "maxUpDownAngle") configFile >> maxUpDownAngle;
    else if (key == "defaultPlot3DPath") configFile >> defaultPlot3DPath;
}

        configFile.close();
    }
};
