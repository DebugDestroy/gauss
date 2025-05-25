#pragma once
#include <string>

namespace core {

class Config {
public:
    int fieldWidth, fieldHeight;
    double defaultCenterX, defaultCenterY, defaultSigmaX, defaultSigmaY, defaultHeight;
    std::string FiltrationLogLevelInterface, FiltrationLogLevelControl;
    std::string logFileNameInterface, logFileNameControl;
    std::string defaultHelp;
    std::string defaultGnuplot, defaultPlotMetedata, defaultPlotVoronoi, defaultPlotDelaunay, defaultPlotPath, defaultWrite, defaultRead, defaultWriteModeImage, defaultBinMode, defaultPlot3DPath;
    int defaultThreshold, defaultNoisy, defaultKlaster, defaultKlasterKern;
    double defaultstartPointX, defaultstartPointY, defaultendPointX, defaultendPointY;
    double vehicleRadius, maxSideAngle, maxUpDownAngle;

    explicit Config(const std::string& filename);
};

}
