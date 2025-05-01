#pragma once
#include <string>
#include "io/BmpHandler.hpp"  // Для BmpWriteMode
#include "algorithms/ComponentCalculator.hpp"  // Для ThresholdMode

struct DispatcherParams {
    std::string s;
    int A, B;
    double h, x, y, sx, sy;
    std::string filename;
    BmpWriteMode bmp_mode;
    int slice;
    ThresholdMode bin_mode;
    int noisy;
    int k, kk;
    double pointA_x, pointA_y, pointB_x, pointB_y;
};
