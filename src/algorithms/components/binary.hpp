#pragma once

#include <vector>
#include <memory>
#include "core/logger.hpp"
#include "algorithms/gauss/pole.hpp"
#include "core/constants.hpp"

namespace algorithms::components {

enum class ThresholdMode {
    All,
    Peaks,
    Valleys
};

class Binarizer {
private:
    core::Logger& logger;

public:
    Binarizer(core::Logger& lg);

    void bin(std::vector<std::vector<double>>& CopyPole,
                    int slice,
                    std::unique_ptr<algorithms::gauss::Pole>& p,
                    ThresholdMode mode);
};

}
