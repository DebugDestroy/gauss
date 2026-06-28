#pragma once

#include <vector>
#include <memory>
#include "core/logger.hpp"

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

    void bin(std::vector<std::vector<double>>& binaryMap,
                    int slice,
                    const std::vector<std::vector<double>>& field,
                    ThresholdMode mode);
};

}
