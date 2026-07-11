#pragma once

#include <vector>
#include <memory>
#include "core/logger.hpp"

namespace algorithms::components {

class Binarizer {
private:
    core::Logger& logger;

public:
    Binarizer(core::Logger& lg);

    void bin(std::vector<std::vector<double>>& binaryMap,
                    int slice,
                    const std::vector<std::vector<double>>& field);
};

}
