#pragma once

#include <vector>
#include "core/constants.hpp"
#include "core/logger.hpp"
#include "algorithms/components/components_analysis.hpp"

namespace algorithms::components {
class Copier {
private:
    core::Logger& logger;

public:
    explicit Copier(core::Logger& lg);
    void removeNoise(std::vector<std::vector<double>>& field, const std::vector<Component>& components);
};
}
