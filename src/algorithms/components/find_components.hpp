#pragma once

#include <vector>
#include <memory>
#include "core/constants.hpp"
#include "core/logger.hpp"
#include "algorithms/gauss/pole.hpp"
#include "algorithms/components/components_analysis.hpp"

namespace algorithms::components {
class ComponentCalculator {
private:
    core::Logger& logger;

    int incrementAndCollect(std::vector<std::vector<double>>& componenta, 
                            std::vector<std::vector<double>>& CopyPole, 
                            int x, int y, int i, int& pixelCount);

public:
    explicit ComponentCalculator(core::Logger& lg);
    
    void wave(int noisy, 
              std::vector<Component>& componenti, 
              std::vector<std::vector<double>>& CopyPole, 
              std::unique_ptr<algorithms::gauss::Pole>& p);
};
}
