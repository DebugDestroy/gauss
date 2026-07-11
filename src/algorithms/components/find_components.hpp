#pragma once

#include <vector>
#include <memory>
#include "core/constants.hpp"
#include "core/logger.hpp"
#include "algorithms/components/components_analysis.hpp"
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::components {
class ComponentCalculator {
private:
    core::Logger& logger;
    
    
    void collectComponent(const std::vector<std::vector<double>>& binaryMap,
                          std::vector<std::vector<bool>>& visited,
                          std::size_t startX,
                          std::size_t startY,
                          std::vector<algorithms::geometry::Pixel>& pixels);

public:
    explicit ComponentCalculator(core::Logger& lg);
    
    void wave(std::size_t noisy, 
              std::vector<Component>& componenti, 
              std::vector<std::vector<double>>& binaryMap, 
              std::vector<std::vector<double>>& field);
};
}
