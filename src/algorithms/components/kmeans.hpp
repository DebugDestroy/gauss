#pragma once

#include <vector>
#include <string>

#include "core/logger.hpp"

namespace algorithms::components {
class KMeans {
private:
    core::Logger& logger;

    void logCenters(const std::vector<std::vector<double>>& clusterCenters, const std::string& prefix = "") const;
    void logPoint(const std::vector<double>& point, const std::string& prefix = "") const;
    double distance(const std::vector<double>& a, const std::vector<double>& b) const;

public:
    struct ClusterResult {
        std::vector<int> labels;
        std::vector<std::vector<double>> centers;
    };

    std::vector<int> labels;
    std::vector<std::vector<double>> centers;

    explicit KMeans(core::Logger& lg);

    void initializeCenters(const std::vector<std::vector<double>>& data, std::size_t k);
    ClusterResult cluster(const std::vector<std::vector<double>>& data, std::size_t k);
    ClusterResult kmeansWithKernels(const std::vector<std::vector<double>>& data, std::size_t k, std::size_t kernelSize);
};
}
