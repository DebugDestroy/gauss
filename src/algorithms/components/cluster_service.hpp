#pragma once
#include <vector>
#include "core/config.hpp"
#include "core/logger.hpp"
#include "algorithms/components/components_analysis.hpp"
#include "algorithms/components/kmeans.hpp"
#include "visualization/color.hpp"
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::components {
class ClusterService {
private:
    core::Logger& logger;
    KMeans& kMeans;
    visualization::ColorGenerator& colorGenerator;
    std::vector<std::vector<double>> kMeansData;

public:
    ClusterService(core::Logger& log, KMeans& kmeans, visualization::ColorGenerator& colorGen);
    
    std::vector<std::array<int, 3>> getClusterColors(int clusterCount);
    std::vector<algorithms::geometry::PointD> getClusterCenters(const std::vector<Component>& components, const core::Config& config);
    void prepareKMeansData(const std::vector<std::vector<double>>& copyPole);
    void applyClusterResults(const KMeans::ClusterResult& result, std::vector<std::vector<double>>& pixelMatrix);
    const std::vector<std::vector<double>>& getKMeansData() const;
};
}
