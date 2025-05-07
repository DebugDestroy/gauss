#pragma once
#include <vector>

// Локальные заголовки
#include "algorithms/KMeans.hpp"
#include "core/Logger.hpp"
#include "core/Config.hpp"
#include "services/Geometry.hpp" // для PointD
#include "algorithms/Component.hpp" // для Component
#include "services/ColorGenerator.hpp"

class ClusterService {
private:
    Logger& logger;
    KMeans kMeans;
    ColorGenerator colorGenerator;
    std::vector<std::vector<double>> kMeansData;

public:
    ClusterService(Logger& log) : logger(log), kMeans(log), colorGenerator() {}
    
    std::vector<std::array<int, 3>> getClusterColors(int clusterCount) {
        return colorGenerator.generateColors(clusterCount, logger);
    }

    std::vector<PointD> getClusterCenters(const std::vector<Component>& components, const Config& config) {
        std::vector<PointD> centers;
        for (const auto& component : components) {
            if (std::isnan(component.center_x) || std::isnan(component.center_y)) {
                logger.logMessage(LogLevel::Warning, "Skipping invalid cluster center (NaN)");
                continue;
            }
            if (component.center_x >= 0 && component.center_x < config.fieldWidth &&
                component.center_y >= 0 && component.center_y < config.fieldHeight) {
                centers.emplace_back(component.center_x, component.center_y);
            } else {
                logger.logMessage(LogLevel::Error, 
                    std::string("Invalid cluster center: (") + std::to_string(component.center_x) + 
                    ", " + std::to_string(component.center_y) + ")");
            }
        }
        return centers;
    }

    void prepareKMeansData(const std::vector<std::vector<double>>& copyPole) {
        kMeansData.clear();
        for (size_t y = 0; y < copyPole.size(); ++y) {
            for (size_t x = 0; x < copyPole[0].size(); ++x) {
                if (copyPole[y][x] > 0) { 
                    kMeansData.push_back({static_cast<double>(x), static_cast<double>(y)});
                }
            }
        }
        logger.logMessage(LogLevel::Debug, 
            std::string("Prepared ") + std::to_string(kMeansData.size()) + " points for K-means");
    }

   void applyClusterResults(const KMeans::ClusterResult& result, 
                           std::vector<std::vector<double>>& pixelMatrix) {
        auto colors = getClusterColors(result.centers.size());
        for (size_t i = 0; i < result.labels.size(); ++i) {
            int label = result.labels[i];
            if (label >= 0 && label < static_cast<int>(colors.size())) {
                int x = static_cast<int>(kMeansData[i][0]);
                int y = static_cast<int>(kMeansData[i][1]);
                const auto& color = colors[label];
                pixelMatrix[y][x] = (color[0] + color[1] + color[2]) / 3.0;
            }
        }
    }


    const std::vector<std::vector<double>>& getKMeansData() const {
        return kMeansData;
    }

    KMeans& getKMeans() {
        return kMeans;
    }
};
