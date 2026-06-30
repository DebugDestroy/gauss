#pragma once
#include "algorithms/geometry/geometry_structures.hpp"
#include "algorithms/path/common/path_metrics.hpp"

#include <vector>
#include <chrono>

namespace statistics {
class StatisticsManager {
public:
    // Начать отсчет
    void startTimer();
    
    // Завершить таймер и посчитать время
    void finishTimer(algorithms::path::PathMetrics& metrics);

    // Вычислить метрики по пути
    void computePathMetrics(
        algorithms::path::PathMetrics& metrics,
        const std::vector<algorithms::geometry::Pixel>& path);
        
    // Сбросить метрики
    void reset(algorithms::path::PathMetrics& metrics);
            
private:
    std::chrono::high_resolution_clock::time_point startTime;
};

} // namespace statistics
