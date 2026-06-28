#include "algorithms/components/binary.hpp"
#include "core/constants.hpp"
#include <string>

namespace algorithms::components {

Binarizer::Binarizer(core::Logger& lg) : logger(lg) {
        logger.trace("[binary] Инициализация бинаризации поля");
    }
    
 void Binarizer::bin(std::vector<std::vector<double>>& binaryMap, 
            int slice, 
            const std::vector<std::vector<double>>& field, 
            ThresholdMode mode) {
        logger.info(std::string("[binary::bin] Бинаризация данных, slice=") + 
                  std::to_string(slice) + ", mode=" + 
                  (mode == ThresholdMode::Peaks ? "Peaks" : 
                   mode == ThresholdMode::Valleys ? "Valleys" : "All"));

        if (field.empty()) {
            logger.error("[binary::bin] Ошибка: данные высот не инициализированы!");
            return;
        }
        
        binaryMap = field;
        int symmetric_slice = 2*core::MID_GRAY - slice;

        for (int x = 0; x < (int)field[0].size(); ++x) {
            for (int y = 0; y < (int)field.size(); ++y) {
                double value = field[y][x];
                switch (mode) {
                    case ThresholdMode::All: {
                        bool is_peak = (slice >= core::MID_GRAY) ? (value > slice) : (value > symmetric_slice);
                        bool is_valley = (slice >= core::MID_GRAY) ? (value < symmetric_slice) : (value < slice);
                        binaryMap[y][x] = (is_peak || is_valley) ? core::WHITE : core::BLACK;
                        break;
                    }
                    case ThresholdMode::Peaks:
                        binaryMap[y][x] = (slice >= core::MID_GRAY) ? (value > slice) : (value > symmetric_slice);
                        break;
                    case ThresholdMode::Valleys:
                        binaryMap[y][x] = (slice >= core::MID_GRAY) ? (value < symmetric_slice) : (value < slice);
                        break;
                }
            } 
        }

        logger.info("[binary::bin] Бинаризация завершена");
    }
}
