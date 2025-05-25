#include "algorithms/components/binary.hpp"
#include <string>

namespace algorithms::components {

Binarizer::Binarizer(core::Logger& lg) : logger(lg) {
        logger.trace("[binary] Инициализация бинаризации поля");
    }
    
 void Binarizer::bin(std::vector<std::vector<double>>& CopyPole, 
            int slice, 
            std::unique_ptr<algorithms::gauss::Pole>& p, 
            ThresholdMode mode) {
        logger.info(std::string("[binary::bin] Бинаризация данных, slice=") + 
                  std::to_string(slice) + ", mode=" + 
                  (mode == ThresholdMode::Peaks ? "Peaks" : 
                   mode == ThresholdMode::Valleys ? "Valleys" : "All"));

        if (p == nullptr) {
            logger.error("[binary::bin] Ошибка: данные высот не инициализированы!");
            return;
        }
        
        CopyPole = p->field;
        int symmetric_slice = 2*core::MID_GRAY - slice;

        for (int x = 0; x < (int)p->field[0].size(); ++x) {
            for (int y = 0; y < (int)p->field.size(); ++y) {
                double value = p->field[y][x];
                switch (mode) {
                    case ThresholdMode::All: {
                        bool is_peak = (slice >= core::MID_GRAY) ? (value > slice) : (value > symmetric_slice);
                        bool is_valley = (slice >= core::MID_GRAY) ? (value < symmetric_slice) : (value < slice);
                        CopyPole[y][x] = (is_peak || is_valley) ? core::WHITE : core::BLACK;
                        break;
                    }
                    case ThresholdMode::Peaks:
                        CopyPole[y][x] = (slice >= core::MID_GRAY) ? (value > slice) : (value > symmetric_slice);
                        break;
                    case ThresholdMode::Valleys:
                        CopyPole[y][x] = (slice >= core::MID_GRAY) ? (value < symmetric_slice) : (value < slice);
                        break;
                }
            } 
        }

        logger.info("[binary::bin] Бинаризация завершена");
    }
}
