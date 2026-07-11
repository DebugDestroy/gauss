#include "algorithms/components/binary.hpp"
#include "core/constants.hpp"
#include <string>

namespace algorithms::components {

Binarizer::Binarizer(core::Logger& lg) : logger(lg) {
        logger.trace("[binary] Инициализация бинаризации поля");
    }
    
void Binarizer::bin(std::vector<std::vector<double>>& binaryMap,
                    int slice,
                    const std::vector<std::vector<double>>& field)
{
    logger.info(std::string("[binary::bin] Бинаризация данных, slice=") +
                std::to_string(slice));

    if (field.empty()) {
        logger.error("[binary::bin] Ошибка: данные высот не инициализированы!");
        return;
    }

    binaryMap.assign(field.size(), std::vector<double>(field[0].size()));

    for (size_t y = 0; y < field.size(); ++y)
    {
        for (size_t x = 0; x < field[0].size(); ++x)
        {
            binaryMap[y][x] =
                (std::abs(field[y][x] - core::MID_GRAY) >= slice)
                    ? core::WHITE
                    : core::BLACK;
        }
    }

    logger.info("[binary::bin] Бинаризация завершена");
}
}
