#include <vector>
#include <algorithm>
#include <string>
#include <cmath>

// Локальные заголовки
#include "algorithms/components/copier_service.hpp"

namespace algorithms::components {

Copier::Copier(core::Logger& lg)
    : logger(lg) {}

void Copier::removeNoise(
        std::vector<std::vector<double>>& field,
        const std::vector<Component>& components)
{
    logger.trace("[Copier::removeNoise] Starting noise removal");

    logger.debug(
        std::string("[Copier::removeNoise] Ширина = ") +
        std::to_string(field[0].size()) + "Длина = " +
        std::to_string(field.size()) +
        " field, " +
        std::to_string(components.size()) +
        " components");

    // Очищаем поле
    std::size_t totalZeroed = 0;

    for (auto& row : field) {
        totalZeroed += row.size();
        std::fill(row.begin(), row.end(), core::BLACK);
    }

    logger.debug(
        "[Copier::removeNoise] Zeroed " +
        std::to_string(totalZeroed) +
        " pixels");

    std::size_t totalCopied = 0;
    std::size_t componentsProcessed = 0;

    for (const auto& comp : components) {

        const auto& compData = comp.componenta;
        std::size_t componentPixels = 0;

        for (const auto& pixel : compData) {

            int x = pixel.x;
            int y = pixel.y;

            if (x < 0 || x >= static_cast<int>(field[0].size()) ||
                y < 0 || y >= static_cast<int>(field.size()))
            {
                continue;
            }

            field[static_cast<std::size_t>(y)][static_cast<std::size_t>(x)] = core::WHITE;

            componentPixels++;
            totalCopied++;
        }

        componentsProcessed++;

        logger.debug(
            "[Copier::removeNoise] Component #" +
            std::to_string(componentsProcessed) +
            ": copied " +
            std::to_string(componentPixels) +
            " pixels");
    }

    logger.debug(
        std::string("[Copier::removeNoise] Completed. Stats:\n") +
        "  Total components processed: " +
        std::to_string(componentsProcessed) + "\n" +
        "  Total pixels copied: " +
        std::to_string(totalCopied) + "\n" +
        "  Zeroed pixels: " +
        std::to_string(totalZeroed - totalCopied));
}

} // namespace algorithms::components
