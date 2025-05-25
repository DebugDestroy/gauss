#include "visualization/color.hpp"

#include <cmath>        // fmod, fabs
#include <ctime>        // time()
#include <cstdlib>      // rand, srand
#include <string>       // std::string
#include <stdexcept>    // std::invalid_argument
#include <algorithm>    // std::clamp

#include "core/constants.hpp"

namespace visualization {

static constexpr double MIN_COLOR_VALUE = 50.0;

static void logColorConversion(double hue, const std::array<int, 3>& rgb, core::Logger& logger) {
    logger.debug("[ColorGenerator] Generated color - HSL: " +
                 std::to_string(hue) + "°, RGB: (" +
                 std::to_string(rgb[0]) + ", " +
                 std::to_string(rgb[1]) + ", " +
                 std::to_string(rgb[2]) + ")");
}

static uint8_t clampColorComponent(double value) {
    return static_cast<uint8_t>(std::clamp(
        value,
        MIN_COLOR_VALUE,
        static_cast<double>(core::WHITE)
    ));
}

std::vector<std::array<int, 3>> ColorGenerator::generateColors(int numColors, core::Logger& logger) {
    logger.trace("[ColorGenerator::generateColors] Starting color generation for " + 
                 std::to_string(numColors) + " colors");

    if (numColors <= 0) {
        logger.error("[ColorGenerator::generateColors] Invalid number of colors: " +
                     std::to_string(numColors));
        throw std::invalid_argument("Number of colors must be greater than 0.");
    }

    std::vector<std::array<int, 3>> colors;
    unsigned int seed = static_cast<unsigned int>(time(0));
    srand(seed);
    logger.debug("[ColorGenerator] Seed initialized: " + std::to_string(seed));

    for (int i = 0; i < numColors; ++i) {
        double hue = static_cast<double>(i) / numColors * 360.0 + rand() % 30;
        hue = fmod(hue, 360.0);

        // HSL to RGB conversion
        double r, g, b;
        double C = 1.0;
        double X = C * (1 - fabs(fmod((hue / 60.0), 2) - 1));
        double m = 0.0;

        if (hue < 60) {
            r = C; g = X; b = 0;
        } else if (hue < 120) {
            r = X; g = C; b = 0;
        } else if (hue < 180) {
            r = 0; g = C; b = X;
        } else if (hue < 240) {
            r = 0; g = X; b = C;
        } else if (hue < 300) {
            r = X; g = 0; b = C;
        } else {
            r = C; g = 0; b = X;
        }

        std::array<int, 3> rgb = {
            clampColorComponent((r + m) * core::WHITE),
            clampColorComponent((g + m) * core::WHITE),
            clampColorComponent((b + m) * core::WHITE)
        };

        colors.push_back(rgb);
        logColorConversion(hue, rgb, logger);
    }

    logger.info("[ColorGenerator::generateColors] Successfully generated " +
                std::to_string(colors.size()) + " distinct colors");

    return colors;
}

}
