#include "visualization/color.hpp"

#include <cmath>        // fmod, fabs
#include <string>       // std::string
#include <stdexcept>    // std::invalid_argument
#include <algorithm>    // std::clamp

#include "core/constants.hpp"

namespace visualization {

ColorGenerator::ColorGenerator(core::Logger& lg, std::mt19937& generator) : logger(lg), gen(generator) {}

void ColorGenerator::logColorConversion(double hue, const std::array<int, 3>& rgb) const {
    logger.debug("[ColorGenerator] Generated color - HSL: " +
                 std::to_string(hue) + "°, RGB: (" +
                 std::to_string(rgb[0]) + ", " +
                 std::to_string(rgb[1]) + ", " +
                 std::to_string(rgb[2]) + ")");
}

uint8_t ColorGenerator::clampColorComponent(double value) {
    return static_cast<uint8_t>(std::clamp(
        value,
        0.0,
        static_cast<double>(core::WHITE)
    ));
}

std::vector<std::array<int, 3>> ColorGenerator::generateColors(std::size_t numColors) {
    logger.trace("[ColorGenerator::generateColors] Starting color generation for " + 
                 std::to_string(numColors) + " colors");

    if (numColors <= 0) {
        logger.error("[ColorGenerator::generateColors] Invalid number of colors: " +
                     std::to_string(numColors));
        return {};
    }

    std::vector<std::array<int, 3>> colors;
    std::uniform_real_distribution<double> hueNoise(0.0, 30.0);
    
    for (std::size_t i = 0; i < numColors; ++i) {
        double hue =
              static_cast<double>(i) / static_cast<double>(numColors) * 360.0 +
              hueNoise(gen);
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
        logColorConversion(hue, rgb);
    }

    logger.debug("[ColorGenerator::generateColors] Successfully generated " +
                std::to_string(colors.size()) + " distinct colors");

    return colors;
}

}
