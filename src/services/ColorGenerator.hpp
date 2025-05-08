#pragma once
#include <vector>       // Для std::vector
#include <array>        // Для std::array
#include <string>       // Для std::string
#include <cmath>        // Для fmod, fabs
#include <ctime>        // Для time()
#include <cstdlib>      // Для srand(), rand()
#include <algorithm>    // Для std::min, std::max
#include <stdexcept>    // Для std::invalid_argument

// Локальные заголовки
#include "core/Constants.hpp"  // Подключаем константы

struct Color{ // структура для хранения rgb цвета
    uint8_t r, g, b;
};

class ColorGenerator {
private:
    static constexpr double MIN_COLOR_VALUE = 50.0;
    
    static void logColorConversion(double hue, const std::array<int, 3>& rgb, Logger& logger) {
        logger.debug(std::string("[ColorGenerator] Generated color - HSL: ") + 
                   std::to_string(hue) + "°, RGB: (" +
                   std::to_string(rgb[0]) + ", " +
                   std::to_string(rgb[1]) + ", " +
                   std::to_string(rgb[2]) + ")");
    }

    static uint8_t clampColorComponent(double value) {
        return static_cast<uint8_t>(std::clamp(
            value,
            MIN_COLOR_VALUE,
            static_cast<double>(Constants::WHITE)
        ));
    }

public:
    static std::vector<std::array<int, 3>> generateColors(int numColors, Logger& logger) {
        logger.trace(std::string("[ColorGenerator::generateColors] Starting color generation for ") + 
                     std::to_string(numColors) + " colors");

        if (numColors <= 0) {
            logger.error(std::string("[ColorGenerator::generateColors] Invalid number of colors: ") + 
                        std::to_string(numColors));
            throw std::invalid_argument("Number of colors must be greater than 0.");
        }

        std::vector<std::array<int, 3>> colors;
        unsigned int seed = static_cast<unsigned int>(time(0));
        srand(seed);
        logger.debug(std::string("[ColorGenerator] Seed initialized: ") + std::to_string(seed));

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

            // Scale to MIN_COLOR_VALUE-Constants::WHITE range
            std::array<int, 3> rgb = {
                clampColorComponent((r + m) * Constants::WHITE),
                clampColorComponent((g + m) * Constants::WHITE),
                clampColorComponent((b + m) * Constants::WHITE)
            };

            colors.push_back(rgb);
            logColorConversion(hue, rgb, logger);
        }

        logger.info(std::string("[ColorGenerator::generateColors] Successfully generated ") + 
                   std::to_string(colors.size()) + " distinct colors");
        return colors;
    }
};
