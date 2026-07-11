#pragma once

#include <vector>
#include <array>
#include <cstdint>
#include <random>

#include "core/logger.hpp"

namespace visualization {

struct Color {
    uint8_t r, g, b;
};

class ColorGenerator {
private:
core::Logger& logger;
std::mt19937 &gen;

void logColorConversion(double hue, const std::array<int, 3>& rgb) const;
static uint8_t clampColorComponent(double value);

public:
    std::vector<std::array<int, 3>> generateColors(std::size_t numColors);
    
    explicit ColorGenerator(core::Logger& lg, std::mt19937& generator);
};

}
