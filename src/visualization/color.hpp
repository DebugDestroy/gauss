#pragma once

#include <vector>
#include <array>
#include <cstdint>

#include "core/logger.hpp"

namespace visualization {

struct Color {
    uint8_t r, g, b;
};

class ColorGenerator {
public:
    static std::vector<std::array<int, 3>> generateColors(int numColors, core::Logger& logger);
};

}
