#pragma once

#include <vector>
#include <string>
#include "core/logger.hpp"  // Предполагаем, что Logger здесь

namespace algorithms::gauss {

class Pole {
private:
    core::Logger& logger;

    void logResizeOperation(int old_rows, int old_cols, int new_rows, int new_cols) const;

public:
    std::vector<std::vector<double>> field;

    Pole(int A, int B, core::Logger& lg);

    void resize(int A, int B);

    void logFieldStats();
};
}
