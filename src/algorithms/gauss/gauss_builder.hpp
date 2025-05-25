#pragma once

#include <vector>
#include <memory>
#include "core/logger.hpp"
#include "core/constants.hpp"
#include "algorithms/gauss/pole.hpp"

namespace algorithms::gauss {
struct Gaus {
    double h, x0, y0, sigma_x, sigma_y;
    Gaus(double h, double x0, double y0, double sigma_x, double sigma_y);
};

class GaussBuilder {
private:
    core::Logger& logger;
    void logGaussParameters(const Gaus& g) const;

public:
    explicit GaussBuilder(core::Logger& lg);

    void addgauss(double h, double x0, double y0, double sigma_x, double sigma_y, std::vector<Gaus>& gaussi);
    void init(int A, int B, std::unique_ptr<Pole>& p);
    void generate(std::unique_ptr<Pole>& p, std::vector<Gaus>& gaussi);
};
}
