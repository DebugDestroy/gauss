#pragma once

#include <vector>
#include <memory>
#include <random>

#include "core/logger.hpp"
#include "core/constants.hpp"
#include "algorithms/geometry/geometry_structures.hpp"

namespace algorithms::gauss {
struct Gaus {
    double h, x0, y0, sigma_x, sigma_y;
    Gaus(double h_, double x0_, double y0_, double sigma_x_, double sigma_y_);
};

class GaussBuilder {
private:
    core::Logger& logger;
    void logGaussParameters(const Gaus& g) const;
    std::mt19937 &gen;
    
public:
    explicit GaussBuilder(core::Logger& lg, std::mt19937& generator);

    void addgauss(double h, double x0, double y0, double sigma_x, double sigma_y, std::vector<Gaus>& gaussi);
    void addgaussRandom(
    double xmin, double xmax,
    double ymin, double ymax,
    double sx_min, double sx_max,
    double sy_min, double sy_max,
    double h_min, double h_max,
    std::size_t count_min, std::size_t count_max,
    std::vector<Gaus>& gaussi);
    
    void init(int A, int B, std::vector<std::vector<double>>& field);
    void generate(std::vector<std::vector<double>>& field, std::vector<Gaus>& gaussi);
    void saveGaussiansToFile(const std::string& filename, const std::vector<Gaus>& gaussi);
    
   static double heightAt(
    double x,
    double y,
    const std::vector<Gaus>& gaussi);
    
   static double heightAt(
    const algorithms::geometry::PointD& point,
    const std::vector<Gaus>& gaussi);
};
}
