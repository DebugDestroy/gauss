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
    double minX;
    double maxX;
    double minY;
    double maxY;
    double invSigmaX;
    double invSigmaY;
    Gaus(double h_, double x0_, double y0_, double sigma_x_, double sigma_y_);
};

struct GaussGridCell
{
    std::vector<std::size_t> gaussIndices;
};

struct GaussGrid
{
    int cellSize;
    int rows;
    int cols;
    std::vector<GaussGridCell> cells;

    const GaussGridCell& at(int row, int col) const
    {
        return cells[static_cast<std::size_t>(row * cols + col)];
    }

    GaussGridCell& at(int row, int col)
    {
        return cells[static_cast<std::size_t>(row * cols + col)];
    }
};

class GaussBuilder {
private:
    core::Logger& logger;
    void logGaussParameters(const Gaus& g) const;
    std::mt19937 &gen;
    GaussGrid gaussGrid;
    
    bool getCell(
        double x,
        double y,
        int& row,
        int& col) const;
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
    
    void buildGaussGrid(
    int fieldWidth,
    int fieldHeight,
    int cellSize,
    const std::vector<Gaus>& gaussi);
    
    void init(int A, int B, std::vector<std::vector<double>>& field);
    void generate(std::vector<std::vector<double>>& field, std::vector<Gaus>& gaussi);
    void saveGaussiansToFile(const std::string& filename, const std::vector<Gaus>& gaussi);
    
    double heightAt(
    double x,
    double y,
    const std::vector<Gaus>& gaussi) const;
    
    double heightAt(
    const algorithms::geometry::PointD& point,
    const std::vector<Gaus>& gaussi) const;
};
}
