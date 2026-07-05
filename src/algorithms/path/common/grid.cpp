#include "algorithms/path/common/grid.hpp"
#include "algorithms/path/common/conditions.hpp"

namespace algorithms::path::common {

GridBuilder::GridBuilder(core::Logger& lg) : logger(lg) {}

void GridBuilder::buildGrid(Grid &grid,
               int fieldWidth,
               int fieldHeight,
               int cellSize) {
    grid.cellSize = cellSize;

    const int rows = fieldHeight / cellSize;
    const int cols = fieldWidth / cellSize;
    grid.rows = rows;
    grid.cols = cols;
    
    grid.cells.resize(rows * cols);

for (int row = 0; row < rows; ++row) {
    for (int col = 0; col < cols; ++col) {
    
        auto& cell = grid.cells[row * cols + col];
        logger.trace(std::string("Cell col=") + std::to_string(col) + " and row=" + std::to_string(row));
        cell.row = row;
        cell.col = col;
    }
}
}

void GridBuilder::buildNavGrid(Grid& grid,
                              const std::vector<std::vector<double>>& binaryMap,
                              const Conditions& conditions,
                              int noisy)
{
    logger.trace("Building nav grid");

    for (auto& cell : grid.cells)
        cell.traversable = conditions.isCellFree(cell, binaryMap,  grid.cellSize, noisy);

    logger.trace("Nav grid built");
}

GridCell GridBuilder::connectPointToGrid(const Grid& grid,
                                     const algorithms::geometry::Pixel& point)
{
    logger.trace("Connecting point to grid");

    const int row = point.y / grid.cellSize;
    const int col = point.x / grid.cellSize;

    return grid.cells[row * grid.cols + col];
}

} // namespace algorithms::path::common
