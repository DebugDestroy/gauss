#pragma once

#include <vector>

#include "core/logger.hpp"
#include "algorithms/geometry/geometry_structures.hpp" // для Pixel

namespace algorithms::path::common {
    class Conditions;   // forward declare
}

namespace algorithms::path::common {

struct GridCell {
    int row; // Индекс ячейки
    int col;

    bool traversable = true;
};

struct Grid {
    int cellSize;
    int rows; 
    int cols;
    std::vector<GridCell> cells;
};

inline std::vector<GridCell> getNeighbours(const Grid& grid,
                                           const GridCell& cell)
{
    std::vector<GridCell> neighbours;

    for (int dr = -1; dr <= 1; ++dr) {
        for (int dc = -1; dc <= 1; ++dc) {

            if (dr == 0 && dc == 0)
                continue;

            const int row = cell.row + dr;
            const int col = cell.col + dc;

            if (row < 0 || row >= grid.rows ||
                col < 0 || col >= grid.cols)
                continue;

            const auto& neighbour = grid.cells[row * grid.cols + col];

            neighbours.push_back(neighbour);
        }
    }

    return neighbours;
}

class GridBuilder {
private:
    core::Logger& logger;

public:
    GridBuilder (core::Logger& lg);
    void buildGrid(Grid &grid,
               int fieldWidth,
               int fieldHeight,
               int cellSize);
               
    void buildNavGrid(Grid& grid,
                      const std::vector<std::vector<double>>& binaryMap,
                      const Conditions& conditions,
                      int noisy);
   
    GridCell connectPointToGrid(const Grid& grid,
                           const algorithms::geometry::Pixel& point);

};

} // namespace algorithms::path::common
