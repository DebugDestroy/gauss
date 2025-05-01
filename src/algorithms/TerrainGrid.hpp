#pragma once
#include <vector>
#include <utility> // для std::pair
#include <cmath> // для std::atan, M_PI
#include <string>
#include <sstream>

// Локальные заголовки
#include "core/Logger.hpp"
#include "core/Geometry.hpp" // для PointD, Edge
#include "core/Pole.hpp" // для Pole

// Структура для хранения данных ячейки сетки
struct GridCell {
    PointD position;       // Координаты центра ячейки
    double height = 0;     // Высота в этой точке
    double slopeX = 0;     // Наклон по оси X (в радианах)
    double slopeY = 0;     // Наклон по оси Y (в радианах)
    bool visited = false;
};

// Класс для работы с сеткой рельефа
class TerrainGrid {
private:
    Logger& logger;

    void logGridInitialization(int width, int height) const {
        std::ostringstream oss;
        oss << "[TerrainGrid] Grid initialized:\n"
            << "  Dimensions: " << width << "x" << height << "\n"
            << "  Cell spacing: 1.0 units\n"
            << "  Y-axis inverted (top-left origin)";
        logger.debug(oss.str());
    }

    void logSlopeCalculation(int x, int y, double slopeX, double slopeY) const {
        if (logger.getLogLevel() <= LogLevel::Debug) {
            std::ostringstream oss;
            oss << "[TerrainGrid] Slope at (" << x << "," << y << "):\n"
                << "  slopeX: " << slopeX << " (" << atan(slopeX) * 180/M_PI << "°)\n"
                << "  slopeY: " << slopeY << " (" << atan(slopeY) * 180/M_PI << "°)";
            logger.debug(oss.str());
        }
    }

public:
    std::vector<std::vector<GridCell>> cells;

    TerrainGrid(int width, int height, Logger& lg) : logger(lg) {
        logger.trace("[TerrainGrid] Constructor called");
        initialize(width, height);
    }

    void initialize(int width, int height) {
        logger.trace("[TerrainGrid::initialize] Initializing grid");
        
        if (width <= 0 || height <= 0) {
            logger.error(std::string("[TerrainGrid::initialize] Invalid dimensions: ") + 
                       std::to_string(width) + "x" + std::to_string(height));
            throw std::invalid_argument("Grid dimensions must be positive");
        }

        cells.resize(height, std::vector<GridCell>(width));

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                cells[y][x].position = PointD(x + 0.5, (height - y - 1) + 0.5);
                
                if (x == 0 && y == 0) {  // Логируем только первую ячейку для примера
                    logger.trace(std::string("[TerrainGrid] Sample cell [0][0] position: (") + 
                               std::to_string(cells[y][x].position.x) + ", " +
                               std::to_string(cells[y][x].position.y) + ")");
                }
            }
        }

        logGridInitialization(width, height);
        logger.info("[TerrainGrid] Grid initialization completed");
    }
    
    std::pair<double, double> getEdgeSlopes(const Edge& edge) const {
        logger.trace(std::string("[TerrainGrid::getEdgeSlopes] Calculating slopes for edge: (") +
                   std::to_string(edge.a.x) + "," + std::to_string(edge.a.y) + ") -> (" +
                   std::to_string(edge.b.x) + "," + std::to_string(edge.b.y) + ")");

        PointD dir = {edge.b.x - edge.a.x, edge.b.y - edge.a.y};
        double length = std::hypot(dir.x, dir.y);
        
        if (length < 1e-6) {
            logger.debug("[TerrainGrid::getEdgeSlopes] Edge has zero length");
            return {0, 0};
        }

        dir.x /= length;
        dir.y /= length;
        int x = static_cast<int>((edge.a.x + edge.b.x) / 2);
        int y = static_cast<int>((edge.a.y + edge.b.y) / 2);

        if (x <= 0 || y <= 0 || x >= static_cast<int>(cells[0].size()-1) || y >= static_cast<int>(cells.size()-1)) {
            logger.warning(std::string("[TerrainGrid::getEdgeSlopes] Edge center out of bounds: (") +
                         std::to_string(x) + "," + std::to_string(y) + ")");
            return {0, 0};
        }

        double forwardSlope = atan(cells[y][x].slopeX * dir.x + cells[y][x].slopeY * dir.y) * 180.0 / M_PI;
        double sideSlope = atan(cells[y][x].slopeX * (-dir.y) + cells[y][x].slopeY * dir.x) * 180.0 / M_PI;

        logger.debug(std::string("[TerrainGrid::getEdgeSlopes] Calculated slopes:\n") +
                   "  Forward: " + std::to_string(forwardSlope) + "°\n" +
                   "  Side: " + std::to_string(sideSlope) + "°");
        
        return {forwardSlope, sideSlope};
    }

    void calculateSlopes(const Pole& elevationData) {
        logger.trace("[TerrainGrid::calculateSlopes] Starting slope calculations");
        
        if (elevationData.field.empty()) {
            logger.error("[TerrainGrid::calculateSlopes] Empty elevation data provided");
            throw std::invalid_argument("Pole data is empty");
        }

        const int height = cells.size();
        const int width = cells[0].size();

        if (height != static_cast<int>(elevationData.field.size()) || width != static_cast<int>(elevationData.field[0].size())) {
            logger.error(std::string("[TerrainGrid::calculateSlopes] Dimension mismatch: ") +
                       std::to_string(width) + "x" + std::to_string(height) + " vs " +
                       std::to_string(elevationData.field[0].size()) + "x" + 
                       std::to_string(elevationData.field.size()));
            throw std::invalid_argument("Grid and elevation data dimensions mismatch");
        }

        size_t processed_cells = 0;
        for (int y = 1; y < height-1; ++y) {
            for (int x = 1; x < width-1; ++x) {
                cells[y][x].height = elevationData.field[y][x];
                cells[y][x].slopeX = (elevationData.field[y][x+1] - elevationData.field[y][x-1]) / 2.0;
                cells[y][x].slopeY = (elevationData.field[y+1][x] - elevationData.field[y-1][x]) / 2.0;
                processed_cells++;

                if (x == 1 && y == 1) {  // Логируем только первую ячейку для примера
                    logSlopeCalculation(x, y, cells[y][x].slopeX, cells[y][x].slopeY);
                }
            }
        }

        logger.info(std::string("[TerrainGrid::calculateSlopes] Completed\n") +
                  "  Processed cells: " + std::to_string(processed_cells) + "\n" +
                  "  Border cells skipped: " + 
                  std::to_string(width*height - processed_cells));
    }

    void resetVisited() {
        size_t reset_count = 0;
        for (auto& row : cells) {
            for (auto& cell : row) {
                if (cell.visited) {
                    cell.visited = false;
                    reset_count++;
                }
            }
        }
        logger.debug(std::string("[TerrainGrid::resetVisited] Reset ") + 
                    std::to_string(reset_count) + " visited cells");
    }
};
