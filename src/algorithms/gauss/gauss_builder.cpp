#include "algorithms/gauss/gauss_builder.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <fstream>

namespace algorithms::gauss {

    Gaus::Gaus(double h_, double x0_, double y0_, double sigma_x_, double sigma_y_) : h(h_), x0(x0_), y0(y0_), sigma_x(sigma_x_), sigma_y(sigma_y_)
    {
        minX = x0 - 3.0 * sigma_x;
        maxX = x0 + 3.0 * sigma_x;
        minY = y0 - 3.0 * sigma_y;
        maxY = y0 + 3.0 * sigma_y;
        invSigmaX = 1.0 / sigma_x;
        invSigmaY = 1.0 / sigma_y;
    }
    
    void GaussBuilder::logGaussParameters(const Gaus& g) const {
        logger.debug(std::string("[GaussBuilder] Added Gaussian - ") +
                   "h: " + std::to_string(g.h) + 
                   ", x0: " + std::to_string(g.x0) +
                   ", y0: " + std::to_string(g.y0) +
                   ", σx: " + std::to_string(g.sigma_x) +
                   ", σy: " + std::to_string(g.sigma_y));
    }

    GaussBuilder::GaussBuilder(core::Logger& lg, std::mt19937& generator) : logger(lg), gen(generator) {}

    void GaussBuilder::addgauss(double h, double x0, double y0, double sigma_x, double sigma_y, std::vector<Gaus>& gaussi) {
        logger.trace("[GaussBuilder::addgauss] Adding new Gaussian distribution");
        
        gaussi.emplace_back(h, x0, y0, sigma_x, sigma_y);
        logGaussParameters(gaussi.back());
    }
    
void GaussBuilder::addgaussRandom(
    double xmin, double xmax,
    double ymin, double ymax,
    double sx_min, double sx_max,
    double sy_min, double sy_max,
    double h_min, double h_max,
    std::size_t count_min, std::size_t count_max,
    std::vector<Gaus>& gaussi)
{
    logger.info("[GaussBuilder] Автогенерация гауссов");

    std::uniform_real_distribution<> dx(xmin, xmax);
    std::uniform_real_distribution<> dy(ymin, ymax);
    std::uniform_real_distribution<> dsx(sx_min, sx_max);
    std::uniform_real_distribution<> dsy(sy_min, sy_max);
    std::uniform_real_distribution<> dh(h_min, h_max);
    std::uniform_int_distribution<std::size_t> dcount(count_min, count_max);

    std::size_t count = dcount(gen);

    for (std::size_t i = 0; i < count; ++i) {
        double x0 = dx(gen);
        double y0 = dy(gen);
        double sigma_x = dsx(gen);
        double sigma_y = dsy(gen);
        double h = dh(gen);

        addgauss(h, x0, y0, sigma_x, sigma_y, gaussi);
    }

    logger.info("[GaussBuilder] Сгенерировано " + std::to_string(count) + " гауссов");
}

void GaussBuilder::buildGaussGrid(
    int fieldWidth,
    int fieldHeight,
    int cellSize,
    const std::vector<Gaus>& gaussi)
{
    gaussGrid.cellSize = cellSize;

    gaussGrid.cols =
        (fieldWidth + cellSize - 1) / cellSize;

    gaussGrid.rows =
        (fieldHeight + cellSize - 1) / cellSize;

    gaussGrid.cells.clear();
    gaussGrid.cells.resize(
        static_cast<std::size_t>(gaussGrid.rows * gaussGrid.cols));

    for (std::size_t i = 0; i < gaussi.size(); ++i)
    {
        const Gaus& g = gaussi[i];

        int minCol = std::max(
            0,
            static_cast<int>(std::floor(g.minX / cellSize)));

        int maxCol = std::min(
            gaussGrid.cols - 1,
            static_cast<int>(std::floor(g.maxX / cellSize)));

        int minRow = std::max(
            0,
            static_cast<int>(std::floor(g.minY / cellSize)));

        int maxRow = std::min(
            gaussGrid.rows - 1,
            static_cast<int>(std::floor(g.maxY / cellSize)));

        for (int row = minRow; row <= maxRow; ++row)
        {
            for (int col = minCol; col <= maxCol; ++col)
            {
                gaussGrid.at(row, col)
                    .gaussIndices.push_back(i);
            }
        }
    }

    logger.info(
        "[GaussBuilder] Gauss grid built: " +
        std::to_string(gaussGrid.rows) + "x" +
        std::to_string(gaussGrid.cols));
}

bool GaussBuilder::getCell(
    double x,
    double y,
    int& row,
    int& col) const
{
    col = static_cast<int>(x / gaussGrid.cellSize);
    row = static_cast<int>(y / gaussGrid.cellSize);


    if (row < 0 || row >= gaussGrid.rows ||
        col < 0 || col >= gaussGrid.cols)
    {
        return false;
    }

    return true;
}
    
    void GaussBuilder::init(int fieldWidth, int fieldHeight, std::vector<std::vector<double>>& field) {
         logger.trace("[GaussBuilder::init] Initializing field");
         field.assign(static_cast<std::size_t>(fieldHeight), std::vector<double>(static_cast<std::size_t>(fieldWidth), 0));
         logger.info("Ширина = " + std::to_string(fieldWidth) + ", Длина = " +  std::to_string(fieldHeight));
         logger.info("Столбцов = " + std::to_string(field[0].size()) + ", Строк = " +  std::to_string(field.size()));
    }
    
    void GaussBuilder::generate(std::vector<std::vector<double>>& field, std::vector<Gaus>& gaussi) {
        logger.trace("[GaussBuilder::generate] Starting field generation");
        
        if (field.empty()) {
            logger.error("[GaussBuilder::generate] Pole not initialized!");
            return;
        }

        const size_t width = field[0].size();
        const size_t height = field.size();
        logger.debug(std::string("[GaussBuilder::generate] Generating field ") +
                   std::to_string(width) + "x" + std::to_string(height) + 
                   " with " + std::to_string(gaussi.size()) + " Gaussians");

        // Заполнение базовым значением
        for (auto& row : field) {
            std::fill(row.begin(), row.end(), core::MID_GRAY);
        }
        logger.debug(std::string("[GaussBuilder::generate] Base level set to core::MID_GRAY = ") + std::to_string(core::MID_GRAY));

        // Применение гауссиан
        size_t total_pixels = 0;
        size_t clamped_pixels = 0;
        
for (const auto& g : gaussi) {

    logger.trace(std::string("[GaussBuilder::generate] Applying Gaussian at (") +
                 std::to_string(g.x0) + "," +
                 std::to_string(g.y0) + ")");

    if (g.sigma_x <= 0 || g.sigma_y <= 0) {
        logger.warning("[GaussBuilder::generate] Skipping Gaussian with zero sigma");
        continue;
    }

    // Ограничиваем область 3 сигмами

    size_t minX = static_cast<size_t>(std::max(0, static_cast<int>(std::floor(g.minX))));
    size_t maxX = static_cast<size_t>(std::min(static_cast<int>(width - 1),
                                           static_cast<int>(std::ceil(g.maxX))));

    size_t minY = static_cast<size_t>(std::max(0, static_cast<int>(std::floor(g.minY))));
    size_t maxY = static_cast<size_t>(std::min(static_cast<int>(height - 1),
                                           static_cast<int>(std::ceil(g.maxY))));
    
    for (size_t y = minY; y <= maxY; ++y) {
        for (size_t x = minX; x <= maxX; ++x) {

            double dx = (static_cast<double>(x) - g.x0) * g.invSigmaX;
            double dy = (static_cast<double>(y) - g.y0) * g.invSigmaY;

            double value =
                g.h * std::exp(-(dx * dx + dy * dy) * 0.5);

            field[y][x] += value;

            if (field[y][x] < core::BLACK ||
                field[y][x] > core::WHITE)
            {
                ++clamped_pixels;
            }

            ++total_pixels;

            field[y][x] = std::clamp(
                field[y][x],
                static_cast<double>(core::BLACK),
                static_cast<double>(core::WHITE)
            );
        }
    }
}

    if (total_pixels != 0) {
       
        logger.info(std::string("[GaussBuilder::generate] Generation completed\n") +
                  "  Total pixels processed: " + std::to_string(total_pixels) + "\n" +
                  "  Clamped pixels: " + std::to_string(clamped_pixels) + "\n" +
                  "  Clamping ratio: " + 
                  std::to_string(static_cast<double>(clamped_pixels)/static_cast<double>(total_pixels)*100) + "%");
    }
  }
  
void GaussBuilder::saveGaussiansToFile(
    const std::string& filename,
    const std::vector<Gaus>& gaussi)
{
    std::ofstream out(filename);

    if (!out.is_open()) {
        logger.error("[GaussBuilder] Не удалось открыть файл для записи!");
        return;
    }

    for (const auto& g : gaussi) {
        out << "g "
            << g.x0 << " "
            << g.y0 << " "
            << g.sigma_x << " "
            << g.sigma_y << " "
            << g.h << "\n";
    }

    out.close();

    logger.info("[GaussBuilder] Гауссы сохранены в файл: " + filename);
}

double GaussBuilder::heightAt(
    double x,
    double y,
    const std::vector<Gaus>& gaussi) const
{
    double value = core::MID_GRAY;
    
    if (gaussGrid.cells.empty())
    {
        for (const auto& g : gaussi)
        {
            if (x < g.minX || x > g.maxX ||
                y < g.minY || y > g.maxY)
                continue;

            double dx = (x - g.x0) * g.invSigmaX;
            double dy = (y - g.y0) * g.invSigmaY;

            value += g.h *
                std::exp(-(dx*dx + dy*dy)*0.5);
        }

        return std::clamp(
            value,
            static_cast<double>(core::BLACK),
            static_cast<double>(core::WHITE)
        );
    }
    
    int col;
    int row;

    if (!getCell(x, y, row, col))
        return value;

    const auto& cell = gaussGrid.at(row,col);


    for (std::size_t index : cell.gaussIndices)
    {
        const auto& g = gaussi[index];


        if (x < g.minX || x > g.maxX ||
            y < g.minY || y > g.maxY)
        {
            continue;
        }


        const double dx = (x - g.x0) * g.invSigmaX;
        const double dy = (y - g.y0) * g.invSigmaY;


        value += g.h *
            std::exp(-(dx * dx + dy * dy) * 0.5);
    }


    return std::clamp(
        value,
        static_cast<double>(core::BLACK),
        static_cast<double>(core::WHITE));
}

double GaussBuilder::heightAt(
    const algorithms::geometry::PointD& point,
    const std::vector<Gaus>& gaussi) const
{
    return heightAt(point.x, point.y, gaussi);
}
    
}
