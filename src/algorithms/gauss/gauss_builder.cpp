#include "algorithms/gauss/gauss_builder.hpp"
#include <cmath>
#include <algorithm>
#include <string>
#include <random>
#include <fstream>

namespace algorithms::gauss {

    Gaus::Gaus(double h, double x0, double y0, double sigma_x, double sigma_y) : h(h), x0(x0), y0(y0), sigma_x(sigma_x), sigma_y(sigma_y) {}
    
    void GaussBuilder::logGaussParameters(const Gaus& g) const {
        logger.debug(std::string("[GaussBuilder] Added Gaussian - ") +
                   "h: " + std::to_string(g.h) + 
                   ", x0: " + std::to_string(g.x0) +
                   ", y0: " + std::to_string(g.y0) +
                   ", σx: " + std::to_string(g.sigma_x) +
                   ", σy: " + std::to_string(g.sigma_y));
    }

    GaussBuilder::GaussBuilder(core::Logger& lg) : logger(lg) {}

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
    int count_min, int count_max,
    std::vector<Gaus>& gaussi)
{
    logger.info("[GaussBuilder] Автогенерация гауссов");

    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<> dx(xmin, xmax);
    std::uniform_real_distribution<> dy(ymin, ymax);
    std::uniform_real_distribution<> dsx(sx_min, sx_max);
    std::uniform_real_distribution<> dsy(sy_min, sy_max);
    std::uniform_real_distribution<> dh(h_min, h_max);
    std::uniform_int_distribution<> dcount(count_min, count_max);

    int count = dcount(gen);

    for (int i = 0; i < count; ++i) {
        double x0 = dx(gen);
        double y0 = dy(gen);
        double sigma_x = dsx(gen);
        double sigma_y = dsy(gen);
        double h = dh(gen);

        addgauss(h, x0, y0, sigma_x, sigma_y, gaussi);
    }

    logger.info("[GaussBuilder] Сгенерировано " + std::to_string(count) + " гауссов");
}
    
    void GaussBuilder::init(int A, int B, std::unique_ptr<Pole>& p) {
        logger.trace(std::string("[GaussBuilder::init] Initializing field with size ") +
                   std::to_string(A) + "x" + std::to_string(B));
        
        if (!p) {
            logger.debug("[GaussBuilder::init] Creating new Pole object");
            p = std::make_unique<Pole>(A, B, logger);
        } else {
            logger.debug("[GaussBuilder::init] Resizing existing Pole object");
            p->resize(A, B);
        }
        
        logger.info("[GaussBuilder::init] Field initialized successfully");
    }
    
    void GaussBuilder::generate(std::unique_ptr<Pole>& p, std::vector<Gaus>& gaussi) {
        logger.trace("[GaussBuilder::generate] Starting field generation");
        
        if (p == nullptr) {
            logger.error("[GaussBuilder::generate] Pole not initialized!");
            return;
        }

        const size_t width = p->field[0].size();
        const size_t height = p->field.size();
        logger.debug(std::string("[GaussBuilder::generate] Generating field ") +
                   std::to_string(width) + "x" + std::to_string(height) + 
                   " with " + std::to_string(gaussi.size()) + " Gaussians");

        // Заполнение базовым значением
        for (auto& row : p->field) {
            std::fill(row.begin(), row.end(), core::MID_GRAY);
        }
        logger.debug(std::string("[GaussBuilder::generate] Base level set to core::MID_GRAY = ") + std::to_string(core::MID_GRAY));

        // Применение гауссиан
        size_t total_pixels = 0;
        size_t clamped_pixels = 0;
        
        for (const auto& g : gaussi) {
            logger.trace(std::string("[GaussBuilder::generate] Applying Gaussian at (") +
                       std::to_string(g.x0) + "," + std::to_string(g.y0) + ")");
            
            if (g.sigma_x <= 0 || g.sigma_y <= 0) {
        logger.warning("[GaussBuilder::generate] Skipping Gaussian with zero sigma");
        continue;
    }
    
            for (size_t x = 0; x < width; ++x) {
                for (size_t y = 0; y < height; ++y) {
                    double value = g.h * exp(-((pow((x - g.x0)/g.sigma_x, 2) + 
                                             pow((y - g.y0)/g.sigma_y, 2))/2));
                    p->field[y][x] += value;
                    
                    // Логирование clamping
                    if (p->field[y][x] <= core::BLACK || p->field[y][x] >= core::WHITE) {
                        clamped_pixels++;
                    }
                    total_pixels++;
                    
                    p->field[y][x] = std::clamp(p->field[y][x], 
                           static_cast<double>(core::BLACK), 
                           static_cast<double>(core::WHITE));
                }
            }
        }

    if (total_pixels != 0) {
       
        logger.info(std::string("[GaussBuilder::generate] Generation completed\n") +
                  "  Total pixels processed: " + std::to_string(total_pixels) + "\n" +
                  "  Clamped pixels: " + std::to_string(clamped_pixels) + "\n" +
                  "  Clamping ratio: " + 
                  std::to_string((double)clamped_pixels/total_pixels*100) + "%");
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
}
