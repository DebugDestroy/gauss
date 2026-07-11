#include <vector>
#include <string>
#include <sstream>
#include <cmath> // для std::fabs, std::sqrt

// Локальные заголовки
#include "algorithms/components/components_analysis.hpp"

namespace algorithms::components {
    
    void Component::logComponentStats() {
        std::ostringstream oss;
        oss << "Full component stats:\n"
            << "  Bounds: x[" << min_x << ".." << max_x << "], y[" << min_y << ".." << max_y << "]\n"
            << "  Center: (" << center_x << ", " << center_y << ")\n"
            << "  Pixels: " << pixelCount << "\n"
            << "  Eigenvalues: " << eigenvalue1 << ", " << eigenvalue2 << "\n"
            << "  Eigenvector1: (" << eigenvec1_x << ", " << eigenvec1_y << ")\n"
            << "  Eigenvector2: (" << eigenvec2_x << ", " << eigenvec2_y << ")";
            
        logger.trace("[Component] " + oss.str());
    }
    
   // Основной конструктор (для готовой матрицы)
    Component::Component(core::Logger& lg, const std::vector<algorithms::geometry::Pixel>& inputComponenta) : logger(lg), componenta(inputComponenta) {
        logger.trace("[Component] Constructor started");
        calculate_metadata();
        logger.trace(std::string("[Component] Created with pixel count: ") + std::to_string(pixelCount));
        logComponentStats(); // Теперь печатает все поля
    }

void Component::calculate_metadata() {
    int count;
    double sum_x;
    double sum_y;

    double cov_xx;
    double cov_xy;
    double cov_yy;

    double trace;
    double det;
    double discriminant;

    double dx;
    double dy;

    double norm1;
    double norm2;

    // Проверка на пустую компоненту
    if (componenta.empty()) {
        logger.warning("[Component::calculate_metadata] Empty component");

        min_x = min_y = 0;
        max_x = max_y = 0;
        center_x = center_y = 0.0;
        pixelCount = 0;

        eigenvalue1 = eigenvalue2 = 0.0;

        eigenvec1_x = 1.0;
        eigenvec1_y = 0.0;
        eigenvec2_x = 0.0;
        eigenvec2_y = 1.0;

        return;
    }

    logger.trace("[Component::calculate_metadata] Starting");

    // Инициализация
   min_x = componenta[0].x;
   max_x = componenta[0].x;

   min_y = componenta[0].y;
   max_y = componenta[0].y;

    sum_x = 0.0;
    sum_y = 0.0;

    pixelCount = static_cast<int>(componenta.size());
    count = pixelCount;

    // Границы и центр масс
    for (const auto& p : componenta) {
        min_x = std::min(min_x, p.x);
        max_x = std::max(max_x, p.x);

        min_y = std::min(min_y, p.y);
        max_y = std::max(max_y, p.y);

        sum_x += p.x;
        sum_y += p.y;
    }

    center_x = sum_x / count;
    center_y = sum_y / count;

    logger.debug(
        "[Component::calculate_metadata] Center = (" +
        std::to_string(center_x) + ", " +
        std::to_string(center_y) + ")"
    );

    // Ковариационная матрица
    cov_xx = 0.0;
    cov_xy = 0.0;
    cov_yy = 0.0;

    for (const auto& p : componenta) {
        dx = p.x - center_x;
        dy = p.y - center_y;

        cov_xx += dx * dx;
        cov_xy += dx * dy;
        cov_yy += dy * dy;
    }

    cov_xx /= count;
    cov_xy /= count;
    cov_yy /= count;

    logger.debug(
        "[Component::calculate_metadata] Covariance matrix: xx=" +
        std::to_string(cov_xx) +
        ", xy=" + std::to_string(cov_xy) +
        ", yy=" + std::to_string(cov_yy)
    );

    // Собственные значения
    trace = cov_xx + cov_yy;
    det = cov_xx * cov_yy - cov_xy * cov_xy;

    discriminant = trace * trace - 4.0 * det;

    if (discriminant < 0.0) {
        discriminant = 0.0;
    }

    eigenvalue1 = (trace + std::sqrt(discriminant)) / 2.0;
    eigenvalue2 = (trace - std::sqrt(discriminant)) / 2.0;

    // Собственные векторы
    if (std::fabs(cov_xy) > core::EPSILON) {

        eigenvec1_x = cov_xy;
        eigenvec1_y = eigenvalue1 - cov_xx;

        eigenvec2_x = cov_xy;
        eigenvec2_y = eigenvalue2 - cov_xx;

    } else {

        if (cov_xx >= cov_yy) {
            eigenvec1_x = 1.0;
            eigenvec1_y = 0.0;

            eigenvec2_x = 0.0;
            eigenvec2_y = 1.0;
        } else {
            eigenvec1_x = 0.0;
            eigenvec1_y = 1.0;

            eigenvec2_x = 1.0;
            eigenvec2_y = 0.0;
        }
    }

    // Нормировка первого вектора
    norm1 = std::sqrt(
        eigenvec1_x * eigenvec1_x +
        eigenvec1_y * eigenvec1_y
    );

    if (norm1 > core::EPSILON) {
        eigenvec1_x /= norm1;
        eigenvec1_y /= norm1;
    }

    // Нормировка второго вектора
    norm2 = std::sqrt(
        eigenvec2_x * eigenvec2_x +
        eigenvec2_y * eigenvec2_y
    );

    if (norm2 > core::EPSILON) {
        eigenvec2_x /= norm2;
        eigenvec2_y /= norm2;
    }

    logger.trace("[Component::calculate_metadata] Finished");
}
}
