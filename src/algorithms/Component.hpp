#pragma once
#include <vector>
#include <string>
#include <sstream>
#include <cmath> // для std::fabs, std::sqrt

// Локальные заголовки
#include "core/Constants.hpp"  // Подключаем константы
#include "core/Logger.hpp"

enum class ThresholdMode {// Для bin
    All,     // и ямы, и горы (по модулю)
    Peaks,   // только горы (без модуля, значение > slice)
    Valleys  // только ямы (без модуля, значение < -slice)
};

class Component {
private:
     Logger& logger;
    
    // Полный лог всех полей компонента
    void logComponentStats() {
        std::ostringstream oss;
        oss << "Full component stats:\n"
            << "  Bounds: x[" << min_x << ".." << max_x << "], y[" << min_y << ".." << max_y << "]\n"
            << "  Center: (" << center_x << ", " << center_y << ")\n"
            << "  Pixels: " << pixelCount << "\n"
            << "  Eigenvalues: " << eigenvalue1 << ", " << eigenvalue2 << "\n"
            << "  Eigenvector1: (" << eigenvec1_x << ", " << eigenvec1_y << ")\n"
            << "  Eigenvector2: (" << eigenvec2_x << ", " << eigenvec2_y << ")";
            
        logger.info("[Component] " + oss.str());
    }
    
public:
    std::vector<std::vector<double>> componenta;
    int min_x, min_y, max_x, max_y;
    double center_x, center_y;
    double eigenvec1_x, eigenvec1_y;
    double eigenvec2_x, eigenvec2_y;
    double eigenvalue1, eigenvalue2;
    int pixelCount;  // Поле для хранения количества пикселей
    
   // Основной конструктор (для готовой матрицы)
    Component(Logger& lg, const std::vector<std::vector<double>>& inputComponenta, int count) : logger(lg), componenta(inputComponenta), pixelCount(count) {
        logger.trace("[Component] Constructor started");
        calculate_metadata();
        logger.debug(std::string("[Component] Created with pixel count: ") + std::to_string(pixelCount));
        logComponentStats(); // Теперь печатает все поля
    }

    void calculate_metadata() {
    // Находим границы и центр
    min_x = componenta[0].size();
    min_y = componenta.size();
    max_x = 0;
    max_y = 0;
    double sum_x = 0, sum_y = 0;
    int count = 0;
    
    logger.trace("[Component::calculate_metadata] Starting calculate_metadata");
    logger.debug("[Component::calculate_metadata] Scanning component pixels...");
    for (int y = 0; y < static_cast<int>(componenta.size()); ++y) {
        for (int x = 0; x < static_cast<int>(componenta[0].size()); ++x) {
            if (std::fabs(componenta[y][x] - Constants::WHITE) < Constants::EPSILON) {
                min_x = std::min(min_x, x);
                min_y = std::min(min_y, y);
                max_x = std::max(max_x, x);
                max_y = std::max(max_y, y);
                sum_x += x;
                sum_y += y;
                count++;
            }
        }
    }

    if (count == 0) {
        logger.warning("[Component::calculate_metadata] Component has no valid pixels!");
        center_x = center_y = 0;
        return;
    }

    center_x = sum_x / count;
    center_y = sum_y / count;
    
    logger.debug(std::string("[Component::calculate_metadata] Center calculated: (") + 
                    std::to_string(center_x) + std::string(", ") + 
                    std::to_string(center_y) + std::string(")"));

    // Корректный расчет ковариационной матрицы
    double cov_xx = 0, cov_xy = 0, cov_yy = 0;
    logger.trace("[Component::calculate_metadata] Calculating covariance matrix...");
    for (int y = min_y; y <= max_y; ++y) {
        for (int x = min_x; x <= max_x; ++x) {
            if (std::fabs(componenta[y][x] - Constants::WHITE) < Constants::EPSILON) {  // Исправлено условие
                double dx = x - center_x;
                double dy = y - center_y;
                cov_xx += dx * dx;
                cov_xy += dx * dy;
                cov_yy += dy * dy;
            }
        }
    }

    cov_xx /= count;
    cov_xy /= count;
    cov_yy /= count;
    
    logger.debug(std::string("[Component::calculate_metadata] Covariance matrix: xx=") + std::to_string(cov_xx) +
                    std::string(", xy=") + std::to_string(cov_xy) +
                    std::string(", yy=") + std::to_string(cov_yy));
                    
    // Собственные значения
    double trace = cov_xx + cov_yy;
    double det = cov_xx * cov_yy - cov_xy * cov_xy;
    eigenvalue1 = (trace + sqrt(trace * trace - 4 * det)) / 2;
    eigenvalue2 = (trace - sqrt(trace * trace - 4 * det)) / 2;
    
    logger.debug(std::string("[Component::calculate_metadata] Eigenvalues calculated: λ1=") + 
                    std::to_string(eigenvalue1) + 
                    std::string(", λ2=") + std::to_string(eigenvalue2));
                    
    // Корректный расчет собственных векторов
    if (fabs(cov_xy) > Constants::EPSILON) {
        logger.trace("[Component::calculate_metadata] Calculating eigenvectors for non-diagonal matrix");
        // Первый собственный вектор (для eigenvalue1)
        eigenvec1_x = cov_xy;
        eigenvec1_y = eigenvalue1 - cov_xx;
        
        // Второй собственный вектор (для eigenvalue2)
        eigenvec2_x = cov_xy;
        eigenvec2_y = eigenvalue2 - cov_xx;
    } else {
        logger.trace("[Component::calculate_metadata] Using identity matrix for eigenvectors (diagonal case)");
        eigenvec1_x = 1; eigenvec1_y = 0;
        eigenvec2_x = 0; eigenvec2_y = 1;
    }

    // Нормализация векторов
    double norm1 = sqrt(eigenvec1_x * eigenvec1_x + eigenvec1_y * eigenvec1_y);
    double norm2 = sqrt(eigenvec2_x * eigenvec2_x + eigenvec2_y * eigenvec2_y);
    
    eigenvec1_x /= norm1; eigenvec1_y /= norm1;
    eigenvec2_x /= norm2; eigenvec2_y /= norm2;
    
    logger.debug(std::string("[Component::calculate_metadata] Normalized eigenvectors: ") +
                    "v1=(" + std::to_string(eigenvec1_x) + "," + std::to_string(eigenvec1_y) + "), " +
                    "v2=(" + std::to_string(eigenvec2_x) + "," + std::to_string(eigenvec2_y) + ")");
        
        logger.trace("[Component::calculate_metadata] calculate_metadata completed");
}
};
