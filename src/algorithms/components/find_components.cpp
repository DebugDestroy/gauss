#include <vector>
#include <algorithm> // для std::remove_if, std::find, std::all_of
#include <cmath> // для std::isnan, std::sqrt
#include <numeric> // для std::iota
#include <string>
#include <sstream>
#include <limits> // для std::numeric_limits

// Локальные заголовки
#include "algorithms/components/find_components.hpp"

namespace algorithms::components {

    int ComponentCalculator::incrementAndCollect(std::vector<std::vector<double>>& componenta, 
                          std::vector<std::vector<double>>& CopyPole, 
                          int x, int y, int i, int& pixelCount) {
        logger.trace("[ComponentCalculator::incrementAndCollect] Вход: "
                "координаты=(" + std::to_string(x) + "," + std::to_string(y) + "), "
                "глубина_рекурсии=" + std::to_string(i) + ", "
                "текущий_счетчик=" + std::to_string(pixelCount));

        if (x < 1 || y < 1 || x > (int)componenta[0].size() - 2 || 
            y > (int)componenta.size() - 2 || std::fabs(CopyPole[y][x] - core::WHITE) > core::EPSILON) {
            logger.trace("[ComponentCalculator::incrementAndCollect] Выход за границы или неподходящее значение");
            return pixelCount;
        }

        if (std::fabs(CopyPole[y][x] - core::WHITE) < core::EPSILON) {
            CopyPole[y][x] = core::BLACK;
            pixelCount++;
            componenta[y][x] = core::WHITE;
            
            logger.trace("[ComponentCalculator::incrementAndCollect] Обработан пиксель: "
                   "координаты=(" + std::to_string(x) + "," + std::to_string(y) + "), "
                   "новый_счетчик=" + std::to_string(pixelCount));

            // Рекурсивные вызовы
            incrementAndCollect(componenta, CopyPole, x + 1, y, i + 1, pixelCount);
            incrementAndCollect(componenta, CopyPole, x - 1, y, i + 1, pixelCount);
            incrementAndCollect(componenta, CopyPole, x, y + 1, i + 1, pixelCount);
            incrementAndCollect(componenta, CopyPole, x, y - 1, i + 1, pixelCount);
        }
        
        return pixelCount;
    }

    ComponentCalculator::ComponentCalculator(core::Logger& lg) : logger(lg) {
        logger.trace("[ComponentCalculator] Инициализация калькулятора компонент");
    }
    
    void ComponentCalculator::wave(int noisy, 
             std::vector<Component>& componenti, 
             std::vector<std::vector<double>>& CopyPole, 
             std::unique_ptr<algorithms::gauss::Pole>& p) {
        logger.info(std::string("[ComponentCalculator::wave] Начало волнового алгоритма, порог=") + 
                  std::to_string(noisy));

        if (p == nullptr) {
            logger.error("[ComponentCalculator::wave] Ошибка: данные высот не инициализированы!");
            return;
        }

        int rows = p->field.size();
        int cols = (rows > 0) ? p->field[0].size() : 0;
        std::vector<Component> noiseComponents;

        logger.debug(std::string("Размер данных: ") + std::to_string(cols) + "x" + std::to_string(rows));

        for (int y = 2; y < rows - 2; ++y) {
            for (int x = 2; x < cols - 2; ++x) {
                if (std::fabs(CopyPole[y][x] - core::WHITE) < core::EPSILON) {
                    std::vector<std::vector<double>> componentData(rows, std::vector<double>(cols, core::BLACK));
                    int pixelCount = 0;
                    incrementAndCollect(componentData, CopyPole, x, y, 0, pixelCount);

                    if (pixelCount >= noisy) {
                        Component component(logger, componentData, pixelCount);
                        componenti.push_back(component);
                        
                        logger.debug(std::string("[ComponentCalculator::wave] Значимая компонента: ") +
                                   "pixels=" + std::to_string(pixelCount) +
                                   ", center=(" + std::to_string(component.center_x) + "," + 
                                   std::to_string(component.center_y) + ")" +
                                   ", size=(" + std::to_string(component.max_x - component.min_x) + 
                                   "x" + std::to_string(component.max_y - component.min_y) + ")");
                    } else {
                        Component noiseComponent(logger, componentData, pixelCount);
                        noiseComponents.push_back(noiseComponent);
                        
                        logger.trace(std::string("[ComponentCalculator::wave] Шумовая компонента: ") +
                                    "pixels=" + std::to_string(pixelCount) +
                                    ", center=(" + std::to_string(noiseComponent.center_x) + "," + 
                                    std::to_string(noiseComponent.center_y) + ")");

                        // Удаляем шум из основного поля
                        for (int i = 0; i < rows; ++i) {
                            for (int j = 0; j < cols; ++j) {
                                if (std::fabs(componentData[i][j] - core::WHITE) < core::EPSILON) {
                                    p->field[i][j] = core::MID_GRAY;
                                }
                            }
                        }
                    }
                }
            }
        }
        
        logger.info("[ComponentCalculator::wave] Обработка завершена");
        logger.debug(std::string("Найдено значимых компонент: ") + std::to_string(componenti.size()) +
                   ", шумовых компонент: " + std::to_string(noiseComponents.size()));
    }
}
