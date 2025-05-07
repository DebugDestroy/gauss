#pragma once
#include <vector>
#include <algorithm>
#include <string>

// Локальные заголовки
#include "core/Constants.hpp"  // Подключаем константы
#include "core/Logger.hpp"
#include "algorithms/Component.hpp"

class Copier {
private:
    Logger& logger;

public:
    Copier(Logger& lg) : logger(lg) {}

    void removeNoise(std::vector<std::vector<double>>& field, const std::vector<Component>& components) {
        logger.trace("[Copier::removeNoise] Starting noise removal");
        
        // Логируем исходные параметры
             logger.debug(std::string("[Copier::removeNoise] Input parameters: ") +
             std::to_string(field.size()) + "x" + 
             std::to_string(field.empty() ? 0 : field[0].size()) + " field, " +
             std::to_string(components.size()) + " components");
             
        // Обнуляем поле
        size_t total_zeroed = 0;
        for (auto& row : field) {
            total_zeroed += row.size();
            std::fill(row.begin(), row.end(), Constants::BLACK);
        }
        logger.debug(std::string("[Copier::removeNoise] Zeroed ") + std::to_string(total_zeroed) + std::string(" pixels"));

        // Копируем только значимые компоненты
        size_t total_copied = 0;
        size_t components_processed = 0;
        
        for (const auto& comp : components) {
            const auto& compData = comp.componenta;
            size_t component_pixels = 0;
            
            for (size_t i = 0; i < compData.size(); ++i) {
                for (size_t j = 0; j < compData[i].size(); ++j) {
                    if (std::fabs(compData[i][j] - Constants::WHITE) < Constants::EPSILON) {
                        field[i][j] = Constants::WHITE;
                        component_pixels++;
                        total_copied++;
                    }
                }
            }
            
            components_processed++;
            logger.debug(std::string("[Copier::removeNoise] Component #") + std::to_string(components_processed) + 
                       std::string(": copied ") + std::to_string(component_pixels) + std::string(" pixels"));
        }

        // Итоговая статистика
        logger.info(std::string("[Copier::removeNoise] Completed. Stats:\n") +
                  std::string("  Total components processed: ") + std::to_string(components_processed) + std::string("\n") +
                  std::string("  Total pixels copied: ") + std::to_string(total_copied) + std::string("\n") +
                  std::string("  Zeroed pixels: ") + std::to_string(total_zeroed - total_copied));
    }
};
