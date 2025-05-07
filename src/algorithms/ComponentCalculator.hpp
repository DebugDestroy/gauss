#pragma once
#include <vector>
#include <algorithm> // для std::remove_if, std::find, std::all_of
#include <cmath> // для std::isnan, std::sqrt
#include <numeric> // для std::iota
#include <string>
#include <sstream>
#include <limits> // для std::numeric_limits

// Локальные заголовки
#include "core/Constants.hpp"  // Подключаем константы
#include "core/Logger.hpp"
#include "services/Pole.hpp" // для std::unique_ptr<Pole>
#include "algorithms/Component.hpp" // для Component, ThresholdMode

class ComponentCalculator {
private:
    Logger& logger;

    int incrementAndCollect(std::vector<std::vector<double>>& componenta, 
                          std::vector<std::vector<double>>& CopyPole, 
                          int x, int y, int i, int& pixelCount) {
        logger.trace("[ComponentCalculator::incrementAndCollect] Вход: "
                "координаты=(" + std::to_string(x) + "," + std::to_string(y) + "), "
                "глубина_рекурсии=" + std::to_string(i) + ", "
                "текущий_счетчик=" + std::to_string(pixelCount));

        if (x < 1 || y < 1 || x > (int)componenta[0].size() - 2 || 
            y > (int)componenta.size() - 2 || std::fabs(CopyPole[y][x] - Constants::WHITE) > Constants::EPSILON) {
            logger.trace("[ComponentCalculator::incrementAndCollect] Выход за границы или неподходящее значение");
            return pixelCount;
        }

        if (std::fabs(CopyPole[y][x] - Constants::WHITE) < Constants::EPSILON) {
            CopyPole[y][x] = Constants::BLACK;
            pixelCount++;
            componenta[y][x] = Constants::WHITE;
            
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
public:
    ComponentCalculator(Logger& lg) : logger(lg) {
        logger.trace("[ComponentCalculator] Инициализация калькулятора компонент");
    }

    void bin(std::vector<std::vector<double>>& CopyPole, 
            int slice, 
            std::unique_ptr<Pole>& p, 
            ThresholdMode mode) {
        logger.info(std::string("[ComponentCalculator::bin] Бинаризация данных, slice=") + 
                  std::to_string(slice) + ", mode=" + 
                  (mode == ThresholdMode::Peaks ? "Peaks" : 
                   mode == ThresholdMode::Valleys ? "Valleys" : "All"));

        if (p == nullptr) {
            logger.error("[ComponentCalculator::bin] Ошибка: данные высот не инициализированы!");
            return;
        }
        
        CopyPole = p->field;
        int symmetric_slice = 2*Constants::MID_GRAY - slice;

        for (int x = 0; x < (int)p->field[0].size(); ++x) {
            for (int y = 0; y < (int)p->field.size(); ++y) {
                double value = p->field[y][x];
                switch (mode) {
                    case ThresholdMode::All: {
                        bool is_peak = (slice >= Constants::MID_GRAY) ? (value > slice) : (value > symmetric_slice);
                        bool is_valley = (slice >= Constants::MID_GRAY) ? (value < symmetric_slice) : (value < slice);
                        CopyPole[y][x] = (is_peak || is_valley) ? Constants::WHITE : Constants::BLACK;
                        break;
                    }
                    case ThresholdMode::Peaks:
                        CopyPole[y][x] = (slice >= Constants::MID_GRAY) ? (value > slice) : (value > symmetric_slice);
                        break;
                    case ThresholdMode::Valleys:
                        CopyPole[y][x] = (slice >= Constants::MID_GRAY) ? (value < symmetric_slice) : (value < slice);
                        break;
                }
            } 
        }

        logger.info("[ComponentCalculator::bin] Бинаризация завершена");
    }
    
    void wave(int noisy, 
             std::vector<Component>& componenti, 
             std::vector<std::vector<double>>& CopyPole, 
             std::unique_ptr<Pole>& p) {
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
                if (std::fabs(CopyPole[y][x] - Constants::WHITE) < Constants::EPSILON) {
                    std::vector<std::vector<double>> componentData(rows, std::vector<double>(cols, Constants::BLACK));
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
                                if (std::fabs(componentData[i][j] - Constants::WHITE) < Constants::EPSILON) {
                                    p->field[i][j] = Constants::MID_GRAY;
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
};
