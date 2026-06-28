#include <vector>
#include <algorithm> // для std::remove_if, std::find, std::all_of
#include <cmath> // для std::isnan, std::sqrt
#include <numeric> // для std::iota
#include <string>
#include <sstream>
#include <stack>
#include <limits> // для std::numeric_limits

// Локальные заголовки
#include "algorithms/components/find_components.hpp"

namespace algorithms::components {

void ComponentCalculator::collectComponent(
        const std::vector<std::vector<double>>& binaryMap,
        std::vector<std::vector<bool>>& visited,
        int startX,
        int startY,
        std::vector<algorithms::geometry::Pixel>& pixels)
{
    std::stack<algorithms::geometry::Pixel> stack;

    stack.push({startX, startY});

    while (!stack.empty()) {

        algorithms::geometry::Pixel current = stack.top();
        stack.pop();

        int x = current.x;
        int y = current.y;

        if (x < 0 || x >= static_cast<int>(binaryMap[0].size()) ||
            y < 0 || y >= static_cast<int>(binaryMap.size()))
        {
            continue;
        }

        if (visited[y][x]) {
            continue;
        }

        if (std::fabs(binaryMap[y][x] - core::WHITE) > core::EPSILON) {
            continue;
        }

        visited[y][x] = true;

        pixels.push_back({x, y});

        stack.push({x + 1, y});
        stack.push({x - 1, y});
        stack.push({x, y + 1});
        stack.push({x, y - 1});
    }
}

    ComponentCalculator::ComponentCalculator(core::Logger& lg) : logger(lg) {
        logger.trace("[ComponentCalculator] Инициализация калькулятора компонент");
    }
    
void ComponentCalculator::wave(
        int noisy,
        std::vector<Component>& componenti,
        std::vector<std::vector<double>>& binaryMap,
        std::vector<std::vector<double>>& field)
{
    logger.info("[ComponentCalculator::wave] Начало волнового алгоритма, порог=" +
                std::to_string(noisy));

    if (field.empty()) {
        logger.error("[ComponentCalculator::wave] Поле не инициализировано");
        return;
    }

    componenti.clear();

    const int rows = field.size();
    const int cols = field[0].size();

    std::vector<std::vector<bool>> visited(
        rows,
        std::vector<bool>(cols, false));

    for (int y = 0; y < rows; ++y) {
        for (int x = 0; x < cols; ++x) {

            if (visited[y][x]) {
                continue;
            }

            if (std::fabs(binaryMap[y][x] - core::WHITE) > core::EPSILON) {
                continue;
            }

            std::vector<algorithms::geometry::Pixel> componentPixels;

            collectComponent(binaryMap,
                             visited,
                             x,
                             y,
                             componentPixels);

            const int pixelCount =
                static_cast<int>(componentPixels.size());

            if (pixelCount >= noisy) {

                componenti.emplace_back(logger,
                                        componentPixels);

            } else {

                // удаляем шум из поля
                for (const auto& p : componentPixels) {

                    int px = p.x;
                    int py = p.y;

                    field[py][px] = core::MID_GRAY;
                    binaryMap[py][px] = core::BLACK;
                }
            }
        }
    }

    logger.info("[ComponentCalculator::wave] Обработка завершена");

    logger.debug("Найдено значимых компонент: " +
                 std::to_string(componenti.size()));
}
}
