#pragma once
#include <vector>
#include <numeric> // для std::iota
#include <random> // для std::mt19937, std::random_device
#include <algorithm> // для std::shuffle
#include <limits> // для std::numeric_limits
#include <string>
#include <sstream>
#include <array>

// Локальные заголовки
#include "core/Logger.hpp"
#include "utils/ColorGenerator.hpp" // для ColorGenerator

class KMeans {
private:
    Logger& logger;

    void logCenters(const std::vector<std::vector<double>>& centers, const std::string& prefix = "") const {
        std::ostringstream oss;
        oss << prefix << "Центры кластеров (" << centers.size() << "):";
        for (size_t i = 0; i < centers.size(); ++i) {
            oss << "\n  Кластер " << i << ": (" << centers[i][0] << ", " << centers[i][1] << ")";
        }
        logger.debug(oss.str());
    }

    void logPoint(const std::vector<double>& point, const std::string& prefix = "") const {
        logger.trace(prefix + "Точка (" + std::to_string(point[0]) + ", " + std::to_string(point[1]) + ")");
    }

public:
    std::vector<int> labels;
    std::vector<std::vector<double>> centers; // Центры кластеров

    struct ClusterResult {
        std::vector<int> labels;
        std::vector<std::array<int, 3>> colors;
        std::vector<std::vector<double>> centers;
    };

    KMeans(Logger& lg) : logger(lg) {
        logger.trace("[KMeans] Инициализация алгоритма K-средних");
    }

    void initializeCenters(const std::vector<std::vector<double>>& data, int k) {
        logger.info("[KMeans::initializeCenters] Инициализация центров кластеров");
        logger.debug(std::string("Количество кластеров: ") + std::to_string(k) + 
                   ", количество точек: " + std::to_string(data.size()));

        if (data.empty() || k <= 0) {
            logger.error("[KMeans::initializeCenters] Ошибка: пустые данные или k <= 0");
            return;
        }
        
        centers.clear();
        std::vector<size_t> indices(data.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});
        
        for (int i = 0; i < k && i < static_cast<int>(data.size()); ++i) {
            centers.push_back(data[indices[i]]);
            logPoint(data[indices[i]], "Выбран начальный центр: ");
        }

        logCenters(centers, "Инициализированы центры: ");
    }

    double distance(const std::vector<double>& a, const std::vector<double>& b) const {
        double dist = std::hypot(a[0]-b[0], a[1]-b[1]);
        logger.trace(std::string("[KMeans::distance] Расстояние между (") + 
                   std::to_string(a[0]) + "," + std::to_string(a[1]) + ") и (" +
                   std::to_string(b[0]) + "," + std::to_string(b[1]) + ") = " +
                   std::to_string(dist));
        return dist;
    }

    ClusterResult cluster(const std::vector<std::vector<double>>& data, int k) {
        logger.info("[KMeans::cluster] Начало кластеризации");
        logger.debug(std::string("Количество кластеров: ") + std::to_string(k) + 
                   ", точек: " + std::to_string(data.size()));

        ClusterResult result;
        initializeCenters(data, k);
        bool changed;
        int iteration = 0;
        
        do {
            iteration++;
            logger.debug(std::string("[KMeans::cluster] Итерация ") + std::to_string(iteration));
            changed = false;
            labels.assign(data.size(), -1);

            // Назначение меток
            for (size_t i = 0; i < data.size(); ++i) {
                double minDist = std::numeric_limits<double>::max();
                for (size_t j = 0; j < centers.size(); ++j) {
                    double dist = distance(data[i], centers[j]);
                    if (dist < minDist) {
                        minDist = dist;
                        labels[i] = j;
                    }
                }
                logger.trace(std::string("[KMeans::cluster] Точка ") + std::to_string(i) + 
                            " отнесена к кластеру " + std::to_string(labels[i]));
            }

            // Обновление центров
            std::vector<std::vector<double>> newCenters(centers.size(), {0.0, 0.0});
            std::vector<int> counts(centers.size(), 0);

            for (size_t i = 0; i < data.size(); ++i) {
                int label = labels[i];
                if (label >= 0) {
                    newCenters[label][0] += data[i][0];
                    newCenters[label][1] += data[i][1];
                    counts[label]++;
                }
            }

            // Проверка изменений
            for (size_t j = 0; j < centers.size(); ++j) {
                if (counts[j] > 0) {
                    newCenters[j][0] /= counts[j];
                    newCenters[j][1] /= counts[j];
                    double centerDist = distance(newCenters[j], centers[j]);
                    if (centerDist > 1e-6) {
                        changed = true;
                    }
                    logger.debug(std::string("[KMeans::cluster] Кластер ") + std::to_string(j) + 
                               ": точек=" + std::to_string(counts[j]) + 
                               ", смещение центра=" + std::to_string(centerDist));
                } else {
                    logger.warning(std::string("[KMeans::cluster] Кластер ") + std::to_string(j) + 
                                 " остался без точек!");
                }
            }
            
            centers = newCenters;
            logCenters(centers, "Обновленные центры: ");

        } while (changed);

        logger.info(std::string("[KMeans::cluster] Кластеризация завершена за ") + 
                   std::to_string(iteration) + " итераций");

        result.labels = labels;
        result.centers = centers;
        result.colors = ColorGenerator::generateColors(k, logger);
        
        logger.debug(std::string("[KMeans::cluster] Сгенерировано ") + 
                    std::to_string(result.colors.size()) + " цветов для кластеров");
        logCenters(result.centers, "Финальные центры кластеров: ");
        
        return result;
    }

    ClusterResult kmeansWithKernels(const std::vector<std::vector<double>>& data, int k, int kernelSize) {
        logger.info("[KMeans::kmeansWithKernels] Начало кластеризации с ядрами");
        logger.debug(std::string("Параметры: k=") + std::to_string(k) + 
                   ", kernelSize=" + std::to_string(kernelSize) + 
                   ", точек=" + std::to_string(data.size()));

        ClusterResult result;
        if (data.empty() || k <= 0 || kernelSize <= 0) {
            logger.error("[KMeans::kmeansWithKernels] Ошибка: пустые данные или неверные параметры");
            return result;
        }

        // Шаг 1: Выборка ядерных центров
        std::vector<std::vector<double>> kernelCenters;
        std::vector<size_t> indices(data.size());
        std::iota(indices.begin(), indices.end(), 0);
        
        size_t samplesToTake = std::min(k * kernelSize, static_cast<int>(data.size()));
        if (samplesToTake == 0) {
            logger.error("[KMeans::kmeansWithKernels] Ошибка: недостаточно точек для выборки");
            return result;
        }

        std::shuffle(indices.begin(), indices.end(), std::mt19937{std::random_device{}()});
        for (size_t i = 0; i < samplesToTake; ++i) {
            kernelCenters.push_back(data[indices[i]]);
        }

        logger.debug(std::string("[KMeans::kmeansWithKernels] Выбрано ") + 
                    std::to_string(kernelCenters.size()) + " ядерных точек");

        // Шаг 2: Первичная кластеризация
        initializeCenters(kernelCenters, k * kernelSize);
        auto kernelResult = cluster(kernelCenters, k * kernelSize);

        // Шаг 3: Агрегация центров
        std::vector<std::vector<double>> aggregatedCenters(k, {0.0, 0.0});
        std::vector<int> counts(k, 0);

        for (size_t i = 0; i < kernelResult.labels.size(); ++i) {
            int label = kernelResult.labels[i];
            if (label >= 0 && label < k) {
                aggregatedCenters[label][0] += kernelCenters[i][0];
                aggregatedCenters[label][1] += kernelCenters[i][1];
                counts[label]++;
            }
        }

        // Шаг 4: Нормализация или замена пустых кластеров
        for (int i = 0; i < k; ++i) {
            if (counts[i] > 0) {
                aggregatedCenters[i][0] /= counts[i];
                aggregatedCenters[i][1] /= counts[i];
                logger.debug(std::string("[KMeans::kmeansWithKernels] Агрегированный центр ") + 
                            std::to_string(i) + " на основе " + 
                            std::to_string(counts[i]) + " точек");
            } else {
                aggregatedCenters[i] = data[std::rand() % data.size()];
                logger.warning(std::string("[KMeans::kmeansWithKernels] Кластер ") + 
                             std::to_string(i) + " был пуст, заменен случайной точкой");
            }
        }

        // Финальная кластеризация
        centers = aggregatedCenters;
        result = cluster(data, k);
        if (static_cast<int>(result.centers.size()) != k) {
            logger.error("[KMeans::kmeansWithKernels] Несоответствие количества кластеров!");
        }
        result.colors = ColorGenerator::generateColors(k, logger);

        logger.info("[KMeans::kmeansWithKernels] Кластеризация с ядрами завершена");
        logCenters(result.centers, "Финальные центры: ");
        
        return result;
    }
};
