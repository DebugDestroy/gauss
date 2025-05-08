#pragma once
#include <vector>      // Для std::vector
#include <string>      // Для std::string
#include <stdexcept>   // Для std::invalid_argument

class Pole {
private:
    Logger& logger;

    void logResizeOperation(int old_rows, int old_cols, int new_rows, int new_cols) const {
        logger.debug(std::string("[Pole] Resizing field from ") +
                    std::to_string(old_rows) + std::string("x") + std::to_string(old_cols) + 
                    std::string(" to ") + std::to_string(new_rows) + std::string("x") + std::to_string(new_cols));
    }

public:
    std::vector<std::vector<double>> field;

    Pole(int A, int B, Logger& lg) : logger(lg) {
        logger.trace(std::string("[Pole] Constructor called with dimensions ") +
                   std::to_string(A) + std::string("x") + std::to_string(B));
        
        if (A <= 0 || B <= 0) {
            logger.error(std::string("[Pole] Invalid dimensions: ") +
                       std::to_string(A) + std::string("x") + std::to_string(B));
            throw std::invalid_argument("Field dimensions must be positive");
        }

        
        field.resize(A, std::vector<double>(B, 0));
        
        logger.info(std::string("[Pole] Field initialized successfully. Dimensions: ") +
                   std::to_string(field.size()) + std::string("x") + 
                   std::to_string(field.empty() ? 0 : field[0].size()));
    }
    
    void resize(int A, int B) {
        logger.trace("[Pole::resize] Starting resize operation");
        
        if (A <= 0 || B <= 0) {
            logger.error(std::string("[Pole::resize] Invalid target dimensions: ") +
                       std::to_string(A) + std::string("x") + std::to_string(B));
            throw std::invalid_argument("Field dimensions must be positive");
        }

        const size_t old_rows = field.size();
        const size_t old_cols = old_rows > 0 ? field[0].size() : 0;
        
        logResizeOperation(old_rows, old_cols, A, B);
        
        // Сохраняем данные при уменьшении размера
        if (A < static_cast<int>(old_rows)) {
            logger.debug(std::string("[Pole::resize] Truncating rows from ") + 
                       std::to_string(old_rows) + std::string(" to ") + std::to_string(A));
        }
        
        field.resize(A);
        
        for (auto& row : field) {
            if (B < static_cast<int>(row.size())) {
                logger.debug(std::string("[Pole::resize] Truncating columns from ") +
                           std::to_string(row.size()) + std::string(" to ") + std::to_string(B));
            }
            row.resize(B, 0);
        }

        logger.info(std::string("[Pole::resize] Resize completed. New dimensions: ") +
                   std::to_string(field.size()) + std::string("x") + 
                   std::to_string(field.empty() ? 0 : field[0].size()));
    }

    void logFieldStats() {
        logger.debug(std::string("[Pole] Current field statistics:\n") +
             std::string("  Dimensions: ") + std::to_string(field.size()) + "x" + 
             std::to_string(field.empty() ? 0 : field[0].size()) + "\n" +
             std::string("  Memory usage: ~") + 
             std::to_string(field.size() * (field.empty() ? 0 : field[0].size()) * sizeof(double)) + 
             " bytes");
    }
};
