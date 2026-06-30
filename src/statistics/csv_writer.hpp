#pragma once

#include "core/logger.hpp"
#include "algorithms/path/common/path_metrics.hpp"

#include <fstream>
#include <string>

namespace statistics {

class CsvWriter {
public:
    CsvWriter(core::Logger& lg);
    ~CsvWriter();

    // Записать одну строку метрик
    void writeMetrics(const algorithms::path::PathMetrics& metrics,
                      const std::string& filename);

private:
    void writeHeader();
    bool open(const std::string& filename);
    void close();

private:
    core::Logger& logger;
    std::ofstream file;
    std::string currentFile;
    bool headerWritten = false;
};

} // namespace statistics
