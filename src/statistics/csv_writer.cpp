#include "statistics/csv_writer.hpp"

#include <filesystem>

namespace statistics {

CsvWriter::CsvWriter(core::Logger& lg) : logger(lg) {}

CsvWriter::~CsvWriter() {
    logger.trace("CsvWriter: destructor called");
    close();
}

bool CsvWriter::open(const std::string& filename) {
    logger.info("CsvWriter: open requested for file = " + filename);
    
    if (file.is_open() && currentFile == filename) {
        logger.debug("CsvWriter: file already open and matches current file");
        return true;
    }

    close();

    currentFile = filename;

    bool fileExists = std::filesystem::exists(filename);
    logger.debug("CsvWriter: file exists = " + std::string(fileExists ? "true" : "false"));
    
    file.open(filename, std::ios::out | std::ios::app);
    if (!file.is_open()) {
        logger.error("CsvWriter: ERROR - failed to open file: " + filename);
        return false;
    }

    headerWritten = fileExists && std::filesystem::file_size(filename) > 0;
    
    logger.debug("CsvWriter: file opened successfully, headerWritten = " +
                 std::string(headerWritten ? "true" : "false"));
                 
    return true;
}

void CsvWriter::writeHeader() {
    logger.debug("CsvWriter: writeHeader called");
    
    if (headerWritten) {
        logger.debug("CsvWriter: header already written, skipping");
        return;
    }
    
    file << "environment,algorithm,executionTimeMs,pathNodes,expandedNodes,euclideanLength,pixelLength,pathFound\n";
    headerWritten = true;
    
    logger.info("CsvWriter: header written");
}

void CsvWriter::writeMetrics(const algorithms::path::PathMetrics& metrics,
                             const std::string& filename) {
    if (!open(filename)) {
        logger.error("CsvWriter: ERROR - open() failed, metrics not written");
        return;
    }
    
    writeHeader();

    file
        << metrics.environment << ","
        << metrics.algorithmName << ","
        << metrics.executionTimeMs << ","
        << metrics.pathNodes << ","
        << metrics.expandedNodes << ","
        << metrics.euclideanLength << ","
        << metrics.pixelLength << ","
        << metrics.pathFound
        << "\n";

    file.flush();
    logger.info("CsvWriter: metrics written successfully");
}

void CsvWriter::close() {
     logger.debug("CsvWriter: close called");
     
     if (file.is_open()) {
         file.close();
         logger.debug("CsvWriter: file closed");
     } else {
         logger.debug("CsvWriter: close called but file was not open");
     }

    currentFile.clear();
    headerWritten = false;
}

} // namespace statistics
