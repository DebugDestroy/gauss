#include "core/logger.hpp"
#include <iostream>
#include <chrono>
#include <ctime>
#include <iomanip>

namespace core {

Logger::Logger(const std::string& fileName, const std::string& module, const std::string& configLevel)
    : moduleName(module), currentLogLevel(stringToLevel(configLevel)) {
    if (currentLogLevel != LogLevel::Off) {
        logFile.open(fileName, std::ios::out | std::ios::app);
        if (!logFile.is_open()) {
            std::cerr << "Failed to open log file: " << fileName << std::endl;
        }
    }
}

Logger::~Logger() {
    if (logFile.is_open()) {
        logFile.close();
    }
}

LogLevel Logger::stringToLevel(const std::string& levelStr) {
    if (levelStr == "TRACE")    return LogLevel::Trace;
    if (levelStr == "DEBUG")    return LogLevel::Debug;
    if (levelStr == "INFO")     return LogLevel::Info;
    if (levelStr == "WARNING")  return LogLevel::Warning;
    if (levelStr == "ERROR")    return LogLevel::Error;
    if (levelStr == "CRITICAL") return LogLevel::Critical;
    if (levelStr == "OFF")      return LogLevel::Off;
    return LogLevel::Info;
}

const char* Logger::levelToString(LogLevel level) {
    switch(level) {
        case LogLevel::Trace:    return "TRACE";
        case LogLevel::Debug:    return "DEBUG";
        case LogLevel::Info:     return "INFO";
        case LogLevel::Warning:  return "WARNING";
        case LogLevel::Error:    return "ERROR";
        case LogLevel::Critical: return "CRITICAL";
        default:                 return "UNKNOWN";
    }
}

void Logger::logMessage(LogLevel level, const std::string& message) {
    if (level < currentLogLevel || !logFile.is_open()) return;

    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    logFile << "[" << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << "] "
            << "[" << levelToString(level) << "] "
            << "[" << moduleName << "] "
            << message << std::endl;
}

void Logger::trace(const std::string& message)    { logMessage(LogLevel::Trace, message); }
void Logger::debug(const std::string& message)    { logMessage(LogLevel::Debug, message); }
void Logger::info(const std::string& message)     { logMessage(LogLevel::Info, message); }
void Logger::warning(const std::string& message)  { logMessage(LogLevel::Warning, message); }
void Logger::error(const std::string& message)    { logMessage(LogLevel::Error, message); }
void Logger::critical(const std::string& message) { logMessage(LogLevel::Critical, message); }

}
