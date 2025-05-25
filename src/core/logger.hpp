#pragma once

#include <string>
#include <fstream>

namespace core {

enum class LogLevel {
    Trace,
    Debug,
    Info,
    Warning,
    Error,
    Critical,
    Off
};

class Logger {
private:
    std::string moduleName;
    LogLevel currentLogLevel;
    std::ofstream logFile;

    static const char* levelToString(LogLevel level);
    static LogLevel stringToLevel(const std::string& levelStr);

public:
    Logger(const std::string& fileName, const std::string& module, const std::string& configLevel);
    ~Logger();

    LogLevel getLogLevel() const;
    
    void logMessage(LogLevel level, const std::string& message);
    void trace(const std::string& message);
    void debug(const std::string& message);
    void info(const std::string& message);
    void warning(const std::string& message);
    void error(const std::string& message);
    void critical(const std::string& message);
};

}
