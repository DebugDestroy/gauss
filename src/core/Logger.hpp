#pragma once
#include <string>
#include <fstream>
#include <ctime>
#include <chrono>
#include <iomanip>

// Уровни логирования от самого детального до самого важного
enum class LogLevel {
    Trace,    // Максимальная детализация (трассировка выполнения)
    Debug,    // Отладочная информация
    Info,     // Информационные сообщения
    Warning,  // Предупреждения
    Error,    // Ошибки (не в алгоритме)
    Critical, // Критические ошибки (в алгоритме)
    Off       // Логирование отключено
};

class Logger {
private:
    std::string moduleName;
    LogLevel currentLogLevel;
    std::ofstream logFile;
    
    // Преобразуем уровень логирования в строку
    static const char* levelToString(LogLevel level) {
        switch(level) {
            case LogLevel::Trace:    return "TRACE";
            case LogLevel::Debug:   return "DEBUG";
            case LogLevel::Info:    return "INFO";
            case LogLevel::Warning: return "WARNING";
            case LogLevel::Error:   return "ERROR";
            case LogLevel::Critical: return "CRITICAL";
            default:                return "UNKNOWN";
        }
    }
    
    // Преобразуем строку из конфига в LogLevel
    static LogLevel stringToLevel(const std::string& levelStr) {
        if (levelStr == "TRACE")    return LogLevel::Trace;
        if (levelStr == "DEBUG")    return LogLevel::Debug;
        if (levelStr == "INFO")     return LogLevel::Info;
        if (levelStr == "WARNING")  return LogLevel::Warning;
        if (levelStr == "ERROR")    return LogLevel::Error;
        if (levelStr == "CRITICAL") return LogLevel::Critical;
        if (levelStr == "OFF")      return LogLevel::Off;
        
        // По умолчанию INFO, если неизвестный уровень
        return LogLevel::Info;
    }

public:
    Logger(const std::string& fileName, const std::string& module, const std::string& configLevel)
        : moduleName(module), currentLogLevel(stringToLevel(configLevel)) {
        if (currentLogLevel != LogLevel::Off) {
            logFile.open(fileName, std::ios::out | std::ios::app);
            if (!logFile.is_open()) {
                std::cerr << "Failed to open log file: " << fileName << std::endl;
            }
        }
    }

    ~Logger() {
        if (logFile.is_open()) {
            logFile.close();
        }
    }
    
    LogLevel getLogLevel() const {
        return currentLogLevel;
    }
    
    void logMessage(LogLevel level, const std::string& message) {
        // Проверяем, нужно ли логировать это сообщение
        if (level < currentLogLevel || !logFile.is_open()) {
            return;
        }
        
        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        
        logFile << "[" << std::put_time(std::localtime(&now_c), "%Y-%m-%d %H:%M:%S") << "] "
                << "[" << levelToString(level) << "] "
                << "[" << moduleName << "] "
                << message << std::endl;
    }
    
    // Удобные методы для каждого уровня
    void trace(const std::string& message) {
        logMessage(LogLevel::Trace, message);
    }
    
    void debug(const std::string& message) {
        logMessage(LogLevel::Debug, message);
    }
    
    void info(const std::string& message) {
        logMessage(LogLevel::Info, message);
    }
    
    void warning(const std::string& message) {
        logMessage(LogLevel::Warning, message);
    }
    
    void error(const std::string& message) {
        logMessage(LogLevel::Error, message);
    }
    
    void critical(const std::string& message) {
        logMessage(LogLevel::Critical, message);
    }
};
