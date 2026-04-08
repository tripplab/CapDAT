#include "logger.hpp"

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>

/**
 * @brief Construct a logger with default INFO verbosity.
 *
 * Terminal output is enabled by default. File logging remains disabled until
 * setLogFile() is called successfully.
 */
Logger::Logger()
    : verbosity_(LogLevel::INFO) {
}

/**
 * @brief Set the active terminal verbosity threshold.
 */
void Logger::setVerbosity(LogLevel level) {
    verbosity_ = level;
}

/**
 * @brief Return the current terminal verbosity threshold.
 */
LogLevel Logger::verbosity() const {
    return verbosity_;
}

/**
 * @brief Enable logging to a file.
 *
 * The file is opened in output mode, replacing previous contents if the file
 * already exists.
 *
 * @throws std::runtime_error if the file cannot be opened.
 */
void Logger::setLogFile(const std::string& path) {
    log_file_.open(path, std::ios::out);

    if (!log_file_.is_open()) {
        throw std::runtime_error("Failed to open log file: " + path);
    }
}

/**
 * @brief Return whether file logging is currently enabled and open.
 */
bool Logger::hasLogFile() const {
    return log_file_.is_open();
}

/**
 * @brief Write an INFO-level message.
 */
void Logger::info(const std::string& message) {
    log(LogLevel::INFO, message);
}

/**
 * @brief Write a WARNING-level message.
 */
void Logger::warning(const std::string& message) {
    log(LogLevel::WARNING, message);
}

/**
 * @brief Write an ERROR-level message.
 */
void Logger::error(const std::string& message) {
    log(LogLevel::ERROR, message);
}

/**
 * @brief Write a DEBUG-level message.
 */
void Logger::debug(const std::string& message) {
    log(LogLevel::DEBUG, message);
}

/**
 * @brief Core internal logging function.
 *
 * The implementation uses a simple threshold rule for terminal output:
 * a message is printed if its severity is less than or equal to the configured
 * verbosity threshold according to the enum ordering.
 *
 * File output, when enabled, records all messages regardless of terminal
 * filtering.
 */
void Logger::log(LogLevel level, const std::string& message) {
    const std::string formatted =
        "[" + currentTimestamp() + "] "
        "[" + levelToString(level) + "] "
        + message;

    if (static_cast<int>(level) <= static_cast<int>(verbosity_)) {
        if (level == LogLevel::ERROR || level == LogLevel::WARNING) {
            std::cerr << formatted << '\n';
        } else {
            std::cout << formatted << '\n';
        }
    }

    if (log_file_.is_open()) {
        log_file_ << formatted << '\n';
        log_file_.flush();
    }
}

/**
 * @brief Convert a LogLevel to a printable label.
 */
std::string Logger::levelToString(LogLevel level) const {
    switch (level) {
        case LogLevel::ERROR:
            return "ERROR";
        case LogLevel::WARNING:
            return "WARNING";
        case LogLevel::INFO:
            return "INFO";
        case LogLevel::DEBUG:
            return "DEBUG";
        default:
            return "UNKNOWN";
    }
}

/**
 * @brief Return a human-readable timestamp string for log messages.
 *
 * The format is:
 * YYYY-MM-DD HH:MM:SS
 */
std::string Logger::currentTimestamp() const {
    const auto now = std::chrono::system_clock::now();
    const std::time_t now_c = std::chrono::system_clock::to_time_t(now);

    std::tm tm_snapshot{};
#if defined(_WIN32)
    localtime_s(&tm_snapshot, &now_c);
#else
    localtime_r(&now_c, &tm_snapshot);
#endif

    std::ostringstream oss;
    oss << std::put_time(&tm_snapshot, "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

// NOTE ON VERBOSITY FILTERING:
//
// The logger uses a simple numeric-threshold policy based on the ordering of
// LogLevel. This keeps v01 easy to understand and avoids unnecessary logging
// complexity in the foundation release.
//
// If later workflows need category-based filtering, multiple sinks, or more
// nuanced verbosity behavior, the filtering policy can be expanded without
// radically changing the public API.

// NOTE ON FILE FLUSHING:
//
// The implementation flushes the file stream after each message. This is a
// deliberate reliability-first choice for v01 because it improves traceability
// in the event of crashes or interrupted runs.
//
// The tradeoff is a bit more I/O overhead. If future profiling shows that log
// flushing materially affects throughput in larger workflows, we could relax
// this policy or make it configurable.

// NOTE ON TIMESTAMP SOURCE:
//
// currentTimestamp() uses std::chrono::system_clock because the goal here is
// human-readable logging, not monotonic elapsed-time measurement. This is
// appropriate for log records, whereas runtime duration measurement belongs in
// Timer and should rely on steady_clock instead.
