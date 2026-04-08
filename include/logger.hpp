#ifndef CAPDAT_LOGGER_HPP
#define CAPDAT_LOGGER_HPP

#include <fstream>
#include <string>

/**
 * @brief Severity levels supported by the CapDAT logger.
 *
 * The numeric ordering is chosen so that lower-severity-detail levels can be
 * filtered using a simple threshold policy:
 *
 * - ERROR   : fatal or near-fatal issues
 * - WARNING : recoverable issues
 * - INFO    : normal runtime progress messages
 * - DEBUG   : verbose diagnostic information
 */
enum class LogLevel {
    ERROR = 0,
    WARNING = 1,
    INFO = 2,
    DEBUG = 3
};

/**
 * @brief Lightweight logger for terminal and optional file output.
 *
 * In v01, the logger is intended to support:
 * - concise user-facing runtime messages,
 * - optional traceability through a log file,
 * - simple verbosity control,
 * - consistent message formatting.
 *
 * The design is intentionally small and dependency-free. More advanced features
 * such as structured logging, sinks, thread safety, or log rotation are out of
 * scope for the foundation release.
 */
class Logger {
public:
    /**
     * @brief Construct a logger with default INFO verbosity and terminal output enabled.
     */
    Logger();

    /**
     * @brief Set the active terminal verbosity threshold.
     *
     * Messages with severity less verbose than the threshold may be suppressed
     * from terminal output depending on the implementation policy.
     *
     * @param level Desired verbosity threshold.
     */
    void setVerbosity(LogLevel level);

    /**
     * @brief Return the current terminal verbosity threshold.
     *
     * @return Current verbosity level.
     */
    [[nodiscard]] LogLevel verbosity() const;

    /**
     * @brief Enable logging to a file.
     *
     * @param path Output log file path.
     *
     * The implementation should open the file and keep it ready for subsequent
     * log writes. Failure to open should be handled clearly by the .cpp file.
     */
    void setLogFile(const std::string& path);

    /**
     * @brief Return whether file logging is currently enabled and open.
     *
     * @return True if a log file stream is active.
     */
    [[nodiscard]] bool hasLogFile() const;

    /**
     * @brief Write an INFO-level message.
     *
     * @param message Text to log.
     */
    void info(const std::string& message);

    /**
     * @brief Write a WARNING-level message.
     *
     * @param message Text to log.
     */
    void warning(const std::string& message);

    /**
     * @brief Write an ERROR-level message.
     *
     * @param message Text to log.
     */
    void error(const std::string& message);

    /**
     * @brief Write a DEBUG-level message.
     *
     * @param message Text to log.
     */
    void debug(const std::string& message);

private:
    /**
     * @brief Core internal logging function.
     *
     * @param level Severity level of the message.
     * @param message Message text.
     */
    void log(LogLevel level, const std::string& message);

    /**
     * @brief Convert a LogLevel to its printable label.
     *
     * Expected examples include "INFO", "WARNING", etc.
     *
     * @param level Severity level.
     * @return Printable level label.
     */
    [[nodiscard]] std::string levelToString(LogLevel level) const;

    /**
     * @brief Return a formatted timestamp string for log messages.
     *
     * In v01, this should remain lightweight and human-readable.
     *
     * @return Timestamp string.
     */
    [[nodiscard]] std::string currentTimestamp() const;

private:
    // -------------------------------------------------------------------------
    // Logger state
    // -------------------------------------------------------------------------

    LogLevel verbosity_ = LogLevel::INFO;
    std::ofstream log_file_;
};

// NOTE ON LOGGER SCOPE:
//
// In v01, Logger is intentionally small and concrete. The project description
// asks for traceable runtime feedback and optional log-file output, but does
// not require a large logging framework.
//
// A future revision could introduce richer features such as thread-safe
// logging, structured fields, sink abstraction, or message categories if batch
// workflows or parallel analyses later justify that complexity. For the
// foundation release, a simple concrete logger is the most maintainable choice.

// NOTE ON LOGLEVEL ORDERING:
//
// The enum ordering is chosen to make threshold-based filtering easy to
// implement. This is a practical convention rather than a deep architectural
// commitment.
//
// If later logging behavior becomes more sophisticated, the filtering policy
// could be revised without changing the public API very much.

// NOTE ON <fstream> IN THE HEADER:
//
// The header stores std::ofstream directly as a member because Logger owns the
// file stream lifecycle in v01. This keeps ownership explicit and avoids
// introducing extra indirection too early.
//
// If compile-time dependency reduction ever becomes important, the stream could
// later be hidden behind a pimpl or another internal abstraction. For v01, the
// direct approach is clearer and simpler.

#endif // CAPDAT_LOGGER_HPP
