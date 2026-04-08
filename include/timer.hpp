#ifndef CAPDAT_TIMER_HPP
#define CAPDAT_TIMER_HPP

#include <chrono>

/**
 * @brief Lightweight wall-clock timer for runtime measurement.
 *
 * In v01, Timer is intended to measure coarse execution intervals such as:
 * - total program runtime,
 * - time spent reading input,
 * - time spent parsing and hierarchy construction.
 *
 * The class is intentionally simple and dependency-free. It is not meant to be
 * a profiling framework, and it does not attempt to measure CPU time, thread
 * time, or fine-grained instrumentation beyond basic elapsed wall-clock time.
 */
class Timer {
public:
    /**
     * @brief Default constructor.
     *
     * The timer starts in a non-running state.
     */
    Timer() = default;

    /**
     * @brief Start or restart the timer.
     *
     * Calling start() resets the start point and marks the timer as running.
     */
    void start();

    /**
     * @brief Stop the timer.
     *
     * Calling stop() records the stop point. Repeated calls should be harmless
     * in v01, though the precise behavior is implementation-defined.
     */
    void stop();

    /**
     * @brief Return the elapsed wall-clock time in seconds.
     *
     * If the timer is still running, the implementation may either compute
     * elapsed time up to "now" or use the last recorded stop point if already
     * stopped. The .cpp implementation should define this clearly.
     *
     * @return Elapsed time in seconds.
     */
    [[nodiscard]] double elapsedSeconds() const;

    /**
     * @brief Return whether the timer is currently running.
     *
     * @return True if start() has been called and stop() has not yet finalized
     * the running interval.
     */
    [[nodiscard]] bool isRunning() const;

private:
    using Clock = std::chrono::steady_clock;
    using TimePoint = Clock::time_point;

private:
    // -------------------------------------------------------------------------
    // Timer state
    // -------------------------------------------------------------------------

    TimePoint start_time_{};
    TimePoint end_time_{};
    bool running_ = false;
};

// NOTE ON CLOCK CHOICE:
//
// std::chrono::steady_clock is used because it is monotonic and therefore more
// appropriate for elapsed-time measurement than a wall-clock source that may
// jump due to system-time adjustments.
//
// For v01, this is the most sensible default for reporting runtime in a stable
// and predictable way.

// NOTE ON TIMER SCOPE:
//
// The Timer class is intentionally minimal. The project description asks for a
// lightweight utility to report elapsed runtime, not a full profiler.
//
// A future revision could add scoped timers, named timers, nested timing, or
// aggregate profiling support if later analytical modules become more complex.
// For the foundation release, a small explicit timer is easier to audit and
// sufficient for the required reporting.

#endif // CAPDAT_TIMER_HPP
