#include "timer.hpp"

/**
 * @brief Start or restart the timer.
 *
 * Calling start() resets the start time and marks the timer as running.
 */
void Timer::start() {
    start_time_ = Clock::now();
    end_time_ = TimePoint{};
    running_ = true;
}

/**
 * @brief Stop the timer.
 *
 * If the timer is currently running, this records the stop time and marks the
 * timer as no longer running. If the timer was already stopped, this call is
 * harmless in v01.
 */
void Timer::stop() {
    if (running_) {
        end_time_ = Clock::now();
        running_ = false;
    }
}

/**
 * @brief Return the elapsed wall-clock time in seconds.
 *
 * Behavior:
 * - If the timer is currently running, elapsed time is measured from the
 *   stored start time to the current time.
 * - If the timer has been stopped, elapsed time is measured from the stored
 *   start time to the stored stop time.
 * - If the timer was never started, the elapsed time is reported as 0.0.
 *
 * @return Elapsed time in seconds.
 */
double Timer::elapsedSeconds() const {
    // If the timer was never started, return zero.
    if (start_time_ == TimePoint{}) {
        return 0.0;
    }

    const TimePoint end_point = running_ ? Clock::now() : end_time_;
    const std::chrono::duration<double> elapsed = end_point - start_time_;
    return elapsed.count();
}

/**
 * @brief Return whether the timer is currently running.
 */
bool Timer::isRunning() const {
    return running_;
}

// NOTE ON DEFAULT-INITIALIZED TIME POINTS:
//
// The implementation uses a default-constructed TimePoint as a simple sentinel
// to detect the "never started" state. This is acceptable in v01 because it
// keeps the class small and avoids introducing extra state just for startup
// bookkeeping.
//
// A future revision could replace this with an explicit "has_started_" flag if
// the timer API grows more complex or if clearer state semantics become useful.

// NOTE ON ELAPSED-TIME POLICY:
//
// elapsedSeconds() is defined to work both while the timer is running and after
// it has been stopped. This makes the class more convenient for lightweight
// runtime reporting and intermediate progress diagnostics.
//
// For the foundation release, this is more practical than forcing callers to
// manage separate querying modes.
