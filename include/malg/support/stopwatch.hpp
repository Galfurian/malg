/// @file stopwatch.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Stopwatch for taking elapsed times.

#pragma once

#include <iostream>
#include <utility>
#include <iomanip>
#include <chrono>

class StopWatch {
private:
    std::chrono::high_resolution_clock::time_point _last_time_point;
    std::chrono::duration<double> _partial_duration;
    std::chrono::duration<double> _total_duration;
    unsigned _laps;

public:
    StopWatch()
        : _last_time_point(std::chrono::high_resolution_clock::now()),
          _partial_duration(std::chrono::duration<double>::zero()),
          _laps(0)
    {
    }

    inline void reset()
    {
        _last_time_point  = std::chrono::high_resolution_clock::now();
        _partial_duration = std::chrono::duration<double>::zero();
        _total_duration   = std::chrono::duration<double>::zero();
        _laps             = 0;
    }

    inline void start()
    {
        _last_time_point = std::chrono::high_resolution_clock::now();
    }

    inline void stop()
    {
        auto now          = std::chrono::high_resolution_clock::now();
        _partial_duration = now - _last_time_point;
        _last_time_point  = now;
        _total_duration += _partial_duration;
        ++_laps;
    }

    auto partial_duration() const
    {
        return _partial_duration;
    }

    auto total_duration() const
    {
        return _total_duration;
    }

    template <typename T = std::chrono::seconds>
    auto partial_elapsed() const
    {
        return std::chrono::duration_cast<T>(_partial_duration).count();
    }

    template <typename T = std::chrono::seconds>
    auto total_elapsed() const
    {
        return std::chrono::duration_cast<T>(_total_duration).count();
    }

    StopWatch mean() const
    {
        StopWatch sw(*this);
        sw._total_duration /= _laps;
        return sw;
    }

    friend std::ostream &operator<<(std::ostream &lhs, const StopWatch &rhs)
    {
        // Get the duration.
        auto duration = rhs.total_duration();
        // Get the time.
        auto h = std::chrono::duration_cast<std::chrono::hours>(duration);
        if (h.count()) {
            lhs << std::setw(3) << h.count() << "H ";
            duration -= h;
        }
        auto m = std::chrono::duration_cast<std::chrono::minutes>(duration);
        if (m.count()) {
            lhs << std::setw(3) << m.count() << "M ";
            duration -= m;
        }
        auto s = std::chrono::duration_cast<std::chrono::seconds>(duration);
        if (s.count()) {
            lhs << std::setw(3) << s.count() << "s ";
            duration -= s;
        }
        auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(duration);
        if (ms.count()) {
            lhs << std::setw(3) << ms.count() << "m ";
            duration -= ms;
        }
        auto us = std::chrono::duration_cast<std::chrono::microseconds>(duration);
        if (us.count()) {
            lhs << std::setw(3) << us.count() << "u ";
            duration -= us;
        }
        auto ns = std::chrono::duration_cast<std::chrono::nanoseconds>(duration);
        if (ns.count()) {
            lhs << std::setw(3) << ns.count() << "n ";
            duration -= ns;
        }
        return lhs;
    }
};