/// @file solver.hpp
/// @brief Simplification of the code available at:
/// https://github.com/headmyshoulder/odeint-v2

#pragma once

#include "it_algebra.hpp"
#include "stepper_euler.hpp"
#include "stepper_rk4.hpp"

namespace malg::control::solver
{

template <class Stepper, class System, class Observer>
constexpr inline void integrate_one_step(
    Stepper &stepper,
    System &system,
    typename Stepper::container_type &state,
    typename Stepper::time_type &time,
    typename Stepper::time_type time_delta,
    Observer observer) noexcept
{
    stepper.do_step(system, state, time, time_delta);
    observer(state, time);
    time += time_delta;
}

template <class Stepper, class System, class Observer>
constexpr inline unsigned integrate_const(
    Stepper &stepper,
    System &system,
    typename Stepper::container_type &state,
    typename Stepper::time_type start_time,
    typename Stepper::time_type end_time,
    typename Stepper::time_type time_delta,
    Observer observer) noexcept
{
    unsigned iteration = 0;
    observer(state, start_time);
    while (start_time < end_time) {
        integrate_one_step(stepper, system, state, start_time, time_delta, observer);
        ++iteration;
    }
    return iteration;
}

} // namespace malg::control::solver