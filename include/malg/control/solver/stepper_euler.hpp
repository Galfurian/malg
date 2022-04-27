/// @file stepper_euler.hpp
/// @brief Simplification of the code available at:
/// https://github.com/headmyshoulder/odeint-v2

#pragma once

namespace malg::control::solver
{

template <class State, class Time>
class stepper_euler {
public:
    using order_type     = unsigned short;
    using time_type      = Time;
    using container_type = State;
    using value_type     = typename State::value_type;

    stepper_euler(const State &x_init)
        : m_dxdt(x_init)
    {
        // Nothing to do.
    }

    constexpr inline order_type order_step() const
    {
        return 1;
    }

    /// @brief Performs one step with the knowledge of dxdt(t)
    /// @param system
    /// @param x
    /// @param dxdt
    /// @param t
    /// @param dt
    template <class System>
    constexpr inline void do_step(System &system, State &x, const State &dxdt, Time t, Time dt) noexcept
    {
        it_algebra::increment(x.begin(), x.end(), dxdt.begin(), dt);
    }

    /// @brief Performs one step.
    /// @param system
    /// @param x
    /// @param t
    /// @param dt
    template <class System>
    constexpr inline void do_step(System &system, State &x, Time t, Time dt) noexcept
    {
        system(x, m_dxdt, t);
        this->do_step(system, x, m_dxdt, t, dt);
    }

private:
    State m_dxdt;
};

} // namespace malg::control::solver