/// @file stepper_rk4.hpp
/// @brief Simplification of the code available at:
/// https://github.com/headmyshoulder/odeint-v2

#pragma once

namespace malg::control::solver
{

template <class State, class Time>
class stepper_rk4 {
public:
    using order_type     = unsigned short;
    using time_type      = Time;
    using container_type = State;
    using value_type     = typename State::value_type;

    stepper_rk4(const State &x_init)
        : m_dxdt(x_init),
          m_dxt(x_init),
          m_dxm(x_init),
          m_dxh(x_init),
          m_xt(x_init)
    {
        // Nothing to do.
    }

    constexpr inline order_type order_step() const
    {
        return 4;
    }

    template <class System>
    constexpr inline void do_step(System &system, State &x, Time t, Time dt)
    {
        // dt * dxdt = k1
        system(x, m_dxdt, t);
        this->do_step(system, x, m_dxdt, t, dt);
    }

private:
    template <class System>
    constexpr inline void do_step(System &system, State &x, const State &dxdt, const Time t, const Time dt)
    {
        const Time val1 = static_cast<Time>(1.0);
        const Time dh   = static_cast<Time>(0.5) * dt;
        const Time th   = t + dh;
        const Time dt6  = dt / static_cast<Time>(6.0);
        const Time dt3  = dt / static_cast<Time>(3.0);

        // dt * dxdt = k1 (computed before calling this function)

        // xt = x + dh * dxdt
        it_algebra::scale_sum(m_xt.begin(), m_xt.end(), val1, x.begin(), dh, dxdt.begin());

        // dt * m_dxt = k2
        system(m_xt, m_dxt, th);

        // xt = x + dh*m_dxt
        it_algebra::scale_sum(m_xt.begin(), m_xt.end(), val1, x.begin(), dh, m_dxt.begin());

        // dt * m_dxm = k3
        system(m_xt, m_dxm, th);

        // xt = x + dt*m_dxm
        it_algebra::scale_sum(m_xt.begin(), m_xt.end(), val1, x.begin(), dt, m_dxm.begin());

        // dt * m_dxh = k4
        system(m_xt, m_dxh, t + dt);

        // x += dt/6 * ( m_dxdt + m_dxt + val2*m_dxm )
        it_algebra::scale_sum_inplace(x.begin(), x.end(), dt6, dxdt.begin(), dt3, m_dxt.begin(), dt3, m_dxm.begin(), dt6, m_dxh.begin());
    }

    State m_dxdt;
    State m_dxt;
    State m_dxm;
    State m_dxh;
    State m_xt;
};

} // namespace malg::control::solver