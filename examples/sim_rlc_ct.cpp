/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include "malg/control/solver/solver.hpp"
#include "malg/control/control.hpp"
#include "malg/io.hpp"

using Time   = double;
using State  = malg::Vector<double>;
using Input  = malg::Vector<double>;
using Output = malg::Vector<double>;

class Model {
public:
    /// Resistance value.
    double R0;
    /// Inductance value.
    double L0;
    /// Capacitance value.
    double C0;
    /// Continous time state-space model.
    malg::control::StateSpace<double> sys;
    /// Input matrix.
    Input u;
    /// Output matrix.
    Output y;

    /// @brief Construct a new RLC.
    /// @param _R0 resistance value.
    /// @param _L0 inductance value.
    /// @param _C0 capacitance value.
    Model(double _R0 = 100, double _L0 = 0.01, double _C0 = 0.001)
        : R0(_R0),
          L0(_L0),
          C0(_C0),
          sys(),
          u(1),
          y(1)
    {
        sys.A = { { 0., 1. }, { -1. / (L0 * C0), -R0 / L0 } };
        sys.B = { { 0. }, { 1. / (L0 * C0) } };
        sys.C = { { 1., 0. } };
        sys.D = { { 0. } };
    }

    /// @brief DC motor behaviour.
    /// @param x the current state.
    /// @param dxdt the final state.
    /// @param t the current time.
    inline constexpr void operator()(const State &x, State &dxdt, Time t) noexcept
    {
        dxdt = malg::dot(sys.A, x) + malg::dot(sys.B, u);
        y    = malg::dot(sys.C, x) + malg::dot(sys.D, u);
    }
};

/// @brief The dc motor itself.
struct ObserverPrint {
    void operator()(const State &x, const Time &t)
    {
        std::cout << std::fixed << std::setprecision(4) << t << " " << x[0] << " " << x[1] << "\n";
    }
};

int main(int, char *[])
{
    Model rlc;
    State x{ .0, .0 };
    Time time_start = 0.0;
    Time time_end   = 1.0;
    Time dt         = 0.0001;

    malg::control::solver::stepper_euler<State, Time> stepper(x);
    ObserverPrint observer;

    std::cout << std::fixed << std::setprecision(4);
    for (Time t = time_start; t < time_end; t += dt) {
        rlc.u[0] = 1.5 * std::sin(2 * M_PI * t);
        stepper.do_step(rlc, x, t, dt);
        std::cout << t << " " << rlc.u[0] << " -> " << rlc.y[0] << "\n";
    }

    return 0;
}