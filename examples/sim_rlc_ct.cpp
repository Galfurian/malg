/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include "malg/control/solver/solver.hpp"
#include "malg/control/control.hpp"
#include "malg/io.hpp"

#include <stopwatch/stopwatch.hpp>

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
    inline void operator()(const State &x, State &dxdt, Time t)
    {
        u[0] = 1.5 * std::sin(2 * M_PI * t);

        dxdt = malg::dot(sys.A, x) + malg::dot(sys.B, u);
    }
};

/// @brief The dc motor itself.
struct ObserverPrint {
    void operator()(const State &x, Time t)
    {
        std::cout << std::fixed << std::setprecision(4) << t << " " << x[0] << " " << x[1] << "\n";
    }
};

/// @brief The dc motor itself.
struct NoObserver {
    void operator()(const State &, Time)
    {
        // Nothing to do.
    }
};

template <typename T>
constexpr inline T compute_samples(Time time_start, Time time_end, Time time_delta, Time sampling = 1.0)
{
    return static_cast<T>(((time_end - time_start) / time_delta) * sampling);
}

int main(int, char *[])
{
    {
        Model system;
        State state{ .0, .0 };
        Time time_start    = 0.0;
        Time time_end      = 1.0;
        Time time_delta    = 0.0001;
        const auto samples = compute_samples<std::size_t>(time_start, time_end, time_delta);
        unsigned steps     = 0;

        malg::control::solver::stepper_euler<State, Time> stepper(state.size());
        NoObserver observer;
        stopwatch::Stopwatch sw;

        std::cout << std::fixed << std::setprecision(5);
        std::cout << "Fixed step simulation.\n";
        std::cout << "Total time points " << samples << "\n";
        std::cout << "Starting simulation.\n";

        sw.start();
        steps = malg::control::solver::integrate_const(stepper, observer, system, state, time_start, time_end, time_delta);
        sw.stop();

        std::cout << "Terminating simulation.\n";
        std::cout << "Final state " << state << "\n";
        std::cout << "Elapsed time " << sw << "\n";
        std::cout << "Integration steps " << steps << "\n\n";
    }
    {
        Model system;
        State state{ .0, .0 };
        Time time_start    = 0.0;
        Time time_end      = 1.0;
        Time time_delta    = 0.0001;
        const auto samples = compute_samples<std::size_t>(time_start, time_end, time_delta);
        unsigned steps     = 0;

        malg::control::solver::stepper_adaptive_rk4<State, Time, 2> stepper(state.size());
        NoObserver observer;
        stopwatch::Stopwatch sw;

        std::cout << std::fixed << std::setprecision(5);
        std::cout << "Variable step simulation.\n";
        std::cout << "Fixed step total time points " << samples << "\n";
        std::cout << "Starting simulation.\n";

        sw.start();
        steps = malg::control::solver::integrate_adaptive(stepper, observer, system, state, time_start, time_end, time_delta);
        sw.stop();

        std::cout << "Terminating simulation.\n";
        std::cout << "Final state " << state << "\n";
        std::cout << "Elapsed time " << sw << "\n";
        std::cout << "Integration steps " << steps << "\n\n";
    }
    return 0;
}