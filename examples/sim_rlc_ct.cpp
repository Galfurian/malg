/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include "malg/control/control.hpp"
#include "malg/io.hpp"

#include <stopwatch/stopwatch.hpp>

#include <solver/stepper/stepper_adaptive.hpp>
#include <solver/stepper/stepper_euler.hpp>
#include <solver/stepper/stepper_rk4.hpp>
#include <solver/detail/observer.hpp>
#include <solver/solver.hpp>

#ifdef MALG_ENABLE_PLOT
#include <matplot/matplot.h>
#endif

using Time     = double;
using Variable = double;
using State    = malg::Vector<Variable>;
using Input    = malg::Vector<Variable>;
using Output   = malg::Vector<Variable>;

class Model {
public:
    /// Resistance value.
    Variable R0;
    /// Inductance value.
    Variable L0;
    /// Capacitance value.
    Variable C0;
    /// Continous time state-space model.
    malg::control::StateSpace<Variable> sys;
    /// Input matrix.
    Input u;
    /// Output matrix.
    Output y;

    /// @brief Construct a new RLC.
    /// @param _R0 resistance value.
    /// @param _L0 inductance value.
    /// @param _C0 capacitance value.
    Model(Variable _R0 = 1, Variable _L0 = 0.5, Variable _C0 = 0.5)
        : R0(_R0),
          L0(_L0),
          C0(_C0),
          sys(),
          u(1),
          y(1)
    {
        sys.A = {
            { -R0 / L0, 1. / L0 },
            { 1. / C0, 0 }
        };
        sys.B = {
            { 1. / L0 },
            { 0. }
        };
        sys.C = {
            { 1., 0. }
        };
        sys.D = {
            { 0. }
        };
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

template <std::size_t DECIMATION = 0>
struct ObserverSave : public solver::detail::DecimationObserver<DECIMATION> {
    std::vector<Variable> time, x0, x1;
    ObserverSave() = default;
    inline void operator()(const State &x, const Time &t) noexcept
    {
        if (this->observe()) {
            time.emplace_back(t);
            x0.emplace_back(x[0]);
            x1.emplace_back(x[1]);
        }
    }
};

template <typename T>
constexpr inline T compute_samples(Time time_start, Time time_end, Time time_delta, Time sampling = 1.0)
{
    return static_cast<T>(((time_end - time_start) / time_delta) * sampling);
}

// matplot::line_spec::marker_style next_marker(unsigned & marker)
// {
//     return static_cast<matplot::line_spec::marker_style>(marker++);
// }

int main(int, char *[])
{
    Model model;

    // Initial and runtime states.
    const State x0{ .0, .0 };
    State x;

    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 1.0;
    const Time time_delta = 0.0001;
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta);

    // Setup the solvers.
    const auto Error      = solver::ErrorFormula::Mixed;
    const auto Iterations = 3;
    using Euler           = solver::stepper_euler<State, Time>;
    using Rk4             = solver::stepper_rk4<State, Time>;
    using AdaptiveEuler   = solver::stepper_adaptive<Euler, Iterations, Error>;
    using AdaptiveRk4     = solver::stepper_adaptive<Rk4, Iterations, Error>;

    AdaptiveEuler adaptive_euler(time_delta);
    AdaptiveRk4 adaptive_rk4(time_delta);
    Euler euler;
    Rk4 rk4;

#ifdef MALG_ENABLE_PLOT
    ObserverSave obs_adaptive_euler;
    ObserverSave obs_adaptive_rk4;
    ObserverSave obs_euler;
    ObserverSave obs_rk4;
#else
    solver::detail::NoObserver obs_adaptive_euler;
    solver::detail::NoObserver obs_adaptive_rk4;
    solver::detail::NoObserver obs_euler;
    solver::detail::NoObserver obs_rk4;
#endif

    stopwatch::Stopwatch sw;

    std::cout << std::fixed;
    std::cout << "Total time points with fixed integration step " << samples << "\n\n";

    std::cout << "Simulating with `Adaptive Euler`...\n";
    x = x0;
    sw.start();
    solver::integrate_adaptive(adaptive_euler, obs_adaptive_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Adaptive RK4`...\n";
    x = x0;
    sw.start();
    solver::integrate_adaptive(adaptive_rk4, obs_adaptive_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Euler`...\n";
    x = x0;
    sw.start();
    solver::integrate_fixed(euler, obs_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `RK4`...\n";
    x = x0;
    sw.start();
    solver::integrate_fixed(rk4, obs_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Adaptive Euler took " << std::setw(12) << adaptive_euler.steps() << " steps, for a total of " << sw[0] << "\n";
    std::cout << "    Adaptive RK4   took " << std::setw(12) << adaptive_rk4.steps() << " steps, for a total of " << sw[1] << "\n";
    std::cout << "    Euler          took " << std::setw(12) << euler.steps() << " steps, for a total of " << sw[2] << "\n";
    std::cout << "    RK4            took " << std::setw(12) << rk4.steps() << " steps, for a total of " << sw[3] << "\n";

#ifdef MALG_ENABLE_PLOT
    auto colors = matplot::palette::accent(4);
    auto color  = colors.begin();

    matplot::hold(matplot::on);

    color = colors.begin();
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.x0, 16)->color(matplot::to_array(*color++)).marker_style("o");
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.x0, 32)->color(matplot::to_array(*color++)).marker_style("d");
    matplot::plot(obs_euler.time, obs_euler.x0)->color(matplot::to_array(*color++));
    matplot::plot(obs_rk4.time, obs_rk4.x0)->color(matplot::to_array(*color++));

    color = colors.begin();
    matplot::scatter(obs_adaptive_euler.time, obs_adaptive_euler.x1, 16)->color(matplot::to_array(*color++)).marker_style("o");
    matplot::scatter(obs_adaptive_rk4.time, obs_adaptive_rk4.x1, 32)->color(matplot::to_array(*color++)).marker_style("d");
    matplot::plot(obs_euler.time, obs_euler.x1)->color(matplot::to_array(*color++));
    matplot::plot(obs_rk4.time, obs_rk4.x1)->color(matplot::to_array(*color++));

    matplot::legend(
        { "Adaptive Euler.x0",
          "Adaptive RK4.x0",
          "Euler.x0",
          "RK4.x0",
          "Adaptive Euler.x1",
          "Adaptive RK4.x1",
          "Euler.x1",
          "RK4.x1" });
    matplot::show();
#endif
    return 0;
}