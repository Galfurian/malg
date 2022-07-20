/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include "malg/control/control.hpp"
#include "malg/io.hpp"

#include <stopwatch/stopwatch.hpp>
#include <solver/detail/observer.hpp>
#include <solver/stepper/stepper_adaptive_euler.hpp>
#include <solver/stepper/stepper_adaptive_rk4.hpp>
#include <solver/stepper/stepper_euler.hpp>
#include <solver/stepper/stepper_rk4.hpp>
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

template <std::size_t DECIMATION>
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
    Model model;
    const State x0{ .0, .0 };
    State x;
    const Time time_start = 0.0;
    const Time time_end   = 1.0;
    const Time time_delta = 0.0001;
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta);
    std::size_t steps     = 0;

    solver::stepper_adaptive_euler<State, Time> adaptive_euler;
    solver::stepper_adaptive_rk4<State, Time> adaptive_rk4;
    solver::stepper_euler<State, Time> euler;
    solver::stepper_rk4<State, Time> rk4;

#ifdef MALG_ENABLE_PLOT
    const auto downsamples = compute_samples<std::size_t>(time_start, time_end, time_delta, time_delta);

    ObserverSave<0> observer_adaptive_euler;
    ObserverSave<0> observer_adaptive_rk4;
    ObserverSave<downsamples> observer_euler;
    ObserverSave<downsamples> observer_rk4;
#else
    NoObserver observer_adaptive_euler;
    NoObserver observer_adaptive_rk4;
    NoObserver observer_euler;
    NoObserver observer_rk4;
#endif

    stopwatch::Stopwatch sw;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Total time points " << samples << "\n";

    x = x0;
    std::cout << "Starting simulation (Adaptive Euler).\n";
    sw.start();
    steps = solver::integrate_adaptive(adaptive_euler, observer_adaptive_euler, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = x0;
    std::cout << "Starting simulation (Adaptive RK4).\n";
    sw.start();
    steps = solver::integrate_adaptive(adaptive_rk4, observer_adaptive_rk4, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = x0;
    std::cout << "Starting simulation (Euler).\n";
    sw.start();
    steps = solver::integrate_fixed(euler, observer_euler, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

    x = x0;
    std::cout << "Starting simulation (RK4).\n";
    sw.start();
    steps = solver::integrate_fixed(rk4, observer_rk4, model, x, time_start, time_end, time_delta);
    sw.stop();
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

#ifdef MALG_ENABLE_PLOT
    auto colors      = matplot::palette::accent(8);
    auto color_index = 0u;
    matplot::line_handle lh;
    matplot::hold(matplot::on);
    lh = matplot::plot(observer_adaptive_euler.time, observer_adaptive_euler.x0);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_adaptive_rk4.time, observer_adaptive_rk4.x0);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_euler.time, observer_euler.x0);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_rk4.time, observer_rk4.x0);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_adaptive_euler.time, observer_adaptive_euler.x1);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_adaptive_rk4.time, observer_adaptive_rk4.x1);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_euler.time, observer_euler.x1);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
    lh = matplot::plot(observer_rk4.time, observer_rk4.x1);
    lh->line_width(3);
    lh->color(matplot::to_array(colors[color_index++]));
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