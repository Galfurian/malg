/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include "malg/control/control.hpp"
#include "malg/io.hpp"

#include <timelib/stopwatch.hpp>

#include <chainsaw/detail/observer.hpp>
#include <chainsaw/solver.hpp>
#include <chainsaw/stepper/stepper_adaptive.hpp>
#include <chainsaw/stepper/stepper_euler.hpp>
#include <chainsaw/stepper/stepper_rk4.hpp>

#ifdef ENABLE_PLOT
#include <gpcpp/gnuplot.hpp>
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
    explicit Model()
        : Model(1, 0.5, 0.5)
    {
        // Nothing to do.
    }

    /// @brief Construct a new RLC.
    /// @param _R0 resistance value.
    /// @param _L0 inductance value.
    /// @param _C0 capacitance value.
    Model(Variable _R0, Variable _L0, Variable _C0)
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

template <std::size_t DECIMATION = 0>
struct ObserverSave : public chainsaw::detail::ObserverDecimate<State, Time, DECIMATION> {
    inline void operator()(const State &x, const Time &t) override
    {
        if (this->observe()) {
            time.emplace_back(t);
            x0.emplace_back(x[0]);
            x1.emplace_back(x[1]);
        }
    }
    std::vector<Variable> time, x0, x1;
};

template <typename T>
inline T compute_samples(Time time_start, Time time_end, Time time_delta, Time sampling)
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
    const auto samples    = compute_samples<std::size_t>(time_start, time_end, time_delta, 1.0);

    // Setup the solvers.
    const auto Error      = chainsaw::ErrorFormula::Mixed;
    const auto Iterations = 3;
    using Euler           = chainsaw::stepper_euler<State, Time>;
    using Rk4             = chainsaw::stepper_rk4<State, Time>;
    using AdaptiveEuler   = chainsaw::stepper_adaptive<Euler, Iterations, Error>;
    using AdaptiveRk4     = chainsaw::stepper_adaptive<Rk4, Iterations, Error>;

    AdaptiveEuler adaptive_euler;
    adaptive_euler.set_tollerance(1e-02);
    adaptive_euler.set_min_delta(1e-09);
    adaptive_euler.set_max_delta(1e-02);
    AdaptiveRk4 adaptive_rk4;
    adaptive_rk4.set_tollerance(1e-02);
    adaptive_rk4.set_min_delta(1e-09);
    adaptive_rk4.set_max_delta(1e-02);
    Euler euler;
    Rk4 rk4;

#ifdef ENABLE_PLOT
    using Observer = ObserverSave<0>;
#else
    using Observer = chainsaw::detail::ObserverPrint<State, Time, 0>;
#endif
    Observer obs_adaptive_euler;
    Observer obs_adaptive_rk4;
    Observer obs_euler;
    Observer obs_rk4;

    timelib::Stopwatch sw;

    std::cout << std::fixed;
    std::cout << "Total time points with fixed integration step " << samples << "\n\n";

    std::cout << "Simulating with `Adaptive Euler`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_adaptive(adaptive_euler, obs_adaptive_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Adaptive RK4`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_adaptive(adaptive_rk4, obs_adaptive_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `Euler`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_fixed(euler, obs_euler, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "Simulating with `RK4`...\n";
    x = x0;
    sw.start();
    chainsaw::integrate_fixed(rk4, obs_rk4, model, x, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Adaptive Euler took " << std::setw(12) << adaptive_euler.steps() << " steps, for a total of " << sw[0] << "\n";
    std::cout << "    Adaptive RK4   took " << std::setw(12) << adaptive_rk4.steps() << " steps, for a total of " << sw[1] << "\n";
    std::cout << "    Euler          took " << std::setw(12) << euler.steps() << " steps, for a total of " << sw[2] << "\n";
    std::cout << "    RK4            took " << std::setw(12) << rk4.steps() << " steps, for a total of " << sw[3] << "\n";

#ifdef ENABLE_PLOT

    // Create a Gnuplot instance.
    gpcpp::Gnuplot gnuplot;

    // Set up the plot with grid, labels, and line widths
    gnuplot.set_title("Comparison of Adaptive Methods and RK4")
        .set_xlabel("Time (s)") // X-axis label
        .set_ylabel("Values")   // Y-axis label
        .set_grid()             // Show grid.
        .set_legend();          // Enable legend.

    // Plot Adaptive Euler x0 with marker size 16
    gnuplot.set_plot_type(gpcpp::plot_type_t::points) // Points style for scatter plot
        .set_line_type(gpcpp::line_type_t::solid)     // Solid line style (if needed)
        .set_line_width(3)                            // Line width
        .set_point_size(1)                            // Marker size
        .plot_xy(obs_adaptive_euler.time, obs_adaptive_euler.x0, "Adaptive Euler.x0");

    // Plot Adaptive RK4 x0 with marker size 32
    gnuplot.set_plot_type(gpcpp::plot_type_t::points) // Points style for scatter plot
        .set_line_type(gpcpp::line_type_t::solid)     // Solid line style (if needed)
        .set_line_width(3)                            // Line width
        .set_point_size(2)                            // Marker size
        .plot_xy(obs_adaptive_rk4.time, obs_adaptive_rk4.x0, "Adaptive RK4.x0");

    // Plot Euler x0 with line style
    gnuplot.set_plot_type(gpcpp::plot_type_t::lines) // Lines style
        .set_line_type(gpcpp::line_type_t::solid)    // Solid line style
        .set_line_width(3)                           // Line width
        .plot_xy(obs_euler.time, obs_euler.x0, "Euler.x0");

    // Plot RK4 x0 with line style
    gnuplot.set_plot_type(gpcpp::plot_type_t::lines) // Lines style
        .set_line_type(gpcpp::line_type_t::solid)    // Solid line style
        .set_line_width(3)                           // Line width
        .plot_xy(obs_rk4.time, obs_rk4.x0, "RK4.x0");

    // Plot Adaptive Euler x1 with marker size 16
    gnuplot.set_plot_type(gpcpp::plot_type_t::points) // Points style for scatter plot
        .set_line_type(gpcpp::line_type_t::solid)     // Solid line style (if needed)
        .set_line_width(3)                            // Line width
        .set_point_size(1)                            // Marker size
        .plot_xy(obs_adaptive_euler.time, obs_adaptive_euler.x1, "Adaptive Euler.x1");

    // Plot Adaptive RK4 x1 with marker size 32
    gnuplot.set_plot_type(gpcpp::plot_type_t::points) // Points style for scatter plot
        .set_line_type(gpcpp::line_type_t::solid)     // Solid line style (if needed)
        .set_line_width(3)                            // Line width
        .set_point_size(2)                            // Marker size
        .plot_xy(obs_adaptive_rk4.time, obs_adaptive_rk4.x1, "Adaptive RK4.x1");

    // Plot Euler x1 with line style
    gnuplot.set_plot_type(gpcpp::plot_type_t::lines) // Lines style
        .set_line_type(gpcpp::line_type_t::solid)    // Solid line style
        .set_line_width(3)                           // Line width
        .plot_xy(obs_euler.time, obs_euler.x1, "Euler.x1");

    // Plot RK4 x1 with line style
    gnuplot.set_plot_type(gpcpp::plot_type_t::lines) // Lines style
        .set_line_type(gpcpp::line_type_t::solid)    // Solid line style
        .set_line_width(3)                           // Line width
        .plot_xy(obs_rk4.time, obs_rk4.x1, "RK4.x1");

    // Show the plot
    gnuplot.show();

#endif
    return 0;
}