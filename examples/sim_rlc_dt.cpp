/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include "malg/control/control.hpp"
#include "malg/io.hpp"

#include <timelib/stopwatch.hpp>

#ifdef ENABLE_PLOT
#include <gpcpp/gnuplot.hpp>
#endif

using Time     = double;
using Variable = double;
using State    = malg::Vector<double>;
using Input    = malg::Vector<double>;
using Output   = malg::Vector<double>;

template <typename T>
inline T compute_samples(Time time_start, Time time_end, Time time_delta, Time sampling)
{
    return static_cast<T>(((time_end - time_start) / time_delta) * sampling);
}

int main(int, char *[])
{
    // Define the model's variables.
    Variable R0 = 1;
    Variable L0 = 0.5;
    Variable C0 = 0.5;

    // Define the time variables.
    const Time time_start     = 0.0;
    const Time time_end       = 2.0;
    const Time time_delta     = 0.001;
    const std::size_t samples = compute_samples<std::size_t>(time_start, time_end, time_delta, 1.0);
    std::size_t steps         = 0;

    // Define the state-space model.
    malg::control::StateSpace<Variable> sys;
    sys.A = {
        { -R0 / L0, -1. / L0 },
        { 1. / C0, 0 }
    };
    sys.B = {
        { 1. / L0 },
        { 0. }
    };
    sys.C = {
        { 0., 1. }
    };
    sys.D = {
        { 0. }
    };

    // Discretize the systme.
    auto dsys = malg::control::c2d(sys, time_delta);

#ifdef ENABLE_PLOT
    std::vector<Variable> time, u0, x0, x1;
#endif

    // == Simulate ============================================================
    Input u{ .0 };     // Input matrix.
    State x{ .0, .0 }; // State matrix.
    Output y{ .0 };    // Output matrix.

    timelib::Stopwatch sw;

    std::cout << std::fixed << std::setprecision(5);
    std::cout << "Total time points " << samples << "\n";

    std::cout << "Starting simulation.\n";
    sw.start();

    for (Time t = time_start; t < time_end; t += time_delta, ++steps) {
        u[0] = 3.3 * std::sin(2 * M_PI * t);

        x = malg::dot(dsys.A, x) + malg::dot(dsys.B, u);
        y = malg::dot(dsys.C, x) + malg::dot(dsys.D, u);

#ifdef ENABLE_PLOT
        time.emplace_back(t);
        u0.emplace_back(u[0]);
        x0.emplace_back(x[0]);
        x1.emplace_back(x[1]);
#else
        std::cout << t << " " << u[0] << " -> " << y[0] << "\n";
#endif
    }
    std::cout << "Terminating simulation.\n";
    std::cout << "Elapsed time " << sw << "\n";
    std::cout << "Integration steps " << steps << "\n\n";

#ifdef ENABLE_PLOT
    // Create a Gnuplot instance.
    gpcpp::Gnuplot gnuplot;

    // Set up the plot with grid, labels, and line widths
    gnuplot.set_title("Plot of u0, x0, x1")
        .set_xlabel("Time (s)") // X-axis label
        .set_ylabel("Values")   // Y-axis label
        .set_grid()
        .set_legend();

    // Plot u0 with line width 3
    gnuplot.set_line_width(3)                     // Line width
        .set_plot_type(gpcpp::plot_type_t::lines) // Line style
        .set_line_type(gpcpp::line_type_t::solid) // Solid line style
        .plot_xy(time, u0, "u0");

    // Plot x0 with line width 3
    gnuplot.set_line_width(3)                     // Line width
        .set_plot_type(gpcpp::plot_type_t::lines) // Line style
        .set_line_type(gpcpp::line_type_t::solid) // Solid line style
        .plot_xy(time, x0, "x0");

    // Plot x1 with line width 3
    gnuplot.set_line_width(3)                     // Line width
        .set_plot_type(gpcpp::plot_type_t::lines) // Line style
        .set_line_type(gpcpp::line_type_t::solid) // Solid line style
        .plot_xy(time, x1, "x1");

    gnuplot.show();
#endif
    return 0;
}