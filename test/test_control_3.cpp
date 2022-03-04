#include "malg/control/control.hpp"
#include "malg/io.hpp"

#include <matplot/matplot.h>

int main(int argc, char *argv[])
{
    // Timesteps.
    std::vector<double> vts = matplot::linspace(0.1, 1, 4);
    // Total simulation time.
    double total_time = 20;
    // Define the state-space model.
    malg::control::StateSpace<double> sys{
        { { -3, -1.5 }, { 5, 0 } },
        { { 1 }, { 0 } },
        { { 0.5, 1.5 } },
        { { 0 } }
    };

    matplot::figure();
    matplot::title("Discretization Test");
    matplot::xlabel("Time (s)");
    matplot::legend(true);
    matplot::hold(true);

    // Plot the input that all the discrete systems receive.
    {
        std::vector<double> vt, vu;
        for (double t = 0.; t < total_time; t += vts.front()) {
            vt.emplace_back(t);
            vu.emplace_back(((1 * std::sin(2 * M_PI * t * 0.1)) < 0) ? 1 : 0);
        }
        auto p = matplot::plot(vt, vu);
        p->line_width(2);
        p->display_name("Input");
    }

    // Simulate the discretized systems.
    for (auto ts : vts) {
        // Discretize the systme.
        auto dsys = malg::control::c2d(sys, ts);
        // Input matrix.
        malg::Matrix<double> u{ { .0 } };
        // State matrix.
        malg::Matrix<double> x{ { .0 }, { .0 } };
        // Output matrix.
        malg::Matrix<double> y{ { .0 } };

        // Vector for plot purposes.
        std::vector<double> vt, vy;

        // Simulation.
        for (double t = .0; t < total_time; t += ts) {
            // Generate the input.
            u(0, 0) = ((1 * std::sin(2 * M_PI * t * 0.1)) < 0) ? 1 : 0;

            // Simulate a step of the system.
            y = dsys.C * x + dsys.D * u;
            x = dsys.A * x + dsys.B * u;

            // Save the values for plot purposes.
            vt.emplace_back(t);
            vy.emplace_back(y(0, 0));
        }

        // Plot this specific discrete system.
        auto p = matplot::plot(vt, vy);
        p->line_width(2);
        p->display_name("TS(" + std::to_string(ts) + ")");
    }

    // Show the plot.
    matplot::show();
    return 0;
}
