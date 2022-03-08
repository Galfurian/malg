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

        // Simulation.
        for (double t = .0; t < total_time; t += ts) {
            // Generate the input.
            u(0, 0) = ((1 * std::sin(2 * M_PI * t * 0.1)) < 0) ? 1 : 0;

            // Simulate a step of the system.
            y = dsys.C * x + dsys.D * u;
            x = dsys.A * x + dsys.B * u;

        }

    }

    return 0;
}
