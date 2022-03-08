#include "malg/control/control.hpp"
#include "malg/io.hpp"

#include <matplot/matplot.h>

int main(int argc, char *argv[])
{
    // Define the model's variables.
    double R0 = 100;
    double L0 = 0.01;
    double C0 = 0.001;
    // Define the time-step.
    double ts = 1e-03;

    // Define the state-space model.
    malg::control::StateSpace<double> sys{
        { { 0., 1. }, { -1. / (L0 * C0), -R0 / L0 } },
        { { 0. }, { 1. / (L0 * C0) } },
        { { 1., 0. } },
        { { 0. } }
    };

    // Discretize the systme.
    auto dsys = malg::control::c2d(sys, ts);

    // == Simulate ============================================================
    malg::Matrix<double> u{ { .0 } };         // Input matrix.
    malg::Matrix<double> x{ { .0 }, { .0 } }; // State matrix.
    malg::Matrix<double> y{ { .0 } };         // Output matrix.

    // Simulation
    for (double t = .0, tt = 4; t < tt; t += ts) {
        u(0, 0) = 1.5 * std::sin(2 * M_PI * t);

        y = dsys.C * x + dsys.D * u;
        x = dsys.A * x + dsys.B * u;
    }

    return 0;
}
