/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include "malg/control/control.hpp"
#include "malg/io.hpp"

using Time   = double;
using State  = malg::Vector<double>;
using Input  = malg::Vector<double>;
using Output = malg::Vector<double>;

int main(int, char *[])
{
    // Define the model's variables.
    double R0 = 100;
    double L0 = 0.01;
    double C0 = 0.001;

    // Define the time variables.
    Time time_start = 0.0;
    Time time_end   = 1.0;
    Time dt         = 0.0001;

    // Define the state-space model.
    malg::control::StateSpace<double> sys{
        { { 0., 1. }, { -1. / (L0 * C0), -R0 / L0 } },
        { { 0. }, { 1. / (L0 * C0) } },
        { { 1., 0. } },
        { { 0. } }
    };

    // Discretize the systme.
    auto dsys = malg::control::c2d(sys, dt);

    // == Simulate ============================================================
    Input u{ .0 };     // Input matrix.
    State x{ .0, .0 }; // State matrix.
    Output y{ .0 };    // Output matrix.

    // Simulation
    std::cout << std::fixed << std::setprecision(4);
    for (Time t = time_start; t < time_end; t += dt) {
        u[0] = 1.5 * std::sin(2 * M_PI * t);

        x = malg::dot(dsys.A, x) + malg::dot(dsys.B, u);
        y = malg::dot(dsys.C, x) + malg::dot(dsys.D, u);

        std::cout << t << " " << u[0] << " -> " << y[0] << "\n";
    }

    return 0;
}