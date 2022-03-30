#include "malg/control/control.hpp"
#include "malg/io.hpp"

int main(int argc, char *argv[])
{
    // Define the model's variables.
    double R0 = 100;
    double L0 = 0.01;
    double C0 = 0.001;
    // Define the time-step.
    double ts  = 1e-03;
    unsigned T = 1000;

    // Define the state-space model.
    malg::control::StateSpace<double> sys{
        { { 0., 1. }, { -1. / (L0 * C0), -R0 / L0 } },
        { { 0. }, { 1. / (L0 * C0) } },
        { { 1., 0. } },
        { { 0. } }
    };

    std::cout << malg::to_matlab("A", sys.A) << "\n";
    std::cout << malg::to_matlab("B", sys.B) << "\n";
    std::cout << malg::to_matlab("C", sys.C) << "\n";
    std::cout << malg::to_matlab("D", sys.D) << "\n";

    // Discretize the systme.
    auto dsys = malg::control::c2d(sys, ts);

    // == Simulate ============================================================
    // Input matrix.
    malg::Matrix<double> u(1, T);
    // State matrix.
    malg::Matrix<double> x(2, T);
    // Output matrix.
    malg::Matrix<double> y(1, T);

    malg::utility::col(u, 0) = {
        { 1 }
    };

    std::cout << std::fixed;
    // Simulation
    for (unsigned i = 0; i < (T - 1); ++i) {
        malg::utility::col(y, i)     = dsys.C * malg::utility::col(x, i) + dsys.D * malg::utility::col(u, i);
        malg::utility::col(x, i + 1) = dsys.A * malg::utility::col(x, i) + dsys.B * malg::utility::col(u, i);
    }

    std::cout << y << "\n";

    return 0;
}
