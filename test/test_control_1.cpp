#include "malg/control/control.hpp"
#include "malg/io.hpp"

int main(int argc, char *argv[])
{
    unsigned size = 4;
    // Create the state-space model.
    malg::control::StateSpace<double> sys{
        malg::utility::rand_matrix<double>(size, size, -10, 10),
        malg::utility::rand_matrix<double>(size, 1, -10, 10),
        malg::utility::rand_matrix<double>(size, size, -10, 10),
        malg::utility::zeros<double>(size, 1),
    };
    // Perform the discretization.
    auto dsys = malg::control::c2d(sys, 1e-03);
    // Select the poles.
    malg::Matrix<double> poles{ { -1, -6, -5, -3 } };
    // Apply ackerman and compute the gain.
    auto K = malg::control::acker(sys.A, sys.B, poles);
    std::cout << sys << "\n";
    std::cout << dsys << "\n";
    std::cout << K << "\n";
    return 0;
}
