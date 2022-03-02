#include "malg/control/control.hpp"
#include "malg/io.hpp"

int main(int argc, char *argv[])
{
    unsigned size = 3;
    // Create the state-space model.
    malg::control::StateSpace<double> sys{
        .A = malg::utility::rand_matrix<double>(size, size, -10, 10),
        .B = malg::utility::rand_matrix<double>(size, 1, -10, 10),
        .C = malg::utility::rand_matrix<double>(size, size, -10, 10),
        .D = malg::utility::zeros<double>(size, 1),
    };
    // Perform the discretization.
    auto dsys = malg::control::c2d(sys, 1e-03);
    // Select the poles.
    malg::Matrix<double> poles{ { -1, -6, -5 } };
    // Apply ackerman and compute the gain.
    auto K = malg::control::acker(dsys.A, dsys.B, poles);

    std::cout << sys << "\n";
    std::cout << dsys << "\n";
    std::cout << K << "\n";
    return 0;
}
