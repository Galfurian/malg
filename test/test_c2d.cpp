#include "malg/control/control.hpp"
#include "malg/io.hpp"


int main(int argc, char *argv[])
{
    unsigned size = 10;
    // Create the state-space model.
    malg::control::StateSpace<double> sys{
        .A = malg::utility::rand_matrix<double>(size, size, -10, 10),
        .B = malg::utility::rand_matrix<double>(size, 1, -10, 10),
        .C = malg::utility::rand_matrix<double>(size, size, -10, 10),
        .D = malg::utility::zeros<double>(size, 1),
    };
    // Perform the discretization.
    auto dsys = malg::control::c2d(sys, 1e-03);
    std::cout << sys << "\n";
    std::cout << dsys << "\n";
    return 0;
}
