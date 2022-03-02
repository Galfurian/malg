/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include "malg/control/control.hpp"
#include "malg/io.hpp"

int main(int argc, char *argv[])
{
    auto vec_out = malg::utility::rand_vector<double>(5, -10, 10);
    std::cout << malg::dump_vector("vec_out", vec_out) << "\n";
    std::ofstream vec_out_file("vector.malg");
    if (vec_out_file.is_open()) {
        vec_out_file << vec_out;
        vec_out_file.close();
        std::ifstream vec_in_file("vector.malg");
        if (vec_in_file.is_open()) {
            malg::Vector<double> vec_in;
            vec_in_file >> vec_in;
            std::cout << malg::dump_vector("vec_in", vec_in) << "\n";
        }
    }

    auto mat_out = malg::utility::rand_matrix<double>(5, 5, -10, 10);
    std::cout << malg::dump_matrix("mat_out", mat_out) << "\n";
    std::ofstream mat_out_file("matrix.malg");
    if (mat_out_file.is_open()) {
        mat_out_file << mat_out;
        mat_out_file.close();
        std::ifstream mat_in_file("matrix.malg");
        if (mat_in_file.is_open()) {
            malg::Matrix<double> mat_in;
            mat_in_file >> mat_in;
            std::cout << malg::dump_matrix("mat_in", mat_in) << "\n";
        }
    }

    return 0;
}