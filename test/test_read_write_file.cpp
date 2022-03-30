/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include "malg/control/control.hpp"
#include "malg/io.hpp"

int test_vec()
{
    malg::Vector<double> data_in, data_out;
    // Generate the vector.
    data_out = malg::utility::rand_vector<double>(5, -10, 10);
    // Save the vector.
    std::ofstream data_out_file("vector.malg");
    if (!data_out_file.is_open())
        return 1;
    data_out_file << std::setprecision(9);
    data_out_file << data_out;
    data_out_file.close();
    // Read the vector.
    std::ifstream data_in_file("vector.malg");
    if (!data_in_file.is_open())
        return 1;
    data_in_file >> data_in;
    data_in_file.close();
    // Set the tollerance.
    feq::tolerance() = 1e-06;
    // Check correctness.
    if (!malg::all(data_in == data_out))
        return 1;
    return 0;
}

int test_matrix()
{
    malg::Matrix<double> data_in, data_out;
    // Generate the matrix.
    data_out = malg::utility::rand_matrix<double>(5, 5, -10, 10);
    // Save the matrix.
    std::ofstream data_out_file("matrix.malg");
    if (!data_out_file.is_open())
        return 1;
    data_out_file << std::setprecision(9);
    data_out_file << data_out;
    data_out_file.close();
    // Read the matrix.
    std::ifstream data_in_file("matrix.malg");
    if (!data_in_file.is_open())
        return 1;
    data_in_file >> data_in;
    data_in_file.close();
    // Set the tollerance.
    feq::tolerance() = 1e-06;
    // Check correctness.
    if (!malg::all(data_in == data_out))
        return 1;
    return 0;
}

int main(int argc, char *argv[])
{
    if (test_vec())
        return 1;
    if (test_matrix())
        return 1;
    return 0;
}