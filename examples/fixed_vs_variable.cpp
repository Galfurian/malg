/// @file test_rlc.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Discretizes a continuous-time state-space description of an RLC
/// circut and simulates it.

#include <stopwatch/stopwatch.hpp>

#include <malg/matrix_array.hpp>
#include <malg/utility.hpp>
#include <malg/matrix.hpp>
#include <malg/linalg.hpp>
#include <malg/math.hpp>
#include <malg/io.hpp>

int main(int, char *[])
{
    const std::size_t samples                      = 500;
    const std::size_t rows                         = 100;
    const std::size_t cols                         = 100;
    const malg::Matrix<double> a                   = malg::utility::rand_matrix<double>(rows, cols, -10, 10);
    const malg::Matrix<double> b                   = malg::utility::rand_matrix<double>(rows, cols, -10, 10);
    const malg::MatrixArray<double, rows, cols> _a = a;
    const malg::MatrixArray<double, rows, cols> _b = b;

    //malg::Matrix<double> c;
    //malg::MatrixArray<double, rows, cols> _c;

    stopwatch::Stopwatch sw;

    std::cout << "[Dynamic Size] Performing multiplication.\n";
    stopwatch::ntimes<samples>(sw, [&]() { a * b; });
    std::cout << "[Dynamic Size] Done, elapsed time : " << sw.mean() << "\n";

    std::cout << "[Fixed Size] Performing multiplication.\n";
    stopwatch::ntimes<samples>(sw, [&]() { _a * _b; });
    std::cout << "[Fixed Size] Done, elapsed time : " << sw.mean() << "\n";

    return 0;
}
