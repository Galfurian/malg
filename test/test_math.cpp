#include "malg/matrix.hpp"
#include "malg/linalg.hpp"
#include "malg/math.hpp"
#include "malg/io.hpp"

#include <cassert>

int main(int argc, char *argv[])
{
    malg::Matrix<double> a = malg::utility::rand_matrix<double>(100, 100, -5, 5);
    std::cout << malg::to_matlab("a", a) << "\n";
    std::cout << malg::linalg::inverse(a) << "\n";
    return 0;
}
