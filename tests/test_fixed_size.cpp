/// @file test_fixed_size.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include "malg/fixed_size/math.hpp"
#include "malg/fixed_size/view.hpp"
#include "malg/fixed_size/io.hpp"

int main(int, char *[])
{
    malg::Vector<double, 5> v1 = { 1, 2, 3, 4, 5 };
    malg::Vector<double, 5> v2 = { 5, 4, 3, 2, 1 };

    malg::Matrix<double, 3, 3> m1 = { { 1, 2, 3 },
                                      { 4, 5, 6 },
                                      { 7, 8, 9 } };

    malg::Matrix<double, 3, 3> m2 = { { 9, 8, 7 },
                                      { 6, 5, 4 },
                                      { 3, 2, 1 } };

    std::cout << v1 << "\n";
    std::cout << v2 << "\n";
    std::cout << m1 << "\n";
    std::cout << m2 << "\n";

    std::cout << (v1 + v2) << "\n";
    std::cout << (v1 + 1.5) << "\n";
    std::cout << (1.5 + v1) << "\n";
    std::cout << (v1 == v2) << "\n";

    std::cout << (m1 + m2) << "\n";
    std::cout << (m1 + 1.5) << "\n";
    std::cout << (1.5 + m1) << "\n";
    std::cout << (m1 == m2) << "\n";

    for (std::size_t i = 0; i < (3 + 3); ++i) {
        std::cout << m1[i] << "\n";
    }

    // auto w1 = malg::view<0, 2, 0, 2>(m1);
}