/// @file test_complex.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include "malg/complex_math.hpp"
#include "malg/io.hpp"
#include "malg/utility.hpp"

#include <iostream>

template <typename T1, typename T2>
void test_complex_1()
{
    std::complex<T1> t1(8, 6);
    std::complex<T2> t2(4, 3);
    std::cout << t1 << " + " << t2 << " = " << (t1 + t2) << " | " << t2 << " + " << t1 << " = " << (t2 + t1) << "\n";
    std::cout << t1 << " - " << t2 << " = " << (t1 - t2) << " | " << t2 << " - " << t1 << " = " << (t2 - t1) << "\n";
    std::cout << t1 << " * " << t2 << " = " << (t1 * t2) << " | " << t2 << " * " << t1 << " = " << (t2 * t1) << "\n";
    std::cout << t1 << " / " << t2 << " = " << (t1 / t2) << " | " << t2 << " / " << t1 << " = " << (t2 / t1) << "\n";
}

template <typename T1, typename T2>
void test_complex_2()
{
    std::complex<T1> t1(8, 2);
    T2 t2(4);
    std::cout << t1 << " + " << t2 << " = " << (t1 + t2) << " | " << t2 << " + " << t1 << " = " << (t2 + t1) << "\n";
    std::cout << t1 << " - " << t2 << " = " << (t1 - t2) << " | " << t2 << " - " << t1 << " = " << (t2 - t1) << "\n";
    std::cout << t1 << " * " << t2 << " = " << (t1 * t2) << " | " << t2 << " * " << t1 << " = " << (t2 * t1) << "\n";
    std::cout << t1 << " / " << t2 << " = " << (t1 / t2) << " | " << t2 << " / " << t1 << " = " << (t2 / t1) << "\n";
}

template <typename T1, typename T2>
void test_complex_3()
{
    auto t1 = malg::utility::rand_matrix<std::complex<T1>>(4, 4, -10, 10);
    auto t2 = malg::utility::rand_matrix<std::complex<T2>>(4, 4, -10, 10);
    auto t3 = malg::utility::rand_matrix<T2>(4, 4, -10, 10);
    std::cout << malg::to_matlab("t1", t1) << "\n";
    std::cout << malg::to_matlab("t2", t2) << "\n";
    std::cout << malg::to_matlab("t3", t3) << "\n";
    std::cout << t1 << " + " << t2 << " = " << (t1 + t2) << " | " << t2 << " + " << t1 << " = " << (t2 + t1) << "\n";
    std::cout << t1 << " - " << t2 << " = " << (t1 - t2) << " | " << t2 << " - " << t1 << " = " << (t2 - t1) << "\n";
    std::cout << t1 << " * " << t2 << " = " << (t1 * t2) << " | " << t2 << " * " << t1 << " = " << (t2 * t1) << "\n";
}

int main(int, char *[])
{
    test_complex_1<double, double>();
    test_complex_1<double, float>();
    test_complex_1<float, double>();
    test_complex_1<float, float>();

    test_complex_2<double, float>();
    test_complex_2<double, int>();

    test_complex_3<double, double>();

    return 0;
}