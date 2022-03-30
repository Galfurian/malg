/// @file test_complex.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#include "malg/complex_math.hpp"

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

int main(int argc, char *argv[])
{
    test_complex_1<float, double>();
    test_complex_2<float, double>();
    test_complex_2<double, int>();
    return 0;
}