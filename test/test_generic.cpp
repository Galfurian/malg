#include "malg/control/control.hpp"
#include "malg/io.hpp"

#include <cassert>

int test_sum()
{
    auto a = malg::utility::rand_matrix<int>(5, 5, -100, 100);
    auto b = malg::utility::rand_matrix<int>(5, 5, -100, 100);
    std::cout << malg::to_cpp("a", a) << ";\n";
    std::cout << malg::to_cpp("b", b) << ";\n";
    std::cout << malg::to_cpp("sum", a + b) << ";\n";
    std::cout << malg::to_cpp("sum", a + b) << ";\n";
    return 0;
}

int main(int argc, char *argv[])
{
    test_sum();
    std::cout << malg::control::polyreduce(malg::Vector<int>({ 0, 0, 1, 2, 3 })) << " == 1 2 3 \n";
    std::cout << malg::control::polyreduce(malg::Vector<int>({ 1, 2, 3, 0, 0 })) << " == 1 2 3 0 0 \n";
    std::cout << malg::control::polyreduce(malg::Vector<int>({ 1, 0, 3 })) << " == 1 0 3 \n";
    std::cout << malg::control::polyreduce(malg::Vector<int>({ 0, 0, 0 })) << " == 0 \n";
    return 0;
}
