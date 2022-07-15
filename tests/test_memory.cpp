/// @file test_complex.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief

#define MEM_TRACE

#include "malg/matrix.hpp"
#include "malg/io.hpp"

#include <cstring>

void test_memory_1()
{
    {
        std::cout << "\n";
        std::cout << "# " << std::string(40, '=') << "\n";
        std::cout << "# Create the matrix.\n";
        malg::Matrix<int> a{
            { 10, 1, 2, 6, -3 },
            { -4, 4, -1, -4, 7 },
            { -4, -2, 0, 9, -9 },
            { -4, -8, -4, -1, 7 },
            { -4, -9, -6, -6, 8 },
            { -5, 0, -1, 8, 10 }
        };
        std::cout << "# Move the matrix.\n";
        malg::Matrix<int> b = std::move(a);
        std::cout << "\n";
        std::cout << "# Deallocating...\n";
    }
    std::cout << "# END\n";
}

void test_memory_2()
{
    {
        std::cout << "\n";
        std::cout << "# " << std::string(40, '=') << "\n";
        std::cout << "# Create the matrix.\n";
        malg::Matrix<int> a{
            { 1, 2, 3, 4 },
            { 5, 6, 7, 8 },
            { 9, 10, 11, 12 }
        };
        std::cout << "# This should NOT resize the matrix.\n";
        a.resize(3, 4);
        std::cout << "\n";
        std::cout << "# Deallocating...\n";
    }
    std::cout << "# END\n";
}

void test_memory_3()
{
    {
        std::cout << "\n";
        std::cout << "# " << std::string(40, '=') << "\n";
        std::cout << "# Create the matrix.\n";
        malg::Matrix<int> a{
            { 1, 2, 3, 4 },
            { 5, 6, 7, 8 },
            { 9, 10, 11, 12 }
        };
        std::cout << "# This should only require a realloc to bigger size.\n";
#ifdef ROW_MAJOR
        a.resize(4, 4);
#else
        a.resize(3, 5);
#endif

        std::cout << "# This should only require a realloc to smaller size.\n";
#ifdef ROW_MAJOR
        a.resize(3, 4);
#else
        a.resize(3, 3);
#endif
        std::cout << "\n";
        std::cout << "# Deallocating...\n";
    }
    std::cout << "# END\n";
}

void test_memory_4()
{
    {
        std::cout << "\n";
        std::cout << "# " << std::string(40, '=') << "\n";
        std::cout << "# Create the matrix.\n";
        malg::Matrix<int> a(3, 3, { 1, 2, 3, 4, 5, 6, 7, 8, 9 });
        std::cout << a << "\n\n";
        std::cout << "# This should require a malloc for bigger size, and deallocation of the smaller one.\n";
        a.resize(4, 5);
        std::cout << a << "\n";
        std::cout << "\n";
        std::cout << "# Deallocating...\n";
    }
    std::cout << "# END\n";
}

int main(int, char *[])
{
    test_memory_1();
    test_memory_2();
    test_memory_3();
    test_memory_4();
    return 0;
}