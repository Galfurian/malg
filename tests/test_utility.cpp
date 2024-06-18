#include "malg/control/control.hpp"
#include "malg/io.hpp"

int test_linspace()
{
    {
        auto result    = malg::utility::linspace<int>(0, 40, 5);
        auto reference = malg::Vector<int>({ 0, 10, 20, 30, 40 });
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    {
        auto result    = malg::utility::linspace<double>(0, 1, 5);
        auto reference = malg::Vector<double>({ 0, 0.25, 0.5, 0.75, 1 });
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    return 0;
}

int test_stack()
{
    malg::Matrix<int> a{
        { 10, 1 },
        { -4, 4 },
    };
    malg::Matrix<int> b{
        { 9, -9 },
        { -1, 7 },
    };
    {
        malg::Matrix<int> reference = malg::Matrix<int>{
            { 10, 1 },
            { -4, 4 },
            { 9, -9 },
            { -1, 7 }
        };
        malg::Matrix<int> result = malg::utility::vstack(a, b);
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    {
        malg::Matrix<int> reference = malg::Matrix<int>{
            { 10, 1 },
            { -4, 4 },
            { 0, 0 },
            { 0, 0 }
        };
        malg::Matrix<int> result = malg::utility::vstack(a, malg::utility::zeros<int>(2, 2));
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    {
        malg::Matrix<int> reference = malg::Matrix<int>{
            { 10, 1, 9, -9 },
            { -4, 4, -1, 7 },
        };
        malg::Matrix<int> result = malg::utility::hstack(a, b);
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    {
        malg::Matrix<int> reference = malg::Matrix<int>{
            { 10, 1, 0, 0 },
            { -4, 4, 0, 0 },
        };
        malg::Matrix<int> result = malg::utility::hstack(a, malg::utility::zeros<int>(2, 2));
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    return 0;
}

int test_transform()
{
    {
        const auto input     = malg::Matrix<int>{ { 10 }, { -4 }, { 9 }, { -1 } };
        const auto reference = malg::Vector<int>{ 10, -4, 9, -1 };
        const auto result    = malg::utility::to_vector(input);
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    {
        const auto input     = malg::Matrix<int>{ { 10, -4, 9, -1 } };
        const auto reference = malg::Vector<int>{ 10, -4, 9, -1 };
        const auto result    = malg::utility::to_vector(input);
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    {
        const auto input     = malg::Vector<int>{ 10, -4, 9, -1 };
        const auto reference = malg::Matrix<int>{ { 10 }, { -4 }, { 9 }, { -1 } };
        const auto result    = malg::utility::to_matrix(input, false);
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    {
        const auto input     = malg::Vector<int>{ 10, -4, 9, -1 };
        const auto reference = malg::Matrix<int>{ { 10, -4, 9, -1 } };
        const auto result    = malg::utility::to_matrix(input, true);
        if (!malg::all(result == reference)) {
            std::cerr << result << " != " << reference << "\n";
            return 1;
        }
    }
    return 0;
}

int main(int, char *[])
{
    if (test_linspace()) {
        return 1;
    }
    if (test_stack()) {
        return 1;
    }
    if (test_transform()) {
        return 1;
    }
    return 0;
}
