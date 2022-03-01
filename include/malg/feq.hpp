/// @file feq.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Functions for checking equality between floating point values.

#pragma once

#include <limits>
#include <cmath>

template <typename T1, typename T2>
inline bool __approximately_equal(T1 a, T2 b, double epsilon = std::numeric_limits<double>::epsilon())
{
    return std::fabs(a - b) <= epsilon * std::fmax(std::fabs(a), std::fabs(b));
}
