/// @file feq.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Functions for checking equality between floating point values.

#pragma once

#include <cmath>
#include <limits>

namespace malg::feq
{

static inline auto &tolerance()
{
    static double tol = std::numeric_limits<double>::epsilon();
    return tol;
}

/// @brief Checks if the two floating point values are equal.
/// @param a the first value.
/// @param b the second value.
/// @returns true if they are approximately equal.
/// @returns false otherwise.
template <typename T1, typename T2>
inline bool approximately_equal(T1 a, T2 b)
{
    return std::fabs(a - b) <= tolerance() * std::fmax(std::fabs(a), std::fabs(b));
}

} // namespace feq
