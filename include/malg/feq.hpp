/// @file feq.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Functions for checking equality between floating point values.

#pragma once

#include <limits>
#include <cmath>

namespace feq
{

static inline auto &tolerance()
{
    static double tol = std::numeric_limits<double>::epsilon();
    return tol;
}

template <typename T1, typename T2>
inline bool approximately_equal(T1 a, T2 b)
{
    return std::fabs(a - b) <= tolerance() * std::fmax(std::fabs(a), std::fabs(b));
}

} // namespace feq
