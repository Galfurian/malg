/// @file complex_math.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Math operators between complex values.

#pragma once

#include "malg/type_traits.hpp"

template <typename T1, typename T2>
inline constexpr auto operator+(const std::complex<T1> &__x, const std::complex<T2> &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r += __y;
    return __r;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2>, T2>>
inline constexpr auto operator+(const std::complex<T1> &__x, const T2 &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r += __y;
    return __r;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline constexpr auto operator+(const T1 &__x, const std::complex<T2> &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r += __y;
    return __r;
}

template <typename T1, typename T2>
inline constexpr auto operator-(const std::complex<T1> &__x, const std::complex<T2> &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r -= __y;
    return __r;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2>, T2>>
inline constexpr auto operator-(const std::complex<T1> &__x, const T2 &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r -= __y;
    return __r;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline constexpr auto operator-(const T1 &__x, const std::complex<T2> &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r -= __y;
    return __r;
}

template <typename T1, typename T2>
inline constexpr auto operator*(const std::complex<T1> &__x, const std::complex<T2> &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r *= __y;
    return __r;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2>, T2>>
inline constexpr auto operator*(const std::complex<T1> &__x, const T2 &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r *= __y;
    return __r;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline constexpr auto operator*(const T1 &__x, const std::complex<T2> &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r *= __y;
    return __r;
}

template <typename T1, typename T2>
inline constexpr auto operator/(const std::complex<T1> &__x, const std::complex<T2> &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r /= __y;
    return __r;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2>, T2>>
inline constexpr auto operator/(const std::complex<T1> &__x, const T2 &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r /= __y;
    return __r;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline constexpr auto operator/(const T1 &__x, const std::complex<T2> &__y)
{
    using T             = malg::extract_common_type_t<T1, T2>;
    std::complex<T> __r = __x;
    __r /= __y;
    return __r;
}
