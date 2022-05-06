/// @file complex_math.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Math operators between complex values.

#pragma once

#include "malg/type_traits.hpp"

/// @brief Sums two complex values.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2>
inline constexpr auto operator+(const std::complex<T1> &lhs, const std::complex<T2> &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp += rhs;
    return tmp;
}

/// @brief Sums a complex value and a numeric value.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2>, T2>>
inline constexpr auto operator+(const std::complex<T1> &lhs, const T2 &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp += rhs;
    return tmp;
}

/// @brief Sums a numeric value and a complex value.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline constexpr auto operator+(const T1 &lhs, const std::complex<T2> &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp += rhs;
    return tmp;
}

/// @brief Subtracts a complex value from another one.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2>
inline constexpr auto operator-(const std::complex<T1> &lhs, const std::complex<T2> &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp -= rhs;
    return tmp;
}

/// @brief Subtracts a numeric value from a complex value.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2>, T2>>
inline constexpr auto operator-(const std::complex<T1> &lhs, const T2 &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp -= rhs;
    return tmp;
}

/// @brief Subtracts a complex value from a numeric value.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline constexpr auto operator-(const T1 &lhs, const std::complex<T2> &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp -= rhs;
    return tmp;
}

/// @brief Multiplies two complex values.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2>
inline constexpr auto operator*(const std::complex<T1> &lhs, const std::complex<T2> &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp *= rhs;
    return tmp;
}

/// @brief Multiplies a complex value and a numeric value.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2>, T2>>
inline constexpr auto operator*(const std::complex<T1> &lhs, const T2 &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp *= rhs;
    return tmp;
}

/// @brief Multiplies a numeric value and a complex value.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline constexpr auto operator*(const T1 &lhs, const std::complex<T2> &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp *= rhs;
    return tmp;
}

/// @brief Divides two complex values.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2>
inline constexpr auto operator/(const std::complex<T1> &lhs, const std::complex<T2> &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp /= rhs;
    return tmp;
}

/// @brief Divides a complex value by a numeric value.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2>, T2>>
inline constexpr auto operator/(const std::complex<T1> &lhs, const T2 &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp /= rhs;
    return tmp;
}

/// @brief Divides a numeric value by a complex value.
/// @param lhs the fist value.
/// @param rhs the fist value.
/// @returns the result of the operation.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline constexpr auto operator/(const T1 &lhs, const std::complex<T2> &rhs)
{
    std::complex<std::common_type_t<T1, T2>> tmp = lhs;
    tmp /= rhs;
    return tmp;
}
