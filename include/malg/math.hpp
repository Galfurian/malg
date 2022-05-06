/// @file math.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Basic mathematical functions.

#pragma once

#include "malg/complex_math.hpp"
#include "malg/matrix_base.hpp"
#include "malg/matrix.hpp"
#include "malg/vector.hpp"
#include "malg/feq.hpp"

#include <type_traits>

/// @brief Sums two matrices.
/// @param a first matrix.
/// @param b second matrix.
/// @returns the resulting matrix.
template <typename T1, typename T2>
auto operator+(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    // Check the sizes.
    assert((a.rows() == b.rows()) && "Matrices has different number of rows.");
    assert((a.cols() == b.cols()) && "Matrices has different number of colmuns.");
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Create the result matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Compute the result.
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = a[i] + b[i];
    return result;
}

/// @brief Sums a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @returns the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] + b;
    return result;
}

/// @brief Sums two matrices.
/// @param a first matrix.
/// @param b second matrix.
/// @returns matrix a.
template <typename T1, typename T2>
auto &operator+=(malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    // Check the sizes.
    assert((a.rows() == b.rows()) && "Matrices has different number of rows.");
    assert((a.cols() == b.cols()) && "Matrices has different number of colmuns.");
    // Compute the result.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] += b[i];
    return a;
}

/// @brief Sums a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @returns matrix a.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+=(malg::MatrixBase<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] += b;
    return a;
}

/// @brief Substraction between two matrices.
/// @param a first matrix.
/// @param b second matrix.
/// @returns the resulting matrix.
template <typename T1, typename T2>
auto operator-(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    // Check the sizes.
    assert((a.rows() == b.rows()) && "Matrices has different number of rows.");
    assert((a.cols() == b.cols()) && "Matrices has different number of colmuns.");
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Create the result matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Compute the result.
    for (std::size_t i = 0; i < result.size(); ++i)
        result[i] = a[i] - b[i];
    return result;
}

/// @brief Substraction between a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @returns the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] - b;
    return result;
}

/// @brief Substraction between two matrices.
/// @param a first matrix.
/// @param b second matrix.
/// @returns matrix a.
template <typename T1, typename T2>
auto &operator-=(malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    // Check the sizes.
    assert((a.rows() == b.rows()) && "Matrices has different number of rows.");
    assert((a.cols() == b.cols()) && "Matrices has different number of colmuns.");
    // Compute the result.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] -= b[i];
    return a;
}

/// @brief Substraction between a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @returns matrix a.
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-=(malg::MatrixBase<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] -= b;
    return a;
}

/// @brief Multiplication between two matrices.
/// @param a first *xM matrix.
/// @param b second Mx* matrix.
/// @returns the resulting matrix.
template <typename T1, typename T2>
auto operator*(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    // Check the sizes.
    assert((a.cols() == b.rows()) &&
           "For matrix multiplication, the number of columns in the first"
           "matrix must be equal to the number of rows in the second matrix.");
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), b.cols());
    // Perform the operation.
    for (std::size_t r = 0; r < a.rows(); r++)
        for (std::size_t c = 0; c < b.cols(); c++)
            for (std::size_t k = 0; k < b.rows(); k++)
                result(r, c) += a(r, k) * b(k, c);
    return result;
}

/// @brief Multiplication between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @returns the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator*(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] * b;
    return result;
}

/// @brief Multiplication between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @returns the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T1> || malg::is_complex_v<T1>, T1>>
inline auto operator*(const T1 &a, const malg::MatrixBase<T2> &b)
{
    return b * a;
}

/// @brief Multiplication between two matrices.
/// @param a first *xM matrix.
/// @param b second Mx* matrix.
/// @returns the resulting matrix.
template <typename T1, typename T2>
auto &operator*=(malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    // Check the sizes.
    assert((a.cols() == b.rows()) &&
           "For matrix multiplication, the number of columns in the first"
           "matrix must be equal to the number of rows in the second matrix.");
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the operation.
    for (std::size_t r = 0; r < a.rows(); r++)
        for (std::size_t c = 0; c < b.cols(); c++)
            for (std::size_t k = 0; k < b.rows(); k++)
                result(r, c) += a(r, k) * b(k, c);
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = result[i];
    return a;
}

/// @brief Multiplication between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @returns the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto &operator*=(malg::MatrixBase<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = a[i] * b;
    return a;
}

/// @brief Division between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @returns the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator/(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] / b;
    return result;
}

/// @brief Division between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @returns the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto &operator/=(malg::MatrixBase<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = a[i] / b;
    return a;
}

/// @brief Equality comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @returns matrix with logical values: true if the element is the same, false
/// otherwise.
template <typename T1, typename T2>
auto operator==(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = feq::approximately_equal(a[i], b[i]);
        else if constexpr (malg::is_complex_v<T1> && malg::is_complex_v<T2>)
            result[i] = feq::approximately_equal(std::real(a[i]), std::real(b[i])) &&
                        feq::approximately_equal(std::imag(a[i]), std::imag(b[i]));
        else
            result[i] = a[i] == b[i];
    }
    return result;
}

/// @brief Equality comparison operator.
/// @param a the first vector.
/// @param b the first vector.
/// @returns vector with logical values: true if the element is the same, false
/// otherwise.
template <typename T1, typename T2>
auto operator==(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    malg::Vector<bool> result(a.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = feq::approximately_equal(a[i], b[i]);
        else if constexpr (malg::is_complex_v<T1> && malg::is_complex_v<T2>)
            result[i] = feq::approximately_equal(std::real(a[i]), std::real(b[i])) &&
                        feq::approximately_equal(std::imag(a[i]), std::imag(b[i]));
        else
            result[i] = a[i] == b[i];
    }
    return result;
}

/// @brief Equality comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @returns matrix with logical values: true if the element is the same, false
/// otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator==(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = feq::approximately_equal(a[i], b);
        else if constexpr (malg::is_complex_v<T1> && malg::is_complex_v<T2>)
            result[i] = feq::approximately_equal(std::real(a[i]), std::real(b)) &&
                        feq::approximately_equal(std::imag(a[i]), std::imag(b));
        else
            result[i] = a[i] == b;
    }
    return result;
}

/// @brief Equality comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @returns matrix with logical values: true if the element is the same, false
/// otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator==(const malg::Vector<T1> &a, const T2 &b)
{
    malg::Vector<bool> result(a.size(), false);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = feq::approximately_equal(a[i], b);
        else if constexpr (malg::is_complex_v<T1> && malg::is_complex_v<T2>)
            result[i] = feq::approximately_equal(std::real(a[i]), std::real(b)) &&
                        feq::approximately_equal(std::imag(a[i]), std::imag(b));
        else
            result[i] = a[i] == b;
    }
    return result;
}

/// @brief Inequality comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @returns matrix with logical values: true if the element is different, false
/// otherwise.
template <typename T1, typename T2>
auto operator!=(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = !feq::approximately_equal(a[i], b[i]);
        else if constexpr (malg::is_complex_v<T1> && malg::is_complex_v<T2>)
            result[i] = !feq::approximately_equal(std::real(a[i]), std::real(b[i])) ||
                        !feq::approximately_equal(std::imag(a[i]), std::imag(b[i]));
        else
            result[i] = a[i] != b[i];
    }
    return result;
}

/// @brief Inequality comparison operator.
/// @param a the first vector.
/// @param b the first vector.
/// @returns vector with logical values: true if the element is different, false
/// otherwise.
template <typename T1, typename T2>
auto operator!=(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    malg::Vector<bool> result(a.size(), false);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = !feq::approximately_equal(a[i], b[i]);
        else if constexpr (malg::is_complex_v<T1> && malg::is_complex_v<T2>)
            result[i] = !feq::approximately_equal(std::real(a[i]), std::real(b[i])) ||
                        !feq::approximately_equal(std::imag(a[i]), std::imag(b[i]));
        else
            result[i] = a[i] != b[i];
    }
    return result;
}

/// @brief Inequality comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @returns matrix with logical values: true if the element is different, false
/// otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator!=(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = !feq::approximately_equal(a[i], b);
        else if constexpr (malg::is_complex_v<T1> && malg::is_complex_v<T2>)
            result[i] = !feq::approximately_equal(std::real(a[i]), std::real(b)) ||
                        !feq::approximately_equal(std::imag(a[i]), std::imag(b));
        else
            result[i] = a[i] != b;
    }
    return result;
}

/// @brief Inequality comparison operator.
/// @param a the vector.
/// @param b the scalar.
/// @returns vector with logical values: true if the element is different, false
/// otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator!=(const malg::Vector<T1> &a, const T2 &b)
{
    malg::Vector<bool> result(a.size(), false);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = !feq::approximately_equal(a[i], b);
        else if constexpr (malg::is_complex_v<T1> && malg::is_complex_v<T2>)
            result[i] = !feq::approximately_equal(std::real(a[i]), std::real(b)) ||
                        !feq::approximately_equal(std::imag(a[i]), std::imag(b));
        else
            result[i] = a[i] != b;
    }
    return result;
}

/// @brief Greather-than comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @returns matrix with logical values: true if the first element greather-than
/// the second, false otherwise.
template <typename T1, typename T2>
auto operator>(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] > b[i];
    return result;
}

/// @brief Greather-than equal comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @returns matrix with logical values: true if the first element greather-than
/// equal the second, false otherwise.
template <typename T1, typename T2>
auto operator>=(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] >= b[i];
    return result;
}

/// @brief Greather-than comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @returns matrix with logical values: true if the first element greather-than
/// the second, false otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator>(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] > b;
    return result;
}

/// @brief Greather-than equal comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @returns matrix with logical values: true if the first element greather-than
/// equal the second, false otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator>=(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] >= b;
    return result;
}

/// @brief Lesser-than comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @returns matrix with logical values: true if the first element lesser-than
/// the second, false otherwise.
template <typename T1, typename T2>
auto operator<(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] < b[i];
    return result;
}

/// @brief Lesser-than equal comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @returns matrix with logical values: true if the first element lesser-than
/// equal the second, false otherwise.
template <typename T1, typename T2>
auto operator<=(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] <= b[i];
    return result;
}

/// @brief Lesser-than comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @returns matrix with logical values: true if the first element lesser-than
/// the second, false otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator<(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] < b;
    return result;
}

/// @brief Lesser-than equal comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @returns matrix with logical values: true if the first element lesser-than
/// equal the second, false otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator<=(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] <= b;
    return result;
}

// ============================================================================
// VECTOR
// ============================================================================

// ====================================
// SUM
// ====================================

/// @brief Element-wise sum of two vector.
/// @param a first vector.
/// @param b second vector.
/// @returns a vector containing the sum of the two vectors.
template <typename T1, typename T2>
inline constexpr auto operator+(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] + b[i];
    return result;
}

/// @brief Sums a scalar to a vector.
/// @param a the vector.
/// @param b the scalar.
/// @returns a vector containing the sum of the scalar and the vector.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] + b;
    return result;
}

/// @brief Element-wise sum of two vector, stores the solution in the left-hand side one.
/// @param a first vector.
/// @param b second vector.
/// @returns a reference to the left-hand side vector, containing the sum of the two vectors.
template <typename T1, typename T2>
inline auto &operator+=(malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = a[i] + b[i];
    return a;
}

/// @brief Sums a scalar to a vector, stores the solution in the left-hand side one.
/// @param a the vector.
/// @param b the scalar.
/// @returns a reference to the left-hand side vector, containing the sum of the scalar and the vector.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+=(const malg::Vector<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = a[i] + b;
    return a;
}

// ====================================
// SUB
// ====================================

/// @brief Element-wise subtraction of two vector.
/// @param a first vector.
/// @param b second vector.
/// @returns a vector containing the subtraction of the second vector from the first.
template <typename T1, typename T2>
inline auto operator-(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] - b[i];
    return result;
}

/// @brief Subtracts a scalar from a vector.
/// @param a the vector.
/// @param b the scalar.
/// @returns a vector containing the subtraction of the scalar from the vector.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] - b;
    return result;
}

/// @brief Element-wise subtraction of two vector, stores the solution in the left-hand side one.
/// @param a first vector.
/// @param b second vector.
/// @returns a reference to the left-hand side vector, containing the subtraction of the second vector from the first.
template <typename T1, typename T2>
inline auto &operator-=(malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = a[i] - b[i];
    return a;
}

/// @brief Subtracts a scalar from a vector, stores the solution in the left-hand side one.
/// @param a the vector.
/// @param b the scalar.
/// @returns a reference to the left-hand side vector, containing the subtraction of the scalar from the vector.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-=(const malg::Vector<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = a[i] - b;
    return a;
}

// ====================================
// MUL
// ====================================

/// @brief Element-wise product of two vectors.
/// @param a first vector.
/// @param b second vector.
/// @returns a vector containing the element-wise multiplicaiton of the two vectors.
template <typename T1, typename T2>
inline auto operator*(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] * b[i];
    return result;
}

/// @brief Multiplies a scalar and a vector.
/// @param a the vector.
/// @param b the scalar.
/// @returns a vector containing the multiplication of the scalar for the vector.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator*(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] * b;
    return result;
}

/// @brief Multiplies a scalar and a vector.
/// @param a the vector.
/// @param b the scalar.
/// @returns a vector containing the multiplication of the scalar for the vector.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline auto operator*(const T1 &a, const malg::Vector<T2> &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a * b[i];
    return result;
}

/// @brief Multiplies a scalar and a vector, stores the solution in the input vector.
/// @param a the vector.
/// @param b the scalar.
/// @returns a reference to the input vector, containing the multiplication of the scalar and the vector.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator*=(const malg::Vector<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = a[i] * b;
    return a;
}

// ====================================
// DIV
// ====================================

/// @brief Divides the elements of a vector by a scalar.
/// @param a the vector.
/// @param b the scalar.
/// @returns a vector containing the division of the elements of the vector by the scalar.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator/(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = a[i] / b;
    return result;
}

/// @brief Divides the elements of a vector by a scalar, stores the solution in the input vector.
/// @param a the vector.
/// @param b the scalar.
/// @returns a reference to the input vector, containing the division of the elements of the vector by the scalar.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator/=(const malg::Vector<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        a[i] = a[i] / b;
    return a;
}

// ========================================================
// SPECIFIC FUNCTIONS
// ========================================================

namespace malg
{

/// @brief Computes the dot product between two vector.
/// @param a the first vector.
/// @param b the second vector.
/// @returns the scalar value resulting from the dot product.
template <typename T1, typename T2>
inline constexpr auto dot(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    data_type_t result = data_type_t(0.);
    // Perform the operation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result += a[i] * b[i];
    return result;
}

/// @brief Computes the dot product between a matrix and a vector.
/// @param A the matrix.
/// @param b the vector.
/// @returns a vector with the same **rows** of *A*, resulting from the dot product.
template <typename T1, typename T2>
inline constexpr auto dot(const malg::MatrixBase<T1> &A, const malg::Vector<T2> &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Create the result vector.
    malg::Vector<data_type_t> result(A.rows());
    // Perform the computation.
    for (std::size_t r = 0; r < A.rows(); ++r)
        for (std::size_t c = 0; c < A.cols(); ++c)
            result[r] += A(r, c) * b[c];
    return result;
}

/// @brief Calls the function **fun** on the elements of **a** and stores the result in a vector.
/// @param a the vector.
/// @param fun the function.
/// @returns a vector containing the results of calling **fun** on the elements of **a**.
template <typename T, typename Function>
inline auto element_wise_function(const malg::Vector<T> &a, Function fun)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T, decltype(fun(std::declval<T>()))>>;
    // Create the resulting vector.
    malg::Vector<data_type_t> result(a);
    // Apply the function to all the elements.
    std::for_each(result.begin(), result.end(), fun);
    // Return the vector.
    return result;
}

/// @brief Calls the function **fun** on the elements of **A** and stores the result in a matrix.
/// @param A the matrix.
/// @param fun the function.
/// @returns a matrix containing the results of calling **fun** on the elements of **A**.
template <typename T, typename Function>
inline auto element_wise_function(const malg::MatrixBase<T> &A, Function fun)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T, decltype(fun(std::declval<T>()))>>;
    // Create the resulting matrix.
    malg::Matrix<data_type_t> result(A.rows(), A.cols());
    // Apply the function to all the elements.
    std::for_each(result.begin(), result.end(), fun);
    // Return the matrix.
    return result;
}

/// @brief Calls the binary function **fun** on the elements of **a** and **b** and stores the result in a vector.
/// @param a the first vector.
/// @param b the second vector.
/// @param fun the function.
/// @returns a vector containing the results of calling **fun** on the elements of **a** and **b**.
template <typename T1, typename T2, typename Function>
inline auto element_wise_binary_function(const malg::Vector<T1> &a, const malg::Vector<T2> &b, Function fun)
{
    assert(a.size() == b.size());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_3_t<T1, T2, decltype(fun(std::declval<T1>(), std::declval<T2>()))>>;
    // Create the resulting vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the computation.
    for (std::size_t i = 0; i < a.size(); ++i)
        result[i] = fun(a[i], b[i]);
    // Return the vector.
    return result;
}

/// @brief Calls the function **fun** on the elements of **A** and **A** and stores the result in a matrix.
/// @param A the first matrix.
/// @param B the second matrix.
/// @param fun the function.
/// @returns a matrix containing the results of calling **fun** on the elements of **A** and **B**.
template <typename T1, typename T2, typename Function>
inline auto element_wise_binary_function(const malg::MatrixBase<T1> &A, const malg::MatrixBase<T2> &B, Function fun)
{
    assert(A.rows() == B.rows());
    assert(A.cols() == B.cols());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_3_t<T1, T2, decltype(fun(std::declval<T1>(), std::declval<T2>()))>>;
    // Create the resulting matrix.
    malg::Matrix<data_type_t> result(A.rows(), A.cols());
    // Perform the computation.
    for (std::size_t i = 0; i < A.size(); ++i)
        result[i] = fun(A[i], B[i]);
    // Return the matrix.
    return result;
}

/// @brief Extracts the real part of the values of a.
/// @param a the matrix.
/// @returns the real part of the values of a.
template <typename T>
auto real(const malg::MatrixBase<std::complex<T>> &a)
{
    return element_wise_function(a, [](const std::complex<T> &value) { return value.real(); });
}

/// @brief Extracts the real part of the values of a.
/// @param a the matrix.
/// @returns the real part of the values of a.
template <typename T>
auto real(const malg::Vector<std::complex<T>> &a)
{
    return element_wise_function(a, [](const std::complex<T> &value) { return value.real(); });
}

/// @brief Extracts the imaginary part of the values of a.
/// @param a the matrix.
/// @returns imaginary part of the values of a.
template <typename T>
auto imag(const malg::MatrixBase<std::complex<T>> &a)
{
    return element_wise_function(a, [](const std::complex<T> &value) { return value.imag(); });
}

/// @brief Extracts the imaginary part of the values of a.
/// @param a the matrix.
/// @returns imaginary part of the values of a.
template <typename T>
auto imag(const malg::Vector<std::complex<T>> &a)
{
    return element_wise_function(a, [](const std::complex<T> &value) { return value.imag(); });
}

/// @brief Transforms each element of a to its absolute value.
/// @param a the input matrix.
/// @returns the same matrix but its absolute values.
template <typename T>
auto abs(const malg::MatrixBase<T> &a)
{
    return element_wise_function(a, [](const T &value) { return std::abs(value); });
}

/// @brief Transforms each element of a to its absolute value.
/// @param a the input vector.
/// @returns the same vector but its absolute values.
template <typename T>
auto abs(const malg::Vector<T> &a)
{
    return element_wise_function(a, [](const T &value) { return std::abs(value); });
}

/// @brief Returns a matrix containing the sign of the values in the input matrix.
/// @param a the input matrix.
/// @returns a matrix containing -1 and +1 based on the signs of the values of a.
template <typename T>
auto sign(const malg::MatrixBase<T> &a)
{
    if constexpr (malg::is_complex_v<T>)
        return element_wise_function(a, [](const T &value) { return T(std::real(value) > 0 ? +1 : -1, std::imag(value) > 0 ? +1 : -1); });
    else
        return element_wise_function(a, [](const T &value) { return value > 0 ? +1 : -1; });
}

/// @brief Returns a matrix containing the sign of the values in the input matrix.
/// @param a the input matrix.
/// @returns a matrix containing -1 and +1 based on the signs of the values of a.
template <typename T>
auto sign(const malg::Vector<T> &a)
{
    if constexpr (malg::is_complex_v<T>)
        return element_wise_function(a, [](const T &value) { return T(std::real(value) > 0 ? +1 : -1, std::imag(value) > 0 ? +1 : -1); });
    else
        return element_wise_function(a, [](const T &value) { return value > 0 ? +1 : -1; });
}

/// @brief Vector projection of **a** onto **b**.
/// @param a the first vector.
/// @param b the second vector.
/// @returns the vector projection.
template <typename T1, typename T2>
inline auto projection(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    return b * (malg::dot(a, b) / malg::dot(b, b));
}

/// @brief Computes the vector length (or magnitude) by using the Pythagorean theorem, of a given column **c** of matrix **A**.
/// @param A the matrix.
/// @param c the column for which we compute the magnitude.
/// @param r_start the starting row.
/// @param r_end the ending row.
/// @returns the mangnitude of the column vector.
template <typename T>
inline auto vector_length(malg::MatrixBase<T> &A, std::size_t c, std::size_t r_start = 0, std::size_t r_end = -1)
{
    std::remove_const_t<malg::is_complex_t<T>> accum = 0;
    for (std::size_t r = r_start; r < std::min(r_end, A.rows()); r++) {
        if constexpr (malg::is_complex_v<T>)
            accum += std::norm(A(r, c));
        else
            accum += A(r, c) * A(r, c);
    }
    return std::sqrt(accum);
}

/// @brief Computes the linear combination of matrices.
/// @param A the first matrix.
/// @param a the first scalar.
/// @param B the second matrix.
/// @param b the second scalar.
/// @returns a matrix containing the linear combination.
template <typename T1, typename T2, typename T3, typename T4>
auto linear_combination(const malg::MatrixBase<T1> &A, const T2 &a, const malg::MatrixBase<T3> &B, const T4 &b)
{
    return element_wise_binary_function(A, B, [&](const T1 lhs, const T1 rhs) { lhs *a + rhs *b; });
}

/// @brief Sums the element of the vector.
/// @param v the input vector.
/// @returns the sum of the elements.
template <typename T>
inline auto sum(const Vector<T> &v)
{
    std::remove_const_t<T> s = 0;
    for (std::size_t i = 0; i < v.size(); ++i)
        s += v[i];
    return s;
}

/// @brief Compute the trace of A, i.e., the sum of the elements along the main diagonal.
/// @param A the input matrix.
/// @returns the sum of the diagonal elements.
template <typename T>
inline auto trace(const MatrixBase<T> &A)
{
#if 0
    return sum(diag(A));
#else
    std::remove_const_t<T> result = 0;
    for (std::size_t i = 0; i < A.rows(); ++i)
        result += A(i, i);
    return result;
#endif
}

/// @brief Returns the minimum value inside the vector, and its position.
/// @param v the input vector.
/// @returns a pair containing the value and its position inside the vector.
template <typename T>
inline auto min(const Vector<T> &v)
{
    using data_type_t = std::remove_const_t<T>;
    if (v.empty())
        return std::make_pair(data_type_t(0.), std::size_t(0));
    std::size_t min_val_pos = 0;
    for (std::size_t i = 1; i < v.size(); ++i)
        if (v[i] < v[min_val_pos])
            min_val_pos = i;
    return std::make_pair(v[min_val_pos], min_val_pos);
}

/// @brief Returns a vector with the minimum value for each column of the input matrix.
/// @param a the input matrix.
/// @returns the minimum values of the columns of a.
template <typename T>
inline auto min(const Matrix<T> &a)
{
    using data_type_t = std::remove_const_t<T>;
    if (a.empty())
        return Vector<data_type_t>();
    Vector<data_type_t> values(a.cols(), data_type_t(0.));
    for (std::size_t c = 0, r, min_val_pos; c < a.cols(); ++c) {
        for (min_val_pos = 0, r = 1; r < a.rows(); ++r)
            if (a(r, c) < a(min_val_pos, c))
                min_val_pos = r;
        values[c] = a(min_val_pos, c);
    }
    return values;
}

/// @brief Returns the maximum value inside the vector, and its position.
/// @param v the input vector.
/// @returns a pair containing the value and its position inside the vector.
template <typename T>
inline auto max(const Vector<T> &v)
{
    using data_type_t = std::remove_const_t<T>;
    if (v.empty())
        return std::make_pair(data_type_t(0.), std::size_t(0));
    std::size_t max_val_pos = 0;
    for (std::size_t i = 1; i < v.size(); ++i)
        if (v[i] > v[max_val_pos])
            max_val_pos = i;
    return std::make_pair(v[max_val_pos], max_val_pos);
}

/// @brief Returns a vector with the maximum value for each column of the input matrix.
/// @param a the input matrix.
/// @returns the maximum values of the columns of a.
template <typename T>
inline auto max(const Matrix<T> &a)
{
    using data_type_t = std::remove_const_t<T>;
    if (a.empty())
        return Vector<data_type_t>();
    Vector<data_type_t> values(a.cols(), data_type_t(0.));
    for (std::size_t c = 0, r, max_val_pos; c < a.cols(); ++c) {
        for (max_val_pos = 0, r = 1; r < a.rows(); ++r)
            if (a(r, c) > a(max_val_pos, c))
                max_val_pos = r;
        values[c] = a(max_val_pos, c);
    }
    return values;
}

/// @brief Checks if all elements are non-zero (or true).
/// @param a the input matrix.
/// @returns true if all elements are non-zero (or true).
/// @returns false if even one element is zero (or false).
template <typename T>
inline bool all(const MatrixBase<T> &a)
{
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] == 0)
            return false;
    return true;
}

/// @brief Checks if all elements are non-zero (or true).
/// @param v the input vector.
/// @returns true if all elements are non-zero (or true).
/// @returns false if even one element is zero (or false).
template <typename T>
inline bool all(const Vector<T> &v)
{
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] == 0)
            return false;
    return true;
}

/// @brief Checks if at least one element is non-zero (or true).
/// @param a the input matrix.
/// @returns true if at least one element is non-zero (or true).
/// @returns false if all elements are zero (or false).
template <typename T>
inline bool any(const MatrixBase<T> &a)
{
    for (std::size_t i = 0; i < a.size(); ++i)
        if (a[i] != 0)
            return true;
    return false;
}

/// @brief Checks if at least one element is non-zero (or true).
/// @param v the input vector.
/// @returns true if at least one element is non-zero (or true).
/// @returns false if all elements are zero (or false).
template <typename T>
inline bool any(const Vector<T> &v)
{
    for (std::size_t i = 0; i < v.size(); ++i)
        if (v[i] != 0)
            return true;
    return false;
}

/// @brief Computes the square norm of the vector.
/// @param v the vector.
/// @return the square norm of the vector.
template <typename T>
inline auto square_norm(const malg::Vector<T> &v)
{
    std::remove_const_t<malg::is_complex_t<T>> accum = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if constexpr (malg::is_complex_v<T>)
            accum += std::norm(v[i]);
        else
            accum += v[i] * v[i];
    }
    return std::sqrt(accum);
}

/// @brief The Frobenius norm of a matrix.
/// @param A the input matrix.
/// @returns the norm.
template <typename T>
inline auto square_norm(const malg::MatrixBase<T> &A)
{
    std::remove_const_t<malg::is_complex_t<T>> accum = 0;
    // Compute the sum of squares of the elements of the given matrix.
    for (std::size_t r = 0; r < A.rows(); ++r) {
        for (std::size_t c = 0; c < A.cols(); ++c) {
            if constexpr (malg::is_complex_v<T>)
                accum += std::norm(A(r, c));
            else
                accum += A(r, c) * A(r, c);
        }
    }
    // Return the square root of the sum of squares.
    return std::sqrt(accum);
}

/// @brief Computes the infinity norm of a vector, i.e., largest magnitude among each element of a vector.
/// @param v the input vector.
/// @returns the infinity norm.
template <typename T>
inline auto infinity_norm(const malg::Vector<T> &v)
{
    using data_type_t = std::remove_const_t<T>;
    if (v.empty())
        return data_type_t(0.);
    data_type_t max = std::abs(v[0]), tmp;
    for (std::size_t i = 1; i < v.size(); ++i) {
        tmp = std::abs(v[i]);
        if (max < tmp) {
            max = tmp;
        }
    }
    return max;
}

/// @brief Computes the infinity norm of a matrix, i.e., largest infinity norm among the rows of the matrix.
/// @param A the input matrix.
/// @returns the infinity norm.
template <typename T>
inline auto infinity_norm(const malg::MatrixBase<T> &A)
{
    using data_type_t = std::remove_const_t<malg::is_complex_t<T>>;
    if (A.empty())
        return data_type_t(0.);
    data_type_t max{}, accum{};
    for (std::size_t r = 0; r < A.rows(); ++r) {
        accum = 0.;
        for (std::size_t c = 0; c < A.cols(); ++c)
            accum += std::abs(A(r, c));
        max = std::max(max, accum);
    }
    return max;
}

/// @brief The Euclidean norm of the lower leading diagonal of a square matrix.
/// @param A the input matrix.
/// @returns the square root of the sum of all the squares.
template <typename T>
auto sub_norm(const malg::MatrixBase<T> &A)
{
    using data_type_t = std::remove_const_t<malg::is_complex_t<T>>;
    if (A.empty())
        return data_type_t(0.);
    assert(malg::utility::is_square(A));
    data_type_t accum = 0;
    for (std::size_t r = 1; r < A.rows(); ++r) {
        for (std::size_t c = 0; c < r; ++c) {
            if constexpr (malg::is_complex_v<T>)
                accum += std::norm(A(r, c));
            else
                accum += A(r, c) * A(r, c);
        }
    }
    return std::sqrt(accum);
}

/// @brief Scale down matrix A by a power of 2, such that norm(A) < 1.
/// @param A the input matrix.
/// @returns the square root of the sum of all the squares.
template <typename T>
auto log2_ceil(const malg::MatrixBase<T> &A)
{
    using data_type_t      = std::remove_const_t<malg::is_complex_t<T>>;
    std::size_t iterations = 0;
    data_type_t scale      = 1.0;
    const auto norm        = malg::infinity_norm(A);
    while ((norm * scale) > 1.0) {
        scale *= 0.5;
        ++iterations;
    }
    return std::make_pair(iterations, scale);
}

} // namespace malg
