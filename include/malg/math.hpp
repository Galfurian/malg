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
/// @return the resulting matrix.
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
    for (unsigned i = 0; i < result.size(); ++i)
        result[i] = a[i] + b[i];
    return result;
}

/// @brief Sums a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] + b;
    return result;
}

/// @brief Sums two matrices.
/// @param a first matrix.
/// @param b second matrix.
/// @return matrix a.
template <typename T1, typename T2>
auto &operator+=(malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    // Check the sizes.
    assert((a.rows() == b.rows()) && "Matrices has different number of rows.");
    assert((a.cols() == b.cols()) && "Matrices has different number of colmuns.");
    // Compute the result.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] += b[i];
    return a;
}

/// @brief Sums a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @return matrix a.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+=(malg::MatrixBase<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] += b;
    return a;
}

/// @brief Substraction between two matrices.
/// @param a first matrix.
/// @param b second matrix.
/// @return the resulting matrix.
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
    for (unsigned i = 0; i < result.size(); ++i)
        result[i] = a[i] - b[i];
    return result;
}

/// @brief Substraction between a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] - b;
    return result;
}

/// @brief Substraction between two matrices.
/// @param a first matrix.
/// @param b second matrix.
/// @return matrix a.
template <typename T1, typename T2>
auto &operator-=(malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    // Check the sizes.
    assert((a.rows() == b.rows()) && "Matrices has different number of rows.");
    assert((a.cols() == b.cols()) && "Matrices has different number of colmuns.");
    // Compute the result.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] -= b[i];
    return a;
}

/// @brief Substraction between a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @return matrix a.
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-=(malg::MatrixBase<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] -= b;
    return a;
}

/// @brief Multiplication between two matrices.
/// @param a first *xM matrix.
/// @param b second Mx* matrix.
/// @return the resulting matrix.
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
    for (unsigned r = 0; r < a.rows(); r++)
        for (unsigned c = 0; c < b.cols(); c++)
            for (unsigned k = 0; k < b.rows(); k++)
                result(r, c) += a(r, k) * b(k, c);
    return result;
}

/// @brief Multiplication between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator*(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] * b;
    return result;
}

/// @brief Multiplication between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T1> || malg::is_complex_v<T1>, T1>>
inline auto operator*(const T1 &a, const malg::MatrixBase<T2> &b)
{
    return b * a;
}

/// @brief Multiplication between two matrices.
/// @param a first *xM matrix.
/// @param b second Mx* matrix.
/// @return the resulting matrix.
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
    for (unsigned r = 0; r < a.rows(); r++)
        for (unsigned c = 0; c < b.cols(); c++)
            for (unsigned k = 0; k < b.rows(); k++)
                result(r, c) += a(r, k) * b(k, c);
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = result[i];
    return a;
}

/// @brief Multiplication between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto &operator*=(malg::MatrixBase<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] * b;
    return a;
}

/// @brief Division between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator/(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] / b;
    return result;
}

/// @brief Division between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto &operator/=(malg::MatrixBase<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] / b;
    return a;
}

/// @brief Equality comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @return matrix with logical values: true if the element is the same, false
/// otherwise.
template <typename T1, typename T2>
auto operator==(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
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
/// @return vector with logical values: true if the element is the same, false
/// otherwise.
template <typename T1, typename T2>
auto operator==(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    malg::Vector<bool> result(a.size());
    for (unsigned i = 0; i < a.size(); ++i) {
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
/// @return matrix with logical values: true if the element is the same, false
/// otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator==(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
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
/// @return matrix with logical values: true if the element is the same, false
/// otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator==(const malg::Vector<T1> &a, const T2 &b)
{
    malg::Vector<bool> result(a.size(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
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
/// @return matrix with logical values: true if the element is different, false
/// otherwise.
template <typename T1, typename T2>
auto operator!=(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
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
/// @return vector with logical values: true if the element is different, false
/// otherwise.
template <typename T1, typename T2>
auto operator!=(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    malg::Vector<bool> result(a.size(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
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
/// @return matrix with logical values: true if the element is different, false
/// otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator!=(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
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
/// @return vector with logical values: true if the element is different, false
/// otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator!=(const malg::Vector<T1> &a, const T2 &b)
{
    malg::Vector<bool> result(a.size(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
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
/// @return matrix with logical values: true if the first element greather-than
/// the second, false otherwise.
template <typename T1, typename T2>
auto operator>(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] > b[i];
    return result;
}

/// @brief Greather-than equal comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @return matrix with logical values: true if the first element greather-than
/// equal the second, false otherwise.
template <typename T1, typename T2>
auto operator>=(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] >= b[i];
    return result;
}

/// @brief Greather-than comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @return matrix with logical values: true if the first element greather-than
/// the second, false otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator>(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] > b;
    return result;
}

/// @brief Greather-than equal comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @return matrix with logical values: true if the first element greather-than
/// equal the second, false otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator>=(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] >= b;
    return result;
}

/// @brief Lesser-than comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @return matrix with logical values: true if the first element lesser-than
/// the second, false otherwise.
template <typename T1, typename T2>
auto operator<(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] < b[i];
    return result;
}

/// @brief Lesser-than equal comparison operator.
/// @param a the first matrix.
/// @param b the first matrix.
/// @return matrix with logical values: true if the first element lesser-than
/// equal the second, false otherwise.
template <typename T1, typename T2>
auto operator<=(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] <= b[i];
    return result;
}

/// @brief Lesser-than comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @return matrix with logical values: true if the first element lesser-than
/// the second, false otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator<(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] < b;
    return result;
}

/// @brief Lesser-than equal comparison operator.
/// @param a the matrix.
/// @param b the scalar.
/// @return matrix with logical values: true if the first element lesser-than
/// equal the second, false otherwise.
template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator<=(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] <= b;
    return result;
}

// ========================================================
// VECTOR
// ========================================================

template <typename T1, typename T2>
inline auto operator+(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] + b[i];
    return result;
}

template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] + b;
    return result;
}

template <typename T1, typename T2>
inline auto &operator+=(malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] + b[i];
    return a;
}

template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+=(const malg::Vector<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] + b;
    return a;
}

template <typename T1, typename T2>
inline auto operator-(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] - b[i];
    return result;
}

template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] - b;
    return result;
}

template <typename T1, typename T2>
inline auto &operator-=(malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] - b[i];
    return a;
}

template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-=(const malg::Vector<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] - b;
    return a;
}

template <typename T1, typename T2>
inline auto operator*(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    data_type_t result = data_type_t(0.);
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result += a[i] * b[i];
    return result;
}

template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator*(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] * b;
    return result;
}

template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline auto operator*(const T1 &a, const malg::Vector<T2> &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    data_type_t result = data_type_t(0.);
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result += a * b[i];
    return result;
}

template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator*=(const malg::Vector<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] * b;
    return a;
}

template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator/(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Declare the output vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] / b;
    return result;
}

template <typename T1, typename T2, typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator/=(const malg::Vector<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] / b;
    return a;
}

// ========================================================
// SPECIFIC FUNCTIONS
// ========================================================

namespace malg
{

template <typename T1, typename T2>
inline auto dot(const malg::MatrixBase<T1> &a, const malg::Vector<T2> &b)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T1, T2>>;
    // Create the result vector.
    malg::Vector<data_type_t> result(a.rows());
    // Perform the computation.
    for (unsigned r = 0; r < a.rows(); ++r)
        for (unsigned c = 0; c < a.cols(); ++c)
            result[r] += a(r, c) * b[c];
    return result;
}

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

template <typename T, typename Function>
inline auto element_wise_function(const malg::MatrixBase<T> &a, Function fun)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T, decltype(fun(std::declval<T>()))>>;
    // Create the resulting matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Apply the function to all the elements.
    std::for_each(result.begin(), result.end(), fun);
    // Return the matrix.
    return result;
}

template <typename T1, typename T2, typename Function>
inline auto element_wise_binary_function(const malg::Vector<T1> &a, const malg::Vector<T2> &b, Function fun)
{
    assert(a.size() == b.size());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_3_t<T1, T2, decltype(fun(std::declval<T1>(), std::declval<T2>()))>>;
    // Create the resulting vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the computation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = fun(a[i], b[i]);
    // Return the vector.
    return result;
}

template <typename T1, typename T2, typename Function>
inline auto element_wise_binary_function(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b, Function fun)
{
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_3_t<T1, T2, decltype(fun(std::declval<T1>(), std::declval<T2>()))>>;
    // Create the resulting matrix.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the computation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = fun(a[i], b[i]);
    // Return the matrix.
    return result;
}

template <typename T1, typename T2>
inline auto element_wise_product(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    return element_wise_binary_function(a, b, [](const T1 lhs, const T2 rhs) { return lhs * rhs; });
}

template <typename T1, typename T2>
inline auto element_wise_product(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    return element_wise_binary_function(a, b, [](const T1 lhs, const T2 rhs) { return lhs * rhs; });
}

template <typename T1, typename T2>
inline auto element_wise_div(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    return element_wise_binary_function(a, b, [](const T1 lhs, const T2 rhs) { return lhs / rhs; });
}

template <typename T1, typename T2>
inline auto element_wise_div(const malg::MatrixBase<T1> &a, const malg::MatrixBase<T2> &b)
{
    return element_wise_binary_function(a, b, [](const T1 lhs, const T2 rhs) { return lhs / rhs; });
}

/// @brief Extracts the real part of the values of a.
/// @param a the matrix.
/// @return the real part of the values of a.
template <typename T>
auto real(const malg::MatrixBase<std::complex<T>> &a)
{
    return element_wise_function(a, [](const std::complex<T> &value) { return value.real(); });
}

/// @brief Extracts the real part of the values of a.
/// @param a the matrix.
/// @return the real part of the values of a.
template <typename T>
auto real(const malg::Vector<std::complex<T>> &a)
{
    return element_wise_function(a, [](const std::complex<T> &value) { return value.real(); });
}

/// @brief Extracts the imaginary part of the values of a.
/// @param a the matrix.
/// @return imaginary part of the values of a.
template <typename T>
auto imag(const malg::MatrixBase<std::complex<T>> &a)
{
    return element_wise_function(a, [](const std::complex<T> &value) { return value.imag(); });
}

/// @brief Extracts the imaginary part of the values of a.
/// @param a the matrix.
/// @return imaginary part of the values of a.
template <typename T>
auto imag(const malg::Vector<std::complex<T>> &a)
{
    return element_wise_function(a, [](const std::complex<T> &value) { return value.imag(); });
}

/// @brief Transforms each element of a to its absolute value.
/// @param a the input matrix.
/// @return the same matrix but its absolute values.
template <typename T>
auto abs(const malg::MatrixBase<T> &a)
{
    return element_wise_function(a, [](const T &value) { return std::abs(value); });
}

/// @brief Transforms each element of a to its absolute value.
/// @param a the input vector.
/// @return the same vector but its absolute values.
template <typename T>
auto abs(const malg::Vector<T> &a)
{
    return element_wise_function(a, [](const T &value) { return std::abs(value); });
}

/// @brief Returns a matrix containing the sign of the values in the input matrix.
/// @param a the input matrix.
/// @return a matrix containing -1 and +1 based on the signs of the values of a.
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
/// @return a matrix containing -1 and +1 based on the signs of the values of a.
template <typename T>
auto sign(const malg::Vector<T> &a)
{
    if constexpr (malg::is_complex_v<T>)
        return element_wise_function(a, [](const T &value) { return T(std::real(value) > 0 ? +1 : -1, std::imag(value) > 0 ? +1 : -1); });
    else
        return element_wise_function(a, [](const T &value) { return value > 0 ? +1 : -1; });
}

template <typename T1, typename T2>
inline auto projection(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    return b * (malg::dot(a, b) / malg::dot(b, b));
}

template <typename T>
inline auto norm(const malg::Vector<T> &v)
{
    std::remove_const_t<is_complex_t<T>> accum = 0;
    for (unsigned i = 0; i < v.size(); ++i) {
        if constexpr (malg::is_complex_v<T>)
            accum += std::norm(v[i]);
        else
            accum += v[i] * v[i];
    }
    return std::sqrt(accum);
}

/// @brief The Euclidean norm of a matrix, the square root of the sum of all the
///        squares.
/// @param A the input matrix.
/// @return the norm.
template <typename T>
inline auto norm(const malg::MatrixBase<T> &A)
{
    std::remove_const_t<is_complex_t<T>> accum = 0;
    for (unsigned r = 0; r < A.rows(); ++r) {
        for (unsigned c = 0; c < A.cols(); ++c) {
            if constexpr (malg::is_complex_v<T>)
                accum += std::norm(A(r, c));
            else
                accum += A(r, c) * A(r, c);
        }
    }
    return std::sqrt(accum);
}

/// @brief The Euclidean norm of the lower leading diagonal of a square matrix.
/// @param A the input matrix.
/// @return the square root of the sum of all the squares.
template <typename T>
auto sub_norm(const malg::MatrixBase<T> &A)
{
    assert(malg::utility::is_square(A));
    std::remove_const_t<is_complex_t<T>> accum = 0;
    for (unsigned r = 1; r < A.rows(); ++r) {
        for (unsigned c = 0; c < r; ++c) {
            if constexpr (malg::is_complex_v<T>)
                accum += std::norm(A(r, c));
            else
                accum += A(r, c) * A(r, c);
        }
    }
    return std::sqrt(accum);
}

template <typename T>
inline auto vector_length(malg::MatrixBase<T> &A, unsigned c, unsigned r_start = 0)
{
    std::remove_const_t<is_complex_t<T>> accum = 0;
    for (unsigned r = r_start; r < A.rows(); r++) {
        if constexpr (malg::is_complex_v<T>)
            accum += std::norm(A(r, c));
        else
            accum += A(r, c) * A(r, c);
    }
    return std::sqrt(accum);
}

// Linear combination of matrices.
template <typename T1, typename T2, typename T3, typename T4>
auto linear_combination(
    const malg::MatrixBase<T1> &A,
    const T2 &a,
    const malg::MatrixBase<T3> &B,
    const T4 &b)
{
    assert((B.cols() == A.rows()) && (B.rows() == A.cols()));
    auto C = A * a;
    for (unsigned r = 0; r < A.rows(); ++r)
        for (unsigned c = 0; c < A.cols(); ++c)
            C(r, c) += b * B(r, c);
    return C;
}

/// @brief Extracts the diagonal elements from the matrix.
/// @param matrix the matrix.
/// @return the diagonal elements.
template <typename T>
inline auto sum(const Vector<T> &v)
{
    std::remove_const_t<T> s = 0;
    for (unsigned i = 0; i < v.size(); ++i)
        s += v[i];
    return s;
}

/// @brief Compute the trace of A, i.e., the sum of the elements along the main diagonal.
/// @param A the input matrix.
/// @return the sum of the diagonal elements.
template <typename T>
inline auto trace(const MatrixBase<T> &A)
{
#if 0
    return sum(diag(A));
#else
    std::remove_const_t<T> result = 0;
    for (unsigned i = 0; i < A.rows(); ++i)
        result += A(i, i);
    return result;
#endif
}

/// @brief Returns the minimum value inside the vector, and its position.
/// @param v the input vector.
/// @return a pair containing the value and its position inside the vector.
template <typename T>
inline auto min(const Vector<T> &v)
{
    using data_type_t = std::remove_const_t<T>;
    if (v.empty())
        return std::make_pair(data_type_t(0.), unsigned(0));
    unsigned min_val_pos = 0;
    for (unsigned i = 1; i < v.size(); ++i)
        if (v[i] < v[min_val_pos])
            min_val_pos = i;
    return std::make_pair(v[min_val_pos], min_val_pos);
}

/// @brief Returns a vector with the minimum value for each column of the input matrix.
/// @param a the input matrix.
/// @return the minimum values of the columns of a.
template <typename T>
inline auto min(const Matrix<T> &a)
{
    using data_type_t = std::remove_const_t<T>;
    if (a.empty())
        return Vector<data_type_t>();
    Vector<data_type_t> values(a.cols(), data_type_t(0.));
    for (unsigned c = 0, r, min_val_pos; c < a.cols(); ++c) {
        for (min_val_pos = 0, r = 1; r < a.rows(); ++r)
            if (a(r, c) < a(min_val_pos, c))
                min_val_pos = r;
        values[c] = a(min_val_pos, c);
    }
    return values;
}

/// @brief Returns the maximum value inside the vector, and its position.
/// @param v the input vector.
/// @return a pair containing the value and its position inside the vector.
template <typename T>
inline auto max(const Vector<T> &v)
{
    using data_type_t = std::remove_const_t<T>;
    if (v.empty())
        return std::make_pair(data_type_t(0.), unsigned(0));
    unsigned max_val_pos = 0;
    for (unsigned i = 1; i < v.size(); ++i)
        if (v[i] > v[max_val_pos])
            max_val_pos = i;
    return std::make_pair(v[max_val_pos], max_val_pos);
}

/// @brief Returns a vector with the maximum value for each column of the input matrix.
/// @param a the input matrix.
/// @return the maximum values of the columns of a.
template <typename T>
inline auto max(const Matrix<T> &a)
{
    using data_type_t = std::remove_const_t<T>;
    if (a.empty())
        return Vector<data_type_t>();
    Vector<data_type_t> values(a.cols(), data_type_t(0.));
    for (unsigned c = 0, r, max_val_pos; c < a.cols(); ++c) {
        for (max_val_pos = 0, r = 1; r < a.rows(); ++r)
            if (a(r, c) > a(max_val_pos, c))
                max_val_pos = r;
        values[c] = a(max_val_pos, c);
    }
    return values;
}

/// @brief Checks if all elements are non-zero (or true).
/// @param a the input matrix.
/// @return true if all elements are non-zero (or true).
/// @return false if even one element is zero (or false).
template <typename T>
inline bool all(const MatrixBase<T> &a)
{
    for (unsigned i = 0; i < a.size(); ++i)
        if (a[i] == 0)
            return false;
    return true;
}

/// @brief Checks if all elements are non-zero (or true).
/// @param a the input vector.
/// @return true if all elements are non-zero (or true).
/// @return false if even one element is zero (or false).
template <typename T>
inline bool all(const Vector<T> &v)
{
    for (unsigned i = 0; i < v.size(); ++i)
        if (v[i] == 0)
            return false;
    return true;
}

/// @brief Checks if at least one element is non-zero (or true).
/// @param a the input matrix.
/// @return true if at least one element is non-zero (or true).
/// @return false if all elements are zero (or false).
template <typename T>
inline bool any(const MatrixBase<T> &a)
{
    for (unsigned i = 0; i < a.size(); ++i)
        if (a[i] != 0)
            return true;
    return false;
}

/// @brief Checks if at least one element is non-zero (or true).
/// @param a the input vector.
/// @return true if at least one element is non-zero (or true).
/// @return false if all elements are zero (or false).
template <typename T>
inline bool any(const Vector<T> &v)
{
    for (unsigned i = 0; i < v.size(); ++i)
        if (v[i] != 0)
            return true;
    return false;
}

} // namespace malg
