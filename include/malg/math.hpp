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
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the result matrix.
    malg::Matrix<T> result(a.rows(), a.cols(), T(0.));
    // Compute the result.
    for (unsigned i = 0; i < result.size(); ++i)
        result[i] = a[i] + b[i];
    return result;
}

/// @brief Sums a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output matrix.
    malg::Matrix<T> result(a.rows(), a.cols(), T(0.));
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the result matrix.
    malg::Matrix<T> result(a.rows(), a.cols(), T(0.));
    // Compute the result.
    for (unsigned i = 0; i < result.size(); ++i)
        result[i] = a[i] - b[i];
    return result;
}

/// @brief Substraction between a matrix and a scalar.
/// @param a first matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output matrix.
    malg::Matrix<T> result(a.rows(), a.cols(), T(0.));
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
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output matrix.
    malg::Matrix<T> result(a.rows(), b.cols(), T(0.));
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator*(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output matrix.
    malg::Matrix<T> result(a.rows(), a.cols(), T(0.));
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
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output matrix.
    malg::Matrix<T> result(a.rows(), a.cols(), T(0.));
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator/(const malg::MatrixBase<T1> &a, const T2 &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output matrix.
    malg::Matrix<T> result(a.rows(), a.cols(), T(0.));
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] / b;
    return result;
}

/// @brief Division between a matrix and a scalar.
/// @param a first NxM matrix.
/// @param b scalar value.
/// @return the resulting matrix.
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
auto operator==(const malg::MatrixBase<T1> &a, const malg::Matrix<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = __approximately_equal(a[i], b[i]);
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator==(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = __approximately_equal(a[i], b);
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
auto operator!=(const malg::MatrixBase<T1> &a, const malg::Matrix<T2> &b)
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = !__approximately_equal(a[i], b[i]);
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
auto operator!=(const malg::MatrixBase<T1> &a, const T2 &b)
{
    malg::Matrix<bool> result(a.rows(), a.cols(), false);
    for (unsigned i = 0; i < a.size(); ++i) {
        if constexpr (std::is_floating_point_v<T1> || std::is_floating_point_v<T2>)
            result[i] = !__approximately_equal(a[i], b);
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output vector.
    malg::Vector<T> result(a.size(), T(0.));
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] + b[i];
    return result;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator+(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output vector.
    malg::Vector<T> result(a.size(), T(0.));
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

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output vector.
    malg::Vector<T> result(a.size(), T(0.));
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] - b[i];
    return result;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator-(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output vector.
    malg::Vector<T> result(a.size(), T(0.));
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

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output vector.
    malg::Vector<T> result(a.size(), T(0.));
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] * b[i];
    return result;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator*(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output vector.
    malg::Vector<T> result(a.size(), T(0.));
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] * b;
    return result;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T1>, T1>>
inline auto operator*(const T1 &a, const malg::Vector<T2> &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output vector.
    malg::Vector<T> result(b.size(), T(0.));
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a * b[i];
    return result;
}

template <typename T1, typename T2>
inline auto &operator*=(malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] * b[i];
    return a;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator*=(const malg::Vector<T1> &a, const T2 &b)
{
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] * b;
    return a;
}

template <typename T1, typename T2>
inline auto operator/(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output vector.
    malg::Vector<T> result(a.size(), T(0.));
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] / b[i];
    return result;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
inline auto operator/(const malg::Vector<T1> &a, const T2 &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Declare the output vector.
    malg::Vector<T> result(a.size(), T(0.));
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = a[i] / b;
    return result;
}

template <typename T1, typename T2>
inline auto &operator/=(malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    assert(a.size() == b.size());
    // Perform the operation.
    for (unsigned i = 0; i < a.size(); ++i)
        a[i] = a[i] / b[i];
    return a;
}

template <
    typename T1,
    typename T2,
    typename = typename std::enable_if_t<std::is_arithmetic_v<T2> || malg::is_complex_v<T2>, T2>>
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
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the result vector.
    malg::Vector<T> result(a.rows());
    // Perform the computation.
    for (unsigned r = 0; r < a.rows(); ++r)
        for (unsigned c = 0; c < a.cols(); ++c)
            result[r] += a(r, c) * b[c];
    return result;
}

template <typename T, typename TF>
inline auto element_wise_function(const malg::Vector<T> &a, TF fun)
{
    using data_type_t = std::remove_const_t<T>;
    // Create the resulting vector.
    malg::Vector<data_type_t> result(a.size());
    // Perform the computation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = fun(a[i]);
    return result;
}

template <typename T, typename TF>
inline auto element_wise_function(const malg::Matrix<T> &a, TF fun)
{
    using data_type_t = std::remove_const_t<T>;
    // Create the resulting vector.
    malg::Matrix<data_type_t> result(a.rows(), a.cols());
    // Perform the computation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = fun(a[i]);
    return result;
}

template <typename T1, typename T2, typename TF>
inline auto element_wise_binary_function(const malg::Vector<T1> &a, const malg::Vector<T2> &b, TF fun)
{
    assert(a.size() == b.size());
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the resulting vector.
    malg::Vector<T> result(a.size());
    // Perform the computation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = fun(a[i], b[i]);
    return result;
}

template <typename T1, typename T2, typename TF>
inline auto element_wise_binary_function(const malg::Matrix<T1> &a, const malg::Matrix<T2> &b, TF fun)
{
    assert(a.rows() == b.rows());
    assert(a.cols() == b.cols());
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the resulting matrix.
    malg::Matrix<T> result(a.rows(), a.cols());
    // Perform the computation.
    for (unsigned i = 0; i < a.size(); ++i)
        result[i] = fun(a[i], b[i]);
    return result;
}

template <typename T1, typename T2>
inline auto element_wise_product(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    return element_wise_binary_function(a, b, [](const T1 lhs, const T2 rhs) { return lhs * rhs; });
}

template <typename T1, typename T2>
inline auto element_wise_product(const malg::Matrix<T1> &a, const malg::Matrix<T2> &b)
{
    return element_wise_binary_function(a, b, [](const T1 lhs, const T2 rhs) { return lhs * rhs; });
}

template <typename T1, typename T2>
inline auto element_wise_div(const malg::Vector<T1> &a, const malg::Vector<T2> &b)
{
    return element_wise_binary_function(a, b, [](const T1 lhs, const T2 rhs) { return lhs / rhs; });
}

template <typename T1, typename T2>
inline auto element_wise_div(const malg::Matrix<T1> &a, const malg::Matrix<T2> &b)
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
    std::remove_const_t<extract_value_t<T>> accum = 0;
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
    std::remove_const_t<extract_value_t<T>> accum = 0;
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
    std::remove_const_t<extract_value_t<T>> accum = 0;
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
    std::remove_const_t<extract_value_t<T>> accum = 0;
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

} // namespace malg
