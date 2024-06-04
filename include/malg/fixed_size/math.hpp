/// @file math.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Basic mathematical functions.

#pragma once

#include "malg/fixed_size/vector.hpp"
#include "malg/fixed_size/matrix.hpp"
#include "malg/type_traits.hpp"

#include <type_traits>
#include <complex>
#include <utility>

namespace malg::details
{

template <typename T, typename Function, std::size_t... N>
[[nodiscard]] constexpr auto apply_unary_helper(
    const malg::VectorBase<T, sizeof...(N)> &l,
    Function func,
    std::integer_sequence<std::size_t, N...> s)
{
    using return_type_t = malg::Vector<decltype(func(l[0])), s.size()>;
    return return_type_t{ func(l[N])... };
}

template <class T, class U, class F, std::size_t... N>
[[nodiscard]] constexpr auto apply_binary_helper(
    const malg::VectorBase<T, sizeof...(N)> &l,
    const malg::VectorBase<U, sizeof...(N)> &r,
    F func,
    std::integer_sequence<std::size_t, N...> s)
{
    using return_type_t = malg::Vector<decltype(func(l[0], r[0])), s.size()>;
    return return_type_t{ func(l[N], r[N])... };
}

template <class T, class U, class F, std::size_t... N>
[[nodiscard]] constexpr auto apply_binary_helper(
    const U &l,
    const malg::VectorBase<T, sizeof...(N)> &r,
    F func,
    std::integer_sequence<std::size_t, N...> s)
{
    using return_type_t = malg::Vector<decltype(func(l, r[0])), s.size()>;
    return return_type_t{ func(l, r[N])... };
}

template <class T, class U, class F, std::size_t... N>
[[nodiscard]] constexpr auto apply_binary_helper(
    const malg::VectorBase<T, sizeof...(N)> &l,
    const U &r,
    F func,
    std::integer_sequence<std::size_t, N...> s)
{
    using return_type_t = malg::Vector<decltype(func(l[0], r)), s.size()>;
    return return_type_t{ func(l[N], r)... };
}

template <typename T, typename Function, std::size_t... N>
constexpr void apply_unary_helper(
    malg::VectorBase<T, sizeof...(N)> &l,
    Function func,
    std::integer_sequence<std::size_t, N...> s)
{
    ((l[N] = func(l[N])), ...);
}

template <class T, class U, class F, std::size_t... N>
constexpr void apply_binary_helper(
    malg::VectorBase<T, sizeof...(N)> &l,
    const malg::VectorBase<U, sizeof...(N)> &r,
    F func,
    std::integer_sequence<std::size_t, N...> s)
{
    ((l[N] = func(l[N], r[N])), ...);
}

template <class T, class U, class F, std::size_t... N>
constexpr void apply_binary_helper(
    malg::VectorBase<T, sizeof...(N)> &l,
    const U &r,
    F func,
    std::integer_sequence<std::size_t, N...> s)
{
    ((l[N] = func(l[N], r)), ...);
}

template <std::size_t N, class F, class T>
constexpr auto apply_unary(F func, T &&arg)
{
    return apply_unary_helper(arg, func, std::make_integer_sequence<std::size_t, N>{});
}

template <std::size_t N, class F, class T1, class T2>
constexpr auto apply_binary(F func, const T1 &arg1, const T2 &arg2)
{
    return apply_binary_helper(arg1, arg2, func, std::make_integer_sequence<std::size_t, N>{});
}

template <std::size_t N, class F, class T>
constexpr void apply_unary(F func, T &arg)
{
    apply_unary_helper(arg, func, std::make_integer_sequence<std::size_t, N>{});
}

template <std::size_t N, class F, class T1, class T2>
constexpr void apply_binary(F func, T1 &arg1, const T2 &arg2)
{
    apply_binary_helper(arg1, arg2, func, std::make_integer_sequence<std::size_t, N>{});
}

struct negate {
    template <class T>
    [[nodiscard]] constexpr auto operator()(const T &l) const noexcept
    {
        return -l;
    }
};

struct real {
    template <class T>
    [[nodiscard]] constexpr auto operator()(const std::complex<T> &l) const noexcept
    {
        return l.real();
    }
};

struct imag {
    template <class T>
    [[nodiscard]] constexpr auto operator()(const std::complex<T> &l) const noexcept
    {
        return l.imag();
    }
};

struct abs {
    template <class T>
    [[nodiscard]] constexpr auto operator()(const T &l) const noexcept
    {
        return std::abs(l);
    }
};

struct sign {
    template <class T>
    [[nodiscard]] constexpr auto operator()(const T &l) const noexcept
    {
        return l > 0 ? +1 : -1;
    }
};

struct plus {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l + r;
    }
};

struct minus {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l - r;
    }
};

struct multiplies {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l * r;
    }
};

struct divides {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l / r;
    }
};

struct equal {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l == r;
    }
};

struct not_equal_to {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l != r;
    }
};

struct less {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l < r;
    }
};

struct less_equal {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l <= r;
    }
};

struct greater {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l > r;
    }
};

struct greater_equal {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l >= r;
    }
};

} // namespace malg::details

#define DEFINE_UNARY_VECTOR_OPERATIONS(OP, FUNC)                     \
    template <class T, std::size_t N>                                \
    [[nodiscard]] constexpr auto OP(const malg::VectorBase<T, N> &l) \
    {                                                                \
        return malg::details::apply_unary<N>(FUNC{}, l);             \
    }

#define DEFINE_BINARY_VECTOR_OPERATIONS(OP, FUNC)                                                                     \
    template <class T, class U, std::size_t N>                                                                        \
    [[nodiscard]] constexpr auto OP(const malg::VectorBase<T, N> &l, const malg::VectorBase<U, N> &r)                 \
    {                                                                                                                 \
        return malg::details::apply_binary<N>(FUNC{}, l, r);                                                          \
    }                                                                                                                 \
    template <class T, class U, std::size_t N, typename = typename std::enable_if_t<std::is_arithmetic<U>::value, U>> \
    [[nodiscard]] constexpr auto OP(const U &l, const malg::VectorBase<T, N> &r)                                      \
    {                                                                                                                 \
        return malg::details::apply_binary<N>(FUNC{}, l, r);                                                          \
    }                                                                                                                 \
    template <class T, class U, std::size_t N, typename = typename std::enable_if_t<std::is_arithmetic<U>::value, U>> \
    [[nodiscard]] constexpr auto OP(const malg::VectorBase<T, N> &l, const U &r)                                      \
    {                                                                                                                 \
        return malg::details::apply_binary<N>(FUNC{}, l, r);                                                          \
    }

#define DEFINE_INPLACE_VECTOR_OPERATIONS(OP, FUNC)                                                                    \
    template <class T, class U, std::size_t N>                                                                        \
    [[nodiscard]] constexpr auto OP(malg::VectorBase<T, N> &l, const malg::VectorBase<U, N> &r)                       \
    {                                                                                                                 \
        return malg::details::apply_binary<N>(FUNC{}, l, r);                                                          \
    }                                                                                                                 \
    template <class T, class U, std::size_t N, typename = typename std::enable_if_t<std::is_arithmetic<U>::value, U>> \
    [[nodiscard]] constexpr auto OP(malg::VectorBase<T, N> &l, const U &r)                                            \
    {                                                                                                                 \
        return malg::details::apply_binary<N>(FUNC{}, l, r);                                                          \
    }

DEFINE_UNARY_VECTOR_OPERATIONS(operator-, malg::details::negate)
DEFINE_BINARY_VECTOR_OPERATIONS(operator+, malg::details::plus)
DEFINE_BINARY_VECTOR_OPERATIONS(operator-, malg::details::minus)
DEFINE_BINARY_VECTOR_OPERATIONS(operator*, malg::details::multiplies)
DEFINE_BINARY_VECTOR_OPERATIONS(operator/, malg::details::divides)

DEFINE_INPLACE_VECTOR_OPERATIONS(operator+=, malg::details::plus)
DEFINE_INPLACE_VECTOR_OPERATIONS(operator-=, malg::details::minus)
DEFINE_INPLACE_VECTOR_OPERATIONS(operator*=, malg::details::multiplies)
DEFINE_INPLACE_VECTOR_OPERATIONS(operator/=, malg::details::divides)

DEFINE_BINARY_VECTOR_OPERATIONS(operator==, malg::details::equal)
DEFINE_BINARY_VECTOR_OPERATIONS(operator!=, malg::details::not_equal_to)
DEFINE_BINARY_VECTOR_OPERATIONS(operator<, malg::details::less)
DEFINE_BINARY_VECTOR_OPERATIONS(operator<=, malg::details::less_equal)
DEFINE_BINARY_VECTOR_OPERATIONS(operator>, malg::details::greater)
DEFINE_BINARY_VECTOR_OPERATIONS(operator>=, malg::details::greater_equal)

namespace malg
{

/// @brief Computes the dot product between two vector.
/// @param a the first vector.
/// @param b the second vector.
/// @returns the scalar value resulting from the dot product.
template <class T, class U, std::size_t N>
inline auto dot(const malg::Vector<T, N> &a, const malg::Vector<U, N> &b)
{
    malg::Vector<decltype((a[0] * b[0])), N> result = 0;
    // Perform the operation.
    for (std::size_t i = 0; i < N; ++i) {
        result += a[i] * b[i];
    }
    return result;
}

/// @brief Computes the dot product between a matrix and a vector.
/// @param A the matrix.
/// @param b the vector.
/// @returns a vector with the same **rows** of *A*, resulting from the dot product.
template <class T, class U, std::size_t N1, std::size_t N2>
inline auto dot(const malg::Matrix<T, N1, N2> &A, const malg::Vector<U, N1> &b)
{
    malg::Vector<decltype((A[0][0] * b[0])), N1> result = 0;
    // Perform the computation.
    for (std::size_t r = 0; r < N1; ++r) {
        for (std::size_t c = 0; c < N2; ++c) {
            result[r] += A[r][c] * b[c];
        }
    }
    return result;
}

/// @brief Extracts the real part of the values of a.
/// @param a the matrix.
/// @returns the real part of the values of a.
DEFINE_UNARY_VECTOR_OPERATIONS(real, malg::details::real)

/// @brief Extracts the imaginary part of the values of a.
/// @param a the matrix.
/// @returns imaginary part of the values of a.
DEFINE_UNARY_VECTOR_OPERATIONS(imag, malg::details::imag)

/// @brief Transforms each element of a to its absolute value.
/// @param a the input matrix.
/// @returns the same matrix but its absolute values.
DEFINE_UNARY_VECTOR_OPERATIONS(abs, malg::details::abs)

/// @brief Returns a matrix containing the sign of the values in the input matrix.
/// @param a the input matrix.
/// @returns a matrix containing -1 and +1 based on the signs of the values of a.
DEFINE_UNARY_VECTOR_OPERATIONS(sign, malg::details::sign)

/// @brief Vector projection of **a** onto **b**.
/// @param a the first vector.
/// @param b the second vector.
/// @returns the vector projection.
template <typename T, typename U, std::size_t N>
inline auto projection(const malg::Vector<T, N> &a, const malg::Vector<U, N> &b)
{
    return b * (malg::dot(a, b) / malg::dot(b, b));
}

/// @brief Computes the vector length (or magnitude) by using the Pythagorean theorem, of a given column **c** of matrix **A**.
/// @param A the matrix.
/// @param c the column for which we compute the magnitude.
/// @param r_start the starting row.
/// @param r_end the ending row.
/// @returns the mangnitude of the column vector.
template <typename T, std::size_t N1, std::size_t N2>
inline auto vector_length(malg::Matrix<T, N1, N2> &A, std::size_t c, std::size_t r_start, std::size_t r_end)
{
    decltype((A[0][0] * A[0][0])) result = 0;
    for (std::size_t r = r_start; r < std::min(r_end, N1); r++) {
        if (malg::is_complex_v<T>) {
            result += std::norm(A[r][c]);
        } else {
            result += A[r][c] * A[r][c];
        }
    }
    return std::sqrt(result);
}

/// @brief Computes the linear combination of matrices.
/// @param A the first matrix.
/// @param a the first scalar.
/// @param B the second matrix.
/// @param b the second scalar.
/// @returns a matrix containing the linear combination.
template <typename T1, typename T2, typename T3, typename T4, std::size_t N1, std::size_t N2>
auto linear_combination(const malg::Matrix<T1, N1, N2> &A, const T2 &a, const malg::Matrix<T3, N1, N2> &B, const T4 &b)
{
    malg::Matrix<decltype((A[0][0] * a + B[0][0] * b)), N1, N2> result;
    for (std::size_t i = 0; i < A.size(); ++i) {
        result(i) = A(i) * a + B(i) * b;
    }
}

/// @brief Sums the element of the vector.
/// @param v the input vector.
/// @returns the sum of the elements.
template <typename T, std::size_t N>
inline auto sum(const malg::Vector<T, N> &v)
{
    std::remove_const_t<T> result = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        result += v[i];
    }
    return result;
}

/// @brief Compute the trace of A, i.e., the sum of the elements along the main diagonal.
/// @param A the input matrix.
/// @returns the sum of the diagonal elements.
template <typename T, std::size_t N>
inline auto trace(const malg::Matrix<T, N, N> &A)
{
    std::remove_const_t<T> result = 0;
    for (std::size_t i = 0; i < N; ++i) {
        result += A[i][i];
    }
}

/// @brief Returns the minimum value inside the vector, and its position.
/// @param v the input vector.
/// @returns a pair containing the value and its position inside the vector.
template <typename T, std::size_t N>
inline auto min(const malg::Vector<T, N> &v)
{
    using data_type_t = std::remove_const_t<T>;
    if (v.empty()) {
        return std::make_pair(data_type_t(0.), static_cast<std::size_t>(0));
    }
    std::size_t min_val_pos = 0;
    for (std::size_t i = 1; i < v.size(); ++i) {
        if (v[i] < v[min_val_pos]) {
            min_val_pos = i;
        }
    }
    return std::make_pair(v[min_val_pos], min_val_pos);
}

/// @brief Returns a vector with the minimum value for each column of the input matrix.
/// @param a the input matrix.
/// @returns the minimum values of the columns of a.
template <typename T, std::size_t N1, std::size_t N2>
inline auto min(malg::Matrix<T, N1, N2> &A)
{
    using data_type_t = std::remove_const_t<T>;
    malg::Vector<data_type_t, N2> result;
    if (!A.empty()) {
        for (std::size_t c = 0, r = 0, min_val_pos; c < N2; ++c) {
            for (min_val_pos = 0, r = 1; r < N1; ++r) {
                if (A[r][c] < A[min_val_pos][c]) {
                    min_val_pos = r;
                }
            }
            result[c] = A[min_val_pos][c];
        }
    }
    return result;
}

/// @brief Returns the maximum value inside the vector, and its position.
/// @param v the input vector.
/// @returns a pair containing the value and its position inside the vector.
template <typename T, std::size_t N>
inline auto max(const malg::Vector<T, N> &v)
{
    using data_type_t = std::remove_const_t<T>;
    if (v.empty()) {
        return std::make_pair(data_type_t(0.), static_cast<std::size_t>(0));
    }
    std::size_t max_val_pos = 0;
    for (std::size_t i = 1; i < v.size(); ++i) {
        if (v[i] > v[max_val_pos]) {
            max_val_pos = i;
        }
    }
    return std::make_pair(v[max_val_pos], max_val_pos);
}

/// @brief Returns a vector with the maximum value for each column of the input matrix.
/// @param a the input matrix.
/// @returns the maximum values of the columns of a.
template <typename T, std::size_t N1, std::size_t N2>
inline auto max(malg::Matrix<T, N1, N2> &A)
{
    using data_type_t = std::remove_const_t<T>;
    malg::Vector<data_type_t, N2> result;
    if (!A.empty()) {
        for (std::size_t c = 0, r = 0, max_val_pos; c < N2; ++c) {
            for (max_val_pos = 0, r = 1; r < N1; ++r) {
                if (A[r][c] > A[max_val_pos][c]) {
                    max_val_pos = r;
                }
            }
            result[c] = A[max_val_pos][c];
        }
    }
    return result;
}

/// @brief Checks if all elements are non-zero (or true).
/// @param a the input matrix.
/// @returns true if all elements are non-zero (or true).
/// @returns false if even one element is zero (or false).
template <typename T, std::size_t N1, std::size_t N2>
inline auto all(const malg::Matrix<T, N1, N2> &a)
{
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a(i) == 0) {
            return false;
        }
    }
    return true;
}

/// @brief Checks if all elements are non-zero (or true).
/// @param v the input vector.
/// @returns true if all elements are non-zero (or true).
/// @returns false if even one element is zero (or false).
template <typename T, std::size_t N>
inline auto all(const malg::Vector<T, N> &v)
{
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (v[i] == 0) {
            return false;
        }
    }
    return true;
}

/// @brief Checks if at least one element is non-zero (or true).
/// @param a the input matrix.
/// @returns true if at least one element is non-zero (or true).
/// @returns false if all elements are zero (or false).
template <typename T, std::size_t N1, std::size_t N2>
inline auto any(const malg::Matrix<T, N1, N2> &a)
{
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a(i) != 0) {
            return true;
        }
    }
    return false;
}

/// @brief Checks if at least one element is non-zero (or true).
/// @param v the input vector.
/// @returns true if at least one element is non-zero (or true).
/// @returns false if all elements are zero (or false).
template <typename T, std::size_t N>
inline auto any(const malg::Vector<T, N> &v)
{
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (v[i] != 0) {
            return true;
        }
    }
    return false;
}

/// @brief The Frobenius norm of a matrix.
/// @param A the input matrix.
/// @returns the norm.
template <typename T, std::size_t N1, std::size_t N2>
inline auto square_norm(const malg::Matrix<T, N1, N2> &A)
{
    std::remove_const_t<malg::is_complex_t<T>> accum = 0;
    // Compute the sum of squares of the elements of the given matrix.
    for (std::size_t r = 0; r < N1; ++r) {
        for (std::size_t c = 0; c < N2; ++c) {
            if (malg::is_complex_v<T>) {
                accum += std::norm(A[r][c]);
            } else {
                accum += A[r][c] * A[r][c];
            }
        }
    }
    // Return the square root of the sum of squares.
    return std::sqrt(accum);
}

/// @brief Computes the square norm of the vector.
/// @param v the vector.
/// @return the square norm of the vector.
template <typename T, std::size_t N>
inline auto square_norm(const malg::Vector<T, N> &v)
{
    std::remove_const_t<malg::is_complex_t<T>> accum = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        if (malg::is_complex_v<T>) {
            accum += std::norm(v[i]);
        } else {
            accum += v[i] * v[i];
        }
    }
    return std::sqrt(accum);
}

/// @brief Computes the infinity norm of a matrix, i.e., largest infinity norm among the rows of the matrix.
/// @param A the input matrix.
/// @returns the infinity norm.
template <typename T, std::size_t N1, std::size_t N2>
inline auto infinity_norm(const malg::Matrix<T, N1, N2> &A)
{
    using data_type_t = std::remove_const_t<malg::is_complex_t<T>>;
    if (A.empty()) {
        return data_type_t(0.);
    }
    data_type_t max{}, accum{};
    for (std::size_t r = 0; r < N1; ++r) {
        accum = 0.;
        for (std::size_t c = 0; c < N2; ++c) {
            accum += std::abs(A[r][c]);
        }
        max = std::max(max, accum);
    }
    return max;
}

/// @brief Computes the infinity norm of a vector, i.e., largest magnitude among each element of a vector.
/// @param v the input vector.
/// @returns the infinity norm.
template <typename T, std::size_t N>
inline auto infinity_norm(const malg::Vector<T, N> &v)
{
    using data_type_t = std::remove_const_t<T>;
    if (v.empty()) {
        return data_type_t(0.);
    }
    data_type_t max = std::abs(v[0]), tmp;
    for (std::size_t i = 1; i < v.size(); ++i) {
        tmp = std::abs(v[i]);
        if (max < tmp) {
            max = tmp;
        }
    }
    return max;
}

/// @brief The Euclidean norm of the lower leading diagonal of a square matrix.
/// @param A the input matrix.
/// @returns the square root of the sum of all the squares.
template <typename T, std::size_t N>
inline auto sub_norm(const malg::Matrix<T, N, N> &A)
{
    using data_type_t = std::remove_const_t<malg::is_complex_t<T>>;
    if (A.empty()) {
        return data_type_t(0.);
    }
    data_type_t accum = 0;
    for (std::size_t r = 1; r < N; ++r) {
        for (std::size_t c = 0; c < r; ++c) {
            if (malg::is_complex_v<T>) {
                accum += std::norm(A[r][c]);
            } else {
                accum += A[r][c] * A[r][c];
            }
        }
    }
    return std::sqrt(accum);
}

/// @brief Scale down matrix A by a power of 2, such that norm(A) < 1.
/// @param A the input matrix.
/// @returns the square root of the sum of all the squares.
template <typename T, std::size_t N1, std::size_t N2>
inline auto log2_ceil(const malg::Matrix<T, N1, N2> &A)
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
