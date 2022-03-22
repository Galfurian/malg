/// @file utility.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Generic utility and interface functions.

#pragma once

#include "malg/matrix_base.hpp"
#include "malg/type_traits.hpp"
#include "malg/matrix.hpp"
#include "malg/vector.hpp"
#include "malg/view.hpp"

#include <algorithm>
#include <complex>
#include <cassert>
#include <random>

namespace malg::utility
{

/// @brief Find the basis for representing numbers on the computer.
/// @return the basis.
constexpr inline auto basis() noexcept
{
    double x, eins, b;
    x = eins = b = 1;
    while ((x + eins) - x == eins)
        x *= 2;
    while ((x + b) == x)
        b *= 2;
    return static_cast<int>((x + b) - x);
}

/// @brief Creates a matrix with the given value on the diagonal
/// @param rows the number of rows.
/// @param cols the number of columns.
/// @param value value for the diagonal.
/// @return the generated matrix.
template <typename T>
constexpr inline auto eye(unsigned rows, unsigned cols, T value = T(1.)) noexcept
{
    using data_type = std::remove_const_t<T>;
    assert((rows > 0) && "You must provide a number of rows greather than zero");
    assert((cols > 0) && "You must provide a number of columns greather than zero");
    Matrix<data_type> m(rows, cols, data_type(0.));
    for (unsigned i = 0; i < std::min(rows, cols); ++i)
        m(i, i) = value;
    return m;
}

/// @brief Creates an identity matrix.
/// @param N size of the identity matrix.
/// @param value value for the diagonal.
/// @return the identity matrix.
template <typename T>
constexpr inline auto identity(unsigned N, T value = T(1.)) noexcept
{
    return eye(N, N, value);
}

/// @brief Creates a matrix with all zeros.
/// @param rows the number of rows.
/// @param columns the number of columns.
/// @return the newly created matrix.
template <typename T>
constexpr inline auto zeros(unsigned rows, unsigned cols) noexcept
{
    using data_type = std::remove_const_t<T>;
    assert((rows > 0) && "You must provide a number of rows greather than zero");
    assert((cols > 0) && "You must provide a number of columns greather than zero");
    return Matrix<data_type>(rows, cols, 0);
}

/// @brief Creates a vector with all zeros.
/// @param size the size of the vector.
/// @return the newly created vector.
template <typename T>
constexpr inline auto zeros(unsigned size) noexcept
{
    using data_type = std::remove_const_t<T>;
    assert((size > 0) && "You must provide a size greather than zero");
    return Vector<data_type>(size, 0);
}

/// @brief Creates a matrix with all ones.
/// @param rows the number of rows.
/// @param columns the number of columns.
/// @return the newly created matrix.
template <typename T>
constexpr inline auto ones(unsigned rows, unsigned cols) noexcept
{
    using data_type = std::remove_const_t<T>;
    assert((rows > 0) && "You must provide a number of rows greather than zero");
    assert((cols > 0) && "You must provide a number of columns greather than zero");
    return Matrix<data_type>(rows, cols, 1);
}

/// @brief Creates a vector with all ones.
/// @param size the size of the vector.
/// @return the newly created vector.
template <typename T>
constexpr inline auto ones(unsigned size) noexcept
{
    using data_type = std::remove_const_t<T>;
    assert((size > 0) && "You must provide a size greather than zero");
    return Vector<data_type>(size, 1);
}

/// @brief Extracts the diagonal elements from the matrix.
/// @param matrix the matrix.
/// @return the diagonal elements.
template <typename T>
constexpr inline auto diag(const MatrixBase<T> &matrix) noexcept
{
    using data_type = std::remove_const_t<T>;
    unsigned size   = std::min(matrix.rows(), matrix.cols());
    Vector<data_type> result(size, data_type(0.));
    for (unsigned i = 0; i < size; ++i)
        result[i] = matrix(i, i);
    return result;
}

/// @brief Extracts the diagonal elements from the matrix.
/// @param matrix the matrix.
/// @return the diagonal elements.
template <typename T>
constexpr inline auto diag(const Vector<T> &v) noexcept
{
    using data_type = std::remove_const_t<T>;
    Matrix<data_type> matrix(v.size(), v.size(), data_type(0.));
    for (unsigned i = 0; i < v.size(); ++i)
        matrix(i, i) = v[i];
    return matrix;
}

/// @brief Extracts the diagonal elements from the matrix.
/// @param matrix the matrix.
/// @return the diagonal elements.
template <typename T>
constexpr inline auto to_matrix(const Vector<T> &v, bool row_matrix) noexcept
{
    using data_type = std::remove_const_t<T>;
    Matrix<data_type> matrix(row_matrix ? 1u : v.size(), row_matrix ? v.size() : 1u, data_type(0));
    for (unsigned i = 0; i < v.size(); ++i)
        matrix[i] = v[i];
    return matrix;
}

template <typename T1, typename T2>
constexpr inline auto vstack(const MatrixBase<T1> &A, const MatrixBase<T2> &B) noexcept
{
    assert(A.cols() == B.cols());
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the output matrix.
    Matrix<T> C = A;
    // Resize the matrix by adding rows for B.
    C.resize(A.rows() + B.rows(), A.cols());
    // Append B to C.
    for (unsigned r = 0, c; r < B.rows(); ++r)
        for (c = 0; c < B.cols(); ++c)
            C(r + A.rows(), c) = B(r, c);
    return C;
}

template <typename T1, typename T2>
constexpr inline auto vstack(const MatrixBase<T1> &A, const Vector<T2> &b) noexcept
{
    assert(A.cols() == b.size());
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the output matrix.
    Matrix<T> C = A;
    // Resize the matrix by adding the row for b.
    C.resize(A.rows() + 1, A.cols());
    // Append b to C.
    for (unsigned i = 0; i < b.size(); ++i)
        C(A.rows(), i) = b[i];
    return C;
}

template <typename T1, typename T2>
constexpr inline auto hstack(const MatrixBase<T1> &A, const MatrixBase<T2> &B) noexcept
{
    assert(A.rows() == B.rows());
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the output matrix.
    Matrix<T> C = A;
    // Resize the matrix by adding cols for B.
    C.resize(A.rows(), A.cols() + B.cols());
    // Append B to C.
    for (unsigned c = 0, r; c < B.cols(); ++c)
        for (r = 0; r < B.rows(); ++r)
            C(r, c + A.cols()) = B(r, c);
    return C;
}

template <typename T1, typename T2>
constexpr inline auto hstack(const MatrixBase<T1> &A, const Vector<T2> &b) noexcept
{
    assert(A.rows() == b.size());
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the output matrix.
    Matrix<T> C = A;
    // Resize the matrix by adding the row for b.
    C.resize(A.rows(), A.cols() + 1);
    // Append b to C.
    for (unsigned i = 0; i < b.size(); ++i)
        C(i, A.cols()) = b[i];
    return C;
}

template <typename T>
constexpr inline auto extract(const MatrixBase<T> &matrix,
                    unsigned start_row    = 0,
                    unsigned end_row      = -1,
                    unsigned start_column = 0,
                    unsigned end_column   = -1) noexcept
{
    using data_type = std::remove_const_t<T>;
    end_row         = std::min(matrix.rows(), end_row);
    end_column      = std::min(matrix.cols(), end_column);
    assert(start_row < end_row);
    assert(start_column < end_column);
    Matrix<data_type> result(end_row - start_row, end_column - start_column, data_type(0.));
    for (unsigned r = start_row; r < end_row; ++r)
        for (unsigned c = start_column; c < end_column; ++c)
            result(r - start_row, c - start_column) = matrix(r, c);
    return result;
}

/// @brief Extracts a column from the matrix.
/// @param matrix the matrix.
/// @param column the column to extract.
/// @return the extracted column as a vector.
template <typename T>
constexpr inline auto extract_column(const MatrixBase<T> &matrix,
                           unsigned column,
                           unsigned start_row = 0,
                           unsigned end_row   = -1) noexcept
{
    using data_type = std::remove_const_t<T>;
    end_row         = std::min(matrix.rows(), end_row);
    assert(start_row < end_row);
    assert(column < matrix.cols());
    Vector<data_type> result(end_row - start_row, data_type(0.));
    for (unsigned r = start_row; r < end_row; ++r)
        result[r - start_row] = matrix(r, column);
    return result;
}

/// @brief Extracts a row from the matrix.
/// @param matrix the matrix.
/// @param row the row to extract.
/// @return the extracted row as a vector.
template <typename T>
constexpr inline auto extract_row(const MatrixBase<T> &matrix,
                        unsigned row,
                        unsigned start_column = 0,
                        unsigned end_column   = -1) noexcept
{
    using data_type = std::remove_const_t<T>;
    end_column      = std::min(matrix.cols(), end_column);
    assert(row < matrix.rows());
    assert(start_column < end_column);
    Vector<data_type> result(end_column - start_column, data_type(0.));
    for (unsigned c = start_column; c < end_column; ++c)
        result[c - start_column] = matrix(row, c);
    return result;
}

/// @brief Swaps the values of the two rows.
/// @param matrix the matrix.
/// @param i the first row.
/// @param j the second row.
template <typename T>
constexpr inline void swap_rows(MatrixBase<T> &matrix, unsigned i, unsigned j, unsigned start_column = 0, unsigned end_column = -1) noexcept
{
    end_column = std::min(end_column, matrix.cols());
    for (unsigned c = start_column; c < end_column; ++c)
        std::swap(matrix(i, c), matrix(j, c));
}

/// @brief Swaps the values of the two columns.
/// @param matrix the matrix.
/// @param i the first column.
/// @param j the second column.
template <typename T>
constexpr inline void swap_cols(MatrixBase<T> &matrix, unsigned i, unsigned j, unsigned start_row = 0, unsigned end_row = -1) noexcept
{
    end_row = std::min(matrix.rows(), end_row);
    for (unsigned r = start_row; r < end_row; ++r)
        std::swap(matrix(r, i), matrix(r, j));
}

template <typename T>
constexpr inline auto &remove_row(Matrix<T> &matrix, unsigned row) noexcept
{
    assert(row < matrix.rows());
    for (unsigned r = row; r < (matrix.rows() - 1); ++r) {
        for (unsigned c = 0; c < matrix.cols(); ++c) {
            std::swap(matrix(r, c), matrix(r + 1, c));
        }
    }
    matrix.resize(matrix.rows() - 1, matrix.cols());
    return matrix;
}

template <typename T>
constexpr inline auto &remove_column(Matrix<T> &matrix, unsigned column) noexcept
{
    assert(column < matrix.cols());
    for (unsigned c = column; c < (matrix.cols() - 1); ++c) {
        for (unsigned r = 0; r < matrix.rows(); ++r) {
            std::swap(matrix(r, c), matrix(r, c + 1));
        }
    }
    matrix.resize(matrix.rows(), matrix.cols() - 1);
    return matrix;
}

/// @brief Checks if the matrix is square.
/// @param matrix the matrix.
/// @return true if the matrix is square.
/// @return false if the matrix is not square.
template <typename T>
constexpr inline auto is_square(const MatrixBase<T> &matrix) noexcept
{
    return matrix.rows() == matrix.cols();
}

/// @brief Checks if the matrix is a column vector.
/// @param matrix the matrix.
/// @return true if the matrix is a column vector.
/// @return false if the matrix is not a column vector.
template <typename T>
constexpr inline auto is_column_vector(const MatrixBase<T> &matrix) noexcept
{
    return matrix.cols() == 1;
}

/// @brief Checks if the matrix is a row vector.
/// @param matrix the matrix.
/// @return true if the matrix is a row vector.
/// @return false if the matrix is not a row vector.
template <typename T>
constexpr inline auto is_row_vector(const MatrixBase<T> &matrix) noexcept
{
    return matrix.rows() == 1;
}

/// @brief Checks if the matrix is symmetric.
/// @param matrix the matrix.
/// @return true if the matrix is symmetric.
/// @return false if the matrix is not symmetric.
template <typename T>
constexpr inline auto is_symmetric(const MatrixBase<T> &matrix) noexcept
{
    if (!is_square(matrix))
        return false;
    for (unsigned i = 0; i < matrix.rows(); ++i)
        for (unsigned j = 0; j < matrix.cols(); ++j)
            if (matrix(i, j) != matrix(j, i))
                return false;
    return true;
}

/// @brief Checks if the matrix is skew symmetric.
/// @param matrix the matrix.
/// @return true if the matrix is skew symmetric.
/// @return false if the matrix is not skew symmetric.
template <typename T>
constexpr inline auto is_skew_symmetric(const MatrixBase<T> &matrix) noexcept
{
    if (!is_square(matrix))
        return false;
    for (unsigned i = 0; i < matrix.rows(); ++i)
        for (unsigned j = i + 1; j < matrix.cols(); ++j)
            if (matrix(i, j) != -matrix(j, i))
                return false;
    return true;
}

/// @brief Checks if the matrix is diagonal.
/// @param matrix the matrix.
/// @return true if the matrix is diagonal.
/// @return false if the matrix is not diagonal.
template <typename T>
constexpr inline auto is_diagonal(const MatrixBase<T> &matrix) noexcept
{
    if (!is_square(matrix))
        return false;
    for (unsigned i = 0; i < matrix.rows(); ++i)
        for (unsigned j = 0; j < matrix.cols(); ++j)
            if ((i != j) && (matrix(i, j) != 0))
                return false;
    return true;
}

/// @brief Checks if the matrix contains all zeros.
/// @param matrix the matrix.
/// @return true if the matrix contains all zeros.
/// @return false if the matrix does not contain all zeros.
template <typename T>
constexpr inline auto is_zero(const MatrixBase<T> &matrix) noexcept
{
    return std::all_of(matrix.begin(), matrix.end(), [](const T &x) { return x == 0; });
}

/// @brief Checks if the matrix contains the same value everywhere.
/// @param matrix the matrix.
/// @return true if the matrix contains the same value everywhere.
/// @return false if the matrix does not contain the same value everywhere.
template <typename T>
constexpr inline auto is_constant(const MatrixBase<T> &matrix) noexcept
{
    return std::adjacent_find(matrix.begin(), matrix.end(), std::not_equal_to{}) == matrix.end();
}

/// @brief Checks if the matrix contains the same value everywhere.
/// @param a the first matrix.
/// @param b the second matrix.
/// @param tolerance the tollerated difference between the two matrices.
/// @return true if the two arrays are equal within the given tolerance.
/// @return false otherwise.
template <typename T>
constexpr inline auto all_close(const Vector<T> &a, const Vector<T> &b, double tolerance = 0.0001) noexcept
{
    assert(a.size() == b.size());
    for (unsigned i = 0; i < a.size(); ++i)
        if (std::abs(a[i] - b[i]) > tolerance)
            return false;
    return true;
}

/// @brief Checks if the matrix contains the same value everywhere.
/// @param a the first matrix.
/// @param b the second matrix.
/// @param tolerance the tollerated difference between the two matrices.
/// @return true if the two arrays are equal within the given tolerance.
/// @return false otherwise.
template <typename T>
constexpr inline auto all_close(const MatrixBase<T> &a, const MatrixBase<T> &b, double tolerance = 0.0001) noexcept
{
    assert((a.rows() == b.rows()) && (a.cols() == b.cols()));
    for (unsigned i = 0; i < a.size(); ++i)
        if (std::abs(a[i] - b[i]) > tolerance)
            return false;
    return true;
}

/// @brief Checks if the matrix is the identity matrix.
/// @param matrix the matrix.
/// @return true if the matrix is the identity matrix.
/// @return false if the matrix is not the identity matrix.
template <typename T>
constexpr inline auto is_identity(const MatrixBase<T> &matrix) noexcept
{
    if (!is_square(matrix))
        return false;
    return (matrix == utility::identity<T>(matrix.cols()));
}

/// @brief Checks if the matrix is orthogonal.
/// @param matrix the matrix.
/// @return true if the matrix is orthogonal.
/// @return false if the matrix is not orthogonal.
template <typename T>
constexpr inline auto is_orthogonal(const MatrixBase<T> &matrix) noexcept
{
    if (!is_square(matrix))
        return false;
    return (matrix * matrix.transpose() == utility::identity<T>(matrix.cols()));
}

/// @brief Checks if the matrix is invertible.
/// @param matrix the matrix.
/// @return true if the matrix is invertible.
/// @return false if the matrix is not invertible.
template <typename T>
constexpr inline auto is_invertible(const MatrixBase<T> &matrix) noexcept
{
    return matrix.determinant() != 0;
}

/// @brief Checks if the matrix is linearly dependent.
/// @param matrix the matrix.
/// @return true if the matrix is linearly dependent.
/// @return false if the matrix is not linearly dependent.
template <typename T>
constexpr inline auto is_linearly_dependent(const MatrixBase<T> &matrix) noexcept
{
    return matrix.determinant() == 0;
}

/// @brief Checks if the matrix is linearly independent.
/// @param matrix the matrix.
/// @return true if the matrix is linearly independent.
/// @return false if the matrix is not linearly independent.
template <typename T>
constexpr inline auto is_linearly_independent(const MatrixBase<T> &matrix) noexcept
{
    return !matrix.is_linearly_dependent();
}

/// @brief Checks if the matrix is a lower triangular one.
/// @param matrix the matrix.
/// @return true if the matrix is a lower triangular one.
/// @return false if the matrix is not a lower triangular one.
template <typename T>
constexpr inline auto is_lower_triangular(const MatrixBase<T> &matrix) noexcept
{
    if (!is_square(matrix))
        return false;
    for (unsigned r = 0, c; r < matrix.rows(); ++r)
        for (c = r + 1; c < matrix.cols(); ++c)
            if (matrix(r, c))
                return false;
    return true;
}

/// @brief Checks if the matrix is a upper triangular one.
/// @param matrix the matrix.
/// @return true if the matrix is a upper triangular one.
/// @return false if the matrix is not a upper triangular one.
template <typename T>
constexpr inline auto is_upper_triangular(const MatrixBase<T> &matrix) noexcept
{
    if (!is_square(matrix))
        return false;
    for (unsigned r = 0, c; r < matrix.rows(); ++r)
        for (c = 0; c < r; ++c)
            if (matrix(r, c))
                return false;
    return true;
}

template <typename T>
constexpr inline auto &reshape(malg::Matrix<T> &m, unsigned rows, unsigned cols) noexcept
{
    return m.reshape(rows, cols);
}

template <typename T>
constexpr inline auto &resize(malg::Matrix<T> &m, unsigned rows, unsigned cols) noexcept
{
    return m.resize(rows, cols);
}

template <typename T>
constexpr inline T inner_product(T *first1, T *last1, T *first2, T init) noexcept
{
    while (first1 != last1) {
        init = std::move(init) + (*first1) * (*first2);
        ++first1, ++first2;
    }
    return init;
}

template <typename MatrixType>
constexpr inline auto view(MatrixType &matrix, unsigned start_row = 0, unsigned end_row = -1, unsigned start_col = 0, unsigned end_col = -1) noexcept
{
    return View(matrix, start_row, end_row, start_col, end_col);
}

template <typename MatrixType>
constexpr inline auto row(MatrixType &matrix, unsigned row, unsigned start_col = 0, unsigned end_col = -1) noexcept
{
    return View(matrix, row, row + 1, start_col, end_col);
}

template <typename MatrixType>
constexpr inline auto col(MatrixType &matrix, unsigned col, unsigned start_row = 0, unsigned end_row = -1) noexcept
{
    return View(matrix, start_row, end_row, col, col + 1);
}

template <class T>
constexpr inline T accumulate(T *first, T *last, T init) noexcept
{
    while (first != last) {
        init = std::move(init) + (*first);
        ++first;
    }
    return init;
}

/// @brief Find the basis for representing numbers on the computer.
/// @tparam T
/// @return int
template <typename T>
static constexpr inline auto basis() noexcept
{
    T x = 1, eins = 1, b = 1;
    while ((x + eins) - x == eins)
        x *= 2;
    while ((x + b) == x)
        b *= 2;
    return (int)((x + b) - x);
}

template <typename T, typename T2>
constexpr inline auto rand_matrix(unsigned rows, unsigned cols, T2 min, T2 max) noexcept
{
    if constexpr (malg::is_complex_v<T>) {
        using data_type_t    = typename T::value_type;
        using distribution_t = std::conditional_t<
            std::is_floating_point_v<data_type_t>,
            std::uniform_real_distribution<data_type_t>,
            std::uniform_int_distribution<data_type_t>>;
        distribution_t dist(min, max);
        std::random_device dev;
        std::mt19937 rng(dev());
        Matrix<std::complex<data_type_t>> m(rows, cols, std::complex<data_type_t>(0., 0.));
        for (unsigned i = 0; i < m.size(); ++i)
            m[i] = std::complex<data_type_t>(dist(rng), dist(rng));
        return m;
    } else {
        using distribution_t = std::conditional_t<
            std::is_floating_point_v<T>,
            std::uniform_real_distribution<T>,
            std::uniform_int_distribution<T>>;
        distribution_t dist(min, max);
        std::random_device dev;
        std::mt19937 rng(dev());
        Matrix<T> m(rows, cols, static_cast<T>(0));
        for (unsigned i = 0; i < m.size(); ++i)
            m[i] = dist(rng);
        return m;
    }
}

template <typename T, typename T2>
constexpr inline auto rand_vector(unsigned size, T2 min, T2 max) noexcept
{
    if constexpr (malg::is_complex_v<T>) {
        using data_type_t    = typename T::value_type;
        using distribution_t = std::conditional_t<
            std::is_floating_point_v<data_type_t>,
            std::uniform_real_distribution<data_type_t>,
            std::uniform_int_distribution<data_type_t>>;
        distribution_t dist(min, max);
        std::random_device dev;
        std::mt19937 rng(dev());
        Vector<T> v(size, std::complex<data_type_t>(0., 0.));
        for (unsigned i = 0; i < v.size(); ++i)
            v[i] = std::complex<data_type_t>(dist(rng), dist(rng));
        return v;
    } else {
        using distribution_t = std::conditional_t<
            std::is_floating_point_v<T>,
            std::uniform_real_distribution<T>,
            std::uniform_int_distribution<T>>;
        distribution_t dist(min, max);
        std::random_device dev;
        std::mt19937 rng(dev());
        Vector<T> v(size, T(0.));
        for (unsigned i = 0; i < v.size(); ++i)
            v[i] = dist(rng);
        return v;
    }
}

} // namespace malg::utility