/// @file utility.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Generic utility and interface functions.

#pragma once

#include "malg/matrix_base.hpp"
#include "malg/feq.hpp"
#include "malg/matrix.hpp"
#include "malg/type_traits.hpp"
#include "malg/vector.hpp"
#include "malg/view.hpp"

#include <algorithm>
#include <complex>
#include <cstdint>
#include <limits>
#include <random>

namespace malg::utility
{

/// @brief Find the basis for representing numbers on the computer.
/// @returns the basis.
template <typename T>
inline auto basis()
{
    T x = 1, eins = 1, b = 1;
    while ((x + eins) - x == eins) {
        x *= 2;
    }
    while ((x + b) == x) {
        b *= 2;
    }
    return static_cast<int>((x + b) - x);
}

/// @brief Creates a matrix with the given value on the diagonal
/// @param rows the number of rows.
/// @param cols the number of columns.
/// @param value value for the diagonal.
/// @returns the generated matrix.
template <typename T>
inline auto eye(std::size_t rows, std::size_t cols, T value)
{
    if (rows == 0) {
        throw std::invalid_argument("You must provide a number of rows greather than zero.");
    }
    if (cols == 0) {
        throw std::invalid_argument("You must provide a number of columns greather than zero.");
    }
    using data_type = std::remove_const_t<T>;
    Matrix<data_type> m(rows, cols, data_type(0.));
    for (std::size_t i = 0; i < std::min(rows, cols); ++i) {
        m(i, i) = value;
    }
    return m;
}

/// @brief Creates an identity matrix.
/// @param N size of the identity matrix.
/// @param value value for the diagonal.
/// @returns the identity matrix.
template <typename T>
inline auto identity(std::size_t N, T value)
{
    return utility::eye<T>(N, N, value);
}

/// @brief Creates a matrix with all zeros.
/// @param rows the number of rows.
/// @param cols the number of columns.
/// @returns the newly created matrix.
template <typename T>
inline auto zeros(std::size_t rows, std::size_t cols)
{
    if (rows == 0) {
        throw std::invalid_argument("You must provide a number of rows greather than zero.");
    }
    if (cols == 0) {
        throw std::invalid_argument("You must provide a number of columns greather than zero.");
    }
    using data_type = std::remove_const_t<T>;
    return Matrix<data_type>(rows, cols, 0);
}

/// @brief Creates a vector with all zeros.
/// @param size the size of the vector.
/// @returns the newly created vector.
template <typename T>
inline auto zeros(std::size_t size)
{
    if (size == 0) {
        throw std::invalid_argument("You must provide a size greather than zero.");
    }
    using data_type = std::remove_const_t<T>;
    return Vector<data_type>(size, 0);
}

/// @brief Creates a matrix with all ones.
/// @param rows the number of rows.
/// @param cols the number of columns.
/// @returns the newly created matrix.
template <typename T>
inline auto ones(std::size_t rows, std::size_t cols)
{
    if (rows == 0) {
        throw std::invalid_argument("You must provide a number of rows greather than zero.");
    }
    if (cols == 0) {
        throw std::invalid_argument("You must provide a number of columns greather than zero.");
    }
    using data_type = std::remove_const_t<T>;
    return Matrix<data_type>(rows, cols, 1);
}

/// @brief Creates a vector with all ones.
/// @param size the size of the vector.
/// @returns the newly created vector.
template <typename T>
inline auto ones(std::size_t size)
{
    if (size == 0) {
        throw std::invalid_argument("You must provide a size greather than zero.");
    }
    using data_type = std::remove_const_t<T>;
    return Vector<data_type>(size, 1);
}

/// @brief Extracts the diagonal elements from the matrix.
/// @param A the matrix.
/// @returns the diagonal elements.
template <typename T>
inline auto diag(const MatrixBase<T> &A)
{
    using data_type  = std::remove_const_t<T>;
    std::size_t size = std::min(A.rows(), A.cols());
    Vector<data_type> result(size, data_type(0.));
    for (std::size_t i = 0; i < size; ++i) {
        result[i] = A(i, i);
    }
    return result;
}

/// @brief Generates a diagonal matrix from the elements of the vector..
/// @param a the vector.
/// @returns the diagonal matrix.
template <typename T>
inline auto diag(const Vector<T> &a)
{
    using data_type = std::remove_const_t<T>;
    Matrix<data_type> matrix(a.size(), a.size(), data_type(0.));
    for (std::size_t i = 0; i < a.size(); ++i) {
        matrix(i, i) = a[i];
    }
    return matrix;
}

/// @brief Transforms a vector into a matrix.
/// @param a the input vector.
/// @param row_matrix if the output matrix should be a row matrix (true), or a column matrix (false).
/// @returns the generated matrix.
template <typename T>
inline auto to_matrix(const Vector<T> &a, bool row_matrix)
{
    using data_type = std::remove_const_t<T>;
    Matrix<data_type> matrix(row_matrix ? 1u : a.size(), row_matrix ? a.size() : 1u, data_type(0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        matrix[i] = a[i];
    }
    return matrix;
}

/// @brief Transforms a matrix into a vector.
/// @param A the input matrix.
/// @returns the vector.
template <typename T>
inline auto ravel(const MatrixBase<T> &A)
{
    using data_type = std::remove_const_t<T>;
    Vector<data_type> a(A.size());
    for (std::size_t i = 0; i < a.size(); ++i) {
        a[i] = A[i];
    }
    return a;
}

/// @brief Vertically stacks two matrices.
/// @param A the first (*xN) matrix.
/// @param B the second (*xN) matrix.
/// @returns the two matrices stacked into one.
template <typename T1, typename T2>
inline auto vstack(const MatrixBase<T1> &A, const MatrixBase<T2> &B)
{
    if (A.cols() != B.cols()) {
        throw std::invalid_argument("Vstack requires matrices with the same number of colmuns.");
    }
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the output matrix.
    Matrix<T> C = A;
    // Resize the matrix by adding rows for B.
    C.resize(A.rows() + B.rows(), A.cols());
    // Append B to C.
    for (std::size_t r = 0, c = 0; r < B.rows(); ++r) {
        for (c = 0; c < B.cols(); ++c) {
            C(r + A.rows(), c) = B(r, c);
        }
    }
    return C;
}

/// @brief Vertically stacks a matrix and a vector.
/// @param A the matrix of size (*xN).
/// @param b the vector of size N.
/// @returns the result of the stacking.
template <typename T1, typename T2>
inline auto vstack(const MatrixBase<T1> &A, const Vector<T2> &b)
{
    if (A.cols() != b.size()) {
        throw std::invalid_argument("The matrix must have the same number of columns as the size of the vector.");
    }
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the output matrix.
    Matrix<T> C = A;
    // Resize the matrix by adding the row for b.
    C.resize(A.rows() + 1, A.cols());
    // Append b to C.
    for (std::size_t i = 0; i < b.size(); ++i) {
        C(A.rows(), i) = b[i];
    }
    return C;
}

/// @brief Horizontally stacks two matrices.
/// @param A the first (Mx*) matrix.
/// @param B the second (Mx*) matrix.
/// @returns the two matrices stacked into one.
template <typename T1, typename T2>
inline auto hstack(const MatrixBase<T1> &A, const MatrixBase<T2> &B)
{
    if (A.rows() != B.rows()) {
        throw std::invalid_argument("Hstack requires matrices with the same number of rows.");
    }
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the output matrix.
    Matrix<T> C(A);
    // Resize the matrix by adding cols for B.
    C.resize(A.rows(), A.cols() + B.cols());
    // Append B to C.
    for (std::size_t c = 0, r = 0; c < B.cols(); ++c) {
        for (r = 0; r < B.rows(); ++r) {
            C(r, c + A.cols()) = B(r, c);
        }
    }
    return C;
}

/// @brief Horizontally stacks a matrix and a vector.
/// @param A the matrix of size (Mx*).
/// @param b the vector of size A.
/// @returns the result of the stacking.
template <typename T1, typename T2>
inline auto hstack(const MatrixBase<T1> &A, const Vector<T2> &b)
{
    if (A.rows() != b.size()) {
        throw std::invalid_argument("The matrix must have the same number of rows as the size of the vector.");
    }
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the output matrix.
    Matrix<T> C = A;
    // Resize the matrix by adding the row for b.
    C.resize(A.rows(), A.cols() + 1);
    // Append b to C.
    for (std::size_t i = 0; i < b.size(); ++i) {
        C(i, A.cols()) = b[i];
    }
    return C;
}

/// @brief Horizontally stacks a matrix and a vector.
/// @param A the matrix of size (Mx*).
/// @param b the vector of size A.
/// @returns the result of the stacking.
template <typename T1, typename T2>
inline auto concatenate(const Vector<T1> &a, const Vector<T2> &b)
{
    // Select the right type.
    using T = malg::extract_common_type_t<T1, T2>;
    // Create the output matrix.
    Vector<T> c(a);
    // Resize the matrix by adding the row for b.
    c.resize(a.size() + b.size());
    // Append b to C.
    for (std::size_t i = 0; i < b.size(); ++i) {
        c[a.size() + i] = b[i];
    }
    return c;
}

/// @brief Extracs a sub-matrix.
/// @param matrix the input matrix.
/// @param start_row the starting row.
/// @param end_row the ending row.
/// @param start_column the starting column.
/// @param end_column the ending column.
/// @returns the extracted sub-matrix.
template <typename T>
inline auto extract(const MatrixBase<T> &matrix,
                    std::size_t start_row    = 0,
                    std::size_t end_row      = std::numeric_limits<std::size_t>::max(),
                    std::size_t start_column = 0,
                    std::size_t end_column   = std::numeric_limits<std::size_t>::max())
{
    if (start_row >= end_row) {
        throw std::invalid_argument("The starting row must be lower than the ending row.");
    }
    if (start_column >= end_column) {
        throw std::invalid_argument("The starting column must be lower than the ending column.");
    }
    using data_type = std::remove_const_t<T>;
    end_row         = std::min(matrix.rows(), end_row);
    end_column      = std::min(matrix.cols(), end_column);
    Matrix<data_type> result(end_row - start_row, end_column - start_column, data_type(0.));
    for (std::size_t r = start_row; r < end_row; ++r) {
        for (std::size_t c = start_column; c < end_column; ++c) {
            result(r - start_row, c - start_column) = matrix(r, c);
        }
    }
    return result;
}

/// @brief Extracts a column from the matrix.
/// @param matrix the matrix.
/// @param column the column to extract.
/// @param start_row the starting row.
/// @param end_row the ending row.
/// @returns the extracted column as a vector.
template <typename T>
inline auto extract_column(const MatrixBase<T> &matrix,
                           std::size_t column,
                           std::size_t start_row = 0,
                           std::size_t end_row   = std::numeric_limits<std::size_t>::max())
{
    if (start_row >= end_row) {
        throw std::invalid_argument("The starting row must be lower than the ending row.");
    }
    if (column >= matrix.cols()) {
        throw std::invalid_argument("The selected column is outsize the matrix.");
    }
    using data_type = std::remove_const_t<T>;
    end_row         = std::min(matrix.rows(), end_row);
    Vector<data_type> result(end_row - start_row, data_type(0.));
    for (std::size_t r = start_row; r < end_row; ++r) {
        result[r - start_row] = matrix(r, column);
    }
    return result;
}

/// @brief Extracts a row from the matrix.
/// @param matrix the matrix.
/// @param row the row to extract.
/// @param start_column the starting column.
/// @param end_column the ending column.
/// @returns the extracted row as a vector.
template <typename T>
inline auto extract_row(const MatrixBase<T> &matrix,
                        std::size_t row,
                        std::size_t start_column = 0,
                        std::size_t end_column   = std::numeric_limits<std::size_t>::max())
{
    if (row >= matrix.rows()) {
        throw std::invalid_argument("The selected row is outsize the matrix.");
    }
    if (start_column >= end_column) {
        throw std::invalid_argument("The starting column must be lower than the ending column.");
    }
    using data_type = std::remove_const_t<T>;
    end_column      = std::min(matrix.cols(), end_column);
    Vector<data_type> result(end_column - start_column, data_type(0.));
    for (std::size_t c = start_column; c < end_column; ++c) {
        result[c - start_column] = matrix(row, c);
    }
    return result;
}

/// @brief Swaps the values of the two rows.
/// @param matrix the matrix.
/// @param i the first row.
/// @param j the second row.
/// @param start_column the starting column.
/// @param end_column the ending column.
template <typename T>
inline void swap_rows(MatrixBase<T> &matrix, std::size_t i, std::size_t j, std::size_t start_column = 0, std::size_t end_column = std::numeric_limits<std::size_t>::max())
{
    end_column = std::min(end_column, matrix.cols());
    for (std::size_t c = start_column; c < end_column; ++c) {
        std::swap(matrix(i, c), matrix(j, c));
    }
}

/// @brief Swaps the values of the two columns.
/// @param matrix the matrix.
/// @param i the first column.
/// @param j the second column.
/// @param start_row the starting row.
/// @param end_row the ending row.
template <typename T>
inline void swap_cols(MatrixBase<T> &matrix, std::size_t i, std::size_t j, std::size_t start_row = 0, std::size_t end_row = std::numeric_limits<std::size_t>::max())
{
    end_row = std::min(matrix.rows(), end_row);
    for (std::size_t r = start_row; r < end_row; ++r) {
        std::swap(matrix(r, i), matrix(r, j));
    }
}

/// @brief Removes the given row from the matrix.
/// @param matrix the input matrix.
/// @param row the row that must be removed.
/// @returns a reference to the input matrix with the row removed.
template <typename T>
inline auto &remove_row(Matrix<T> &matrix, std::size_t row)
{
    if (row >= matrix.rows()) {
        throw std::invalid_argument("The selected row is outsize the matrix.");
    }
    for (std::size_t r = row; r < (matrix.rows() - 1); ++r) {
        for (std::size_t c = 0; c < matrix.cols(); ++c) {
            std::swap(matrix(r, c), matrix(r + 1, c));
        }
    }
    matrix.resize(matrix.rows() - 1, matrix.cols());
    return matrix;
}

/// @brief Removes the given column from the matrix.
/// @param matrix the input matrix.
/// @param column the column that must be removed.
/// @returns a reference to the input matrix with the column removed.
template <typename T>
inline auto &remove_column(Matrix<T> &matrix, std::size_t column)
{
    if (column >= matrix.cols()) {
        throw std::invalid_argument("The selected column is outsize the matrix.");
    }
    for (std::size_t c = column; c < (matrix.cols() - 1); ++c) {
        for (std::size_t r = 0; r < matrix.rows(); ++r) {
            std::swap(matrix(r, c), matrix(r, c + 1));
        }
    }
    matrix.resize(matrix.rows(), matrix.cols() - 1);
    return matrix;
}

/// @brief Checks if the matrix is square.
/// @param matrix the matrix.
/// @returns true if the matrix is square.
/// @returns false if the matrix is not square.
template <typename T>
inline auto is_square(const MatrixBase<T> &matrix)
{
    return matrix.rows() == matrix.cols();
}

/// @brief Checks if the matrix is a column vector.
/// @param matrix the matrix.
/// @returns true if the matrix is a column vector.
/// @returns false if the matrix is not a column vector.
template <typename T>
inline auto is_column_vector(const MatrixBase<T> &matrix)
{
    return matrix.cols() == 1;
}

/// @brief Checks if the matrix is a row vector.
/// @param matrix the matrix.
/// @returns true if the matrix is a row vector.
/// @returns false if the matrix is not a row vector.
template <typename T>
inline auto is_row_vector(const MatrixBase<T> &matrix)
{
    return matrix.rows() == 1;
}

/// @brief Checks if the matrix is symmetric.
/// @param matrix the matrix.
/// @returns true if the matrix is symmetric.
/// @returns false if the matrix is not symmetric.
template <typename T>
inline auto is_symmetric(const MatrixBase<T> &matrix)
{
    if (!is_square(matrix)) {
        return false;
    }
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        for (std::size_t j = 0; j < matrix.cols(); ++j) {
            if (matrix(i, j) != matrix(j, i)) {
                return false;
            }
        }
    }
    return true;
}

/// @brief Checks if the matrix is skew symmetric.
/// @param matrix the matrix.
/// @returns true if the matrix is skew symmetric.
/// @returns false if the matrix is not skew symmetric.
template <typename T>
inline auto is_skew_symmetric(const MatrixBase<T> &matrix)
{
    if (!is_square(matrix)) {
        return false;
    }
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        for (std::size_t j = i + 1; j < matrix.cols(); ++j) {
            if (matrix(i, j) != -matrix(j, i)) {
                return false;
            }
        }
    }
    return true;
}

/// @brief Checks if the matrix is diagonal.
/// @param matrix the matrix.
/// @returns true if the matrix is diagonal.
/// @returns false if the matrix is not diagonal.
template <typename T>
inline auto is_diagonal(const MatrixBase<T> &matrix)
{
    if (!is_square(matrix)) {
        return false;
    }
    for (std::size_t i = 0; i < matrix.rows(); ++i) {
        for (std::size_t j = 0; j < matrix.cols(); ++j) {
            if ((i != j) && (matrix(i, j) != 0)) {
                return false;
            }
        }
    }
    return true;
}

/// @brief Checks if the matrix contains all zeros.
/// @param matrix the matrix.
/// @returns true if the matrix contains all zeros.
/// @returns false if the matrix does not contain all zeros.
template <typename T>
inline auto is_zero(const MatrixBase<T> &matrix)
{
    return std::all_of(matrix.begin(), matrix.end(), [](const T &x) { return x == 0; });
}

/// @brief Checks if the matrix contains the same value everywhere.
/// @param matrix the matrix.
/// @returns true if the matrix contains the same value everywhere.
/// @returns false if the matrix does not contain the same value everywhere.
template <typename T>
inline auto is_constant(const MatrixBase<T> &matrix)
{
    return std::adjacent_find(matrix.begin(), matrix.end(), std::not_equal_to{}) == matrix.end();
}

/// @brief Checks if the matrix contains the same value everywhere.
/// @param a the first matrix.
/// @param b the second matrix.
/// @param tolerance the tollerated difference between the two matrices.
/// @returns true if the two arrays are equal within the given tolerance.
/// @returns false otherwise.
template <typename T>
inline auto all_close(const Vector<T> &a, const Vector<T> &b, double tolerance = 0.0001)
{
    // Check the sizes.
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors has different sizes.");
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) > tolerance) {
            return false;
        }
    }
    return true;
}

/// @brief Checks if the matrix contains the same value everywhere.
/// @param A the first matrix.
/// @param B the second matrix.
/// @param tolerance the tollerated difference between the two matrices.
/// @returns true if the two arrays are equal within the given tolerance.
/// @returns false otherwise.
template <typename T>
inline auto all_close(const MatrixBase<T> &A, const MatrixBase<T> &B, double tolerance = 0.0001)
{
    // Check the sizes.
    if (A.rows() != B.rows()) {
        throw std::invalid_argument("Matrices has different number of rows.");
    }
    if (A.cols() != B.cols()) {
        throw std::invalid_argument("Matrices has different number of colmuns.");
    }
    for (std::size_t i = 0; i < A.size(); ++i) {
        if (std::abs(A[i] - B[i]) > tolerance) {
            return false;
        }
    }
    return true;
}

/// @brief Checks if the matrix is the identity matrix.
/// @param matrix the matrix.
/// @returns true if the matrix is the identity matrix.
/// @returns false if the matrix is not the identity matrix.
template <typename T>
inline auto is_identity(const MatrixBase<T> &matrix)
{
    if (!is_square(matrix)) {
        return false;
    }
    return (matrix == utility::identity<T>(matrix.cols()));
}

/// @brief Checks if the matrix is orthogonal.
/// @param matrix the matrix.
/// @returns true if the matrix is orthogonal.
/// @returns false if the matrix is not orthogonal.
template <typename T>
inline auto is_orthogonal(const MatrixBase<T> &matrix)
{
    if (!is_square(matrix)) {
        return false;
    }
    return (matrix * matrix.transpose() == utility::identity<T>(matrix.cols()));
}

/// @brief Checks if the matrix is invertible.
/// @param matrix the matrix.
/// @returns true if the matrix is invertible.
/// @returns false if the matrix is not invertible.
template <typename T>
inline auto is_invertible(const MatrixBase<T> &matrix)
{
    return matrix.determinant() != 0;
}

/// @brief Checks if the matrix is linearly dependent.
/// @param matrix the matrix.
/// @returns true if the matrix is linearly dependent.
/// @returns false if the matrix is not linearly dependent.
template <typename T>
inline auto is_linearly_dependent(const MatrixBase<T> &matrix)
{
    return matrix.determinant() == 0;
}

/// @brief Checks if the matrix is linearly independent.
/// @param matrix the matrix.
/// @returns true if the matrix is linearly independent.
/// @returns false if the matrix is not linearly independent.
template <typename T>
inline auto is_linearly_independent(const MatrixBase<T> &matrix)
{
    return !matrix.is_linearly_dependent();
}

/// @brief Checks if the matrix is a lower triangular one.
/// @param matrix the matrix.
/// @returns true if the matrix is a lower triangular one.
/// @returns false if the matrix is not a lower triangular one.
template <typename T>
inline auto is_lower_triangular(const MatrixBase<T> &matrix)
{
    if (!is_square(matrix)) {
        return false;
    }
    for (std::size_t r = 0, c = 0; r < matrix.rows(); ++r) {
        for (c = r + 1; c < matrix.cols(); ++c) {
            if (matrix(r, c)) {
                return false;
            }
        }
    }
    return true;
}

/// @brief Checks if the matrix is a upper triangular one.
/// @param matrix the matrix.
/// @returns true if the matrix is a upper triangular one.
/// @returns false if the matrix is not a upper triangular one.
template <typename T>
inline auto is_upper_triangular(const MatrixBase<T> &matrix)
{
    if (!is_square(matrix)) {
        return false;
    }
    for (std::size_t r = 0, c = 0; r < matrix.rows(); ++r) {
        for (c = 0; c < r; ++c) {
            if (matrix(r, c)) {
                return false;
            }
        }
    }
    return true;
}

/// @brief Changes the shape of the matrix **A**, but not the overall number of elements.
/// @param A the matrix.
/// @param rows the new number of rows.
/// @param cols the new number of columns.
/// @return the matrix **A** reshaped.
template <typename T>
inline auto &reshape(malg::Matrix<T> &A, std::size_t rows, std::size_t cols)
{
    return A.reshape(rows, cols);
}

/// @brief Changes the shape of the matrix **A**, without preserving the overall number of elements.
/// @param A the matrix.
/// @param rows the new number of rows.
/// @param cols the new number of columns.
/// @return the matrix **A** resized.
template <typename T>
inline auto &resize(malg::Matrix<T> &A, std::size_t rows, std::size_t cols)
{
    return A.resize(rows, cols);
}

/// @brief Computes the inner product of two lists by using iterators.
/// @param first1 the starting iterator of the **first list**.
/// @param last1 the final iterator of the **first list**.
/// @param first2 the starting iterator of the **second list**.
/// @param init the initial value.
/// @return the inner product of the two lists.
template <typename T>
inline T inner_product(T *first1, T *last1, T *first2, T init)
{
    while (first1 != last1) {
        init = std::move(init) + (*first1) * (*first2);
        ++first1, ++first2;
    }
    return init;
}

/// @brief Sums the lement of the given list.
/// @param first the starting iterator of the **list**.
/// @param last the final iterator of the **list**.
/// @param init the initial value.
/// @return the sums of the element of the list.
template <class T>
inline T accumulate(T *first, T *last, T init)
{
    while (first != last) {
        init = std::move(init) + (*first);
        ++first;
    }
    return init;
}

/// @brief Generates num points between min and max and return as vector.
/// @tparam T The type of the vector.
/// @param min The minimum value.
/// @param max The maximum value.
/// @param num The number of elements.
/// @return The generated vector.
template <typename T>
[[nodiscard]] inline auto linspace(T min, T max, std::size_t num = 100)
{
    malg::Vector<T> result(num, 0);
    if (num > 0) {
        if (num >= 2) {
            for (std::size_t i = 0; i < num - 1UL; ++i) {
                result[i] = min + (i * (max - min)) / std::floor(num - 1);
            }
        }
        result[num - 1] = max;
    }
    return result;
}

/// @brief Return a row vector with 'num' elements logarithmically spaced from 10^first to 10^last.
/// @tparam T The type of the vector.
/// @param first The first exponent.
/// @param last The last exponent.
/// @param num The number of elements.
/// @param base The base.
/// @return The generated vector.
template <typename T>
inline std::vector<T> logspace(T first, T last, std::size_t num = 50, T base = 10)
{
    T current_value = first, step = (last - first) / (num - 1);
    malg::Vector<T> result(num, 1);
    for (std::size_t i = 0L; i < num; ++i) {
        result[i] = std::pow(base, current_value);
        current_value += step;
    }
    return result;
}

/// @brief Generates a **view** for the given matrix **A**.
/// @param A the matrix.
/// @param start_row the starating row.
/// @param end_row the final row.
/// @param start_col the starating column.
/// @param end_col the final column.
/// @return the generated view.
template <typename MatrixType>
inline auto view(MatrixType &A, std::size_t start_row = 0, std::size_t end_row = std::numeric_limits<std::size_t>::max(), std::size_t start_col = 0, std::size_t end_col = std::numeric_limits<std::size_t>::max())
{
    return View(&A, start_row, end_row, start_col, end_col);
}

/// @brief Generates a **row view** for the given matrix **A**.
/// @param A the matrix.
/// @param row the row.
/// @param start_col the starating column.
/// @param end_col the final column.
/// @return the generated view.
template <typename MatrixType>
inline auto row(MatrixType &A, std::size_t row, std::size_t start_col = 0, std::size_t end_col = std::numeric_limits<std::size_t>::max())
{
    return View(&A, row, row + 1, start_col, end_col);
}

/// @brief Generates a **column view** for the given matrix **A**.
/// @param A the matrix.
/// @param col the column.
/// @param start_row the starating row.
/// @param end_row the final row.
/// @return the generated view.
template <typename MatrixType>
inline auto col(MatrixType &A, std::size_t col, std::size_t start_row = 0, std::size_t end_row = std::numeric_limits<std::size_t>::max())
{
    return View(&A, start_row, end_row, col, col + 1);
}

/// @brief Return indices that are non-zero in the flattened version of a.
/// @param A the input matrix.
/// @return an array containing the indices of the non-zero elements of A.
template <class T>
inline auto flatnonzero(const Matrix<T> &A)
{
    std::vector<std::size_t> nonzero_indices{};
    for (std::size_t i = 0; i < A.size(); ++i) {
        if (!malg::feq::approximately_equal(A[i], 0.)) {
            nonzero_indices.emplace_back(i);
        }
    }
    return nonzero_indices;
}

/// @brief Return a vector of indices of non-zero elements of matrix A,
/// 		as a row if A is a row vector, or as a column otherwise.
/// @tparam T The type of the matrix.
/// @param A		  The matrix to analyse.
/// @param n		  If provided, return only the first n indices.
/// @param first_last If true it returns the fist n, otherwise the last n.
/// @return A vector of indices.
template <class T>
inline auto find(const Matrix<T> &A, std::size_t n = std::numeric_limits<std::size_t>::max(), bool first_last = true)
{
    // Find all the indexes of nonzeros elements.
    std::vector<size_t> index = utility::flatnonzero(A);
    if (n == 0) {
        return Matrix<std::size_t>(0UL, 0UL);
    }
    if ((n > 0) && (n < index.size())) {
        std::size_t start = static_cast<std::size_t>(*index.begin());
        if (first_last) {
            index = std::vector<size_t>(start + 0UL, start + n);
        } else {
            index = std::vector<size_t>(start + index.size() - n, start + index.size());
        }
    }
    // Prepare the array for the indexes. If the input is a row-vector, return the indexes as a row-vector.
    if (utility::is_row_vector(A)) {
        return Matrix<std::size_t>(1UL, index.size(), index);
    }
    return Matrix<std::size_t>(index.size(), 1UL, index);
}

/// @brief Compares the two matrices.
/// @param A the first matrix.
/// @param B the second matrix.
/// @return
///      0 if they are the same,
///     -i if the i-th element of A is smaller than the i-th element of B,
///      i if the i-th element of A is greather than the i-th element of B.
template <typename T>
inline int64_t compare(const malg::MatrixBase<T> &A, const malg::MatrixBase<T> &B)
{
    // Check the sizes.
    if (A.size() != B.size()) {
        throw std::invalid_argument("Matrices has different size.");
    }
    for (std::size_t i = 0; i < A.size(); ++i) {
        if (A[i] < B[i]) {
            return -static_cast<int64_t>(i);
        }
        if (A[i] > B[i]) {
            return static_cast<int64_t>(i);
        }
    }
    return 0;
}

/// @brief Compares the two vectors.
/// @param a the first vector.
/// @param b the second vector.
/// @return
///      0 if they are the same,
///     -i if the i-th element of a is smaller than the i-th element of b,
///      i if the i-th element of a is greather than the i-th element of b.
template <typename T>
inline int64_t compare(const malg::Vector<T> &a, const malg::Vector<T> &b)
{
    // Check the sizes.
    if (a.size() != b.size()) {
        throw std::invalid_argument("Matrices has different size.");
    }
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i] < b[i]) {
            return -static_cast<int64_t>(i);
        }
        if (a[i] > b[i]) {
            return static_cast<int64_t>(i);
        }
    }
    return 0;
}

/// @brief Generates a random matrix.
/// @param rows the number of rows.
/// @param cols the number of columns.
/// @param min the minimum value.
/// @param max the maximum value.
/// @return the random matrix.
template <typename T1, typename T2>
inline auto rand_matrix(std::size_t rows, std::size_t cols, T2 min, T2 max)
{
    if constexpr (malg::is_complex_v<T1>) {
        // First, we extract the underlying datatype for the complex value.
        using data_type_t = typename T1::value_type;
        // Then, we define the type of for the distribution.
        using distribution_t = std::conditional_t<
            std::is_floating_point_v<data_type_t>,
            std::uniform_real_distribution<data_type_t>,
            std::uniform_int_distribution<data_type_t>>;
        // Create the distribution.
        distribution_t dist(min, max);
        // Create the interface to the random number generator.
        std::random_device dev;
        // Create the random number generator.
        std::mt19937 rng(dev());
        // Initialize the matrix.
        Matrix<std::complex<data_type_t>> m(rows, cols, std::complex<data_type_t>(0., 0.));
        // Generate the matrix.
        for (std::size_t i = 0; i < m.size(); ++i) {
            m[i] = std::complex<data_type_t>(dist(rng), dist(rng));
        }
        return m;
    } else {
        // We define the type of for the distribution.
        using distribution_t = std::conditional_t<
            std::is_floating_point_v<T1>,
            std::uniform_real_distribution<T1>,
            std::uniform_int_distribution<T1>>;
        // Create the distribution.
        distribution_t dist(min, max);
        // Create the interface to the random number generator.
        std::random_device dev;
        // Create the random number generator.
        std::mt19937 rng(dev());
        // Initialize the matrix.
        Matrix<T1> m(rows, cols, static_cast<T1>(0));
        // Generate the matrix.
        for (std::size_t i = 0; i < m.size(); ++i) {
            m[i] = dist(rng);
        }
        return m;
    }
}

/// @brief Generates a random vector.
/// @param size the size of the vector.
/// @param min the minimum value.
/// @param max the maximum value.
/// @return the random vector.
template <typename T1, typename T2>
inline auto rand_vector(std::size_t size, T2 min, T2 max)
{
    if constexpr (malg::is_complex_v<T1>) {
        // First, we extract the underlying datatype for the complex value.
        using data_type_t = typename T1::value_type;
        // Then, we define the type of for the distribution.
        using distribution_t = std::conditional_t<
            std::is_floating_point_v<data_type_t>,
            std::uniform_real_distribution<data_type_t>,
            std::uniform_int_distribution<data_type_t>>;
        // Create the distribution.
        distribution_t dist(min, max);
        // Create the interface to the random number generator.
        std::random_device dev;
        // Create the random number generator.
        std::mt19937 rng(dev());
        // Initialize the vector.
        Vector<T1> v(size, std::complex<data_type_t>(0., 0.));
        // Generate the vector.
        for (std::size_t i = 0; i < v.size(); ++i) {
            v[i] = std::complex<data_type_t>(dist(rng), dist(rng));
        }
        return v;
    } else {
        // We define the type of for the distribution.
        using distribution_t = std::conditional_t<
            std::is_floating_point_v<T1>,
            std::uniform_real_distribution<T1>,
            std::uniform_int_distribution<T1>>;
        // Create the distribution.
        distribution_t dist(min, max);
        // Create the interface to the random number generator.
        std::random_device dev;
        // Create the random number generator.
        std::mt19937 rng(dev());
        // Initialize the vector.
        Vector<T1> v(size, T1(0.));
        // Generate the vector.
        for (std::size_t i = 0; i < v.size(); ++i) {
            v[i] = dist(rng);
        }
        return v;
    }
}

} // namespace malg::utility
