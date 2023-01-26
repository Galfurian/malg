/// @file linalg.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Linear algebra algorithms.

#pragma once

#include "malg/utility.hpp"
#include "malg/vector.hpp"
#include "malg/matrix.hpp"
#include "malg/math.hpp"

#include <algorithm>
#include <complex>

namespace malg::linalg
{

/// @brief Computes the transpose of the matrix.
/// @param A the matrix.
/// @returns the matrix transposed.
template <typename T>
inline auto transpose(const MatrixBase<T> &A)
{
    Matrix<std::remove_const_t<T>> result(A.cols(), A.rows(), 0.0);
    for (std::size_t r = 0, c = 0; r < A.rows(); ++r)
        for (c = 0; c < A.cols(); ++c)
            result(c, r) = A(r, c);
    return result;
}

/// @brief Computes the factorial of the input value.
/// @param n the input value.
/// @returns the factorial.
template <typename T>
inline constexpr T factorial(T n)
{
    return (n <= 0) ? 1 : n * factorial<T>(n - 1);
}

/// @brief Computes the power of the input matrix.
/// @param A the input matrix.
/// @param N the power to compute.
/// @returns the result of the operation.
template <typename T>
inline auto powm(const MatrixBase<T> &A, std::size_t N)
{
    auto R = utility::eye<std::remove_const_t<T>>(A.rows(), A.cols());
    for (std::size_t k = 0; k < N; ++k)
        R *= A;
    return R;
}

/// @brief Computes the hermitian transpose of the complex matrix.
/// @param A the matrix.
/// @returns the matrix transposed.
template <typename T>
inline auto hermitian_transpose(const MatrixBase<std::complex<T>> &A)
{
    Matrix<std::remove_const_t<T>> result(A.cols(), A.rows(), 0.0);
    for (std::size_t r = 0, c = 0; r < A.rows(); ++r)
        for (c = 0; c < A.cols(); ++c)
            result(c, r) = std::conj(A(r, c));
    return result;
}

/// @brief Computes the transpose and takes the complex conjugate of each entry.
/// @param matrix the matrix.
/// @returns the matrix transposed.
template <typename T>
inline auto conjugate_transpose(const MatrixBase<std::complex<T>> &matrix)
{
    Matrix<std::complex<std::remove_const_t<T>>> result(matrix.cols(), matrix.rows(), 0.);
    for (std::size_t r = 0, c = 0; r < matrix.rows(); ++r)
        for (c = 0; c < matrix.cols(); ++c)
            result(c, r) = std::conj(matrix(r, c));
    return result;
}

/// @brief Function to get cofactor of the matrix.
/// @param matrix the input matrix.
/// @param p the row that must be removed.
/// @param q the column that must be removed.
/// @returns An [N-1, N-1] matrix, generated by removing row p and column q from
/// the input matrix.
template <typename T>
inline auto cofactor(const MatrixBase<T> &matrix, std::size_t p, std::size_t q)
{
    // Create the output matrix.
    Matrix<std::remove_const_t<T>> output(matrix.rows() - 1, matrix.cols() - 1, 0);
    // Create the indexing variables.
    std::size_t i, j, row, col;
    // Looping for each element of the matrix.
    for (i = 0, j = 0, row = 0, col = 0; row < matrix.rows(); ++row) {
        for (col = 0; col < matrix.cols(); ++col) {
            // Copying only those element which are not in given row and column.
            if ((row != p) && (col != q)) {
                output(i, j++) = matrix(row, col);
                // When the row is filled, increase row index and reset col index.
                if (j == (matrix.cols() - 1))
                    j = 0, ++i;
            }
        }
    }
    return output;
}

/// @brief Use Gaussian Elimination to create an Upper Diagonalized matrix then
/// multiplying the diagonal to get the determinant.
/// @param matrix the input matrix.
/// @returns the determinant of the matrix.
template <typename CT>
inline auto determinant(const MatrixBase<CT> &matrix)
{
    if (!malg::utility::is_square(matrix))
        throw std::invalid_argument("The input matrix is not square.");
    // Create an alias for the non-constant type.
    using T = std::remove_const_t<CT>;
    // Get the size of the matrix.
    const std::size_t N = matrix.cols();
    // Fast exit with a 1x1.
    if (N == 1) {
        return matrix(0, 0);
    }
    // Fast exit with a 2x2.
    if (N == 2) {
        return matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0);
    }
    // Fast exit with a 3x3.
    if (N == 3) {
        return matrix(0, 0) * (matrix(1, 1) * matrix(2, 2) - matrix(2, 1) * matrix(1, 2)) -
               matrix(0, 1) * (matrix(1, 0) * matrix(2, 2) - matrix(1, 2) * matrix(2, 0)) +
               matrix(0, 2) * (matrix(1, 0) * matrix(2, 1) - matrix(1, 1) * matrix(2, 0));
    }
    // Fast exit with a 4x4.
    if (N == 4) {
        return matrix(0, 3) * matrix(1, 2) * matrix(2, 1) * matrix(3, 0) - matrix(0, 2) * matrix(1, 3) * matrix(2, 1) * matrix(3, 0) -
               matrix(0, 3) * matrix(1, 1) * matrix(2, 2) * matrix(3, 0) + matrix(0, 1) * matrix(1, 3) * matrix(2, 2) * matrix(3, 0) +
               matrix(0, 2) * matrix(1, 1) * matrix(2, 3) * matrix(3, 0) - matrix(0, 1) * matrix(1, 2) * matrix(2, 3) * matrix(3, 0) -
               matrix(0, 3) * matrix(1, 2) * matrix(2, 0) * matrix(3, 1) + matrix(0, 2) * matrix(1, 3) * matrix(2, 0) * matrix(3, 1) +
               matrix(0, 3) * matrix(1, 0) * matrix(2, 2) * matrix(3, 1) - matrix(0, 0) * matrix(1, 3) * matrix(2, 2) * matrix(3, 1) -
               matrix(0, 2) * matrix(1, 0) * matrix(2, 3) * matrix(3, 1) + matrix(0, 0) * matrix(1, 2) * matrix(2, 3) * matrix(3, 1) +
               matrix(0, 3) * matrix(1, 1) * matrix(2, 0) * matrix(3, 2) - matrix(0, 1) * matrix(1, 3) * matrix(2, 0) * matrix(3, 2) -
               matrix(0, 3) * matrix(1, 0) * matrix(2, 1) * matrix(3, 2) + matrix(0, 0) * matrix(1, 3) * matrix(2, 1) * matrix(3, 2) +
               matrix(0, 1) * matrix(1, 0) * matrix(2, 3) * matrix(3, 2) - matrix(0, 0) * matrix(1, 1) * matrix(2, 3) * matrix(3, 2) -
               matrix(0, 2) * matrix(1, 1) * matrix(2, 0) * matrix(3, 3) + matrix(0, 1) * matrix(1, 2) * matrix(2, 0) * matrix(3, 3) +
               matrix(0, 2) * matrix(1, 0) * matrix(2, 1) * matrix(3, 3) - matrix(0, 0) * matrix(1, 2) * matrix(2, 1) * matrix(3, 3) -
               matrix(0, 1) * matrix(1, 0) * matrix(2, 2) * matrix(3, 3) + matrix(0, 0) * matrix(1, 1) * matrix(2, 2) * matrix(3, 3);
    }
    // We need to create a temporary.
    Matrix<T> A(matrix);
    // Create the indexing variables.
    std::size_t c, r = 0, k;
    // Initialize the determinant, and create both pivot and ratio variable.
    T det = static_cast<T>(1.), pivot, ratio;
    // We convert the temporary to upper triangular form.
    for (c = 0; c < N; ++c) {
        // If we have a negative value on the diagonal, we need to move it
        // somewhere else.
        if (A(c, c) == 0.) {
            // Right now, I'm trying to find a place below the current
            k = c + 1;
            while ((k < A.rows()) && (A(k, c) == 0.))
                k++;
            // If we did not find a non-zero value, we have a singular matrix.
            if (k == A.rows())
                break;
            // Swap the rows.
            utility::swap_rows(A, c, k);
            // Every time we swap rows, we need to change the sign to the
            // determinant.
            det *= -1;
        }
        // Store the pivot.
        pivot = A(c, c);
        for (r = c + 1; r < N; ++r) {
            ratio = A(r, c) / pivot;
            for (k = c; k < N; ++k) {
                A(r, k) -= ratio * A(c, k);
            }
        }
        // Multiply the pivot for the determinant variable.
        det *= pivot;
    }
    return det;
}

/// @brief Computes the adjoint of this matrix.
/// @param matrix the input matrix.
/// @returns the adjoint of the matrix.
template <typename T>
inline auto adjoint(const MatrixBase<T> &matrix)
{
    if (!malg::utility::is_square(matrix))
        throw std::invalid_argument("The input matrix is not square.");
    // Get the size of the matrix.
    const std::size_t N = matrix.cols();
    // Return 1.
    if (N == 1)
        return Matrix<std::remove_const_t<T>>(1, 1, 1);
    // Prepare the output matrix.
    Matrix<std::remove_const_t<T>> adj(N, N, 0);
    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {
            // Get cofactor of A[i][j]
            auto support = CofactorView(&matrix, i, j);
            // Sign of adj[j][i] positive if sum of row and column indexes is
            // even. Interchanging rows and columns to get the transpose of the
            // cofactor matrix.
            adj(j, i) = (((i + j) % 2 == 0) ? 1 : -1) * linalg::determinant(support);
        }
    }
    return adj;
}

/// @brief Computes the inverse of this matrix.
/// @param matrix the input matrix.
/// @returns the inverse of the matrix.
template <typename T>
inline auto inverse(const MatrixBase<T> &matrix)
{
    if (!malg::utility::is_square(matrix))
        throw std::invalid_argument("The input matrix is not square.");
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T, double>>;
    // Compute the determinant.
    data_type_t det = linalg::determinant(matrix);
    // If determinant is zero, the matrix is singular.
    if (det == 0.) {
        std::cerr << "Matrix is singular.\n";
        return Matrix<data_type_t>();
    }
    // Find adjoint of the matrix.
    auto adjoint = linalg::adjoint(matrix);
    // Create a matrix for the result.
    Matrix<data_type_t> inv(matrix.rows(), matrix.cols(), 0);
    // Find Inverse using formula "inv(A) = adj(A)/det(A)".
    for (std::size_t r = 0; r < matrix.rows(); ++r)
        for (std::size_t c = 0; c < matrix.cols(); ++c)
            inv(r, c) = adjoint(r, c) / det;
    return inv;
}

/// @brief Computes the QR decomposition of the input matrix.
/// @param A the input matrix.
/// @returns a pair of matrices (Q, R).
template <typename T>
inline auto qr_decomposition(const Matrix<T> &A)
{
    // Select the right type.
    using data_type_t = std::remove_const_t<malg::extract_common_type_t<T, double>>;
    // This function is required to build the householder.
    static auto make_householder = [](const Vector<data_type_t> &a) {
        // Find prependicular vector to mirror.
        auto u = a / (a[0] + std::copysign(malg::square_norm(a), a[0]));
        u[0]   = 1;
        // Finding Householder projection.
        return utility::eye(a.size(), a.size(), data_type_t(1.)) -
               (utility::to_matrix(u, false) * utility::to_matrix(u, true)) * (2. / (u * u));
    };

    // Initialize matrix R.
    Matrix<data_type_t> R(A);
    // Initialize matrix Q.
    Matrix<data_type_t> Q = utility::eye<data_type_t>(R.rows(), R.rows(), data_type_t(1.));
    // Create hoseholder.
    Matrix<data_type_t> H(R.rows(), R.rows(), data_type_t(0.));

    // Compute Q and R.
    for (std::size_t i = 0; i < (R.cols() - (R.cols() == R.rows())); ++i) {
        // Initialize holder.
        for (std::size_t r = 0; r < H.rows(); ++r)
            for (std::size_t c = 0; c < H.cols(); ++c)
                H(r, c) = (r == c) ? 1 : 0;
        // Calculate Householder matrix i: rows and i: columns from R i: rows and ith column
        View(H, i, -1, i, -1) = make_householder(utility::extract_column(R, i, i, -1));
        // Update Q and R.
        Q = Q * H;
        R = H * R;
    }
    // Clean up the lower triangular elements of R.
    for (std::size_t r = 1; r < R.rows(); ++r)
        for (std::size_t c = 0; c < std::min(r, R.cols()); ++c)
            R(r, c) = 0;
    return std::make_pair(std::move(Q), std::move(R));
}

/// @brief Computes the LU decomposition of the input matrix.
/// @param A the input matrix.
/// @returns a pair of matrices (L, U).
template <typename T>
inline auto lu_decomposition(const MatrixBase<T> &A)
{
    if (!malg::utility::is_square(A))
        throw std::invalid_argument("The input matrix is not square.");
    std::size_t N = A.rows();
    Matrix<T> l(N, N, 0);
    Matrix<T> u(N, N, 0);

    std::size_t i, j, k;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (j < i)
                l(j, i) = 0;
            else {
                l(j, i) = A(j, i);
                for (k = 0; k < i; k++) {
                    l(j, i) = l(j, i) - l(j, k) * u(k, i);
                }
            }
        }
        for (j = 0; j < N; j++) {
            if (j < i)
                u(i, j) = 0;
            else if (j == i)
                u(i, j) = 1;
            else {
                u(i, j) = A(i, j) / l(i, i);
                for (k = 0; k < i; k++) {
                    u(i, j) = u(i, j) - ((l(i, k) * u(k, j)) / l(i, i));
                }
            }
        }
    }
    return std::make_pair(l, u);
}

/// @brief Solves the linear system of equations `Ax=b`, basically it finds `x`.
/// @param A the equation matrix (m*n).
/// @param b the vector (m).
/// @returns the vector `x`.
template <typename T>
auto solve(const Matrix<T> &A, const Vector<T> &b)
{
    std::size_t rows = A.rows(), cols = A.cols();
    if (rows == cols) {
        return malg::dot(linalg::inverse(A), b);
    }
    // Overdetermined system.
    if (rows > cols) {
        auto [_q, _r] = linalg::qr_decomposition(A);
        return malg::dot(
            linalg::inverse(View(_r, 0, cols, 0, cols)),
            malg::dot(
                linalg::transpose(View(_q, 0, rows, 0, cols)),
                b));
    }
    // Underdetermined system.
    auto At = linalg::transpose(A);
    return malg::dot(At * linalg::inverse(A * At), b);
}

/// @brief Divides the two matrices.
/// @param A the first matrix.
/// @param B the second matrix.
/// @returns the result of the division.
template <typename T1, typename T2>
inline auto div(const MatrixBase<T1> &A, const MatrixBase<T2> &B)
{
    return A * linalg::inverse(B);
}

/// @brief Computes the exponential of the input matrix.
/// @param A the matrix.
/// @returns the resulting exponential.
template <typename T>
inline auto exp(const MatrixBase<T> &A)
{
    Matrix<std::remove_const_t<T>> result(A.cols(), A.rows(), 0.0);
    for (std::size_t r = 0, c = 0; r < A.rows(); ++r)
        for (c = 0; c < A.cols(); ++c)
            result(c, r) = std::exp(A(r, c));
    return result;
}

/// @brief Copmutes the exponenation of the input matrix.
/// @param A the input matrix.
/// @param accuracy the desired accuracy.
/// @returns the result of the exponential.
template <typename T>
inline auto expm(const MatrixBase<T> &A, double accuracy = 0.00001)
{
    // NOTE: use simplified form of "scale and square"
    // Trade faster coversion in power series for a couple of additional square operations
    // TODO: better/faster algorithm than power series?

    // Scale down matrix A by a power of 2, such that norm(A) < 1.
    const auto [iterations, scale] = malg::log2_ceil(A);
    // Apply the scaling.
    const auto scaled_a = A * scale;

    // Compute power series for e^(A/(2^iterations))
    // init (k = 0)
    const auto batch_size      = A.rows() * A.rows();
    const auto square_accuracy = accuracy * accuracy * scale * scale;
    auto mtk                   = malg::utility::eye<T>(A.rows(), A.cols()); // scaled_a to the power k
    auto ret                   = malg::utility::eye<T>(A.rows(), A.cols()); // sum of power seriees
    auto fac_inv               = double{ 1.0 };                             // inverse faculty
    auto rel_square_diff       = square_accuracy + 1.0;

    for (std::size_t batch_start_idx = 1; (rel_square_diff > square_accuracy && fac_inv != 0.0); batch_start_idx += batch_size) {
        auto local_accum = malg::Matrix<T>(A.rows(), A.cols());
        for (std::size_t i = 0; i < batch_size; ++i) {
            const double k = static_cast<double>(batch_start_idx + i);
            fac_inv        = fac_inv * (1.0 / k);
            if (feq::approximately_equal(fac_inv, 0.0)) {
                break;
            }
            mtk = mtk * scaled_a;
            local_accum += mtk * fac_inv;
        }
        ret += local_accum;
        // Caclulate relative change in this iteration
        // TODO: properly guard against division by zero
        const malg::Matrix<T> rel_error = malg::linalg::div(local_accum * local_accum, ret * ret + accuracy);
        rel_square_diff                 = malg::square_norm(rel_error);
    };
    // raise the result
    for (std::size_t k = 0; k < iterations; ++k) {
        ret = ret * ret;
    }
    return ret;
}

} // namespace malg::linalg
