/// @file control.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Control algorithms and data-structures.

#pragma once

#include "malg/matrix.hpp"
#include "malg/linalg.hpp"

#include <tuple>

namespace malg::control
{

/// @brief Continuous-time state space model.
template <typename T>
class StateSpace {
public:
    /// System matrix.
    Matrix<T> A;
    /// Input matrix.
    Matrix<T> B;
    /// Output matrix.
    Matrix<T> C;
    /// Feedforward matrix.
    Matrix<T> D;

    /// @brief Construct a new State Space object.
    StateSpace()
        : A(),
          B(),
          C(),
          D()
    {
        // Nothing to do.
    }

    /// @brief Construct a new State Space object, with the given matrices.
    StateSpace(const Matrix<T> &_A, const Matrix<T> &_B, const Matrix<T> &_C, const Matrix<T> &_D)
        : A(_A),
          B(_B),
          C(_C),
          D(_D)
    {
        // Nothing to do.
    }
};

/// @brief Discrete-time state space model.
template <typename T>
class DiscreteStateSpace : public StateSpace<T> {
public:
    /// Sample time used for the discretization.
    T sample_time;

    /// @brief Construct a new Discrete State Space object.
    DiscreteStateSpace()
        : StateSpace<T>(),
          sample_time(.0)
    {
        // Nothing to do.
    }

    /// @brief Construct a new Discrete State Space object, with the given matrices.
    DiscreteStateSpace(const Matrix<T> &_A, const Matrix<T> &_B, const Matrix<T> &_C, const Matrix<T> &_D, T _sample_time)
        : StateSpace<T>(_A, _B, _C, _D),
          sample_time(_sample_time)
    {
        // Nothing to do.
    }
};

/// @brief Discretize the state space model, with the given time-step.
/// @param sys the continuous-time state space model.
/// @param sample_time the sample time step.
/// @returns the discretized state space.
template <typename T>
inline auto c2d(const StateSpace<T> &sys, T sample_time)
{
    // Create the discrete time state-space model.
    DiscreteStateSpace<T> dsys;
    // Compute the discretized matrices.
    dsys.A = linalg::expm(sys.A * sample_time);
    dsys.B = linalg::inverse(sys.A) * (dsys.A - utility::identity<T>(sys.A.rows())) * sys.B;
    dsys.C = sys.C;
    dsys.D = sys.D;
    // Set the sampling period.
    dsys.sample_time = sample_time;
    // Check the dimensions.
    assert(dsys.A.rows() == dsys.A.cols());
    assert(dsys.B.rows() == dsys.A.rows());
    assert(dsys.C.cols() == dsys.A.rows());
    return dsys;
}

/// @brief Simulate one step with the given discretized system.
/// @param sys the system to simulate.
/// @param x the current state.
/// @param u the current input.
/// @returns a tuple containing the next state, and the output.
template <typename T>
inline std::pair<Vector<T>, Vector<T>> simulate_step(
    const DiscreteStateSpace<T> &sys,
    const Vector<T> &x,
    const Vector<T> &u)
{
    return std::make_pair(
        dot(sys.A, x) + dot(sys.B, u),
        dot(sys.C, x) + dot(sys.D, u));
}

/// @brief Simulate the discrete time state-space model.
/// @param sys the discrete time state-space model.
/// @param u the inputs to provide.
/// @param simulated_time the total simulated time.
/// @param init the initial state.
/// @returns auto
template <typename T>
inline auto simulate(
    const DiscreteStateSpace<T> &sys,
    const std::vector<Vector<T>> &u,
    T simulated_time,
    Vector<T> init)
{
    // Compute the simulation steps.
    auto steps = std::min(u.size(), static_cast<size_t>(simulated_time / sys.sample_time));
    // Prepare a vector for the time vector.
    std::vector<T> t(steps + 1, 0.);
    // Prepare a vector for the state vector, specifically just the initial value.
    std::vector<Vector<T>> x(1, init);
    // Prepare a vector for the output.
    std::vector<Vector<T>> y;
    // Perform the simulation.
    double rtime = 0.;
    for (std::size_t i = 0; i < steps; ++i, rtime += sys.sample_time) {
        t[i]      = rtime;
        auto step = simulate_step(sys, x.back(), u[i]);
        x.emplace_back(step.first);
        y.emplace_back(step.second);
    }
    return std::make_tuple(t, x, y);
}

/// @brief Controllabilty matrix.
/// @param A State matrix of the system (NxN).
/// @param B Input matrix of the system (NxQ).
/// @returns Controllability matrix.
/// @details CM = [ B    A*B    A^2*B    ...    A^(N-1)*B ]
template <typename T>
inline auto ctrb(const MatrixBase<T> &A, const MatrixBase<T> &B)
{
    if (!utility::is_square(A)) {
        throw std::runtime_error("ctrb: matrix A must be square.");
    }
    if (A.rows() != B.rows()) {
        throw std::runtime_error("ctrb: A and B matrices dimensions doesn't match.");
    }
    // Create the controllability matrix.
    Matrix<T> result = B;
    // Construct the controllability matrix.
    for (std::size_t i = 1; i < A.rows(); ++i)
        result = utility::hstack(result, linalg::powm(A, i) * B);
    return result;
}

/// @brief Observability matrix.
/// @param A State matrix of the system (NxN).
/// @param C Output matrix of the system (PxN).
/// @returns Observability matrix.
/// @details
///      | C         |
///      | C*A       |
/// OM = | C*A^2     |
///      | ...       |
///      | C*A^(n-1) |
template <typename T>
inline auto obsv(const MatrixBase<T> &A, const MatrixBase<T> &C)
{
    if (!utility::is_square(A)) {
        throw std::runtime_error("ctrb: matrix A must be square.");
    }
    if (A.rows() != C.cols()) {
        throw std::runtime_error("ctrb: A and C matrices dimensions doesn't match.");
    }
    // Create the observability matrix.
    Matrix<T> result = C;
    // Construct the observability matrix.
    for (std::size_t i = 1; i < A.rows(); ++i)
        result = utility::vstack(result, C * linalg::powm(A, i));
    return result;
}

/// @brief Computes the coefficients of the polynomial whose roots are the elements of a.
/// @param a the input vector.
/// @returns coefficients of the polynomial.
template <typename T>
inline auto poly(const Vector<T> &a)
{
    Vector<T> c(a.size() + 1);
    c[0] = 1;
    for (std::size_t j = 0; j < a.size(); ++j)
        for (std::size_t i = j + 1; i >= 1; --i)
            c[i] -= a[j] * c[i - 1];
    return c;
}

/// @brief Reduce a polynomial coefficient vector to a minimum number of terms
/// by stripping off any leading zeros.
/// @param a the input coefficient vector.
/// @returns the reduced vector.
template <typename T>
inline auto polyreduce(const Vector<T> &a)
{
    std::size_t i, j;
    for (j = 0; j < a.size(); ++j)
        if (a[j] != 0)
            break;
    if (j == a.size())
        return Vector<T>({ 0 });
    Vector<T> c(a.size() - j);
    for (i = j; i < a.size(); ++i)
        c[i - j] = a[i];
    return c;
}

/// @brief Computes the coefficients of the polynomial.
/// @param A Could be a row/column vector, or a matrix.
/// @returns coefficients of the polynomial.
template <typename T>
inline auto poly(const MatrixBase<T> &A)
{
    // An empty matrix.
    if ((A.size() == 0))
        return Vector<T>{ 1 };

    // A row/column vector.
    if ((A.rows() == 1) || (A.cols() == 1)) {
        Vector<T> c(A.size() + 1);
        c[0] = 1;
        for (std::size_t j = 0; j < A.size(); ++j)
            for (std::size_t i = j + 1; i >= 1; --i)
                c[i] -= A[j] * c[i - 1];
        return c;
    }
    // If we are dealing with a matrix, A must be square.
    if (!utility::is_square(A)) {
        throw std::runtime_error("poly: Matrix A must be square.");
    }
    auto n = A.rows();
    Vector<T> c(n + 1);
    Matrix<T> B = A;
    auto last_c = trace(B);
    c[0]        = 1;
    c[1]        = -last_c;
    auto I      = utility::identity<T>(B.rows());
    for (std::size_t m = 2; m < (n + 1); ++m) {
        B      = A * (B - (I * last_c));
        last_c = trace(B) / static_cast<double>(m);
        c[m]   = -last_c;
    }
    return c;
}

/// @brief Pole placement using Ackermann method.
/// @param A state matrix of the system.
/// @param B input matrix of the system.
/// @param poles for generating the desired closed-loop behaviour.
/// @returns Gains such that A - BK has given eigenvalues.
template <typename T>
inline auto acker(const MatrixBase<T> &A, const MatrixBase<T> &B, const MatrixBase<T> &poles)
{
    if (!utility::is_square(A) || !utility::is_column_vector(B) || (A.rows() != B.rows())) {
        std::cerr << "acker: matrix A and vector B not conformal.\n";
        exit(1);
    }
    // poles could be complex.
    if (!utility::is_row_vector(poles) || poles.empty()) {
        std::cerr << "acker: poles must be a vector.\n";
        exit(1);
    }
    // Make sure the system is controllable.
    auto ct = ctrb<T>(A, B);
    if (linalg::determinant(ct) == 0) {
        std::cerr << "acker: system not reachable, pole placement invalid.\n";
        exit(1);
    }
    // Compute the desired characteristic polynomial.
    auto p = poly(poles);
    // Place the poles using Ackermann's method.
    auto Ap = utility::zeros<T>(A.rows(), A.cols());
    for (std::size_t i = 0; i < (A.rows() + 1); ++i) {
        Ap += linalg::powm(A, (A.cols() - i)) * p[i];
    }
    // Prepare the selection matrix.
    auto selection = utility::zeros<T>(1, ct.rows());
    // We are interested only in the last row.
    selection[ct.rows() - 1] = 1;
    // Compute the final matrix and select only the last row.
    return selection * linalg::inverse(ct) * Ap;
}

/// @brief Generates an empty input vector.
/// @param sys The system for which we want to generate the input.
/// @returns The empty input vector.
template <typename T>
inline auto generate_input_zero(const StateSpace<T> &sys)
{
    return Vector<T>(sys.B.cols());
}

/// @brief Generates an empty state vector.
/// @param sys The system for which we want to generate the state.
/// @returns The empty state vector.
template <typename T>
inline auto generate_state_zero(const StateSpace<T> &sys)
{
    return Vector<T>(sys.A.rows());
}

} // namespace malg::control
