/// @file eigen.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Eigenvalues extraction algorithm.
/// @details
/// The original code is posted here:
///     http://cplusplus.com/forum/beginner/220486/2/
/// the author was inspired by these work:
///     http://people.inf.ethz.ch/arbenz/ewp/Lnotes/chapter4.pdf
///     https://www.math.kth.se/na/SF2524/matber15/qrmethod.pdf
///     https://www.ams.org/journals/mcom/2002-71-240/S0025-5718-01-01387-4/S0025-5718-01-01387-4.pdf

#include "malg/utility.hpp"
#include "malg/matrix.hpp"
#include "malg/linalg.hpp"
#include "malg/math.hpp"

#include <complex>
#include <cassert>

namespace malg::eigen
{

// Wilkinson shift in QR algorithm
template <typename T>
inline auto wilkinson_shift(const malg::Matrix<std::complex<T>> &A)
{
    std::complex<T> s(0.0);
    unsigned N = A.rows() - 1;
    if (N > 0) {
        // Bottom-right elements.
        auto a = A(N - 1, N - 1), b = A(N - 1, N), c = A(N, N - 1), d = A(N, N);
        auto delta = std::sqrt((a + d) * (a + d) - 4.0 * (a * d - b * c));
        auto s1    = 0.5 * (a + d + delta);
        auto s2    = 0.5 * (a + d - delta);
        s          = (std::norm(s1 - d) < std::norm(s2 - d) ? s1 : s2);
    }
    return s;
}

// Reduce A to hessenberg form A = P H P-1 where P is unitary, H is hessenberg
//                             i.e. P-1 A P = H
// A hessenberg Matrix is upper triangular plus single non-zero diagonal below main diagonal
template <typename T>
inline auto hessenberg(const malg::Matrix<std::complex<T>> &A)
{
    int N  = A.rows();
    auto P = malg::utility::identity<std::complex<T>>(N);
    auto H = A;
    for (int k = 0; k < N - 2; k++) // k is the working column
    {
        // X vector, based on the elements from k+1 down in the k-th column
        auto xlength = malg::vector_length(H, k, k + 1);

        // U vector ( normalise X - rho.|x|.e_k )
        malg::Vector<std::complex<T>> U(N, 0);
        std::complex<T> rho(1.0, 0.);
        auto xk  = H(k + 1, k);
        auto axk = std::abs(xk);
        if (axk > std::numeric_limits<T>::epsilon())
            rho = -xk / axk;
        U[k + 1]     = xk - rho * xlength;
        auto ulength = std::norm(U[k + 1]);
        for (int i = k + 2; i < N; i++) {
            U[i] = H(i, k);
            ulength += std::norm(U[i]);
        }
        ulength = std::max(std::sqrt(ulength), std::numeric_limits<T>::epsilon());
        for (int i = k + 1; i < N; i++)
            U[i] /= ulength;

        // Householder Matrix: P = I - 2 U U*T
        auto PK = malg::utility::identity<std::complex<T>>(N);
        for (int r = k + 1; r < N; ++r)
            for (int c = k + 1; c < N; ++c)
                PK(r, c) -= 2.0 * U[r] * std::conj(U[c]);

        // Transform as PK*T H PK.   Note: PK is unitary, so PK*T = P
        H = PK * (H * PK);
        P = P * PK;
    }
    return std::make_pair(P, H);
}

// Factorises a hessenberg malg::Matrix<std::complex<T>> A as QR, where Q is unitary and R is upper triangular
// Uses N-1 Givens rotations
template <typename T>
auto qr_factorise_givens(const malg::Matrix<std::complex<T>> &A)
{
    malg::Matrix<std::complex<T>> Q = malg::utility::identity<std::complex<T>>(A.rows());
    malg::Matrix<std::complex<T>> R(A);
    for (unsigned i = 1, j = 0; i < A.rows(); ++i, j = i - 1) {
        // i : The row number
        // j : aiming to zero the element one place below the diagonal.
        if (std::abs(R(i, j)) < std::numeric_limits<T>::epsilon())
            continue;
        // Form the Givens malg::Matrix<std::complex<T>>
        std::complex<T> c = R(j, j);
        std::complex<T> s = -std::conj(R(i, j));
        double length     = std::sqrt(std::norm(c) + std::norm(s));
        c /= length;
        s /= length;
        std::complex<T> cstar            = std::conj(c); //  G*T = ( c* -s )     G = (  c  s  )     <--- j
        std::complex<T> sstar            = std::conj(s); //        ( s*  c )         ( -s* c* )     <--- i
        malg::Matrix<std::complex<T>> RR = R;
        malg::Matrix<std::complex<T>> QQ = Q;
        for (unsigned m = 0; m < A.cols(); m++) {
            R(j, m) = (cstar * RR(j, m)) - (s * RR(i, m));
            R(i, m) = (sstar * RR(j, m)) + (c * RR(i, m));
            Q(m, j) = (c * QQ(m, j)) - (sstar * QQ(m, i));
            Q(m, i) = (s * QQ(m, j)) + (cstar * QQ(m, i));
        }
    }
    return std::make_pair(Q, R);
}

// Apply the QR algorithm to the matrix A.
//
// Multi-stage:
//    - transform to a hessenberg matrix
//    - apply QR factorisation based on Givens rotations
//    - uses (single) Wilkinson shift - double-shift version in development
//
// Should give a Shur decomposition A = P T P-1 where P is unitary, T is upper triangular
//                             i.e. P-1 A P = T
// Eigenvalues of A should be the diagonal elements of T
// If A is hermitian T would be diagonal and the eigenvectors would be the columns of P
template <typename T>
auto QRHessenberg(const malg::Matrix<std::complex<T>> &A)
{
    const int ITERMAX      = 10000;
    const double TOLERANCE = 1e-14;

    int N = A.rows();

    malg::Matrix<std::complex<T>> Told(N, N, 0.);
    malg::Matrix<std::complex<T>> I = malg::utility::identity<std::complex<T>>(N);

    // Stage 1: transform to hessenberg Matrix ( HM = Hessenberg Matrix, UT = unitary transformation )
    auto [HM, UT] = hessenberg(A);

    // Stage 2: apply QR factorisation (using Givens rotations)
    int iter        = 1;
    double residual = 1.0;
    while ((residual > TOLERANCE) && (iter < ITERMAX)) {
        Told = UT;

        // Spectral shift
        std::complex<T> mu = wilkinson_shift(UT);
        if (abs(mu) < std::numeric_limits<T>::epsilon())
            mu = 1.0; // prevent unitary matrices causing a problem
        UT = malg::linear_combination(UT, 1.0, I, -mu);

        // Basic QR algorithm by Givens rotation
        auto [Q, R] = qr_factorise_givens(UT);

        UT = R * Q;
        HM = HM * Q;

        // Reverse shift.
        UT = malg::linear_combination(UT, 1.0, I, mu);

        // Calculate residuals
        // 1. Change on iteration.
        residual = malg::norm(malg::linear_combination(UT, 1.0, Told, -1.0));
        // 2. Below-diagonal elements.
        residual += malg::sub_norm(UT);
        iter++;
        //std::cout << "\nIteration: " << iter << "   Residual: " << residual << std::endl;
    }
    //std::cout << "\nQR iterations: " << iter << "   Residual: " << residual << std::endl;
    if (residual > TOLERANCE) {
        std::cout << "***** WARNING ***** QR algorithm not converged\n";
    }
    return std::make_pair(HM, UT);
}

// Find the eigenvectors of upper-triangular matrix T; returns them as column vectors of matrix E
// The eigenvalues are necessarily the diagonal elements of T
// NOTE: if there are repeated eigenvalues, then THERE MAY NOT BE N EIGENVECTORS
template <typename T>
auto eigenvectorUpper(const malg::Matrix<std::complex<T>> &M)
{
    bool fullset = true;

    int N = M.rows();

    // Columns of E will hold the eigenvectors.
    malg::Matrix<std::complex<T>> E(N, N, 0.);

    malg::Matrix<std::complex<T>> TT = M;
    for (int L = N - 1; L >= 0; L--) // find Lth eigenvector, working from the bottom
    {
        bool ok = true;

        malg::Vector<std::complex<T>> V(N, 0.0);
        std::complex<T> lambda = M(L, L);
        // TT = M - lambda I
        // Solve TT.V = 0
        // free choice of this component
        // back-substitute for other components
        for (int k = 0; k < N; k++)
            TT(k, k) = M(k, k) - lambda;

        V[L] = 1.0;

        for (int i = L - 1; i >= 0; i--) {
            V[i] = 0.0;
            for (int j = i + 1; j <= L; j++) V[i] -= TT(i, j) * V[j];
            if (std::abs(TT(i, i)) < std::numeric_limits<T>::epsilon()) // problem with repeated eigenvalues
            {
                if (std::abs(V[i]) > std::numeric_limits<T>::epsilon())
                    ok = false; // incomplete set; use the lower-L one only
                V[i] = 0.0;
            } else {
                V[i] = V[i] / TT(i, i);
            }
        }

        if (ok) {
            // Normalise
            double length = malg::norm(V);
            for (int i = 0; i <= L; i++)
                E(i, L) = V[i] / length;
        } else {
            fullset = false;
            for (int i = 0; i <= L; i++) E(i, L) = 0.0;
        }
    }

    if (!fullset) {
        std::cout << "\n***** WARNING ***** Can't find N independent eigenvectors\n";
        std::cout << "   Some will be set to zero\n";
    }

    return std::make_pair(E, fullset);
}

} // namespace malg::eigen
