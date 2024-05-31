/// @file io.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Input output functions for matrices.

#pragma once

#include "malg/fixed_size/vector.hpp"
#include "malg/fixed_size/matrix.hpp"
#include "malg/fixed_size/view.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

namespace malg::details
{

/// @brief Returns the longhest value inside the data.
/// @param data the input data.
/// @param precision the desired precision.
/// @returns the number of characters for the longhest value.
template <typename T>
inline auto get_longhest_value(T &data, std::streamsize precision)
{
    std::stringstream ss;
    ss.precision(precision);
    int longhest = 0, length;
    for (std::size_t i = 0; i < data.size(); ++i) {
        ss.str("");
        ss << data[i];
        length = static_cast<int>(ss.str().length());
        if (longhest < length) {
            longhest = length;
        }
    }
    return longhest;
}

} // namespace malg::details

/// @brief Stream operator for Vector.
template <typename T, std::size_t N>
std::ostream &operator<<(std::ostream &lhs, const malg::Vector<T, N> &rhs)
{
    std::size_t i;
    // Get the longhest value.
    auto length = malg::details::get_longhest_value(rhs, lhs.precision());
    // Print the vector.
    for (i = 0; i < rhs.size(); ++i) {
        lhs << std::setw(length) << rhs[i];
        if (i < (rhs.size() - 1)) {
            lhs << " ";
        }
    }
    return lhs;
}

/// @brief Stream operator for MatrixBase.
template <typename T, std::size_t N1, std::size_t N2>
std::ostream &operator<<(std::ostream &lhs, const malg::Matrix<T, N1, N2> &rhs)
{
    std::size_t r = 0, c = 0;
    // Get the longhest value.
    auto length = malg::details::get_longhest_value(rhs, lhs.precision());
    // Print the matrix.
    for (r = 0; r < N1; ++r) {
        for (c = 0; c < N2; ++c) {
            lhs << std::setw(length) << rhs[r][c];
            if (c < (N2 - 1)) {
                lhs << " ";
            }
        }
        if (r < (N1 - 1)) {
            lhs << "\n";
        }
    }
    return lhs;
}

/// @brief Stream operator for VectorView.
template <typename T, std::size_t N, std::size_t Start = 0, std::size_t End = N>
std::ostream &operator<<(std::ostream &lhs, const malg::VectorView<T, N, Start, End> &rhs)
{
    std::size_t i;
    // Get the longhest value.
    auto length = malg::details::get_longhest_value(rhs, lhs.precision());
    // Print the vector.
    for (i = 0; i < rhs.size(); ++i) {
        lhs << std::setw(length) << rhs[i];
        if (i < (rhs.size() - 1)) {
            lhs << " ";
        }
    }
    return lhs;
}

/// @brief Stream operator for MatrixBase.
template <typename T, std::size_t N1, std::size_t N2, std::size_t RowStart, std::size_t RowEnd, std::size_t ColStart, std::size_t ColEnd>
std::ostream &operator<<(std::ostream &lhs, const malg::MatrixView<T, N1, N2, RowStart, RowEnd, ColStart, ColEnd> &rhs)
{
    std::size_t r = 0, c = 0;
    // Get the longhest value.
    auto length = malg::details::get_longhest_value(rhs, lhs.precision());
    // Print the matrix.
    for (r = 0; r < (RowEnd - RowStart); ++r) {
        for (c = 0; c < (ColEnd - ColStart); ++c) {
            lhs << std::setw(length) << rhs[r][c];
            if (c < ((ColEnd - ColStart) - 1)) {
                lhs << " ";
            }
        }
        if (r < ((RowEnd - RowStart) - 1)) {
            lhs << "\n";
        }
    }
    return lhs;
}

/// @brief Stream operator for writing Vector on file.
template <typename T, std::size_t N>
std::ofstream &operator<<(std::ofstream &lhs, const malg::Vector<T, N> &rhs)
{
    std::size_t i;
    // Get the longhest value.
    auto length = malg::details::get_longhest_value(rhs, lhs.precision());
    // Print the vector size.
    lhs << "V " << rhs.size() << "\n";
    // Print the vector.
    for (i = 0; i < rhs.size(); ++i) {
        lhs << std::setw(length) << rhs[i];
        if (i < (rhs.size() - 1)) {
            lhs << " ";
        }
    }
    lhs << "\n";
    return lhs;
}

/// @brief Stream operator for writing MatrixBase on file.
template <typename T, std::size_t N1, std::size_t N2>
std::ofstream &operator<<(std::ofstream &lhs, const malg::Matrix<T, N1, N2> &rhs)
{
    std::size_t r = 0, c = 0;
    // Get the longhest value.
    auto length = malg::details::get_longhest_value(rhs, lhs.precision());
    // Print the matrix size.
    lhs << "M " << N1 << " " << N2 << "\n";
    // Print the matrix.
    for (r = 0; r < N1; ++r) {
        for (c = 0; c < N2; ++c) {
            lhs << std::setw(length) << rhs[r][c];
            if (c < (N2 - 1)) {
                lhs << " ";
            }
        }
        if (r < (N1 - 1)) {
            lhs << "\n";
        }
    }
    lhs << "\n";
    return lhs;
}

/// @brief Stream operator for reading a Vector from file.
template <typename T, std::size_t N>
std::ifstream &operator>>(std::ifstream &lhs, malg::Vector<T, N> &rhs)
{
    // Read the type.
    char type;
    lhs >> type;
    if (type != 'V') {
        std::cerr << "The file does not contain a valid vector (not a vector).\n";
        return lhs;
    }
    // Read the size.
    std::size_t size;
    lhs >> size;
    if (size != N) {
        std::cerr << "The file does not contain a valid vector (wrong size).\n";
        return lhs;
    }
    // Read the vector.
    for (std::size_t i = 0; i < N; ++i) {
        lhs >> rhs[i];
    }
    return lhs;
}

/// @brief Stream operator for reading a Matrix from file.
template <typename T, std::size_t N1, std::size_t N2>
std::ifstream &operator>>(std::ifstream &lhs, malg::Matrix<T, N1, N2> &rhs)
{
    // Read the type.
    char type;
    lhs >> type;
    if (type != 'M') {
        std::cerr << "The file does not contain a valid matrix.\n";
        return lhs;
    }
    // Read the rows.
    std::size_t rows;
    lhs >> rows;
    if (rows != N1) {
        std::cerr << "The file does not contain a valid vector (wrong rows).\n";
        return lhs;
    }
    // Read the cols.
    std::size_t cols;
    lhs >> cols;
    if (cols != N2) {
        std::cerr << "The file does not contain a valid vector (wrong cols).\n";
        return lhs;
    }
    // Read the matrix.
    for (std::size_t r = 0; r < N1; ++r) {
        for (std::size_t c = 0; c < N2; ++c) {
            lhs >> rhs[r][c];
        }
    }
    return lhs;
}
