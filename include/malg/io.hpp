/// @file io.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Input output functions for matrices.

#pragma once

#include "malg/matrix_base.hpp"
#include "malg/vector.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

template <typename T>
inline unsigned __get_longhest_value(T &ds, unsigned precision = 6)
{
    std::stringstream ss;
    ss.precision(precision);
    unsigned longhest = 0, length;
    for (unsigned i = 0; i < ds.size(); ++i) {
        ss.str("");
        ss << ds[i];
        length = ss.str().length();
        if (longhest < length)
            longhest = length;
    }
    return longhest;
}

template <typename T>
inline std::ostream &operator<<(std::ostream &lhs, const malg::Vector<T> &rhs)
{
    // Get the longhest value.
    unsigned i, length = __get_longhest_value(rhs, lhs.precision());
    // Print the vector.
    for (i = 0; i < rhs.size(); ++i) {
        lhs << std::setw(length) << rhs[i];
        if (i < (rhs.size() - 1))
            lhs << " ";
    }
    return lhs;
}

template <typename T>
inline std::ostream &operator<<(std::ostream &lhs, const malg::MatrixBase<T> &rhs)
{
    // Get the longhest value.
    unsigned r, c, length = __get_longhest_value(rhs, lhs.precision());
    // Print the matrix.
    for (r = 0; r < rhs.rows(); ++r) {
        for (c = 0; c < rhs.cols(); ++c) {
            lhs << std::setw(length) << rhs(r, c);
            if (c < (rhs.cols() - 1))
                lhs << " ";
        }
        if (r < (rhs.rows() - 1))
            lhs << "\n";
    }
    return lhs;
}

template <typename T>
inline std::ofstream &operator<<(std::ofstream &lhs, const malg::Vector<T> &rhs)
{
    // Get the longhest value.
    unsigned i, length = __get_longhest_value(rhs, lhs.precision());
    // Print the vector size.
    lhs << rhs.size() << "\n";
    // Print the vector.
    for (i = 0; i < rhs.size(); ++i) {
        lhs << std::setw(length) << rhs[i];
        if (i < (rhs.size() - 1))
            lhs << " ";
    }
    lhs << "\n";
    return lhs;
}

template <typename T>
inline std::ofstream &operator<<(std::ofstream &lhs, const malg::MatrixBase<T> &rhs)
{
    // Get the longhest value.
    unsigned r, c, length = __get_longhest_value(rhs, lhs.precision());
    // Print the matrix size.
    lhs << rhs.rows() << " " << rhs.cols() << "\n";
    // Print the matrix.
    for (r = 0; r < rhs.rows(); ++r) {
        for (c = 0; c < rhs.cols(); ++c) {
            lhs << std::setw(length) << rhs(r, c);
            if (c < (rhs.cols() - 1))
                lhs << " ";
        }
        if (r < (rhs.rows() - 1))
            lhs << "\n";
    }
    lhs << "\n";
    return lhs;
}

template <typename T>
inline std::ifstream &operator>>(std::ifstream &lhs, malg::Vector<T> &rhs)
{
    unsigned size, i;
    // Read the size.
    lhs >> size;
    // Prepare the vector.
    rhs.resize(size);
    // Read the vector.
    for (i = 0; i < rhs.size(); ++i)
        lhs >> rhs[i];
    return lhs;
}

template <typename T>
inline std::ifstream &operator>>(std::ifstream &lhs, malg::Matrix<T> &rhs)
{
    unsigned rows, cols, r, c;
    // Read the size.
    lhs >> rows;
    lhs >> cols;
    // Prepare the matrix.
    rhs.resize(rows, cols);
    // Read the matrix.
    for (r = 0; r < rhs.rows(); ++r)
        for (c = 0; c < rhs.cols(); ++c)
            lhs >> rhs(r, c);
    return lhs;
}

/// @brief Output stream function.
/// @param lhs the stream.
/// @param rhs the state space model.
/// @return the original stream.
template <typename T>
inline std::ostream &operator<<(std::ostream &lhs, const malg::control::StateSpace<T> &rhs)
{
    lhs << "A =\n"
        << rhs.A << "\n"
        << "B =\n"
        << rhs.B << "\n"
        << "C =\n"
        << rhs.C << "\n"
        << "D =\n"
        << rhs.D << "\n"
        << "Continuous-time state-space model.\n";
    return lhs;
}

/// @brief Output stream function.
/// @param lhs the stream.
/// @param rhs the state space model.
/// @return the original stream.
template <typename T>
inline std::ostream &operator<<(std::ostream &lhs, const malg::control::DiscreteStateSpace<T> &rhs)
{
    lhs << "A =\n"
        << rhs.A << "\n"
        << "B =\n"
        << rhs.B << "\n"
        << "C =\n"
        << rhs.C << "\n"
        << "D =\n"
        << rhs.D << "\n"
        << "Sample time: " << rhs.sample_time << " seconds\n"
        << "Discrete-time state-space model.\n";
    return lhs;
}

namespace malg
{

template <typename T>
inline std::string dump_matrix(const char *name, const malg::MatrixBase<T> &m)
{
    std::stringstream ss;
    ss << name << " = \n"
       << m;
    return ss.str();
}

template <typename T>
inline std::string dump_vector(const char *name, const malg::Vector<T> &v)
{
    std::stringstream ss;
    ss << name << " = " << v;
    return ss.str();
}

template <typename T>
inline std::string to_matlab(const char *name, const malg::Vector<T> &v)
{
    std::stringstream ss;
    ss << name << " = [";
    for (unsigned i = 0; i < v.size(); ++i) {
        ss << v[i];
        if (i < (v.size() - 1))
            ss << " ";
    }
    ss << "]";
    return ss.str();
}

template <typename T>
inline std::string to_matlab(const char *name, const malg::MatrixBase<T> &m)
{
    std::stringstream ss;
    ss << name << " = [";
    for (unsigned r = 0, c; r < m.rows(); ++r) {
        ss << "[";
        for (c = 0; c < m.cols(); ++c) {
            if constexpr (malg::is_complex_v<T>) {
                ss << m(r, c).real() << std::showpos << m(r, c).imag() << "i" << std::noshowpos;
            } else {
                ss << m(r, c);
            }
            if (c < (m.cols() - 1))
                ss << " ";
        }
        ss << "];";
    }
    ss << "]";
    return ss.str();
}

template <typename T>
inline std::string to_cpp(const char *name, const malg::Vector<T> &v)
{
    std::stringstream ss;
    ss << "malg::Vector<double> " << name << " = {";
    for (unsigned i = 0; i < v.size(); ++i) {
        ss << v[i];
        if (i < (v.size() - 1))
            ss << ", ";
    }
    ss << "};\n";
    return ss.str();
}

template <typename T>
inline std::string to_cpp(const char *name, const malg::MatrixBase<T> &m)
{
    std::stringstream ss;
    ss << "malg::Matrix<double> " << name << " = {";
    for (unsigned r = 0, c; r < m.rows(); ++r) {
        ss << "{";
        for (c = 0; c < m.cols(); ++c) {
            ss << m(r, c);
            if (c < (m.cols() - 1))
                ss << ",";
        }
        ss << "}";
        if (r < (m.rows() - 1))
            ss << ",";
    }
    ss << "}";
    return ss.str();
}

} // namespace malg