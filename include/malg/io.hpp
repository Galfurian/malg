/// @file io.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Input output functions for matrices.

#pragma once

#include "malg/control/control.hpp"
#include "malg/matrix_base.hpp"
#include "malg/vector.hpp"
#include "malg/view.hpp"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>

/// @brief Returns the longhest value inside the data.
/// @param data the input data.
/// @param precision the desired precision.
/// @returns the number of characters for the longhest value.
template <typename T>
inline auto __get_longhest_value(T &data, std::streamsize precision = 6)
{
    std::stringstream ss;
    ss.precision(precision);
    int longhest = 0, length;
    for (std::size_t i = 0; i < data.size(); ++i) {
        ss.str("");
        ss << data[i];
        length = static_cast<int>(ss.str().length());
        if (longhest < length)
            longhest = length;
    }
    return longhest;
}

/// @brief Stream operator for Vector.
template <typename T>
inline std::ostream &operator<<(std::ostream &lhs, const malg::Vector<T> &rhs)
{
    std::size_t i;
    // Get the longhest value.
    auto length = __get_longhest_value(rhs, lhs.precision());
    // Print the vector.
    for (i = 0; i < rhs.size(); ++i) {
        lhs << std::setw(length) << rhs[i];
        if (i < (rhs.size() - 1))
            lhs << " ";
    }
    return lhs;
}

/// @brief Stream operator for MatrixBase.
template <typename T>
inline std::ostream &operator<<(std::ostream &lhs, const malg::MatrixBase<T> &rhs)
{
    std::size_t r, c;
    // Get the longhest value.
    auto length = __get_longhest_value(rhs, lhs.precision());
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

/// @brief Stream operator for writing Vector on file.
template <typename T>
inline std::ofstream &operator<<(std::ofstream &lhs, const malg::Vector<T> &rhs)
{
    std::size_t i;
    // Get the longhest value.
    auto length = __get_longhest_value(rhs, lhs.precision());
    // Print the vector size.
    lhs << "V " << rhs.size() << "\n";
    // Print the vector.
    for (i = 0; i < rhs.size(); ++i) {
        lhs << std::setw(length) << rhs[i];
        if (i < (rhs.size() - 1))
            lhs << " ";
    }
    lhs << "\n";
    return lhs;
}

/// @brief Stream operator for writing MatrixBase on file.
template <typename T>
inline std::ofstream &operator<<(std::ofstream &lhs, const malg::MatrixBase<T> &rhs)
{
    std::size_t r, c;
    // Get the longhest value.
    auto length = __get_longhest_value(rhs, lhs.precision());
    // Print the matrix size.
    lhs << "M " << rhs.rows() << " " << rhs.cols() << "\n";
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

/// @brief Stream operator for Range.
inline std::ostream &operator<<(std::ostream &lhs, const malg::Range &rhs)
{
    lhs << "[" << rhs.start << ", " << rhs.stop << "]";
    return lhs;
}

/// @brief Stream operator for reading a Vector from file.
template <typename T>
inline std::ifstream &operator>>(std::ifstream &lhs, malg::Vector<T> &rhs)
{
    std::size_t size, i;
    char type;
    // Read the type.
    lhs >> type;
    if (type != 'V') {
        std::cerr << "The file does not contain a valid vector.\n";
        return lhs;
    }
    // Read the size.
    lhs >> size;
    // Prepare the vector.
    rhs.resize(size);
    // Read the vector.
    for (i = 0; i < rhs.size(); ++i)
        lhs >> rhs[i];
    return lhs;
}

/// @brief Stream operator for reading a Matrix from file.
template <typename T>
inline std::ifstream &operator>>(std::ifstream &lhs, malg::Matrix<T> &rhs)
{
    std::size_t rows, cols, r, c;
    char type;
    // Read the type.
    lhs >> type;
    if (type != 'M') {
        std::cerr << "The file does not contain a valid matrix.\n";
        return lhs;
    }
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
/// @returns the original stream.
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
/// @returns the original stream.
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

/// @brief Prepares a string for printing the vector.
/// @param name the name to show.
/// @param v the vector to display.
/// @returns the output string.
template <typename T>
inline std::string dump_vector(const char *name, const malg::Vector<T> &v)
{
    std::stringstream ss;
    ss << name << " = " << v;
    return ss.str();
}

/// @brief Prepares a string for printing the matrix.
/// @param name the name to show.
/// @param m the matrix to display.
/// @returns the output string.
template <typename T>
inline std::string dump_matrix(const char *name, const malg::MatrixBase<T> &m)
{
    std::stringstream ss;
    ss << name << " = \n"
       << m;
    return ss.str();
}

/// @brief Prepares a matlab declaration of the input vector.
/// @param name the name to show.
/// @param v the vector to transform.
/// @returns the output string.
template <typename T>
inline std::string to_matlab(const char *name, const malg::Vector<T> &v)
{
    std::stringstream ss;
    ss << name << " = [";
    for (std::size_t i = 0; i < v.size(); ++i) {
        ss << v[i];
        if (i < (v.size() - 1))
            ss << " ";
    }
    ss << "]";
    return ss.str();
}

/// @brief Prepares a matlab declaration of the input matrix.
/// @param name the name to show.
/// @param m the matrix to transform.
/// @returns the output string.
template <typename T>
inline std::string to_matlab(const char *name, const malg::MatrixBase<T> &m)
{
    std::stringstream ss;
    ss << name << " = [";
    for (std::size_t r = 0, c; r < m.rows(); ++r) {
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

/// @brief Prepares a c++ declaration of the input vector.
/// @param name the name to show.
/// @param type the declaration type to use.
/// @param v the vector to transform.
/// @returns the output string.
template <typename T>
inline std::string to_cpp(const char *name, const char *type, const malg::Vector<T> &v)
{
    std::stringstream ss;
    ss << "malg::Vector<" << type << "> " << name << " = {";
    for (std::size_t i = 0; i < v->size(); ++i) {
        if constexpr (malg::is_complex_v<T>) {
            ss << v->operator[](i).real() << std::showpos << v->operator[](i).imag() << "i" << std::noshowpos;
        } else {
            ss << v->operator[](i);
        }
        if (i < (v->size() - 1))
            ss << ", ";
    }
    ss << "}";
    return ss.str();
}

/// @brief Prepares a c++ declaration of the input matrix.
/// @param name the name to show.
/// @param type the declaration type to use.
/// @param m the matrix to transform.
/// @returns the output string.
template <typename T>
inline std::string to_cpp(const char *name, const char *type, const malg::MatrixBase<T> &m)
{
    std::stringstream ss;
    ss << "malg::Matrix<" << type << "> " << name << " = {";
    for (std::size_t r = 0, c; r < m->rows(); ++r) {
        ss << "{";
        for (c = 0; c < m->cols(); ++c) {
            if constexpr (malg::is_complex_v<T>) {
                ss << m->operator()(r, c).real() << std::showpos << m->operator()(r, c).imag() << "i" << std::noshowpos;
            } else {
                ss << m->operator()(r, c);
            }
            if (c < (m->cols() - 1))
                ss << ",";
        }
        ss << "}";
        if (r < (m->rows() - 1))
            ss << ",";
    }
    ss << "}";
    return ss.str();
}

} // namespace malg