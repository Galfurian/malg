/// @file matrix_base.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Base matrix abstract class.

#pragma once

#include <cstddef>

/// @brief If NOT commented, the matrices works in row-major. Otherwise, they
/// work in column-major.
#define ROW_MAJOR

namespace malg
{

/// @brief The base class for matrix-like structures.
template <typename T>
class MatrixBase {
protected:
    /// Rows of the matrix.
    std::size_t _rows;
    /// Columns of the matrix.
    std::size_t _cols;

public:
    /// The data types of the element of the matrix.
    using value_type = T;
    /// The data types for iterating the element of the vector.
    using iterator = T *;
    /// The data types for iterating the element of the vector.
    using const_iterator = const T *;

    /// @brief Construct a new Matrix Base object.
    MatrixBase()
        : _rows(0),
          _cols(0)
    {
    }

    /// @brief Construct a new Matrix Base object.
    /// @param rows
    /// @param cols
    MatrixBase(std::size_t rows, std::size_t cols)
        : _rows(rows),
          _cols(cols)
    {
        // Nothing to do.
    }

    /// @brief Destroy the Matrix object.
    virtual ~MatrixBase() = default;

    /// @brief Get the number of rows of the matrix.
    /// @returns the number of rows.
    inline virtual std::size_t rows() const
    {
        return _rows;
    }

    /// @brief Get the number of columns of the matrix.
    /// @returns the number of columns.
    inline virtual std::size_t cols() const
    {
        return _cols;
    }

    /// @brief Returns the total size of the matrix.
    /// @returns the total size of the matrix.
    inline virtual std::size_t size() const
    {
        return _rows * _cols;
    }

    /// @brief Checks if the matrix is empty.
    /// @returns true if it is empty.
    /// @returns false otherwise.
    inline virtual bool empty() const
    {
        return (_rows == 0) || (_cols == 0);
    }

    /// @brief Returns a pointer to the internal data.
    /// @return the pointer.
    virtual iterator data() = 0;

    /// @brief Returns a constant pointer to the internal data.
    /// @return the constant pointer.
    virtual const_iterator data() const = 0;

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    virtual T &operator[](std::size_t pos) = 0;

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    virtual const T &operator[](std::size_t pos) const = 0;

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the reference to the accessed item.
    virtual T &operator()(std::size_t row, std::size_t col) = 0;

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the const reference to the accessed item.
    virtual const T &operator()(std::size_t row, std::size_t col) const = 0;

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the reference to the accessed item.
    virtual T &at(std::size_t row, std::size_t col) = 0;

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the const reference to the accessed item.
    virtual const T &at(std::size_t row, std::size_t col) const = 0;
};

} // namespace malg
