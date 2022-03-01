/// @file matrix_base.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Base matrix abstract class.

#pragma once

/// @brief If NOT commented, the matrices works in row-major. Otherwise, they
/// work in column-major.
#define ROW_MAJOR

namespace malg
{

template <typename T>
class MatrixBase {
public:
    /// Rows of the matrix.
    unsigned _rows;
    /// Columns of the matrix.
    unsigned _cols;

public:
    using this_type = MatrixBase<T>;
    using base_type = T;

    /// @brief Construct a new Matrix Base object.
    MatrixBase()
        : _rows(0),
          _cols(0)
    {
    }

    /// @brief Construct a new Matrix Base object.
    /// @param rows
    /// @param cols
    MatrixBase(unsigned rows, unsigned cols)
        : _rows(rows),
          _cols(cols)
    {
        // Nothing to do.
    }

    /// @brief Destroy the Matrix object.
    virtual ~MatrixBase() = default;

    /// @brief Get the number of rows of the matrix.
    /// @return the number of rows.
    inline virtual unsigned rows() const
    {
        return _rows;
    }

    /// @brief Get the number of columns of the matrix.
    /// @return the number of columns.
    inline virtual unsigned cols() const
    {
        return _cols;
    }

    /// @brief Returns the total size of the matrix.
    /// @return the total size of the matrix.
    inline virtual unsigned size() const
    {
        return _rows * _cols;
    }

    /// @brief Checks if the matrix is empty.
    /// @return true if it is empty.
    /// @return false otherwise.
    inline virtual bool empty() const
    {
        return (_rows == 0) || (_cols == 0);
    }

    virtual T *data() = 0;

    virtual const T *data() const = 0;

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    virtual T &operator[](unsigned pos) = 0;

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    virtual const T &operator[](unsigned pos) const = 0;

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    virtual T &operator()(unsigned row, unsigned col) = 0;

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    virtual const T &operator()(unsigned row, unsigned col) const = 0;

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    virtual T &at(unsigned row, unsigned col) = 0;

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    virtual const T &at(unsigned row, unsigned col) const = 0;
};

} // namespace malg
