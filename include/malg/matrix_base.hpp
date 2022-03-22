/// @file matrix_base.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Base matrix abstract class.

#pragma once

/// @brief If NOT commented, the matrices works in row-major. Otherwise, they
/// work in column-major.
#define ROW_MAJOR

namespace malg
{

using size_type_t = unsigned;

template <typename T>
class MatrixBase {
public:
    /// Rows of the matrix.
    size_type_t _rows;
    /// Columns of the matrix.
    size_type_t _cols;

public:
    using this_type = MatrixBase<T>;
    using base_type = T;

    /// @brief Construct a new Matrix Base object.
    constexpr MatrixBase() noexcept
        : _rows(0),
          _cols(0)
    {
    }

    /// @brief Construct a new Matrix Base object.
    /// @param rows
    /// @param cols
    constexpr MatrixBase(size_type_t rows, size_type_t cols) noexcept
        : _rows(rows),
          _cols(cols)
    {
        // Nothing to do.
    }

    /// @brief Destroy the Matrix object.
    virtual ~MatrixBase() = default;

    /// @brief Get the number of rows of the matrix.
    /// @return the number of rows.
    constexpr inline virtual size_type_t rows() const noexcept
    {
        return _rows;
    }

    /// @brief Get the number of columns of the matrix.
    /// @return the number of columns.
    constexpr inline virtual size_type_t cols() const noexcept
    {
        return _cols;
    }

    /// @brief Returns the total size of the matrix.
    /// @return the total size of the matrix.
    constexpr inline virtual size_type_t size() const noexcept
    {
        return _rows * _cols;
    }

    /// @brief Checks if the matrix is empty.
    /// @return true if it is empty.
    /// @return false otherwise.
    constexpr inline virtual bool empty() const noexcept
    {
        return (_rows == 0) || (_cols == 0);
    }

    virtual T *data() = 0;

    virtual const T *data() const = 0;

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    virtual T &operator[](size_type_t pos) = 0;

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    virtual const T &operator[](size_type_t pos) const = 0;

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    virtual T &operator()(size_type_t row, size_type_t col) = 0;

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    virtual const T &operator()(size_type_t row, size_type_t col) const = 0;

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    virtual T &at(size_type_t row, size_type_t col) = 0;

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    virtual const T &at(size_type_t row, size_type_t col) const = 0;
};

} // namespace malg
