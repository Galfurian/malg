/// @file matrix.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief The matrix class.

#pragma once

#include "malg/type_traits.hpp"
#include "malg/matrix_base.hpp"
#include "malg/vector.hpp"

#include <initializer_list>
#include <algorithm>
#include <vector>

namespace malg
{

struct Range;

template <typename MatrixType>
class View;

/// @brief MatrixArray structure.
template <typename T, std::size_t ROWS, std::size_t COLS>
class MatrixArray : public MatrixBase<T> {
public:
    /// The internal data.
    std::array<T, ROWS * COLS> _data;

public:
    /// The data types of the element of the matrix.
    using value_type = T;

    /// @brief Construct a new MatrixArray object.
    /// @param rows
    /// @param cols
    /// @param initial
    constexpr MatrixArray() noexcept
        : MatrixBase<T>(ROWS, COLS),
          _data()
    {
        // Nothing to do.
    }

    /// @brief Construct a new MatrixArray object.
    /// @param rows
    /// @param cols
    /// @param initial
    constexpr MatrixArray(const T &initial) noexcept
        : MatrixBase<T>(ROWS, COLS),
          _data()
    {
        for (std::size_t i = 0; i < _data.size(); ++i)
            _data[i] = initial;
    }

    /// @brief Construct a new MatrixArray object.
    /// @param data
    constexpr MatrixArray(const std::initializer_list<std::initializer_list<T>> &data)
        : MatrixBase<T>(ROWS, COLS),
          _data()
    {
        // Check the number of rows.
        if (data.size() != ROWS)
            throw std::invalid_argument("The input list has an invalid number of rows.");
        // Check the number of columns.
        if (data.begin()->size() != COLS)
            throw std::invalid_argument("The input list has an invalid number of columns.");
        // Get an interator for the data.
        auto r_it = data.begin();
        for (std::size_t r = 0, c = 0; r < this->rows(); ++r, ++r_it) {
            auto c_it = (*r_it).begin();
            for (c = 0; c < this->cols(); ++c, ++c_it) {
                this->at(r, c) = *c_it;
            }
        }
    }

    /// @brief Construct a new MatrixArray object.
    /// @param rhs
    constexpr MatrixArray(const MatrixBase<T> &rhs)
        : MatrixBase<T>(ROWS, COLS),
          _data()
    {
        // Check the number of rows.
        if (rhs.rows() != ROWS)
            throw std::invalid_argument("The input list has an invalid number of rows.");
        // Check the number of columns.
        if (rhs.cols() != COLS)
            throw std::invalid_argument("The input list has an invalid number of columns.");
        // Copy the data.
        for (std::size_t r = 0; r < ROWS; ++r)
            for (std::size_t c = 0; c < COLS; ++c)
                this->at(r, c) = rhs(r, c);
    }

    /// @brief Construct a new MatrixArray object.
    template <typename T2>
    constexpr MatrixArray(const MatrixBase<T2> &rhs)
        : MatrixBase<T>(ROWS, COLS),
          _data()
    {
        // Check the number of rows.
        if (rhs.rows() != ROWS)
            throw std::invalid_argument("The input list has an invalid number of rows.");
        // Check the number of columns.
        if (rhs.cols() != COLS)
            throw std::invalid_argument("The input list has an invalid number of columns.");
        // Copy the data.
        for (std::size_t r = 0; r < ROWS; ++r)
            for (std::size_t c = 0; c < COLS; ++c)
                this->at(r, c) = rhs(r, c);
    }

    /// @brief Destroy the MatrixArray object.
    virtual ~MatrixArray() noexcept override
    {
        // Nothing to do.
    }

    /// @brief Returns a pointer to the internal data.
    /// @return the pointer.
    constexpr inline T *data() noexcept override
    {
        return _data.data();
    }

    /// @brief Returns a constant pointer to the internal data.
    /// @return the constant pointer.
    constexpr inline const T *data() const noexcept override
    {
        return _data.data();
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    constexpr auto &operator=(const MatrixBase<T> &rhs)
    {
        if (&rhs == this)
            return *this;
        // Check the number of rows.
        if (rhs.rows() != ROWS)
            throw std::invalid_argument("The input list has an invalid number of rows.");
        // Check the number of columns.
        if (rhs.cols() != COLS)
            throw std::invalid_argument("The input list has an invalid number of columns.");
        // Copy the data.
        for (std::size_t r = 0; r < ROWS; ++r)
            for (std::size_t c = 0; c < COLS; ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    template <typename T2>
    constexpr auto &operator=(const MatrixBase<T2> &rhs)
    {
        // Check the number of rows.
        if (rhs.rows() != ROWS)
            throw std::invalid_argument("The input list has an invalid number of rows.");
        // Check the number of columns.
        if (rhs.cols() != COLS)
            throw std::invalid_argument("The input list has an invalid number of columns.");
        // Copy the data.
        for (std::size_t r = 0; r < ROWS; ++r)
            for (std::size_t c = 0; c < COLS; ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Construct a new MatrixArray object.
    /// @param data
    constexpr auto &operator=(const std::initializer_list<std::initializer_list<T>> &data)
    {
        // Check the number of rows.
        if (data.size() != ROWS)
            throw std::invalid_argument("The input list has an invalid number of rows.");
        // Check the number of columns.
        if (data.begin()->size() != COLS)
            throw std::invalid_argument("The input list has an invalid number of columns.");
        // Get the number of rows.
        this->_rows = data.size();
        // Get the number of columns.
        this->_cols = data.begin()->size();
        // Get an interator for the data.
        auto r_it = data.begin();
        for (std::size_t r = 0, c = 0; r < this->rows(); ++r, ++r_it) {
            auto c_it = (*r_it).begin();
            for (c = 0; c < this->cols(); ++c, ++c_it) {
                this->at(r, c) = *c_it;
            }
        }
        return *this;
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    constexpr inline T &operator[](std::size_t pos) noexcept override
    {
        return _data[pos];
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    constexpr inline const T &operator[](std::size_t pos) const noexcept override
    {
        return _data[pos];
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the reference to the accessed item.
    constexpr inline T &operator()(std::size_t row, std::size_t col) noexcept override
    {
#ifdef ROW_MAJOR
        return _data[row * COLS + col];
#else
        return _data[col * ROWS + row];
#endif
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the const reference to the accessed item.
    constexpr inline const T &operator()(std::size_t row, std::size_t col) const noexcept override
    {
#ifdef ROW_MAJOR
        return _data[row * COLS + col];
#else
        return _data[col * ROWS + row];
#endif
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the reference to the accessed item.
    constexpr inline T &at(std::size_t row, std::size_t col) noexcept override
    {
#ifdef ROW_MAJOR
        return _data[row * COLS + col];
#else
        return _data[col * ROWS + row];
#endif
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the const reference to the accessed item.
    constexpr inline const T &at(std::size_t row, std::size_t col) const noexcept override
    {
#ifdef ROW_MAJOR
        return _data[row * COLS + col];
#else
        return _data[col * ROWS + row];
#endif
    }

    /// @brief A constant iterator poiting to the beginning of the internal data.
    /// @return the iterator.
    constexpr inline auto begin() const noexcept
    {
        return _data.begin();
    }

    /// @brief Iterator poiting to the beginning of the internal data.
    /// @return the iterator.
    constexpr inline auto begin() noexcept
    {
        return _data.begin();
    }

    /// @brief A constant iterator poiting to the end of the internal data.
    /// @return the iterator.
    constexpr inline auto end() const noexcept
    {
        return _data.end();
    }

    /// @brief Iterator poiting to the end of the internal data.
    /// @return the iterator.
    constexpr inline auto end() noexcept
    {
        return _data.end();
    }
};

} // namespace malg
