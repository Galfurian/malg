/// @file view.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief View classes, which stores nothing but provide structured access to
/// real matrices.

#pragma once

#include "malg/matrix.hpp"

#include <iostream>

namespace malg
{

template <typename MatrixType>
class View : public MatrixBase<typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::base_type, typename MatrixType::base_type>> {
public:
    using T = typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::base_type, typename MatrixType::base_type>;

    MatrixType &_matrix;
    size_type_t _start_row, _end_row;
    size_type_t _start_col, _end_col;

public:
    constexpr View(MatrixType &matrix,
         size_type_t start_row = 0,
         size_type_t end_row   = -1,
         size_type_t start_col = 0,
         size_type_t end_col   = -1) noexcept
        : _matrix(matrix),
          _start_row(start_row),
          _end_row(end_row),
          _start_col(start_col),
          _end_col(end_col)
    {
        assert(start_row < end_row);
        assert(start_col < end_col);
    }

    /// @brief Get the number of rows of the matrix.
    /// @return the number of rows.
    constexpr inline size_type_t rows() const noexcept override
    {
        return std::min(_end_row, _matrix.rows()) - _start_row;
    }

    /// @brief Get the number of columns of the matrix.
    /// @return the number of columns.
    constexpr inline size_type_t cols() const noexcept override
    {
        return std::min(_end_col, _matrix.cols()) - _start_col;
    }

    /// @brief Returns the total size of the matrix.
    /// @return the total size of the matrix.
    constexpr inline size_type_t size() const noexcept override
    {
        return this->rows() * this->cols();
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    template <typename T2>
    constexpr inline auto &operator=(const MatrixBase<T2> &rhs) noexcept
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (size_type_t r = 0; r < rhs.rows(); ++r)
            for (size_type_t c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    constexpr inline auto &operator=(const View<MatrixType> &rhs) noexcept
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (size_type_t r = 0; r < rhs.rows(); ++r)
            for (size_type_t c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    template <typename T2>
    constexpr inline auto &operator=(const View<T2> &rhs) noexcept
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (size_type_t r = 0; r < rhs.rows(); ++r)
            for (size_type_t c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    constexpr inline T *data() noexcept override
    {
        return _matrix.data();
    }

    constexpr inline const T *data() const noexcept override
    {
        return _matrix.data();
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    constexpr inline T &operator[](size_type_t pos) noexcept override
    {
#ifdef ROW_MAJOR
        return this->at((pos / this->rows()), (pos % this->rows()));
#else
        return this->at((pos % this->rows()), (pos / this->rows()));
#endif
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    constexpr inline const T &operator[](size_type_t pos) const noexcept override
    {
#ifdef ROW_MAJOR
        return this->at((pos / this->rows()), (pos % this->rows()));
#else
        return this->at((pos % this->rows()), (pos / this->rows()));
#endif
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    constexpr inline T &operator()(size_type_t row, size_type_t col) noexcept override
    {
        return this->at(row, col);
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    constexpr inline const T &operator()(size_type_t row, size_type_t col) const noexcept override
    {
        return this->at(row, col);
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    constexpr inline T &at(size_type_t row, size_type_t col) noexcept override
    {
        return _matrix.at(_start_row + row, _start_col + col);
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    constexpr inline const T &at(size_type_t row, size_type_t col) const noexcept override
    {
        return _matrix.at(_start_row + row, _start_col + col);
    }
};

template <typename MatrixType>
class CofactorView : public MatrixBase<typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::base_type, typename MatrixType::base_type>> {
public:
    using T = typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::base_type, typename MatrixType::base_type>;

    MatrixType &_matrix;
    size_type_t _excluded_row;
    size_type_t _excluded_col;

public:
    constexpr CofactorView(MatrixType &matrix,
                 size_type_t excluded_row,
                 size_type_t excluded_col) noexcept
        : _matrix(matrix),
          _excluded_row(excluded_row),
          _excluded_col(excluded_col)
    {
        // Nothing to do.
    }

    /// @brief Get the number of rows of the matrix.
    /// @return the number of rows.
    constexpr inline size_type_t rows() const noexcept override
    {
        return _matrix.rows() - 1;
    }

    /// @brief Get the number of columns of the matrix.
    /// @return the number of columns.
    constexpr inline size_type_t cols() const noexcept override
    {
        return _matrix.cols() - 1;
    }

    /// @brief Returns the total size of the matrix.
    /// @return the total size of the matrix.
    constexpr inline size_type_t size() const noexcept override
    {
        return this->rows() * this->cols();
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    template <typename T2>
    constexpr inline auto &operator=(const MatrixBase<T2> &rhs) noexcept
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (size_type_t r = 0; r < rhs.rows(); ++r)
            for (size_type_t c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    constexpr inline T *data() noexcept override
    {
        return _matrix.data();
    }

    constexpr inline const T *data() const noexcept override
    {
        return _matrix.data();
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    constexpr inline T &operator[](size_type_t pos) noexcept override
    {
#ifdef ROW_MAJOR
        return this->at((pos / this->rows()), (pos % this->rows()));
#else
        return this->at((pos % this->cols()), (pos / this->cols()));
#endif
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    constexpr inline const T &operator[](size_type_t pos) const noexcept override
    {
#ifdef ROW_MAJOR
        return this->at((pos / this->rows()), (pos % this->rows()));
#else
        return this->at((pos % this->cols()), (pos / this->cols()));
#endif
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    constexpr inline T &operator()(size_type_t row, size_type_t col) noexcept override
    {
        return this->at(row, col);
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    constexpr inline const T &operator()(size_type_t row, size_type_t col) const noexcept override
    {
        return this->at(row, col);
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    constexpr inline T &at(size_type_t row, size_type_t col) noexcept override
    {
        return _matrix.at(row + (row >= _excluded_row), col + (col >= _excluded_col));
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    constexpr inline const T &at(size_type_t row, size_type_t col) const noexcept override
    {
        return _matrix.at(row + (row >= _excluded_row), col + (col >= _excluded_col));
    }
};

} // namespace malg
