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
    unsigned _start_row, _end_row;
    unsigned _start_col, _end_col;

public:
    View(MatrixType &matrix,
         unsigned start_row = 0,
         unsigned end_row   = -1,
         unsigned start_col = 0,
         unsigned end_col   = -1)
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
    inline unsigned rows() const override
    {
        return std::min(_end_row, _matrix.rows()) - _start_row;
    }

    /// @brief Get the number of columns of the matrix.
    /// @return the number of columns.
    inline unsigned cols() const override
    {
        return std::min(_end_col, _matrix.cols()) - _start_col;
    }

    /// @brief Returns the total size of the matrix.
    /// @return the total size of the matrix.
    inline unsigned size() const override
    {
        return this->rows() * this->cols();
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    template <typename T2>
    auto &operator=(const MatrixBase<T2> &rhs)
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (unsigned r = 0; r < rhs.rows(); ++r)
            for (unsigned c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    auto &operator=(const View<MatrixType> &rhs)
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (unsigned r = 0; r < rhs.rows(); ++r)
            for (unsigned c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    template <typename T2>
    auto &operator=(const View<T2> &rhs)
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (unsigned r = 0; r < rhs.rows(); ++r)
            for (unsigned c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    inline T *data() override
    {
        return _matrix.data();
    }

    inline const T *data() const override
    {
        return _matrix.data();
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    inline T &operator[](unsigned pos) override
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
    inline const T &operator[](unsigned pos) const override
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
    inline T &operator()(unsigned row, unsigned col) override
    {
        return this->at(row, col);
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    inline const T &operator()(unsigned row, unsigned col) const override
    {
        return this->at(row, col);
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    inline T &at(unsigned row, unsigned col) override
    {
        return _matrix.at(_start_row + row, _start_col + col);
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    inline const T &at(unsigned row, unsigned col) const override
    {
        return _matrix.at(_start_row + row, _start_col + col);
    }
};

template <typename MatrixType>
class CofactorView : public MatrixBase<typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::base_type, typename MatrixType::base_type>> {
public:
    using T = typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::base_type, typename MatrixType::base_type>;

    MatrixType &_matrix;
    unsigned _excluded_row;
    unsigned _excluded_col;

public:
    CofactorView(MatrixType &matrix,
                 unsigned excluded_row,
                 unsigned excluded_col)
        : _matrix(matrix),
          _excluded_row(excluded_row),
          _excluded_col(excluded_col)
    {
        // Nothing to do.
    }

    /// @brief Get the number of rows of the matrix.
    /// @return the number of rows.
    inline unsigned rows() const override
    {
        return _matrix.rows() - 1;
    }

    /// @brief Get the number of columns of the matrix.
    /// @return the number of columns.
    inline unsigned cols() const override
    {
        return _matrix.cols() - 1;
    }

    /// @brief Returns the total size of the matrix.
    /// @return the total size of the matrix.
    inline unsigned size() const override
    {
        return this->rows() * this->cols();
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    template <typename T2>
    auto &operator=(const MatrixBase<T2> &rhs)
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (unsigned r = 0; r < rhs.rows(); ++r)
            for (unsigned c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    inline T *data() override
    {
        return _matrix.data();
    }

    inline const T *data() const override
    {
        return _matrix.data();
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    inline T &operator[](unsigned pos) override
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
    inline const T &operator[](unsigned pos) const override
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
    inline T &operator()(unsigned row, unsigned col) override
    {
        return this->at(row, col);
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    inline const T &operator()(unsigned row, unsigned col) const override
    {
        return this->at(row, col);
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    inline T &at(unsigned row, unsigned col) override
    {
        return _matrix.at(row + (row >= _excluded_row), col + (col >= _excluded_col));
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    inline const T &at(unsigned row, unsigned col) const override
    {
        return _matrix.at(row + (row >= _excluded_row), col + (col >= _excluded_col));
    }
};

} // namespace malg
