/// @file matrix2d.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief The matrix class.

#pragma once

#include "malg/type_traits.hpp"
#include "malg/matrix_base.hpp"
#include "malg/vector.hpp"

#include <initializer_list>
#include <cassert>

namespace malg
{

class Range;

template <typename MatrixType>
class View;

template <typename T>
class Matrix : public MatrixBase<T> {
public:
    Vector<T> _data;

public:
    using value_type = T;

    /// @brief Construct a new Matrix object.
    constexpr Matrix() noexcept
        : MatrixBase<T>(),
          _data()
    {
        // Nothing to do.
    }

    /// @brief Construct a new Matrix object.
    /// @param rows
    /// @param cols
    /// @param initial
    constexpr Matrix(size_type_t rows, size_type_t cols, const T &initial = T(0)) noexcept
        : MatrixBase<T>(rows, cols),
          _data(rows * cols, initial)
    {
        // Nothing to do.
    }

    /// @brief Construct a new Matrix object.
    /// @param rows
    /// @param cols
    /// @param data
    constexpr Matrix(size_type_t rows, size_type_t cols, const std::initializer_list<T> &data) noexcept
        : MatrixBase<T>(rows, cols),
          _data(rows * cols)
    {
        auto it = data.begin();
        for (size_type_t r = 0, c; r < this->rows(); ++r)
            for (c = 0; c < this->cols(); ++c, ++it)
                this->at(r, c) = *it;
    }

    /// @brief Construct a new Matrix object.
    /// @param data
    constexpr Matrix(const std::initializer_list<std::initializer_list<T>> &data) noexcept
        : MatrixBase<T>(),
          _data()
    {
        // Get the number of rows.
        this->_rows = data.size();
        assert(this->rows() > 0);
        // Get the number of columns.
        this->_cols = (*data.begin()).size();
        assert(this->cols() > 0);
        // Initialize the vector.
        _data = Vector<T>(this->rows() * this->cols());
        // Get an interator for the data.
        auto r_it = data.begin();
        for (size_type_t r = 0, c = 0; r < this->rows(); ++r, ++r_it) {
            auto c_it = (*r_it).begin();
            for (c = 0; c < this->cols(); ++c, ++c_it) {
                this->at(r, c) = *c_it;
            }
        }
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs
    constexpr Matrix(const MatrixBase<T> &rhs) noexcept
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (size_type_t i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    template <typename T2>
    constexpr Matrix(const MatrixBase<T2> &rhs) noexcept
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (size_type_t i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs
    constexpr Matrix(const Matrix<T> &rhs) noexcept
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (size_type_t i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    template <typename T2>
    constexpr Matrix(const Matrix<T2> &rhs) noexcept
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (size_type_t i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs
    constexpr Matrix(Matrix<T> &&rhs) noexcept
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(std::move(rhs._data))
    {
        // Nothing to do.
    }

    /// @brief Construct a new Matrix object.
    template <typename T2>
    constexpr Matrix(Matrix<T2> &&rhs) noexcept
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(std::move(rhs._data))
    {
        // Nothing to do.
    }

    /// @brief Destroy the Matrix object.
    virtual ~Matrix() noexcept override
    {
        // Nothing to do.
    }

    constexpr inline T *data() noexcept override
    {
        return _data.data();
    }

    constexpr inline const T *data() const noexcept override
    {
        return _data.data();
    }

    constexpr auto &reshape(size_type_t rows, size_type_t cols) noexcept
    {
        assert(this->size() == (rows * cols));
        this->_rows = rows;
        this->_cols = cols;
        return *this;
    }

    constexpr auto &resize(size_type_t rows, size_type_t cols) noexcept
    {
        malg::Matrix<T> m(rows, cols, static_cast<T>(0));
        for (size_type_t r = 0; r < std::min(this->rows(), rows); ++r)
            for (size_type_t c = 0; c < std::min(this->cols(), cols); ++c)
                m(r, c) = this->at(r, c);
        // Move the data.
        _data = std::move(m._data);
        // Set the new rows and columns.
        this->_rows = rows;
        this->_cols = cols;
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    constexpr auto &operator=(const Matrix<T> &rhs) noexcept
    {
        if (&rhs == this)
            return *this;
        if (rhs.size() > 0) {
            // Initialize the new vector of data.
            _data = Vector<T>(rhs.rows() * rhs.cols());
            // Set the new rows and columns.
            this->_rows = rhs.rows();
            this->_cols = rhs.cols();
            // Copy the content.
            for (size_type_t r = 0; r < this->rows(); ++r)
                for (size_type_t c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    constexpr auto &operator=(const MatrixBase<T> &rhs) noexcept
    {
        if (&rhs == this)
            return *this;
        if (rhs.size() > 0) {
            // Initialize the new vector of data.
            _data = Vector<T>(rhs.rows() * rhs.cols());
            // Set the new rows and columns.
            this->_rows = rhs.rows();
            this->_cols = rhs.cols();
            // Copy the content.
            for (size_type_t r = 0; r < this->rows(); ++r)
                for (size_type_t c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    template <typename T2>
    constexpr auto &operator=(const MatrixBase<T2> &rhs) noexcept
    {
        if (rhs.size() > 0) {
            // Initialize the new vector of data.
            _data = Vector<T>(rhs.rows() * rhs.cols());
            // Set the new rows and columns.
            this->_rows = rhs.rows();
            this->_cols = rhs.cols();
            // Copy the content.
            for (size_type_t r = 0; r < this->rows(); ++r)
                for (size_type_t c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    template <typename T2>
    constexpr auto &operator=(const Matrix<T2> &rhs) noexcept
    {
        if constexpr (std::is_same_v<T, T2>)
            if (&rhs == this)
                return *this;
        if (rhs.size() > 0) {
            // Initialize the new vector of data.
            _data = Vector<T>(rhs.rows() * rhs.cols());
            // Set the new rows and columns.
            this->_rows = rhs.rows();
            this->_cols = rhs.cols();
            // Copy the content.
            for (size_type_t r = 0; r < this->rows(); ++r)
                for (size_type_t c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    constexpr auto &operator=(Matrix<T> &&rhs) noexcept
    {
        if ((&rhs != this) && (rhs.size() > 0)) {
            // Copy the content.
            _data = std::move(rhs._data);
            // Set the new rows and columns.
            this->_rows = rhs.rows();
            this->_cols = rhs.cols();
        }
        return *this;
    }

    /// @brief Construct a new Matrix object.
    /// @param data
    constexpr auto &operator=(const std::initializer_list<std::initializer_list<T>> &data) noexcept
    {
        // Get the number of rows.
        this->_rows = data.size();
        assert(this->rows() > 0);
        // Get the number of columns.
        this->_cols = (*data.begin()).size();
        assert(this->cols() > 0);
        // Initialize the vector.
        _data = Vector<T>(this->rows() * this->cols());
        // Get an interator for the data.
        auto r_it = data.begin();
        for (size_type_t r = 0, c = 0; r < this->rows(); ++r, ++r_it) {
            auto c_it = (*r_it).begin();
            for (c = 0; c < this->cols(); ++c, ++c_it) {
                this->at(r, c) = *c_it;
            }
        }
        return *this;
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    constexpr inline T &operator[](size_type_t pos) noexcept override
    {
        return _data[pos];
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    constexpr inline const T &operator[](size_type_t pos) const noexcept override
    {
        return _data[pos];
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

    /// @brief Operator for generating a view of the matrix.
    /// @param row the range of rows to extract.
    /// @param col the range of columns to extract.
    /// @return the desired view.
    constexpr View<Matrix<T>> operator()(Range row, Range col) noexcept;

    /// @brief Operator for generating a view of the matrix.
    /// @param row the range of rows to extract.
    /// @param col the range of columns to extract.
    /// @return the desired view.
    constexpr View<const Matrix<T>> operator()(Range row, Range col) const noexcept;

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the reference to the accessed item.
    constexpr inline T &at(size_type_t row, size_type_t col) noexcept override
    {
#ifdef ROW_MAJOR
        return _data[(row * this->cols()) + col];
#else
        return _data[(col * this->rows()) + row];
#endif
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @return the const reference to the accessed item.
    constexpr inline const T &at(size_type_t row, size_type_t col) const noexcept override
    {
#ifdef ROW_MAJOR
        return _data[(row * this->cols()) + col];
#else
        return _data[(col * this->rows()) + row];
#endif
    }

    constexpr inline auto begin() const noexcept
    {
        return _data.begin();
    }

    constexpr inline auto begin() noexcept
    {
        return _data.begin();
    }

    constexpr inline auto end() const noexcept
    {
        return _data.end();
    }

    constexpr inline auto end() noexcept
    {
        return _data.end();
    }
};

} // namespace malg
