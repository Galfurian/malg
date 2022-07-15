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

/// @brief Matrix structure.
template <typename T>
class Matrix : public MatrixBase<T> {
public:
    /// The internal data.
    Vector<T> _data;

public:
    /// The data types of the element of the matrix.
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
    constexpr Matrix(std::size_t rows, std::size_t cols, const T &initial = T(0)) noexcept
        : MatrixBase<T>(rows, cols),
          _data(rows * cols, initial)
    {
        // Nothing to do.
    }

    /// @brief Construct a new Matrix object.
    /// @param rows
    /// @param cols
    /// @param data
    constexpr Matrix(std::size_t rows, std::size_t cols, const std::initializer_list<T> &data) noexcept
        : MatrixBase<T>(rows, cols),
          _data(rows * cols)
    {
        auto it = data.begin();
        for (std::size_t r = 0, c; r < this->rows(); ++r)
            for (c = 0; c < this->cols(); ++c, ++it)
                this->at(r, c) = *it;
    }

    /// @brief Construct a new Matrix object.
    /// @param rows
    /// @param cols
    /// @param data
    constexpr Matrix(std::size_t rows, std::size_t cols, const std::vector<T> &data) noexcept
        : MatrixBase<T>(rows, cols),
          _data(rows * cols)
    {
        auto it = data.begin();
        for (std::size_t r = 0, c; r < this->rows(); ++r)
            for (c = 0; c < this->cols(); ++c, ++it)
                this->at(r, c) = *it;
    }

    /// @brief Construct a new Matrix object.
    /// @param data
    constexpr Matrix(const std::initializer_list<std::initializer_list<T>> &data)
        : MatrixBase<T>(),
          _data()
    {
        // Check the number of rows.
        if (data.size() <= 0)
            throw std::invalid_argument("The input list has no rows.");
        // Check the number of columns.
        if (data.begin()->size() <= 0)
            throw std::invalid_argument("The input list has no columns.");
        // Get the number of rows.
        this->_rows = data.size();
        // Get the number of columns.
        this->_cols = data.begin()->size();
        // Initialize the vector.
        _data = Vector<T>(this->rows() * this->cols());
        // Get an interator for the data.
        auto r_it = data.begin();
        for (std::size_t r = 0, c = 0; r < this->rows(); ++r, ++r_it) {
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
        for (std::size_t i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    template <typename T2>
    constexpr Matrix(const MatrixBase<T2> &rhs) noexcept
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (std::size_t i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs
    constexpr Matrix(const Matrix<T> &rhs) noexcept
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (std::size_t i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    template <typename T2>
    constexpr Matrix(const Matrix<T2> &rhs) noexcept
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (std::size_t i = 0; i < (this->rows() * this->cols()); ++i)
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

    /// @brief Changes the shape of the matrix, but not the overall number of elements.
    /// @param rows the new number of rows.
    /// @param cols the new number of columns.
    /// @return a reference to this matrix.
    constexpr auto &reshape(std::size_t rows, std::size_t cols)
    {
        if (this->size() != (rows * cols))
            throw std::invalid_argument("The new shape has a different number of elements.");
        this->_rows = rows;
        this->_cols = cols;
        return *this;
    }

    /// @brief Changes the shape of the matrix, without preserving the overall number of elements.
    /// @param rows the new number of rows.
    /// @param cols the new number of columns.
    /// @return a reference to this matrix.
    constexpr auto &resize(std::size_t rows, std::size_t cols) noexcept
    {
        if ((this->rows() == rows) && (this->cols() == cols))
            return *this;
#ifdef ROW_MAJOR
        if (this->cols() == cols) {
#else
        if (this->rows() == rows) {
#endif
            // Resize the data.
            _data.resize(rows * cols);
        } else {
            malg::Matrix<T> m(rows, cols, static_cast<T>(0));
            for (std::size_t r = 0; r < std::min(this->rows(), rows); ++r)
                for (std::size_t c = 0; c < std::min(this->cols(), cols); ++c)
                    m(r, c) = this->at(r, c);
            // Move the data.
            _data = std::move(m._data);
        }
        // Set the new rows and columns.
        this->_rows = rows;
        this->_cols = cols;
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
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
            for (std::size_t r = 0; r < this->rows(); ++r)
                for (std::size_t c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
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
            for (std::size_t r = 0; r < this->rows(); ++r)
                for (std::size_t c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
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
            for (std::size_t r = 0; r < this->rows(); ++r)
                for (std::size_t c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
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
            for (std::size_t r = 0; r < this->rows(); ++r)
                for (std::size_t c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
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
        // Check the number of rows.
        if (data.size() <= 0)
            return *this;
        // Check the number of columns.
        if (data.begin()->size() <= 0)
            return *this;
        // Get the number of rows.
        this->_rows = data.size();
        // Get the number of columns.
        this->_cols = data.begin()->size();
        // Initialize the vector.
        _data = Vector<T>(this->rows() * this->cols());
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
        return this->at(row, col);
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the const reference to the accessed item.
    constexpr inline const T &operator()(std::size_t row, std::size_t col) const noexcept override
    {
        return this->at(row, col);
    }

    /// @brief Operator for generating a view of the matrix.
    /// @param row the range of rows to extract.
    /// @param col the range of columns to extract.
    /// @returns the desired view.
    constexpr View<Matrix<T>> operator()(Range row, Range col) noexcept;

    /// @brief Operator for generating a view of the matrix.
    /// @param row the range of rows to extract.
    /// @param col the range of columns to extract.
    /// @returns the desired view.
    constexpr View<const Matrix<T>> operator()(Range row, Range col) const noexcept;

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the reference to the accessed item.
    constexpr inline T &at(std::size_t row, std::size_t col) noexcept override
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
    /// @returns the const reference to the accessed item.
    constexpr inline const T &at(std::size_t row, std::size_t col) const noexcept override
    {
#ifdef ROW_MAJOR
        return _data[(row * this->cols()) + col];
#else
        return _data[(col * this->rows()) + row];
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
