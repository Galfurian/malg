/// @file matrix.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief The matrix class.

#pragma once

#include "malg/type_traits.hpp"
#include "malg/matrix_base.hpp"
#include "malg/vector.hpp"

#include <algorithm>
#include <initializer_list>
#include <vector>

namespace malg
{

struct Range;

template <typename MatrixType>
class View;

/// @brief Matrix structure.
template <typename T>
class Matrix : public MatrixBase<T> {
private:
    /// The internal data.
    Vector<T> _data;
    /// Rows of the matrix.
    using MatrixBase<T>::_rows;
    /// Columns of the matrix.
    using MatrixBase<T>::_cols;

public:
    /// The data types of the element of the matrix.
    using value_type = T;
    /// The data types for iterating the element of the vector.
    using iterator = typename MatrixBase<T>::iterator;
    /// The data types for iterating the element of the vector.
    using const_iterator = typename MatrixBase<T>::const_iterator;

    /// @brief Construct a new Matrix object.
    Matrix()
        : MatrixBase<T>(),
          _data()
    {
        // Nothing to do.
    }

    /// @brief Construct a new Matrix object.
    /// @param rows
    /// @param cols
    /// @param initial
    Matrix(std::size_t rows, std::size_t cols, const T &initial = T(0))
        : MatrixBase<T>(rows, cols),
          _data(rows * cols, initial)
    {
        // Nothing to do.
    }

    /// @brief Construct a new Matrix object.
    /// @param rows
    /// @param cols
    /// @param data
    Matrix(std::size_t rows, std::size_t cols, const std::initializer_list<T> &data)
        : MatrixBase<T>(rows, cols),
          _data(rows * cols)
    {
        auto it = data.begin();
        for (std::size_t r = 0, c = 0; r < _rows; ++r) {
            for (c = 0; c < _cols; ++c, ++it) {
                this->at(r, c) = *it;
            }
        }
    }

    /// @brief Construct a new Matrix object.
    /// @param rows
    /// @param cols
    /// @param data
    Matrix(std::size_t rows, std::size_t cols, const std::vector<T> &data)
        : MatrixBase<T>(rows, cols),
          _data(rows * cols)
    {
        auto it = data.begin();
        for (std::size_t r = 0, c = 0; r < _rows; ++r) {
            for (c = 0; c < _cols; ++c, ++it) {
                this->at(r, c) = *it;
            }
        }
    }

    /// @brief Construct a new Matrix object.
    /// @param data
    Matrix(const std::initializer_list<std::initializer_list<T>> &data)
        : MatrixBase<T>(),
          _data()
    {
        // Check the number of rows.
        if (data.size() <= 0) {
            throw std::invalid_argument("The input list has no rows.");
        }
        // Check the number of columns.
        if (data.begin()->size() <= 0) {
            throw std::invalid_argument("The input list has no columns.");
        }
        // Get the number of rows.
        this->_rows = data.size();
        // Get the number of columns.
        this->_cols = data.begin()->size();
        // Initialize the vector.
        _data = Vector<T>(_rows * _cols);
        // Get an interator for the data.
        auto r_it = data.begin();
        for (std::size_t r = 0, c = 0; r < _rows; ++r, ++r_it) {
            auto c_it = (*r_it).begin();
            for (c = 0; c < _cols; ++c, ++c_it) {
                this->at(r, c) = *c_it;
            }
        }
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs
    explicit Matrix(const MatrixBase<T> &rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (std::size_t i = 0; i < (_rows * _cols); ++i) {
            _data[i] = rhs[i];
        }
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs the other matrix.
    template <typename T2>
    explicit Matrix(const MatrixBase<T2> &rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (std::size_t i = 0; i < (_rows * _cols); ++i) {
            _data[i] = rhs[i];
        }
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs the other matrix.
    Matrix(const Matrix<T> &rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (std::size_t i = 0; i < (_rows * _cols); ++i) {
            _data[i] = rhs[i];
        }
    }

    /// @brief Construct a new Matrix object.
    /// @tparam T2 the type of the other matrix.
    /// @param rhs the other matrix.
    template <typename T2>
    explicit Matrix(const Matrix<T2> &rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (std::size_t i = 0; i < (_rows * _cols); ++i) {
            _data[i] = rhs[i];
        }
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs
    Matrix(Matrix<T> &&rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(std::move(rhs._data))
    {
        // Nothing to do.
    }

    /// @brief Construct a new Matrix object.
    /// @tparam T2 the type of the other matrix.
    /// @param rhs the other matrix.
    template <typename T2>
    explicit Matrix(Matrix<T2> &&rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(std::move(rhs._data))
    {
        // Nothing to do.
    }

    /// @brief Destroy the Matrix object.
    virtual ~Matrix() override
    {
        // Nothing to do.
    }

    /// @brief Returns a pointer to the internal data.
    /// @return the pointer.
    inline iterator data() override
    {
        return _data.data();
    }

    /// @brief Returns a constant pointer to the internal data.
    /// @return the constant pointer.
    inline const_iterator data() const override
    {
        return _data.data();
    }

    /// @brief Changes the shape of the matrix, but not the overall number of elements.
    /// @param rows the new number of rows.
    /// @param cols the new number of columns.
    /// @return a reference to this matrix.
    auto &reshape(std::size_t rows, std::size_t cols)
    {
        if (this->size() != (rows * cols)) {
            throw std::invalid_argument("The new shape has a different number of elements.");
        }
        this->_rows = rows;
        this->_cols = cols;
        return *this;
    }

    /// @brief Changes the shape of the matrix, without preserving the overall number of elements.
    /// @param rows the new number of rows.
    /// @param cols the new number of columns.
    /// @return a reference to this matrix.
    auto &resize(std::size_t rows, std::size_t cols)
    {
        if ((_rows == rows) && (_cols == cols)) {
            return *this;
        }
#ifdef ROW_MAJOR
        if (_cols == cols) {
#else
        if (_rows == rows) {
#endif
            // Resize the data.
            _data.resize(rows * cols);
        } else {
            malg::Matrix<T> m(rows, cols, static_cast<T>(0));
            for (std::size_t r = 0; r < std::min(_rows, rows); ++r) {
                for (std::size_t c = 0; c < std::min(_cols, cols); ++c) {
                    m(r, c) = this->at(r, c);
                }
            }
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
    auto &operator=(const Matrix<T> &rhs)
    {
        if (&rhs == this) {
            return *this;
        }
        if (rhs.size() > 0) {
            // Initialize the new vector of data.
            _data = Vector<T>(rhs.rows() * rhs.cols());
            // Set the new rows and columns.
            this->_rows = rhs.rows();
            this->_cols = rhs.cols();
            // Copy the content.
            for (std::size_t r = 0; r < _rows; ++r) {
                for (std::size_t c = 0; c < _cols; ++c) {
                    this->at(r, c) = rhs(r, c);
                }
            }
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    auto &operator=(const MatrixBase<T> &rhs)
    {
        if (&rhs == this) {
            return *this;
        }
        if (rhs.size() > 0) {
            // Initialize the new vector of data.
            _data = Vector<T>(rhs.rows() * rhs.cols());
            // Set the new rows and columns.
            this->_rows = rhs.rows();
            this->_cols = rhs.cols();
            // Copy the content.
            for (std::size_t r = 0; r < _rows; ++r) {
                for (std::size_t c = 0; c < _cols; ++c) {
                    this->at(r, c) = rhs(r, c);
                }
            }
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    template <typename T2>
    auto &operator=(const MatrixBase<T2> &rhs)
    {
        if (rhs.size() > 0) {
            // Initialize the new vector of data.
            _data = Vector<T>(rhs.rows() * rhs.cols());
            // Set the new rows and columns.
            this->_rows = rhs.rows();
            this->_cols = rhs.cols();
            // Copy the content.
            for (std::size_t r = 0; r < _rows; ++r) {
                for (std::size_t c = 0; c < _cols; ++c) {
                    this->at(r, c) = rhs(r, c);
                }
            }
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    template <typename T2>
    Matrix<T> &operator=(const Matrix<T2> &rhs)
    {
        if (&rhs != this) {
            if (!std::is_same_v<T, T2>) {
                if (rhs.size() > 0) {
                    // Initialize the new vector of data.
                    _data = Vector<T>(rhs.rows() * rhs.cols());
                    // Set the new rows and columns.
                    this->_rows = rhs.rows();
                    this->_cols = rhs.cols();
                    // Copy the content.
                    for (std::size_t r = 0; r < _rows; ++r) {
                        for (std::size_t c = 0; c < _cols; ++c) {
                            this->at(r, c) = rhs(r, c);
                        }
                    }
                }
            }
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    auto &operator=(Matrix<T> &&rhs)
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

    /// @brief Assign operator.
    /// @param data the data to use.
    /// @return a reference to this.
    auto &operator=(const std::initializer_list<std::initializer_list<T>> &data)
    {
        // Check the number of rows.
        if (data.size() <= 0) {
            return *this;
        }
        // Check the number of columns.
        if (data.begin()->size() <= 0) {
            return *this;
        }
        // Get the number of rows.
        this->_rows = data.size();
        // Get the number of columns.
        this->_cols = data.begin()->size();
        // Initialize the vector.
        _data = Vector<T>(_rows * _cols);
        // Get an interator for the data.
        auto r_it = data.begin();
        for (std::size_t r = 0, c = 0; r < _rows; ++r, ++r_it) {
            auto c_it = (*r_it).begin();
            for (c = 0; c < _cols; ++c, ++c_it) {
                this->at(r, c) = *c_it;
            }
        }
        return *this;
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    inline T &operator[](std::size_t pos) override
    {
        return _data[pos];
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    inline const T &operator[](std::size_t pos) const override
    {
        return _data[pos];
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the reference to the accessed item.
    inline T &operator()(std::size_t row, std::size_t col) override
    {
        return this->at(row, col);
    }

    /// @brief Operator for accessing the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the const reference to the accessed item.
    inline const T &operator()(std::size_t row, std::size_t col) const override
    {
        return this->at(row, col);
    }

    /// @brief Operator for generating a view of the matrix.
    /// @param row the range of rows to extract.
    /// @param col the range of columns to extract.
    /// @returns the desired view.
    View<Matrix<T>> operator()(Range row, Range col);

    /// @brief Operator for generating a view of the matrix.
    /// @param row the range of rows to extract.
    /// @param col the range of columns to extract.
    /// @returns the desired view.
    View<const Matrix<T>> operator()(Range row, Range col) const;

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the reference to the accessed item.
    inline T &at(std::size_t row, std::size_t col) override
    {
#ifdef ROW_MAJOR
        return _data[(row * _cols) + col];
#else
        return _data[(col * _rows) + row];
#endif
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the const reference to the accessed item.
    inline const T &at(std::size_t row, std::size_t col) const override
    {
#ifdef ROW_MAJOR
        return _data[(row * _cols) + col];
#else
        return _data[(col * _rows) + row];
#endif
    }

    /// @brief A constant iterator poiting to the beginning of the internal data.
    /// @return the iterator.
    inline auto begin() const
    {
        return _data.begin();
    }

    /// @brief Iterator poiting to the beginning of the internal data.
    /// @return the iterator.
    inline auto begin()
    {
        return _data.begin();
    }

    /// @brief A constant iterator poiting to the end of the internal data.
    /// @return the iterator.
    inline auto end() const
    {
        return _data.end();
    }

    /// @brief Iterator poiting to the end of the internal data.
    /// @return the iterator.
    inline auto end()
    {
        return _data.end();
    }
};

} // namespace malg
