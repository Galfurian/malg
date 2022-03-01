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

template <typename T>
class Matrix : public MatrixBase<T> {
public:
    Vector<T> _data;

public:
    using this_type = Matrix<T>;
    using base_type = T;

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
    Matrix(unsigned rows, unsigned cols, const T &initial)
        : MatrixBase<T>(rows, cols),
          _data(rows * cols, initial)
    {
        // Nothing to do.
    }

    /// @brief Construct a new Matrix object.
    /// @param rows
    /// @param cols
    /// @param data
    Matrix(unsigned rows, unsigned cols, const std::initializer_list<T> &data)
        : MatrixBase<T>(rows, cols),
          _data(rows * cols)
    {
        auto it = data.begin();
        for (unsigned r = 0, c; r < this->rows(); ++r)
            for (c = 0; c < this->cols(); ++c, ++it)
                this->at(r, c) = *it;
    }

    /// @brief Construct a new Matrix object.
    /// @param data
    Matrix(const std::initializer_list<std::initializer_list<T>> &data)
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
        for (unsigned r = 0, c = 0; r < this->rows(); ++r, ++r_it) {
            auto c_it = (*r_it).begin();
            for (c = 0; c < this->cols(); ++c, ++c_it) {
                this->at(r, c) = *c_it;
            }
        }
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs
    Matrix(const MatrixBase<T> &rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (unsigned i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    template <typename T2>
    Matrix(const MatrixBase<T2> &rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (unsigned i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    /// @param rhs
    Matrix(const Matrix<T> &rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (unsigned i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
    }

    /// @brief Construct a new Matrix object.
    template <typename T2>
    Matrix(const Matrix<T2> &rhs)
        : MatrixBase<T>(rhs.rows(), rhs.cols()),
          _data(rhs.rows() * rhs.cols())
    {
        for (unsigned i = 0; i < (this->rows() * this->cols()); ++i)
            _data[i] = rhs[i];
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
    template <typename T2>
    Matrix(Matrix<T2> &&rhs)
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

    inline T *data() override
    {
        return _data.data();
    }

    inline const T *data() const override
    {
        return _data.data();
    }

    auto &reshape(unsigned rows, unsigned cols)
    {
        assert(this->size() == (rows * cols));
        this->_rows = rows;
        this->_cols = cols;
        return *this;
    }

    auto &resize(unsigned rows, unsigned cols)
    {
        malg::Matrix<T> m(rows, cols, static_cast<T>(0));
        for (unsigned r = 0; r < std::min(this->rows(), rows); ++r)
            for (unsigned c = 0; c < std::min(this->cols(), cols); ++c)
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
    auto &operator=(const Matrix<T> &rhs)
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
            for (unsigned r = 0; r < this->rows(); ++r)
                for (unsigned c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    auto &operator=(const MatrixBase<T> &rhs)
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
            for (unsigned r = 0; r < this->rows(); ++r)
                for (unsigned c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
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
            for (unsigned r = 0; r < this->rows(); ++r)
                for (unsigned c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
    template <typename T2>
    auto &operator=(const Matrix<T2> &rhs)
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
            for (unsigned r = 0; r < this->rows(); ++r)
                for (unsigned c = 0; c < this->cols(); ++c)
                    this->at(r, c) = rhs(r, c);
        }
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @return a reference to this matrix.
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

    /// @brief Construct a new Matrix object.
    /// @param data
    auto &operator=(const std::initializer_list<std::initializer_list<T>> &data)
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
        for (unsigned r = 0, c = 0; r < this->rows(); ++r, ++r_it) {
            auto c_it = (*r_it).begin();
            for (c = 0; c < this->cols(); ++c, ++c_it) {
                this->at(r, c) = *c_it;
            }
        }
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    inline T &operator[](unsigned pos) override
    {
        return _data[pos];
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @return the reference to the accessed item.
    inline const T &operator[](unsigned pos) const override
    {
        return _data[pos];
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
    inline const T &at(unsigned row, unsigned col) const override
    {
#ifdef ROW_MAJOR
        return _data[(row * this->cols()) + col];
#else
        return _data[(col * this->rows()) + row];
#endif
    }

    inline auto begin() const
    {
        return _data.begin();
    }

    inline auto begin()
    {
        return _data.begin();
    }

    inline auto end() const
    {
        return _data.end();
    }

    inline auto end()
    {
        return _data.end();
    }
};

} // namespace malg
