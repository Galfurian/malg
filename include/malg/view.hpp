/// @file view.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief View classes, which stores nothing but provide structured access to
/// real matrices.

#pragma once

#include "malg/matrix.hpp"

#include <iostream>

namespace malg
{

/// @brief Pairs of indices that identify a range.
struct Range {
    /// The starting index.
    std::size_t start;
    /// The final index.
    std::size_t stop;

    /// @brief Returns the [0, inf] index.
    /// @return the built range.
    static inline Range all()
    {
        return { 0, static_cast<std::size_t>(-1) };
    }
};

/// @brief A view of a matrix-like structure.
template <typename MatrixType>
class View : public MatrixBase<typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::value_type, typename MatrixType::value_type>> {
public:
    /// The datatype of the referenced matrix.
    using T = typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::value_type, typename MatrixType::value_type>;
    /// Pointer to the matrix.
    MatrixType *_matrix;
    /// The range of rows.
    Range _row;
    /// The range of columns.
    Range _col;

public:
    /// @brief Creates a **View** for the given **matrix**.
    /// @param matrix the matrix.
    /// @param start_row the starating row.
    /// @param end_row the final row.
    /// @param start_col the starating column.
    /// @param end_col the final column.
    View(MatrixType *matrix, std::size_t start_row = 0, std::size_t end_row = std::numeric_limits<std::size_t>::max(), std::size_t start_col = 0, std::size_t end_col = std::numeric_limits<std::size_t>::max())
        : _matrix(matrix),
          _row{ start_row, end_row },
          _col{ start_col, end_col }
    {
        assert(_row.start < _row.stop);
        assert(_col.start < _col.stop);
    }

    /// @brief Creates a **View** for the given **matrix**.
    /// @param matrix the matrix.
    /// @param row the range for the rows.
    /// @param col the range for the columns.
    View(MatrixType *matrix, Range row = Range::all(), Range col = Range::all())
        : _matrix(matrix),
          _row{ row },
          _col{ col }
    {
        assert(_row.start < _row.stop);
        assert(_col.start < _col.stop);
    }

    /// @brief Get the number of rows of the matrix.
    /// @returns the number of rows.
    inline std::size_t rows() const override
    {
        return std::min(_row.stop, _matrix->rows()) - _row.start;
    }

    /// @brief Get the number of columns of the matrix.
    /// @returns the number of columns.
    inline std::size_t cols() const override
    {
        return std::min(_col.stop, _matrix->cols()) - _col.start;
    }

    /// @brief Returns the total size of the matrix.
    /// @returns the total size of the matrix.
    inline std::size_t size() const override
    {
        return this->rows() * this->cols();
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    template <typename T2>
    inline auto &operator=(const MatrixBase<T2> &rhs)
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (std::size_t r = 0; r < rhs.rows(); ++r)
            for (std::size_t c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    inline auto &operator=(const View<MatrixType> &rhs)
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (std::size_t r = 0; r < rhs.rows(); ++r)
            for (std::size_t c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    template <typename T2>
    inline auto &operator=(const View<T2> &rhs)
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (std::size_t r = 0; r < rhs.rows(); ++r)
            for (std::size_t c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Construct a new Matrix object.
    /// @param data the data to use.
    /// @return a reference to this.
    auto &operator=(const std::initializer_list<std::initializer_list<T>> &data)
    {
        // Check the number of rows.
        if (0 == data.size()) {
            std::cerr << "View: first dimension of initializer list is empty.\n";
            return *this;
        }
        if (this->rows() != data.size()) {
            std::cerr << "View: first dimension of initializer list is incompatible (" << this->rows() << " vs " << data.size() << ").\n";
            return *this;
        }
        // Check the number of columns.
        if (0 == data.begin()->size()) {
            std::cerr << "View: second dimension of initializer list is empty.\n";
            return *this;
        }
        if (this->cols() != data.begin()->size()) {
            std::cerr << "View: second dimension of initializer list is incompatible (" << this->cols() << " vs " << data.begin()->size() << ").\n";
            return *this;
        }
        // Get an interator for the data.
        auto r_it = data.begin();
        decltype(r_it->begin()) c_it;
        for (std::size_t r = 0, c = 0; r < this->rows(); ++r, ++r_it) {
            c_it = r_it->begin();
            for (c = 0; c < this->cols(); ++c, ++c_it) {
                this->at(r, c) = *c_it;
            }
        }
        return *this;
    }

    /// @brief Returns a pointer to the internal data.
    /// @return the pointer.
    inline T *data() override
    {
        return _matrix->data();
    }

    /// @brief Returns a constant pointer to the internal data.
    /// @return the constant pointer.
    inline const T *data() const override
    {
        return _matrix->data();
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    inline T &operator[](std::size_t pos) override
    {
#ifdef ROW_MAJOR
        return this->at((pos / this->cols()), (pos % this->cols()));
#else
        return this->at((pos % this->rows()), (pos / this->rows()));
#endif
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    inline const T &operator[](std::size_t pos) const override
    {
#ifdef ROW_MAJOR
        return this->at((pos / this->cols()), (pos % this->cols()));
#else
        return this->at((pos % this->rows()), (pos / this->rows()));
#endif
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

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the reference to the accessed item.
    inline T &at(std::size_t row, std::size_t col) override
    {
        return _matrix->at(_row.start + row, _col.start + col);
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the const reference to the accessed item.
    inline const T &at(std::size_t row, std::size_t col) const override
    {
        return _matrix->at(_row.start + row, _col.start + col);
    }
};

template <typename T>
View<const Matrix<T>> Matrix<T>::operator()(Range row, Range col) const
{
    return View(this, row, col);
}

template <typename T>
View<Matrix<T>> Matrix<T>::operator()(Range row, Range col)
{
    return View(this, row, col);
}

/// @brief A view of another matrix-like structure, without a specific row and column.
template <typename MatrixType>
class CofactorView : public MatrixBase<typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::value_type, typename MatrixType::value_type>> {
public:
    /// The datatype of the referenced matrix.
    using T = typename std::conditional_t<std::is_const_v<MatrixType>, const typename MatrixType::value_type, typename MatrixType::value_type>;
    /// Pointer to the matrix.
    MatrixType *_matrix;
    /// Index of the excluded row.
    std::size_t _excluded_row;
    /// Index of the excluded column.
    std::size_t _excluded_col;

public:
    /// @brief Creates a new cofactor view.
    /// @param matrix the matrix of which to create the view.
    /// @param excluded_row the row to exclude.
    /// @param excluded_col the column to exclude.
    CofactorView(MatrixType *matrix,
                 std::size_t excluded_row,
                 std::size_t excluded_col)
        : _matrix(matrix),
          _excluded_row(excluded_row),
          _excluded_col(excluded_col)
    {
        // Nothing to do.
    }

    /// @brief Get the number of rows of the matrix.
    /// @returns the number of rows.
    inline std::size_t rows() const override
    {
        return _matrix->rows() - 1;
    }

    /// @brief Get the number of columns of the matrix.
    /// @returns the number of columns.
    inline std::size_t cols() const override
    {
        return _matrix->cols() - 1;
    }

    /// @brief Returns the total size of the matrix.
    /// @returns the total size of the matrix.
    inline std::size_t size() const override
    {
        return this->rows() * this->cols();
    }

    /// @brief Assign operator.
    /// @param rhs the other matrix.
    /// @returns a reference to this matrix.
    template <typename T2>
    inline auto &operator=(const MatrixBase<T2> &rhs)
    {
        assert(this->rows() == rhs.rows());
        assert(this->cols() == rhs.cols());
        for (std::size_t r = 0; r < rhs.rows(); ++r)
            for (std::size_t c = 0; c < rhs.cols(); ++c)
                this->at(r, c) = rhs(r, c);
        return *this;
    }

    /// @brief Returns a pointer to the content.
    /// @returns the pointer.
    inline T *data() override
    {
        return _matrix->data();
    }

    /// @brief Returns a constant pointer to the content.
    /// @returns the constant pointer.
    inline const T *data() const override
    {
        return _matrix->data();
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    inline T &operator[](std::size_t pos) override
    {
#ifdef ROW_MAJOR
        return this->at((pos / this->rows()), (pos % this->rows()));
#else
        return this->at((pos % this->cols()), (pos / this->cols()));
#endif
    }

    /// @brief Operator for accessing the matrix linearly.
    /// @param pos the liner position.
    /// @returns the reference to the accessed item.
    inline const T &operator[](std::size_t pos) const override
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

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the reference to the accessed item.
    inline T &at(std::size_t row, std::size_t col) override
    {
        return _matrix->at(row + (row >= _excluded_row), col + (col >= _excluded_col));
    }

    /// @brief Alternative function to access the matrix.
    /// @param row the accessed row.
    /// @param col the accessed column.
    /// @returns the const reference to the accessed item.
    inline const T &at(std::size_t row, std::size_t col) const override
    {
        return _matrix->at(row + (row >= _excluded_row), col + (col >= _excluded_col));
    }
};

} // namespace malg
