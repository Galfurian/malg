/// @file view.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief VectorView classes, which stores nothing but provide structured access to
/// real matrices.

#pragma once

#include "malg/fixed_size/matrix.hpp"

namespace malg
{

template <typename T, std::size_t N, std::size_t Start = 0, std::size_t End = N>
class VectorView {
public:
    using value_type             = T;
    using pointer                = value_type *;
    using const_pointer          = const value_type *;
    using reference              = value_type &;
    using const_reference        = const value_type &;
    using iterator               = value_type *;
    using const_iterator         = const value_type *;
    using size_type              = std::size_t;
    using difference_type        = std::ptrdiff_t;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    // Constructor takes a reference to the original data
    VectorView(VectorBase<T, N> &vec)
        : data(vec)
    {
        static_assert(Start < End && End <= N, "Invalid view range");
    }

    std::size_t size() const
    {
        return End - Start;
    }

    // Element access.
    reference operator[](size_type idx)
    {
        if (idx >= End - Start) {
            throw std::out_of_range("Index out of range");
        }
        return data[Start + idx];
    }

    const_reference operator[](size_type idx) const
    {
        if (idx >= End - Start) {
            throw std::out_of_range("Index out of range");
        }
        return data[Start + idx];
    }

    [[nodiscard]] constexpr iterator begin() noexcept
    {
        return iterator(&data[Start]);
    }

    [[nodiscard]] constexpr const_iterator begin() const noexcept
    {
        return const_iterator(&data[Start]);
    }

    [[nodiscard]] constexpr iterator end() noexcept
    {
        return iterator(&data[End]);
    }

    [[nodiscard]] constexpr const_iterator end() const noexcept
    {
        return const_iterator(&data[End]);
    }

    [[nodiscard]] constexpr reverse_iterator rbegin() noexcept
    {
        return reverse_iterator(this->end());
    }

    [[nodiscard]] constexpr const_reverse_iterator rbegin() const noexcept
    {
        return const_reverse_iterator(this->end());
    }

    [[nodiscard]] constexpr reverse_iterator rend() noexcept
    {
        return reverse_iterator(this->begin());
    }

    [[nodiscard]] constexpr const_reverse_iterator rend() const noexcept
    {
        return const_reverse_iterator(this->begin());
    }

    [[nodiscard]] constexpr const_iterator cbegin() const noexcept
    {
        return const_iterator(&data[Start]);
    }

    [[nodiscard]] constexpr const_iterator cend() const noexcept
    {
        return const_iterator(&data[End]);
    }

    [[nodiscard]] constexpr const_reverse_iterator crbegin() const noexcept
    {
        return const_reverse_iterator(this->end());
    }

    [[nodiscard]] constexpr const_reverse_iterator crend() const noexcept
    {
        return const_reverse_iterator(this->begin());
    }

private:
    VectorBase<T, N> &data;
};

template <typename T, std::size_t N1, std::size_t N2, std::size_t RowStart, std::size_t RowEnd, std::size_t ColStart, std::size_t ColEnd>
class MatrixView {
public:
    using value_type             = T;
    using pointer                = value_type *;
    using const_pointer          = const value_type *;
    using reference              = value_type &;
    using const_reference        = const value_type &;
    using iterator               = value_type *;
    using const_iterator         = const value_type *;
    using size_type              = std::size_t;
    using difference_type        = std::ptrdiff_t;
    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;

    // Constructor takes a reference to the original data
    MatrixView(Matrix<T, N1, N2> &mat)
        : data(mat)
    {
        static_assert(RowStart < RowEnd && RowEnd <= N1 && ColStart < ColEnd && ColEnd <= N2, "Invalid view range");
    }

    std::size_t size() const
    {
        return RowEnd - RowStart;
    }

    VectorView<T, N2, ColStart, ColEnd> operator[](size_type idx)
    {
        if (idx >= RowEnd - RowStart) {
            throw std::out_of_range("Row index out of range");
        }
        return VectorView<T, N2, ColStart, ColEnd>(data[RowStart + idx]);
    }

    const VectorView<T, N2, ColStart, ColEnd> operator[](size_type idx) const
    {
        if (idx >= RowEnd - RowStart) {
            throw std::out_of_range("Row index out of range");
        }
        return VectorView<T, N2, ColStart, ColEnd>(data[RowStart + idx]);
    }

private:
    Matrix<T, N1, N2> &data;
};

template <std::size_t Start, std::size_t End, typename T, std::size_t N>
[[nodiscard]] constexpr auto view(VectorBase<T, N> &v)
{
    return VectorView<T, N, Start, End>(v);
}

template <std::size_t RowStart, std::size_t RowEnd, std::size_t ColStart, std::size_t ColEnd, typename T, std::size_t N1, std::size_t N2>
[[nodiscard]] constexpr auto view(Matrix<T, N1, N2> &v)
{
    return MatrixView<T, N1, N2, RowStart, RowEnd, ColStart, ColEnd>(v);
}

} // namespace malg
