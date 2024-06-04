/// @file vector.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief The vector class, used by the matrix class.

#pragma once

#include <cstddef>
#include <iterator>

namespace malg
{

/// @brief The vector class.
template <typename T, std::size_t N>
struct VectorBase {
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

    inline constexpr size_type size() const noexcept
    {
        return N;
    }

    inline constexpr size_type max_size() const noexcept
    {
        return N;
    }

    inline constexpr bool empty() const noexcept
    {
        return N == 0;
    }

    // Element access.
    virtual reference operator[](size_type idx) = 0;

    virtual const_reference operator[](size_type idx) const = 0;
};

/// @brief The vector class.
template <typename T, std::size_t N>
struct Vector : public VectorBase<T, N> {
    using typename VectorBase<T, N>::value_type;
    using typename VectorBase<T, N>::pointer;
    using typename VectorBase<T, N>::const_pointer;
    using typename VectorBase<T, N>::reference;
    using typename VectorBase<T, N>::const_reference;
    using typename VectorBase<T, N>::iterator;
    using typename VectorBase<T, N>::const_iterator;
    using typename VectorBase<T, N>::size_type;
    using typename VectorBase<T, N>::difference_type;
    using typename VectorBase<T, N>::reverse_iterator;
    using typename VectorBase<T, N>::const_reverse_iterator;

    // Default constructor
    Vector() = default;

    // Constructor using initializer list
    Vector(std::initializer_list<T> init_list)
    {
        if (init_list.size() != N) {
            throw std::length_error("Initializer list size does not match Vector size.");
        }
        std::copy(init_list.begin(), init_list.end(), _data);
    }

    // Copy constructor
    Vector(const Vector &other)
    {
        std::copy(other._data, other._data + N, _data);
    }

    // Copy assignment operator
    Vector &operator=(const Vector &other)
    {
        if (this != &other) {
            std::copy(other._data, other._data + N, _data);
        }
        return *this;
    }

    // Move constructor
    Vector(Vector &&other) noexcept
    {
        std::move(other._data, other._data + N, _data);
    }

    // Move assignment operator
    Vector &operator=(Vector &&other) noexcept
    {
        if (this != &other) {
            std::move(other._data, other._data + N, _data);
        }
        return *this;
    }

    constexpr void fill(const T &value)
    {
        std::fill_n(this->begin(), this->size(), value);
    }

    constexpr void swap(malg::Vector<T, N> &other) noexcept
    {
        std::swap_ranges(this->begin(), this->end(), other.begin());
    }

    [[nodiscard]] constexpr iterator begin() noexcept
    {
        return iterator(_data);
    }

    [[nodiscard]] constexpr const_iterator begin() const noexcept
    {
        return const_iterator(_data);
    }

    [[nodiscard]] constexpr iterator end() noexcept
    {
        return iterator(_data + N);
    }

    [[nodiscard]] constexpr const_iterator end() const noexcept
    {
        return const_iterator(_data + N);
    }

    [[nodiscard]] constexpr reverse_iterator rbegin() noexcept
    {
        return reverse_iterator(end());
    }

    [[nodiscard]] constexpr const_reverse_iterator rbegin() const noexcept
    {
        return const_reverse_iterator(end());
    }

    [[nodiscard]] constexpr reverse_iterator rend() noexcept
    {
        return reverse_iterator(begin());
    }

    [[nodiscard]] constexpr const_reverse_iterator rend() const noexcept
    {
        return const_reverse_iterator(begin());
    }

    [[nodiscard]] constexpr const_iterator cbegin() const noexcept
    {
        return const_iterator(_data);
    }

    [[nodiscard]] constexpr const_iterator cend() const noexcept
    {
        return const_iterator(_data + N);
    }

    [[nodiscard]] constexpr const_reverse_iterator crbegin() const noexcept
    {
        return const_reverse_iterator(this->end());
    }

    [[nodiscard]] constexpr const_reverse_iterator crend() const noexcept
    {
        return const_reverse_iterator(this->begin());
    }

    // Element access.
    reference operator[](size_type idx) override
    {
        if (idx >= N) {
            throw std::out_of_range("Index out of range");
        }
        return _data[idx];
    }

    const_reference operator[](size_type idx) const override
    {
        if (idx >= N) {
            throw std::out_of_range("Index out of range");
        }
        return _data[idx];
    }

    constexpr reference at(size_type n)
    {
        if (n >= N) {
            std::__throw_out_of_range_fmt("Vector::at: n (%u) >= N (%u)", n, N);
        }
        return _data[n];
    }

    constexpr const_reference at(size_type n) const
    {
        if (n >= N) {
            std::__throw_out_of_range_fmt("Vector::at: n (%u) >= N (%u)", n, N);
        }
        return _data[n];
    }

    [[nodiscard]] constexpr reference front() noexcept
    {
        if (N == 0) {
            std::__throw_length_error("Vector::front: vector is empty");
        }
        return _data[0];
    }

    constexpr const_reference front() const noexcept
    {
        if (N == 0) {
            std::__throw_length_error("Vector::front: vector is empty");
        }
        return _data[0];
    }

    [[nodiscard]] constexpr reference back() noexcept
    {
        if (N == 0) {
            std::__throw_length_error("Vector::back: vector is empty");
        }
        return _data[N - 1];
    }

    constexpr const_reference back() const noexcept
    {
        if (N == 0) {
            std::__throw_length_error("Vector::back: vector is empty");
        }
        return _data[N - 1];
    }

    [[nodiscard]] constexpr pointer data() noexcept
    {
        return &_data;
    }

    [[nodiscard]] constexpr const_pointer data() const noexcept
    {
        return &_data;
    }

protected:
    T _data[N];
};

} // namespace malg
