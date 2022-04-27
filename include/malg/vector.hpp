/// @file vector.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief The vector class, used by the matrix class.

#pragma once

#include <initializer_list>
#include <cstdlib>
#include <cassert>
#include <utility>

namespace malg
{

using size_type_t = unsigned;

template <typename T>
class Vector {
public:
    using value_type = T;

    /// @brief Construct a new Vector object.
    constexpr Vector() noexcept
        : _data(nullptr),
          _size(0)
    {
        // Nothing to do.
    }

    /// @brief Construct a new Vector object.
    /// @param size dimension of the vector.
    constexpr Vector(size_type_t size) noexcept
        : _data(nullptr),
          _size(0)
    {
        this->allocate(size);
    }

    /// @brief Construct a new Vector object.
    /// @param size dimension of the vector.
    /// @param value the vector will be initialized with it.
    constexpr Vector(size_type_t size, T value) noexcept
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(size))
            for (size_type_t i = 0; i < _size; ++i)
                _data[i] = value;
    }

    /// @brief Construct a new Vector object.
    /// @param size dimension of the vector.
    /// @param data the initialization data.
    constexpr Vector(size_type_t size, const std::initializer_list<T> &data) noexcept
        : _data(nullptr),
          _size(0)
    {
        assert(size == data.size());
        if (this->allocate(size)) {
            auto it = data.begin();
            for (size_type_t i = 0; i < _size; ++i)
                _data[i] = *it++;
        }
    }

    /// @brief Construct a new Vector object.
    /// @param values list of values.
    constexpr Vector(const std::initializer_list<T> &values) noexcept
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(values.size())) {
            auto it = values.begin();
            for (size_type_t i = 0; i < _size; ++i)
                _data[i] = *it++;
        }
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    constexpr Vector(const Vector<T> &other) noexcept
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(other.size()))
            for (size_type_t i = 0; i < _size; ++i)
                _data[i] = other[i];
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    template <typename T2>
    constexpr Vector(const Vector<T2> &other) noexcept
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(other.size()))
            for (size_type_t i = 0; i < _size; ++i)
                _data[i] = other[i];
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    constexpr Vector(Vector<T> &&other) noexcept
        : _data(std::move(other._data)),
          _size(other._size)
    {
        other._data = nullptr;
        other._size = 0;
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    template <typename T2>
    constexpr Vector(Vector<T2> &other) noexcept
        : _data(std::move(other._data)),
          _size(other._size)
    {
        other._data = nullptr;
        other._size = 0;
    }

    /// @brief Destroy the Vector object.
    constexpr ~Vector() noexcept
    {
        this->deallocate();
    }

    /// @brief Gives access to the data.
    /// @return pointer to the data.
    inline constexpr const auto data() const noexcept
    {
        return _data;
    }

    /// @brief Gives access to the data.
    /// @return pointer to the data.
    inline constexpr auto data() noexcept
    {
        return _data;
    }

    /// @brief Returns the dimension of the vector.
    /// @return the dimension of the vector.
    inline constexpr auto size() const noexcept
    {
        return _size;
    }

    /// @brief Checks if the vector is empty.
    /// @return true if it is empty.
    /// @return false otherwise.
    inline constexpr auto empty() const noexcept
    {
        return _size == 0;
    }

    inline constexpr auto &resize(const Vector &other) noexcept
    {
        this->allocate(other.size());
        return *this;
    }

    inline constexpr auto &resize(size_type_t size) noexcept
    {
        this->allocate(size);
        return *this;
    }

    /// @brief Operator for accessing the vector.
    /// @param pos the position.
    /// @return the reference to the accessed item.
    inline constexpr auto &operator[](size_type_t pos) noexcept
    {
        return _data[pos];
    }

    /// @brief Operator for accessing the vector.
    /// @param pos the position.
    /// @return the reference to the accessed item.
    inline constexpr const auto &operator[](size_type_t pos) const noexcept
    {
        return _data[pos];
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @return this vector.
    inline constexpr auto &operator=(const Vector<T> &other) noexcept
    {
        if (this->allocate(other.size()))
            for (size_type_t i = 0; i < _size; ++i)
                _data[i] = other[i];
        return *this;
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @return this vector.
    template <typename T2>
    inline constexpr auto &operator=(const Vector<T2> &other) noexcept
    {
        if (this->allocate(other.size()))
            for (size_type_t i = 0; i < _size; ++i)
                _data[i] = other[i];
        return *this;
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @return this vector.
    inline constexpr auto &operator=(Vector<T> &&other) noexcept
    {
        this->deallocate();
        _data       = std::move(other._data);
        _size       = other._size;
        other._data = nullptr;
        other._size = 0;
        return *this;
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @return this vector.
    template <typename T2>
    inline constexpr auto &operator=(Vector<T2> &&other) noexcept
    {
        this->deallocate();
        _data       = std::move(other._data);
        _size       = other._size;
        other._data = nullptr;
        other._size = 0;
        return *this;
    }

    /// @brief Provides a pointer to the begining of the vector.
    /// @return the pointer.
    inline constexpr const auto begin() const noexcept
    {
        return _data;
    }

    /// @brief Provides a pointer to the begining of the vector.
    /// @return the pointer.
    inline constexpr auto begin() noexcept
    {
        return _data;
    }

    /// @brief Provides a pointer to the end of the vector.
    /// @return the pointer.
    inline constexpr const auto end() const noexcept
    {
        return &_data[_size];
    }

    /// @brief Provides a pointer to the end of the vector.
    /// @return the pointer.
    inline constexpr auto end() noexcept
    {
        return &_data[_size];
    }

private:
    inline constexpr auto deallocate() noexcept
    {
        if (_size > 0) {
            std::free(_data);
            _data = nullptr;
            _size = 0;
            return true;
        }
        return false;
    }

    inline constexpr auto allocate(size_type_t size) noexcept
    {
        if (size == 0)
            return this->deallocate();
        // If the size is the same, we don't need to act.
        if (size == _size)
            return true;
        // Re-allocate the data.
        _data = static_cast<T *>(std::realloc(_data, sizeof(T) * size));
        // If we are allocating more size, we clean the memory.
        if (size > _size) {
            for (size_type_t i = _size; i < size; ++i)
                _data[i] = T(0.);
        }
        // If it is the first time we allocate, clean all the data.
        else if (_size == 0) {
            for (size_type_t i = 0; i < size; ++i)
                _data[i] = T(0.);
        }
        // Set the size.
        _size = size;
        return true;
    }

    /// The data.
    T *_data;
    /// The size of the vector.
    size_type_t _size;
};

} // namespace malg