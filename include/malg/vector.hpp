/// @file vector.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief The vector class, used by the matrix class.

#pragma once

#include <cassert>
#include <cstdlib>
#include <initializer_list>
#include <utility>

//#define MEM_TRACE
#ifdef MEM_TRACE
#include <iostream>
#endif

namespace malg
{

/// @brief The vector class.
template <typename T>
class Vector {
public:
    /// The data types of the element of the vector.
    using value_type = T;
    /// The data types for iterating the element of the vector.
    using iterator = T *;
    /// The data types for iterating the element of the vector.
    using const_iterator = const T *;

    /// @brief Construct a new Vector object.
    Vector()
        : _data(nullptr),
          _size(0)
    {
        // Nothing to do.
    }

    /// @brief Construct a new Vector object.
    /// @param size dimension of the vector.
    explicit Vector(std::size_t size)
        : _data(nullptr),
          _size(0)
    {
        this->allocate(size);
    }

    /// @brief Construct a new Vector object.
    /// @param size dimension of the vector.
    /// @param value the vector will be initialized with it.
    Vector(std::size_t size, T value)
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(size)) {
            for (std::size_t i = 0; i < _size; ++i) {
                _data[i] = value;
            }
        }
    }

    /// @brief Construct a new Vector object.
    /// @param size dimension of the vector.
    /// @param data the initialization data.
    Vector(std::size_t size, const std::initializer_list<T> &data)
        : _data(nullptr),
          _size(0)
    {
        assert(size == data.size());
        if (this->allocate(size)) {
            auto it = data.begin();
            for (std::size_t i = 0; i < _size; ++i) {
                _data[i] = *it++;
            }
        }
    }

    /// @brief Construct a new Vector object.
    /// @param values list of values.
    Vector(const std::initializer_list<T> &values)
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(values.size())) {
            auto it = values.begin();
            for (std::size_t i = 0; i < _size; ++i) {
                _data[i] = *it++;
            }
        }
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    Vector(const Vector<T> &other)
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(other.size())) {
            for (std::size_t i = 0; i < _size; ++i) {
                _data[i] = other[i];
            }
        }
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    template <typename T2>
    explicit Vector(const Vector<T2> &other)
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(other.size())) {
            for (std::size_t i = 0; i < _size; ++i) {
                _data[i] = other[i];
            }
        }
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    Vector(Vector<T> &&other)
        : _data(),
          _size()
    {
        if (this != &other) {
            this->deallocate();
            std::swap(_data, other._data);
            std::swap(_size, other._size);
        }
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    template <typename T2>
    explicit Vector(Vector<T2> &&other)
        : _data(),
          _size()
    {
        if (this != &other) {
            this->deallocate();
            std::swap(_data, other._data);
            std::swap(_size, other._size);
        }
    }

    /// @brief Destroy the Vector object.
    ~Vector()
    {
        this->deallocate();
    }

    /// @brief Gives access to the data.
    /// @returns pointer to the data.
    inline const_iterator data() const
    {
        return _data;
    }

    /// @brief Gives access to the data.
    /// @returns pointer to the data.
    inline iterator data()
    {
        return _data;
    }

    /// @brief Returns the dimension of the vector.
    /// @returns the dimension of the vector.
    inline std::size_t size() const
    {
        return _size;
    }

    /// @brief Checks if the vector is empty.
    /// @returns true if it is empty.
    /// @returns false otherwise.
    inline bool empty() const
    {
        return _size == 0;
    }

    /// @brief Resizes the vector based on the size of the **other**.
    /// @param other the other vector.
    /// @return a reference to this vector.
    inline auto &resize(const Vector &other)
    {
        this->allocate(other.size());
        return *this;
    }

    /// @brief Resizes the vector based on the given **size**.
    /// @param size the new size.
    /// @return a reference to this vector.
    inline auto &resize(std::size_t size)
    {
        this->allocate(size);
        return *this;
    }

    /// @brief Operator for accessing the vector.
    /// @param pos the position.
    /// @returns the reference to the accessed item.
    inline auto &operator[](std::size_t pos)
    {
        return _data[pos];
    }

    /// @brief Operator for accessing the vector.
    /// @param pos the position.
    /// @returns the reference to the accessed item.
    inline const T &operator[](std::size_t pos) const
    {
        return _data[pos];
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @returns this vector.
    inline auto &operator=(const Vector<T> &other)
    {
        if (this != &other) {
            if (this->allocate(other.size())) {
                for (std::size_t i = 0; i < _size; ++i) {
                    _data[i] = other[i];
                }
            }
        }
        return *this;
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @returns this vector.
    template <typename T2>
    inline auto &operator=(const Vector<T2> &other)
    {
        if (this->allocate(other.size())) {
            for (std::size_t i = 0; i < _size; ++i) {
                _data[i] = other[i];
            }
        }
        return *this;
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @returns this vector.
    inline auto &operator=(Vector<T> &&other)
    {
        if (this != &other) {
            this->deallocate();
            std::swap(_data, other._data);
            std::swap(_size, other._size);
        }
        return *this;
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @returns this vector.
    template <typename T2>
    inline auto &operator=(Vector<T2> &&other)
    {
        if (this != &other) {
            this->deallocate();
            std::swap(_data, other._data);
            std::swap(_size, other._size);
        }
        return *this;
    }

    /// @brief Provides a pointer to the begining of the vector.
    /// @returns the pointer.
    inline const_iterator begin() const
    {
        return _data;
    }

    /// @brief Provides a pointer to the begining of the vector.
    /// @returns the pointer.
    inline const_iterator cbegin() const
    {
        return _data;
    }

    /// @brief Provides a pointer to the begining of the vector.
    /// @returns the pointer.
    inline iterator begin()
    {
        return _data;
    }

    /// @brief Provides a pointer to the end of the vector.
    /// @returns the pointer.
    inline const_iterator end() const
    {
        return _data + _size;
    }

    /// @brief Provides a pointer to the end of the vector.
    /// @returns the pointer.
    inline const_iterator cend() const
    {
        return _data + _size;
    }

    /// @brief Provides a pointer to the end of the vector.
    /// @returns the pointer.
    inline iterator end()
    {
        return _data + _size;
    }

private:
    inline auto deallocate()
    {
        if (_size > 0) {
#ifdef MEM_TRACE
            std::cout << "Deallocationg " << (sizeof(T) * _size) << " bytes (" << _data << ")\n";
#endif
            std::free(_data);
            _data = nullptr;
            _size = 0;
            return true;
        }
        return false;
    }

    inline auto allocate(std::size_t size)
    {
        if (size == 0) {
            return this->deallocate();
        }
        // If the size is the same, we don't need to act.
        if (size == _size) {
            return true;
        }
        // Re-allocate the data.
        _data = static_cast<T *>(std::realloc(_data, sizeof(T) * size));
#ifdef MEM_TRACE
        if (_size == 0)
            std::cout << "Allocating    " << (sizeof(T) * size) << " bytes (" << _data << ")\n";
        else
            std::cout << "Re-allocating " << (sizeof(T) * _size) << " bytes to " << (sizeof(T) * size) << " bytes (" << _data << ")\n";
#endif
        // If we are allocating more size, we clean the memory.
        if (size > _size) {
            for (std::size_t i = _size; i < size; ++i) {
                _data[i] = T(0.);
            }
        }
        // If it is the first time we allocate, clean all the data.
        else if (_size == 0) {
            for (std::size_t i = 0; i < size; ++i) {
                _data[i] = T(0.);
            }
        }
        // Set the size.
        _size = size;
        return true;
    }

    /// The data.
    T *_data;
    /// The size of the vector.
    std::size_t _size;
};

} // namespace malg