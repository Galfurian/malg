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

template <typename T>
class Vector {
public:
    /// @brief Construct a new Vector object.
    Vector()
        : _data(nullptr),
          _size(0)
    {
        // Nothing to do.
    }

    /// @brief Construct a new Vector object.
    /// @param size dimension of the vector.
    Vector(unsigned size)
        : _data(nullptr),
          _size(0)
    {
        this->allocate(size);
    }

    /// @brief Construct a new Vector object.
    /// @param size dimension of the vector.
    /// @param value the vector will be initialized with it.
    Vector(unsigned size, T value)
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(size))
            for (unsigned i = 0; i < _size; ++i)
                _data[i] = value;
    }

    /// @brief Construct a new Vector object.
    /// @param size dimension of the vector.
    /// @param data the initialization data.
    Vector(unsigned size, const std::initializer_list<T> &data)
        : _data(nullptr),
          _size(0)
    {
        assert(size == data.size());
        if (this->allocate(size)) {
            auto it = data.begin();
            for (unsigned i = 0; i < _size; ++i)
                _data[i] = *it++;
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
            for (unsigned i = 0; i < _size; ++i)
                _data[i] = *it++;
        }
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    Vector(const Vector<T> &other)
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(other.size()))
            for (unsigned i = 0; i < _size; ++i)
                _data[i] = other[i];
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    template <typename T2>
    Vector(const Vector<T2> &other)
        : _data(nullptr),
          _size(0)
    {
        if (this->allocate(other.size()))
            for (unsigned i = 0; i < _size; ++i)
                _data[i] = other[i];
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    Vector(Vector<T> &&other)
        : _data(std::move(other._data)),
          _size(other._size)
    {
        other._data = nullptr;
        other._size = 0;
    }

    /// @brief Copy construct a new Vector object.
    /// @param other the other vector.
    template <typename T2>
    Vector(Vector<T2> &other)
        : _data(std::move(other._data)),
          _size(other._size)
    {
        other._data = nullptr;
        other._size = 0;
    }

    /// @brief Destroy the Vector object.
    ~Vector()
    {
        this->deallocate();
    }

    /// @brief Gives access to the data.
    /// @return pointer to the data.
    inline const T *data() const
    {
        return _data;
    }

    /// @brief Gives access to the data.
    /// @return pointer to the data.
    inline T *data()
    {
        return _data;
    }

    /// @brief Returns the dimension of the vector.
    /// @return the dimension of the vector.
    inline unsigned size() const
    {
        return _size;
    }

    /// @brief Checks if the vector is empty.
    /// @return true if it is empty.
    /// @return false otherwise.
    inline virtual bool empty() const
    {
        return _size == 0;
    }

    auto &resize(unsigned size)
    {
        this->allocate(size);
        return *this;
    }

    /// @brief Operator for accessing the vector.
    /// @param pos the position.
    /// @return the reference to the accessed item.
    T &operator[](unsigned pos)
    {
        assert(pos < _size);
        return _data[pos];
    }

    /// @brief Operator for accessing the vector.
    /// @param pos the position.
    /// @return the reference to the accessed item.
    const T &operator[](unsigned pos) const
    {
        assert(pos < _size);
        return _data[pos];
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @return this vector.
    Vector<T> &operator=(const Vector<T> &other)
    {
        if (this->allocate(other.size()))
            for (unsigned i = 0; i < _size; ++i)
                _data[i] = other[i];
        return *this;
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @return this vector.
    template <typename T2>
    Vector<T> &operator=(const Vector<T2> &other)
    {
        if (this->allocate(other.size()))
            for (unsigned i = 0; i < _size; ++i)
                _data[i] = other[i];
        return *this;
    }

    /// @brief Assignment operator.
    /// @param other the other vector.
    /// @return this vector.
    Vector<T> &operator=(Vector<T> &&other)
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
    Vector<T> &operator=(Vector<T2> &&other)
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
    inline const T *begin() const
    {
        return _data;
    }

    /// @brief Provides a pointer to the begining of the vector.
    /// @return the pointer.
    inline T *begin()
    {
        return _data;
    }

    /// @brief Provides a pointer to the end of the vector.
    /// @return the pointer.
    inline const T *end() const
    {
        return &_data[_size];
    }

    /// @brief Provides a pointer to the end of the vector.
    /// @return the pointer.
    inline T *end()
    {
        return &_data[_size];
    }

private:
    inline bool deallocate()
    {
        if (_size > 0) {
            std::free(_data);
            _data = nullptr;
            _size = 0;
            return true;
        }
        return false;
    }

    inline bool allocate(unsigned size)
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
            for (unsigned i = _size; i < size; ++i)
                _data[i] = T(0.);
        }
        // If it is the first time we allocate, clean all the data.
        else if (_size == 0) {
            for (unsigned i = 0; i < size; ++i)
                _data[i] = T(0.);
        }
        // Set the size.
        _size = size;
        return true;
    }

    /// The data.
    T *_data;
    /// The size of the vector.
    unsigned _size;
};

} // namespace malg