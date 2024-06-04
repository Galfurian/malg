/// @file matrix.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief The matrix class.

#pragma once

#include "malg/fixed_size/vector.hpp"

namespace malg
{

template <class T, std::size_t N1, std::size_t N2 = N1>
class Matrix : public malg::Vector<malg::Vector<T, N2>, N1> {
public:
    // Constructor using initializer list
    Matrix(std::initializer_list<std::initializer_list<T>> init_list)
    {
        if (init_list.size() != N1) {
            throw std::length_error("Initializer list size does not match Vector size.");
        }
        std::size_t r = 0;
        for (const auto &it : init_list) {
            if (it.size() != N2) {
                throw std::length_error("Initializer list size does not match Vector size.");
            }
            std::copy(it.begin(), it.end(), this->_data[r++].begin());
        }
    }

    inline constexpr std::size_t size() const noexcept
    {
        return N1 * N2;
    }

    inline constexpr bool empty() const noexcept
    {
        return (N1 == 0) || (N2 == 0);
    }

    inline T &operator()(std::size_t index)
    {
#ifdef ROW_MAJOR
        std::size_t row = index / N2;
        std::size_t col = index % N2;
#else
        std::size_t col = index / N1;
        std::size_t row = index % N1;
#endif
        return this->_data[row][col];
    }

    inline const T &operator()(std::size_t index) const
    {
#ifdef ROW_MAJOR
        std::size_t row = index / N2;
        std::size_t col = index % N2;
#else
        std::size_t col = index / N1;
        std::size_t row = index % N1;
#endif
        return this->_data[row][col];
    }
};

} // namespace malg
