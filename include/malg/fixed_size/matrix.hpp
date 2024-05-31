/// @file matrix.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief The matrix class.

#pragma once

#include "malg/fixed_size/vector.hpp"

namespace malg
{

template <class T, std::size_t N1, std::size_t N2 = N1>
using Matrix = malg::Vector<malg::Vector<T, N2>, N1>;

} // namespace malg
