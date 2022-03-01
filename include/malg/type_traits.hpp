/// @file complex_type_traits.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Helper functions for dealing with template types.

#pragma once

#include <type_traits>
#include <complex>

namespace malg
{

template <typename T>
struct is_complex : public std::false_type {
};

template <typename T>
struct is_complex<std::complex<T>> : public std::true_type {
    using type = T;
};

template <typename T>
inline constexpr bool is_complex_v = is_complex<T>::value;

template <typename T>
struct extract_value : public std::false_type {
    using type = T;
};

template <typename T>
struct extract_value<std::complex<T>> : public std::true_type {
    using type = T;
};

template <typename T>
using extract_value_t = extract_value<T>::type;

template <
    typename T1,
    typename T2>
using extract_common_type_t =
    std::remove_const_t<
        std::conditional_t<
            is_complex_v<T1> || is_complex_v<T2>,
            std::complex<std::common_type_t<extract_value_t<T1>, extract_value_t<T2>>>,
            std::common_type_t<T1, T2>>>;

} // namespace malg