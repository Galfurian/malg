/// @file type_traits.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Helper functions for dealing with template types.

#pragma once

#include <complex>
#include <type_traits>

namespace malg
{

/// @brief This structure has a twofold purpose, it can be used to check if a
/// type is a complex, or to extract the internal type of a complex.
/// @tparam T a normal type.
template <typename T>
struct is_complex : public std::false_type {
    using type = T; ///< An alias for T.
};

/// @brief This structure has a twofold purpose, it can be used to check if a
/// type is a complex, or to extract the internal type of a complex.
/// @tparam T the internal type of a complex.
template <typename T>
struct is_complex<std::complex<T>> : public std::true_type {
    using type = T; ///< An alias for T.
};

/// @brief Alias template for checking if a type is complex.F
/// @tparam T the type we are analyzing.
/// @details
/// is_complex_v< double               > -> false
/// is_complex_v< std::complex<double> > -> true
template <typename T>
inline const bool is_complex_v = is_complex<T>::value;

/// @brief Alias template for extracting the internal type of a complex.
/// @tparam T the type we are analyzing.
/// @details
/// is_complex_t< double               > -> double
/// is_complex_t< std::complex<double> > -> double
template <typename T>
using is_complex_t = typename is_complex<T>::type;

/// @brief Given two types it extracts the common type between them.
/// @tparam T1 the first type.
/// @tparam T2 the second type.
/// @details
/// Basically, if either one of the two types is a complex, the output will be a complex.
/// extract_common_type_t< int                 , double              > -> double
/// extract_common_type_t< int                 , std::complex<float> > -> std::complex<float>
/// extract_common_type_t< std::complex<double>, std::complex<float> > -> std::complex<double>
template <typename T1, typename T2>
using extract_common_type_t = std::conditional_t<is_complex_v<T1> || is_complex_v<T2>, std::complex<std::common_type_t<is_complex_t<T1>, is_complex_t<T2>>>, std::common_type_t<T1, T2>>;

/// @brief Given three types it extracts the common type between them.
/// @tparam T1 the first type.
/// @tparam T2 the second type.
/// @details
/// extract_common_type_t< int                 , double              > -> double
/// extract_common_type_t< int                 , std::complex<float> > -> std::complex<float>
/// extract_common_type_t< std::complex<double>, std::complex<float> > -> std::complex<double>
template <typename T1, typename T2, typename T3>
using extract_common_type_3_t = malg::extract_common_type_t<malg::extract_common_type_t<T1, T2>, T3>;

} // namespace malg