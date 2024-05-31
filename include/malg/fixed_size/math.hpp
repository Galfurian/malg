/// @file math.hpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)
/// @brief Basic mathematical functions.

#pragma once

#include "malg/fixed_size/vector.hpp"
#include "malg/fixed_size/matrix.hpp"

#include <type_traits>
#include <utility>

namespace malg::details
{

template <typename T, typename Function, std::size_t... N>
[[nodiscard]] constexpr auto apply_unary_helper(
    const malg::VectorBase<T, sizeof...(N)> &l,
    Function func,
    std::integer_sequence<std::size_t, N...> s)
{
    using return_type_t = malg::Vector<decltype(func(l[0])), s.size()>;
    return return_type_t{ func(l[N])... };
}

template <class T, class U, class F, std::size_t... N>
[[nodiscard]] constexpr auto apply_binary_helper(
    const malg::VectorBase<T, sizeof...(N)> &l,
    const malg::VectorBase<U, sizeof...(N)> &r,
    F func,
    std::integer_sequence<std::size_t, N...> s)
{
    using return_type_t = malg::Vector<decltype(func(l[0], r[0])), s.size()>;
    return return_type_t{ func(l[N], r[N])... };
}

template <class T, class U, class F, std::size_t... N>
[[nodiscard]] constexpr auto apply_binary_helper(
    const U &l,
    const malg::VectorBase<T, sizeof...(N)> &r,
    F func,
    std::integer_sequence<std::size_t, N...> s)
{
    using return_type_t = malg::Vector<decltype(func(l, r[0])), s.size()>;
    return return_type_t{ func(l, r[N])... };
}

template <class T, class U, class F, std::size_t... N>
[[nodiscard]] constexpr auto apply_binary_helper(
    const malg::VectorBase<T, sizeof...(N)> &l,
    const U &r,
    F func,
    std::integer_sequence<std::size_t, N...> s)
{
    using return_type_t = malg::Vector<decltype(func(l[0], r)), s.size()>;
    return return_type_t{ func(l[N], r)... };
}

template <std::size_t N, class F, class T>
constexpr auto apply_unary(F func, T &&arg)
{
    return apply_unary_helper(arg, func, std::make_integer_sequence<std::size_t, N>{});
}

template <std::size_t N, class F, class T1, class T2>
constexpr auto apply_binary(F func, const T1 &arg1, const T2 &arg2)
{
    return apply_binary_helper(arg1, arg2, func, std::make_integer_sequence<std::size_t, N>{});
}

struct negate {
    template <class T>
    [[nodiscard]] constexpr auto operator()(const T &l) const noexcept
    {
        return -l;
    }
};

struct plus {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l + r;
    }
};

struct minus {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l - r;
    }
};

struct multiplies {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l * r;
    }
};

struct divides {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l / r;
    }
};

struct equal {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l == r;
    }
};

struct not_equal_to {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l != r;
    }
};

struct less {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l < r;
    }
};

struct less_equal {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l <= r;
    }
};

struct greater {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l > r;
    }
};

struct greater_equal {
    template <class T, class U>
    [[nodiscard]] constexpr auto operator()(const T &l, const U &r) const noexcept
    {
        return l >= r;
    }
};

} // namespace malg::details

#define DEFINE_UNARY_VECTOR_OPERATIONS(OP, FUNC)                     \
    template <class T, std::size_t N>                                \
    [[nodiscard]] constexpr auto OP(const malg::VectorBase<T, N> &l) \
    {                                                                \
        return malg::details::apply_unary<N>(FUNC{}, l);             \
    }

#define DEFINE_BINARY_VECTOR_OPERATIONS(OP, FUNC)                                                                     \
    template <class T, class U, std::size_t N>                                                                        \
    [[nodiscard]] constexpr auto OP(const malg::VectorBase<T, N> &l, const malg::VectorBase<U, N> &r)                 \
    {                                                                                                                 \
        return malg::details::apply_binary<N>(FUNC{}, l, r);                                                          \
    }                                                                                                                 \
    template <class T, class U, std::size_t N, typename = typename std::enable_if_t<std::is_arithmetic<U>::value, U>> \
    [[nodiscard]] constexpr auto OP(const U &l, const malg::VectorBase<T, N> &r)                                      \
    {                                                                                                                 \
        return malg::details::apply_binary<N>(FUNC{}, l, r);                                                          \
    }                                                                                                                 \
    template <class T, class U, std::size_t N, typename = typename std::enable_if_t<std::is_arithmetic<U>::value, U>> \
    [[nodiscard]] constexpr auto OP(const malg::VectorBase<T, N> &l, const U &r)                                      \
    {                                                                                                                 \
        return malg::details::apply_binary<N>(FUNC{}, l, r);                                                          \
    }

DEFINE_UNARY_VECTOR_OPERATIONS(operator-, malg::details::negate)
DEFINE_BINARY_VECTOR_OPERATIONS(operator+, malg::details::plus)
DEFINE_BINARY_VECTOR_OPERATIONS(operator-, malg::details::minus)
DEFINE_BINARY_VECTOR_OPERATIONS(operator*, malg::details::multiplies)
DEFINE_BINARY_VECTOR_OPERATIONS(operator/, malg::details::divides)
DEFINE_BINARY_VECTOR_OPERATIONS(operator==, malg::details::equal)
DEFINE_BINARY_VECTOR_OPERATIONS(operator!=, malg::details::not_equal_to)
DEFINE_BINARY_VECTOR_OPERATIONS(operator<, malg::details::less)
DEFINE_BINARY_VECTOR_OPERATIONS(operator<=, malg::details::less_equal)
DEFINE_BINARY_VECTOR_OPERATIONS(operator>, malg::details::greater)
DEFINE_BINARY_VECTOR_OPERATIONS(operator>=, malg::details::greater_equal)
