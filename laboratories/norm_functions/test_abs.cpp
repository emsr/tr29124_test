/**
 *
 */

#include <cmath>
#include <initializer_list>
#include <algorithm>
#include <iostream>

template<typename _Tp>
  constexpr bool
  lessabs(const _Tp& a, const _Tp& b)
  { return std::abs(a) < std::abs(b); }

template<typename _Tp>
  constexpr _Tp
  minabs(const _Tp& a, const _Tp& b)
  { return std::abs(std::min(a, b, lessabs<_Tp>)); }

template<typename _Tp>
  constexpr _Tp
  minabs(std::initializer_list<_Tp> il)
  { return std::abs(std::min(il, lessabs<_Tp>)); }

template<typename _Tp>
  constexpr _Tp
  maxabs(const _Tp& a, const _Tp& b)
  { return std::abs(std::max(a, b, lessabs<_Tp>)); }

template<typename _Tp>
  constexpr _Tp
  maxabs(std::initializer_list<_Tp> il)
  { return std::abs(std::max(il, lessabs<_Tp>)); }

template<typename _Tp>
  constexpr std::pair<_Tp, _Tp>
  minmaxabs(const _Tp& a, const _Tp& b)
  {
    auto [c, d] = std::minmax(a, b, lessabs<_Tp>);
    return {std::abs(c), std::abs(d)};
  }

template<typename _Tp>
  constexpr std::pair<_Tp, _Tp>
  minmaxabs(std::initializer_list<_Tp> il)
  {
    auto [c, d] = std::minmax(il, lessabs<_Tp>);
    return {std::abs(c), std::abs(d)};
  }

int
main()
{
  auto a = minabs({-4, 3, 2});
  std::cout << "minabs({-4, 3, 2}) = " << a << '\n';

  auto b = maxabs({-4, 3, 2});
  std::cout << "maxabs({-4, 3, 2}) = " << b << '\n';

  auto [c, d] = minmaxabs({-4, 3, 2});
  std::cout << "minmaxabs({-4, 3, 2}) = " << c << ' ' << d << '\n';
}
