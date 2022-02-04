/**
 *
 */

#include <cmath>
#include <initializer_list>
#include <algorithm>
#include <iostream>

template<typename Tp>
  constexpr bool
  lessabs(const Tp& a, const Tp& b)
  { return std::abs(a) < std::abs(b); }

template<typename Tp>
  constexpr Tp
  minabs(const Tp& a, const Tp& b)
  { return std::abs(std::min(a, b, lessabs<Tp>)); }

template<typename Tp>
  constexpr Tp
  minabs(std::initializer_list<Tp> il)
  { return std::abs(std::min(il, lessabs<Tp>)); }

template<typename Tp>
  constexpr Tp
  maxabs(const Tp& a, const Tp& b)
  { return std::abs(std::max(a, b, lessabs<Tp>)); }

template<typename Tp>
  constexpr Tp
  maxabs(std::initializer_list<Tp> il)
  { return std::abs(std::max(il, lessabs<Tp>)); }

template<typename Tp>
  constexpr std::pair<Tp, Tp>
  minmaxabs(const Tp& a, const Tp& b)
  {
    auto [c, d] = std::minmax(a, b, lessabs<Tp>);
    return {std::abs(c), std::abs(d)};
  }

template<typename Tp>
  constexpr std::pair<Tp, Tp>
  minmaxabs(std::initializer_list<Tp> il)
  {
    auto [c, d] = std::minmax(il, lessabs<Tp>);
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
