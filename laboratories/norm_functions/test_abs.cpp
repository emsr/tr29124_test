/**
 *
 */

#include <cmath>
#include <initializer_list>
#include <algorithm>
#include <iostream>

template<typename _Tp>
  constexpr bool
  lessabs(const _Tp& __a, const _Tp& __b)
  { return std::abs(__a) < std::abs(__b); }

template<typename _Tp>
  constexpr _Tp
  minabs(const _Tp& __a, const _Tp& __b)
  { return std::abs(std::min(__a, __b, lessabs<_Tp>)); }

template<typename _Tp>
  constexpr _Tp
  minabs(std::initializer_list<_Tp> __il)
  { return std::abs(std::min(__il, lessabs<_Tp>)); }

template<typename _Tp>
  constexpr _Tp
  maxabs(const _Tp& __a, const _Tp& __b)
  { return std::abs(std::max(__a, __b, lessabs<_Tp>)); }

template<typename _Tp>
  constexpr _Tp
  maxabs(std::initializer_list<_Tp> __il)
  { return std::abs(std::max(__il, lessabs<_Tp>)); }

template<typename _Tp>
  constexpr std::pair<_Tp, _Tp>
  minmaxabs(const _Tp& __a, const _Tp& __b)
  {
    auto [__c, __d] = std::minmax(__a, __b, lessabs<_Tp>);
    return {std::abs(__c), std::abs(__d)};
  }

template<typename _Tp>
  constexpr std::pair<_Tp, _Tp>
  minmaxabs(std::initializer_list<_Tp> __il)
  {
    auto [__c, __d] = std::minmax(__il, lessabs<_Tp>);
    return {std::abs(__c), std::abs(__d)};
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
