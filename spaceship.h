#ifndef SPACESHIP_H
#define SPACESHIP_H 1

#include <type_traits>

/**
 * Implement the spaceship operator <==> that returns
 *  +1 if a < b
 *   0 if a == b
 *  -1 if a > b
 */
template<typename _Int>
  struct _Spaceship
  {
    using __type = _Int;

    constexpr int
    operator()(__type __a, __type __b)
    { return __a < __b ? +1 : (__b < __a ? -1 : 0); }
  };

/**
 * We need, not an operator, but a class with three specializations of the spaceship operator.
 */
template<typename _Int, int _Sign>
  struct _SpaceshipType : public std::integral_constant<int, _Sign>
  {};

template<typename _Int>
  using _SpaceLess = _SpaceshipType<_Int, +1>;
template<typename _Int>
  using _SpaceEqual = _SpaceshipType<_Int, 0>;
template<typename _Int>
  using _SpaceGreater = _SpaceshipType<_Int, -1>;

#endif // SPACESHIP_H
