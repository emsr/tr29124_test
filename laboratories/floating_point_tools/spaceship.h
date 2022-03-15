#ifndef SPACESHIP_H
#define SPACESHIP_H 1

#include <type_traits>

/**
 * Implement the spaceship operator <==> that returns
 *  +1 if a < b
 *   0 if a == b
 *  -1 if a > b
 */
template<typename Int>
  struct Spaceship
  {
    using type = Int;

    constexpr int
    operator()(type a, type b)
    { return a < b ? +1 : (b < a ? -1 : 0); }
  };

/**
 * We need, not an operator, but a class with three specializations of the spaceship operator.
 */
template<typename Int, int _Sign>
  struct SpaceshipType : public std::integral_constant<int, _Sign>
  {};

template<typename Int>
  using SpaceLess = SpaceshipType<Int, +1>;
template<typename Int>
  using SpaceEqual = SpaceshipType<Int, 0>;
template<typename Int>
  using SpaceGreater = SpaceshipType<Int, -1>;

#endif // SPACESHIP_H
