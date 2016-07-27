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
    using __stype = std::make_signed_t<_Int>;

    constexpr __stype
    operator()(__type __a, __type __b)
    {
      return __a < __b ? __stype(-1)
		       : (__b < __a ? __stype(+1)
				    : __stype(0));
    }
  };

template<typename _Int>
  using _Spaceship_t = typename _Spaceship<_Int>::__stype;

/**
 * We need, not an operator, but a class with three specializations of the spaceship operator.
 */
template<typename _Int, _Spaceship_t<_Int> _Sign>
  struct _SpaceshipType
  {
  };

#endif // SPACESHIP_H
