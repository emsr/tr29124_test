#ifndef SOLVER_LOW_DEGREE_H
#define SOLVER_LOW_DEGREE_H 1

#include <experimental/array>

#include "solution.h"

namespace __gnu_cxx
{


  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 2>
    __quadratic(const _Iter& __coef);

  template<typename _Real>
    inline std::array<solution_t<_Real>, 2>
    __quadratic(_Real __c0, _Real __c1, _Real __c2)
    {
      using std::experimental::make_array;
      return __quadratic<_Real>(make_array(__c0, __c1, __c2));
    }


  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 3>
    __cubic(const _Iter& __coef);

  template<typename _Real>
    inline std::array<solution_t<_Real>, 3>
    __cubic(_Real __c0, _Real __c1, _Real __c2, _Real __c3)
    {
      using std::experimental::make_array;
      return __cubic<_Real>(make_array(__c0, __c1, __c2, __c3));
    }


  template<typename _Real, typename _Iter>
    std::array<solution_t<_Real>, 4>
    __quartic(const _Iter& __coef);

  template<typename _Real>
    inline std::array<solution_t<_Real>, 4>
    __quartic(_Real __c0, _Real __c1, _Real __c2, _Real __c3, _Real __c4)
    {
      using std::experimental::make_array;
      return __quartic<_Real>(make_array(__c0, __c1, __c2, __c3, __c4));
    }


} // namespace __gnu_cxx


#include "solver_low_degree.tcc"


#endif // SOLVER_LOW_DEGREE_H
