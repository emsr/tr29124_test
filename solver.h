#ifndef SOLVER_H
#define SOLVER_H 1

#include <experimental/array>
#include <complex>
#include <variant>
#include <iosfwd>

namespace __gnu_cxx
{

  template<typename _Real>
    using solution_t
	    = std::variant<std::monostate, _Real, std::complex<_Real>>;

  //template<typename _Real>
  //  operator bool(const solution_t<_Real>& __x)
  //  { return __x.index() != 0; }

  template<typename _Real>
    solution_t<_Real>
    real(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return solution_t<_Real>{};
      else if (__x.index() == 1)
	return solution_t<_Real>(std::get<1>(__x));
      else
	return solution_t<_Real>(std::real(std::get<2>(__x)));
    }

  template<typename _Real>
    solution_t<_Real>
    imag(const solution_t<_Real>& __x)
    {
      if (__x.index() == 0)
	return solution_t<_Real>{};
      else if (__x.index() == 1)
	return solution_t<_Real>(_Real{0});
      else
	return solution_t<_Real>(std::imag(std::get<2>(__x)));
    }

  /**
   * Lexicographic order of solutions as complex numbers.
   * Null solutions compare as less than except to another null solution.
   */
  template<typename _Real>
    bool
    operator<(const solution_t<_Real>& __a, const solution_t<_Real>& __b)
    {
      if (__a.index() == 0 && __b.index() == 0)
	return false;
      else if (__a.index() == 0)
	return true;
      else if (__b.index() == 0)
	return false;
      else
	{
	  const auto __rea = real(__a);
	  const auto __reb = real(__b);
	  if (std::get<1>(__rea) < std::get<1>(__reb))
	    return true;
	  else if (__rea == __reb)
	    return std::get<1>(imag(__a)) < std::get<1>(imag(__b));
	  else
	    return false;
	}
    }


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

  template<typename _Char, typename _Real>
    std::basic_ostream<_Char>&
    operator<<(std::basic_ostream<_Char>& __out,
	       const __gnu_cxx::solution_t<_Real>& __sln)
    {
      const auto __idx = __sln.index();
      if (__idx == 0)
	__out << "null";
      else if (__idx == 1)
	__out << std::get<1>(__sln);
      else if (__idx == 2)
	__out << std::get<2>(__sln);
      return __out;
    }

#include "solver.tcc"


#endif  //  SOLVER_H
