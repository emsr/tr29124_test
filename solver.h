#ifndef SOLVER_H
#define SOLVER_H 1

#include <experimental/array>
#include <complex>
#include <variant>
#include <iosfwd>

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
	  if (__rea < __reb)
	    return true;
	  else if (__rea == __reb)
	    return imag(__a) < imag(__b);
	  else
	    return false;
	}
    }

  template<typename _Char, typename _Real>
    std::basic_ostream<_Char>&
    operator<<(std::basic_ostream<_Char>& __out,
	       const solution_t<_Real>& __sln)
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

  template<typename _Real>
    std::array<solution_t<_Real>, 2>
    quadratic(const std::array<_Real, 3>& coef);

  template<typename _Real>
    std::array<solution_t<_Real>, 3>
    cubic(const std::array<_Real, 4>& coef);

  template<typename _Real>
    std::array<solution_t<_Real>, 4>
    quartic(const std::array<_Real, 5>& coef);


#include "solver.tcc"


#endif  //  SOLVER_H
