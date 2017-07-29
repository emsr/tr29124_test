#ifndef SOLVER_H
#define SOLVER_H 1

#include <experimental/array>
#include <complex>
#include <variant>
#include <iosfwd>

  template<typename _Real>
    using solution_t
	    = std::variant<std::monostate, _Real, std::complex<_Real>>;

  template<typename _Char, typename _Real>
    std::basic_ostream<_Char>&
    operator<<(std::basic_ostream<_Char> &__out,
	       const solution_t<_Real> &__sln)
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
    quadratic(const std::array<_Real, 3>& asub);

  template<typename _Real>
    std::array<solution_t<_Real>, 3>
    cubic(const std::array<_Real, 4>& asub);

  template<typename _Real>
    std::array<solution_t<_Real>, 4>
    quartic(const std::array<_Real, 5>& asub);


#include "solver.tcc"


#endif  //  SOLVER_H
