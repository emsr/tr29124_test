#ifndef SOLVER_H
#define SOLVER_H 1

#include <experimental/array>
#include <complex>
#include <variant>

  template<typename _Real>
    using solution_t
	    = std::variant<std::monostate, _Real, std::complex<_Real>>;

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
