/*
$HOME/bin/bin/g++ -std=c++17 -I. -o variant_help variant_help.cpp
./variant_help
*/

#include <complex>
#include <type_traits>
#include <iostream>
#include <variant>

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

int
main()
{
  std::cout << "is_move_constructible_v<double> = " << std::is_move_constructible_v<double> << '\n';
  std::cout << "is_swappable_v<double>          = " << std::is_swappable_v<double> << '\n';

  std::cout << "is_move_constructible_v<complex<double>> = " << std::is_move_constructible_v<std::complex<double>> << '\n';
  std::cout << "is_swappable_v<complex<double>>          = " << std::is_swappable_v<std::complex<double>> << '\n';

  std::cout << "is_move_constructible_v<monostate> = " << std::is_move_constructible_v<std::monostate> << '\n';
  std::cout << "is_swappable_v<monostate>          = " << std::is_swappable_v<std::monostate> << '\n';

  std::cout << "is_move_constructible_v<solution_t<double>> = " << std::is_move_constructible_v<solution_t<double>> << '\n';
  std::cout << "is_swappable_v<solution_t<double>>          = " << std::is_swappable_v<solution_t<double>> << '\n';
}
