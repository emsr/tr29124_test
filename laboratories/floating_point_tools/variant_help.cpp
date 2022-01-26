/**
 *
 */

#include <complex>
#include <type_traits>
#include <iostream>
#include <variant>

  template<typename _Real>
    using solution_t
	    = std::variant<std::monostate, _Real, std::complex<_Real>>;

  //template<typename _Real>
  //  operator bool(const solution_t<_Real>& x)
  //  { return x.index() != 0; }

  template<typename _Real>
    solution_t<_Real>
    real(const solution_t<_Real>& x)
    {
      if (x.index() == 0)
	return solution_t<_Real>{};
      else if (x.index() == 1)
	return solution_t<_Real>(std::get<1>(x));
      else
	return solution_t<_Real>(std::real(std::get<2>(x)));
    }

  template<typename _Real>
    solution_t<_Real>
    imag(const solution_t<_Real>& x)
    {
      if (x.index() == 0)
	return solution_t<_Real>{};
      else if (x.index() == 1)
	return solution_t<_Real>(_Real{0});
      else
	return solution_t<_Real>(std::imag(std::get<2>(x)));
    }

  /**
   * Lexicographic order of solutions as complex numbers.
   * Null solutions compare as less than except to another null solution.
   */
  template<typename _Real>
    bool
    operator<(const solution_t<_Real>& a, const solution_t<_Real>& b)
    {
      if (a.index() == 0 && b.index() == 0)
	return false;
      else if (a.index() == 0)
	return true;
      else if (b.index() == 0)
	return false;
      else
	{
	  const auto rea = real(a);
	  const auto reb = real(b);
	  if (std::get<1>(rea) < std::get<1>(reb))
	    return true;
	  else if (rea == reb)
	    return std::get<1>(imag(a)) < std::get<1>(imag(b));
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
