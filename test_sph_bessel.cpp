/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -I. -o test_sph_bessel test_sph_bessel.cpp wrap_boost.cpp -lquadmath
./test_sph_bessel > test_sph_bessel.txt

$HOME/bin/bin/g++ -std=gnu++17 -I. -o test_sph_bessel test_sph_bessel.cpp wrap_boost.cpp -lquadmath
./test_sph_bessel > test_sph_bessel.txt

g++ -std=gnu++17 -DNO_LOGBQ -I. -o test_sph_bessel test_sph_bessel.cpp wrap_boost.cpp -lquadmath
./test_sph_bessel > test_sph_bessel.txt
*/

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <limits>

template<typename _Tp>
  void
  test_sph_bessel(_Tp proto = _Tp{})
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    const auto _S_pi = __gnu_cxx::__const_pi<_Real>();

    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    auto fname = [](std::string_view front, int n, std::string_view back)
		 {
		   std::ostringstream out;
		   out << front << n << back;
		   return out.str();
		 };

    std::cout << '\n';
    for (int n = 0; n <= 50; ++n)
      {
	std::cout << ' ' << std::setw(width) << "x";
	std::cout << ' ' << std::setw(width) << fname("j_", n, "(x)");
	std::cout << ' ' << std::setw(width) << fname("j'_", n, "(x)");
	std::cout << ' ' << std::setw(width) << fname("n_", n, "(x)");
	std::cout << ' ' << std::setw(width) << fname("n'_", n, "(x)");
	std::cout << ' ' << std::setw(width) << "x^2 W[j,n]";
	std::cout << ' ' << std::setw(width) << fname("i_", n, "(x)");
	std::cout << ' ' << std::setw(width) << fname("i'_", n, "(x)");
	std::cout << ' ' << std::setw(width) << fname("k'_", n, "(x)");
	std::cout << ' ' << std::setw(width) << fname("k_", n, "(x)");
	std::cout << ' ' << std::setw(width) << "-(2x^2/pi) W[i,k]";
	std::cout << '\n';
	for (int i = 0; i <= 100; ++i)
	  {
	    auto x = i * _Tp{0.1Q};
	    auto Wjn = x * x;
	    auto Wik = -_Tp{2} * Wjn / _S_pi;
	    std::cout << ' ' << std::setw(width) << x;
	    try
	      {
		auto Sph = std::__detail::__sph_bessel_jn(n, x);
		std::cout << ' ' << std::setw(width) << Sph.__j_value;
		std::cout << ' ' << std::setw(width) << Sph.__j_deriv;
		std::cout << ' ' << std::setw(width) << Sph.__n_value;
		std::cout << ' ' << std::setw(width) << Sph.__n_deriv;
		std::cout << ' ' << std::setw(width) << Wjn * Sph.__Wronskian();
	      }
	    catch (std::exception& err)
	      {
		std::cerr << err.what() << '\n';
	      }
	    try
	      {
		auto Sph = std::__detail::__sph_bessel_ik(n, x);
		std::cout << ' ' << std::setw(width) << Sph.__i_value;
		std::cout << ' ' << std::setw(width) << Sph.__i_deriv;
		std::cout << ' ' << std::setw(width) << Sph.__k_value;
		std::cout << ' ' << std::setw(width) << Sph.__k_deriv;
		std::cout << ' ' << std::setw(width) << Wik * Sph.__Wronskian();
	      }
	    catch (std::exception& err)
	      {
		std::cerr << err.what() << '\n';
	      }
	    std::cout << '\n';
	  }
      }
  }

template<typename _Tp>
  void
  test_std_bessel(_Tp proto = _Tp{})
  {
    using _Val = _Tp;
    using _Real = std::__detail::__num_traits_t<_Val>;
    const auto _S_pi = __gnu_cxx::__const_pi<_Real>();

    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    auto fname = [](std::string_view front, int n, std::string_view back)
		 {
		   std::ostringstream out;
		   out << front << n << back;
		   return out.str();
		 };

    std::cout << '\n';
    for (int n = 0; n <= 50; ++n)
      {
	std::cout << ' ' << std::setw(width) << "x";
	std::cout << ' ' << std::setw(width) << fname("j_", n, "(x)");
	std::cout << ' ' << std::setw(width) << fname("n_", n, "(x)");
	for (int i = 0; i <= 100; ++i)
	  {
	    auto x = i * _Tp{0.1Q};
	    std::cout << ' ' << std::setw(width) << x;
	    try
	      {
		auto jx = std::sph_bessel(n, x);
		auto nx = std::sph_neumann(n, x);
		std::cout << ' ' << std::setw(width) << jx;
		std::cout << ' ' << std::setw(width) << nx;
	      }
	    catch (std::exception& err)
	      {
		std::cerr << err.what() << '\n';
	      }
	    std::cout << '\n';
	  }
      }
  }

int
main()
{
  test_sph_bessel<float>();
  test_sph_bessel<double>();
  test_sph_bessel<long double>();
  //test_sph_bessel<__float128>();

  test_std_bessel<float>();
  test_std_bessel<double>();
  test_std_bessel<long double>();
  //test_std_bessel<__float128>();
}

