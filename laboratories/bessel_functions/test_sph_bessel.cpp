/**
 *
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>
#include <stdexcept>
#include <limits>

#include <emsr/float128_io.h>

template<typename _Tp>
  void
  test_sph_bessel(_Tp proto = _Tp{})
  {
    const auto s_pi = emsr::pi_v<_Tp>;

    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    auto fname = [](std::string_view front, int n, std::string_view back)
		 {
		   std::ostringstream out;
		   out << front << n << back;
		   return out.str();
		 };

    for (int n = 0; n <= 50; ++n)
      {
	std::cout << "\n\n  n = " << n << '\n';
	std::cout << ' ' << std::setw(w) << "x";
	std::cout << ' ' << std::setw(w) << fname("j_", n, "(x)");
	std::cout << ' ' << std::setw(w) << fname("j'_", n, "(x)");
	std::cout << ' ' << std::setw(w) << fname("n_", n, "(x)");
	std::cout << ' ' << std::setw(w) << fname("n'_", n, "(x)");
	std::cout << ' ' << std::setw(w) << "x^2 W[j,n]";
	std::cout << ' ' << std::setw(w) << fname("i_", n, "(x)");
	std::cout << ' ' << std::setw(w) << fname("i'_", n, "(x)");
	std::cout << ' ' << std::setw(w) << fname("k'_", n, "(x)");
	std::cout << ' ' << std::setw(w) << fname("k_", n, "(x)");
	std::cout << ' ' << std::setw(w) << "-(2x^2/pi) W[i,k]";
	std::cout << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (int i = 0; i <= 100; ++i)
	  {
	    auto x = i * del;
	    auto Wjn = x * x;
	    auto Wik = -_Tp{2} * Wjn / s_pi;
	    std::cout << ' ' << std::setw(w) << x;
	    try
	      {
		auto Sph = emsr::detail::sph_bessel_jn(n, x);
		std::cout << ' ' << std::setw(w) << Sph.j_value;
		std::cout << ' ' << std::setw(w) << Sph.j_deriv;
		std::cout << ' ' << std::setw(w) << Sph.n_value;
		std::cout << ' ' << std::setw(w) << Sph.n_deriv;
		std::cout << ' ' << std::setw(w) << Wjn * Sph.Wronskian();
	      }
	    catch (std::exception& err)
	      {
		std::cerr << err.what() << '\n';
	      }
	    try
	      {
		auto Sph = emsr::detail::sph_bessel_ik(n, x);
		std::cout << ' ' << std::setw(w) << Sph.i_value;
		std::cout << ' ' << std::setw(w) << Sph.i_deriv;
		std::cout << ' ' << std::setw(w) << Sph.k_value;
		std::cout << ' ' << std::setw(w) << Sph.k_deriv;
		std::cout << ' ' << std::setw(w) << Wik * Sph.Wronskian();
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
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto w = 8 + std::cout.precision();

    auto fname = [](std::string_view front, int n, std::string_view back)
		 {
		   std::ostringstream out;
		   out << front << n << back;
		   return out.str();
		 };

    for (int n = 0; n <= 50; ++n)
      {
	std::cout << "\n\n  n = " << n << '\n';
	std::cout << ' ' << std::setw(w) << "x";
	std::cout << ' ' << std::setw(w) << fname("j_", n, "(x)");
	std::cout << ' ' << std::setw(w) << fname("n_", n, "(x)");
	const auto del = _Tp{1} / _Tp{10};
	for (int i = 0; i <= 100; ++i)
	  {
	    auto x = i * del;
	    std::cout << ' ' << std::setw(w) << x;
	    try
	      {
		auto jx = emsr::sph_bessel(n, x);
		auto nx = emsr::sph_neumann(n, x);
		std::cout << ' ' << std::setw(w) << jx;
		std::cout << ' ' << std::setw(w) << nx;
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

