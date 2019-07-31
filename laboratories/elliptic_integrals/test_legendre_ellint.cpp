/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <wrap_gsl.h>

template<typename _Tp>
  void
  test_comp_ellint_1()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    const auto w = 8 + std::cout.precision();

    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i < 100; ++i)
      {
	_Tp k = i * (0.01L);
	const auto K_gcc = std::comp_ellint_1(k);
	const auto K_gsl = gsl::comp_ellint_1(k);
	std::cout << std::setw(5) << k
                  << ' ' << std::setw(w) << K_gcc
                  << ' ' << std::setw(w) << K_gsl
                  << ' ' << std::setw(w) << K_gcc - K_gsl
                  << '\n';
      }
  }

template<typename _Tp>
  void
  test_comp_ellint_2()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    const auto w = 8 + std::cout.precision();

    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i < 100; ++i)
      {
	_Tp k = i * (0.01L);
	const auto E_gcc = std::comp_ellint_2(k);
	const auto E_gsl = gsl::comp_ellint_2(k);
	std::cout << std::setw(5) << k
                  << ' ' << std::setw(w) << E_gcc
                  << ' ' << std::setw(w) << E_gsl
                  << ' ' << std::setw(w) << E_gcc - E_gsl
                  << '\n';
      }
  }

template<typename _Tp>
  void
  test_ellint_1()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    const auto w = 8 + std::cout.precision();
    const auto _S_pi_2 = __gnu_cxx::numbers::__pi_half_v<_Tp>;

    const std::vector<_Tp> ks({0.0L, 0.5L, 0.75L, 0.9L, 0.95L, 0.99L});
    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i <= 90; ++i)
      {
	const auto phi = i * _S_pi_2 / _Tp{90};
	std::cout << ' ' << std::setw(w) << phi;
	for (const auto k : ks)
	  {
	    const auto F_gcc = std::ellint_1(k, phi);
	    const auto F_gsl = gsl::ellint_1(k, phi);
	    std::cout << ' ' << std::setw(w) << F_gcc
                      << ' ' << std::setw(w) << F_gsl
                      << ' ' << std::setw(w) << F_gcc - F_gsl;
	  }
	std::cout << '\n';
      }
  }

template<typename _Tp>
  void
  test_ellint_2()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    const auto w = 8 + std::cout.precision();
    const auto _S_pi_2 = __gnu_cxx::numbers::__pi_half_v<_Tp>;

    const std::vector<_Tp> ks({0.0L, 0.5L, 0.75L, 0.9L, 0.95L, 0.99L});
    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i <= 90; ++i)
      {
	const auto phi = i * _S_pi_2 / _Tp{90};
	std::cout << ' ' << std::setw(w) << phi;
	for (const auto k : ks)
	  {
	    const auto E_gcc = std::ellint_2(k, phi);
	    const auto E_gsl = gsl::ellint_2(k, phi);
	    std::cout << ' ' << std::setw(w) << E_gcc
                      << ' ' << std::setw(w) << E_gsl
                      << ' ' << std::setw(w) << E_gcc - E_gsl;
	  }
	std::cout << '\n';
      }
  }

int
main()
{
  test_comp_ellint_1<double>();
  test_ellint_1<double>();

  test_comp_ellint_2<double>();
  test_ellint_2<double>();
}
