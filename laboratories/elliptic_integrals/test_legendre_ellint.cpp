/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>

#include <emsr/numeric_limits.h>
#include <emsr/math_constants.h>
#include <emsr/sf_ellint.h>

#include <wrap_gsl.h>

template<typename _Tp>
  void
  test_comp_ellint_1()
  {
    std::cout.precision(emsr::digits10<_Tp>());
    const auto w = 8 + std::cout.precision();

    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i < 100; ++i)
      {
	_Tp k = i * (0.01L);
	const auto K_gcc = emsr::comp_ellint_1(k);
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
    std::cout.precision(emsr::digits10<_Tp>());
    const auto w = 8 + std::cout.precision();

    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i < 100; ++i)
      {
	_Tp k = i * (0.01L);
	const auto E_gcc = emsr::comp_ellint_2(k);
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
    std::cout.precision(emsr::digits10<_Tp>());
    const auto w = 8 + std::cout.precision();
    const auto s_pi_2 = emsr::pi_v<_Tp> / _Tp{2};

    const std::vector<_Tp> ks({0.0L, 0.5L, 0.75L, 0.9L, 0.95L, 0.99L});
    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i <= 90; ++i)
      {
	const auto phi = i * s_pi_2 / _Tp{90};
	std::cout << ' ' << std::setw(w) << phi;
	for (const auto k : ks)
	  {
	    const auto F_gcc = emsr::ellint_1(k, phi);
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
    std::cout.precision(emsr::digits10<_Tp>());
    const auto w = 8 + std::cout.precision();
    const auto s_pi_2 = emsr::pi_v<_Tp> / _Tp{2};

    const std::vector<_Tp> ks({0.0L, 0.5L, 0.75L, 0.9L, 0.95L, 0.99L});
    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i <= 90; ++i)
      {
	const auto phi = i * s_pi_2 / _Tp{90};
	std::cout << ' ' << std::setw(w) << phi;
	for (const auto k : ks)
	  {
	    const auto E_gcc = emsr::ellint_2(k, phi);
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
