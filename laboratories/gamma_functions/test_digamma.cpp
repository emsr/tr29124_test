/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <wrap_gsl.h>

template<typename _Tp>
  void
  test_digamma(_Tp __proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(__proto));
    auto w = 8 + std::cout.precision();

    auto x_start = -9.9375L;
    auto x_stop = +10.0625L;
    const unsigned int max = 801;
    const auto delta = (x_stop - x_start) / (max - 1);

    std::cout << '\n' << '\n';
    for (unsigned int i = 0; i < max; ++i)
      {
	_Tp x = x_start + i * delta;
	_Tp y_gcc = __gnu_cxx::digamma(x);
	_Tp y_gsl = gsl::digamma(x);
	std::cout << std::setw(5) << x
                  << ' ' << std::setw(w) << y_gcc
                  << ' ' << std::setw(w) << y_gsl
                  << ' ' << std::setw(w) << y_gcc - y_gsl
                  << '\n';
      }

    std::cout << '\n' << '\n';
    for (unsigned int i = 1; i <= 200; ++i)
      {
	_Tp x = i * 0.25;
	_Tp y_gcc = __gnu_cxx::digamma(x);
	_Tp y_gsl = gsl::digamma(x);
	std::cout << std::setw(5) << x
                  << ' ' << std::setw(w) << y_gcc
                  << ' ' << std::setw(w) << y_gsl
                  << ' ' << std::setw(w) << y_gcc - y_gsl;
	if (i % 4 == 0)
	  std::cout << ' ' << std::setw(w) << __gnu_cxx::harmonic<_Tp>(i / 4);
	std::cout << '\n';
      }
  }

int
main()
{
  test_digamma<double>();

  return 0;
}


