/**
 *
 */

#include <limits>
#include <iostream>
#include <iomanip>
#include <cmath>

#include <emsr/float128_io.h>
#include <emsr/sf_bessel.h>

#include <wrap_boost.h>

template<typename Tp>
  void
  RunSphHankel1(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10<Tp>(proto));
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    auto max_max_abs_frac = Tp{-1};
    for (int inu = 0; inu <= +500; ++inu)
      {
	auto max_abs_frac = Tp{-1};
	auto nu = inu * 0.1;
	std::cout << '\n';
	for (int iz = 1; iz <= +1000; ++iz)
	  {
	    auto z = iz * 0.01;
	    auto h1s = emsr::sph_hankel_1(nu, z);
	    auto h1b = beast::sph_hankel_1(nu, z);
	    auto abs_frac = std::abs((h1s - h1b) / h1b);
	    std::cout << ' ' << std::setw(width) << nu
		      << ' ' << std::setw(width) << z
		      << ' ' << std::setw(2 * width) << h1s
		      << ' ' << std::setw(2 * width) << h1b
		      << ' ' << std::setw(width) << abs_frac
		      << '\n';
	    if (abs_frac > max_abs_frac)
	      max_abs_frac = abs_frac;
	  }
	std::cout << " max_abs_frac = " << std::setw(width) << max_abs_frac << '\n';
	if (max_abs_frac > max_max_abs_frac)
	  max_max_abs_frac = max_abs_frac;
      }
    std::cout << " max_max_abs_frac = " << std::setw(width) << max_max_abs_frac << '\n';
  }

int
main()
{
  RunSphHankel1<float>();
  RunSphHankel1<double>();
  RunSphHankel1<long double>();
#ifdef EMSR_HAVE_FLOAT128
  RunSphHankel1<__float128>();
#endif
}
