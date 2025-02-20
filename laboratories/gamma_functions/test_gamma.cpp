/**
 *
 */

#include <limits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <complex>

#include <emsr/special_functions.h>
#include <emsr/fp_type_util.h>
#include <emsr/float128_io.h>

#include <wrap_boost.h>

template<typename Tp, typename _Gamma>
  void
  test_gamma(_Gamma gamma)
  {
    using Val = Tp;
    using Real = emsr::num_traits_t<Val>;

    std::cout.precision(std::numeric_limits<Real>::digits10);
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << '\n'
	      << ' ' << std::setw(width) << "a"
	      << ' ' << std::setw(width) << "log_gamma"
	      << ' ' << std::setw(width) << "std::lgamma"
	      << ' ' << std::setw(width) << "log_gamma"
	      << ' ' << std::setw(width) << "delta_std"
	      << ' ' << std::setw(width) << "delta_boost"
	      << '\n';
    int i_min = -200;
    for (int i = i_min; i <= +500; ++i)
      {
	auto a = Tp{0.10L} * i;
	auto gam = gamma(a - Tp{1});
	auto gam0 = std::lgamma(a);
	auto glgam = emsr::detail::log_gamma(a);
	auto blgam = beast::lgamma(a);
	std::cout << ' ' << std::setw(width) << a
		  << ' ' << std::setw(width) << gam
		  << ' ' << std::setw(width) << gam0
		  << ' ' << std::setw(width) << glgam
		  << ' ' << std::setw(width) << (gam - gam0) / std::abs(gam0)
		  << ' ' << std::setw(width) << (glgam - blgam) / std::abs(blgam)
		  << '\n';
      }
  }

int
main()
{

  std::cout << "\n\nLanczos Algorithm\n\n";
  test_gamma<double>(emsr::detail::lanczos_log_gamma1p<double>);

  std::cout << "\n\nSpouge Algorithm\n\n";
  test_gamma<double>(emsr::detail::spouge_log_gamma1p<double>);
}
