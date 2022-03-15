/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <emsr/float128_io.h>
#include <emsr/sf_trig.h>

template<typename Tp>
  void
  run_sin_cosh_pi(Tp proto = Tp{})
  {
    const Tp s_pi = emsr::pi_v<Tp>;

    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 10 + std::cout.precision();

    std::cout << '\n';
    std::cout << std::setw(width) << " x"
	      << std::setw(width) << " sinh_pi (GCC)"
	      << std::setw(width) << " delta sinh(pi x)"
	      << std::setw(width) << " cosh_pi (GCC)"
	      << std::setw(width) << " delta cosh(pi x)"
	      << std::setw(width) << " tanh_pi (GCC)"
	      << std::setw(width) << " delta tanh(pi x)"
	      << '\n';
    std::cout << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << std::setw(width) << " ==============="
	      << '\n';
    const auto del = Tp{1}/Tp{10};
    for (int i = -1600; i <= +1600; ++i)
      {
	auto x = del * i;
	auto sinh_pi_g = emsr::sinh_pi(x);
	auto cosh_pi_g = emsr::cosh_pi(x);
	auto tanh_pi_g = emsr::tanh_pi(x);
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << sinh_pi_g
		  << ' ' << std::setw(width) << sinh_pi_g - std::sinh(s_pi * x)
		  << ' ' << std::setw(width) << cosh_pi_g
		  << ' ' << std::setw(width) << cosh_pi_g - std::cosh(s_pi * x)
		  << ' ' << std::setw(width) << tanh_pi_g
		  << ' ' << std::setw(width) << tanh_pi_g - std::tanh(s_pi * x)
		  << '\n';
      }
    std::cout << '\n';
  }


int
main()
{
  std::cout << "\nfloat\n=====\n\n";
  run_sin_cosh_pi<float>();

  std::cout << "\ndouble\n=====\n\n";
  run_sin_cosh_pi<double>();

  std::cout << "\nlong double\n=====\n\n";
  run_sin_cosh_pi<long double>();

  //std::cout << "\n__float128\n=====\n\n";
  //run_sin_cosh_pi<__float128>();
}
