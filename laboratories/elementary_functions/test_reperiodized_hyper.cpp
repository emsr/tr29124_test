/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

#include <ext/float128_io.h>

template<typename _Tp>
  void
  run_sin_cosh_pi(_Tp proto = _Tp{})
  {
    const _Tp _S_pi = __gnu_cxx::math::__pi_v<_Tp>;

    std::cout.precision(__gnu_cxx::__digits10(proto));
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
    const auto del = _Tp{1}/_Tp{10};
    for (int i = -1600; i <= +1600; ++i)
      {
	auto x = del * i;
	auto sinh_pi_g = __gnu_cxx::sinh_pi(x);
	auto cosh_pi_g = __gnu_cxx::cosh_pi(x);
	auto tanh_pi_g = __gnu_cxx::tanh_pi(x);
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << sinh_pi_g
		  << ' ' << std::setw(width) << sinh_pi_g - std::sinh(_S_pi * x)
		  << ' ' << std::setw(width) << cosh_pi_g
		  << ' ' << std::setw(width) << cosh_pi_g - std::cosh(_S_pi * x)
		  << ' ' << std::setw(width) << tanh_pi_g
		  << ' ' << std::setw(width) << tanh_pi_g - std::tanh(_S_pi * x)
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
