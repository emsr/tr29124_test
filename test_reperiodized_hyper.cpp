/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_reperiodized_hyper test_reperiodized_hyper.cpp wrap_boost.cpp -lquadmath
./test_reperiodized_hyper > test_reperiodized_hyper.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -DNO_SINH_COSH_PI -I. -o test_reperiodized_hyper test_reperiodized_hyper.cpp wrap_boost.cpp -lgsl -lgslcblas -lquadmath
./test_reperiodized_hyper > test_reperiodized_hyper.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>
#include <bits/float128_io.h>

#include <ext/cmath>

template<typename _Tp>
  void
  run_sin_cosh_pi(_Tp proto = _Tp{})
  {
    const _Tp _S_pi = __gnu_cxx::__const_pi<_Tp>(proto);

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 10 + std::cout.precision();

    std::cout << std::endl;
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
    for (int i = -1600; i <= +1600; ++i)
      {
	auto x = _Tp(0.1Q * i);
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
    std::cout << std::endl;
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

  std::cout << "\n__float128\n=====\n\n";
  run_sin_cosh_pi<__float128>();
}
