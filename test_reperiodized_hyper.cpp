/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -I. -o test_reperiodized_hyper test_reperiodized_hyper.cpp wrap_boost.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_reperiodized_hyper > test_reperiodized_hyper.txt

g++ -std=gnu++17 -g -DNO_SINH_COSH_PI -I. -o test_reperiodized_hyper test_reperiodized_hyper.cpp wrap_boost.cpp -lgsl -lgslcblas -lquadmath
./test_reperiodized_hyper > test_reperiodized_hyper.txt
*/

#include <iostream>
#include <iomanip>
#include <limits>

#include <ext/cmath>

template<typename _Tp>
  void
  run_sin_cosh_pi()
  {
    constexpr _Tp _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
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
  std::cout << "\ndouble\n=====\n\n";
  run_sin_cosh_pi<double>();

  std::cout << "\nlong double\n=====\n\n";
  run_sin_cosh_pi<long double>();
}
