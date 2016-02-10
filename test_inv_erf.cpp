// $HOME/bin_specfun/bin/g++ -D__STDCPP_WANT_MATH_SPEC_FUNCS__ -o test_inv_erf test_inv_erf.cpp -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_inv_erf > test_inv_erf.txt

// AAOF pp. 408-409.

#include <vector>
#include <iostream>
#include <limits>
#include <cmath>
#include "float128.h"

template<typename _Tp>
  void
  test_inv_erf()
  {
    //  Build the series coefficients.
    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_pi =  _Tp{3.1415926535897932384626433832795029Q};
    const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797Q};

    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto width = 6 + std::cout.precision();

    const int n_max = 250;
    std::vector<_Tp> a;
    a.push_back(1);
    for (int n = 1; n < n_max; ++n)
      {
	auto atemp = _Tp{0};
	for (int k = 1; k <= n; ++k)
          atemp += _Tp(2 * (k - 1) + 1) * a[k - 1]
		 * _Tp(2 * (n - k) + 1) * a[n - k]
		 / _Tp(k * (2 * k - 1));
	atemp /= _Tp(2 * n + 1);
	a.push_back(atemp);
      }
    for (auto aa : a)
      std::cout << ' ' << aa << '\n';

    std::cout << "\n\n"
	      << std::setw(width) << "x"
	      << std::setw(width) << "inv_erf(x)"
	      << std::setw(width) << "erf(inv_erf(x))"
	      << std::setw(width) << "erf(inv_erf(x)) - x"
	      << '\n';
    for (int i = -200; i <= 200; ++i)
      {
	auto x = i * _Tp{0.01Q};
	auto chi = _S_sqrt_pi * x / _Tp{2};
	auto chi2 = chi * chi;
	auto chip = chi;
	auto inverf = _Tp{0};
	for (int k = 0; k < n_max; ++k)
	  {
	    auto term = a[k] * chip;
	    inverf += term;
	    if (std::abs(term) / std::abs(inverf) < _S_eps)
	      break;
	    chip *= chi2;
	  }
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << inverf
		  << ' ' << std::setw(width) << std::erf(inverf)
		  << ' ' << std::setw(width) << std::erf(inverf) - x
		  << '\n';
      }

    std::cout << "\n\n"
	      << std::setw(width) << "x"
	      << std::setw(width) << "erf(x)"
	      << std::setw(width) << "inv_erf(erf(x))"
	      << std::setw(width) << "inv_erf(erf(x)) - x"
	      << '\n';
    for (int i = -200; i <= 200; ++i)
      {
	auto x = i * _Tp{0.01Q};
	auto erf = std::erf(x);
	auto chi = _S_sqrt_pi * erf / _Tp{2};
	auto chi2 = chi * chi;
	auto chip = chi;
	auto inverf = _Tp{0};
	for (int k = 0; k < n_max; ++k)
	  {
	    auto term = a[k] * chip;
	    inverf += term;
	    if (std::abs(term) / std::abs(inverf) < _S_eps)
	      break;
	    chip *= chi2;
	  }
	std::cout << ' ' << std::setw(width) << x
		  << ' ' << std::setw(width) << erf
		  << ' ' << inverf
		  << ' ' << inverf - x
		  << '\n';
      }
  }

int
main()
{
  std::cout << "\n\n  float\n";
  std::cout << "  =====\n";
  test_inv_erf<float>();

  std::cout << "\n\n  double\n";
  std::cout << "  ======\n";
  test_inv_erf<double>();

  std::cout << "\n\n  long double\n";
  std::cout << "  ===========\n";
  test_inv_erf<long double>();

  std::cout << "\n\n  __float128\n";
  std::cout << "  ==========\n";
  test_inv_erf<__float128>();
}
