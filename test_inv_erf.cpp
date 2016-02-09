// g++ -std=c++14 -o test_inv_erf test_inv_erf.cpp

// ./test_inv_erf

// AAOF pp. 408-409.

#include<vector>
#include<iostream>
#include<limits>
#include<cmath>

int
main()
{
  //  Build the series coefficients.
  using _Tp = double;
  const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
  const auto _S_pi =  _Tp{3.1415926535897932384626433832795029L};
  const auto _S_sqrt_pi = _Tp{1.772453850905516027298167483341145182797L};

  std::cout.precision(std::numeric_limits<_Tp>::digits10);

  std::vector<_Tp> a;
  a.push_back(1);
  for (int j = 1; j < 100; ++j)
    {
      auto atemp = _Tp{0};
      for (int k = 1; k <= j; ++k)
        atemp += a[k - 1] * a[j - k] / _Tp(k) / _Tp(2 * k - 1);
      atemp /= _Tp(2 * j + 1);
      a.push_back(atemp);
    }
  for (auto aa : a)
    std::cout << ' ' << aa << '\n';

  for (int i = -200; i <= 200; ++i)
    {
      auto x = i * 0.01;
      auto chi = _S_sqrt_pi * x / _Tp{2};
      auto chi2 = chi * chi;
      auto inverf = _Tp{0};
      auto chip = chi;
      for (int k = 0; k < 100; ++k)
	{
	  auto term = a[k] * chip;
	  inverf += term;
	  if (std::abs(term) < _S_eps)
	    break;
	  chip *= chi2;
	}
      std::cout << ' ' << x << ' ' << inverf << ' ' << erf(inverf) << '\n';
    }
}
