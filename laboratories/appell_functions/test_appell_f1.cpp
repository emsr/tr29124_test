/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

template<typename _Tp>
  _Tp
  __appell_f1_series(_Tp alpha, _Tp betax, _Tp betay, _Tp gamma, _Tp x, _Tp y)
  {
    const _Tp _S_eps = std::numeric_limits<_Tp>::epsilon();
    const int N = 100;
    auto f1 = _Tp{0};
    auto termx = _Tp{1};
    for (int m = 0; m < N; ++m)
      {
	auto termy = _Tp{1};
	for (int n = 0; n < N; ++n)
	  {
	    const auto term = termx * termy;
	    f1 += term;
	    if (std::abs(term) < _S_eps * std::abs(f1))
	      break;
	    termy *= (alpha + n) * (betay + n) * y / (gamma + n) / (n + 1);
	  }
	termx *= (alpha + m) * (betax + m) * x / (gamma + m) / (m + 1);
      }
    return f1;
  }

template<typename _Tp>
  _Tp
  __appell_f2_series(_Tp alpha, _Tp beta, _Tp betap, _Tp gamma, _Tp gammap, _Tp x, _Tp y)
  {
  }

template<typename _Tp>
  _Tp
  __appell_f3_series(_Tp alpha, _Tp alphap, _Tp beta, _Tp betap, _Tp gamma, _Tp x, _Tp y)
  {
  }

template<typename _Tp>
  _Tp
  __appell_f4_series(_Tp alpha, _Tp beta, _Tp gamma, _Tp gammap, _Tp x, _Tp y)
  {
  }

template<typename _Tp>
  void
  test_appell_f1(_Tp alpha, _Tp betax, _Tp betay, _Tp gamma)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 8 + std::cout.precision();

    const auto del = _Tp{0.01};
    for (int m = -99; m <= 99; ++m)
      {
        auto x = m * del;
	for (int n = -99; n <= 99; ++n)
	  {
	    auto y = n * del;
	    auto f1 = __appell_f1_series(alpha, betax, betay, gamma, x, y);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << y
		      << ' ' << std::setw(w) << f1
		      << '\n';
	  }
	std::cout << '\n';
      }
  }

int
main()
{
  test_appell_f1(1.0, 1.5, -0.5, 0.5);
}

