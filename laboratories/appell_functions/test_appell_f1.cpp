/**
 *
 */

#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>

template<typename Tp>
  Tp
  appell_f1_series(Tp alpha, Tp betax, Tp betay, Tp gamma, Tp x, Tp y)
  {
    const Tp _S_eps = std::numeric_limits<Tp>::epsilon();
    const int N = 100;
    auto f1 = Tp{0};
    auto termx = Tp{1};
    for (int m = 0; m < N; ++m)
      {
	auto termy = Tp{1};
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

template<typename Tp>
  Tp
  appell_f2_series(Tp alpha, Tp beta, Tp betap, Tp gamma, Tp gammap, Tp x, Tp y)
  {
  }

template<typename Tp>
  Tp
  appell_f3_series(Tp alpha, Tp alphap, Tp beta, Tp betap, Tp gamma, Tp x, Tp y)
  {
  }

template<typename Tp>
  Tp
  appell_f4_series(Tp alpha, Tp beta, Tp gamma, Tp gammap, Tp x, Tp y)
  {
  }

template<typename Tp>
  void
  test_appell_f1(Tp alpha, Tp betax, Tp betay, Tp gamma)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();

    const auto del = Tp{0.01};
    for (int m = -99; m <= 99; ++m)
      {
        auto x = m * del;
	for (int n = -99; n <= 99; ++n)
	  {
	    auto y = n * del;
	    auto f1 = appell_f1_series(alpha, betax, betay, gamma, x, y);
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

