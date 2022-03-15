/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>

#include <emsr/math_constants.h>
#include <emsr/numeric_limits.h>

/**
 * Integer powers of 2.
 */
template<typename Tp>
  constexpr Tp
  pow2(int n)
  {
    if (n == 0)
      return Tp{1};
    else
      {
	const long p = 0x1L << std::abs(n);
	return n > 0 ? Tp(p) : Tp{1} / Tp(p);
      }
  }

/**
 * 
 */
template<typename Tp>
  constexpr Tp
  exp_bailey(Tp x, int _J = 8)
  {
    constexpr auto s_N = (std::numeric_limits<Tp>::max_digits10 + 1) / 2;
    constexpr auto s_log_2 = emsr::ln2_v<Tp>;
    if (_J <= 0)
      _J = s_N - 1;
    const int n = x / s_log_2;
    const auto r = (x - n * s_log_2) / pow2<Tp>(_J);
    auto term = Tp{1};
    auto expr = term;
    for (int k = 1; k < s_N; ++k)
      expr += (term *= r / Tp(k));
    for (int i = 0; i < _J; ++i)
      expr *= expr;
    return expr * pow2<Tp>(n);
  }

template<typename Tp>
  void
  test_exp(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = 8 + std::cout.precision();

    for (int i = -100; i <= +100; ++i)
      {
	Tp x = Tp{0.1L} * i;
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << exp_bailey(x)
		  << ' ' << std::setw(w) << exp_bailey(x, 0)
		  << ' ' << std::setw(w) << std::exp(x)
		  << '\n';
      }
  }

int
main()
{
  std::cout << "\nTesting float...\n";
  test_exp(1.0F);

  std::cout << "\nTesting double...\n";
  test_exp(1.0);

  std::cout << "\nTesting long double...\n";
  test_exp(1.0L);
}
