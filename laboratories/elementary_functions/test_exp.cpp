/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>

/**
 * Integer powers of 2.
 */
template<typename _Tp>
  constexpr _Tp
  pow2(int n)
  {
    if (n == 0)
      return _Tp{1};
    else
      {
	const long p = 0x1L << std::abs(n);
	return n > 0 ? _Tp(p) : _Tp{1} / _Tp(p);
      }
  }

/**
 * 
 */
template<typename _Tp>
  constexpr _Tp
  exp_bailey(_Tp x, int _J = 8)
  {
    constexpr auto _S_N = (std::numeric_limits<_Tp>::max_digits10 + 1) / 2;
    constexpr auto _S_log_2 = emsr::ln2_v<_Tp>;
    if (_J <= 0)
      _J = _S_N - 1;
    const int n = x / _S_log_2;
    const auto r = (x - n * _S_log_2) / pow2<_Tp>(_J);
    auto term = _Tp{1};
    auto expr = term;
    for (int k = 1; k < _S_N; ++k)
      expr += (term *= r / _Tp(k));
    for (int i = 0; i < _J; ++i)
      expr *= expr;
    return expr * pow2<_Tp>(n);
  }

template<typename _Tp>
  void
  test_exp(_Tp proto = _Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    auto w = 8 + std::cout.precision();

    for (int i = -100; i <= +100; ++i)
      {
	_Tp x = _Tp{0.1L} * i;
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
