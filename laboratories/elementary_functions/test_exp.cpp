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
  __pow2(int __n)
  {
    if (__n == 0)
      return _Tp{1};
    else
      {
	const long __p = 0x1L << std::abs(__n);
	return __n > 0 ? _Tp(__p) : _Tp{1} / _Tp(__p);
      }
  }

/**
 * 
 */
template<typename _Tp>
  constexpr _Tp
  __exp_bailey(_Tp __x, int _J = 8)
  {
    constexpr auto _S_N = (std::numeric_limits<_Tp>::max_digits10 + 1) / 2;
    constexpr auto _S_log_2 = __gnu_cxx::numbers::__ln_2_v<_Tp>;
    if (_J <= 0)
      _J = _S_N - 1;
    const int __n = __x / _S_log_2;
    const auto __r = (__x - __n * _S_log_2) / __pow2<_Tp>(_J);
    auto __term = _Tp{1};
    auto __expr = __term;
    for (int __k = 1; __k < _S_N; ++__k)
      __expr += (__term *= __r / _Tp(__k));
    for (int __i = 0; __i < _J; ++__i)
      __expr *= __expr;
    return __expr * __pow2<_Tp>(__n);
  }

template<typename _Tp>
  void
  test_exp(_Tp __proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(__proto));
    auto w = 8 + std::cout.precision();

    for (int i = -100; i <= +100; ++i)
      {
	_Tp x = _Tp{0.1L} * i;
	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << __exp_bailey(x)
		  << ' ' << std::setw(w) << __exp_bailey(x, 0)
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
