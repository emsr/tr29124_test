/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_lambert_w test_lambert_w.cpp -lquadmath
./test_lambert_w > test_lambert_w.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_lambert_w test_lambert_w.cpp -lquadmath
./test_lambert_w > test_lambert_w.txt
*/

#include <ext/cmath>

  /**
   * This is the third-order Halley's root finding algorithm for Lambert W.
   * The radius of convergence is 1/e but it staggers on pretty well up to above 3.
   */
  template<typename _Tp>
    _Tp
    __lambert_w_series(_Tp __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__z);
      const auto _S_max_iter = 1000u;

      auto _W = __z * (_Tp{1} - __z);
      auto __term = -__z * __z;
      for (auto __k = 3u; __k < _S_max_iter; ++__k)
	{
	  __term *= -__z * std::pow(_Tp(__k) / _Tp(__k - 1), __k - 2);
	  _W += __term;
	  if (std::abs(__term) < _S_eps * std::abs(_W))
	    break;
	}
      return _W;
    }

  /**
   * This is the second-order Newton root finding algorithm for Lambert W.
   */
  template<typename _Tp>
    _Tp
    __lambert_w_newton(_Tp __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__z);
      const auto _S_max_iter = 1000u;

      auto __wk = _Tp{1};
      for (auto __k = 0u; __k < _S_max_iter; ++__k)
	{
          const auto __expwk = std::exp(__wk);
          const auto __wexpwk = __wk * __expwk;
	  const auto __wkp1 = __wk - (__wexpwk - __z)
				   / (_Tp{1} + __wk) / __expwk;
	  const auto __del = std::abs(__wkp1 - __wk);
	  __wk = __wkp1;
	  if (__del < _S_eps)
	    break;
	}
      return __wk;
    }

  /**
   * This is the third-order Halley root finding algorithm for Lambert W.
   */
  template<typename _Tp>
    _Tp
    __lambert_w_halley(_Tp __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__z);
      const auto _S_max_iter = 1000u;

      auto __wk = _Tp{1};
      for (auto __k = 0u; __k < _S_max_iter; ++__k)
	{
          const auto __expwk = std::exp(__wk);
	  const auto __fact = __wk * __expwk - __z;
          const auto __wkp1 = __wk - __fact
		      / ((__wk + 1) * __expwk - (__wk + 2) * __fact / (2 * __wk + 2));
	  const auto __del = std::abs(__wkp1 - __wk);
	  __wk = __wkp1;
	  if (__del < _S_eps)
	    break;
	}
      return __wk;
    }


  /**
   * This is the fifth-order Schroder's update for Lambert W.
   */
  template<typename _Tp>
    _Tp
    __lambert_w_schroder(_Tp __z, _Tp _W)
    {
      const auto __y = __z * std::exp(-_W);
      const auto __f0 = _W - __y;
      const auto __f1 = _Tp{1} + __y;
      const auto __f2 = __y;
      const auto __f11 = __f1 * __f1;
      const auto __f0y = __f0 * __y;
      const auto __f00y = __f0 * __f0y;
      return _W - 4 *__f0 * (6 * __f1 * (__f11 + __f0y) + __f00y)
		/ (__f11 * (24 * __f11 + 36 * __f0y)
		 + 6 * __f00y * (14 * __y + 8 + __f0));
    }


/**
 * This is the asymptotic log series for @f$ W_0(z) @f$
 * as @f$ z \rightarrow \infty @f$.
 * @f[
 *    W_0(z) = \xi - ln(\xi) + \frac{ln(\xi)}{\xi}
 *           - \frac{ln(\xi)}{\xi^2} + \frac{(ln(\xi))^2}{\xi^2}
 * @f]
 * where @f$ \xi = ln(z) @f$
 */
template<typename _Tp>
  _Tp
  __lambert_w_0_log_series(_Tp __z)
  {
    const auto __xi = std::log(__z);
    const auto __lnxi = std::log(__xi);
    return __xi - __lnxi * (_Tp{1} - (_Tp{1} / __xi)
		* (_Tp{1} - (_Tp{1} / __xi) * (_Tp{1} - __lnxi / _Tp{2})));
  }


/**
 * This is the asymptotic log series for @f$ W_1(z) @f$
 * as @f$ z \rightarrow 0- @f$.
 * @f[
 *    W_1(z) = -\eta - ln(\eta) - \frac{ln(\eta)}{\eta}
 *           - \frac{ln(\eta)}{\eta^2} - \frac{(ln(\eta))^2}{\eta^2}
 * @f]
 * where @f$ \eta = ln(-1/z) @f$
 */
template<typename _Tp>
  _Tp
  __lambert_w_1_log_series(_Tp __z)
  {
    const auto __eta = std::log(-_Tp{1} / __z);
    const auto __lneta = std::log(__eta);
    return -__eta - __lneta * (_Tp{1} + (_Tp{1} / __eta)
		* (_Tp{1} + (_Tp{1} / __eta) * (_Tp{1} + __lneta / _Tp{2})));
  }


template<typename _Tp>
  void
  test_lambert_w(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const int N = 100;
    const auto _S_e = __gnu_cxx::__const_e(proto);
    const auto _S_1_e = _Tp{1} / _S_e;
    const auto del = (_S_e + _S_1_e) / N;

    for (int i = 0; i <= N; ++i)
      {
	auto z = -_S_1_e  + del * i;
        auto W_newton = __lambert_w_newton(z);
        auto W_halley = __lambert_w_halley(z);
        auto W_series = __lambert_w_series(z);
	std::cout << ' ' << std::setw(w) << z
		  << ' ' << std::setw(w) << W_newton
		  << ' ' << std::setw(w) << W_halley
		  << ' ' << std::setw(w) << W_series
		  << '\n';
      }
  }

int
main()
{
  test_lambert_w(1.0);
}
