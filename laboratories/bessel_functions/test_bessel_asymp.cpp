 /*
$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_bessel_asymp test_bessel_asymp.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin/lib64:$LD_LIBRARY_PATH ./test_bessel_asymp

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_bessel_asymp test_bessel_asymp.cpp -lquadmath
PATH=$HOME/bin/lib64:$PATH ./test_bessel_asymp
*/

#include <iostream>
#include <sstream>
#include <ext/cmath>

bool VERBOSE = false;

  template <typename _Tp>
    void
    __cyl_bessel_jn_asymp_old(_Tp __nu, _Tp __x,
			      _Tp & __Jnu, _Tp & __Nnu)
    {
      const auto __2nu = 2 * __nu;
      const auto __x8 = 8 * __x;
      auto __k = 1;
      auto __k2m1 = 1;
      auto __P = _Tp(1);
      auto __Q = (__2nu - __k2m1) * (__2nu + __k2m1) / __x8;
      ++__k;
      const auto __eps = std::numeric_limits<_Tp>::epsilon();
      auto __t = _Tp(1);
      do
        {
          __k2m1 += 2;
          __t *= -(__2nu - __k2m1) * (__2nu + __k2m1) / (__k * __x8);
          auto __convP = std::abs(__t) < __eps * std::abs(__P);
          __P += __t;
          ++__k;

          __k2m1 += 2;
          __t *= (__2nu - __k2m1) * (__2nu + __k2m1) / (__k * __x8);
          auto __convQ = std::abs(__t) < __eps * std::abs(__Q);
          __Q += __t;
          ++__k;

          if (__convP && __convQ && __k > (__nu / _Tp(2)))
            break;
        }
      while (__k < 50 * __nu);

      auto __chi = __x - (__nu + _Tp(0.5L))
        	       * __gnu_cxx::__math_constants<_Tp>::__pi_half;
      auto __c = std::cos(__chi);
      auto __s = std::sin(__chi);

      auto __coef = std::sqrt(_Tp(2)
        	  / (__gnu_cxx::__math_constants<_Tp>::__pi * __x));
      __Jnu = __coef * (__c * __P - __s * __Q);
      __Nnu = __coef * (__s * __P + __c * __Q);

      return;
    }

  /**
   * @brief This routine computes the asymptotic cylindrical Bessel
   * 	    and Neumann functions of order nu: @f$ J_{\nu}(x) @f$,
   * 	    @f$ N_{\nu}(x) @f$.  Use this for @f$ x >> nu^2 + 1 @f$.
   *
   * @f[
   *   J_{\nu}(z) = \left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k}(\nu)}{z^{2k}}
   *  - \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * and
   * @f[
   *   N_{\nu}(z) = \left(\frac{2}{\pi z}\right)^{1/2} \left(
   *    \sin(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k}(\nu)}{z^{2k}}
   *  + \cos(\omega)\sum_{k=0}^{\infty}(-1)^k\frac{a_{2k+1}(\nu)}{z^{2k+1}}
   *    \right)
   * @f]
   * where @f$ \omega = z - \nu\pi/2 - \pi/4 @f$ and
   * @f[
   *   a_{k}(\nu) = \frac{(4\nu^2 - 1^2)(4\nu^2 - 3^2)...(4\nu^2 - (2k-1)^2)}
   *                     {8^k k!}
   * @f]
   * There sums work everywhere but on the negative real axis:
   * @f$ |ph(z)| < \pi - \delta @f$.
   *
   * References:
   *  (1) Handbook of Mathematical Functions,
   * 	  ed. Milton Abramowitz and Irene A. Stegun,
   * 	  Dover Publications,
   * 	  Section 9 p. 364, Equations 9.2.5-9.2.10
   *
   * @param  __nu  The order of the Bessel functions.
   * @param  __x   The argument of the Bessel functions.
   * @return A struct containing the cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tnu, typename _Tp>
    __gnu_cxx::__cyl_bessel_t<_Tnu, _Tp, _Tp>
    __cyl_bessel_jn_asymp(_Tnu __nu, _Tp __x)
    {
      using __bess_t = __gnu_cxx::__cyl_bessel_t<_Tp, _Tp, _Tp>;

      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__x);
      const auto __2nu = _Tp{2} * __nu;
      const auto __4nu2 = __2nu * __2nu;
      const auto __8x = _Tp{8} * __x;
      auto __k = 0;
      auto __bk_xk = _Tp{1};
      auto _Rsum = __bk_xk;
      auto __ak_xk = _Tp{1};
if (VERBOSE) std::cout << ' ' << std::setw(20) << __ak_xk << '\n';
      auto _Psum = __ak_xk;
      auto __convP = false;
      ++__k;
      auto __2km1 = 1;
      __bk_xk *= (__4nu2 + __2km1 * (__2km1 + 2)) / __8x;
      auto _Ssum = __bk_xk;
      __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / __8x;
if (VERBOSE) std::cout << ' ' << std::setw(20) << __ak_xk << '\n';
      auto _Qsum = __ak_xk;
      auto __convQ = false;
      auto __ak_xk_prev = std::abs(__ak_xk);
      do
	{
	  ++__k;
	  __2km1 += 2;
	  __bk_xk = -(__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk / (__k * __8x);
	  _Rsum += __bk_xk;
	  __ak_xk *= -(__2nu - __2km1) * (__2nu + __2km1) / (__k * __8x);
if (VERBOSE) std::cout << ' ' << std::setw(20) << __ak_xk << '\n';
	  if (__k > __nu / _Tp{2} && std::abs(__ak_xk) > __ak_xk_prev)
	    break;
	  _Psum += __ak_xk;
	  __ak_xk_prev = std::abs(__ak_xk);
	  __convP = std::abs(__ak_xk) < _S_eps * std::abs(_Psum);

	  ++__k;
	  __2km1 += 2;
	  __bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk / (__k * __8x);
	  _Ssum += __bk_xk;
	  __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / (__k * __8x);
if (VERBOSE) std::cout << ' ' << std::setw(20) << __ak_xk << '\n';
	  if (__k > __nu / _Tp{2} && std::abs(__ak_xk) > __ak_xk_prev)
	    break;
	  _Qsum += __ak_xk;
	  __ak_xk_prev = std::abs(__ak_xk);
	  __convQ = std::abs(__ak_xk) < _S_eps * std::abs(_Qsum);

	  if (__convP && __convQ)
	    break;
	}
      while (__k < _Tp{20} * __nu);

      const auto __omega = __x - (__nu + _Tp{0.5L}) * _S_pi_2;
      const auto __c = std::cos(__omega);
      const auto __s = std::sin(__omega);

      const auto __coef = std::sqrt(_Tp{2} / (_S_pi * __x));
      return __bess_t{__nu, __x,
		__coef * (__c * _Psum - __s * _Qsum),
		-__coef * (__s * _Rsum + __c * _Ssum),
		__coef * (__s * _Psum + __c * _Qsum),
		__coef * (__c * _Rsum - __s * _Ssum)};
    }

  /**
   * @brief This routine computes the asymptotic modified cylindrical
   * 	    Bessel and functions of order nu: @f$ I_{\nu}(x) @f$,
   * 	    @f$ N_{\nu}(x) @f$.  Use this for @f$ x >> nu^2 + 1 @f$.
   *
   * References:
   *  (1) Handbook of Mathematical Functions,
   * 	  ed. Milton Abramowitz and Irene A. Stegun,
   * 	  Dover Publications,
   * 	  Section 9 p. 364, Equations 9.2.5-9.2.10
   *
   * @param  __nu  The order of the Bessel functions.
   * @param  __x   The argument of the Bessel functions.
   * @return A struct containing the modified cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tnu, typename _Tp>
    __gnu_cxx::__cyl_mod_bessel_t<_Tnu, _Tp, _Tp>
    __cyl_bessel_ik_asymp(_Tnu __nu, _Tp __x)
    {
      using __bess_t = __gnu_cxx::__cyl_mod_bessel_t<_Tp, _Tp, _Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(__x);
      const auto __2nu = _Tp{2} * __nu;
      const auto __4nu2 = __2nu * __2nu;
      const auto __8x = _Tp{8} * __x;
      auto __k = 0;
      auto __bk_xk = _Tp{1};
      auto _Rsum = __bk_xk;
      auto __ak_xk = _Tp{1};
if (VERBOSE) std::cout << ' ' << std::setw(20) << __ak_xk << '\n';
      auto _Psum = __ak_xk;
      auto __convP = false;
      ++__k;
      auto __2km1 = 1;
      __bk_xk *= (__4nu2 + __2km1 * (__2km1 + 2)) / __8x;
      auto _Ssum = __bk_xk;
      __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / __8x;
if (VERBOSE) std::cout << ' ' << std::setw(20) << __ak_xk << '\n';
      auto _Qsum = __ak_xk;
      auto __convQ = false;
      auto __ak_xk_prev = std::abs(__ak_xk);
      do
	{
	  ++__k;
	  __2km1 += 2;
	  __bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk / (__k * __8x);
	  _Rsum += __bk_xk;
	  __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / (__k * __8x);
if (VERBOSE) std::cout << ' ' << std::setw(20) << __ak_xk << '\n';
	  if (__k > __nu / _Tp{2} && std::abs(__ak_xk) > __ak_xk_prev)
	    break;
	  _Psum += __ak_xk;
	  __ak_xk_prev = std::abs(__ak_xk);
	  __convP = std::abs(__ak_xk) < _S_eps * std::abs(_Psum);

	  ++__k;
	  __2km1 += 2;
	  __bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * __ak_xk / (__k * __8x);
	  _Ssum += __bk_xk;
	  __ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / (__k * __8x);
if (VERBOSE) std::cout << ' ' << std::setw(20) << __ak_xk << '\n';
	  if (__k > __nu / _Tp{2} && std::abs(__ak_xk) > __ak_xk_prev)
	    break;
	  _Qsum += __ak_xk;
	  __ak_xk_prev = std::abs(__ak_xk);
	  __convQ = std::abs(__ak_xk) < _S_eps * std::abs(_Qsum);

	  if (__convP && __convQ)
	    break;
	}
      while (__k < _Tp{20} * __nu);

      const auto __coef = std::sqrt(_Tp{1} / (_Tp{2} * _S_pi * __x));
      return __bess_t{__nu, __x,
		      __coef * std::exp(__x) * (_Psum - _Qsum),
		      __coef * std::exp(__x) * (_Rsum - _Ssum),
		      _S_pi * __coef * std::exp(-__x) * (_Psum + _Qsum),
		      -_S_pi * __coef * std::exp(-__x) * (_Rsum + _Ssum)};
    }

template<typename _Tnu, typename _Tp>
  void
  test_bessel_asymp(_Tnu nu, _Tp x, bool use_internal, bool use_internal_old)
  {
    const auto _S_pi_2 = __gnu_cxx::__const_pi_half(x);
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = 8 + std::cout.precision();

    if (use_internal)
      {
	VERBOSE = true;
	auto Bess = std::__detail::__cyl_bessel_jn(nu, x);
	auto BessAsym = __cyl_bessel_jn_asymp(nu, x);
	std::cout << '\n';
	std::cout << "nu = " << std::setw(w) << nu << '\n';
	std::cout << "x  = " << std::setw(w) << x << '\n';
	std::cout << "Jnu  = " << std::setw(w) << Bess.__J_value << '\n';
	std::cout << "Nnu  = " << std::setw(w) << Bess.__N_value << '\n';
	std::cout << "Jnua = " << std::setw(w) << BessAsym.__J_value << '\n';
	std::cout << "Nnua = " << std::setw(w) << BessAsym.__N_value << '\n';
	std::cout << "Jnua - Jnu = " << std::setw(w) << BessAsym.__J_value - Bess.__J_value << '\n';
	std::cout << "Nnua - Nnu = " << std::setw(w) << BessAsym.__N_value - Bess.__N_value << '\n';
	std::cout << "Jpnu  = " << std::setw(w) << Bess.__J_deriv << '\n';
	std::cout << "Npnu  = " << std::setw(w) << Bess.__N_deriv << '\n';
	std::cout << "Jpnua = " << std::setw(w) << BessAsym.__J_deriv << '\n';
	std::cout << "Npnua = " << std::setw(w) << BessAsym.__N_deriv << '\n';
	std::cout << "Jpnua - Jpnu = " << std::setw(w) << BessAsym.__J_deriv - Bess.__J_deriv << '\n';
	std::cout << "Npnua - Npnu = " << std::setw(w) << BessAsym.__N_deriv - Bess.__N_deriv << '\n';
	std::cout << "pi x Wronski / 2 = " << std::setw(w) << _S_pi_2 * x * BessAsym.__Wronskian() << '\n';
	VERBOSE = false;
      }

    std::cout << '\n';
    std::cout << ' ' << std::setw(w) << "x"
	      << ' ' << std::setw(w) << "Jnu"
	      << ' ' << std::setw(w) << "Nnu"
	      << ' ' << std::setw(w) << "Jnua"
	      << ' ' << std::setw(w) << "Nnua"
	      << ' ' << std::setw(w) << "Jnua - Jnu"
	      << ' ' << std::setw(w) << "Nnua - Nnu";
    if (!use_internal_old)
      {
	std::cout << ' ' << std::setw(w) << "Jpnu"
		  << ' ' << std::setw(w) << "Npnu"
		  << ' ' << std::setw(w) << "Jpnua"
		  << ' ' << std::setw(w) << "Npnua"
		  << ' ' << std::setw(w) << "Jpnua - Jpnu"
		  << ' ' << std::setw(w) << "Npnua - Npnu"
		  << ' ' << std::setw(w) << "pi x Wronski / 2";
      }
    std::cout << '\n';

    do
      {
	__gnu_cxx::__cyl_bessel_t<_Tp, _Tp, _Tp> Bess;
	try
	{
	  Bess = std::__detail::__cyl_bessel_jn(nu, x);
	}
	catch (std::exception& e)
	{
	  std::cout << '\n' << "Couldn't run main Bessel function." << '\n';
	}

	__gnu_cxx::__cyl_bessel_t<_Tp, _Tp, _Tp> BessAsym;
	if (use_internal)
	  BessAsym = __cyl_bessel_jn_asymp(nu, x);
	else if (use_internal_old)
	  {
            _Tp Jnua = 0.0, Nnua = 0.0, Jpnua = 0.0, Npnua = 0.0;
            __cyl_bessel_jn_asymp_old(nu, x, Jnua, Nnua);
            BessAsym.__nu_arg = nu;
            BessAsym.__x_arg = x;
            BessAsym.__J_value = Jnua;
            BessAsym.__J_deriv = Jpnua;
            BessAsym.__N_value = Nnua;
            BessAsym.__N_deriv = Npnua;
	  }
	else
	  BessAsym = std::__detail::__cyl_bessel_jn_asymp(nu, x);

	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << Bess.__J_value
		  << ' ' << std::setw(w) << Bess.__N_value
		  << ' ' << std::setw(w) << BessAsym.__J_value
		  << ' ' << std::setw(w) << BessAsym.__N_value
		  << ' ' << std::setw(w) << BessAsym.__J_value - Bess.__J_value
		  << ' ' << std::setw(w) << BessAsym.__N_value - Bess.__N_value;
	if (!use_internal_old)
	  {
	    std::cout << ' ' << std::setw(w) << Bess.__J_deriv
		      << ' ' << std::setw(w) << Bess.__N_deriv
		      << ' ' << std::setw(w) << BessAsym.__J_deriv
		      << ' ' << std::setw(w) << BessAsym.__N_deriv
		      << ' ' << std::setw(w) << BessAsym.__J_deriv - Bess.__J_deriv
		      << ' ' << std::setw(w) << BessAsym.__N_deriv - Bess.__N_deriv
		      << ' ' << std::setw(w) << _S_pi_2 * x * BessAsym.__Wronskian();
	  }
	std::cout << '\n';

	x += 1000.0;
      }
    while (x <= 100000.0);
  }

int
main(int n_app_args, char ** app_arg)
{
  long double nu = 20.0L;
  long double x = 1000.0L;

  bool use_internal = false;
  bool use_internal_old = false;
  if (n_app_args > 1)
    {
      int i_internal = 0;
      std::istringstream in(app_arg[1]);
      in >> i_internal;
      if (i_internal > 0)
	use_internal = true;
      if (i_internal < 0)
	use_internal_old = true;
    }
  if (n_app_args > 2)
    {
      std::istringstream in(app_arg[2]);
      in >> nu;
    }
  if (n_app_args > 3)
    {
      std::istringstream in(app_arg[3]);
      in >> x;
    }

  std::cout << '\n';
  std::cout << "use internal function     = " << std::boolalpha << use_internal << '\n';

  std::cout << '\n';
  std::cout << "use old internal function = " << std::boolalpha << use_internal_old << '\n';

  std::cout << '\n';
  std::cout << "nu = " << nu << '\n';

  std::cout << '\n';
  std::cout << "x  = " << x << '\n';

  test_bessel_asymp(float(nu), float(x), use_internal, use_internal_old);

  test_bessel_asymp(double(nu), double(x), use_internal, use_internal_old);

  test_bessel_asymp(nu, x, use_internal, use_internal_old);

  return 0;
}
