/**
 *
 */

#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include <emsr/math_constants.h>
#include <emsr/sf_bessel.h>


bool VERBOSE = false;

  template <typename Tp>
    void
    cyl_bessel_jn_asymp_old(Tp nu, Tp x,
			      Tp & Jnu, Tp & Nnu)
    {
      const auto __2nu = 2 * nu;
      const auto x8 = 8 * x;
      auto k = 1;
      auto k2m1 = 1;
      auto P = Tp(1);
      auto Q = (__2nu - k2m1) * (__2nu + k2m1) / x8;
      ++k;
      const auto eps = std::numeric_limits<Tp>::epsilon();
      auto t = Tp(1);
      do
        {
          k2m1 += 2;
          t *= -(__2nu - k2m1) * (__2nu + k2m1) / (k * x8);
          auto convP = std::abs(t) < eps * std::abs(P);
          P += t;
          ++k;

          k2m1 += 2;
          t *= (__2nu - k2m1) * (__2nu + k2m1) / (k * x8);
          auto convQ = std::abs(t) < eps * std::abs(Q);
          Q += t;
          ++k;

          if (convP && convQ && k > (nu / Tp(2)))
            break;
        }
      while (k < 50 * nu);

      auto chi = x - (nu + Tp(0.5L))
        	       * emsr::pi_v<Tp> / Tp{2};
      auto c = std::cos(chi);
      auto s = std::sin(chi);

      auto coef = std::sqrt(Tp(2)
        	  / (emsr::pi_v<Tp> * x));
      Jnu = coef * (c * P - s * Q);
      Nnu = coef * (s * P + c * Q);

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
   * @param  nu  The order of the Bessel functions.
   * @param  x   The argument of the Bessel functions.
   * @return A struct containing the cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tnu, typename Tp>
    emsr::cyl_bessel_t<_Tnu, Tp, Tp>
    cyl_bessel_jn_asymp(_Tnu nu, Tp x)
    {
      using bess_t = emsr::cyl_bessel_t<Tp, Tp, Tp>;

      const auto s_eps = emsr::epsilon(x);
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = s_pi / Tp{2};
      const auto __2nu = Tp{2} * nu;
      const auto __4nu2 = __2nu * __2nu;
      const auto __8x = Tp{8} * x;
      auto k = 0;
      auto bk_xk = Tp{1};
      auto _Rsum = bk_xk;
      auto ak_xk = Tp{1};
if (VERBOSE) std::cout << ' ' << std::setw(20) << ak_xk << '\n';
      auto _Psum = ak_xk;
      auto convP = false;
      ++k;
      auto __2km1 = 1;
      bk_xk *= (__4nu2 + __2km1 * (__2km1 + 2)) / __8x;
      auto _Ssum = bk_xk;
      ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / __8x;
if (VERBOSE) std::cout << ' ' << std::setw(20) << ak_xk << '\n';
      auto _Qsum = ak_xk;
      auto convQ = false;
      auto ak_xk_prev = std::abs(ak_xk);
      do
	{
	  ++k;
	  __2km1 += 2;
	  bk_xk = -(__4nu2 + __2km1 * (__2km1 + 2)) * ak_xk / (k * __8x);
	  _Rsum += bk_xk;
	  ak_xk *= -(__2nu - __2km1) * (__2nu + __2km1) / (k * __8x);
if (VERBOSE) std::cout << ' ' << std::setw(20) << ak_xk << '\n';
	  if (k > nu / Tp{2} && std::abs(ak_xk) > ak_xk_prev)
	    break;
	  _Psum += ak_xk;
	  ak_xk_prev = std::abs(ak_xk);
	  convP = std::abs(ak_xk) < s_eps * std::abs(_Psum);

	  ++k;
	  __2km1 += 2;
	  bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * ak_xk / (k * __8x);
	  _Ssum += bk_xk;
	  ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / (k * __8x);
if (VERBOSE) std::cout << ' ' << std::setw(20) << ak_xk << '\n';
	  if (k > nu / Tp{2} && std::abs(ak_xk) > ak_xk_prev)
	    break;
	  _Qsum += ak_xk;
	  ak_xk_prev = std::abs(ak_xk);
	  convQ = std::abs(ak_xk) < s_eps * std::abs(_Qsum);

	  if (convP && convQ)
	    break;
	}
      while (k < Tp{20} * nu);

      const auto omega = x - (nu + Tp{0.5L}) * s_pi_2;
      const auto c = std::cos(omega);
      const auto s = std::sin(omega);

      const auto coef = std::sqrt(Tp{2} / (s_pi * x));
      return bess_t{nu, x,
		coef * (c * _Psum - s * _Qsum),
		-coef * (s * _Rsum + c * _Ssum),
		coef * (s * _Psum + c * _Qsum),
		coef * (c * _Rsum - s * _Ssum)};
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
   * @param  nu  The order of the Bessel functions.
   * @param  x   The argument of the Bessel functions.
   * @return A struct containing the modified cylindrical Bessel functions
   *         of the first and second kinds and their derivatives.
   */
  template<typename _Tnu, typename Tp>
    emsr::cyl_mod_bessel_t<_Tnu, Tp, Tp>
    cyl_bessel_ik_asymp(_Tnu nu, Tp x)
    {
      using bess_t = emsr::cyl_mod_bessel_t<Tp, Tp, Tp>;
      const auto s_eps = emsr::epsilon(x);
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_pi_2 = s_pi / Tp{2};
      const auto __2nu = Tp{2} * nu;
      const auto __4nu2 = __2nu * __2nu;
      const auto __8x = Tp{8} * x;
      auto k = 0;
      auto bk_xk = Tp{1};
      auto _Rsum = bk_xk;
      auto ak_xk = Tp{1};
if (VERBOSE) std::cout << ' ' << std::setw(20) << ak_xk << '\n';
      auto _Psum = ak_xk;
      auto convP = false;
      ++k;
      auto __2km1 = 1;
      bk_xk *= (__4nu2 + __2km1 * (__2km1 + 2)) / __8x;
      auto _Ssum = bk_xk;
      ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / __8x;
if (VERBOSE) std::cout << ' ' << std::setw(20) << ak_xk << '\n';
      auto _Qsum = ak_xk;
      auto convQ = false;
      auto ak_xk_prev = std::abs(ak_xk);
      do
	{
	  ++k;
	  __2km1 += 2;
	  bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * ak_xk / (k * __8x);
	  _Rsum += bk_xk;
	  ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / (k * __8x);
if (VERBOSE) std::cout << ' ' << std::setw(20) << ak_xk << '\n';
	  if (k > nu / Tp{2} && std::abs(ak_xk) > ak_xk_prev)
	    break;
	  _Psum += ak_xk;
	  ak_xk_prev = std::abs(ak_xk);
	  convP = std::abs(ak_xk) < s_eps * std::abs(_Psum);

	  ++k;
	  __2km1 += 2;
	  bk_xk = (__4nu2 + __2km1 * (__2km1 + 2)) * ak_xk / (k * __8x);
	  _Ssum += bk_xk;
	  ak_xk *= (__2nu - __2km1) * (__2nu + __2km1) / (k * __8x);
if (VERBOSE) std::cout << ' ' << std::setw(20) << ak_xk << '\n';
	  if (k > nu / Tp{2} && std::abs(ak_xk) > ak_xk_prev)
	    break;
	  _Qsum += ak_xk;
	  ak_xk_prev = std::abs(ak_xk);
	  convQ = std::abs(ak_xk) < s_eps * std::abs(_Qsum);

	  if (convP && convQ)
	    break;
	}
      while (k < Tp{20} * nu);

      const auto coef = std::sqrt(Tp{1} / (Tp{2} * s_pi * x));
      return bess_t{nu, x,
		      coef * std::exp(x) * (_Psum - _Qsum),
		      coef * std::exp(x) * (_Rsum - _Ssum),
		      s_pi * coef * std::exp(-x) * (_Psum + _Qsum),
		      -s_pi * coef * std::exp(-x) * (_Rsum + _Ssum)};
    }

template<typename _Tnu, typename Tp>
  void
  test_bessel_asymp(_Tnu nu, Tp x, bool use_internal, bool use_internal_old)
  {
    const auto s_pi_2 = emsr::pi_v<Tp> / Tp{2};
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = 8 + std::cout.precision();

    if (use_internal)
      {
	VERBOSE = true;
	auto Bess = emsr::detail::cyl_bessel_jn(nu, x);
	auto BessAsym = cyl_bessel_jn_asymp(nu, x);
	std::cout << '\n';
	std::cout << "nu = " << std::setw(w) << nu << '\n';
	std::cout << "x  = " << std::setw(w) << x << '\n';
	std::cout << "Jnu  = " << std::setw(w) << Bess.J_value << '\n';
	std::cout << "Nnu  = " << std::setw(w) << Bess.N_value << '\n';
	std::cout << "Jnua = " << std::setw(w) << BessAsym.J_value << '\n';
	std::cout << "Nnua = " << std::setw(w) << BessAsym.N_value << '\n';
	std::cout << "Jnua - Jnu = " << std::setw(w) << BessAsym.J_value - Bess.J_value << '\n';
	std::cout << "Nnua - Nnu = " << std::setw(w) << BessAsym.N_value - Bess.N_value << '\n';
	std::cout << "Jpnu  = " << std::setw(w) << Bess.J_deriv << '\n';
	std::cout << "Npnu  = " << std::setw(w) << Bess.N_deriv << '\n';
	std::cout << "Jpnua = " << std::setw(w) << BessAsym.J_deriv << '\n';
	std::cout << "Npnua = " << std::setw(w) << BessAsym.N_deriv << '\n';
	std::cout << "Jpnua - Jpnu = " << std::setw(w) << BessAsym.J_deriv - Bess.J_deriv << '\n';
	std::cout << "Npnua - Npnu = " << std::setw(w) << BessAsym.N_deriv - Bess.N_deriv << '\n';
	std::cout << "pi x Wronski / 2 = " << std::setw(w) << s_pi_2 * x * BessAsym.Wronskian() << '\n';
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
	emsr::cyl_bessel_t<Tp, Tp, Tp> Bess;
	try
	{
	  Bess = emsr::detail::cyl_bessel_jn(nu, x);
	}
	catch (std::exception& e)
	{
	  std::cout << '\n' << "Couldn't run main Bessel function." << '\n';
	}

	emsr::cyl_bessel_t<Tp, Tp, Tp> BessAsym;
	if (use_internal)
	  BessAsym = cyl_bessel_jn_asymp(nu, x);
	else if (use_internal_old)
	  {
            Tp Jnua = 0.0, Nnua = 0.0, Jpnua = 0.0, Npnua = 0.0;
            cyl_bessel_jn_asymp_old(nu, x, Jnua, Nnua);
            BessAsym.nu_arg = nu;
            BessAsym.x_arg = x;
            BessAsym.J_value = Jnua;
            BessAsym.J_deriv = Jpnua;
            BessAsym.N_value = Nnua;
            BessAsym.N_deriv = Npnua;
	  }
	else
	  BessAsym = emsr::detail::cyl_bessel_jn_asymp(nu, x);

	std::cout << ' ' << std::setw(w) << x
		  << ' ' << std::setw(w) << Bess.J_value
		  << ' ' << std::setw(w) << Bess.N_value
		  << ' ' << std::setw(w) << BessAsym.J_value
		  << ' ' << std::setw(w) << BessAsym.N_value
		  << ' ' << std::setw(w) << BessAsym.J_value - Bess.J_value
		  << ' ' << std::setw(w) << BessAsym.N_value - Bess.N_value;
	if (!use_internal_old)
	  {
	    std::cout << ' ' << std::setw(w) << Bess.J_deriv
		      << ' ' << std::setw(w) << Bess.N_deriv
		      << ' ' << std::setw(w) << BessAsym.J_deriv
		      << ' ' << std::setw(w) << BessAsym.N_deriv
		      << ' ' << std::setw(w) << BessAsym.J_deriv - Bess.J_deriv
		      << ' ' << std::setw(w) << BessAsym.N_deriv - Bess.N_deriv
		      << ' ' << std::setw(w) << s_pi_2 * x * BessAsym.Wronskian();
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
