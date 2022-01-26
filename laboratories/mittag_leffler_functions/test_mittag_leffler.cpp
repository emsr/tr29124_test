/**
 *
 */

#include <iostream>
#include <iomanip>
#include <algorithm> // max({...})
#include <cmath>

#include <emsr/numeric_limits.h>
#include <emsr/float128_io.h>
#include <emsr/special_functions.h>

#include <emsr/integration.h>

  /* Monotone integrand for the Mittag-Leffler function. */
  template<typename _Tp>
    std::complex<_Tp>
    mittag_leffler_K(_Tp alpha, _Tp beta, _Tp chi,
		       const std::complex<_Tp>& z)
    {
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto chip1 = std::pow(chi, _Tp{1} / alpha);
      const auto chip2 = std::pow(chi, (_Tp{1} - beta) / alpha);
      return chip2
	   * std::exp(-chip1)
	   * (chi * std::sin(_S_pi * (_Tp{1} - beta))
	      - z * std::sin(_S_pi * (_Tp{1} - beta + alpha)))
	   / (chi * chi - _Tp{2} * chi * z * std::cos(alpha * _S_pi)
		 + z * z)
	   / _S_pi / alpha;
    }

  /* Monotone integral for the Mittag-Leffler function. */
  template<typename _Tp>
    std::complex<_Tp>
    mittag_leffler_K_integral(_Tp alpha, _Tp beta,
				_Tp chi_min, _Tp chi_max,
				const std::complex<_Tp>& z)
    {
      const auto _S_eps = emsr::epsilon(chi_min);
      auto func = [alpha, beta, z](_Tp chi)
		    -> std::complex<_Tp>
		    { return mittag_leffler_K(alpha, beta, chi, z); };

      const auto epsabs = _Tp{100} * _S_eps;
      const auto epsrel = _Tp{0};
      auto ws = emsr::cquad_workspace<_Tp, std::complex<_Tp>>();

      auto quad
	= emsr::cquad_integrate(ws, func, chi_min, chi_max,
				     epsabs, epsrel);

      return quad.result;
    }

  /* Oscillatory integrand for the Mittag-Leffler function. */
  template<typename _Tp>
    std::complex<_Tp>
    mittag_leffler_P(_Tp alpha, _Tp beta, _Tp epsilon, _Tp phi,
		       const std::complex<_Tp>& z)
    {
      const auto _S_i = std::complex<_Tp>{0, 1};
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto epsp1 = std::pow(epsilon, _Tp{1} / alpha);
      const auto rat = _Tp{1} + (_Tp{1} - beta) / alpha;
      const auto epsp2 = std::pow(epsilon, rat);
      const auto omega = phi * rat + epsp1 * std::sin(phi / alpha);
      return epsp2
	   * std::exp(epsp1 * std::cos(phi / alpha))
	   * std::polar(_Tp{1}, omega)
	   / (epsilon * _S_i - z)
	   / _Tp{2} / _S_pi / alpha;
    }

  /* Oscillatory integral for the Mittag-Leffler function. */
  template<typename _Tp>
    std::complex<_Tp>
    mittag_leffler_P_integral(_Tp alpha, _Tp beta, _Tp epsilon,
				_Tp phi_min, _Tp phi_max,
				const std::complex<_Tp>& z)
    {
      const auto _S_eps = emsr::epsilon(phi_min);
      auto func = [alpha, beta, epsilon, z](_Tp phi)
		    -> std::complex<_Tp>
		    {
		      return mittag_leffler_P(alpha, beta,
						epsilon, phi, z);
		    };

      const auto epsabs = _Tp{100} * _S_eps;
      const auto epsrel = _Tp{0};
      auto ws = emsr::cquad_workspace<_Tp, std::complex<_Tp>>();

      auto quad
	= emsr::cquad_integrate(ws, func, phi_min, phi_max,
				     epsabs, epsrel);

      return quad.result;
    }


  /**
   * Compute the Mittag-Leffer function:
   * @f[
   *   E_{\alpha,\beta}(z) = \sum_{k=0}^{\infty}\frac{z^k}{\beta + \alpha k},
   *   \mbox{  } \alpha > 0, \beta \elem \complex, z \elem \complex
   * @f]
   *
   * @see COMPUTATION OF THE MITTAG-LEFFLER FUNCTION @f$ E_{\alpha,\beta}(z) @f$
   * AND ITS DERIVATIVE, Rudolf Gorenflo, Joulia Loutchko & Yuri Luchko
   */
  template<typename _Tp>
    std::complex<_Tp>
    mittag_leffler(_Tp alpha, _Tp beta, const std::complex<_Tp>& z)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_eps = emsr::epsilon(alpha);
      const auto _S_2pi = emsr::tau_v<_Tp>;
      const auto _S_pi = emsr::pi_v<_Tp>;

      const auto az = std::abs(z);
      if (alpha > _Tp{1})
	{
          unsigned int k0 = _Tp{1} + std::floor(alpha);
	  const auto alpha0 = alpha / k0;
	  const auto rho0 = std::pow(z, _Tp{1} / _Tp(k0));
	  const auto lamb = _S_2pi / _Tp(k0);

	  auto E = _Cmplx{0};
	  for (auto k = 0u; k < k0; ++k)
	    {
	      auto zk = rho0 * std::polar(_Tp{1}, lamb * _Tp(k));
	      E += mittag_leffler(alpha0, beta, zk);
	    }
	  return E / _Tp(k0);
	}
      else if (az < _S_eps)
	return emsr::detail::gamma_reciprocal(beta);
      else if (az < _Tp{1})
	{
	  unsigned int k0 = std::max(std::ceil((_Tp{1} - beta) / alpha),
				std::ceil(std::log(_S_eps * (_Tp{1} - az))
					    / std::log(az)));
	  auto E = _Cmplx{0};
	  auto zk = _Cmplx{1};
	  for (auto k = 0u; k <= k0; ++k)
	    {
	      const auto arg = beta + alpha * k;
	      const auto term = zk
			* emsr::detail::gamma_reciprocal(arg);
	      E += term;
	      if (std::abs(term) < _S_eps)
		break;
	      zk *= z;
	    }
	  return E;
	}
      else if (az > std::floor(_Tp{10} + _Tp{5} * alpha))
	{
	  unsigned int k0 = std::floor(-std::log(_S_eps) / std::log(az));
	  auto E = _Cmplx{0};
	  auto zk = _Cmplx{1};
	  for (auto k = 1u; k <= k0; ++k)
	    {
	      zk /= z;
	      E += zk * emsr::detail::gamma_reciprocal(beta - alpha * k);
	    }
	  if (std::arg(z)
	      < _S_pi * (alpha / _Tp{4} + std::min(_Tp{1}, alpha) / _Tp{2}))
	    {
	      const auto zp1 = std::pow(z, _Tp{1} / alpha);
	      const auto zp2 = std::pow(z, (_Tp{1} - beta) / alpha);
	      const auto extra = zp2 * std::exp(zp1) / alpha;
	      return extra - E;
	    }
	  else
	    return -E;
	}
      else
	{
	  auto chi0 = _Tp{0};
	  if (beta >= _Tp{0})
	    chi0 = std::max({_Tp{1}, _Tp{2} * az,
			std::pow(-std::log(_S_pi * _S_eps / _Tp{6}), alpha)});
	  else
	    {
	      const auto abeta = std::abs(beta);
	      chi0 = std::max({std::pow(_Tp{1} + abeta, alpha),
			 _Tp{2} * az, 
		std::pow(-_Tp{2}
			  * std::log(_S_pi * _S_eps
			  / (_Tp{6} * (abeta + _Tp{2})
			   * std::pow(_Tp{2} * abeta, abeta))),
			 alpha)});
	    }

	  const auto absarz = std::abs(std::arg(z));
	  if (absarz > alpha * _S_pi + _S_eps)
	    {
	      if (beta <= _Tp{1})
		return mittag_leffler_K_integral(alpha, beta,
						   _Tp{0}, chi0, z);
	      else
		{
		  const auto api = _S_pi * alpha;
		  return mittag_leffler_K_integral(alpha, beta,
						     _Tp{1}, chi0, z)
		       + mittag_leffler_P_integral(alpha, beta, _Tp{1},
						     -api, api, z);
		}
	    }
	  else if (absarz < alpha * _S_pi - _S_eps)
	    {
	      const auto zp1 = std::pow(z, _Tp{1} / alpha);
	      const auto zp2 = std::pow(z, (_Tp{1} - beta) / alpha);
	      const auto extra = zp2 * std::exp(zp1) / alpha;
	      if (beta <= _Tp{1})
		return mittag_leffler_K_integral(alpha, beta,
						   _Tp{0}, chi0, z)
			+ extra;
	      else
		{
		  const auto lo = az / _Tp{2};
		  const auto api = _S_pi * alpha;
		  return mittag_leffler_K_integral(alpha, beta,
						     lo, chi0, z)
		       + mittag_leffler_P_integral(alpha, beta, lo,
							-api, api, z)
		       + extra;
		}
	    }
	  else
	    {
	      const auto lo = (az + _Tp{1}) / _Tp{2};
	      const auto api = _S_pi * alpha;
	      return mittag_leffler_K_integral(alpha, beta,
						 lo, chi0, z)
		   + mittag_leffler_P_integral(alpha, beta, lo,
						 -api, api, z);
	    }
	}
    }

  /**
   * Compute the derivative of the Mittag-Leffer function:
   * @f[
   *   E_{\alpha,\beta}(z) = \sum_{k=0}^{\infty}
   *                       \frac{z^k}{\Gamma(\beta + \alpha k)},
   *   \mbox{  } \alpha > 0, \beta \elem \complex, z \elem \complex
   * @f]
   *
   * @see COMPUTATION OF THE MITTAG-LEFFLER FUNCTION @f$ E_{\alpha,\beta}(z) @f$
   * AND ITS DERIVATIVE, Rudolf Gorenflo, Joulia Loutchko & Yuri Luchko
   */
  template<typename _Tp>
    std::complex<_Tp>
    mittag_leffler_deriv(_Tp alpha, _Tp beta,
			   const std::complex<_Tp>& z)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_eps = emsr::epsilon(alpha);

      const auto az = std::abs(z);
      if (az < _Tp{1})
	{
	  auto k1 = _Tp{0};
	  if (alpha > _Tp{1})
	    k1 = _Tp{1} + (_Tp{2} - alpha - beta) / (alpha - _Tp{1});
	  else
	    {
	      const auto D = _Tp{1}
			     + alpha * (alpha - _Tp{4} * beta + _Tp{6});
	      const auto omega = alpha + beta - _Tp{3} / _Tp{2};
	      const auto rat = _Tp{1} + (_Tp{3} - alpha - beta) / alpha;
	      if (D <= _Tp{0})
		k1 = rat;
	      else
		k1 = std::max(rat,
			_Tp{1}
			+ (_Tp{1} - _Tp{2} * omega * alpha + std::sqrt(D))
				  / (2 * alpha * alpha));
	    }
	  k1 = std::ceil(k1);
	  unsigned int k0 = std::max(k1,
				 std::ceil(std::log(_S_eps * (_Tp{1} - az))
					 / std::log(az)));
	  auto Ep = _Cmplx{0};
	  auto zk = _Cmplx{1};
	  for (auto k = 0u; k <= k0; ++k)
	    {
	      Ep += _Tp(k + 1) * zk
		    * emsr::detail::gamma_reciprocal(beta
						      + alpha * _Tp(k + 1));
	      zk *= z;
	    }
	  return Ep;
	}
      else
	return (mittag_leffler(alpha, beta - _Tp{1}, z)
	      - (beta - _Tp{1}) * mittag_leffler(alpha, beta, z))
	     / alpha / z;
    }

template<typename _Tp>
  void
  test_mittag_leffler(_Tp proto = _Tp{})
  {
    using namespace std::complex_literals;
    using _Cmplx = std::complex<_Tp>;

    std::cout.precision(emsr::digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    // Figure 1
    {
      const auto alpha = _Tp{1} / _Tp{4};
      const auto beta = _Tp{1};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 100; ++i)
	{
	  auto t = i * del;
	  auto ml_val = mittag_leffler(alpha, beta, _Cmplx(-t, 0));
	  auto ml_der = -mittag_leffler_deriv(alpha, beta, _Cmplx(-t, 0));
	  std::cout << std::setw(width) << t
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::real(ml_der)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 2
    {
      const auto alpha = _Tp{7} / _Tp{4};
      const auto beta = _Tp{1};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto t = i * del;
	  auto ml_val = mittag_leffler(alpha, beta, _Cmplx(-t, 0));
	  auto ml_der = -mittag_leffler_deriv(alpha, beta, _Cmplx(-t, 0));
	  std::cout << std::setw(width) << t
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::real(ml_der)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 3
    {
      const auto alpha = _Tp{9} / _Tp{4};
      const auto beta = _Tp{1};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 1000; ++i)
	{
	  auto t = i * del;
	  auto ml_val = mittag_leffler(alpha, beta, _Cmplx(-t, 0));
	  auto ml_der = -mittag_leffler_deriv(alpha, beta, _Cmplx(-t, 0));
	  std::cout << std::setw(width) << t
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::real(ml_der)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 4
    {
      const auto alpha = _Tp{3} / _Tp{4};
      const auto beta = _Tp{1};
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto phase = alpha * _S_pi / _Tp{4};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 50; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 5
    {
      const auto alpha = _Tp{3} / _Tp{4};
      const auto beta = _Tp{1};
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto phase = alpha * _S_pi / _Tp{2};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 6
    {
      const auto alpha = _Tp{3} / _Tp{4};
      const auto beta = _Tp{1};
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto phase = _Tp{3} * alpha * _S_pi / _Tp{4};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 7
    {
      const auto alpha = _Tp{3} / _Tp{4};
      const auto beta = _Tp{1};
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto phase = _S_pi;
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 200; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 8
    {
      const auto alpha = _Tp{5} / _Tp{4};
      const auto beta = _Tp{1};
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto phase = alpha * _S_pi / _Tp{4};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 100; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 9
    {
      const auto alpha = _Tp{5} / _Tp{4};
      const auto beta = _Tp{1};
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto phase = alpha * _S_pi / _Tp{2};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 10
    {
      const auto alpha = _Tp{5} / _Tp{4};
      const auto beta = _Tp{1};
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto phase = _Tp{3} * alpha * _S_pi / _Tp{4};
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }

    // Figure 11
    {
      const auto alpha = _Tp{5} / _Tp{4};
      const auto beta = _Tp{1};
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto phase =  _S_pi;
      std::cout << '\n';
      std::cout << '\n';
      const auto del = _Tp{1} / _Tp{10};
      for (int i = 0; i <= 1000; ++i)
	{
	  auto z = std::polar(i * del, phase);
	  auto ml_val = mittag_leffler(alpha, beta, z);
	  std::cout << std::setw(width) << std::abs(z)
		    << std::setw(width) << std::abs(ml_val)
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::imag(ml_val)
		    << '\n';
	}
      std::cout << std::flush;
    }
  }

int
main()
{
  test_mittag_leffler(1.0);
}
