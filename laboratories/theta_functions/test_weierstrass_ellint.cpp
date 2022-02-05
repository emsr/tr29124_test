/**
 *
 */

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>

#include <emsr/solver_low_degree.h>
#include <emsr/fp_type_util.h>
#include <emsr/numeric_limits.h>
#include <emsr/sf_theta.h>

  /**
   * Return the elliptic modular function by product expansion:
   * @f[
   *    \lambda(q) = 16 q \prod_{n=1}^{\infty}
   *             \left(\frac{(1 + q^{2n})}{(1 + q^{2n-1})}\right)^8
   * @f]
   *
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   */
  template<typename Tp>
    std::complex<Tp>
    elliptic_modular_prod(const std::complex<Tp>& q)
    {
      using Real = emsr::num_traits_t<Tp>;
      using Cmplx = std::complex<Real>;
      constexpr std::size_t s_max_iter = 50;
      const auto s_eps = emsr::epsilon(std::abs(q));

      auto qqn =  Cmplx{1};
      auto prod = Cmplx{1};
      for (std::size_t i = 1; i < s_max_iter; ++i)
	{
	  qqn *= q;
	  auto fact = Real{1} / (Real{1} + qqn);
	  qqn *= q;
	  fact *= (Real{1} + qqn);
	  fact *= fact;
	  fact *= fact;
	  fact *= fact;
	  prod *= fact;
	  if (std::abs(qqn) < s_eps)
	    break;
	}
      return Real{16} * q * prod;
    }

  /**
   * Return the Dedekind eta function by product expansion:
   * @f[
   *    \eta(q) = q^{1/12} \prod_{n=1}^{\infty}(1 - q^{2n})
   * @f]
   *
   * @param q The elliptic nome, @f$ |q| < 1 @f$.
   */
  template<typename Tp>
    std::complex<Tp>
    dedekind_eta_prod(const std::complex<Tp>& q)
    {
      using Real = emsr::num_traits_t<Tp>;
      using Cmplx = std::complex<Real>;
      constexpr std::size_t s_max_iter = 50;
      const auto s_eps = emsr::epsilon(std::abs(q));

      const auto qq = q * q;
      auto qqn = Cmplx{1};
      auto prod = Cmplx{1};
      for (std::size_t i = 1; i < s_max_iter; ++i)
	{
	  qqn *= qq;
	  prod *= (Real{1} - qqn);
	  if (std::abs(qqn) < s_eps)
	    break;
	}
      return std::pow(q, Real{1} / Real{12}) * prod;
    }

  /**
   * Return the Weierstrass elliptic function.
   *
   * @param omega1 The first lattice frequency.
   * @param omega3 The third lattice frequency.
   * @param z The argument.
   */
  template<typename _Tp1, typename _Tp3 = std::complex<_Tp1>,
	  typename Tp = typename emsr::detail::jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome>
    // @todo Include the arg type Tp.
    typename emsr::detail::jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome
    weierstrass_p(_Tp1 omega1, _Tp3 omega3, Tp z)
    {
      // @todo Include the arg type Tp.
      //using _Type = typename emsr::detail::jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome;
      using Real = emsr::num_traits_t<Tp>;
      using Cmplx = std::complex<Real>;
      const auto s_pi = emsr::pi_v<Tp>;
      const auto s_i = Cmplx{0, 1};

      const auto lattice = emsr::detail::jacobi_lattice_t(omega1, omega3);
      const auto theta0 = emsr::detail::jacobi_theta_0_t(lattice);
      const auto roots = emsr::detail::weierstrass_roots_t(theta0, omega1);

      const auto tau = omega3 / omega1;
      if (std::imag(tau) <= Real{0})
	throw std::domain_error("Im(omega3/omega1) must be positive.");
      const auto q = std::exp(s_i * s_pi * tau);

      const Cmplx arg = s_pi * z / Real{2} / omega1;
      const auto theta2 = emsr::detail::jacobi_theta_2(q, arg);
      const auto numer = s_pi * theta0.th3 * theta0.th4 * theta2;
      const auto theta1 = emsr::detail::jacobi_theta_1(q, arg);
      const auto denom = Tp{2} * omega1 * theta1;
      const auto rat = numer / denom;

      return roots.e1 + rat * rat;
    }

  /**
   * 
   */
  template<typename _Tp1, typename _Tp3 = std::complex<_Tp1>,
	  typename Tp = typename emsr::detail::jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome>
    // @todo Include the arg type Tp.
    typename emsr::detail::jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome
    weierstrass_zeta(_Tp1 omega1, _Tp3 omega3, Tp z)
    {
      const auto fancy_P = weierstrass_p(omega1, omega3, z);
      const auto lattice = emsr::detail::jacobi_lattice_t(omega1, omega3);
      const auto root = emsr::detail::weierstrass_roots_t(lattice);
      const auto carlson_RG = emsr::detail::ellint_rg(fancy_P - root.e1,
							 fancy_P - root.e2,
							 fancy_P - root.e3);
      return Tp{2} * carlson_RG - z * fancy_P;
    }

  /**
   * 
   */
  template<typename _Tp1, typename _Tp3 = std::complex<_Tp1>,
	  typename Tp = typename emsr::detail::jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome>
    // @todo Include the arg type Tp.
    typename emsr::detail::jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome
    weierstrass_sigma(_Tp1 omega1, _Tp3 omega3, Tp z)
    {
      // @todo Include the arg type Tp.
      using Real = emsr::num_traits_t<Tp>;
      const auto s_pi = emsr::pi_v<Tp>;
      const auto lattice = emsr::detail::jacobi_lattice_t(omega1, omega3);
      const auto theta0 = emsr::detail::jacobi_theta_0_t(lattice);
      const auto omega_1 = lattice.omega_1();
      const auto eta_1 = theta0.eta_1;
      const auto fac = s_pi / (Real{2} * omega_1);
      return std::exp(eta_1 * z * z / (Real{2} * omega_1))
	   * emsr::detail::jacobi_theta_1(lattice.ellnome(), fac * z)
	   / fac / theta0.th1p;
    }


/**/
template<typename Tp>
  void
  test_weierstrass_ellint()
  {
    using Real = emsr::num_traits_t<Tp>;
    using Cmplx = std::complex<Real>;

    std::cout.precision(emsr::digits10<Real>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    std::cout << '\n';
    const auto omega1 = Cmplx{2, 0};
    const auto omega3 = Cmplx{0, 3};
    const auto del = Tp{0.0625};
    for (int ir = -100; ir <= +100; ++ir)
      {
	std::cout << '\n';
	for (int ii = -100; ii <= +100; ++ii)
	  {
	    //const auto z = 0.1 * ir + 0.1i * ii;
	    // The above fails because of C99 complex.
	    // Which needs to die.
	    const auto z = Cmplx(del * ir, del * ii);
	    const auto fancyP = weierstrass_p(omega1, omega3, z);
	    std::cout << ' ' << std::setw(w) << std::real(z)
		      << ' ' << std::setw(w) << std::imag(z)
		      << ' ' << std::setw(w) << std::real(fancyP)
		      << ' ' << std::setw(w) << std::imag(fancyP)
		      << ' ' << std::setw(w) << std::abs(fancyP)
		      << '\n';
	  }
      }
  }


/**/
template<typename Tp>
  void
  test_weierstrass_zeta()
  {
    using Real = emsr::num_traits_t<Tp>;
    using Cmplx = std::complex<Real>;

    std::cout.precision(emsr::digits10<Real>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    std::cout << '\n';
    const auto omega1 = Cmplx{2, 0};
    const auto omega3 = Cmplx{0, 3};
    const auto del = Tp{0.0625};
    for (int ir = -100; ir <= +100; ++ir)
      {
	std::cout << '\n';
	for (int ii = -100; ii <= +100; ++ii)
	  {
	    //const auto z = 0.1 * ir + 0.1i * ii;
	    // The above fails because of C99 complex.
	    // Which needs to die.
	    const auto z = Cmplx(del * ir, del * ii);
	    const auto zeta = weierstrass_zeta(omega1, omega3, z);
	    std::cout << ' ' << std::setw(w) << std::real(z)
		      << ' ' << std::setw(w) << std::imag(z)
		      << ' ' << std::setw(w) << std::real(zeta)
		      << ' ' << std::setw(w) << std::imag(zeta)
		      << ' ' << std::setw(w) << std::abs(zeta)
		      << '\n';
	  }
      }
  }


/**/
template<typename Tp>
  void
  test_weierstrass_sigma()
  {
    using Real = emsr::num_traits_t<Tp>;
    using Cmplx = std::complex<Real>;

    std::cout.precision(emsr::digits10<Real>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    std::cout << '\n';
    const auto omega1 = Cmplx{2, 0};
    const auto omega3 = Cmplx{0, 3};
    const auto del = Tp{0.0625};
    for (int ir = -100; ir <= +100; ++ir)
      {
	std::cout << '\n';
	for (int ii = -100; ii <= +100; ++ii)
	  {
	    //const auto z = 0.1 * ir + 0.1i * ii;
	    // The above fails because of C99 complex.
	    // Which needs to die.
	    const auto z = Cmplx(del * ir, del * ii);
	    const auto sigma = weierstrass_sigma(omega1, omega3, z);
	    std::cout << ' ' << std::setw(w) << std::real(z)
		      << ' ' << std::setw(w) << std::imag(z)
		      << ' ' << std::setw(w) << std::real(sigma)
		      << ' ' << std::setw(w) << std::imag(sigma)
		      << ' ' << std::setw(w) << std::abs(sigma)
		      << '\n';
	  }
      }
  }

/**/
template<typename Tp>
  void
  test_weierstrass_invariants()
  {
    using Real = emsr::num_traits_t<Tp>;
    using Cmplx = std::complex<Real>;
    const auto s_pi = emsr::pi_v<Tp>;

    std::cout.precision(emsr::digits10<Real>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    bool polar = false;
    std::cout << '\n';
    for (int ir = 1; ir < 100; ++ir)
      {
	const auto r = Real{0.01L} * ir;
	std::cout << '\n';
	for (int iphi = 0; iphi <= 360; ++iphi)
	  {
	    const auto phi = s_pi * iphi / Real{180};
	    const auto q = std::polar(r, phi);
	    const auto lambda = emsr::detail::jacobi_lattice_t<Cmplx, Cmplx>(q);
	    const auto [g2, g3] = emsr::detail::weierstrass_invariants_t(lambda);
	    // The solvers only deal with real coefficients at this time.
	    //const auto [e1, e2, e3] = emsr::cubic(-g3, -g2, Cmplx{0}, Cmplx{4});
	    if (polar)
	      std::cout << ' ' << std::setw(w) << r
			<< ' ' << std::setw(w) << phi;
	    else
	      std::cout << ' ' << std::setw(w) << q.real()
			<< ' ' << std::setw(w) << q.imag();
	    std::cout << ' ' << std::setw(w) << std::real(g2)
		      << ' ' << std::setw(w) << std::imag(g2)
		      << ' ' << std::setw(w) << std::abs(g2)
		      << ' ' << std::setw(w) << std::real(g3)
		      << ' ' << std::setw(w) << std::imag(g3)
		      << ' ' << std::setw(w) << std::abs(g3)
		      << '\n';
	  }
      }
  }

int
main()
{
  test_weierstrass_ellint<double>();
  test_weierstrass_invariants<double>();
  test_weierstrass_zeta<double>();
  test_weierstrass_sigma<double>();
}
