/**
 *
 */

#include <cmath>
#include <complex>
#include <iostream>
#include <iomanip>
#include <emsr/solver_low_degree.h>

  /**
   * Return the elliptic modular function by product expansion:
   * @f[
   *    \lambda(q) = 16 q \prod_{n=1}^{\infty}
   *             \left(\frac{(1 + q^{2n})}{(1 + q^{2n-1})}\right)^8
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __elliptic_modular_prod(const std::complex<_Tp>& __q)
    {
      using _Real = emsr::num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      constexpr std::size_t _S_max_iter = 50;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__q));

      auto __qqn =  _Cmplx{1};
      auto __prod = _Cmplx{1};
      for (std::size_t __i = 1; __i < _S_max_iter; ++__i)
	{
	  __qqn *= __q;
	  auto __fact = _Real{1} / (_Real{1} + __qqn);
	  __qqn *= __q;
	  __fact *= (_Real{1} + __qqn);
	  __fact *= __fact;
	  __fact *= __fact;
	  __fact *= __fact;
	  __prod *= __fact;
	  if (std::abs(__qqn) < _S_eps)
	    break;
	}
      return _Real{16} * __q * __prod;
    }

  /**
   * Return the Dedekind eta function by product expansion:
   * @f[
   *    \eta(q) = q^{1/12} \prod_{n=1}^{\infty}(1 - q^{2n})
   * @f]
   *
   * @param __q The elliptic nome, @f$ |q| < 1 @f$.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __dedekind_eta_prod(const std::complex<_Tp>& __q)
    {
      using _Real = emsr::num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      constexpr std::size_t _S_max_iter = 50;
      const auto _S_eps = __gnu_cxx::__epsilon(std::abs(__q));

      const auto __qq = __q * __q;
      auto __qqn = _Cmplx{1};
      auto __prod = _Cmplx{1};
      for (std::size_t __i = 1; __i < _S_max_iter; ++__i)
	{
	  __qqn *= __qq;
	  __prod *= (_Real{1} - __qqn);
	  if (std::abs(__qqn) < _S_eps)
	    break;
	}
      return std::pow(__q, _Real{1} / _Real{12}) * __prod;
    }

  /**
   * Return the Weierstrass elliptic function.
   *
   * @param __omega1 The first lattice frequency.
   * @param __omega3 The third lattice frequency.
   * @param __z The argument.
   */
  template<typename _Tp1, typename _Tp3 = std::complex<_Tp1>,
	  typename _Tp = typename std::__detail::__jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome>
    // @todo Include the arg type _Tp.
    typename std::__detail::__jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome
    __weierstrass_p(_Tp1 __omega1, _Tp3 __omega3, _Tp __z)
    {
      // @todo Include the arg type _Tp.
      //using _Type = typename std::__detail::__jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome;
      using _Real = emsr::num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto _S_i = _Cmplx{0, 1};

      const auto __lattice = std::__detail::__jacobi_lattice_t(__omega1, __omega3);
      const auto __theta0 = std::__detail::__jacobi_theta_0_t(__lattice);
      const auto __roots = std::__detail::__weierstrass_roots_t(__theta0, __omega1);

      const auto __tau = __omega3 / __omega1;
      if (std::imag(__tau) <= _Real{0})
	std::__throw_domain_error("Im(omega3/omega1) must be positive.");
      const auto __q = std::exp(_S_i * _S_pi * __tau);

      const _Cmplx __arg = _S_pi * __z / _Real{2} / __omega1;
      const auto __theta2 = std::__detail::__jacobi_theta_2(__q, __arg);
      const auto __numer = _S_pi * __theta0.th3 * __theta0.th4 * __theta2;
      const auto __theta1 = std::__detail::__jacobi_theta_1(__q, __arg);
      const auto __denom = _Tp{2} * __omega1 * __theta1;
      const auto __rat = __numer / __denom;

      return __roots.__e1 + __rat * __rat;
    }

  /**
   * 
   */
  template<typename _Tp1, typename _Tp3 = std::complex<_Tp1>,
	  typename _Tp = typename std::__detail::__jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome>
    // @todo Include the arg type _Tp.
    typename std::__detail::__jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome
    __weierstrass_zeta(_Tp1 __omega1, _Tp3 __omega3, _Tp __z)
    {
      const auto __fancy_P = __weierstrass_p(__omega1, __omega3, __z);
      const auto __lattice = std::__detail::__jacobi_lattice_t(__omega1, __omega3);
      const auto __root = std::__detail::__weierstrass_roots_t(__lattice);
      const auto __carlson_RG = std::__detail::__ellint_rg(__fancy_P - __root.__e1,
							 __fancy_P - __root.__e2,
							 __fancy_P - __root.__e3);
      return _Tp{2} * __carlson_RG - __z * __fancy_P;
    }

  /**
   * 
   */
  template<typename _Tp1, typename _Tp3 = std::complex<_Tp1>,
	  typename _Tp = typename std::__detail::__jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome>
    // @todo Include the arg type _Tp.
    typename std::__detail::__jacobi_lattice_t<_Tp1, _Tp3>::_Tp_Nome
    __weierstrass_sigma(_Tp1 __omega1, _Tp3 __omega3, _Tp __z)
    {
      // @todo Include the arg type _Tp.
      using _Real = emsr::num_traits_t<_Tp>;
      const auto _S_pi = emsr::pi_v<_Tp>;
      const auto __lattice = std::__detail::__jacobi_lattice_t(__omega1, __omega3);
      const auto __theta0 = std::__detail::__jacobi_theta_0_t(__lattice);
      const auto __omega_1 = __lattice.__omega_1();
      const auto __eta_1 = __theta0.eta_1;
      const auto __fac = _S_pi / (_Real{2} * __omega_1);
      return std::exp(__eta_1 * __z * __z / (_Real{2} * __omega_1))
	   * std::__detail::__jacobi_theta_1(__lattice.__ellnome(), __fac * __z)
	   / __fac / __theta0.th1p;
    }


/**/
template<typename _Tp>
  void
  test_weierstrass_ellint()
  {
    using _Real = emsr::num_traits_t<_Tp>;
    using _Cmplx = std::complex<_Real>;

    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    std::cout << '\n';
    const auto omega1 = _Cmplx{2, 0};
    const auto omega3 = _Cmplx{0, 3};
    const auto del = _Tp{0.0625};
    for (int ir = -100; ir <= +100; ++ir)
      {
	std::cout << '\n';
	for (int ii = -100; ii <= +100; ++ii)
	  {
	    //const auto z = 0.1 * ir + 0.1i * ii;
	    // The above fails because of C99 complex.
	    // Which needs to die.
	    const auto z = _Cmplx(del * ir, del * ii);
	    const auto fancyP = __weierstrass_p(omega1, omega3, z);
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
template<typename _Tp>
  void
  test_weierstrass_zeta()
  {
    using _Real = emsr::num_traits_t<_Tp>;
    using _Cmplx = std::complex<_Real>;

    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    std::cout << '\n';
    const auto omega1 = _Cmplx{2, 0};
    const auto omega3 = _Cmplx{0, 3};
    const auto del = _Tp{0.0625};
    for (int ir = -100; ir <= +100; ++ir)
      {
	std::cout << '\n';
	for (int ii = -100; ii <= +100; ++ii)
	  {
	    //const auto z = 0.1 * ir + 0.1i * ii;
	    // The above fails because of C99 complex.
	    // Which needs to die.
	    const auto z = _Cmplx(del * ir, del * ii);
	    const auto zeta = __weierstrass_zeta(omega1, omega3, z);
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
template<typename _Tp>
  void
  test_weierstrass_sigma()
  {
    using _Real = emsr::num_traits_t<_Tp>;
    using _Cmplx = std::complex<_Real>;

    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    std::cout << '\n';
    const auto omega1 = _Cmplx{2, 0};
    const auto omega3 = _Cmplx{0, 3};
    const auto del = _Tp{0.0625};
    for (int ir = -100; ir <= +100; ++ir)
      {
	std::cout << '\n';
	for (int ii = -100; ii <= +100; ++ii)
	  {
	    //const auto z = 0.1 * ir + 0.1i * ii;
	    // The above fails because of C99 complex.
	    // Which needs to die.
	    const auto z = _Cmplx(del * ir, del * ii);
	    const auto sigma = __weierstrass_sigma(omega1, omega3, z);
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
template<typename _Tp>
  void
  test_weierstrass_invariants()
  {
    using _Real = emsr::num_traits_t<_Tp>;
    using _Cmplx = std::complex<_Real>;
    const auto _S_pi = emsr::pi_v<_Tp>;

    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    bool polar = false;
    std::cout << '\n';
    for (int ir = 1; ir < 100; ++ir)
      {
	const auto r = _Real{0.01L} * ir;
	std::cout << '\n';
	for (int iphi = 0; iphi <= 360; ++iphi)
	  {
	    const auto phi = _S_pi * iphi / _Real{180};
	    const auto q = std::polar(r, phi);
	    const auto lambda = std::__detail::__jacobi_lattice_t<_Cmplx, _Cmplx>(q);
	    const auto [g2, g3] = std::__detail::__weierstrass_invariants_t(lambda);
	    // The solvers only deal with real coefficients at this time.
	    //const auto [e1, e2, e3] = __gnu_cxx::__cubic(-g3, -g2, _Cmplx{0}, _Cmplx{4});
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
