/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_weierstrass_ellint test_weierstrass_ellint.cpp -lquadmath -Lwrappers/debug -lwrap_boost
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_weierstrass_ellint > test_weierstrass_ellint.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_weierstrass_ellint test_weierstrass_ellint.cpp -lquadmath -Lwrappers/debug -lwrap_boost
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_weierstrass_ellint > test_weierstrass_ellint.txt
*/

#include <ext/cmath>
#include <complex>
#include <iostream>
#include <iomanip>

  /**
   * Return the elliptic modular function.
   * @f[
   *    \lambda(q) = 16 q \prod_{n=1}^{\infty}
   *             \left(\frac{(1 + q^{2n})}{(1 + q^{2n-1})}\right)^8
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __elliptic_modular_prod(const std::complex<_Tp>& __q)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
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
   * Return the Dedekind eta function.
   * @f[
   *    \eta(q) = q^{1/12} \prod_{n=1}^{\infty}(1 - q^{2n})
   * @f]
   */
  template<typename _Tp>
    std::complex<_Tp>
    __dedekind_eta_prod(const std::complex<_Tp>& __q)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
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
   */
  template<typename _Tp>
    std::complex<_Tp>
    __weierstrass_p(std::complex<_Tp> __omega1, std::complex<_Tp> __omega3,
		    std::complex<_Tp> __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;

      const auto _S_pi = __gnu_cxx::__const_pi<_Real>();
      const auto _S_i = _Cmplx{0, 1};

      const auto __tau = __omega3 / __omega1;
      if (std::imag(__tau) <= _Tp{0})
	std::__throw_domain_error("Im(omega3/omega1) must be positive.");
      const auto __q = std::exp(_S_i * _S_pi * __tau);

      const auto __theta2 = std::__detail::__jacobi_theta_2(__q, _Cmplx{0});
      const auto __theta2p2 = __theta2 * __theta2;
      const auto __theta2p4 = __theta2p2 * __theta2p2;
      const auto __theta4 = std::__detail::__jacobi_theta_4(__q, _Cmplx{0});
      const auto __theta4p2 = __theta4 * __theta4;
      const auto __theta4p4 = __theta4p2 * __theta4p2;
      const auto __e1 = _S_pi * _S_pi * (__theta2p4 + _Tp{2} * __theta4p4)
		      / (_Tp{12} * __omega1 * __omega1);

      const _Cmplx __arg = _S_pi * __z / _Tp{2} / __omega1;
      const auto __theta3 = std::__detail::__jacobi_theta_3(__q, _Cmplx{0});
      const auto __theta2a = std::__detail::__jacobi_theta_2(__q, __arg);
      const auto __numer = _S_pi * __theta3 * __theta4 * __theta2a;
      const auto __theta1a = std::__detail::__jacobi_theta_1(__q, __arg);
      const auto __denom = _Tp{2} * __omega1 * __theta1a;
      const auto __rat = __numer / __denom;

      return __e1 + __rat * __rat;
    }

/**/
template<typename _Tp>
  void
  test_weierstrass_ellint()
  {
    using _Real = std::__detail::__num_traits_t<_Tp>;
    using _Cmplx = std::complex<_Real>;

    std::cout.precision(__gnu_cxx::__digits10<_Real>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

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

int
main()
{
  test_weierstrass_ellint<double>();
}
