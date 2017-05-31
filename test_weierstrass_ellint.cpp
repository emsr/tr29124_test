/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_weierstrass_ellint test_weierstrass_ellint.cpp -lquadmath -Lwrappers/debug -lwrap_boost
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_weierstrass_ellint > test_weierstrass_ellint.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_weierstrass_ellint test_weierstrass_ellint.cpp -lquadmath -Lwrappers/debug -lwrap_boost
./test_weierstrass_ellint > test_weierstrass_ellint.txt
*/

#include <ext/cmath>
#include <complex>
#include <iostream>
#include <iomanip>

  /**
   * Return the Weierstrass elliptic function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    weierstrass_p(std::complex<_Tp> __omega1, std::complex<_Tp> __omega3,
		  std::complex<_Tp> __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_pi = __gnu_cxx::__const_pi<_Real>();
      const auto __tau = __omega3 / __omega1;
      // Confine z to the fundamental parallelogram defined by pi and tau*pi.
      const auto __zi = std::fmod(std::imag(__z), _S_pi * std::imag(__tau));
      const auto __zr = std::fmod(std::real(__z), _S_pi);
      const auto __q = std::__detail::__polar_pi(_Tp{1}, __tau);
      const auto __theta2 = std::__detail::__theta_2(__q, _Cmplx{0});
      const auto __theta2p2 = __theta2 * __theta2;
      const auto __theta2p4 = __theta2p2 * __theta2p2;
      const auto __theta4 = std::__detail::__theta_4(__q, _Cmplx{0});
      const auto __theta4p2 = __theta4 * __theta4;
      const auto __theta4p4 = __theta4p2 * __theta4p2;
      const auto __e1 = _S_pi * _S_pi * (__theta2p4 + _Tp{2} * __theta4p4)
		      / (_Tp{12} * __omega1 * __omega1);
      const _Cmplx __arg = _S_pi * __z / _Tp{2} / __omega1;
      const auto __theta3 = std::__detail::__theta_3(__q, _Cmplx{0});
      const auto __theta2a = std::__detail::__theta_2(__q, __arg);
      const auto __numer = _S_pi * __theta3 * __theta4 * __theta2a;
      const auto __theta1a = std::__detail::__theta_1(__q, __arg);
      const auto __denom = _Tp{2} * __omega1 * __theta1a;
      const auto __rat = __numer / __denom;
      return __e1 * __rat * __rat;
    }

/**/
template<typename _Tp>
  void
  test_weierstrass_ellint()
  {
    std::cout.precision(__gnu_cxx::__digits10<_Tp>());
    auto w = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto omega1 = _Tp{2};
    const auto omega3 = _Tp{3};
    const auto del = _Tp{0.0625};
    for (int ir = -100; ir <= +100; ++ir)
      for (int ii = -100; ii <= +100; ++ii)
	{
	  //const auto z = 0.1 * ir + 0.1i * ii;
	  // The above fails because of C99 complex.
	  // Which needs to die.
	  const auto z = std::complex<_Tp>(del * ir, del * ii);
	  const auto fancyP = weierstrass_p(omega1, omega3, z);
	  std::cout << ' ' << std::setw(w) << std::real(z)
		    << ' ' << std::setw(w) << std::imag(z)
		    << ' ' << std::setw(w) << std::real(fancyP)
		    << ' ' << std::setw(w) << std::imag(fancyP)
		    << '\n';
	}
  }

int
main()
{
  test_weierstrass_ellint<double>();
}
