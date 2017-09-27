/*
$HOME/bin/bin/g++ -std=c++17 -g -Wall -Wextra -o test_meixner_pollaczek test_meixner_pollaczek.cpp
./test_meixner_pollaczek > test_meixner_pollaczek.txt
*/

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>

template<typename _Tp>
  struct
  __meixner_pollaczek_t
  {
    std::complex<_Tp> __value;
    std::complex<_Tp> __factor;
  };

/**
 * Compute the Meixner-Pollaczek polynomial by recursion:
 * @f[
 *    0 = (n + 1) P_{n+1}(x)
 *      - 2 [x \sin \phi + (n + \lambda) \cos \phi] P_n(x)
 *      + (n + 2\lambda - 1) P_{n-1}(x) = 0
 * @f]
 * where @f$ P_n(x) = P^{(\lambda)}_n(x;\phi) @f$
 */
template<typename _Tp, typename _TpX>
  __meixner_pollaczek_t<_Tp>
  __meixner_pollaczek_recur(int n, _Tp lambda, _Tp phi, _TpX x)
  {
    std::complex<_Tp> Pnm1 = _Tp{1};
    if (n == 0)
      return {Pnm1, _Tp{1}};

    constexpr std::complex<_Tp> i{0, 1};
    const auto ix = i * _Tp(x);
    const auto sinphi = std::sin(phi);
    const auto cosphi = std::cos(phi);

    auto fact = std::complex<_Tp>{1};

    std::complex<_Tp> Pn = _Tp{2} * lambda * std::polar(_Tp{1}, phi)
			 - _Tp{2} * i * (lambda + ix) * std::sin(phi);

    if (n == 1)
      return {Pn, fact};

    auto Pnp1 = (_Tp{2} * (_Tp(x) * sinphi + (lambda + _Tp{1}) * cosphi) * Pn
		 - _Tp{2} * lambda * Pnm1) / _Tp{2};

    for (int k = 2; k < n; ++k)
      {
	Pnm1 = Pn;
	Pn = Pnp1;
	Pnp1 = (_Tp{2} * (_Tp(x) * sinphi + (lambda + _Tp(k)) * cosphi) * Pn
		 - (_Tp{2} * lambda + _Tp(k - 1)) * Pnm1) / _Tp(k + 1);
      }

    return {Pnp1, fact};
  }

/**
 * Compute the Meixner-Pollaczek polynomial defined by
 * @f[
 *    P^{(\lambda)}_n(x;\phi) = \frac{(2\lambda)_n}{n!}
 *             {}_2F_1(-n, \lambda + ix; 2\lambda; 1 - e^{-i2\phi})
 * @f]
 */
template<typename _Tp, typename _TpX>
  __meixner_pollaczek_t<_Tp>
  __meixner_pollaczek(int n, _Tp lambda, _Tp phi, _TpX x)
  {
    if (std::isnan(lambda))
      return {std::complex<_Tp>(lambda), std::complex<_Tp>{}};
    else if (std::isnan(phi))
      return {std::complex<_Tp>(phi), std::complex<_Tp>{}};
    else if (std::isnan(x))
      return {std::complex<_Tp>(x), std::complex<_Tp>{}};
    else
      return __meixner_pollaczek_recur(n, lambda, phi, x);
  }

template<typename _Tp>
  void
  test_meixner_pollaczek(int n_max, _Tp lambda, _Tp phi)
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = -320; i <= 320; ++i)
	  {
	    auto x = i * _Tp{0.01L};
	    auto P = __meixner_pollaczek(n, lambda, phi, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << std::real(P.__value)
		      << ' ' << std::setw(w) << std::imag(P.__value)
		      << '\n';
	  }
      }
  }

int
main()
{
  constexpr long double pi = 3.1415926535897932384626433832795029L;
  test_meixner_pollaczek(10, 1.0f, float(pi) / 2.0f);
}
