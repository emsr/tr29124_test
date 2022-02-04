/**
 *
 */

#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
#include <complex>

template<typename Tp>
  struct
  meixner_pollaczek_t
  {
    std::complex<Tp> value;
    std::complex<Tp> factor;
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
template<typename Tp, typename _TpX>
  meixner_pollaczek_t<Tp>
  meixner_pollaczek_recur(int n, Tp lambda, Tp phi, _TpX x)
  {
    std::complex<Tp> Pnm1 = Tp{1};
    if (n == 0)
      return {Pnm1, Tp{1}};

    constexpr std::complex<Tp> i{0, 1};
    const auto ix = i * Tp(x);
    const auto sinphi = std::sin(phi);
    const auto cosphi = std::cos(phi);

    auto fact = std::complex<Tp>{1};

    std::complex<Tp> Pn = Tp{2} * lambda * std::polar(Tp{1}, phi)
			 - Tp{2} * i * (lambda + ix) * std::sin(phi);

    if (n == 1)
      return {Pn, fact};

    auto Pnp1 = (Tp{2} * (Tp(x) * sinphi + (lambda + Tp{1}) * cosphi) * Pn
		 - Tp{2} * lambda * Pnm1) / Tp{2};

    for (int k = 2; k < n; ++k)
      {
	Pnm1 = Pn;
	Pn = Pnp1;
	Pnp1 = (Tp{2} * (Tp(x) * sinphi + (lambda + Tp(k)) * cosphi) * Pn
		 - (Tp{2} * lambda + Tp(k - 1)) * Pnm1) / Tp(k + 1);
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
template<typename Tp, typename _TpX>
  meixner_pollaczek_t<Tp>
  meixner_pollaczek(int n, Tp lambda, Tp phi, _TpX x)
  {
    if (std::isnan(lambda))
      return {std::complex<Tp>(lambda), std::complex<Tp>{}};
    else if (std::isnan(phi))
      return {std::complex<Tp>(phi), std::complex<Tp>{}};
    else if (std::isnan(x))
      return {std::complex<Tp>(x), std::complex<Tp>{}};
    else
      return meixner_pollaczek_recur(n, lambda, phi, x);
  }

template<typename Tp>
  void
  test_meixner_pollaczek(int n_max, Tp lambda, Tp phi)
  {
    std::cout.precision(std::numeric_limits<Tp>::digits10);
    auto w = std::cout.precision() + 8;

    for (int n = 0; n <= n_max; ++n)
      {
	std::cout << '\n' << '\n' << " n = " << n << '\n';
	for (int i = -320; i <= 320; ++i)
	  {
	    auto x = i * Tp{0.01L};
	    auto P = meixner_pollaczek(n, lambda, phi, x);
	    std::cout << ' ' << std::setw(w) << x
		      << ' ' << std::setw(w) << std::real(P.value)
		      << ' ' << std::setw(w) << std::imag(P.value)
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
