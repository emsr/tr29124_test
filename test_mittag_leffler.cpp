/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_mittag_leffler test_mittag_leffler.cpp -lquadmath
./test_mittag_leffler > test_mittag_leffler.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_mittag_leffler test_mittag_leffler.cpp -lquadmath
./test_mittag_leffler > test_mittag_leffler.txt
*/

#include <iostream>
#include <iomanip>
#include <complex>
#include <ext/cmath>
#include <bits/numeric_limits.h>
#include <bits/float128_io.h>
#include <bits/specfun.h>

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
    __mittag_leffler(_Tp __alpha, _Tp __beta, std::complex<_Tp> __z)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__alpha);
      const auto _S_2pi = __gnu_cxx::__const_2_pi(__alpha);
      const auto _S_pi = __gnu_cxx::__const_pi(__alpha);

      const auto __az = std::abs(__z);
      if (__alpha > _Tp{1})
	{
          unsigned int __k0 = 1u + std::floor(__alpha);
	  auto __alpha0 = __alpha / __k0;
	  auto __rho0 = std::pow(__z, _Tp{1} / _Tp(__k0));

	  auto __E = _Cmplx{0};
	  for (auto __k = 0u; __k < __k0; ++__k)
	    {
	      auto __zk = __rho0
	      		* std::polar(_Tp{1}, _S_2pi * _Tp(__k) / _Tp(__k0));
	      __E += __mittag_leffler(__alpha0, __beta, __zk);
	    }
	  return __E / _Tp(__k0);
	}
      else if (__az < _S_eps)
	return _Tp{1} * std::__detail::__gamma_reciprocal(_Cmplx(__beta));
      else if (__az < _Tp{1})
	{
	  unsigned int __k0 = std::max(std::ceil((_Tp{1} - __beta) / __alpha),
				std::ceil(std::log(_S_eps * (_Tp{1} - __az))
					    / std::log(__az)));
	  auto __E = _Cmplx{0};
	  auto __zk = _Cmplx{1};
	  for (auto __k = 0u; __k <= __k0; ++__k)
	    {
	      //__E += __zk / __gnu_cxx::tgamma(_Cmplx(__beta + __alpha * __k));
	      __E += __zk * std::__detail::__gamma_reciprocal(_Cmplx(__beta + __alpha * __k));
	      __zk *= __z;
	    }
	  return __E;
	}
      else if (__az > std::floor(_Tp{10} + _Tp{5} * __alpha))
	{
	  unsigned int __k0 = std::floor(-std::log(_S_eps) / std::log(__az));
	  auto __E = _Cmplx{0};
	  auto __zk = _Cmplx{1};
	  for (auto __k = 1u; __k <= __k0; ++__k)
	    {
	      __zk /= __z;
	      __E += __zk * std::__detail::__gamma_reciprocal(_Cmplx(__beta - __alpha * __k));
	    }
	  if (std::arg(__z)
	      < _S_pi * (__alpha / _Tp{4} + std::min(_Tp{1}, __alpha) / _Tp{2}))
	    {
	      auto __zpoa = std::pow(__z, _Tp{1} / __alpha);
	      return std::pow(__zpoa, _Tp{1} - __beta)
		   * std::exp(__zpoa) / __alpha;
		   - __E;
	    }
	  else
	    return -__E;
	}
      else
	{
	  return 0;
	}
    }

  template<typename _Tp>
    std::complex<_Tp>
    __mittag_leffler_deriv(_Tp __alpha, _Tp __beta, std::complex<_Tp> __z)
    {
      using _Cmplx = std::complex<_Tp>;
      const auto _S_eps = __gnu_cxx::__epsilon(__alpha);

      const auto __az = std::abs(__z);
      if (__az < _Tp{1})
	{
	  unsigned int __k1 = 0u;
	  if (__alpha > _Tp{1})
	    __k1 = 1 + (_Tp{2} - __alpha - __beta) / (__alpha - _Tp{1});
	  else
	    {
	      auto __D = _Tp{1}
		       + __alpha * (__alpha - 4 * __beta + 6);
	      auto __omega = __alpha + __beta - _Tp{1.5Q};
	      if (__D <= _Tp{0})
		__k1 = 1 + (3 - __alpha - __beta) / __alpha;
	      else
		__k1 = std::max(1 + (3 - __alpha - __beta) / __alpha,
				1 + (1 - 2 * __omega * __alpha + std::sqrt(__D))
				  / (2 * __alpha * __alpha));
	    }
	  unsigned int __k0 = std::max(__k1,
		static_cast<unsigned int>(std::log(_S_eps * (_Tp{1} - __az))
					/ std::log(__az)));
	  auto __Ep = _Cmplx{0};
	  auto __zk = _Cmplx{1};
	  for (auto __k = 0u; __k <= __k0; ++__k)
	    {
	      //__Ep += _Tp(__k + 1) * __zk
		//    / __gnu_cxx::tgamma(_Cmplx(__beta + __alpha * (__k + 1)));
	      __Ep += _Tp(__k + 1) * __zk
		     * std::__detail::__gamma_reciprocal(_Cmplx(__beta + __alpha * (__k + 1)));
	      __zk *= __z;
	    }
	  return __Ep;
	}
      else
	return (__mittag_leffler(__alpha, __beta - 1, __z)
	      - (__beta - 1) * __mittag_leffler(__alpha, __beta, __z))
	     / __alpha / __z;
    }

template<typename _Tp>
  void
  test_mittag_leffler(_Tp proto = _Tp{})
  {
    using namespace std::literals::complex_literals;
    using _Cmplx = std::complex<_Tp>;

    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    {
      auto __alpha = _Tp{0.25Q};
      auto __beta = _Tp{1.00Q};
      std::cout << '\n';
      for (int i = 0; i < 100; ++i)
	{
	  auto t = i * _Tp{0.1Q};
	  auto ml_val = __mittag_leffler(__alpha, __beta, _Cmplx(-t, 0));
	  auto ml_der = -__mittag_leffler_deriv(__alpha, __beta, _Cmplx(-t, 0));
	  std::cout << std::setw(width) << t
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::real(ml_der)
		    << '\n';
	}
    }

    {
      auto __alpha = _Tp{1.75Q};
      auto __beta = _Tp{1.00Q};
      std::cout << '\n';
      for (int i = 0; i < 500; ++i)
	{
	  auto t = i * _Tp{0.1Q};
	  auto ml_val = __mittag_leffler(__alpha, __beta, _Cmplx(-t, 0));
	  auto ml_der = -__mittag_leffler_deriv(__alpha, __beta, _Cmplx(-t, 0));
	  std::cout << std::setw(width) << t
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::real(ml_der)
		    << '\n';
	}
    }

    {
      auto __alpha = _Tp{2.25Q};
      auto __beta = _Tp{1.00Q};
      std::cout << '\n';
      for (int i = 0; i < 1000; ++i)
	{
	  auto t = i * _Tp{0.1Q};
	  auto ml_val = __mittag_leffler(__alpha, __beta, _Cmplx(-t, 0));
	  auto ml_der = -__mittag_leffler_deriv(__alpha, __beta, _Cmplx(-t, 0));
	  std::cout << std::setw(width) << t
		    << std::setw(width) << std::real(ml_val)
		    << std::setw(width) << std::real(ml_der)
		    << '\n';
	}
    }
  }

int
main()
{
  test_mittag_leffler(1.0);
}
