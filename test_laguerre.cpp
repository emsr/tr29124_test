/*
$HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_laguerre test_laguerre.cpp -lquadmath
./test_laguerre > test_laguerre.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_laguerre test_laguerre.cpp -lquadmath
./test_laguerre > test_laguerre.txt
*/

#include <iostream>
#include <iomanip>
#include <tr1/cmath>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>

  /**
   * Compute the Laguerre polynomial ratio by continued fraction:
   * @f[
   *    \frac{L_n^{(\alpha)}(x)}{L'_n^{(\alpha)(x)} = \frac{x}{n-}
   *       \frac{(n+\alpha)n}{2n+\alpha-1-x-}
   *       \frac{(n-1+\alpha)(n-1)}{2n+\alpha-3-x-}
   *       \frac{(n-2+\alpha)(n-2)}{2n+\alpha-5-x-} ...
   *       \frac{1+\alpha}{1+\alpha-x}
   * @f]
   *
   * @see RICHARD J. MATHAR, GAUSS-LAGUERRE AND GAUSS-HERMITE QUADRATURE
   * ON 64, 96 AND 128 NODES
   *
   * @see T. S. Shao, T. C. Chen, and R. M. Frank,
   * Tables of zeros and Gaussian weights of certain associated Laguerre
   * polynomials and the related generalized Hermite polynomials,
   * Math. Comp. 18 (1964), no. 88, 598{616. MR 0166397 (29 #3674)
   */
  template<typename _Tp>
    _Tp
    __laguerre_ratio(unsigned int __n, _Tp __x)
    {
      auto __frac = _Tp{0};
      for (auto __i = __n - 1; __i >= 0; --__i)
	__frac = _Tp(__n - __i) * _Tp(__n - __i)
		/ (_Tp(2 * (__n - __i) - 1) - __x - __frac);
      __frac = __x / (__n - __frac);
      return __frac;
    }

  /**
   * Compute the Laguerre polynomial ratio by continued fraction:
   * @f[
   *    \frac{L_n^{(\alpha)}(x)}{L'_n^{(\alpha)(x)} = \frac{x}{n-}
   *       \frac{(n+\alpha)n}{2n+\alpha-1-x-}
   *       \frac{(n-1+\alpha)(n-1)}{2n+\alpha-3-x-}
   *       \frac{(n-2+\alpha)(n-2)}{2n+\alpha-5-x-} \cdots
   *       \frac{1+\alpha}{1+\alpha-x}
   * @f]
   *
   * @see RICHARD J. MATHAR, GAUSS-LAGUERRE AND GAUSS-HERMITE QUADRATURE
   * ON 64, 96 AND 128 NODES
   *
   * @see T. S. Shao, T. C. Chen, and R. M. Frank,
   * Tables of zeros and Gaussian weights of certain associated Laguerre
   * polynomials and the related generalized Hermite polynomials,
   * Math. Comp. 18 (1964), no. 88, 598{616. MR 0166397 (29 #3674)
   */
  template<typename _Tp, typename _Ta>
    _Tp
    __laguerre_ratio(unsigned int __n, _Ta __alpha, _Tp __x)
    {
      auto __frac = _Tp{0};
      for (auto __i = __n - 1; __i >= 0; --__i)
	__frac = _Tp(__n - __i + __alpha) * _Tp(__n - __i)
		/ (_Tp(2 * (__n - __i) - 1 + __alpha) - __x - __frac);
      __frac = __x / (__n - __frac);
      return __frac;
    }

  /**
   * Return a vector of zeros of the Laguerre polynomial of degree n.
   */
  template<typename _Tp>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __laguerre_zeros(unsigned int __n, _Tp __proto)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__proto);
      const unsigned int _S_maxit = 1000u;

      std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __pt(__n);

      auto __z = _Tp{0};
      auto __w = _Tp{0};
      for (auto __i = 1u; __i <= __n; ++__i)
	{
	  // Clever approximations for roots.
	  if (__i == 1)
	    __z += (3.0) / (1.0 + 2.4 * __n);
	  else if (__i == 2)
	    __z += (15.0) / (1.0 + 2.5 * __n);
	  else
	    {
	      auto __ai = __i - 2;
	      __z += ((1.0 + 2.55 * __ai) / (1.9 * __ai))
		   * (__z - __pt[__i - 3].__zero);
	    }
	  // Iterate TTRR for polynomial values
	  for (auto __its = 1u; __its <= _S_maxit; ++__its)
	    {
	      auto __L2 = _Tp{0};
	      auto __L1 = _Tp{1};
	      for (auto __j = 1u; __j <= __n; ++__j)
		{
		  auto __L3 = __L2;
		  __L2 = __L1;
		  __L1 = ((_Tp(2 * __j - 1) - __z) * __L2
		       - (_Tp(__j - 1)) * __L3) / _Tp(__j);
		}
	      // Derivative.
	      auto __Lp = (_Tp(__n) * __L1 - _Tp(__n) * __L2) / __z;
	      // Newton's rule for root.
	      auto __z1 = __z;
	      __z = __z1 - __L1 / __Lp;
	      if (std::abs(__z - __z1) <= _S_eps)
		{
		  __w = _Tp{-1} / (__Lp * __n * __L2);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__laguerre_zeros: "
					 "Too many iterations");
	   }
	  __pt[__i - 1].__zero = __z;
	  __pt[__i - 1].__weight = __w;
	}
      return __pt;
    }

  /**
   * Return an array of abscissae and weights for the Gauss-Laguerre rule.
   */
  template<typename _Tp, typename _Tn>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __laguerre_zeros(unsigned int __n, _Tn __alpha, _Tp __proto)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__proto);
      const unsigned int _S_maxit = 1000;

      std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __pt(__n);

      for (auto __i = 1u; __i <= __n; ++__i)
	{
	  auto __z = _Tp{0};
	  auto __w = _Tp{0};
	  // Clever approximations for roots.
	  if (__i == 1)
	    __z += (1.0 + __alpha)
		 * (3.0 + 0.92 * __alpha) / (1.0 + 2.4 * __n + 1.8 * __alpha);
	  else if (__i == 2)
	    __z += (15.0 + 6.25 * __alpha) / (1.0 + 2.5 * __n + 0.9 * __alpha);
	  else
	    {
	      auto __ai = __i - 2;
	      __z += ((1.0 + 2.55 * __ai) / (1.9 * __ai)
		     + 1.26 * __ai * __alpha / (1.0 + 3.5 * __ai))
		   * (__z - __pt[__i - 3].__zero) / (1.0 + 0.3 * __alpha);
	    }
	  // Iterate TTRR for polynomial values
	  for (auto __its = 1u; __its <= _S_maxit; ++__its)
	    {
	      auto __L2 = _Tp{0};
	      auto __L1 = _Tp{1};
	      for (auto __j = 1u; __j <= __n; ++__j)
		{
		  auto __L3 = __L2;
		  __L2 = __L1;
		  __L1 = ((_Tp(2 * __j - 1 + __alpha) - __z) * __L2
			- (_Tp(__j - 1 + __alpha)) * __L3) / _Tp(__j);
		}
	      // Derivative.
	      auto __Lp = (_Tp(__n) * __L1 - _Tp(__n + __alpha) * __L2) / __z;
	      // Newton's rule for root.
	      auto __z1 = __z;
	      __z = __z1 - __L1 / __Lp;
	      if (std::abs(__z - __z1) <= _S_eps)
		{
		  auto __exparg = std::lgamma(_Tp(__alpha + __n))
				- std::lgamma(_Tp(__n));
		  __w = -std::exp(__exparg) / (__Lp * __n * __L2);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__laguerre_zeros: "
					 "Too many iterations");
	   }
	  __pt[__i - 1].__zero = __z;
	  __pt[__i - 1].__weight = __w;
	}
    return __pt;
  }

template<typename _Tp>
  void
  test_laguerre(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    for (int n = 0; n <= 50; ++n)
      {
	auto pt = __laguerre_zeros(n, proto);
	std::cout << "\nl = " << std::setw(4) << n << ":\n";
	for (auto [z, w] : pt)
	  std::cout << ' ' << std::setw(width) << z
		    << ' ' << std::setw(width) << w
		    << '\n';
      }
  }

int
main()
{
  test_laguerre(1.0F);

  test_laguerre(1.0);

  test_laguerre(1.0L);

  test_laguerre(1.0Q);

  return 0;
}

