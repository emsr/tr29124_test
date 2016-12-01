/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_hermite test_hermite.cpp -lquadmath
./test_hermite > test_hermite.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_hermite test_hermite.cpp -lquadmath
./test_hermite > test_hermite.txt
*/

#include <sstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <limits>

#include <bits/float128_io.h>
#include "new_hermite.tcc"
#include "LentzContinuedFraction.tcc"

#include "quadrature/integration.h"

  /**
   * Compute the Hermite polynomial ratio by continued fraction:
   * @f[
   *    \frac{H_n(x)}{H_{n-1}(x)} = 2n\frac{H_n(x)}{H'_n(x)}
   *       = b_0 + K_{k=1}^{n-1}\left(\frac{-2(n-k)}{2x}\right)
   *      b_k = 2x, \mbox{ } k >= 0 \mbox{ } a_k = -2(n-k) \mbox{ } k >= 1
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
    __hermite_ratio(unsigned int __n, _Tp __x)
    {
      auto __a = [__n](std::size_t __k, _Tp) { return _Tp(-2 * (__n - __k)); };
      using _AFun = decltype(__a);

      auto __b = [__x](std::size_t, _Tp) { return _Tp{2} * __x; };
      using _BFun = decltype(__b);

      auto __w = [](std::size_t, _Tp) { return _Tp{0}; };
      using _TailFun = decltype(__w);

      using _CFrac = _LentzContinuedFraction<_Tp, _AFun, _BFun, _TailFun>;
      _CFrac __Hrat(__a, __b, __w);

      return __Hrat();
    }

  template<typename _Tp>
    std::vector<_Tp> 
    __hermite_zeros(unsigned int __n, _Tp __proto = _Tp{})
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__proto);
      const unsigned int _S_maxit = 1000u;
      const auto _S_pim4 = _Tp{0.7511255444649424828587030047762276930510L};

      std::vector<_Tp> __zero(__n);
      std::vector<_Tp> __weight(__n);

      auto __m = (__n + 1) / 2;
      _Tp __z;
      _Tp __w;
      for (auto __i = 1u; __i <= __m; ++__i)
	{
	  if (__i == 1)
	    __z = std::sqrt(_Tp(2 * __n + 1))
		- 1.85575 * std::pow(_Tp(2 * __n + 1), -0.166667);
	  else if (__i == 2)
	    __z -= 1.14 * std::pow(_Tp(__n), 0.426) / __z;
	  else if (__i == 3)
	    __z = 1.86 * __z - 0.86 * __zero[0];
	  else if (__i == 4)
	    __z = 1.91 * __z - 0.91 * __zero[1];
	  else
	    __z = 2.0 * __z - __zero[__i - 3];
	  for (auto __its = 1u; __its <= _S_maxit; ++__its)
	    {
	      auto __H1 = _S_pim4;
	      auto __H2 = _Tp{0};
	      for (auto __j = 1u; __j <= __n; ++__j)
		{
		  auto __H3 = __H2;
		  __H2 = __H1;
		  __H1 = __z * std::sqrt(_Tp{2} / __j) * __H2
		       - std::sqrt(_Tp(__j - 1) / _Tp(__j)) * __H3;
		}
	      auto __Hp = std::sqrt(_Tp(2 * __n)) * __H2;
	      auto __z1 = __z;
	      __z = __z1 - __H1 / __Hp;
	      if (std::abs(__z - __z1) <= _S_eps)
		{
		  __w = 2.0 / (__Hp * __Hp);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__hermite_zeros: "
					 "Too many iterations");
	    }
	  __zero[__i - 1] = __z;
	  __zero[__n - __i] = -__z;
	  __weight[__i] = __w;
	  __weight[__n - __i] = __w;
	}

      return __zero;
    }

template<typename _Tp>
  void
  test_hermite(_Tp proto = _Tp{})
  {
    const auto _S_pi = __gnu_cxx::__const_pi(proto);

    auto fname = [](std::string_view front, int n, std::string_view back)
		 {
		   std::ostringstream out;
		   out << front << n << back;
		   return out.str();
		 };

    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout.flags(std::ios::showpoint);
    auto width = 8 + std::cout.precision();

    std::cout << "\f\n\n  Table of integer sqrt\n";
    std::cout << "  =====================\n";
    for (int n = 0; n <= 100; ++n)
      std::cout << "  " << std::setw(4) << n
		<< "  " << std::setw(width) << std::sqrt(_Tp(n)) << '\n';

    std::cout << "\f\n\n  Table of factorial sqrt\n";
    std::cout << "  =======================\n";
    auto fact = _Tp{1};
    for (int n = 0; n <= 100; ++n)
      {
	std::cout << "  " << std::setw(4) << n
		  << "  " << std::setw(width) << fact << '\n';
	fact *= std::sqrt(_Tp(n + 1));
      }

    std::cout << "\f\n\n  Examine asymptotic transition region\n";
    std::cout << "  ====================================\n";
    for (int n = 4; n <= 200; ++n)
      {
	auto xt = std::sqrt(_Tp{2} * n);
	auto del = _Tp{0.5Q} * xt / _Tp{200};
	std::cout << "  " << std::setw(6) << "n   = " << std::setw(4) << n << '\n';
	std::cout << "  " << std::setw(6) << "x_t = " << std::setw(width) << xt << '\n';
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("Hr_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Ha_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta";
	std::cout << "  " << std::setw(width) << "delta-ratio";
	std::cout << "  " << std::setw(width) << fname("Ha2_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta2";
	std::cout << "  " << std::setw(width) << "delta2-ratio";
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	for (int i = 0; i <= 201; ++i)
          {
            auto x = _Tp{0.75Q} * xt + i * del;
            auto h = __poly_hermite_recursion(n, x);
            auto ht = __poly_hermite_asymp(n, x);
            auto h2 = __poly_hermite_asymp2(n, x);
	    if (std::abs(x - xt) < del)
	      std::cout << ">>";
	    else if (std::abs(x - _Tp{0.95Q} * xt) < del)
	      std::cout << "> ";
	    else if (std::abs(x - _Tp{1.05Q} * xt) < del)
	      std::cout << "> ";
            else
              std::cout << "  ";
	    std::cout << std::setw(width) << x
		      << "  " << std::setw(width) << h
		      << "  " << std::setw(width) << ht
		      << "  " << std::setw(width) << ht - h
		      << "  " << std::setw(width) << (ht - h) / h
		      << "  " << std::setw(width) << h2
		      << "  " << std::setw(width) << h2 - h
		      << "  " << std::setw(width) << (h2 - h) / h
		      << '\n';
          }
      }

    std::cout << "\f\n\n  Examine asymptotic transition region\n";
    std::cout << "  ====================================\n";
    for (int n : {1000, 2000, 5000, 10000})
      {
	auto xt = std::sqrt(_Tp{2} * n);
	auto del = _Tp{0.5Q} * xt / _Tp{200};
	std::cout << "  " << std::setw(6) << "n   = " << std::setw(4) << n << '\n';
	std::cout << "  " << std::setw(6) << "x_t = " << std::setw(width) << xt << '\n';
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("Hr_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Ha_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta";
	std::cout << "  " << std::setw(width) << "delta-ratio";
	std::cout << "  " << std::setw(width) << fname("Ha2_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta2";
	std::cout << "  " << std::setw(width) << "delta2-ratio";
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	for (int i = 1; i <= 201; ++i)
          {
            auto x = _Tp{0.75Q} * xt + i * del;
            auto h = __poly_hermite_recursion(n, x);
            auto ht = __poly_hermite_asymp(n, x);
            auto h2 = __poly_hermite_asymp2(n, x);
	    if (std::abs(x - xt) < del)
	      std::cout << ">>";
	    else if (std::abs(x - _Tp{0.95Q} * xt) < del)
	      std::cout << "> ";
	    else if (std::abs(x - _Tp{1.05Q} * xt) < del)
	      std::cout << "> ";
            else
              std::cout << "  ";
	    std::cout << std::setw(width) << x
		      << "  " << std::setw(width) << h
		      << "  " << std::setw(width) << ht
		      << "  " << std::setw(width) << ht - h
		      << "  " << std::setw(width) << (ht - h) / h
		      << "  " << std::setw(width) << h2
		      << "  " << std::setw(width) << h2 - h
		      << "  " << std::setw(width) << (h2 - h) / h
		      << '\n';
          }
      }

    std::cout << "\f\n\n  Compare recursion and asymptotic\n";
    std::cout << "  ================================\n";
    for (int n = 0; n <= 50; ++n)
      {
	auto xt = std::sqrt(_Tp{2} * n);
	std::cout << "  " << std::setw(6) << "n   = " << std::setw(4) << n << '\n';;
	std::cout << "  " << std::setw(6) << "x_t = " << std::setw(width) << xt << '\n';;
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("Hr_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Ha_", n, "(x)");
	std::cout << "  " << std::setw(width) << "delta";
	std::cout << "  " << std::setw(width) << "delta-ratio";
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	auto del = _Tp{0.1Q};
	for (int i = 0; i <= 100; ++i)
          {
            auto x = i * del;
            auto h = __poly_hermite_recursion(n, x);
            auto ht = __poly_hermite_asymp(n, x);
	    if (std::abs(x - xt) < del)
	      std::cout << ">>";
	    else if (std::abs(x - _Tp{0.95Q} * xt) < del)
	      std::cout << "> ";
	    else if (std::abs(x - _Tp{1.05Q} * xt) < del)
	      std::cout << "> ";
            else
              std::cout << "  ";
            std::cout << std::setw(width) << x
		      << "  " << std::setw(width) << h
		      << "  " << std::setw(width) << ht
		      << "  " << std::setw(width) << ht - h
		      << "  " << std::setw(width) << (ht - h) / h
		      << '\n';
          }
      }

    std::cout << "\f\n\n  Compare normalizations\n";
    std::cout << "  ======================\n";
    auto power = _Tp{1};
    for (int n = 0; n <= 50; ++n)
      {
	// sqrt(factorial(n) * 2**n sqrt(pi))
	auto norm = std::exp(_Tp{0.5Q} * std::lgamma(n + 1)) * std::sqrt(power * std::sqrt(_S_pi));
	power *= _Tp{2};
	std::cout << "  " << std::setw(width) << "n = " << n << '\n';
	std::cout << "  " << std::setw(width) << "norm = " << norm << '\n';
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << fname("H_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Hn_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("He_", n, "(x)");
	std::cout << "  " << std::setw(width) << fname("Hen_", n, "(x)");
	std::cout << '\n';
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << "  " << std::setw(width) << "------------";
	std::cout << '\n';
	for (int i = 0; i <= 100; ++i)
          {
            auto x = i * _Tp{0.1Q};
            auto h = __poly_hermite_recursion(n, x);
            auto hn = __poly_hermite_norm_recursion(n, x);
            auto he = __poly_prob_hermite_recursion(n, x);
            auto hen = __poly_prob_hermite_norm_recursion(n, x);
            std::cout << "  " << std::setw(width) << x;
            std::cout << "  " << std::setw(width) << h;
            std::cout << "  " << std::setw(width) << hn;
            std::cout << "  " << std::setw(width) << he;
            std::cout << "  " << std::setw(width) << hen;
            std::cout << '\n';
          }
      }

    for (int n = 0; n <= 50; ++n)
      {
	auto zero = __hermite_zeros(n, proto);
	std::cout << "\nl = " << std::setw(4) << n << ":\n";
	for (auto z : zero)
	  std::cout << ' ' << std::setw(width) << z << '\n';
      }

    return;
  }

template<typename _Tp>
  void
  test_hermite_norm(_Tp proto)
  {
    return;
  }


int
main()
{
  test_hermite(1.0F);

  test_hermite(1.0);

  test_hermite(1.0L);

  test_hermite(1.0Q);

  return 0;
}
