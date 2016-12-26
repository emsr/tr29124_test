/*
$HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_legendre test_legendre.cpp -lquadmath
./test_legendre > test_legendre.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_legendre test_legendre.cpp -lquadmath
./test_legendre > test_legendre.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <stdexcept>
#include <limits>
#include <vector>

#include "legendre.tcc"
#include "bits/sf_legendre.tcc"

  /**
   * Build a list of zeros and weights for the Gauss-Legendre integration rule
   * for the Legendre polynomial of degree @c l.
   */
  template<typename _Tp>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __legendre_zeros(unsigned int __l, _Tp proto = _Tp{})
    {
      const auto _S_eps = __gnu_cxx::__epsilon(proto);
      const auto _S_pi = __gnu_cxx::__const_pi(proto);
      const unsigned int _S_maxit = 1000u;

      std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __pt(__l);

      auto __m = __l / 2;

      // Treat the central zero for odd order specially.
      // Be careful to avoid overflow of the factorials.
      // An alternative would be to proceed with the recursion
      // for large order.
      if (__l & 1)
	{
	  if (__l < std::__detail::_S_num_factorials<_Tp>)
	    {
	      auto __lm = __l - 1;
	      auto __lmfact = std::__detail::__factorial<_Tp>(__lm);
	      auto __mm = __lm / 2;
	      auto __mmfact = std::__detail::__factorial<_Tp>(__mm);
	      auto __Plm1 = (__lm & 1 ? -1 : 1) * __lmfact / __mmfact / __mmfact
			    / std::pow(_Tp{2}, __lm);
	      auto __Ppl = __l * __Plm1;
	      __pt[__m].__zero = _Tp{0};
	      __pt[__m].__weight = _Tp{2} / __Ppl / __Ppl;
	    }
	  else
	    {
	      auto __lm = __l - 1;
	      auto __lmfact = std::__detail::__log_factorial<_Tp>(__lm);
	      auto __mm = __lm / 2;
	      auto __mmfact = std::__detail::__log_factorial<_Tp>(__mm);
	      auto __Plm1 = (__lm & 1 ? -1 : 1)
			  * std::exp(__lmfact - 2 * __mmfact)
			  / std::pow(_Tp{2}, __lm);
	      auto __Ppl = __l * __Plm1;
	      __pt[__m].__zero = _Tp{0};
	      __pt[__m].__weight = _Tp{2} / __Ppl / __Ppl;
	    }
	}

      for (auto __i = 1u; __i <= __m; ++__i)
	{
	  // Clever approximation of root.
	  auto __z = std::cos(_S_pi * (__i - _Tp{1} / _Tp{4})
				    / (__l + _Tp{1} / _Tp{2}));
	  auto __z1 = __z;
	  auto __w = _Tp{0};
	  for (auto __its = 0u; __its < _S_maxit; ++__its)
	    {
	      // Compute __P, __P1, and __P2 the Legendre polynomials of order
	      // l, l-1, l-2 respectively by iterating through the recursion
	      // relation for the Legendre polynomials.
	      // Compute __Pp the derivative of the Legendre polynomial of order l.
	      auto __P1 = _Tp{0};
	      auto __P = _Tp{1};
	      for  (auto __k = 1u; __k <= __l; ++__k)
		{
		  auto __P2 = __P1;
		  __P1 = __P;
		  // Recursion for Legendre polynomials.
		  __P = ((_Tp{2} * __k - _Tp{1}) * __z * __P1
		      - (__k - _Tp{1}) * __P2) / __k;
		}
	      // Recursion for the derivative of The Legendre polynomial.
	      auto __Pp = __l * (__z * __P - __P1) / (__z * __z - _Tp{1});
	      __z1 = __z;
	      // Converge on root by Newton's method.
	      __z = __z1 - __P / __Pp;
	      if (std::abs(__z - __z1) < _S_eps)
		{
		  __w = _Tp{2} / ((_Tp{1} - __z * __z) * __Pp * __Pp);
		  break;
		}
	      if (__its > _S_maxit)
		std::__throw_logic_error("__legendre_zeros: "
					 "Too many iterations");
	    }

	  __pt[__i - 1].__zero = -__z;
	  __pt[__l - __i].__zero = __z;
	  __pt[__i - 1].__weight = __w;
	  __pt[__l - __i].__weight = __w;
	}

      return __pt;
    }

template<typename _Tp>
  void
  test_legendre(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    for (int l = 0; l <= 50; ++l)
      {
	std::cout << "  " << std::setw(width) << "x";
	std::cout << "  " << std::setw(width) << "P_" << l << "(x)";
	std::cout << '\n';
	for (int i = -100; i <= 100; ++i)
	  {
	    auto x = i * _Tp{0.01Q};
	    _Tp P, P_l, P_l0;
	    try
	      {
		P = __legendre_p(l, x);
		P_l = std::__detail::__poly_legendre_p(l, x);
		P_l0 = std::__detail::__assoc_legendre_p(l, 0, x);
	      }
	    catch (std::exception & err)
	      {
		std::cerr << err.what() << '\n';
	      }
	    std::cout << "  " << std::setw(width) << x;
	    std::cout << "  " << std::setw(width) << P;
	    std::cout << "  " << std::setw(width) << P_l;
	    std::cout << "  " << std::setw(width) << P_l0;
	    std::cout << "  " << std::setw(width) << P_l - P;
	    std::cout << '\n';
	  }
      }

    for (unsigned int l = 0; l <= 3; ++l)
      {
	for (unsigned int m = l; m <= l; --m)
	  {
	    std::cout << "  " << std::setw(width) << "x";
	    std::cout << "  " << std::setw(width) << "P_" << l << "_" << m << "(x)";
	    std::cout << '\n';
	    for (int i = -100; i <= 100; ++i)
	      {
		auto x = i * _Tp{0.01Q};
		_Tp P_lm;
		try
		  {
		    P_lm = std::__detail::__assoc_legendre_p(l, m, x);

		  }
		catch (std::exception & err)
		  {
		    std::cerr << err.what() << '\n';
		  }
		std::cout << "  " << std::setw(width) << x;
		std::cout << "  " << std::setw(width) << P_lm;
		std::cout << '\n';
	      }
	  }
      }

    for (int l = 0; l <= 1024; ++l)
      {
	auto pt = __legendre_zeros(l, proto);
	std::cout << "\nl = " << std::setw(4) << l << ":\n";
	for (auto [z, w] : pt)
	  std::cout << ' ' << std::setw(width) << z
		    << ' ' << std::setw(width) << w
		    << '\n';
      }
  }

int
main()
{
  test_legendre(1.0F);

  test_legendre(1.0);

  test_legendre(1.0L);

  test_legendre(1.0Q);

  return 0;
}

