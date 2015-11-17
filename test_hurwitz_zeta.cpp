// $HOME/bin/bin/g++ -o test_hurwitz_zeta test_hurwitz_zeta.cpp

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>

  template<typename _Tp>
    _Tp
    __hurwitz_zeta_euler_maclaurin(_Tp __s, _Tp __a)
    {
      constexpr auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      constexpr int _S_N = 10 + std::numeric_limits<_Tp>::digits10 / 2;
      constexpr int _S_jmax = 12;
      const auto __pmax  = std::pow(_Tp(_S_N) + __a, -__s);
      auto __sfact = __s;
      auto __pfact = __pmax / (_S_N + __a);
      auto __ans = __pmax * ((_S_N + __a) / (__s - _Tp{1}) + _Tp{0.5L});

      //  Coefficients for Euler-Maclaurin summation:
      //  B_{2j}/(2j)!
      static constexpr _Tp
      _S_hzeta_c[15]
      {
	1.0000000000000000000000000000L,
	8.3333333333333333333333333333e-02L,
       -1.3888888888888888888888888889e-03L,
	3.3068783068783068783068783069e-05L,
       -8.2671957671957671957671957672e-07L,
	2.0876756987868098979210090321e-08L,
       -5.2841901386874931848476822022e-10L,
	1.3382536530684678832826980975e-11L,
       -3.3896802963225828668301953912e-13L,
	8.5860620562778445641359054504e-15L,
       -2.1748686985580618730415164239e-16L,
	5.5090028283602295152026526089e-18L,
       -1.3954464685812523340707686264e-19L,
	3.5347070396294674716932299778e-21L,
       -8.9535174270375468504026113181e-23L
      };

      for(auto __k = 0; __k < _S_N; ++__k)
        __ans += std::pow(__k + __a, -__s);

      for(auto __j = 0; __j <= _S_jmax; ++__j)
        {
	  auto __delta = _S_hzeta_c[__j + 1] * __sfact * __pfact;
	  __ans += __delta;
	  if(std::abs(__delta / __ans) < _Tp{0.5L} * _S_eps)
	    break;
	  __sfact *= (__s + _Tp(2 * __j + 1)) * (__s + _Tp(2 * __j + 2));
	  __pfact /= (_S_N + __a) * (_S_N + __a);
        }

      return __ans;
    }

  template<typename _Tp>
    _Tp
    __hurwitz_zeta_glob(_Tp __s, _Tp __a)
    {
  std::cout.precision(16);
      constexpr auto _S_eps = 1000 * std::numeric_limits<_Tp>::epsilon();
      //  Max before overflow?
      constexpr auto _S_max = std::numeric_limits<_Tp>::max();
      constexpr auto _S_inf = std::numeric_limits<_Tp>::infinity();

      if (__s == +_Tp{0})
	return _S_inf;

      constexpr unsigned int _S_maxit = 10000;
       //  Zeroth order contribution already calculated.
      auto __zeta = _Tp{0.5L};
      for (unsigned int __n = 1; __n < _S_maxit; ++__n)
	{
	  bool __punt = false;
	  auto __term = _Tp{1}; // Again, the zeroth order.
	  auto __bincoeff = _Tp{1};
	  for (unsigned int __k = 1; __k <= __n; ++__k)
	    {
	      __bincoeff *= -_Tp(__n - __k + 1) / _Tp(__k);
	      if (std::abs(__bincoeff) > _S_max)
	      {
		__punt = true;
		break;
	      }
	      __term += __bincoeff * std::pow(_Tp(__k + __a), _Tp{1} - __s);
std::cout << "        " << __k << ' ' << __bincoeff << ' ' << __term << '\n';
	    }
	  if (__punt)
	    break;
	  __term /= _Tp(__n + 1);
	  __zeta += __term;
	  if (std::abs(__term / __zeta) < _S_eps)
	    break;
std::cout << "    " << __n << ' ' << __term << ' ' << __zeta << '\n';
	}

      __zeta /= __s - _Tp{1};

      return __zeta;
    }

int
main()
{
  std::cout.precision(16);
  for (auto ia = 0; ia < 100; ++ia)
    {
      long double a = 0.1 * ia;
      std::cout << ' ' << a << '\n';
      for (auto is = 0; is < 100; ++is)
	{
	  long double s = 0.1 * is;
	  std::cout << ' ' << std::setw(4) << s
                    << ' ' << std::setw(20) << __hurwitz_zeta_euler_maclaurin(s, a)
                    << ' ' << std::setw(20) << __hurwitz_zeta_glob(s, a) << '\n';
	}
    }
}

