/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_anger_weber test_anger_weber.cpp -lquadmath
./test_anger_weber > test_anger_weber.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_anger_weber test_anger_weber.cpp -lquadmath
./test_anger_weber > test_anger_weber.txt
*/

#include <ext/cmath>

  template<typename _Tp>
    struct _AngerWeberState
    {
      _Tp __nu;
      _Tp __z;
      _Tp __Jbold;
      _Tp __Ebold;
    };

  template<typename _Tp>
    _AngerWeberState<_Tp>
    __anger_weber_sum(_Tp __nu, _Tp __z)
    {
      //using _Val = _Tp;
      //using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_eps = __gnu_cxx::__epsilon(__z);

      const auto __z2 = __z / _Tp{2};
      auto _GamArg11 = _Tp{1} + __nu / _Tp{2};
      auto _GamArg12 = _Tp{1} - __nu / _Tp{2};
      auto _GamArg21 = _Tp{3} / _Tp{2} + __nu / _Tp{2};
      auto _GamArg22 = _Tp{3} / _Tp{2} - __nu / _Tp{2};
      auto _Gam11 = std::tgamma(_GamArg11);
      auto _Gam12 = std::tgamma(_GamArg12);
      auto _Gam21 = std::tgamma(_GamArg21);
      auto _Gam22 = std::tgamma(_GamArg22);
      auto __term1 = _Tp{1} / (_Gam11 * _Gam12);
      auto _S1 = __term1;
      auto __term2 = __z2 / (_Gam21 * _Gam22);
      auto _S2 = __term2;
      for (auto __k = 1u; __k < 10000u; ++__k)
	{
	  __term1 *= -__z2 / _GamArg11 * __z2 / _GamArg12;
	  _S1 += __term1;
	  _GamArg11 += _Tp{1};
	  _GamArg12 += _Tp{1};

	  __term2 *= -__z2 / _GamArg21 * __z2 / _GamArg22;
	  _S2 += __term2;
	  _GamArg21 += _Tp{1};
	  _GamArg22 += _Tp{1};

	  if (std::abs(__term1) < _S_eps * std::abs(_S1)
	   && std::abs(__term2) < _S_eps * std::abs(_S2))
	    break;
	}
      //auto [__sin, __cos] = __sincos_pi(__nu / _Tp{2});
      auto __ph = std::__detail::__sincos_pi(__nu / _Tp{2});
      return _AngerWeberState<_Tp>{__nu, __z,
				   __ph.__cos_value * _S1
				 + __ph.__sin_value * _S2,
				   __ph.__sin_value * _S1
				 - __ph.__cos_value * _S2};
    }

  /**
   * Compute Anger and Weber functions for fixed order @f$ \nu @f$
   * and large agument @f$ |z| @f$.
   *
   * @see http://dlmf.nist.gov/11.11#i
   */
  template<typename _Tp>
    _AngerWeberState<_Tp>
    __anger_weber_asymp_arg(_Tp __nu, _Tp __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__z));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      const auto _S_max_iter = 1000u;
      const auto __z2 = __z * __z;

      auto __F_z2k = _Tp{1};
      auto __G_z2k = _Tp{1};

      auto __Fsum = __F_z2k;
      auto __Gsum = __G_z2k;
      for (auto __k = 1u; __k < _S_max_iter; ++__k)
	{
	  __F_z2k *= (__nu - _Tp(2 * __k - 1)) * (__nu + _Tp(2 * __k - 1))
		   / __z2;
	  __Fsum += __F_z2k;
	  __G_z2k *= (__nu - _Tp(2 * __k)) * (__nu + _Tp(2 * __k))
		   / __z2;
	  __Gsum += __G_z2k;
	}

      auto __ph = std::__detail::__sincos_pi(__nu / _Tp{2});
      auto _Bess = __cyl_bessel(__nu, __z);
      return _AngerWeberState<_Tp>{__nu, __z,
				   _Bess._J_value
				     + __ph.__sin_value
					* (__Fsum + __nu * __Gsum / __z)
					/ _S_pi / __z,
				  -_Bess._N_value
				     - (_Tp{1} + __ph.__cos_value) * __Fsum
					/ _S_pi / __z
				     - (_Tp{1} - __ph.__cos_value) * __Gsum
					* __nu / _S_pi / __z / __z};
    }

  /**
   * Compute Anger and Weber functions for large order @f$ |\nu| @f$
   * and fixed agument @f$ z @f$.
   *
   * @see http://dlmf.nist.gov/11.11#ii
   */
  template<typename _Tp>
    _AngerWeberState<_Tp>
    __anger_weber_asymp_order(_Tp __nu, _Tp __z)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      const auto __sinnp = __gnu_cxx::sin_pi(__nu);
      const auto __sinnpd2 = __gnu_cxx::sin_pi(__nu / _Tp{2});
      const auto __cosnpd2 = __gnu_cxx::cos_pi(__nu / _Tp{2});
      const auto __nufact = __nu * __z / (__nu * __nu - _Tp{1});
      return _AngerWeberState<_Tp>{__nu, __z,
				   __sinnp * (_Tp{1} - __nufact) / __nu / _S_pi,
				   _Tp{2} * (__sinnpd2 + __nufact * __cosnpd2)
					  / __nu / _S_pi};
    }

  /**
   * Compute Anger and Weber functions for large order @f$ \nu @f$
   * and fixed ratio @f$ z/\nu @f$.
   *
   * @see http://dlmf.nist.gov/11.11#iii
   */
  template<typename _Tp>
    _AngerWeberState<_Tp>
    __anger_weber_asymp_uniform(_Tp __nu, _Tp __z)
    {
    }

  /* Use the reciprocal gamma function... Fails.  WTF. */
  template<typename _Tp>
    _AngerWeberState<_Tp>
    __anger_weber_sum_recip(_Tp __nu, _Tp __z)
    {
      //using _Val = _Tp;
      //using _Real = std::__detail::__num_traits_t<_Val>;
      const auto _S_eps = __gnu_cxx::__epsilon(__z);

      const auto __z2 = __z / _Tp{2};
      auto _GamArg11 = _Tp{1} + __nu / _Tp{2};
      auto _GamArg12 = _Tp{1} - __nu / _Tp{2};
      auto _GamArg21 = _Tp{3} / _Tp{2} + __nu / _Tp{2};
      auto _GamArg22 = _Tp{3} / _Tp{2} - __nu / _Tp{2};
      auto _Gam11 = std::__detail::__gamma_reciprocal(_GamArg11);
      auto _Gam12 = std::__detail::__gamma_reciprocal(_GamArg12);
      auto _Gam21 = std::__detail::__gamma_reciprocal(_GamArg21);
      auto _Gam22 = std::__detail::__gamma_reciprocal(_GamArg22);
      auto __term1 = _Gam11 * _Gam12;
      auto _S1 = __term1;
      auto __term2 = __z2 * _Gam21 * _Gam22;
      auto _S2 = __term2;
      for (auto __k = 1u; __k < 10000u; ++__k)
	{
	  __term1 *= -__z2 / _GamArg11 * __z2 / _GamArg12;
	  _S1 += __term1;
	  _GamArg11 += _Tp{1};
	  _GamArg12 += _Tp{1};

	  __term2 *= -__z2 / _GamArg21 * __z2 / _GamArg22;
	  _S2 += __term2;
	  _GamArg21 += _Tp{1};
	  _GamArg22 += _Tp{1};

	  if (std::abs(__term1) < _S_eps * std::abs(_S1)
	   && std::abs(__term2) < _S_eps * std::abs(_S2))
	    break;
	}
      //auto [__sin, __cos] = __sincos_pi(__nu / _Tp{2});
      auto __ph = std::__detail::__sincos_pi(__nu / _Tp{2});
      return _AngerWeberState<_Tp>{__nu, __z,
				   __ph.__cos_value * _S1
				 + __ph.__sin_value * _S2,
				   __ph.__sin_value * _S1
				 - __ph.__cos_value * _S2};
    }


template<typename _Tp>
  void
  test_anger_weber(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    for (auto nu : {_Tp{0}, _Tp{0.5Q}, _Tp{1}, _Tp{1.5Q}, _Tp{2}, _Tp{5}})
      {
	std::cout << "\n\n nu = " << std::setw(4) << nu << '\n';
	std::cout << ' ' << std::setw(4) << "z"
		  << ' ' << std::setw(width) << "Jbold"
		  << ' ' << std::setw(width) << "Ebold"
		  << '\n';
	std::cout << ' ' << std::setw(4) << "-"
		  << ' ' << std::setw(width) << "-----"
		  << ' ' << std::setw(width) << "-----"
		  << '\n';
	for (int k = -80; k <= 80; ++k)
	  {
	    auto z = _Tp{0.1Q} * k;
	    auto AW = __anger_weber_sum(nu, z);
	    std::cout << ' ' << std::setw(4) << z
		      << ' ' << std::setw(width) << AW.__Jbold
		      << ' ' << std::setw(width) << AW.__Ebold
		      << '\n';
	  }
      }
  }

int
main()
{
  test_anger_weber(1.0);
}
