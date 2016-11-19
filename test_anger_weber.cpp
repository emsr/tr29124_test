/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_anger_weber test_anger_weber.cpp -lquadmath
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

      const auto __z2 = __z / _Tp{2};
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));
      auto __z2k = _Tp{1};
      auto _GamArg11 = _Tp{1} + __nu / _Tp{2};
      auto _GamArg12 = _Tp{1} - __nu / _Tp{2};
      auto _GamArg21 = _Tp{3} / _Tp{2} + __nu / _Tp{2};
      auto _GamArg22 = _Tp{3} / _Tp{2} - __nu / _Tp{2};
      auto _Gam11 = std::tgamma(_GamArg11);
      auto _Gam12 = std::tgamma(_GamArg12);
      auto _Gam21 = std::tgamma(_GamArg21);
      auto _Gam22 = std::tgamma(_GamArg22);
      auto __term1 = _Tp{1} / (_Gam11 * _Gam12);
      auto __term2 = _Tp{1} / (_Gam21 * _Gam22);
      auto _S1 = __term1;
      auto _S2 = __term2;
      for (int __k = 1; __k < 100; ++__k)
	{
	  _GamArg11 += _Tp{1};
	  _GamArg12 += _Tp{1};
	  __z2k *= __z2;
	  __term1 *= _Tp{-1} * __z2k / (_GamArg11 * _GamArg12);
	  _S1 += __term1;

	  _GamArg21 += _Tp{1};
	  _GamArg22 += _Tp{1};
	  __z2k *= __z2;
	  __term1 *= _Tp{-1} * __z2k / (_GamArg21 * _GamArg22);
	  _S2 += __term2;
	}
      //auto [__sin, __cos] = __sincos_pi(__nu / _Tp{2});
      auto __ph = __sincos_pi(__nu / _Tp{2});
      return _AngerWeberState<_Tp>{__nu, __z,
				   __cos * _S1, __ph.__sin_value * _S2,
				   __sin * _S1, -__ph.__cos_value * _S2};
    }


int
main()
{
}
