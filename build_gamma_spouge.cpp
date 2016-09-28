/*
$HOME/bin_tr29124/bin/g++ -I. -o build_gamma_spouge build_gamma_spouge.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./build_gamma_spouge > build_gamma_spouge.txt

$HOME/bin/bin/g++ -std=gnu++14 -DNO_LOGBQ -I. -o build_gamma_spouge build_gamma_spouge.cpp -lquadmath
./build_gamma_spouge > build_gamma_spouge.txt
*/

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include <bits/float128.h>

  template<typename _Tp>
    void
    spouge()
    {
      std::cout.precision(std::numeric_limits<_Tp>::digits10);
      std::cout << std::showpoint << std::scientific;
      auto width = 8 + std::cout.precision();

      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_pi  = __gnu_cxx::__math_constants<_Tp>::__pi;
      const auto _S_2pi = _Tp{2} * _S_pi;
      auto __a = _Tp{1};
      const auto __fact = _Tp{1} / std::sqrt(_S_2pi);
      while (_S_eps <= __fact * std::pow(_S_2pi, -__a) / std::sqrt(__a))
	{
          std::cout << "err = " << __fact * std::pow(_S_2pi, -__a) / std::sqrt(__a) << '\n';
          __a += _Tp{1};
	}
      std::cout << "a = " << __a << '\n';

      std::vector<_Tp> __c;
      auto __factc = _Tp{1};
      __c.push_back(__factc * std::sqrt(__a - 1) * std::exp(__a - 1));
      std::cout << "c_0 = " << __c.back() << '\n';
      auto __sum = __factc * std::exp(__a - 1);
      for (int __k = 1; __k < std::ceil(__a); ++__k)
	{
	  __factc *= -_Tp{1} / _Tp(__k);
	  auto __ak = _Tp(__a - __k - 1);
	  __c.push_back(__factc * std::pow(__ak, _Tp(__k + 0.5Q)) * std::exp(__ak));
	  std::cout << "c_" << __k << " = " << __c.back() << '\n';
	}

      auto __log_gamma_spouge =
	[=](_Tp __z)
	-> _Tp
	{
	  // Reflection is right but auto and use of functions won't compile.
	  //if (__z <= -__a)
	  //  return std::log(_S_pi) - std::log(std::sin(_S_pi * __z)) - __log_gamma_spouge(_Tp{1} - __z);
	  //else
	    {
	      auto __sum = std::sqrt(_S_2pi);
	      for (int __k = 0; __k < __c.size(); ++__k)
		__sum += __c[__k] / (__z + __k + 1);
	      return std::log(__sum)
		   + (__z + _Tp{0.5Q}) * std::log(__z + __a)
		   - (__z + __a) - std::log(__z);
	    }
	};

      std::cout << '\n'
		<< ' ' << std::setw(width) << "z"
		<< ' ' << std::setw(width) << "spouge"
		<< ' ' << std::setw(width) << "lgamma"
		<< ' ' << std::setw(width) << "delta"
		<< '\n';
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = _Tp{0.01Q} * i;
	  std::cout << ' ' << std::setw(width) << z
		    << ' ' << std::setw(width) << __log_gamma_spouge(z)
		    << ' ' << std::setw(width) << std::lgamma(z)
		    << ' ' << std::setw(width) << __log_gamma_spouge(z) - std::lgamma(z) << '\n';
	}

      //  Try to invert using Newton...
      auto __log_gamma_inv
      {
	[=](_Tp __y)
	-> _Tp
	{
	  constexpr auto _S_log_10 = __gnu_cxx::__math_constants<_Tp>::__log10_e;
	  constexpr auto _S_ln_2 = __gnu_cxx::__math_constants<_Tp>::__ln_2;
	  constexpr auto _S_ln_pi = __gnu_cxx::__math_constants<_Tp>::__ln_pi;
	  constexpr auto _S_log_sqrt_2pi = (_S_ln_2 + _S_ln_pi) / _Tp{2};
	  auto __x = __y * _S_log_10 - _S_log_sqrt_2pi;
	  auto __x0 = __x;
	  for (int __i = 0; __i < 100; ++__i)
	    __x = (__x0 + __x) / std::log(__x);
	  return __x;
	}
      };

      std::cout << '\n'
		<< ' ' << std::setw(width) << "z"
		<< ' ' << std::setw(width) << "spouge"
		<< ' ' << std::setw(width) << "lgamma"
		<< ' ' << std::setw(width) << "delta"
		<< ' ' << std::setw(width) << "inverse"
		<< ' ' << std::setw(width) << "delta_inv"
		<< '\n';
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = _Tp{0.01Q} * i;
	  _Tp x, y;
	  std::cout << ' ' << std::setw(width) << z
		    << ' ' << std::setw(width) << (y = __log_gamma_spouge(z))
		    << ' ' << std::setw(width) << std::lgamma(z)
		    << ' ' << std::setw(width) << __log_gamma_spouge(z) - std::lgamma(z)
		    << ' ' << std::setw(width) << (x = __log_gamma_inv(y))
		    << ' ' << std::setw(width) << x - z
		    << '\n';
	}
    }

  /**
   * A struct for Spouge algorithm Chebyshev arrays of coefficients.
   */
  template<typename _Tp>
    class _GammaSpouge
    {
    };

  template<>
    class _GammaSpouge<float>
    {
    public:
    private:
      static constexpr std::array<float, 7>
      _S_cheby
      {
	2.901419e+03F,
	-5.929168e+03F,
	4.148274e+03F,
	-1.164761e+03F,
	1.174135e+02F,
	-2.786588e+00F,
	3.775392e-03F,
      };
    };

  template<>
    class _GammaSpouge<double>
    {
    public:
    private:
      static constexpr std::array<double, 18>
      _S_cheby
      {
	2.785716565770350e+08,
	-1.693088166941517e+09,
	4.549688586500031e+09,
	-7.121728036151557e+09,
	7.202572947273274e+09,
	-4.935548868770376e+09,
	2.338187776097503e+09,
	-7.678102458920741e+08,
	1.727524819329867e+08,
	-2.595321377008346e+07,
	2.494811203993971e+06,
	-1.437252641338402e+05,
	4.490767356961276e+03,
	-6.505596924745029e+01,
	3.362323142416327e-01,
	-3.817361443986454e-04,
	3.273137866873352e-08,
	-7.642333165976788e-15,
      };
    };

  template<>
    class _GammaSpouge<long double>
    {
    public:
    private:
      static constexpr std::array<long double, 22>
      _S_cheby
      {
	1.681473171108908244e+10L,
	-1.269150315503303974e+11L,
	4.339449429013039995e+11L,
	-8.893680202692714895e+11L,
	1.218472425867950986e+12L,
	-1.178403473259353616e+12L,
	8.282455311246278274e+11L,
	-4.292112878930625978e+11L,
	1.646988347276488710e+11L,
	-4.661514921989111004e+10L,
	9.619972564515443397e+09L,
	-1.419382551781042824e+09L,
	1.454145470816386107e+08L,
	-9.923020719435758179e+06L,
	4.253557563919127284e+05L,
	-1.053371059784341875e+04L,
	1.332425479537961437e+02L,
	-7.118343974029489132e-01L,
	1.172051640057979518e-03L,
	-3.323940885824119041e-07L,
	4.503801674404338524e-12L,
	-5.320477002211632680e-20L,
      };
    };


#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  template<>
    class _GammaSpouge<__float128>
    {
    public:
    private:
      static constexpr std::array<__float128, 40>
      _S_cheby
      {
	 1.488707141702962349642653219904701e+18Q,
	-2.109024888057172629401888926607933e+19Q,
	 1.417814156004491800291266261243181e+20Q,
	-6.017978256631943800541949222618531e+20Q,
	 1.810303371559858936645834474459794e+21Q,
	-4.106765117367170106628989406117407e+21Q,
	 7.299495565278381248758472830920901e+21Q,
	-1.042658085606638732437326916809544e+22Q,
	 1.218085331069918351049721248713456e+22Q,
	-1.178425136120009375043789561151657e+22Q,
	 9.524552911654596104822221145923981e+21Q,
	-6.470848515795635226546940452396693e+21Q,
	 3.710085524031367371144056238201665e+21Q,
	-1.799208570324458717498947049552296e+21Q,
	 7.385190850868692826830617550673786e+20Q,
	-2.564085902302497902610094974644620e+20Q,
	 7.515077439232198453710859556732366e+19Q,
	-1.853282816124970754507275065198273e+19Q,
	 3.827869971220176147395256318872661e+18Q,
	-6.582119247090799001499961431273749e+17Q,
	 9.351471019350666459087977506154173e+16Q,
	-1.087561099332606693611359806619008e+16Q,
	 1.023682895113200510083891964414548e+15Q,
	-7.692392155137210713863610487671446e+13Q,
	 4.538847752507923706006970038411543e+12Q,
	-2.061098065296257001203592871224284e+11Q,
	 7.029050529288024691234709461995537e+09Q,
	-1.746978910485928115412705762596253e+08Q,
	 3.048305382188066788914727298031581e+06Q,
	-3.562744619669248206793629956902640e+04Q,
	 2.625939669969919276517969162685481e+02Q,
	-1.127892258247142779807411800956670e+00Q,
	 2.538684230850237447901013114737673e-03Q,
	-2.583211865372192597501314906164821e-06Q,
	 9.590040155479072110068965236664741e-10Q,
	-9.347057516377501180809930331336725e-14Q,
	 1.386220962446620618073769872331415e-18Q,
	-1.138134595144440607620794028391791e-24Q,
	 5.491908932091769529804906638081618e-33Q,
	-1.332629445370080503706686280760692e-46Q,
      };
    };
#endif

int
main()
{
  std::cout << "\n\nSpouge Algorithm\n\n";
  std::cout << "\nspouge<float>\n";
  spouge<float>();
  std::cout << "\nspouge<double>\n";
  spouge<double>();
  std::cout << "\nspouge<long double>\n";
  spouge<long double>();
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\nspouge<__float128>\n";
  spouge<__float128>();
#endif
}
