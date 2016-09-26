/*
$HOME/bin_tr29124/bin/g++ -I. -o test_gamma test_gamma.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./test_gamma > test_gamma.txt

$HOME/bin/bin/g++ -std=gnu++14 -DNO_LOGBQ -I. -o test_gamma test_gamma.cpp -lquadmath
./test_gamma > test_gamma.txt
*/

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include <bits/float128.h>

   //  Checbyshev coefficient matrix.
   constexpr int
   __cheby(unsigned int __n, unsigned int __k)
   {
     if (__k > __n)
       return 0;
     else if (__n == 1)
       return 1;
     else if (__n == 2)
       {
	 if (__k == 1)
	   return 0;
	 else if (__k == 2)
	   return 1;
       }
     else
       {
	 if (__k == 1)
	   return -__cheby(__n - 2, 1);
	 else if (__k == __n)
	   return 2 * __cheby(__n - 1, __k - 1);
	 else
	   return 2 * __cheby(__n - 1, __k - 1) - __cheby(__n - 2, __k);
       }
   }

  template<typename _Tp>
    _Tp
    __p(unsigned int __k, _Tp __g)
    {
      const auto _S_pi  = __gnu_cxx::__math_constants<_Tp>::__pi;
      auto __fact = std::sqrt(_Tp{2} / _S_pi);
      auto __sum = __cheby(2 * __k + 1, 1) * __fact
		 * std::exp(_Tp(__g + 0.5L))
		 / std::sqrt(_Tp(__g + 0.5L));
      for (int __a = 1; __a <= __k; ++__a)
	{
	  __fact *= _Tp(2 * __a - 1) / 2;
	  __sum += __cheby(2 * __k + 1, 2 * __a + 1) * __fact
		 * std::pow(_Tp(__a + __g + 0.5L), -_Tp(__a + 0.5L))
		 * std::exp(_Tp(__a + __g + 0.5L));
	}
      return __sum;
    }

  template<typename _Tp>
    void
    lanczos()
    {
      std::cout.precision(std::numeric_limits<_Tp>::digits10);
      std::cout << std::showpoint << std::scientific;

      // From Pugh..
      int __n_old = 0;
      int __n = -2 - 0.3 * std::log(std::numeric_limits<_Tp>::epsilon());
      std::cout << "n = " << __n << '\n';

      auto __g = __n - _Tp{0.5L};
      std::cout << "g = " << __g << '\n';
      while (__n != __n_old)
	{
	  std::cout << '\n';
	  std::vector<_Tp> __a;
	  for (unsigned int k = 1; k <= __n; ++k)
	    {
	      for (unsigned int j = 1; j <= k; ++j)
        	std::cout << "  C(" << std::setw(2) << k
			    << ", " << std::setw(2) << j
			    << ") = " << std::setw(4) << __cheby(k, j);
              std::cout << '\n';
	    }

	  std::cout << '\n';
	  auto __prev = std::numeric_limits<_Tp>::max();
	  for (unsigned int __k = 0; __k <= __n + 5; ++__k)
	    {
	      auto __curr = __p(__k, __g);
	      if (std::abs(__curr) > std::abs(__prev))
		{
		  __n_old = __n;
		  __n = __k;
		  __g = __n - _Tp{0.5L};
		  break;
		}
	      __prev = __curr;
	      std::cout << "  p(" << __k << ", " << __g << ") = " << __curr << '\n';
	    }
	  std::cout << "n = " << __n << '\n';
	}

      auto __log_gamma_lanczos =
	[=](_Tp __z)
	-> _Tp
	{
	  constexpr auto _S_ln_2 = __gnu_cxx::__math_constants<_Tp>::__ln_2;
	  constexpr auto _S_ln_pi = __gnu_cxx::__math_constants<_Tp>::__ln_pi;
	  constexpr auto _S_log_sqrt_2pi = (_S_ln_2 + _S_ln_pi) / _Tp{2};
	  auto __fact = _Tp{1};
	  auto __sum = _Tp{0.5L} * __p(0, __g);
	  for (unsigned int __k = 1; __k < __n; ++__k)
	    {
	      __fact *= (__z - __k + 1) / (__z + __k);
	      __sum += __fact * __p(__k, __g);
	    }
	  return _S_log_sqrt_2pi + std::log(__sum)
	       + (__z + 0.5L) * std::log(__z + __g + 0.5L)
	       - (__z + __g + 0.5L) - std::log(__z);
	};

      std::cout << '\n';
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = _Tp{0.01Q} * i;
	  std::cout << ' ' << z
		    << ' ' << __log_gamma_lanczos(z)
		    << ' ' << std::lgamma(z)
		    << ' ' << __log_gamma_lanczos(z) - std::lgamma(z) << '\n';
	}
    }

  template<typename _Tp>
    void
    spouge()
    {
      std::cout.precision(std::numeric_limits<_Tp>::digits10);
      std::cout << std::showpoint << std::scientific;

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
	  __c.push_back(__factc * std::pow(__ak, _Tp(__k + 0.5L)) * std::exp(__ak));
	  std::cout << "c_" << __k << " = " << __c.back() << '\n';
	}

      auto __log_gamma_spouge =
	[=](_Tp __z)
	-> _Tp
	{
	  // Reflection is right but auto and use of functions won't compile.
	  //if (__z <= -a)
	    //return std::log(_S_pi) - std::log(std::sin(_S_pi * __z)) - __log_gamma_spouge(_Tp{1} - __z);
	  //else
	    {
	      auto __sum = std::sqrt(_S_2pi);
	      for (int __k = 0; __k < __c.size(); ++__k)
		__sum += __c[__k] / (__z + __k + 1);
	      return std::log(__sum)
		   + (__z + _Tp{0.5L}) * std::log(__z + __a)
		   - (__z + __a) - std::log(__z);
	    }
	};

      for (int i = 0; i <= 500; ++i)
	{
	  auto z = _Tp{0.01L} * i;
	  std::cout << ' ' << z
		    << ' ' << __log_gamma_spouge(z)
		    << ' ' << std::lgamma(z)
		    << ' ' << __log_gamma_spouge(z) - std::lgamma(z) << '\n';
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

      for (int i = 0; i <= 500; ++i)
	{
	  auto z = _Tp{0.01L} * i;
	  _Tp x, y;
	  std::cout << ' ' << z
		    << ' ' << (y = __log_gamma_spouge(z))
		    << ' ' << std::lgamma(z)
		    << ' ' << __log_gamma_spouge(z) - std::lgamma(z)
		    << ' ' << (x = __log_gamma_inv(y))
		    << ' ' << x - y
		    << '\n';
	}
    }


  /**
   * A struct for Lanczos algorithm Chebyshev arrays of coefficients.
   */
  template<typename _Tp>
    class _GammaLanczos
    {
    };

  template<>
    class _GammaLanczos<float>
    {
    public:
    private:
      static constexpr float _S_g = 6.5F;
      static constexpr std::array<float, 7>
      _S_cheby
      {
	3.307139e+02F,
	-2.255998e+02F,
	6.989520e+01F,
	-9.058929e+00F,
	4.110107e-01F,
	-4.150391e-03F,
	-3.417969e-03F,
      };
    };

  template<>
    class _GammaLanczos<double>
    {
    public:
    private:
      static constexpr double _S_g = 9.5;
      static constexpr std::array<double, 10>
      _S_cheby
      {
	5.557569219204146e+03,
	-4.248114953727554e+03,
	1.881719608233706e+03,
	-4.705537221412237e+02,
	6.325224688788239e+01,
	-4.206901076213398e+00,
	1.202512485324405e-01,
	-1.141081476816908e-03,
	2.055079676210880e-06,
	1.280568540096283e-09,
      };
    };

  template<>
    class _GammaLanczos<long double>
    {
    public:
    private:
      static constexpr long double _S_g = 10.5L;
      static constexpr std::array<long double, 11>
      _S_cheby
      {
	1.440399692024250728e+04L,
	-1.128006201837065341e+04L,
	5.384108670160999829e+03L,
	-1.536234184127325861e+03L,
	2.528551924697309561e+02L,
	-2.265389090278717887e+01L,
	1.006663776178612579e+00L,
	-1.900805731354182626e-02L,
	1.150508317664389324e-04L,
	-1.208915136885480024e-07L,
	-1.518856151960790157e-10L,
      };
    };

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  template<>
    class _GammaLanczos<__float128>
    {
    public:
    private:
      static constexpr __float128 _S_g = 13.5Q;
      static constexpr std::array<__float128, 14>
      _S_cheby
      {
	 2.564476893267270739326759539239521e+05Q,
	-2.115503710351143292058877626137631e+05Q,
	 1.184019145386031178030546551084647e+05Q,
	-4.454725318050137764345651484051269e+04Q,
	 1.108444304734911694218772523380670e+04Q,
	-1.779039702565864879298077142685014e+03Q,
	 1.775917720477127714789670669062174e+02Q,
	-1.046153541164985987873520383376722e+01Q,
	 3.366596234159295518095091028749941e-01Q,
	-5.260241043537793269590387229607226e-03Q,
	 3.290526994052737279641248387319708e-05Q,
	-5.791840871423216610203896349432783e-08Q,
	 1.318192361117106013253744235732287e-11Q,
	-3.352799410216973507737605805778536e-17Q,
      };
    };
#endif


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
  std::cout << "\n\nLanczos Algorithm\n\n";
  std::cout << "\nlanczos<float>\n";
  lanczos<float>();
  std::cout << "\nlanczos<double>\n";
  lanczos<double>();
  std::cout << "\nlanczos<long double>\n";
  lanczos<long double>();
#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)
  std::cout << "\nlanczos<__float128>\n";
  lanczos<__float128>();
#endif

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
