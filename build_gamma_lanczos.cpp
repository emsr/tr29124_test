/*
$HOME/bin_tr29124/bin/g++ -I. -o build_gamma_lanczos build_gamma_lanczos.cpp -lquadmath
LD_LIBRARY_PATH=$HOME/bin_tr29124/lib64:$LD_LIBRARY_PATH ./build_gamma_lanczos > build_gamma_lanczos.txt

$HOME/bin/bin/g++ -std=gnu++14 -DNO_LOGBQ -I. -o build_gamma_lanczos build_gamma_lanczos.cpp -lquadmath
./build_gamma_lanczos > build_gamma_lanczos.txt
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
		 * std::exp(_Tp(__g + 0.5Q))
		 / std::sqrt(_Tp(__g + 0.5Q));
      for (int __a = 1; __a <= __k; ++__a)
	{
	  __fact *= _Tp(2 * __a - 1) / 2;
	  __sum += __cheby(2 * __k + 1, 2 * __a + 1) * __fact
		 * std::pow(_Tp(__a + __g + 0.5Q), -_Tp(__a + 0.5Q))
		 * std::exp(_Tp(__a + __g + 0.5Q));
	}
      return __sum;
    }

  template<typename _Tp>
    void
    lanczos()
    {
      std::cout.precision(std::numeric_limits<_Tp>::digits10);
      std::cout << std::showpoint << std::scientific;
      auto width = 8 + std::cout.precision();

      // From Pugh..
      int __n_old = 0;
      int __n = -2 - 0.3 * std::log(std::numeric_limits<_Tp>::epsilon());
      std::cout << "n_Pugh = " << __n << '\n';

      auto __g = __n - _Tp{0.5Q};
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
		  __g = __n - _Tp{0.5Q};
		  break;
		}
	      __prev = __curr;
	      std::cout << "  p(" << __k << ", " << __g << ") = " << __curr << '\n';
	    }
	  std::cout << "n = " << __n << '\n';
	}

      auto __log_gamma1p_lanczos =
	[=](_Tp __z)
	-> _Tp
	{
	  constexpr auto _S_ln_2 = __gnu_cxx::__math_constants<_Tp>::__ln_2;
	  constexpr auto _S_ln_pi = __gnu_cxx::__math_constants<_Tp>::__ln_pi;
	  constexpr auto _S_log_sqrt_2pi = (_S_ln_2 + _S_ln_pi) / _Tp{2};
	  auto __fact = _Tp{1};
	  auto __sum = _Tp{0.5Q} * __p(0, __g);
	  for (unsigned int __k = 1; __k < __n; ++__k)
	    {
	      __fact *= (__z - __k + 1) / (__z + __k);
	      __sum += __fact * __p(__k, __g);
	    }
	  return _S_log_sqrt_2pi + std::log(__sum)
	       + (__z + _Tp{0.5Q}) * std::log(__z + __g + _Tp{0.5Q})
	       - (__z + __g + _Tp{0.5Q});
	};

      std::cout << '\n'
		<< ' ' << std::setw(width) << "z"
		<< ' ' << std::setw(width) << "lanczos"
		<< ' ' << std::setw(width) << "lgamma"
		<< ' ' << std::setw(width) << "delta"
		<< '\n';
      for (int i = 0; i <= 500; ++i)
	{
	  auto z = _Tp{0.01Q} * i;
	  std::cout << ' ' << std::setw(width) << z
		    << ' ' << std::setw(width) << __log_gamma1p_lanczos(z - _Tp{1})
		    << ' ' << std::setw(width) << std::lgamma(z)
		    << ' ' << std::setw(width) << __log_gamma1p_lanczos(z - _Tp{1}) - std::lgamma(z)
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
}
