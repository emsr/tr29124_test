// $HOME/bin_specfun/bin/g++ -o test_gamma test_gamma.cpp -lquadmath

// LD_LIBRARY_PATH=$HOME/bin_specfun/lib64:$LD_LIBRARY_PATH ./test_gamma > test_gamma.txt

#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <ext/cmath>
#include "float128.h"

//  Checbyshev coefficient matrix.
int
__c(unsigned int __n, unsigned int __k)
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
	return -__c(__n - 2, 1);
      else if (__k == __n)
	return 2 * __c(__n - 1, __k - 1);
      else
	return 2 * __c(__n - 1, __k - 1) - __c(__n - 2, __k);
    }
}

template<typename _Tp>
  _Tp
  __p(unsigned int __k, _Tp __g)
  {
    const auto _S_pi  = _Tp{3.1415926535897932384626433832795029L};
    auto __fact = std::sqrt(_Tp{2} / _S_pi);
    auto __sum = __c(2 * __k + 1, 1) * __fact
	       * std::exp(__g + 0.5)
	       / std::sqrt(__g + 0.5);
    for (int __a = 1; __a <= __k; ++__a)
      {
	__fact *= _Tp(2 * __a - 1) / 2;
	__sum += __c(2 * __k + 1, 2 * __a + 1) * __fact
	       * std::pow(__a + __g + 0.5L, -_Tp(__a + 0.5L))
	       * std::exp(__a + __g + 0.5L);
      }
    return __sum;
  }

template<typename _Tp>
  void
  lanczos()
  {
    std::cout.precision(std::numeric_limits<_Tp>::digits10);

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
			  << ") = " << std::setw(4) << __c(k, j);
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

    constexpr auto _S_log_sqrt_2pi = 9.189385332046727417803297364056176398620e-1L;

    auto
    __log_gamma_lanczos = [=](_Tp __z) -> _Tp
    {
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

    const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
    const auto _S_pi  = _Tp{3.1415926535897932384626433832795029L};
    const auto _S_2pi = _Tp{6.283185307179586476925286766559005768391L};
    auto a = _Tp{1};
    const auto __fact = _Tp{1} / std::sqrt(_S_2pi);
    while (_S_eps <= __fact * std::pow(_S_2pi, -a) / std::sqrt(a))
      {
        std::cout << "err = " << __fact * std::pow(_S_2pi, -a) / std::sqrt(a) << '\n';
        a += _Tp{1};
      }
    std::cout << "a = " << a << '\n';

    std::vector<_Tp> c;
    auto __factc = _Tp{1};
    c.push_back(__factc * std::sqrt(a - 1) * std::exp(a - 1));
    std::cout << "c_0 = " << c.back() << '\n';
    auto __sum = __factc * std::exp(a - 1);
    for (int __k = 1; __k < std::ceil(a); ++__k)
      {
	__factc *= -_Tp{1} / _Tp(__k);
	auto __ak = _Tp(a - __k - 1);
	c.push_back(__factc * std::pow(__ak, _Tp(__k + 0.5L)) * std::exp(__ak));
	std::cout << "c_" << __k << " = " << c.back() << '\n';
      }

    auto
    __log_gamma_spouge = [=](_Tp __z) -> _Tp
    {
      // Reflection is right but auto and use of functions won't compile.
      //if (__z <= -a)
	//return std::log(_S_pi) - std::log(std::sin(_S_pi * __z)) - __log_gamma_spouge(_Tp{1} - __z);
      //else
	{
	  auto __sum = std::sqrt(_S_2pi);
	  for (int __k = 0; __k < c.size(); ++__k)
	    __sum += c[__k] / (__z + __k + 1);
	  return std::log(__sum)
	       + (__z + _Tp{0.5L}) * std::log(__z + a)
	       - (__z + a) - std::log(__z);
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
    const auto _S_log_10 = std::log(10.0);
    const auto _S_log_sqrt2pi = 0.5 * std::log(_S_2pi);
    auto inv_log_gamma
    {
      [=](_Tp y)
      -> _Tp
      {
	auto x = y * _S_log_10 - _S_log_sqrt2pi;
	auto x0 = x;
	for (int i = 0; i < 100; ++i)
	  x = (x0 + x) / std::log(x);
	return x;
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
		  << ' ' << (x = inv_log_gamma(y))
		  << ' ' << x - y
		  << '\n';
      }
  }

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
