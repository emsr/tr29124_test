/*
LD_LIBRARY_PATH=..:$LD_LIBRARY_PATH $HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I.. -o build_clenshaw_curtis build_clenshaw_curtis.cpp -lquadmath -L.. -lwburkhardt
LD_LIBRARY_PATH=..:$LD_LIBRARY_PATH ./build_clenshaw_curtis > build_clenshaw_curtis.txt
*/

#include <vector>
#include <ext/cmath>
#include <bits/specfun_state.h>
#include <iostream>
#include <iomanip>
#include "wrap_burkhardt.h"

  /**
   * Return a vector of zeros of the Chebyshev function of the second kind
   * of order @f$ n @f$, @f$ U_n(x) @f$.
   * The zeros are given by:
   * @f[
   *   x_k = \cos\left(\frac{k\pi}{n+1}\right), k \elem {1, ..., n}
   * @f]
   */
  template<typename _Tp>
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
    __chebyshev_u_zeros(unsigned int __n)
    {
      const auto _S_pi = __gnu_cxx::__const_pi<_Tp>();
      std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __pt(__n);
      for (unsigned int __k = 1; __k <= __n; ++__k)
	{
	  auto __arg = _Tp(__k) / _Tp(__n + 1);
	  auto __half = __gnu_cxx::__fp_is_equal<_Tp>(__arg, _Tp{0.5L});
	  auto __z = (__half ? _Tp{0} : __gnu_cxx::cos_pi(__arg));
	  auto __w = _S_pi * (_Tp{1} - __z * __z) / _Tp(__n + 1);
	  __pt[__k - 1].__zero = __z;
	  __pt[__k - 1].__weight = __w;
	}
      return __pt;
    }

/**
 * 
 *
 * @f[
 *    w_k = frac{c_k}{n}\left[1-\sum_{j=1}^{[n/2]}\frac{b_k}{4j^2-1}
 *            \cos\left(\frac{2jk\pi}{n}\right)\right]
 *    \mbox{   } k = 0, 1, ..., n
 * @f]
 * 
 * @see Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules
 */
template<typename _Tp, std::size_t _Num>
  std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
  build_clenshaw_curtis_sum()
  {
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __out(_Num + 1);
    if (_Num == 1)
      {
	__out[0].__zero = _Tp{0};
	__out[0].__weight = _Tp{2};
	return __out;
      }
    else
      {
	const auto _S_pi = __gnu_cxx::__const_pi<_Tp>();
	auto uz = __chebyshev_u_zeros<_Tp>(_Num - 1);
	__out[0].__zero = _Tp{+1};
	__out[0].__weight = _Tp{1} / (_Num * _Num - 1 + _Num % 2);
	for (auto __k = 1u; __k <= uz.size(); ++__k)
	  {
	    __out[__k].__zero = uz[__k - 1].__zero;

	    auto __sum = _Tp{0};
	    for (auto __j = 1u; __j <= _Num / 2; ++__j)
	      {
		auto __b = _Tp(__j == _Num / 2 ? 1 : 2);
		__sum += __b * std::cos(2 * __j * __k * _S_pi / _Num)
		       / _Tp(4 * __j * __j - 1);
	      }
	    auto __w = _Tp{2} * (_Tp{1} - __sum) / _Tp{_Num};
	    __out[__k].__weight = __w;
	  }
	__out[_Num].__zero = _Tp{-1};
	__out[_Num].__weight = __out[0].__weight;
	return __out;
      }
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{c_k}{n}\left[1-\sum_{j=1}^{[n/2]}\frac{b_k}{4j^2-1}
 *            \cos\left(\frac{2jk\pi}{n}\right)\right]
 *    \mbox{   } k = 0, 1, ..., n
 * @f]
 * 
 * @see Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules
 */
template<typename _Tp, std::size_t _Num>
  std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
  build_clenshaw_curtis_fft()
  {
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __out(_Num + 1);
    return __out;
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{2}{n}\left[1-2\sum_{j=1}^{[n/2]}\frac{1}{4j^2-1}
 *            \cos\left(j\frac{(2k+1)\pi}{n}\right)\right]
 *    \mbox{   } k = 0, 1, ..., n-1
 * @f]
 */
template<typename _Tp, std::size_t _Num>
  std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
  build_fejer_1_sum()
  {
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __out(_Num + 1);
    if (_Num == 1)
      {
	__out[0].__zero = _Tp{0};
	__out[0].__weight = _Tp{2};
	return __out;
      }
    else
      {
	const auto _S_pi = __gnu_cxx::__const_pi<_Tp>();
	auto uz = __chebyshev_u_zeros<_Tp>(_Num - 1);
	__out[0].__zero = _Tp{+1};
	__out[0].__weight = _Tp{1} / (_Num * _Num - 1 + _Num % 2);
	for (auto __k = 1u; __k <= uz.size(); ++__k)
	  {
	    __out[__k].__zero = uz[__k - 1].__zero;

	    auto __sum = _Tp{0};
	    for (auto __j = 1u; __j <= _Num / 2; ++__j)
	      __sum += _Tp{2} * std::cos(__j * (2 * __k + 1) * _S_pi / _Num)
		     / _Tp(4 * __j * __j - 1);
	    auto __w = _Tp{2} * (_Tp{1} - __sum) / _Tp{_Num};
	    __out[__k].__weight = __w;
	  }
	__out[_Num].__zero = _Tp{-1};
	__out[_Num].__weight = __out[0].__weight;
	return __out;
      }
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{2}{n}\left[1-2\sum_{j=1}^{[n/2]}\frac{1}{4j^2-1}
 *            \cos\left(j\frac{(2k+1)\pi}{n}\right)\right]
 *    \mbox{   } k = 0, 1, ..., n-1
 * @f]
 */
template<typename _Tp, std::size_t _Num>
  std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
  build_fejer_1_fft()
  {
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __out(_Num);
    return __out;
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{4}{n}\sin\left(\frac{k\pi}{n}\right)
 *        \sum_{j=1}^{[n/2]}\frac{1}{2j-1}\sin left((2j-1)\frac{k\pi}{n}\right)
 *    \mbox{   } k = 0, 1, ..., n
 * @f]
 */
template<typename _Tp, std::size_t _Num>
  std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
  build_fejer_2_sum()
  {
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __out(_Num + 1);
    if (_Num == 0)
      {
	__out[0].__zero = _Tp{0};
	__out[0].__weight = _Tp{2};
	return __out;
      }
    else if (_Num == 1)
      {
	__out[0].__zero = _Tp{-0.5L};
	__out[0].__weight = _Tp{1};
	__out[1].__zero = _Tp{+0.5L};
	__out[1].__weight = _Tp{1};
	return __out;
      }
    else
      {
	const auto _S_pi = __gnu_cxx::__const_pi<_Tp>();
	auto uz = __chebyshev_u_zeros<_Tp>(_Num - 1);
	__out[0].__zero = _Tp{+1};
	__out[0].__weight = _Tp{1} / (_Num * _Num - 1 + _Num % 2);
	for (auto __k = 1u; __k <= uz.size(); ++__k)
	  {
	    __out[__k].__zero = uz[__k - 1].__zero;

	    auto __sum = _Tp{0};
	    for (auto __j = 1u; __j <= _Num / 2; ++__j)
	      __sum += std::sin((2 * __j - 1) * __k * _S_pi / _Num)
		     / _Tp(2 * __j - 1);
	    auto __w = _Tp{4} * std::sin(__k * _S_pi / _Num) * __sum
		     / _Tp{_Num};
	    __out[__k].__weight = __w;
	  }
	__out[_Num].__zero = _Tp{-1};
	__out[_Num].__weight = __out[0].__weight;
	return __out;
      }
  }

/**
 * 
 *
 * @f[
 *    w_k = frac{4}{n}\sin\left(j\frac{k\pi}{n}\right)
 *        \sum_{j=1}^{[n/2]}\frac{1}{2j-1}\sin left((2j-1)\frac{k\pi}{n}\right)
 *    \mbox{   } k = 0, 1, ..., n
 * @f]
 */
template<typename _Tp, std::size_t _Num>
  std::vector<__gnu_cxx::__quadrature_point_t<_Tp>>
  build_fejer_2_fft()
  {
    std::vector<__gnu_cxx::__quadrature_point_t<_Tp>> __out(_Num + 1);
    return __out;
  }

int
main()
{
  std::cout.precision(__gnu_cxx::__digits10<long double>());
  auto w = 8 + std::cout.precision();

  auto cc24b = burkhardt::clenshaw_curtis_rule(25);

  auto cc24 = build_clenshaw_curtis_sum<long double, 24>();
  std::cout << "\nClenshaw-Curtis 24\n";
  int i = 0;
  for (const auto& cc : cc24)
    {
      std::cout << std::setw(w) << cc.__zero << ' '
		<< std::setw(w) << cc.__weight << ' '
		<< std::setw(w) << cc.__zero - cc24b[i].__zero << ' '
		<< std::setw(w) << cc.__weight - cc24b[i].__weight
		<< '\n';
      ++i;
    }

  auto cc48b = burkhardt::clenshaw_curtis_rule(49);

  std::cout << "\nClenshaw-Curtis 48\n";
  auto cc48 = build_clenshaw_curtis_sum<long double, 48>();
  i = 0;
  for (const auto& cc : cc48)
    {
      std::cout << std::setw(w) << cc.__zero << ' '
		<< std::setw(w) << cc.__weight << ' '
		<< std::setw(w) << cc.__zero - cc48b[i].__zero << ' '
		<< std::setw(w) << cc.__weight - cc48b[i].__weight
		<< '\n';
      ++i;
    }

  auto f1_24b = burkhardt::fejer_1_rule(25);

  auto f1_24 = build_fejer_1_sum<long double, 24>();
  std::cout << "\nFejer1 24\n";
  i = 0;
  for (const auto& f1 : f1_24)
    {
      std::cout << std::setw(w) << f1.__zero << ' '
		<< std::setw(w) << f1.__weight << ' '
		<< std::setw(w) << f1.__zero - f1_24b[i].__zero << ' '
		<< std::setw(w) << f1.__weight - f1_24b[i].__weight
		<< '\n';
      ++i;
    }

  auto f2_24b = burkhardt::fejer_2_rule(25);

  auto f2_24 = build_fejer_2_sum<long double, 24>();
  std::cout << "\nFejer2 24\n";
  i = 0;
  for (const auto& f2 : f2_24)
    {
      std::cout << std::setw(w) << f2.__zero << ' '
		<< std::setw(w) << f2.__weight << ' '
		<< std::setw(w) << f2.__zero - f2_24b[i].__zero << ' '
		<< std::setw(w) << f2.__weight - f2_24b[i].__weight
		<< '\n';
      ++i;
    }
}
