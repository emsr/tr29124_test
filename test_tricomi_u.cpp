/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_tricomi_u test_tricomi_u.cpp -lquadmath -L. -lwgsl -lburkhardt
./test_tricomi_u > test_tricomi_u.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -I. -o test_tricomi_u test_tricomi_u.cpp -lquadmath -L. -lwgsl -lburkhardt
./test_tricomi_u > test_tricomi_u.txt
*/

#include <ext/cmath>
#include <limits>
#include <iostream>
#include <iomanip>

//namespace std
//{
//namespace __detail
//{

  /**
   * @brief  Return the Tricomi confluent hypergeometric function
   * @f[
   *   U(a,c,x) = \frac{\Gamma(1-c)}{\Gamma(a-c+1)} {}_1F_1(a;c;x)
   *       + \frac{\Gamma(c-1)}{\Gamma(a)} x^{1-c} {}_1F_1(a-c+1;2-c;x)
   * @f]
   * @param  __a  The @a numerator parameter.
   * @param  __c  The @a denominator parameter.
   * @param  __x  The argument of the confluent hypergeometric function.
   * @return  The Tricomi confluent hypergeometric function.
   */
  template<typename _Tp>
    _Tp
    __tricomi_u_naive(_Tp __a, _Tp __c, _Tp __x)
    {
      auto __U1 = _Tp{};
      auto __b = __a - __c + _Tp{1};
      auto __ib = __gnu_cxx::__fp_is_integer(__b);
      if (!__ib || (__ib && __ib() > 0))
	__U1 = std::__detail::__gamma(_Tp{1} - __c)
	     * std::__detail::__conf_hyperg(__a, __c, __x)
	     / std::__detail::__gamma(__b);

      auto __U2 = _Tp{};
      auto __ia = __gnu_cxx::__fp_is_integer(__a);
      if (!__ia || (__ia && __ia() > 0))
	__U2 = std::__detail::__gamma(__c - _Tp{1})
	     * std::pow(__x, _Tp{1} - __c)
	     * std::__detail::__conf_hyperg(__b, _Tp{2} - __c, __x)
	     / std::__detail::__gamma(__a);

      return __U1 + __U2;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __tricomi_u_asymp(_Tp __a, _Tp __c, _Tp __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__z);
      const unsigned int _S_max_iter = 100000u;
      auto __b = __a - __c + _Tp{1};
      auto __term = _Tp{1};
      auto _Usum = _Tp{1};
      for (auto __k = 1u; __k < _S_max_iter; ++__k)
	{
	  __term *= -(__a + _Tp(__k - 1)) * (__b + _Tp(__k - 1))
		  / _Tp(__k) / __z;
	  _Usum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(_Usum))
	    break;
	}
      return std::pow(__z, -__a) * _Usum;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __tricomi_u_c_pos_int(_Tp __a, int __m, _Tp __z)
    {
      const auto _S_eps = __gnu_cxx::__epsilon(__z);
      const unsigned int _S_max_iter = 100000u;
//std::cout << '\n';

      auto __term1 = _Tp{1};
      auto _U1 = __term1;
      for (auto __k = 1; __k <= __m - 2; ++__k)
	{
	  __term1 *= (__a + _Tp(-__m + __k)) * __z
		   / _Tp(1 - __m + __k) / _Tp(__k);
	  _U1 += __term1;
//std::cout << "_U1 = " << _U1 << '\n';
	}
      _U1 *= std::__detail::__factorial<_Tp>(__m - 2)
	   / std::__detail::__gamma(__a)
	   / std::pow(__z, _Tp(__m - 1));

      const auto __b = __a + _Tp(1 - __m);
      auto __psi2 = std::__detail::__psi(__a)
		  - std::__detail::__psi<_Tp>(1)
		  - std::__detail::__psi<_Tp>(__m)
		  + std::log(__z);
      auto __fact2 = _Tp{1};
      auto _U2 = __psi2;
      for (auto __k = 1u; __k < _S_max_iter; ++__k)
	{
	  __psi2 += _Tp{1} / (__a + _Tp(__k - 1))
		  - _Tp{1} / _Tp(__k)
		  - _Tp{1} / _Tp(__m + __k - 1);
	  __fact2 *= (__a + _Tp(__k)) * __z / _Tp(__m + __k) / _Tp(__k);
	  auto __term2 = __psi2 * __fact2;
	  _U2 += __term2;
//std::cout << "_U2 = " << _U2 << '\n';
	  if (std::abs(__term2) < _S_eps * std::abs(_U2))
	    break;
	}
      _U2 *= (__m & 1 ? -1 : +1)
	   / std::__detail::__factorial<_Tp>(__m - 1)
	   / std::__detail::__gamma(__b);

      return _U1 + _U2;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __tricomi_u_c_nonpos_int(_Tp __a, _Tp __m, _Tp __z)
    {
      auto __b = __a + _Tp(1 - __m);
      return std::pow(__z, _Tp(1 - __m))
	   * __tricomi_u_c_pos_int(__b, 2 - __m, __z);
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __tricomi_u_ac_int(_Tp __n, _Tp __m, _Tp __z)
    {
      const unsigned int _S_max_iter = 100000u;
      const auto _S_eps = __gnu_cxx::__epsilon(__z);
      auto __term = _Tp{1};
      auto _Usum = _Tp{1};
      for (auto __k = 1u; __k < -__n; ++__k)
	{
	  __term *= _Tp(__n + __k - 1) / _Tp(__m + __k - 1) / _Tp(__k) * __z;
	  _Usum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(_Usum))
	    break;
	}
      return (__n & 1 ? -1 : +1) * __gnu_cxx::pochhammer(__m, -__n) * _Usum;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __tricomi_u_series(_Tp __a, _Tp __c, _Tp __z)
    {
      return 0;
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __tricomi_u(_Tp __a, _Tp __c, _Tp __z)
    {
      auto __aint = __gnu_cxx::__fp_is_integer(__a);
      auto __cint = __gnu_cxx::__fp_is_integer(__c);
      if (__cint)
	{
	  if (__cint() > 0)
	    return __tricomi_u_c_pos_int(__a, __cint(), __z);
	  else
	    return __tricomi_u_c_nonpos_int(__a, __cint(), __z);
	}
      else
	return __tricomi_u_series(__a, __c, __z);
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __whittaker_m(_Tp __kappa, _Tp __mu, _Tp __z)
    {
      return std::exp(-__z / 2) * std::pow(__z, 0.5 + __mu)
	   * __gnu_cxx::conf_hyperg(0.5 + __mu - __kappa, 1 + 2 * __mu, __z);
    }

  /**
   * 
   */
  template<typename _Tp>
    _Tp
    __whittaker_w(_Tp __kappa, _Tp __mu, _Tp __z)
    {
      return std::exp(-__z / 2) * std::pow(__z, 0.5 + __mu)
	   * __tricomi_u(0.5 + __mu - __kappa, 1 + 2 * __mu, __z);
    }

//} // namespace __detail
//} // namespace std

template<typename _Tp>
  void
  test_tricomi_u(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    auto a = _Tp{1.2Q};
    auto c = _Tp{0.2Q};
    for (int i = 0; i < +200; ++i)
    {
      auto z = _Tp{0.1Q} * i;
      std::cout << ' ' << std::setw(6) << z
		<< ' ' << std::setw(width) << __tricomi_u_naive(a, c, z)
		<< '\n';
    }

    std::cout << "\nInteger c = m\n";
    auto z = _Tp{0.5};
    std::cout << " a = " << std::setw(6) << a << '\n';
    std::cout << " z = " << std::setw(6) << z << '\n';
    for (auto m = 1u; m <= +20; ++m)
    {
      std::cout << ' ' << std::setw(6) << m
		<< ' ' << std::setw(width) << __tricomi_u_naive(a, _Tp(m), z)
		<< ' ' << std::setw(width) << __tricomi_u_c_pos_int(a, m, z)
		<< '\n';
    }
  }

int
main()
{
  test_tricomi_u(1.0);
}
