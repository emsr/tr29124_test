/*
$HOME/bin_tr29124/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_erfc test_erfc.cpp -lquadmath
./test_erfc > test_erfc.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_erfc test_erfc.cpp -lquadmath
./test_erfc > test_erfc.txt
*/

#include <iostream>
#include <iomanip>
#include <cmath>
#include <LentzContinuedFraction.tcc>

  /**
   * Return the complementary error function by series:
   * @f[
   *    erfc(z) = e^{-z^2}\sum_{k=0}^{\infty} \frac{(-z)^k}{\Gamma(1+k/2)}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __erfc_series_1(_Tp __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      const auto _S_max_iter = 1000;
      const auto __rez = std::real(__z);
      const auto _S_eps = __gnu_cxx::__epsilon(__rez);
      const auto _S_g3d2 = __gnu_cxx::__const_root_pi(__rez) / _Real{2};
      auto __mzk = _Tp{1};
      _Real __irgam[2]{_Real{1}, _S_g3d2};
      auto __sum = _Tp{0};
      for (int __k = 0; __k < _S_max_iter; ++__k)
	{
	  const auto __term = __mzk * __irgam[__k & 1];
	  __sum += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__sum))
	    break;
	  __mzk *= -__z;
	  __irgam[__k & 1] /= __k;
	}
      return std::exp(-__z * __z) * __sum;
    }

  /**
   * Return the complementary error function by series:
   * @f[
   *    erfc(z) = 
   * @f]
   * @see On the Evaluation of Integrals Related to the Error Function
   * of Real and Complex Variable with High Relative Precision, Masatake Mori,
   * Publ. RIMS, Kyoto Univ. 19 (1983), 1081-1094
   */
  template<typename _Tp>
    _Tp
    __erfc_series(_Tp __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      const auto _S_max_iter = 1000;
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__z));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__z));

      const auto __h = _Real{1} / _Real{2};

      const auto __z2 = __z * __z;
      auto __sum = _Tp{0};
      auto __a = std::exp(__h * __h);
      auto __b = _Real{1};
      const auto __c = _Real{1} / __a / __a;
      for (auto __k = 1; __k < _S_max_iter; ++__k)
      {
	__a *= __c;
	__b *= __a;
	const auto __kh = __k * __h;
	const auto __kh2 = __kh * __kh;
	const auto __term = __b / (__kh2 + __z2);
	__sum += __term;
	if (std::abs(__term) < _S_eps * std::abs(__sum))
	  break;
      }
      auto __erfc = _Tp{2} * __z * __h * std::exp(-__z2)
		  * (_Tp{1} / __z2 / _Tp{2} + __sum) / _S_pi;

      auto __rez = std::real(__z);
      auto __imz = std::imag(__z);
      if (__rez < _Real{0})
	__erfc += _Real{2};
      if (std::abs(__rez) < _S_pi / __h)
	{
	  auto __sz = _Tp(__rez < _Tp{0} ? -1 : +1) * __z;
	  auto __Rcorr = _Tp{2} / std::expm1(_Tp{2} * _S_pi * __sz / __h);
	  if (__rez >= _Real{0})
	    {
	      if (std::abs(__imz) < _S_pi / _Tp{2} - __rez)
		__erfc -= __Rcorr;
	    }
	  else
	    {
	      if (std::abs(__imz) < _S_pi / _Tp{2} + __rez)
		__erfc += __Rcorr;
	    }
	}

      return __erfc;
    }

  /**
   * Use the continued fraction to evaluate the complementary
   * error function.
   */
  template<typename _Tp>
    _Tp
    __erfc_cont_frac(_Tp __z)
    {
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(std::real(__z));
      const auto __a = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 1 ? _Tp{1} : _Tp(__k - 1) / _Tp{2}; };
      using _AFun = decltype(__a);
      const auto __b = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 0 ? _Tp{0} : __k & 1 ? __z * __z : _Tp{1}; };
      using _BFun = decltype(__b);
      const auto __w = [__b](size_t __k, _Tp __z)
		 ->_Tp
		 { return __b(__k, __z)
		  / (std::sqrt(_Tp{1} + _Tp(2 * __k) / __z / __z) - _Tp{1})
		  / _Tp{2}; };
      using _WFun = decltype(__w);
      using _CFrac = _LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>;
      const _CFrac __cf(__a, __b, __w);
      return std::exp(-__z * __z) * __cf(__z) / _S_sqrt_pi;
    }

  /**
   * Use the continued fraction to evaluate the complementary
   * error function.
   * @f[
   *    w(z) = 
   * @f]
   */
  template<typename _Tp>
    _Tp
    __fadeeva_cont_frac(_Tp __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_i = _Cmplx{0, 1};
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(std::real(__z));
      const auto __a = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 1 ? _Tp{1} : _Tp(__k - 1) / _Tp{2}; };
      using _AFun = decltype(__a);
      const auto __b = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 0 ? _Tp{0} : __k & 1 ? __z * __z : _Tp{1}; };
      using _BFun = decltype(__b);
      const auto __w = [__b](size_t __k, _Tp __z)
		 ->_Tp
		 { return __b(__k, __z)
		  / (std::sqrt(_Tp{1} + _Tp(2 * __k) / __z / __z) - _Tp{1})
		  / _Tp{2}; };
      using _WFun = decltype(__w);
      using _CFrac = _LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>;
      const _CFrac __cf(__a, __b, __w);
      return _S_i * __cf(__z) / _S_sqrt_pi;
    }

  /**
   * Use the even continued fraction to evaluate the complementary
   * error function.
   * @f[
   *    erfc(z) = 
   * @f]
   */
  template<typename _Tp>
    _Tp
    __erfc_cont_frac_even(_Tp __z)
    {
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(std::real(__z));
      const auto __a = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 1
			 ? _Tp{2} * __z
			 : __k & 1 ? __z * __z : _Tp{1}; };
      using _AFun = decltype(__a);
      const auto __b = [](size_t __k, _Tp __z)
		 ->_Tp
		 { return __k == 0
			? _Tp{0}
			: __k & 1 ? _Tp{2} * __z * __z
				  : _Tp{1}; };
      using _BFun = decltype(__b);
      const auto __w = [__b](size_t, _Tp)->_Tp{ return _Tp{0}; };
      using _WFun = decltype(__w);
      using _CFrac = _LentzContinuedFraction<_Tp, _AFun, _BFun, _WFun>;
      const _CFrac __cf(__a, __b, __w);
      return std::exp(-__z * __z) * __cf(__z) / _S_sqrt_pi;
    }

  /**
   * Return the Fadeeva function.
   * @f[
   *    w(z) = e^{-z^2}erfc(-iz)
   * @f]
   */
  template<typename _Tp>
    _Tp
    __fadeeva(_Tp __z)
    {
      using _Real = std::__detail::__num_traits_t<_Tp>;
      using _Cmplx = std::complex<_Real>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::real(__z));
      const auto _S_i = _Cmplx{0, 1};
      if (__isnan(__z))
	return _Cmplx{_S_NaN, _S_NaN};
      else if (std::real(__z) < _Real{0})
	return _Real{2} * std::exp(-__z * __z) - __fadeeva(-__z);
      else if (std::abs(__z) < _Real{15})
	return 0;//FIXME
      else
	return __fadeeva_cont_frac(__z);
    }

template<typename _Tp>
  void
  test_erfc(_Tp proto = _Tp{})
  {
    using namespace std::literals::complex_literals;

    std::cout.precision(__gnu_cxx::__digits10(proto));
    auto width = std::cout.precision() + 8;
    std::cout << std::showpoint << std::scientific;

    const auto del = _Tp{1} / _Tp{100};
    for (int i = -200; i <= 1000; ++i)
      {
	auto x = del * i;
	auto __erfcs = __erfc_series(x);
	auto __erfccf = 0;//__erfc_cont_frac(x);
	std::cout << ' ' << x
		  << ' ' << std::setw(width) << std::erfc(x)
		  << ' ' << std::setw(width) << __erfcs
		  << ' ' << std::setw(width) << __erfccf
		  << '\n';
      }
  }


int
main()
{
  test_erfc(1.0);
}

