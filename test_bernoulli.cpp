/*
$HOME/bin_specfun/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_bernoulli test_bernoulli.cpp -lquadmath
./test_bernoulli > test_bernoulli.txt

$HOME/bin/bin/g++ -std=gnu++17 -g -Wall -Wextra -Wno-psabi -I. -o test_bernoulli test_bernoulli.cpp -lquadmath
./test_bernoulli > test_bernoulli.txt
*/ 

#include <cmath>
//#include <complex>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <vector>
// Why no math stuff without this?
#include <ext/math_const.h>
#include <bits/sf_gamma.tcc>


  template<typename _RealTp, typename _IntTp,
	   _IntTp _Num = 1, _IntTp _Den = 1>
    constexpr _RealTp
    __frac = _RealTp(_Num) / _RealTp(_Den);

  template<typename _RealTp,
	   unsigned long long _Num = 1,
           unsigned long long _Den = 1>
    constexpr _RealTp
    frac = __frac<_RealTp, unsigned long long, _Num, _Den>;

  /**
   * Return the Bernoulli number from lookup or by series expansion.
   */
  template<typename _Tp>
    _Tp
    __bernoulli_series(unsigned int __n)
    {
      static constexpr std::size_t _S_len = 24;
      static constexpr _Tp
      _S_num[_S_len]
      {
	 frac<_Tp>,
	-frac<_Tp, 1ull, 2ull>,
	 frac<_Tp, 1ull, 6ull>,            _Tp(0ull),
	-frac<_Tp, 1ull, 30ull>,           _Tp(0ull),
	 frac<_Tp, 1ull, 42ull>,           _Tp(0ull),
	-frac<_Tp, 1ull, 30ull>,           _Tp(0ull),
	 frac<_Tp, 5ull, 66ull>,           _Tp(0ull),
	-frac<_Tp, 691ull, 2730ull>,       _Tp(0ull),
	 frac<_Tp, 7ull, 6ull>,            _Tp(0ull),
	-frac<_Tp, 3617ull, 510ull>,       _Tp(0ull),
	 frac<_Tp, 43867ull, 798ull>,      _Tp(0ull),
	-frac<_Tp, 174611ull, 330ull>,     _Tp(0ull),
	 frac<_Tp, 854513ull, 138ull>,     _Tp(0ull)
      };

      if (__n == 0)
	return frac<_Tp>;
      else if (__n == 1)
	return -frac<_Tp, 1, 2>;
      else if (__n % 2 == 1) // Take care of the rest of the odd ones.
	return _Tp(0);
      else if (__n < _S_len) // return small evens that are painful for the series.
	return _S_num[__n];
      else
	{
	  _Tp __fact = _Tp(1);
	  if ((__n / 2) % 2 == 0)
	    __fact *= -1;
	  for (auto __k = 1u; __k <= __n; ++__k)
	    __fact *= __k / __gnu_cxx::__const_pi<_Tp>();
	  __fact *= _Tp(2);

	  _Tp __sum = _Tp(0);
	  for (auto __i = 1u; __i < 1000; ++__i)
	    {
	      _Tp __term = std::pow(_Tp(__i), -_Tp(__n));
	      if (__term < __gnu_cxx::__epsilon<_Tp>())
        	break;
	      __sum += __term;
	    }

	  return __fact * __sum;
	}
    }

  /**
   * Return the Bernoulli polynomial @f$ B_n(x) @f$ of order n at argument x.
   *
   * The values at 0 and 1 are equal to the corresponding Bernoulli number:
   * @f[
   *   B_n(0) = B_n(1) = B_n
   * @f]
   *
   * The derivative is proportional to the previous polynomial:
   * @f[
   *   B_n'(x) = n * B_{n-1}(x)
   * @f]
   *
   * The series expansion is:
   * @f[
   *   B_n(x) = \sum_{k=0}^{n} B_k binom{n}{k} x^{n-k}
   * @f]
   *
   * A useful argument promotion is:
   * @f[
   *   B_n(x+1) - B_n(x) = n * x^{n-1}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __bernoulli(unsigned int __n, _Tp __x)
    {
      auto _B_n = std::__detail::__bernoulli<_Tp>(0);
      auto __bincoef = _Tp{1};
      for (auto __k = 1u; __k <= __n; ++__k)
      {
	__bincoef *= _Tp(__n + 1 - __k) / _Tp(__k);
	_B_n = __x * _B_n + __bincoef * std::__detail::__bernoulli<_Tp>(__k);
      }

      return _B_n;
    }

  /**
   * Return the Euler number from lookup or by series expansion.
   *
   * The Euler numbers are given by the recursive sum:
   * @f[
   *   E_n = B_n(1) = B_n
   * @f]
   * where @f$ E_0 = 1 @f$, @f$ E_1 = 0 @f$, @f$ E_2 = -1 @f$
   *
   * @todo Find a way to predict the maximum Euler number for a type.
   */
  template<typename _Tp>
    _Tp
    __euler_series(unsigned int __n)
    {
      static constexpr std::size_t _S_len = 22;
      static constexpr _Tp
      _S_num[_S_len]
      {
	 1ll, 0,
	-1ll, 0ll,
	 5ll, 0ll,
	-61ll, 0ll,
	 1385ll, 0ll,
	-50521ll, 0ll,
	 2702765ll, 0ll,
	-199360981ll, 0ll,
	 19391512145ll, 0ll,
	-2404879675441ll, 0ll,
	 370371188237525ll, 0ll,
	//-69348874393137901ll, 0ll,
      };

      if (__n == 0)
	return _Tp{1};
      else if (__n & 1)
	return _Tp{0};
      else if (__n == 2)
        return _Tp{-1};
      else if (__n < _S_len)
	return _S_num[__n];
      else
	{
	  std::vector<_Tp> _En(__n + 1);
	  _En[0] = _Tp{1};
	  _En[1] = _Tp{0};
	  _En[2] = _Tp{-1};

	  for (auto __i = 3u; __i <= __n; ++__i)
	    {
	      _En[__i] = 0;

	      if (__i % 2 == 0)
		{
		  for (auto __j = 2u; __j <= __i; __j += 2u)
		    _En[__i] -= std::__detail::__bincoef<_Tp>(__i, __j)
			      * _En[__i - __j];
		}
	    }
	  return _En[__n];
	}
    }

  /**
   * @brief This returns Euler number @f$ E_n @f$.
   *
   * @param __n the order n of the Euler number.
   * @return  The Euler number of order n.
   */
  template<typename _Tp>
    _Tp //_GLIBCXX14_CONSTEXPR
    __euler(unsigned int __n)
    { return __euler_series<_Tp>(__n); }

  /**
   * Return the Euler polynomial @f$ E_n(x) @f$ of order n at argument x.
   *
   * The derivative is proportional to the previous polynomial:
   * @f[
   *   E_n'(x) = n E_{n-1}(x)
   * @f]
   *
   * @f[
   *   E_n(1/2) = \frac{E_n}{2^n}, \mbox{ where } E_n
   *             \mbox{ is the n-th Euler number.}
   * @f]
   *
   */
  template<typename _Tp>
    _Tp
    __euler(unsigned int __n, _Tp __x)
    {
      auto __bx1 = __bernoulli(__n + 1, __x );
      auto __bx2 = __bernoulli(__n + 1, _Tp{0.5L} * __x );

      auto _E_n = _Tp{2} * (__bx1 - __bx2 * std::pow(_Tp{2}, _Tp(__n + 1)))
		/ _Tp(__n + 1);

      return _E_n;
    }

  /**
   * Return the Stirling number of the second kind by series expansion.
   * The series is:
   * @f[
   *   \sigma_n^{(m)} = \sum_{k=0}^{m}\frac{(-1)^{m-k}k^n}{(m-k)!k!}
   * @f]
   *
   * @todo Find a way to predict the maximum Stirling number for a type.
   */
  template<typename _Tp>
    _Tp
    __stirling_2_series(unsigned int __n, unsigned int __m)
    {
      if (__m > std::__detail::_S_num_factorials<_Tp>)
	{
	  auto _S2 = _Tp{0};
	  for (auto __k = 0u; __k <= __m; ++__k)
	    {
	      auto __lf1 = std::__detail::__log_factorial<_Tp>(__k);
	      auto __lf2 = std::__detail::__log_factorial<_Tp>(__m - __k);
	      _S2 += ((__m - __k) & 1 ? _Tp{-1} : _Tp{1})
		   * std::exp(__n * std::log(__k) - __lf1 - __lf2);
	    }
	  return _S2;
	}
      else
	{
	  auto _S2 = _Tp{0};
	  for (auto __k = 0u; __k <= __m; ++__k)
	    {
	      _S2 += ((__m - __k) & 1 ? _Tp{-1} : _Tp{1})
		   * std::pow(__k, __n)
		   / std::__detail::__factorial<_Tp>(__k)
		   / std::__detail::__factorial<_Tp>(__m - __k);
	    }
	  // @todo Only round if the sum is less than
	  // the maximum representable integer.
	  // Find or make a tool for this.
	  return std::nearbyint(_S2);
	}
    }

  /**
   * Return the Stirling number of the second kind by recursion.
   * The recursion is
   * @f[
   *   \sigma_n^{(m)} = m \sigma_{n-1}^{(m)} + \sigma_{n-1}^{(m-1)}
   * @f]
   * with starting values
   * @f[
   *   \sigma_0^{(0\rightarrow m)} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   \sigma_{0\rightarrow n}^{(0)} = {1, 0, 0, ..., 0}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __stirling_2_recur(unsigned int __n, unsigned int __m)
    {
      if (__n == 0)
	return _Tp(__m == 0);
      else if (__m == 0)
	return _Tp(__n == 0);
      else
	{
	  std::vector<_Tp> __sigold(__m + 1), __signew(__m + 1);
	  __sigold[1] = _Tp{1};
	  if (__n == 1)
	    return __sigold[__m];
	  for (auto __in = 1u; __in <= __n; ++__in)
	    {
	      __signew[1] = __sigold[1];
	      for (auto __im = 2u; __im <= __m; ++__im)
		__signew[__im] = __im * __sigold[__im] + __sigold[__im - 1];
	      std::swap(__sigold, __signew);
	    }
	  return __signew[__m];
	}
    }

  /**
   * Return the Stirling number of the second kind from lookup
   * or by series expansion.
   *
   * @todo Look into asymptotic solutions.
   */
  template<typename _Tp>
    _Tp
    __stirling_2(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return _Tp{0};
      else if (__m == __n)
	return _Tp{1};
      else if (__m == 0 && __n >= 10)
	return _Tp{0};
      else
	return __stirling_2_recur<_Tp>(__n, __m);
    }

  /**
   * Return the Stirling number of the second kind.
   *
   * @todo Look into asymptotic solutions.
   */
  template<typename _Tp>
    _Tp
    __log_stirling_2(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return -std::numeric_limits<_Tp>::infinity();
      else if (__m == __n)
	return _Tp{0};
      else if (__m == 0 && __n >= 1)
	return -std::numeric_limits<_Tp>::infinity();
      else
	return std::log(__stirling_2<_Tp>(__n, __m));
    }

  /**
   * Return the Stirling number of the first kind by series expansion.
   * N.B. This seems to be a total disaster.
   */
  template<typename _Tp>
    _Tp
    __stirling_1_series(unsigned int __n, unsigned int __m)
    {
      if (2 * __n - __m > std::__detail::_S_num_factorials<_Tp> / 2)
	{
	  auto _S1 = _Tp{0};
	  for (auto __k = 0u; __k <= __n - __m; ++__k)
	    {
	      auto __nmk = __n - __m + __k;
	      auto __lbc1 = std::__detail::__log_bincoef<_Tp>(__n - 1 + __k, __nmk);
	      auto __slbc1 = std::__detail::__log_bincoef_sign<_Tp>(__n - 1 + __k, __nmk);
	      auto __lbc2 = std::__detail::__log_bincoef<_Tp>(2 * __n - __m, __nmk);
	      auto __slbc2 = std::__detail::__log_bincoef_sign<_Tp>(2 * __n - __m, __nmk);
	      _S1 += (__k & 1 ? _Tp{-1} : _Tp{1}) * __slbc1 * __slbc2
		   * std::exp(__lbc1 + __lbc2 + __log_stirling_2<_Tp>(__nmk, __k));
	    }
	  return _S1;
	}
      else
	{
	  auto _S1 = _Tp{0};
	  for (auto __k = 0u; __k <= __n - __m; ++__k)
	    {
	      auto __nmk = __n - __m + __k;
	      _S1 += (__k & 1 ? _Tp{-1} : _Tp{1})
		   * std::__detail::__bincoef<_Tp>(__n - 1 + __k, __nmk)
		   * std::__detail::__bincoef<_Tp>(2 * __n - __m, __nmk)
		   * __stirling_2<_Tp>(__nmk, __k);
	    }
	  // @todo Only round if the sum is less than
	  // the maximum representable integer.
	  // Find or make a tool for this.
	  return std::nearbyint(_S1);
	}
    }

  /**
   * Return the Stirling number of the first kind by recursion.
   * The recursion is
   * @f[
   *   S_{n+1}^{(m)} = S_n^{(m-1)} - n S_n^{(m)} \mbox{ or }
   *   S_n^{(m)} = S_{n-1}^{(m-1)} - (n-1) S_{n-1}^{(m)}
   * @f]
   * with starting values
   * @f[
   *   S_0^{(0\rightarrow m)} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   S_{0\rightarrow n}^{(0)} = {1, 0, 0, ..., 0}
   * @f]
   */
  template<typename _Tp>
    _Tp
    __stirling_1_recur(unsigned int __n, unsigned int __m)
    {
      if (__n == 0)
	return _Tp(__m == 0);
      else if (__m == 0)
	return _Tp(__n == 0);
      else
	{
	  std::vector<_Tp> _Sold(__m + 1), _Snew(__m + 1);
	  _Sold[1] = _Tp{1};
	  if (__n == 1)
	    return _Sold[__m];
	  for (auto __in = 1u; __in <= __n; ++__in)
	    {
	      for (auto __im = 1u; __im <= __m; ++__im)
		_Snew[__im] = _Sold[__im - 1] - __in * _Sold[__im];
	      std::swap(_Sold, _Snew);
	    }
	  return _Snew[__m];
	}
    }

  /**
   * Return the Stirling number of the first kind from lookup
   * or by series expansion if terms of Stirling numbers of the second kind.
   *
   * @todo Look into asymptotic solutions.
   */
  template<typename _Tp>
    _Tp
    __stirling_1(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return _Tp{0};
      else if (__m == __n)
	return _Tp{1};
      else if (__m == 0 && __n >= 1)
	return _Tp{0};
      else
        return __stirling_1_recur<_Tp>(__n, __m);
    }

  /**
   * Return the logarithm of the absolute value of Stirling number
   * of the first kind.
   */
  template<typename _Tp>
    _Tp
    __log_stirling_1(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return -std::numeric_limits<_Tp>::infinity();
      else if (__m == __n)
	return _Tp{0};
      else if (__m == 0 && __n >= 1)
	return -std::numeric_limits<_Tp>::infinity();
      else
	return std::log(std::abs(__stirling_1<_Tp>(__n, __m)));
    }

  /**
   * Return the sign of the exponent of the logarithm of the Stirling number
   * of the first kind.
   */
  template<typename _Tp>
    _Tp
    __log_stirling_1_sign(unsigned int __n, unsigned int __m)
    { return (__n + __m) & 1 ? _Tp{-1} : _Tp{+1}; }


template<typename _Tp>
  void
  test_bernoulli(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "\n Bernoulli numbers\n";
    for (auto n = 0u; n <= 200; ++n)
      std::cout << ' ' << std::setw(4) << n
		<< ' ' << std::setw(width) << __bernoulli_series<_Tp>(n)
		<< '\n';

    std::cout << "\n Bernoulli polynomials\n";
    for (auto n = 0u; n <= 50; ++n)
      {
	std::cout << '\n' << ' ' << std::setw(4) << n << '\n';
	for (auto i = 0u; i <= 50; ++i)
	  {
	    auto x = _Tp{0.1Q} * i;
	    auto _B_n = __bernoulli(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << _B_n
		      << '\n';
	  }
      }

    std::cout << "\n Euler numbers\n";
    for (auto n = 0u; n <= 200; ++n)
      std::cout << ' ' << std::setw(4) << n
		<< ' ' << std::setw(width) << __euler<_Tp>(n)
		<< '\n';

    std::cout << "\n Euler polynomials\n";
    for (auto n = 0u; n <= 50; ++n)
      {
	std::cout << '\n' << ' ' << std::setw(4) << n << '\n';
	for (auto i = 0u; i <= 50; ++i)
	  {
	    auto x = _Tp{0.1Q} * i;
	    auto _E_n = __euler(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << _E_n
		      << '\n';
	  }
      }

    std::cout << "\n Stirling numbers of the second kind";
    for (auto n = 0u; n <= 100; ++n)
      {
	std::cout << '\n';
	for (auto m = 0u; m <= n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		    << ' ' << std::setw(width) << __stirling_2_series<_Tp>(n, m)
		    << ' ' << std::setw(width) << __stirling_2_recur<_Tp>(n, m)
		    << '\n';
      }

    std::cout << "\n Stirling numbers of the first kind";
    for (auto n = 0u; n <= 100; ++n)
      {
	std::cout << '\n';
	for (auto m = 0u; m <= n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		    << ' ' << std::setw(width) << __stirling_1_series<_Tp>(n, m)
		    << ' ' << std::setw(width) << __stirling_1_recur<_Tp>(n, m)
		    << '\n';
      }
/*
 */
  }

int
main()
{
  //test_bernoulli(0.0F);

  test_bernoulli(0.0);

  test_bernoulli(0.0L);

  test_bernoulli(0.0Q);

  return 0;
}

