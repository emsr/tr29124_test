/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <vector>
#include <ext/float128_math.h>
#include <ext/float128_io.h>

#include <wrap_burkhardt.h>


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
   * According to Modern Computer Arithmetic a stable recursion
   * for the Bernoulli numbers is
   * @f[
   *    \sum_{k=0}^{n} \frac{C_k}{(2n+1-2k)!4^{n-k}} = \frac{1}{(2n)! 4^n}
   * @f]
   * where @f$ C_k = B_{2k} / (2k)! @f$ is the scaled Bernoulli number.
   */
  template<typename _Tp>
    _Tp
    __bernoulli_scaled_recur(unsigned int __n)
    {
      return _Tp{0};
    }


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
	    __fact *= __k / emsr::pi_v<_Tp>;
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
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  auto _B_n = std::__detail::__bernoulli<_Tp>(0);
	  auto __binomial = _Tp{1};
	  for (auto __k = 1u; __k <= __n; ++__k)
	    {
	      __binomial *= _Tp(__n + 1 - __k) / _Tp(__k);
	      _B_n = __x * _B_n + __binomial
		   * std::__detail::__bernoulli<_Tp>(__k);
	    }
	  return _B_n;
	}
    }

  /**
   * Return the periodic Bernoulli polynomial @f$ \tilde{B}_n(x) @f$
   * of order n at argument x.
   * @f[
   *   \tilde{B}_n(x) = B_n(x), 0 <= x < 1
   * @f]
   * @f[
   *   \tilde{B}_n(x+1) = \tilde{B}_n(x), x >= 1
   * @f]
   */
  template<typename _Tp>
    _Tp
    __bernoulli_period(unsigned int __n, _Tp __x)
    {
      int __p = int(__x);
      return __bernoulli(__n, __x - _Tp(__p));
    }

  /**
   * Return the generalized Bernoulli polynomial @f$ B^{(a)}_n(x) @f$
   * of order n, degree a, at argument x.
   *
   * A recursion is:
   * @f[
   *   aB^{(a+1)}_n(x) = (a-n)B^{(a)}_n(x) + k(x-a)B^{(a)}_{n-1}(x)
   * @f]
   * Starting with @f$ B^{(a)}_0(x) = 1 @f$, @f$  @f$.
   */
  template<typename _Tp>
    _Tp
    __bernoulli(unsigned int __n, _Tp __a, _Tp __x)
    {
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  /* Not so easy after all :-( */;
	}
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
      using std::__detail::__binomial;
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
		    _En[__i] -= __binomial<_Tp>(__i, __j)
			      * _En[__i - __j];
		}
	    }
	  return _En[__n];
	}
    }

  /**
   * @brief This returns Euler number @f$ E_n @f$
   * with an asymptotic formula:
   * @f[
   *    E_{2n} = (-1)^n8\sqrt{\frac{n}{\pi}}
   *   \left( \frac{4n}{\pi e} \frac{480n^2+9}{480n^2-1} \right)^{2n}
   * @f]
   *
   * @param __n the order n of the Euler number.
   * @return  The Euler number of order n.
   */
  template<typename _Tp>
    inline _Tp
    __euler_asymp(unsigned int __n)
    {
      if (__n & 1)
	return _Tp{0};
      else
	{
	  const auto _S_e = emsr::e_v<_Tp>;
	  const auto _S_pi = emsr::pi_v<_Tp>;
	  const auto __n2 = _Tp(__n * __n);
	  return emsr::parity<_Tp>(__n / 2) * _Tp{8}
		* std::sqrt(_Tp(__n) / _S_pi)
		* std::pow(_Tp(4 * __n) / (_S_pi * _S_e)
			 * (__n2 + _Tp{9} / _Tp{480})
			 / (__n2 - _Tp{1} / _Tp{480}), _Tp(2 * __n));
	}
    }

  /**
   * @brief This returns Euler number @f$ E_n @f$.
   *
   * @param __n the order n of the Euler number.
   * @return  The Euler number of order n.
   */
  template<typename _Tp>
    inline _Tp
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
   */
  template<typename _Tp>
    _Tp
    __euler(unsigned int __n, _Tp __x)
    {
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else
	{
	  auto __bx1 = __bernoulli(__n + 1, __x );
	  auto __bx2 = __bernoulli(__n + 1, _Tp{0.5L} * __x );

	  auto _E_n = _Tp{2} * (__bx1 - __bx2 * std::pow(_Tp{2}, _Tp(__n + 1)))
		    / _Tp(__n + 1);

	  return _E_n;
	}
    }

  /**
   * Return the periodic Euler polynomial @f$ \tilde{E}_n(x) @f$
   * of order n at argument x.
   *
   * @f[
   *   \tilde{E}_n(x) = E_n(x), 0 <= x < 1
   * @f]
   * @f[
   *   \tilde{E}_n(x+1) = -\tilde{E}_n(x), x >= 1
   * @f]
   */
  template<typename _Tp>
    _Tp
    __euler_period(unsigned int __n, _Tp __x)
    {
      int __p = int(__x);
      return ((__p & 1) ? -1 : +1) * __euler(__n, __x - _Tp(__p));
    }

  /**
   * Return the Stirling number of the second kind by series expansion.
   * The series is:
   * @f[
   *   \sigma_n^{(m)} = \sum_{k=0}^{m}\frac{(-1)^{m-k}k^n}{(m-k)!k!}
   * @f]
   *
   * @todo Find a way to predict the maximum Stirling number supported
   *       for a given type.
   */
  template<typename _Tp>
    _Tp
    __stirling_2_series(unsigned int __n, unsigned int __m)
    {
      using std::__detail::__log_factorial;
      using std::__detail::__factorial;
      if (__m > std::__detail::_S_num_factorials<_Tp>)
	{
	  auto _S2 = _Tp{0};
	  for (auto __k = 0u; __k <= __m; ++__k)
	    {
	      auto __lf1 = __log_factorial<_Tp>(__k);
	      auto __lf2 = __log_factorial<_Tp>(__m - __k);
	      _S2 += (((__m - __k) & 1) ? _Tp{-1} : _Tp{1})
		   * std::exp(__n * std::log(__k) - __lf1 - __lf2);
	    }
	  return _S2;
	}
      else
	{
	  auto _S2 = _Tp{0};
	  for (auto __k = 0u; __k <= __m; ++__k)
	    {
	      _S2 += (((__m - __k) & 1) ? _Tp{-1} : _Tp{1})
		   * std::pow(__k, __n)
		   / __factorial<_Tp>(__k)
		   / __factorial<_Tp>(__m - __k);
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
   *   \newcommand{\stirling}[2]{\genfrac{\{}{\}}{0pt}{0}{#1}{#2}}
   *   \stirling{n}{m} = m \stirling{n-1}{m} + \stirling{n-1}{m-1}
   * @f]
   * with starting values
   * @f[
   *   \newcommand{\stirling}[2]{\genfrac{\{}{\}}{0pt}{0}{#1}{#2}}
   *   \stirling{0}{0\rightarrow m} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   \newcommand{\stirling}[2]{\genfrac{\{}{\}}{0pt}{0}{#1}{#2}}
   *   \stirling{0\rightarrow n}{0} = {1, 0, 0, ..., 0}
   * @f]
   *
   * The Stirling number of the second kind is denoted by other symbols
   * in the literature: 
   * @f$ \sigma_n^{(m)} @f$, @f$ \textit{S}_n^{(m)} @f$ and others.
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
   * The series is:
   * @f[
   *   \sigma_n^{(m)} = \sum_{k=0}^{m}\frac{(-1)^{m-k}k^n}{(m-k)!k!}
   * @f]
   *
   * @todo Find asymptotic expressions for the Stirling numbers.
   */
  template<typename _Tp>
    _Tp
    __stirling_2(unsigned int __n, unsigned int __m)
    {
      if (__m > __n)
	return _Tp{0};
      else if (__m == __n)
	return _Tp{1};
      else if (__m == 0 && __n >= 1)
	return _Tp{0};
      else
	return __stirling_2_recur<_Tp>(__n, __m);
    }

  /**
   * Return the Stirling number of the second kind.
   *
   * @todo Find asymptotic expressions for the Stirling numbers.
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
   * Maybe accumulate the positive and negative terms separately
   * and add them att the end.
   */
  template<typename _Tp>
    _Tp
    __stirling_1_series(unsigned int __n, unsigned int __m)
    {
      using std::__detail::__log_binomial;
      using std::__detail::__log_binomial_sign;
      using std::__detail::__binomial;
      if (2 * __n - __m > std::__detail::_S_num_factorials<_Tp> / 2)
	{
	  auto _S1 = _Tp{0};
	  for (auto __k = 0u; __k <= __n - __m; ++__k)
	    {
	      const auto __nmpk = __n - __m + __k;
	      const auto __nmmk = __n - __m - __k;
	      const auto __lbc1 = __log_binomial<_Tp>(__n - 1 + __k, __nmpk);
	      const auto __slbc1 = __log_binomial_sign<_Tp>(__n - 1 + __k, __nmpk);
	      const auto __lbc2 = __log_binomial<_Tp>(2 * __n - __m, __nmmk);
	      const auto __slbc2 = __log_binomial_sign<_Tp>(2 * __n - __m, __nmmk);
	      _S1 += emsr::parity<_Tp>(__k) * __slbc1 * __slbc2
		   * std::exp(__lbc1 + __lbc2 + __log_stirling_2<_Tp>(__nmpk, __k));
	    }
	  return _S1;
	}
      else
	{
	  auto _S1 = _Tp{0};
	  for (auto __k = 0u; __k <= __n - __m; ++__k)
	    {
	      const auto __nmpk = __n - __m + __k;
	      const auto __nmmk = __n - __m - __k;
	      _S1 += emsr::parity<_Tp>(__k)
		   * __binomial<_Tp>(__n - 1 + __k, __nmpk)
		   * __binomial<_Tp>(2 * __n - __m, __nmmk)
		   * __stirling_2<_Tp>(__nmpk, __k);
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
   * Return the Stirling number of the first kind.
   *
   * The Stirling numbers of the first kind are the coefficients of
   * the Pocchammer polynomials:
   * @f[
   *   (x)_n = \sum_{k=0}^{n} S_n^{(k)} x^k
   * @f]
   *
   * The recursion is
   * @f[
   *   S_{n+1}^{(m)} = S_n^{(m-1)} - n S_n^{(m)} \mbox{ or }
   * @f]
   * with starting values
   * @f[
   *   S_0^{(0\rightarrow m)} = {1, 0, 0, ..., 0}
   * @f]
   * and
   * @f[
   *   S_{0\rightarrow n}^{(0)} = {1, 0, 0, ..., 0}
   * @f]
   *
   * @todo Find asymptotic expressions for the Stirling numbers.
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
    inline _Tp
    __log_stirling_1_sign(unsigned int __n, unsigned int __m)
    { return ((__n + __m) & 1) ? _Tp{-1} : _Tp{+1}; }

  /**
   * Return the Eulerian number of the first kind by recursion.
   * The recursion is
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    _Tp
    __eulerian_1_recur(unsigned int __n, unsigned int __m)
    {
      if (__m == 0)
	return _Tp{1};
      else if (__m >= __n)
	return _Tp{0};
      else if (__m == __n - 1)
	return _Tp{1};
      else if (__n - __m - 1 < __m) // Symmetry.
	return __eulerian_1_recur<_Tp>(__n, __n - __m - 1);
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<_Tp> _Aold(__m + 1), _Anew(__m + 1);
	  _Aold[0] = _Tp{1};
	  _Anew[0] = _Tp{1};
	  _Anew[1] = _Tp{1};
	  for (auto __in = 3u; __in <= __n; ++__in)
	    {
	      std::swap(_Aold, _Anew);
	      for (auto __im = 1u; __im <= __m; ++__im)
		_Anew[__im] = (__in - __im) * _Aold[__im - 1]
			    + (__im + 1) * _Aold[__im];
	    }
	  return _Anew[__m];
	}
    }

  /**
   * Return the Eulerian number of the first kind.
   * The Eulerian numbers are defined by recursion:
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    inline _Tp
    __eulerian_1(unsigned int __n, unsigned int __m)
    { return __eulerian_1_recur<_Tp>(__n, __m); }

  /**
   * Return a vector Eulerian numbers of the first kind by recursion.
   * The recursion is
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    std::vector<_Tp>
    __eulerian_1_recur(unsigned int __n)
    {
      if (__n == 0)
	return std::vector<_Tp>(1, _Tp{1});
      //else if (__m == __n - 1)
	//return _Tp{1};
      //else if (__n - __m - 1 < __m) // Symmetry.
	//return __eulerian_1_recur<_Tp>(__n, __n - __m - 1);
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<_Tp> _Aold(__n + 1), _Anew(__n + 1);
	  _Aold[0] = _Anew[0] = _Tp{1};
	  _Anew[1] = _Tp{1};
	  for (auto __in = 3u; __in <= __n; ++__in)
	    {
	      std::swap(_Aold, _Anew);
	      for (auto __im = 1u; __im <= __n; ++__im)
		_Anew[__im] = (__in - __im) * _Aold[__im - 1]
			    + (__im + 1) * _Aold[__im];
	    }
	  return _Anew;
	}
    }

  /**
   * Return a vector Eulerian numbers of the first kind.
   * The Eulerian numbers are defined by recursion:
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    inline std::vector<_Tp>
    __eulerian_1(unsigned int __n)
    { return __eulerian_1_recur<_Tp>(__n); }

  /**
   * Return the Eulerian number of the second kind by recursion.
   * The recursion is:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    _Tp
    __eulerian_2_recur(unsigned int __n, unsigned int __m)
    {
      if (__m == 0)
	return _Tp{1};
      else if (__m >= __n)
	return _Tp{0};
      else if (__n == 0)
	return _Tp{1};
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<_Tp> _Aold(__m + 1), _Anew(__m + 1);
	  _Aold[0] = _Tp{1};
	  _Anew[0] = _Tp{1};
	  _Anew[1] = _Tp{2};
	  for (auto __in = 3u; __in <= __n; ++__in)
	    {
	      std::swap(_Aold, _Anew);
	      for (auto __im = 1u; __im <= __m; ++__im)
		_Anew[__im] = (2 * __in - __im - 1) * _Aold[__im - 1]
			    + (__im + 1) * _Aold[__im];
	    }
	  return _Anew[__m];
	}
    }

  /**
   * Return the Eulerian number of the second kind.
   * The Eulerian numbers of the second kind are defined by recursion:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    inline _Tp
    __eulerian_2(unsigned int __n, unsigned int __m)
    { return __eulerian_2_recur<_Tp>(__n, __m); }

  /**
   * Return the Eulerian number of the second kind by recursion.
   * The recursion is:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    std::vector<_Tp>
    __eulerian_2_recur(unsigned int __n)
    {
      if (__n == 0)
	return std::vector<_Tp>(1, _Tp{1});
      //else if (__m >= __n)
	//return _Tp{0};
      //else if (__n == 0)
	//return _Tp{1};
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<_Tp> _Aold(__n + 1), _Anew(__n + 1);
	  _Aold[0] = _Anew[0] = _Tp{1};
	  _Anew[1] = _Tp{2};
	  for (auto __in = 3u; __in <= __n; ++__in)
	    {
	      std::swap(_Aold, _Anew);
	      for (auto __im = 1u; __im <= __n; ++__im)
		_Anew[__im] = (2 * __in - __im - 1) * _Aold[__im - 1]
			    + (__im + 1) * _Aold[__im];
	    }
	  return _Anew;
	}
    }

  /**
   * Return the Eulerian number of the second kind.
   * The Eulerian numbers of the second kind are defined by recursion:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename _Tp>
    inline std::vector<_Tp>
    __eulerian_2(unsigned int __n)
    { return __eulerian_2_recur<_Tp>(__n); }

  /**
   * Return the Lah number by downward recurrence.
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename _Tp>
    _Tp
    __lah_recur(unsigned int __n, unsigned int __k)
    {
      if (__k > __n)
	return _Tp{0};
      else if (__n == 0)
	return (__k == 0 ? _Tp{1} : _Tp{0});
      else
	{
	  _Tp _Lnn = 1;
	  for (unsigned int __i = 1u; __i <= __n - __k; ++__i)
	    _Lnn *= _Tp(__n - __i + 1) * _Tp(__n - __i) / _Tp(__i);
	  return _Lnn;
	}
    }

  /**
   * Return a vector of Lah numbers by downward recurrence.
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename _Tp>
    std::vector<_Tp>
    __lah_recur(unsigned int __n)
    {
      if (__n == 0)
	return std::vector<_Tp>(1, _Tp{1});
      else
	{
	  std::vector<_Tp> _L(__n + 1);
	  _Tp _Lnn = 1;
	  _L[__n] = _Lnn;
	  for (unsigned int __i = 1u; __i <= __n; ++__i)
	    {
	      _Lnn *= _Tp(__n - __i + 1) * _Tp(__n - __i) / _Tp(__i);
	      _L[__n - __i] = _Lnn;
	    }
	  return _L;
	}
    }

  /**
   * Return a vector of Lah numbers.
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename _Tp>
    inline std::vector<_Tp>
    __lah(unsigned int __n)
    { return __lah_recur<_Tp>(__n); }

  /**
   * Return a vector of the Bell numbers by summation.
   * @f[
   *   B(n) = \sum_{k=0}^{n}S_n^{(k)}
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename _Tp>
    std::vector<_Tp>
    __bell_series(unsigned int __n)
    {
      std::vector<_Tp> __bell(__n + 1);
      __bell[0] = _Tp{1};

      /// @todo Test for blowup in Bell number summation.
      for (unsigned int __i = 1; __i <= __n; ++__i)
	for (unsigned int __j = 1; __j <= __i; ++__j)
	  __bell[__i] += __bell[__i - __j]
		       * __gnu_cxx::binomial<_Tp>(__i - 1, __j - 1);

      return __bell;
    }

  /**
   * Return a vector of the Bell numbers.
   */
  template<typename _Tp>
    inline std::vector<_Tp>
    __bell(unsigned int __n)
    { return __bell_series<_Tp>(__n); }

  /**
   * Evaluate the Bell polynomial
   * @f[
   *   B(n) = \sum_{k=0}^{n}S_n^{(k)}x^k
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename _Tp, typename _Up>
    inline _Up
    __bell(unsigned int __n, _Up __x)
    {
      const auto _Sn = __stirling_2<_Tp>(__n);
      auto __bell = _Sn[__n];
      for (unsigned int __i = 1; __i < __n; ++__i)
	__bell = _Sn[__n - __i] + __x * __bell;
      return __bell;
    }

template<typename _Tp>
  void
  test_bernoulli(_Tp proto = _Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "\n\n Bell numbers\n";
    for (auto n = 1u; n <= 20; ++n)
      {
	auto bell = __bell_series<long long>(n);
        auto bellv = burkhardt::bell(n);
	std::cout << '\n';
	for (auto k = 0u; k <= n; ++k)
	  std::cout << ' ' << std::setw(4) << n
	       	    << ' ' << std::setw(4) << k
		    << ' ' << std::setw(width) << bell[k]
		    << ' ' << std::setw(width) << bellv[k]
		    << '\n';
      }

    std::cout << "\n\n Lah numbers\n";
    for (auto n = 1u; n <= 20; ++n)
      {
	const auto lahv = __lah<_Tp>(n);
	std::cout << '\n';
	for (auto k = 0u; k <= n; ++k)
	  std::cout << ' ' << std::setw(4) << n
	       	    << ' ' << std::setw(4) << k
		    << ' ' << std::setw(width) << __lah_recur<_Tp>(n, k)
		    << ' ' << std::setw(width) << lahv[k]
		    << '\n';
      }

    std::cout << "\n\n Bernoulli numbers\n";
    for (auto n = 0u; n <= 200; ++n)
      std::cout << ' ' << std::setw(4) << n
		<< ' ' << std::setw(width) << __bernoulli_series<_Tp>(n)
		<< '\n';

    std::cout << "\n\n Bernoulli polynomials\n";
    for (auto n = 0u; n <= 50; ++n)
      {
	std::cout << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (auto i = 0u; i <= 50; ++i)
	  {
	    auto x = del * i;
	    auto _B_n = __bernoulli(n, x);
	    auto _tildeB_n = __bernoulli_period(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << _B_n
		      << ' ' << std::setw(width) << _tildeB_n
		      << '\n';
	  }
      }

    std::cout << "\n\n Euler numbers\n";
    for (auto n = 0u; n <= 200; ++n)
      std::cout << ' ' << std::setw(4) << n
		<< ' ' << std::setw(width) << __euler<_Tp>(n)
		<< '\n';

    std::cout << "\n\n Euler polynomials\n";
    for (auto n = 0u; n <= 50; ++n)
      {
	std::cout << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = _Tp{1} / _Tp{10};
	for (auto i = 0u; i <= 50; ++i)
	  {
	    auto x = del * i;
	    auto _E_n = __euler(n, x);
	    auto _tildeE_n = __euler_period(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << _E_n
		      << ' ' << std::setw(width) << _tildeE_n
		      << '\n';
	  }
      }

    std::cout << "\n\n Stirling numbers of the second kind";
    for (auto n = 0u; n <= 100; ++n)
      {
        const auto s2v = std::__detail::__stirling_2<_Tp>(n);
	std::cout << '\n';
	for (auto m = 0u; m <= n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		    << ' ' << std::setw(width) << __stirling_2_series<_Tp>(n, m)
		    << ' ' << std::setw(width) << __stirling_2_recur<_Tp>(n, m)
		    << ' ' << std::setw(width) << s2v[m]
		    << '\n';
      }

    std::cout << "\n\n Stirling numbers of the first kind";
    for (auto n = 0u; n <= 100; ++n)
      {
        const auto s1v = std::__detail::__stirling_1<_Tp>(n);
	std::cout << '\n';
	for (auto m = 0u; m <= n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		    << ' ' << std::setw(width) << __stirling_1_series<_Tp>(n, m)
		    << ' ' << std::setw(width) << __stirling_1_recur<_Tp>(n, m)
		    << ' ' << std::setw(width) << s1v[m]
		    << '\n';
      }

    std::cout << "\n\n Eulerian numbers of the first kind";
    for (auto n = 1u; n <= 10; ++n)
      {
        const auto e1v = __eulerian_1<_Tp>(n);
	std::cout << '\n';
	for (auto m = 0u; m < n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		  //  << ' ' << std::setw(width) << __eulerian_1_series<_Tp>(n, m)
		    << ' ' << std::setw(width) << __eulerian_1_recur<_Tp>(n, m)
		    << ' ' << std::setw(width) << burkhardt::eulerian_1(n, m)
		    << ' ' << std::setw(width) << e1v[m]
		    << '\n';
      }

    std::cout << "\n\n Eulerian numbers of the second kind";
    for (auto n = 1u; n <= 10; ++n)
      {
        const auto e2v = __eulerian_2<_Tp>(n);
	std::cout << '\n';
	for (auto m = 0u; m < n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		  //  << ' ' << std::setw(width) << __eulerian_2_series<_Tp>(n, m)
		    << ' ' << std::setw(width) << __eulerian_2_recur<_Tp>(n, m)
		  //  << ' ' << std::setw(width) << burkhardt::eulerian_2(n, m)
		    << ' ' << std::setw(width) << e2v[m]
		    << '\n';
      }
  }

int
main()
{
  //test_bernoulli(0.0F);

  test_bernoulli(0.0);

  test_bernoulli(0.0L);

  //test_bernoulli(0.0Q);

  return 0;
}

