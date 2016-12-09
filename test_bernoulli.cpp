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
      static constexpr std::size_t __len = 24;
      static constexpr _Tp
      __num[__len]
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

      if (__n == 1)
	return -frac<_Tp, 1, 2>;

      //  Take care of the rest of the odd ones.
      if (__n % 2 == 1)
	return _Tp(0);

      //  Take care of some small evens that are painful for the series.
      if (__n < __len)
	return __num[__n];


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
      static constexpr std::size_t __len = 0;
      static constexpr _Tp
      __num[__len]
      {
      };

      if (__n == 0)
	return _Tp{1};

      if (__n == 1)
	return _Tp{0};

      if (__n == 2)
        return _Tp{-1};

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
/*
*/
  }

int
main()
{
  test_bernoulli(0.0);

  test_bernoulli(0.0L);

  test_bernoulli(0.0Q);

  return 0;
}

