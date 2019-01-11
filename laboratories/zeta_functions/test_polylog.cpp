/*
$HOME/bin_specfun/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_polylog test_polylog.cpp -lquadmath -Lwrappers/debug -lwrap_cephes
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_polylog > test_polylog.txt

LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH $HOME/bin_specfun/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_polylog test_polylog.cpp -lquadmath -Lwrappers/debug -lwrap_cephes

$HOME/bin/bin/g++ -std=gnu++2a -g -Wall -Wextra -Wno-psabi -I. -o test_polylog test_polylog.cpp -lquadmath -Lwrappers/debug -lwrap_cephes
PATH=wrappers/debug:$PATH ./test_polylog > test_polylog.txt
LD_LIBRARY_PATH=wrappers/debug:$LD_LIBRARY_PATH ./test_polylog > test_polylog.txt
*/

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <bits/float128_io.h>

#include "wrap_cephes.h"

  /**
   * This class manages the termination of series.
   * Termination conditions involve both a maximum iteration count
   * and a relative precision.
   */
  template<typename _Tp>
    class _Terminator
    {
    private:

      using _Real = std::__detail::__num_traits_t<_Tp>;
      const std::size_t _M_max_iter;
      std::size_t _M_curr_iter;
      _Real _M_toler;

    public:

      _Terminator(std::size_t __max_iter, _Real __mul = _Real{1})
      : _M_max_iter(__max_iter), _M_curr_iter{0},
	_M_toler(std::abs(__mul) * std::numeric_limits<_Real>::epsilon())
      { }

      /// Return the current number of terms summed.
      std::size_t
      num_terms() const
      { return this->_M_curr_iter; }

      /// Detect if the sum should terminate either because the incoming term
      /// is small enough or the maximum number of terms has been reached.
      bool
      operator()(_Tp __term, _Tp __sum)
      {
	if (this->_M_curr_iter >= this->_M_max_iter
	    || ++this->_M_curr_iter == this->_M_max_iter)
	  return true;
	else if (std::abs(__term) < this->_M_toler * std::abs(__sum))
	  return true;
	else
	  return false;
      }
    };

  /**
   * This class manages the termination of asymptotic series.
   * In particular, this termination watches for the growth of the sequence
   * of terms to stop the series.
   *
   * Termination conditions involve both a maximum iteration count
   * and a relative precision.
   */
  template<typename _Tp>
    class _AsympTerminator
    {
    private:

      using _Real = std::__detail::__num_traits_t<_Tp>;
      const std::size_t _M_max_iter;
      std::size_t _M_curr_iter;
      _Real _M_toler;
      _Real _M_prev_term = std::numeric_limits<_Real>::max();
      bool _M_stop_asymp = false;

    public:

      _AsympTerminator(std::size_t __max_iter, _Real __mul = _Real{1})
      : _M_max_iter(__max_iter), _M_curr_iter{0},
	_M_toler(std::abs(__mul) * std::numeric_limits<_Real>::epsilon())
      { }

      /// Filter a term before applying it to the sum.
      _Tp
      operator<<(_Tp __term)
      {
	if (std::abs(__term) > this->_M_prev_term)
	  {
	    this->_M_stop_asymp = true;
	    return _Tp{0};
	  }
	else
	  return __term;
      }

      /// Return the current number of terms summed.
      std::size_t
      num_terms() const
      { return this->_M_curr_iter; }

      /// Detect if the sum should terminate either because the incoming term
      /// is small enough or because the terms are starting to grow or
      //  the maximum number of terms has been reached.
      bool
      operator()(_Tp __term, _Tp __sum)
      {
	if (this->_M_stop_asymp)
	  return true;
	else
	  {
	    const auto __aterm = std::abs(__term);
	    this->_M_stop_asymp = (__aterm > this->_M_prev_term);
	    this->_M_prev_term = __aterm;
	    if (this->_M_curr_iter >= this->_M_max_iter
	    || ++this->_M_curr_iter == this->_M_max_iter)
	      return true;
	    else if (__aterm < this->_M_toler * std::abs(__sum))
	      return true;
	    else if (this->_M_stop_asymp)
	      return true;
	    else
	      return false;
	  }
      }
    };

  /**
   * Compute the polylogarithm for nonpositive integer order.
   * @f[
   *    Li_{-n}(z) = \sum_{k=0}^{n} S(n+1,k+1) \left(\frac{z}{1-z}\right)^{k+1}
   *           \mbox{    } n = 0,1,2, ...
   * @f]
   */
  template<typename _Tp>
    _Tp
    __polylog_nonpos_int_1(int __n, _Tp __x)
    {
      if (__x == _Tp{0})
	return _Tp{0};
      else if (__x == _Tp{1})
	return std::numeric_limits<_Tp>::infinity();
      else
	{
	  const auto __arg = __x / (_Tp{1} - __x);
	  auto __fact = __arg;
	  auto __term = __arg * __gnu_cxx::stirling_2<_Tp>(1 - __n, 1);
	  auto __sum = __term;
	  for (int __k = 1; __k <= -__n; ++__k)
	    {
	      __fact *= __k * __arg;
	      __term = __fact * __gnu_cxx::stirling_2<_Tp>(1 - __n, 1 + __k);
	      __sum += __term;
	    }
	  return __sum;
	}
    }

  /**
   * Compute the polylogarithm for negative integer order.
   * @f[
   *    Li_{-n}(z) = (-1)^{n+1}
   *                \sum_{k=0}^{n} S(n+1,k+1) \left(\frac{-1}{1-z}\right)^{k+1}
   *           \mbox{    } n = 0,1,2, ...
   * @f]
   * where @f$ S(n,k) @f$ are the Sterling numbers of the second kind.
   */
  template<typename _Tp>
    _Tp
    __polylog_nonpos_int_2(int __n, _Tp __x)
    {
      if (__x == _Tp{0})
	return _Tp{0};
      else if (__x == _Tp{1})
	return std::numeric_limits<_Tp>::infinity();
      else
	{
	  const auto __arg = _Tp{-1} / (_Tp{1} - __x);
	  auto __fact = __arg;
	  auto __term = __fact * __gnu_cxx::stirling_2<_Tp>(1 - __n, 1);
	  auto __sum = __term;
	  for (int __k = 1; __k <= -__n; ++__k)
	    {
	      __fact *= __k * __arg;
	      __term = __fact * __gnu_cxx::stirling_2<_Tp>(1 - __n, 1 + __k);
	      __sum += __term;
	    }
	  __sum *= __gnu_cxx::__parity<_Tp>(1 - __n);
	  return __sum;
	}
    }

  /**
   * Compute the polylogarithm for negative integer order.
   * @f[
   *    Li_{-n}(z) = \frac{1}{(1-z)^{n+1}}\sum_{k=0}^{n-1}
   *                 \sum_{k=0}^{n-1} \left< n \over k \right> z^{n-k}
   *           \mbox{    } n = 1,2, ...
   * @f]
   */
  template<typename _Tp>
    _Tp
    __polylog_nonpos_int_3(int __n, _Tp __x)
    {
      if (__x == _Tp{0})
	return _Tp{0};
      else if (__x == _Tp{1})
	return std::numeric_limits<_Tp>::infinity();
      else
	{
	  auto __fact = _Tp{1} / __x;
	  auto __term = __fact * __gnu_cxx::eulerian_1<_Tp>(-__n, 0);
	  auto __sum = __term;
	  for (int __k = 1; __k < -__n; ++__k)
	    {
	      __fact /= __x;
	      __term = __fact * __gnu_cxx::eulerian_1<_Tp>(-__n, __k);
	      __sum += __term;
	    }
	  __sum *= std::pow(__x / (_Tp{1} - __x), _Tp(1 - __n));
	  return __sum;
	}
    }

  /**
   * Compute the polylogarithm for negative integer order.
   * @f[
   *   Li_{-p}(e^w) = p!(-w)^{-(p+1)}
   *     - \sum_{k=0}^{\infty} \frac{B_{p+2k+q+1}}{(p+2k+q+1)!}
   *                           \frac{(p+2k+q)!}{(2k+q)!}w^{2k+q}
   * @f]
   * where @f$ q = (p+1) mod 2 @f$.
   *
   * @param __n the negative integer index @f$ n = -p @f$.
   * @param __w the argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_neg(int __n, std::complex<_Tp> __w)
    {
      const auto _S_inf = std::numeric_limits<_Tp>::infinity();
      if (__gnu_cxx::__fp_is_zero(__w))
	return std::complex<_Tp>{0};
      else if (__gnu_cxx::__fp_is_equal(__w, _Tp{1}))
	return std::complex<_Tp>{_S_inf, _Tp{0}};
      else
	{
	  const int __p = -__n;
	  const int __pp = 1 + __p;
	  const int __q = (__p & 1) ? 0 : 1;
	  const auto __w2 = __w * __w;
	  auto __wp = (__p & 1) ? std::complex<_Tp>{1} : __w;
	  unsigned int __2k = __q;
	  auto __gam = std::__detail::__factorial<_Tp>(__p + __2k);
	  const auto __pfact = std::__detail::__factorial<_Tp>(__p);
	  auto __res = __pfact * std::pow(-__w, _Tp(-__pp));
	  auto __sum = std::complex<_Tp>{};
	  constexpr unsigned int __maxit = 300;
	  _Terminator<std::complex<_Tp>> __done(__maxit);
	  while (true)
	    {
	      const auto __id = (__p + __2k + 1) / 2;
	      if (__id == std::__detail::_Num_Euler_Maclaurin_zeta)
		break;
	      auto __term = __gam * __wp
		* _Tp(std::__detail::_S_Euler_Maclaurin_zeta[__id]);
	      __sum += __term;
	      if (__done(__term, __sum))
		break;
	      __gam *= _Tp(__p + __2k + 1) / _Tp(__2k + 1)
		     * _Tp(__p + __2k + 2) / _Tp(__2k + 2);
	      __wp *= __w2;
	      __2k += 2;
	    }
	  __res -= __sum;
	  return __res;
	}
    }

  /**
   * This function treats the cases of negative integer index @f$ s = -p @f$
   * which are even.
   *
   * In that case the sine occuring in the expansion occasionally
   * takes on the value zero.
   * We use that to provide an optimized series for p = 2n:
   *
   * @f[
   *   Li_{-p}(e^w) = \Gamma(1+p) (-w)^{-p-1} - A_p(w) - \sigma B_p(w)
   * @f]
   * with
   * @f[
   *   A_p(w) = 2 (2\pi)^{-p-1} \frac{p!}{(2\pi)^{-p/2}}
   *           \left(1 + \frac{w^2}{4\pi^2}\right)^{-(p+1)/2}
   *          \cos\left[(1 + p)ArcTan\left(\frac{2\pi}{w}\right)\right]
   * @f]
   * and 
   * @f[
   *   B_p(w) = - 2 (2\pi)^{-p-1} \sum_{k = 0}^{\infty} 
   *           \frac{(2k + 1 + p)!}{(2k + 1)!}
   *           (-1)^k \left(\frac{w}{2\pi}\right)^{2k+1} [\zeta(2 + 2k + p) - 1]
   * @f]
   * This is suitable for @f$ |w| < 2 \pi @f$
   * The original series is (This might be worthwhile if we use
   * the already present table of the Bernoullis)
   * @f[
   *   Li_{-p}(e^w) = p! (-w)^{-p-1}
   *      - \sigma (2\pi)^p /\pi \sum_{k = 0}^{\infty}
   *       \frac{(2k + 1 + p)!}{(2k + 1)!}
   *       (-1)^k \left(\frac{w}{2\pi}\right)^{2k+1} \zeta(2 + 2k + p)
   * @f]
   * Where @f$ p = 4k (\sigma = 1) @f$ or @f$ p = 4k + 2 (\sigma = -1) @f$.
   *
   * @param __p the integral index @f$ p = 4k @f$ or @f$ p = 4k + 2 @f$.
   * @param __w The complex argument w
   * @return the value of the Polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_neg_even(unsigned int __p, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__w));
      const int __pp = 1 + __p;
      const int __sigma = __p % 4 == 0 ? +1 : -1;
      const auto __lnp = std::__detail::__log_factorial<_Tp>(__p);
      auto __res = std::exp(__lnp - _Tp(__pp) * std::log(-__w));
      auto __wup = __w / _S_2pi;
      auto __wq = -__wup * __wup;
      auto __pref = _Tp{2} * std::pow(_S_2pi, _Tp(-__pp));
      // Subtract the expression A_p(w)
      __pref *= __sigma;
      __res -= __pref
	     * std::exp(__lnp - _Tp{0.5L} * __pp * std::log(_Tp{1} - __wq))
	     * std::cos(_Tp(__pp) * std::atan(_Tp{1} / __wup));
      unsigned int __2k = 0;
      auto __gam = std::__detail::__factorial<_Tp>(1 + __p);
      std::complex<_Tp> __sum;
      constexpr unsigned int __maxit = 400;
      _Terminator<std::complex<_Tp>> __done(__maxit);
      while (true)
	{
	  auto __term = __gam * __wup
		      * std::__detail::__riemann_zeta_m_1<_Tp>(__2k + 2 + __p);
	  __sum += __term;
	  if (__done(__term, __sum))
	    break;
	  __gam *= _Tp(__2k + 2 + __p) / _Tp(__2k + 2)
		 * _Tp(__2k + 3 + __p) / _Tp(__2k + 3);
	  __2k += 2;
	  __wup *= __wq;
	}
      __res -= __pref * __sum;
      return __res;
    }

  /**
   * This function treats the cases of negative integer index @f$ s = -p @f$
   * which are odd.
   *
   * In that case the sine occuring in the expansion occasionally vanishes.
   * We use that to provide an optimized series for @f$ p = 1 + 2k @f$:
   * In the template parameter sigma we transport whether
   * @f$ p = 4k + 1 (\sigma = 1) @f$ or @f$ p = 4k + 3  (\sigma = -1) @f$.
   *
   * @f[
   *   Li_{-p}(e^w) = \Gamma(1+p) (-w)^{-p-1} + \sigma A_p(w) - \sigma B_p(w)
   * @f]
   * with
   * @f[
   *   A_p(w) = 2 (2\pi)^{-p-1} p!
   *          \left(1 + \frac{w^2}{4\pi^2}\right)^{-(p + 1)/2}
   *           \cos((1 + p) ArcTan(2 \pi / w))
   * @f]
   * and 
   * @f[
   *   B_p(w) = 2 (2\pi)^{-p-1}
   *            \sum_{k=0}^{\infty}\frac{(2k + p)!}{(2k)!}
   *      \left(\frac{-w^2}{4\pi^2}\right)^k [\zeta(2k + p + 1) - 1]
   * @f]
   * This is suitable for @f$ |w| < 2 \pi @f$.
   * The original series is (This might be worthwhile if we use
   * the already present table of the Bernoullis)
   * @f[
   *   Li_{-p}(e^w) = p! (-w)^{-p-1}
   *      - 2\sigma (2\pi)^{-p-1} \sum_{k = 0}^{\infty}
   *       \frac{(2k + p)!}{(2k)!}
   *       (-1)^k \left(\frac{w}{2\pi}\right)^{2k} \zeta(1 + 2k + p)
   * @f]
   *
   * @param __p the integral index @f$ p = 4k + 1 @f$ or @f$ p = 4k + 3 @f$.
   * @param __w The complex argument w.
   * @return The value of the Polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_neg_odd(unsigned int __p, std::complex<_Tp> __w)
    {
      const auto _S_2pi = __gnu_cxx::__const_2_pi(std::real(__w));
      const int __pp = 1 + __p;
      const int __sigma = __p % 4 == 1 ? +1 : -1;
      const auto __lnp = std::__detail::__log_factorial<_Tp>(__p);
      auto __res = std::exp(__lnp - _Tp(__pp) * std::log(-__w));
      auto __wup = __w / _S_2pi;
      auto __wq = -__wup * __wup;
      auto __pref = _Tp{2} * std::pow(_S_2pi, _Tp(-__pp));
      __pref *= __sigma;
      // Add the expression A_p(w)
      __res += __pref
	     * std::exp(__lnp - _Tp{0.5L} * __pp * std::log(_Tp{1} - __wq))
	     * std::cos(_Tp(__pp) * std::atan(_Tp{1} / __wup));
      auto __gam = std::__detail::__factorial<_Tp>(__p);
      unsigned int __2k = 0;
      std::complex<_Tp> __sum;
      constexpr unsigned int __maxit = 400;
      _Terminator<std::complex<_Tp>> __done(__maxit);
      while (true)
	{
	  auto __term = __gam * __wup
		      * std::__detail::__riemann_zeta_m_1<_Tp>(__2k + 1 + __p);
	  __sum += __term;
	  if (__done(__term, __sum))
	    break;
	  __gam *= _Tp(__2k + 1 + __p) / _Tp(__2k + 1)
		 * _Tp(__2k + 2 + __p) / _Tp(__2k + 2);
	  __2k += 2;
	  __wup *= __wq;
	}
      __res -= __pref * __sum;
      return __res;
  }

  /**
   * This function treats the cases of negative integer index s
   * and branches accordingly
   *
   * @param __s the integer index s.
   * @param __w The Argument w
   * @return The value of the Polylogarithm evaluated by a suitable function.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_neg_old(int __s, std::complex<_Tp> __w)
    { // negative integer __s
      const auto __p = -__s;
      switch (__p % 2)
      {
      case 0:
	return __polylog_exp_neg_even<_Tp>(__p, __w);
      case 1:
	return __polylog_exp_neg_odd<_Tp>(__p, __w);
      default: // We shouldn't need this.
	return std::complex<_Tp>{};
      }
    }

  /**
   * This function implements the asymptotic series for the polylog.
   * It is given by
   * @f[
   *    2 \sum_{k=0}^{\infty} \zeta(2k) w^{s-2k}/\Gamma(s-2k+1)
   *       -i \pi w^{s-1}/\Gamma(s)
   * @f]
   * for @f$ Re(w) >> 1 @f$
   *
   * Don't check this against Mathematica 8.
   * For real u the imaginary part of the polylog is given by
   * @f$ Im(Li_s(e^u)) = - \pi u^{s-1}/\Gamma(s) @f$.
   * Check this relation for any benchmark that you use.
   *
   * @param __s the real index s.
   * @param __w the large complex argument w.
   * @return the value of the polylogarithm.
   */
  template<typename _Tp>
    std::complex<_Tp>
    __polylog_exp_asymp(_Tp __s, std::complex<_Tp> __w)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__s);
      // wgamma = w^{s-1} / Gamma(s)
      auto __wgamma = std::pow(__w, __s - _Tp{1}) * std::__detail::__gamma_reciprocal(__s);
      auto __res = std::complex<_Tp>(_Tp{0}, -_S_pi) * __wgamma;
      // wgamma = w^s / Gamma(s+1)
      __wgamma *= __w / __s;
      constexpr unsigned int __maxiter = 100;
      _AsympTerminator<std::complex<_Tp>> __done(__maxiter);
      // zeta(0) w^s / Gamma(s + 1)
      std::complex<_Tp> __oldterm = -_Tp{0.5L} * __wgamma;
      __res += _Tp{2} * __oldterm;
      std::complex<_Tp> __term;
      auto __wq = _Tp{1} / (__w * __w);
      int __k = 1;
      while (true)
	{
	  __wgamma *= __wq * (__s + _Tp(1 - 2 * __k)) * (__s + _Tp(2 - 2 * __k));
	  __term = std::__detail::__riemann_zeta<_Tp>(2 * __k) * __wgamma;
	  __res += __done << _Tp{2} * __term;
	  if (__done(_Tp{2} * __term, __res))
	    break;
	  __oldterm = __term;
	  ++__k;
	}
      return __res;
    }

/**
 * Compute the polylogarithm for negative integer order.
 */
template<typename Tp>
  void
  test_polylog_neg_int(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    //std::cout << "\nTest negative integer order\n";
    for (auto n : {-1, -2, -3, -4, -5})
      {
	std::cout << "\n\nNegative integer degree: n = " << n << '\n';
	const auto del = Tp{1} / Tp{20};
	for (int i = -200; i <= 20; ++i)
	  {
	    auto x = del * i;
	    auto Ls_rat1 = __polylog_nonpos_int_1(n, x);
	    auto Ls_rat2 = __polylog_nonpos_int_2(n, x);
	    auto Ls_rat3 = __polylog_nonpos_int_3(n, x);
	    auto Ls_nint = x == Tp{0}
			 ? Tp{0}
			 : [n, x]() -> Tp
			   { auto w = std::log(std::complex<Tp>(x));
			     return std::real(std::__detail::__polylog_exp_neg(n, w)); }();
	    auto Ls_gnu = x == Tp{0}
			? Tp{0}
			: [n, x]() -> Tp
			  { auto w = std::log(std::complex<Tp>(x));
			    return std::real(__polylog_exp_neg_old(n, w)); }();
	    std::cout << ' ' << n
		      << ' ' << x
		      << ' ' << Ls_rat1
		      << ' ' << Ls_rat2
		      << ' ' << Ls_rat3
		      << ' ' << Ls_nint
		      << ' ' << Ls_gnu
		      << '\n';
	  }
      }
    std::cout << '\n' << std::flush;
  }

template<typename Tp>
  void
  test_polylog_0(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    std::cout << "\n\nTest against algebraic function for Li_0\n";
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i <= 10; ++i)
      {
	auto x = del * i;
	auto Li0 = __gnu_cxx::__fp_is_zero(Tp{1} - x)
		 ? std::numeric_limits<Tp>::quiet_NaN()
		 : x / (Tp{1} - x);
	auto Lis_gnu = __gnu_cxx::polylog(Tp{0}, x);
	std::cout << ' ' << x
		  << ' ' << Li0
		  << ' ' << Lis_gnu
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }

template<typename Tp>
  void
  test_polylog_1(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    std::cout << "\n\nTest against algebraic function for Li_1\n";
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i <= 10; ++i)
      {
	auto x = del * i;
	auto Li1 = __gnu_cxx::__fp_is_zero(Tp{1} - x)
		 ? std::numeric_limits<Tp>::quiet_NaN()
		 : -std::log(Tp{1} - x);
	auto Lis_gnu = __gnu_cxx::polylog(Tp{1}, x);
	std::cout << ' ' << x
		  << ' ' << Li1
		  << ' ' << Lis_gnu
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }

template<typename Tp>
  void
  test_polylog_dilog(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    std::cout << "\n\nTest against local dilog\n";
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i <= 10; ++i)
      {
	auto x = del * i;
	auto Ls_dilog = __gnu_cxx::dilog(x);
	auto Ls_gnu = __gnu_cxx::polylog(Tp(2), x);
	std::cout << ' ' << x
		  << ' ' << Ls_dilog
		  << ' ' << Ls_gnu
		  << '\n';
      }
    std::cout << '\n' << std::flush;
  }

template<typename Tp>
  void
  test_polylog_cephes(Tp proto = Tp{})
  {
    std::cout.precision(__gnu_cxx::__digits10(proto));
    std::cout << std::scientific;

    //std::cout << "\nTest against Cephes for integer order\n";
    for (auto n : {0, 1, 2, 3, 4, 5})
      {
	std::cout << "\n\nNon-negative integer degree: n = " << n << '\n';
	const auto del = Tp{1} / Tp{10};
	for (int i = -200; i <= 10; ++i)
	  {
	    auto x = del * i;
	    auto Ls_ceph = cephes::polylog(n, x);
	    auto Ls_gnu = __gnu_cxx::polylog(Tp(n), x);
	    std::cout << ' ' << n
		      << ' ' << x
		      << ' ' << Ls_ceph
		      << ' ' << Ls_gnu
		      << '\n';
	  }
      }
    std::cout << '\n' << std::flush;
  }

template<typename Tp>
  void
  TestPolyLog(Tp proto = Tp{})
  {
    const auto _S_pi = __gnu_cxx::__const_pi(proto);
    const auto _S_2pi = __gnu_cxx::__const_2_pi(proto);

    std::cout.precision(__gnu_cxx::__digits10(proto) - 1);
    std::cout << std::scientific;
    const auto w = std::cout.precision() + 8;

    std::cout << '\n';

    //std::size_t n = 5000;

    // this part of code was for performance testing.
    // the old implementation takes about 2.8s on my core2 and the new one 0.8s
    //     for(std::size_t i = 0; i < n; ++i)
    //     {
    //       Tp x = Tp{10}* static_cast<Tp>(i)/n + 1.1;
    // //      std::cout << std::scientific<<x<<' '<<
    //       std::__detail::__riemann_zeta(x)
    // //      std::tr1::__detail::__riemann_zeta(x)
    //       ;//between 1 and 10 riemann_zeta_glob is called
    // //      <<'\n';
    //     }

    // Something that didn't work in the original implementation
    std::cout << std::riemann_zeta(std::complex<Tp>{5.1, 0.5}) << '\n';
    std::cout << __gnu_cxx::hurwitz_zeta(Tp{5.1}, Tp{0.5}) << '\n';
    std::cout << __gnu_cxx::hurwitz_zeta(Tp{5.1}, std::complex<Tp>{0.5}) << '\n';
    std::cout << std::__detail::__hurwitz_zeta_polylog(Tp{5.1}, std::complex<Tp>{0.5}) << '\n';
    std::cout << std::__detail::__polylog_exp(Tp{2.5}, std::complex<Tp>(Tp{15}, Tp{1})) << '\n';

    for(std::size_t k = 0; k < 32; ++k)
    {
      std::cout << "=======  " << k << "  ==========" << '\n';
      auto w = std::complex<Tp>(Tp{0}, _S_2pi * k / Tp{32});
      std::cout << std::__detail::__polylog_exp(Tp{4}, w) << '\n';
      std::cout << std::__detail::__polylog_exp(-Tp{4}, w) << '\n';
      std::cout << std::__detail::__polylog_exp(Tp{2.6}, w) << '\n';
      std::cout << std::__detail::__polylog_exp(Tp{-2.6}, w) << '\n';
    }
    std::cout << '\n';

    std::cout << '\n' << std::__detail::__polylog_exp(Tp{2.6}, std::complex<Tp>(_S_pi, _S_pi)) << '\n';

    for(std::size_t k = 0; k < 10; ++k)
    {
      auto w = std::complex<Tp>(-_S_pi / 2 - _S_pi / 5, 0);
      std::cout << std::__detail::__polylog_exp(-Tp{4}, w) << '\n';
      std::cout << std::__detail::__polylog_exp_sum(-Tp{4}, w) << '\n';
    }
    std::cout << '\n' << std::flush;

    std::cout << '\n' << std::__detail::__polylog_exp_neg(Tp{-50.5}, std::complex<Tp>(Tp{1}, Tp{1})) << '\n';
    std::cout << '\n' << std::__detail::__polylog_exp_neg(Tp{-5}, std::complex<Tp>(Tp{1}, Tp{1})) << '\n';
    std::cout << '\n' << std::__detail::__polylog_exp_pos(Tp{2.3}, std::complex<Tp>(Tp{1}, Tp{1})) << '\n';
    //Don't trust Mathematica for small s
    std::cout << '\n' << std::__detail::__polylog_exp_asymp(Tp{60.4}, std::complex<Tp>(Tp{30}, Tp{0})) << '\n';

    // auto l = 2;
    // auto p = std::atan(l);
    // Tp alpha[] = {0.5, 1, 1.5, 4};
    // std::ofstream data("el20.txt");
    // for(std::size_t a = 0; a < sizeof(alpha) / sizeof(Tp); ++a)
    // {
    //   for(int s = -1 ; s <= 1; s += 2)
    //   {
    //     for(auto k = -_S_pi; k < _S_pi; k += Tp{0.002})
    //       data << k << ' ' << std::sqrt(Tp{1} + l * l) * real(std::exp(std::complex<Tp>(0, -s * p)) / (std::exp(std::complex<Tp>(0, k)) - std::exp(-alpha[a]))) << '\n';
    //     data << "&" << '\n';
    //   }
    // }

    const auto del01 = Tp{1} / Tp{100};
    const auto del05 = Tp{1} / Tp{20};

    std::ofstream test("test_polylog.dat");
    test.precision(std::cout.precision());
    for (auto s = Tp{2.5}; s < Tp{3.5}; s += del01)
      test << s << ' ' << std::setw(w) << std::real(std::__detail::__polylog(s, Tp{2})) - Tp{2} << '\n';
    test << '\n' << std::flush;

    std::cout << '\n' << std::__detail::__polylog(Tp{3.1}, Tp{2}) << '\n';
    std::cout << '\n' << std::__detail::__polylog_exp_pos(Tp{3.1}, std::complex<Tp>(std::log(Tp{2}))) << '\n';

    std::cout << "\nTest function 1 [PolyLog_Exp_pos(k,exp(i2pik)]:\n";
    for (std::size_t k = 3; k < 8; ++k)
      for (Tp x = 0; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_pos(k, std::polar(Tp{1}, _S_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest function 2 [PolyLog_Exp_pos(k,x)]:\n";
    for (std::size_t k = 3; k < 8; ++k)
      for (Tp x = 0; x < 6.28; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_pos(k, x)
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest function 3 [PolyLog_Exp_neg(s<0, exp(i2pik)]:\n";
    for (Tp k = -Tp{8}; k < Tp{0}; k += Tp{1} / Tp{13})
      for(Tp x = 0; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_neg(k, std::polar(Tp{1}, _S_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest function 4 + 5 [PolyLog_Exp_neg(k<0, exp(i2pik]:\n";
    for (int k = -40; k < 0; ++k)
      for (Tp x = 0; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_neg(k, std::polar(Tp{1}, _S_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest series 6 [PolyLog_Exp_pos(s, exp(i2pix)]:\n";
    for (Tp k = Tp{1} / Tp{7}; k < Tp{13}; k += Tp{1} / Tp{11})
      for (Tp x = Tp{0}; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_pos(k, std::polar(Tp{1}, _S_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest series 7 [PolyLog_Exp_asym(k, 100 exp(i2pix))]:\n";
    for (Tp k = Tp{-13}; k < Tp{13}; k += Tp{1} / Tp{11})
      for (Tp x = Tp{0}; x < Tp{1}; x += del01)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_asymp(k, Tp{100} * std::polar(Tp{1}, _S_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest series 8 [PolyLog_Exp_negative_real_part(k, x)]:\n";
    for (Tp k = -Tp{13}; k < Tp{13}; k += Tp{1} / Tp{11})
      for (Tp x = -Tp{7} / Tp{10} * _S_pi; x > -_S_2pi; x -= del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << std::__detail::__polylog_exp_sum(k, x)
		  << '\n';
    std::cout << '\n' << std::flush;
  }

int
main()
{
  test_polylog_0(1.0);

  test_polylog_1(1.0);

  test_polylog_dilog(1.0);

  test_polylog_cephes(1.0);

  test_polylog_neg_int(1.0);

  TestPolyLog<double>();

  TestPolyLog<float>();

  TestPolyLog<long double>();

  // This works but it takes forever.
  TestPolyLog<__float128>();

  return 0;
}
