/**
 *
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <complex>

#include <emsr/float128_math.h>
#include <emsr/float128_io.h>
#include <emsr/fp_type_util.h>
#include <emsr/sf_stirling.h>
#include <emsr/sf_euler.h>
#include <emsr/sf_gamma.h>
#include <emsr/sf_zeta.h>
#include <emsr/sf_polylog.h>

#include <wrap_cephes.h>

  /**
   * This class manages the termination of series.
   * Termination conditions involve both a maximum iteration count
   * and a relative precision.
   */
  template<typename Tp>
    class _Terminator
    {
    private:

      using Real = emsr::num_traits_t<Tp>;
      const std::size_t _M_max_iter;
      std::size_t _M_curr_iter;
      Real _M_toler;

    public:

      _Terminator(std::size_t max_iter, Real mul = Real{1})
      : _M_max_iter(max_iter), _M_curr_iter{0},
	_M_toler(std::abs(mul) * std::numeric_limits<Real>::epsilon())
      { }

      /// Return the current number of terms summed.
      std::size_t
      num_terms() const
      { return this->_M_curr_iter; }

      /// Detect if the sum should terminate either because the incoming term
      /// is small enough or the maximum number of terms has been reached.
      bool
      operator()(Tp term, Tp sum)
      {
	if (this->_M_curr_iter >= this->_M_max_iter
	    || ++this->_M_curr_iter == this->_M_max_iter)
	  return true;
	else if (std::abs(term) < this->_M_toler * std::abs(sum))
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
  template<typename Tp>
    class _AsympTerminator
    {
    private:

      using Real = emsr::num_traits_t<Tp>;
      const std::size_t _M_max_iter;
      std::size_t _M_curr_iter;
      Real _M_toler;
      Real _M_prev_term = std::numeric_limits<Real>::max();
      bool _M_stop_asymp = false;

    public:

      _AsympTerminator(std::size_t max_iter, Real mul = Real{1})
      : _M_max_iter(max_iter), _M_curr_iter{0},
	_M_toler(std::abs(mul) * std::numeric_limits<Real>::epsilon())
      { }

      /// Filter a term before applying it to the sum.
      Tp
      operator<<(Tp term)
      {
	if (std::abs(term) > this->_M_prev_term)
	  {
	    this->_M_stop_asymp = true;
	    return Tp{0};
	  }
	else
	  return term;
      }

      /// Return the current number of terms summed.
      std::size_t
      num_terms() const
      { return this->_M_curr_iter; }

      /// Detect if the sum should terminate either because the incoming term
      /// is small enough or because the terms are starting to grow or
      //  the maximum number of terms has been reached.
      bool
      operator()(Tp term, Tp sum)
      {
	if (this->_M_stop_asymp)
	  return true;
	else
	  {
	    const auto aterm = std::abs(term);
	    this->_M_stop_asymp = (aterm > this->_M_prev_term);
	    this->_M_prev_term = aterm;
	    if (this->_M_curr_iter >= this->_M_max_iter
	    || ++this->_M_curr_iter == this->_M_max_iter)
	      return true;
	    else if (aterm < this->_M_toler * std::abs(sum))
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
  template<typename Tp>
    Tp
    polylog_nonpos_int_1(int n, Tp x)
    {
      if (x == Tp{0})
	return Tp{0};
      else if (x == Tp{1})
	return std::numeric_limits<Tp>::infinity();
      else
	{
	  const auto arg = x / (Tp{1} - x);
	  auto fact = arg;
	  auto term = arg * emsr::stirling_2<Tp>(1 - n, 1);
	  auto sum = term;
	  for (int k = 1; k <= -n; ++k)
	    {
	      fact *= k * arg;
	      term = fact * emsr::stirling_2<Tp>(1 - n, 1 + k);
	      sum += term;
	    }
	  return sum;
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
  template<typename Tp>
    Tp
    polylog_nonpos_int_2(int n, Tp x)
    {
      if (x == Tp{0})
	return Tp{0};
      else if (x == Tp{1})
	return std::numeric_limits<Tp>::infinity();
      else
	{
	  const auto arg = Tp{-1} / (Tp{1} - x);
	  auto fact = arg;
	  auto term = fact * emsr::stirling_2<Tp>(1 - n, 1);
	  auto sum = term;
	  for (int k = 1; k <= -n; ++k)
	    {
	      fact *= k * arg;
	      term = fact * emsr::stirling_2<Tp>(1 - n, 1 + k);
	      sum += term;
	    }
	  sum *= emsr::parity<Tp>(1 - n);
	  return sum;
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
  template<typename Tp>
    Tp
    polylog_nonpos_int_3(int n, Tp x)
    {
      if (x == Tp{0})
	return Tp{0};
      else if (x == Tp{1})
	return std::numeric_limits<Tp>::infinity();
      else
	{
	  auto fact = Tp{1} / x;
	  auto term = fact * emsr::eulerian_1<Tp>(-n, 0);
	  auto sum = term;
	  for (int k = 1; k < -n; ++k)
	    {
	      fact /= x;
	      term = fact * emsr::eulerian_1<Tp>(-n, k);
	      sum += term;
	    }
	  sum *= std::pow(x / (Tp{1} - x), Tp(1 - n));
	  return sum;
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
   * @param n the negative integer index @f$ n = -p @f$.
   * @param w the argument w.
   * @return the value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg(int n, std::complex<Tp> w)
    {
      const auto s_inf = std::numeric_limits<Tp>::infinity();
      if (emsr::fp_is_zero(w))
	return std::complex<Tp>{0};
      else if (emsr::fp_is_equal(w, Tp{1}))
	return std::complex<Tp>{s_inf, Tp{0}};
      else
	{
	  const int p = -n;
	  const int pp = 1 + p;
	  const int q = (p & 1) ? 0 : 1;
	  const auto w2 = w * w;
	  auto wp = (p & 1) ? std::complex<Tp>{1} : w;
	  unsigned int __2k = q;
	  auto gam = emsr::detail::factorial<Tp>(p + __2k);
	  const auto pfact = emsr::detail::factorial<Tp>(p);
	  auto res = pfact * std::pow(-w, Tp(-pp));
	  auto sum = std::complex<Tp>{};
	  constexpr unsigned int maxit = 300;
	  _Terminator<std::complex<Tp>> done(maxit);
	  while (true)
	    {
	      const auto id = (p + __2k + 1) / 2;
	      if (id == emsr::detail::Num_Euler_Maclaurin_zeta)
		break;
	      auto term = gam * wp
		* Tp(emsr::detail::Euler_Maclaurin_zeta[id]);
	      sum += term;
	      if (done(term, sum))
		break;
	      gam *= Tp(p + __2k + 1) / Tp(__2k + 1)
		     * Tp(p + __2k + 2) / Tp(__2k + 2);
	      wp *= w2;
	      __2k += 2;
	    }
	  res -= sum;
	  return res;
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
   * @param p the integral index @f$ p = 4k @f$ or @f$ p = 4k + 2 @f$.
   * @param w The complex argument w
   * @return the value of the Polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg_even(unsigned int p, const std::complex<Tp>& w)
    {
      const auto s_2pi = emsr::tau_v<Tp>;
      const int pp = 1 + p;
      const int sigma = p % 4 == 0 ? +1 : -1;
      const auto lnp = emsr::detail::log_factorial<Tp>(p);
      auto res = std::exp(lnp - Tp(pp) * std::log(-w));
      auto wup = w / s_2pi;
      auto wq = -wup * wup;
      auto pref = Tp{2} * std::pow(s_2pi, Tp(-pp));
      // Subtract the expression A_p(w)
      pref *= sigma;
      res -= pref
	     * std::exp(lnp - Tp{0.5L} * pp * std::log(Tp{1} - wq))
	     * std::cos(Tp(pp) * std::atan(Tp{1} / wup));
      unsigned int __2k = 0;
      auto gam = emsr::detail::factorial<Tp>(1 + p);
      std::complex<Tp> sum;
      constexpr unsigned int maxit = 400;
      _Terminator<std::complex<Tp>> done(maxit);
      while (true)
	{
	  auto term = gam * wup
		      * emsr::detail::riemann_zeta_m_1<Tp>(__2k + 2 + p);
	  sum += term;
	  if (done(term, sum))
	    break;
	  gam *= Tp(__2k + 2 + p) / Tp(__2k + 2)
		 * Tp(__2k + 3 + p) / Tp(__2k + 3);
	  __2k += 2;
	  wup *= wq;
	}
      res -= pref * sum;
      return res;
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
   * @param p the integral index @f$ p = 4k + 1 @f$ or @f$ p = 4k + 3 @f$.
   * @param w The complex argument w.
   * @return The value of the Polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg_odd(unsigned int p, const std::complex<Tp>& w)
    {
      const auto s_2pi = emsr::tau_v<Tp>;
      const int pp = 1 + p;
      const int sigma = p % 4 == 1 ? +1 : -1;
      const auto lnp = emsr::detail::log_factorial<Tp>(p);
      auto res = std::exp(lnp - Tp(pp) * std::log(-w));
      auto wup = w / s_2pi;
      auto wq = -wup * wup;
      auto pref = Tp{2} * std::pow(s_2pi, Tp(-pp));
      pref *= sigma;
      // Add the expression A_p(w)
      res += pref
	     * std::exp(lnp - Tp{0.5L} * pp * std::log(Tp{1} - wq))
	     * std::cos(Tp(pp) * std::atan(Tp{1} / wup));
      auto gam = emsr::detail::factorial<Tp>(p);
      unsigned int __2k = 0;
      std::complex<Tp> sum;
      constexpr unsigned int maxit = 400;
      _Terminator<std::complex<Tp>> done(maxit);
      while (true)
	{
	  auto term = gam * wup
		      * emsr::detail::riemann_zeta_m_1<Tp>(__2k + 1 + p);
	  sum += term;
	  if (done(term, sum))
	    break;
	  gam *= Tp(__2k + 1 + p) / Tp(__2k + 1)
		 * Tp(__2k + 2 + p) / Tp(__2k + 2);
	  __2k += 2;
	  wup *= wq;
	}
      res -= pref * sum;
      return res;
  }

  /**
   * This function treats the cases of negative integer index s
   * and branches accordingly
   *
   * @param s the integer index s.
   * @param w The Argument w
   * @return The value of the Polylogarithm evaluated by a suitable function.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_neg_old(int s, std::complex<Tp> w)
    { // negative integer s
      const auto p = -s;
      switch (p % 2)
      {
      case 0:
	return polylog_exp_neg_even<Tp>(p, w);
      case 1:
	return polylog_exp_neg_odd<Tp>(p, w);
      default: // We shouldn't need this.
	return std::complex<Tp>{};
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
   * @param s the real index s.
   * @param w the large complex argument w.
   * @return the value of the polylogarithm.
   */
  template<typename Tp>
    std::complex<Tp>
    polylog_exp_asymp(Tp s, std::complex<Tp> w)
    {
      const auto s_pi = emsr::pi_v<Tp>;
      // wgamma = w^{s-1} / Gamma(s)
      auto wgamma = std::pow(w, s - Tp{1}) * emsr::detail::gamma_reciprocal(s);
      auto res = std::complex<Tp>(Tp{0}, -s_pi) * wgamma;
      // wgamma = w^s / Gamma(s+1)
      wgamma *= w / s;
      constexpr unsigned int maxiter = 100;
      _AsympTerminator<std::complex<Tp>> done(maxiter);
      // zeta(0) w^s / Gamma(s + 1)
      std::complex<Tp> oldterm = -Tp{0.5L} * wgamma;
      res += Tp{2} * oldterm;
      std::complex<Tp> term;
      auto wq = Tp{1} / (w * w);
      int k = 1;
      while (true)
	{
	  wgamma *= wq * (s + Tp(1 - 2 * k)) * (s + Tp(2 - 2 * k));
	  term = emsr::detail::riemann_zeta<Tp>(2 * k) * wgamma;
	  res += done << Tp{2} * term;
	  if (done(Tp{2} * term, res))
	    break;
	  oldterm = term;
	  ++k;
	}
      return res;
    }

/**
 * Compute the polylogarithm for negative integer order.
 */
template<typename Tp>
  void
  test_polylog_neg_int(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::scientific;

    //std::cout << "\nTest negative integer order\n";
    for (auto n : {-1, -2, -3, -4, -5})
      {
	std::cout << "\n\nNegative integer degree: n = " << n << '\n';
	const auto del = Tp{1} / Tp{20};
	for (int i = -200; i <= 20; ++i)
	  {
	    auto x = del * i;
	    auto Ls_rat1 = polylog_nonpos_int_1(n, x);
	    auto Ls_rat2 = polylog_nonpos_int_2(n, x);
	    auto Ls_rat3 = polylog_nonpos_int_3(n, x);
	    auto Ls_nint = x == Tp{0}
			 ? Tp{0}
			 : [n, x]() -> Tp
			   { auto w = std::log(std::complex<Tp>(x));
			     return std::real(emsr::detail::polylog_exp_neg(n, w)); }();
	    auto Ls_gnu = x == Tp{0}
			? Tp{0}
			: [n, x]() -> Tp
			  { auto w = std::log(std::complex<Tp>(x));
			    return std::real(polylog_exp_neg_old(n, w)); }();
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
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::scientific;

    std::cout << "\n\nTest against algebraic function for Li_0\n";
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i <= 10; ++i)
      {
	auto x = del * i;
	auto Li0 = emsr::fp_is_zero(Tp{1} - x)
		 ? std::numeric_limits<Tp>::quiet_NaN()
		 : x / (Tp{1} - x);
	auto Lis_gnu = emsr::polylog(Tp{0}, x);
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
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::scientific;

    std::cout << "\n\nTest against algebraic function for Li_1\n";
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i <= 10; ++i)
      {
	auto x = del * i;
	auto Li1 = emsr::fp_is_zero(Tp{1} - x)
		 ? std::numeric_limits<Tp>::quiet_NaN()
		 : -std::log(Tp{1} - x);
	auto Lis_gnu = emsr::polylog(Tp{1}, x);
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
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::scientific;

    std::cout << "\n\nTest against local dilog\n";
    const auto del = Tp{1} / Tp{10};
    for (int i = -200; i <= 10; ++i)
      {
	auto x = del * i;
	auto Ls_dilog = emsr::dilog(x);
	auto Ls_gnu = emsr::polylog(Tp(2), x);
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
    std::cout.precision(emsr::digits10(proto));
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
	    auto Ls_gnu = emsr::polylog(Tp(n), x);
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
    const auto s_pi = emsr::pi_v<Tp>;
    const auto s_2pi = emsr::tau_v<Tp>;

    std::cout.precision(emsr::digits10(proto) - 1);
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
    //       emsr::detail::riemann_zeta(x)
    // //      std::tr1::detail::riemann_zeta(x)
    //       ;//between 1 and 10 riemann_zeta_glob is called
    // //      <<'\n';
    //     }

    // Something that didn't work in the original implementation
    std::cout << emsr::riemann_zeta(std::complex<Tp>{5.1, 0.5}) << '\n';
    std::cout << emsr::hurwitz_zeta(Tp{5.1}, Tp{0.5}) << '\n';
    std::cout << emsr::hurwitz_zeta(Tp{5.1}, std::complex<Tp>{0.5}) << '\n';
    std::cout << emsr::detail::hurwitz_zeta_polylog(Tp{5.1}, std::complex<Tp>{0.5}) << '\n';
    std::cout << emsr::detail::polylog_exp(Tp{2.5}, std::complex<Tp>(Tp{15}, Tp{1})) << '\n';

    for(std::size_t k = 0; k < 32; ++k)
      {
	auto w = std::complex<Tp>(Tp{0}, s_2pi * k / Tp{32});
	std::cout << "=======  " << k << "  ==========" << '\n';
	std::cout << emsr::detail::polylog_exp(Tp{4}, w) << '\n';
	std::cout << emsr::detail::polylog_exp(Tp{-4}, w) << '\n';
	std::cout << emsr::detail::polylog_exp(Tp{2.6}, w) << '\n';
	std::cout << emsr::detail::polylog_exp(Tp{-2.6}, w) << '\n';
      }
    std::cout << '\n';

    std::cout << '\n' << emsr::detail::polylog_exp(Tp{2.6}, std::complex<Tp>(s_pi, s_pi)) << '\n';

    for(std::size_t k = 0; k < 10; ++k)
      {
	auto w = std::complex<Tp>(-s_pi / 2 - s_pi / 5, 0);
	std::cout << emsr::detail::polylog_exp(Tp{-4}, w) << '\n';
	std::cout << emsr::detail::polylog_exp_sum(Tp{-4}, w) << '\n';
      }
    std::cout << '\n' << std::flush;

    std::cout << '\n' << emsr::detail::polylog_exp_neg(Tp{-50.5}, std::complex<Tp>(Tp{1}, Tp{1})) << '\n';
    std::cout << '\n' << emsr::detail::polylog_exp_neg(Tp{-5}, std::complex<Tp>(Tp{1}, Tp{1})) << '\n';
    std::cout << '\n' << emsr::detail::polylog_exp_pos(Tp{2.3}, std::complex<Tp>(Tp{1}, Tp{1})) << '\n';
    //Don't trust Mathematica for small s
    std::cout << '\n' << emsr::detail::polylog_exp_asymp(Tp{60.4}, std::complex<Tp>(Tp{30}, Tp{0})) << '\n';

    auto l = Tp{2};
    auto p = std::atan(l);
    std::ofstream data("polylog_el20.txt");
    for(auto alpha : {Tp{0.5}, Tp{1}, Tp{1.5}, Tp{4}})
      {
	for(int s : {-1, 1})
	  {
	    for(auto k = -s_pi; k < s_pi; k += Tp{0.002})
	      {
		data << k << ' ' << std::sqrt(Tp{1} + l * l)
				* real(std::exp(std::complex<Tp>(0, -s * p))
				/ (std::exp(std::complex<Tp>(0, k)) - std::exp(-alpha))) << '\n';
		data << "&" << '\n';
	      }
	  }
      }

    const auto del01 = Tp{1} / Tp{100};
    const auto del05 = Tp{1} / Tp{20};

    std::ofstream test("test_polylog.dat");
    test.precision(std::cout.precision());
    for (auto s = Tp{2.5}; s < Tp{3.5}; s += del01)
      test << s << ' ' << std::setw(w) << std::real(emsr::detail::polylog(s, Tp{2})) - Tp{2} << '\n';
    test << '\n' << std::flush;

    std::cout << '\n' << emsr::detail::polylog(Tp{3.1}, Tp{2}) << '\n';
    std::cout << '\n' << emsr::detail::polylog_exp_pos(Tp{3.1}, std::complex<Tp>(std::log(Tp{2}))) << '\n';

    std::cout << "\nTest function 1 [PolyLog_Exp_pos(k,exp(i2pik)]:\n";
    for (std::size_t k = 3; k < 8; ++k)
      for (Tp x = 0; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << emsr::detail::polylog_exp_pos(k, std::polar(Tp{1}, s_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest function 2 [PolyLog_Exp_pos(k,x)]:\n";
    for (std::size_t k = 3; k < 8; ++k)
      for (Tp x = 0; x < 6.28; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << emsr::detail::polylog_exp_pos(k, x)
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest function 3 [PolyLog_Exp_neg(s<0, exp(i2pik)]:\n";
    for (Tp k = Tp{-8}; k < Tp{0}; k += Tp{1} / Tp{13})
      for(Tp x = 0; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << emsr::detail::polylog_exp_neg(k, std::polar(Tp{1}, s_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest function 4 + 5 [PolyLog_Exp_neg(k<0, exp(i2pik]:\n";
    for (int k = -40; k < 0; ++k)
      for (Tp x = 0; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << emsr::detail::polylog_exp_neg(k, std::polar(Tp{1}, s_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest series 6 [PolyLog_Exp_pos(s, exp(i2pix)]:\n";
    for (Tp k = Tp{1} / Tp{7}; k < Tp{13}; k += Tp{1} / Tp{11})
      for (Tp x = Tp{0}; x < Tp{1}; x += del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << emsr::detail::polylog_exp_pos(k, std::polar(Tp{1}, s_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest series 7 [PolyLog_Exp_asym(k, 100 exp(i2pix))]:\n";
    for (Tp k = Tp{-13}; k < Tp{13}; k += Tp{1} / Tp{11})
      for (Tp x = Tp{0}; x < Tp{1}; x += del01)
	std::cout << k
		  << ' ' << x
		  << ' ' << emsr::detail::polylog_exp_asymp(k, Tp{100} * std::polar(Tp{1}, s_2pi * x))
		  << '\n';
    std::cout << '\n' << std::flush;

    std::cout << "\nTest series 8 [PolyLog_Exp_negative_real_part(k, x)]:\n";
    for (Tp k = Tp{-13}; k < Tp{13}; k += Tp{1} / Tp{11})
      for (Tp x = Tp{-7} / Tp{10} * s_pi; x > -s_2pi; x -= del05)
	std::cout << k
		  << ' ' << x
		  << ' ' << emsr::detail::polylog_exp_sum(k, x)
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
  //TestPolyLog<__float128>();

  return 0;
}
