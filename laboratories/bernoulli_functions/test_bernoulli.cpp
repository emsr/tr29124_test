/**
 *
 */

#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <limits>
#include <vector>

#include <emsr/float128_math.h>
#include <emsr/float128_io.h>
#include <emsr/math_constants.h>
#include <emsr/numeric_limits.h>
#include <emsr/special_functions.h>

#include <wrap_burkhardt.h>


  //template<typename RealTp, typename IntTp,
//	   IntTp Num = 1, IntTp Den = 1>
//    constexpr RealTp
//    frac = RealTp(Num) / RealTp(Den);

  template<typename RealTp,
	   unsigned long long Num = 1,
           unsigned long long Den = 1>
    constexpr RealTp
    frac = RealTp(Num) / RealTp(Den);

  /**
   * According to Modern Computer Arithmetic a stable recursion
   * for the Bernoulli numbers is
   * @f[
   *    \sum_{k=0}^{n} \frac{C_k}{(2n+1-2k)!4^{n-k}} = \frac{1}{(2n)! 4^n}
   * @f]
   * where @f$ C_k = B_{2k} / (2k)! @f$ is the scaled Bernoulli number.
   */
  template<typename Tp>
    Tp
    bernoulli_scaled_recur(unsigned int n)
    {
      return Tp{0};
    }


  /**
   * Return the Bernoulli number from lookup or by series expansion.
   */
  template<typename Tp>
    Tp
    bernoulli_series(unsigned int n)
    {
      static constexpr std::size_t s_len = 24;
      static constexpr Tp
      s_num[s_len]
      {
	 frac<Tp>,
	-frac<Tp, 1ull, 2ull>,
	 frac<Tp, 1ull, 6ull>,            Tp(0ull),
	-frac<Tp, 1ull, 30ull>,           Tp(0ull),
	 frac<Tp, 1ull, 42ull>,           Tp(0ull),
	-frac<Tp, 1ull, 30ull>,           Tp(0ull),
	 frac<Tp, 5ull, 66ull>,           Tp(0ull),
	-frac<Tp, 691ull, 2730ull>,       Tp(0ull),
	 frac<Tp, 7ull, 6ull>,            Tp(0ull),
	-frac<Tp, 3617ull, 510ull>,       Tp(0ull),
	 frac<Tp, 43867ull, 798ull>,      Tp(0ull),
	-frac<Tp, 174611ull, 330ull>,     Tp(0ull),
	 frac<Tp, 854513ull, 138ull>,     Tp(0ull)
      };

      if (n == 0)
	return frac<Tp>;
      else if (n == 1)
	return -frac<Tp, 1, 2>;
      else if (n % 2 == 1) // Take care of the rest of the odd ones.
	return Tp(0);
      else if (n < s_len) // return small evens that are painful for the series.
	return s_num[n];
      else
	{
	  Tp fact = Tp(1);
	  if ((n / 2) % 2 == 0)
	    fact *= -1;
	  for (auto k = 1u; k <= n; ++k)
	    fact *= k / emsr::pi_v<Tp>;
	  fact *= Tp(2);

	  Tp sum = Tp(0);
	  for (auto i = 1u; i < 1000; ++i)
	    {
	      Tp term = std::pow(Tp(i), -Tp(n));
	      if (term < emsr::epsilon<Tp>())
        	break;
	      sum += term;
	    }

	  return fact * sum;
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
  template<typename Tp>
    Tp
    bernoulli(unsigned int n, Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	{
	  auto B_n = emsr::detail::bernoulli<Tp>(0);
	  auto binomial = Tp{1};
	  for (auto k = 1u; k <= n; ++k)
	    {
	      binomial *= Tp(n + 1 - k) / Tp(k);
	      B_n = x * B_n + binomial
		   * emsr::detail::bernoulli<Tp>(k);
	    }
	  return B_n;
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
  template<typename Tp>
    Tp
    bernoulli_period(unsigned int n, Tp x)
    {
      int p = int(x);
      return bernoulli(n, x - Tp(p));
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
  template<typename Tp>
    Tp
    bernoulli(unsigned int n, Tp a, Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
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
  template<typename Tp>
    Tp
    euler_series(unsigned int n)
    {
      using emsr::detail::binomial;
      static constexpr std::size_t s_len = 22;
      static constexpr Tp
      s_num[s_len]
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

      if (n == 0)
	return Tp{1};
      else if (n & 1)
	return Tp{0};
      else if (n == 2)
        return Tp{-1};
      else if (n < s_len)
	return s_num[n];
      else
	{
	  std::vector<Tp> En(n + 1);
	  En[0] = Tp{1};
	  En[1] = Tp{0};
	  En[2] = Tp{-1};

	  for (auto i = 3u; i <= n; ++i)
	    {
	      En[i] = 0;

	      if (i % 2 == 0)
		{
		  for (auto j = 2u; j <= i; j += 2u)
		    En[i] -= binomial<Tp>(i, j)
			      * En[i - j];
		}
	    }
	  return En[n];
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
   * @param n the order n of the Euler number.
   * @return  The Euler number of order n.
   */
  template<typename Tp>
    inline Tp
    euler_asymp(unsigned int n)
    {
      if (n & 1)
	return Tp{0};
      else
	{
	  const auto s_e = emsr::e_v<Tp>;
	  const auto s_pi = emsr::pi_v<Tp>;
	  const auto n2 = Tp(n * n);
	  return emsr::parity<Tp>(n / 2) * Tp{8}
		* std::sqrt(Tp(n) / s_pi)
		* std::pow(Tp(4 * n) / (s_pi * s_e)
			 * (n2 + Tp{9} / Tp{480})
			 / (n2 - Tp{1} / Tp{480}), Tp(2 * n));
	}
    }

  /**
   * @brief This returns Euler number @f$ E_n @f$.
   *
   * @param n the order n of the Euler number.
   * @return  The Euler number of order n.
   */
  template<typename Tp>
    inline Tp
    euler(unsigned int n)
    { return euler_series<Tp>(n); }

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
  template<typename Tp>
    Tp
    euler(unsigned int n, Tp x)
    {
      if (std::isnan(x))
	return std::numeric_limits<Tp>::quiet_NaN();
      else
	{
	  auto bx1 = bernoulli(n + 1, x );
	  auto bx2 = bernoulli(n + 1, Tp{0.5L} * x );

	  auto E_n = Tp{2} * (bx1 - bx2 * std::pow(Tp{2}, Tp(n + 1)))
		    / Tp(n + 1);

	  return E_n;
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
  template<typename Tp>
    Tp
    euler_period(unsigned int n, Tp x)
    {
      int p = int(x);
      return ((p & 1) ? -1 : +1) * euler(n, x - Tp(p));
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
  template<typename Tp>
    Tp
    stirling_2_series(unsigned int n, unsigned int m)
    {
      using emsr::detail::log_factorial;
      using emsr::detail::factorial;
      if (m > emsr::detail::s_num_factorials<Tp>)
	{
	  auto S2 = Tp{0};
	  for (auto k = 0u; k <= m; ++k)
	    {
	      auto lf1 = log_factorial<Tp>(k);
	      auto lf2 = log_factorial<Tp>(m - k);
	      S2 += (((m - k) & 1) ? Tp{-1} : Tp{1})
		   * std::exp(n * std::log(k) - lf1 - lf2);
	    }
	  return S2;
	}
      else
	{
	  auto S2 = Tp{0};
	  for (auto k = 0u; k <= m; ++k)
	    {
	      S2 += (((m - k) & 1) ? Tp{-1} : Tp{1})
		   * std::pow(k, n)
		   / factorial<Tp>(k)
		   / factorial<Tp>(m - k);
	    }
	  // @todo Only round if the sum is less than
	  // the maximum representable integer.
	  // Find or make a tool for this.
	  return std::nearbyint(S2);
	}
      // Why this warning?
      return Tp{0};
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
  template<typename Tp>
    Tp
    stirling_2_recur(unsigned int n, unsigned int m)
    {
      if (n == 0)
	return Tp(m == 0);
      else if (m == 0)
	return Tp(n == 0);
      else
	{
	  std::vector<Tp> sigold(m + 1), signew(m + 1);
	  sigold[1] = Tp{1};
	  if (n == 1)
	    return sigold[m];
	  for (auto in = 1u; in <= n; ++in)
	    {
	      signew[1] = sigold[1];
	      for (auto im = 2u; im <= m; ++im)
		signew[im] = im * sigold[im] + sigold[im - 1];
	      std::swap(sigold, signew);
	    }
	  return signew[m];
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
  template<typename Tp>
    Tp
    stirling_2(unsigned int n, unsigned int m)
    {
      if (m > n)
	return Tp{0};
      else if (m == n)
	return Tp{1};
      else if (m == 0 && n >= 1)
	return Tp{0};
      else
	return stirling_2_recur<Tp>(n, m);
    }

  /**
   * Return the Stirling number of the second kind.
   *
   * @todo Find asymptotic expressions for the Stirling numbers.
   */
  template<typename Tp>
    Tp
    log_stirling_2(unsigned int n, unsigned int m)
    {
      if (m > n)
	return -std::numeric_limits<Tp>::infinity();
      else if (m == n)
	return Tp{0};
      else if (m == 0 && n >= 1)
	return -std::numeric_limits<Tp>::infinity();
      else
	return std::log(stirling_2<Tp>(n, m));
    }

  /**
   * Return the Stirling number of the first kind by series expansion.
   * N.B. This seems to be a total disaster.
   * Maybe accumulate the positive and negative terms separately
   * and add them att the end.
   */
  template<typename Tp>
    Tp
    stirling_1_series(unsigned int n, unsigned int m)
    {
      using emsr::detail::log_binomial;
      using emsr::detail::log_binomial_sign;
      using emsr::detail::binomial;
      if (2 * n - m > emsr::detail::s_num_factorials<Tp> / 2)
	{
	  auto S1 = Tp{0};
	  for (auto k = 0u; k <= n - m; ++k)
	    {
	      const auto nmpk = n - m + k;
	      const auto nmmk = n - m - k;
	      const auto lbc1 = log_binomial<Tp>(n - 1 + k, nmpk);
	      const auto slbc1 = log_binomial_sign<Tp>(n - 1 + k, nmpk);
	      const auto lbc2 = log_binomial<Tp>(2 * n - m, nmmk);
	      const auto slbc2 = log_binomial_sign<Tp>(2 * n - m, nmmk);
	      S1 += emsr::parity<Tp>(k) * slbc1 * slbc2
		   * std::exp(lbc1 + lbc2 + log_stirling_2<Tp>(nmpk, k));
	    }
	  return S1;
	}
      else
	{
	  auto S1 = Tp{0};
	  for (auto k = 0u; k <= n - m; ++k)
	    {
	      const auto nmpk = n - m + k;
	      const auto nmmk = n - m - k;
	      S1 += emsr::parity<Tp>(k)
		   * binomial<Tp>(n - 1 + k, nmpk)
		   * binomial<Tp>(2 * n - m, nmmk)
		   * stirling_2<Tp>(nmpk, k);
	    }
	  // @todo Only round if the sum is less than
	  // the maximum representable integer.
	  // Find or make a tool for this.
	  return std::nearbyint(S1);
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
  template<typename Tp>
    Tp
    stirling_1_recur(unsigned int n, unsigned int m)
    {
      if (n == 0)
	return Tp(m == 0);
      else if (m == 0)
	return Tp(n == 0);
      else
	{
	  std::vector<Tp> Sold(m + 1), Snew(m + 1);
	  Sold[1] = Tp{1};
	  if (n == 1)
	    return Sold[m];
	  for (auto in = 1u; in <= n; ++in)
	    {
	      for (auto im = 1u; im <= m; ++im)
		Snew[im] = Sold[im - 1] - in * Sold[im];
	      std::swap(Sold, Snew);
	    }
	  return Snew[m];
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
  template<typename Tp>
    Tp
    stirling_1(unsigned int n, unsigned int m)
    {
      if (m > n)
	return Tp{0};
      else if (m == n)
	return Tp{1};
      else if (m == 0 && n >= 1)
	return Tp{0};
      else
        return stirling_1_recur<Tp>(n, m);
    }

  /**
   * Return the logarithm of the absolute value of Stirling number
   * of the first kind.
   */
  template<typename Tp>
    Tp
    log_stirling_1(unsigned int n, unsigned int m)
    {
      if (m > n)
	return -std::numeric_limits<Tp>::infinity();
      else if (m == n)
	return Tp{0};
      else if (m == 0 && n >= 1)
	return -std::numeric_limits<Tp>::infinity();
      else
	return std::log(std::abs(stirling_1<Tp>(n, m)));
    }

  /**
   * Return the sign of the exponent of the logarithm of the Stirling number
   * of the first kind.
   */
  template<typename Tp>
    inline Tp
    log_stirling_1_sign(unsigned int n, unsigned int m)
    { return ((n + m) & 1) ? Tp{-1} : Tp{+1}; }

  /**
   * Return the Eulerian number of the first kind by recursion.
   * The recursion is
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    Tp
    eulerian_1_recur(unsigned int n, unsigned int m)
    {
      if (m == 0)
	return Tp{1};
      else if (m >= n)
	return Tp{0};
      else if (m == n - 1)
	return Tp{1};
      else if (n - m - 1 < m) // Symmetry.
	return eulerian_1_recur<Tp>(n, n - m - 1);
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<Tp> Aold(m + 1), Anew(m + 1);
	  Aold[0] = Tp{1};
	  Anew[0] = Tp{1};
	  Anew[1] = Tp{1};
	  for (auto in = 3u; in <= n; ++in)
	    {
	      std::swap(Aold, Anew);
	      for (auto im = 1u; im <= m; ++im)
		Anew[im] = (in - im) * Aold[im - 1]
			    + (im + 1) * Aold[im];
	    }
	  return Anew[m];
	}
    }

  /**
   * Return the Eulerian number of the first kind.
   * The Eulerian numbers are defined by recursion:
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    inline Tp
    eulerian_1(unsigned int n, unsigned int m)
    { return eulerian_1_recur<Tp>(n, m); }

  /**
   * Return a vector Eulerian numbers of the first kind by recursion.
   * The recursion is
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    std::vector<Tp>
    eulerian_1_recur(unsigned int n)
    {
      if (n == 0)
	return std::vector<Tp>(1, Tp{1});
      //else if (m == n - 1)
	//return Tp{1};
      //else if (n - m - 1 < m) // Symmetry.
	//return eulerian_1_recur<Tp>(n, n - m - 1);
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<Tp> Aold(n + 1), Anew(n + 1);
	  Aold[0] = Anew[0] = Tp{1};
	  Anew[1] = Tp{1};
	  for (auto in = 3u; in <= n; ++in)
	    {
	      std::swap(Aold, Anew);
	      for (auto im = 1u; im <= n; ++im)
		Anew[im] = (in - im) * Aold[im - 1]
			    + (im + 1) * Aold[im];
	    }
	  return Anew;
	}
    }

  /**
   * Return a vector Eulerian numbers of the first kind.
   * The Eulerian numbers are defined by recursion:
   * @f[
   *   A(n,m) = (n-m)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    inline std::vector<Tp>
    eulerian_1(unsigned int n)
    { return eulerian_1_recur<Tp>(n); }

  /**
   * Return the Eulerian number of the second kind by recursion.
   * The recursion is:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    Tp
    eulerian_2_recur(unsigned int n, unsigned int m)
    {
      if (m == 0)
	return Tp{1};
      else if (m >= n)
	return Tp{0};
      else if (n == 0)
	return Tp{1};
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<Tp> Aold(m + 1), Anew(m + 1);
	  Aold[0] = Tp{1};
	  Anew[0] = Tp{1};
	  Anew[1] = Tp{2};
	  for (auto in = 3u; in <= n; ++in)
	    {
	      std::swap(Aold, Anew);
	      for (auto im = 1u; im <= m; ++im)
		Anew[im] = (2 * in - im - 1) * Aold[im - 1]
			    + (im + 1) * Aold[im];
	    }
	  return Anew[m];
	}
    }

  /**
   * Return the Eulerian number of the second kind.
   * The Eulerian numbers of the second kind are defined by recursion:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    inline Tp
    eulerian_2(unsigned int n, unsigned int m)
    { return eulerian_2_recur<Tp>(n, m); }

  /**
   * Return the Eulerian number of the second kind by recursion.
   * The recursion is:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    std::vector<Tp>
    eulerian_2_recur(unsigned int n)
    {
      if (n == 0)
	return std::vector<Tp>(1, Tp{1});
      //else if (m >= n)
	//return Tp{0};
      //else if (n == 0)
	//return Tp{1};
      else
	{
	  // Start recursion with n == 2 (already returned above).
	  std::vector<Tp> Aold(n + 1), Anew(n + 1);
	  Aold[0] = Anew[0] = Tp{1};
	  Anew[1] = Tp{2};
	  for (auto in = 3u; in <= n; ++in)
	    {
	      std::swap(Aold, Anew);
	      for (auto im = 1u; im <= n; ++im)
		Anew[im] = (2 * in - im - 1) * Aold[im - 1]
			    + (im + 1) * Aold[im];
	    }
	  return Anew;
	}
    }

  /**
   * Return the Eulerian number of the second kind.
   * The Eulerian numbers of the second kind are defined by recursion:
   * @f[
   *   A(n,m) = (2n-m-1)A(n-1,m-1) + (m+1)A(n-1,m) \mbox{ for } n > 0
   * @f]
   */
  template<typename Tp>
    inline std::vector<Tp>
    eulerian_2(unsigned int n)
    { return eulerian_2_recur<Tp>(n); }

  /**
   * Return the Lah number by downward recurrence.
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename Tp>
    Tp
    lah_recur(unsigned int n, unsigned int k)
    {
      if (k > n)
	return Tp{0};
      else if (n == 0)
	return (k == 0 ? Tp{1} : Tp{0});
      else
	{
	  Tp Lnn = 1;
	  for (unsigned int i = 1u; i <= n - k; ++i)
	    Lnn *= Tp(n - i + 1) * Tp(n - i) / Tp(i);
	  return Lnn;
	}
    }

  /**
   * Return a vector of Lah numbers by downward recurrence.
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename Tp>
    std::vector<Tp>
    lah_recur(unsigned int n)
    {
      if (n == 0)
	return std::vector<Tp>(1, Tp{1});
      else
	{
	  std::vector<Tp> L(n + 1);
	  Tp Lnn = 1;
	  L[n] = Lnn;
	  for (unsigned int i = 1u; i <= n; ++i)
	    {
	      Lnn *= Tp(n - i + 1) * Tp(n - i) / Tp(i);
	      L[n - i] = Lnn;
	    }
	  return L;
	}
    }

  /**
   * Return a vector of Lah numbers.
   * @f[
   *   L(n,k-1) = \frac{k(k-1)}{n-k+1}L(n,k);  L(n,n) = 1
   * @f]
   */
  template<typename Tp>
    inline std::vector<Tp>
    lah(unsigned int n)
    { return lah_recur<Tp>(n); }

  /**
   * Return a vector of the Bell numbers by summation.
   * @f[
   *   B(n) = \sum_{k=0}^{n}S_n^{(k)}
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename Tp>
    std::vector<Tp>
    bell_series(unsigned int n)
    {
      std::vector<Tp> bell(n + 1);
      bell[0] = Tp{1};

      /// @todo Test for blowup in Bell number summation.
      for (unsigned int i = 1; i <= n; ++i)
	for (unsigned int j = 1; j <= i; ++j)
	  bell[i] += bell[i - j]
		       * emsr::binomial<Tp>(i - 1, j - 1);

      return bell;
    }

  /**
   * Return a vector of the Bell numbers.
   */
  template<typename Tp>
    inline std::vector<Tp>
    bell(unsigned int n)
    { return bell_series<Tp>(n); }

  /**
   * Evaluate the Bell polynomial
   * @f[
   *   B(n) = \sum_{k=0}^{n}S_n^{(k)}x^k
   * @f]
   * where @f$ S_n^{(k)} @f$ are the Stirling numbers of the second kind.
   */
  template<typename Tp, typename Up>
    inline Up
    bell(unsigned int n, Up x)
    {
      const auto Sn = stirling_2<Tp>(n);
      auto bell = Sn[n];
      for (unsigned int i = 1; i < n; ++i)
	bell = Sn[n - i] + x * bell;
      return bell;
    }

template<typename Tp>
  void
  test_bernoulli(Tp proto = Tp{})
  {
    std::cout.precision(emsr::digits10(proto));
    std::cout << std::showpoint << std::scientific;
    auto width = 8 + std::cout.precision();

    std::cout << "\n\n Bell numbers\n";
    for (auto n = 1u; n <= 20; ++n)
      {
	auto bell = bell_series<long long>(n);
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
	const auto lahv = lah<Tp>(n);
	std::cout << '\n';
	for (auto k = 0u; k <= n; ++k)
	  std::cout << ' ' << std::setw(4) << n
	       	    << ' ' << std::setw(4) << k
		    << ' ' << std::setw(width) << lah_recur<Tp>(n, k)
		    << ' ' << std::setw(width) << lahv[k]
		    << '\n';
      }

    std::cout << "\n\n Bernoulli numbers\n";
    for (auto n = 0u; n <= 200; ++n)
      std::cout << ' ' << std::setw(4) << n
		<< ' ' << std::setw(width) << bernoulli_series<Tp>(n)
		<< '\n';

    std::cout << "\n\n Bernoulli polynomials\n";
    for (auto n = 0u; n <= 50; ++n)
      {
	std::cout << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = Tp{1} / Tp{10};
	for (auto i = 0u; i <= 50; ++i)
	  {
	    auto x = del * i;
	    auto B_n = bernoulli(n, x);
	    auto tildeB_n = bernoulli_period(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << B_n
		      << ' ' << std::setw(width) << tildeB_n
		      << '\n';
	  }
      }

    std::cout << "\n\n Euler numbers\n";
    for (auto n = 0u; n <= 200; ++n)
      std::cout << ' ' << std::setw(4) << n
		<< ' ' << std::setw(width) << euler<Tp>(n)
		<< '\n';

    std::cout << "\n\n Euler polynomials\n";
    for (auto n = 0u; n <= 50; ++n)
      {
	std::cout << '\n' << ' ' << std::setw(4) << n << '\n';
	const auto del = Tp{1} / Tp{10};
	for (auto i = 0u; i <= 50; ++i)
	  {
	    auto x = del * i;
	    auto E_n = euler(n, x);
	    auto tildeE_n = euler_period(n, x);
	    std::cout << ' ' << std::setw(width) << x
		      << ' ' << std::setw(width) << E_n
		      << ' ' << std::setw(width) << tildeE_n
		      << '\n';
	  }
      }

    std::cout << "\n\n Stirling numbers of the second kind";
    for (auto n = 0u; n <= 100; ++n)
      {
        const auto s2v = emsr::detail::stirling_2<Tp>(n);
	std::cout << '\n';
	for (auto m = 0u; m <= n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		    << ' ' << std::setw(width) << stirling_2_series<Tp>(n, m)
		    << ' ' << std::setw(width) << stirling_2_recur<Tp>(n, m)
		    << ' ' << std::setw(width) << s2v[m]
		    << '\n';
      }

    std::cout << "\n\n Stirling numbers of the first kind";
    for (auto n = 0u; n <= 100; ++n)
      {
        const auto s1v = emsr::detail::stirling_1<Tp>(n);
	std::cout << '\n';
	for (auto m = 0u; m <= n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		    << ' ' << std::setw(width) << stirling_1_series<Tp>(n, m)
		    << ' ' << std::setw(width) << stirling_1_recur<Tp>(n, m)
		    << ' ' << std::setw(width) << s1v[m]
		    << '\n';
      }

    std::cout << "\n\n Eulerian numbers of the first kind";
    for (auto n = 1u; n <= 10; ++n)
      {
        const auto e1v = eulerian_1<Tp>(n);
	std::cout << '\n';
	for (auto m = 0u; m < n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		  //  << ' ' << std::setw(width) << eulerian_1_series<Tp>(n, m)
		    << ' ' << std::setw(width) << eulerian_1_recur<Tp>(n, m)
		    << ' ' << std::setw(width) << burkhardt::eulerian_1(n, m)
		    << ' ' << std::setw(width) << e1v[m]
		    << '\n';
      }

    std::cout << "\n\n Eulerian numbers of the second kind";
    for (auto n = 1u; n <= 10; ++n)
      {
        const auto e2v = eulerian_2<Tp>(n);
	std::cout << '\n';
	for (auto m = 0u; m < n; ++m)
	  std::cout << ' ' << std::setw(4) << n
		    << ' ' << std::setw(4) << m
		  //  << ' ' << std::setw(width) << eulerian_2_series<Tp>(n, m)
		    << ' ' << std::setw(width) << eulerian_2_recur<Tp>(n, m)
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

