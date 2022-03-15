#ifndef CHEBYSHEV_TCC
#define CHEBYSHEV_TCC 1

#include <utility>
#include <iterator>
#include <algorithm> // For find_if.
#include <numbers>

#include <emsr/math_constants.h>

/**
 * @file chebyshev.tcc Implementation for C++ Chebyshev methods.
 */

namespace emsr
{

  /**
   * @brief  Construct a Chebyshev fit of a function.
   *         Given a function @c func, the lower limit @c a, the upper limit @c b,
   *         and the number of points @c n this routine computes the coefficients
   *         of a Chebyshev polynomial expansion such that
   * @f[
   *   f(x) = k\sum_{0}^{n-1} c_k T_k(x) - c_0/2.
   * @f]
   * Use this routine with moderately large n of 30 or 50.
   * Then truncate the series to a smaller number of terms to satisfy
   * the accuracy requirements.
   */
  template<typename Tp>
    Chebyshev<Tp>::Chebyshev(Tp a, Tp b, unsigned n, Tp func(Tp))
    : m_lower(a),
      m_upper(b),
      m_coef(n)
    {
      constexpr Tp _S_pi = emsr::pi_v<Tp>;

      if (this->upper() == this->lower())
	throw std::domain_error("_Chebyshev: input domain has zero length");

      std::vector<Tp> f(n);

      const auto bma = (m_upper - m_lower) / Tp{2};
      const auto bpa = (m_upper + m_lower) / Tp{2};
      for (unsigned int k = 0; k < n; ++k)
	{
	  Tp y = std::cos(_S_pi * (k + Tp(0.5L)) / n);
	  f[k] = func(bpa + y * bma);
	}

      const auto fac = Tp{2} / n;
      for (unsigned int j = 0; j < n; ++j)
	{
	  auto sum = Tp{};
	  for (unsigned int k = 0; k < n; ++k)
	    sum += f[k] * std::cos(_S_pi * j * (k + Tp(0.5L)) / n);
	  m_coef[j] = fac * sum;
	}
      return;
    }


  /**
   * Constructs a Chebyshev expansion from a polynomial.
   *
   * A polynomial valid in the range x=a to x=b is specified by coefficients d[0..m-1].
   * An equivalent array of Chebyshev coefficients, c[0..m-1], is output.
   * The index mfew will be set to the index of the first nonzero Chebyshev coefficient smaller than err.
   * Then the polynomial is approximated by the first mfew cefficients c[0..mfew-1].
   */
  template<typename Tp>
    template<typename Up>
    Chebyshev<Tp>::Chebyshev(Up a, Up b, const emsr::Polynomial<Up>& poly)
    : m_lower(a),
      m_upper(b),
      m_coef(poly.order())
    {
      if (this->upper() == this->lower())
	throw std::domain_error("Chebyshev: input domain has zero length");

      //  Map the interval of the polynomial from x:[a,b] to y:[-1,+1].
      poly_shift_coeff((-Tp{2} - this->upper() - this->lower())
		     / (this->upper() - this->lower()),
		       (+Tp{2} - this->upper() - this->lower())
		     / (this->upper() - this->lower()), poly);

      Tp pow = Tp{1};
      this->m_coef[0] = Tp{2} * poly[0];
      for (int k = 1; k < poly.size(); ++k)
	{
	  this->m_coef[k] = Tp{};
	  auto fac = poly[k] / pow;
	  int jm = k;
	  int jp = 1;
	  for (int j = k; j >= 0; j -= 2, --jm, ++jp)
	    {
	      //  Increment this and lower orders of Chebyshev
	      //  with the combinatorial coefficient times d[k].
	      this->m_coef[j] += fac;
	      fac *= Tp(jm) / Tp(jp);
	    }
	  pow += pow;
	}

      //  Map the interval of the polynomial from y:[-1,+1] to x:[a,b].
      poly_shift_coeff(this->lower(), this->upper(), poly);
    }


  /**
   * Chebyshev evaluation.
   */
  template<typename Tp>
    Tp
    Chebyshev<Tp>::operator()(Tp x) const
    {
      if ((x - this->lower()) * (x - this->upper()) > Tp{})
	throw std::domain_error("Chebyshev::operator(): argument out of range");

      if (this->m_coef.size() == 0)
	return Tp{};

      auto d = Tp{};
      auto dd = Tp{};
      auto y = (Tp{2} * x - this->lower() - this->upper())
	       / (this->upper() - this->lower());
      auto y2 = Tp{2} * y;
      for (unsigned int j = this->m_coef.size() - 1; j >= 1; --j)
	dd = std::exchange(d, y2 * d - dd + this->m_coef[j]);

      return y * d - dd + this->m_coef[0] / Tp{2};
    }


  /**
   * This routine returns a Chebyshev fit of the derivative of this Chebyshev fit.
   */
  template<typename Tp>
    Chebyshev<Tp>
    Chebyshev<Tp>::derivative() const
    {
      if (this->m_coef.size() < 3)
	return Chebyshev<Tp>{};

      unsigned int n = this->m_coef.size();
      std::vector<Tp> cder(n);
      cder[n - 1] = Tp{};
      cder[n - 2] = Tp(2 * (n - 1)) * this->m_coef[n - 1];
      for (unsigned int j = n - 3; j >= 0; --j)
	cder[j] = cder[j + 2] + Tp(2 * (j + 1)) * this->m_coef[j + 1];

      const auto con = Tp{2} / (this->upper() - this->lower());
      for (unsigned int j = 0; j < n; ++j)
	cder[j] *= con;

      return Chebyshev<Tp>(this->lower(), this->upper(), std::move(cder));
    }


  /**
   * This routine returns the Chebyshev fit of the integral of this Chebyshev fit.
   * The constant of integration is set so that the integral vanishes at a.
   */
  template<typename Tp>
    Chebyshev<Tp>
    Chebyshev<Tp>::integral() const
    {
      unsigned int n = this->m_coef.size();
      std::vector<Tp> cint(n);
      auto sum = Tp{};
      auto fact = Tp{1};
      const auto con = (this->upper() - this->lower()) / Tp{4};
      for (unsigned int j = 1; j < n - 1; ++j)
	{
	  //  Accumulate the constant of integration in sum.
	  //  fact will equal +/-1 to alternate the sign of the terms.
	  cint[j] = con
		* (this->m_coef[j - 1] - this->m_coef[j + 1]) / j;
	  sum += fact * cint[j];
	  fact = -fact;
	}
      cint[n - 1] = con * this->m_coef[n - 2] / (n - 1);
      sum += fact * cint[n - 1];

      //  Set the constant of integration.
      cint[0] = Tp{2} * sum;

      return Chebyshev<Tp>(this->lower(), this->upper(), cint);
    }


  /**
   * This routine returns the array d[0..n-1], of coefficients of a polynomial
   * expansion which is equivalent to the Chebyshev fit.
   */
  template<typename Tp>
    emsr::Polynomial<Tp>
    Chebyshev<Tp>::to_polynomial() const
    {
      unsigned int n = this->m_coef.size();

      std::vector<Tp> d(n);
      std::vector<Tp> dd(n);

      d[0] = this->m_coef[n - 1];
      for (unsigned int j = n - 2; j >= 1; --j)
	{
	  for (unsigned int k = n - j; k >= 1; --k)
	    dd[k] = std::exchange(d[k], Tp{2} * d[k - 1] - dd[k]);
	  dd[0] = std::exchange(d[0], -dd[0] + this->m_coef[j]);
	}
      for (unsigned int j = n - 1; j >= 1; --j)
	d[j] = d[j - 1] - dd[j];
      d[0] = -dd[0] + this->m_coef[0] / Tp{2};

      //  Map the interval [-1,+1] to [a,b].
      poly_shift_coeff(this->lower(), this->upper(), d);

      return emsr::Polynomial<Tp>(d);
    }


  /**
   * Truncate the Chebyshev series so that the first nonzero neglected term
   * is less than eps.  Skipping zero coefficients takes care of series
   * with alternating zero and nonzero terms.
   */
  template<typename Tp>
    template<typename Up>
      void
      Chebyshev<Tp>::truncate(Up eps)
      {
	// @todo Uniform container erasure would be better.
	auto beg = std::find_if(std::begin(this->m_coef),
				  std::end(this->m_coef),
			[eps](Tp c)
			{ return c != Tp{} && std::abs(c) < eps; });
	this->m_coef.erase(beg, std::end(this->m_coef));
      }


  /**
   * Write the Chebyshev expansion to an output stream.
   */
  template<typename Tp,
	   typename CharT, typename Traits = std::char_traits<CharT>>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& os,
	       const Chebyshev<Tp>& cheb)
    {
      auto width = std::numeric_limits<Tp>::max_digits10;
      auto prec = os.precision(width);
      os << '{';
      auto begin = std::begin(cheb.m_coef);
      auto end = std::end(cheb.m_coef);
      if (begin != end)
	{
	  os << *begin;
	  for (++begin; begin != end; ++begin)
	    os << ',' << *begin;
	}
      os << '}';
      os.precision(prec);
      return os;
    }

  /**
   * Economize a polynomial poly over a range [a, b] by returning
   * a new lower-order polynomial so that evaluations of the new polynomial
   * will equal those of the old polynomial within the specified tolerance eps
   * over the requested range.
   */
  template<typename Tp>
    emsr::Polynomial<Tp>
    economize_polynomial(Tp a, Tp b, const emsr::Polynomial<Tp>& poly,
			 Tp eps)
    {
      if (poly.order() < 3)
	throw std::domain_error("poly_economize: input arrays too small");

      //  Convert the polynomial into a Chebyshev series.
      Chebyshev<Tp> cpoly(a, b, poly);

      //  Truncate the Chebyshev series to the first term less than eps.
      cpoly.truncate(eps);

      //  Convert the truncated Chebyshev series back into a polynomial.
      return cpoly.to_polynomial();
    }

} // namespace emsr

#endif // CHEBYSHEV_TCC
