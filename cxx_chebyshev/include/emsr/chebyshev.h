#ifndef CHEBYSHEV_H
#define CHEBYSHEV_H 1

#include <initializer_list>
#include <vector>
#include <iosfwd>

#include <emsr/polynomial.h>

/**
 * @file chebyshev.h Interface for C++ Chebyshev methods.
 */

namespace emsr
{

  /**
   * @brief  _Chebyshev represents a Chebyshev fit of a function.
   * Given a function @c func, the lower limit @c a, the upper limit @c b,
   * and the number of points @c n _Chebyshev contains the coefficients
   * of a Chebyshev polynomial expansion such that
   * @f[
   *    f(x) = k\sum_{0}^{n-1} c_k T_k(x) - c_0/2.
   * @f]
   * Use these constructors with moderately large n of 30 or 50.
   * Then truncate the Chebyshev fit to a smaller number of terms to satisfy
   * the accuracy requirements.
   */
  template<typename Tp>
    class Chebyshev
    {
    public:

      using value_type = Tp;

      Chebyshev()
      : m_lower{-1}, m_upper{+1},
	m_coef{}
      { }

      Chebyshev(Tp a, Tp b, unsigned n, Tp func(Tp));

      Chebyshev(Tp a, Tp b, std::initializer_list<Tp> il)
      : m_lower{a}, m_upper{b},
	m_coef{il}
      { }

      Chebyshev(Tp a, Tp b, const std::vector<Tp>& coeff)
      : m_lower{a}, m_upper{b},
	m_coef{coeff}
      { }

      Chebyshev(Tp a, Tp b, std::vector<Tp>&& coeff)
      : m_lower{a}, m_upper{b},
	m_coef{std::move(coeff)}
      { }

      template<typename Up>
        Chebyshev(Up a, Up b, const emsr::Polynomial<Up>& poly);

      Chebyshev derivative() const;

      Chebyshev integral() const;

      std::size_t
      order() const
      { return this->m_coef.size(); }

      value_type
      lower() const
      { return this->m_lower; }

      value_type
      upper() const
      { return this->m_upper; }

      emsr::Polynomial<Tp>
      to_polynomial() const;

      value_type operator()(value_type x) const;

      template<typename Up>
        void truncate(Up eps = std::numeric_limits<Up>::epsilon());

      template<typename Tp1, typename CharT, typename Traits>
	friend std::basic_ostream<CharT, Traits>&
	operator<<(std::basic_ostream<CharT, Traits>& os,
		   const Chebyshev<Tp1>& cheb);

    private:
      Tp m_lower;
      Tp m_upper;
      std::vector<Tp> m_coef;
    };

} // namespace emsr

#include <emsr/chebyshev.tcc>

#endif // CHEBYSHEV_H
