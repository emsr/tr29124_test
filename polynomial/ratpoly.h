
#include <iostream>

#include "polynomial.h"

#ifndef _EXT_RATPOLY_HPP
#define _EXT_RATPOLY_HPP 1

namespace __gnu_cxx
{

  /**
   *
   */
  template<typename _Tp>
    class _RationalPolynomial
    {
    public:
      /**
       *  Typedefs.
       */
      using polynomial_type = _Polynomial<_Tp>;
      using value_type = typename polynomial_type::value_type;
// The above should be _Polynomial<_Tp>to be consistent with rational.h
/* These might not make sense.
      using reference = typename polynomial_type::reference;
      using const_reference = typename polynomial_type::const_reference;
      using pointer = typename polynomial_type::pointer;
      using const_pointer = typename polynomial_type::const_pointer;
      using iterator = typename polynomial_type::iterator;
      using const_iterator = typename polynomial_type::const_iterator;
      using reverse_iterator = typename polynomial_type::reverse_iterator;
      using const_reverse_iterator = typename polynomial_type::const_reverse_iterator;
*/
      using size_type = typename polynomial_type::size_type;
      using difference_type = typename polynomial_type::difference_type;

      /**
       *  Create a zero degree polynomial with value zero.
       */
      _RationalPolynomial()
      : _M_num(), _M_den()
      { }

      /**
       *  Copy ctor.
       */
      _RationalPolynomial(const _RationalPolynomial&) = default;

      _RationalPolynomial(const _Polynomial<_Tp>& __num,
      			  const _Polynomial<_Tp>& __den)
      : _M_num(__num), _M_den(__den)
      { }

      /**
       *  Evaluate the polynomial at the input point.
       */
      value_type
      operator()(value_type __x) const
      { return this->_M_num(__x) / this->_M_den(__x); }

      /**
       *  Unary plus.
       */
      _RationalPolynomial
      operator+() const
      { return *this; }

      /**
       *  Unary minus.
       */
      _RationalPolynomial
      operator-() const
      { return _RationalPolynomial(*this) *= value_type(-1); }

      /**
       *  Copy assignment.
       */
      _RationalPolynomial&
      operator=(const _RationalPolynomial&) = default;

      /**
       *  Add a rational polynomial to this rational polynomial.
       */
      _RationalPolynomial&
      operator+=(const _RationalPolynomial& __x)
      {
        this->numer() = this->numer() * __x.denom() + this->denom() * __x.numer();
        this->denom() *= __x.denom();
	return *this;
      }

      /**
       *  Subtract a rational polynomial from this rational polynomial.
       */
      _RationalPolynomial&
      operator-=(const _RationalPolynomial& __x)
      {
        this->numer() = this->numer() * __x.denom() - this->denom() * __x.numer();
        this->denom() *= __x.denom();
	return *this;
      }

      /**
       *  Multiply this rational polynomial by a rational polynomial.
       */
      _RationalPolynomial&
      operator*=(const _RationalPolynomial& __x)
      {
	this->numer() *= __x.numer();
	this->denom() *= __x.denom();
	return *this;
      }

      /**
       *  Divide this rational polynomial by a rational polynomial.
       */
      _RationalPolynomial&
      operator/=(const _RationalPolynomial& __x)
      {
	this->numer() *= __x.denom();
	this->denom() *= __x.numer();
	return *this;
      }

      const _Polynomial<value_type>&
      numer() const
      { return this->_M_num; }

      _Polynomial<value_type>&
      numer()
      { return this->_M_num; }

      const _Polynomial<value_type>&
      denom() const
      { return this->_M_den; }

      _Polynomial<value_type>&
      denom()
      { return this->_M_den; }

    private:

      _Polynomial<value_type> _M_num;
      _Polynomial<value_type> _M_den;
    };

  /**
   *  Write a polynomial to a stream.
   *  The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename _Tp>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& __os, const _RationalPolynomial<_Tp>& __poly)
    {
      __os << __poly.numer() << "/" << __poly.denom();
      return __os;
    }

  /**
   *  Read a polynomial from a stream.
   *  The input format can be a plain scalar (zero degree polynomial)
   *  or a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename _Tp>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& __is, _RationalPolynomial<_Tp>& __poly)
    {
      _Polynomial<_Tp> __numer, __denom;
      __is >> __numer;
      if (!__is.fail())
	{
	  CharT __ch;
	  __is >> __ch;
	  if (__ch != '/')
	    __is.setstate(std::ios_base::failbit);
	  else
	    {
	      __is >> __denom;
	      if (!__is.fail())
		{
		  __poly.numer() = __numer;
		  __poly.denom() = __denom;
		}
	    }
	}
      return __is;
    }

} // namespace __gnu_cxx

#endif // _EXT_RATPOLY_HPP
