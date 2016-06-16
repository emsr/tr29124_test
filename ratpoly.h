#include <bits/polynomial.h>

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
       *  @todo Should we grab these from _M_coeff (i.e. std::vector<_Tp>)?
       */
      using polynomial_type = typename _Polynomial<_Tp>;
      using value_type = typename polynomial_type::value_type;
      using reference = typename polynomial_type::reference;
      using const_reference = typename polynomial_type::const_reference;
      using pointer = typename polynomial_type::pointer;
      using const_pointer = typename polynomial_type::const_pointer;
      using iterator = typename polynomial_type::iterator;
      using const_iterator = typename polynomial_type::const_iterator;
      using reverse_iterator = typename polynomial_type::reverse_iterator;
      using const_reverse_iterator = typename polynomial_type::const_reverse_iterator;
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
      _RationalPolynomial(const _Polynomial&) = default;

      _RationalPolynomial(const _Polynomial<_Up>& __num,
      			  const _Polynomial<_Up>& __den)
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
        this->_M_num = this->_M_num * __x._M_den + this->_M_den * __x._M_num;
        this->_M_den *= __x._M_den;
	return *this;
      }

      /**
       *  Subtract a rational polynomial from this rational polynomial.
       */
      _RationalPolynomial&
      operator-=(const _RationalPolynomial& __x)
      {
        this->_M_num = this->_M_num * __x._M_den - this->_M_den * __x._M_num;
        this->_M_den *= __x._M_den;
	return *this;
      }

      /**
       *  Multiply this rational polynomial by a rational polynomial.
       */
      _RationalPolynomial&
      operator*=(const _RationalPolynomial& __x)
      {
	this->_M_num *= __x._M_num;
	this->_M_den *= __x._M_den;
	return *this;
      }

      /**
       *  Divide this rational polynomial by a rational polynomial.
       */
      _RationalPolynomial&
      operator/=(const _RationalPolynomial& __x)
      {
	this->_M_num *= __x._M_den;
	this->_M_den *= __x._M_num;
	return *this;
      }

    private:
      _Polynomial<value_type> _M_num;
      _Polynomial<value_type> _M_den;
    }
#endif // _EXT_RATPOLY_HPP
