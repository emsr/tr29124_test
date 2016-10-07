// Math extensions -*- C++ -*-

// Copyright (C) 2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file ext/polynomial.h
 *  This file is a GNU extension to the Standard C++ Library.
 */

#ifndef _EXT_POLYNOMIAL_H
#define _EXT_POLYNOMIAL_H 1

#pragma GCC system_header

#if __cplusplus < 201103L
# include <bits/c++0x_warning.h>
#else

#include <initializer_list>
#include <vector>
#include <iosfwd>
#include <limits>
#include <array>
#include <utility> // For exchange

namespace std {
  template<typename _Tp>
    class complex;
}

/**
 *  detail: Do we want this to always have a size of at least one? a_0 = _Tp{}?  YES.
 *  detail: Should I punt on the initial power?  YES.
 *
 *  If high degree coefficients are zero, should I resize down? YES (or provide another word for order).
 *  How to access coefficients (bikeshed)?
 *    poly[i];
 *    coefficient(i);
 *    operator[](int i);
 *    begin(), end()?
 *    const _Tp* coefficients(); // Access for C, Fortran.
 *  How to set individual coefficients?
 *    poly[i] = c;
 *    coefficient(i, c);
 *    coefficient(i) = c;
 *  How to handle division?
 *    operator/ and throw out remainder?
 *    operator% to return the remainder?
 *    std::pair<> div(const _Polynomial& __a, const _Polynomial& __b) or remquo.
 *    void divmod(const _Polynomial& __a, const _Polynomial& __b,
 *                _Polynomial& __q, _Polynomial& __r);
 *  Should factory methods like derivative and integral be globals?
 *  I could have members:
 *    _Polynomial& integrate(_Tp c);
 *    _Polynomial& differentiate();
 *  Largest coefficient:
 *    Enforce coefficient of largest power be nonzero?
 *    Return an 'effective' order? Larest nonzero coefficient?
 *    Monic polynomial has largest coefficient as 1.  Subclass?
 */
namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *
   */
  template<typename _Tp>
    class _Polynomial
    {
    public:
      /**
       *  Typedefs.
       *  @todo Should we grab these from _M_coeff (i.e. std::vector<_Tp>)?
       */
      using value_type = typename std::vector<_Tp>::value_type;
      using reference = typename std::vector<value_type>::reference;
      using const_reference = typename std::vector<value_type>::const_reference;
      using pointer = typename std::vector<value_type>::pointer;
      using const_pointer = typename std::vector<value_type>::const_pointer;
      using iterator = typename std::vector<value_type>::iterator;
      using const_iterator = typename std::vector<value_type>::const_iterator;
      using reverse_iterator = typename std::vector<value_type>::reverse_iterator;
      using const_reverse_iterator = typename std::vector<value_type>::const_reverse_iterator;
      using size_type = typename std::vector<_Tp>::size_type;
      using difference_type = typename std::vector<_Tp>::difference_type;

      /**
       *  Create a zero degree polynomial with value zero.
       */
      _Polynomial()
      : _M_coeff(1)
      { }

      /**
       *  Copy ctor.
       */
      _Polynomial(const _Polynomial&) = default;

      template<typename _Up>
	_Polynomial(const _Polynomial<_Up>& __poly)    
	: _M_coeff{}
	{
          for (const auto __c : __poly)
	    this->_M_coeff.push_back(static_cast<value_type>(__c));
	}

      /**
       *  Create a polynomial of just one constant term.
       */
      explicit
      _Polynomial(value_type __a, size_type __degree = 0)
      : _M_coeff(__degree + 1)
      { this->_M_coeff[__degree] = __a; }

      /**
       *  Create a polynomial from an initializer list of coefficients.
       */
      _Polynomial(std::initializer_list<value_type> __ila)
      : _M_coeff(__ila)
      { }

      /**
       *  Create a polynomial from an input iterator range of coefficients.
       */
      template<typename InIter,
	       typename = std::_RequireInputIter<InIter>>
	_Polynomial(const InIter& __abegin, const InIter& __aend)
	: _M_coeff(__abegin, __aend)
	{ }

      /**
       *  Use Lagrange interpolation to construct a polynomial passing through
       *  the data points.  The degree will be one less than the number of points.
       */
      template<typename InIter,
	       typename = std::_RequireInputIter<InIter>>
	_Polynomial(const InIter& __xbegin, const InIter& __xend,
		    const InIter& __ybegin)
	: _M_coeff()
	{
	  std::vector<_Polynomial<value_type>> __numer;
	  std::vector<_Polynomial<value_type>> __denom;
	  for (auto __xi = __xbegin; __xi != __xend; ++__xi)
	    {
	      for (auto __xj = __xi + 1; __xj != __xend; ++__xj)
		__denom.push_back(value_type(*__xj) - value_type(*__xi));
	      __numer.push_back({-value_type(*__xi), value_type{1}});
	    }
	}

      /**
       *  Swap the polynomial with another polynomial.
       */
      void
      swap(_Polynomial& __poly)
      { this->_M_coeff.swap(__poly._M_coeff); }

      /**
       *  Evaluate the polynomial at the input point.
       */
      value_type
      operator()(value_type __x) const
      {
	if (this->degree() > 0)
	  {
	    value_type __poly(this->coefficient(this->degree()));
	    for (int __i = this->degree() - 1; __i >= 0; --__i)
	      __poly = __poly * __x + this->coefficient(__i);
	    return __poly;
	  }
	else
	  return value_type{};
      }

      /**
       *  Evaluate the polynomial at the input point.
       */
      template<typename _Up>
	auto
	operator()(_Up __x) const
	-> decltype(value_type{} * _Up{})
	{
	  if (this->degree() > 0)
	    {
	      auto __poly(_Up{1} * this->coefficient(this->degree()));
	      for (int __i = this->degree() - 1; __i >= 0; --__i)
		__poly = __poly * __x + this->coefficient(__i);
	      return __poly;
	    }
	  else
	    return value_type{} * _Up{};
	}

      /**
       *  Evaluate the polynomial using a modification of Horner's rule which
       *  exploits the fact that the polynomial coefficients are all real.
       *
       *  The algorithm is discussed in detail in:
       *  Knuth, D. E., The Art of Computer Programming: Seminumerical
       *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
       * 
       *  If n is the degree of the polynomial,
       *  n - 3 multiplies and 4 * n - 6 additions are saved.
       */
      template<typename _Up>
	auto
	operator()(const std::complex<_Up>& __z) const
	-> decltype(value_type{} * std::complex<_Up>{});

      /**
       *  Evaluate the polynomial at a range of input points.
       *  The output is written to the output iterator which
       *  must be large enough to contain the results.
       *  The next available output iterator is returned.
       */
      template<typename InIter, typename OutIter,
	       typename = std::_RequireInputIter<InIter>>
	OutIter
	operator()(const InIter& __xbegin, const InIter& __xend,
        	   OutIter& __pbegin) const
	{
	  for (; __xbegin != __xend; ++__xbegin)
	    __pbegin++ = (*this)(__xbegin++);
	  return __pbegin;
	}

      //  Could/should this be done by output iterator range?
      template<size_type N>
	void
	eval(value_type __x, std::array<value_type, N>& __arr);

      /**
       *  Evaluate the polynomial and its derivatives at the point x.
       *  The values are placed in the output range starting with the
       *  polynomial value and continuing through higher derivatives.
       */
      template<typename OutIter>
	void
	eval(value_type __x, OutIter __b, OutIter __e);

      /**
       *  Evaluate the even part of the polynomial at the input point.
       */
      value_type
      even(value_type __x) const;

      /**
       *  Evaluate the odd part of the polynomial at the input point.
       */
      value_type
      odd(value_type __x) const;

      /**
       *  Evaluate the even part of the polynomial using a modification
       *  of Horner's rule which exploits the fact that the polynomial
       *  coefficients are all real.
       *
       *  The algorithm is discussed in detail in:
       *  Knuth, D. E., The Art of Computer Programming: Seminumerical
       *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
       * 
       *  If n is the degree of the polynomial,
       *  n - 3 multiplies and 4 * n - 6 additions are saved.
       */
      template<typename _Up>
	auto
	even(const std::complex<_Up>& __z) const
	-> decltype(value_type{} * std::complex<_Up>{});

      /**
       *  Evaluate the odd part of the polynomial using a modification
       *  of Horner's rule which exploits the fact that the polynomial
       *  coefficients are all real.
       *
       *  The algorithm is discussed in detail in:
       *  Knuth, D. E., The Art of Computer Programming: Seminumerical
       *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
       * 
       *  If n is the degree of the polynomial,
       *  n - 3 multiplies and 4 * n - 6 additions are saved.
       */
      template<typename _Up>
	auto
	odd(const std::complex<_Up>& __z) const
	-> decltype(value_type{} * std::complex<_Up>{});

      /**
       *  Return the derivative of the polynomial.
       */
      _Polynomial
      derivative() const
      {
	_Polynomial __res(value_type{},
			  this->degree() > 0UL ? this->degree() - 1 : 0UL);
	for (size_type __n = this->degree(), __i = 1; __i <= __n; ++__i)
	  __res._M_coeff[__i - 1] = __i * _M_coeff[__i];
	return __res;
      }

      /**
       *  Return the integral of the polynomial with given integration constant.
       */
      _Polynomial
      integral(value_type __c = value_type{}) const
      {
	_Polynomial __res(value_type{}, this->degree() + 1);
	__res._M_coeff[0] = __c;
	for (size_type __n = this->degree(), __i = 0; __i <= __n; ++__i)
	  __res._M_coeff[__i + 1] = _M_coeff[__i] / value_type(__i + 1);
	return __res;
      }

      /**
       *  Unary plus.
       */
      _Polynomial
      operator+() const
      { return *this; }

      /**
       *  Unary minus.
       */
      _Polynomial
      operator-() const
      { return _Polynomial(*this) *= value_type(-1); }

      /**
       *  Assign from a scalar.
       *  The result is a zero degree polynomial equal to the scalar.
       */
      _Polynomial&
      operator=(const value_type& __x)
      {
	_M_coeff = {__x};
	return *this;
      }

      /**
       *  Copy assignment.
       */
      _Polynomial&
      operator=(const _Polynomial&) = default;

      template<typename _Up>
	_Polynomial&
	operator=(const _Polynomial<_Up>& __poly)
	{
	  if (&__poly != this)
	    {
	      this->_M_coeff.clear();
	      for (const auto __c : __poly)
		this->_M_coeff = static_cast<value_type>(__c);
	      return *this;
	    }
	}

      /**
       *  Assign from an initialiser list.
       */
      _Polynomial&
      operator=(std::initializer_list<value_type> __ila)
      {
	_M_coeff = __ila;
	return *this;
      }

      /**
       *  Add a scalar to the polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator+=(const _Up& __x)
	{
	  degree(this->degree()); // Resize if necessary.
	  _M_coeff[0] += static_cast<value_type>(__x);
	  return *this;
	}

      /**
       *  Subtract a scalar from the polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator-=(const _Up& __x)
	{
	  degree(this->degree()); // Resize if necessary.
	  _M_coeff[0] -= static_cast<value_type>(__x);
	  return *this;
	}

      /**
       *  Multiply the polynomial by a scalar.
       */
      template<typename _Up>
	_Polynomial&
	operator*=(const _Up& __x)
	{
	  degree(this->degree()); // Resize if necessary.
	  for (size_type __i = 0; __i < _M_coeff.size(); ++__i)
	    _M_coeff[__i] *= static_cast<value_type>(__x);
	  return *this;
	}

      /**
       *  Divide the polynomial by a scalar.
       */
      template<typename _Up>
	_Polynomial&
	operator/=(const _Up& __x)
	{
	  for (size_type __i = 0; __i < _M_coeff.size(); ++__i)
	    this->_M_coeff[__i] /= static_cast<value_type>(__x);
	  return *this;
	}

      /**
       *  Take the modulus of the polynomial relative to a scalar.
       *  The result is always null.
       */
      template<typename _Up>
	_Polynomial&
	operator%=(const _Up&)
	{
	  degree(0UL); // Resize.
	  this->_M_coeff[0] = value_type{};
	  return *this;
	}

      /**
       *  Add another polynomial to the polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator+=(const _Polynomial<_Up>& __poly)
	{
	  this->degree(std::max(this->degree(), __poly.degree()));
	  for (size_type __n = __poly.degree(), __i = 0; __i <= __n; ++__i)
	    this->_M_coeff[__i] += static_cast<value_type>(__poly._M_coeff[__i]);
	  return *this;
	}

      /**
       *  Subtract another polynomial from the polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator-=(const _Polynomial<_Up>& __poly)
	{
	  this->degree(std::max(this->degree(), __poly.degree())); // Resize if necessary.
	  for (size_type __n = __poly.degree(), __i = 0; __i <= __n; ++__i)
	    this->_M_coeff[__i] -= static_cast<value_type>(__poly._M_coeff[__i]);
	  return *this;
	}

      /**
       *  Multiply the polynomial by another polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator*=(const _Polynomial<_Up>& __poly);

      /**
       *  Divide the polynomial by another polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator/=(const _Polynomial<_Up>& __poly)
	{
	  _Polynomial<value_type >__quo, __rem;
	  divmod(*this, __poly, __quo, __rem);
	  *this = __quo;
	  return *this;
	}

      /**
       *  Take the modulus of (modulate?) the polynomial relative to another polynomial.
       */
      template<typename _Up>
	_Polynomial&
	operator%=(const _Polynomial<_Up>& __poly)
	{
	  _Polynomial<value_type >__quo, __rem;
	  divmod(*this, __poly, __quo, __rem);
	  *this = __rem;
	  return *this;
	}

      /**
       *  Return the degree or the power of the largest coefficient.
       */
      size_type
      degree() const
      { return (this->_M_coeff.size() > 0 ? this->_M_coeff.size() - 1 : 0); }

      /**
       *  Set the degree or the power of the largest coefficient.
       */
      void
      degree(size_type __degree)
      { this->_M_coeff.resize(__degree + 1UL); }

      /**
       *  Return the size of the coefficient sequence.
       */
      size_type
      size() const
      { return this->_M_coeff.size(); }


      /**
       *  Return the @c ith coefficient with range checking.
       */
      value_type
      coefficient(size_type __i) const
      { return this->_M_coeff.at(__i); }

      /**
       *  Set coefficient @c i to @c val with range checking.
       */
      void
      coefficient(size_type __i, value_type __val)
      { this->_M_coeff.at(__i) = __val; }

      /**
       *  Return a @c const pointer to the coefficient sequence.
       */
      const value_type*
      coefficients() const
      { this->_M_coeff.data(); }

      /**
       *  Return coefficient @c i.
       */
      value_type
      operator[](size_type __i) const
      { return this->_M_coeff[__i]; }

      /**
       *  Return coefficient @c i as an lvalue.
       */
      reference
      operator[](size_type __i)
      { return this->_M_coeff[__i]; }

      /**
       *  Return an iterator to the beginning of the coefficient sequence.
       */
      iterator
      begin()
      { return this->_M_coeff.begin(); }

      /**
       *  Return an iterator to one past the end of the coefficient sequence.
       */
      iterator
      end()
      { return this->_M_coeff.end(); }

      /**
       *  Return a @c const iterator the beginning
       *  of the coefficient sequence.
       */
      const_iterator
      begin() const
      { return this->_M_coeff.begin(); }

      /**
       *  Return a @c const iterator to one past the end
       *  of the coefficient sequence.
       */
      const_iterator
      end() const
      { return this->_M_coeff.end(); }

      /**
       *  Return a @c const iterator the beginning
       *  of the coefficient sequence.
       */
      const_iterator
      cbegin() const
      { return this->_M_coeff.cbegin(); }

      /**
       *  Return a @c const iterator to one past the end
       *  of the coefficient sequence.
       */
      const_iterator
      cend() const
      { return this->_M_coeff.cend(); }

      reverse_iterator
      rbegin()
      { return this->_M_coeff.rbegin(); }

      reverse_iterator
      rend()
      { return this->_M_coeff.rend(); }

      const_reverse_iterator
      rbegin() const
      { return this->_M_coeff.rbegin(); }

      const_reverse_iterator
      rend() const
      { return this->_M_coeff.rend(); }

      const_reverse_iterator
      crbegin() const
      { return this->_M_coeff.crbegin(); }

      const_reverse_iterator
      crend() const
      { return this->_M_coeff.crend(); }

      template<typename CharT, typename Traits, typename _Tp1>
	friend std::basic_istream<CharT, Traits>&
	operator>>(std::basic_istream<CharT, Traits>&, _Polynomial<_Tp1>&);

      template<typename _Tp1>
	friend bool
	operator==(const _Polynomial<_Tp1>& __pa,
		   const _Polynomial<_Tp1>& __pb);

    private:

      std::vector<value_type> _M_coeff;
    };

  /**
   *  Return the sum of a polynomial with a scalar.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() + _Up())>
    operator+(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() + _Up())>(__poly) += __x; }

  /**
   *  Return the difference of a polynomial with a scalar.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() - _Up())>
    operator-(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() - _Up())>(__poly) -= __x; }

  /**
   *  Return the product of a polynomial with a scalar.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() * _Up())>
    operator*(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() * _Up())>(__poly) *= __x; }

  /**
   *  Return the quotient of a polynomial with a scalar.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator/(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() / _Up())>(__poly) /= __x; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator%(const _Polynomial<_Tp>& __poly, const _Up& __x)
    { return _Polynomial<decltype(_Tp() / _Up())>(__poly) %= __x; }

  /**
   *  Return the sum of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() + _Up())>
    operator+(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() + _Up())>(__pa) += __pb; }

  /**
   *  Return the difference of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() - _Up())>
    operator-(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() - _Up())>(__pa) -= __pb; }

  /**
   *  Return the product of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() * _Up())>
    operator*(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() * _Up())>(__pa) *= __pb; }

  /**
   *  Return the quotient of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator/(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() / _Up())>(__pa) /= __pb; }

  /**
   *  Return the modulus or remainder of one polynomial relative to another one.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator%(const _Polynomial<_Tp>& __pa, const _Polynomial<_Up>& __pb)
    { return _Polynomial<decltype(_Tp() / _Up())>(__pa) %= __pb; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() + _Up())>
    operator+(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() + _Up())>(__x) += __poly; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() - _Up())>
    operator-(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() - _Up())>(__x) -= __poly; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() * _Up())>
    operator*(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() * _Up())>(__x) *= __poly; }

  /**
   *  Return the quotient of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator/(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() / _Up())>(__x) /= __poly; }

  /**
   *  Return the modulus or remainder of one polynomial relative to another one.
   */
  template<typename _Tp, typename _Up>
    inline _Polynomial<decltype(_Tp() / _Up())>
    operator%(const _Tp& __x, const _Polynomial<_Up>& __poly)
    { return _Polynomial<decltype(_Tp() / _Up())>(__x) %= __poly; }

  /**
   *  Divide two polynomials returning the quotient and remainder.
   */
  template<typename _Tp>
    void
    divmod(const _Polynomial<_Tp>& __pa, const _Polynomial<_Tp>& __pb,
           _Polynomial<_Tp>& __quo, _Polynomial<_Tp>& __rem);

  /**
   *  Write a polynomial to a stream.
   *  The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename _Tp>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& __os,
	       const _Polynomial<_Tp>& __poly);

  /**
   *  Read a polynomial from a stream.
   *  The input format can be a plain scalar (zero degree polynomial)
   *  or a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename _Tp>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& __is,
	       _Polynomial<_Tp>& __poly);

  /**
   *  Return true if two polynomials are equal.
   */
  template<typename _Tp>
    inline bool
    operator==(const _Polynomial<_Tp>& __pa, const _Polynomial<_Tp>& __pb)
    { return __pa._M_coeff == __pb._M_coeff; }

  /**
   *  Return false if two polynomials are equal.
   */
  template<typename _Tp>
    inline bool
    operator!=(const _Polynomial<_Tp>& __pa, const _Polynomial<_Tp>& __pb)
    { return !(__pa == __pb); }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#include <ext/polynomial.tcc>

#endif // C++11

#endif // _EXT_POLYNOMIAL_H

