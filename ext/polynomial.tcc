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

/** @file ext/polynomial.tcc
 *  This file is a GNU extension to the Standard C++ Library.
 */

#ifndef _EXT_POLYNOMIAL_TCC
#define _EXT_POLYNOMIAL_TCC 1

#pragma GCC system_header

#if __cplusplus < 201103L
# include <bits/c++0x_warning.h>
#else

#include <ios>
#include <complex>

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

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
  template<typename _Tp>
    template<typename _Up>
      auto
      _Polynomial<_Tp>::operator()(const std::complex<_Up>& __z) const
      -> decltype(_Polynomial<_Tp>::value_type{} * std::complex<_Up>{})
      {
	const auto __r = _Tp{2} * std::real(__z);
	const auto __s = std::norm(__z);
	size_type __n = this->degree();
	auto __aa = this->coefficient(__n);
	auto __bb = this->coefficient(__n - 1);
	for (size_type __j = 2; __j <= __n; ++__j)
	  __bb = this->coefficient(__n - __j)
	       - __s * std::exchange(__aa, __bb + __r * __aa);
	return __aa * __z + __bb;
      };

    //  Could/should this be done by output iterator range?
  template<typename _Tp>
    template<typename _Polynomial<_Tp>::size_type N>
      void
      _Polynomial<_Tp>::eval(typename _Polynomial<_Tp>::value_type __x,
			     std::array<_Polynomial<_Tp>::value_type, N>& __arr)
      {
	if (__arr.size() > 0)
	  {
	    __arr.fill(value_type{});
	    const size_type __sz = _M_coeff.size();
	    __arr[0] = this->coefficient(__sz - 1);
            for (int __i = __sz - 2; __i >= 0; --__i)
	      {
		int __nn = std::min(__arr.size() - 1, __sz - 1 - __i);
		for (int __j = __nn; __j >= 1; --__j)
		  __arr[__j] = __arr[__j] * __x + __arr[__j - 1];
		__arr[0] = __arr[0] * __x + this->coefficient(__i);
	      }
	    //  Now put in the factorials.
	    value_type __fact = value_type(1);
	    for (size_type __n = __arr.size(), __i = 2; __i < __n; ++__i)
	      {
		__fact *= value_type(__i);
		__arr[__i] *= __fact;
	      }
	  }
      }

  /**
   *  Evaluate the polynomial and its derivatives at the point x.
   *  The values are placed in the output range starting with the
   *  polynomial value and continuing through higher derivatives.
   */
  template<typename _Tp>
    template<typename OutIter>
      void
      _Polynomial<_Tp>::eval(typename _Polynomial<_Tp>::value_type __x,
			     OutIter __b, OutIter __e)
      {
	if(__b != __e)
	  {
	    std::fill(__b, __e, value_type{});
	    const size_type __sz = _M_coeff.size();
	    *__b = _M_coeff[__sz - 1];
            for (int __i = __sz - 2; __i >= 0; --__i)
	      {
		for (auto __it = std::reverse_iterator<OutIter>(__e);
		     __it != std::reverse_iterator<OutIter>(__b) - 1; ++__it)
		  *__it = *__it * __x + *(__it + 1);
		*__b = *__b * __x + _M_coeff[__i];
	      }
	    //  Now put in the factorials.
	    int __i = 0;
	    value_type __fact = value_type(++__i);
	    for (auto __it = __b + 1; __it != __e; ++__it)
	      {
		__fact *= value_type(__i);
		*__it *= __fact;
		++__i;
	      }
	  }
      }

  /**
   *  Evaluate the even part of the polynomial at the input point.
   */
  template<typename _Tp>
    typename _Polynomial<_Tp>::value_type
    _Polynomial<_Tp>::even(typename _Polynomial<_Tp>::value_type __x) const
    {
      if (this->degree() > 0)
	{
	  auto __odd = this->degree() % 2;
	  value_type __poly(this->coefficient(this->degree() - __odd));
	  for (int __i = this->degree() - __odd - 2; __i >= 0; __i -= 2)
	    __poly = __poly * __x * __x + this->coefficient(__i);
	  return __poly;
	}
      else
	return value_type{};
    }

  /**
   *  Evaluate the odd part of the polynomial at the input point.
   */
  template<typename _Tp>
    typename _Polynomial<_Tp>::value_type
    _Polynomial<_Tp>::odd(typename _Polynomial<_Tp>::value_type __x) const
    {
      if (this->degree() > 0)
	{
	  auto __even = (this->degree() % 2 == 0 ? 1 : 0);
	  value_type __poly(this->coefficient(this->degree() - __even));
	  for (int __i = this->degree() - __even - 2; __i >= 0; __i -= 2)
	    __poly = __poly * __x * __x + this->coefficient(__i);
	  return __poly * __x;
	}
      else
	return value_type{};
    }

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
  template<typename _Tp>
    template<typename _Up>
      auto
      _Polynomial<_Tp>::even(const std::complex<_Up>& __z) const
      -> decltype(typename _Polynomial<_Tp>::value_type{} * std::complex<_Up>{})
      {
	if (this->degree() > 0)
	  {
	    const auto __zz = __z * __z;
	    const auto __r = _Tp{2} * std::real(__zz);
	    const auto __s = std::norm(__zz);
	    auto __odd = this->degree() % 2;
	    size_type __n = this->degree() - __odd;
	    auto __aa = this->coefficient(__n);
	    auto __bb = this->coefficient(__n - 2);
	    for (size_type __j = 4; __j <= __n; __j += 2)
	      __bb = this->coefficient(__n - __j)
		   - __s * std::exchange(__aa, __bb + __r * __aa);
	    return __aa * __zz + __bb;
	  }
	else
	  return decltype(value_type{} * std::complex<_Up>{}){};
      };

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
  template<typename _Tp>
    template<typename _Up>
      auto
      _Polynomial<_Tp>::odd(const std::complex<_Up>& __z) const
      -> decltype(_Polynomial<_Tp>::value_type{} * std::complex<_Up>{})
      {
	if (this->degree() > 0)
	  {
	    const auto __zz = __z * __z;
	    const auto __r = _Tp{2} * std::real(__zz);
	    const auto __s = std::norm(__zz);
	    auto __even = (this->degree() % 2 == 0 ? 1 : 0);
	    size_type __n = this->degree() - __even;
	    auto __aa = this->coefficient(__n);
	    auto __bb = this->coefficient(__n - 2);
	    for (size_type __j = 4; __j <= __n; __j += 2)
	      __bb = this->coefficient(__n - __j)
		   - __s * std::exchange(__aa, __bb + __r * __aa);
	    return __z * (__aa * __zz + __bb);
	  }
	else
	  return decltype(value_type{} * std::complex<_Up>{}){};
      };

    /**
     *  Multiply the polynomial by another polynomial.
     */
  template<typename _Tp>
    template<typename _Up>
      _Polynomial<_Tp>&
      _Polynomial<_Tp>::operator*=(const _Polynomial<_Up>& __poly)
      {
	//  Test for zero size polys and do special processing?
	const size_type __m = this->degree();
	const size_type __n = __poly.degree();
	std::vector<value_type> __new_coeff(__m + __n + 1);
	for (size_type __i = 0; __i <= __m; ++__i)
	  for (size_type __j = 0; __j <= __n; ++__j)
	    __new_coeff[__i + __j] += this->_M_coeff[__i]
				* static_cast<value_type>(__poly._M_coeff[__j]);
	this->_M_coeff = __new_coeff;
	return *this;
      }

  /**
   *  Divide two polynomials returning the quotient and remainder.
   */
  template<typename _Tp>
    void
    divmod(const _Polynomial<_Tp>& __pa, const _Polynomial<_Tp>& __pb,
           _Polynomial<_Tp>& __quo, _Polynomial<_Tp>& __rem)
    {
      __rem = __pa;
      __quo = _Polynomial<_Tp>(_Tp(), __pa.degree());
      const std::size_t __na = __pa.degree();
      const std::size_t __nb = __pb.degree();
      if (__nb <= __na)
	{
	  for (int __k = __na - __nb; __k >= 0; --__k)
	    {
	      __quo.coefficient(__k, __rem.coefficient(__nb + __k)
				   / __pb.coefficient(__nb));
	      for (int __j = __nb + __k - 1; __j >= __k; --__j)
		__rem.coefficient(__j, __rem.coefficient(__j)
				     - __quo.coefficient(__k)
				     * __pb.coefficient(__j - __k));
	    }
	  for (int __j = __nb; __j <= __na; ++__j)
	    __rem.coefficient(__j, _Tp());
	}
    }

  /**
   *  Write a polynomial to a stream.
   *  The format is a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename _Tp>
    std::basic_ostream<CharT, Traits>&
    operator<<(std::basic_ostream<CharT, Traits>& __os,
	       const _Polynomial<_Tp>& __poly)
    {
      int __old_prec = __os.precision(std::numeric_limits<_Tp>::max_digits10);
      __os << "(";
      for (size_t __i = 0; __i < __poly.degree(); ++__i)
        __os << __poly.coefficient(__i) << ",";
      __os << __poly.coefficient(__poly.degree());
      __os << ")";
      __os.precision(__old_prec);
      return __os;
    }

  /**
   *  Read a polynomial from a stream.
   *  The input format can be a plain scalar (zero degree polynomial)
   *  or a parenthesized comma-delimited list of coefficients.
   */
  template<typename CharT, typename Traits, typename _Tp>
    std::basic_istream<CharT, Traits>&
    operator>>(std::basic_istream<CharT, Traits>& __is,
	       _Polynomial<_Tp>& __poly)
    {
      _Tp __x;
      CharT __ch;
      __is >> __ch;
      if (__ch == '(')
	{
	  do
	    {
	      __is >> __x >> __ch;
	      __poly._M_coeff.push_back(__x);
	    }
	  while (__ch == ',');
	  if (__ch != ')')
	    __is.setstate(std::ios_base::failbit);
	}
      else
	{
	  __is.putback(__ch);
	  __is >> __x;
	  __poly = __x;
	}
      return __is;
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __gnu_cxx

#endif // C++11

#endif // _EXT_POLYNOMIAL_TCC

