
#include <initializer_list>
#include <vector>
#include <iosfwd>
#include <limits>
#include <array>
#include <ios>
#include <cmath>
#include <complex>

/**
 *  detail: Do we want this to always have a size of at least one? a_0 = _Tp()?  NO.
 *  detail: Should I punt on the initial power?  YES.
 *
 *  If high degree coefficients are zero, should I resize down?
 *  How to access coefficients (bikeshed)?
 *    poly[i];
 *    coefficient(i);
 *    operator[](int i);
 *    begin(), end()?
 *  How to set individual coefficients?
 *    poly[i] = c;
 *    coefficient(i, c);
 *    coefficient(i) = c;
 *  How to handle division?
 *    operator/ and throw out remainder?
 *    operator% to return the remainder?
 *    std::pair<> div(const _Polynomial& __a, const _Polynomial& __b)
 *    void divmod(const _Polynomial& __a, const _Polynomial& __b,
 *                _Polynomial& __q, _Polynomial& __r);
 *  Should factory methods like derivative and integral be globals?
 *  I could have members:
 *    _Polynomial& integrate(_Tp c);
 *    _Polynomial& differentiate();
 */
namespace __gnu_cxx
{

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
    using reference = typename std::vector<_Tp>::reference;
    using const_reference = typename std::vector<_Tp>::const_reference;
    using pointer = typename std::vector<_Tp>::pointer;
    using const_pointer = typename std::vector<_Tp>::const_pointer;
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
     *  .
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
	}
      else
	return value_type{};
    }

    /**
     *  Evaluate the polynomial at the input point.
     */
    template<typename _Tp2>
      auto
      operator()(_Tp2 __x) const
      -> decltype(value_type{} * _Tp2())
      {
	if (this->degree() > 0)
	  {
	    auto __poly(this->coefficient(this->degree()) * _Tp2(1));
	    for (int __i = this->degree() - 1; __i >= 0; --__i)
	      __poly = __poly * __x + this->coefficient(__i);
	  }
	else
	  return value_type{};
      }

    /**
     *  The following polynomial evaluations are done using
     *  a modified of Horner's rule which exploits the fact that
     *  the polynomial coefficients are all real.
     *  The algorithm is discussed in detail in:
     *  Knuth, D. E., The Art of Computer Programming: Seminumerical
     *  Algorithms (Vol. 2) Third Ed., Addison-Wesley, pp 486-488, 1998.
     * 
     *  If n is the degree of the polynomial, n - 3 multiplies are
     *  saved and 4 * n - 6 additions are saved.
     */
    template<typename _Tp2>
      auto
      operator()(std::complex<_Tp2> __z) const
      -> decltype(value_type{} * std::complex<_Tp2>{})
      {
	const auto __r = _Tp{2} * std::real(__z);
	const auto __s = std::norm(__z);
	auto __aa = this->coefficient(this->degree());
	auto __bb = this->coefficient(this->degree() - 1);
	for (int __j = 1; __j <= this->degree(); ++__j)
	  {
	    auto __cc  = __s * __aa;
	    __aa = __bb + __r * __aa;
	    __bb = this->coefficient(this->degree() - __j) - __cc;
	  }
	return __aa * __z + __bb;
      };

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
      eval(value_type __x, std::array<value_type, N>& __arr)
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
	    for (size_t __i = 2; __i < __arr.size(); ++__i)
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
    template<typename OutIter>
      void
      eval(value_type __x, OutIter __b, OutIter __e)
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
     *  Return the derivative of the polynomial.
     */
    _Polynomial
    derivative() const
    {
      _Polynomial __res(value_type{}, (this->degree() > 0UL ? this->degree() - 1 : 0UL));
      for (size_type __i = 1; __i <= this->degree(); ++__i)
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
      for (size_type __i = 0; __i <= this->degree(); ++__i)
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
	this->degree(std::max(this->degree(), __poly.degree())); // Resize if necessary.
	for (size_type __i = 0; __i <= __poly.degree(); ++__i)
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
	for (size_type __i = 0; __i <= __poly.degree(); ++__i)
	  this->_M_coeff[__i] -= static_cast<value_type>(__poly._M_coeff[__i]);
	return *this;
      }

    /**
     *  Multiply the polynomial by another polynomial.
     */
    template<typename _Up>
      _Polynomial&
      operator*=(const _Polynomial<_Up>& __poly)
      {
	//  Test for zero size polys and do special processing?
	std::vector<value_type> __new_coeff(this->degree() + __poly.degree() + 1);
	for (size_type __i = 0; __i < this->_M_coeff.size(); ++__i)
	  for (size_type __j = 0; __j < __poly._M_coeff.size(); ++__j)
	    __new_coeff[__i + __j] += this->_M_coeff[__i]
				* static_cast<value_type>(__poly._M_coeff[__j]);
	this->_M_coeff = __new_coeff;
	return *this;
      }

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

    value_type
    coefficient(size_type __i) const
    { return (this->_M_coeff.size() > __i ? this->_M_coeff[__i] : value_type{}); }

    void
    coefficient(size_type __i, value_type __val)
    { this->_M_coeff.at(__i) = __val; }

    value_type
    operator[](size_type __i) const
    { return this->_M_coeff[__i]; }

    reference
    operator[](size_type __i)
    { return this->_M_coeff[__i]; }

    iterator
    begin()
    { return this->_M_coeff.begin(); }

    iterator
    end()
    { return this->_M_coeff.end(); }

    const_iterator
    begin() const
    { return this->_M_coeff.begin(); }

    const_iterator
    end() const
    { return this->_M_coeff.end(); }

    const_iterator
    cbegin() const
    { return this->_M_coeff.cbegin(); }

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
      operator==(const _Polynomial<_Tp1>& __pa, const _Polynomial<_Tp1>& __pb);

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
    operator<<(std::basic_ostream<CharT, Traits>& __os, const _Polynomial<_Tp>& __poly)
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
    operator>>(std::basic_istream<CharT, Traits>& __is, _Polynomial<_Tp>& __poly)
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

} // namespace __gnu_cxx

/*
//  Here is another idea from boost.
namespace __gnu_cxx
{

  template<typename _Tp, size_t _Num, size_t _K>
    constexpr _Tp
    _Polynomial_eval_help(const std::array<_Tp, _Num>& __arr, _Tp __x)
    { return _Polynomial_eval_help<_Tp, _Num, _K - 1>(__arr, __x) * __x + __arr[_K]; }

  template<typename _Tp, size_t _Num>
    constexpr _Tp
    _Polynomial_eval_help<_Tp, _Num, 0>(const std::array<_Tp, _Num>& __arr, _Tp __x)
    { return _Tp(1); }

  template<typename _Tp, size_t _Num>
    constexpr _Tp
    _Polynomial_eval(const std::array<_Tp, _Num>& __arr, _Tp __x)
    { return _Polynomial_eval_help<_Tp, _Num, _Num - 1>(__arr, __x); }

  template<typename _Tp>
    constexpr _Tp
    _Polynomial_eval<_Tp, 0>(const std::array<_Tp, 0>& __arr, _Tp __x)
    { return _Tp(0); }

} // namespace __gnu_cxx
*/
