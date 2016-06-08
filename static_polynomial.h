
#include <initializer_list>
#include <vector>
#include <iosfwd>
#include <limits>
#include <array>
#include <ios>
#include <cmath>
#include <complex>

/**
 *  detail: Do we want this to always have a size of at least one? a_0 = _Tp()?  YES.
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
 *    std::pair<> div(const _StaticPolynomial& __a, const _StaticPolynomial& __b)
 *    void divmod(const _StaticPolynomial& __a, const _StaticPolynomial& __b,
 *                _StaticPolynomial& __q, _StaticPolynomial& __r);
 *  Should factory methods like derivative and integral be globals?
 *  I could have members:
 *    _StaticPolynomial& integrate(_Tp c);
 *    _StaticPolynomial& differentiate();
 */
namespace __gnu_cxx
{

  /**
   *
   */
  template<typename _Tp, std::size_t _Num>
    class _StaticPolynomial
    {
    public:
      /**
       *  Typedefs.
       *  @todo Should we grab these from _M_coeff (i.e. std::array<_Tp, _Num>)?
       */
      using value_type = typename std::array<_Tp, _Num>::value_type;
      using reference = typename std::array<_Tp, _Num>::reference;
      using const_reference = typename std::array<_Tp, _Num>::const_reference;
      using pointer = typename std::array<_Tp, _Num>::pointer;
      using const_pointer = typename std::array<_Tp, _Num>::const_pointer;
      using iterator = typename std::array<value_type, _Num>::iterator;
      using const_iterator = typename std::array<value_type, _Num>::const_iterator;
      using reverse_iterator = typename std::array<value_type, _Num>::reverse_iterator;
      using const_reverse_iterator = typename std::array<value_type, _Num>::const_reverse_iterator;
      using size_type = typename std::array<_Tp, _Num>::size_type;
      using difference_type = typename std::array<_Tp, _Num>::difference_type;

      /**
       *  Create a zero degree polynomial with value zero.
       */
      constexpr
      _StaticPolynomial()
      : _M_coeff{}
      { }

      /**
       *  Copy ctor.
       */
      constexpr _StaticPolynomial(const _StaticPolynomial&) = default;

      template<typename _Up>
	constexpr
	_StaticPolynomial(const _StaticPolynomial<_Up, _Num>& __poly)    
	: _M_coeff{}
	{
          for (auto __i = 0ULL; __i < _Num; ++__i)
	    this->_M_coeff[__i] = static_cast<value_type>(__poly._M_coeff[__i]));
	}

      /**
       *  Create a polynomial of just one constant term.
       */
      constexpr explicit
      _StaticPolynomial(value_type __a, size_type __degree = 0)
      : _M_coeff(__degree + 1)
      { this->_M_coeff[__degree] = __a; }

      /**
       *  Create a polynomial from an initializer list of coefficients.
       */
      constexpr
      _StaticPolynomial(std::initializer_list<value_type> __ila)
      : _M_coeff(__ila)
      { }

      /**
       *  Create a polynomial from an input iterator range of coefficients.
       */
      template<typename InIter,
	       typename = std::_RequireInputIter<InIter>>
	constexpr
	_StaticPolynomial(const InIter& __abegin, const InIter& __aend)
	: _M_coeff(__abegin, __aend)
	{ }

      /**
       *  Swap the polynomial with another polynomial.
       */
      void
      swap(_StaticPolynomial& __poly)
      { this->_M_coeff.swap(__poly._M_coeff); }

      /**
       *  Evaluate the polynomial at the input point.
       */
      constexpr value_type
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
	constexpr auto
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
	constexpr auto
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
	constexpr OutIter
	operator()(const InIter& __xbegin, const InIter& __xend,
        	   OutIter& __pbegin) const
	{
	  for (; __xbegin != __xend; ++__xbegin)
	    __pbegin++ = (*this)(__xbegin++);
	  return __pbegin;
	}

      //  Could/should this be done by output iterator range?
      template<size_type N>
	constexpr void
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
	constexpr void
	eval(value_type __x, OutIter __b, OutIter __e)
	{
	  if(__b != __e)
	    {
	      std::fill(__b, __e, value_type{});
	      constexpr size_type __sz = _M_coeff.size();
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
      _StaticPolynomial
      derivative() const
      {
	_StaticPolynomial __res(value_type{}, (this->degree() > 0UL ? this->degree() - 1 : 0UL));
	for (size_type __i = 1; __i <= this->degree(); ++__i)
	  __res._M_coeff[__i - 1] = __i * _M_coeff[__i];
	return __res;
      }

      /**
       *  Return the integral of the polynomial with given integration constant.
       */
      _StaticPolynomial
      integral(value_type __c = value_type{}) const
      {
	_StaticPolynomial __res(value_type{}, this->degree() + 1);
	__res._M_coeff[0] = __c;
	for (size_type __i = 0; __i <= this->degree(); ++__i)
	  __res._M_coeff[__i + 1] = _M_coeff[__i] / value_type(__i + 1);
	return __res;
      }
 
      /**
       *  Unary plus.
       */
      _StaticPolynomial&
      operator+() const
      { return *this; }

      /**
       *  Unary minus.
       */
      _StaticPolynomial&
      operator-() const
      { return _StaticPolynomial(*this) *= value_type(-1); }

      /**
       *  Copy assignment.
       */
      _StaticPolynomial&
      operator=(const _StaticPolynomial&) = default;

      template<typename _Up>
	_StaticPolynomial&
	operator=(const _StaticPolynomial<_Up, _Num>& __poly)
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
      constexpr _StaticPolynomial&
      operator=(std::initializer_list<value_type> __ila)
      {
	for (size_type __i = 0;
	     __i <= std::min(this->degree(), __ila.size()); ++__i)
	  this->_M_coeff[__i] = __ila[__i];
	return *this;
      }

      /**
       *  Add a scalar to the polynomial.
       */
      template<typename _Up>
	constexpr _StaticPolynomial&
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
	constexpr _StaticPolynomial&
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
	constexpr _StaticPolynomial&
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
	constexpr _StaticPolynomial&
	operator/=(const _Up& __x)
	{
	  for (size_type __i = 0; __i < _M_coeff.size(); ++__i)
	    this->_M_coeff[__i] /= static_cast<value_type>(__x);
	  return *this;
	}

      /**
       *  Add another polynomial to the polynomial.
       */
      template<typename _Up>
	_StaticPolynomial&
	operator+=(const _StaticPolynomial<_Up, _Num>& __poly)
	{
	  for (size_type __i = 0; __i <= __poly.degree(); ++__i)
	    this->_M_coeff[__i] += static_cast<value_type>(__poly._M_coeff[__i]);
	  return *this;
	}

      /**
       *  Subtract another polynomial from the polynomial.
       */
      template<typename _Up>
	_StaticPolynomial&
	operator-=(const _StaticPolynomial<_Up, _Num>& __poly)
	{
	  for (size_type __i = 0; __i <= __poly.degree(); ++__i)
	    this->_M_coeff[__i] -= static_cast<value_type>(__poly._M_coeff[__i]);
	  return *this;
	}

      /**
       *  Return the degree or the power of the largest coefficient.
       */
      constexpr size_type
      degree() const
      { return (this->_M_coeff.size() > 0 ? this->_M_coeff.size() - 1 : 0); }

      constexpr value_type
      coefficient(size_type __i) const
      { return (this->_M_coeff.size() > __i ? this->_M_coeff[__i] : value_type{}); }

      void
      coefficient(size_type __i, value_type __val)
      { this->_M_coeff.at(__i) = __val; }

      constexpr value_type
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

      template<typename _Tp1>
	friend bool
	operator==(const _StaticPolynomial<_Tp1, _Num>& __pa,
		   const _StaticPolynomial<_Tp1, _Num>& __pb);

    private:

      std::array<value_type, _Num> _M_coeff;
    };

  /**
   *  Return the sum of a polynomial with a scalar.
   */
  template<typename _Tp, std::size_t _Num, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() + _Up()), _Num>
    operator+(const _StaticPolynomial<_Tp, _Num>& __poly, const _Up& __x)
    { return _StaticPolynomial<decltype(_Tp() + _Up()), _Num>(__poly) += __x; }

  /**
   *  Return the difference of a polynomial with a scalar.
   */
  template<typename _Tp, std::size_t _Num, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() - _Up()), _Num>
    operator-(const _StaticPolynomial<_Tp, _Num>& __poly, const _Up& __x)
    { return _StaticPolynomial<decltype(_Tp() - _Up())>(__poly) -= __x; }

  /**
   *  Return the product of a polynomial with a scalar.
   */
  template<typename _Tp, std::size_t _Num, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() * _Up()), _Num>
    operator*(const _StaticPolynomial<_Tp, _Num>& __poly, const _Up& __x)
    { return _StaticPolynomial<decltype(_Tp() * _Up())>(__poly) *= __x; }

  /**
   *  Return the quotient of a polynomial with a scalar.
   */
  template<typename _Tp, std::size_t _Num, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() / _Up()), _Num>
    operator/(const _StaticPolynomial<_Tp, _Num>& __poly, const _Up& __x)
    { return _StaticPolynomial<decltype(_Tp() / _Up())>(__poly) /= __x; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() / _Up())>
    operator%(const _StaticPolynomial<_Tp>& __poly, const _Up& __x)
    { return _StaticPolynomial<decltype(_Tp() / _Up())>(__poly) %= __x; }

  /**
   *  Return the sum of two polynomials.
   */
  template<typename _Tp, std::size_t _NumT, typename _Up, std::size_t _NumU>
    inline _StaticPolynomial<decltype(_Tp() + _Up()), std::max(_NumT, _NumU)>
    operator+(const _StaticPolynomial<_Tp, _NumT>& __pa,
	      const _StaticPolynomial<_Up, _NumU>& __pb)
    {
      return _StaticPolynomial<decltype(_Tp() + _Up()),
			       std::max(_NumT, _NumU)>(__pa) += __pb;
    }

  /**
   *  Return the difference of two polynomials.
   */
  template<typename _Tp, std::size_t _NumT, typename _Up, std::size_t _NumU>
    inline _StaticPolynomial<decltype(_Tp() - _Up()), std::max(_NumT, _NumU)>
    operator-(const _StaticPolynomial<_Tp, _NumT>& __pa,
	      const _StaticPolynomial<_Up, _NumU>& __pb)
    {
      return _StaticPolynomial<decltype(_Tp() - _Up()),
			       std::max(_NumT, _NumU)>(__pa) -= __pb;
    }

  /**
   *  Return the product of two polynomials.
   */
  template<typename _Tp, std::size_t _NumT, typename _Up, std::size_t _NumU>
    inline _StaticPolynomial<decltype(_Tp() * _Up()), _NumT * _NumU>
    operator*(const _StaticPolynomial<_Tp, _NumT>& __pa,
	      const _StaticPolynomial<_Up, _NumU>& __pb)
    {
      return _StaticPolynomial<decltype(_Tp() * _Up(),
			       _NumT * _NumU)>(__pa) *= __pb;
    }

  /**
   *  Return the quotient of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() / _Up())>
    operator/(const _StaticPolynomial<_Tp>& __pa, const _StaticPolynomial<_Up>& __pb)
    { return _StaticPolynomial<decltype(_Tp() / _Up())>(__pa) /= __pb; }

  /**
   *  Return the modulus or remainder of one polynomial relative to another one.
   */
  template<typename _Tp, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() / _Up())>
    operator%(const _StaticPolynomial<_Tp>& __pa, const _StaticPolynomial<_Up>& __pb)
    { return _StaticPolynomial<decltype(_Tp() / _Up())>(__pa) %= __pb; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() + _Up())>
    operator+(const _Tp& __x, const _StaticPolynomial<_Up>& __poly)
    { return _StaticPolynomial<decltype(_Tp() + _Up())>(__x) += __poly; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() - _Up())>
    operator-(const _Tp& __x, const _StaticPolynomial<_Up>& __poly)
    { return _StaticPolynomial<decltype(_Tp() - _Up())>(__x) -= __poly; }

  /**
   *
   */
  template<typename _Tp, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() * _Up())>
    operator*(const _Tp& __x, const _StaticPolynomial<_Up>& __poly)
    { return _StaticPolynomial<decltype(_Tp() * _Up())>(__x) *= __poly; }

  /**
   *  Return the quotient of two polynomials.
   */
  template<typename _Tp, typename _Up>
    inline _StaticPolynomial<decltype(_Tp() / _Up())>
    operator/(const _Tp& __x, const _StaticPolynomial<_Up>& __poly)
    { return _StaticPolynomial<decltype(_Tp() / _Up())>(__x) /= __poly; }

  /**
   *  Return true if two polynomials are equal.
   */
  template<typename _Tp, std::size_t _NumA, std::size_t _NumB>
    inline bool
    operator==(const _StaticPolynomial<_Tp, _NumA>&,
	       const _StaticPolynomial<_Tp, _NumB>&)
    { return false; }

  template<typename _Tp, std::size_t _Num>
    inline bool
    operator==(const _StaticPolynomial<_Tp, _Num>& __pa,
	       const _StaticPolynomial<_Tp, _Num>& __pb)
    { return __pa._M_coeff == __pb._M_coeff; }

  /**
   *  Return false if two polynomials are equal.
   */
  template<typename _Tp>
    inline bool
    operator!=(const _StaticPolynomial<_Tp>& __pa,
	       const _StaticPolynomial<_Tp>& __pb)
    { return !(__pa == __pb); }

} // namespace __gnu_cxx
