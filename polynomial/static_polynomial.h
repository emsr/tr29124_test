
#include <limits>
//#include <experimental/array> // For make_array/to_array.
#include <array>
#include <complex>

namespace __gnu_cxx
{

  /**
   * This is a constant size polynomial.
   * It is really meant to just evaluate canned polynomial literals.
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
	    this->_M_coeff[__i] = static_cast<value_type>(__poly._M_coeff[__i]);
	}

      //  Constructor from C-type array.
      template<typename _Up>
	constexpr
	_StaticPolynomial(const _Up (&__arr)[_Num])
	: _M_coeff{}
	{
          for (auto __i = 0ULL; __i < _Num; ++__i)
	    this->_M_coeff[__i] = static_cast<value_type>(__arr[__i]);
	}

      //  Constructor from initializer_list array.
      constexpr
      _StaticPolynomial(std::initializer_list<_Tp> __il)
      : _M_coeff{}
      {
	//static_assert(__il.size() == _Num, "");
	std::size_t __i = 0;
	for (auto&& __coeff : __il)
	  this->_M_coeff[__i++] = __coeff;
      }

      /**
       *  Create a polynomial - actually a monomial - of just one term.
       */
      constexpr explicit
      _StaticPolynomial(value_type __a, size_type __degree = 0)
      : _M_coeff(__degree + 1)
      {
        static_assert(__degree < _Num, "_StaticPolynomial: degree out of range");
        this->_M_coeff[__degree] = __a;
      }

      /**
       *  Create a polynomial from an argument list of coefficients.
      constexpr
      _StaticPolynomial(value_type&& __aa0, value_type&&... __aa)
      : _M_coeff(std::experimental::make_array(__aa0, __aa...))
      { }
       */

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
	    return __poly;
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
	      return __poly;
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
       *  Evaluate the even part of the polynomial at the input point.
       */
      constexpr value_type
      even(value_type __x) const
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
      constexpr value_type
      odd(value_type __x) const
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
      template<typename _Tp2>
	auto
	even(std::complex<_Tp2> __z) const
	-> decltype(value_type{} * std::complex<_Tp2>{})
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
	    return decltype(value_type{} * std::complex<_Tp2>{}){};
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
      template<typename _Tp2>
	auto
	odd(std::complex<_Tp2> __z) const
	-> decltype(value_type{} * std::complex<_Tp2>{})
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
	    return decltype(value_type{} * std::complex<_Tp2>{}){};
	};

      /**
       *  Return the derivative of the polynomial.
       */
      std::enable_if_t<0 < _Num, _StaticPolynomial<_Tp, _Num - 1>>
      derivative() const
      {
	_StaticPolynomial<_Tp, _Num - 1> __res;
	for (size_type __i = 1; __i <= this->degree(); ++__i)
	  __res._M_coeff[__i - 1] = __i * _M_coeff[__i];
	return __res;
      }

      /**
       *  Return the integral of the polynomial with given integration constant.
       */
      _StaticPolynomial<_Tp, _Num + 1>
      integral(value_type __c = value_type{}) const
      {
	_StaticPolynomial<_Tp, _Num + 1> __res;
	__res._M_coeff[0] = __c;
	for (size_type __i = 0; __i <= this->degree(); ++__i)
	  __res._M_coeff[__i + 1] = _M_coeff[__i] / value_type(__i + 1);
	return __res;
      }

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
  template<typename _Tp, std::size_t _NumA, std::size_t _NumB>
    inline bool
    operator!=(const _StaticPolynomial<_Tp, _NumA>& __pa,
	       const _StaticPolynomial<_Tp, _NumB>& __pb)
    { return true; }

  /**
   *  Return false if two polynomials are equal.
   */
  template<typename _Tp, std::size_t _Num>
    inline bool
    operator!=(const _StaticPolynomial<_Tp, _Num>& __pa,
	       const _StaticPolynomial<_Tp, _Num>& __pb)
    { return !(__pa == __pb); }

} // namespace __gnu_cxx
