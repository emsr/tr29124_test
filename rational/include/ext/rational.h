
#ifndef _RATIONAL_HPP
#define _RATIONAL_HPP 1

#include <iostream>
#include <ios>
#include <stdexcept>
#include <limits>
#include <ratio>
#include <cassert>
#include <experimental/numeric>

// operator%=, operator%? = Meaningless.
// std::hash... We should do that.
// floor/ceil... Like what went into chrono.

namespace __gnu_cxx
{

  class _Bad_Rational
  : public std::domain_error
  {
  public:
    explicit
    _Bad_Rational()
    : std::domain_error("_Bad_Rational: zero denominator")
    { }
  };

  template<typename _IntTp>
    class _Rational;

  template<typename _IntTp>
    _Rational<_IntTp>
    abs(const _Rational<_IntTp>&);

  template<typename _IntTp = std::intmax_t>
    class _Rational
    {
      static_assert(std::numeric_limits<_IntTp>::is_specialized, "");

    public:
      typedef _IntTp value_type;

      constexpr
      _Rational()
      : _M_num{0}, _M_den{1}
      { }

      constexpr explicit
      _Rational(value_type __num)
      : _M_num(__num), _M_den(1)
      { }

      template<intmax_t _Num, intmax_t _Den>
	constexpr explicit
	_Rational(std::ratio<_Num, _Den>)
	: _M_num(std::ratio<_Num, _Den>::num),
	  _M_den(std::ratio<_Num, _Den>::den)
	{ }

      _Rational(value_type __num, value_type __den)
      : _M_num(__num), _M_den(__den)
      { this->_M_normalize(); }

      // Default copy constructor and assignment are fine

      // Add assignment from value_type
      _Rational&
      operator=(value_type __num)
      { return this->assign(__num, 1); }

      // Assign in place
      _Rational&
      assign(value_type, value_type);

      // Access to representation
      constexpr value_type
      num() const
      { return this->_M_num; }

      constexpr value_type
      den() const
      { return this->_M_den; }

      // Arithmetic assignment operators
      _Rational& operator+=(_Rational);
      _Rational& operator-=(_Rational);
      _Rational& operator*=(_Rational);
      _Rational& operator/=(_Rational);

      _Rational& operator+=(value_type);
      _Rational& operator-=(value_type);
      _Rational& operator*=(value_type);
      _Rational& operator/=(value_type);

      // Increment and decrement
      const _Rational& operator++();
      const _Rational& operator--();

      // Operator not
      bool
      operator!() const
      { return !this->_M_num; }

      // Boolean conversion
      operator
      bool() const
      { return this->operator!(); }

      // Comparison operators
      bool operator<(const _Rational&) const;
      bool operator==(const _Rational&) const;

      bool operator<(value_type) const;
      bool operator>(value_type) const;
      bool operator==(value_type) const;

    private:

      value_type _M_num;
      value_type _M_den;

      // Representation note: Fractions are kept in normalized form at all
      // times. normalized form is defined as gcd(num,den) == 1 and den > 0.
      // In particular, note that the implementation of abs() below relies
      // on den always being positive.
      bool _M_valid() const;
      void _M_normalize();
    };

  // Assign in place
  template<typename _IntTp>
    inline _Rational<_IntTp>&
    _Rational<_IntTp>::assign(value_type __num, value_type __den)
    {
      this->_M_num = __num;
      this->_M_den = __den;
      this->_M_normalize();
      return *this;
    }

  // Unary plus and minus
  template<typename _IntTp>
    inline _Rational<_IntTp>
    operator+(const _Rational<_IntTp>& __rat)
    { return __rat; }

  template<typename _IntTp>
    inline _Rational<_IntTp>
    operator-(const _Rational<_IntTp>& __rat)
    { return _Rational<_IntTp>(-__rat.num(), __rat.den()); }

  // Arithmetic assignment operators
  template<typename _IntTp>
    _Rational<_IntTp>&
    _Rational<_IntTp>::operator+=(_Rational<_IntTp> __rat)
    {
      // This calculation avoids overflow, and minimises the number of expensive
      // calculations. Thanks to Nickolay Mladenov for this algorithm.
      //
      // Proof:
      // We have to compute a/b + c/d, where gcd(a,b)=1 and gcd(b,c)=1.
      // Let g = gcd(b,d), and b = b1*g, d=d1*g. Then gcd(b1,d1)=1
      //
      // The result is (a*d1 + c*b1) / (b1*d1*g).
      // Now we have to normalize this ratio.
      // Let's assume h | gcd((a*d1 + c*b1), (b1*d1*g)), and h > 1
      // If h | b1 then gcd(h,d1)=1 and hence h|(a*d1+c*b1) => h|a.
      // But since gcd(a,b1)=1 we have h=1.
      // Similarly h|d1 leads to h=1.
      // So we have that h | gcd((a*d1 + c*b1) , (b1*d1*g)) => h|g
      // Finally we have gcd((a*d1 + c*b1), (b1*d1*g)) = gcd((a*d1 + c*b1), g)
      // Which proves that instead of normalizing the result, it is better to
      // divide num and den by gcd((a*d1 + c*b1), g)

      value_type __gcd = std::experimental::gcd(this->_M_den, __rat.den());
      this->_M_den /= __gcd;  // = b1 from the calculations above
      this->_M_num = this->_M_num * (__rat.den() / __gcd) + __rat.num() * this->_M_den;
      __gcd = std::experimental::gcd(this->_M_num, __gcd);
      this->_M_num /= __gcd;
      this->_M_den *= __rat.den() / __gcd;

      return *this;
    }

  template<typename _IntTp>
    _Rational<_IntTp>&
    _Rational<_IntTp>::operator-=(_Rational<_IntTp> __rat)
    {
      // This calculation avoids overflow, and minimises the number of expensive
      // calculations. It corresponds exactly to the += case above
      value_type __gcd = std::experimental::gcd(this->_M_den, __rat.den());
      this->_M_den /= __gcd;
      this->_M_num = this->_M_num * (__rat.den() / __gcd) - __rat.num() * this->_M_den;
      __gcd = std::experimental::gcd(this->_M_num, __gcd);
      this->_M_num /= __gcd;
      this->_M_den *= __rat.den() / __gcd;

      return *this;
    }

  template<typename _IntTp>
    _Rational<_IntTp>&
    _Rational<_IntTp>::operator*=(_Rational<_IntTp> __rat)
    {
      // Avoid overflow and preserve normalization
      value_type __gcd1 = std::experimental::gcd(this->_M_num, __rat.den());
      value_type __gcd2 = std::experimental::gcd(__rat.num(), this->_M_den);
      this->_M_num = (this->_M_num / __gcd1) * (__rat.num() / __gcd2);
      this->_M_den = (this->_M_den / __gcd2) * (__rat.den() / __gcd1);
      return *this;
    }

  template<typename _IntTp>
    _Rational<_IntTp>&
    _Rational<_IntTp>::operator/=(_Rational<_IntTp> __rat)
    {
      constexpr value_type zero{0};

      // Trap division by zero
      if (__rat.num() == zero)
	throw _Bad_Rational();
      if (this->_M_num == zero)
	return *this;

      // Avoid overflow and preserve normalization
      value_type __gcd1 = std::experimental::gcd(this->_M_num, __rat.num());
      value_type __gcd2 = std::experimental::gcd(__rat.den(), this->_M_den);
      this->_M_num = (this->_M_num / __gcd1) * (__rat.den() / __gcd2);
      this->_M_den = (this->_M_den / __gcd2) * (__rat.num() / __gcd1);

      if (this->_M_den < zero)
	{
          this->_M_num = -this->_M_num;
          this->_M_den = -this->_M_den;
	}
      return *this;
    }

  // Mixed-mode operators
  template<typename _IntTp>
    inline _Rational<_IntTp>&
    _Rational<_IntTp>::operator+=(value_type __i)
    { return this->operator+=(_Rational<_IntTp>(__i)); }

  template<typename _IntTp>
    inline _Rational<_IntTp>&
    _Rational<_IntTp>::operator-=(value_type __i)
    { return this->operator-=(_Rational<_IntTp>(__i)); }

  template<typename _IntTp>
    inline _Rational<_IntTp>&
    _Rational<_IntTp>::operator*=(value_type __i)
    { return this->operator*=(_Rational<_IntTp>(__i)); }

  template<typename _IntTp>
    inline _Rational<_IntTp>&
    _Rational<_IntTp>::operator/=(value_type __i)
    { return this->operator/=(_Rational<_IntTp>(__i)); }

  // Increment and decrement
  template<typename _IntTp>
    inline const _Rational<_IntTp>&
    _Rational<_IntTp>::operator++()
    {
      // This can never denormalise the fraction
      this->_M_num += this->_M_den;
      return *this;
    }

  template<typename _IntTp>
    inline const _Rational<_IntTp>&
    _Rational<_IntTp>::operator--()
    {
      // This can never denormalise the fraction
      this->_M_num -= this->_M_den;
      return *this;
    }

  template<typename _IntTp>
    bool
    _Rational<_IntTp>::operator<(const _Rational<_IntTp>& r) const
    {
      constexpr value_type zero{0};

      assert(this->_M_den > zero);
      assert(r.den() > zero);

      // Determine relative order by expanding each value to its simple continued
      // fraction representation using the Euclidian GCD algorithm.
      struct
      { value_type n, d, q, r; }
      ts{
	this->_M_num, this->_M_den,
	static_cast<value_type>(this->_M_num / this->_M_den),
	static_cast<value_type>(this->_M_num % this->_M_den)
      },
      rs
      {
	r.num(), r.den(),
	static_cast<value_type>(r.num() / r.den()),
	static_cast<value_type>(r.num() % r.den())
      };
      unsigned reverse = 0u;

      // Normalize negative moduli by repeatedly adding the (positive) denominator
      // and decrementing the quotient.  Later cycles should have all positive
      // values, so this only has to be done for the first cycle.  (The rules of
      // C++ require a nonnegative quotient & remainder for a nonnegative dividend
      // & positive divisor.)
      while (ts.r < zero)
	{
	  ts.r += ts.d;
	  --ts.q;
	}
      while (rs.r < zero)
	{
	  rs.r += rs.d;
	  --rs.q;
	}

      // Loop through and compare each variable's continued-fraction components
      while (true)
	{
	  // The quotients of the current cycle are the continued-fraction
	  // components.  Comparing two c.f. is comparing their sequences,
	  // stopping at the first difference.
	  if (ts.q != rs.q)
	    {
	      // Since reciprocation changes the relative order of two variables,
	      // and c.f. use reciprocals, the less/greater-than test reverses
	      // after each index.  (Start w/ non-reversed @ whole-number place.)
	      return reverse ? ts.q > rs.q : ts.q < rs.q;
	    }

	  // Prepare the next cycle
	  reverse ^= 1u;

	  if ((ts.r == zero) || (rs.r == zero))
            {
              // At least one variable's c.f. expansion has ended
              break;
            }

	  ts.n = ts.d;         ts.d = ts.r;
	  ts.q = ts.n / ts.d;  ts.r = ts.n % ts.d;
	  rs.n = rs.d;         rs.d = rs.r;
	  rs.q = rs.n / rs.d;  rs.r = rs.n % rs.d;
	}

      // Compare infinity-valued components for otherwise equal sequences
      if (ts.r == rs.r)
	{
          // Both remainders are zero, so the next (and subsequent) c.f.
          // components for both sequences are infinity.  Therefore, the sequences
          // and their corresponding values are equal.
          return false;
	}
      else
	{
          // Exactly one of the remainders is zero, so all following c.f.
          // components of that variable are infinity, while the other variable
          // has a finite next c.f. component.  So that other variable has the
          // lesser value (modulo the reversal flag!).
          return (ts.r != zero) != static_cast<bool>(reverse);
	}
    }

  template<typename _IntTp>
    bool
    _Rational<_IntTp>::operator<(value_type __i) const
    {
      constexpr value_type zero{0};

      // Break value into mixed-fraction form, w/ always-nonnegative remainder
      assert(this->_M_den > zero);
      value_type __quo = this->_M_num / this->_M_den,
		 __rem = this->_M_num % this->_M_den;
      while (__rem < zero)
	{
	  __rem += this->_M_den;
	  --__quo;
	}

      // Compare with just the quotient, since the remainder always bumps the
      // value up.  [Since q = floor(n/d), and if n/d < i then q < i, if n/d == i
      // then q == i, if n/d == i + r/d then q == i, and if n/d >= i + 1 then
      // q >= i + 1 > i; therefore n/d < i iff q < i.]
      return __quo < __i;
    }

  template<typename _IntTp>
    bool
    _Rational<_IntTp>::operator>(value_type __i) const
    {
      // Trap equality first
      if (this->_M_num == __i && this->_M_den == _IntTp(1))
	return false;

      // Otherwise, we can use operator<
      return !this->operator<(__i);
    }

  template<typename _IntTp>
    inline bool
    _Rational<_IntTp>::operator==(const _Rational<_IntTp>& __rat) const
    { return ((this->_M_num == __rat.num) && (this->_M_den == __rat.den())); }

  template<typename _IntTp>
    inline bool
    _Rational<_IntTp>::operator==(value_type __i) const
    { return ((this->_M_den == _IntTp{1}) && (this->_M_num == __i)); }

  // Invariant check
  template<typename _IntTp>
    inline bool
    _Rational<_IntTp>::_M_valid() const
    {
      return this->_M_den > value_type{0}
          && std::experimental::gcd(this->_M_num, this->_M_den) == value_type{1};
    }

  // Normalisation
  template<typename _IntTp>
    void
    _Rational<_IntTp>::_M_normalize()
    {
      constexpr value_type zero{0};

      if (this->_M_den == zero)
	throw _Bad_Rational();

      // Handle the case of zero separately, to avoid division by zero
      if (this->_M_num == zero)
	{
          this->_M_den = value_type{1};
          return;
	}

      value_type __gcd = std::experimental::gcd(this->_M_num, this->_M_den);

      this->_M_num /= __gcd;
      this->_M_den /= __gcd;

      // Ensure that the denominator is positive
      if (this->_M_den < zero)
	{
          this->_M_num = -this->_M_num;
          this->_M_den = -this->_M_den;
	}

      assert(this->_M_valid());
    }

  // Type conversion
  template<typename _Tp, typename _IntTp>
    inline _Tp
    _Rational_cast(const _Rational<_IntTp>& __src)
    {
      return static_cast<_Tp>(__src.num())
           / static_cast<_Tp>(__src.den());
    }

  // Symmetric math operators.
  template<typename _IntTp>
    _Rational<_IntTp>
    operator+(const _Rational<_IntTp>& __r1, const _Rational<_IntTp>& __r2)
    { return _Rational<_IntTp>(__r1) += __r2; }

  template<typename _IntTp, typename _RealTp>
    _RealTp
    operator+(_RealTp __r1, const _Rational<_IntTp>& __r2)
    { return __r1 + _Rational_cast<_RealTp>(__r2); }

  template<typename _IntTp, typename _RealTp>
    _RealTp
    operator+(const _Rational<_IntTp>& __r1, _RealTp __r2)
    { return _Rational_cast<_RealTp>(__r1) + __r2; }

  template<typename _IntTp>
    _Rational<_IntTp>
    operator-(const _Rational<_IntTp>& __r1, const _Rational<_IntTp>& __r2)
    { return _Rational<_IntTp>(__r1) -= __r2; }

  template<typename _IntTp, typename _RealTp>
    _RealTp
    operator-(_RealTp __r1, const _Rational<_IntTp>& __r2)
    { return __r1 - _Rational_cast<_RealTp>(__r2); }

  template<typename _IntTp, typename _RealTp>
    _RealTp
    operator-(const _Rational<_IntTp>& __r1, _RealTp __r2)
    { return _Rational_cast<_RealTp>(__r1) - __r2; }

  template<typename _IntTp>
    _Rational<_IntTp>
    operator*(const _Rational<_IntTp>& __r1, const _Rational<_IntTp>& __r2)
    { return _Rational<_IntTp>(__r1) *= __r2; }

  template<typename _IntTp, typename _RealTp>
    _RealTp
    operator*(_RealTp __r1, const _Rational<_IntTp>& __r2)
    { return __r1 * _Rational_cast<_RealTp>(__r2); }

  template<typename _IntTp, typename _RealTp>
    _RealTp
    operator*(const _Rational<_IntTp>& __r1, _RealTp __r2)
    { return _Rational_cast<_RealTp>(__r1) * __r2; }

  template<typename _IntTp>
    _Rational<_IntTp>
    operator/(const _Rational<_IntTp>& __r1, const _Rational<_IntTp>& __r2)
    { return _Rational<_IntTp>(__r1) /= __r2; }

  template<typename _IntTp, typename _RealTp>
    _RealTp
    operator/(_RealTp __r1, const _Rational<_IntTp>& __r2)
    { return __r1 / _Rational_cast<_RealTp>(__r2); }

  template<typename _IntTp, typename _RealTp>
    _RealTp
    operator/(const _Rational<_IntTp>& __r1, _RealTp __r2)
    { return _Rational_cast<_RealTp>(__r1) / __r2; }

  namespace detail
  {

    // A utility class to reset the format flags for an istream at end
    // of scope, even in case of exceptions
    struct resetter
    {
      resetter(std::istream& __is)
      : _M_is(__is), _M_fl(__is.flags())
      { }
      ~resetter()
      { _M_is.flags(_M_fl); }

      std::istream& _M_is;
      std::ios_base::fmtflags _M_fl;
    };

  }

  // Input and output
  template<typename _IntTp>
    std::istream&
    operator>>(std::istream& __is, _Rational<_IntTp>& __rat)
    {
      using __value_type = typename _Rational<_IntTp>::value_type;

      auto __num = __value_type{0}, __den = __value_type{1};
      char __c = 0;
      detail::resetter sentry(__is);

      __is >> __num;
      __c = __is.get();

      if (__c != '/')
	__is.clear(std::ios_base::badbit);  // old GNU c++ lib has no ios_base

  #if !defined(__GNUC__) || (defined(__GNUC__) && (__GNUC__ >= 3))
      __is >> std::noskipws;
  #else
      __is.unsetf(ios::skipws); // compiles, but seems to have no effect.
  #endif
      __is >> __den;

      if (__is)
	__rat.assign(__num, __den);

      return __is;
    }

  // Add manipulators for output format?
  template<typename _IntTp>
    std::ostream&
    operator<<(std::ostream& __os, const _Rational<_IntTp>& __rat)
    {
      __os << __rat.num() << '/' << __rat.den();
      return __os;
    }

  // Do not use any abs() defined on _IntTp - it isn't worth it, given the
  // difficulties involved (Koenig lookup required, there may not *be* an abs()
  // defined, etc etc).
  template<typename _IntTp>
    inline _Rational<_IntTp>
    abs(const _Rational<_IntTp>& __rat)
    {
      if (__rat.num() >= _IntTp{0})
	return __rat;
      else
	return _Rational<_IntTp>(-__rat.num(), __rat.den());
    }

} // namespace __gnu_cxx

#endif // _RATIONAL_H