
#ifndef _RATIONAL_H
#define _RATIONAL_H 1

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

#include <ext/rational.tcc>

#endif // _RATIONAL_H
