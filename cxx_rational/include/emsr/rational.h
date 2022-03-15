
#ifndef RATIONAL_H
#define RATIONAL_H 1

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

namespace emsr
{

  class BadRational
  : public std::domain_error
  {
  public:
    explicit
    BadRational()
    : std::domain_error("BadRational: zero denominator")
    { }
  };

  template<typename IntTp>
    class Rational;

  template<typename IntTp>
    Rational<IntTp>
    abs(const Rational<IntTp>&);

  template<typename IntTp = std::intmax_t>
    class Rational
    {
      static_assert(std::numeric_limits<IntTp>::is_specialized, "");

    public:
      typedef IntTp value_type;

      constexpr
      Rational()
      : m_num{0}, m_den{1}
      { }

      constexpr explicit
      Rational(value_type num)
      : m_num(num), m_den(1)
      { }

      template<intmax_t Num, intmax_t Den>
	constexpr explicit
	Rational(std::ratio<Num, Den>)
	: m_num(std::ratio<Num, Den>::num),
	  m_den(std::ratio<Num, Den>::den)
	{ }

      Rational(value_type num, value_type den)
      : m_num(num), m_den(den)
      { this->m_normalize(); }

      // Default copy constructor and assignment are fine

      // Add assignment from value_type
      Rational&
      operator=(value_type num)
      { return this->assign(num, 1); }

      // Assign in place
      Rational&
      assign(value_type, value_type);

      // Access to representation
      constexpr value_type
      num() const
      { return this->m_num; }

      constexpr value_type
      den() const
      { return this->m_den; }

      // Arithmetic assignment operators
      Rational& operator+=(Rational);
      Rational& operator-=(Rational);
      Rational& operator*=(Rational);
      Rational& operator/=(Rational);

      Rational& operator+=(value_type);
      Rational& operator-=(value_type);
      Rational& operator*=(value_type);
      Rational& operator/=(value_type);

      // Increment and decrement
      const Rational& operator++();
      const Rational& operator--();

      // Operator not
      bool
      operator!() const
      { return !this->m_num; }

      // Boolean conversion
      operator
      bool() const
      { return this->operator!(); }

      // Comparison operators
      bool operator<(const Rational&) const;
      bool operator==(const Rational&) const;

      bool operator<(value_type) const;
      bool operator>(value_type) const;
      bool operator==(value_type) const;

    private:

      value_type m_num;
      value_type m_den;

      // Representation note: Fractions are kept in normalized form at all
      // times. normalized form is defined as gcd(num,den) == 1 and den > 0.
      // In particular, note that the implementation of abs() below relies
      // on den always being positive.
      bool m_valid() const;
      void m_normalize();
    };

  // Assign in place
  template<typename IntTp>
    inline Rational<IntTp>&
    Rational<IntTp>::assign(value_type num, value_type den)
    {
      this->m_num = num;
      this->m_den = den;
      this->m_normalize();
      return *this;
    }

  // Unary plus and minus
  template<typename IntTp>
    inline Rational<IntTp>
    operator+(const Rational<IntTp>& rat)
    { return rat; }

  template<typename IntTp>
    inline Rational<IntTp>
    operator-(const Rational<IntTp>& rat)
    { return Rational<IntTp>(-rat.num(), rat.den()); }

  // Mixed-mode operators
  template<typename IntTp>
    inline Rational<IntTp>&
    Rational<IntTp>::operator+=(value_type i)
    { return this->operator+=(Rational<IntTp>(i)); }

  template<typename IntTp>
    inline Rational<IntTp>&
    Rational<IntTp>::operator-=(value_type i)
    { return this->operator-=(Rational<IntTp>(i)); }

  template<typename IntTp>
    inline Rational<IntTp>&
    Rational<IntTp>::operator*=(value_type i)
    { return this->operator*=(Rational<IntTp>(i)); }

  template<typename IntTp>
    inline Rational<IntTp>&
    Rational<IntTp>::operator/=(value_type i)
    { return this->operator/=(Rational<IntTp>(i)); }

  // Increment and decrement
  template<typename IntTp>
    inline const Rational<IntTp>&
    Rational<IntTp>::operator++()
    {
      // This can never denormalise the fraction
      this->m_num += this->m_den;
      return *this;
    }

  template<typename IntTp>
    inline const Rational<IntTp>&
    Rational<IntTp>::operator--()
    {
      // This can never denormalise the fraction
      this->m_num -= this->m_den;
      return *this;
    }

  template<typename IntTp>
    bool
    Rational<IntTp>::operator>(value_type i) const
    {
      // Trap equality first
      if (this->m_num == i && this->m_den == IntTp(1))
	return false;

      // Otherwise, we can use operator<
      return !this->operator<(i);
    }

  template<typename IntTp>
    inline bool
    Rational<IntTp>::operator==(const Rational<IntTp>& rat) const
    { return ((this->m_num == rat.num) && (this->m_den == rat.den())); }

  template<typename IntTp>
    inline bool
    Rational<IntTp>::operator==(value_type i) const
    { return ((this->m_den == IntTp{1}) && (this->m_num == i)); }

  // Invariant check
  template<typename IntTp>
    inline bool
    Rational<IntTp>::m_valid() const
    {
      return this->m_den > value_type{0}
          && std::experimental::gcd(this->m_num, this->m_den) == value_type{1};
    }

  // Type conversion
  template<typename Tp, typename IntTp>
    inline Tp
    Rational_cast(const Rational<IntTp>& src)
    {
      return static_cast<Tp>(src.num())
           / static_cast<Tp>(src.den());
    }

  // Symmetric math operators.
  template<typename IntTp>
    Rational<IntTp>
    operator+(const Rational<IntTp>& r1, const Rational<IntTp>& r2)
    { return Rational<IntTp>(r1) += r2; }

  template<typename IntTp, typename RealTp>
    RealTp
    operator+(RealTp r1, const Rational<IntTp>& r2)
    { return r1 + Rational_cast<RealTp>(r2); }

  template<typename IntTp, typename RealTp>
    RealTp
    operator+(const Rational<IntTp>& r1, RealTp r2)
    { return Rational_cast<RealTp>(r1) + r2; }

  template<typename IntTp>
    Rational<IntTp>
    operator-(const Rational<IntTp>& r1, const Rational<IntTp>& r2)
    { return Rational<IntTp>(r1) -= r2; }

  template<typename IntTp, typename RealTp>
    RealTp
    operator-(RealTp r1, const Rational<IntTp>& r2)
    { return r1 - Rational_cast<RealTp>(r2); }

  template<typename IntTp, typename RealTp>
    RealTp
    operator-(const Rational<IntTp>& r1, RealTp r2)
    { return Rational_cast<RealTp>(r1) - r2; }

  template<typename IntTp>
    Rational<IntTp>
    operator*(const Rational<IntTp>& r1, const Rational<IntTp>& r2)
    { return Rational<IntTp>(r1) *= r2; }

  template<typename IntTp, typename RealTp>
    RealTp
    operator*(RealTp r1, const Rational<IntTp>& r2)
    { return r1 * Rational_cast<RealTp>(r2); }

  template<typename IntTp, typename RealTp>
    RealTp
    operator*(const Rational<IntTp>& r1, RealTp r2)
    { return Rational_cast<RealTp>(r1) * r2; }

  template<typename IntTp>
    Rational<IntTp>
    operator/(const Rational<IntTp>& r1, const Rational<IntTp>& r2)
    { return Rational<IntTp>(r1) /= r2; }

  template<typename IntTp, typename RealTp>
    RealTp
    operator/(RealTp r1, const Rational<IntTp>& r2)
    { return r1 / Rational_cast<RealTp>(r2); }

  template<typename IntTp, typename RealTp>
    RealTp
    operator/(const Rational<IntTp>& r1, RealTp r2)
    { return Rational_cast<RealTp>(r1) / r2; }

  namespace detail
  {

    // A utility class to reset the format flags for an istream at end
    // of scope, even in case of exceptions
    struct resetter
    {
      resetter(std::istream& is)
      : m_is(is), m_fl(is.flags())
      { }
      ~resetter()
      { m_is.flags(m_fl); }

      std::istream& m_is;
      std::ios_base::fmtflags m_fl;
    };

  }

  // Add manipulators for output format?
  template<typename IntTp>
    std::ostream&
    operator<<(std::ostream& os, const Rational<IntTp>& rat)
    {
      os << rat.num() << '/' << rat.den();
      return os;
    }

  // Do not use any abs() defined on IntTp - it isn't worth it, given the
  // difficulties involved (Koenig lookup required, there may not *be* an abs()
  // defined, etc etc).
  template<typename IntTp>
    inline Rational<IntTp>
    abs(const Rational<IntTp>& rat)
    {
      if (rat.num() >= IntTp{0})
	return rat;
      else
	return Rational<IntTp>(-rat.num(), rat.den());
    }

} // namespace emsr

#include <emsr/rational.tcc>

#endif // RATIONALH
