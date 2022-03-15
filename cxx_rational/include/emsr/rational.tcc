
#ifndef RATIONAL_TCC
#define RATIONAL_TCC 1

namespace emsr
{

  // Arithmetic assignment operators
  template<typename IntTp>
    Rational<IntTp>&
    Rational<IntTp>::operator+=(Rational<IntTp> rat)
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

      value_type gcd = std::experimental::gcd(this->m_den, rat.den());
      this->m_den /= gcd;  // = b1 from the calculations above
      this->m_num = this->m_num * (rat.den() / gcd) + rat.num() * this->m_den;
      gcd = std::experimental::gcd(this->m_num, gcd);
      this->m_num /= gcd;
      this->m_den *= rat.den() / gcd;

      return *this;
    }

  template<typename IntTp>
    Rational<IntTp>&
    Rational<IntTp>::operator-=(Rational<IntTp> rat)
    {
      // This calculation avoids overflow, and minimises the number of expensive
      // calculations. It corresponds exactly to the += case above
      value_type gcd = std::experimental::gcd(this->m_den, rat.den());
      this->m_den /= gcd;
      this->m_num = this->m_num * (rat.den() / gcd) - rat.num() * this->m_den;
      gcd = std::experimental::gcd(this->m_num, gcd);
      this->m_num /= gcd;
      this->m_den *= rat.den() / gcd;

      return *this;
    }

  template<typename IntTp>
    Rational<IntTp>&
    Rational<IntTp>::operator*=(Rational<IntTp> rat)
    {
      // Avoid overflow and preserve normalization
      value_type gcd1 = std::experimental::gcd(this->m_num, rat.den());
      value_type gcd2 = std::experimental::gcd(rat.num(), this->m_den);
      this->m_num = (this->m_num / gcd1) * (rat.num() / gcd2);
      this->m_den = (this->m_den / gcd2) * (rat.den() / gcd1);
      return *this;
    }

  template<typename IntTp>
    Rational<IntTp>&
    Rational<IntTp>::operator/=(Rational<IntTp> rat)
    {
      constexpr value_type zero{0};

      // Trap division by zero
      if (rat.num() == zero)
	throw BadRational();
      if (this->m_num == zero)
	return *this;

      // Avoid overflow and preserve normalization
      value_type gcd1 = std::experimental::gcd(this->m_num, rat.num());
      value_type gcd2 = std::experimental::gcd(rat.den(), this->m_den);
      this->m_num = (this->m_num / gcd1) * (rat.den() / gcd2);
      this->m_den = (this->m_den / gcd2) * (rat.num() / gcd1);

      if (this->m_den < zero)
	{
          this->m_num = -this->m_num;
          this->m_den = -this->m_den;
	}
      return *this;
    }

  template<typename IntTp>
    bool
    Rational<IntTp>::operator<(const Rational<IntTp>& r) const
    {
      constexpr value_type zero{0};

      assert(this->m_den > zero);
      assert(r.den() > zero);

      // Determine relative order by expanding each value to its simple continued
      // fraction representation using the Euclidian GCD algorithm.
      struct
      { value_type n, d, q, r; }
      ts{
	this->m_num, this->m_den,
	static_cast<value_type>(this->m_num / this->m_den),
	static_cast<value_type>(this->m_num % this->m_den)
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

  template<typename IntTp>
    bool
    Rational<IntTp>::operator<(value_type i) const
    {
      constexpr value_type zero{0};

      // Break value into mixed-fraction form, w/ always-nonnegative remainder
      assert(this->m_den > zero);
      value_type quo = this->m_num / this->m_den,
		 rem = this->m_num % this->m_den;
      while (rem < zero)
	{
	  rem += this->m_den;
	  --quo;
	}

      // Compare with just the quotient, since the remainder always bumps the
      // value up.  [Since q = floor(n/d), and if n/d < i then q < i, if n/d == i
      // then q == i, if n/d == i + r/d then q == i, and if n/d >= i + 1 then
      // q >= i + 1 > i; therefore n/d < i iff q < i.]
      return quo < i;
    }

  // Normalisation
  template<typename IntTp>
    void
    Rational<IntTp>::m_normalize()
    {
      constexpr value_type zero{0};

      if (this->m_den == zero)
	throw BadRational();

      // Handle the case of zero separately, to avoid division by zero
      if (this->m_num == zero)
	{
          this->m_den = value_type{1};
          return;
	}

      value_type gcd = std::experimental::gcd(this->m_num, this->m_den);

      this->m_num /= gcd;
      this->m_den /= gcd;

      // Ensure that the denominator is positive
      if (this->m_den < zero)
	{
          this->m_num = -this->m_num;
          this->m_den = -this->m_den;
	}

      assert(this->m_valid());
    }

  // Input and output
  template<typename IntTp>
    std::istream&
    operator>>(std::istream& is, Rational<IntTp>& rat)
    {
      using value_type = typename Rational<IntTp>::value_type;

      auto num = value_type{0}, den = value_type{1};
      char c = 0;
      detail::resetter sentry(is);

      is >> num;
      c = is.get();

      if (c != '/')
	is.clear(std::ios_base::badbit);  // old GNU c++ lib has no ios_base

  #if !defined(__GNUC__) || (defined(__GNUC__) && (__GNUC__ >= 3))
      is >> std::noskipws;
  #else
      is.unsetf(ios::skipws); // compiles, but seems to have no effect.
  #endif
      is >> den;

      if (is)
	rat.assign(num, den);

      return is;
    }

} // namespace emsr

#endif // RATIONAL_TCC
