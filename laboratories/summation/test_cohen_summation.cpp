/*
g++ -g -o test_cohen_summation test_cohen_summation.cpp
./test_cohen_summation > test_cohen_summation.txt
*/

#include <iostream>
#include <iomanip>
#include <vector>

#include <array>
#include <limits>
#include <cmath>

// m_d[0] is never used in the wild (or shouldn't be).

  /**
   * This is a Cohen summation accelerator for alternating series:
   * @f[
   *   S = \sum_{k=0}^\infty a_k
   * @f]
   *
   * This summation computes a maximum number of terms based either on
   * @f[
   *   n > \frac{\log(4/\epsilon)}{\log(3 + \sqrt{8})}
   * @f]
   * where @f$ \epsilon @f$ is machine epsilon.
   * Or when the three-term recurrence relation on integers fails:
   * @f[
   *   d_{n+1} = 6 d_n - d_{n-1}
   * @f]
   * with starting values @f$ d_{-1} = 17 @f$ and @f$ d_0 = 3 @f$
   * The series will not be marked converged until that number of terms is reached
   * and will not accept more after.
   *
   * @see Henri Cohen, Fernando Rodriguez Villegas, and Don Zagier,
   *      Convergence Acceleration of Alternating Series,
   *      Experimental Mathematics 9:1, page 3
   */
  template<typename Tp>
    class CohenSum
    {
    public:

      using value_type = Tp;

      ///  Default constructor.
      CohenSum()
      : m_sum{}, m_term{}, m_d{17, 3}, m_b{-1}, m_c{1}, m_num_terms{0}, m_converged{false}
      {
        constexpr auto eps = std::numeric_limits<value_type>::epsilon();
        constexpr auto m = std::numeric_limits<long long>::max() / 6;
        this->m_max_num_terms = std::log(4 / eps) / std::log(fact);
        for (unsigned int n = 0; n <= this->m_max_num_terms; ++n)
        {
          if (this->m_d[1] > m)
          {
            this->m_max_num_terms = n;
std::cout << "BREAK\n";
            break;
          }
          const auto d = 6 * this->m_d[1] - this->m_d[0];
          this->m_d[0] = this->m_d[1];
          this->m_d[1] = d;
        }
        this->m_c = -this->m_d[1];
std::cout << "    n = " << std::setw(22) << this->m_max_num_terms << "; d = " << std::setw(22) << this->m_d[1] << "; b = " << std::setw(22) << this->m_b << "; c = " << std::setw(22) << this->m_c << '\n';
      }

      /// Add a new term to the sum.
      CohenSum&
      operator+=(value_type term)
      {
        if (this->m_converged)
          return *this;
        this->m_c = this->m_b - this->m_c;
        this->m_term = this->m_c * term;
        this->m_sum += this->m_term;
        const std::ptrdiff_t n = this->m_max_num_terms;
        const std::ptrdiff_t k = this->m_num_terms;
        this->m_b *= value_type{2} * (k + n) * (k - n) / (2 * k + 1) / (k + 1);
        ++this->m_num_terms;
        if (this->m_num_terms >= this->m_max_num_terms)
          this->m_converged = true;
std::cout << "    n = " << std::setw(2) << n << "; k = " << std::setw(2) << k;
std::cout << "    b = " << std::setw(22) << this->m_b << "; c = " << std::setw(22) << this->m_c << "; s = " << std::setw(22) << this->m_sum << '\n';
        return *this;
      }

      /// Subtract a new term from the sum.
      CohenSum&
      operator-=(value_type term)
      { return this->operator+=(-term); }

      /// Return true if the sum converged.
      bool
      converged() const
      { return this->m_converged; }

      /// Return false if the sum converged.
      operator
      bool() const
      { return !this->converged(); }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->m_sum / this->m_d[1]; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->m_num_terms; }

      /// Return the maximum number of terms contributing to the sum.
      std::size_t
      max_num_terms() const
      { return this->m_max_num_terms; }

      /// Return the current last term contributing to the sum.
      value_type
      term() const
      { return this->m_term; }

      static constexpr auto sqrt8 = value_type{2.828427124746190097603377448419396157138L};
      static constexpr auto fact = value_type{3} + sqrt8;

      ///  Reset the sum to it's initial state.
      CohenSum&
      reset()
      {
	this->m_sum = value_type{};
	this->m_term = value_type{};
	//this->m_d = {17, 3};
        //this->m_max_num_terms;
	this->m_num_terms = 0;
	this->m_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      CohenSum&
      reset(value_type first_term)
      {
	this->reset();
	this->operator+=(first_term);
	return *this;
      }

    //private:

      value_type m_sum;
      value_type m_term;
      std::array<long long, 2> m_d;
      value_type m_b;
      value_type m_c;
      std::size_t m_max_num_terms;
      std::size_t m_num_terms;
      bool m_converged;
    };

int
main()
{
  constexpr auto m = std::numeric_limits<long long>::max() / 6;
  constexpr auto w = 1 + std::numeric_limits<long long>::digits10;
  std::cout.precision(w);

  CohenSum<long double> CS;

  long double d_fp = 1;
  long double cc = CS.fact;
  long double fact = 1;
  CS.m_d = {17, 3};
  for (int k = 0; k < 50; ++k)
  {
    if (CS.m_d[1] > m)
      break;

    const auto d = 6 * CS.m_d[1] - CS.m_d[0];
    CS.m_d[0] = CS.m_d[1];
    CS.m_d[1] = d;
    std::cout << std::setw(2) << k << " d = " << std::setw(w) << CS.m_d[1] << std::setw(6 + w) << d_fp << '\n';
    fact *= cc;
    d_fp = (fact + 1 / fact) / 2;
  }
  std::cout << "   m = " << std::setw(w) << std::numeric_limits<long long>::max() << '\n';
  std::cout << "   n = " << std::setw(w) << CS.m_max_num_terms << '\n';
  std::cout << "   d = " << std::setw(w) << CS.m_d[1] << '\n';
  constexpr auto pipid12 = 0.8224670334241132182362075833230125946103L;
  long double naive = 0, sign = 1;
  for (unsigned int k = 0; k < CS.max_num_terms(); ++k)
  {
    CS += 1.0L / ((k + 1) * (k + 1));
    naive += sign / ((k + 1) * (k + 1));
    sign = -sign;
  }
  std::cout << "   CS    = " << std::setw(w) << CS() << '\n';
  std::cout << "   naive = " << std::setw(w) << naive << '\n';
  std::cout << "   true  = " << std::setw(w) << pipid12 << '\n';

  // https://www.johndcook.com/blog/2020/08/06/cohen-acceleration/
  for (int n = 1; n <= 25; ++n)
  {
    long double d = std::pow(3.0L + std::sqrt(8.0L), n);
    d = (d + 1 / d) / 2;
    std::cout << '\n';
    std::cout << "  n = " << std::setw(2) << n << '\n';
    std::cout << "  d     = " << std::setw(w) << d << '\n';

    std::vector<long double> a(n);
    long double naive = 0;
    long double sign = 1;
    for (int k = 0; k < n; ++k)
    {
      a[k] = 1.0L / (k + 1) / (k + 1);
      naive += sign * a[k];
      sign = -sign;
    }

    long double b = -1;
    long double c = -d;
    long double s = 0;
    for (int k = 0; k < n; ++k)
    {
      c = b - c;
      s += c * a[k];
      b *= (k + n) * (k - n) / ((k + 0.5) * (k + 1));
      std::cout << "    b = " << std::setw(w+2) << b << "; c = " << std::setw(w+2) << c << "; s = " << std::setw(w+2) << s << '\n';
    }

    std::cout << "  naive = " << std::setw(w) << naive << '\n';
    std::cout << "  sum   = " << std::setw(w) << s / d << '\n';
    std::cout << "  err   = " << std::setw(w) << s / d - pipid12 << '\n';
  }
}
