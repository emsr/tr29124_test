
#include <limits>
#include <cmath>

namespace emsr
{

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
      : m_sum{}, m_term{}, m_d{}, m_b{-1}, m_c{1}, m_num_terms{0}, m_converged{false}
      {
        constexpr auto eps = std::numeric_limits<value_type>::epsilon();
        constexpr auto m = std::numeric_limits<long long>::max() / 6;
        this->m_max_num_terms = std::log(4 / eps) / std::log(fact);
        long long d_prev = 17;
        long long d_curr = 3;
        for (unsigned int n = 0; n <= this->m_max_num_terms; ++n)
        {
          if (d_curr > m)
          {
            this->m_max_num_terms = n;
            break;
          }
          const auto d_next = 6 * d_curr - d_prev;
          d_prev = d_curr;
          d_curr = d_next;
        }
        this->m_d = d_curr;
        this->m_c = -this->m_d;
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
      { return this->m_sum / this->m_d; }

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
	// Keep m_d
        // Keep m_max_num_terms;
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

    private:

      value_type m_sum;
      value_type m_term;
      long long m_d;
      value_type m_b;
      value_type m_c;
      std::size_t m_max_num_terms;
      std::size_t m_num_terms;
      bool m_converged;
    };

} // namespace emsr
