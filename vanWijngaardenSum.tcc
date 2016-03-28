
#include <vector>

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  /**
   *
   */
  template<typename _Tp>
    class _VanWijngaardenSum
    {
    public:

      using value_type = _Tp;

      ///  Default constructor.
      _VanWijngaardenSum() = default;

      ///  Constructor taking the first term.
      explicit _VanWijngaardenSum(value_type __first_term)
      : _VanWijngaardenSum{}
      { operator+=(__first_term); }

      /// Add a new term to the sum.
      _VanWijngaardenSum&
      operator+=(value_type __term)
      {
	if (__isnan(__term))
	  throw std::runtime_error("_VanWijngaardenSum: bad term");
	if (std::abs(__term) == std::numeric_limits<value_type>::infinity())
	  throw std::runtime_error("_VanWijngaardenSum: infinite term");

	++this->_M_num_terms;

	if (this->_M_delta.size() == 0)
	  {
	    this->_M_delta.push_back(__term);
	    this->_M_sum += value_type{0.5L} * this->_M_delta.back();
	  }
	else
	  {
	    auto __temp = this->_M_delta[0];
	    this->_M_delta[0] = __term;
	    auto __n = this->_M_delta.size();
	    for (auto __j = 0;
		 __j < __n - 1; ++__j)
	      __temp = std::exchange(this->_M_delta[__j + 1],
				value_type{0.5L} * (this->_M_delta[__j] + __temp));
	    auto __next = value_type{0.5L} * (this->_M_delta.back() + __temp);
	    if (std::abs(__next) < std::abs(this->_M_delta.back()))
	      {
		this->_M_delta.push_back(__next);
		this->_M_sum += value_type{0.5L} * this->_M_delta.back();
	      }
	    else
	      this->_M_sum += __next;
	  }
	return *this;
      }

      /// Subtract a new term from the sum.
      _VanWijngaardenSum&
      operator-=(value_type __term)
      { return operator+=(-__term); }

      /// Return true if the sum converged.
      operator
      bool() const
      { return !this->_M_converged; }

      /// Return the current value of the sum.
      value_type
      operator()() const
      { return this->_M_sum; }

      /// Return the current number of terms contributing to the sum.
      std::size_t
      num_terms() const
      { return this->_M_num_terms; }

      ///  Reset the sum to it's initial state.
      _VanWijngaardenSum&
      reset()
      {
	this->_M_sum = value_type{};
	this->_M_delta.clear();
	this->_M_num_terms = 0;
	this->_M_converged = false;
	return *this;
      }

      ///  Restart the sum with the first new term.
      _VanWijngaardenSum&
      reset(value_type __first_term)
      {
	this->reset();
	this->operator+=(__first_term);
	return *this;
      }

    private:

      value_type _M_sum;
      std::vector<value_type> _M_delta;
      std::size_t _M_num_terms;
      bool _M_converged;
    };

//
// A functor for a vanWijnGaarden compressor must have
// _Tp operator()(int) that returns a term in the original defining series.
//
template<typename _Tp>
  class __lerch_term
  {
  public:

    using value_type = _Tp;

    __lerch_term(value_type __z, value_type __s, value_type __a)
    : _M_z{__z}, _M_s{__s}, _M_a{__a}
    { }

    value_type
    operator()(std::size_t __i) const
    {
      return std::pow(_M_z, value_type(__i))
	   / std::pow(_M_a + value_type(__i), _M_s);
    }

  private:

    value_type _M_z;
    value_type _M_s;
    value_type _M_a;
  };

  /**
   *  This performs a series compression on a monotone series - converting
   *  it to an alternating series - for the regular van Wijngaarden sum.
   *  ADL for ctors anyone?  I'd like to put a lambda in here*
   */
  template<typename _TermFn>
    class _VanWijngaardenCompressor
    {
    public:

      _VanWijngaardenCompressor(_TermFn __term_fn)
      : _M_term_fn{__term_fn}
      { }

      auto
      operator[](std::size_t __j) const
      {
	using value_type = decltype(this->_M_term_fn(__j));
	constexpr auto _S_min = std::numeric_limits<value_type>::min();
	constexpr auto _S_eps = std::numeric_limits<value_type>::epsilon();
	// Maximum number of iterations before 2^k overflow.
	constexpr auto __k_max = std::numeric_limits<std::size_t>::digits;

	auto __sum = value_type{};
	auto __two2k = std::size_t{1};
	for (auto __k = std::size_t{0}; __k < __k_max; __k += std::size_t{2})
	  {
	    // Index for the term in the original series.
	    auto __i = std::size_t{0};
	    if (__builtin_mul_overflow(__two2k, __j + 1, &__i))
	      throw std::runtime_error("_VanWijngaardenCompressor: "
				       "index overflow");	  
	    --__i;

	    // Increment the sum.
	    auto __term = __two2k * this->_M_term_fn(__i);
	    __sum += __term;

	    // Stop summation if either the sum is zero
	    // or if |term / sum| is below requested accuracy.
	    if (std::abs(__sum) <= _S_min
	     || std::abs(__term / __sum) < value_type{1.0e-2} * _S_eps)
	      break;

	    if (__builtin_mul_overflow(__two2k, std::size_t{2}, &__two2k))
	      throw std::runtime_error("_VanWijngaardenCompressor: "
				       "index overflow");
	  }

	auto __sign = (__j % 2 == 1 ? -1 : +1);
	return __sign * __sum;
      }

    private:

      _TermFn _M_term_fn;
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std
