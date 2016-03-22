
#include <vector>

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

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
      {
	operator+=(__first_term);
      }

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
	this->_M_sum = __sum;
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

/*
// This is an example of the conversion of a non-alternating series to an alternating one
// for the Lerch transcendent.
// Think of a way to generalize this.
// Also, why not throw in the (-1)^j?  DONE..
// I think in place of z, s, a you hand in a functor that, given j, returns the b_j coefficient.
template<typename _Tp>
  _Tp
  aj(_Tp __z, _Tp __s, _Tp __a, unsigned long long j, _Tp __eps)
  {
    constexpr auto _S_min = std::numeric_limits<_Tp>::min();
    // Maximum number of iterations before 2^k overflow.
    constexpr auto __k_max = std::numeric_limits<unsigned long long>::digits;

    auto __sum = _Tp{}
    auto __two2k = 1ULL;
    for (auto __k = 0ULL; __k < __k_max; ++__k)
      {
	// Index for the term in the original series.	auto __ind = 0ULL;
	if (!__builtin_mul_overflow(__two2k, __j + 1, &__ind))
	  throw std::runtime_error("aj: integer overflow");	  
	--__ind;

	// Increment the sum.
	auto __z2ind = std::pow(__z, __ind);
	auto __bjk = __two2k * __z2ind / std::pow(__a + __ind, __s);
	__sum += bjk;

	// Stop summation if either sum is zero or |term/sum|
	// is below requested accuracy.
	if (std::abs(__sum) <= _S_min || std::abs(__bjk / __sum) < 1.0e-2 * __eps)
	  break;

	++__k;
	if (!__builtin_mul_overflow(__two2k, 2, &__two2k))
	  throw std::runtime_error("aj: integer overflow");
      }

    auto __sign = (__j % 2 == 1 ? -1 : +1);
    return __sign * __sum;
  }
*/
_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std
