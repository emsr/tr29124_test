
#include <vector>

namespace std _GLIBCXX_VISIBILITY(default)
{
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION

  template<typename _Tp>
    class __vanWijngaardenSum
    {
    public:
      __vanWijngaardenSum() = default;

      explicit __vanWijngaardenSum(_Tp __sum)
      : _M_sum(__sum), _M_delta{}
      { }

      __vanWijngaardenSum&
      operator+=(_Tp __term)
      {
	if (__isnan(__term))
	  throw std::runtime_error("__vanWijngaardenSum: bad term");
	if (std::fabs(__term) == std::numeric_limits<_Tp>::infinity())
	  throw std::runtime_error("__vanWijngaardenSum: infinite term");

	if (this->_M_delta.size() == 0)
	  {
	    this->_M_delta.push_back(__term);
	    this->_M_sum += _Tp{0.5L} * this->_M_delta.back();
	  }
	else
	  {
	    auto __temp = this->_M_delta[0];
	    this->_M_delta[0] = __term;
	    auto __n = this->_M_delta.size();
	    for (auto __j = 0;
		 __j < __n - 1; ++__j)
	      __temp = std::exchange(this->_M_delta[__j + 1],
				_Tp{0.5L} * (this->_M_delta[__j] + __temp));
	    auto __next = _Tp{0.5L} * (this->_M_delta.back() + __temp);
	    if (std::fabs(__next) < std::fabs(this->_M_delta.back()))
	      {
		this->_M_delta.push_back(__next);
		this->_M_sum += _Tp{0.5L} * this->_M_delta.back();
	      }
	    else
	      this->_M_sum += __next;
	  }
	return *this;
      }

      __vanWijngaardenSum&
      operator-=(_Tp __term)
      { return operator+=(-__term); }

      _Tp
      operator()() const
      { return this->_M_sum; }

      __vanWijngaardenSum&
      reset(_Tp __sum = _Tp{})
      {
	this->_M_sum = __sum;
	this->_M_delta.clear();
	return *this;
      }

    private:
      _Tp _M_sum;
      std::vector<_Tp> _M_delta;
    };

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std
