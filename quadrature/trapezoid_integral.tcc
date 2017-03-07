
#include <cmath>

template<typename _Tp, typename _Func>
  _Tp
  trapezoid_integral<_Tp, _Func>::operator()()
  {
    auto __sum_prev = -std::numeric_limits<_Tp>::max() / 1000;
    for (std::size_t __j = 0; __j < _S_max_iter; ++__j)
      {
	auto __sum = this->step();
	if (std::abs(__sum - __sum_prev) < _M_err * std::abs(__sum_prev))
	  return __sum;
	if (std::abs(__sum) < _M_err && std::abs(__sum_prev) < _M_err && __j > 6 )
	  return __sum;
	__sum_prev = __sum;
      }
  }

template<typename _Tp, typename _Func>
  _Tp
  trapezoid_integral<_Tp, _Func>::step()
  {
    if (_M_iter == 0)
      {
	_M_iter = 1;
        _M_sum = (_M_b - _M_a) * (_M_fun(_M_a) + _M_fun(_M_b)) / _Tp{2};
        _M_pow2 = 1;
      }
    else
      {
        const auto __del = (_M_b - _M_a) / _M_pow2;
        auto __x = _M_a + __del / _Tp{2};
	auto __sum = _Tp{0};
        for (std::size_t __j = 1; __j <= _M_pow2; ++__j, __x += __del)
	  __sum += _M_fun(__x);
        _M_sum = (__sum + (_M_b - _M_a) * __sum / _M_pow2) / _Tp{2};
        _M_pow2 *= 2;
      }
    return _M_sum;
  }
