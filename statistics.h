#ifndef STATISTICS_H
#define STATISTICS_H 1

#include <bits/complex_util.h>

/**
 *  Incremental computation of statistics.
 */
template<typename _Tp>
  struct _Statistics
  {
    _Statistics&
    operator<<(_Tp __diff)
    {
      ++_M_count;
      auto __old_mean = _M_mean;
      _M_mean = (_M_type(__diff) + _M_type(_M_count - 1) * _M_mean) / _M_type(_M_count);
      auto __del_mean = _M_mean - __old_mean;
      auto __del_diff = _M_type(__diff) - _M_mean;
      if (_M_count > 1)
	_M_variance = (_M_type(_M_count - 2) * _M_variance * _M_variance
		    + _M_type(_M_count - 1) * __del_mean * __del_mean
		    + __del_diff * __del_diff) / _M_type(_M_count - 1);
      if (__diff < _M_min)
	_M_min = __diff;
      if (__diff > _M_max)
	_M_max = __diff;

      return *this;
    }

    static constexpr bool _M_is_complex = __gnu_cxx::is_complex_v<_Tp>;

    using _M_type = std::conditional_t<__gnu_cxx::is_complex_v<_Tp>,
				       std::complex<long double>, long double>;

    _Tp
    count() const
    { return _Tp(_M_count); }

    _Tp
    mean() const
    { return _Tp(_M_mean); }

    _Tp
    variance() const
    { return _Tp(_M_variance); }

    _Tp
    std_deviation() const
    { return _Tp(std::sqrt(_M_variance)); }

    _Tp
    min() const
    { return _Tp(_M_min); }

    _Tp
    max() const
    { return _Tp(_M_max); }

    std::size_t _M_count = 0;
    _M_type _M_mean = 0;
    _M_type _M_variance = 0;
    _M_type _M_min = std::numeric_limits<long double>::max();
    _M_type _M_max = std::numeric_limits<long double>::lowest();
  };

#endif // STATISTICS_H
