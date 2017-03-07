
#include <cstddef>
#include <limits>

template<typename _Tp, typename _Func>
  class trapezoid_integral
  {
  public:

    trapezoid_integral(_Func __fun, _Tp __a, _Tp __b, _Tp __err)
    : _M_fun(__fun), _M_a(__a), _M_b(__b), _M_err(__err), _M_sum()
    { }

    _Tp operator()();

  private:

    _Tp step();

    _Func _M_fun;
    _Tp _M_a;
    _Tp _M_b;
    _Tp _M_err;
    _Tp _M_sum;
    std::size_t _M_iter = 0;
    std::size_t _M_pow2 = 0;
    std::size_t _S_max_iter = std::numeric_limits<std::size_t>::digits - 1;
  };

#include "trapezoid_integral.tcc"
