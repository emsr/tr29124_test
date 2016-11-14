
#include <stdexcept>
#include <bits/numeric_limits.h>

/**
 * You need a functor belching a, b
 * an iterator pair for a and an iterator for b
 *  
 */
template<typename _Tp, typename _AFun, typename _BFun, typename _TailFun>
  class _ForwardContinuedFraction
  {
  private:

    _AFun _M_num;
    _BFun _M_den;
    _TailFun _M_tail;

  public:

    using _ARet = decltype(_M_num(0ull, _Tp{}));
    using _BRet = decltype(_M_den(0ull, _Tp{}));
    using _Ret = decltype(_ARet(0, _Tp{}) / _BRet(0, _Tp{}));

    _ForwardContinuedFraction(_AFun __a, _BFun __b, _TailFun __w)
    : _M_num(__a),
      _M_den(__b),
      _M_tail(__w)
    { }

    _Ret
    operator()(_Tp __x) const
    {
      const auto _S_fp_min = __gnu_cxx::__min(__x);
      const auto _S_eps = __gnu_cxx::__epsilon(__x);
      const int _S_max_iter = 1000;

      auto __b = _M_den(1, __x);
      auto _Den2 = __b;
      auto __a = _M_num(1, __x);
      auto _Num2 = __a;
      __b = _M_den(2, __x);
      auto _Num1 = __b * _Num2;
      __a = _M_num(2, __x);
      auto _Den1 = __a + __b * _Den2;
      std::size_t __k = 3;
      while (true)
	{
	  __a = _M_num(__i, __x);
	  __b = _M_den(__i, __x);
	  auto _Num = __b * _Num1 + __a * _Num2;
	  auto _Den = __b * _Den1 + __a * _Den2;
	  _Num2 = _Num1;
	  _Den2 = _Den1;
	  _Num1 = _Num;
	  _Den1 = _Den;
	}
       += _M_den(0, __x);
    }
  };

