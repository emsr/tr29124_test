
#include <stdexcept>
#include <bits/numeric_limits.h>

/**
 * You need a functor belching a, b
 * an iterator pair for a and an iterator for b
 *  
 */
template<typename _Tp, typename _AFun, typename _BFun, typename _TailFun>
  class _LentzContinuedFraction
  {
  private:

    _AFun _M_num;
    _BFun _M_den;
    _TailFun _M_tail;

  public:

    using _ARet = decltype(_M_num(0ull, _Tp{}));
    using _BRet = decltype(_M_den(0ull, _Tp{}));
    using _Ret = decltype(_ARet(0, _Tp{}) / _BRet(0, _Tp{}));

    _LentzContinuedFraction(_AFun __a, _BFun __b, _TailFun __w)
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
      _Ret __c(_Tp{1} / _S_fp_min);
      auto __d(_Tp{1} / __b);
      auto __h(__d);
      std::size_t __i = 1;
      while (true)
	{
	  auto __a = _M_num(__i, __x);
	  __b = _M_den(__i, __x);
	  __d = _Tp{1} / (__a * __d + __b);
	  __c = __b + __a / __c;
	  auto __del = __c * __d;
	  __h *= __del;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    break;
	  if (__i > _S_max_iter)
	    throw std::runtime_error("_LentzContinuedFraction: " "continued fraction evaluation failed");
	    //std::__throw_runtime_error(__N"_LentzContinuedFraction: "
		//		   "continued fraction evaluation failed"));
	  ++__i;
	}
      __h += _M_den(0, __x);
      return __h;
    }
  };

