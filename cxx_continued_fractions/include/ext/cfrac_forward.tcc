#ifndef CFRAC_FORWARD_TCC
#define CFRAC_FORWARD_TCC 1

#include <stdexcept>
#include <bits/numeric_limits.h>

/**
 * You need a functors returning the numerator and denominator a, b
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
    using _Val = emsr::fp_promote_t<_Tp, _ARet, _BRet>;
    using _Real = emsr::num_traits_t<_Val>;

    const _Real _S_fp_min = std::numeric_limits<_Real>::min();

    _Real _M_rel_error = std::numeric_limits<_Real>::epsilon();
    int _M_max_iter = 1000;

    _ForwardContinuedFraction(_AFun __a, _BFun __b, _TailFun __w)
    : _M_num(__a),
      _M_den(__b),
      _M_tail(__w)
    { }

    _Ret
    operator()(_Tp __x) const
    {
      auto __b = _M_den(0, __x);
      auto _Den2 = __b;
      auto __a = _M_num(1, __x);
      auto _Num2 = __a;
      __b = _M_den(2, __x);
      auto _Num1 = __b * _Num2;
      __a = _M_num(2, __x);
      auto _Den1 = __a + __b * _Den2;
      std::size_t __k = 2;
      while (true)
	{
	  __a = _M_num(__k, __x);
	  __b = _M_den(__k, __x);
	  auto _Num = __b * _Num1 + __a * _Num2;
	  auto _Den = __b * _Den1 + __a * _Den2;
	  _Num2 = _Num1;
	  _Den2 = _Den1;
	  _Num1 = _Num;
	  _Den1 = _Den;
	  ++__k;
	  if (__k > _M_max_iter)
	    throw std::runtime_error("_ForwardContinuedFraction: "
				     "continued fraction evaluation failed");
	}
    }
  };

#endif // CFRAC_FORWARD_TCC
