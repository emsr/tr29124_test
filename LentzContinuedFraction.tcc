
/**
 * You need a functor belching a, b
 * an iterator pair for a and an iterator for b
 *  
 */
template<typename _Tp>
  class _LentzContinuedFraction
  {
    _Tp
    operator()
      std::complex<_Tp> __b(_Tp{1}, __t);
      std::complex<_Tp> __c(_Tp{1} / _S_fp_min);
      auto __d(_Tp{1} / __b);
      auto __h(__d);
      int __i = 2;
      while (true)
	{
	  auto __a = -_Tp(__i - 1) * _Tp(__i - 1);
	  __b += _Tp{2};
	  __d = _Tp{1} / (__a * __d + __b);
	  __c = __b + __a / __c;
	  auto __del = __c * __d;
	  __h *= __del;
	  if (std::abs(__del - _Tp{1}) < _S_eps)
	    break;
	  if (__i > _S_max_iter)
	    std::__throw_runtime_error(__N"_LentzContinuedFraction: "
				   "continued fraction evaluation failed"));
	  ++__i;
	}
      __h *= std::polar(_Tp{1}, -__t);
      _Ci = -__h.real();
      _Si = _S_pi_2 + __h.imag();
  };

