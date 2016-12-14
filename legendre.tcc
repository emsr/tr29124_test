

  template <typename _Tp>
    _Tp
    __legendre_p(unsigned int __l, _Tp __x)
    {
      if (__l < 0)
	throw std::domain_error("__legendre_p: Bad order");
      if (__x < _Tp{-1} || __x > _Tp{+1})
	throw std::domain_error("__legendre_p: Bad argument");

      auto _P_lm2 = _Tp{1};
      if (__l == 0)
	return _P_lm2;

      auto _P_lm1 = __x;
      if (__l == 1)
	return _P_lm1;

      if (__x == _Tp{+1})
	return _Tp{+1};

      if (__x == _Tp{-1})
	return (__l % 2 == 1 ? _Tp{-1} : _Tp{+1});

      _Tp _P_l;
      for (auto __ll = 2u; __ll <= __l; ++__ll)
	{
          _P_l = (_Tp(2 * __ll - 1) * __x * _P_lm1
        	- _Tp(__ll - 1) * _P_lm2) / _Tp(__ll);
          _P_lm2 = _P_lm1;
          _P_lm1 = _P_l;
	}

      return _P_l;
    }

/**
 * Return the Lobatto polynomial of order @f$ l @f$ and argument @f$ x @f$.
 * The Lobatto polynomials are closely related to Legendre polynomials
 * but have the distinction that nodes computed from them include the endpoints
 * at [-1,+1].
 * @f[
 *   Lo_l(x) = (1 - x^2)P'_l(x) = l \left[ P_{l - 1}(x) - x P_l(x) \right]
 * @f[
 */
  template <typename _Tp>
    _Tp
    __lobatto(unsigned int __l, _Tp __x)
    {
      if (__x == _Tp{-1} || __x == _Tp{+1})
	return _Tp{0};
      else
	{
	  auto _Lo_lm2 = _Tp{0};
	  if (__l == 0)
	    return _Lo_lm2;

	  auto _Lo_lm1 = _Tp{1} - __x * __x;
	  if (__l == 1)
	    return _Lo_lm1;

	  auto _P_lm2 = _Tp{1};
	  auto _P_lm1 = __x;
	  _Tp _P_l;
	  auto _Pp_lm1 = _Tp{0};
	  auto _Pp_l = _Tp{1};
	  _Tp _Lo_l;
	  _Tp _Lop_l;
	  for (auto __ll = 2u; __ll <= __l; ++__ll)
	    {
	      _P_l = (_Tp(2 * __ll - 1) * __x * _P_lm1
		   - _Tp(__ll - 1) * _P_lm2) / _Tp(__ll);
	      // Recursion for the derivative of the Legendre polynomial.
	      _Pp_lm1 = _Pp_l;
	      _Lo_l = -__l * (__x * _P_l - _P_lm1);
	      _Pp_l = _Lo_l / (_Tp{1} - __x * __x);
	      _Lop_l = __l * (_Pp_lm1 - _P_l - __x * _Pp_l);
	      _P_lm2 = _P_lm1;
	      _P_lm1 = _P_l;
	    }
	  return _Lo_l;
	}
    }
