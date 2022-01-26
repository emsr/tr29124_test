

  template <typename _Tp>
    _Tp
    legendre_p(unsigned int l, _Tp x)
    {
      if (l < 0)
	throw std::domain_error("legendre_p: Bad order");
      if (x < _Tp{-1} || x > _Tp{+1})
	throw std::domain_error("legendre_p: Bad argument");

      auto _P_lm2 = _Tp{1};
      if (l == 0)
	return _P_lm2;

      auto _P_lm1 = x;
      if (l == 1)
	return _P_lm1;

      if (x == _Tp{+1})
	return _Tp{+1};

      if (x == _Tp{-1})
	return (l % 2 == 1 ? _Tp{-1} : _Tp{+1});

      _Tp _P_l;
      for (auto ll = 2u; ll <= l; ++ll)
	{
          _P_l = (_Tp(2 * ll - 1) * x * _P_lm1
        	- _Tp(ll - 1) * _P_lm2) / _Tp(ll);
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
    std::pair<_Tp, _Tp>
    lobatto(unsigned int l, _Tp x)
    {
      if (x == _Tp{-1} || x == _Tp{+1})
	return std::make_pair(_Tp{0}, _Tp{0}); // FIXME deriv!
      else
	{
	  auto _Lo_lm2 = _Tp{0};
	  if (l == 0)
	    return std::make_pair(_Lo_lm2, _Tp{0});

	  auto _Lo_lm1 = _Tp{1} - x * x;
	  if (l == 1)
	    return std::make_pair(_Lo_lm1, -_Tp{2} * x);

	  auto _P_lm2 = _Tp{1};
	  auto _P_lm1 = x;
	  _Tp _P_l;
	  auto _Pp_lm1 = _Tp{0};
	  auto _Pp_l = _Tp{1};
	  _Tp _Lo_l;
	  _Tp _Lop_l;
	  for (auto ll = 2u; ll <= l; ++ll)
	    {
	      _P_l = (_Tp(2 * ll - 1) * x * _P_lm1
		   - _Tp(ll - 1) * _P_lm2) / _Tp(ll);
	      // Recursion for the derivative of the Legendre polynomial.
	      _Pp_lm1 = _Pp_l;
	      _Lo_l = -l * (x * _P_l - _P_lm1);
	      _Pp_l = _Lo_l / (_Tp{1} - x * x);
	      _Lop_l = l * (_Pp_lm1 - _P_l - x * _Pp_l);
	      _P_lm2 = _P_lm1;
	      _P_lm1 = _P_l;
	    }
	  return std::make_pair(_Lo_l, _Lop_l);
	}
    }
