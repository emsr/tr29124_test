

  template<typename _Tp>
    _Tp
    get_zeta(_Tp __y)
    {
      static constexpr auto _S_2d3   = _Tp{0.6666666666666666666666666666666666666667Q};
      static constexpr auto _S_lncon = _Tp{0.2703100720721095879853420769762327577152Q}; // -(2/3)ln(2/3)
      if (__y == _Tp{0})
	return std::numeric_limits<_Tp>::infinity();
      else if (__y <= _Tp{1})
	{
	  auto __w = std::sqrt((_Tp{1} + __y) * (_Tp{1} - __y));
	  // Compute xi = ln(1 + (1 - y^2)^(1/2)) - ln(y) - (1 - y^2)^(1/2) = (2/3)(zeta)^(3/2)
	  // using default branch of logarithm and square root.
	  auto __xi = std::log(_Tp{1} + __w) - std::log(__y) - __w;

	  auto __logxi = std::log(__xi);

	  // Compute ln(zeta), zeta.
	  auto __logzeta = _S_2d3 * __logxi + _S_lncon;
	  auto __zeta = std::exp(__logzeta);
	  return __zeta;
	}
      else
	{
	  auto __w = std::sqrt((__y + _Tp{1}) * (__y - _Tp{1}));
	  // Compute xi = (y^2 - 1)^(1/2) - arcsec(y) = (2/3)(-zeta)^(3/2)
	  // using default branch of logarithm and square root.
	  auto __xi = __w - std::acos(_Tp{1} / __y);

	  auto __logxi = std::log(__xi);

	  // Compute ln(-zeta), zeta.
	  auto __logmzeta = _S_2d3 * __logxi + _S_lncon;
	  auto __zeta = -std::exp(__logmzeta);
	  return __zeta;
	}
    }
