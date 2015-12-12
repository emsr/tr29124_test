/**  @param[in]  T  
 *
 *  @param[out]  _Ai
 *  @param[out]  _Aip  
 *  @param[out]  _Bi
 *  @param[out]  _Bip  
 *  @param[out]  __w1
 *  @param[out]  __w1p  
 *  @param[out]  __w2
 *  @param[out]  __w2p  
 */
  airy(std::complex<double> __t,
       std::complex<double>& _Ai, std::complex<double>& _Aip,
       std::complex<double>& _Bi, std::complex<double>& _Bip,
       std::complex<double>& __w1, std::complex<double>& __w1p,
       std::complex<double>& __w2, std::complex<double>& __w2p)
  {
    using __cmplx = std::complex<double>;

    __cmplx _F,_G,_U,_V,_Fp,_Gp,_Up,_Vp,b1,b2,b3,b4,b5,
            b6,b7,b8,b9,b10,_Ai1,_Aip1,_Ai2,_Aip2,i,x,x1,x2,x3,x4;
    bool once{false};
    constexpr double _S_pi{3.1415926535e0};
    constexpr double _S_gamma1d3{2.588194037928068e-01},
    constexpr double _S_gamma2d3{3.550280538878172e-01},
    constexpr double _F_k[9]
    {
      1.0 / 6.0,
      4.0 / 720.0,
      28.0 / 362880.0,
      28.0 / 47900160.0,
      280.0 * 13.0 / 130767.0e7,
      280.0 * 208.0 / 640237.0e10,
      280.0 * 208.0 * 19.0 / 5.1090942e19,
      280.0 * 208.0 * 418.0 / 6.2044840173e23,
      280.0 * 208.0 * 418.0 * 25.0 / 1.0888869e28
    };
    constexpr double _Fp_k[9];
    constexpr double _G_k[9]
     {
      2.0 / 24.0,
      10.0 / 5040.0,
      80.0 / 36288.0e2,
      880.0 / 622702.0e4,
      880.0 * 14.0 / 209228.0e8,
      880.0 * 238.0 / 121645.0e12,
      880.0 * 4780.0 / 1.1240007278e21,
      880.0 * 4780.0 * 23.0D0 / 1.551121e25,
      880.0 * 4780.0 * 598.0D0 / 3.048883446e29
    };
    constexpr double _Gp_k[9];
    constexpr auto _S_i = __cmplx(0.0, 1.0);
    constexpr double _S_big = 3.5;
    double cn[20], dn[20];
    if (once)
    {
      for (int __n = 1; __n <= 9; ++__n)
	{
	  _Fp_k[__n] = 3 * __n * _F_k[__n];
          _Gp_k[__n] = (3 * __n + 1) * _G_k[__n];
	}
      const auto _S_sqrtpi = std::sqrt(_S_pi);
      cn[1] = 15.0 / 216.0;
      dn[1] = -(7.0 / 5.0) * cn[1];
      for (int __n = 2; __n <= 11; ++__n)
	{
	  auto __m = __n - 1;
	  cn[__n] = cn[__m] * (6 * __n - 5) * (6 * __n - 3) * (6 * __n - 1)
		/ (216 * __n * (2 * __n - 1));
	  dn[__n] = -cn[__n] * (6 * __n + 1) / (6 * __n - 1);
	}
      once = true;
    }
    if (std::abs(__t) <= _S_big)
      {
	auto _F = __cmplx{1};
	auto _G = t;
	for (int __n = 1; __n <= 9; ++__n)
	  {
	    if (std::abs(__t) < 1.0e-10)
	      break;
	    auto __xx = std::log10(std::abs(t)) * (3 * __n + 1) + std::log10(_G_k[__n]);
	    if (__xx < -36.0)
	      break;
	    _F += std::pow(_F_k[__n] * __t, 3 * __n);
	    _G += std::pow(_G_k[__n] * __t, 3 * __n + 1);
	  }
	auto _U = std::sqrt(3.0 * _S_pi) * (_S_gamma2d3 * _F + _S_gamma1d3 * _G);
	auto _V = _S_sqrtpi * (_S_gamma2d3 * _F - _S_gamma1d3 * _G);
	__w1 = _U - _S_i * _V;
	__w2 = _U + _S_i * _V;
	_Bi = _U / _S_sqrtpi;
	_Ai = _V / _S_sqrtpi;
	_Fp = __cmplx{};
	_Gp = __cmplx{1};
	for (int __n = 1; __n <= 9; ++__n)
	  {
	    if (std::abs(__t) < 1.0e-10)
	      break;
	    auto __xx = std::log10(std::abs(__t)) * 3 * __n + std::log10(_Gp_k[__n]);
	    if (__xx < -36.0)
	      break;
	    _Fp += std::pow(_Fp_k[__n] * __t, 3 * __n - 1);
	    _Gp += std::pow(_Gp_k[__n] * __t, 3 * __n);
	  }
	_Up = std::sqrt(3.0 * _S_pi) * (_S_gamma2d3 * _Fp + _S_gamma1d3 * _Gp);
	_Vp = _S_sqrtpi * (_S_gamma2d3 * _Fp - _S_gamma1d3 * _Gp);
	__w1p = _Up - _S_i * _Vp;
	__w2p = _Up + _S_i * _Vp;
	_Bip = _Up / _S_sqrtpi;
	_Aip = _Vp / _S_sqrtpi;
	return;
      }
    else // |t| > 3.5
      {
	if (std::real(__t) > 0.0)
	  {
	    auto __zeta0 = std::pow((2.0 / 3.0) * __t, 1.5);
	    auto __qrt0 = std::pow(__t, -0.25);
	    auto __iqrt0 = std::pow(__t, 0.25);
	    auto __ezeta0 = std::exp(-__zeta0);
	    _Ai = 1.0;
	    _Aip = 1.0;
	    auto __sign0 = -1;
	    auto __izeta0 = 1 / __zeta0;
	    for (int n = 1; n <= 11; ++n)
	      {
		_Ai += __sign0 * cn[__n] * __izeta0;
		_Aip += __sign0 * dn[__n] * __izeta0;
		__sign0 *= -1;
		__izeta0 /= __zeta0;
	      }
	    _Ai *= 0.5 * __qrt0 * __ezeta0 / _S_sqrtpi;
	    _Aip *= -0.5 * __iqrt0 * __ezeta0 / _S_sqrtpi;
	    auto __t1 = __t * std::exp(+2.0 * _S_pi * _S_i / 3.0);
	    auto __t2 = __t * std::exp(-2.0 * _S_pi * _S_i / 3.0);
	    auto __zeta1 = std::pow((2.0 / 3.0) * __t1, 1.5);
	    auto __zeta2 = std::pow((2.0 / 3.0) * __t2, 1.5);
	    auto __qrt1 = std::pow(__t1, -0.25);
	    auto __qrt2 = std::pow(__t2, -0.25);
	    auto __iqrt1 = std::pow(__t1, 0.25);
	    auto __iqrt2 = std::pow(__t2, 0.25);
	    auto __ezeta1 = std::exp(-__zeta1);
	    auto __ezeta2 = std::exp(-__zeta2);
	    _Ai1 = 1.0;
	    _Aip1 = 1.0;
	    _Ai2 = _Ai1;
	    _Aip2 = _Aip1;
	    auto __sign = -1;
	    auto __izeta1 = 1 / __zeta1;
	    auto __izeta2 = 1 / __zeta2;
	    for (int __n = 1; __n <= 11; ++__n)
	      {
		_Ai1 += __sign * cn[__n] * __izeta1;
		_Ai2 += __sign * cn[__n] * __izeta2;
		_Aip1 += __sign * dn[__n] * __izeta1;
		_Aip2 += __sign * dn[__n] * __izeta2;
		__sign *= -1;
		__izeta1 /= __zeta1;
		__izeta2 /= __zeta2;
	      }
	    _Ai1 *= 0.5 * __ezeta1 * __qrt1 / _S_sqrtpi;
	    _Ai2 *= 0.5 * __ezeta2 * __qrt2 / _S_sqrtpi;
	    _Aip1 *= -0.5 * __iqrt1 * __ezeta1 / _S_sqrtpi;
	    _Aip2 *= -0.5 * __iqrt2 * __ezeta2 / _S_sqrtpi;
	    _Bi = std::exp(+_S_i * _S_pi / 6.0) * _Ai1
		+ std::exp(-_S_i * _S_pi / 6.0) * _Ai2;
	    _Bip = std::exp(+_S_i * 5.0 * _S_pi / 6.0) * _Aip1
		 + std::exp(-_S_i * 5.0 * _S_pi / 6.0) * _Aip2;
	    __w1 = _S_sqrtpi * (_Bi - _S_i * _Ai);
	    __w2 = _S_sqrtpi * (_Bi + _S_i * _Ai);
	    __w1p = _S_sqrtpi * (_Bip - _S_i * _Aip);
	    __w2p = _S_sqrtpi * (_Bip + _S_i * _Aip);
	    return;
	  }
	else // Argument t is on or below real axis.
	  {
	    auto __zeta = std::pow(-(2.0 / 3.0) * __t, 1.5);
	    auto __qrt0 = std::pow(-__t, -0.25);
	    auto __iqrt0 = std::pow(-__t, +0.25);
	    auto __ezeta0 = std::exp(-_S_i * (__zeta + (_S_pi / 4.0)));
	    auto __iezeta0 = std::exp(+_S_i * (__zeta + (_S_pi / 4.0)));
	    __w1 = __cmplx{1};
	    __w2 = __cmplx{1};
	    __w1p = +_S_i;
	    __w2p = -_S_i;
	    auto __ipn = +_S_i;
	    auto __imn = -_S_i;
	    auto __ixn = 1 / __zeta;
	    for (int __n = 1; __n <= 11; ++__n)
	      {
		__w1 += __ipn * cn[__n] * __ixn;
		__w2 += __imn * cn[__n] * __ixn;
		__w1p += +_S_i * __ipn * dn[__n] * __ixn;
		__w2p += -_S_i * __imn * dn[__n] * __ixn;
		__ipn *= +_S_i;
		__imn *= -_S_i;
		__ixn /= __zeta;
	      }
	    __w1 *= __qrt0 * __ezeta0;
	    __w2 *= __qrt0 * __iezeta0;
	    __w1p *= __iqrt0 * __ezeta0;
	    __w2p *= __iqrt0 * __iezeta0;
	    _Bi = (__w1 + __w2) / (2.0 * _S_sqrtpi);
	    _Ai = (__w2 - __w1) / (2.0 * _S_i * _S_sqrtpi);
	    _Bip = (__w1p + __w2p) / (2.0 * _S_sqrtpi);
	    _Aip = (__w2p - __w1p) / (2.0 * _S_i * _S_sqrtpi);
	    return;
	  }
      }
  }
