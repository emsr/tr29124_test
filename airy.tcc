#include <complex>
#include <ext/cmath>

/**
 *  @param[in]  t      The input argument.
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
template<typename _Tp>
  void
  airy(std::complex<_Tp> __t,
       std::complex<_Tp>& _Ai, std::complex<_Tp>& _Aip,
       std::complex<_Tp>& _Bi, std::complex<_Tp>& _Bip,
       std::complex<_Tp>& __w1, std::complex<_Tp>& __w1p,
       std::complex<_Tp>& __w2, std::complex<_Tp>& __w2p)
  {
    using __cmplx = std::complex<_Tp>;

    constexpr _Tp _S_pi = __gnu_cxx::__math_constants<_Tp>::__pi;
    constexpr _Tp _S_sqrt_pi = __gnu_cxx::__math_constants<_Tp>::__root_pi;
    constexpr _Tp _S_eps = __gnu_cxx::__math_constants<_Tp>::__eps;
    constexpr _Tp _S_log10min = std::log10(__gnu_cxx::__math_constants<_Tp>::__min);
    constexpr _Tp _S_Ai0{3.550280538878172392600631860041831763980e-1};
    constexpr _Tp _S_Aip0{2.588194037928067984051835601892039634793e-1};
    constexpr _Tp _S_Bi0{6.149266274460007351509223690936135535960e-1};
    constexpr _Tp _S_Bip0{8.868776642045783582807775119976424596506e-1};
    constexpr auto _S_i = __cmplx(_Tp{0}, _Tp{1});
    constexpr _Tp _S_big = _Tp{3.5};
    constexpr _Tp _F_k[9]
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
    /*constexpr*/ _Tp _Fp_k[9];
    if (true)
      for (int __k = 0; __k < 9; ++__k)
	_Fp_k[__k] = 3 * (__k + 1) * _F_k[__k];
    constexpr _Tp _G_k[9]
     {
      2.0 / 24.0,
      10.0 / 5040.0,
      80.0 / 36288.0e2,
      880.0 / 622702.0e4,
      880.0 * 14.0 / 209228.0e8,
      880.0 * 238.0 / 121645.0e12,
      880.0 * 4780.0 / 1.1240007278e21,
      880.0 * 4780.0 * 23.0 / 1.551121e25,
      880.0 * 4780.0 * 598.0 / 3.048883446e29
    };
    /*constexpr*/ _Tp _Gp_k[9];
    if (true)
      for (int __k = 0; __k < 9; ++__k)
	_Gp_k[__k] = (3 * (__k + 1) + 1) * _G_k[__k];

    _Tp _S_cn[20], _S_dn[20];
    if (true)
    {
      _S_cn[0] = _Tp{15} / _Tp{216};
      _S_dn[0] = -(_Tp{7} / _Tp{5}) * _S_cn[0];
      for (int __n = 2; __n <= 11; ++__n)
	{
	  _S_cn[__n - 1] = _S_cn[__n - 2]
			 * (6 * __n - 5) * (6 * __n - 3) * (6 * __n - 1)
			 / (216 * __n * (2 * __n - 1));
	  _S_dn[__n - 1] = -_S_cn[__n - 1] * (6 * __n + 1) / (6 * __n - 1);
	}
    }

    if (std::abs(__t) <= _S_big)
      {
	const auto __log10t = std::log10(std::abs(__t));
	const auto __ttt = __t * __t * __t;

	auto __term = __cmplx{1};
	auto _F = __cmplx{1};
	auto _G = __t;
	for (int __n = 0; __n < 9; ++__n)
	  {
	    if (std::abs(__t) < _S_eps)
	      break;
	    auto __xx = __log10t * (3 * (__n + 1) + 1)
		      + std::log10(_G_k[__n]);
	    if (__xx < _S_log10min)
	      break;
	    __term *= __ttt;
	    _F += _F_k[__n] * __term;
	    _G += _G_k[__n] * __term * __t;
	  }
	auto _U = std::sqrt(_Tp{3} * _S_pi)
		* (_S_Ai0 * _F + _S_Aip0 * _G);
	auto _V = _S_sqrt_pi * (_S_Ai0 * _F - _S_Aip0 * _G);
	__w1 = _U - _S_i * _V;
	__w2 = _U + _S_i * _V;
	_Bi = _U / _S_sqrt_pi;
	_Ai = _V / _S_sqrt_pi;

	__term = __cmplx{1};
	auto _Fp = __cmplx{0};
	auto _Gp = __cmplx{1};
	for (int __n = 0; __n < 9; ++__n)
	  {
	    if (std::abs(__t) < _S_eps)
	      break;
	    auto __xx = __log10t * 3 * (__n + 1)
		      + std::log10(_Gp_k[__n]);
	    if (__xx < _S_log10min)
	      break;
	    __term *= __ttt;
	    _Fp += _Fp_k[__n] * __term / __t;
	    _Gp += _Gp_k[__n] * __term;
	  }
	auto _Up = std::sqrt(_Tp{3} * _S_pi)
		 * (_S_Ai0 * _Fp + _S_Aip0 * _Gp);
	auto _Vp = _S_sqrt_pi * (_S_Ai0 * _Fp - _S_Aip0 * _Gp);
	__w1p = _Up - _S_i * _Vp;
	__w2p = _Up + _S_i * _Vp;
	_Bip = _Up / _S_sqrt_pi;
	_Aip = _Vp / _S_sqrt_pi;

	return;
      }
    else // |t| > 3.5
      {
	if (std::real(__t) > _Tp{0})
	  {
	    auto __zeta0 = (_Tp{2} / _Tp{3}) * std::pow(__t, 1.5);
	    auto __mqrt0 = std::pow(__t, -0.25);
	    auto __pqrt0 = std::pow(__t, +0.25);
	    auto __ezeta0 = std::exp(-__zeta0);
	    _Ai = __cmplx{1};
	    _Aip = __cmplx{1};
	    auto __fact0 = -_Tp{1} / __zeta0;
	    auto __izeta0 = __cmplx{1};
	    for (int __n = 0; __n < 11; ++__n)
	      {
		__izeta0 *= __fact0;
		_Ai += _S_cn[__n] * __izeta0;
		_Aip += _S_dn[__n] * __izeta0;
	      }
	    _Ai *= 0.5 * __mqrt0 * __ezeta0 / _S_sqrt_pi;
	    _Aip *= -0.5 * __pqrt0 * __ezeta0 / _S_sqrt_pi;

	    auto __t1 = __t * std::exp(+_Tp{2} * _S_pi * _S_i / _Tp{3});
	    auto __t2 = __t * std::exp(-_Tp{2} * _S_pi * _S_i / _Tp{3});
	    auto __zeta1 = (_Tp{2} / _Tp{3}) * std::pow(__t1, 1.5);
	    auto __zeta2 = (_Tp{2} / _Tp{3}) * std::pow(__t2, 1.5);
	    auto __mqrt1 = std::pow(__t1, -0.25);
	    auto __mqrt2 = std::pow(__t2, -0.25);
	    auto __pqrt1 = std::pow(__t1, +0.25);
	    auto __pqrt2 = std::pow(__t2, +0.25);
	    auto __ezeta1 = std::exp(-__zeta1);
	    auto __ezeta2 = std::exp(-__zeta2);
	    auto _Ai1 = __cmplx{1};
	    auto _Ai1p = __cmplx{1};
	    auto _Ai2 = _Ai1;
	    auto _Ai2p = _Ai1p;
	    auto __sign = 1;
	    auto __izeta1 = __cmplx{1};
	    auto __izeta2 = __cmplx{1};
	    for (int __n = 0; __n < 11; ++__n)
	      {
		__sign *= -1;
		__izeta1 /= __zeta1;
		__izeta2 /= __zeta2;
		_Ai1 += __sign * _S_cn[__n] * __izeta1;
		_Ai2 += __sign * _S_cn[__n] * __izeta2;
		_Ai1p += __sign * _S_dn[__n] * __izeta1;
		_Ai2p += __sign * _S_dn[__n] * __izeta2;
	      }
	    _Ai1 *= 0.5 * __ezeta1 * __mqrt1 / _S_sqrt_pi;
	    _Ai2 *= 0.5 * __ezeta2 * __mqrt2 / _S_sqrt_pi;
	    _Ai1p *= -0.5 * __pqrt1 * __ezeta1 / _S_sqrt_pi;
	    _Ai2p *= -0.5 * __pqrt2 * __ezeta2 / _S_sqrt_pi;

	    _Bi = std::exp(+_S_i * _S_pi / _Tp{6}) * _Ai1
		+ std::exp(-_S_i * _S_pi / _Tp{6}) * _Ai2;
	    _Bip = std::exp(+_S_i * _Tp{5} * _S_pi / _Tp{6}) * _Ai1p
		 + std::exp(-_S_i * _Tp{5} * _S_pi / _Tp{6}) * _Ai2p;

	    __w1 = _S_sqrt_pi * (_Bi - _S_i * _Ai);
	    __w2 = _S_sqrt_pi * (_Bi + _S_i * _Ai);
	    __w1p = _S_sqrt_pi * (_Bip - _S_i * _Aip);
	    __w2p = _S_sqrt_pi * (_Bip + _S_i * _Aip);

	    return;
	  }
	else // Argument t is on or left of the imaginary axis.
	  {
	    auto __zeta = (_Tp{2} / _Tp{3}) * std::pow(-__t, 1.5);
	    auto __mqrt = std::pow(-__t, -0.25);
	    auto __pqrt = std::pow(-__t, +0.25);
	    auto __mezeta = std::exp(-_S_i * (__zeta + (_S_pi / _Tp{4})));
	    auto __pezeta = std::exp(+_S_i * (__zeta + (_S_pi / _Tp{4})));
	    __w1 = __cmplx{1};
	    __w2 = __cmplx{1};
	    __w1p = +_S_i;
	    __w2p = -_S_i;
	    auto __ipn = __cmplx{1};
	    auto __imn = __cmplx{1};
	    auto __ixn = __cmplx{1};
	    for (int __n = 0; __n < 11; ++__n)
	      {
		__ipn *= +_S_i;
		__imn *= -_S_i;
		__ixn /= __zeta;
		__w1 += __ipn * _S_cn[__n] * __ixn;
		__w2 += __imn * _S_cn[__n] * __ixn;
		__w1p += +_S_i * __ipn * _S_dn[__n] * __ixn;
		__w2p += -_S_i * __imn * _S_dn[__n] * __ixn;
	      }
	    __w1 *= __mqrt * __mezeta;
	    __w2 *= __mqrt * __pezeta;
	    __w1p *= __pqrt * __mezeta;
	    __w2p *= __pqrt * __pezeta;

	    _Bi = (__w1 + __w2) / (_Tp{2} * _S_sqrt_pi);
	    _Ai = (__w2 - __w1) / (_Tp{2} * _S_i * _S_sqrt_pi);
	    _Bip = (__w1p + __w2p) / (_Tp{2} * _S_sqrt_pi);
	    _Aip = (__w2p - __w1p) / (_Tp{2} * _S_i * _S_sqrt_pi);

	    return;
	  }
      }
  }
