#include <complex>
#include <ext/cmath>
#include "numeric_limits.h"

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
    constexpr _Tp _S_eps = __gnu_cxx::__epsilon(_Tp{});
    constexpr _Tp _S_log10min = __gnu_cxx::__log10_min(_Tp{});
    constexpr _Tp _S_Ai0{3.550280538878172392600631860041831763980e-1};
    constexpr _Tp _S_Aip0{2.588194037928067984051835601892039634793e-1};
    //constexpr _Tp _S_Bi0{6.149266274460007351509223690936135535960e-1};
    //constexpr _Tp _S_Bip0{8.868776642045783582807775119976424596506e-1};
    constexpr auto _S_i = __cmplx(_Tp{0}, _Tp{1});
    constexpr _Tp _S_big = _Tp{3.5};
/*
    constexpr _Tp _F_k[9]
    {
      1.0 / 6.0,       //  1 / 3!==6
      4.0 / 720.0,     //  1*4 / 6!==720
      28.0 / 362880.0, //  1*4*7 / 9!==362880
      28.0 / 47900160.0,  //  (1*4*7*10 / 10) / (12! / 10)==47900160
      280.0 * 13.0 / 130767.0e7,  //  1*4*7*10*13 / (15!==1307674368000)
      280.0 * 208.0 / 640237.0e10,  //  1*4*7*10*13*16 / 18!==6402373705728000
      280.0 * 208.0 * 19.0 / 5.1090942e19,  //  1*4*7*10*13*16*19 / 21!==51090942171709440000
      280.0 * 208.0 * 418.0 / 6.2044840173e23,  //  1*4*7*10*13*16*19*22 / 24!==620448401733239439360000
      280.0 * 208.0 * 418.0 * 25.0 / 1.0888869e28  //  1*4*7*10*13*16*19*22*25 / 27!==10888869450418352160768000000
    };
    constexpr _Tp _Fp_k[9];
    if (true)
      for (int __k = 0; __k < 9; ++__k)
	_Fp_k[__k] = 3 * (__k + 1) * _F_k[__k];
    constexpr _Tp _G_k[9]
     {
      2.0 / 24.0,  //  2 / 4!
      10.0 / 5040.0,  //  2*5 / 7!
      80.0 / 36288.0e2,  //  2*5*8 / 10!==3628800
      880.0 / 622702.0e4,  //  2*5*8*11 / 13!==6227020800
      880.0 * 14.0 / 209228.0e8,  //  2*5*8*11*14 / 16!==20922789888000
      880.0 * 238.0 / 121645.0e12,  //  2*5*8*11*14*17 / 19!==121645100408832000
//            4760.0 was 4780.0!
      880.0 * 4760.0 / 1.1240007278e21,  //  2*5*8*11*14*17*20 / 22!==1124000727777607680000
      880.0 * 4760.0 * 23.0 / 1.551121e25,  //  2*5*8*11*14*17*20*23 / 25!==15511210043330985984000000
      880.0 * 4760.0 * 598.0 / 3.048883446e29  //  2*5*8*11*14*17*20*23*26 / 28!==304888344611713860501504000000
    };
    constexpr _Tp _Gp_k[9];
    if (true)
      for (int __k = 0; __k < 9; ++__k)
	_Gp_k[__k] = (3 * (__k + 1) + 1) * _G_k[__k];
*/
    constexpr int _N_FG = 12;
    constexpr _Tp
    _F_k[_N_FG]
    {
      0.166666666666666666671L,
      0.0055555555555555555554L,
      7.71604938271604938293e-05L,
      5.84549195660306771413e-07L,
      2.78356759838241319725e-09L,
      9.09662613850461829133e-12L,
      2.16586336631062340282e-14L,
      3.9236655186786655847e-17L,
      5.58926712062487975073e-20L,
      6.42444496623549396612e-23L,
      6.08375470287452080124e-26L,
      4.82837674831311174686e-29L
    };
    constexpr _Tp
    _Fp_k[_N_FG]
    {
      1.0L,
      0.0499999999999999999973L,
      0.000925925925925925925952L,
      8.76823793490460157129e-06L,
      5.01042167708834375494e-08L,
      1.91029148908596984119e-10L,
      5.19807207914549616664e-13L,
      1.05938969004323970788e-15L,
      1.67678013618746392524e-18L,
      2.1200668388577130088e-21L,
      2.19015169303482748847e-24L,
      1.88306693184211358127e-27L
    };
    constexpr _Tp
    _G_k[_N_FG]
    {
      0.0833333333333333333356L,
      0.00198412698412698412696L,
      2.20458553791887125215e-05L,
      1.41319585764030208476e-07L,
      5.88831607350125868627e-10L,
      1.72172984605299961594e-12L,
      3.72668797846969613851e-15L,
      6.21114663078282689725e-18L,
      8.21580242167040594874e-21L,
      8.83419615233376983785e-24L,
      7.87361510903188042579e-27L,
      5.9111224542281384578e-30L
    };
    constexpr _Tp
    _Gp_k[_N_FG]
    {
      0.583333333333333333369L,
      0.0198412698412698412696L,
      0.000286596119929453262791L,
      2.26111337222448333562e-06L,
      1.11878005396523915037e-08L,
      3.78780566131659915511e-11L,
      9.31671994617424034639e-14L,
      1.73912105661919153125e-16L,
      2.54689875071782584419e-19L,
      3.00362669179348174486e-22L,
      2.9132375903417957575e-25L,
      2.36444898169125538315e-28L
    };

    constexpr int _N_cd = 20;
    _Tp _S_cn[_N_cd], _S_dn[_N_cd];
    if (true)
    {
      _S_cn[0] = _Tp{15} / _Tp{216};
      _S_dn[0] = -(_Tp{7} / _Tp{5}) * _S_cn[0];
      for (int __n = 2; __n < 11; ++__n)
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
	for (int __n = 0; __n < _N_FG; ++__n)
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
	auto _UU = std::sqrt(_Tp{3} * _S_pi)
		* (_S_Ai0 * _F + _S_Aip0 * _G);
	auto _VV = _S_sqrt_pi * (_S_Ai0 * _F - _S_Aip0 * _G);
	_Bi = _UU / _S_sqrt_pi;
	_Ai = _VV / _S_sqrt_pi;
	__w1 = _UU - _S_i * _VV;
	__w2 = _UU + _S_i * _VV;

	__term = __cmplx{1};
	auto _Fp = __cmplx{0};
	auto _Gp = __cmplx{1};
	for (int __n = 0; __n < _N_FG; ++__n)
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
	auto _UUp = std::sqrt(_Tp{3} * _S_pi)
		  * (_S_Ai0 * _Fp + _S_Aip0 * _Gp);
	auto _VVp = _S_sqrt_pi * (_S_Ai0 * _Fp - _S_Aip0 * _Gp);
	_Bip = _UUp / _S_sqrt_pi;
	_Aip = _VVp / _S_sqrt_pi;
	__w1p = _UUp - _S_i * _VVp;
	__w2p = _UUp + _S_i * _VVp;

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
