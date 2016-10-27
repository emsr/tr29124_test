// Special functions -*- C++ -*-

// Copyright (C) 2006-2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/sf_hermite.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

#ifndef _GLIBCXX_BITS_NEW_HERMITE_TCC
#define _GLIBCXX_BITS_NEW_HERMITE_TCC 1

#pragma GCC system_header

  /**
   *   @brief This routine returns the Hermite polynomial
   *          of order n: @f$ H_n(x) @f$ by recursion on n.
   *
   *   The Hermite polynomial is defined by:
   *   @f[
   *     H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   *   @f]
   *
   *   @param __n The order of the Hermite polynomial.
   *   @param __x The argument of the Hermite polynomial.
   *   @return The value of the Hermite polynomial of order n
   *           and argument x.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_recursion(unsigned int __n, _Tp __x)
    {
      // Compute H_0.
      auto __H_nm2 = _Tp{1};
      if (__n == 0)
	return __H_nm2;

      // Compute H_1.
      auto __H_nm1 = _Tp{2} * __x;
      if (__n == 1)
	return __H_nm1;

      // Compute H_n.
      _Tp __H_n;
      for (unsigned int __i = 2; __i <= __n; ++__i)
	{
	  __H_n = _Tp{2} * (__x * __H_nm1 - _Tp(__i - 1) * __H_nm2);
	  __H_nm2 = __H_nm1;
	  __H_nm1 = __H_n;
	}

      return __H_n;
    }


  /**
   *   @brief This routine returns the Hermite polynomial
   *          of large order n: @f$ H_n(x) @f$.  We assume here
   *          that x >= 0.
   *
   *   The Hermite polynomial is defined by:
   *   @f[
   *     H_n(x) = (-1)^n e^{x^2} \frac{d^n}{dx^n} e^{-x^2}
   *   @f]
   *
   *  @see "Asymptotic analysis of the Hermite polynomials
   *        from their differential-difference equation", 
   *        Diego Dominici, arXiv:math/0601078v1 [math.CA] 4 Jan 2006
   *   @param __n The order of the Hermite polynomial.
   *   @param __x The argument of the Hermite polynomial.
   *   @return The value of the Hermite polynomial of order n
   *           and argument x.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_asymp(unsigned int __n, _Tp __x)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      const auto _S_sqrt_2 = __gnu_cxx::__const_root_2(__x);
      const auto _S_sqrt_2pi = __gnu_cxx::__const_root_2_pi(__x);
      // __x >= 0 in this routine.
      const auto __xturn = std::sqrt(_Tp(2 * __n));
      if (std::abs(__x - __xturn) < _Tp{0.05L} * __xturn)
	{
	  // Transition region x ~ sqrt(2n).
	  const auto __n_2 = _Tp(__n) / _Tp{2};
	  const auto __n6th = std::pow(_Tp(__n), _Tp{1} / _Tp{6});
	  const auto __exparg = _Tp(__n) * std::log(__xturn) - _Tp{3} * __n_2
			      + __xturn * __x;
	  const auto __airyarg = _S_sqrt_2 * (__x - __xturn) * __n6th;
	  auto _Airy = std::__detail::__airy(__airyarg);
	  return _S_sqrt_2pi * __n6th * std::exp(__exparg) * _Airy.__Ai_value;
	}
      else if (__x < __xturn)
	{
	  // Oscillatory region |x| < sqrt(2n).
	  const auto __theta = std::asin(__x / __xturn);
	  const auto __2theta = _Tp{2} * __theta;
	  const auto __n_2 = _Tp(__n) / _Tp{2};
	  const auto __exparg = __n_2 * (_Tp{2} * std::log(__xturn)
					- std::cos(__2theta));
	  const auto __arg = __theta / _Tp{2}
			   + __n_2 * (std::sin(__2theta) + __2theta - _S_pi);
	  return std::sqrt(_Tp{2} / std::cos(__theta))
	       * std::exp(__exparg) * std::cos(__arg);
	}
      else
	{
	  // Exponential region |x| > sqrt(2n).
	  const auto __sigma = std::sqrt((__x - __xturn) * (__x + __xturn));
	  const auto __exparg = _Tp{0.5L} * (__x * (__x - __sigma) - __n)
			      + __n * std::log(__sigma + __x);
	  return std::exp(__exparg)
	       * std::sqrt(_Tp{0.5L} * (_Tp{1} + __x / __sigma));
	}
    }


  template<typename _Tp>
    _Tp
    __poly_hermite_asymp2(unsigned int __n, _Tp __x)
    {
      const auto _S_pi = __gnu_cxx::__const_pi(__x);
      const auto _S_sqrt_2 = __gnu_cxx::__const_root_2(__x);
      const auto _S_sqrt_pi = __gnu_cxx::__const_root_pi(__x);
      const auto _S_Ai0 = _Tp{-2.3381074104597670384891972524467L};
      bool __n_odd = __n % 2 == 1;
      auto __z = std::abs(__x);
      auto __z_turn = std::sqrt(_Tp(2 * __n + 1));
      auto __delta = _S_Ai0 / _S_sqrt_2 / std::pow(__n, _Tp{1} / _Tp{6});
      auto __sign = __n_odd && (__x < _Tp{0}) ? _Tp{-1} : _Tp{1};
      auto __f = _Tp{1};
      for (int __j = 1; __j <= __n; ++__j)
        __f *= std::sqrt(_Tp(__j));

      if (__z < __z_turn + __delta)
	{
	  auto __phi = std::acos(__z / __z_turn);
	  return __f * __sign
		     * (__n_odd ? _S_sqrt_2 : _Tp{1})
		     * std::pow(_Tp{2}, _Tp(__n) / _Tp{2})
		     * std::pow(_Tp{2} / _Tp(__n), _Tp{0.25})
		     / std::sqrt(_S_sqrt_pi * std::sin(__phi))
		     * std::sin(_S_pi * _Tp{0.75} + (__n / _Tp{2} + _Tp{0.25})
				   * (std::sin(_Tp{2} * __phi) - _Tp{2} * __phi))
		     * std::exp(__z * __z / _Tp{2});
	}
      else if (__z > __z_turn - __delta)
	{
	  auto __phi = std::acosh(__z / __z_turn);
	  return __f * __sign
		     * (__n_odd ? _Tp{1} : _Tp{1} / _S_sqrt_2)
		     * std::pow(_Tp{2}, _Tp(__n) / _Tp{2})
		     * std::pow(_Tp(__n), _Tp{-0.25})
		     / std::sqrt(_S_sqrt_2 * _S_sqrt_pi * std::sinh(__phi))
		     * std::exp((_Tp(__n) / _Tp{2} + _Tp{0.25})
			       * (_Tp{2} * __phi - std::sinh(_Tp{2} * __phi)))
		     * std::exp(__z * __z / _Tp{2});
	}
      else
	{
	  auto __arg = (__z - __z_turn)
		     * _S_sqrt_2 * std::pow(_Tp(__n), _Tp{1} / _Tp{6});
	  auto __airy = std::__detail::__airy(__arg);
	  return __f * __sign
		     * (__n_odd ? _S_sqrt_2 : _Tp{1})
		     * std::sqrt(_S_sqrt_pi * _S_sqrt_2)
		     * std::pow(_Tp{2}, _Tp(__n) / _Tp{2})
		     * std::pow(_Tp(__n), _Tp{-1} / _Tp{12})
		     * __airy.__Ai_value * std::exp(__z * __z / _Tp{2});
	}
    }

  struct __sqrt_factorial_t
  {
    unsigned int __n
    long double __value;
  };

  constexpr std::size_t _S_sqrt_fact_len = 101;
  constexpr __sqrt_factorial_t
  _S_sqrt_fact[__sqrt_fact_len]
  {
    {  0,  0.0L},
    {  1,  1.0L},
    {  2,  1.4142135623730950488016887242097L},
    {  3,  1.73205080756887729352744634150587L},
    {  4,  2.0L},
    {  5,  2.23606797749978969640917366873128L},
    {  6,  2.44948974278317809819728407470589L},
    {  7,  2.64575131106459059050161575363926L},
    {  8,  2.8284271247461900976033774484194L},
    {  9,  3.0L},
    { 10,  3.16227766016837933199889354443272L},
    { 11,  3.31662479035539984911493273667069L},
    { 12,  3.46410161513775458705489268301174L},
    { 13,  3.6055512754639892931192212674705L},
    { 14,  3.74165738677394138558374873231655L},
    { 15,  3.8729833462074168851792653997824L},
    { 16,  4.0L},
    { 17,  4.12310562561766054982140985597408L},
    { 18,  4.24264068711928514640506617262909L},
    { 19,  4.35889894354067355223698198385962L},
    { 20,  4.47213595499957939281834733746255L},
    { 21,  4.58257569495584000658804719372801L},
    { 22,  4.69041575982342955456563011354447L},
    { 23,  4.79583152331271954159743806416269L},
    { 24,  4.89897948556635619639456814941178L},
    { 25,  5.0L},
    { 26,  5.09901951359278483002822410902278L},
    { 27,  5.19615242270663188058233902451762L},
    { 28,  5.29150262212918118100323150727852L},
    { 29,  5.38516480713450403125071049154033L},
    { 30,  5.47722557505166113456969782800802L},
    { 31,  5.56776436283002192211947129891855L},
    { 32,  5.65685424949238019520675489683879L},
    { 33,  5.74456264653802865985061146821893L},
    { 34,  5.83095189484530047087415287754558L},
    { 35,  5.91607978309961604256732829156162L},
    { 36,  6.0L},
    { 37,  6.08276253029821968899968424520207L},
    { 38,  6.16441400296897645025019238145424L},
    { 39,  6.24499799839839820584689312093979L},
    { 40,  6.32455532033675866399778708886544L},
    { 41,  6.40312423743284868648821767462181L},
    { 42,  6.480740698407860230965967436088L},
    { 43,  6.557438524302000652344109997636L},
    { 44,  6.63324958071079969822986547334137L},
    { 45,  6.70820393249936908922752100619383L},
    { 46,  6.78232998312526813906455632662597L},
    { 47,  6.85565460040104412493587144908485L},
    { 48,  6.92820323027550917410978536602349L},
    { 49,  7.0L},
    { 50,  7.07106781186547524400844362104849L},
    { 51,  7.14142842854284999799939981136727L},
    { 52,  7.21110255092797858623844253494099L},
    { 53,  7.28010988928051827109730249152703L},
    { 54,  7.34846922834953429459185222411767L},
    { 55,  7.41619848709566294871139744080071L},
    { 56,  7.4833147735478827711674974646331L},
    { 57,  7.54983443527074969723668480694612L},
    { 58,  7.61577310586390828566141102715832L},
    { 59,  7.68114574786860817576968702173137L},
    { 60,  7.7459666924148337703585307995648L},
    { 61,  7.8102496759066543941297227357591L},
    { 62,  7.87400787401181101968503444881201L},
    { 63,  7.93725393319377177150484726091778L},
    { 64,  8.0L},
    { 65,  8.06225774829854965236661323030377L},
    { 66,  8.12403840463596036045988356826604L},
    { 67,  8.18535277187244996995370372473393L},
    { 68,  8.24621125123532109964281971194815L},
    { 69,  8.30662386291807485258426274490749L},
    { 70,  8.36660026534075547978172025785187L},
    { 71,  8.42614977317635863063413990620274L},
    { 72,  8.48528137423857029281013234525819L},
    { 73,  8.54400374531753116787164832623971L},
    { 74,  8.60232526704262677172947353504971L},
    { 75,  8.66025403784438646763723170752936L},
    { 76,  8.71779788708134710447396396771923L},
    { 77,  8.77496438739212206040638830741631L},
    { 78,  8.83176086632784685476404272695925L},
    { 79,  8.88819441731558885009144167540873L},
    { 80,  8.94427190999915878563669467492511L},
    { 81,  9.0L},
    { 82,  9.05538513813741662657380816698407L},
    { 83,  9.11043357914429888194562610468867L},
    { 84,  9.16515138991168001317609438745602L},
    { 85,  9.21954445729288731000227428176279L},
    { 86,  9.27361849549570375251641607399017L},
    { 87,  9.32737905308881504555447554232056L},
    { 88,  9.38083151964685910913126022708893L},
    { 89,  9.43398113205660381132066037762264L},
    { 90,  9.48683298050513799599668063329816L},
    { 91,  9.53939201416945649152621586023227L},
    { 92,  9.59166304662543908319487612832539L},
    { 93,  9.64365076099295499576003104743266L},
    { 94,  9.69535971483265802814888115084531L},
    { 95,  9.7467943448089639068384131998996L},
    { 96,  9.79795897113271239278913629882357L},
    { 97,  9.84885780179610472174621141491763L},
    { 98,  9.89949493661166534161182106946789L},
    { 99,  9.94987437106619954734479821001206L},
    {100, 10.0L},
  };

  template<typename _Tp>
    _Tp
    __sqrt_factorial(unsigned int __n)
    {
      if (__n < _S_sqrt_fact_len)
	return _S_sqrt_fact[__n];
      else
	return 0.0;
    };

  template<typename _Tp>
    _Tp
    __poly_hermite_norm_factor(unsigned int __n, _Tp __x)
    {
      const auto _S_root_pi = __gnu_cxx::__const_root_pi(__x);
      const auto _S_root4_pi = std::sqrt(_S_root_pi);
      const auto _S_root_2 = __gnu_cxx::__const_root_2(__x);
      return _S_root4_pi
	   * std::pow(_S_root_2, _Tp(__n)) * __sqrt_factorial(__n);
    }

  /**
   *  @brief  Compute the normalized Hermite polynomial by recursion.
   *  @todo  Tabulate sqrt(int) or even sqrt(i)/sqrt(i+1) and add a helper.
   */
  template<typename _Tp>
    _Tp
    __poly_hermite_norm_recursion(unsigned int __n, _Tp __x)
    {
      const auto _S_inv_root4_pi
	= _Tp{7.511255444649424828587030047762276930510e-1L};
      const auto _S_root_2 = __gnu_cxx::__const_root_2(__x);

      auto __Hn_nm2 = _S_inv_root4_pi;
      if (__n == 0)
	return __Hn_nm2;
      auto __Hn_nm1 = _S_root_2 * __x * __Hn_nm2;
      if (__n == 1)
	return __Hn_nm1;

      _Tp __Hn_n;
      for (unsigned int __i = 2; __i < __n; ++__i)
	{
	  __Hn_n = __x * std::sqrt(_Tp{2} / _Tp(__i)) * __Hn_nm1
		 - std::sqrt(_Tp(__i - 1) / _Tp(__i)) * __Hn_nm2;
	  __Hn_nm2 = __Hn_nm1;
	  __Hn_nm1 = __Hn_n;
	}

      return __Hn_n;
    }

  /**
   * @brief This routine returns the Probabilists Hermite polynomial
   * 	    of order n: @f$ He_n(x) @f$ by recursion on n.
   *
   * The Hermite polynomial is defined by:
   * @f[
   *   He_n(x) = (-1)^n e^{x^2/2} \frac{d^n}{dx^n} e^{-x^2/2}
   * @f]
   * or
   * @f[
   *   He_n(x) = \frac{1}{2^{-n/2}}H_n\left(\frac{x}{\sqrt{2}\right)
   * @f]
   * where @f$ H_n(x) @f$ is the Physicists Hermite function.
   *
   * @param __n The order of the Hermite polynomial.
   * @param __x The argument of the Hermite polynomial.
   * @return The value of the Hermite polynomial of order n
   * 	     and argument x.
   */
  template<typename _Tp>
    _Tp
    __poly_prob_hermite_recursion(unsigned int __n, _Tp __x)
    {
      // Compute He_0.
      auto __He_nm2 = _Tp{1};
      if (__n == 0)
	return __He_nm2;

      // Compute He_1.
      auto __He_nm1 = __x;
      if (__n == 1)
	return __He_nm1;

      // Compute He_n.
      _Tp __He_n;
      for (unsigned int __i = 2; __i <= __n; ++__i)
	{
	  __He_n = __x * __He_nm1 - _Tp(__i - 1) * __He_nm2;
	  __He_nm2 = __He_nm1;
	  __He_nm1 = __He_n;
	}

      return __He_n;
    }

  /**
   *  @brief  Compute the normalized probabalists Hermite polynomial by recursion.
   */
  template<typename _Tp>
    _Tp
    __poly_prob_hermite_norm_recursion(unsigned int __n, _Tp __x)
    {
      const auto _S_inv_root4_pi
	= _Tp{7.511255444649424828587030047762276930510e-1L};
      const auto _S_root_2 = __gnu_cxx::__const_root_2(__x);

      // Compute Hen_0.
      auto __Hen_nm2 = _S_inv_root4_pi;
      if (__n == 0)
	return __Hen_nm2;

      // Compute Hen_1.
      auto __Hen_nm1 = __x;
      if (__n == 1)
	return __Hen_nm1;

      // Compute Hen_n.
      _Tp __Hen_n;
      for (unsigned int __i = 2; __i <= __n; ++__i)
	{
	  __Hen_n = (__x * __Hen_nm1 - _Tp(__i - 1) * __Hen_nm2)
		  / std::sqrt(_Tp(__i));
	  __Hen_nm2 = __Hen_nm1;
	  __Hen_nm1 = __Hen_n;
	}

      return __Hen_n;
    }

  template<typename _Tp>
    _Tp
    __poly_prob_hermite_norm_factor(unsigned int __n, _Tp __x)
    {
      const auto _S_root_pi = __gnu_cxx::__const_root_pi(__x);
      const auto _S_root_2 = __gnu_cxx::__const_root_2(__x);
      const auto _S_root4_2pi = std::sqrt(_S_root_2 * _S_root_pi);
      return _S_root4_2pi * __sqrt_factorial(__n);
    }

#endif // _GLIBCXX_BITS_NEW_HERMITE_TCC
