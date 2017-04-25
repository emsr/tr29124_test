// Special functions -*- C++ -*-

// Copyright (C) 2006-2017 Free Software Foundation, Inc.
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

/** @file bits/sf_zeta.tcc
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

//
// ISO C++ 14882 TR29124: Mathematical Special Functions
//

// Written by Edward Smith-Rowland and Florian Goth
//
// References:
// (1) Handbook of Mathematical Functions,
//     Ed. by Milton Abramowitz and Irene A. Stegun,
//     Dover Publications, New-York, Section 5, pp. 807-808.
// (2) The Gnu Scientific Library, http://www.gnu.org/software/gsl
// (3) Gamma, Exploring Euler's Constant, Julian Havil,
//     Princeton, 2003.
// (4) David C. Wood, "The Computation of Polylogarithms."

#ifndef _GLIBCXX_BITS_SF_ZETA_TCC
#define _GLIBCXX_BITS_SF_ZETA_TCC 1

#pragma GCC system_header

#include <ext/math_const.h>

namespace std _GLIBCXX_VISIBILITY(default)
{
// Implementation-space details.
namespace __detail
{
_GLIBCXX_BEGIN_NAMESPACE_VERSION


  /**
   * Coefficients for Euler-Maclaurin summation of zeta functions.
   * @f[
   *    B_{2j} / (2j)!
   * @f]
   * where @f$ B_k @f$ are the Bernoulli numbers.
   */
  constexpr size_t _Num_Euler_Maclaurin_zeta = 100;
  constexpr long double
  _S_Euler_Maclaurin_zeta[_Num_Euler_Maclaurin_zeta]
  {
    1.00000000000000000000000000000000000L,
    8.33333333333333333333333333333333293e-02L,
   -1.38888888888888888888888888888888875e-03L,
    3.30687830687830687830687830687830609e-05L,
   -8.26719576719576719576719576719576597e-07L,
    2.08767569878680989792100903212014296e-08L,
   -5.28419013868749318484768220217955604e-10L,
    1.33825365306846788328269809751291227e-11L,
   -3.38968029632258286683019539124944218e-13L,
    8.58606205627784456413590545042562615e-15L,
   -2.17486869855806187304151642386591768e-16L,
    5.50900282836022951520265260890225438e-18L,
   -1.39544646858125233407076862640635480e-19L,
    3.53470703962946747169322997780379902e-21L,
   -8.95351742703754684639940801672890898e-23L,
    2.26795245233768305922449726817928506e-24L,
   -5.74479066887220244232839527972348697e-26L,
    1.45517247561486490107622443104134417e-27L,
   -3.68599494066531017606286927671534186e-29L,
    9.33673425709504466636710365024250844e-31L,
   -2.36502241570062993304902926977940878e-32L,
    5.99067176248213430064218240871649208e-34L,
   -1.51745488446829026064464819837699250e-35L,
    3.84375812545418822940606216740290214e-37L,
   -9.73635307264669102780496423687655647e-39L,
    2.46624704420068095513732412234574675e-40L,
   -6.24707674182074368796151949260113716e-42L,
    1.58240302446449142838660289637807111e-43L,
   -4.00827368594893596494573716493578672e-45L,
    1.01530758555695563022273865751378251e-46L,
   -2.57180415824187174746079460115444631e-48L,
    6.51445603523381492510893884688687796e-50L,
   -1.65013099068965245381972311645983560e-51L,
    4.17983062853947589044505904302589394e-53L,
   -1.05876346677029087587739692698831756e-54L,
    2.68187919126077066314325024787533130e-56L,
   -6.79327935110742120171687695308170512e-58L,
    1.72075776166814048850302218744398370e-59L,
   -4.35873032934889383811051833724115685e-61L,
    1.10407929036846667370868730253333297e-62L,
   -2.79666551337813450363217687544853825e-64L,
    7.08403650167947018923360554701085951e-66L,
   -1.79440740828922406419836715043711635e-67L,
    4.54528706361109610084319503460215356e-69L,
   -1.15133466319820517965514570874684656e-70L,
    2.91636477109236135051215065347539762e-72L,
   -7.38723826349733755172097357215101415e-74L,
    1.87120931176379530341655669167307802e-75L,
   -4.73982855776179939823365705874837915e-77L,
    1.20061259933545065010289119344900538e-78L,
   -3.04118724151429237818089382050570806e-80L,
    7.70341727470510626032951412805395999e-82L,
   -1.95129839090988306787181681770634078e-83L,
    4.94269656515946146653024823540482418e-85L,
   -1.25199966591718479000903037235065000e-86L,
    3.17135220176351545507490909160914066e-88L,
   -8.03312897073533444702820185587950603e-90L,
    2.03481533916614656707738578184657923e-91L,
   -5.15424746644747384952223952829841139e-93L,
    1.30558613521494672211652811590429162e-94L,
   -3.30708831417509124211476473245870569e-96L,
    8.37695256004909128671576487001515732e-98L,
   -2.12190687174971376532740302523869392e-99L,
    5.37485289561228024546639030950850747e-101L,
   -1.36146614321720693646766878619458764e-102L,
    3.44863402799339902711847245019887448e-104L,
   -8.73549204163835504185857147126132564e-106L,
    2.21272598339254970369646016351296266e-107L,
   -5.60490039283722415865004549966830368e-109L,
    1.41973785499917876418113219390120402e-110L,
   -3.59623799825876265563506189711913868e-112L,
    9.10937726607823184392152343960921788e-114L,
   -2.30743221710912328319788632553096580e-115L,
    5.84479408529900198086678715736051924e-117L,
   -1.48050363717057449175295580276805879e-118L,
    3.75015952262271968009868219553161640e-120L,
   -9.49926504199295827433414434212738749e-122L,
    2.40619194446751986855735721071770656e-123L,
   -6.09495539710268473564828955326838575e-125L,
    1.54387023770424714174247152273407761e-126L,
   -3.91066899685929231015496868644244724e-128L,
    9.90584028987942974230798861783873238e-130L,
   -2.50917865785635531226030684021698888e-131L,
    6.35582378960245978851062983907498289e-133L,
   -1.60994897346162503541124270943910111e-134L,
    4.07804838987246223591940353095558419e-136L,
   -1.03298172453151901860990367375763718e-137L,
    2.61657327526092016567608044513104821e-139L,
   -6.62785753340862278242754283133858879e-141L,
    1.67885592574436730038021560419761003e-142L,
   -4.25259173903430060933988427590406140e-144L,
    1.07719407136645728953530432238003752e-145L,
   -2.72856445808395795763235897809610861e-147L,
    6.91153451343701367068428806267885056e-149L,
   -1.75071214421576824424282262560183841e-150L,
    4.43460566672391957275357495883452114e-152L,
   -1.12329873784871496672434722468912147e-153L,
    2.84534894256939708475369589042384123e-155L,
   -7.20735306841513677412535382503699916e-157L,
    1.82564385955014175253212078464905862e-158L
  };

  /**
   *  @brief Compute the dilogarithm function @f$ Li_2(x) @f$
   *         by summation for x <= 1.
   *
   *  The dilogarithm function is defined by:
   *   @f[
   *     Li_2(x) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} \mbox{ for } s > 1
   *   @f]
   *  For |x| near 1 use the reflection formulae:
   *   @f[
   *     Li_2(-x) + Li_2(1-x) = \frac{\pi^2}{6} - \ln(x) \ln(1-x)
   *   @f]
   *   @f[
   *     Li_2(-x) - Li_2(1-x) - \frac{1}{2}Li_2(1-x^2)
   *         = -\frac{\pi^2}{12} - \ln(x) \ln(1-x)
   *   @f]
   *  For x < -1 use the reflection formula:
   *   @f[
   *     Li_2(1-x) - Li_2(1-\frac{1}{1-x}) - \frac{1}{2}(\ln(x))^2
   *   @f]
   */
  template<typename _Tp>
    _Tp
    __dilog(_Tp __x)
    {
      constexpr unsigned long long _S_maxit = 100000ULL;
      const auto _S_eps = 10 * __gnu_cxx::__epsilon(std::real(__x));
      const auto _S_pipio6 = __gnu_cxx::__const_pi_sqr_div_6(std::real(__x));
      if (__isnan(__x))
	return __gnu_cxx::__quiet_NaN(std::real(__x));
      else if (__x > _Tp{+1})
	std::__throw_range_error(__N("dilog: argument greater than one"));
      else if (__x < _Tp{-1})
	{
	  auto __lnfact = std::log(_Tp{1} - __x);
	  return -__dilog(_Tp{1} - _Tp{1} / (_Tp{1} - __x))
		 - _Tp{0.5L} * __lnfact * __lnfact;
	}
      else if (__x == _Tp{0})
	return _Tp{0};
      else if (__x == _Tp{1})
	return _S_pipio6;
      else if (__x == -_Tp{1})
	return -_Tp{0.5L} * _S_pipio6;
      else if (__x > _Tp{0.5L})
	return _S_pipio6 - std::log(__x) * std::log(_Tp{1} - __x)
	     - __dilog(_Tp{1} - __x);
      else if (__x < -_Tp{0.5L})
	return -_Tp{0.5L} * _S_pipio6 - std::log(_Tp{1} - __x) * std::log(-__x)
	     + __dilog(_Tp{1} + __x) -_Tp{0.5L} * __dilog(_Tp{1} - __x * __x);
      else
	{
	  _Tp __sum = 0;
	  _Tp __fact = 1;
	  for (auto __i = 1ULL; __i < _S_maxit; ++__i)
	    {
	      __fact *= __x;
	      auto __term = __fact / (__i * __i);
	      __sum += __term;
	      if (std::abs(__term) < _S_eps * std::abs(__sum))
		break;
	      if (__i + 1 == _S_maxit)
		std::__throw_runtime_error(__N("__dilog: sum failed"));
	    }
	  return __sum;
	}
    }


  /**
   * @brief  Compute the Riemann zeta function @f$ \zeta(s) @f$
   * 	     by summation for s > 1.
   *
   * The Riemann zeta function is defined by:
   *  @f[
   * 	\zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
   *  @f]
   * For s < 1 use the reflection formula:
   *  @f[
   * 	\zeta(s) = (2\pi)^s \Gamma(1-s) \zeta(1-s) / \pi
   *  @f]
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta_sum(_Tp __s)
    {
      using _Val = _Tp;
      using _Real = __num_traits_t<_Val>;
      // A user shouldn't get to this.
      if (std::real(__s) < _Real{1})
	std::__throw_domain_error(__N("__riemann_zeta_sum: "
				      "Bad argument in zeta sum."));
      else if (std::real(__s) > _Real{1})
	{
	  constexpr unsigned int _S_max_iter = 10000;
	  const auto _S_eps = __gnu_cxx::__epsilon(std::real(__s));
	  auto __zeta = _Val{1};
	  for (unsigned int __k = 2; __k < _S_max_iter; ++__k)
	    {
	      auto __term = std::pow(_Val(__k), -__s);
	      __zeta += __term;
	      if (std::abs(__term) < _S_eps * std::abs(__zeta)
		  || std::abs(__term) < _S_eps
		     && std::abs(__zeta) < _Real{100} * _S_eps)
		break;
	    }
	  return __zeta;
	}
      else
	{
	  const auto _S_pi = __gnu_cxx::__const_pi(std::real(__s));
	  auto __zeta = std::pow(_Real{2} * _S_pi, __s)
		      * __sin_pi(_Real{0.5L} * __s) * __gamma(_Val{1} - __s)
		      * __riemann_zeta_sum(_Val{1} - __s) / _S_pi;
	  return __zeta;
	}
    }


  /**
   * @brief  Evaluate the Riemann zeta function @f$ \zeta(s) @f$
   * 	     by an alternate series for s > 0.
   *
   * This is a specialization of the code for the Hurwitz zeta function.
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta_euler_maclaurin(_Tp __s)
    {
      using _Val = _Tp;
      using _Real = __num_traits_t<_Val>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__s));
      const auto _S_N = 10 + __gnu_cxx::__digits10(std::real(__s)) / _Tp{2};
      const auto _S_jmax = _Num_Euler_Maclaurin_zeta - 1;

      const auto __pmax  = std::pow(_Val{_S_N + 1}, -__s);
      const auto __denom = _Val{_S_N + 1} * _Val{_S_N + 1};
      auto __ans = __pmax * (_Val{_S_N + 1} / (__s - _Val{1}) + _Val{0.5L});
      for (auto __k = 0; __k < _S_N; ++__k)
        __ans += std::pow(_Val(__k + 1), -__s);

      auto __fact = __pmax * __s / _Val{_S_N + 1};
      auto __delta_prev = __gnu_cxx::__max(__s);
      for (auto __j = 0; __j < _S_jmax; ++__j)
        {
	  auto __delta = _S_Euler_Maclaurin_zeta[__j + 1] * __fact;
	  if (std::abs(__delta) > __delta_prev)
	    break;
	  __delta_prev = std::abs(__delta);
	  __ans += __delta;
	  if (std::abs(__delta) < _Real{0.5L} * _S_eps * std::abs(__ans)
	      || std::abs(__delta) < _S_eps
		 && std::abs(__ans) < _Real{100} * _S_eps)
	    break;
	  __fact *= (__s + _Val(2 * __j + 1)) * (__s + _Val(2 * __j + 2))
		  / __denom;
        }

      return __ans;
    }


  /**
   * @brief  Evaluate the Riemann zeta function by series for all s != 1.
   * 	     Convergence is great until largish negative numbers.
   * 	     Then the convergence of the > 0 sum gets better.
   *
   * The series is:
   * @f[
   * 	\zeta(s) = \frac{1}{1-2^{1-s}}
   * 		   \sum_{n=0}^{\infty} \frac{1}{2^{n+1}}
   * 		   \sum_{k=0}^{n} (-1)^k \frac{n!}{(n-k)!k!} (k+1)^{-s}
   * @f]
   * Havil 2003, p. 206.
   *
   * The Riemann zeta function is defined by:
   * @f[
   * 	\zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for s > 1
   * @f]
   * For s < 1 use the reflection formula:
   * @f[
   * 	\zeta(s) = (2\pi)^s \Gamma(1-s) \zeta(1-s) / \pi
   * @f]
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta_m_1_glob(_Tp __s)
    {
      using _Val = _Tp;
      using _Real = __num_traits_t<_Val>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__s));
      //  Max e exponent before overflow.
      const auto __max_binom
		 = std::exp(__gnu_cxx::__max_exponent10(std::real(__s))
			    * std::log(_Real{10}) - _Real{1});

      auto __zeta_m_1 = _Val{0};
      // This for loop starts at 1 because we already calculated the
      // value of the zeroeth order in __zeta_m_1 above
      const unsigned int __maxit = 10000;
      auto __num = _Real{0.25L};
      for (unsigned int __i = 1; __i < __maxit; ++__i)
	{
	  bool __punt = false;
	  auto __binom = _Real{1};
	  auto __term = _Val{0};
	  // This for loop starts at 1 because we already calculated the value
	  // of the zeroeth order in __term above.
	  for (unsigned int __j = 1; __j <= __i; ++__j)
	    {
	      __binom *= -_Real(__i - __j + 1) / _Real(__j);
	      if(std::abs(__binom) > __max_binom )
		{
		  // This only gets hit for x << 0.
		  __punt = true;
		  break;
		}
	      __term += __binom * std::pow(_Val(1 + __j), -__s);
	    }
	  if (__punt)
	    break;
	  __term *= __num;
	  __zeta_m_1 += __term;
	  if (std::abs(__term) < _S_eps * std::abs(__zeta_m_1)
	      || std::abs(__term) < _S_eps
		 && std::abs(__zeta_m_1) < _Real{100} * _S_eps)
	    break;
	  __num *= _Real{0.5L};
	}
      const auto __pow2 = std::pow(_Val{2}, _Val{1} - __s);
      __zeta_m_1 += __pow2;
      __zeta_m_1 /= _Val{1} - __pow2;
      return __zeta_m_1;
    }

  template<typename _Tp>
    _Tp
    __riemann_zeta_glob(_Tp __s)
    {
      using _Val = _Tp;
      using _Real = __num_traits_t<_Val>;

      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__s));
      const auto _S_pi_2 = __gnu_cxx::__const_pi_half(std::real(__s));

      //  This series works until the binomial coefficient blows up
      //  so use reflection.
      if (std::real(__s) < _Real{0})
	{
	  if (__gnu_cxx::__fp_is_even_integer(__s))
	    return _Val{0};
	  else
	    {
	      auto __zeta = __riemann_zeta_glob(_Val{1} - __s);
	      __zeta *= std::pow(_Real{2} * _S_pi, __s)
		     * __sin_pi(_Real{0.5L} * __s)
		     * __gamma(_Val{1} - __s) / _S_pi;
	      return __zeta;
	    }
	}
      else
	return _Tp{1} + __riemann_zeta_m_1_glob(__s);
    }


  /**
   * @brief  Compute the Riemann zeta function @f$ \zeta(s) @f$
   * 	 using the product over prime factors.
   *
   * @f[
   * 	\zeta(s) = \Pi_{i=1}^\infty \frac{1}{1 - p_i^{-s}}
   * @f]
   * where @f$ {p_i} @f$ are the prime numbers.
   *
   * The Riemann zeta function is defined by:
   * @f[
   *    \renewcommand\Re{\operatorname{Re}}
   *    \renewcommand\Im{\operatorname{Im}}
   * 	\zeta(s) = \sum_{k=1}^{\infty} \frac{1}{k^{s}} for \Re{s} > 1
   * @f]
   * For \Re(s) < 1 use the reflection formula:
   * @f[
   * 	\zeta(s) = (2\pi)^s \Gamma(1-s) \zeta(1-s) / \pi
   * @f]
   *
   * @param __s The argument
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta_product(_Tp __s)
    {
      using _Val = _Tp;
      using _Real = __num_traits_t<_Val>;
      extern const unsigned long __prime_list[];

      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__s));
      constexpr unsigned long
        _S_num_primes = sizeof(unsigned long) != 8 ? 256 : 256 + 48;

      auto __zeta = _Val{1};
      for (unsigned long __i = 0; __i < _S_num_primes; ++__i)
	{
	  const auto __fact = _Val{1}
			    - std::pow(_Real(__prime_list[__i]), -__s);
	  __zeta *= __fact;
	  if (std::abs(_Tp{1} - __fact) < _S_eps) // Assume zeta near 1.
	    break;
	}

      __zeta = _Tp{1} / __zeta;

      return __zeta;
    }


  /**
   * Table of zeta(n) - 1 from 0 - 120.
   * MPFR @ 128 bits.
   */
  constexpr size_t
  _S_num_zetam1 = 121;

  constexpr long double
  _S_zetam1[_S_num_zetam1]
  {
    -1.5L,                                          //   0
    std::numeric_limits<long double>::infinity(),   //   1
    6.449340668482264364724151666460251892177e-1L,  //   2
    2.020569031595942853997381615114499907647e-1L,  //   3
    8.232323371113819151600369654116790277462e-2L,  //   4
    3.692775514336992633136548645703416805713e-2L,  //   5
    1.734306198444913971451792979092052790186e-2L,  //   6
    8.349277381922826839797549849796759599843e-3L,  //   7
    4.077356197944339378685238508652465258950e-3L,  //   8
    2.008392826082214417852769232412060485604e-3L,  //   9
    9.945751278180853371459589003190170060214e-4L,  //  10
    4.941886041194645587022825264699364686068e-4L,  //  11
    2.460865533080482986379980477396709604160e-4L,  //  12
    1.227133475784891467518365263573957142749e-4L,  //  13
    6.124813505870482925854510513533374748177e-5L,  //  14
    3.058823630702049355172851064506258762801e-5L,  //  15
    1.528225940865187173257148763672202323739e-5L,  //  16
    7.637197637899762273600293563029213088257e-6L,  //  17
    3.817293264999839856461644621939730454694e-6L,  //  18
    1.908212716553938925656957795101353258569e-6L,  //  19
    9.539620338727961131520386834493459437919e-7L,  //  20
    4.769329867878064631167196043730459664471e-7L,  //  21
    2.384505027277329900036481867529949350419e-7L,  //  22
    1.192199259653110730677887188823263872549e-7L,  //  23
    5.960818905125947961244020793580122750393e-8L,  //  24
    2.980350351465228018606370506936601184471e-8L,  //  25
    1.490155482836504123465850663069862886482e-8L,  //  26
    7.450711789835429491981004170604119454712e-9L,  //  27
    3.725334024788457054819204018402423232885e-9L,  //  28
    1.862659723513049006403909945416948061669e-9L,  //  29
    9.313274324196681828717647350212198135677e-10L, //  30
    4.656629065033784072989233251220071062704e-10L, //  31
    2.328311833676505492001455975940495024831e-10L, //  32
    1.164155017270051977592973835456309516528e-10L, //  33
    5.820772087902700889243685989106305417368e-11L, //  34
    2.910385044497099686929425227884046410669e-11L, //  35
    1.455192189104198423592963224531842098334e-11L, //  36
    7.275959835057481014520869012338059265263e-12L, //  37
    3.637979547378651190237236355873273513051e-12L, //  38
    1.818989650307065947584832100730085030987e-12L, //  39
    9.094947840263889282533118386949087534482e-13L, //  40
    4.547473783042154026799112029488570339961e-13L, //  41
    2.273736845824652515226821577978691217250e-13L, //  42
    1.136868407680227849349104838025906441861e-13L, //  43
    5.684341987627585609277182967524068526363e-14L, //  44
    2.842170976889301855455073704942662033022e-14L, //  45
    1.421085482803160676983430714173953721447e-14L, //  46
    7.105427395210852712877354479956799457540e-15L, //  47
    3.552713691337113673298469534059343240771e-15L, //  48
    1.776356843579120327473349014400279865980e-15L, //  49
    8.881784210930815903096091386391386649172e-16L, //  50
    4.440892103143813364197770940268122986877e-16L, //  51
    2.220446050798041983999320094204660286072e-16L, //  52
    1.110223025141066133720544569921388238976e-16L, //  53
    5.551115124845481243723736590509454214979e-17L, //  54
    2.775557562136124172581632453854098189256e-17L, //  55
    1.387778780972523276283909490650087159020e-17L, //  56
    6.938893904544153697446085326249613606421e-18L, //  57
    3.469446952165922624744271496109153952849e-18L, //  58
    1.734723476047576572048972969937766807693e-18L, //  59
    8.673617380119933728342055067345929347336e-19L, //  60
    4.336808690020650487497023565906367637200e-19L, //  61
    2.168404344997219785013910168326102593523e-19L, //  62
    1.084202172494241406301271116539546929420e-19L, //  63
    5.421010862456645410918700404413613405660e-20L, //  64
    2.710505431223468831954621311921825336782e-20L, //  65
    1.355252715610116458148523399711726681995e-20L, //  66
    6.776263578045189097995298741000309894844e-21L, //  67
    3.388131789020796818085703100408435571778e-21L, //  68
    1.694065894509799165406492747048108912984e-21L, //  69
    8.470329472546998348246992605151870123760e-22L, //  70
    4.235164736272833347862270482171914722722e-22L, //  71
    2.117582368136194731844209444015663667353e-22L, //  72
    1.058791184068023385226500150767838272960e-22L, //  73
    5.293955920339870323813912246795908957429e-23L, //  74
    2.646977960169852961134116619883673563755e-23L, //  75
    1.323488980084899080309451049270121075322e-23L, //  76
    6.617444900424404067355245869046478332807e-24L, //  77
    3.308722450212171588946956563227359069812e-24L, //  78
    1.654361225106075646229923736818740547512e-24L, //  79
    8.271806125530344403671108096295678003592e-25L, //  80
    4.135903062765160926009383852215164090474e-25L, //  81
    2.067951531382576704395963965944918517449e-25L, //  82
    1.033975765691287099328403715492352137455e-25L, //  83
    5.169878828456431320410159441971315309917e-26L, //  84
    2.584939414228214268127816150081035315909e-26L, //  85
    1.292469707114106670038085128629065184730e-26L, //  86
    6.462348535570531803437454412518556478869e-27L, //  87
    3.231174267785265386134631538638949625204e-27L, //  88
    1.615587133892632521206406698623221009248e-27L, //  89
    8.077935669463162033155494313137423014210e-28L, //  90
    4.038967834731580825616620004023205489476e-28L, //  91
    2.019483917365790349161443228820936759716e-28L, //  92
    1.009741958682895153362818597460967711233e-28L, //  93
    5.048709793414475696133364781858743133105e-29L, //  94
    2.524354896707237824529247247973459828381e-29L, //  95
    1.262177448353618904396004317475305931953e-29L, //  96
    6.310887241768094495295138721816048328696e-30L, //  97
    3.155443620884047239436836171504799139404e-30L, //  98
    1.577721810442023616297279256834389142642e-30L, //  99
    7.888609052210118067801840968499904004972e-31L, // 100
    3.944304526105059027058642826413931148366e-31L, // 101
    1.972152263052529513529321413206965574183e-31L, // 102
    9.860761315262647567646607066034827870915e-32L, // 103
    4.930380657631323783823303533017413935458e-32L, // 104
    2.465190328815661891911651766508706967729e-32L, // 105
    1.232595164407830945955825883254353483864e-32L, // 106
    6.162975822039154729779129416271767419322e-33L, // 107
    3.081487911019577364889564708135883709661e-33L, // 108
    1.540743955509788682444782354067941854830e-33L, // 109
    7.888609052210118067801840968499904004972e-31L, // 110
    3.851859888774471706111955885169854637076e-34L, // 111
    1.925929944387235853055977942584927318538e-34L, // 112
    9.629649721936179265279889712924636592691e-35L, // 113
    4.814824860968089632639944856462318296345e-35L, // 114
    2.407412430484044816319972428231159148173e-35L, // 115
    1.203706215242022408159986214115579574086e-35L, // 116
    6.018531076210112040799931070577897870432e-36L, // 117
    3.009265538105056020399965535288948935216e-36L, // 118
    1.504632769052528010199982767644474467608e-36L, // 119
    7.523163845262640050999913838222372338039e-37L, // 120
  };


  /**
   * @brief  Return the Riemann zeta function @f$ \zeta(s) - 1 @f$
   *
   *
   * @param __s The argument @f$ s != 1 @f$
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta_m_1(_Tp __s)
    {
      using _Real = __num_traits_t<_Tp>;
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__s));
      if (__s == _Real{1})
	return __gnu_cxx::__infinity(std::real(__s));

      auto __n = __gnu_cxx::__fp_is_integer(__s);
      if (__n && __n() >= 0 && __n() < _S_num_zetam1)
	return _Tp(_S_zetam1[__n()]);
      else if (std::real(__s) > _Real{0})
	return __riemann_zeta_m_1_glob(__s);
      else // Re[s] < 0
	return std::pow(_Real{2} * _S_pi, __s)
	     * __sin_pi(_Real{0.5L} * __s)
	     * __gamma(_Tp{1} - __s)
	     * (_Real{1} + __riemann_zeta_m_1(_Tp{1} - __s)) / _S_pi
	     - _Real{1};
    }


  /**
   * @brief  Return the Riemann zeta function @f$ \zeta(s) @f$.
   *
   * The Riemann zeta function is defined by:
   * @f[
   * 	\zeta(s) = \sum_{k=1}^{\infty} k^{-s} \mbox{ for } \Re(s) > 1 \\
   * 		   \frac{(2\pi)^s}{\pi} \sin(\frac{\pi s}{2})
   * 		   \Gamma(1 - s) \zeta(1 - s) \mbox{ for } \Re(s) < 1
   * @f]
   *
   * @param __s The argument
   */
  template<typename _Tp>
    _Tp
    __riemann_zeta(_Tp __s)
    {
      using _Val = _Tp;
      using _Real = __num_traits_t<_Val>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::real(__s));
      const auto _S_inf = __gnu_cxx::__infinity(std::real(__s));
      const auto _S_pi = __gnu_cxx::__const_pi(std::real(__s));
      if (__isnan(__s))
	return _S_NaN;
      else if (__s == _Val{1})
	return _S_inf;
      else
	{
	  const auto __S_max = __gnu_cxx::__digits(std::real(__s));
	  const auto __p = __gnu_cxx::__fp_is_integer(__s);
	  if (__p && __p() >= 0)
	    return _Real{1} + __riemann_zeta_m_1(_Real(__p()));
	  else if (__p && __p() < 0)
	    {
	      const auto __r = _Real(__p());
	      _Tp __zeta = std::pow(_Real{2} * _S_pi, __r)
			 * __sin_pi(_Real{0.5L} * __r)
			 * __gamma(_Real{1} - __r)
			 * (_Real{1} + __riemann_zeta_m_1(_Real{1} - __r))
			 / _S_pi;
	      return __zeta;
	    }
	  else if (std::real(__s) < -_Real{19})
	    {
	      auto __zeta = __riemann_zeta_product(_Val{1} - __s);
	      __zeta *= std::pow(_Real{2} * _S_pi, __s)
		     * __sin_pi(_Real{0.5L} * __s)
		     * std::exp(__log_gamma(_Val{1} - __s))
		     / _S_pi;
	      return __zeta;
	    }
	  else if (std::real(__s) < __S_max)
	    {
	      /// @todo Global double sum or MacLaurin series in riemann_zeta?
	      bool __glob = true;
	      if (__glob)
		return __riemann_zeta_glob(__s);
	      else
		return __riemann_zeta_sum(__s);
	    }
	  else
	    return _Val{1} + std::pow(_Val{2}, -__s);
	}
    }

  /**
   * @brief  Return the Hurwitz zeta function @f$ \zeta(s,a) @f$
   * 	     for all s != 1 and a > -1.
   * @see An efficient algorithm for accelerating the convergence
   * 	  of oscillatory series, useful for computing the
   * 	  polylogarithm and Hurwitz zeta functions, Linas Vep"0160tas
   *
   * @param __s The argument @f$ s != 1 @f$
   * @param __a The scale parameter @f$ a > -1 @f$
   */
  template<typename _Tp>
    _Tp
    __hurwitz_zeta_euler_maclaurin(_Tp __s, _Tp __a)
    {
      using _Val = _Tp;
      using _Real = __num_traits_t<_Val>;
      const auto _S_eps = __gnu_cxx::__epsilon(std::real(__s));
      const int _S_N = 10 + __gnu_cxx::__digits10(std::real(__s)) / 2;
      const int _S_jmax = _Num_Euler_Maclaurin_zeta - 1;

      const auto __pmax  = std::pow(_Val{_S_N} + __a, -__s);
      auto __ans = __pmax
		 * ((_Val{_S_N} + __a) / (__s - _Val{1}) + _Real{0.5L});
      for(auto __k = 0; __k < _S_N; ++__k)
	__ans += std::pow(_Tp(__k) + __a, -__s);

      auto __sfact = __s;
      auto __pfact = __pmax / (_Val(_S_N) + __a);
      for(auto __j = 0; __j < _S_jmax; ++__j)
	{
	  auto __delta = _Real(_S_Euler_Maclaurin_zeta[__j + 1])
			* __sfact * __pfact;
	  __ans += __delta;
	  if (std::abs(__delta) < _Real{0.5L} * _S_eps * std::abs(__ans))
	    break;
	  __sfact *= (__s + _Val(2 * __j + 1)) * (__s + _Val(2 * __j + 2));
	  __pfact /= (_Val{_S_N} + __a) * (_Val{_S_N} + __a);
	}

      return __ans;
    }

  /**
   * @brief  Return the Hurwitz zeta function @f$ \zeta(s,a) @f$
   * for all s != 1 and a > -1.
   *
   * The Hurwitz zeta function is defined by:
   * @f[
   *   \zeta(s,a) = \sum_{n=0}^{\infty} \frac{1}{(n+a)^s}
   * @f]
   * The Riemann zeta function is a special case:
   * @f[
   *   \zeta(s) = \zeta(s,1)
   * @f]
   *
   * @param __s The argument @f$ s != 1 @f$
   * @param __a The scale parameter @f$ a > -1 @f$
   */
  template<typename _Tp>
    _Tp
    __hurwitz_zeta(_Tp __s, _Tp __a)
    {
      using _Val = _Tp;
      using _Real = __num_traits_t<_Val>;
      const auto _S_NaN = __gnu_cxx::__quiet_NaN(std::real(__s));
      const auto _S_inf = __gnu_cxx::__infinity(std::real(__s));
      if (__isnan(__s) || __isnan(__a))
	return _S_NaN;
      else if (__a == _Real{1})
	{
	  if (__s == _Real{1})
	    return _S_inf;
	  else
	    return __riemann_zeta(__s);
	}
      else
        return __hurwitz_zeta_euler_maclaurin(__s, __a);
    }

  /**
   * Return the Debye function.
   * The Debye functions are related to the incomplete Riemann zeta function:
   * @f[
   *    \zeta_x(s) = \frac{1}{\Gamma(s)}\int_{0}^{x}\frac{t^{s-1}}{e^t-1}dt
   *          = \sum_{k=1}^{\infty}\frac{P(s,kx)}{k^s}
   * @f]
   * @f[
   *    Z_x(s) = \frac{1}{\Gamma(s)}\int_{x}^{\infty}\frac{t^{s-1}}{e^t-1}dt
   *          = \sum_{k=1}^{\infty}\frac{Q(s,kx)}{k^s}
   * @f]
   * where @f$ P(a,x), Q(a,x) @f$ is the incomplete gamma function ratios.
   * The Debye functions are:
   * @f[
   *    D_n(x) = \frac{n}{x^n}\int_{0}^{x}\frac{t^n}{e^t-1}dt
   *           = \Gamma(n+1)\zeta_x(n+1)
   * @f]
   * and
   * @f[
   *    \int_{0}^{x}\frac{t^n}{e^t-1}dt = \Gamma(n+1)\zeta_x(n+1)
   * @f]
   *
   * @todo: We should return both the Debye function and it's complement.
   */
  template<typename _Tp>
    _Tp
    __debye(unsigned int __n, _Tp __x)
    {
      if (std::isnan(__x))
	return std::numeric_limits<_Tp>::quiet_NaN();
      else if (__n < 1)
	std::__throw_domain_error("__debye: Degree n must be positive.");
      else if (__x >= _Tp{3})
	{
	  // For values up to 4.80 the list of zeta functions
	  // and the sum up to k < K are huge enough to gain
	  // numeric stability in the sum

	  // n!zeta(n) is the integral for x=inf, Abramowitz & Stegun 27.1.3
	  auto __sum = _Tp{0};
	  if (std::__detail::_S_num_factorials<_Tp>)
	    __sum += std::__detail::__factorial<_Tp>(__n)
		   * std::__detail::__riemann_zeta<_Tp>(__n + 1);
	  else
	    return __gnu_cxx::__infinity(__x);

	  /**
	   * Compute the Debye function:
	   * @f[
	   *    D_n(x) = 1 - \sum_{k = 1}^{\infty} e^{-kx}
	   *       \frac{n}{k}\sum_{m=0}^{n}\frac{n!}{(n-m)!}frac{1}{(kx)^m}
	   * @f]
	   * Abramowitz & Stegun 27.1.2
	   */
	  const std::size_t _S_max_iter = 100;
	  auto __term = _Tp{0};
	  for(unsigned int __k = 1; __k < _S_max_iter; ++__k)
	    {
	      const auto __xk = __x * __k;
	      auto __ksum = _Tp{1} / __xk;
	      auto __kterm = _Tp(__n) * __ksum / __xk;  // n / (xk)^2
	      for (unsigned int __s = 1; __s <= __n; ++__s)
		__ksum += std::exchange(__kterm,
					_Tp(__n - __s) * __kterm / __xk);

	      __term -= std::exp(-__xk) * __ksum * std::pow(__x, _Tp(__n + 1));
	    }
	  __sum += __term;
	  return _Tp(__n) * __sum / std::pow(__x, _Tp(__n));
	}
      else
	{
	  /**
	   * Compute the Debye function:
	   * @f[
	   *    D_n(x) = 1 - \frac{n x}{2(n+1)}
	   *       + n \sum_{k = 1}^{\infty} \frac{B_{2k} x^{2k}}{(2k + n)(2k)!}
	   * @f]
           * for @f$ |x| < 2\pi @f$.
	   * Abramowitz-Stegun 27.1.1
	   */
	  const auto _S_eps = __gnu_cxx::__epsilon(__x);
	  const std::size_t _S_max_iter = 200;
	  const auto _S_1_2pi = __gnu_cxx::__const_one_div_2_pi(__x);
	  const auto __x2pi = __x * _S_1_2pi;
	  const auto __x2pi2 = __x2pi * __x2pi;
	  auto __x2pi2k = __x2pi2;
	  auto __sum = _Tp{0};
	  for(unsigned int __k = 1; __k < _S_max_iter; ++__k)
	    {
	      const auto __term = _Tp{2}
				* std::__detail::__riemann_zeta<_Tp>(2 * __k)
				* __x2pi2k / _Tp(2 * __k + __n);
	      __sum += __term;
	      if (std::abs(__term) < _S_eps * std::abs(__sum))
        	break;
	      __x2pi2k *= -__x2pi2;
	    }
	  __sum *= _Tp(__n);
	  __sum += _Tp{1} - _Tp(__n) * __x / _Tp(2 * (1 + __n));
	  return __sum;
	}
    }

_GLIBCXX_END_NAMESPACE_VERSION
} // namespace __detail
} // namespace std

#endif // _GLIBCXX_BITS_SF_ZETA_TCC
