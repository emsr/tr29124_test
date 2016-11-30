/* integration/tests.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "tests.h"

/* These are the test functions from table 4.1 of the QUADPACK book */

/* f1(x) = x^alpha * log(1/x) */
/* integ(f1,x,0,1) = 1/(alpha + 1)^2 */

_Tp
f1(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *) params;
  return std::pow(x,alpha) * std::log(1/x);
}

/* f2(x) = 4^-alpha / ((x-pi/4)^2 + 16^-alpha) */
/* integ(f2,x,0,1) = arctan((4-pi)4^(alpha-1)) + arctan(pi 4^(alpha-1)) */

_Tp
f2(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *) params;
  return std::pow(_Tp{4},-alpha) / (std::pow((x-M_PI/_Tp{4}),_Tp{2}) + std::pow(_Tp{16},-alpha));
}

/* f3(x) = cos(2^alpha * sin(x)) */
/* integ(f3,x,0,pi) = pi J_0(2^alpha) */

_Tp
f3(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *) params;
  return std::cos(std::pow(_Tp{2},alpha) * std::sin(x));
}

/* Functions 4, 5 and 6 are duplicates of functions  1, 2 and 3 */
/* ....                                                         */

/* f7(x) = |x - 1/3|^alpha */
/* integ(f7,x,0,1) = ((2/3)^(alpha+1) + (1/3)^(alpha+1))/(alpha + 1) */

_Tp
f7(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *) params;
  return std::pow(std::abs(x - (_Tp{1}/_Tp{3})),alpha);
}

/* f8(x) = |x - pi/4|^alpha */
/* integ(f8,x,0,1) = 
   ((1 - pi/4)^(alpha+1) + (pi/4)^(alpha+1))/(alpha + 1) */

_Tp
f8(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *) params;
  return std::pow(std::abs(x - (M_PI/_Tp{4})),alpha);
}

/* f9(x) = sqrt(1 - x^2) / (x + 1 + 2^-alpha) */
/* integ(f9,x,-1,1) = pi/sqrt((1+2^-alpha)^2-1) */

_Tp
f9(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *) params;
  return 1 / ((x + 1 + std::pow(_Tp{2},-alpha)) * std::sqrt(1-x*x));
}

/* f10(x) = std::sin(x)^(alpha - 1) */
/* integ(f10,x,0,pi/2) = 2^(alpha-2) ((Gamma(alpha/2))^2)/Gamma(alpha) */

_Tp
f10(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *)params;
  return std::pow(std::sin(x), alpha-1);
}

/* f11(x) = log(1/x)^(alpha - 1) */
/* integ(f11,x,0,1) = Gamma(alpha) */

_Tp
f11(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *)params;
  return std::pow(std::log(1/x), alpha-1);
}

/* f12(x) = std::exp(20*(x-1)) * sin(2^alpha * x) */
/* integ(f12,x,0,1) = 
   (20 sin(2^alpha) - 2^alpha cos(2^alpha) + 2^alpha exp(-20))
   /(400 + 4^alpha) */

_Tp
f12(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *)params;
  return std::exp(20 * (x - 1)) * std::sin(std::pow(_Tp{2}, alpha) * x);
}

/* f13(x) = cos(2^alpha * x)/sqrt(x(1 - x)) */
/* integ(f13,x,0,1) = pi cos(2^(alpha-1)) J_0(2^(alpha-1))  */

_Tp
f13(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *)params;
  return std::cos(std::pow(_Tp{2}, alpha) * x) / std::sqrt(x * (1 - x));
}

_Tp
f14(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *)params;
  return std::exp(-std::pow(_Tp{2}, -alpha) * x) * std::cos(x) / std::sqrt(x);
}

_Tp
f15(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *)params;
  return x * x * std::exp(-std::pow(_Tp{2}, -alpha) * x);
}

_Tp
f16(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *)params;
  if (x==0 && alpha == 1)
    return 1;  /* make the function continuous in x */
  if (x==0 && alpha > 1)
    return 0;   /* avoid problems with pow(0,1) */
  return std::pow(x, alpha - 1) / std::pow((1 + 10 * x), _Tp{2});
}

_Tp
f17(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *) params;
  return std::pow(_Tp{2}, -alpha)
  	/ (((x - 1) * (x - 1) + std::pow(_Tp{4}, -alpha)) * (x - 2));
}

/* f454(x) = x^3 log|(x^2-1)(x^2-2)| */
/* integ(f454,x,0,inf) = 61 log(2) + (77/4) log(7) - 27 */

_Tp
f454(_Tp x, void* params)
{
  _Tp x2 = x * x;
  _Tp x3 = x * x2;
  params = 0;
  return x3 * std::log(std::abs((x2 - _Tp{1}) * (x2 - _Tp{2})));
}

/* f455(x) = log(x)/(1+100*x^2) */
/* integ(f455,x,0,inf) = -log(10)/20 */

_Tp
f455(_Tp x, void* params)
{
  params = 0;
  return std::log(x) / (_Tp{1} + _Tp{100} * x * x);
}

/* f456(x) = log(x) */
/* integ(f456*std::sin(10 pi x),x,0,1) = -(gamma + log(10pi) - Ci(10pi))/(10pi) */

_Tp
f456(_Tp x, void* params)
{
  params = 0;
  if (x == _Tp{0})
    return 0;
  return std::log(x);
}

/* f457(x) = 1/sqrt(x) */
/* integ(f457*cos(pi x / 2),x,0,+inf) = 1 */

_Tp
f457(_Tp x, void* params)
{
  params = 0;
  if (x == _Tp{0})
    return 0;
  return 1 / std::sqrt(x);
}

/* f458(x) = 1/(1 + log(x)^2)^2 */
/* integ(log(x) f458(x),x,0,1) = (Ci(1) std::sin(1) + (pi/2 - Si(1)) cos(1))/pi 
                               = -0.1892752 */

_Tp
f458(_Tp x, void* params)
{
  params = 0;

  if (x == _Tp{0}) 
    return 0;
  else 
    {
      _Tp u = log(x);
      _Tp v = 1 + u * u;
      return _Tp{1} / (v * v);
    }
}

/* f459(x) = 1/(5 x^3 + 6) */
/* integ(f459/(x-0),x,-1,5) = log(125/631)/18 */

_Tp
f459(_Tp x, void* params)
{
  params = 0;
  return _Tp{1} / (_Tp{5} * x * x * x + _Tp{6});
}

/* myfn1(x) = exp(-x - x^2) */
/* integ(myfn1,x,-inf,inf) = sqrt(pi) exp(-1/4) */

_Tp
myfn1(_Tp x, void* params)
{
  params = 0;
  return std::exp(-x - x * x);
}

/* myfn2(x) = exp(alpha*x) */
/* integ(myfn2,x,-inf,b) = exp(alpha*b)/alpha */

_Tp
myfn2(_Tp x, void* params)
{
  _Tp alpha = *(_Tp *)params;
  return std::exp(alpha * x);
}


/* f_monomial = constant * x^degree */
_Tp
f_monomial(_Tp x, void* params)
{
  struct monomial_params* p = (struct monomial_params *)params;

  return p->constant * gsl_pow_int(x, p->degree);
}

/* integ(f_monomial,x,a b)=constant*(b^(degree+1)-a^(degree+1))/(degree+1) */
_Tp
integ_f_monomial(_Tp a, _Tp b, struct monomial_params* p)
{
  const int degreep1 = p->degree + 1;
  const _Tp bnp1 = gsl_pow_int(b, degreep1);
  const _Tp anp1 = gsl_pow_int(a, degreep1);
  return (p->constant / degreep1)*(bnp1 - anp1);
}

/* f(x) = std::sin(x) */
_Tp
f_sin(_Tp x, void* params)
{
    return std::sin(x);
}

/* integ(f_sin,x,a,b) */
_Tp
integ_f_sin(_Tp a, _Tp b)
{
    return -std::cos(b) + std::cos(a);
}


/* The test functions. */
_Tp
cqf1(_Tp x, void* params)
{
  return std::exp(x);
}
    
_Tp
cqf2(_Tp x, void* params)
{
  return x >= 0.3;
}

_Tp
cqf3(_Tp x, void* params)
{
  return std::sqrt(x);
}


_Tp
cqf4(_Tp x, void* params)
{
  return (_Tp{23} / _Tp{25}) * std::cosh(x) - std::cos(x);
}

_Tp
cqf5(_Tp x, void* params)
{
  _Tp x2 = x*x;
  return _Tp{1} / (x2 * (x2 + 1) + _Tp{0.9});
}

_Tp
cqf6(_Tp x, void* params)
{
  return x * std::sqrt(x);
}

_Tp
cqf7(_Tp x, void* params)
{
  return _Tp{1} / std::sqrt(x);
}

_Tp
cqf8(_Tp x, void* params)
{
  _Tp x2 = x*x;
  return _Tp{1} / (_Tp{1} + x2*x2);
}

_Tp
cqf9(_Tp x, void* params)
{
  return _Tp{2} / (_Tp{2} + std::sin(_Tp{10} * M_PI * x));
}

_Tp
cqf10(_Tp x, void* params)
{
  return _Tp{1} / (_Tp{1} + x);
}

_Tp
cqf11(_Tp x, void* params)
{
  return _Tp{1} / (_Tp{1} + std::exp(x));
}

_Tp
cqf12(_Tp x, void* params)
{
  return x / (std::exp(x) - _Tp{1});
}

_Tp
cqf13(_Tp x, void* params)
{
  return std::sin(_Tp{100} * M_PI * x) / (M_PI * x);
}

_Tp
cqf14(_Tp x, void* params)
{
  return std::sqrt(_Tp{50}) * std::exp(-_Tp{50} * M_PI * x * x);
}

_Tp
cqf15(_Tp x, void* params)
{
  return _Tp{25} * std::exp(-_Tp{25} * x);
}

_Tp
cqf16(_Tp x, void* params)
{
  return _Tp{50} / M_PI * (_Tp{2500} * x * x + _Tp{1});
}

_Tp
cqf17(_Tp x, void* params)
{
  _Tp t1 = _Tp{50} * M_PI * x, t2;
  t2 = std::sin(t1) / t1;
  return _Tp{50} * t2 * t2;
}

_Tp
cqf18(_Tp x, void* params)
{
  return std::cos(std::cos(x)
		+ _Tp{3} * std::sin(x)
		+ _Tp{2} * std::cos(_Tp{2} * x)
		+ _Tp{3} * std::sin(_Tp{2} * x)
		+ _Tp{3} * std::cos(_Tp{3} * x));
}

_Tp
cqf19(_Tp x, void* params)
{
  return std::log(x);
}

_Tp
cqf20(_Tp x, void* params)
{
  return _Tp{1} / (x * x + _Tp{1.005});
}

_Tp
cqf21(_Tp x, void* params)
{
  return _Tp{1} / std::cosh(_Tp{10} * (x - 0.2) * _Tp{2}) + 
    _Tp{1} / std::cosh(_Tp{100} * (x - 0.4) * _Tp{4}) + 
    _Tp{1} / std::cosh(_Tp{1000} * (x - 0.6) * _Tp{8});
}

_Tp
cqf22(_Tp x, void* params)
{
  return _Tp{4} * M_PI * M_PI * x
	* std::sin(_Tp{20} * M_PI * x) * std::cos(_Tp{2} * M_PI * x);
}

_Tp
cqf23(_Tp x, void* params)
{
  _Tp t = _Tp{230} * x - _Tp{30};
  return _Tp{1} / (_Tp{1} + t * t);
}

_Tp
cqf24(_Tp x, void* params)
{
  return std::floor(std::exp(x));
}

_Tp
cqf25(_Tp x, void* params)
{
  return (x < _Tp{1}) * (x + _Tp{1})
     + (_Tp{1} <= x && x <= _Tp{3}) * (_Tp{3} - x)
     + (x > _Tp{3}) * _Tp{2};
}
