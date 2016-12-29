/* quadrature/testcase.h
 *
 * Copyright (C) 2016
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

#ifndef QUADRATURE_TESTCASE_H
#define QUADRATURE_TESTCASE_H 1

template<typename _Tp>
  struct monomial
  {
    int degree;
    _Tp constant;

    monomial(int deg, _Tp c)
    : degree(deg),
      constant(c)
    { }

    _Tp
    operator()(_Tp x) const
    { return constant * std::pow(x, degree); }

    monomial
    integral() const
    { return monomial(degree + 1, constant / (degree + 1)); }

    monomial
    derivative() const
    {
      auto deg = std::max(0, degree - 1);
      return monomial(deg, deg * constant);
    }
  };

template<typename _Tp>
  _Tp
  integrate(const monomial<_Tp>& mon, _Tp a, _Tp b)
  {
    auto integ = mon.integral();
    return integ(b) - integ(a);
  }


#include "testcase.tcc"

#endif // QUADRATURE_TESTCASE_H
