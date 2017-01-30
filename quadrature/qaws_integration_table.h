/* quadrature/qaws_integration_table.h
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Brian Gough
 * Copyright (C) 2016-2017 Free Software Foundation, Inc.
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

#ifndef QAWS_INTEGRATION_TABLE_H
#define QAWS_INTEGRATION_TABLE_H 1

#include <array>

namespace __gnu_test
{

  /**
   * This structure manages integration of functions
   * with optional singular factors
   * @f[
   *   log(x - \alpha) log(\beta - x)
   * @f]
   */
  template<typename _Tp>
    struct qaws_integration_table
    {
      _Tp alpha;
      _Tp beta;
      int mu;
      int nu;
      std::array<_Tp, 25> ri;
      std::array<_Tp, 25> rj;
      std::array<_Tp, 25> rg;
      std::array<_Tp, 25> rh;

      qaws_integration_table(_Tp alpha_in, _Tp beta_in, int mu_in, int nu_in);
      void set(_Tp alpha_in, _Tp beta_in, int mu_in, int nu_in);
      void initialise();
    };

} // namespace __gnu_test

#include "qaws_integration_table.tcc"

#endif // QAWS_INTEGRATION_TABLE_H
