// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2016 Free Software Foundation, Inc.
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
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.
//
// Ported from GSL by Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// This file implements an oscillatory integrand table for use
// in integration schemes
// Based upon gsl-2.1/integration/qmomof.c

#ifndef OSCILLATORY_INTEGRATION_TABLE_TCC
#define OSCILLATORY_INTEGRATION_TABLE_TCC 1

namespace __gnu_test
{
  template<typename _Tp>
    int
    _S_tridiag(size_t n, _Tp *c, _Tp *d, _Tp *e, _Tp *b);
 
  /**
   * Compute Chebyshev moments at level @c level.
   */
  template<typename _Tp>
    void
    oscillatory_integration_table<_Tp>::
    compute_moments(_Tp par, std::size_t __level)
    {
      _Tp v[28], d[25], d1[25], d2[25];

      const size_t noeq = 25;

      const auto par2 = par * par;
      const auto par4 = par2 * par2;
      const auto par22 = par2 + 2.0;

      const auto sinpar = std::sin(par);
      const auto cospar = std::cos(par);

      //
      // Compute the Chebyschev moments with respect to cosine.
      //

      auto ac = 8 * cospar;
      auto as = 24 * par * sinpar;

      v[0] = 2 * sinpar / par;
      v[1] = (8 * cospar + (2 * par2 - 8) * sinpar / par) / par2;
      v[2] = (32 * (par2 - 12) * cospar
              + (2 * ((par2 - 80) * par2 + 192) * sinpar) / par) / par4;

      if (std::abs(par) <= 24)
	{
	  // Compute the moments as the solution of a boundary value
          // problem using the asyptotic expansion as an endpoint.
	  _Tp an = 6;
	  for (auto k = 0u; k < noeq - 1; ++k)
            {
              auto an2 = an * an;
              d[k] = -2 * (an2 - 4) * (par22 - 2 * an2);
              d2[k] = (an - 1) * (an - 2) * par2;
              d1[k + 1] = (an + 3) * (an + 4) * par2;
              v[k + 3] = as - (an2 - 4) * ac;
              an += _Tp{2};
            }

	  auto an2 = an * an;

	  d[noeq - 1] = -2 * (an2 - 4) * (par22 - 2 * an2);
	  v[noeq + 2] = as - (an2 - 4) * ac;
	  v[3] = v[3] - 56 * par2 * v[2];

	  auto ass = par * sinpar;
	  auto asap = (((((210 * par2 - 1) * cospar - (105 * par2 - 63) * ass) / an2
                    - (1 - 15 * par2) * cospar + 15 * ass) / an2 
        	   - cospar + 3 * ass) / an2 
        	  - cospar) / an2;
	  v[noeq + 2] = v[noeq + 2] - 2 * asap * par2 * (an - 1) * (an - 2);

	  _S_tridiag(noeq, d1, d, d2, v + 3);
	}
      else
	{
	  // Compute the moments by forward recursion.
	  auto an = _Tp{4};
	  for (auto k = 3u; k < 13u; ++k)
            {
              _Tp an2 = an * an;
              v[k] = ((an2 - 4) * (2 * (par22 - 2 * an2) * v[k - 1] - ac)
                      + as - par2 * (an + 1) * (an + 2) * v[k - 2]) 
        	/ (par2 * (an - 1) * (an - 2));
              an += _Tp{2};
            }
	}


      for (auto i = 0u; i < 13u; ++i)
	this->chebmo[25 * __level + 2 * i] = v[i];

      //
      // Compute the Chebyschev moments with respect to sine.
      //

      v[0] = 2 * (sinpar - par * cospar) / par2;
      v[1] = (18 - 48 / par2) * sinpar / par2 + (-2 + 48 / par2) * cospar / par;

      ac = -24 * par * cospar;
      as = -8 * sinpar;

      if (std::abs(par) <= 24)
	{
	  // Compute the moments as the solution of a boundary value
          // problem using the asyptotic expansion as an endpoint.
	  auto an = _Tp{5};
	  for (auto k = 0u; k < noeq - 1; ++k)
            {
              auto an2 = an * an;
              d[k] = -2 * (an2 - 4) * (par22 - 2 * an2);
              d2[k] = (an - 1) * (an - 2) * par2;
              d1[k + 1] = (an + 3) * (an + 4) * par2;
              v[k + 2] = ac + (an2 - 4) * as;
              an += _Tp{2};
            }
	  auto an2 = an * an;

	  d[noeq - 1] = -2 * (an2 - 4) * (par22 - 2 * an2);
	  v[noeq + 1] = ac + (an2 - 4) * as;
	  v[2] = v[2] - 42 * par2 * v[1];

	  auto ass = par * cospar;
	  auto asap = (((((105 * par2 - 63) * ass - (210 * par2 - 1) * sinpar) / an2
                    + (15 * par2 - 1) * sinpar
                    - 15 * ass) / an2 - sinpar - 3 * ass) / an2 - sinpar) / an2;
	  v[noeq + 1] = v[noeq + 1] - 2 * asap * par2 * (an - 1) * (an - 2);

	  _S_tridiag(noeq, d1, d, d2, v + 2);
	}
      else
	{
	  // Compute the moments by forward recursion.
	  auto an = _Tp{3};
	  for (auto k = 2u; k < 12u; ++k)
            {
              _Tp an2 = an * an;
              v[k] = ((an2 - _Tp{4}) * (_Tp{2} * (par22 - _Tp{2} * an2) * v[k - 1] + as)
                      + ac - par2 * (an + _Tp{1}) * (an + _Tp{2}) * v[k - 2]) 
        	/ (par2 * (an - 1) * (an - _Tp{2}));
              an += _Tp{2};
            }
	}

      for (auto i = 0u; i < 12u; ++i)
	this->chebmo[25 * __level + 2 * i + 1] = v[i];
    }

  /**
   * Solve a tridiagonal system A x = b:
   *
   *   c[1 .. n - 1]   subdiagonal of the matrix A
   *   d[0 .. n - 1]   diagonal of the matrix A
   *   e[0 .. n - 2]   superdiagonal of the matrix A
   *
   *   b[0 .. n - 1]   right hand side, replaced by the solution vector x
   */
  template<typename _Tp>
    int
    _S_tridiag(size_t n, _Tp *c, _Tp *d, _Tp *e, _Tp *b)
    {
      /* Solves a tridiagonal matrix
       */

      size_t k;

      c[0] = d[0];

      if (n == 0)
	return 0;

      if (n == 1)
	{
	  b[0] = b[0] / d[0];
	  return 0;
	}

      d[0] = e[0];
      e[0] = 0;
      e[n - 1] = 0;

      for (auto k = 0u; k < n - 1; ++k)
	{
	  size_t k1 = k + 1;

	  if (std::abs(c[k1]) >= std::abs(c[k]))
            {
              std::swap(c[k1], c[k]);
              std::swap(d[k1], d[k]);
	      std::swap(e[k1], e[k]);
              std::swap(b[k1], b[k]);
            }

	  if (c[k] == 0)
            return 1;

	  {
            double t = -c[k1] / c[k];

            c[k1] = d[k1] + t * d[k];
            d[k1] = e[k1] + t * e[k];
            e[k1] = 0;
            b[k1] = b[k1] + t * b[k];
	  }
	}

      if (c[n - 1] == 0)
	return 1;

      b[n - 1] = b[n - 1] / c[n - 1];

      b[n - 2] = (b[n - 2] - d[n - 2] * b[n - 1]) / c[n - 2];

      for (k = n ; k > 2; --k)
	{
	  size_t kb = k - 3;
	  b[kb] = (b[kb] - d[kb] * b[kb + 1] - e[kb] * b[kb + 2]) / c[kb];
	}

      return 0;
    }

} // namespace __gnu_test

#endif // OSCILLATORY_INTEGRATION_TABLE_TCC
