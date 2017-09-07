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

#ifndef GAUSS_LAGUERRE_INTEGRATE_H
#define GAUSS_LAGUERRE_INTEGRATE_H 1

#include <vector>

namespace __gnu_ext
{

template<typename _Tp, typename _FuncTp>
  _Tp
  gauss_laguerre_integrate(const _FuncTp& __func,
			   unsigned int __n, _Tp __alpha)
  {
    if(__n == 0)
      std::__throw_domain_error("gauss_laguerre_integrate: "
    				"laguerre order must be greater than 0");
    else
     {
	auto __rule = std::__detail::__laguerre_zeros(__n, __alpha);
	auto __sum = _Tp{0};
	for (const auto& __pt : __rule)
	  {
	    auto __x = __pt.__zero;
	    auto __w = __pt.__weight;
	    __sum += __w * __func(__x);
	  }
	return __sum;
      }
  }

} // namespace __gnu_ext

#endif // GAUSS_LAGUERRE_INTEGRATE_H
