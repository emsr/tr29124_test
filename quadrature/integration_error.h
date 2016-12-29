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
// Ported from GSL by Jason Dick
// Originally written by Brian Gaugh
//
//Used to interpret the error estimate of the itnegration routines
//Based upon gsl-1.9/integration/err.c

#ifndef ERR_H
#define ERR_H 1

#include <cmath>
#include <limits>
#include <string_view>

namespace __gnu_test
{

  template<typename _Tp>
    class _IntegrationError : public std::runtime_error
    {
      _Tp _M_result;
      _Tp _M_abserr;
      int _M_status;

    public:

      _IntegrationError(const char* __what, int __status,
			_Tp __result = _Tp{0}, _Tp __abserr = _Tp{0})
      : std::runtime_error(__what),
	_M_result(__result),
	_M_abserr(__abserr),
	_M_status(__status)
      { }

      int
      status() const
      { return _M_status; }

      _Tp
      result() const
      { return this->_M_result; }

      _Tp
      abserr() const
      { return this->_M_abserr; }
    };

  /**
   * Throws appropriate error if errcode nonzero
   */
  template<typename _Tp>
    void
    __throw__IntegrationError(const char* __what, int __status,
			      _Tp __result = _Tp{0}, _Tp __abserr = _Tp{0})
    {
      _GLIBCXX_THROW_OR_ABORT(
	_IntegrationError(__what, __status, __result, __abserr));
    }

  /**
   * Throws appropriate error if errcode nonzero
   */
  template<typename _Tp>
    void
    __check_error(std::string_view __func, int __errcode,
		  _Tp __result = _Tp{0}, _Tp __abserr = _Tp{0})
    {
      std::ostringstream msg;
      msg << __func << ": ";

      if (__errcode > 2)
	--__errcode;
      switch(__errcode)
	{
	case 0:
	  return;
	case 1:
	  msg << "Number of iterations was insufficient";
	  break;
	case 2:
	  msg << "Cannot reach tolerance because of roundoff error";
	  break;
	case 3:
	  msg << "Bad integrand behavior found in the integration interval";
	  break;
	case 4:
	  msg << "Roundoff error detected in the extrapolation";
	  break;
	case 5:
	  msg << "Integral is divergent, or slowly convergent";
	  break;
        case 6:
	  msg << "Maximum number of subdivisions reached";
	  break;
	default:
	  msg << "Could not integrate function";
	}

      __throw__IntegrationError(msg.str().c_str(), __errcode,
				__result, __abserr);
    }

  template<typename _Tp>
    _Tp
    __rescale_error(_Tp __err,
		    const _Tp __result_abs, const _Tp __result_asc)
    {
      const auto _S_eps = std::numeric_limits<_Tp>::epsilon();
      const auto _S_min = std::numeric_limits<_Tp>::min();

      __err = std::abs(__err);
      if (__result_asc != 0 && __err != _Tp{0})
	{
	  auto __scale = std::pow((_Tp{200} * __err / __result_asc), _Tp{1.5});

	  if (__scale < _Tp{1})
	    __err = __result_asc * __scale;
	  else
	    __err = __result_asc;
	}
      if (__result_abs > _S_min / (_Tp{50} * _S_eps))
	{
	  auto __min_err = _Tp{50} * _S_eps * __result_abs;

	  if (__min_err > __err)
	    __err = __min_err;
	}

      return __err;
    }

} //namespace

#endif // ERR_H
