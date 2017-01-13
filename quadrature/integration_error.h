// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2017 Free Software Foundation, Inc.
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
// Used to interpret the error estimate of the itnegration routines
// Based on gsl/integration/err.c

#ifndef INTEGRATION_ERROR_H
#define INTEGRATION_ERROR_H 1

#include <cmath>
#include <limits>
#include <string_view>

namespace __gnu_test
{

  enum
  {
    NO_ERROR,
    MAX_ITER_ERROR,
    ROUNDOFF_ERROR,
    SINGULAR_ERROR,
    EXTRAP_ROUNDOFF_ERROR,
    DIVERGENCE_ERROR,
    MAX_SUBDIV_ERROR,
    TOLERANCE_ERROR,
    UNKNOWN_ERROR
  };

  template<typename _Tp>
    class _IntegrationError : public std::runtime_error
    {
      _Tp _M_result;
      _Tp _M_abserr;
      int _M_errcode;

    public:

      _IntegrationError(const char* __what, int __errcode,
			_Tp __result = std::numeric_limits<_Tp>::quiet_NaN(),
			_Tp __abserr = std::numeric_limits<_Tp>::quiet_NaN())
      : std::runtime_error(__what),
	_M_result(__result),
	_M_abserr(__abserr),
	_M_errcode(__errcode)
      { }

      int
      error_code() const
      { return _M_errcode; }

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
    __throw__IntegrationError(const char* __what, int __errcode,
			      _Tp __result = std::numeric_limits<_Tp>::quiet_NaN(),
			      _Tp __abserr = std::numeric_limits<_Tp>::quiet_NaN())
    {
      _GLIBCXX_THROW_OR_ABORT(
	_IntegrationError(__what, __errcode, __result, __abserr));
    }

  /**
   * Throws appropriate error if errcode nonzero
   */
  template<typename _Tp>
    void
    __check_error(std::string_view __func, int __errcode,
		  _Tp __result = std::numeric_limits<_Tp>::quiet_NaN(),
		  _Tp __abserr = std::numeric_limits<_Tp>::quiet_NaN())
    {
      std::ostringstream msg;
      msg << __func << ": ";

      if (__errcode > 2)
	--__errcode;
      switch(__errcode)
	{
	case NO_ERROR:
	  return;
	case MAX_ITER_ERROR:
	  msg << "Number of iterations was insufficient";
	  break;
	case ROUNDOFF_ERROR:
	  msg << "Cannot reach tolerance because of roundoff error";
	  break;
	case SINGULAR_ERROR:
	  msg << "Bad integrand behavior found in the integration interval";
	  break;
	case EXTRAP_ROUNDOFF_ERROR:
	  msg << "Roundoff error detected in the extrapolation";
	  break;
	case DIVERGENCE_ERROR:
	  msg << "Integral is divergent, or slowly convergent";
	  break;
	case MAX_SUBDIV_ERROR:
	  msg << "Maximum number of subdivisions reached";
	  break;
        case TOLERANCE_ERROR:
	  msg << "Cannot reach tolerance with maximum order rule";
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

} // namespace __gnu_test

#endif // INTEGRATION_ERROR_H
