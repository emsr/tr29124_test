// TR29124 math special functions -*- C++ -*-

// Copyright (C) 2015 Free Software Foundation, Inc.
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

/** @file bits/hankel.h
 *  This is an internal header file, included by other library headers.
 *  You should not attempt to use it directly.
 */

#ifndef _GLIBCXX_BITS_HANKEL_H
#define _GLIBCXX_BITS_HANKEL_H 1

#pragma GCC system_header

#include <cmath>
//#include <bits/sf_airy.tcc>
//#include <bits/sf_hankel.tcc>
#include "sf_airy.tcc"
#include "sf_hankel_new.tcc"

namespace __gnu_cxx _GLIBCXX_VISIBILITY(default)
{

  inline std::complex<float>
  ccyl_hankel_h1f(std::complex<float> __nu,
		  std::complex<float> __z)
  { return std::__detail::__cyl_hankel_h1<float>(__nu, __z); }

  inline std::complex<long double>
  ccyl_hankel_h1l(std::complex<long double> __nu,
		  std::complex<long double> __z)
  { return std::__detail::__cyl_hankel_h1<long double>(__nu, __z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type>
    ccyl_hankel_h1(std::complex<_Tpnu> __nu, std::complex<_Tp> __z)
    {
      using __type = typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type;
      return std::__detail::__cyl_hankel_h1<__type>(__nu, __z);
    }


  inline std::complex<float>
  ccyl_hankel_h2f(std::complex<float> __nu,
		  std::complex<float> __z)
  { return std::__detail::__cyl_hankel_h2<float>(__nu, __z); }

  inline std::complex<long double>
  ccyl_hankel_h2l(std::complex<long double> __nu,
		  std::complex<long double> __z)
  { return std::__detail::__cyl_hankel_h2<long double>(__nu, __z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type>
    ccyl_hankel_h2(std::complex<_Tpnu> __nu, std::complex<_Tp> __z)
    {
      using __type = typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type;
      return std::__detail::__cyl_hankel_h2<__type>(__nu, __z);
    }


  inline std::complex<float>
  ccyl_besself(std::complex<float> __nu,
	       std::complex<float> __z)
  { return std::__detail::__cyl_bessel<float>(__nu, __z); }

  inline std::complex<long double>
  ccyl_bessell(std::complex<long double> __nu,
	       std::complex<long double> __z)
  { return std::__detail::__cyl_bessel<long double>(__nu, __z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type>
    ccyl_bessel(std::complex<_Tpnu> __nu, std::complex<_Tp> __z)
    {
      using __type = typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type;
      return std::__detail::__cyl_bessel<__type>(__nu, __z);
    }


  inline std::complex<float>
  ccyl_neumannf(std::complex<float> __nu,
	        std::complex<float> __z)
  { return std::__detail::__cyl_neumann<float>(__nu, __z); }

  inline std::complex<long double>
  ccyl_neumannl(std::complex<long double> __nu,
	        std::complex<long double> __z)
  { return std::__detail::__cyl_neumann<long double>(__nu, __z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type>
    ccyl_neumann(std::complex<_Tpnu> __nu, std::complex<_Tp> __z)
    {
      using __type = typename __gnu_cxx::__promote_2<_Tpnu, _Tp>::__type;
      return std::__detail::__cyl_neumann<__type>(__nu, __z);
    }


  inline std::complex<float>
  csph_hankel_h1f(int __n, std::complex<float> __z)
  {return std::__detail::__sph_hankel_h1<float>(__n, __z); }

  inline std::complex<long double>
  csph_hankel_h1l(int __n, std::complex<long double> __z)
  { return std::__detail::__sph_hankel_h1<long double>(__n, __z); }

  template<typename _Tp>
    inline std::complex<typename __gnu_cxx::__promote<_Tp>::__type>
    csph_hankel_h1(int __n, std::complex<_Tp> __z)
    {
      using __type = typename __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__sph_hankel_h1<__type>(__n, __z);
    }


  inline std::complex<float>
  csph_hankel_h2f(int __n, std::complex<float> __z)
  { return std::__detail::__sph_hankel_h2<float>(__n, __z); }

  inline std::complex<long double>
  csph_hankel_h2l(int __n, std::complex<long double> __z)
  { return std::__detail::__sph_hankel_h2<long double>(__n, __z); }

  template<typename _Tp>
    inline std::complex<typename __gnu_cxx::__promote<_Tp>::__type>
    csph_hankel_h2(int __n, std::complex<_Tp> __z)
    {
      using __type = typename __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__sph_hankel_h2<__type>(__n, __z);
    }


  inline std::complex<float>
  csph_besself(int __n, std::complex<float> __z)
  { return std::__detail::__sph_bessel<float>(__n, __z); }

  inline std::complex<long double>
  csph_bessell(int __n, std::complex<long double> __z)
  { return std::__detail::__sph_bessel<long double>(__n, __z); }

  template<typename _Tp>
    inline std::complex<typename __gnu_cxx::__promote<_Tp>::__type>
    csph_bessel(int __n, std::complex<_Tp> __z)
    {
      using __type = typename __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__sph_bessel<__type>(__n, __z);
    }


  inline std::complex<float>
  csph_neumannf(int __n, std::complex<float> __z)
  { return std::__detail::__sph_neumann<float>(__n, __z); }

  inline std::complex<long double>
  csph_neumannl(int __n, std::complex<long double> __z)
  { return std::__detail::__sph_neumann<long double>(__n, __z); }

  template<typename _Tp>
    inline std::complex<typename __gnu_cxx::__promote<_Tp>::__type>
    csph_neumann(int __n, std::complex<_Tp> __z)
    {
      using __type = typename __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__sph_neumann<__type>(__n, __z);
    }


  inline std::complex<float>
  cairy_aif(std::complex<float> __z)
  { return std::__detail::__airy_ai<float>(__z); }

  inline std::complex<long double>
  cairy_ail(std::complex<long double> __z)
  { return std::__detail::__airy_ai<long double>(__z); }

  template<typename _Tp>
    inline std::complex<typename __gnu_cxx::__promote<_Tp>::__type>
    cairy_ai(std::complex<_Tp> __z)
    {
      using __type = typename __gnu_cxx::__promote<_Tp>::__type;
      return std::__detail::__airy_ai<__type>(__z);
    }
//
//
//  inline std::complex<float>
//  cairy_bif(std::complex<float> __z)
//  { return std::__detail::__airy_bi<float>(__z); }
//
//  inline std::complex<long double>
//  cairy_bil(std::complex<long double> __z)
//  { return std::__detail::__airy_bi<long double>(__z); }
//
//  template<typename _Tp>
//    inline std::complex<typename __gnu_cxx::__promote<_Tp>::__type>
//    cairy_bi(std::complex<_Tp> __z)
//    {
//      using __type = typename __gnu_cxx::__promote<_Tp>::__type;
//      return std::__detail:__airy_bi<__type>(__z);
//    }

} // namespace __gnu_cxx

#endif // _GLIBCXX_BITS_HANKEL_H
