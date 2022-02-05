
// Copyright (C) 2015-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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

#ifndef HANKEL_H
#define HANKEL_H 1

#include <cmath>
#include <complex>
#include <emsr/complex128.h>
#include <emsr/sf_airy.tcc>
//#include <emsr/sf_hankel.tcc>
#include <emsr/sf_hankel_new.tcc>

namespace emsr
{

  inline std::complex<float>
  ccyl_hankel_h1f(std::complex<float> nu,
		  std::complex<float> z)
  { return emsr::detail::cyl_hankel_h1<float>(nu, z); }

  inline std::complex<long double>
  ccyl_hankel_h1l(std::complex<long double> nu,
		  std::complex<long double> z)
  { return emsr::detail::cyl_hankel_h1<long double>(nu, z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<typename emsr::promote_2<_Tpnu, _Tp>::type>
    ccyl_hankel_h1(std::complex<_Tpnu> nu, std::complex<_Tp> z)
    {
      using type = typename emsr::promote_2<_Tpnu, _Tp>::type;
      return emsr::detail::cyl_hankel_h1<type>(nu, z);
    }


  inline std::complex<float>
  ccyl_hankel_h2f(std::complex<float> nu,
		  std::complex<float> z)
  { return emsr::detail::cyl_hankel_h2<float>(nu, z); }

  inline std::complex<long double>
  ccyl_hankel_h2l(std::complex<long double> nu,
		  std::complex<long double> z)
  { return emsr::detail::cyl_hankel_h2<long double>(nu, z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<typename emsr::promote_2<_Tpnu, _Tp>::type>
    ccyl_hankel_h2(std::complex<_Tpnu> nu, std::complex<_Tp> z)
    {
      using type = typename emsr::promote_2<_Tpnu, _Tp>::type;
      return emsr::detail::cyl_hankel_h2<type>(nu, z);
    }


  inline std::complex<float>
  ccyl_besself(std::complex<float> nu,
	       std::complex<float> z)
  { return emsr::detail::cyl_bessel<float>(nu, z); }

  inline std::complex<long double>
  ccyl_bessell(std::complex<long double> nu,
	       std::complex<long double> z)
  { return emsr::detail::cyl_bessel<long double>(nu, z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<typename emsr::promote_2<_Tpnu, _Tp>::type>
    ccyl_bessel(std::complex<_Tpnu> nu, std::complex<_Tp> z)
    {
      using type = typename emsr::promote_2<_Tpnu, _Tp>::type;
      return emsr::detail::cyl_bessel<type>(nu, z);
    }


  inline std::complex<float>
  ccyl_neumannf(std::complex<float> nu,
	        std::complex<float> z)
  { return emsr::detail::cyl_neumann<float>(nu, z); }

  inline std::complex<long double>
  ccyl_neumannl(std::complex<long double> nu,
	        std::complex<long double> z)
  { return emsr::detail::cyl_neumann<long double>(nu, z); }

  template<typename _Tpnu, typename _Tp>
    inline std::complex<typename emsr::promote_2<_Tpnu, _Tp>::type>
    ccyl_neumann(std::complex<_Tpnu> nu, std::complex<_Tp> z)
    {
      using type = typename emsr::promote_2<_Tpnu, _Tp>::type;
      return emsr::detail::cyl_neumann<type>(nu, z);
    }


  inline std::complex<float>
  csph_hankel_h1f(int n, std::complex<float> z)
  {return emsr::detail::sph_hankel_h1<float>(n, z); }

  inline std::complex<long double>
  csph_hankel_h1l(int n, std::complex<long double> z)
  { return emsr::detail::sph_hankel_h1<long double>(n, z); }

  template<typename _Tp>
    inline std::complex<typename emsr::promote<_Tp>::type>
    csph_hankel_h1(int n, std::complex<_Tp> z)
    {
      using type = typename emsr::promote<_Tp>::type;
      return emsr::detail::sph_hankel_h1<type>(n, z);
    }


  inline std::complex<float>
  csph_hankel_h2f(int n, std::complex<float> z)
  { return emsr::detail::sph_hankel_h2<float>(n, z); }

  inline std::complex<long double>
  csph_hankel_h2l(int n, std::complex<long double> z)
  { return emsr::detail::sph_hankel_h2<long double>(n, z); }

  template<typename _Tp>
    inline std::complex<typename emsr::promote<_Tp>::type>
    csph_hankel_h2(int n, std::complex<_Tp> z)
    {
      using type = typename emsr::promote<_Tp>::type;
      return emsr::detail::sph_hankel_h2<type>(n, z);
    }


  inline std::complex<float>
  csph_besself(int n, std::complex<float> z)
  { return emsr::detail::sph_bessel<float>(n, z); }

  inline std::complex<long double>
  csph_bessell(int n, std::complex<long double> z)
  { return emsr::detail::sph_bessel<long double>(n, z); }

  template<typename _Tp>
    inline std::complex<typename emsr::promote<_Tp>::type>
    csph_bessel(int n, std::complex<_Tp> z)
    {
      using type = typename emsr::promote<_Tp>::type;
      return emsr::detail::sph_bessel<type>(n, z);
    }


  inline std::complex<float>
  csph_neumannf(int n, std::complex<float> z)
  { return emsr::detail::sph_neumann<float>(n, z); }

  inline std::complex<long double>
  csph_neumannl(int n, std::complex<long double> z)
  { return emsr::detail::sph_neumann<long double>(n, z); }

  template<typename _Tp>
    inline std::complex<typename emsr::promote<_Tp>::type>
    csph_neumann(int n, std::complex<_Tp> z)
    {
      using type = typename emsr::promote<_Tp>::type;
      return emsr::detail::sph_neumann<type>(n, z);
    }


  inline std::complex<float>
  cairy_aif(std::complex<float> z)
  { return emsr::detail::airy_ai<float>(z); }

  inline std::complex<long double>
  cairy_ail(std::complex<long double> z)
  { return emsr::detail::airy_ai<long double>(z); }

  template<typename _Tp>
    inline std::complex<typename emsr::promote<_Tp>::type>
    cairy_ai(std::complex<_Tp> z)
    {
      using type = typename emsr::promote<_Tp>::type;
      return emsr::detail::airy_ai<type>(z);
    }
//
//
//  inline std::complex<float>
//  cairy_bif(std::complex<float> z)
//  { return emsr::detail::airy_bi<float>(z); }
//
//  inline std::complex<long double>
//  cairy_bil(std::complex<long double> z)
//  { return emsr::detail::airy_bi<long double>(z); }
//
//  template<typename _Tp>
//    inline std::complex<typename emsr::promote<_Tp>::type>
//    cairy_bi(std::complex<_Tp> z)
//    {
//      using type = typename emsr::promote<_Tp>::type;
//      return emsr::detail::airy_bi<type>(z);
//    }

} // namespace emsr

#endif // HANKEL_H
