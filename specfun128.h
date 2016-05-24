// -*- C++ -*- header.

// Copyright (C) 2016 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// Under Section 7 of GPL version 3, you are granted additional
// permissions described in the GCC Runtime Library Exception, version
// 3.1, as published by the Free Software Foundation.

// You should have received a copy of the GNU General Public License and
// a copy of the GCC Runtime Library Exception along with this program;
// see the files COPYING3 and COPYING.RUNTIME respectively.  If not, see
// <http://www.gnu.org/licenses/>.

/** @file bits/specfun128.h
 *  This is an internal header file, included by other library headers.
 *  Do not attempt to use it directly. @headername{cmath}
 */

#ifndef SPECFUN128_H
#define SPECFUN128_H 1

#pragma GCC system_header

#if !defined(__STRICT_ANSI__) && defined(_GLIBCXX_USE_FLOAT128)

namespace std
{

  inline __float128
  assoc_laguerreq(unsigned int __n, unsigned int __m, __float128 __x)
  { return __detail::__assoc_laguerre<__float128>(__n, __m, __x); }

  inline __float128
  assoc_legendreq(unsigned int __l, unsigned int __m, __float128 __x)
  { return __detail::__assoc_legendre_p<__float128>(__l, __m, __x); }

  inline __float128
  betaq(__float128 __x, __float128 __y)
  { return __detail::__beta<__float128>(__x, __y); }

  inline __float128
  comp_ellint_1q(__float128 __k)
  { return __detail::__comp_ellint_1<__float128>(__k); }

  inline __float128
  comp_ellint_2q(__float128 __k)
  { return __detail::__comp_ellint_2<__float128>(__k); }

  inline __float128
  comp_ellint_3q(__float128 __k, __float128 __nu)
  { return __detail::__comp_ellint_3<__float128>(__k, __nu); }

  inline __float128
  cyl_bessel_iq(__float128 __nu, __float128 __x)
  { return __detail::__cyl_bessel_i<__float128>(__nu, __x); }

  inline __float128
  cyl_bessel_jq(__float128 __nu, __float128 __x)
  { return __detail::__cyl_bessel_j<__float128>(__nu, __x); }

  inline __float128
  cyl_bessel_kq(__float128 __nu, __float128 __x)
  { return __detail::__cyl_bessel_k<__float128>(__nu, __x); }

  inline __float128
  cyl_neumannq(__float128 __nu, __float128 __x)
  { return __detail::__cyl_neumann_n<__float128>(__nu, __x); }

  inline __float128
  ellint_1q(__float128 __k, __float128 __phi)
  { return __detail::__ellint_1<__float128>(__k, __phi); }

  inline __float128
  ellint_2q(__float128 __k, __float128 __phi)
  { return __detail::__ellint_2<__float128>(__k, __phi); }

  inline __float128
  ellint_3q(__float128 __k, __float128 __nu, __float128 __phi)
  { return __detail::__ellint_3<__float128>(__k, __nu, __phi); }

  inline __float128
  expintq(__float128 __x)
  { return __detail::__expint<__float128>(__x); }

  inline __float128
  hermiteq(unsigned int __n, __float128 __x)
  { return __detail::__poly_hermite<__float128>(__n, __x); }

  inline __float128
  laguerreq(unsigned int __n, __float128 __x)
  { return __detail::__laguerre<__float128>(__n, __x); }

  inline __float128
  legendreq(unsigned int __n, __float128 __x)
  { return __detail::__poly_legendre_p<__float128>(__n, __x); }

  inline __float128
  riemann_zetaq(__float128 __s)
  { return __detail::__riemann_zeta<__float128>(__s); }

  inline __float128
  sph_besselq(unsigned int __n, __float128 __x)
  { return __detail::__sph_bessel<__float128>(__n, __x); }

  inline __float128
  sph_legendreq(unsigned int __l, unsigned int __m, __float128 __theta)
  { return __detail::__sph_legendre<__float128>(__l, __m, __theta); }

  inline __float128
  sph_neumannq(unsigned int __n, __float128 __x)
  { return __detail::__sph_neumann<__float128>(__n, __x); }

} // namespace std

namespace __gnu_cxx
{

  inline __float128
  conf_hypergq(__float128 __a, __float128 __c, __float128 __x)
  { return std::__detail::__conf_hyperg<__float128>(__a, __c, __x); }

  inline __float128
  hypergq(__float128 __a, __float128 __b, __float128 __c, __float128 __x)
  { return std::__detail::__hyperg<__float128>(__a, __b, __c, __x); }

  inline __float128
  conf_hyperg_limq(__float128 __c, __float128 __x)
  { return std::__detail::__conf_hyperg_lim<__float128>(__c, __x); }

  inline __float128
  sinc_piq(__float128 __x)
  { return std::__detail::__sinc_pi<__float128>(__x); }

  inline __float128
  sincq(__float128 __x)
  { return std::__detail::__sinc<__float128>(__x); }

  inline __float128
  logintq(__float128 __x)
  { return std::__detail::__logint<__float128>(__x); }
  inline __float128
  sinintq(__float128 __x)
  { return std::__detail::__sincosint<__float128>(__x).first; }

  inline __float128
  cosintq(__float128 __x)
  { return std::__detail::__sincosint<__float128>(__x).second; }

  inline __float128
  sinhintq(__float128 __x)
  { return std::__detail::__sinhint<__float128>(__x); }

  inline __float128
  coshintq(__float128 __x)
  { return std::__detail::__coshint<__float128>(__x); }

  inline __float128
  jacobi_snq(__float128 __k, __float128 __u)
  {
    return std::get<_GLIBCXX_JACOBI_SN>
		(std::__detail::__jacobi_sncndn<__float128>(__k, __u));
  }

  inline __float128
  jacobi_cnq(__float128 __k, __float128 __u)
  {
    return std::get<_GLIBCXX_JACOBI_CN>
		(std::__detail::__jacobi_sncndn<__float128>(__k, __u));
  }

  inline __float128
  jacobi_dnq(__float128 __k, __float128 __u)
  {
    return std::get<_GLIBCXX_JACOBI_DN>
		(std::__detail::__jacobi_sncndn<__float128>(__k, __u));
  }

  inline __float128
  chebyshev_tq(unsigned int __n, __float128 __x)
  { return std::__detail::__chebyshev_t<__float128>(__n, __x); }

  inline __float128
  chebyshev_uq(unsigned int __n, __float128 __x)
  { return std::__detail::__chebyshev_u<__float128>(__n, __x); }

  inline __float128
  chebyshev_vq(unsigned int __n, __float128 __x)
  { return std::__detail::__chebyshev_v<__float128>(__n, __x); }

  inline __float128
  chebyshev_wq(unsigned int __n, __float128 __x)
  { return std::__detail::__chebyshev_w<__float128>(__n, __x); }

  inline __float128
  jacobiq(unsigned __n, __float128 __alpha, __float128 __beta, __float128 __x)
  { return std::__detail::__poly_jacobi<__float128>(__n, __alpha, __beta, __x); }

  inline __float128
  gegenbauerq(unsigned int __n, __float128 __alpha, __float128 __x)
  { return std::__detail::__gegenbauer_poly<__float128>(__n, __alpha, __x); }

  inline __float128
  zernikeq(unsigned int __n, int __m, __float128 __rho, __float128 __phi)
  { return std::__detail::__zernike<__float128>(__n, __m, __rho, __phi); }

  inline __float128
  radpolyq(unsigned int __n, unsigned int __m, __float128 __rho)
  { return std::__detail::__poly_radial_jacobi(__n, __m, __rho); }

  inline __float128
  sinhc_piq(__float128 __x)
  { return std::__detail::__sinhc_pi<__float128>(__x); }

  inline __float128
  sinhcq(__float128 __x)
  { return std::__detail::__sinhc<__float128>(__x); }

  inline std::complex<__float128>
  cyl_hankel_1q(__float128 __nu, __float128 __z)
  { return std::__detail::__cyl_hankel_1<__float128>(__nu, __z); }

  inline std::complex<__float128>
  cyl_hankel_2q(__float128 __nu, __float128 __z)
  { return std::__detail::__cyl_hankel_2<__float128>(__nu, __z); }

  inline std::complex<__float128>
  sph_hankel_1q(unsigned int __n, __float128 __z)
  { return std::__detail::__sph_hankel_1<__float128>(__n, __z); }

  inline std::complex<__float128>
  sph_hankel_2q(unsigned int __n, __float128 __z)
  { return std::__detail::__sph_hankel_2<__float128>(__n, __z); }

  inline __float128
  sph_bessel_iq(unsigned int __n, __float128 __x)
  {
    __float128 __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<__float128>(__n, __x,
        				  __i_n, __k_n, __ip_n, __kp_n);
    return __i_n;
  }

  inline __float128
  sph_bessel_kq(unsigned int __n, __float128 __x)
  {
    __float128 __i_n, __k_n, __ip_n, __kp_n;
    std::__detail::__sph_bessel_ik<__float128>(__n, __x,
        				  __i_n, __k_n, __ip_n, __kp_n);
    return __k_n;
  }

  inline __float128
  airy_aiq(__float128 __x)
  {
    __float128 __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<__float128>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Ai;
  }

  inline __float128
  airy_biq(__float128 __x)
  {
    __float128 __Ai, __Bi, __Aip, __Bip;
    std::__detail::__airy<__float128>(__x, __Ai, __Bi, __Aip, __Bip);
    return __Bi;
  }

  inline __float128
  gamma_uq(__float128 __n, __float128 __x)
  { return std::__detail::__gamma_u<__float128>(__n, __x); }

  inline __float128
  gamma_lq(__float128 __n, __float128 __x)
  { return std::__detail::__gamma_l<__float128>(__n, __x); }

  inline __float128
  digammaq(__float128 __z)
  { return std::__detail::__psi<__float128>(__z); }

  inline __float128
  dilogq(__float128 __x)
  { return std::__detail::__dilog<__float128>(__x); }

  inline __float128
  comp_ellint_rq(__float128 __x, __float128 __y)
  { return std::__detail::__comp_ellint_rf<__float128>(__x, __y); }

  inline __float128
  ellint_rfq(__float128 __x, __float128 __y, __float128 __z)
  { return std::__detail::__ellint_rf<__float128>(__x, __y, __z); }

  inline __float128
  ellint_rcq(__float128 __x, __float128 __y)
  { return std::__detail::__ellint_rc<__float128>(__x, __y); }

  inline __float128
  ellint_rjq(__float128 __x, __float128 __y, __float128 __z, __float128 __p)
  { return std::__detail::__ellint_rj<__float128>(__x, __y, __z, __p); }

  inline __float128
  ellint_rdq(__float128 __x, __float128 __y, __float128 __z)
  { return std::__detail::__ellint_rd<__float128>(__x, __y, __z); }

  inline __float128
  comp_ellint_rg(__float128 __x, __float128 __y)
  { return std::__detail::__comp_ellint_rg<__float128>(__x, __y); }

  inline __float128
  ellint_rgq(__float128 __x, __float128 __y, __float128 __z)
  { return std::__detail::__ellint_rg<__float128>(__x, __y, __z); }

  inline __float128
  hurwitz_zetaq(__float128 __s, __float128 __a)
  { return std::__detail::__hurwitz_zeta<__float128>(__s, __a); }

  inline __float128
  psiq(__float128 __x)
  { return std::__detail::__psi<__float128>(__x); }

  inline __float128
  ibetaq(__float128 __a, __float128 __b, __float128 __x)
  { return std::__detail::__beta_inc<__float128>(__a, __b, __x); }

  inline __float128
  ibetacq(__float128 __a, __float128 __b, __float128 __x)
  { return 1.0Q - ibetaq(__a, __b, __x); }

  inline __float128
  fresnel_sq(__float128 __x)
  { return std::imag(std::__detail::__fresnel<__float128>(__x)); }

  inline __float128
  fresnel_cq(__float128 __x)
  { return std::real(std::__detail::__fresnel<__float128>(__x)); }

  inline __float128
  dawsonq(__float128 __x)
  { return std::__detail::__dawson<__float128>(__x); }

  inline __float128
  expintq(unsigned int __n, __float128 __x)
  { return std::__detail::__expint<__float128>(__n, __x); }

  inline __float128
  expint_enq(unsigned int __n, __float128 __x)
  { return std::__detail::__expint<__float128>(__n, __x); }

  inline __float128
  lpochhammer_uq(__float128 __a, __float128 __n)
  { return std::__detail::__log_pochhammer_u<__float128>(__a, __n); }

  inline __float128
  lpochhammer_lq(__float128 __a, __float128 __n)
  { return std::__detail::__log_pochhammer_l<__float128>(__a, __n); }

  inline __float128
  pochhammer_uq(__float128 __a, __float128 __n)
  { return std::__detail::__pochhammer_u<__float128>(__a, __n); }

  inline __float128
  pochhammer_lq(__float128 __a, __float128 __n)
  { return std::__detail::__pochhammer_l<__float128>(__a, __n); }

  inline __float128
  factorialq(unsigned int __n)
  { return std::__detail::__factorial<__float128>(__n); }

  inline __float128
  double_factorialq(int __n)
  { return std::__detail::__double_factorial<__float128>(__n); }

  inline __float128
  lfactorialq(unsigned int __n)
  { return std::__detail::__log_factorial<__float128>(__n); }

  inline __float128
  ldouble_factorialq(int __n)
  { return std::__detail::__log_double_factorial<__float128>(__n); }

  inline __float128
  bincoefq(unsigned int __n, unsigned int __k)
  { return std::__detail::__bincoef<__float128>(__n, __k); }

  inline __float128
  lbincoefq(unsigned int __n, unsigned int __k)
  { return std::__detail::__log_bincoef<__float128>(__n, __k); }

  inline __float128
  legendre_qq(unsigned int __n, __float128 __x)
  { return std::__detail::__legendre_q<__float128>(__n, __x); }

  inline __float128
  pgammaq(__float128 __a, __float128 __x)
  { return std::__detail::__pgamma<__float128>(__a, __x); }

  inline __float128
  qgammaq(__float128 __a, __float128 __x)
  { return std::__detail::__qgamma<__float128>(__a, __x); }

  inline __float128
  jacobi_zetaq(__float128 __k, __float128 __phi)
  { return std::__detail::__jacobi_zeta<__float128>(__k, __phi); }

  inline __float128
  heuman_lambdaq(__float128 __k, __float128 __phi)
  { return std::__detail::__heuman_lambda<__float128>(__k, __phi); }

  inline __float128
  comp_ellint_dq(__float128 __k)
  { return std::__detail::__comp_ellint_d<__float128>(__k); }

  inline __float128
  ellint_dq(__float128 __k, __float128 __phi)
  { return std::__detail::__ellint_d<__float128>(__k, __phi); }

  inline __float128
  ellint_el1q(__float128 __x, __float128 __k_c)
  { return std::__detail::__ellint_el1<__float128>(__x, __k_c); }

  inline __float128
  ellint_el2q(__float128 __x, __float128 __k_c, __float128 __a, __float128 __b)
  { return std::__detail::__ellint_el2<__float128>(__x, __k_c, __a, __b); }

  inline __float128
  ellint_el3q(__float128 __x, __float128 __k_c, __float128 __p)
  { return std::__detail::__ellint_el3<__float128>(__x, __k_c, __p); }

  inline __float128
  ellint_celq(__float128 __k_c, __float128 __p, __float128 __a, __float128 __b)
  { return std::__detail::__ellint_cel<__float128>(__k_c, __p, __a, __b); }

  inline std::complex<__float128>
  cyl_hankel_1q(std::complex<__float128> __nu, std::complex<__float128> __x)
  { return std::__detail::__cyl_hankel_1<__float128>(__nu, __x); }

  inline std::complex<__float128>
  cyl_hankel_2q(std::complex<__float128> __nu, std::complex<__float128> __x)
  { return std::__detail::__cyl_hankel_2<__float128>(__nu, __x); }

  inline std::complex<__float128>
  sph_hankel_1q(unsigned int __n, std::complex<__float128> __x)
  { return std::__detail::__sph_hankel_1<__float128>(__n, __x); }

  inline std::complex<__float128>
  sph_hankel_2q(unsigned int __n, std::complex<__float128> __x)
  { return std::__detail::__sph_hankel_2<__float128>(__n, __x); }

  inline std::complex<__float128>
  sph_harmonicq(unsigned int __l, unsigned int __m,
		__float128 __theta, __float128 __phi)
  { return std::__detail::__sph_harmonic<__float128>(__l, __m, __theta, __phi); }

} // namespace __gnu_cxx

#endif // __STRICT_ANSI__ && _GLIBCXX_USE_FLOAT128

#endif // SPECFUN128_H
