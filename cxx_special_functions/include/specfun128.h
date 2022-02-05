
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.

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

#ifdef EMSR_HAVE_FLOAT128

namespace emsr
{

  inline __float128
  assoc_laguerreq(unsigned int n, unsigned int m, __float128 x)
  { return emsr::detail::assoc_laguerre<__float128>(n, m, x); }

  inline __float128
  assoc_legendreq(unsigned int l, unsigned int m, __float128 x)
  { return emsr::detail::assoc_legendre_p<__float128>(l, m, x); }

  inline __float128
  betaq(__float128 x, __float128 y)
  { return emsr::detail::beta<__float128>(x, y); }

  inline __float128
  comp_ellint_1q(__float128 k)
  { return emsr::detail::comp_ellint_1<__float128>(k); }

  inline __float128
  comp_ellint_2q(__float128 k)
  { return emsr::detail::comp_ellint_2<__float128>(k); }

  inline __float128
  comp_ellint_3q(__float128 k, __float128 nu)
  { return emsr::detail::comp_ellint_3<__float128>(k, nu); }

  inline __float128
  cyl_bessel_iq(__float128 nu, __float128 x)
  { return emsr::detail::cyl_bessel_i<__float128>(nu, x); }

  inline __float128
  cyl_bessel_jq(__float128 nu, __float128 x)
  { return emsr::detail::cyl_bessel_j<__float128>(nu, x); }

  inline __float128
  cyl_bessel_kq(__float128 nu, __float128 x)
  { return emsr::detail::cyl_bessel_k<__float128>(nu, x); }

  inline __float128
  cyl_neumannq(__float128 nu, __float128 x)
  { return emsr::detail::cyl_neumann_n<__float128>(nu, x); }

  inline __float128
  ellint_1q(__float128 k, __float128 phi)
  { return emsr::detail::ellint_1<__float128>(k, phi); }

  inline __float128
  ellint_2q(__float128 k, __float128 phi)
  { return emsr::detail::ellint_2<__float128>(k, phi); }

  inline __float128
  ellint_3q(__float128 k, __float128 nu, __float128 phi)
  { return emsr::detail::ellint_3<__float128>(k, nu, phi); }

  inline __float128
  expintq(__float128 x)
  { return emsr::detail::expint<__float128>(x); }

  inline __float128
  hermiteq(unsigned int n, __float128 x)
  { return emsr::detail::hermite<__float128>(n, x); }

  inline __float128
  laguerreq(unsigned int n, __float128 x)
  { return emsr::detail::laguerre<__float128>(n, x); }

  inline __float128
  legendreq(unsigned int n, __float128 x)
  { return emsr::detail::legendre_p<__float128>(n, x); }

  inline __float128
  riemann_zetaq(__float128 s)
  { return emsr::detail::riemann_zeta<__float128>(s); }

  inline __float128
  sph_besselq(unsigned int n, __float128 x)
  { return emsr::detail::sph_bessel<__float128>(n, x); }

  inline __float128
  sph_legendreq(unsigned int l, unsigned int m, __float128 theta)
  { return emsr::detail::sph_legendre<__float128>(l, m, theta); }

  inline __float128
  sph_neumannq(unsigned int n, __float128 x)
  { return emsr::detail::sph_neumann<__float128>(n, x); }

  inline __float128
  conf_hypergq(__float128 a, __float128 c, __float128 x)
  { return emsr::detail::conf_hyperg<__float128>(a, c, x); }

  inline __float128
  hypergq(__float128 a, __float128 b, __float128 c, __float128 x)
  { return emsr::detail::hyperg<__float128>(a, b, c, x); }

  inline __float128
  conf_hyperg_limq(__float128 c, __float128 x)
  { return emsr::detail::conf_hyperg_lim<__float128>(c, x); }

  inline __float128
  sinc_piq(__float128 x)
  { return emsr::detail::sinc_pi<__float128>(x); }

  inline __float128
  sincq(__float128 x)
  { return emsr::detail::sinc<__float128>(x); }

  inline __float128
  sinhc_piq(__float128 x)
  { return emsr::detail::sinhc_pi<__float128>(x); }

  inline __float128
  sinhcq(__float128 x)
  { return emsr::detail::sinhc<__float128>(x); }

  inline __float128
  logintq(__float128 x)
  { return emsr::detail::logint<__float128>(x); }
  inline __float128
  sinintq(__float128 x)
  { return emsr::detail::sincosint<__float128>(x).first; }

  inline __float128
  cosintq(__float128 x)
  { return emsr::detail::sincosint<__float128>(x).second; }

  inline __float128
  sinhintq(__float128 x)
  { return emsr::detail::sinhint<__float128>(x); }

  inline __float128
  coshintq(__float128 x)
  { return emsr::detail::coshint<__float128>(x); }

  inline __float128
  jacobi_snq(__float128 k, __float128 u)
  {
    return std::get<_GLIBCXX_JACOBI_SN>
		(emsr::detail::jacobi_ellint<__float128>(k, u));
  }

  inline __float128
  jacobi_cnq(__float128 k, __float128 u)
  {
    return std::get<_GLIBCXX_JACOBI_CN>
		(emsr::detail::jacobi_ellint<__float128>(k, u));
  }

  inline __float128
  jacobi_dnq(__float128 k, __float128 u)
  {
    return std::get<_GLIBCXX_JACOBI_DN>
		(emsr::detail::jacobi_ellint<__float128>(k, u));
  }

  inline __float128
  chebyshev_tq(unsigned int n, __float128 x)
  { return emsr::detail::chebyshev_t<__float128>(n, x); }

  inline __float128
  chebyshev_uq(unsigned int n, __float128 x)
  { return emsr::detail::chebyshev_u<__float128>(n, x); }

  inline __float128
  chebyshev_vq(unsigned int n, __float128 x)
  { return emsr::detail::chebyshev_v<__float128>(n, x); }

  inline __float128
  chebyshev_wq(unsigned int n, __float128 x)
  { return emsr::detail::chebyshev_w<__float128>(n, x); }

  inline __float128
  jacobiq(unsigned n, __float128 alpha, __float128 beta, __float128 x)
  { return emsr::detail::poly_jacobi<__float128>(n, alpha, beta, x); }

  inline __float128
  gegenbauerq(unsigned int n, __float128 alpha, __float128 x)
  { return emsr::detail::gegenbauer_poly<__float128>(n, alpha, x); }

  inline __float128
  zernikeq(unsigned int n, int m, __float128 rho, __float128 phi)
  { return emsr::detail::zernike<__float128>(n, m, rho, phi); }

  inline __float128
  radpolyq(unsigned int n, unsigned int m, __float128 rho)
  { return emsr::detail::poly_radial_jacobi(n, m, rho); }

  inline __float128
  sinhc_piq(__float128 x)
  { return emsr::detail::sinhc_pi<__float128>(x); }

  inline __float128
  sinhcq(__float128 x)
  { return emsr::detail::sinhc<__float128>(x); }

  inline std::complex<__float128>
  cyl_hankel_1q(__float128 nu, __float128 z)
  { return emsr::detail::cyl_hankel_1<__float128>(nu, z); }

  inline std::complex<__float128>
  cyl_hankel_2q(__float128 nu, __float128 z)
  { return emsr::detail::cyl_hankel_2<__float128>(nu, z); }

  inline std::complex<__float128>
  sph_hankel_1q(unsigned int n, __float128 z)
  { return emsr::detail::sph_hankel_1<__float128>(n, z); }

  inline std::complex<__float128>
  sph_hankel_2q(unsigned int n, __float128 z)
  { return emsr::detail::sph_hankel_2<__float128>(n, z); }

  inline __float128
  sph_bessel_iq(unsigned int n, __float128 x)
  {
    __float128 i_n, k_n, ip_n, kp_n;
    emsr::detail::sph_bessel_ik<__float128>(n, x,
        				  i_n, k_n, ip_n, kp_n);
    return __i_n;
  }

  inline __float128
  sph_bessel_kq(unsigned int n, __float128 x)
  {
    __float128 i_n, k_n, ip_n, kp_n;
    emsr::detail::sph_bessel_ik<__float128>(n, x,
        				  i_n, k_n, ip_n, kp_n);
    return __k_n;
  }

  inline __float128
  airy_aiq(__float128 x)
  {
    __float128 Ai, Bi, Aip, Bip;
    emsr::detail::airy<__float128>(x, Ai, Bi, Aip, Bip);
    return __Ai;
  }

  inline __float128
  airy_biq(__float128 x)
  {
    __float128 Ai, Bi, Aip, Bip;
    emsr::detail::airy<__float128>(x, Ai, Bi, Aip, Bip);
    return __Bi;
  }

  inline __float128
  tgammaq(__float128 n, __float128 x)
  { return emsr::detail::tgamma<__float128>(n, x); }

  inline __float128
  tgamma_lowerq(__float128 n, __float128 x)
  { return emsr::detail::tgamma_lower<__float128>(n, x); }

  inline __float128
  dilogq(__float128 x)
  { return emsr::detail::dilog<__float128>(x); }

  inline __float128
  comp_ellint_rq(__float128 x, __float128 y)
  { return emsr::detail::comp_ellint_rf<__float128>(x, y); }

  inline __float128
  ellint_rfq(__float128 x, __float128 y, __float128 z)
  { return emsr::detail::ellint_rf<__float128>(x, y, z); }

  inline __float128
  ellint_rcq(__float128 x, __float128 y)
  { return emsr::detail::ellint_rc<__float128>(x, y); }

  inline __float128
  ellint_rjq(__float128 x, __float128 y, __float128 z, __float128 p)
  { return emsr::detail::ellint_rj<__float128>(x, y, z, p); }

  inline __float128
  ellint_rdq(__float128 x, __float128 y, __float128 z)
  { return emsr::detail::ellint_rd<__float128>(x, y, z); }

  inline __float128
  comp_ellint_rg(__float128 x, __float128 y)
  { return emsr::detail::comp_ellint_rg<__float128>(x, y); }

  inline __float128
  ellint_rgq(__float128 x, __float128 y, __float128 z)
  { return emsr::detail::ellint_rg<__float128>(x, y, z); }

  inline __float128
  hurwitz_zetaq(__float128 s, __float128 a)
  { return emsr::detail::hurwitz_zeta<__float128>(s, a); }

  inline __float128
  digammaq(__float128 x)
  { return emsr::detail::digamma<__float128>(x); }

  inline __float128
  polygammaq(unsigned int m, __float128 x)
  { return emsr::detail::polygamma<__float128>(m, x); }

  inline __float128
  ibetaq(__float128 a, __float128 b, __float128 x)
  { return emsr::detail::beta_inc<__float128>(a, b, x); }

  inline __float128
  ibetacq(__float128 a, __float128 b, __float128 x)
  { return 1.0Q - ibetaq(a, b, x); }

  inline __float128
  fresnel_sq(__float128 x)
  { return std::imag(emsr::detail::fresnel<__float128>(x)); }

  inline __float128
  fresnel_cq(__float128 x)
  { return std::real(emsr::detail::fresnel<__float128>(x)); }

  inline __float128
  dawsonq(__float128 x)
  { return emsr::detail::dawson<__float128>(x); }

  inline __float128
  expintq(unsigned int n, __float128 x)
  { return emsr::detail::expint<__float128>(n, x); }

  inline __float128
  expint_enq(unsigned int n, __float128 x)
  { return emsr::detail::expint<__float128>(n, x); }

  inline __float128
  lrising_factorialq(__float128 a, __float128 n)
  { return emsr::detail::log_rising_factorial<__float128>(a, n); }

  inline __float128
  lfalling_factorialq(__float128 a, __float128 n)
  { return emsr::detail::log_falling_factorial<__float128>(a, n); }

  inline __float128
  rising_factorialq(__float128 a, __float128 n)
  { return emsr::detail::rising_factorial<__float128>(a, n); }

  inline __float128
  falling_factorialq(__float128 a, __float128 n)
  { return emsr::detail::falling_factorial<__float128>(a, n); }

  inline __float128
  factorialq(unsigned int n)
  { return emsr::detail::factorial<__float128>(n); }

  inline __float128
  double_factorialq(int n)
  { return emsr::detail::double_factorial<__float128>(n); }

  inline __float128
  lfactorialq(unsigned int n)
  { return emsr::detail::log_factorial<__float128>(n); }

  inline __float128
  ldouble_factorialq(int n)
  { return emsr::detail::log_double_factorial<__float128>(n); }

  inline __float128
  binomialq(unsigned int n, unsigned int k)
  { return emsr::detail::binomial<__float128>(n, k); }

  inline __float128
  lbinomialq(unsigned int n, unsigned int k)
  { return emsr::detail::log_binomial<__float128>(n, k); }

  inline __float128
  legendre_qq(unsigned int n, __float128 x)
  { return emsr::detail::legendre_q<__float128>(n, x); }

  inline __float128
  gamma_pq(__float128 a, __float128 x)
  { return emsr::detail::gamma_p<__float128>(a, x); }

  inline __float128
  gamma_qq(__float128 a, __float128 x)
  { return emsr::detail::gamma_q<__float128>(a, x); }

  inline __float128
  jacobi_zetaq(__float128 k, __float128 phi)
  { return emsr::detail::jacobi_zeta<__float128>(k, phi); }

  inline __float128
  heuman_lambdaq(__float128 k, __float128 phi)
  { return emsr::detail::heuman_lambda<__float128>(k, phi); }

  inline __float128
  comp_ellint_dq(__float128 k)
  { return emsr::detail::comp_ellint_d<__float128>(k); }

  inline __float128
  ellint_dq(__float128 k, __float128 phi)
  { return emsr::detail::ellint_d<__float128>(k, phi); }

  inline __float128
  ellint_el1q(__float128 x, __float128 k_c)
  { return emsr::detail::ellint_el1<__float128>(x, k_c); }

  inline __float128
  ellint_el2q(__float128 x, __float128 k_c, __float128 a, __float128 b)
  { return emsr::detail::ellint_el2<__float128>(x, k_c, a, b); }

  inline __float128
  ellint_el3q(__float128 x, __float128 k_c, __float128 p)
  { return emsr::detail::ellint_el3<__float128>(x, k_c, p); }

  inline __float128
  ellint_celq(__float128 k_c, __float128 p, __float128 a, __float128 b)
  { return emsr::detail::ellint_cel<__float128>(k_c, p, a, b); }

  inline std::complex<__float128>
  cyl_hankel_1q(std::complex<__float128> nu, std::complex<__float128> x)
  { return emsr::detail::cyl_hankel_1<__float128>(nu, x); }

  inline std::complex<__float128>
  cyl_hankel_2q(std::complex<__float128> nu, std::complex<__float128> x)
  { return emsr::detail::cyl_hankel_2<__float128>(nu, x); }

  inline std::complex<__float128>
  sph_hankel_1q(unsigned int n, std::complex<__float128> x)
  { return emsr::detail::sph_hankel_1<__float128>(n, x); }

  inline std::complex<__float128>
  sph_hankel_2q(unsigned int n, std::complex<__float128> x)
  { return emsr::detail::sph_hankel_2<__float128>(n, x); }

  inline std::complex<__float128>
  sph_harmonicq(unsigned int l, unsigned int m,
		__float128 theta, __float128 phi)
  { return emsr::detail::sph_harmonic<__float128>(l, m, theta, phi); }

  inline __float128
  polylogq(__float128 s, __float128 w)
  { return emsr::detail::polylog<__float128>(s, w); }

  inline std::complex<__float128>
  polylogq(__float128 s, std::complex<__float128> w)
  { return emsr::detail::polylog<__float128>(s, w); }

  inline __float128
  dirichlet_etaq(__float128 s)
  { return emsr::detail::dirichlet_eta<__float128>(s); }

  inline __float128
  dirichlet_betaq(__float128 s)
  { return emsr::detail::dirichlet_beta<__float128>(s); }

  inline __float128
  dirichlet_lambdaq(__float128 s)
  { return emsr::detail::dirichlet_lambda<__float128>(s); }

  inline __float128
  clausen_slq(unsigned int m, __float128 w)
  { return emsr::detail::clausen_sl<__float128>(m, w); }

  inline __float128
  clausen_clq(unsigned int m, __float128 w)
  { return emsr::detail::clausen_cl<__float128>(m, w); }

  inline __float128
  clausenq(unsigned int m, __float128 w)
  { return emsr::detail::clausen<__float128>(m, w); }

  inline std::complex<__float128>
  clausenq(unsigned int m, std::complex<__float128> w)
  { return emsr::detail::clausen<__float128>(m, w); }

  inline __float128
  theta_1q(__float128 nu, __float128 x)
  { return emsr::detail::theta_1<__float128>(nu, x); }

  inline __float128
  theta_2q(__float128 nu, __float128 x)
  { return emsr::detail::theta_2<__float128>(nu, x); }

  inline __float128
  theta_3q(__float128 nu, __float128 x)
  { return emsr::detail::theta_3<__float128>(nu, x); }

  inline __float128
  theta_4q(__float128 nu, __float128 x)
  { return emsr::detail::theta_4<__float128>(nu, x); }

  inline __float128
  ellnomeq(__float128 k)
  { return emsr::detail::ellnome<__float128>(k); }

  inline __float128
  theta_sq(__float128 k, __float128 x)
  { return emsr::detail::theta_s<__float128>(k, x); }

  inline __float128
  theta_dq(__float128 k, __float128 x)
  { return emsr::detail::theta_d<__float128>(k, x); }

  inline __float128
  theta_cq(__float128 k, __float128 x)
  { return emsr::detail::theta_c<__float128>(k, x); }

  inline __float128
  theta_nq(__float128 k, __float128 x)
  { return emsr::detail::theta_n<__float128>(k, x); }

  inline __float128
  owens_tq(__float128 __h, __float128 a)
  { return emsr::detail::owens_t<__float128>(__h, a); }

  inline __float128
  fermi_diracq(__float128 s, __float128 x)
  { return emsr::detail::fermi_dirac<__float128>(s, x); }

  inline __float128
  bose_einsteinq(__float128 s, __float128 x)
  { return emsr::detail::bose_einstein<__float128>(s, x); }

  inline __float128
  sin_piq(__float128 x)
  { return emsr::detail::sin_pi<__float128>(x); }

  inline __float128
  sinh_piq(__float128 x)
  { return emsr::detail::sinh_pi<__float128>(x); }

  inline __float128
  cos_piq(__float128 x)
  { return emsr::detail::cos_pi<__float128>(x); }

  inline __float128
  cosh_piq(__float128 x)
  { return emsr::detail::cosh_pi<__float128>(x); }

  inline __float128
  tan_piq(__float128 x)
  { return emsr::detail::tan_pi<__float128>(x); }

  inline __float128
  tanh_piq(__float128 x)
  { return emsr::detail::tanh_pi<__float128>(x); }

  inline __sincos_t<__float128>
  sincosq(__float128 x)
  { return emsr::detail::sincos<__float128>(x); }

  inline __sincos_t<__float128>
  sincos_piq(__float128 x)
  { return emsr::detail::sincos_pi<__float128>(x); }

} // namespace emsr

#endif // EMSR_HAVE_FLOAT128

#endif // SPECFUN128_H
