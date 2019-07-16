unit SpecFunX;

{Special functions (extended precision)}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


uses
  sfBasic;   {Basic common code}

(*************************************************************************

 DESCRIPTION   :  Special functions (extended precision)

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  See amath_info.txt/references


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     15.01.10  W.Ehrhardt  Initial version: lngammax
 0.11     17.01.10  we          invgammax
 0.12     18.01.10  we          psix
 0.13     18.01.10  we          lnbetax
 0.14     18.01.10  we          betax
 0.15     20.01.10  we          Elliptic integrals (Carlson style)
 0.16     21.01.10  we          Elliptic integrals (Bulirsch style)
 0.17     23.01.10  we          Elliptic integrals (Maple style)
 0.18     23.01.10  we          sfc_signgammax returns extended
 0.19     23.01.10  we          agmx
 0.20     24.01.10  we          sncndnx
 0.21     24.01.10  we          erfx, erfcx
 0.22     24.01.10  we          normstd_pdfx, normstd_cdfx, normstd_invx
 0.23     27.01.10  we          gammax
 0.24     10.02.10  we          zetax, zeta1px
 0.25     14.02.10  we          e1x
 0.26     15.02.10  we          eix, lix
 0.27     18.02.10  we          enx

 0.28     02.03.10  we          chix, shix
 0.29     03.03.10  we          cix, six, ssix
 0.30     06.03.10  we          dawsonx
 0.31     11.03.10  we          erf_invx, erfc_invx
 0.32     16.03.10  we          dilogx

 0.33     06.04.10  we          bessel_j0x/j1x/jnx/y0x/y1x/ynx
 0.34     06.04.10  we          bessel_i0x/i0ex/i1x/i1ex
 0.35     06.04.10  we          bessel_k0x/k0ex/k1x/k1ex
 0.36     08.04.10  we          cinx, cinhx
 0.37     10.04.10  we          facx
 0.38     11.04.10  we          bessel_inx/knx
 0.39     12.04.10  we          ti2x
 0.40     17.04.10  we          gamma1pm1x
 0.41     22.04.10  we          cl2x

 0.42     19.05.10  we          Fresnelx
 0.43     20.05.10  we          LambertWx, LambertW1x

 0.44     11.07.10  we          ibetax, beta distribution functions
 0.45     13.07.10  we          t_pdfx, t_cdfx, t_invx
 0.46     16.07.10  we          f_pdfx, f_cdfx, f_invx

 1.00.00  17.08.10  we          Common version number after split
 1.00.01  18.08.10  we          gammastarx
 1.00.02  24.08.10  we          incgammax, igammapx, igammaqx
 1.00.03  05.09.10  we          inverse normalised incomplete gamma functions
 1.00.04  05.09.10  we          chi2_pdfx/cdfx/invx
 1.00.05  07.09.10  we          gamma_pdfx/cdfx/invx
 1.00.06  10.09.10  we          trigammax
 1.00.07  15.09.10  we          etax, etam1x, zetam1x

 1.01.00  06.10.10  we          dfacx
 1.01.01  08.10.10  we          legendre_px, legendre_qx, legendre_plmx
 1.01.02  10.10.10  we          chebyshev_tx, chebyshev_ux
 1.01.03  11.10.10  we          gegenbauer_cx
 1.01.04  15.10.10  we          jacobi_px
 1.01.05  16.10.10  we          hermite_hx
 1.01.06  17.10.10  we          laguerrex
 1.01.07  19.10.10  we          gamma_delta_ratiox, gamma_ratiox, pochhammerx
 1.01.08  20.10.10  we          laguerre_assx, laguerre_lx
 1.01.09  20.10.10  we          spherical_harmonicx
 1.01.10  22.10.10  we          zernike_rx
 1.01.11  23.10.10  we          zetahx

 1.02.00  01.11.10  we          bessel_jvx, bessel_yvx
 1.02.01  04.11.10  we          bessel_ivx, bessel_kvx
 1.02.02  05.11.10  we          Airy functions Ai, Ai', Bi, Bi'
 1.02.03  06.11.10  we          sph_bessel_jnx, sph_bessel_ynx
 1.02.04  15.11.10  we          bessel_ivex, bessel_kvex

 1.03.00  03.12.10  we          polylogx
 1.03.01  05.12.10  we          bernoullix, zetaintx
 1.03.02  07.12.10  we          debyex
 1.03.03  15.12.10  we          tetragammax, pentagammax, polygammax
 1.03.04  31.12.10  we          binomialx

 1.04.00  12.02.11  we          Goodwin-Staton integral gsix
 1.04.01  13.02.11  we          expint3x
 1.04.02  14.02.11  we          erfix
 1.04.03  17.02.11  we          RiemannRx

 1.05.00  11.04.11  we          Zero order Kelvin functions berx,beix,kerx,keix
 1.05.01  13.04.11  we          cauchy_pdfx/cdfx/invx
 1.05.02  14.04.11  we          normal_pdfx/cdfx/invx
 1.05.03  14.04.11  we          exp_pdfx/cdfx/invx
 1.05.04  15.04.11  we          lognormal_pdfx/cdfx/invx
 1.05.05  16.04.11  we          logistic_pdfx/cdfx/invx
 1.05.06  16.04.11  we          weibull_pdfx/cdfx/invx
 1.05.07  17.04.11  we          laplace_pdfx/cdfx/invx
 1.05.08  17.04.11  we          changed sh,sc to a,b in gamma_pdfx/cdfx/invx
 1.05.09  19.04.11  we          pareto_pdfx/cdfx/invx
 1.05.10  20.04.11  we          uniform_pdfx/cdfx/invx

 1.06.00  08.05.11  we          Struve functions struve_h0x/h1x/l0x,l1x
 1.06.01  13.05.11  we          theta2x,theta3x,theta4x
 1.06.02  14.05.11  we          theta1px
 1.06.03  17.05.11  we          jacobi_thetax
 1.06.04  19.05.11  we          EllipticModulusx
 1.06.05  19.05.11  we          replaced misleading kc by k in EllipticCKx and EllipticCEx
 1.06.06  20.05.11  we          EllipticNomex

 1.07.00  01.06.11  we          einx
 1.07.01  04.06.11  we          ellint_1x, ellint_2x
 1.07.02  05.06.11  we          comp_ellint_1x/2x/3x
 1.07.03  06.06.11  we          jacobi_amx
 1.07.04  07.06.11  we          jacobi_zetax
 1.07.05  08.06.11  we          ellint_3x
 1.07.06  10.06.11  we          erfcex
 1.07.07  17.06.11  we          jacobi_arcsnx/cnx/dnx
 1.07.08  18.06.11  we          jacobi_snx/cnx/dnx
 1.07.09  23.06.11  we          heuman_lambdax

 1.08.00  03.08.11  we          invgamma renamed to rgamma

 1.09.00  05.10.11  we          lngamma1px

 1.10.00  10.12.11  we          triangular_pdfx/cdfx/invx
 1.10.01  14.12.11  we          binomial_cdfx/pmfx
 1.10.02  15.12.11  we          poisson_cdfx/pmfx
 1.10.03  16.12.11  we          lnfacx
 1.10.04  18.12.11  we          negbinom_cdfx/pmfx
 1.10.05  21.12.11  we          primezetax
 1.10.06  27.12.11  we          hypergeo_cdfx/pmfx
 1.10.07  28.12.11  we          sph_bessel_inx, sph_bessel_inex
 1.10.08  29.12.11  we          sph_bessel_knx, sph_bessel_knex

 1.11.00  05.02.12  we          LerchPhix
 1.11.01  06.02.12  we          DirichletBetax
 1.11.02  06.02.12  we          DirichletLambdax
 1.11.03  07.02.12  we          LegendreChix
 1.11.04  08.02.12  we          polylogrx

 1.12.00  22.03.12  we          poch1x
 1.12.01  11.04.12  we          igammax (non-normalised upper incomplete gamma)
 1.12.02  22.04.12  we          igammalx (non-normalised lower incomplete gamma)
 1.12.03  23.04.12  we          lngammasx
 1.12.04  27.04.12  we          legendre_qlmx
 1.12.05  16.05.12  we          toroidal_qlmx
 1.12.06  20.05.12  we          toroidal_plmx

 1.13.00  10.06.12  we          fermi_diracx
 1.13.01  13.06.12  we          fermi_dirac_m05x/p05x
 1.13.02  15.06.12  we          dawson2x
 1.13.03  16.06.12  we          inerfcx
 1.13.04  23.06.12  we          geix
 1.13.05  26.06.12  we          erfgx
 1.13.06  30.06.12  we          igammatx

 1.16.00  14.03.13  we          Remaining Jacobi elliptic functions
 1.16.01  25.03.13  we          Remaining inverse Jacobi elliptic functions
 1.16.02  28.03.13  we          trilogx

 1.17.00  11.04.13  we          fermi_dirac_p15x
 1.17.01  13.04.13  we          k changed to longint in poisson_pmfx
 1.17.02  14.04.13  we          rayleigh_pdfx/cdfx/invx
 1.17.03  17.04.13  we          beta3x
 1.17.04  19.04.13  we          maxwell_pdfx/cdfx/invx
 1.17.05  20.04.13  we          evt1_pdfx/cdfx/invx
 1.17.06  01.05.13  we          FresnelCx/FresnelSx

 1.18.00  09.05.13  we          Airy/Scorer functions Gi/Hi
 1.18.01  19.05.13  we          hyperg_2F1x
 1.18.02  22.05.13  we          hyperg_2F1rx
 1.18.03  24.05.13  we          hyperg_1F1x, hyperg_1F1rx
 1.18.03  30.05.13  we          hyperg_ux

 1.19.00  16.06.13  we          Whittaker functions
 1.19.01  27.06.13  we          hyperg_0F1x
 1.19.02  01.07.13  we          hyperg_0F1rx

 1.20.00  15.08.13  we          kumaraswamy_pdfx/cdfx/invx
 1.20.01  18.08.13  we          erf_px, erf_qx, erf_zx

 1.21.00  11.09.13  we          Bernoulli polynomials bernpolyx
 1.21.01  17.09.13  we          chebyshev_vx/wx
 1.21.02  24.09.13  we          cosintx, sinintx
 1.21.03  25.09.13  we          BatemanGx

 1.22.00  19.10.13  we          moyal_pdfx/cdfx/invx
 1.22.01  23.10.13  we          EllipticKimx,EllipticECimx

 1.23.00  26.12.13  we          fermi_dirac_p25x

 1.25.00  01.05.14  we          harmonicx
 1.25.01  03.05.14  we          harmonic2x
 1.25.02  04.05.14  we          zipf_pmfx/cdfx

 1.26.00  24.05.14  we          li_invx
 1.26.01  30.05.14  we          lobachevsky_cx/sx
 1.26.02  02.06.14  we          fibpolyx, lucpolyx
 1.26.03  05.06.14  we          levy_pdfx/cdfx/invx

 1.27.00  22.06.14  we          lngamma_invx
 1.27.01  25.06.14  we          psi_invx
 1.27.02  02.07.14  we          ibeta_invx
 1.27.03  12.07.14  we          CylinderDx, CylinderUx, HermiteHx

 1.28.00  17.08.14  we          catalanx
 1.28.01  19.08.14  we          CylinderVx

 1.31.00  23.12.14  we          invgamma_pdfx/cdfx/invx
 1.31.01  26.12.14  we          logseries_cdfx/pmfx: logarithmic (series) distribution
 1.31.02  03.01.15  we          wald_cdfx/pdfx/invx
 1.31.03  09.01.15  we          Euler numbers

 1.32.00  25.01.15  we          ei_invx
 1.32.01  28.03.15  we          comp_ellint_dx, ellint_dx

 1.33.00  04.06.15  we          keplerx
 1.33.01  10.06.15  we          struve_hx, struve_lx
 1.33.02  19.06.15  we          airy_aisx, airy_bisx
 1.33.03  22.06.15  we          ell_rgx

 1.36.00  26.09.15  we          lemniscate functions sin_lemnx, cos_lemnx

 1.38.00  12.05.16  we          comp_ellint_bx, ellint_bx
 1.38.01  18.06.16  we          Bring radical bringx
 1.38.02  25.06.16  we          Wright omegax

 1.39.00  23.08.16  we          inv_gammax
 1.39.01  28.09.16  we          fibfunx
 1.39.02  30.09.16  we          chebyshev_f1x
 1.39.03  07.10.16  we          hermite_hex
 1.39.04  08.10.16  we          kolmogorov_cdfx/invx

 1.40.00  25.05.17  we          euler_qx
 1.40.01  21.06.17  we          detaix

 1.41.00  08.07.17  we          Fresnel(FG)x
 1.41.01  16.07.17  we          Bessel lambda
 1.41.02  19.07.17  we          erfh/2x
 1.41.03  22.07.17  we          chi_pdfx/cdfx/invx

 1.42.00  17.08.17  we          eisx2x

 1.43.00  01.12.17  we          erfi_invx
 1.43.01  06.12.17  we          nakagami_pdfx/cdfx/invx
 1.43.02  08.12.17  we          MarcumQx
 1.43.03  13.12.17  we          E1sx(x) = exp(x)*E1(x)
 1.43.04  14.12.17  we          Eisx(x) = exp(-x)*Ei(x)

 1.44.00  10.01.18  we          RiemannR_invx
 1.44.01  18.01.18  we          Owen's T function OwenTx
 1.44.02  23.01.18  we          erfce_invx
 1.44.03  30.01.18  we          einsteinx
 1.44.04  02.02.18  we          transportx
 1.44.05  09.02.18  we          fermi_dirac_rx

 1.45.00  20.02.18  we          LangevinLx
 1.45.01  21.02.18  we          LangevinL_invx
 1.45.02  26.02.18  we          tix (inverse tangent integral of real order)
 1.45.03  05.03.18  we          Bose-Einstein integral
 1.45.04  11.03.18  we          Rogers-Ramanujan continued fraction rrcfx

 1.46.00  26.03.18  we          bessel_[i,j,k,y]0_intx
 1.46.01  02.04.18  we          EllipticPiCimx
 1.46.02  05.04.18  we          Mathematica style complete elliptic integrals
 1.46.03  12.04.18  we          Neville theta functions
 1.46.04  15.04.18  we          arcclx/slx

 1.47.00  23.04.18  we          emlambdax
 1.47.01  27.04.18  we          Weierstrass elliptic functions
 1.47.02  03.05.18  we          wpg_invx
 1.47.03  04.05.18  we          M_EllipticFx, M_EllipticEx
 1.47.04  06.05.18  we          KleinJx

 1.48.00  14.05.18  we          wpe_derx / wpg_derx
 1.48.01  23.05.18  we          M_EllipticPix
 1.48.02  07.06.18  we          lnbinomialx

 1.49.00  29.07.18  we          exprelnx
 1.49.01  03.08.18  we          hyperg_2F0x

 1.50.00  21.08.18  we          eulerpolyx
 1.50.01  26.08.18  we          besselpolyx
 1.50.02  28.08.18  we          igammap_derx
 1.50.03  01.09.18  we          lnBarnesGx
 1.50.04  06.09.18  we          igammanx(n,x)
 1.50.05  07.09.18  we          eibetax

 1.51.00  16.09.18  we          CoulombCLx
 1.51.01  17.09.18  we          CoulombSLx
 1.51.02  28.09.18  we          CoulombFFPx/Fx/GGpx
 1.51.03  01.10.18  we          expnx
 1.51.04  02.10.18  we          removed igammanx

 1.52.00  12.10.18  we          psistarx
 1.52.01  17.10.18  we          SynchFx, SynchGx

 1.53.00  20.11.18  we          Derivatives of Kelvin functions

***************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2009-2018 Wolfgang Ehrhardt

 This software is provided 'as-is', without any express or implied warranty.
 In no event will the authors be held liable for any damages arising from
 the use of this software.

 Permission is granted to anyone to use this software for any purpose,
 including commercial applications, and to alter it and redistribute it
 freely, subject to the following restrictions:

 1. The origin of this software must not be misrepresented; you must not
    claim that you wrote the original software. If you use this software in
    a product, an acknowledgment in the product documentation would be
    appreciated but is not required.

 2. Altered source versions must be plainly marked as such, and must not be
    misrepresented as being the original software.

 3. This notice may not be removed or altered from any source distribution.
----------------------------------------------------------------------------*)

const
  SpecFunX_Version = SF_BaseVersion+'-X';

{#Z+}
{---------------------------------------------------------------------------}
{--------------------------- Bessel functions ------------------------------}
{---------------------------------------------------------------------------}
{#Z-}
function airy_aix(x: extended): extended;
  {-Return the Airy function Ai(x)}

function airy_aipx(x: extended): extended;
  {-Return the Airy function Ai'(x)}

function airy_aisx(x: extended): extended;
  {-Return the scaled Airy function Ai(x) if x <= 0, Ai(x)*exp(2/3*x^1.5) for x > 0}

function airy_bix(x: extended): extended;
  {-Return the Airy function Bi(x)}

function airy_bipx(x: extended): extended;
  {-Return the Airy function Bi'(x)}

function airy_bisx(x: extended): extended;
  {-Return the scaled Airy function Bi(x) if x <= 0, Bi(x)*exp(-2/3*x^1.5) for x > 0}

function airy_gix(x: extended): extended;
  {-Return the Airy/Scorer function Gi(x) = 1/Pi*integral(sin(x*t+t^3/3), t=0..INF)}

function airy_hix(x: extended): extended;
  {-Return the Airy/Scorer function Hi(x) = 1/Pi*integral(exp(x*t-t^3/3), t=0..INF)}


function bessel_i0x(x: extended): extended;
  {-Return I0(x), the modified Bessel function of the 1st kind, order zero}

function bessel_i0ex(x: extended): extended;
  {-Return I0(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order zero}

function bessel_i1x(x: extended): extended;
  {-Return I1(x), the modified Bessel function of the 1st kind, order one}

function bessel_i1ex(x: extended): extended;
  {-Return I1(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order one}

function bessel_inx(n: integer; x: extended): extended;
  {-Return I_n(x), the modified Bessel function of the 1st kind, order n.}

function bessel_ivx(v, x: extended): extended;
  {-Return I_v(x), the modified Bessel function of the 1st kind, order v.}

function bessel_ivex(v, x: extended): extended;
  {-Return I_v(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order v.}

function bessel_j0x(x: extended): extended;
  {-Return J0(x), the Bessel function of the 1st kind, order zero}

function bessel_j1x(x: extended): extended;
  {-Return J1(x), the Bessel function of the 1st kind, order one}

function bessel_jnx(n: integer; x: extended): extended;
  {-Return J_n(x), the Bessel function of the 1st kind, order n; not suitable for large n or x.}

function bessel_jvx(v, x: extended): extended;
  {-Return J_v(x), the Bessel function of the 1st kind, order v; not suitable for large v.}

function bessel_k0x(x: extended): extended;
  {-Return K0(x), the modified Bessel function of the 2nd kind, order zero, x>0}

function bessel_k0ex(x: extended): extended;
  {-Return K0(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order zero, x>0}

function bessel_k1x(x: extended): extended;
  {-Return K1(x), the modified Bessel function of the 2nd kind, order one, x>0}

function bessel_k1ex(x: extended): extended;
  {-Return K1(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order one, x>0}

function bessel_knx(n: integer; x: extended): extended;
  {-Return K_n(x), the modified Bessel function of the 2nd kind, order n, x>0, not suitable for large n}

function bessel_kvx(v, x: extended): extended;
  {-Return K_v(x), the modified Bessel function of the 2nd kind, order v, x>0}

function bessel_kvex(v, x: extended): extended;
  {-Return K_v(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order v, x>0}

function bessel_y0x(x: extended): extended;
  {-Return Y0(x), the Bessel function of the 2nd kind, order zero; x>0}

function bessel_y1x(x: extended): extended;
  {-Return Y1(x), the Bessel function of the 2nd kind, order one; x>0}

function bessel_ynx(n: integer; x: extended): extended;
  {-Return Y_n(x), the Bessel function of the 2nd kind, order n, x>0, not suitable for large n or x}

function bessel_yvx(v, x: extended): extended;
  {-Return Y_v(x), the Bessel function of the 2nd kind, order v; x > 0; not suitable for large v.}

function bessel_lambdax(v,x: extended): extended;
  {-Compute lambda(v,x) = Gamma(v+1)*J(v,x)/(0.5x)^v for v,x >= 0}


function bessel_i0_intx(u: extended): extended;
  {-Return the integral int(bessel_i0(x), x = 0..u)}

function bessel_j0_intx(u: extended): extended;
  {-Return the integral int(bessel_j0(x), x = 0..u)}

function bessel_k0_intx(u: extended): extended;
  {-Return the integral int(bessel_k0(x), x = 0..u), u >= 0}

function bessel_y0_intx(u: extended): extended;
  {-Return the integral int(bessel_y0(x), x = 0..u), u >= 0}


function kelvin_berx(x: extended): extended;
  {-Return the Kelvin function ber(x)}

function kelvin_beix(x: extended): extended;
  {-Return the Kelvin function bei(x)}

function kelvin_kerx(x: extended): extended;
  {-Return the Kelvin function ker(x), x > 0}

function kelvin_keix(x: extended): extended;
  {-Return the Kelvin function kei(x), x >= 0}

procedure kelvin_kerkeix(x: extended; var kr, ki: extended);
  {-Return the Kelvin functions kr=ker(x), ki=kei(x), x > 0}

procedure kelvin_berbeix(x: extended; var br, bi: extended);
  {-Return the Kelvin functions br=ber(x), bi=bei(x)}

procedure kelvin_derx(x: extended; var berp, beip, kerp, keip: extended);
  {-Return the derivatives of the zero order Kelvin functions, x >= 0}

function kelvin_berpx(x: extended): extended;
  {-Return the Kelvin function ber'(x), x >= 0}

function kelvin_beipx(x: extended): extended;
  {-Return the Kelvin function bei'(x), x >= 0}

function kelvin_kerpx(x: extended): extended;
  {-Return the Kelvin function ker'(x), x > 0}

function kelvin_keipx(x: extended): extended;
  {-Return the Kelvin function kei'(x), x >= 0}


function sph_bessel_jnx(n: integer; x: extended): extended;
  {-Return j_n(x), the spherical Bessel function of the 1st kind, order n}

function sph_bessel_ynx(n: integer; x: extended): extended;
  {-Return y_n(x), the spherical Bessel function of the 2nd kind, order n >=0 , x<>0}

function sph_bessel_inx(n: integer; x: extended): extended;
  {-Return i_n(x), the modified spherical Bessel function of the 1st/2nd kind, order n}

function sph_bessel_inex(n: integer; x: extended): extended;
  {-Return i_n(x)*exp(-|x|), the exponentially scaled modified spherical Bessel function of the 1st/2nd kind, order n}

function sph_bessel_knx(n: integer; x: extended): extended;
  {-Return k_n(x), the modified spherical Bessel function of the 3rd kind, order n, x>0}

function sph_bessel_knex(n: integer; x: extended): extended;
  {-Return k_n(x)*exp(x), the exponentially scaled modified spherical Bessel function of the 3rd kind, order n, x>0}


function struve_h0x(x: extended): extended;
 {-Return H0(x), the Struve function of order 0}

function struve_h1x(x: extended): extended;
  {-Return H1(x), the Struve function of order 1}

function struve_l0x(x: extended): extended;
  {-Return L0(x), the modified Struve function of order 0}

function struve_l1x(x: extended): extended;
  {-Return L1(x), the modified Struve function of order 1}

function struve_hx(v, x: extended): extended;
  {-Return H_v(x), the Struve function of order v, x < 0 only if v is an integer.}

function struve_lx(v, x: extended): extended;
  {-Return L_v(x), the modified Struve function of order v, x < 0 only if v is an integer.}


function CoulombCLx(L: integer; eta: extended): extended;
  {-Return the normalizing constant CL for Coulomb wave function, L >= 0}

function CoulombSLx(L: integer; eta: extended): extended;
  {-Return the Coulomb phase shift sigma_L(eta) for L >= 0}

function CoulombFx(L: integer; eta, x: extended): extended;
  {-Return the regular Coulomb wave functions FL(eta,x) for L >= 0, x > 0}

procedure CoulombFFpx(L: integer; eta, x: extended; var fc,fcp: extended; var ifail: integer);
  {-Return the regular Coulomb wave functions fc=FL(eta,x) and fcp=FL'(eta,x)}
  { for L >= 0, x > 0. No error if ifail=0, failure if ifail < 0, if ifail=1 }
  { then the very inaccurate JWKB approximation has been computed.}

procedure CoulombGGpx(L: integer; eta, x: extended; var gc,gcp: extended; var ifail: integer);
  {-Return the irregular Coulomb wave functions gc=GL(eta,x) and gcp=GL'(eta,x)}
  { for  L >= 0, x > 0. No error if ifail=0, failure if ifail < 0, if ifail=1  }
  { then the very inaccurate JWKB approximation has been computed.}


function SynchFx(x: extended): extended;
  {-Return the first synchrotron function F(x) = integral(x*BesselK(5/3,t), t=x..INF) for x >= 0}

function SynchGx(x: extended): extended;
  {-Return the second synchrotron function G(x) = x*BesselK(2/3,x) for x >= 0}


{#Z+}
{---------------------------------------------------------------------------}
{---------------------- Gamma function and related -------------------------}
{---------------------------------------------------------------------------}
{#Z-}
const
  MAXGAMX  = 1755.455;                  {max. argument for gammax}
  MAXLGMX  = 1.0484814683901952E4928;   {max. argument for lngammax}
  MAXDFACX = 3209;                      {max. argument for dfacx}

function gammax(x: extended): extended;
  {-Return gamma(x), x <= MAXGAMX; invalid if x is a non-positive integer}

function gamma1pm1x(x: extended): extended;
  {-Return gamma(1+x)-1 with increased accuracy for x near 0}

function gammastarx(x: extended): extended;
  {-Return Temme's gammastar(x) = gamma(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x)), x>0.}
  { For large x the asymptotic expansion is gammastar(x) = 1 + 1/12x + O(1/x^2)}

function gamma_delta_ratiox(x,d: extended): extended;
  {-Return gamma(x)/gamma(x+d), accurate even for |d| << |x|}

function gamma_ratiox(x,y: extended): extended;
  {-Return gamma(x)/gamma(y)}

function pochhammerx(a,x: extended): extended;
  {-Return the Pochhammer symbol gamma(a+x)/gamma(a)}

function poch1x(a,x: extended): extended;
  {-Return (pochhammer(a,x)-1)/x, psi(a) if x=0; accurate even for small |x|}

function binomialx(n,k: integer): extended;
  {-Return the binomial coefficient 'n choose k'}

function facx(n: integer): extended;
  {-Return the factorial n!, n<MAXGAM-1; INF if n<0}

function dfacx(n: integer): extended;
  {-Return the double factorial n!!, n<=MAXDFACX; INF for even n<0}

function lnbinomialx(n,k: longint): extended;
  {-Return ln(binomial(n,k)), n >= k >= 0}

function lnfacx(n: longint): extended;
  {-Return ln(n!), INF if n<0}

function lngammax(x: extended): extended;
  {-Return ln(|gammax(x)|), |x| <= MAXLGM, invalid if x is a non-positive integer.}
  { Function signgammax can be used if the sign of gammax(x) is needed.}

function lngammasx(x: extended; var s: integer): extended;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, s=-1,1 is the sign of gamma}

function lngamma1px(x: extended): extended;
  {-Return ln(|gamma(1+x)|) with increased accuracy for x near 0}

function lngamma_invx(y: extended): extended;
  {-Inverse of lngamma: return x with lngamma(x) = y, y >= -0.12142, x > 1.4616}

function inv_gammax(y: extended): extended;
  {-Inverse of gamma: return x with gamma(x) = y, y >= 0.8857421875}

function signgammax(x: extended): extended;
  {-Return sign(gamma(x)), useless for 0 or negative integer}

function rgammax(x: extended): extended;
  {-Return the reciprocal gamma function rgamma = 1/gamma(x)}

procedure incgammax(a,x: extended; var p,q: extended);
  {-Return the normalised incomplete gamma functions P and Q, a>=0, x>=0}
  { P(a,x) = integral(exp(-t)*t^(a-1), t=0..x  )/gamma(a)}
  { Q(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf)/gamma(a)}

function igammapx(a,x: extended): extended;
  {-Return the normalised lower incomplete gamma function P(a,x), a>=0, x>=0}
  { P(a,x) = integral(exp(-t)*t^(a-1), t=0..x)/gamma(a)}

function igammaqx(a,x: extended): extended;
  {-Return the normalised upper incomplete gamma function Q(a,x), a>=0, x>=0}
  { Q(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf)/gamma(a)}

function igammax(a,x: extended): extended;
  {-Return the non-normalised upper incomplete gamma function}
  { GAMMA(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf). If x<0 }
  { the real part is returned and a must be <> -1, -2, ...   }

function igammalx(a,x: extended): extended;
  {-Return the non-normalised lower incomplete gamma function}
  { gamma(a,x) = integral(exp(-t)*t^(a-1), t=0..x); a<>0,-1,-2,..}

function igammatx(a,x: extended): extended;
  {-Return Tricomi's entire incomplete gamma function igammatx(a,x)}
  { = igammal(a,x)/gamma(a)/x^a = P(a,x)/x^a }

procedure incgamma_invx(a,p,q: extended; var x: extended; var ierr: integer);
  {-Return the inverse normalised incomplete gamma function, i.e. calculate}
  { x with P(a,x)=p and Q(a,x)=q. Input parameter a>0, p>=0, q>0 and p+q=1.}
  { ierr is >= 0 for success, < 0 for input errors or iterations failures. }

function igamma_invx(a,p,q: extended): extended;
  {-Return the inverse normalised incomplete gamma function, i.e. calculate}
  { x with P(a,x)=p and Q(a,x)=q. Input parameter a>0, p>=0, q>0 and p+q=1.}

function igammap_invx(a,p: extended): extended;
  {-Inverse incomplete gamma: return x with P(a,x)=p, a>=0, 0<=p<1}

function igammaq_invx(a,q: extended): extended;
  {-Inverse complemented incomplete gamma: return x with Q(a,x)=q, a>=0, 0<q<=1}

function igammap_derx(a,x: extended): extended;
  {-Return the partial derivative with respect to x of the normalised}
  { lower incomplete gamma function P(a,x), x >= 0, a <> 0,-1,-2 ...}

function psix(x: extended): extended;
  {-Return the psi (digamma) function of x, INF if x is a non-positive integer}

function psistarx(x: extended): extended;
  {-Return psi(x) - ln(x), x > 0}

function psi_invx(y: extended): extended;
  {-Inverse of psi, return x with psi(x)=y, y <= ln_MaxExt}

function trigammax(x: extended): extended;
  {-Return the trigamma function of x, INF if x is a negative integer}

function tetragammax(x: extended): extended;
  {-Return the tetragamma function psi''(x), NAN/RTE if x is a negative integer}

function pentagammax(x: extended): extended;
  {-Return the pentagamma function psi'''(x), INF if x is a negative integer}

function polygammax(n: integer; x: extended): extended;
  {-Return the polygamma function: n'th derivative of psi; n>=0, x>0 for n>MAXGAMX}
  { Note: The accuracy may be reduced for n >= MAXGAMX due to ln/exp operations.}

function lnBarnesGx(x: extended): extended;
  {-Return ln(BarnesG(x)), real part for x < 0}

function BatemanGx(x: extended): extended;
  {-Return the Bateman function G(x) = psi((x+1)/2) - psi(x/2), x<>0,-1,-2,...}

function lnbetax(x,y: extended): extended;
  {-Return the logarithm of |beta(x,y)|=|gamma(x)*gamma(y)/gamma(x+y)|}

function betax(x,y: extended): extended;
  {-Return the function beta(x,y)=gamma(x)*gamma(y)/gamma(x+y)}

function ibetax(a, b, x: extended): extended;
  {-Return the normalised incomplete beta function, a>0, b>0, 0 <= x <= 1}
  { ibetax = integral(t^(a-1)*(1-t)^(b-1) / betax(a,b), t=0..x)}

function ibeta_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the normalised incomplete beta function}
  { with a > 0, b > 0, and 0 <= y <= 1.}

function beta3x(a, b, x: extended): extended;
  {-Return the non-normalised incomplete beta function B_x(a,b)}
  { for 0<=x<=1, B_x = integral(t^(a-1)*(1-t)^(b-1), t=0..x).  }

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Zeta functions and polylogarithms   -------------------}
{---------------------------------------------------------------------------}
{#Z-}
function zetax(s: extended): extended;
  {-Return the Riemann zeta function at s, s<>1}

function zeta1px(x: extended): extended;
  {-Return the Riemann zeta function at 1+x, x<>0}

function zetaintx(n: integer): extended;
  {-Return zeta(n) for integer arguments, n<>1}

function zetam1x(s: extended): extended;
  {-Return Riemann zeta(s)-1, s<>1}

function zetahx(s,a: extended): extended;
  {-Return the Hurwitz zeta function zetah(s,a) = sum(1/(i+a)^s, i=0..INF), s<>1, a>0}

function LerchPhix(z,s,a: extended): extended;
  {-Return the Lerch transcendent Phi(z,s,a) = sum(z^n/(n+a)^s, n=0..INF), z <= 1, s >= -1, a >= 0}

function DirichletBetax(s: extended): extended;
  {-Return the Dirichlet beta function sum((-1)^n/(2n+1)^s, n=0..INF)}

function DirichletLambdax(s: extended): extended;
  {-Return the Dirichlet lambda function sum(1/(2n+1)^s, n=0..INF), s<>1}

function LegendreChix(s, x: extended): extended;
  {-Return Legendre's Chi-function chi(s,x); s>=0, |x|<=1, x<>1 if s<=1}

function etax(s: extended): extended;
  {-Return the Dirichlet eta function}

function etaintx(n: integer): extended;
  {-Return the Dirichlet function eta(n) for integer arguments}

function etam1x(s: extended): extended;
  {-Return Dirichlet eta(s)-1}

function primezetax(x: extended): extended;
  {-Return the prime zeta function P(x) = sum(1/p^x, p prime), x > 1/5; }
  { for x<1 the real part of P(x) is returned.}

function harmonicx(x: extended): extended;
  {-Return the harmonic number function H(x) = psi(x+1) + EulerGamma}

function harmonic2x(x,r: extended): extended;
  {-Return the generalized harmonic function H(x,r) = zeta(r)-zetah(r,x+1); x >= -1}

function cl2x(x: extended): extended;
  {-Return the Clausen function: integral(-ln(2*|sin(t/2)|),t=0..x) = Im(Li_2(exp(ix)))}

function ti2x(x: extended): extended;
  {-Return the inverse tangent integral, ti2(x) = integral(arctan(t)/t, t=0..x)}

function tix(s,x: extended): extended;
  {-Return the inverse tangent integral of order s >= 0}

function lobachevsky_cx(x: extended): extended;
  {-Return the Lobachevski function L(x) = integral(-ln(|cos(t)|), t=0..x)}

function lobachevsky_sx(x: extended): extended;
  {-Return the Lobachevski function Lambda(x) = integral(-ln(|2sin(t)|), t=0..x)}

function dilogx(x: extended): extended;
  {-Return dilog(x) = Re(Li_2(x)), Li_2(x) = -integral(ln(1-t)/t, t=0..x)}

function trilogx(x: extended): extended;
  {-Return the trilogarithm function trilog(x) = Re(Li_3(x))}

function polylogx(n: integer; x: extended): extended;
  {-Return the polylogarithm Li_n(x) of integer order; real part for n>0,x>1}

function polylogrx(s, x: extended): extended;
  {-Return the polylogarithm Li_s(x) of real order s >= -1, x <= 256;}
  { s > 0 if x > 256; real part if x > 1}

function bose_einsteinx(s,x: extended): extended;
  {-Return the Bose-Einstein integral of real order s >= -1}

function fermi_diracx(n: integer; x: extended): extended;
  {-Return the integer order Fermi-Dirac integral F_n(x) = 1/n!*integral(t^n/(exp(t-x)+1), t=0..INF)}

function fermi_dirac_rx(s,x: extended): extended;
  {-Return the Fermi-Dirac integral of real order s >= -1}

{#Z+}
{Obsolete functions, use fermi_dirac_rx with s = -1/2, 1/2, 3/2, 5/2 }
{#Z-}

function fermi_dirac_m05x(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}

function fermi_dirac_p05x(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}

function fermi_dirac_p15x(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}

function fermi_dirac_p25x(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Legendre style) -------------------}
{---------------------------------------------------------------------------}
{#Z-}
{---------------------------------------------------------------------------}
function comp_ellint_1x(k: extended): extended;
  {-Return the complete elliptic integral of the 1st kind, real part if |k|>1}

function comp_ellint_2x(k: extended): extended;
  {-Return the complete elliptic integral of the 2nd kind, real part if |k|>1}

function comp_ellint_3x(nu,k: extended): extended;
  {-Return the complete elliptic integral of the 3rd kind, |k|<>1, nu<>1; real part for |k|>1}

function comp_ellint_bx(k: extended): extended;
  {-Return the complete elliptic integral B(k) = (E(k) - kc^2*K(k))/k^2, real part for |k|>1}

function comp_ellint_dx(k: extended): extended;
  {-Return the complete elliptic integral D(k) = (K(k) - E(k))/k^2, real part for |k|>1}

function ellint_1x(phi,k: extended): extended;
  {-Return the Legendre elliptic integral F(phi,k) of the 1st kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}

function ellint_2x(phi,k: extended): extended;
  {-Return the Legendre elliptic integral E(phi,k) of the 2nd kind}
  { = integral(sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}

function ellint_3x(phi,nu,k: extended): extended;
  {-Return the Legendre elliptic integral PI(phi,nu,k) of the 3rd kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2)/(1-nu*sin(x)^2),x=0..phi) with  }
  { |k*sin(phi)|<=1, returns Cauchy principal value if nu*sin(phi)^2>1}

function ellint_bx(phi,k: extended): extended;
  {-Return the Legendre elliptic integral B(phi,k) = (E(phi,k) - kc^2*F(phi,k))/k^2}
  { = integral(cos(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }

function ellint_dx(phi,k: extended): extended;
  {-Return the Legendre elliptic integral D(phi,k) = (F(phi,k) - E(phi,k))/k^2 }
  { = integral(sin(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }

function heuman_lambdax(phi,k: extended): extended;
  {-Return Heuman's function Lambda_0(phi,k) = F(phi,k')/K(k') + 2/Pi*K(k)*Z(phi,k'), |k|<=1}

function jacobi_zetax(phi,k: extended): extended;
  {-Return the Jacobi Zeta function Z(phi,k) = E(phi,k) - E(k)/K(k)*F(phi,k), |k|<=1}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Carlson style) --------------------}
{---------------------------------------------------------------------------}
{#Z-}
function ell_rcx(x,y: extended): extended;
  {-Return Carlson's degenerate elliptic integral RC; x>=0, y<>0}

function ell_rfx(x,y,z: extended): extended;
  {-Return Carlson's elliptic integral of the 1st kind; x,y,z >=0, at most one =0}

function ell_rdx(x,y,z: extended): extended;
  {-Return Carlson's elliptic integral of the 2nd kind; z>0; x,y >=0, at most one =0}

function ell_rgx(x,y,z: extended): extended;
  {-Return Carlson's completely symmetric elliptic integral of the 2nd kind; x,y,z >= 0}

function ell_rjx(x,y,z,r: extended): extended;
  {-Return Carlson's elliptic integral of the 3rd kind; r<>0; x,y,z >=0, at most one =0}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Bulirsch style) -------------------}
{---------------------------------------------------------------------------}
{#Z-}
function cel1x(kc: extended): extended;
  {-Return Bulirsch's complete elliptic integral of the 1st kind, kc<>0}

function cel2x(kc, a, b: extended): extended;
  {-Return Bulirsch's complete elliptic integral of the 2nd kind, kc<>0}

function celx(kc, p, a, b: extended): extended;
  {-Return Bulirsch's general complete elliptic integral, kc<>0, Cauchy principle value if p<0}

function el1x(x,kc: extended): extended;
  {-Return Bulirsch's incomplete elliptic integral of the 1st kind}

function el2x(x,kc,a,b: extended): extended;
  {-Return Bulirsch's incomplete elliptic integral of the 2nd kind, kc<>0}

function el3x(x,kc,p: extended): extended;
  {-Return Bulirsch's incomplete elliptic integral of the 3rd kind, 1+p*x^2<>0}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Maple V style) --------------------}
{---------------------------------------------------------------------------}
{#Z-}
function EllipticFx(z,k: extended): extended;
  {-Return the incomplete elliptic integral of the 1st kind; |z|<=1, |k*z|<1}

function EllipticKx(k: extended): extended;
  {-Return the complete elliptic integral of the 1st kind, |k| < 1}

function EllipticKimx(k: extended): extended;
  {-Return K(i*k), the complete elliptic integral of the 1st kind with}
  { imaginary modulus = integral(1/sqrt(1-x^2)/sqrt(1+k^2*x^2),x=0..1)}

function EllipticCKx(k: extended): extended;
  {-Return the complementary complete elliptic integral of the 1st kind, k<>0}

function EllipticEx(z,k: extended): extended;
  {-Return the incomplete elliptic integrals of the 2nd kind, |z|<=1, |k*z| <= 1}

function EllipticECx(k: extended): extended;
  {-Return the complete elliptic integral of the 2nd kind, real part if |k|>1}

function EllipticECimx(k: extended): extended;
  {-Return E(i*k), the complete elliptic integral of the 2nd kind with}
  { imaginary modulus = integral(sqrt(1+k^2*x^2)/sqrt(1-x^2),x=0..1)  }

function EllipticCEx(k: extended): extended;
  {-Return the complementary complete elliptic integral of the 2nd kind}

function EllipticPix(z,nu,k: extended): extended;
  {-Return the incomplete elliptic integral of the 3rd kind, |z|<=1, |k*z|<1}

function EllipticPiCx(nu,k: extended): extended;
  {-Return the complete elliptic integral of the 3rd kind, |k|<>1, nu<>1; real part for |k|>1}

function EllipticPiCimx(nu,k: extended): extended;
  {-Return Pi(nu, k*i), the complete elliptic integral of the 3rd kind}
  { with imaginary modulus, nu <> 1, real part if nu > 1}

function EllipticCPix(nu,k: extended): extended;
  {-Return the complementary complete elliptic integral of the 3rd kind, k<>0, nu<>1}

{#Z+}
{---------------------------------------------------------------------------}
{---------------- Elliptic integrals (Mathematica style) -------------------}
{---------------------------------------------------------------------------}
{#Z-}
{---------------------------------------------------------------------------}
function M_EllipticKx(m: extended): extended;
  {-Return the complete elliptic integral of the 1st kind,}
  { K(m) = integral(dx/sqrt(1-m*sin(x)^2),x=0..Pi/2), real part for m>1}

function M_EllipticECx(m: extended): extended;
  {-Return the complete elliptic integral of the 2nd kind,}
  { E(m) = integral(sqrt(1-m*sin(x)^2),x=0..Pi/2), real part for m>1}

function M_EllipticPiCx(n,m: extended): extended;
  {-Return the complete elliptic integral of the 3rd kind, n<>1, m<>1, real part}
  { for m>1; Pi(n|m) = integral(1/(1-n*sin(x)^2/sqrt(1-m*sin(x)^2)), x=0..Pi/2)}

function M_EllipticFx(phi,m: extended): extended;
  {-Return the incomplete elliptic integral of the 1st kind,}
  { F(phi,m) = integral(dx/sqrt(1-m*sin(x)^2),x=0..phi), m*sin(phi)^2 <= 1}

function M_EllipticEx(phi,m: extended): extended;
  {-Return the incomplete elliptic integral of the 2nd kind,}
  { E(phi,m) = integral(sqrt(1-m*sin(x)^2),x=0..phi), m*sin(phi)^2 <= 1}

function M_EllipticPix(n,phi,m: extended): extended;
  {-Return the incomplete elliptic integral Pi(n,phi,m) of the 3rd kind}
  { = integral(1/sqrt(1-m*sin(x)^2)/(1-n*sin(x)^2),x=0..phi), with n<>1}
  { m*sin(phi)^2 <= 1, real part if complex.}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Jacobi elliptic and theta functions -------------------}
{---------------------------------------------------------------------------}
{#Z-}

function EllipticModulusx(q: extended): extended;
  {-Return the elliptic modulus k(q) = theta_2(q)^2/theta_3(q)^2, 0 <= q <= 1}

function EllipticNomex(k: extended): extended;
  {-Return the elliptic nome q(k) = exp(-Pi*EllipticCK(k)/EllipticK(k)), |k| < 1}

function jacobi_amx(x,k: extended): extended;
  {-Return the Jacobi amplitude am(x,k)}

function jacobi_arccnx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arccn(x,k), |x| <= 1, x >= sqrt(1 - 1/k^2) if k >= 1}

function jacobi_arccdx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arccd(x,k); |x| <= 1 if |k| < 1; |x| >= 1 if |k| > 1 }

function jacobi_arccsx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arccs(x,k), |x| >= sqrt(k^2-1) for |k|>1}

function jacobi_arcdcx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcdc(x,k); |x| >= 1 if |k| < 1; |x| <= 1 if |k| > 1 }

function jacobi_arcdnx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcdn(x,k), 0 <= x <= 1, x^2 + k^2 > 1 if |k| < 1;  |x| <= 1 if |k| > 1}

function jacobi_arcdsx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcds(x,k), x^2 + k^2 >= 1}

function jacobi_arcncx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcnc(x,k), x >= 1, x^2 <= k^2/(k^2-1) for |k|>1}

function jacobi_arcndx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcnd(x,k), x >= 1, x^2 <= k^2/(1-k^2) if k < 1}

function jacobi_arcnsx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcns(x,k), |x| >= 1, |x| >= k if k>=1}

function jacobi_arcscx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcsc(x,k), |x| <= 1/sqrt(k^2-1) for |k|>1}

function jacobi_arcsdx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcsd(x,k), x^2*(1-k^2) <= 1}

function jacobi_arcsnx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcsn(x,k), |x| <= 1 and |x*k| <= 1}

function jacobi_snx(x,k: extended): extended;
  {-Return the Jacobi elliptic function sn(x,k)}

function jacobi_cnx(x,k: extended): extended;
  {-Return the Jacobi elliptic function cn(x,k)}

function jacobi_dnx(x,k: extended): extended;
  {-Return the Jacobi elliptic function dn(x,k)}

function jacobi_ncx(x,k: extended): extended;
  {-Return the Jacobi elliptic function nc(x,k)}

function jacobi_scx(x,k: extended): extended;
  {-Return the Jacobi elliptic function sc(x,k)}

function jacobi_dcx(x,k: extended): extended;
  {-Return the Jacobi elliptic function dc(x,k)}

function jacobi_ndx(x,k: extended): extended;
  {-Return the Jacobi elliptic function nd(x,k)}

function jacobi_sdx(x,k: extended): extended;
  {-Return the Jacobi elliptic function sd(x,k)}

function jacobi_cdx(x,k: extended): extended;
  {-Return the Jacobi elliptic function cd(x,k)}

function jacobi_nsx(x,k: extended): extended;
  {-Return the Jacobi elliptic function ns(x,k)}

function jacobi_csx(x,k: extended): extended;
  {-Return the Jacobi elliptic function cs(x,k)}

function jacobi_dsx(x,k: extended): extended;
  {-Return the Jacobi elliptic function ds(x,k)}

function jacobi_thetax(n: integer; x,q: extended): extended;
  {-Return the Jacobi theta function theta_n(x,q), n=1..4, 0 <= q < 1}

procedure sncndnx(x,mc: extended; var sn,cn,dn: extended);
  {-Return the Jacobi elliptic functions sn,cn,dn for argument x and}
  { complementary parameter mc.}

function theta1px(q: extended): extended;
  {-Return the derivative  theta1p(q) := d/dx(theta_1(x,q)) at x=0,}
  { = 2*q^(1/4)*sum((-1)^n*(2n+1)*q^(n*(n+1)),n=0..Inf), 0 <= q < 1}

function theta2x(q: extended): extended;
  {-Return Jacobi theta_2(q) = 2*q^(1/4)*sum(q^(n*(n+1)),n=0..Inf) 0 <= q < 1}

function theta3x(q: extended): extended;
  {-Return Jacobi theta_3(q) = 1 + 2*sum(q^(n*n)),n=1..Inf); |q| < 1}

function theta4x(q: extended): extended;
  {-Return Jacobi theta_4(q) = 1 + 2*sum((-1)^n*q^(n*(n+1)),n=1..Inf); |q| < 1}


function ntheta_sx(x, k: extended): extended;
  {-Return the Neville theta_s function, |k| <= 1}

function ntheta_cx(x, k: extended): extended;
  {-Return the Neville theta_c function, |k| <= 1}

function ntheta_dx(x, k: extended): extended;
  {-Return the Neville theta_d function, |k| <= 1}

function ntheta_nx(x, k: extended): extended;
  {-Return the Neville theta_n function, |k| <= 1}


procedure sincos_lemnx(x: extended; var sl,cl: extended);
  {-Return the lemniscate functions sl = sin_lemn(x), cl = cos_lemn(x)}

function sin_lemnx(x: extended): extended;
  {-Return the lemniscate sine function sin_lemn(x)}

function cos_lemnx(x: extended): extended;
  {-Return the lemniscate cosine function cos_lemn(x)}

function arcclx(x: extended): extended;
  {-Return the inverse lemniscate cosine function, |x| <= 1}

function arcslx(x: extended): extended;
  {-Return the inverse lemniscate sine function, |x| <= 1}

{#Z+}
{---------------------------------------------------------------------------}
{--------------- Weierstrass elliptic and modular functions ----------------}
{---------------------------------------------------------------------------}
{#Z-}
function detaix(x: extended): extended;
  {-Return Dedekind eta(i*x), x >= 0}

function emlambdax(y: extended): extended;
  {-Return the elliptic modular function lambda(iy), y >= 0}

function KleinJx(y: extended): extended;
  {-Return Klein's complete invariant J(iy), y>0}

function wplx(x: extended): extended;
  {-Return the Weierstrass function wp(x,1,0)=wpe(x,1/2,0), basic lemniscatic case}

function wpex(x,e1,e2: extended): extended;
  {-Return the Weierstrass function P(x,e1,e2) from the lattice roots e1 < e2}

function wpe_derx(x,e1,e2: extended): extended;
  {-Return Weierstrass P'(x,e1,e2) from the lattice roots e1 < e2}

function wpe_imx(y,e1,e2: extended): extended;
  {-Return the Weierstrass function P(iy,e1,e2) from the lattice roots e1 < e2}

function wpe_invx(y,e1,e2: extended): extended;
  {-Return the smallest positive x with wpe(x)=y, y >= e1}

function wpgx(x,g2,g3: extended): extended;
  {-Return the Weierstrass function P(x,e1,e2) from lattice invariants g2, g3}

function wpg_derx(x,g2,g3: extended): extended;
  {-Return Weierstrass P'(x,e1,e2) from lattice invariants g2, g3}

function wpg_imx(y,g2,g3: extended): extended;
  {-Return the Weierstrass function P(iy, g2, g3)}

function wpg_invx(y,g2,g3: extended): extended;
  {-Return the smallest positive x with wpg(x,g2,g3)=y, y >= e2}

{#Z+}
{---------------------------------------------------------------------------}
{----------------------- Error function and related ------------------------}
{---------------------------------------------------------------------------}
{#Z-}
function dawsonx(x: extended): extended;
  {-Return Dawson's integral: dawson(x) = exp(-x^2)*integral(exp(t^2), t=0..x)}

function dawson2x(p,x: extended): extended;
  {-Return the generalized Dawson integral F(p,x) = exp(-x^p)*integral(exp(t^p), t=0..x); x,p >= 0}

function erfx(x: extended): extended;
  {-Return the error function erf(x) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..x)}

function erfcx(x: extended): extended;
  {-Return the complementary error function erfc(x) = 1-erf(x)}

function erfcex(x: extended): extended;
  {-Return the exponentially scaled complementary error function erfce(x) = exp(x^2)*erfc(x)}

function inerfcx(n: integer; x: extended): extended;
  {-Return the repeated integrals of erfc, n >= -1; scaled with exp(x^2) for x>0}

function erfgx(p,x: extended): extended;
  {-Return the generalized error function integral(exp(-t^p), t=0..x); x,p >= 0}

function erfix(x: extended): extended;
  {-Return the imaginary error function erfi(x) = erf(ix)/i}

function erfhx(x,h: extended): extended;
  {-Accurately compute erf(x+h) - erf(x-h)}

function erf2x(x1,x2: extended):extended;
  {-Accurately compute erf(x2) - erf(x1)}

function erf_invx(x: extended): extended;
  {-Return the inverse function of erf, erf(erf_inv(x)) = x, -1 < x < 1}

function erfc_invx(x: extended): extended;
  {-Return the inverse function of erfc, erfc(erfc_inv(x)) = x, 0 < x < 2}

function erfce_invx(x: extended): extended;
  {-Return the functional inverse of erfce, erfce(erfce_inv(x)) = x, x > 0}

function erfi_invx(y: extended): extended;
  {-Return the inverse imaginary error function: erfi(erf_inv(y))=y}

function erf_px(x: extended): extended;
  {-Return the probability function erf_p = integral(exp(-t^2/2)/sqrt(2*Pi), t=-Inf..x)}

function erf_qx(x: extended): extended;
  {-Return the probability function erf_q = integral(exp(-t^2/2)/sqrt(2*Pi), t=x..Inf)}

function erf_zx(x: extended): extended;
  {-Return the probability function erf_z = exp(-x^2/2)/sqrt(2*Pi)}

function FresnelCx(x: extended): extended;
  {-Return the Fresnel integral C(x)=integral(cos(Pi/2*t^2),t=0..x)}

function FresnelSx(x: extended): extended;
  {-Return the Fresnel integral S(x)=integral(sin(Pi/2*t^2),t=0..x)}

procedure Fresnelx(x: extended; var s,c: extended);
  {-Return the Fresnel integrals S(x)=integral(sin(Pi/2*t^2),t=0..x) and C(x)=integral(cos(Pi/2*t^2),t=0..x)}

procedure FresnelFGx(x: extended; var f,g: extended);
  {-Return the Fresnel auxiliary functions f,g}

function FresnelFx(x: extended): extended;
  {-Return the Fresnel auxiliary function f}

function FresnelGx(x: extended): extended;
  {-Return the Fresnel auxiliary function g}

function gsix(x: extended): extended;
  {-Return the Goodwin-Staton integral = integral(exp(-t*t)/(t+x), t=0..Inf), x <> 0}

function OwenTx(h,a: extended): extended;
  {-Return Owen's T function T(h,a)}

function expint3x(x: extended): extended;
  {-Return the integral(exp(-t^3), t=0..x), x >= 0}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Exponential integrals and related ---------------------}
{---------------------------------------------------------------------------}
{#Z-}
function chix(x: extended): extended;
  {-Return the hyperbolic cosine integral = EulerGamma + ln(|x|) + integral((cosh(t)-1)/t, t=0..|x|)}

function cix(x: extended): extended;
  {-Return the cosine integral, ci(x) = EulerGamma + ln(|x|) + integral((cos(t)-1)/t, t=0..|x|)}

function cinx(x: extended): extended;
  {-Return the entire cosine integral, cin(x) = integral((1-cos(t))/t, t=0..x)}

function cinhx(x: extended): extended;
  {-Return the entire hyperbolic cosine integral, cinh(x) = integral((cosh(t)-1)/t, t=0..x)}

function e1x(x: extended): extended;
  {-Return the exponential integral E1(x) = integral(exp(-x*t)/t, t=1..Inf), x <> 0}

function e1sx(x: extended): extended;
  {-Return E1s(x) = exp(x)*E1(x), x <> 0}

function eix(x: extended): extended;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}

function einx(x: extended): extended;
  {-Return the entire exponential integral ein(x) = integral((1-exp(-t))/t, t=0..x)}

function eisx(x: extended): extended;
  {-Return Eis(x) = exp(-x)*Ei(x), x <> 0}

function eisx2x(x: extended): extended;
  {-Return exp(-x^2)*Ei(x^2), x <> 0}

function ei_invx(x: extended): extended;
  {-Return the functional inverse of Ei(x), ei_inv(ei(x))=x}

function enx(n: longint; x: extended): extended;
  {-Return the exponential integral E_n(x) = integral(exp(-x*t)/t^n, t=1..Inf), x > 0}

function eibetax(n: integer; x: extended): extended;
  {-Return the exponential integral beta(n,x) = int(t^n*exp(-x*t), t=-1..1), n >= 0}

function geix(p,x: extended): extended;
  {-Return the generalized exponential integral E_p(x) = integral(exp(-x*t)/t^p, t=1..Inf), x >= 0}

function lix(x: extended): extended;
  {-Return the logarithmic integral li(x) = PV-integral(1/ln(t), t=0..x), x >= 0, x <> 1}

function li_invx(x: extended): extended;
  {-Return the functional inverse of li(x), li(li_inv(x))=x}

function shix(x: extended): extended;
  {-Return the hyperbolic sine integral, shi(x) = integral(sinh(t)/t, t=0..x)}

function six(x: extended): extended;
  {-Return the sine integral, si(x) = integral(sin(t)/t, t=0..x)}

function ssix(x: extended): extended;
  {-Return the shifted sine integral, ssix(x) = six(x) - Pi/2}

{#Z+}
{---------------------------------------------------------------------------}
{------------------ Orthogonal polynomials and related ---------------------}
{---------------------------------------------------------------------------}
{#Z-}
function chebyshev_tx(n: integer; x: extended): extended;
  {-Return T_n(x), the Chebyshev polynomial of the first kind, degree n}

function chebyshev_ux(n: integer; x: extended): extended;
  {-Return U_n(x), the Chebyshev polynomial of the second kind, degree n}

function chebyshev_vx(n: integer; x: extended): extended;
  {-Return V_n(x), the Chebyshev polynomial of the third kind, degree n>=0}

function chebyshev_wx(n: integer; x: extended): extended;
  {-Return W_n(x), the Chebyshev polynomial of the fourth kind, degree n>=0}

function chebyshev_f1x(v,x: extended): extended;
  {-Return T_v(x), the Chebyshev function the first kind, real part for x<-1}

function gegenbauer_cx(n: integer; a,x: extended): extended;
  {-Return Cn(a,x), the nth Gegenbauer (ultraspherical) polynomial with}
  { parameter a. The degree n must be non-negative; a should be > -0.5 }
  { When a = 0,   C0(0,x) = 1,  and   Cn(0,x) = 2/n*Tn(x)   for n <> 0.}

function hermite_hx(n: integer; x: extended): extended;
  {-Return Hn(x), the nth Hermite polynomial, degree n >= 0}

function hermite_hex(n: integer; x: extended): extended;
  {-Return He_n(x), the nth "probabilists'" Hermite polynomial, degree n >= 0}

function jacobi_px(n: integer; a,b,x: extended): extended;
  {-Return Pn(a,b,x), the nth Jacobi polynomial with parameters a,b. Degree n}
  { must be >= 0; a,b should be > -1 (a+b must not be an integer < -1).}

function laguerrex(n: integer; a,x: extended): extended;
  {-Return Ln(a,x), the nth generalized Laguerre polynomial with parameter a;}
  { degree n must be >= 0. x >=0 and a > -1 are the standard ranges.}

function laguerre_assx(n,m: integer; x: extended): extended;
  {-Return the associated Laguerre polynomial Ln(m,x); n,m >= 0}

function laguerre_lx(n: integer; x: extended): extended;
  {-Return the nth Laguerre polynomial Ln(0,x); n >= 0}

function legendre_px(l: integer; x: extended): extended;
  {-Return P_l(x), the Legendre polynomial/function P_l, degree l}

function legendre_qx(l: integer; x: extended): extended;
  {-Return Q_l(x), the Legendre function of the 2nd kind, degree l >=0, |x| <> 1}

function legendre_plmx(l,m: integer; x: extended): extended;
  {-Return the associated Legendre polynomial P_lm(x)}

function legendre_qlmx(l,m: integer; x: extended): extended;
  {-Return Q(l,m,x), the associated Legendre function of the second kind; l >= 0, l+m >= 0, |x|<>1}

procedure spherical_harmonicx(l, m: integer; theta, phi: extended; var yr,yi: extended);
  {-Return Re and Im of the spherical harmonic function Y_lm(theta,phi)}

function toroidal_plmx(l,m: integer; x: extended): extended;
  {-Return the toroidal harmonic function P(l-0.5,m,x); l,m=0,1; x >= 1}

function toroidal_qlmx(l,m: integer; x: extended): extended;
  {-Return the toroidal harmonic function Q(l-0.5,m,x); l=0,1; x > 1}

function zernike_rx(n,m: integer; r: extended): extended;
  {-Return the Zernike radial polynomial Rnm(r), r >= 0, n >= m >= 0, n-m even}

{#Z+}
{---------------------------------------------------------------------------}
{---------------------- Statistical distributions --------------------------}
{---------------------------------------------------------------------------}
{#Z-}
function beta_pdfx(a, b, x: extended): extended;
  {-Return the probability density function of the beta distribution with}
  { parameters a and b: sfc_beta_pdf = x^(a-1)*(1-x)^(b-1) / beta(a,b)}

function beta_cdfx(a, b, x: extended): extended;
  {-Return the cumulative beta distribution function, a>0, b>0}

function beta_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the beta distribution function. a>0, b>0;}
  { 0 <= y <= 1. Given y the function finds x such that beta_cdf(a, b, x) = y}

function binomial_cdfx(p: extended; n, k: longint): extended;
  {-Return the cumulative binomial distribution function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}

function binomial_pmfx(p: extended; n, k: longint): extended;
  {-Return the binomial distribution probability mass function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}

function cauchy_pdfx(a, b, x: extended): extended;
  {-Return the Cauchy probability density function with location a }
  { and scale b > 0, 1/(Pi*b*(1+((x-a)/b)^2))}

function cauchy_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Cauchy distribution function with location a}
  { and scale b > 0, = 1/2 + arctan((x-a)/b)/Pi}

function cauchy_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Cauchy distribution function}
  { with location a and scale b > 0}

function chi_pdfx(nu: longint; x: extended): extended;
  {-Return the probability density function of the chi distribution, nu>0}

function chi_cdfx(nu: longint; x: extended): extended;
  {-Return the cumulative chi distribution with nu>0 degrees of freedom, x >= 0}

function chi_invx(nu: longint; p: extended): extended;
  {-Return the functional inverse of the chi distribution, nu>0, 0 <= p < 1}

function chi2_pdfx(nu: longint; x: extended): extended;
  {-Return the probability density function of the chi-square distribution, nu>0}

function chi2_cdfx(nu: longint; x: extended): extended;
  {-Return the cumulative chi-square distribution with nu>0 degrees of freedom, x >= 0}

function chi2_invx(nu: longint; p: extended): extended;
  {-Return the functional inverse of the chi-square distribution, nu>0, 0 <= p < 1}

function evt1_pdfx(a, b, x: extended): extended;
  {-Return the probability density function of the Extreme Value Type I distribution}
  { with location a and scale b > 0, result = exp(-(x-a)/b)/b * exp(-exp(-(x-a)/b)) }

function evt1_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Extreme Value Type I distribution function}
  { with location a and scale b > 0;  result = exp(-exp(-(x-a)/b)). }

function evt1_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Extreme Value Type I distribution}
  { function with location a and scale b > 0;  result = a - b*ln(ln(-y)). }

function exp_pdfx(a, alpha, x: extended): extended;
  {-Return the exponential probability density function with location a }
  { and rate alpha > 0, = alpha*exp(-alpha*(x-a)) if x >= a, 0 if x < a.}

function exp_cdfx(a, alpha, x: extended): extended;
  {-Return the cumulative exponential distribution function with location a}
  { and rate alpha > 0, = 1 - exp(-alpha*(x-a)) if x >= a, 0 if x < a.}

function exp_invx(a, alpha, y: extended): extended;
  {-Return the functional inverse of the exponential distribution function with}
  { location a and rate alpha > 0}

function f_pdfx(nu1, nu2: longint; x: extended): extended;
  {-Return the probability density function of the F distribution; x >= 0, nu1, nu2 > 0}

function f_cdfx(nu1, nu2: longint; x: extended): extended;
  {-Return the cumulative F distribution function; x >= 0, nu1, nu2 > 0}

function f_invx(nu1, nu2: longint; y: extended): extended;
  {-Return the functional inverse of the F distribution, nu1, nu2 > 0, 0 <= y <= 1}

function gamma_pdfx(a, b, x: extended): extended;
  {-Return the probability density function of a gamma distribution with shape}
  { a>0, scale b>0: gamma_pdf = x^(a-1)*exp(-x/b)/gamma(a)/b^a, x>0}

function gamma_cdfx(a, b, x: extended): extended;
  {-Return the cumulative gamma distribution function, shape a>0, scale b>0}

function gamma_invx(a, b, p: extended): extended;
  {-Return the functional inverse of the gamma distribution function, shape a>0,}
  { scale b>0, 0 <= p <= 1, i.e. finds x such that gamma_cdf(a, b, x) = p}

function hypergeo_pmfx(n1,n2,n,k: longint): extended;
  {-Return the hypergeometric distribution probability mass function; n,n1,n2 >= 0, n <= n1+n2;}
  { i.e. the probability that among n randomly chosen samples from a container}
  { with n1 type1 objects and n2 type2 objects are exactly k type1 objects:}

function invgamma_pdfx(a, b, x: extended): extended;
  {-Return the probability density function of an inverse gamma distribution}
  { with shape a>0, scale b>0: result = (b/x)^a/x*exp(-b/x)/Gamma(a), x >= 0}

function invgamma_cdfx(a, b, x: extended): extended;
  {-Return the cumulative inverse gamma distribution function, shape a>0, scale}
  { b>0: result = Gamma(a,b/x)/Gamma(a) = Q(a,b/x) = igammaq(a,b/x), x >= 0}

function invgamma_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the inverse gamma distribution function, shape}
  { a>0, scale b>0, 0 <= y <= 1, i.e. find x such that invgamma_cdfx(a, b, x) = y  }

function hypergeo_cdfx(n1,n2,n,k: longint): extended;
  {-Return the cumulative hypergeometric distribution function; n,n1,n2 >= 0, n <= n1+n2}

function kolmogorov_cdfx(x: extended): extended;
  {-Return the limiting form for the cumulative Kolmogorov distribution function}

function kolmogorov_invx(y: extended): extended;
  {-Return the functional inverse of the Kolmogorov distribution}

function kumaraswamy_pdfx(a, b, x: extended): extended;
  {-Return the Kumaraswamy probability density function with shape}
  { parameters a,b>0, 0<=x<=1; result = a*b*x^(a-1)*(1-x^a)^(b-1) }

function kumaraswamy_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Kumaraswamy distribution function with}
  { shape parameters a,b > 0, 0 <= x <= 1; result = 1-(1-x^a)^b}

function kumaraswamy_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Kumaraswamy distribution}
  { with shape parameters a,b > 0; result = [1-(1-y)^(1/b)]^(1/a)}

function laplace_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Laplace distribution function with location a and scale b > 0}

function laplace_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Laplace distribution with location a and scale b > 0}

function laplace_pdfx(a, b, x: extended): extended;
  {-Return the Laplace probability density function with location a}
  { and scale b > 0, result = exp(-abs(x-a)/b) / (2*b)}

function levy_pdfx(a, b, x: extended): extended;
  {-Return the Levy probability density function with}
  { location a and scale parameter b > 0}

function levy_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Levy distribution function with}
  { location a and scale parameter b > 0}

function levy_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Levy distribution}
  { with location a and scale parameter b > 0}

function logistic_pdfx(a, b, x: extended): extended;
  {-Return the logistic probability density function with location a}
  { and scale parameter b > 0, exp(-(x-a)/b)/b/(1+exp(-(x-a)/b))^2}

function logistic_cdfx(a, b, x: extended): extended;
  {-Return the cumulative logistic distribution function with}
  { location a and scale parameter b > 0}

function logistic_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the logistic distribution}
  { with location a and scale parameter b > 0}

function lognormal_pdfx(a, b, x: extended): extended;
  {-Return the log-normal probability density function with}
  { location a and scale parameter b > 0, zero for x <= 0.}

function lognormal_cdfx(a, b, x: extended): extended;
  {-Return the cumulative log-normal distribution function with}
  { location a and scale parameter b > 0, zero for x <= 0.}

function lognormal_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the log-normal distribution}
  { with location a and scale parameter b > 0, 0 < y < 1.}

function logseries_pmfx(a: extended; k: longint): extended;
  {-Return the logarithmic (series) probability mass function}
  { with shape 0 < a < 1, k > 0; result = -a^k/(k*ln(1-a))   }

function logseries_cdfx(a: extended; k: longint): extended;
  {-Return the cumulative logarithmic (series) distribution function with shape 0 < a < 1, k > 0}

function maxwell_pdfx(b, x: extended): extended;
  {-Return the Maxwell probability density function with scale b > 0, x >= 0}

function maxwell_cdfx(b, x: extended): extended;
  {-Return the cumulative Maxwell distribution function with scale b > 0, x >= 0}

function maxwell_invx(b, y: extended): extended;
  {-Return the functional inverse of the Maxwell distribution with scale b > 0}

function moyal_pdfx(a, b, x: extended): extended;
  {-Return the Moyal probability density function with}
  { location a and scale parameter b > 0}

function moyal_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Moyal distribution function with}
  { location a and scale parameter b > 0}

function moyal_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Moyal distribution}
  { with location a and scale parameter b > 0}

function nakagami_pdfx(m, w, x: extended): extended;
  {-Return the probability density function of a gamma distribution with}
  { shape m>0, spread b>0, x>=0: nakagami_pdf = 2x*gamma_pdf(m,w/m,x^2) }

function nakagami_cdfx(m, w, x: extended): extended;
  {-Return the cumulative Nakagami distribution function, shape m>0, spread w>0}

function nakagami_invx(m, w, p: extended): extended;
  {-Return the functional inverse of the Nakagami distribution function, shape m>0,}
  { spread w>0, 0 <= p <= 1, i.e. find x such that nakagami_cdf(m, w, x) = p}

function negbinom_cdfx(p,r: extended; k: longint): extended;
  {-Return the cumulative negative binomial distribution function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}

function negbinom_pmfx(p,r: extended; k: longint): extended;
  {-Return the negative binomial distribution probability mass function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}

function normal_pdfx(mu, sd, x: extended): extended;
  {-Return the normal (Gaussian) probability density function with mean mu}
  { and standard deviation sd>0, exp(-0.5*(x-mu)^2/sd^2) / sqrt(2*Pi*sd^2)}

function normal_cdfx(mu, sd, x: extended): extended;
  {-Return the normal (Gaussian) distribution function}
  { with mean mu and standard deviation sd > 0        }

function normal_invx(mu, sd, y: extended): extended;
  {-Return the functional inverse of the normal (Gaussian) distribution}
  { with mean mu and standard deviation sd > 0, 0 < y < 1.}

function normstd_pdfx(x: extended): extended;
  {-Return the std. normal probability density function exp(-x^2/2)/sqrt(2*Pi)}

function normstd_cdfx(x: extended): extended;
  {-Return the standard normal distribution function}

function normstd_invx(y: extended): extended;
  {-Return the inverse standard normal distribution function, 0 < y < 1.}
  { For x=normstd_inv(y) and y from (0,1), normstd_cdf(x) = y}

function pareto_pdfx(k, a, x: extended): extended;
  {-Return the Pareto probability density function with minimum value k > 0}
  { and shape a, x >= a > 0, result = (a/x)*(k/x)^a}

function pareto_cdfx(k, a, x: extended): extended;
  {-Return the cumulative Pareto distribution function minimum value k > 0}
  { and shape a, x >= a > 0, result = 1-(k/x)^a}

function pareto_invx(k, a, y: extended): extended;
  {-Return the functional inverse of the Pareto distribution with minimum}
  { value k > 0 and shape a, x >= a > 0, result = k/(1-x)^(1/a)}

function poisson_cdfx(mu: extended; k: longint): extended;
  {-Return the cumulative Poisson distribution function with mean mu >= 0}

function poisson_pmfx(mu: extended; k: longint): extended;
  {-Return the Poisson distribution probability mass function with mean mu >= 0}

function rayleigh_pdfx(b, x: extended): extended;
  {-Return the Rayleigh probability density function with}
  { scale b > 0, x >= 0; result = x*exp(-0.5*(x/b)^2)/b^2}

function rayleigh_cdfx(b, x: extended): extended;
  {-Return the cumulative Rayleigh distribution function with scale b > 0, x >= 0}

function rayleigh_invx(b, y: extended): extended;
  {-Return the functional inverse of the Rayleigh distribution with scale b > 0}

function t_pdfx(nu: longint; x: extended): extended;
  {-Return the probability density function of Student's t distribution, nu>0}

function t_cdfx(nu: longint; t: extended): extended;
  {-Return the cumulative Student t distribution with nu>0 degrees of freedom}

function t_invx(nu: longint; p: extended): extended;
  {-Return the functional inverse of Student's t distribution, nu>0, 0 <= p <= 1}

function triangular_pdfx(a, b, c, x: extended): extended;
  {-Return the triangular probability density function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}

function triangular_cdfx(a, b, c, x: extended): extended;
  {-Return the cumulative triangular distribution function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}

function triangular_invx(a, b, c, y: extended): extended;
  {-Return the functional inverse of the triangular distribution with}
  { lower limit a, upper limit b, mode c; a<b, a <= c <= b, 0 <= y <= 1}

function uniform_pdfx(a, b, x: extended): extended;
  {-Return the uniform probability density function on [a,b], a<b}

function uniform_cdfx(a, b, x: extended): extended;
  {-Return the cumulative uniform distribution function on [a,b], a<b}

function uniform_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the uniform distribution on [a,b], a<b}

function wald_pdfx(mu, b, x: extended): extended;
  {-Return the Wald (inverse Gaussian) probability density}
  { function with mean mu > 0, scale b > 0 for x >= 0     }

function wald_cdfx(mu, b, x: extended): extended;
  {-Return the Wald (inverse Gaussian) distribution  }
  { function with mean mu > 0, scale b > 0 for x >= 0}

function wald_invx(mu, b, y: extended): extended;
  {-Return the functional inverse of the Wald (inverse Gaussian)}
  { distribution with mean mu > 0, scale b > 0, 0 <= y < 1.     }

function weibull_pdfx(a, b, x: extended): extended;
  {-Return the Weibull probability density function with shape a > 0}
  { and scale b > 0, result = a*x^(a-1)*exp(-(x/b)^a)/ b^a, x > 0}

function weibull_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Weibull distribution function with}
  { shape parameter a > 0 and scale parameter b > 0}

function weibull_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Weibull distribution}
  { with shape parameter a > 0 and scale parameter b > 0}

function zipf_pmfx(r: extended; k: longint): extended;
  {-Return the Zipf distribution probability mass function k^(-(r+1))/zeta(r+1), r>0, k>0}

function zipf_cdfx(r: extended; k: longint): extended;
  {-Return the cumulative Zipf distribution function H(k,r+1)/zeta(r+1), r>0, k>0}

{#Z+}
{---------------------------------------------------------------------------}
{-------------------- Hypergeometric functions -----------------------------}
{---------------------------------------------------------------------------}
{#Z-}
function hyperg_1F1x(a,b,x: extended): extended;
  {-Return the confluent hypergeometric function 1F1(a,b,x); Kummer's function M(a,b,x)}

function hyperg_1F1rx(a,b,x: extended): extended;
  {-Return the regularized Kummer hypergeometric function 1F1(a,b,x)/Gamma(b)}

function hyperg_ux(a,b,x: extended): extended;
  {-Return Tricomi's confluent hypergeometric function U(a,b,x). If}
  { x<0, then a must be an integer and a<0 or 1+a-b an integer < 0.}

function hyperg_2F1x(a,b,c,x: extended): extended;
  {-Return the Gauss hypergeometric function 2F1(a,b;c;x)}

function hyperg_2F1rx(a,b,c,x: extended): extended;
  {-Return the regularized Gauss hypergeometric function 2F1(a,b,c,x)/Gamma(c)}

function hyperg_2F0x(a,b,x: extended): extended;
  {-Return 2F0(a,b,x), if x>0 then a or b must be a negative integer}

function WhittakerMx(k,m,x: extended): extended;
  {-Return the Whittaker M function = exp(-x/2)*x^(0.5+m) * 1F1(m-k-0.5,2m+1,x)}

function WhittakerWx(k,m,x: extended): extended;
  {-Return the Whittaker W function = exp(-x/2)*x^(0.5+m) * U(m-k-0.5,2m+1,x)}

function hyperg_0F1x(b,x: extended): extended;
  {-Return the confluent hypergeometric limit function 0F1(;b;x)}

function hyperg_0F1rx(b,x: extended): extended;
  {-Return the regularized confluent hypergeometric limit function 0F1(;b;x)/Gamma(b)}

function CylinderDx(v,x: extended): extended;
  {-Return Whittaker's parabolic cylinder function D_v(x)}

function CylinderUx(a,x: extended): extended;
  {-Return the parabolic cylinder function U(a,x)}

function CylinderVx(a,x: extended): extended;
  {-Return the parabolic cylinder function V(a,x) with 2a integer}

function HermiteHx(v,x: extended): extended;
  {-Return the Hermite function H_v(x) of degree v}

{#Z+}
{---------------------------------------------------------------------------}
{--------------------------- Other functions -------------------------------}
{---------------------------------------------------------------------------}
{#Z-}
function agmx(x,y: extended): extended;
  {-Return the arithmetic-geometric mean of |x| and |y|; |x|,|y| < sqrt(MaxExtended)}

function bernoullix(n: integer): extended;
  {-Return the nth Bernoulli number, 0 if n<0 or odd n >= 3}

function bernpolyx(n: integer; x: extended): extended;
  {-Return the Bernoulli polynomial B_n(x), 0 <= n <= MaxBernoulli}

function besselpolyx(n: integer; x: extended): extended;
  {-Return yn(x), the nth Bessel polynomial}

function bringx(x: extended): extended;
  {-Return the Bring radical b := BR(x) with b^5 + b + x = 0}

function catalanx(x: extended): extended;
  {-Return the Catalan function C(x) = binomial(2x,x)/(x+1)}

function cosintx(n: integer; x: extended): extended;
  {-Return cosint(n, x) = integral(cos(t)^n, t=0..x), n>=0}

function debyex(n: integer; x: extended): extended;
  {-Return the Debye function D(n,x) = n/x^n*integral(t^n/(exp(t)-1),t=0..x) of order n>0, x>=0}

function einsteinx(n: integer; x: extended): extended;
  {-Return the Einstein function E_n, n=1..4, x > 0 for n=3,4}

function eulerx(n: integer): extended;
  {-Return the nth Euler number, 0 if n<0 or odd n}

function eulerpolyx(n: integer; x: extended): extended;
  {-Return the Euler polynomial E_n(x), 0 <= n < MaxBernoulli}

function euler_qx(q: extended): extended;
  {-Return the EulerQ function product(1-q^n, n=1..Inf)}

function expnx(n: integer; x: extended): extended;
  {-Return the truncated exponential sum function e_n = sum(x^k/k!, k=0..n), 0 <= n < MAXGAM-1}

function exprelnx(n: longint; x: extended): extended;
  {-Return the relative exponential = (e^x-sum(x^k/k!, k=0..n-1)*n!/x^n}

function fibfunx(v,x: extended): extended;
  {-Return the general Fibonacci function F_v(x)}

function fibpolyx(n: integer; x: extended): extended;
  {-Return the Fibonacci polynomial F_n(x)}

function keplerx(M,e: extended): extended;
  {-Solve Kepler's equation, result x is the eccentric anomaly from the mean anomaly M and the }
  { eccentricity e >= 0; x - e*sin(x) = M, x + x^3/3 = M, or e*sinh(x) - x = M for e <1, =1, >1}

function lucpolyx(n: integer; x: extended): extended;
  {-Return the Lucas polynomial L_n(x)}

function LambertWx(x: extended): extended;
  {-Return the Lambert W function (principal branch), x >= -1/e}

function LambertW1x(x: extended): extended;
  {-Return the Lambert W function (-1 branch), -1/e <= x < 0}

function LangevinLx(x: extended): extended;
  {-Return the Langevin function L(x) = coth(x) - 1/x, L(0) = 0}

function LangevinL_invx(x: extended): extended;
  {-Return the functional inverse of the Langevin function, |x| < 1}

function MarcumQx(m: integer; a,b: extended): extended;
  {-Return the generalized Marcum Q function Q(m,a,b), a,b >= 0}

function omegax(x: extended): extended;
  {-Return the Wright omega function, i.e. the solution w of w + ln(w) = x}

function RiemannRx(x: extended): extended;
  {-Return the Riemann prime counting function R(x), x >= 1/1024}

function RiemannR_invx(x: extended): extended;
  {-Return the functional inverse of R(x), R(RiemannR_inv(x))=x, x >= 1.125}

function rrcfx(q: extended): extended;
  {-Return the Rogers-Ramanujan continued fraction for |q| < 1}

function sinintx(n: integer; x: extended): extended;
  {-Return sinint(n, x) = integral(sin(t)^n, t=0..x), n>=0}

function transportx(n: integer; x: extended): extended;
  {-Return the transport integral J_n(x) for x >= 0, n >= 2,}
  { J_n(x) = integral(t^n*exp(t)/(exp(t)-1)^2, t=0..x)      }

implementation

uses
  AMath,
  sfGamma,   {Gamma function and related}
  sfGamma2,  {inverse/incomplete Gamma/Beta}
  sfZeta,    {Zeta functions and polylogarithms}
  sfErf,     {Error function and related}
  sfBessel,  {Bessel functions}
  sfBessl2,  {Bessel functions: Airy, Kelvin, Coulomb}
  sfEllInt,  {Elliptic integrals}
  sfExpInt,  {Exponential integrals and related}
  sfSDist,   {Statistical distributions}
  sfPoly,    {Orthogonal polynomials and related}
  sfHyperG,  {Hypergeometric functions}
  sfMisc;    {Other functions}

{---------------------------------------------------------------------------}
{--------------------------- Bessel functions ------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function bessel_i0x(x: extended): extended;
  {-Return I0(x), the modified Bessel function of the 1st kind, order zero}
begin
  bessel_i0x := sfc_i0(x);
end;


{---------------------------------------------------------------------------}
function bessel_i0ex(x: extended): extended;
  {-Return I0(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order zero}
begin
  bessel_i0ex := sfc_i0e(x);
end;


{---------------------------------------------------------------------------}
function bessel_i1x(x: extended): extended;
  {-Return I1(x), the modified Bessel function of the 1st kind, order one}
begin
  bessel_i1x := sfc_i1(x);
end;


{---------------------------------------------------------------------------}
function bessel_i1ex(x: extended): extended;
  {-Return I1(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order one}
begin
  bessel_i1ex := sfc_i1e(x);
end;


{---------------------------------------------------------------------------}
function bessel_inx(n: integer; x: extended): extended;
  {-Return I_n(x), the modified Bessel function of the 1st kind, order n.}
begin
  bessel_inx := sfc_in(n,x);
end;


{---------------------------------------------------------------------------}
function bessel_ivx(v, x: extended): extended;
  {-Return I_v(x), the modified Bessel function of the 1st kind, order v.}
begin
  bessel_ivx := sfc_iv(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_ivex(v, x: extended): extended;
  {-Return I_v(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order v.}
begin
  bessel_ivex := sfc_ive(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_j0x(x: extended): extended;
  {-Return J0(x), the Bessel function of the 1st kind, order zero}
begin
  bessel_j0x := sfc_j0(x);
end;


{---------------------------------------------------------------------------}
function bessel_j1x(x: extended): extended;
  {-Return J1(x), the Bessel function of the 1st kind, order one}
begin
  bessel_j1x := sfc_j1(x);
end;


{---------------------------------------------------------------------------}
function bessel_jnx(n: integer; x: extended): extended;
  {-Return J_n(x), the Bessel function of the 1st kind, order n; not suitable for large n or x.}
begin
  bessel_jnx := sfc_jn(n,x);
end;


{---------------------------------------------------------------------------}
function bessel_jvx(v, x: extended): extended;
  {-Return J_v(x), the Bessel function of the 1st kind, order v; not suitable for large v.}
begin
  bessel_jvx := sfc_jv(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_k0x(x: extended): extended;
  {-Return K0(x), the modified Bessel function of the 2nd kind, order zero, x>0}
begin
  bessel_k0x := sfc_k0(x);
end;


{---------------------------------------------------------------------------}
function bessel_k0ex(x: extended): extended;
  {-Return K0(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order zero, x>0}
begin
  bessel_k0ex := sfc_k0e(x);
end;


{---------------------------------------------------------------------------}
function bessel_k1x(x: extended): extended;
  {-Return K1(x), the modified Bessel function of the 2nd kind, order one, x>0}
begin
  bessel_k1x := sfc_k1(x);
end;


{---------------------------------------------------------------------------}
function bessel_k1ex(x: extended): extended;
  {-Return K1(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order one, x>0}
begin
  bessel_k1ex := sfc_k1e(x);
end;


{---------------------------------------------------------------------------}
function bessel_knx(n: integer; x: extended): extended;
  {-Return K_n(x), the modified Bessel function of the 2nd kind, order n, x>0, not suitable for large n}
begin
  bessel_knx := sfc_kn(n,x);
end;


{---------------------------------------------------------------------------}
function bessel_kvx(v, x: extended): extended;
  {-Return K_v(x), the modified Bessel function of the 2nd kind, order v, x>0}
begin
  bessel_kvx := sfc_kv(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_kvex(v, x: extended): extended;
  {-Return K_v(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order v, x>0}
begin
  bessel_kvex := sfc_kve(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_y0x(x: extended): extended;
  {-Return Y0(x), the Bessel function of the 2nd kind, order zero; x>0}
begin
  bessel_y0x := sfc_y0(x);
end;


{---------------------------------------------------------------------------}
function bessel_y1x(x: extended): extended;
  {-Return Y1(x), the Bessel function of the 2nd kind, order one; x>0}
begin
  bessel_y1x := sfc_y1(x);
end;


{---------------------------------------------------------------------------}
function bessel_ynx(n: integer; x: extended): extended;
  {-Return Y_n(x), the Bessel function of the 2nd kind, order n, x>0, not suitable for large n or x}
begin
  bessel_ynx := sfc_yn(n,x);
end;


{---------------------------------------------------------------------------}
function bessel_yvx(v, x: extended): extended;
  {-Return Y_v(x), the Bessel function of the 2nd kind, order v; x > 0; not suitable for large v.}
begin
  bessel_yvx := sfc_yv(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_lambdax(v,x: extended): extended;
  {-Compute lambda(v,x) = Gamma(v+1)*J(v,x)/(0.5x)^v for v,x >= 0}
begin
  bessel_lambdax := sfc_blam(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_i0_intx(u: extended): extended;
  {-Return the integral int(bessel_i0(x), x = 0..u)}
begin
  bessel_i0_intx := sfc_i0int(u);
end;


{---------------------------------------------------------------------------}
function bessel_j0_intx(u: extended): extended;
  {-Return the integral int(bessel_j0(x), x = 0..u)}
begin
  bessel_j0_intx := sfc_j0int(u);
end;


{---------------------------------------------------------------------------}
function bessel_k0_intx(u: extended): extended;
  {-Return the integral int(bessel_k0(x), x = 0..u), u >= 0}
begin
  bessel_k0_intx := sfc_k0int(u);
end;


{---------------------------------------------------------------------------}
function bessel_y0_intx(u: extended): extended;
  {-Return the integral int(bessel_y0(x), x = 0..u), u >= 0}
begin
  bessel_y0_intx := sfc_y0int(u);
end;


{---------------------------------------------------------------------------}
function airy_aix(x: extended): extended;
  {-Return the Airy function Ai(x)}
begin
  airy_aix := sfc_airy_ai(x);
end;


{---------------------------------------------------------------------------}
function airy_aipx(x: extended): extended;
  {-Return the Airy function Ai'(x)}
begin
  airy_aipx := sfc_airy_aip(x);
end;


{---------------------------------------------------------------------------}
function airy_aisx(x: extended): extended;
  {-Return the scaled Airy function Ai(x) if x <= 0, Ai(x)*exp(2/3*x^1.5) for x > 0}
begin
  airy_aisx := sfc_airy_ais(x);
end;


{---------------------------------------------------------------------------}
function airy_bix(x: extended): extended;
  {-Return the Airy function Bi(x)}
begin
  airy_bix := sfc_airy_bi(x);
end;


{---------------------------------------------------------------------------}
function airy_bipx(x: extended): extended;
  {-Return the Airy function Bi'(x)}
begin
  airy_bipx := sfc_airy_bip(x);
end;


{---------------------------------------------------------------------------}
function airy_bisx(x: extended): extended;
  {-Return the scaled Airy function Bi(x) if x <= 0, Bi(x)*exp(-2/3*x^1.5) for x > 0}
begin
  airy_bisx := sfc_airy_bis(x);
end;


{---------------------------------------------------------------------------}
function airy_gix(x: extended): extended;
  {-Return the Airy/Scorer function Gi(x) = 1/Pi*integral(sin(x*t+t^3/3), t=0..INF)}
begin
  airy_gix := sfc_airy_gi(x);
end;


{---------------------------------------------------------------------------}
function airy_hix(x: extended): extended;
  {-Return the Airy/Scorer function Hi(x) = 1/Pi*integral(exp(x*t-t^3/3), t=0..INF)}
begin
  airy_hix := sfc_airy_hi(x);
end;


{---------------------------------------------------------------------------}
function kelvin_berx(x: extended): extended;
  {-Return the Kelvin function ber(x)}
begin
  kelvin_berx := sfc_ber(x);
end;


{---------------------------------------------------------------------------}
function kelvin_beix(x: extended): extended;
  {-Return the Kelvin function bei(x)}
begin
  kelvin_beix := sfc_bei(x);
end;


{---------------------------------------------------------------------------}
function kelvin_kerx(x: extended): extended;
  {-Return the Kelvin function ker(x), x > 0}
begin
  kelvin_kerx := sfc_ker(x);
end;


{---------------------------------------------------------------------------}
function kelvin_keix(x: extended): extended;
  {-Return the Kelvin function kei(x), x >= 0}
begin
  kelvin_keix := sfc_kei(x);
end;


{---------------------------------------------------------------------------}
procedure kelvin_kerkeix(x: extended; var kr, ki: extended);
  {-Return the Kelvin functions kr=ker(x), ki=kei(x), x > 0}
begin
  sfc_ker_kei(x,kr,ki);
end;


{---------------------------------------------------------------------------}
procedure kelvin_berbeix(x: extended; var br, bi: extended);
  {-Return the Kelvin functions br=ber(x), bi=bei(x)}
begin
  sfc_ber_bei(x,br,bi);
end;


{---------------------------------------------------------------------------}
procedure kelvin_derx(x: extended; var berp, beip, kerp, keip: extended);
  {-Return the derivatives of the zero order Kelvin functions, x >= 0}
begin
  sfc_kelvin_der(x, berp, beip, kerp, keip);
end;


{---------------------------------------------------------------------------}
function kelvin_berpx(x: extended): extended;
  {-Return the Kelvin function ber'(x), x >= 0}
begin
  kelvin_berpx := sfc_berp(x);
end;


{---------------------------------------------------------------------------}
function kelvin_beipx(x: extended): extended;
  {-Return the Kelvin function bei'(x), x >= 0}
begin
  kelvin_beipx := sfc_beip(x);
end;


{---------------------------------------------------------------------------}
function kelvin_kerpx(x: extended): extended;
  {-Return the Kelvin function ker'(x), x > 0}
begin
  kelvin_kerpx := sfc_kerp(x);
end;


{---------------------------------------------------------------------------}
function kelvin_keipx(x: extended): extended;
  {-Return the Kelvin function kei'(x), x >= 0}
begin
  kelvin_keipx := sfc_keip(x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_jnx(n: integer; x: extended): extended;
  {-Return j_n(x), the spherical Bessel function of the 1st kind, order n}
begin
  sph_bessel_jnx := sfc_sph_jn(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_ynx(n: integer; x: extended): extended;
  {-Return y_n(x), the spherical Bessel function of the 2nd kind, order n >=0 , x<>0}
begin
  sph_bessel_ynx := sfc_sph_yn(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_inx(n: integer; x: extended): extended;
  {-Return i_n(x), the modified spherical Bessel function of the 1st/2nd kind, order n}
begin
  sph_bessel_inx := sfc_sph_in(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_inex(n: integer; x: extended): extended;
  {-Return i_n(x)*exp(-|x|), the exponentially scaled modified spherical Bessel function of the 1st/2nd kind, order n}
begin
  sph_bessel_inex := sfc_sph_ine(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_knx(n: integer; x: extended): extended;
  {-Return k_n(x), the modified spherical Bessel function of the 3rd kind, order n, x>0}
begin
  sph_bessel_knx := sfc_sph_kn(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_knex(n: integer; x: extended): extended;
  {-Return k_n(x)*exp(x), the exponentially scaled modified spherical Bessel function of the 3rd kind, order n, x>0}
begin
  sph_bessel_knex := sfc_sph_kne(n,x);
end;


{---------------------------------------------------------------------------}
function struve_h0x(x: extended): extended;
 {-Return H0(x), the Struve function of order 0}
begin
  struve_h0x := sfc_struve_h0(x);
end;


{---------------------------------------------------------------------------}
function struve_h1x(x: extended): extended;
  {-Return H1(x), the Struve function of order 1}
begin
  struve_h1x := sfc_struve_h1(x);
end;


{---------------------------------------------------------------------------}
function struve_l0x(x: extended): extended;
  {-Return L0(x), the modified Struve function of order 0}
begin
  struve_l0x := sfc_struve_l0(x);
end;


{---------------------------------------------------------------------------}
function struve_l1x(x: extended): extended;
  {-Return L1(x), the modified Struve function of order 1}
begin
  struve_l1x := sfc_struve_l1(x);
end;


{---------------------------------------------------------------------------}
function struve_hx(v, x: extended): extended;
  {-Return H_v(x), the Struve function of order v, x < 0 only if v is an integer.}
begin
  struve_hx := sfc_struve_h(v,x);
end;


{---------------------------------------------------------------------------}
function struve_lx(v, x: extended): extended;
  {-Return L_v(x), the modified Struve function of order v, x < 0 only if v is an integer.}
begin
  struve_lx := sfc_struve_l(v,x);
end;


{---------------------------------------------------------------------------}
function CoulombCLx(L: integer; eta: extended): extended;
  {-Return the normalizing constant CL for Coulomb wave function, L >= 0}
begin
  CoulombCLx := sfc_coulcl(L,eta);
end;


{---------------------------------------------------------------------------}
function CoulombSLx(L: integer; eta: extended): extended;
  {-Return the Coulomb phase shift sigma_L(eta) for L >= 0}
begin
  CoulombSLx := sfc_cshift(L,eta);
end;


{---------------------------------------------------------------------------}
function CoulombFx(L: integer; eta, x: extended): extended;
  {-Return the regular Coulomb wave functions FL(eta,x) for L >= 0, x > 0}
begin
  CoulombFx := sfc_coul_f(L,eta,x);
end;


{---------------------------------------------------------------------------}
procedure CoulombFFpx(L: integer; eta, x: extended; var fc,fcp: extended; var ifail: integer);
  {-Return the regular Coulomb wave functions fc=FL(eta,x) and fcp=FL'(eta,x)}
  { for L >= 0, x > 0. No error if ifail=0, failure if ifail < 0, if ifail=1 }
  { then the very inaccurate JWKB approximation has been computed.}
begin
  sfc_coul_ffp(L,eta,x, fc,fcp, ifail);
end;


{---------------------------------------------------------------------------}
procedure CoulombGGpx(L: integer; eta, x: extended; var gc,gcp: extended; var ifail: integer);
  {-Return the irregular Coulomb wave functions gc=GL(eta,x) and gcp=GL'(eta,x)}
  { for  L >= 0, x > 0. No error if ifail=0, failure if ifail < 0, if ifail=1  }
  { then the very inaccurate JWKB approximation has been computed.}
begin
  sfc_coul_ggp(L,eta,x, gc,gcp, ifail);
end;


{---------------------------------------------------------------------------}
function SynchFx(x: extended): extended;
  {-Return the first synchrotron function F(x) = integral(x*BesselK(5/3,t), t=x..INF) for x >= 0}
begin
  SynchFx := sfc_synchf(x);
end;


{---------------------------------------------------------------------------}
function SynchGx(x: extended): extended;
  {-Return the second synchrotron function G(x) = x*BesselK(2/3,x) for x >= 0}
begin
  SynchGx := sfc_synchg(x);
end;


{---------------------------------------------------------------------------}
{---------------------- Gamma function and related -------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function signgammax(x: extended): extended;
  {-Return sign(gamma(x)), useless for 0 or negative integer}
begin
  signgammax := sfc_signgamma(x);
end;


{---------------------------------------------------------------------------}
function gammax(x: extended): extended;
  {-Return gamma(x), x <= MAXGAMX; invalid if x is a non-positive integer}
begin
  gammax := sfc_gamma(x);
end;


{---------------------------------------------------------------------------}
function gamma1pm1x(x: extended): extended;
  {-Return gamma(1+x)-1 with increased accuracy for x near 0}
begin
  gamma1pm1x := sfc_gamma1pm1(x);
end;


{---------------------------------------------------------------------------}
function gammastarx(x: extended): extended;
  {-Return Temme's gammastar(x) = gamma(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x)), x>0.}
  { For large x the asymptotic expansion is gammastar(x) = 1 + 1/12x + O(1/x^2)}
begin
  gammastarx := sfc_gstar(x);
end;


{---------------------------------------------------------------------------}
function gamma_delta_ratiox(x,d: extended): extended;
  {-Return gamma(x)/gamma(x+d), accurate even for |d| << |x|}
begin
  gamma_delta_ratiox := sfc_gamma_delta_ratio(x,d);
end;


{---------------------------------------------------------------------------}
function gamma_ratiox(x,y: extended): extended;
  {-Return gamma(x)/gamma(y)}
begin
  gamma_ratiox := sfc_gamma_ratio(x,y);
end;


{---------------------------------------------------------------------------}
function pochhammerx(a,x: extended): extended;
  {-Return the Pochhammer symbol gamma(a+x)/gamma(a)}
begin
  pochhammerx := sfc_pochhammer(a,x);
end;


{---------------------------------------------------------------------------}
function poch1x(a,x: extended): extended;
  {-Return (pochhammer(a,x)-1)/x, psi(a) if x=0; accurate even for small |x|}
begin
  poch1x := sfc_poch1(a,x);
end;


{---------------------------------------------------------------------------}
function binomialx(n,k: integer): extended;
  {-Return the binomial coefficient 'n choose k'}
begin
  binomialx := sfc_binomial(n,k,false);
end;


{---------------------------------------------------------------------------}
function lnbinomialx(n,k: longint): extended;
  {-Return ln(binomial(n,k)), n >= k >= 0}
begin
  lnbinomialx := sfc_lnbinomial(n,k);
end;


{---------------------------------------------------------------------------}
function facx(n: integer): extended;
  {-Return the factorial n!, n<MAXGAM-1; INF if n<0}
begin
  facx := sfc_fac(n);
end;


{---------------------------------------------------------------------------}
function dfacx(n: integer): extended;
  {-Return the double factorial n!!, n<=MAXDFACX; INF for even n<0}
begin
  dfacx := sfc_dfac(n);
end;


{---------------------------------------------------------------------------}
function lnfacx(n: longint): extended;
  {-Return ln(n!), INF if n<0}
begin
  lnfacx := sfc_lnfac(n);
end;


{---------------------------------------------------------------------------}
function lngammax(x: extended): extended;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, invalid if x is a non-positive integer}
  { Function signgamma can be used if the sign of gamma(x) is needed.}
begin
  if IsNanOrInf(x) then lngammax := x
  else if abs(x) > MAXLGMX then lngammax := PosInf_x
  else lngammax := sfc_lngamma(x);
end;


{---------------------------------------------------------------------------}
function lngammasx(x: extended; var s: integer): extended;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, s=-1,1 is the sign of gamma}
begin
  lngammasx := sfc_lngammas(x,s);
end;

{---------------------------------------------------------------------------}
function lngamma1px(x: extended): extended;
  {-Return ln(|gamma(1+x)|) with increased accuracy for x near 0}
begin
  if IsNanOrInf(x) then lngamma1px := x
  else if abs(x) > MAXLGMX then lngamma1px := PosInf_x
  else lngamma1px := sfc_lngamma1p(x);
end;


{---------------------------------------------------------------------------}
function rgammax(x: extended): extended;
  {-Return the reciprocal gamma function rgamma = 1/gamma(x)}
begin
  if x>=1755.5 then rgammax := 0.0
  else rgammax := sfc_rgamma(x);
end;


{---------------------------------------------------------------------------}
procedure incgammax(a,x: extended; var p,q: extended);
  {-Return the normalised incomplete gamma functions P and Q, a>=0, x>=0}
  { P(a,x) = integral(exp(-t)*t^(a-1), t=0..x  )/gamma(a)}
  { Q(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf)/gamma(a)}
begin
  sfc_incgamma(a,x,p,q);
end;


{---------------------------------------------------------------------------}
function igammapx(a,x: extended): extended;
  {-Return the normalised lower incomplete gamma function P(a,x), a>=0, x>=0}
  { P(a,x) = integral(exp(-t)*t^(a-1), t=0..x)/gamma(a)}
begin
  igammapx := sfc_igammap(a,x);
end;


{---------------------------------------------------------------------------}
function igammaqx(a,x: extended): extended;
  {-Return the normalised upper incomplete gamma function Q(a,x), a>=0, x>=0}
  { Q(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf)/gamma(a)}
begin
  igammaqx := sfc_igammaq(a,x);
end;


{---------------------------------------------------------------------------}
function igammax(a,x: extended): extended;
  {-Return the non-normalised upper incomplete gamma function}
  { GAMMA(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf). If x<0 }
  { the real part is returned and a must be <> -1, -2, ...   }
begin
  igammax := sfc_igamma(a,x);
end;


{---------------------------------------------------------------------------}
function igammalx(a,x: extended): extended;
  {-Return the non-normalised lower incomplete gamma function}
  { gamma(a,x) = integral(exp(-t)*t^(a-1), t=0..x); a<>0,-1,-2,..}
begin
  igammalx := sfc_igammal(a,x);
end;


{---------------------------------------------------------------------------}
function igammatx(a,x: extended): extended;
  {-Return Tricomi's entire incomplete gamma function igammatx(a,x)}
  { = igammal(a,x)/gamma(a)/x^a = P(a,x)/x^a }
begin
  igammatx := sfc_git(a,x);
end;


{---------------------------------------------------------------------------}
procedure incgamma_invx(a,p,q: extended; var x: extended; var ierr: integer);
  {-Return the inverse normalised incomplete gamma function, i.e. calculate}
  { x with P(a,x)=p and Q(a,x)=q. Input parameter a>0, p>=0, q>0 and p+q=1.}
  { ierr is >= 0 for success, < 0 for input errors or iterations failures. }
begin
  sfc_incgamma_inv(a,p,q,x,ierr);
end;


{---------------------------------------------------------------------------}
function igamma_invx(a,p,q: extended): extended;
  {-Return the inverse normalised incomplete gamma function, i.e. calculate}
  { x with P(a,x)=p and Q(a,x)=q. Input parameter a>0, p>=0, q>0 and p+q=1.}
begin
  igamma_invx := sfc_igamma_inv(a,p,q);
end;


{---------------------------------------------------------------------------}
function igammap_invx(a,p: extended): extended;
  {-Inverse incomplete gamma: return x with P(a,x)=p, a>=0, 0<=p<1}
begin
  igammap_invx := sfc_igammap_inv(a,p);
end;


{---------------------------------------------------------------------------}
function igammap_derx(a,x: extended): extended;
  {-Return the partial derivative with respect to x of the normalised}
  { lower incomplete gamma function P(a,x), x >= 0, a <> 0,-1,-2 ...}
begin
  igammap_derx := sfc_igammap_der(a,x);
end;


{---------------------------------------------------------------------------}
function igammaq_invx(a,q: extended): extended;
  {-Inverse complemented incomplete gamma: return x with Q(a,x)=q, a>=0, 0<q<=1}
begin
  igammaq_invx := sfc_igammaq_inv(a,q);
end;


{---------------------------------------------------------------------------}
function lngamma_invx(y: extended): extended;
  {-Inverse of lngamma: return x with lngamma(x) = y, y >= -0.12142, x > 1.4616}
begin
  lngamma_invx := sfc_ilng(y);
end;


{---------------------------------------------------------------------------}
function inv_gammax(y: extended): extended;
  {-Inverse of gamma: return x with gamma(x) = y, y >= 0.8857421875}
begin
  inv_gammax := sfc_invgam(y);
end;


{---------------------------------------------------------------------------}
function psix(x: extended): extended;
  {-Return the psi (digamma) function of x, INF if x is a non-positive integer}
begin
  if x=PosInf_x then psix := x
  else psix := sfc_psi(x);
end;


{---------------------------------------------------------------------------}
function psistarx(x: extended): extended;
  {-Return psi(x) - ln(x), x > 0}
begin
  psistarx := sfc_psistar(x);
end;


{---------------------------------------------------------------------------}
function psi_invx(y: extended): extended;
  {-Inverse of psi, return x with psi(x)=y, y <= ln_MaxExt}
begin
  psi_invx := sfc_ipsi(y);
end;


{---------------------------------------------------------------------------}
function trigammax(x: extended): extended;
  {-Return the trigamma function of x, INF if x is a negative integer}
begin
  trigammax := sfc_trigamma(x);
end;


{---------------------------------------------------------------------------}
function tetragammax(x: extended): extended;
  {-Return the tetragamma function psi''(x), NAN/RTE if x is a negative integer}
begin
  tetragammax := sfc_tetragamma(x);
end;


{---------------------------------------------------------------------------}
function pentagammax(x: extended): extended;
  {-Return the pentagamma function psi'''(x), INF if x is a negative integer}
begin
  pentagammax := sfc_pentagamma(x);
end;


{---------------------------------------------------------------------------}
function polygammax(n: integer; x: extended): extended;
  {-Return the polygamma function: n'th derivative of psi; n>=0, x>0 for n>MAXGAMX}
  { Note: The accuracy may be reduced for n >= MAXGAMX due to ln/exp operations.}
begin
  polygammax := sfc_polygamma(n,x);
end;


{---------------------------------------------------------------------------}
function lnBarnesGx(x: extended): extended;
  {-Return ln(BarnesG(x)), real part for x < 0}
begin
  lnBarnesGx := sfc_lnbg(x);
end;


{---------------------------------------------------------------------------}
function BatemanGx(x: extended): extended;
  {-Return the Bateman function G(x) = psi((x+1)/2) - psi(x/2), x<>0,-1,-2,...}
begin
  BatemanGx := sfc_batemang(x);
end;


{---------------------------------------------------------------------------}
function lnbetax(x,y: extended): extended;
  {-Return the logarithm of |beta(x,y)|=|gamma(x)*gamma(y)/gamma(x+y)|}
begin
  lnbetax := sfc_lnbeta(x,y);
end;


{---------------------------------------------------------------------------}
function betax(x,y: extended): extended;
  {-Return the function beta(x,y)=gamma(x)*gamma(y)/gamma(x+y)}
begin
  betax := sfc_beta(x,y);
end;


{---------------------------------------------------------------------------}
function ibetax(a, b, x: extended): extended;
  {-Return the normalised incomplete beta function, a>0, b>0, 0 <= x <= 1}
  { ibetax = integral(t^(a-1)*(1-t)^(b-1) / betax(a,b), t=0..x)}
begin
  ibetax := sfc_ibeta(a,b,x);
end;


{---------------------------------------------------------------------------}
function ibeta_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the normalised incomplete beta function}
  { with a > 0, b > 0, and 0 <= y <= 1.}
begin
  ibeta_invx := sfc_ibeta_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function beta3x(a, b, x: extended): extended;
  {-Return the non-normalised incomplete beta function B_x(a,b)}
  { for 0<=x<=1, B_x = integral(t^(a-1)*(1-t)^(b-1), t=0..x).  }
begin
  beta3x := sfc_nnbeta(a,b,x);
end;


{---------------------------------------------------------------------------}
{------------------- Zeta functions and polylogarithms   -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function zetax(s: extended): extended;
  {-Return the Riemann zeta function at s, s<>1}
begin
  zetax := sfc_zeta(s);
end;


{---------------------------------------------------------------------------}
function zetaintx(n: integer): extended;
  {-Return zeta(n) for integer arguments, n<>1}
begin
  zetaintx := sfc_zetaint(n);
end;


{---------------------------------------------------------------------------}
function zeta1px(x: extended): extended;
  {-Return the Riemann zeta function at 1+x, x<>0}
begin
  zeta1px := sfc_zeta1p(x);
end;


{---------------------------------------------------------------------------}
function etax(s: extended): extended;
  {-Return the Dirichlet eta function}
begin
  etax := sfc_eta(s);
end;


{---------------------------------------------------------------------------}
function etaintx(n: integer): extended;
  {-Return the Dirichlet function eta(n) for integer arguments}
begin
  etaintx := sfc_etaint(n);
end;


{---------------------------------------------------------------------------}
function etam1x(s: extended): extended;
  {-Return Dirichlet eta(s)-1}
begin
  etam1x := sfc_etam1(s);
end;


{---------------------------------------------------------------------------}
function zetam1x(s: extended): extended;
  {-Return Riemann zeta(s)-1, s<>1}
begin
  zetam1x := sfc_zetam1(s);
end;


{---------------------------------------------------------------------------}
function zetahx(s,a: extended): extended;
  {-Return the Hurwitz zeta function zetah(s,a) = sum(1/(i+a)^s, i=0..INF), s<>1, a>0}
begin
  zetahx := sfc_zetah(s,a);
end;


{---------------------------------------------------------------------------}
function LerchPhix(z,s,a: extended): extended;
  {-Return the Lerch transcendent Phi(z,s,a) = sum(z^n/(n+a)^s, n=0..INF), z <= 1, s >= -1, a >= 0}
begin
  LerchPhix := sfc_lerch(z,s,a);
end;


{---------------------------------------------------------------------------}
function DirichletBetax(s: extended): extended;
  {-Return the Dirichlet beta function sum((-1)^n/(2n+1)^s, n=0..INF)}
begin
  DirichletBetax := sfc_dbeta(s);
end;


{---------------------------------------------------------------------------}
function DirichletLambdax(s: extended): extended;
  {-Return the Dirichlet lambda function sum(1/(2n+1)^s, n=0..INF), s<>1}
begin
  DirichletLambdax := sfc_dlambda(s);
end;


{---------------------------------------------------------------------------}
function LegendreChix(s, x: extended): extended;
  {-Return Legendre's Chi-function chi(s,x); s>=0; real part if |x| > 1}
begin
  LegendreChix := sfc_lchi(s, x);
end;


{---------------------------------------------------------------------------}
function primezetax(x: extended): extended;
  {-Return the prime zeta function P(x) = sum(1/p^x, p prime), x > 1/5; }
  { for x<1 the real part of P(x) is returned.}
begin
  primezetax := sfc_pz(x);
end;


{---------------------------------------------------------------------------}
function harmonicx(x: extended): extended;
  {-Return the harmonic number function H(x) = psi(x+1) + EulerGamma}
begin
  harmonicx := sfc_harmonic(x);
end;


{---------------------------------------------------------------------------}
function harmonic2x(x,r: extended): extended;
  {-Return the generalized harmonic function H(x,r) = zeta(r)-zetah(r,x+1); x >= -1}
begin
  harmonic2x := sfc_harmonic2(x,r);
end;


{---------------------------------------------------------------------------}
function cl2x(x: extended): extended;
  {-Return the Clausen function: integral(-ln(2*|sin(t/2)|),t=0..x) = Im(Li_2(exp(ix)))}
begin
  cl2x := sfc_cl2(x);
end;


{---------------------------------------------------------------------------}
function ti2x(x: extended): extended;
  {-Return the inverse tangent integral, ti2(x) = integral(arctan(t)/t, t=0..x)}
begin
  ti2x := sfc_ti2(x);
end;


{---------------------------------------------------------------------------}
function tix(s,x: extended): extended;
  {-Return the inverse tangent integral of order s >= 0}
begin
  tix := sfc_ti(s,x);
end;


{---------------------------------------------------------------------------}
function lobachevsky_cx(x: extended): extended;
  {-Return the Lobachevski function L(x) = integral(-ln(|cos(t)|), t=0..x)}
begin
  lobachevsky_cx := sfc_llci(x);
end;


{---------------------------------------------------------------------------}
function lobachevsky_sx(x: extended): extended;
  {-Return the Lobachevski function Lambda(x) = integral(-ln(|2sin(t)|), t=0..x)}
begin
  lobachevsky_sx := sfc_llsi(x);
end;


{---------------------------------------------------------------------------}
function dilogx(x: extended): extended;
  {-Return dilog(x) = Re(Li_2(x)), Li_2(x) = -integral(ln(1-t)/t, t=0..x)}
begin
  dilogx := sfc_dilog(x);
end;


{---------------------------------------------------------------------------}
function trilogx(x: extended): extended;
  {-Return the trilogarithm function trilog(x) = Re(Li_3(x))}
begin
  trilogx := sfc_trilog(x);
end;


{---------------------------------------------------------------------------}
function polylogx(n: integer; x: extended): extended;
  {-Return the polylogarithm Li_n(x) of integer order; real part for n>0,x>1}
begin
  polylogx := sfc_polylog(n,x);
end;


{---------------------------------------------------------------------------}
function polylogrx(s, x: extended): extended;
  {-Return the polylogarithm Li_s(x) of real order s >= -1, x <= 256;}
  { s > 0 if x > 256; real part if x > 1}
begin
  polylogrx := sfc_polylogr(s,x);
end;


{---------------------------------------------------------------------------}
function bose_einsteinx(s,x: extended): extended;
  {-Return the Bose-Einstein integral of real order s >= -1}
begin
  bose_einsteinx := sfc_beint(s,x);
end;


{---------------------------------------------------------------------------}
function fermi_diracx(n: integer; x: extended): extended;
  {-Return the integer order Fermi-Dirac integral F_n(x) = 1/n!*integral(t^n/(exp(t-x)+1), t=0..INF)}
begin
  fermi_diracx := sfc_fermi_dirac(n,x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac_rx(s,x: extended): extended;
  {-Return the Fermi-Dirac integral of real order s >= -1}
begin
  fermi_dirac_rx := sfc_fdr(s,x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac_m05x(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}
begin
  fermi_dirac_m05x := sfc_fdm05(x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac_p05x(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}
begin
  fermi_dirac_p05x := sfc_fdp05(x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac_p15x(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}
begin
  fermi_dirac_p15x := sfc_fdp15(x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac_p25x(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}
begin
  fermi_dirac_p25x := sfc_fdp25(x);
end;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Legendre style) -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function comp_ellint_1x(k: extended): extended;
  {-Return the complete elliptic integral of the 1st kind, real part if |k|>1}
begin
  comp_ellint_1x := sfc_EllipticK(k);
end;


{---------------------------------------------------------------------------}
function comp_ellint_2x(k: extended): extended;
  {-Return the complete elliptic integral of the 2nd kind, real part if |k|>1}
begin
  comp_ellint_2x := sfc_EllipticEC(k);
end;


{---------------------------------------------------------------------------}
function comp_ellint_3x(nu,k: extended): extended;
  {-Return the complete elliptic integral of the 3rd kind, |k|<>1, nu<>1; real part for |k|>1}
begin
  comp_ellint_3x := sfc_EllipticPiC(nu,k);
end;


{---------------------------------------------------------------------------}
function comp_ellint_bx(k: extended): extended;
  {-Return the complete elliptic integral B(k) = (E(k) - kc^2*K(k))/k^2, real part for |k|>1}
begin
  comp_ellint_bx := sfc_cel_b(k);
end;


{---------------------------------------------------------------------------}
function comp_ellint_dx(k: extended): extended;
  {-Return the complete elliptic integral D(k) = (K(k) - E(k))/k^2, real part for |k|>1}
begin
  comp_ellint_dx := sfc_cel_d(k);
end;


{---------------------------------------------------------------------------}
function ellint_1x(phi,k: extended): extended;
  {-Return the Legendre elliptic integral F(phi,k) of the 1st kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}
begin
  ellint_1x := sfc_ellint_1(phi,k);
end;


{---------------------------------------------------------------------------}
function ellint_2x(phi,k: extended): extended;
  {-Return the Legendre elliptic integral E(phi,k) of the 2nd kind}
  { = integral(sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}
begin
  ellint_2x := sfc_ellint_2(phi,k);
end;


{---------------------------------------------------------------------------}
function ellint_3x(phi,nu,k: extended): extended;
  {-Return the Legendre elliptic integral PI(phi,nu,k) of the 3rd kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2)/(1-nu*sin(x)^2),x=0..phi) with  }
  { |k*sin(phi)|<=1, returns Cauchy principal value if nu*sin(phi)^2>1}
begin
  ellint_3x := sfc_ellint_3(phi,nu,k);
end;


{---------------------------------------------------------------------------}
function ellint_bx(phi,k: extended): extended;
  {-Return the Legendre elliptic integral B(phi,k) = (E(phi,k) - kc^2*F(phi,k))/k^2}
  { = integral(cos(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }
begin
  ellint_bx := sfc_ellint_b(phi,k);
end;

{---------------------------------------------------------------------------}
function ellint_dx(phi,k: extended): extended;
  {-Return the Legendre elliptic integral D(phi,k) = (F(phi,k) - E(phi,k))/k^2 }
  { = integral(sin(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }
begin
  ellint_dx := sfc_ellint_d(phi,k);
end;


{---------------------------------------------------------------------------}
function heuman_lambdax(phi,k: extended): extended;
  {-Return Heuman's function Lambda_0(phi,k) = F(phi,k')/K(k') + 2/Pi*K(k)*Z(phi,k'), |k|<=1}
begin
  heuman_lambdax := sfc_hlambda(phi,k);
end;


{---------------------------------------------------------------------------}
function jacobi_zetax(phi,k: extended): extended;
  {-Return the Jacobi Zeta function Z(phi,k) = E(phi,k) - E(k)/K(k)*F(phi,k), |k|<=1}
begin
  jacobi_zetax := sfc_jzeta(phi,k);
end;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Carlson style) --------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function ell_rcx(x,y: extended): extended;
  {-Return Carlson's degenerate elliptic integral RC; x>=0, y<>0}
begin
  ell_rcx := sfc_ell_rc(x,y);
end;


{---------------------------------------------------------------------------}
function ell_rfx(x,y,z: extended): extended;
  {-Return Carlson's elliptic integral of the 1st kind; x,y,z >=0, at most one =0}
begin
  ell_rfx := sfc_ell_rf(x,y,z);
end;


{---------------------------------------------------------------------------}
function ell_rdx(x,y,z: extended): extended;
  {-Return Carlson's elliptic integral of the 2nd kind; z>0; x,y >=0, at most one =0}
begin
  ell_rdx := sfc_ell_rd(x,y,z);
end;


{---------------------------------------------------------------------------}
function ell_rgx(x,y,z: extended): extended;
  {-Return Carlson's completely symmetric elliptic integral of the 2nd kind; x,y,z >= 0}
begin
  ell_rgx := sfc_ell_rg(x,y,z);
end;


{---------------------------------------------------------------------------}
function ell_rjx(x,y,z,r: extended): extended;
  {-Return Carlson's elliptic integral of the 3rd kind; r<>0; x,y,z >=0, at most one =0}
begin
  ell_rjx := sfc_ell_rj(x,y,z,r);
end;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Bulirsch style) -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function cel1x(kc: extended): extended;
  {-Return Bulirsch's complete elliptic integral of the 1st kind, kc<>0}
begin
  cel1x := sfc_cel1(kc);
end;


{---------------------------------------------------------------------------}
function cel2x(kc, a, b: extended): extended;
  {-Return Bulirsch's complete elliptic integral of the 2nd kind, kc<>0}
begin
  cel2x := sfc_cel2(kc,a,b);
end;


{---------------------------------------------------------------------------}
function celx(kc, p, a, b: extended): extended;
  {-Return Bulirsch's general complete elliptic integral, kc<>0, Cauchy principle value if p<0}
begin
  celx := sfc_cel(kc,p,a,b);
end;


{---------------------------------------------------------------------------}
function el1x(x,kc: extended): extended;
  {-Return Bulirsch's incomplete elliptic integral of the 1st kind}
begin
  el1x := sfc_el1(x,kc);
end;


{---------------------------------------------------------------------------}
function el2x(x,kc,a,b: extended): extended;
  {-Return Bulirsch's incomplete elliptic integral of the 2nd kind, kc<>0}
begin
  el2x := sfc_el2(x,kc,a,b);
end;


{---------------------------------------------------------------------------}
function el3x(x,kc,p: extended): extended;
  {-Return Bulirsch's incomplete elliptic integral of the 3rd kind, 1+p*x^2<>0}
begin
  el3x := sfc_el3(x,kc,p);
end;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Maple V style) --------------------}
{---------------------------------------------------------------------------}

function EllipticFx(z,k: extended): extended;
  {-Return the incomplete elliptic integral of the 1st kind; |z|<=1, |k*z|<1}
begin
  EllipticFx := sfc_EllipticF(z,k);
end;


{---------------------------------------------------------------------------}
function EllipticKx(k: extended): extended;
  {-Return the complete elliptic integral of the 1st kind, |k| < 1}
begin
  EllipticKx := sfc_EllipticK(k);
end;


{---------------------------------------------------------------------------}
function EllipticKimx(k: extended): extended;
  {-Return K(i*k), the complete elliptic integral of the 1st kind with}
  { imaginary modulus = integral(1/sqrt(1-x^2)/sqrt(1+k^2*x^2),x=0..1)}
begin
  EllipticKimx := sfc_EllipticKim(k);
end;


{---------------------------------------------------------------------------}
function EllipticCKx(k: extended): extended;
  {-Return the complementary complete elliptic integral of the 1st kind, k<>0}
begin
  EllipticCKx := sfc_EllipticCK(k);
end;


{---------------------------------------------------------------------------}
function EllipticEx(z,k: extended): extended;
  {-Return the incomplete elliptic integrals of the 2nd kind, |z|<=1, |k*z| <= 1}
begin
  EllipticEx := sfc_EllipticE(z,k);
end;


{---------------------------------------------------------------------------}
function EllipticECx(k: extended): extended;
  {-Return the complete elliptic integral of the 2nd kind, real part for |k|>1}
begin
  EllipticECx := sfc_EllipticEC(k);
end;


{---------------------------------------------------------------------------}
function EllipticECimx(k: extended): extended;
  {-Return E(i*k), the complete elliptic integral of the 2nd kind with}
  { imaginary modulus = integral(sqrt(1+k^2*x^2)/sqrt(1-x^2),x=0..1)  }
begin
  EllipticECimx := sfc_EllipticECim(k);
end;


{---------------------------------------------------------------------------}
function EllipticCEx(k: extended): extended;
  {-Return the complementary complete elliptic integral of the 2nd kind}
begin
  EllipticCEx := sfc_EllipticCE(k);
end;


{---------------------------------------------------------------------------}
function EllipticPix(z,nu,k: extended): extended;
  {-Return the incomplete elliptic integral of the 3rd kind, |z|<=1, |k*z|<1}
begin
  EllipticPix := sfc_EllipticPi(z,nu,k);
end;


{---------------------------------------------------------------------------}
function EllipticPiCx(nu,k: extended): extended;
  {-Return the complete elliptic integral of the 3rd kind, |k|<>1, nu<>1; real part for |k|>1}
begin
  EllipticPiCx := sfc_EllipticPiC(nu,k);
end;


{---------------------------------------------------------------------------}
function EllipticPiCimx(nu,k: extended): extended;
  {-Return Pi(nu, k*i), the complete elliptic integral of the 3rd kind}
  { with imaginary modulus, nu <> 1, real part if nu > 1}
begin
  EllipticPiCimx := sfc_EllipticPiCim(nu,k);
end;


{---------------------------------------------------------------------------}
function EllipticCPix(nu,k: extended): extended;
  {-Return the complementary complete elliptic integral of the 3rd kind, k<>0, nu<>1}
begin
  EllipticCPix := sfc_EllipticCPi(nu,k);
end;


{---------------------------------------------------------------------------}
{---------------- Elliptic integrals (Mathematica style) -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function M_EllipticKx(m: extended): extended;
  {-Return the complete elliptic integral of the 1st kind,}
  { K(m) = integral(dx/sqrt(1-m*sin(x)^2),x=0..Pi/2), real part for m>1}
begin
  M_EllipticKx := sfc_M_EllipticK(m);
end;


{---------------------------------------------------------------------------}
function M_EllipticECx(m: extended): extended;
  {-Return the complete elliptic integral of the 2nd kind,}
  { E(m) = integral(sqrt(1-m*sin(x)^2),x=0..Pi/2), real part for m>1}
begin
  M_EllipticECx := sfc_M_EllipticEC(m);
end;


{---------------------------------------------------------------------------}
function M_EllipticPiCx(n,m: extended): extended;
  {-Return the complete elliptic integral of the 3rd kind, n<>1, m<>1, real part}
  { for m>1; Pi(n|m) = integral(1/(1-n*sin(x)^2/sqrt(1-m*sin(x)^2)), x=0..Pi/2)}
begin
  M_EllipticPiCx := sfc_M_EllipticPiC(n,m);
end;


{---------------------------------------------------------------------------}
function M_EllipticFx(phi,m: extended): extended;
  {-Return the incomplete elliptic integral of the 1st kind,}
  { F(phi,m) = integral(dx/sqrt(1-m*sin(x)^2),x=0..phi), m*sin(phi)^2 <= 1}
begin
  M_EllipticFx := sfc_M_EllipticF(phi,m);
end;


{---------------------------------------------------------------------------}
function M_EllipticEx(phi,m: extended): extended;
  {-Return the incomplete elliptic integral of the 2nd kind,}
  { E(phi,m) = integral(sqrt(1-m*sin(x)^2),x=0..phi), m*sin(phi)^2 <= 1}
begin
  M_EllipticEx := sfc_M_EllipticE(phi,m);
end;


{---------------------------------------------------------------------------}
function M_EllipticPix(n,phi,m: extended): extended;
  {-Return the incomplete elliptic integral Pi(n,phi,m) of the 3rd kind}
  { = integral(1/sqrt(1-m*sin(x)^2)/(1-n*sin(x)^2),x=0..phi), with n<>1}
  { m*sin(phi)^2 <= 1, real part if complex.}
begin
  M_EllipticPix := sfc_M_EllipticPi(n,phi,m);
end;

{---------------------------------------------------------------------------}
{------------------- Jacobi elliptic and theta functions -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function EllipticModulusx(q: extended): extended;
  {-Return the elliptic modulus k(q) = theta_2(q)^2/theta_3(q)^2, 0 <= q <= 1}
begin
  EllipticModulusx := sfc_ellmod(q);
end;


{---------------------------------------------------------------------------}
function EllipticNomex(k: extended): extended;
  {-Return the elliptic nome q(k) = exp(-Pi*EllipticCK(k)/EllipticK(k)), |k| < 1}
begin
  EllipticNomex := sfc_nome(k);
end;


{---------------------------------------------------------------------------}
function jacobi_amx(x,k: extended): extended;
  {-Return the Jacobi amplitude am(x,k)}
begin
  jacobi_amx := sfc_jam(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_thetax(n: integer; x,q: extended): extended;
  {-Return the Jacobi theta function theta_n(x,q), n=1..4, 0 <= q < 1}
begin
  jacobi_thetax := sfc_jtheta(n,x,q);
end;


{---------------------------------------------------------------------------}
procedure sncndnx(x,mc: extended; var sn,cn,dn: extended);
  {-Return the Jacobi elliptic functions sn,cn,dn for argument x and}
  { complementary parameter mc.}
begin
  sfc_sncndn(x,mc,sn,cn,dn);
end;


{---------------------------------------------------------------------------}
function jacobi_snx(x,k: extended): extended;
  {-Return the Jacobi elliptic function sn(x,k)}
begin
  jacobi_snx := sfc_jacobi_sn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_cnx(x,k: extended): extended;
  {-Return the Jacobi elliptic function cn(x,k)}
begin
  jacobi_cnx := sfc_jacobi_cn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_dnx(x,k: extended): extended;
  {-Return the Jacobi elliptic function dn(x,k)}
begin
  jacobi_dnx := sfc_jacobi_dn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_ncx(x,k: extended): extended;
  {-Return the Jacobi elliptic function nc(x,k)}
begin
  jacobi_ncx := sfc_jacobi_nc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_scx(x,k: extended): extended;
  {-Return the Jacobi elliptic function sc(x,k)}
begin
  jacobi_scx := sfc_jacobi_sc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_dcx(x,k: extended): extended;
  {-Return the Jacobi elliptic function dc(x,k)}
begin
  jacobi_dcx := sfc_jacobi_dc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_ndx(x,k: extended): extended;
  {-Return the Jacobi elliptic function nd(x,k)}
begin
  jacobi_ndx := sfc_jacobi_nd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_sdx(x,k: extended): extended;
  {-Return the Jacobi elliptic function sd(x,k)}
begin
  jacobi_sdx := sfc_jacobi_sd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_cdx(x,k: extended): extended;
  {-Return the Jacobi elliptic function cd(x,k)}
begin
  jacobi_cdx := sfc_jacobi_cd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_nsx(x,k: extended): extended;
  {-Return the Jacobi elliptic function ns(x,k)}
begin
  jacobi_nsx := sfc_jacobi_ns(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_csx(x,k: extended): extended;
  {-Return the Jacobi elliptic function cs(x,k)}
begin
  jacobi_csx := sfc_jacobi_cs(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_dsx(x,k: extended): extended;
  {-Return the Jacobi elliptic function ds(x,k)}
begin
  jacobi_dsx := sfc_jacobi_ds(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arccnx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arccn(x,k), |x| <= 1, x >= sqrt(1 - 1/k^2) if k >= 1}
begin
  jacobi_arccnx := sfc_jacobi_arccn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arccdx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arccd(x,k); |x| <= 1 if |k| < 1; |x| >= 1 if |k| > 1 }
begin
  jacobi_arccdx := sfc_jacobi_arccd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arccsx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arccs(x,k), |x| >= sqrt(k^2-1) for |k|>1}
begin
  jacobi_arccsx := sfc_jacobi_arccs(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcdcx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcdc(x,k); |x| >= 1 if |k| < 1; |x| <= 1 if |k| > 1 }
begin
  jacobi_arcdcx := sfc_jacobi_arcdc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcdnx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcdn(x,k), 0 <= x <= 1, x^2 + k^2 > 1 if |k| < 1;  |x| <= 1 if |k| > 1}
begin
  jacobi_arcdnx := sfc_jacobi_arcdn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcdsx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcds(x,k), x^2 + k^2 >= 1}
begin
  jacobi_arcdsx := sfc_jacobi_arcds(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcncx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcnc(x,k), x >= 1, x^2 <= k^2/(k^2-1) for |k|>1}
begin
  jacobi_arcncx := sfc_jacobi_arcnc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcndx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcnd(x,k), x >= 1, x^2 <= k^2/(1-k^2) if k < 1}
begin
  jacobi_arcndx := sfc_jacobi_arcnd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcnsx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcns(x,k), |x| >= 1, |x| >= k if k>=1}
begin
  jacobi_arcnsx := sfc_jacobi_arcns(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcscx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcsc(x,k), |x| <= 1/sqrt(k^2-1) for |k|>1}
begin
  jacobi_arcscx := sfc_jacobi_arcsc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcsdx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcsd(x,k), x^2*(1-k^2) <= 1}
begin
  jacobi_arcsdx := sfc_jacobi_arcsd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcsnx(x,k: extended): extended;
  {-Return the inverse Jacobi elliptic function arcsn(x,k), |x| <= 1 and |x*k| <= 1}
begin
  jacobi_arcsnx := sfc_jacobi_arcsn(x,k);
end;


{---------------------------------------------------------------------------}
function theta1px(q: extended): extended;
  {-Return the derivative  theta1p(q) := d/dx(theta_1(x,q)) at x=0,}
  { = 2*q^(1/4)*sum((-1)^n*(2n+1)*q^(n*(n+1)),n=0..Inf), 0 <= q < 1}
begin
  theta1px := sfc_theta1p(q);
end;


{---------------------------------------------------------------------------}
function theta2x(q: extended): extended;
  {-Return Jacobi theta_2(q) = 2*q^(1/4)*sum(q^(n*(n+1)),n=0..Inf) 0 <= q < 1}
begin
  theta2x := sfc_theta2(q);
end;


{---------------------------------------------------------------------------}
function theta3x(q: extended): extended;
  {-Return Jacobi theta_3(q) = 1 + 2*sum(q^(n*n)),n=1..Inf); |q| < 1}
begin
  theta3x := sfc_theta3(q);
end;


{---------------------------------------------------------------------------}
function theta4x(q: extended): extended;
  {-Return Jacobi theta_4(q) = 1 + 2*sum((-1)^n*q^(n*(n+1)),n=1..Inf); |q| < 1}
begin
  theta4x := sfc_theta4(q);
end;


{---------------------------------------------------------------------------}
function ntheta_sx(x, k: extended): extended;
  {-Return the Neville theta_s function, |k| <= 1}
begin
  ntheta_sx := sfc_ntheta_s(x, k);
end;


{---------------------------------------------------------------------------}
function ntheta_cx(x, k: extended): extended;
  {-Return the Neville theta_c function, |k| <= 1}
begin
  ntheta_cx := sfc_ntheta_c(x, k);
end;


{---------------------------------------------------------------------------}
function ntheta_dx(x, k: extended): extended;
  {-Return the Neville theta_d function, |k| <= 1}
begin
  ntheta_dx := sfc_ntheta_d(x, k);
end;


{---------------------------------------------------------------------------}
function ntheta_nx(x, k: extended): extended;
  {-Return the Neville theta_n function, |k| <= 1}
begin
  ntheta_nx := sfc_ntheta_n(x, k);
end;


{---------------------------------------------------------------------------}
procedure sincos_lemnx(x: extended; var sl,cl: extended);
  {-Return the lemniscate functions sl = sin_lemn(x), cl = cos_lemn(x)}
begin
  sfc_lemn(x,sl,cl);
end;


{---------------------------------------------------------------------------}
function sin_lemnx(x: extended): extended;
  {-Return the lemniscate sine function sin_lemn(x)}
var
  sl,cl: extended;
begin
  sfc_lemn(x,sl,cl);
  sin_lemnx := sl;
end;


{---------------------------------------------------------------------------}
function cos_lemnx(x: extended): extended;
  {-Return the lemniscate cosine function cos_lemn(x)}
var
  sl,cl: extended;
begin
  sfc_lemn(x,sl,cl);
  cos_lemnx := cl;
end;


{---------------------------------------------------------------------------}
function arcclx(x: extended): extended;
  {-Return the inverse lemniscate cosine function, |x| <= 1}
begin
  arcclx := sfc_arccl(x);
end;


{---------------------------------------------------------------------------}
function arcslx(x: extended): extended;
  {-Return the inverse lemniscate cosine function, |x| <= 1}
begin
  arcslx := sfc_arcsl(x);
end;


{---------------------------------------------------------------------------}
{--------------- Weierstrass elliptic and modular functions ----------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function detaix(x: extended): extended;
  {-Return Dedekind eta(i*x), x >= 0}
begin
  detaix := sfc_detai(x);
end;


{---------------------------------------------------------------------------}
function emlambdax(y: extended): extended;
  {-Return the elliptic modular function lambda(iy), y >= 0}
begin
  emlambdax := sfc_emlambda(y);
end;


{---------------------------------------------------------------------------}
function KleinJx(y: extended): extended;
  {-Return Klein's complete invariant J(iy), y>0}
begin
  KleinJx := sfc_KleinJ(y);
end;


{---------------------------------------------------------------------------}
function wplx(x: extended): extended;
  {-Return the Weierstrass function wp(x,1,0)=wpe(x,1/2,0), basic lemniscatic case}
begin
  wplx := sfc_wpl(x);
end;

{---------------------------------------------------------------------------}
function wpex(x,e1,e2: extended): extended;
  {-Return the Weierstrass function P(x,e1,e2) from the lattice roots e1 < e2}
begin
  wpex := sfc_wpe(x,e1,e2);
end;


{---------------------------------------------------------------------------}
function wpe_derx(x,e1,e2: extended): extended;
  {-Return Weierstrass P'(x,e1,e2) from the lattice roots e1 < e2}
begin
  wpe_derx := sfc_wpe_der(x,e1,e2);
end;


{---------------------------------------------------------------------------}
function wpe_imx(y,e1,e2: extended): extended;
  {-Return the Weierstrass function P(iy,e1,e2) from the lattice roots e1 < e2}
begin
  wpe_imx := sfc_wpe_im(y,e1,e2);
end;


{---------------------------------------------------------------------------}
function wpe_invx(y,e1,e2: extended): extended;
  {-Return the smallest positive x with wpe(x)=y, y >= e1}
begin
  wpe_invx := sfc_wpe_inv(y,e1,e2);
end;


{---------------------------------------------------------------------------}
function wpgx(x,g2,g3: extended): extended;
  {-Return the Weierstrass function P(x,e1,e2) from lattice invariants g2, g3}
begin
  wpgx := sfc_wpg(x,g2,g3);
end;


{---------------------------------------------------------------------------}
function wpg_derx(x,g2,g3: extended): extended;
  {-Return Weierstrass P'(x,e1,e2) from lattice invariants g2, g3}
begin
  wpg_derx := sfc_wpg_der(x,g2,g3);
end;


{---------------------------------------------------------------------------}
function wpg_imx(y,g2,g3: extended): extended;
  {-Return the Weierstrass function P(iy, g2, g3)}
begin
  wpg_imx := sfc_wpg_im(y,g2,g3);
end;


{---------------------------------------------------------------------------}
function wpg_invx(y,g2,g3: extended): extended;
  {-Return the smallest positive x with wpg(x,g2,g3)=y, y >= e2}
begin
  wpg_invx := sfc_wpg_inv(y,g2,g3);
end;


{---------------------------------------------------------------------------}
{----------------------- Error function and related ------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function dawsonx(x: extended): extended;
  {-Return Dawson's integral: dawson(x) = exp(-x^2)*integral(exp(t^2), t=0..x)}
begin
  dawsonx := sfc_dawson(x);
end;


{---------------------------------------------------------------------------}
function dawson2x(p,x: extended): extended;
  {-Return the generalized Dawson integral F(p,x) = exp(-x^p)*integral(exp(t^p), t=0..x); x,p >= 0}
begin
  dawson2x := sfc_gendaw(p,x);
end;


{---------------------------------------------------------------------------}
function erfx(x: extended): extended;
  {-Return the error function erf(x) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..x)}
begin
  erfx := sfc_erf(x);
end;


{---------------------------------------------------------------------------}
function erfcx(x: extended): extended;
  {-Return the complementary error function erfc(x) = 1-erf(x)}
begin
  erfcx := sfc_erfc(x);
end;


{---------------------------------------------------------------------------}
function erfcex(x: extended): extended;
  {-Return the exponentially scaled complementary error function erfce(x) = exp(x^2)*erfc(x)}
begin
  erfcex := sfc_erfce(x);
end;


{---------------------------------------------------------------------------}
function inerfcx(n: integer; x: extended): extended;
  {-Return the repeated integrals of erfc, n >= -1; scaled with exp(x^2) for x>0}
begin
  inerfcx := sfc_inerfc(n,x);
end;


{---------------------------------------------------------------------------}
function erfgx(p,x: extended): extended;
  {-Return the generalized error function integral(exp(-t^p), t=0..x); x,p >= 0}
begin
  erfgx := sfc_erfg(p,x);
end;


{---------------------------------------------------------------------------}
function erfix(x: extended): extended;
  {-Return the imaginary error function erfi(x) = erf(ix)/i}
begin
  erfix := sfc_erfi(x);
end;


{---------------------------------------------------------------------------}
function erfhx(x,h: extended): extended;
  {-Accurately compute erf(x+h) - erf(x-h)}
begin
  erfhx := sfc_erfh(x,h);
end;


{---------------------------------------------------------------------------}
function erf2x(x1,x2: extended):extended;
  {-Accurately compute erf(x2) - erf(x1)}
begin
  erf2x := sfc_erf2(x1,x2);
end;


{---------------------------------------------------------------------------}
function erf_invx(x: extended): extended;
  {-Return the inverse function of erf, erf(erf_inv(x)) = x, -1 < x < 1}
begin
  erf_invx := sfc_erf_inv(x);
end;


{---------------------------------------------------------------------------}
function erfi_invx(y: extended): extended;
  {-Return the inverse imaginary error function: erfi(erf_inv(y))=y}
begin
  erfi_invx := sfc_erfi_inv(y);
end;


{---------------------------------------------------------------------------}
function erfc_invx(x: extended): extended;
  {-Return the inverse function of erfc, erfc(erfc_inv(x)) = x, 0 < x < 2}
begin
  erfc_invx := sfc_erfc_inv(x);
end;


{---------------------------------------------------------------------------}
function erfce_invx(x: extended): extended;
  {-Return the functional inverse of erfce, erfce(erfce_inv(x)) = x, x > 0}
begin
  erfce_invx := sfc_erfce_inv(x);
end;


{---------------------------------------------------------------------------}
function erf_px(x: extended): extended;
  {-Return the probability function erf_p = integral(exp(-t^2/2)/sqrt(2*Pi), t=-Inf..x)}
begin
  erf_px := sfc_erf_p(x);
end;


{---------------------------------------------------------------------------}
function erf_qx(x: extended): extended;
  {-Return the probability function erf_q = integral(exp(-t^2/2)/sqrt(2*Pi), t=x..Inf)}
begin
  erf_qx := sfc_erf_q(x);
end;


{---------------------------------------------------------------------------}
function erf_zx(x: extended): extended;
  {-Return the probability function erf_z = exp(-x^2/2)/sqrt(2*Pi)}
begin
  erf_zx := sfc_erf_z(x);
end;


{---------------------------------------------------------------------------}
function expint3x(x: extended): extended;
  {-Return the integral(exp(-t^3), t=0..x), x >= 0}
begin
  expint3x := sfc_expint3(x);
end;


{---------------------------------------------------------------------------}
procedure Fresnelx(x: extended; var s,c: extended);
  {-Return the Fresnel integrals S(x)=integral(sin(Pi/2*t^2),t=0..x) and C(x)=integral(cos(Pi/2*t^2),t=0..x)}
begin
  sfc_fresnel(x,s,c);
end;


{---------------------------------------------------------------------------}
function FresnelCx(x: extended): extended;
  {-Return the Fresnel integral C(x)=integral(cos(Pi/2*t^2),t=0..x)}
var
  s,c: extended;
begin
  sfc_fresnel(x,s,c);
  FresnelCx := c;
end;


{---------------------------------------------------------------------------}
function FresnelSx(x: extended): extended;
  {-Return the Fresnel integral S(x)=integral(sin(Pi/2*t^2),t=0..x)}
var
  s,c: extended;
begin
  sfc_fresnel(x,s,c);
  FresnelSx := s;
end;


{---------------------------------------------------------------------------}
procedure FresnelFGx(x: extended; var f,g: extended);
  {-Return the Fresnel auxiliary functions f,g}
begin
  sfc_fresnel_fg(x,f,g);
end;


{---------------------------------------------------------------------------}
function FresnelFx(x: extended): extended;
  {-Return the Fresnel auxiliary function f}
var
  f,g: extended;
begin
  sfc_fresnel_fg(x,f,g);
  FresnelFx := f;
end;


{---------------------------------------------------------------------------}
function FresnelGx(x: extended): extended;
  {-Return the Fresnel auxiliary function g}
var
  f,g: extended;
begin
  sfc_fresnel_fg(x,f,g);
  FresnelGx := g;
end;


{---------------------------------------------------------------------------}
function gsix(x: extended): extended;
  {-Return the Goodwin-Staton integral = integral(exp(-t*t)/(t+x), t=0..Inf), x <> 0}
begin
  gsix := sfc_gsi(x);
end;


{---------------------------------------------------------------------------}
function OwenTx(h,a: extended): extended;
  {-Return Owen's T function T(h,a)}
begin
  OwenTx := sfc_owent(h,a);
end;


{---------------------------------------------------------------------------}
{------------------- Exponential integrals and related ---------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function chix(x: extended): extended;
  {-Return the hyperbolic cosine integral = EulerGamma + ln(|x|) + integral((cosh(t)-1)/t, t=0..|x|)}
begin
  chix := sfc_chi(x);
end;


{---------------------------------------------------------------------------}
function cix(x: extended): extended;
  {-Return the cosine integral, ci(x) = EulerGamma + ln(|x|) + integral((cos(t)-1)/t, t=0..|x|)}
begin
  cix := sfc_ci(x);
end;


{---------------------------------------------------------------------------}
function cinx(x: extended): extended;
  {-Return the entire cosine integral, cin(x) = integral((1-cos(t))/t, t=0..x)}
begin
  cinx := sfc_cin(x);
end;


{---------------------------------------------------------------------------}
function cinhx(x: extended): extended;
  {-Return the entire hyperbolic cosine integral, cinh(x) = integral((cosh(t)-1)/t, t=0..x)}
begin
  cinhx := sfc_cinh(x);
end;


{---------------------------------------------------------------------------}
function e1x(x: extended): extended;
  {-Return the exponential integral E1(x) = integral(exp(-x*t)/t, t=1..Inf), x <> 0}
begin
  e1x := sfc_e1(x);
end;


{---------------------------------------------------------------------------}
function e1sx(x: extended): extended;
  {-Return E1s(x) = exp(x)*E1(x), x <> 0}
begin
  e1sx := sfc_e1s(x);
end;


{---------------------------------------------------------------------------}
function eix(x: extended): extended;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}
begin
  eix := sfc_ei(x);
end;


{---------------------------------------------------------------------------}
function eisx(x: extended): extended;
  {-Return Eis(x) = exp(-x)*Ei(x), x <> 0}
begin
  eisx := sfc_eis(x);
end;



{---------------------------------------------------------------------------}
function einx(x: extended): extended;
  {-Return the entire exponential integral ein(x) = integral((1-exp(-t))/t, t=0..x)}
begin
  einx := sfc_ein(x);
end;


{---------------------------------------------------------------------------}
function eisx2x(x: extended): extended;
  {-Return exp(-x^2)*Ei(x^2)}
begin
  eisx2x := sfc_eisx2(x);
end;


{---------------------------------------------------------------------------}
function ei_invx(x: extended): extended;
  {-Return the functional inverse of Ei(x), ei_inv(ei(x))=x}
begin
  ei_invx := sfc_ei_inv(x);
end;


{---------------------------------------------------------------------------}
function enx(n: longint; x: extended): extended;
  {-Return the exponential integral E_n(x) = integral(exp(-x*t)/t^n, t=1..Inf), x > 0}
begin
  enx := sfc_en(n,x);
end;


{---------------------------------------------------------------------------}
function eibetax(n: integer; x: extended): extended;
  {-Return the exponential integral beta(n,x) = int(t^n*exp(-x*t), t=-1..1), n >= 0}
begin
  eibetax := sfc_eibeta(n,x);
end;


{---------------------------------------------------------------------------}
function geix(p,x: extended): extended;
  {-Return the generalized exponential integral E_p(x) = integral(exp(-x*t)/t^p, t=1..Inf), x >= 0}
begin
  geix := sfc_gei(p,x);
end;


{---------------------------------------------------------------------------}
function lix(x: extended): extended;
  {-Return the logarithmic integral li(x) = PV-integral(1/ln(t), t=0..x), x >= 0, x <> 1}
begin
  lix := sfc_li(x);
end;


{---------------------------------------------------------------------------}
function li_invx(x: extended): extended;
  {-Return the functional inverse of li(x), li(li_inv(x))=x}
begin
  li_invx := sfc_ali(x);
end;


{---------------------------------------------------------------------------}
function shix(x: extended): extended;
  {-Return the hyperbolic sine integral, shi(x) = integral(sinh(t)/t, t=0..x)}
begin
  shix := sfc_shi(x);
end;


{---------------------------------------------------------------------------}
function six(x: extended): extended;
  {-Return the sine integral, si(x) = integral(sin(t)/t, t=0..x)}
begin
  six := sfc_si(x);
end;


{---------------------------------------------------------------------------}
function ssix(x: extended): extended;
  {-Return the shifted sine integral, ssi(x) = six(x) - Pi/2}
begin
  ssix := sfc_ssi(x);
end;


{---------------------------------------------------------------------------}
{---------------------- Statistical distributions --------------------------}
{---------------------------------------------------------------------------}

function beta_pdfx(a, b, x: extended): extended;
  {-Return the probability density function of the beta distribution with}
  { parameters a and b: beta_pdf = x^(a-1)*(1-x)^(b-1) / beta(a,b)}
begin
  beta_pdfx := sfc_beta_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function beta_cdfx(a, b, x: extended): extended;
  {-Return the cumulative beta distribution function, a>0, b>0}
begin
  beta_cdfx := sfc_beta_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function beta_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the beta distribution function. a>0, b>0;}
  { 0 <= y <= 1. Given y the function finds x such that beta_cdf(a, b, x) = y}
begin
  beta_invx := sfc_beta_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function binomial_cdfx(p: extended; n, k: longint): extended;
  {-Return the cumulative binomial distribution function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}
begin
  binomial_cdfx := sfc_binomial_cdf(p,n,k);
end;


{---------------------------------------------------------------------------}
function binomial_pmfx(p: extended; n, k: longint): extended;
  {-Return the binomial distribution probability mass function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}
begin
  binomial_pmfx := sfc_binomial_pmf(p,n,k);
end;


{---------------------------------------------------------------------------}
function cauchy_pdfx(a, b, x: extended): extended;
  {-Return the Cauchy probability density function with location a }
  { and scale b > 0, 1/(Pi*b*(1+((x-a)/b)^2))}
begin
  cauchy_pdfx := sfc_cauchy_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function cauchy_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Cauchy distribution function with location a}
  { and scale b > 0, = 1/2 + arctan((x-a)/b)/Pi}
begin
  cauchy_cdfx := sfc_cauchy_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function cauchy_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Cauchy distribution function}
  { with location a and scale b > 0}
begin
  cauchy_invx := sfc_cauchy_inv(a, b, y);
end;


function chi_pdfx(nu: longint; x: extended): extended;
  {-Return the probability density function of the chi distribution, nu>0}
begin
  chi_pdfx := sfc_chi_pdf(nu,x);
end;


{---------------------------------------------------------------------------}
function chi_cdfx(nu: longint; x: extended): extended;
  {-Return the cumulative chi distribution with nu>0 degrees of freedom, x >= 0}
begin
  chi_cdfx := sfc_chi_cdf(nu,x);
end;


{---------------------------------------------------------------------------}
function chi_invx(nu: longint; p: extended): extended;
  {-Return the functional inverse of the chi distribution, nu>0, 0 <= p < 1}
begin
  chi_invx := sfc_chi_inv(nu,p);
end;


{---------------------------------------------------------------------------}
function chi2_pdfx(nu: longint; x: extended): extended;
  {-Return the probability density function of the chi-square distribution, nu>0}
begin
  chi2_pdfx := sfc_chi2_pdf(nu,x);
end;


{---------------------------------------------------------------------------}
function chi2_cdfx(nu: longint; x: extended): extended;
  {-Return the cumulative chi-square distribution with nu>0 degrees of freedom, x >= 0}
begin
  chi2_cdfx := sfc_chi2_cdf(nu,x);
end;


{---------------------------------------------------------------------------}
function chi2_invx(nu: longint; p: extended): extended;
  {-Return the functional inverse of the chi-square distribution, nu>0, 0 <= p < 1}
begin
  chi2_invx := sfc_chi2_inv(nu,p);
end;


{---------------------------------------------------------------------------}
function evt1_pdfx(a, b, x: extended): extended;
  {-Return the probability density function of the Extreme Value Type I distribution}
  { with location a and scale b > 0, result = exp(-(x-a)/b)/b * exp(-exp(-(x-a)/b)) }
begin
  evt1_pdfx := sfc_evt1_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function evt1_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Extreme Value Type I distribution function}
  { with location a and scale b > 0;  result = exp(-exp(-(x-a)/b)). }
begin
  evt1_cdfx := sfc_evt1_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function evt1_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Extreme Value Type I distribution}
  { function with location a and scale b > 0;  result = a - b*ln(ln(-y)). }
begin
  evt1_invx := sfc_evt1_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function exp_pdfx(a, alpha, x: extended): extended;
  {-Return the exponential probability density function with location a }
  { and rate alpha > 0, = alpha*exp(-alpha*(x-a)) if x >= a, 0 if x < a.}
begin
  exp_pdfx := sfc_exp_pdf(a, alpha, x);
end;


{---------------------------------------------------------------------------}
function exp_cdfx(a, alpha, x: extended): extended;
  {-Return the cumulative exponential distribution function with location a}
  { and rate alpha > 0, = 1 - exp(-alpha*(x-a)) if x >= a, 0 if x < a.}
begin
  exp_cdfx := sfc_exp_cdf(a, alpha, x);
end;


{---------------------------------------------------------------------------}
function exp_invx(a, alpha, y: extended): extended;
  {-Return the functional inverse of the exponential distribution function with}
  { location a and rate alpha > 0}
begin
  exp_invx := sfc_exp_inv(a, alpha, y);
end;


{---------------------------------------------------------------------------}
function f_pdfx(nu1, nu2: longint; x: extended): extended;
  {-Return the probability density function of the F distribution; x >= 0, nu1, nu2 > 0}
begin
  f_pdfx := sfc_f_pdf(nu1, nu2, x);
end;


{---------------------------------------------------------------------------}
function f_cdfx(nu1, nu2: longint; x: extended): extended;
  {-Return the cumulative F distribution function; x >= 0, nu1, nu2 > 0}
begin
  f_cdfx := sfc_f_cdf(nu1, nu2, x);
end;


{---------------------------------------------------------------------------}
function f_invx(nu1, nu2: longint; y: extended): extended;
  {-Return the functional inverse of the F distribution, nu1, nu2 > 0, 0 <= y <= 1}
begin
  f_invx := sfc_f_inv(nu1, nu2, y);
end;


{---------------------------------------------------------------------------}
function gamma_pdfx(a, b, x: extended): extended;
  {-Return the probability density function of a gamma distribution with shape}
  { a>0, scale b>0: gamma_pdf = x^(a-1)*exp(-x/b)/gamma(a)/b^a, x>0}
begin
  gamma_pdfx := sfc_gamma_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function gamma_cdfx(a, b, x: extended): extended;
  {-Return the cumulative gamma distribution function, shape a>0, scale b>0}
begin
  gamma_cdfx := sfc_gamma_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function gamma_invx(a, b, p: extended): extended;
  {-Return the functional inverse of the gamma distribution function, shape a>0,}
  { scale b>0, 0 <= p <= 1, i.e. finds x such that gamma_cdf(a, b, x) = p}
begin
  gamma_invx := sfc_gamma_inv(a, b, p);
end;


{---------------------------------------------------------------------------}
function hypergeo_pmfx(n1,n2,n,k: longint): extended;
  {-Return the hypergeometric distribution probability mass function; n,n1,n2 >= 0, n <= n1+n2;}
  { i.e. the probability that among n randomly chosen samples from a container}
  { with n1 type1 objects and n2 type2 objects are exactly k type1 objects:}
begin
  hypergeo_pmfx := sfc_hypergeo_pmf(n1,n2,n,k);
end;


{---------------------------------------------------------------------------}
function hypergeo_cdfx(n1,n2,n,k: longint): extended;
  {-Return the cumulative hypergeometric distribution function; n,n1,n2 >= 0, n <= n1+n2}
begin
  hypergeo_cdfx := sfc_hypergeo_cdf(n1,n2,n,k);
end;


{---------------------------------------------------------------------------}
function invgamma_pdfx(a, b, x: extended): extended;
  {-Return the probability density function of an inverse gamma distribution}
  { with shape a>0, scale b>0: result = (b/x)^a/x*exp(-b/x)/Gamma(a), x >= 0}
begin
  invgamma_pdfx := sfc_invgamma_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function invgamma_cdfx(a, b, x: extended): extended;
  {-Return the cumulative inverse gamma distribution function, shape a>0, scale}
  { b>0: result = Gamma(a,b/x)/Gamma(a) = Q(a,b/x) = igammaq(a,b/x), x >= 0}
begin
  invgamma_cdfx := sfc_invgamma_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function invgamma_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the inverse gamma distribution function, shape}
  { a>0, scale b>0, 0 <= y <= 1, i.e. find x such that invgamma_cdf(a, b, x) = y  }
begin
  invgamma_invx := sfc_invgamma_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function kolmogorov_cdfx(x: extended): extended;
  {-Return the limiting form for the cumulative Kolmogorov distribution function}
begin
  kolmogorov_cdfx := sfc_ks_cdf(x);
end;


{---------------------------------------------------------------------------}
function kolmogorov_invx(y: extended): extended;
  {-Return the functional inverse of the Kolmogorov distribution}
begin
  kolmogorov_invx := sfc_ks_inv(y);
end;


{---------------------------------------------------------------------------}
function kumaraswamy_pdfx(a, b, x: extended): extended;
  {-Return the Kumaraswamy probability density function with shape}
  { parameters a,b>0, 0<=x<=1; result = a*b*x^(a-1)*(1-x^a)^(b-1) }
begin
  kumaraswamy_pdfx := sfc_kumaraswamy_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function kumaraswamy_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Kumaraswamy distribution function with}
  { shape parameters a,b > 0, 0 <= x <= 1; result = 1-(1-x^a)^b}
begin
  kumaraswamy_cdfx := sfc_kumaraswamy_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function kumaraswamy_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Kumaraswamy distribution}
  { with shape parameters a,b > 0; result = [1-(1-y)^(1/b)]^(1/a)}
begin
  kumaraswamy_invx := sfc_kumaraswamy_inv(a,b,y)
end;


{---------------------------------------------------------------------------}
function laplace_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Laplace distribution function with location a and scale b > 0}
begin
  laplace_cdfx := sfc_laplace_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function laplace_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Laplace distribution with location a and scale b > 0}
begin
  laplace_invx := sfc_laplace_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function laplace_pdfx(a, b, x: extended): extended;
  {-Return the Laplace probability density function with location a}
  { and scale b > 0, result = exp(-abs(x-a)/b) / (2*b)}
begin
  laplace_pdfx := sfc_laplace_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function levy_pdfx(a, b, x: extended): extended;
  {-Return the Levy probability density function with}
  { location a and scale parameter b > 0}
begin
  levy_pdfx := sfc_levy_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function levy_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Levy distribution function with}
  { location a and scale parameter b > 0}
begin
  levy_cdfx := sfc_levy_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function levy_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Levy distribution}
  { with location a and scale parameter b > 0}
begin
  levy_invx := sfc_levy_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function logistic_pdfx(a, b, x: extended): extended;
  {-Return the logistic probability density function with location a}
  { and scale parameter b > 0, exp(-(x-a)/b)/b/(1+exp(-(x-a)/b))^2}
begin
  logistic_pdfx := sfc_logistic_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function logistic_cdfx(a, b, x: extended): extended;
  {-Return the cumulative logistic distribution function with}
  { location a and scale parameter b > 0}
begin
  logistic_cdfx := sfc_logistic_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function logistic_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the logistic distribution}
  { with location a and scale parameter b > 0}
begin
  logistic_invx := sfc_logistic_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function lognormal_pdfx(a, b, x: extended): extended;
  {-Return the log-normal probability density function with}
  { location a and scale parameter b > 0, zero for x <= 0.}
begin
  lognormal_pdfx := sfc_lognormal_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function lognormal_cdfx(a, b, x: extended): extended;
  {-Return the cumulative log-normal distribution function with}
  { location a and scale parameter b > 0, zero for x <= 0.}
begin
  lognormal_cdfx := sfc_lognormal_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function lognormal_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the log-normal distribution}
  { with location a and scale parameter b > 0, 0 < y < 1.}
begin
  lognormal_invx := sfc_lognormal_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function logseries_pmfx(a: extended; k: longint): extended;
  {-Return the logarithmic (series) probability mass function}
  { with shape 0 < a < 1, k > 0; result = -a^k/(k*ln(1-a))   }
begin
  logseries_pmfx := sfc_ls_pmf(a,k);
end;


{---------------------------------------------------------------------------}
function logseries_cdfx(a: extended; k: longint): extended;
  {-Return the cumulative logarithmic (series) distribution function with shape 0 < a < 1, k > 0}
begin
  logseries_cdfx := sfc_ls_cdf(a,k);
end;


{---------------------------------------------------------------------------}
function maxwell_pdfx(b, x: extended): extended;
  {-Return the Maxwell probability density function with scale b > 0, x >= 0}
begin
  maxwell_pdfx := sfc_maxwell_pdf(b,x);
end;


{---------------------------------------------------------------------------}
function maxwell_cdfx(b, x: extended): extended;
  {-Return the cumulative Maxwell distribution function with scale b > 0, x >= 0}
begin
  maxwell_cdfx := sfc_maxwell_cdf(b,x);
end;


{---------------------------------------------------------------------------}
function maxwell_invx(b, y: extended): extended;
  {-Return the functional inverse of the Maxwell distribution with scale b > 0}
begin
  maxwell_invx := sfc_maxwell_inv(b,y);
end;


{---------------------------------------------------------------------------}
function moyal_pdfx(a, b, x: extended): extended;
  {-Return the Moyal probability density function with}
  { location a and scale parameter b > 0}
begin
  moyal_pdfx := sfc_moyal_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function moyal_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Moyal distribution function with}
  { location a and scale parameter b > 0}
begin
  moyal_cdfx := sfc_moyal_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function moyal_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Moyal distribution}
  { with location a and scale parameter b > 0}
begin
  moyal_invx := sfc_moyal_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function nakagami_pdfx(m, w, x: extended): extended;
  {-Return the probability density function of a gamma distribution with}
  { shape m>0, spread b>0, x>=0: nakagami_pdf = 2x*gamma_pdf(m,w/m,x^2) }
begin
  nakagami_pdfx := sfc_nakagami_pdf(m, w, x);
end;


{---------------------------------------------------------------------------}
function nakagami_cdfx(m, w, x: extended): extended;
  {-Return the cumulative Nakagami distribution function, shape m>0, spread w>0}
begin
  nakagami_cdfx := sfc_nakagami_cdf(m, w, x);
end;


{---------------------------------------------------------------------------}
function nakagami_invx(m, w, p: extended): extended;
  {-Return the functional inverse of the Nakagami distribution function, shape m>0,}
  { spread w>0, 0 <= p <= 1, i.e. find x such that nakagami_cdf(m, w, x) = p}
begin
  nakagami_invx := sfc_nakagami_inv(m, w, p);
end;

{---------------------------------------------------------------------------}
function negbinom_cdfx(p,r: extended; k: longint): extended;
  {-Return the cumulative negative binomial distribution function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}
begin
  negbinom_cdfx := sfc_negbinom_cdf(p,r,k);
end;


{---------------------------------------------------------------------------}
function negbinom_pmfx(p,r: extended; k: longint): extended;
  {-Return the negative binomial distribution probability mass function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}
begin
  negbinom_pmfx := sfc_negbinom_pmf(p,r,k);
end;


{---------------------------------------------------------------------------}
function normal_pdfx(mu, sd, x: extended): extended;
  {-Return the normal (Gaussian) probability density function with mean mu}
  { and standard deviation sd>0, exp(-0.5*(x-mu)^2/sd^2) / sqrt(2*Pi*sd^2)}
begin
  normal_pdfx := sfc_normal_pdf(mu,sd,x);
end;


{---------------------------------------------------------------------------}
function normal_cdfx(mu, sd, x: extended): extended;
  {-Return the normal (Gaussian) distribution function}
  { with mean mu and standard deviation sd > 0        }
begin
  normal_cdfx := sfc_normal_cdf(mu,sd,x);
end;


{---------------------------------------------------------------------------}
function normal_invx(mu, sd, y: extended): extended;
  {-Return the functional inverse of the normal (Gaussian) distribution}
  { with mean mu and standard deviation sd > 0, 0 < y < 1.}
begin
  normal_invx := sfc_normal_inv(mu, sd, y);
end;


{---------------------------------------------------------------------------}
function normstd_pdfx(x: extended): extended;
  {-Return the std. normal probability density function exp(-x^2/2)/sqrt(2*Pi)}
begin
  normstd_pdfx := sfc_normstd_pdf(x);
end;


{---------------------------------------------------------------------------}
function normstd_cdfx(x: extended): extended;
  {-Return the standard normal distribution function}
begin
  normstd_cdfx := sfc_normstd_cdf(x);
end;


{---------------------------------------------------------------------------}
function normstd_invx(y: extended): extended;
  {-Return the inverse standard normal distribution function, 0 <= y <= 1.}
  { For x=normstd_inv(y) and y from [0,1], normstd_cdf(x) = y}
begin
  if y <= 0.0 then normstd_invx := -MaxExtended
  else if y >= 1.0 then normstd_invx := MaxExtended
  else normstd_invx := sfc_normstd_inv(y);
end;


{---------------------------------------------------------------------------}
function pareto_pdfx(k, a, x: extended): extended;
  {-Return the Pareto probability density function with minimum value k > 0}
  { and shape a, x >= a > 0, result = (a/x)*(k/x)^a}
begin
  pareto_pdfx := sfc_pareto_pdf(k,a,x);
end;


{---------------------------------------------------------------------------}
function pareto_cdfx(k, a, x: extended): extended;
  {-Return the cumulative Pareto distribution function minimum value k > 0}
  { and shape a, x >= a > 0, result = 1-(k/x)^a}
begin
  pareto_cdfx := sfc_pareto_cdf(k,a,x);
end;


{---------------------------------------------------------------------------}
function pareto_invx(k, a, y: extended): extended;
  {-Return the functional inverse of the Pareto distribution with minimum}
  { value k > 0 and shape a, x >= a > 0, result = k/(1-x)^(1/a)}
begin
  pareto_invx := sfc_pareto_inv(k,a,y);
end;


{---------------------------------------------------------------------------}
function poisson_cdfx(mu: extended; k: longint): extended;
  {-Return the cumulative Poisson distribution function with mean mu >= 0}
begin
  poisson_cdfx := sfc_poisson_cdf(mu,k);
end;


{---------------------------------------------------------------------------}
function poisson_pmfx(mu: extended; k: longint): extended;
  {-Return the Poisson distribution probability mass function with mean mu >= 0}
begin
  poisson_pmfx := sfc_poisson_pmf(mu,k);
end;

{---------------------------------------------------------------------------}
function rayleigh_pdfx(b, x: extended): extended;
  {-Return the Rayleigh probability density function with}
  { scale b > 0, x >= 0; result = x*exp(-0.5*(x/b)^2)/b^2}
begin
  rayleigh_pdfx := sfc_rayleigh_pdf(b,x);
end;


{---------------------------------------------------------------------------}
function rayleigh_cdfx(b, x: extended): extended;
  {-Return the cumulative Rayleigh distribution function with scale b > 0, x >= 0}
begin
  rayleigh_cdfx := sfc_rayleigh_cdf(b,x);
end;


{---------------------------------------------------------------------------}
function rayleigh_invx(b, y: extended): extended;
  {-Return the functional inverse of the Rayleigh distribution with scale b > 0}
begin
  rayleigh_invx := sfc_rayleigh_inv(b,y);
end;


{---------------------------------------------------------------------------}
function t_pdfx(nu: longint; x: extended): extended;
  {-Return the probability density function of Student's t distribution, nu>0}
begin
  t_pdfx := sfc_t_pdf(nu,x);
end;


{---------------------------------------------------------------------------}
function t_cdfx(nu: longint; t: extended): extended;
  {-Return the cumulative Student t distribution with nu>0 degrees of freedom}
begin
  t_cdfx := sfc_t_cdf(nu,t);
end;


{---------------------------------------------------------------------------}
function t_invx(nu: longint; p: extended): extended;
  {-Return the functional inverse of Student's t distribution, nu>0, 0 <= p <= 1}
begin
  t_invx := sfc_t_inv(nu,p);
end;


{---------------------------------------------------------------------------}
function triangular_pdfx(a, b, c, x: extended): extended;
  {-Return the triangular probability density function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}
begin
  triangular_pdfx := sfc_triangular_pdf(a,b,c,x);
end;


{---------------------------------------------------------------------------}
function triangular_cdfx(a, b, c, x: extended): extended;
  {-Return the cumulative triangular distribution function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}
begin
  triangular_cdfx := sfc_triangular_cdf(a,b,c,x);
end;


{---------------------------------------------------------------------------}
function triangular_invx(a, b, c, y: extended): extended;
  {-Return the functional inverse of the triangular distribution with}
  { lower limit a, upper limit b, mode c; a<b, a <= c <= b, 0 <= y <= 1}
begin
  triangular_invx := sfc_triangular_inv(a,b,c,y);
end;


{---------------------------------------------------------------------------}
function uniform_pdfx(a, b, x: extended): extended;
  {-Return the uniform probability density function on [a,b], a<b}
begin
  uniform_pdfx := sfc_uniform_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function uniform_cdfx(a, b, x: extended): extended;
  {-Return the cumulative uniform distribution function on [a,b], a<b}
begin
  uniform_cdfx := sfc_uniform_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function uniform_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the uniform distribution on [a,b], a<b}
begin
  uniform_invx := sfc_uniform_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function weibull_pdfx(a, b, x: extended): extended;
  {-Return the Weibull probability density function with shape a > 0}
  { and scale b > 0, result = a*x^(a-1)*exp(-(x/b)^a)/ b^a, x > 0}
begin
  weibull_pdfx := sfc_weibull_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function weibull_cdfx(a, b, x: extended): extended;
  {-Return the cumulative Weibull distribution function with}
  { shape parameter a > 0 and scale parameter b > 0}
begin
  weibull_cdfx := sfc_weibull_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function weibull_invx(a, b, y: extended): extended;
  {-Return the functional inverse of the Weibull distribution}
  { with shape parameter a > 0 and scale parameter b > 0}
begin
  weibull_invx := sfc_weibull_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function wald_pdfx(mu, b, x: extended): extended;
  {-Return the Wald (inverse Gaussian) probability density}
  { function with mean mu > 0, scale b > 0 for x >= 0     }
begin
  wald_pdfx := sfc_wald_pdf(mu, b, x);
end;


{---------------------------------------------------------------------------}
function wald_cdfx(mu, b, x: extended): extended;
  {-Return the Wald (inverse Gaussian) distribution  }
  { function with mean mu > 0, scale b > 0 for x >= 0}
begin
  wald_cdfx := sfc_wald_cdf(mu, b, x);
end;


{---------------------------------------------------------------------------}
function wald_invx(mu, b, y: extended): extended;
  {-Return the functional inverse of the Wald (inverse Gaussian)}
  { distribution with mean mu > 0, scale b > 0, 0 <= y < 1.     }
begin
  wald_invx := sfc_wald_inv(mu, b, y);
end;


{---------------------------------------------------------------------------}
function zipf_pmfx(r: extended; k: longint): extended;
  {-Return the Zipf distribution probability mass function k^(-(r+1))/zeta(r+1), r>0, k>0}
begin
  zipf_pmfx := sfc_zipf_pmf(r,k);
end;


{---------------------------------------------------------------------------}
function zipf_cdfx(r: extended; k: longint): extended;
  {-Return the cumulative Zipf distribution function H(k,r+1)/zeta(r+1), r>0, k>0}
begin
  zipf_cdfx := sfc_zipf_cdf(r,k);
end;


{---------------------------------------------------------------------------}
{------------------ Orthogonal polynomials and related ---------------------}
{---------------------------------------------------------------------------}

function chebyshev_tx(n: integer; x: extended): extended;
  {-Return Tn(x), the Chebyshev polynomial of the first kind, degree n}
begin
  chebyshev_tx := sfc_chebyshev_t(n,x);
end;


{---------------------------------------------------------------------------}
function chebyshev_ux(n: integer; x: extended): extended;
  {-Return Un(x), the Chebyshev polynomial of the second kind, degree n}
begin
  chebyshev_ux := sfc_chebyshev_u(n,x);
end;


{---------------------------------------------------------------------------}
function chebyshev_vx(n: integer; x: extended): extended;
  {-Return V_n(x), the Chebyshev polynomial of the third kind, degree n>=0}
begin
  chebyshev_vx := sfc_chebyshev_v(n,x);
end;


{---------------------------------------------------------------------------}
function chebyshev_wx(n: integer; x: extended): extended;
  {-Return W_n(x), the Chebyshev polynomial of the fourth kind, degree n>=0}
begin
  chebyshev_wx := sfc_chebyshev_w(n,x);
end;


{---------------------------------------------------------------------------}
function chebyshev_f1x(v,x: extended): extended;
  {-Return T_v(x), the Chebyshev function the first kind, real part for x<-1}
begin
  chebyshev_f1x := sfc_chebf1(v,x);
end;


{---------------------------------------------------------------------------}
function gegenbauer_cx(n: integer; a,x: extended): extended;
  {-Return Cn(a,x), the nth Gegenbauer (ultraspherical) polynomial with}
  { parameter a. The degree n must be non-negative; a should be > -0.5 }
  { When a = 0,   C0(0,x) = 1,  and   Cn(0,x) = 2/n*Tn(x)   for n <> 0.}
begin
  gegenbauer_cx := sfc_gegenbauer_c(n,a,x);
end;


{---------------------------------------------------------------------------}
function hermite_hx(n: integer; x: extended): extended;
  {-Return Hn(x), the nth Hermite polynomial, degree n >= 0}
begin
  hermite_hx := sfc_hermite_h(n,x);
end;


{---------------------------------------------------------------------------}
function hermite_hex(n: integer; x: extended): extended;
  {-Return He_n(x), the nth "probabilists'" Hermite polynomial, degree n >= 0}
begin
  hermite_hex := sfc_hermite_he(n,x);
end;


{---------------------------------------------------------------------------}
function jacobi_px(n: integer; a,b,x: extended): extended;
  {-Return Pn(a,b,x), the nth Jacobi polynomial with parameters a,b. Degree n}
  { must be >= 0; a,b should be > -1 (a+b must not be an integer < -1).}
begin
  jacobi_px := sfc_jacobi_p(n,a,b,x);
end;


{---------------------------------------------------------------------------}
function laguerrex(n: integer; a,x: extended): extended;
  {-Return Ln(a,x), the nth generalized Laguerre polynomial with parameter a;}
  { degree n must be >= 0. x >=0 and a > -1 are the standard ranges.}
begin
  laguerrex := sfc_laguerre(n,a,x);
end;


{---------------------------------------------------------------------------}
function laguerre_assx(n,m: integer; x: extended): extended;
  {-Return the associated Laguerre polynomial Ln(m,x); n,m >= 0}
begin
  laguerre_assx := sfc_laguerre_ass(n,m,x);
end;


{---------------------------------------------------------------------------}
function laguerre_lx(n: integer; x: extended): extended;
  {-Return the nth Laguerre polynomial Ln(0,x); n >= 0}
begin
  laguerre_lx := sfc_laguerre_l(n,x);
end;


{---------------------------------------------------------------------------}
function legendre_px(l: integer; x: extended): extended;
  {-Return P_l(x), the Legendre polynomial/function P_l, degree l}
begin
  legendre_px := sfc_legendre_p(l,x);
end;


{---------------------------------------------------------------------------}
function legendre_qx(l: integer; x: extended): extended;
  {-Return Q_l(x), the Legendre function of the 2nd kind, degree l >=0, |x| <> 1}
begin
 legendre_qx := sfc_legendre_q(l,x);
end;


{---------------------------------------------------------------------------}
function legendre_plmx(l,m: integer; x: extended): extended;
  {-Return the associated Legendre polynomial P_lm(x)}
begin
  legendre_plmx := sfc_legendre_plm(l,m,x);
end;


{---------------------------------------------------------------------------}
function legendre_qlmx(l,m: integer; x: extended): extended;
  {-Return Q(l,m,x), the associated Legendre function of the second kind; l >= 0, l+m >= 0, |x|<>1}
begin
  legendre_qlmx := sfc_legendre_qlm(l,m,x);
end;


{---------------------------------------------------------------------------}
procedure spherical_harmonicx(l, m: integer; theta, phi: extended; var yr,yi: extended);
  {-Return Re and Im of the spherical harmonic function Y_lm(theta,phi)}
begin
  sfc_spherical_harmonic(l,m,theta,phi,yr,yi);
end;


{---------------------------------------------------------------------------}
function toroidal_plmx(l,m: integer; x: extended): extended;
  {-Return the toroidal harmonic function P(l-0.5,m,x); l,m=0,1; x >= 1}
begin
  toroidal_plmx := sfc_thp(l,m,x);
end;


{---------------------------------------------------------------------------}
function toroidal_qlmx(l,m: integer; x: extended): extended;
  {-Return the toroidal harmonic function Q(l-0.5,m,x); l=0,1; x > 1}
begin
  toroidal_qlmx := sfc_thq(l,m,x);
end;


{---------------------------------------------------------------------------}
function zernike_rx(n,m: integer; r: extended): extended;
  {-Return the Zernike radial polynomial Rnm(r), r >= 0, n >= m >= 0, n-m even}
begin
  zernike_rx := sfc_zernike_r(n,m,r);
end;


{---------------------------------------------------------------------------}
{-------------------- Hypergeometric functions -----------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function hyperg_1F1x(a,b,x: extended): extended;
  {-Return the confluent hypergeometric function 1F1(a,b,x); Kummer's function M(a,b,x)}
begin
  hyperg_1F1x := sfc_1f1(a,b,x);
end;


{---------------------------------------------------------------------------}
function hyperg_1F1rx(a,b,x: extended): extended;
  {-Return the regularized Kummer hypergeometric function 1F1(a,b,x)/Gamma(b)}
begin
  hyperg_1F1rx := sfc_1f1r(a,b,x);
end;


{---------------------------------------------------------------------------}
function hyperg_ux(a,b,x: extended): extended;
  {-Return Tricomi's confluent hypergeometric function U(a,b,x). If}
  { x<0, then a must be an integer and a<0 or 1+a-b an integer < 0.}
begin
  hyperg_ux := sfc_chu(a,b,x);
end;


{---------------------------------------------------------------------------}
function hyperg_2F1x(a,b,c,x: extended): extended;
  {-Return the Gauss hypergeometric function 2F1(a,b;c;x)}
begin
  hyperg_2F1x := sfc_2f1(a,b,c,x);
end;


{---------------------------------------------------------------------------}
function hyperg_2F1rx(a,b,c,x: extended): extended;
  {-Return the regularized Gauss hypergeometric function 2F1(a,b,c,x)/Gamma(c)}
begin
  hyperg_2F1rx := sfc_2f1r(a,b,c,x);
end;


{---------------------------------------------------------------------------}
function hyperg_2F0x(a,b,x: extended): extended;
  {-Return 2F0(a,b,x), if x>0 then a or b must be a negative integer}
begin
  hyperg_2F0x := sfc_2f0(a,b,x);
end;


{---------------------------------------------------------------------------}
function WhittakerMx(k,m,x: extended): extended;
  {-Return the Whittaker M function = exp(-x/2)*x^(0.5+m) * 1F1(m-k-0.5,2m+1,x)}
begin
  WhittakerMx := sfc_whitm(k,m,x);
end;


{---------------------------------------------------------------------------}
function WhittakerWx(k,m,x: extended): extended;
  {-Return the Whittaker W function = exp(-x/2)*x^(0.5+m) * U(m-k-0.5,2m+1,x)}
begin
  WhittakerWx := sfc_whitw(k,m,x);
end;


{---------------------------------------------------------------------------}
function hyperg_0F1x(b,x: extended): extended;
  {-Return the confluent hypergeometric limit function 0F1(;b;x)}
begin
  hyperg_0F1x := sfc_0f1(b,x);
end;


{---------------------------------------------------------------------------}
function hyperg_0F1rx(b,x: extended): extended;
  {-Return the regularized confluent hypergeometric limit function 0F1(;b;x)/Gamma(b)}
begin
  hyperg_0F1rx := sfc_0f1r(b,x);
end;


{---------------------------------------------------------------------------}
function CylinderDx(v,x: extended): extended;
  {-Return Whittaker's parabolic cylinder function D_v(x)}
begin
  CylinderDx := sfc_pcfd(v,x);
end;


{---------------------------------------------------------------------------}
function CylinderUx(a,x: extended): extended;
  {-Return the parabolic cylinder function U(a,x)}
begin
  CylinderUx := sfc_pcfu(a,x);
end;


{---------------------------------------------------------------------------}
function CylinderVx(a,x: extended): extended;
  {-Return the parabolic cylinder function V(a,x) with 2a integer}
begin
  CylinderVx := sfc_pcfv(a,x);
end;


{---------------------------------------------------------------------------}
function HermiteHx(v,x: extended): extended;
  {-Return the Hermite function H_v(x) of degree v}
begin
  HermiteHx := sfc_pcfhh(v,x);
end;



{---------------------------------------------------------------------------}
{--------------------------- Other functions -------------------------------}
{---------------------------------------------------------------------------}

function bernoullix(n: integer): extended;
  {-Return the nth Bernoulli number, 0 if n<0 or odd n >= 3}
begin
  bernoullix := sfc_bernoulli(n);
end;


{---------------------------------------------------------------------------}
function bernpolyx(n: integer; x: extended): extended;
  {-Return the Bernoulli polynomial B_n(x), 0 <= n <= MaxBernoulli}
begin
  bernpolyx := sfc_bernpoly(n,x);
end;


{---------------------------------------------------------------------------}
function catalanx(x: extended): extended;
  {-Return the Catalan function C(x) = binomial(2x,x)/(x+1)}
begin
  catalanx := sfc_catf(x);
end;


{---------------------------------------------------------------------------}
function agmx(x,y: extended): extended;
  {-Return the arithmetic-geometric mean of |x| and |y|; |x|,|y| < sqrt(MaxExtended)}
begin
  agmx := sfc_agm(x,y);
end;


{---------------------------------------------------------------------------}
function besselpolyx(n: integer; x: extended): extended;
  {-Return yn(x), the nth Bessel polynomial}
begin
  besselpolyx := sfc_besspoly(n,x);
end;


{---------------------------------------------------------------------------}
function bringx(x: extended): extended;
  {-Return the Bring radical b := BR(x) with b^5 + b + x = 0}
begin
  bringx := sfc_br(x);
end;


{---------------------------------------------------------------------------}
function einsteinx(n: integer; x: extended): extended;
  {-Return the Einstein function E_n, n=1..4, x > 0 for n=3,4}
begin
  einsteinx := sfc_einstein(n,x);
end;


{---------------------------------------------------------------------------}
function debyex(n: integer; x: extended): extended;
  {-Return the Debye function D(n,x) = n/x^n*integral(t^n/(exp(t)-1),t=0..x) of order n>0, x>=0}
begin
  debyex := sfc_debye(n,x);
end;


{---------------------------------------------------------------------------}
function eulerx(n: integer): extended;
  {-Return the nth Euler number, 0 if n<0 or odd n}
begin
  eulerx := sfc_euler(n);
end;


{---------------------------------------------------------------------------}
function eulerpolyx(n: integer; x: extended): extended;
  {-Return the Euler polynomial E_n(x), 0 <= n < MaxBernoulli}
begin
  eulerpolyx := sfc_epoly(n,x);
end;


{---------------------------------------------------------------------------}
function euler_qx(q: extended): extended;
  {-Return the EulerQ function product(1-q^n, n=1..Inf)}
begin
  euler_qx := sfc_eulerq(q);
end;


{---------------------------------------------------------------------------}
function expnx(n: integer; x: extended): extended;
  {-Return the truncated exponential sum function e_n = sum(x^k/k!, k=0..n), 0 <= n < MAXGAM-1}
begin
  expnx := sfc_expn(n,x);
end;


{---------------------------------------------------------------------------}
function exprelnx(n: longint; x: extended): extended;
  {-Return the relative exponential = (e^x-sum(x^k/k!, k=0..n-1)*n!/x^n}
begin
  exprelnx := sfc_exprel(n,x);
end;


{---------------------------------------------------------------------------}
function fibfunx(v,x: extended): extended;
  {-Return the general Fibonacci function F_v(x)}
begin
  fibfunx := sfc_fibf(v,x);
end;


{---------------------------------------------------------------------------}
function fibpolyx(n: integer; x: extended): extended;
  {-Return the Fibonacci polynomial F_n(x)}
begin
  fibpolyx := sfc_fpoly(n,x);
end;


{---------------------------------------------------------------------------}
function lucpolyx(n: integer; x: extended): extended;
  {-Return the Lucas polynomial L_n(x)}
begin
  lucpolyx := sfc_lpoly(n,x);
end;


{---------------------------------------------------------------------------}
function keplerx(M,e: extended): extended;
  {-Solve Kepler's equation, result x is the eccentric anomaly from the mean anomaly M and the }
  { eccentricity e >= 0; x - e*sin(x) = M, x + x^3/3 = M, or e*sinh(x) - x = M for e <1, =1, >1}
begin
  keplerx := sfc_kepler(M, e);
end;


{---------------------------------------------------------------------------}
function LambertWx(x: extended): extended;
  {-Return the Lambert W function (principal branch), x >= -1/e}
begin
  LambertWx := sfc_LambertW(x);
end;


{---------------------------------------------------------------------------}
function LambertW1x(x: extended): extended;
  {-Return the Lambert W function (-1 branch), -1/e <= x < 0}
begin
  LambertW1x := sfc_LambertW1(x);
end;


{---------------------------------------------------------------------------}
function LangevinLx(x: extended): extended;
  {-Return the Langevin function L(x) = coth(x) - 1/x, L(0) = 0}
begin
  LangevinLx := sfc_langevin(x);
end;


{---------------------------------------------------------------------------}
function LangevinL_invx(x: extended): extended;
  {-Return the functional inverse of the Langevin function, |x| < 1}
begin
  LangevinL_invx := sfc_llinv(x);
end;


{---------------------------------------------------------------------------}
function MarcumQx(m: integer; a,b: extended): extended;
  {-Return the generalized Marcum Q function Q(m,a,b), a,b >= 0}
begin
  MarcumQx := sfc_mq(m,a,b);
end;


{---------------------------------------------------------------------------}
function omegax(x: extended): extended;
  {-Return the Wright omega function, i.e. the solution w of w + ln(w) = x}
begin
  omegax := sfc_wo(x);
end;


{---------------------------------------------------------------------------}
function RiemannRx(x: extended): extended;
  {-Return the Riemann prime counting function R(x), x >= 1/1024}
begin
  RiemannRx := sfc_ri(x);
end;


{---------------------------------------------------------------------------}
function RiemannR_invx(x: extended): extended;
  {-Return the functional inverse of R(x), R(RiemannR_inv(x))=x, x >= 1.125}
begin
  RiemannR_invx := sfc_ari(x);
end;


{---------------------------------------------------------------------------}
function rrcfx(q: extended): extended;
  {-Return the Rogers-Ramanujan continued fraction for |q| < 1}
begin
  rrcfx := sfc_rrcf(q);
end;


{---------------------------------------------------------------------------}
function cosintx(n: integer; x: extended): extended;
  {-Return cosint(n, x) = integral(cos(t)^n, t=0..x), n>=0}
begin
  cosintx := sfc_cosint(n,x);
end;


{---------------------------------------------------------------------------}
function sinintx(n: integer; x: extended): extended;
  {-Return sinint(n, x) = integral(sin(t)^n, t=0..x), n>=0}
begin
  sinintx := sfc_sinint(n,x);
end;


{---------------------------------------------------------------------------}
function transportx(n: integer; x: extended): extended;
  {-Return the transport integral J_n(x) for x >= 0, n >= 2,}
  { J_n(x) = integral(t^n*exp(t)/(exp(t)-1)^2, t=0..x)      }
begin
  transportx := sfc_trint(n,x);
end;


end.
