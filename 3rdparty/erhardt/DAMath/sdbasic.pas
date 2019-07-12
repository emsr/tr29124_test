unit sdBasic;

{Basic code for double precision special functions}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}

{$define damtools_intde} {Use DE integration routines from DAMTools}

uses
  DAMath;


(*************************************************************************

 DESCRIPTION   :  Basic code for double precision special functions:
                  Definitions, constants, etc.

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  define damtools_intde to use the intde functions from DAMTools

 REFERENCES    :  Index in damath_info.txt/references

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  05.02.13  W.Ehrhardt  sdBasic:  Initial BP7 version from AMath.sfbasic
 1.00.01  05.02.13  we          sdBasic:  sfd_bernoulli
 1.00.02  05.02.13  we          sdMisc:   sfd_agm, sfd_agm2
 1.00.03  06.02.13  we          sdEllInt: THexDblW constants for sfd_nome
 1.00.04  06.02.13  we          sdEllInt: adjusted some constants for double precision
 1.00.05  06.02.13  we          sdEllInt: changed test 0<x to x<-eps_d  in sfd_ellint_?
 1.00.06  06.02.13  we          sdEllInt: changed test abs(t-Pi_2)
 1.00.07  07.02.13  we          sdMisc:   sfd_LambertW, sfd_LambertW1
 1.00.08  07.02.13  we          sdEllInt: adjusted Carlson eps constants
 1.00.09  07.02.13  we          sdEllInt: sfd_EllipticPiC/CPi always with cel
 1.00.10  07.02.13  we          sdBasic:  sfd_bernoulli with bx16 (halved iteration count)
 1.00.11  07.02.13  we          sdGamma:  stirf, sfd_gaminv_small, sfd_gamma_medium
 1.00.12  08.02.13  we          sdGamma:  sfd_lngcorr, sfd_gstar, lngamma_small, lanczos
 1.00.13  08.02.13  we          sdbasic:  BoFHex as THexDblW
 1.00.14  08.02.13  we          sdGamma:  sfd_fac, sfd_dfac, sfd_binomial, sfd_psi, sfd_poch1
 1.00.15  08.02.13  we          sdZeta:   sfd_dilog, sfd_etaint, sfd_zetaint
 1.00.16  09.02.13  we          sdZeta:   zetap with 53-bit logic from [19]
 1.00.17  09.02.13  we          sdZeta:   etam1pos
 1.00.18  09.02.13  we          sdGamma:  sfd_trigamma/tetragamma/pentagamma/polygamma
 1.00.19  09.02.13  we          sdZeta:   sfd_cl2
 1.00.20  10.02.13  we          sdMisc:   sfd_debye
 1.00.21  10.02.13  we          sdPoly:   adjusted tol in lq_hf
 1.00.22  10.02.13  we          sdGamma:  sfd_lngamma uses FacTab
 1.00.23  10.02.13  we          sdMisc:   sfd_ri
 1.00.24  10.02.13  we          sdExpInt: sfd_ei
 1.00.25  10.02.13  we          sdExpInt: modified sfd_eiex to compute exp(x) if expx < 0
 1.00.26  10.02.13  we          sdExpInt: en_series avoids FPC nonsense
 1.00.27  10.02.13  we          sdExpInt: sfd_ssi, auxfg
 1.00.28  11.02.13  we          sdGamma:  igam_aux, sfd_igprefix
 1.00.29  11.02.13  we          sdErf:    erf_small, sfd_erfc_iusc, sfd_erfi
 1.00.30  11.02.13  we          sdErf:    sfd_fresnel
 1.00.31  12.02.13  we          sdSDist:  sfd_normstd_cdf, sfd_normstd_pdf
 1.00.32  12.02.13  we          sdZeta:   sfd_fdm05, sfd_fdp05
 1.00.33  12.02.13  we          sdBessel: Constant IvMaxX in sfd_iv
 1.00.34  12.02.13  we          sdBessel: sfd_i0
 1.00.35  12.02.13  we          sdBessel: sfd_j0, sfd_y0, bess_m0p0
 1.00.36  13.02.13  we          sdBessel: sfd_j1, sfd_y1, bess_m1p1
 1.00.37  13.02.13  we          sdBessel: sfd_in
 1.00.38  13.02.13  we          sdBessel: fix two near overflow cases in bessel_jy
 1.00.39  14.02.13  we          sdBessel: Airy functions
 1.00.40  14.02.13  we          sdBessel: Kelvin functions
 1.00.41  14.02.13  we          sdBessel: Struve functions
 1.00.42  15.02.13  we          sdGamma:  improved sfd_git
 1.00.43  16.02.13  we          sdZeta:   improved sfd_zetah
 1.00.44  20.02.13  we          SpecFunD: SpecFunD_Version
 1.00.45  01.03.13  we          (some)    Reduced some Chebyshev degrees
 1.00.46  02.02.13  we          sdBessel: Fixed value of IvMaxX in sfd_iv

 1.01.00  15.03.13  we          sdEllInt: Remaining Jacobi elliptic functions
 1.01.01  25.03.13  we          sdEllInt: Remaining inverse Jacobi elliptic functions
 1.01.02  28.03.13  we          sdZeta:   sfd_zetah for s<0
 1.01.03  28.03.13  we          sdZeta:   sfd_trilog
 1.01.04  30.03.13  we          sdZeta:   sfd_zetah for a near 1 and s<0
 1.01.05  30.03.13  we          sdGamma:  improved sfd_binomial

 1.02.00  07.04.13  we          sdSDist:  Improved sfd_chi2_pdf for very small x
 1.02.01  11.04.13  we          sdZeta:   sfd_fdp15
 1.02.02  13.04.13  we          sdSDist:  k changed to longint in sfd_poisson_pmf
 1.02.03  14.04.13  we          sdSDist:  sfd_rayleigh_pdf/cdf/inv
 1.02.04  18.04.13  we          sdGamma:  improved sfd_beta
 1.02.05  18.04.13  we          sdGamma:  ibeta_series with parameter normalised
 1.02.06  18.04.13  we          sdGamma:  sfd_nnbeta
 1.02.07  19.04.13  we          sdSDist:  sfd_maxwell_pdf/cdf/inv
 1.02.08  20.04.13  we          sdSDist:  sfd_evt1_pdf/cdf/inv
 1.02.09  20.04.13  we          sdSDist:  sfd_maxwell_cdf/inv with lower igamma functions
 1.02.10  21.04.13  we          sdGamma:  lanczos_gm05 = lanczos_g - 0.5 in implementation
 1.02.11  27.04.13  we          sdGamma:  polygam_negx for n=11,12
 1.02.12  01.05.13  we          SpecFunD: FresnelC/FresnelS

 1.03.00  09.05.13  we          sdBessel: Airy/Scorer functions sfd_airy_gi/hi
 1.03.01  11.05.13  we          sdGamma:  improved sfd_polygamma for large order: e.g. (3001, 871)
 1.03.02  13.05.13  we          sdPoly:   Change limits to 2e16 in p0ph,p1mh,p1ph
 1.03.03  27.05.13  we          sdHyperG: New Hypergeometric functions unit
 1.03.04  27.05.13  we          sdGamma:  sfd_nnbeta for a<=0 or b<=0
 1.03.05  31.05.13  we          sdHyperG: AMath changes for sfd_chu

 1.04.00  15.06.13  we          sdBessel: Fix typo in h1v_large
 1.04.01  15.06.13  we          sdBessel: sfd_bess_kv2
 1.04.02  15.06.13  we          sdHyperG: Temme's chu and uabx from AMath
 1.04.03  16.06.13  we          sdHyperG: Whittaker functions sfd_whitm/sfd _whitw
 1.04.04  28.06.13  we          sdBessel: Check overflow / return INF in bessel_jy
 1.04.05  28.06.13  we          sdHyperG: sfd_0f1
 1.04.06  01.07.13  we          sdHyperG: sfd_0f1r

 1.05.00  28.07.13  we          sdHyperG: MAXIT=64 in h1f1_tricomi
 1.05.01  28.07.13  we          sdGamma:  sfd_git for x<0
 1.05.02  15.08.13  we          sdSDist:  Removed code fragment for y>1 in sfd_beta_inv
 1.05.03  15.08.13  we          sdSDist:  Check a,b > 0 in sfd_beta_cdf
 1.05.04  15.08.13  we          sdSDist:  sfd_kumaraswamy_pdf/cdf/inv
 1.05.05  16.08.13  we          sdZeta:   sfd_zeta for small s
 1.05.06  16.08.13  we          sdZeta:   sfd_polylog with sfd_zetaint
 1.05.07  15.08.13  we          sdErf:    sfd_erf_z, sfd_erf_p, sfd_erf_q
 1.05.08  18.08.13  we          sdsDist:  normal/normstd_pdf/cdf with erf_z and erf_p

 1.06.00  07.09.13  we          sdZeta:   Improved sfd_zetah/hurwitz_formula
 1.06.01  12.09.13  we          sdZeta:   Bernoulli polynomials sfd_bernpoly
 1.06.02  18.09.13  we          sdPoly:   sfd_chebyshev_v/w
 1.06.03  24.09.13  we          sdMisc:   sfd_cosint, sfd_sinint
 1.06.04  25.09.13  we          sdGamma:  sfd_batemang
 1.06.05  25.09.13  we          (some)    use const one_d to avoid FPC nonsense
 1.06.06  26.09.13  we          sdMisc:   Fixed quasi-periodic code for sfd_cos/sinint
 1.06.07  28.09.13  we          sdEllInt: special_reduce_modpi, used in sfd_ellint_1/2/3 and sfd_hlambda

 1.07.00  19.10.13  we          sdsDist:  sfd_moyal_pdf/cdf/inv
 1.07.01  23.10.13  we          sdEllInt: sfd_EllipticKim, sfd_EllipticECim
 1.07.02  04.11.13  we          sdSDist:  improved sfd_t_pdf for small x^2/nu

 1.08.00  26.12.13  we          sdZeta:   sfd_fdp25

 1.09.00  24.03.14  we          sdHyperG: chu_luke only with $ifdef in dchu
 1.09.01  28.03.14  we          sdBessel: sfd_yn with LnPi from DAMath
 1.09.02  28.03.14  we          sdHyperG: InfOrNan tests in sfd_1f1

 1.10.00  18.04.14  we          sdGamma:  polygam_negx (larger n for x<0}
 1.10.01  02.05.14  we          sdZeta:   sfd_harmonic
 1.10.02  03.05.14  we          sdZeta:   sfd_harmonic2
 1.10.03  04.05.14  we          sdSDist:  sfd_zipf_pmf/cdf

 1.11.00  24.05.14  we          sdMisc:   sfd_ali functional inverse of sfd_li
 1.11.01  29.05.14  we          sdZeta:   Removed redundancy in sfd_harmonic2
 1.11.02  29.05.14  we          sdZeta:   polylogneg with array of single
 1.11.03  30.05.14  we          SpecFunD: Fix zipf_pmf/cdf function type
 1.11.04  31.05.14  we          sdZeta:   Lobachevski sfd_llci, sfd_llsi; improved harm2core
 1.11.05  04.06.14  we          sdMisc:   sfd_fpoly/sfd_lpoly
 1.11.06  05.06.14  we          sdSDist:  sfd_levy_pdf/cdf/inv

 1.12.00  16.06.14  we          Incomplete/inverse functions moved to sdGamma2
 1.12.01  16.06.14  we          sdGamma:  const lanczos_gm05: double
 1.12.02  16.06.14  we          sdGamma2: sfd_ibeta_inv (code from sdSDist)
 1.12.03  16.06.14  we          sdZeta:   Fermi/Dirac, Lobachewsky, harmonic functions moved to sdZeta2
 1.12.04  22.06.14  we          sdGamma2: sfd_ilng (inverse of lngamma)
 1.12.05  25.06.14  we          sdGamma2: sfd_ipsi (inverse of psi)
 1.12.06  03.07.14  we          SpecFunD: ibeta_inv
 1.12.07  13.07.14  we          sdHyperG: sfd_pcfd,sfd_pcfu,sfd_pcfhh

 1.13.00  30.07.14  we          sdBasic:  moebius[1..83]
 1.13.01  30.07.14  we          sdZeta:   Improved and expanded primezeta sfd_pz
 1.13.02  06.08.14  we          sdGamma:  First check IsNaND(x) in sfd_lngammas
 1.13.03  06.08.14  we          sdGamma:  improved sfd_gdr_pos for small x
 1.13.04  11.08.14  we          sdBessel: Iv := PosInf_v if Kv,Kv1=0 in bessel_ik or if x >= IvMaxXH in sfd_iv
 1.13.05  11.08.14  we          sdGamma:  special case x=x+y in sfd_beta
 1.13.06  11.08.14  we          sdHyperG: fix: sfd_pcfd,sfd_pcfu,sfd_pcfhh use double
 1.13.07  17.08.14  we          sdMisc:   Catalan function sfd_catf
 1.13.08  19.08.14  we          sdHyperG: sfd_pcfv
 1.13.09  20.08.14  we          sdGamma:  sfd_gamma_delta_ratio returns 1 if x=x+d

 1.14.00  01.09.14  we          sdHyperG: sfd_chu for x<0
 1.14.01  04.09.14  we          sdHyperG: fix h1f1_tricomi (division by t=0)
 1.14.02  07.09.14  we          sdSDist:  sfd_logistic_inv uses logit function
 1.14.03  06.10.14  we          sdSDist:  sfd_logistic_cdf uses logistic function

 1.15.00  26.10.14  we          sdHyperG: No Kummer transformation in 1F1 if b=0,-1,-2,...
 1.15.01  29.10.14  we          sdHyperG: 1F1 = 1 + a/b*(exp(x)-1) for very small a,b
 1.15.02  07.11.14  we          sdBasic:  Small Bernoulli numbers B2nHex moved to DAMath

 1.16.00  15.12.14  we          sdGamma:  small changes in sfd_batemang
 1.16.01  16.12.14  we          sdGamma2: Check IsNanOrInf in sfd_ilng
 1.16.02  23.12.14  we          sdSDist:  sfd_invgamma_pdf/cdf/inv
 1.16.03  26.12.14  we          sdSDist:  logarithmic (series) distribution
 1.16.04  03.01.15  we          sdSDist:  sfd_wald_cdf/pdf/inv
 1.16.05  05.01.15  we          sdSDist:  Error if k<1 in sfd_ls_cdf/pmf
 1.16.06  05.01.15  we          sdHyperG: Special case b=1 in gen_0f1
 1.16.07  09.01.15  we          sdMisc:   Euler numbers sfd_euler

 1.17.00  25.01.15  we          sdExpInt: sfd_ali from sdMisc
 1.17.01  25.01.15  we          sdExpInt: sfd_ei_inv
 1.17.02  28.03.15  we          sdEllInt: sfd_cel_d, sfd_ellint_d
 1.17.03  01.05.15  we          sdExpInt: fix sfd_ei_inv for x Nan

 1.18.00  18.05.15  we          sdZeta:   Improved sfd_etam1 (adjusted Borwein constants)
 1.18.01  18.05.15  we          sdZeta:   sfd_polylogr and sfd_lerch for s >= -1
 1.18.02  25.05.15  we          sdZeta:   sfd_lerch(1,s,a) for s<>1 and special case z=0
 1.18.03  04.06.15  we          sdBasic:  remove CSInitD
 1.18.04  04.06.15  we          sdMisc:   sfd_kepler
 1.18.05  10.06.15  we          sdBessel: new IJ_series replaces Jv_series
 1.18.06  10.06.15  we          sdBessel: rewrite of sfd_in: use IJ_series and CF1_I
 1.18.07  10.06.15  we          sdBessel: improved bess_m0p0/bess_m1p1 with rem_2pi_sym and Kahan summation
 1.18.08  10.06.15  we          sdBessel: sfd_struve_l/h, sfd_struve_l0/1 use sfd_struve_l
 1.18.09  14.06.15  we          sdBessel: sfd_struve_h/l with real nu >= 0
 1.18.10  19.06.15  we          sdBessel: scaled Airy functions sfd_airy_ais, sfd_airy_bis
 1.18.11  22.06.15  we          sdEllInt: avoid overflow in sfd_ell_rf, better arg checks in Carlson functions
 1.18.12  22.06.15  we          sdEllInt: Carlson sfd_ell_rg
 1.18.13  26.06.15  we          sdZeta2:  Fix: Type double in sfd_fdp25 and harm2core
 1.18.14  26.06.15  we          sdHyperG: Fix: double variables in chu_negx

 1.19.00  08.07.15  we          sdBessel: special cases Yv(+- 1/2, x)
 1.19.01  08.07.15  we          sdBessel: sfd_struve_h/l for v < 0
 1.19.02  08.07.15  we          (most):   removed NOBASM/BASM conditional defines
 1.19.03  09.07.15  we          sdBessel: fix sfd_struve_h/l(v,0) for v<=-1
 1.19.04  19.07.15  we          sdZeta:   eta near s=1
 1.19.05  19.07.15  we          sdEllint: sfd_el1(x, kc) = arctan(x) for kc^2=1
 1.19.06  01.08.15  we          sdHyperG: Avoid some overflows in h2f1_sum, improved h2f1_trans

 1.20.00  25.08.15  we          (most):   Use THREE & SIXX to avoid 'optimization'
 1.20.01  25.08.15  we          sdZeta:   sfd_bernpoly: avoid integer overflow and FPC optimization error
 1.20.02  26.08.15  we          sdBessel: Sqrt2 in kelvin functions (anti-'optimizing')

 1.21.00  26.09.15  we          sdEllInt: lemniscate functions sin_lemn, cos_lemn
 1.21.01  29.09.15  we          sdHyperG: sfd_pcfv(a,x) for a < 0
 1.21.02  11.10.15  we          sdGamma:  polygamma for x<0 uses Euler series for Pi*cot(Pi*x)

 1.22.00  14.02.16  we          sdErf:    Adjusted trivial ranges of erf/erfc/erfce to double

 1.23.00  12.05.16  we          sdEllInt: sfd_cel_b, sfd_ellint_b
 1.23.01  20.06.16  we          sdMisc:   Bring radical sfd_br
 1.23.02  25.06.16  we          sdMisc:   Wright omega sfd_wo
 1.23.03  26.06.16  we          sdMisc:   Nan/Inf handling for sfd_wo
 1.23.04  27.06.16  we          sdMisc:   Nan/Inf handling for sfd_br

 1.24.00  23.08.16  we          sdGamma2: sfd_invgam (inverse of gamma)
 1.24.01  28.09.16  we          sdMisc:   sfd_fibf
 1.24.02  30.09.16  we          sdPoly:   sfd_chebf1
 1.24.03  07.10.16  we          sdPoly:   sfd_hermite_he
 1.24.04  09.10.16  we          sdSDist:  Kolmogorov distribution sfd_ks_cdf/inv

 1.25.00  21.05.17  we          sdEllInt: sfd_EllipticK, sfd_EllipticEC for k>1
 1.25.01  25.05.17  we          sdMisc:   sfd_eulerq
 1.25.02  26.05.17  we          sdMisc:   Improved AE for sfd_eulerq if q ~ -1
 1.25.03  09.06.17  we          sdSDist:  sfd_cauchy_inv uses tanpi
 1.25.04  22.06.17  we          sdMisc:   sfd_detai

 1.26.00  07.07.17  we          sdErf:    new sfd_fresnel_fg, improved sfd_fresnel
 1.26.01  10.07.17  we          sdErf:    sfd_fresnel_fg for x<0
 1.26.02  17.07.17  we          sdBessel: sfd_blam
 1.26.03  19.07.17  we          sdErf:    sfd_erfh, sfd_erf2
 1.26.04  22.07.17  we          sdSDist:  Fix sfd_chi2_pdf for x=0
 1.26.05  22.07.17  we          sdSDist:  sfd_chi_pdf/cdf/inv

 1.27.00  17.08.17  we          sdExpint: sfd_eisx2
 1.27.01  17.08.17  we          sdErf:    sfd_gsi for x < 0

 1.28.00  02.12.17  we          sdErf:    sfd_erfi_inv
 1.28.01  02.12.17  we          sdErf:    sfd_expint3 with sfd_igammal
 1.28.02  02.12.17  we          (some):   Suppress warnings: Local variable does not seem to be initialized
 1.28.03  06.12.17  we          sdSDist:  range checks in sfd_gamma_inv
 1.28.04  06.12.17  we          sdSDist:  sfd_nakagami_pdf/cdf/inv
 1.28.05  07.12.17  we          sdMisc:   change extended to double in sfd_detai
 1.28.06  08.12.17  we          sdMisc:   sfd_mq (Marcum Q)
 1.28.07  12.12.17  we          sdExpint: separate functions e1_small, e1ex_large
 1.28.08  13.12.17  we          sdExpint: sfd_e1s(x) = exp(x)*E1(x)
 1.28.09  14.12.17  we          sdExpint: function ei_asymp, sfd_eisx2 with ei_asymp
 1.28.10  14.12.17  we          sdExpint: sfd_eis(x) = exp(-x)*Ei(x)
 1.28.11  14.12.17  we          sdExpint: e1_small, e1ex_large with adjusted rational approximations
 1.28.12  19.12.17  we          sdSDist:  sfd_ks_inv with modified regula falsi

 1.29.00  04.01.18  we          sdMisc:   Improved sfd_ri
 1.29.01  10.01.18  we          sdMisc:   sfd_ari (inverse of RiemannR)
 1.29.02  18.01.18  we          sdErf:    minor changes in sfd_erf_p/q
 1.29.03  18.01.18  we          sdErf:    sfd_owent
 1.29.04  23.01.18  we          sdErf:    sfd_erfce_inv
 1.29.05  30.01.18  we          sdMisc:   sfd_einstein
 1.29.06  02.02.18  we          sdMisc:   rewrite sfd_debye (always integrate if x>4}
 1.29.07  02.02.18  we          sdMisc:   sfd_trint
 1.29.08  04.02.18  we          sdZeta:   Rewrite polylog, returns real part for x>1,n>1
 1.29.09  04.02.18  we          sdZeta:   sfd_zetaint(1) = NaN
 1.29.10  05.02.18  we          sdZeta:   sfd_polylog(n,x) for x > 1
 1.29.11  06.02.18  we          sdbasic:  sfd_intde/i_p
 1.29.12  06.02.18  we          sdBessel: sfd_struve_h with sfd_intdei_p
 1.29.13  07.02.18  we          sdbasic:  conditional define amtools_intde
 1.29.14  07.02.18  we          sdZeta:   sfd_polylogr for x < -1
 1.29.15  09.02.18  we          sdZeta:   Merged with (old) unit sdZeta2
 1.29.16  09.02.18  we          sdZeta:   sfd_fdr (Fermi-Dirac real order)
 1.29.17  12.02.18  we          sdZeta:   sfd_polylogr(s,x)=x for large s
 1.29.18  14.02.18  we          (some)    fixed some parameter/return types

 1.30.00  21.02.18  we          sdMisc:   sfd_langevin: Langevin function from DAMath
 1.30.01  21.02.18  we          sdMisc:   sfd_llinv: Inverse Langevin function
 1.30.02  24.02.18  we          sdZeta:   sfd_lerch for z < -1
 1.30.03  25.02.18  we          sdZeta:   fix NaN check for r in sfd_harmonic2
 1.30.04  25.02.18  we          sdZeta:   fix overflow for LerchPhi interand
 1.30.05  26.02.18  we          sdZeta:   sfd_ti (inverse tangent integral of real order)
 1.30.06  01.03.18  we          sdZeta:   improved lerch integration
 1.30.07  04.03.18  we          sdZeta:   sfd_polylogr for 1 < x <= 256
 1.30.08  05.03.18  we          sdZeta:   sfd_lchi for |x| > 1
 1.30.09  05.03.18  we          sdZeta:   Bose-Einstein integral sfd_beint
 1.30.10  06.03.18  we          sdZeta:   sfd_polylogr x > 256
 1.30.11  11.03.18  we          sdMisc:   Rogers-Ramanujan continued fraction sfd_rrcf

 1.31.00  26.03.18  we          sdBessel: sfd_[i,j,k,y]0int
 1.31.01  30.03.18  we          sdEllInt: sfd_cel_d and sfd_cel_b for |k| >= 1
 1.31.02  31.03.18  we          sdEllInt: sfd_EllipticPiC for |k| > 1
 1.31.03  02.04.18  we          sdEllInt: sfd_EllipticPiCim
 1.31.04  05.04.18  we          sdEllInt: Fix sfd_EllipticPiC for k^2 >= nu > 1
 1.31.05  05.04.18  we          sdEllInt: Mathematica style complete elliptic integrals
 1.31.06  12.04.18  we          sdEllInt: Neville theta functions
 1.31.07  15.04.18  we          sdEllInt: sfd_arccl/sl

 1.32.00  19.04.18  we          (some):   Fixes for FPC311
 1.32.01  22.04.18  we          sdEllInt: sfd_el2 computed with Carlson
 1.32.02  24.04.18  we          sdEllInt: sfd_emlambda(i*y): elliptic modular function
 1.32.03  24.04.18  we          sdEllInt: sfd_detai from sdmisc
 1.32.04  28.04.18  we          sdEllInt: Weierstrass elliptic functions
 1.32.05  01.05.18  we          sdEllInt: sfd_wpg with JacobiCN, CompAuxWP, equireduce, wpdup
 1.32.06  03.05.18  we          sdEllInt: sfd_wpg_inv
 1.32.07  06.05.18  we          sdEllInt: sfd_KleinJ

 1.33.00  14.05.18  we          sdEllInt: Basic lemniscatic case: separate function for argument
 1.33.01  14.05.18  we          sdEllInt: sfd_wpe_der / sfd_wpg_der
 1.33.02  23.05.18  we          sdEllInt: sfd_M_EllipticF := phi for m=0
 1.33.03  23.05.18  we          sdEllInt: sfd_M_EllipticPi
 1.33.04  08.06.18  we          sdGamma:  sfd_lnbinomial

 1.34.00  30.07.18  we          sdMisc:   sfd_exprel
 1.34.01  03.08.18  we          sdHyperG: Fix tmax update in hyp2f0
 1.34.02  03.08.18  we          sdHyperG: sfd_2f0

 1.35.00  20.08.18  we          sdEllint: IsNanOrInf check sfd_EllipticPiCim
 1.35.01  20.08.18  we          sdErf:    IsNanOrInf check sfd_erfi_inv
 1.35.02  20.08.18  we          sdBasic:  Make MaxBernoulli available for TPH
 1.35.03  21.08.18  we          sdMisc:   sfd_epoly
 1.35.04  23.08.18  we          sdGamma:  Check NaNs in sfd_gamma
 1.35.04  23.08.18  we          specfund: Removed MAXGAM check in gamma and gamma1pm1
 1.35.05  23.08.18  we          sdEllint: Implemented max iteration check in some functions
 1.35.06  23.08.18  we          sdExpInt: NanOrInf check in sfd_en
 1.35.07  23.08.18  we          sdZeta:   max iteration check in liae_pos
 1.35.08  23.08.18  we          sdBessel: NanOrInf check in sfd_jn, sfd_in
 1.35.09  23.08.18  we          sdGamma2: NanOrInf check in sfd_incgamma_ex, sfd_ibeta, sfd_nnbeta
 1.35.10  23.08.18  we                    Move Marcum Q from sdMisc to sdErf
 1.35.11  27.08.18  we          sdMisc:   sfd_besspoly
 1.35.12  28.08.18  we          sdGamma2: sfd_igammap_der
 1.35.13  01.09.18  we          sdGamma2: sdGamma:  sfd_lnbg
 1.35.14  07.09.18  we          sdGamma2: improved sfd_igamma
 1.35.15  07.09.18  we          sdGamma2: sfd_igamman
 1.35.16  07.09.18  we          sdExpInt: sfd_en for n<0 and x>0
 1.35.17  07.09.18  we          sdExpInt: sfd_eibeta

 1.36.00  16.09.18  we          sdBessl2: New unit with Airy/Kelvin/Struve functions from sdBessel
 1.36.01  16.09.18  we          sdBessl2: sfd_coulcl
 1.36.02  17.09.18  we          sdBessl2: sfc_cshift (Coulomb phase shift)
 1.36.03  28.09.18  we          sdBessl2: sfd_coul_ffp, sfd_coul_ggp, sfd_coul_f
 1.36.04  29.09.18  we          sdGamma2: sfd_igamma with Laguerre expansion
 1.36.05  01.10.18  we          sdMisc:   sfd_expn
 1.36.06  02.10.18  we          sdGamma2: removed sfd_igamman, sfd_igamma for x<0, a<>-1,-2,...,
 1.36.07  02.10.18  we          sdGamma2: sfd_igammal for x < 0
 1.36.08  05.10.18  we          sdHyperG: threshold in h1f1_acf

 1.37.00  12.10.18  we          sdGamma:  sfd_psistar
 1.37.01  17.10.18  we          sdBessl2: sfd_synchf, sfd_synchg

 1.38.00  20.11.18  we          sdBessl2: Make skp4hx/sinkp4 global
 1.38.01  20.11.18  we          sdBessl2: Derivatives of Kelvin functions
 1.38.02  21.11.18  we          sdBessl2: Fix NaN results in sfd_berp/beip/kerp/keip

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

(*-------------------------------------------------------------------------
  This Pascal code uses material and ideas from open source and public
  domain libraries, see the file '3rdparty.ama' for the licenses.
---------------------------------------------------------------------------*)


const
  SD_BaseVersion = '1.38.02';

{$ifdef J_OPT}
var
{$else}
const
{$endif}
  RTE_NoConvergence: integer = 234;   {RTE for no convergence, set negative}
                                      {to continue (with inaccurate result)}
  RTE_ArgumentRange: integer = 235;   {RTE for argument(s) out of range, set}
                                      {negative to continue with NAN result.}
                                      {Tested if $R+ option is enabled}


procedure sincosPix2(x: double; var s,c: double);
  {-Return s=sin(Pi/2*x^2), c=cos(Pi/2*x^2); (s,c)=(0,1) for abs(x) >= 2^53}

function sfd_bernoulli(n: integer): double;
  {-Return the nth Bernoulli number, 0 if n<0 or odd n >= 3}


{#Z+}
{The next two functions are renamed copies from the DAMTools unit.}
{If $damtools_intde is defined the DAMtools function are called.  }
{#Z-}

procedure sfd_intde_p(f: TFuncDP; p: pointer; a, b, eps: double; var result, abserr: double;
                      var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over the finite interval (a,b)}
  { using Double Exponential (DE) transformation. Parameters:  }

  {#F}
  { f      - integrand f(x), must be analytic over (a,b)   }
  { p      - untyped pointer to parameters for the function}
  { a      - lower limit of integration                    }
  { b      - upper limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_d/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}


procedure sfd_intdei_p(f: TFuncDP; p: pointer; a, eps: double; var result, abserr: double;
                       var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has no oscillatory factor. Parameters:   }

  {#F}
  { f      - integrand f(x), must be analytic over (a,INF) }
  { p      - untyped pointer to parameters for the function}
  { a      - lower limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_d/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}

const
  MaxBernoulli = 259;  {Maximum index for Bernoulli numbers}

{#Z+}
{---------------------------------------------------------------------------}
{Moebius function æ(n) for small arguments}
const
  n_moeb = 83;
  moebius: array[1..n_moeb] of shortint = (
             1, -1, -1,  0, -1,  1, -1,  0,  0,  1,
            -1,  0, -1,  1,  1,  0, -1,  0, -1,  0,
             1,  1, -1,  0,  0,  1,  0,  0, -1, -1,
            -1,  0,  1,  1,  1,  0, -1,  1,  1,  0,
            -1, -1, -1,  0,  0,  1, -1,  0,  0,  0,
             1,  0, -1,  0,  1,  0,  1,  1, -1,  0,
            -1,  1,  0,  0,  1, -1, -1,  0,  1, -1,
            -1,  0, -1,  1,  0,  0,  1, -1, -1,  0,
             0,  1, -1);


{---------------------------------------------------------------------------}
{Bernoulli numbers. The B(2n), n = 0..MaxB2nSmall are in DAMath}

const
  NBoF  = 20;
  BoFHex: array[0..NBoF] of THexDblW = ({Bernoulli(2k+2)/(2k+2)! }
            ($5555,$5555,$5555,$3FB5),  {+8.3333333333333333336E-2}
            ($6C17,$16C1,$C16C,$BF56),  {-1.3888888888888888888E-3}
            ($1567,$BC01,$566A,$3F01),  {+3.3068783068783068783E-5}
            ($EF0B,$9334,$BD77,$BEAB),  {-8.2671957671957671957E-7}
            ($0EBE,$2BF7,$6A8F,$3E56),  {+2.0876756987868098980E-8}
            ($267F,$D644,$2805,$BE02),  {-5.2841901386874931847E-10}
            ($9162,$C4E0,$6DB2,$3DAD),  {+1.3382536530684678833E-11}
            ($955C,$1F79,$DA4E,$BD57),  {-3.3896802963225828668E-13}
            ($2E9E,$1D65,$5587,$3D03),  {+8.5860620562778445639E-15}
            ($ACF1,$68CA,$57D9,$BCAF),  {-2.1748686985580618731E-16}
            ($376F,$F09C,$67E1,$3C59),  {+5.5090028283602295151E-18}
            ($2B5C,$033A,$97D9,$BC04),  {-1.3954464685812523341E-19}
            ($AD06,$D7C6,$B132,$3BB0),  {+3.5347070396294674718E-21}
            ($1C16,$D59F,$0F72,$BB5B),  {-8.9535174270375468504E-23}
            ($A26D,$A4CC,$EF2D,$3B05),  {+2.2679524523376830603E-24}
            ($E38B,$F96D,$C77D,$BAB1),  {-5.7447906688722024451E-26}
            ($1B62,$DE52,$D299,$3A5C),  {+1.4551724756148649018E-27}
            ($74A7,$6565,$5CDE,$BA07),  {-3.6859949406653101781E-29}
            ($4ADF,$DB3B,$EFE8,$39B2),  {+9.3367342570950446721E-31}
            ($61FF,$9047,$B322,$B95E),  {-2.3650224157006299346E-32}
            ($8464,$F932,$E25F,$3908)); {+5.9906717624821343044E-34}
{#Z-}

{$ifdef debug}
type
  sfd_debug_str = string[255];

const
  NDCTR = 8;

var
  sfd_diagctr: array[0..NDCTR-1] of longint; {General counters for diagnostics}
  sfd_debug_output: boolean;                 {Really doing the outputs if true}

procedure sfd_dump_diagctr;
  {-Writeln diagnostic counters}

procedure sfd_write_debug_str(const msg: sfd_debug_str);
  {-Writeln or Outputdebugstr of msg}
{$endif}


implementation

{$undef non_empty}
{$undef USE_OutputDebugString}

{$ifdef damtools_intde}
  {$define non_empty}
{$endif}

{$ifdef Debug}
  {$ifdef WIN32or64}
    {$ifndef VirtualPascal}
      {$define USE_OutputDebugString}
      {$define non_empty}
    {$endif}
  {$endif}
{$endif}

{$ifdef non_empty}
  uses
  {$ifdef USE_OutputDebugString}
    {$ifdef damtools_intde}
      DAMTools,
    {$endif}
    {$ifdef UNIT_SCOPE}
      winapi.windows;
    {$else}
      windows;
    {$endif}
  {$else}
    {$ifdef damtools_intde}
      DAMTools;
    {$endif}
  {$endif}
{$endif}



{---------------------------------------------------------------------------}
procedure sincosPix2(x: double; var s,c: double);
  {-Return s=sin(Pi/2*x^2), c=cos(Pi/2*x^2); (s,c)=(0,1) for abs(x) >= 2^53}
var
  n,f,g: double;
begin
  {Note that this routine is very sensible to argument changes!     }

  {Demo of sincosPix2 sensibility to argument changes               }
  {Computing sin(Pi/2*x^2) with sincosPix2 and MPArith for x=10000.1}
  {mp_float default bit precision = 240,  decimal precision = 72.2  }
  {Machine eps for double = 2.22044604925E-0016                     }
  {d = x(double)      = +10000.100000000000364                      }
  {x - d              = -3.6379788070917129517E-13                  }
  {rel. diff          = -3.6379424276674362773E-17                  }
  {a = sincosPix2     =  1.57073287395724723E-0002                  }
  {b = sin(Pi/2*x^2)  = +1.5707317311820675753E-2                   }
  {c = sin(Pi/2*d^2)  = +1.5707328739572472151E-2                   }
  {(c - b)            = +1.1427751796397653663E-8                   }
  {(c - b)/b          = +7.2754319337507758102E-7                   }
  {(c - a)            = -1.2497147346576779516E-19                  }
  {(c - a)/b          = -7.9562582829926944803E-18                  }

  if THexDblW(x)[3] and $7FF0 >= $4340 then begin
    {abs(x) >= 2^53}
    s := 0.0;
    c := 1.0;
  end
  else begin
    {c, s depend on frac(x) and int(x) mod 4. This code is based on  }
    {W. Van Snyder: Remark on algorithm 723: Fresnel integrals, 1993.}
    {http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.101.7180}
    {Fortran source from http://netlib.org/toms/723}
    f := abs(x);
    n := int(f);
    f := f - n;
    g := 2.0*frac(f*int(0.5*n));
    if frac(0.5*n)=0.0 then sincosPi(0.5*f*f + g, s, c)
    else begin
      sincosPi((0.5*f*f + f) + g, c, s);
      c := -c;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_bernoulli(n: integer): double;
  {-Return the nth Bernoulli number, 0 if n<0 or odd n >= 3}
var
  bn,p4: double;
  m: integer;
const
  bx16: array[0..8] of THexDblw = (    {Bernoulli(16*n+128)}
          ($C7AA,$9397,$49D6,$D78B),   {-5.2500923086774131347E+113}
          ($DC82,$055D,$312B,$DBFC),   {-1.2806926804084747554E+135}
          ($F740,$CDB0,$7C97,$E095),   {-1.8437723552033869952E+157}
          ($C6C9,$5F0A,$3B07,$E554),   {-1.3116736213556958091E+180}
          ($6528,$135D,$5994,$EA34),   {-3.9876744968232205368E+203}
          ($ED60,$24DA,$6069,$EF33),   {-4.5902296220617917542E+227}
          ($3B1F,$D638,$879E,$F44F),   {-1.8059559586909309354E+252}
          ($2A28,$6304,$1403,$F984),   {-2.2244891682179835734E+277}
          ($3856,$7CD9,$8C92,$FED2));  {-7.9502125045885251673E+302}
begin
   if odd(n) or (n<0) then begin
     if n=1 then sfd_bernoulli := -0.5
     else sfd_bernoulli := 0.0;
   end
   else begin
     m := n div 2;
     if m<=MaxB2nSmall then sfd_bernoulli := double(B2nHex[m])
     else if n>MaxBernoulli then sfd_bernoulli := PosInf_d
     else begin
       {When n is even, B(2n) = -2*(-1)^n*m!/(2Pi)^m*zeta(m) with m=2n. For }
       {large m (e.g. m>63) zeta(m) is very close to 1 and we can derive the}
       {asymptotic recursion formula B(m+1) = -m*(m+1)/(2Pi)^2 * B(m). The  }
       {avg. iteration count is <2, the max.rel. error < 1.9*eps_d for n=100}
       m  := (n - 120) div 16;
       bn := double(bx16[m]);
       m  := 16*m + 128;
       p4 := 4.0*PiSqr;
       if n>m then begin
         while n>m do begin
           inc(m,2);
           bn := bn/p4*m*(1-m);
         end;
       end
       else begin
         while m>n do begin
           bn := bn/m/(1-m)*p4;
           dec(m,2);
         end;
       end;
       sfd_bernoulli := bn;
     end;
  end;
end;


{$ifdef damtools_intde}

{---------------------------------------------------------------------------}
procedure sfd_intde_p(f: TFuncDP; p: pointer; a, b, eps: double; var result, abserr: double;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over the finite interval (a,b)}
  { using Double Exponential (DE) transformation. Parameters:  }
begin
  DAMTools.intde_p(f, p, a, b, eps, result, abserr, neval, ier);
end;


{---------------------------------------------------------------------------}
procedure sfd_intdei_p(f: TFuncDP; p: pointer; a, eps: double; var result, abserr: double;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has no oscillatory factor. Parameters:   }
begin
  DAMTools.intdei_p(f, p, a, eps, result, abserr, neval, ier);
end;

{$else}

{---------------------------------------------------------------------------}
procedure sfd_intde_p(f: TFuncDP; p: pointer; a, b, eps: double; var result, abserr: double;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over the finite interval (a,b)}
  { using Double Exponential (DE) transformation. Parameters:  }

  {#F}
  { f      - integrand f(x), must be analytic over (a,b)   }
  { p      - untyped pointer to parameters for the function}
  { a      - lower limit of integration                    }
  { b      - upper limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_d/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}
const
  mmax = 256;
  efs  = 0.1;
  hoff = 8.5;
var
  m: integer;
var
  epsln,epsh,h0,ehp,ehm,epst,ba,ir,iv,h,iback,
  irback,t,ep,em,xw,xa,wg,fa,fb,err,errt,errh,errd: double;
begin

  result := 0;
  abserr := 0;
  neval  := 0;

  if eps < eps_d/32 then begin
    ier := 1;
    exit;
  end
  else begin
    ier := 0;
    if a=b then exit;
  end;

  epsln := 1.0 - ln(efs*eps);
  epsh  := sqrt(efs*eps);
  h0    := hoff/epsln;
  ehp   := exp(h0);
  ehm   := 1.0/ehp;
  epst  := exp(-ehm*epsln);
  ba    := b - a;
  ir    := f(0.5*(a+b),p);
  ir    := 0.25*ir*ba;
  neval := 1;
  iv    := ir*Pi;
  err   := abs(iv)*epst;
  errh  := 0.0;

  h := 2.0*h0;
  m := 1;
  repeat
    iback  := iv;
    irback := ir;
    t := 0.5*h;
    repeat
      em := exp(t);
      ep := pi_2*em;
      em := pi_2/em;
      repeat
        xw := 1.0/(1.0 + exp(ep-em));
        xa := ba*xw;
        wg := xa*(1.0-xw);
        fa := a+xa;
        fb := b-xa;
        fa := f(fa,p)*wg;
        fb := f(fb,p)*wg;
        inc(neval,2);
        ir := ir + (fa+fb);
        iv := iv + (fa+fb)*(ep+em);
        errt:= (abs(fa) + abs(fb))*(ep+em);
        if m=1 then err := err + errt*epst;
        ep := ep*ehp;
        em := em*ehm;
      until (errt <= err) and (xw <= epsh);
      t := t + h;
    until t >= h0;
    if m=1 then begin
      errh := (err/epst)*epsh*h0;
      errd := 1.0 + 2.0*errh;
    end
    else errd := h*(abs(iv - 2.0*iback) + 4.0*abs(ir - 2.0*irback));
    h := 0.5*h;
    m := m+m;
  until (errd <=errh) or (m > mmax);

  result := h*iv;
  if errd > errh then begin
    ier := 3;
    abserr := errd*m;
  end
  else abserr := 0.5*errh*epsh*m/efs;
end;


{---------------------------------------------------------------------------}
procedure sfd_intdei_p(f: TFuncDP; p: pointer; a, eps: double; var result, abserr: double;
                var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has no oscillatory factor. Parameters:   }

  {#F}
  { f      - integrand f(x), must be analytic over (a,INF) }
  { p      - untyped pointer to parameters for the function}
  { a      - lower limit of integration                    }
  { eps    - relative error requested                      }
  { result - approximation to the integral,  if ier=0 or 3 }
  { abserr - estimate of the absolute error, if ier=0 or 3 }
  { neval  - number of integrand evaluations               }
  { ier    - 0: OK, 1: eps < eps_d/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}
const
  mmax = 256;
  efs  = 0.1;
  hoff = 11.0;
var
  m: integer;
var
  epsln,epsh,h0,ehp,ehm,epst,ir,iv,h,iback,irback,
  t,ep,em,xp,xm,fp,fm,err,errt,errh,errd: double;
begin

  result := 0.0;
  abserr := 0.0;
  neval  := 0;

  if eps < eps_d/32 then begin
    ier := 1;
    exit;
  end
  else ier := 0;

  epsln := 1.0 - ln(efs*eps);
  epsh  := sqrt(efs*eps);
  h0    := hoff/epsln;
  ehp   := exp(h0);
  ehm   := 1.0/ehp;
  epst  := exp(-ehm*epsln);
  ir    := f(a+1.0,p);
  neval := 1;
  iv    := ir*pi_2;
  err   := abs(iv)*epst;
  errh  := 0.0;

  h := 2.0*h0;
  m := 1;
  repeat
    iback  := iv;
    irback := ir;
    t := 0.5*h;
    repeat
      em := exp(t);
      ep := pi_4*em;
      em := pi_4/em;
      repeat
        xp  := exp(ep-em);
        xm  := 1.0/xp;
        inc(neval,2);
        fp  := f(a+xp,p);
        fm  := f(a+xm,p);
        fp  := fp*xp;
        fm  := fm*xm;
        ir  := ir + (fp+fm);
        iv  := iv + (fp+fm)*(ep+em);
        errt:= (abs(fp) + abs(fm))*(ep+em);
        if m=1 then err := err + errt*epst;
        ep  := ep*ehp;
        em  := em*ehm;
      until (errt <= err) and (xm <= epsh);
      t := t + h;
    until t >= h0;
    if m=1 then begin
      errh := (err/epst)*epsh*h0;
      errd := 1.0 + 2.0*errh;
    end
    else errd := h*(abs(iv - 2.0*iback) + 4.0*abs(ir - 2.0*irback));
    h := 0.5*h;
    m := m+m;
  until (errd <= errh) or (m >= mmax);

  result := h*iv;
  if errd > errh then begin
    ier := 3;
    abserr := errd*m;
  end
  else abserr := 0.5*errh*epsh*m/efs;

end;

{$endif}


{$ifdef debug}

{$ifdef USE_OutputDebugString}
{---------------------------------------------------------------------------}
procedure sfd_write_debug_str(const msg: sfd_debug_str);
  {-Writeln or Outputdebugstr of msg}
var
  ax: ansistring;
begin
  if sfd_debug_output then begin
    if IsConsole then writeln(msg)
    else begin
      ax := msg;
      OutputDebugString(pchar({$ifdef D12Plus}string{$endif}(ax)));
    end;
  end;
end;
{$else}
{---------------------------------------------------------------------------}
procedure sfd_write_debug_str(const msg: sfd_debug_str);
  {-Writeln or Outputdebugstr of msg}
begin
  if sfd_debug_output then writeln(msg);
end;
{$endif}


{---------------------------------------------------------------------------}
procedure sfd_dump_diagctr;
  {-Writeln diagnostic counters}
var
  k: integer;
begin
  write('Diag ctr:');
  for k:=0 to NDCTR-1 do write(k:3,':',sfd_diagctr[k]);
  writeln;
end;

begin
  sfd_debug_output := false;
  fillchar(sfd_diagctr, sizeof(sfd_diagctr),0);
{$endif}

end.
