unit sfBasic;

{Basic common code for special functions (extended/double precision)}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}

{$define amtools_intde} {Use DE integration routines from AMTools}


uses
  AMath;


(*************************************************************************

 DESCRIPTION   :  Basic common code for special functions (extended/double precision)
                  Definitions, constants, etc.

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  define amtools_intde to use the intde functions from AMTools

 REFERENCES    :  Index in amath_info.txt/references


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     29.12.09  W.Ehrhardt  Initial version: gamma
 0.11     30.12.09  we          lngamma, signgamma
 0.12     31.12.09  we          lngamma2, fix lngamma
 0.13     01.01.10  we          agm, agm2
 0.14     01.01.10  we          complementary elliptic integrals
 0.15     01.01.10  we          elliptic integrals
 0.16     01.01.10  we          cel, cel1, cel2
 0.17     02.01.10  we          ell_rc
 0.18     03.01.10  we          el1, el2
 0.19     03.01.10  we          ell_rf
 0.20     04.01.10  we          ell_rd
 0.21     04.01.10  we          ell_rj, el3
 0.22     04.01.10  we          Range EllipticF extended to |z|<=1 |k*z|<1
 0.23     04.01.10  we          EllipticE, EllipticEC, for |k|=1
 0.24     04.01.10  we          el1,el2 for kc=0
 0.25     05.01.10  we          Jacobi elliptic functions sn,cn,dn
 0.26     06.01.10  we          erf, erfc
 0.27     06.01.10  we          normstd_pdf, normstd_cdf, normstd_inv
 0.28     08.01.10  we          psi(digamma)
 0.29     09.01.10  we          beta, lnbeta
 0.30     10.01.10  we          rgamma

 0.31     13.01.10  we          unit SFCommon
 0.32     17.01.10  we          removed lngamma2
 0.33     21.01.10  we          adjusted eps for ell_rf (^6 scaling is not optimal)
 0.34     22.01.10  we          use absolute values of arguments in agm/2
 0.35     22.01.10  we          Asymptotic values for small/large args in Maple style
                                complete elliptic integrals of 1st and 2nd kind
 0.36     23.01.10  we          sfc_signgamma returns extended
 0.37     23.01.10  we          removed TP5 code
 0.38     24.01.10  we          sncndn: use sincos, handle case x=0
 0.39     24.01.10  we          Inf handling in sfc_normstd_cdf
 0.40     27.01.10  we          sfc_gamma

 0.41     10.02.10  we          removed Carlson debug code
 0.42     10.02.10  we          sfc_zeta, sfc_zeta1p
 0.43     14.02.10  we          sfc_e1
 0.44     15.02.10  we          sfc_ei, sfc_li
 0.45     16.02.10  we          improved sfc_li(x) for x>exp(6) via sfc_eiex
 0.46     17.02.10  we          sfc_en
 0.47     18.02.10  we          en_cfrac from [7]
 0.48     22.02.10  we          sfc_gamma: Pi/MaxExtended instead of Pi/MaxDouble

 0.49     01.03.10  we          e1(x) = -ei(-x) for x < 0
 0.50     02.03.10  we          sfc_shi
 0.51     03.03.10  we          sfc_ci, sfc_si, sfc_ssi
 0.52     04.03.10  we          sfc_chi
 0.53     06.03.10  we          sfc_dawson
 0.54     11.03.10  we          sfc_erf_inv, sfc_erfc_inv
 0.55     11.03.10  we          CSInitX
 0.56     16.03.10  we          sfc_dilog
 0.57     16.03.10  we          Pi^2/6 as hex for VER5X

 0.58     01.04.10  we          sfc_j0, sfc_y0
 0.59     02.04.10  we          sfc_j1, sfc_y1
 0.60     03.04.10  we          sfc_jn, cfrac_jn
 0.61     04.04.10  we          sfc_yn
 0.62     04.04.10  we          sfc_i0, sfc_i0e
 0.63     04.04.10  we          sfc_i1, sfc_i1e
 0.64     05.04.10  we          sfc_k0, sfc_k0e
 0.65     05.04.10  we          sfc_k1, sfc_k1e
 0.66     07.04.10  we          sfc_cin, sfc_cinh
 0.67     08.04.10  we          fix auxfg for x>sqrt(MaxExtended)
 0.68     10.04.10  we          sfc_fac
 0.69     10.04.10  we          sfc_in, sfc_kn
 0.70     12.04.10  we          sfc_ti2
 0.71     17.04.10  we          improved sfc_lngamma
 0.72     17.04.10  we          sfc_gamma1pm1
 0.73     22.04.10  we          sfc_cl2

 0.74     03.05.10  we          sincosPix2
 0.75     04.05.10  we          sfc_fresnel
 0.76     18.05.10  we          sfc_LambertW, sfc_LambertW1
 0.77     19.05.10  we          Fix sfc_LambertW1 for x denormal
 0.78     19.05.10  we          Typed const RTE_NoConvergence
 0.79     19.05.10  we          Fix for non-BASM in sfc_LambertW1
 0.80     20.05.10  we          Better decomposition em1h/em1l in sfc_LambertW
 0.81     21.05.10  we          Improved accuracy in sfc_normstd_cdf

 0.82     10.07.10  we          Typed const RTE_ArgumentRange
 0.83     10.07.10  we          sfc_ibeta
 0.84     10.07.10  we          sfc_beta_pdf
 0.85     10.07.10  we          sfc_beta_inv, sfc_beta_cdf
 0.86     13.07.10  we          sfc_t_pdf, sfc_t_cdf, sfc_t_inv
 0.87     16.07.10  we          sfc_f_pdf, sfc_f_cdf, sfc_f_inv
 0.88     18.07.10  we          Argument check $ifopt R+

 0.89     17.08.10  we          SFCommon split into several units

 1.00.00  17.08.10  we          Common version number after split
 1.00.01  17.08.10  we          sfGamma: Temme's gstar function
 1.00.02  23.08.10  we          sfGamma: Incomplete Gamma functions
 1.00.03  25.08.10  we          sfc_incgamma: cases when a is integer or half integer
 1.00.04  29.08.10  we          sfc_incgamma: dax calculated via Lanczos sum
 1.00.05  01.09.10  we          sfc_incgamma_ex, sfc_igprefix
 1.00.06  04.09.10  we          sfc_incgamma_inv
 1.00.07  05.09.10  we          sfc_incgamma_inv: use eps_d for p+q=1 check
 1.00.08  05.09.10  we          sfc_incgamma_inv: return Inf if p=1 and q=0
 1.00.09  05.09.10  we          sfSDist: sfc_chi2_pdf/cdf/inv
 1.00.10  07.09.10  we          sfSDist: sfc_gamma_pdf/cdf/inv
 1.00.11  07.09.10  we          avoid infinite loop if x is NAN in sfc_ei/sfc_e1
 1.00.12  08.09.10  we          Improved arg checking and NAN/INF handling
 1.00.13  09.09.10  we          sfc_fresnel handles NANs
 1.00.14  10.09.10  we          sfc_trigamma
 1.00.15  11.09.10  we          ibeta_series with iteration limit and convergence check
 1.00.16  11.09.10  we          sfBessel: Improved arg checking and NAN/INF handling
 1.00.17  12.09.10  we          Extended over/underflow check in sfc_jn
 1.00.18  14.09.10  we          sfc_eta, sfc_etam1
 1.00.19  15.09.10  we          sfc_zetam1

 1.01.00  06.10.10  we          sfGamma: sfc_dfac
 1.01.01  06.10.10  we          sfPoly:  legendre_pq
 1.01.02  06.10.10  we          sfPoly:  sfc_legendre_p, sfc_legendre_q
 1.01.03  07.10.10  we          sfPoly:  legendre_plmf, sfc_legendre_plm
 1.01.04  10.10.10  we          sfPoly:  sfc_chebyshev_t, sfc_chebyshev_u
 1.01.05  10.10.10  we          fix sign in sfc_chebyshev_t if x<0, n>=NTR
 1.01.06  11.10.10  we          sfPoly:  sfc_gegenbauer_c
 1.01.07  15.10.10  we          sfPoly:  sfc_jacobi_p
 1.01.08  16.10.10  we          sfMisc:  Fix sfc_LambertW1 for WIN16
 1.01.09  16.10.10  we          sfPoly:  sfc_hermite_h
 1.01.10  17.10.10  we          sfPoly:  sfc_laguerre
 1.01.11  18.10.10  we          sfGamma: Lanczos sum with and without expg scale
 1.01.12  19.10.10  we          gamma_delta_ratio, gamma_ratio, pochhammer
 1.01.13  20.10.10  we          sfPoly:  sfc_laguerre_ass, sfc_laguerre_l
 1.01.14  20.10.10  we          sfPoly:  sfc_spherical_harmonic
 1.01.15  22.10.10  we          sfPoly:  sfc_zernike_r
 1.01.16  22.10.10  we          sfPoly:  legendre_plmf with sfc_gamma_ratio
 1.01.17  23.10.10  we          sfGamma: sfc_zetah

 1.02.00  31.10.10  we          sfBessel: sfc_yn: improved pre-checks
 1.02.01  01.11.10  we          sfBessel: sfc_jv, sfc_yv, sfc_bess_jyv
 1.02.02  03.11.10  we          sfBessel: sfc_iv, sfc_kv, sfc_bess_ikv
 1.02.03  04.11.10  we          sfBessel: NAN/INF handling for real order functions
 1.02.04  05.11.10  we          sfMisc  : Tolerance 1.0078125*eps_x in LambertW/1
 1.02.05  05.11.10  we          sfBessel: Airy functions Ai, Ai', Bi, Bi'
 1.02.06  06.11.10  we          sfBessel: sfc_sph_jn, sfc_sph_yn
 1.02.07  06.11.10  we          sfBessel: increased MAXIT in CF1_j
 1.02.08  11.11.10  we          sfGamma : sfc_beta uses sfc_pochhammer for integer x+y=0
 1.02.09  11.11.10  we          sfGamma : Fixed missing check in sfc_lnbeta
 1.02.10  14.11.10  we          sfBessel: Fixed sfc_i0e for small x
 1.02.11  15.11.10  we          sfBessel: sfc_ive, sfc_kve

 1.03.00  01.12.10  we          sfBasic : sfc_bernoulli and small B(2n) as interfaced typed const
 1.03.01  03.12.10  we          sfMisc  : sfc_polylog
 1.03.02  05.12.10  we          sfBasic : sfc_zetaint
 1.03.03  07.12.10  we          sfMisc  : sfc_debye
 1.03.04  07.12.10  we          sfGamma : Fixed sfc_psi(+INF) = +INF
 1.03.05  13.12.10  we          sfGamma : sfc_trigamma with Hurwitz zeta
 1.03.06  14.12.10  we          sfGamma : sfc_tetragamma, sfc_pentagamma
 1.03.07  15.12.10  we          sfGamma : sfc_polygamma
 1.03.08  19.12.10  we          sfGamma : sfc_(ln)gamma return NAN/RTE if x is a non-positive integer
 1.03.09  29.12.10  we          sfGamma : improved sfc_dfac for a few small values to return integer
 1.03.10  30.12.10  we          sfGamma : sfc_binomial

 1.04.00  10.02.11  we          sfBessel: Improved J0,J1 for |x| >= 500; Y0,Y1 for x >= 1600
 1.04.01  12.02.11  we          sfErf   : Goodwin-Staton integral sfc_gsi
 1.04.02  13.02.11  we          sfErf   : sfc_expint3
 1.04.03  14.02.11  we          sfErf   : imaginary error function erfi
 1.04.04  17.02.11  we          sfMisc  : sfc_ri

 1.05.00  31.03.11  we          Fix for some units: uses sfBasic for RTE_ArgumentRange ifopt R+
 1.05.01  10.04.11  we          sfBessel: Zero order Kelvin functions ber,bei,ker,kei
 1.05.02  10.04.11  we          sfc_kei(0) = -Pi/4
 1.05.03  13.04.11  we          sfSDist : sfc_cauchy_pdf/cdf/inv
 1.05.04  14.04.11  we          sfSDist : sfc_normal_pdf/cdf/inv
 1.05.05  14.04.11  we          sfSDist : sfc_exp_pdf/cdf/inv
 1.05.06  15.04.11  we          sfSDist : sfc_lognormal_pdf/cdf/inv
 1.05.07  16.04.11  we          sfSDist : sfc_logistic_pdf/cdf/inv
 1.05.08  16.04.11  we          sfSDist : sfc_weibull_pdf/cdf/inv
 1.05.09  17.04.11  we          sfSDist : sfc_laplace_pdf/cdf/inv
 1.05.10  17.04.11  we          sfSDist : changed sh,sc to a,b in sfc_gamma_pdf/cdf/inv
 1.05.11  19.04.11  we          sfSDist : sfc_pareto_pdf/cdf/inv
 1.05.12  20.04.11  we          sfSDist : sfc_uniform_pdf/cdf/inv

 1.06.00  05.05.11  we          sfBessel: sfc_struve_h0
 1.06.01  05.05.11  we          sfBessel: sfc_struve_h1
 1.06.02  06.05.11  we          sfBessel: sfc_struve_l1
 1.06.03  07.05.11  we          sfBessel: sfc_struve_l0
 1.06.04  12.05.11  we          sfc_sncndn moved from sfMisc to sfEllInt
 1.06.05  13.05.11  we          sfEllInt: sfc_theta2/3/4
 1.06.06  14.05.11  we          sfEllInt: sfc_theta1p
 1.06.07  16.05.11  we          sfEllInt: sfc_jtheta
 1.06.08  17.05.11  we          sfEllInt: hack: fix overflow for theta_ts for q near 1
 1.06.09  19.05.11  we          sfEllInt: sfc_ellmod
 1.06.10  19.05.11  we          replaced misleading kc by k in EllipticCK and EllipticCE
 1.06.11  20.05.11  we          sfEllInt: sfc_nome
 1.06.12  22.05.11  we          sfEllInt: improved sfc_nome with Chebyshev for k <= 1/8
 1.06.13  26.05.11  we          sfBessel: fix for very large x in kelvin_large
 1.06.14  26.05.11  we          sfEllInt: improved NAN handling for theta functions

 1.07.00  01.06.11  we          sfExpInt: sfc_ein
 1.07.01  02.06.11  we          sfEllInt: reordered argument check in sfc_ell_rc
 1.07.02  02.06.11  we          sfEllInt: sfc_ellint_1/2
 1.07.03  03.06.11  we          sfEllInt: sfc_ellint_2: case k=1
 1.07.04  03.06.11  we          sfEllInt: sfc_EllipticPiC/CPi with sfc_cel
 1.07.05  04.06.11  we          sfEllInt: inline agm in sfc_EllipticK
 1.07.06  04.06.11  we          sfEllInt: sfc_ellint_2: fix phi close to n*Pi_2
 1.07.07  05.06.11  we          comp_ellint_1/2/3(x)
 1.07.08  06.06.11  we          sfEllInt: sfc_jam
 1.07.09  07.06.11  we          sfEllInt: sfc_jzeta
 1.07.10  08.06.11  we          sfEllInt: sfc_ellint_3
 1.07.11  10.06.11  we          sfEllInt: fix D2/D3 issue in sfc_ellint_3
 1.07.12  10.06.11  we          sfErf:    sfc_erfce
 1.07.13  15.06.11  we          sfEllInt: parameter check in sfc_EllipticK
 1.07.14  16.06.11  we          sfEllInt: sfc_jam for |k| > 1
 1.07.15  17.06.11  we          sfEllInt: sfc_ellint_1 for |k| > 1
 1.07.16  17.06.11  we          sfEllInt: sfc_jacobi_arcsn/cn/dn
 1.07.17  18.06.11  we          sfEllInt: sfc_jacobi_sn/cn/dn
 1.07.18  19.06.11  we          sfEllInt: sfc_ellint_2/3 for |k| > 1
 1.07.19  19.06.11  we          SpecFun:  Fix result for Nan argument in LambertW/1
 1.07.20  19.06.11  we          sfMisc:   Move IsNanOrInf test to top in LambertW1
 1.07.21  23.06.11  we          sfEllInt: sfc_hlambda
 1.07.22  26.06.11  we          sfMisc:   simplified range check in sfc_LambertW/1

 1.08.00  03.07.11  we          sfGamma:  Asymptotic expansion for x>1024 in stirf
 1.08.01  12.07.11  we          sfGamma:  polygam_negx: polygamma for x<0 and 1<=n<=10
 1.08.02  17.07.11  we          sfExpInt: range sfc_li(x) now x>=0, x<>1
 1.08.03  26.07.11  we          sfEllInt: improved sfc_EllipticEC
 1.08.04  26.07.11  we          sfEllInt: nu<>1 in sfc_EllipticCPi,sfc_EllipticPiC
 1.08.05  27.07.11  we          sfEllInt: Reciprocal-Modulus Transformation in sfc_EllipticE
 1.08.06  02.08.11  we          sfEllInt: sfc_el1 uses arcsinh for kc=0
 1.08.07  03.08.11  we          sfGamma:  invgamma renamed to rgamma
 1.08.08  04.08.11  we          sfGamma:  small simplification in sfc_gstar
 1.08.09  06.08.11  we          sfGamma:  sfc_dfac for odd negative n, fac/dfac return INF if undefined
 1.08.10  08.08.11  we          sfGamma:  const ETAEPS, improved etam1pos
 1.08.11  09.08.11  we          sfGamma:  improved sfc_zetah
 1.08.12  13.08.11  we          sfBessel: cases |v|=0,1 in sfc_kve
 1.08.13  13.08.11  we          sfBessel: v=-1 in sfc_iv/e
 1.08.14  15.08.11  we          sfSDist:  Fix sfc_t_cdf for t>0 and large nu
 1.08.15  15.08.11  we          sfErf:    Range check sfc_erf(c)_inv
 1.08.16  15.08.11  we          sfSDist:  sfc_normstd_inv uses sfc_erfc_inv

 1.09.00  14.09.11  we          sfGamma:  improved lngamma_small
 1.09.01  14.09.11  we          sfGamma:  sfc_lngamma1p
 1.09.02  17.09.11  we          sfGamma:  lanczos interfaced
 1.09.03  17.09.11  we          sfGamma:  sfc_ibetaprefix
 1.09.04  18.09.11  we          sfGamma:  improved sfc_ibeta
 1.09.05  06.10.11  we          sfGamma:  sfc_zetah for 0 <= s < 1

 1.10.00  24.11.11  we          sfGamma:  sfc_lnfac
 1.10.01  10.12.11  we          sfSDist:  sfc_triangular_pdf/cdf/inv
 1.10.02  13.12.11  we          sfSDist:  sfc_beta_pdf uses sfc_ibetaprefix
 1.10.03  14.12.11  we          sfSDist:  sfc_binomial_cdf/pmf
 1.10.04  15.12.11  we          sfGamma:  sfc_igprefix for TP5 without Lanczos
 1.10.05  15.12.11  we          sfSDist:  sfc_poisson_cdf/pmf
 1.10.06  16.12.11  we          sfMisc:   sfc_debye with DE integration for x>4 and n>6
 1.10.07  18.12.11  we          sfSDist:  sfc_negbinom_cdf/pmf
 1.10.08  21.12.11  we          sfGamma:  TP5 for sfc_zetah/sfc_zetam1
 1.10.09  21.12.11  we          sfMisc:   sfc_pz
 1.10.10  22.12.11  we          sfGamma:  etam1pos with exp7
 1.10.11  22.12.11  we          sfMisc:   sfc_pz with Cohen's speed-up $ifdef BASM
 1.10.12  26.12.11  we          sfSDist:  sfc_hypergeo_cdf/pmf
 1.10.13  27.12.11  we          sfSDist:  Fix sfc_hypergeo_pmf(0,0,0,k)
 1.10.14  28.12.11  we          sfBessel: sfc_sph_in, sfc_sph_ine
 1.10.15  29.12.11  we          sfBessel: sfc_sph_kn, sfc_sph_kne
 1.10.16  29.12.11  we          sfBessel: Fix bug for small negative x in sfc_sph_ine
 1.10.17  02.01.12  we          sfMisc:   improved sfc_polylog for n<-9 and -4 < x < -11.5

 1.11.00  03.02.12  we          sfZeta:   New unit with zeta functions and polylogarithms
 1.11.01  04.02.12  we          sfZeta:   sfc_lerch
 1.11.02  05.02.12  we          sfZeta:   sfc_lerch for z=-1: use lphi_aj instead of zetah
 1.11.03  06.02.12  we          sfZeta:   sfc_dbeta
 1.11.04  06.02.12  we          sfZeta:   sfc_dlambda
 1.11.05  07.02.12  we          sfZeta:   sfc_lchi
 1.11.06  08.02.12  we          sfZeta:   sfc_polylogr
 1.11.07  15.02.12  we          SpecFun:  use Ext2Dbl with el3, gamma, beta, e1

 1.12.00  25.02.12  we          sfMisc:   Avoid over/underflow for large x in debye_dei
 1.12.01  13.03.12  we          sfZeta:   sfc_lchi(1,x) = arctanh(x)
 1.12.02  20.03.12  we          sfBasic:  BoFHex array: Bernoulli(2k+2)/(2k+2)!
 1.12.03  20.03.12  we          sfZeta:   sfc_zetah with BoFHex
 1.12.04  21.03.12  we          sfGamma:  sfc_gamma_delta_ratio for negative arguments
 1.12.05  21.03.12  we          sfGamma:  improved sfc_pochhammer
 1.12.06  22.03.12  we          sfGamma:  sfc_poch1
 1.12.07  04.04.12  we          sfPoly:   sfc_legendre_p, sfc_legendre_plm for |x| > 1
 1.12.08  09.04.12  we          sfGamma:  fix special case if prefix is zero in sfc_incgamma_ex
 1.12.09  11.04.12  we          sfGamma:  sfc_igamma (non-normalised upper incomplete gamma)
 1.12.10  17.04.12  we          sfPoly:   sfc_legendre_q for |x| > 1
 1.12.11  18.04.12  we          sfGamma:  sfc_lngammas
 1.12.12  22.04.12  we          sfGamma:  sfc_gamma for x < -MAXGAMX
 1.12.13  22.04.12  we          sfGamma:  sfc_igammal (non-normalised lower incomplete gamma)
 1.12.14  27.04.12  we          sfPoly:   sfc_legendre_qlm
 1.12.15  16.05.12  we          sfPoly:   sfc_thq (toroidal harmonics)
 1.12.16  20.05.12  we          sfPoly:   sfc_thp (toroidal harmonics)
 1.12.17  24.05.12  we          sfPoly:   sfc_thq for m < 0
 1.12.18  29.05.12  we          sfBasic:  fix missing bit in B2nHex[0]

 1.13.00  07.06.12  we          sfGamma:  fix sfc_igprefix for a<1 and underflow of exp(-x)
 1.13.01  08.06.12  we          sfGamma:  improved sfc_rgamma
 1.13.02  09.06.12  we          sfGamma:  sfc_taylor
 1.13.03  09.06.12  we          sfZeta:   sfc_polylog uses sfc_taylor
 1.13.04  09.06.12  we          sfZeta:   sfc_fermi_dirac
 1.13.05  10.06.12  we          sfZeta:   fix n>x in fd_asymp_exp_int
 1.13.06  11.06.12  we          sfZeta:   sfc_zeta uses sfc_zetaint
 1.13.07  12.06.12  we          sfZeta:   sfc_etaint
 1.13.08  12.06.12  we          sfBasic:  sfc_diagctr
 1.13.09  13.06.12  we          sfZeta:   sfc_fdm05, sfc_fdp05
 1.13.10  15.06.12  we          sfErf:    sfc_gendaw (generalized Dawson integral)
 1.13.11  16.06.12  we          sfErf:    sfc_inerfc (repeated integrals of erfc)
 1.13.12  20.06.12  we          sfExpInt: corrected bound in sfc_en (11140.0 -> 11390.2)
 1.13.13  22.06.12  we          sfExpInt: uniform asymptotic expansion ep_largep for n>=10000
 1.13.14  23.06.12  we          sfExpInt: generalized exponential integral sfc_gei
 1.13.15  26.06.12  we          sfErf:    sfc_erfg (generalized error function)
 1.13.16  29.06.12  we          sfBasic:  sfc_write_debug_str
 1.13.17  30.06.12  we          sfGamma:  sfc_git (Tricomi's incomplete gamma)
 1.13.18  01.07.12  we          sfGamma:  check loss of accuracy in sfc_igammal

 1.14.00  20.07.12  we          sfGamma:  removed special treatment for a=0 in sfc_poch
 1.14.01  23.07.12  we          sfExpInt: more compact sfc_eiex
 1.14.02  31.01.13  we          sfBasic:  sfc_write_debug_str: fix and D17 adjustment

 1.15.00  15.02.13  we          sfGamma:  improved sfc_git
 1.15.01  15.02.13  we          sfErf:    fix sfc_expint3 for 0.4e-6 .. 6e-6
 1.15.02  16.02.13  we          sfGamma:  fix sfc_trigamma/tetragamma/pentagamma for very large arguments
 1.15.03  16.02.13  we          sfZeta:   improved sfc_zetah
 1.15.04  19.02.13  we          sfBessel: improved quick check in sfc_jn
 1.15.05  20.02.13  we          sfBessel: Yv_series
 1.15.06  20.02.13  we          sfBessel: handle some near overflows in bessel_jy
 1.15.07  20.02.13  we          sfBessel: improved check in sfc_yn

 1.16.00  14.03.13  we          sfEllInt: Remaining Jacobi elliptic functions
 1.16.01  15.03.13  we          sfEllInt: Improved sfc_jacobi_arcsn/cn/dn
 1.16.02  20.03.13  we          sfZeta:   fix sfc_zetah for j>NBoF
 1.16.03  22.03.13  we          sfEllInt: sfc_jacobi_arcsc, sfc_jacobi_arccs
 1.16.04  23.03.13  we          sfEllInt: sfc_jacobi_arcnc, sfc_jacobi_arcns
 1.16.05  23.03.13  we          sfEllInt: sfc_jacobi_arcnd
 1.16.06  24.03.13  we          sfEllInt: sfc_jacobi_arccd, sfc_jacobi_arcdc
 1.16.07  24.03.13  we          sfEllInt: sfc_jacobi_arcsd
 1.16.08  25.03.13  we          sfEllInt: sfc_jacobi_arcds
 1.16.09  27.03.13  we          sfZeta:   sfc_zetah for s<0
 1.16.10  28.03.13  we          sfZeta:   sfc_trilog
 1.16.11  29.03.13  we          sfZeta:   TP5 fix for sfc_zeta
 1.16.12  29.03.13  we          sfZeta:   sfc_zetah for a near 1 and s<0
 1.16.13  30.03.13  we          sfGamma:  improved sfc_binomial

 1.17.00  07.04.13  we          sfSDist:  Improved sfc_chi2_pdf for very small x
 1.17.01  10.04.13  we          sfZeta:   sfc_fdp15
 1.17.02  13.04.13  we          sfSDist:  k changed to longint in sfc_poisson_pmf
 1.17.03  14.04.13  we          sfSDist:  sfc_rayleigh_pdf/cdf/inv
 1.17.04  15.04.13  we          sfGamma:  improved sfc_beta
 1.17.05  16.04.13  we          sfGamma:  ibeta_series with parameter normalised
 1.17.06  17.04.13  we          sfGamma:  sfc_nnbeta
 1.17.07  19.04.13  we          sfSDist:  sfc_maxwell_pdf/cdf/inv
 1.17.08  20.04.13  we          sfSDist:  sfc_evt1_pdf/cdf/inv
 1.17.09  20.04.13  we          sfSDist:  sfc_maxwell_cdf/inv with lower igamma functions
 1.17.10  21.04.13  we          sfGamma:  lanczos_gm05 = lanczos_g - 0.5 in implementation
 1.17.11  27.04.13  we          sfGamma:  polygam_negx for n=11,12
 1.17.12  01.05.13  we          FresnelC/FresnelS

 1.18.00  09.05.13  we          sfBessel: Airy/Scorer functions Gi/Hi
 1.18.01  09.05.13  we          sfBessel: Prevent some wrong compiler optimizations for div by 3
 1.18.02  11.05.13  we          sfGamma:  improved sfc_polygamma for large order: e.g. (3000,200)
 1.18.03  14.05.13  we          sfHyperG: 2F1: major rewrite (additional linear transformations)
 1.18.04  14.05.13  we          sfHyperG: 2F1: fix y=0 in h2f1_sum
 1.18.05  16.05.13  we          sfHyperG: 2F1: removed recurrence on c (to make c-a-b > 0)
 1.18.06  17.05.13  we          sfHyperG: 2F1: h2f1poly with parameter err
 1.18.07  18.05.13  we          sfHyperG: 2F1: no Luke for TP5
 1.18.08  19.05.13  we          (some)    Prevent some wrong compiler optimizations for div by 3
 1.18.09  21.05.13  we          sfHyperG: 2F1: regularized function sfc_2f1r
 1.18.10  21.05.13  we          sfHyperG: 1F1: h1f1_sum_laguerre
 1.18.11  23.05.13  we          sfHyperG: 1F1: h1f1_tricomi only if 2bx > 4ax
 1.18.12  23.05.13  we          sfHyperG: 1F1: fix h1f1_acf if x >= Ln_MaxExt
 1.18.13  24.05.13  we          sfHyperG: 1F1: regularized function sfc_1f1r
 1.18.14  26.05.13  we          sfGamma:  sfc_nnbeta for a<=0 or b<=0
 1.18.15  27.05.13  we          sfHyperG: CHU: handle a or b negative integer
 1.18.16  28.05.13  we          sfHyperG: CHU: special cases a=-1, x=0, a=b
 1.18.17  29.05.13  we          sfHyperG: CHU: case b=2a, allow a or a+1-b negative integer in luke
 1.18.18  30.05.13  we          sfHyperG: CHU: handle b=a+n+1
 1.18.19  31.05.13  we          sfHyperG: CHU: fix some double/overflow cases in dchu

 1.19.00  06.06.13  we          sfBessel: Fix typo in h1v_large
 1.19.01  08.06.13  we          sfBessel: sfc_bess_kv2
 1.19.02  09.06.13  we          sfHyperG: CHU: fix sign in special case a=-1
 1.19.03  12.06.13  we          sfHyperG: CHU: implementation of Temme's chu and uabx
 1.19.04  14.06.13  we          sfHyperG: CHU: use dchu with Kummer if b<0.0 and t>T_MAXA
 1.19.05  16.06.13  we          sfHyperG: Whittaker functions sfc_whitm/sfc_whitw
 1.19.06  27.06.13  we          sfBessel: Check overflow / return INF in bessel_jy
 1.19.07  27.06.13  we          sfHyperG: sfc_0f1
 1.19.08  28.06.13  we          sfHyperG: gen_0f1, sfc_0f1r

 1.20.00  27.07.13  we          sfHyperG: MAXIT=64 in h1f1_tricomi
 1.20.01  27.07.13  we          sfGamma:  sfc_git for x<0
 1.20.02  14.08.13  we          sfSDist:  Removed code fragment for y>1 in sfc_beta_inv
 1.20.03  14.08.13  we          sfSDist:  Check a,b > 0 in sfc_beta_cdf
 1.20.04  15.08.13  we          sfSDist:  sfc_kumaraswamy_pdf/cdf/inv
 1.20.05  16.08.13  we          sfZeta:   sfc_zeta for small s
 1.20.06  16.08.13  we          sfZeta:   sfc_polylog with sfc_zetaint
 1.20.07  17.08.13  we          sfSDist:  Improved sfc_lognormal_pdf with expmx2h
 1.20.08  18.08.13  we          sfErf:    sfc_erf_p, sfc_erf_q, sfc_erf_z
 1.20.09  18.08.13  we          sfSDist:  normal/normstd_pdf/cdf with erf_z and erf_p

 1.21.00  07.09.13  we          sfZeta:   Improved sfc_zetah/hurwitz_formula
 1.21.01  11.09.13  we          sfZeta:   Bernoulli polynomials sfc_bernpoly
 1.21.02  14.09.13  we          sfPoly:   sfc_chebyshev_v
 1.21.03  15.09.13  we          sfPoly:   sfc_chebyshev_w
 1.21.04  24.09.13  we          sfMisc:   sfc_cosint, sfc_sinint
 1.21.05  25.09.13  we          sfGamma:  sfc_batemang
 1.21.06  26.09.13  we          sfMisc:   Fixed quasi-periodic code for sfc_cos/sinint
 1.21.07  27.09.13  we          sfEllInt: Fixed quasi-periodic code for sfc_ellint_1/2/3
 1.21.08  28.09.13  we          sfEllInt: special_reduce_modpi, used in sfc_ellint_1/2/3 and sfc_hlambda

 1.22.00  19.10.13  we          sfSDist:  sfc_moyal_pdf/cdf/inv
 1.22.01  22.10.13  we          sfEllInt: sfc_EllipticKim
 1.22.02  23.10.13  we          sfEllInt: sfc_EllipticECim
 1.22.03  04.11.13  we          sfSDist:  improved sfc_t_pdf for small x^2/nu

 1.23.00  26.12.13  we          sfZeta:   sfc_fdp25
 1.23.01  26.12.13  we          sfBasic:  complete references in amath_info.txt

 1.24.00  11.03.14  we          sfExpInt: Adjust arrays sizes in sfc_ci, si_medium
 1.24.01  24.03.14  we          sfHyperG: chu_luke only with $ifdef in dchu
 1.24.02  26.03.14  we          sfBessel: sfc_yn with LnPi from AMath

 1.25.00  16.04.14  we          sfGamma:  Improved polygam_negx (larger n for x<0}
 1.25.01  01.05.14  we          sfZeta:   sfc_harmonic
 1.25.02  03.05.14  we          sfZeta:   sfc_harmonic2
 1.25.03  04.05.14  we          sfSDist:  sfc_zipf_pmf/cdf

 1.26.00  24.05.14  we          sfMisc:   sfc_ali (functional inverse of sfc_li)
 1.26.01  27.05.14  we          sfZeta:   Removed redundancy in sfc_harmonic2
 1.26.02  28.05.14  we          sfZeta:   sfc_harmonic with nch = 22
 1.26.03  29.05.14  we          sfZeta:   polylogneg with array of single
 1.26.04  30.05.14  we          SpecFun:  Fix zipf_pmf/cdf function type
 1.26.05  30.05.14  we          sfZeta:   Lobachevski sfc_llci, sfc_llsi
 1.26.06  31.05.14  we          sfZeta:   improved harm2core
 1.26.07  01.06.14  we          sfMisc:   Fibonacci polynomials sfc_fpoly
 1.26.08  02.06.14  we          sfMisc:   Lucas polynomials sfc_lpoly
 1.26.09  05.06.14  we          sfSDist:  sfc_levy_pdf/cdf/inv

 1.27.00  15.06.14  we          Inverse/incomplete Gamma/Beta moved to unit sfGamma2
 1.27.01  15.06.14  we          sfGamma:  const lanczos_gm05: extended
 1.27.02  15.06.14  we          sfGamma2: sfc_ibeta_inv (code from sfSDist)
 1.27.03  15.06.14  we          sfZeta:   Fermi/Dirac, Lobachewsky, harmonic functions moved to sfZeta2
 1.27.04  22.06.14  we          sfGamma2: sfc_ilng (inverse of lngamma)
 1.27.05  25.06.14  we          sfGamma2: sfc_ipsi (inverse of psi/digamma)
 1.27.06  02.07.14  we          SpecFun/X ibeta_inv/x
 1.27.07  11.07.14  we          sfHyperG: sfc_pcfd
 1.27.08  11.07.14  we          sfHyperG: sfc_pcfu
 1.27.09  12.07.14  we          sfHyperG: sfc_pcfhh

 1.28.00  29.07.14  we          sfBasic:  moebius[1..95]
 1.28.01  29.07.14  we          sfZeta:   Improved and expanded primezeta sfc_pz
 1.28.02  05.08.14  we          sfZeta:   sfc_zeta for s < -MaxGAMX
 1.28.03  05.08.14  we          sfGamma:  First check IsNaN(x) in sfc_lngamma
 1.28.04  05.08.14  we          sfGamma:  improved sfc_gdr_pos for small x
 1.28.05  08.08.14  we          sfBessel: Iv := PosInf_x if Kv,Kv1=0 in bessel_ik or if x >= IvMaxXH in sfc_iv
 1.28.06  10.08.14  we          sfGamma:  fix sfc_lnbeta for VER50
 1.28.07  11.08.14  we          sfGamma:  special case x=x+y in sfc_beta
 1.28.08  17.08.14  we          sfMisc:   Catalan function sfc_catf
 1.28.09  19.08.14  we          sfHyperG: sfc_pcfv
 1.28.10  20.08.14  we          sfGamma:  sfc_gamma_delta_ratio returns 1 if x=x+d

 1.29.00  01.09.14  we          sfHyperG: sfc_chu for x<0
 1.29.01  04.09.14  we          sfHyperG: fix h1f1_tricomi (division by t=0)
 1.29.02  07.09.14  we          sfSDist:  sfc_logistic_inv uses logit function
 1.29.03  06.10.14  we          sfSDist:  sfc_logistic_cdf uses logistic function

 1.30.00  26.10.14  we          sfHyperG: No Kummer transformation in 1F1 if b=0,-1,-2,...
 1.30.01  28.10.14  we          sfHyperG: 1F1 = 1 + a/b*(exp(x)-1) for very small a,b
 1.30.02  06.11.14  we          sfBasic:  Small Bernoulli numbers B2nHex moved to AMath

 1.31.00  14.12.14  we          sfGamma:  small changes in sfc_batemang
 1.31.01  16.12.14  we          sfGamma2: Check IsNanOrInf in sfc_ilng
 1.31.02  23.12.14  we          sfSDist:  sfc_invgamma_pdf/cdf/inv
 1.31.03  26.12.14  we          sfSDist:  logarithmic (series) distribution
 1.31.04  31.12.14  we          sfSDist:  sfc_wald_cdf/pdf
 1.31.05  01.01.15  we          sfSDist:  sfc_wald_inv
 1.31.06  02.01.15  we          sfSDist:  Improved sfc_wald_inv (restart with mode)
 1.31.07  05.01.15  we          sfHyperG: Special case b=1 in gen_0f1
 1.31.08  05.01.15  we          sfSDist:  Error if k<1 in sfc_ls_cdf/pmf

 1.32.00  25.01.15  we          sfExpint: sfc_ali from sfMisc
 1.32.01  25.01.15  we          sfExpint: sfc_ei_inv
 1.32.02  26.03.15  we          sfEllInt: sfc_cel_d
 1.32.03  27.03.15  we          sfEllInt: sfc_ellint_d
 1.32.04  01.05.15  we          sfExpint: fix sfc_ei_inv for x Nan

 1.33.00  17.05.15  we          sfZeta:   sfc_lerch for s >= -1
 1.33.01  18.05.15  we          sfZeta:   sfc_polylogr for s >= -1
 1.33.02  18.05.15  we          sfZeta:   Improved sfc_etam1 (adjusted Borwein constants)
 1.33.03  24.05.15  we          sfZeta:   sfc_lerch(1,s,a) for s<>1 and special case z=0
 1.33.04  04.06.15  we          sfMisc:   sfc_kepler
 1.33.05  04.06.15  we          sfBasic:  remove CSInitX
 1.33.06  07.06.15  we          sfBessel: new IJ_series replaces Jv_series
 1.33.07  07.06.15  we          sfBessel: rewrite of sfc_in: use IJ_series and CF1_I
 1.33.08  08.06.15  we          sfBessel: improved bess_m0p0/bess_m1p1 with rem_2pi_sym and Kahan summation
 1.33.09  09.06.15  we          sfBessel: sfc_struve_l
 1.33.10  10.06.15  we          sfBessel: sfc_struve_h
 1.33.11  10.06.15  we          sfBessel: sfc_struve_l0/1 use sfc_struve_l
 1.33.12  14.06.15  we          sfBessel: sfc_struve_h/l with real nu >= 0
 1.33.13  19.06.15  we          sfBessel: scaled Airy functions sfc_airy_ais, sfc_airy_bis
 1.33.14  21.06.15  we          sfEllInt: avoid overflow in sfc_ell_rf, better arg checks in Carlson functions
 1.33.15  22.06.15  we          sfEllInt: Carlson sfc_ell_rg

 1.34.00  05.07.15  we          sfBessel: special cases Yv(+- 1/2, x)
 1.34.01  06.07.15  we          sfBessel: sfc_struve_h for v < 0
 1.34.02  07.07.15  we          sfBessel: sfc_struve_l for v < 0
 1.34.03  09.07.15  we          sfBessel: fix sfc_struve_h/l(v,0) for v<=-1
 1.34.04  18.07.15  we          sfZeta:   eta near s=1
 1.34.04  18.07.15  we          sfZeta:   eta near s=1
 1.34.05  19.07.15  we          sfEllint: sfc_el1(x, kc) = arctan(x) for kc^2=1
 1.34.06  01.08.15  we          sfHyperG: Avoid some overflows in h2f1_sum, improved h2f1_trans

 1.35.00  25.08.15  we          (some):   Use THREE & SIXX to avoid 'optimization'
 1.35.01  25.08.15  we          sfZeta:   sfc_bernpoly: avoid integer overflow and FPC optimization error
 1.35.02  26.08.15  we          sfBessel: Sqrt2 in kelvin functions (anti-'optimizing')

 1.36.00  26.09.15  we          sfEllInt: lemniscate functions sin_lemnx, cos_lemnx
 1.36.01  29.09.15  we          sfHyperG: sfc_pcfv(a,x) for a<0
 1.36.02  11.10.15  we          sfGamma:  polygamma for x<0 uses Euler series for Pi*cot(Pi*x)

 1.37.00  14.02.16  we          sfErf:    Constants for trivial ranges of erf/erfc/erfce

 1.38.00  11.05.16  we          sfEllInt: sfc_cel_b
 1.38.01  12.05.16  we          sfEllInt: sfc_ellint_b
 1.38.02  17.06.16  we          sfMisc:   Bring radical sfc_br
 1.38.03  20.06.16  we          sfMisc:   Improved sfc_br
 1.38.04  25.06.16  we          sfMisc:   Wright omega sfc_wo
 1.38.05  26.06.16  we          sfMisc:   Nan/Inf handling for sfc_wo
 1.38.06  27.06.16  we          sfMisc:   Nan/Inf handling for sfc_br

 1.39.00  23.08.16  we          sfGamma2: sfc_invgam (inverse of gamma)
 1.39.01  28.09.16  we          sfMisc:   General Fibonacci function
 1.39.02  29.09.16  we          sfPoly:   sfc_chebf1
 1.39.03  06.10.16  we          sfPoly:   sfc_hermite_he
 1.39.04  08.10.16  we          sfSDist:  sfc_ks_cdf/inv

 1.40.00  19.05.17  we          sfEllInt: sfc_EllipticK  for k>1
 1.40.01  19.05.17  we          sfEllInt: sfc_EllipticEC for k>1
 1.40.02  25.05.17  we          sfMisc:   sfc_eulerq
 1.40.03  26.05.17  we          sfMisc:   Improved AE for sfc_eulerq if q ~ -1
 1.40.04  08.06.17  we          sfSDist:  sfc_cauchy_inv uses tanpi
 1.40.05  20.06.17  we          sfMisc:   sfc_detai
 1.40.06  27.06.17  we          sfMisc:   Nan/Inf handling for sfc_detai
 1.40.07  29.06.17  we          (some):   Code removed for TP5-TP6, TPW1-D1

 1.41.00  07.07.17  we          sfErf:    new sfc_fresnel_fg, improved sfc_fresnel
 1.41.01  10.07.17  we          sfErf:    sfc_fresnel_fg for x<0
 1.41.02  16.07.17  we          sfBessel: Bessel lambda sfc_blam
 1.41.03  18.07.17  we          sfErf:    New sfc_erfh
 1.41.04  19.07.17  we          sfErf:    New sfc_erf2
 1.41.05  22.07.17  we          sfSDist:  Fix sfc_chi2_pdf for x=0
 1.41.06  22.07.17  we          sfSDist:  sfc_chi_pdf/cdf/inv

 1.42.00  16.08.17  we          sfExpInt: sfc_eisx2
 1.42.01  17.08.17  we          sfErf:    sfc_gsi for x < 0

 1.43.00  16.08.17  we          sfErf:    sfc_erfi_inv
 1.43.01  02.12.17  we          sfErf:    sfc_expint3 with sfc_igammal
 1.43.02  02.12.17  we          (some):   Suppress warnings: Local variable does not seem to be initialized
 1.43.03  05.12.17  we          sfSDist:  range checks in sfc_gamma_inv
 1.43.04  06.12.17  we          sfSDist:  sfc_nakagami_pdf/cdf/inv
 1.43.05  06.12.17  we          specfun:  Ext2Dbl for some pdf/inv
 1.43.06  08.12.17  we          sfMisc:   generalized Marcum Q function sfc_mq
 1.43.07  12.12.17  we          sfExpInt: separate functions e1_small, e1ex_large
 1.43.08  13.12.17  we          sfExpInt: sfc_e1s(x) = exp(x)*E1(x)
 1.43.09  14.12.17  we          sfExpInt: function ei_asymp, sfc_eisx2 with ei_asymp
 1.43.10  14.12.17  we          sfExpInt: sfc_eis(x) = exp(-x)*Ei(x)
 1.43.11  18.12.17  we          sfSDist:  sfc_ks_inv with modified regula falsi

 1.44.00  04.01.18  we          sfMisc:   Improved sfc_ri
 1.44.01  10.01.18  we          sfMisc:   sfc_ari
 1.44.02  17.01.18  we          sfErf:    minor changes in sfc_erf_p/q
 1.44.03  18.01.18  we          sfErf:    Owen's T function sfc_owent
 1.44.04  23.01.18  we          sfErf:    sfc_erfce_inv
 1.44.05  29.01.18  we          sfMisc:   sfc_einstein
 1.44.06  01.02.18  we          sfMisc:   rewrite sfc_debye (always integrate if x>4}
 1.44.07  02.02.18  we          sfMisc:   Transport integrals sfc_trint
 1.44.08  04.02.18  we          sfZeta:   Rewrite polylog, returns real part for x>1,n>1
 1.44.09  04.02.18  we          sfZeta:   sfc_zetaint(1) = NaN
 1.44.10  05.02.18  we          sfZeta:   sfc_polylog(n,x) for x > 1
 1.44.11  06.02.18  we          sfbasic:  sfc_intde/i_p
 1.44.12  06.02.18  we          sfBessel: sfc_struve_h with sfc_intdei_p
 1.44.13  07.02.18  we          sfZeta:   sfc_polylogr for x < -1
 1.44.14  07.02.18  we          sfbasic:  conditional define amtools_intde
 1.44.15  08.02.18  we          sfZeta:   Merged with (old) unit sfZeta2
 1.44.16  09.02.18  we          sfZeta:   sfc_fdr (Fermi-Dirac real order)
 1.44.17  12.02.18  we          sfZeta:   sfc_polylogr(s,x)=x for large s

 1.45.00  20.02.18  we          sfMisc:   sfc_langevin: Langevin function from AMath
 1.45.01  21.02.18  we          sfMisc:   sfc_llinv: Inverse Langevin function
 1.45.02  24.02.18  we          sfZeta:   sfc_lerch for z < -1
 1.45.03  25.02.18  we          sfZeta:   fix NaN check for r in sfc_harmonic2
 1.45.04  26.02.18  we          sfZeta:   fix overflow for LerchPhi integrand
 1.45.05  26.02.18  we          sfZeta:   sfc_ti (inverse tangent integral of real order)
 1.45.06  01.03.18  we          sfZeta:   improved lerch integration
 1.45.07  03.03.18  we          sfZeta:   sfc_polylogr for 1 < x <= 256
 1.45.08  05.03.18  we          sfZeta:   sfc_lchi for |x| > 1
 1.45.09  05.03.18  we          sfZeta:   Bose-Einstein integral sfc_beint
 1.45.10  06.03.18  we          sfZeta:   sfc_polylogr for x > 256
 1.45.11  11.03.18  we          sfMisc:   Rogers-Ramanujan continued fraction sfc_rrcf

 1.46.00  21.03.18  we          sfBessel: sfc_j0int
 1.46.01  22.03.18  we          sfBessel: sfc_y0int
 1.46.02  23.03.18  we          sfBessel: sfc_k0int
 1.46.03  24.03.18  we          sfBessel: sfc_i0int
 1.46.04  26.03.18  we          specfun:  Fix nakagami return types
 1.46.05  27.03.18  we          sfBessel: improved sfc_[i,j,y]0int
 1.46.06  30.03.18  we          sfEllInt: sfc_cel_d for |k| >= 1
 1.46.07  30.03.18  we          sfEllInt: sfc_cel_b for |k| >= 1
 1.46.08  31.03.18  we          sfEllInt: sfc_EllipticPiC for |k| > 1
 1.46.09  02.04.18  we          sfEllInt: sfc_EllipticPiCim
 1.46.10  04.04.18  we          sfEllInt: Fix sfc_EllipticPiC for k^2 > nu > 1
 1.46.11  04.04.18  we          sfEllInt: Fix sfc_EllipticPiC for k^2 = nu > 1
 1.46.12  05.04.18  we          sfEllInt: Mathematica style complete elliptic integrals
 1.46.13  05.04.18  we          sfEllInt: Removed conditional define use_bulirsch (same as DAMath)
 1.46.14  12.04.18  we          sfEllInt: Neville theta functions
 1.46.15  15.04.18  we          sfEllInt: sfc_arccl/sl

 1.47.00  19.04.18  we          (some):   Fixes for FPC311
 1.47.01  22.04.18  we          sfEllInt: sfc_el2 computed with Carlson
 1.47.02  23.04.18  we          sfEllInt: sfc_mlambda(i*y): elliptic modular function
 1.47.03  24.04.18  we          sfEllInt: sfc_detai from sfmisc
 1.47.04  25.04.18  we          sfEllInt: sfc_wpl
 1.47.05  26.04.18  we          sfEllInt: sfc_wpe, sfc_wpe_im
 1.47.06  27.04.18  we          sfEllInt: sfc_wpe_inv
 1.47.07  27.04.18  we          sfEllInt: sfc_wpg, sfc_wpg_im
 1.47.08  28.04.18  we          sfEllInt: sfc_wpg with JacobiCN
 1.47.09  30.04.18  we          sfEllInt: CompAuxWP
 1.47.10  30.04.18  we          sfEllInt: equireduce
 1.47.11  01.05.18  we          sfEllInt: sfc_wpg uses wpdup for small x
 1.47.12  03.05.18  we          sfEllInt: sfc_wpg_inv
 1.47.13  04.05.18  we          sfEllInt: changed argument check in sfc_ell_rd
 1.47.14  04.05.18  we          sfEllInt: sfc_M_EllipticF/sfc_M_EllipticE
 1.47.15  05.05.18  we          sfEllInt: sfc_KleinJ

 1.48.00  13.05.18  we          sfEllInt: Basic lemniscatic case: separate function for argument
 1.48.01  14.05.18  we          sfEllInt: sfc_wpe_der / sfc_wpg_der
 1.48.02  20.05.18  we          sfEllInt: sfc_M_EllipticF := phi for m=0
 1.48.03  23.05.18  we          sfEllInt: sfc_M_EllipticPi
 1.48.04  07.06.18  we          sfGamma:  sfc_lnbinomial

 1.49.00  26.07.18  we          sfZeta:   Avoid D25 warning
 1.49.01  29.07.18  we          sfMisc:   sfc_exprel
 1.49.02  02.08.18  we          sfHyperG: Fix tmax update in hyp2f0
 1.49.03  03.08.18  we          sfHyperG: sfc_2f0

 1.50.00  14.08.18  we          sfEllint: IsNanOrInf check sfc_EllipticPiCim
 1.50.01  14.08.18  we          sfErf:    IsNanOrInf check sfc_erfi_inv
 1.50.02  19.08.18  we          sfBasic:  Make MaxBernoulli available for TPH
 1.50.03  21.08.18  we          sfMisc:   sfc_epoly
 1.50.04  23.08.18  we          sfGamma:  Check NaNs in sfc_gamma
 1.50.05  23.08.18  we          sfEllint: Implemented max iteration check in some functions
 1.50.06  23.08.18  we          sfExpInt: NanOrInf check in sfc_en
 1.50.07  23.08.18  we          sfZeta:   max iteration check in liae_pos
 1.50.08  23.08.18  we          sfBessel: NanOrInf check in sfc_jn, sfc_in
 1.50.09  23.08.18  we          sfGamma2: NanOrInf check in sfc_incgamma_ex, sfc_ibeta, sfc_nnbeta
 1.50.10  23.08.18  we                    Move Marcum Q from sfMisc to sfErf
 1.50.11  26.08.18  we          sfMisc:   sfc_besspoly
 1.50.12  28.08.18  we          sfGamma2: sfc_igammap_der
 1.50.13  28.08.18  we          SpecFun:  Fix return type of M_EllipticPiC
 1.50.14  01.09.18  we          sfGamma:  sfc_lnbg
 1.50.15  06.09.18  we          sfGamma2: sfc_igamman
 1.50.16  07.09.18  we          sfExpInt: sfc_en for n<0 and x>0
 1.50.17  07.09.18  we          sfExpInt: sfc_eibeta
 1.50.18  07.09.18  we          sfGamma2: improved sfc_igamma

 1.51.00  16.09.18  we          sfBessl2: New unit with Airy/Kelvin/Struve functions from sfBessel
 1.51.01  16.09.18  we          sfBessl2: sfc_coulcl
 1.51.02  17.09.18  we          sfBessl2: sfc_cshift
 1.51.03  28.09.18  we          sfBessl2: sfc_coul_ffp, sfc_coul_ggp, sfc_coul_f
 1.51.04  29.09.18  we          sfGamma2: sfc_igamma with Laguerre expansion
 1.51.05  01.10.18  we          sfMisc:   sfc_expn
 1.51.06  02.10.18  we          sfGamma2: removed sfc_igamman, sfc_igamma for x<0, a<>-1,-2,...
 1.51.07  02.10.18  we          sfGamma2: sfc_igammal for x < 0

 1.52.00  09.10.18  we          sfEllInt: IsNanOrInf check sfc_EllipticPiCim
 1.52.01  09.10.18  we          sfEllInt: Implemented max iteration check in some functions
 1.52.02  11.10.18  we          sfGamma:  sfc_psistar
 1.52.03  16.10.18  we          sfBessl2: sfc_synchf
 1.52.04  17.10.18  we          sfBessl2: sfc_synchg

 1.53.00  19.11.18  we          sfBessl2: Make skp4hx/sinkp4 global
 1.53.01  19.11.18  we          sfBessl2: sfc_kelvin_der: derivatives of all Kelvin functions
 1.53.02  20.11.18  we          sfBessl2: Separate functions sfc_berp/beip/kerp/keip
 1.53.03  21.11.18  we          sfBessl2: Fix NaN results in sfc_berp/beip/kerp/keip

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
  SF_BaseVersion = '1.53.03';

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

procedure sincosPix2(x: extended; var s,c: extended);
  {-Return s=sin(Pi/2*x^2), c=cos(Pi/2*x^2); (s,c)=(0,1) for abs(x) >= 2^64}

function sfc_bernoulli(n: integer): extended;
  {-Return the nth Bernoulli number, 0 if n<0 or odd n >= 3}

{#Z+}
{The next two functions are renamed copies from the amtools unit.}
{If $amtools_intde is defined the AMtools function are called.   }
{#Z-}

{---------------------------------------------------------------------------}
procedure sfc_intde_p(f: TFuncXP; p: pointer; a, b, eps: extended; var result, abserr: extended;
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
  { ier    - 0: OK, 3: max. iterations                     }
  {#F}

procedure sfc_intdei_p(f: TFuncXP; p: pointer; a, eps: extended; var result, abserr: extended;
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
  { ier    - 0: OK, 1: eps < eps_x/32                      }
  {          3: max. iterations, result/abserr have values }
  {#F}


const
  MaxBernoulli = 2312;  {Maximum index for Bernoulli numbers}

{#Z+}
{---------------------------------------------------------------------------}
{Moebius function æ(n) for small arguments}
const
  n_moeb = 95;
  moebius: array[1..n_moeb] of shortint = (
             1, -1, -1,  0, -1,  1, -1,  0,  0,  1,
            -1,  0, -1,  1,  1,  0, -1,  0, -1,  0,
             1,  1, -1,  0,  0,  1,  0,  0, -1, -1,
            -1,  0,  1,  1,  1,  0, -1,  1,  1,  0,
            -1, -1, -1,  0,  0,  1, -1,  0,  0,  0,
             1,  0, -1,  0,  1,  0,  1,  1, -1,  0,
            -1,  1,  0,  0,  1, -1, -1,  0,  1, -1,
            -1,  0, -1,  1,  0,  0,  1, -1, -1,  0,
             0,  1, -1,  0,  1,  1,  1,  0, -1,  0,
             1,  0,  1,  1,  1);


{---------------------------------------------------------------------------}
{Bernoulli numbers. The B(2n), n = 0..MaxB2nSmall are in AMath}
const
  NBoF  = 20;
  BoFHex: array[0..NBoF] of THexExtW = (      {Bernoulli(2k+2)/(2k+2)! }
            ($AAAB,$AAAA,$AAAA,$AAAA,$3FFB),  {+8.3333333333333333336E-2}
            ($B60B,$0B60,$60B6,$B60B,$BFF5),  {-1.3888888888888888888E-3}
            ($355E,$08AB,$55E0,$8AB3,$3FF0),  {+3.3068783068783068783E-5}
            ($5563,$A778,$BC99,$DDEB,$BFEA),  {-8.2671957671957671957E-7}
            ($ED15,$B875,$795F,$B354,$3FE5),  {+2.0876756987868098980E-8}
            ($F7BA,$2133,$2EB2,$9140,$BFE0),  {-5.2841901386874931847E-10}
            ($0FDC,$048B,$9627,$EB6D,$3FDA),  {+1.3382536530684678833E-11}
            ($E260,$CCAA,$70FB,$BED2,$BFD5),  {-3.3896802963225828668E-13}
            ($ECEE,$2974,$38EB,$9AAC,$3FD0),  {+8.5860620562778445639E-15}
            ($8664,$5567,$CB46,$FABE,$BFCA),  {-2.1748686985580618731E-16}
            ($77D8,$E1BB,$0F84,$CB3F,$3FC5),  {+5.5090028283602295151E-18}
            ($E371,$D15A,$C819,$A4BE,$BFC0),  {-1.3954464685812523341E-19}
            ($3204,$3568,$96BE,$8589,$3FBB),  {+3.5347070396294674718E-21}
            ($B3E6,$F8E0,$96AC,$D87B,$BFB5),  {-8.9535174270375468504E-23}
            ($6AD7,$6513,$6D26,$AF79,$3FB0),  {+2.2679524523376830603E-24}
            ($57E2,$6F1C,$EFCB,$8E3B,$BFAB),  {-5.7447906688722024451E-26}
            ($0D7F,$90DB,$CEF2,$E694,$3FA5),  {+1.4551724756148649018E-27}
            ($3676,$2BA5,$F32B,$BAE6,$BFA0),  {-3.6859949406653101781E-29}
            ($F634,$DA56,$46D9,$977F,$3F9B),  {+9.3367342570950446721E-31}
            ($F768,$3B0F,$1482,$F599,$BF95),  {-2.3650224157006299346E-32}
            ($2185,$9423,$FFC9,$C712,$3F90)); {+5.9906717624821343044E-34}
{#Z-}

{$ifdef debug}
type
  sfc_debug_str = string[255];

const
  NDCTR = 8;
var
  sfc_diagctr: array[0..NDCTR-1] of longint; {General counters for diagnostics}
  sfc_debug_output: boolean;                 {Really doing the outputs if true}

procedure sfc_dump_diagctr;
  {-Writeln diagnostic counters}

procedure sfc_write_debug_str(const msg: sfc_debug_str);
  {-Writeln or Outputdebugstr of msg}
{$endif}


implementation


{$undef non_empty}
{$undef USE_OutputDebugString}

{$ifdef amtools_intde}
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
    {$ifdef amtools_intde}
      AMTools,
    {$endif}
    {$ifdef UNIT_SCOPE}
      winapi.windows;
    {$else}
      windows;
    {$endif}
  {$else}
    {$ifdef amtools_intde}
      AMTools;
    {$endif}
  {$endif}
{$endif}



{---------------------------------------------------------------------------}
procedure sincosPix2(x: extended; var s,c: extended);
  {-Return s=sin(Pi/2*x^2), c=cos(Pi/2*x^2); (s,c)=(0,1) for abs(x) >= 2^64}
var
  n,f,g: extended;
begin
  {Note that this routine is very sensible to argument changes!}

  {Example: for x = 10000.1 the s value is 1.57073173006607621E-2, and the}
  {result for the same extended x in MPArith gives 1.57073173006607619E-2,}
  {the difference is 2.181E-19. If x is represented with MPArith's default}
  {240 bit precision the value is 1.57073173118206758E-2 with a difference}
  {of 1.115991E-11. The differences of the MPArith arguments were absolute}
  {10000.1 - x(extended) = 3.553E-16, relative 3.553E-20}

  if TExtRec(x).xp and $7FFF >= $403F then begin
    {abs(x) >= 2^64}
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
function sfc_bernoulli(n: integer): extended;
  {-Return the nth Bernoulli number, 0 if n<0 or odd n >= 3}
var
  bn,p4: extended;
  m: integer;
const
  bx32: array[0..68] of THexExtw = (        {Bernoulli(32*n+128)}
          ($52FD,$BE3D,$B49C,$DA4E,$C178),  {-5.2500923086774133900E+113}
          ($FF69,$87B9,$BE6D,$ABE4,$C209),  {-1.8437723552033869727E+157}
          ($4263,$EB29,$A09A,$A2CC,$C2A3),  {-3.9876744968232207445E+203}
          ($F6AA,$C1D8,$F6B1,$FC3C,$C344),  {-1.8059559586909309014E+252}
          ($B09F,$C9C2,$93E6,$9464,$C3ED),  {-7.9502125045885252855E+302}
          ($9CD8,$18C3,$6CDC,$9557,$C49B),  {-1.9158673530031573512E+355}
          ($4195,$04C1,$073E,$A49C,$C54E),  {-1.6181135520838065925E+409}
          ($739D,$571F,$11CD,$8B25,$C606),  {-3.3538244754399336890E+464}
          ($6270,$777F,$2738,$86CE,$C6C2),  {-1.2747336393845383643E+521}
          ($7FAB,$7570,$F9E6,$EADB,$C781),  {-6.9702671752326436791E+578}
          ($4BFF,$0F88,$1349,$95D1,$C845),  {-4.4656129117005935721E+637}
          ($77B2,$1772,$C601,$EAC7,$C90B),  {-2.8113943110814931668E+697}
          ($215A,$6E86,$91AA,$C206,$C9D5),  {-1.4934073418970347178E+758}
          ($5C7F,$D540,$F299,$9400,$CAA2),  {-5.8578844571843963827E+819}
          ($9A59,$4576,$14FE,$B949,$CB71),  {-1.5084075750891085972E+882}
          ($F479,$B007,$0BB5,$AB72,$CC43),  {-2.2966905497257900647E+945}
          ($C675,$9DA9,$A628,$D590,$CD17),  {-1.8830664597658071289E+1009}
          ($0CCE,$2F96,$B80B,$A49C,$CDEE),  {-7.6426951845750635249E+1073}
          ($3EA8,$63C6,$E69D,$9180,$CEC7),  {-1.4228762140674037692E+1139}
          ($4F15,$2C2F,$CD6C,$899F,$CFA2),  {-1.1338551314597730180E+1205}
          ($3CE2,$FDE0,$5490,$82C2,$D07F),  {-3.6304766452190330209E+1271}
          ($3E25,$B6E8,$D002,$EB8A,$D15D),  {-4.4077763139644032336E+1338}
          ($C753,$F062,$8250,$BEAA,$D23E),  {-1.9238585974016076227E+1406}
          ($0C77,$16EC,$4BAC,$840D,$D321),  {-2.8737803110109338078E+1474}
          ($4BEB,$3C03,$6AAB,$9587,$D405),  {-1.4036962824375857856E+1543}
          ($55C7,$0248,$8353,$84AE,$D4EB),  {-2.1491058465822269708E+1612}
          ($E368,$D032,$65E9,$B163,$D5D2),  {-9.9151825072869898178E+1681}
          ($BDEF,$D7FE,$04EF,$AC37,$D6BB),  {-1.3287264083350159101E+1752}
          ($7A90,$3C2E,$9735,$EA9B,$D7A5),  {-4.9971936871087436659E+1822}
          ($9FC1,$AC12,$A7FA,$D91F,$D891),  {-5.1070475539922579356E+1893}
          ($C234,$DE06,$003E,$8470,$D97F),  {-1.3759818684504019193E+1965}
          ($D312,$510C,$861B,$CEFB,$DA6D),  {-9.4989226808690414704E+2036}
          ($B32D,$A0BE,$8870,$C9B7,$DB5D),  {-1.6356184195123832511E+2109}
          ($F8D5,$2E47,$4E40,$EF06,$DC4E),  {-6.8487492292184253538E+2181}
          ($4706,$AFF3,$6F23,$A81A,$DD41),  {-6.8082216303292659152E+2254}
          ($269B,$331E,$4461,$892E,$DE35),  {-1.5706145196821570477E+2328}
          ($B8F7,$3321,$712A,$FE3E,$DF29),  {-8.2289799425426414984E+2401}
          ($D598,$89BA,$1E5D,$830E,$E020),  {-9.5930887251893861974E+2475}
          ($5173,$0867,$AB5F,$9368,$E117),  {-2.4402644753692194590E+2550}
          ($E63F,$76FF,$2870,$B191,$E20F),  {-1.3295801865055646274E+2625}
          ($F314,$7CE5,$BF4B,$E10C,$E308),  {-1.5244028786126305655E+2700}
          ($C94F,$DECD,$6EC9,$9389,$E403),  {-3.6161842428643845094E+2775}
          ($5940,$BF69,$792C,$C4E9,$E4FE),  {-1.7464299021830195980E+2851}
          ($A1C0,$F220,$16A1,$83B6,$E5FB),  {-1.6907945786261035527E+2927}
          ($92B8,$8BC5,$F001,$AE03,$E6F8),  {-3.2332881196060927594E+3003}
          ($DBF1,$FA98,$71AD,$DFDC,$E7F6),  {-1.2040771662431166519E+3080}
          ($BEF9,$BB8F,$DE3B,$8A4F,$E8F6),  {-8.6141885195778356538E+3156}
          ($D7D6,$F08B,$2222,$A20A,$E9F6),  {-1.1685695524501785594E+3234}
          ($E709,$C538,$6958,$B1BD,$EAF7),  {-2.9684329116433388662E+3311}
          ($9EFA,$77E4,$442C,$B459,$EBF9),  {-1.3950642936150073193E+3389}
          ($9693,$9589,$9A1E,$A753,$ECFC),  {-1.1989884572796487307E+3467}
          ($1F33,$FD26,$C024,$8C5F,$EE00),  {-1.8635271302518619006E+3545}
          ($BFDB,$D0DD,$0E46,$D2AF,$EF04),  {-5.1817804064721174490E+3623}
          ($36A4,$A1D6,$BAE4,$8BF7,$F00A),  {-2.5511413022051873656E+3702}
          ($04E7,$CCC3,$B596,$A2FF,$F110),  {-2.2016594557163485556E+3781}
          ($0B04,$CD7E,$45F3,$A4C4,$F217),  {-3.2985559030415466710E+3860}
          ($E49A,$0504,$EDB4,$8F39,$F31F),  {-8.4995513646331021388E+3939}
          ($9F07,$D98B,$CC09,$D433,$F427),  {-3.7328650389340701819E+4019}
          ($1417,$52A2,$A454,$84CC,$F531),  {-2.7699188001914063216E+4099}
          ($3D97,$75F1,$0EF2,$8B3C,$F63B),  {-3.4434756013646314132E+4179}
          ($198B,$2EA0,$521E,$F293,$F745),  {-7.1133726432191288073E+4259}
          ($A708,$C185,$31CE,$AE2D,$F851),  {-2.4224631302940113182E+4340}
          ($955A,$427F,$7027,$CC98,$F95D),  {-1.3495915388686739922E+4421}
          ($E7A1,$2BA2,$3BF4,$C31E,$FA6A),  {-1.2208793515089649122E+4502}
          ($B45B,$AA6B,$972C,$95FC,$FB78),  {-1.7804374121649736373E+4583}
          ($CB82,$33DD,$DA5B,$B88F,$FC86),  {-4.1563772250412736023E+4664}
          ($90AA,$6E95,$658B,$B48A,$FD95),  {-1.5426825919798643530E+4746}
          ($C03F,$1D53,$D309,$8B77,$FEA5),  {-9.0434733129342411538E+4827}
          ($3D11,$CC1D,$E67F,$A912,$FFB5)); {-8.3194707606651575345E+4909}
begin
  if odd(n) or (n<0) then begin
    if n=1 then sfc_bernoulli := -0.5
    else sfc_bernoulli := 0.0;
  end
  else begin
    m := n div 2;
    if m<=MaxB2nSmall then sfc_bernoulli := extended(B2nHex[m])
    else if n>MaxBernoulli then sfc_bernoulli := PosInf_x
    else begin
      {When n is even, B(2n) = -2*(-1)^n*m!/(2Pi)^m*zeta(m) with m=2n. For }
      {large m (e.g. m>63) zeta(m) is very close to 1 and we can derive the}
      {asymptotic recursion formula B(m+1) = -m*(m+1)/(2Pi)^2 * B(m). The  }
      {avg. iteration count is <4, the max. rel. error=4.5*eps_x for n=878.}
      m  := (n - 112) div 32;
      bn := extended(bx32[m]);
      m  := 32*m + 128;
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
      sfc_bernoulli := bn;
    end;
  end;
end;


{$ifdef amtools_intde}

{---------------------------------------------------------------------------}
procedure sfc_intde_p(f: TFuncXP; p: pointer; a, b, eps: extended; var result, abserr: extended;
                      var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over the finite interval (a,b)}
  { using Double Exponential (DE) transformation. Parameters:  }
begin
  amtools.intde_p(f, p, a, b, eps, result, abserr, neval, ier);
end;


{---------------------------------------------------------------------------}
procedure sfc_intdei_p(f: TFuncXP; p: pointer; a, eps: extended; var result, abserr: extended;
                       var neval: longint; var ier: integer);
  {-Automatic quadrature of f(x) over (a,INF) using Double Exponential}
  { transformation when f(x) has no oscillatory factor. Parameters:   }
begin
  amtools.intdei_p(f, p, a, eps, result, abserr, neval, ier);
end;

{$else}

{---------------------------------------------------------------------------}
procedure sfc_intde_p(f: TFuncXP; p: pointer; a, b, eps: extended; var result, abserr: extended;
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
  { ier    - 0: OK, 3: max. iterations                     }
  {#F}
const
  mmax = 256;
  efs  = 0.1;
  hoff = 8.5;
var
  m: integer;
var
  epsln,epsh,h0,ehp,ehm,epst,ba,ir,iv,h,iback,
  irback,t,ep,em,xw,xa,wg,fa,fb,err,errt,errh,errd: extended;
begin

  result := 0;
  abserr := 0;
  neval  := 0;
  ier    := 0;
  if a=b then exit;

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
procedure sfc_intdei_p(f: TFuncXP; p: pointer; a, eps: extended; var result, abserr: extended;
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
  { ier    - 0: OK, 1: eps < eps_x/32                      }
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
  t,ep,em,xp,xm,fp,fm,err,errt,errh,errd: extended;
begin

  result := 0.0;
  abserr := 0.0;
  neval  := 0;

  if eps < eps_x/32.0 then begin
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
  ir    := f(a+1.0, p);
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
        fp  := f(a+xp, p);
        fm  := f(a+xm, p);
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
procedure sfc_write_debug_str(const msg: sfc_debug_str);
  {-Writeln or Outputdebugstr of msg}
var
  ax: ansistring;
begin
  if sfc_debug_output then begin
    if IsConsole then writeln(msg)
    else begin
      ax := msg;
      OutputDebugString(pchar({$ifdef D12Plus}string{$endif}(ax)));
    end;
  end;
end;
{$else}
{---------------------------------------------------------------------------}
procedure sfc_write_debug_str(const msg: sfc_debug_str);
  {-Writeln or Outputdebugstr of msg}
begin
  if sfc_debug_output then writeln(msg);
end;
{$endif}


{---------------------------------------------------------------------------}
procedure sfc_dump_diagctr;
  {-Writeln diagnostic counters}
var
  k: integer;
begin
  write('Diag ctr:');
  for k:=0 to NDCTR-1 do write(k:3,':',sfc_diagctr[k]);
  writeln;
end;

begin
  sfc_debug_output := false;
  fillchar(sfc_diagctr, sizeof(sfc_diagctr),0);
{$endif}

end.
