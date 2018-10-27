unit sdBasic;

{Basic code for double precision special functions}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


uses DAMath;


(*************************************************************************

 DESCRIPTION   :  Basic code for double precision special functions:
                  Definitions, constants, etc.

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

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

 1.11.00  24.05.14  we          sdMisc:   sfd_ali functional inverse of sfc_li
 1.11.01  29.05.14  we          sdZeta:   Removed redundancy in sfc_harmonic2
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
 1.12.07  13.07.14  we          sdHyperG: sfc_pcfd,sfc_pcfu,sfc_pcfhh

 1.13.00  30.07.14  we          sdBasic:  moebius[1..83]
 1.13.01  30.07.14  we          sdZeta:   Improved and expanded primezeta sfd_pz
 1.13.02  06.08.14  we          sdGamma:  First check IsNaND(x) in sfd_lngammas
 1.13.03  06.08.14  we          sdGamma:  improved sfd_gdr_pos for small x
 1.13.04  11.08.14  we          sdBessel: Iv := PosInf_v if Kv,Kv1=0 in bessel_ik or if x >= IvMaxXH in sfc_iv
 1.13.05  11.08.14  we          sdGamma:  special case x=x+y in sfd_beta
 1.13.06  11.08.14  we          sdHyperG: fix: sfc_pcfd,sfc_pcfu,sfc_pcfhh use double
 1.13.07  17.08.14  we          sdMisc:   Catalan function sfd_catf
 1.13.08  19.08.14  we          sdHyperG: sfc_pcfv
 1.13.09  20.08.14  we          sdGamma:  sfd_gamma_delta_ratio returns 1 if x=x+d

 1.14.00  01.09.14  we          sdHyperG: sfd_chu for x<0
 1.14.01  04.09.14  we          sdHyperG: fix h1f1_tricomi (division by t=0)
 1.14.02  07.09.14  we          sdSDist:  sfd_logistic_inv uses logit function
 1.14.03  06.10.14  we          sdSDist:  sfd_logistic_cdf uses logistic function

 1.15.00  26.10.14  we          sdHyperG: No Kummer transformation in 1F1 if b=0,-1,-2,...
 1.15.01  29.10.14  we          sdHyperG: 1F1 = 1 + a/b*(exp(x)-1) for very small a,b
 1.15.02  07.11.14  we          sdBasic:  Small Bernoulli numbers B2nHex moved to DAMath

 1.16.00  15.12.14  we          sdGamma:  small changes in sfc_batemang
 1.16.01  16.12.14  we          sdGamma2: Check IsNanOrInf in sfd_ilng
 1.16.02  23.12.14  we          sdSDist:  sfd_invgamma_pdf/cdf/inv
 1.16.03  26.12.14  we          sdSDist:  logarithmic (series) distribution
 1.16.04  03.01.15  we          sdSDist:  sfd_wald_cdf/pdf/inv
 1.16.05  05.01.15  we          sdSDist:  Error if k<1 in sfc_ls_cdf/pmf
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
 1.18.06  10.06.15  we          sdBessel: rewrite of sfc_in: use IJ_series and CF1_I
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

***************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2009-2015 Wolfgang Ehrhardt

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
  SD_BaseVersion = '1.21.02';

const
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
  MaxBernoulli = 259;

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
  {-writeln diagnostic counters}

procedure sfd_write_debug_str(const msg: sfd_debug_str);
  {-Writeln or Outputdebugstr of msg}
{$endif}


implementation


{$ifdef Debug}
  {$ifdef WIN32or64}
    {$ifndef VirtualPascal}
      {$define USE_OutputDebugString}
    {$endif}
  {$endif}
  {$ifdef USE_OutputDebugString}
    {$ifdef UNIT_SCOPE}
      uses winapi.windows;
    {$else}
      uses windows;
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
  {-writeln diagnostic counters}
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
