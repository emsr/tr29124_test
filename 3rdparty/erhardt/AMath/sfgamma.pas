unit sfGamma;

{Common code for special functions: Gamma function and related}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Common code for special functions: Gamma function and related

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in amath_info.txt/references

                  [7] Cephes Mathematical Library, Version 2.8
                      http://www.moshier.net/#Cephes or http://www.netlib.org/cephes/
                 [14] SLATEC Common Mathematical Library, Version 4.1, July 1993
                      (general purpose mathematical and statistical routines written in Fortran 77)
                      http://www.netlib.org/slatec
                 [19] Boost C++ Libraries, Release 1.42.0, 2010.
                      http://www.boost.org/
                 [20] Special functions by Wayne Fullerton,
                      http://www.netlib.org/fn
                      Almost identical to the FNLIB subset of SLATEC [14]
                 [26] N.M. Temme, A Set of Algorithms for the Incomplete Gamma Functions,
                      Probability in the Engineering and Informational Sciences, 8 (1994), pp.291-307,
                      available as http://oai.cwi.nl/oai/asset/10080/10080A.pdf
                 [27] A.R. Didonato, A.H. Morris, Computation of the Incomplete Gamma Function Ratios
                      and their Inverse. ACM TOMS, Vol 12, No 4, Dec 1986, pp.377-393.
                      Fortran source: ACM TOMS 13 (1987), pp.318-319; available from
                      http://netlib.org/toms/654
                 [32] D.E. Knuth: The Art of computer programming; Volume 1, Fundamental Algorithms,
                      3rd ed., 1997; Volume 2, Seminumerical Algorithms, 3rd ed., 1998.
                 [50] A. Erdelyi et al., Higher Transcendental Functions Vol. I-III, California
                      Institute of Technology - Bateman Manuscript Project, 1953-1955,
                      Available via http://en.wikipedia.org/wiki/Bateman_Manuscript_Project
                 [78] V.S. Adamchik, Contributions to the Theory of the Barnes function,
                      2003, https://arxiv.org/pdf/math/0308086v1.pdf


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  17.08.10  W.Ehrhardt  Common version number after split
 1.00.01  17.08.10  we          Temme's gstar function
 1.00.02  23.08.10  we          sfc_incgamma, sfc_igammap, sfc_igammaq
 1.00.03  25.08.10  we          sfc_incgamma: cases when a is integer or half integer
 1.00.04  29.08.10  we          sfc_incgamma: dax calculated via Lanczos sum
 1.00.05  01.09.10  we          sfc_incgamma_ex, sfc_igprefix
 1.00.06  04.09.10  we          sfc_incgamma_inv
 1.00.07  05.09.10  we          sfc_incgamma_inv: use eps_d for p+q=1 check
 1.00.08  05.09.10  we          sfc_incgamma_inv: return Inf if p=1 and q=0
 1.00.09  06.09.10  we          sfc_incgamma_inv: Bug fixed code for DM Eq (34)
 1.00.10  08.09.10  we          Improved arg checking and NAN/INF handling
 1.00.11  10.09.10  we          sfc_trigamma
 1.00.12  11.09.10  we          ibeta_series with iteration limit and convergence check
 1.00.13  14.09.10  we          sfc_eta, sfc_etam1
 1.00.14  15.09.10  we          sfc_zetam1

 1.01.00  06.10.10  we          sfc_dfac
 1.01.01  18.10.10  we          Lanczos sum with and without expg scale
 1.01.02  19.10.10  we          sfc_gamma_delta_ratio, sfc_gamma_ratio, sfc_pochhammer
 1.01.03  23.10.10  we          sfc_zetah

 1.02.00  11.11.10  we          sfc_beta uses sfc_pochhammer for integer x+y=0
 1.02.01  11.11.10  we          Fixed missing check in sfc_lnbeta

 1.03.00  07.12.10  we          Fixed sfc_psi(+INF) = +INF
 1.03.01  13.12.10  we          sfc_trigamma with Hurwitz zeta
 1.03.02  14.12.10  we          sfc_tetragamma, sfc_pentagamma
 1.03.03  15.12.10  we          sfc_polygamma
 1.03.04  19.12.10  we          sfc_(ln)gamma return NAN/RTE if x is a non-positive integer
 1.03.05  29.12.10  we          improved sfc_dfac for a few small values to return integer
 1.03.06  30.12.10  we          sfc_binomial

 1.08.00  03.07.11  we          Asymptotic expansion for x>1024 in stirf
 1.08.01  12.07.11  we          polygam_negx: polygamma for x<0 and 1<=n<=10
 1.08.02  03.08.11  we          invgamma renamed to rgamma
 1.08.03  04.08.11  we          small simplification in sfc_gstar
 1.08.04  06.08.11  we          sfc_dfac for odd negative n, fac/dfac return INF if undefined
 1.08.05  08.08.11  we          const ETAEPS, improved etam1pos
 1.08.06  09.08.11  we          improved sfc_zetah

 1.09.00  14.09.11  we          improved lngamma_small
 1.09.01  14.09.11  we          sfc_lngamma1p
 1.09.02  17.09.11  we          lanczos interfaced
 1.09.03  17.09.11  we          sfc_ibetaprefix
 1.09.04  18.09.11  we          improved sfc_ibeta
 1.09.05  06.10.11  we          sfc_zetah for 0 <= s < 1

 1.10.00  24.11.11  we          sfc_lnfac
 1.10.01  15.12.11  we          sfc_igprefix for TP5 without Lanczos
 1.10.02  21.12.11  we          TP5 for sfc_zetah/sfc_zetam1
 1.10.03  22.12.11  we          etam1pos with exp7

 1.11.00  03.02.12  we          Moved zeta and eta functions to sfZeta

 1.12.00  21.03.12  we          sfc_gamma_delta_ratio for negative arguments
 1.12.01  21.03.12  we          improved sfc_pochhammer
 1.12.02  22.03.12  we          sfc_poch1
 1.12.03  09.04.12  we          fix special case if prefix is zero in sfc_incgamma_ex
 1.12.04  11.04.12  we          sfc_igamma (non-normalised upper incomplete gamma)
 1.12.05  18.04.12  we          sfc_lngammas
 1.12.06  22.04.12  we          sfc_gamma for x < -MAXGAMX
 1.12.07  22.04.12  we          sfc_igammal (non-normalised lower incomplete gamma)

 1.13.00  07.06.12  we          fix sfc_igprefix for a<1 and underflow of exp(-x)
 1.13.01  08.06.12  we          improved sfc_rgamma
 1.13.02  09.06.12  we          sfc_taylor
 1.13.03  28.06.12  we          igaml_mser
 1.13.04  30.06.12  we          sfc_git (Tricomi's incomplete gamma)
 1.13.05  01.07.12  we          check loss of accuracy in sfc_igammal

 1.14.00  20.07.12  we          removed special treatment for a=0 in sfc_poch

 1.15.00  15.02.13  we          improved sfc_git
 1.15.01  16.02.13  we          fix sfc_trigamma/tetragamma/pentagamma for very large arguments

 1.16.00  30.03.13  we          improved sfc_binomial

 1.17.00  15.04.13  we          improved sfc_beta
 1.17.01  16.04.13  we          ibeta_series with parameter normalised
 1.17.02  17.04.13  we          sfc_nnbeta
 1.17.03  21.04.13  we          lanczos_gm05 = lanczos_g - 0.5 in implementation
 1.17.04  27.04.13  we          polygam_negx for n=11,12

 1.18.00  11.05.13  we          improved sfc_polygamma for large order: e.g. (3000,200)
 1.18.01  19.05.13  we          Prevent some wrong compiler optimizations for div by 3
 1.18.02  26.05.13  we          sfc_nnbeta for a<=0 or b<=0

 1.20.00  27.07.13  we          sfc_git for x<0

 1.21.00  25.09.13  we          sfc_batemang

 1.25.00  16.04.14  we          Improved polygam_negx (larger n for x<0}

 1.27.00  15.06.14  we          Incomplete/inverse functions moved to sfGamma2
 1.27.01  15.06.14  we          const lanczos_gm05: extended

 1.28.00  05.08.14  we          First check IsNaN(x) in sfc_lngamma
 1.28.01  05.08.14  we          improved sfc_gdr_pos for small x
 1.28.02  10.08.14  we          fix sfc_lnbeta for VER50
 1.28.03  11.08.14  we          special case x=x+y in sfc_beta
 1.28.04  20.08.14  we          sfc_gamma_delta_ratio returns 1 if x=x+d

 1.31.00  14.12.14  we          small changes in sfc_batemang

 1.35.00  25.08.15  we          Use THREE & SIXX to avoid 'optimization'

 1.36.00  11.10.15  we          polygamma for x<0 uses Euler series for Pi*cot(Pi*x)

 1.40.00  29.06.17  we          Code removed for TP5-TP6, TPW1-D1

 1.48.00  07.06.18  we          sfc_lnbinomial

 1.50.00  23.08.18  we          Check NaNs in sfc_gamma
 1.50.01  01.09.18  we          sfc_lnbg

 1.52.00  11.10.18  we          sfc_psistar

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
  MAXGAMX = 1755.455;  {max. argument for sfc_gamma}
  XLGAE   = 8.0;       {threshold for asymptotic expansion of sfc_lngamma}

function sfc_lngcorr(x: extended): extended;
  {-Return lngamma correction term lngamma(x) - ((x-0.5)*ln(x) - x + ln(sqrt(2*Pi)), x>=8}

function sfc_rgamma(x: extended): extended;
  {-Return the reciprocal gamma function rgamma = 1/gamma(x)}

function sfc_binomial(n,k: integer; dbl: boolean): extended;
  {-Return the binomial coefficient 'n choose k'; set dbl=true, if used in SpecFun}

function sfc_lnbinomial(n,k: longint): extended;
  {-Return ln(binomial(n,k)), n >= k >= 0}

function sfc_fac(n: integer): extended;
  {-Return the factorial n!, n<MAXGAMX-1; INF if n<0}

function sfc_lnfac(n: longint): extended;
  {-Return ln(n!), INF if n<0}

function sfc_dfac(n: integer): extended;
  {-Return the double factorial n!!, n<=3209; INF for even n<0}

function sfc_gamma(x: extended): extended;
  {-Return gamma(x), |x| <= MAXGAMX; NAN/RTE if x is a non-positive integer}

function sfc_gamma_medium(x: extended): extended;
  {-Return gamma(x), |x| <= 13, x negative integer produces div by 0}

function sfc_gaminv_small(x: extended): extended;
  {- 1/gamma(x) for small argument values, |x| < 0.03125}

function sfc_gamma1pm1(x: extended): extended;
  {-Return gamma(1+x)-1 with increased accuracy for x near 0}

function sfc_gamma_delta_ratio(x,d: extended): extended;
  {-Return gamma(x)/gamma(x+d), accurate even for |d| << |x|}

function sfc_gamma_ratio(x,y: extended): extended;
  {-Return gamma(x)/gamma(y)}

function sfc_pochhammer(a,x: extended): extended;
  {-Return the Pochhammer symbol (a)_x = gamma(a+x)/gamma(a).}

function sfc_poch1(a,x: extended): extended;
  {-Return (pochhammer(a,x)-1)/x, psi(a) if x=0; accurate even for small |x|}

function sfc_lngamma(x: extended): extended;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, NAN/RTE if x is a non-positive integer}

function sfc_lngammas(x: extended; var s: integer): extended;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, s=-1,1 is the sign of gamma}

function sfc_lngamma1p(x: extended): extended;
  {-Return ln(|gamma(1+x)|) with increased accuracy for x near 0}

function sfc_signgamma(x: extended): extended;
  {-Return sign(gamma(x)), useless for 0 or negative integer}

function sfc_gstar(x: extended): extended;
  {-Return Temme's gstar(x) = gamma(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x)), x>0.}
  { For large x the asymptotic expansion is gstar(x) = 1 + 1/12x + O(1/x^2)}

function sfc_psi(x: extended): extended;
  {-Return the psi (digamma) function of x, INF if x is a non-positive integer}

function sfc_psistar(x: extended): extended;
  {-Return psi(x) - ln(x), x > 0}

function sfc_trigamma(x: extended): extended;
  {-Return the trigamma function psi'(x), INF if x is a negative integer}

function sfc_tetragamma(x: extended): extended;
  {-Return the tetragamma function psi''(x), NAN/RTE if x is a negative integer}

function sfc_pentagamma(x: extended): extended;
  {-Return the pentagamma function psi'''(x), INF if x is a negative integer}

function sfc_polygamma(n: integer; x: extended): extended;
  {-Return the polygamma function: n'th derivative of psi; n>=0, x>0 for n>MAXGAMX}
  { Note: The accuracy may be reduced for n >= MAXGAMX due to ln/exp operations.}

function sfc_batemang(x: extended): extended;
  {-Return the Bateman function G(x) = psi((x+1)/2) - psi(x/2), x<>0,-1,-2,...}

function sfc_beta(x,y: extended): extended;
  {-Return the function beta(x,y)=gamma(x)*gamma(y)/gamma(x+y)}

function sfc_lnbeta(x,y: extended): extended;
  {-Return the logarithm of |beta(x,y)|=|gamma(x)*gamma(y)/gamma(x+y)|}

function sfc_taylor(x: extended; n: integer): extended;
  {-Return the Taylor coefficient x^n/n!, n>=0}

function sfc_lnbg(x: extended): extended;
  {-Return ln(BarnesG(x)), real part for x < 0}

{#Z+}
const
  lanczos_gm05: extended = 12.64453125 + 3.375000000005456968e-5; {lanczos_g-0.5}
  {lanczos_g = asd(13.144565) = ($0000,$B257,$2363,$D250,$4002) = 13.1445650000000000545696821...}

function lanczos(x: extended; egscale: boolean): extended;
  {-Return the Lanczos sum for x, exp(g) scaled if egscale=true}
{#Z-}


implementation

uses
  AMath,
  memh,
  sfBasic,  {Basic common code}
  sfZeta;   {Zeta functions and polylogarithms}

{---------------------------------------------------------------------------}
function lanczos(x: extended; egscale: boolean): extended;
  {-Return the Lanczos sum for x, exp(g) scaled if egscale=true}
const
  lsh:  array[0..12] of THexExtW = (
          ($6F92,$D992,$7F81,$A01D,$402C),  {+4.4012138428004608955E+13}
          ($CD54,$4004,$20CB,$974E,$402C),  {+4.1590453358593200516E+13}
          ($F967,$596F,$659A,$8311,$402B),  {+1.8013842787117996778E+13}
          ($71D8,$A66C,$D4A5,$899F,$4029),  {+4.7287362634753888969E+12}
          ($96F3,$2C67,$5538,$C317,$4026),  {+8.3791008362840464705E+11}
          ($817F,$4B70,$3C24,$C4AA,$4023),  {+1.0558370727342993449E+11}
          ($125B,$9FAE,$C3E8,$908F,$4020),  {+9.7013636184949994935E+9}
          ($B883,$315C,$CD76,$9C24,$401C),  {+6.5491439754820526409E+8}
          ($29BA,$97D5,$7978,$F5F5,$4017),  {+3.2238322942133565306E+7}
          ($EA59,$87B0,$11C1,$89C2,$4013),  {+1.1285142194970914380E+6}
          ($0534,$EF56,$966A,$D053,$400D),  {+2.6665793784598589447E+4}
          ($8931,$781E,$A7EE,$BEF0,$4007),  {+3.8188012486329268705E+2}
          ($2CB3,$B138,$98FF,$A06C,$4000)); {+2.5066282746310005025}
  lseh: array[0..12] of THexExtW = (
          ($EBF0,$180B,$E131,$A434,$4019),  {+8.6091529534185372176E+7}
          ($A3D3,$B6ED,$E125,$9B2B,$4019),  {+8.1354505178580112428E+7}
          ($78B7,$D786,$C498,$866A,$4018),  {+3.5236626388154619108E+7}
          ($F7DC,$EF2B,$16FC,$8D24,$4016),  {+9.2498149880244712949E+6}
          ($8D43,$C678,$81BB,$C813,$4013),  {+1.6390242166871469602E+6}
          ($D9B5,$7AB8,$B435,$C9B0,$4010),  {+2.0653081576412250327E+5}
          ($A231,$1063,$6764,$9441,$400D),  {+1.8976701935302889156E+4}
          ($059E,$8F9A,$3482,$A022,$4009),  {+1.2810689099125594799E+3}
          ($35CD,$8C87,$6555,$FC3E,$4004),  {+6.3060933434202345361E+1}
          ($950B,$1B89,$3411,$8D47,$4000),  {+2.2074709097925276381}
          ($3037,$2E58,$56F1,$D5A6,$3FFA),  {+5.2160586946135054274E-2}
          ($653E,$8220,$AD06,$C3D1,$3FF4),  {+7.4699038089154483163E-4}
          ($6259,$720D,$001B,$A486,$3FED)); {+4.9031805734598718626E-6}
type
  TX12 = array[0..12] of extended;

  function lratev(const num: TX12; z: extended): extended;
  const
    denom: TX12 = ( 0.0,
                    39916800.0,
                    120543840.0,
                    150917976.0,
                    105258076.0,
                    45995730.0,
                    13339535.0,
                    2637558.0,
                    357423.0,
                    32670.0,
                    1925.0,
                    66.0,
                    1.0);
  var
    s1,s2: extended;
    i: integer;
  begin
    if abs(z)<=1.0 then begin
      s1 := num[12];
      s2 := denom[12];
      for i:=11 downto 0 do begin
        s1 := s1*z + num[i];
        s2 := s2*z + denom[i];
      end;
    end
    else begin
      z  := 1.0/z;
      s1 := num[0];
      s2 := denom[0];
      for i:=1 to 12 do begin
        s1 := s1*z + num[i];
        s2 := s2*z + denom[i];
      end;
    end;
    lratev := s1/s2;
  end;

begin
  {Ref: Boost [19], lanczos.hpp, struct lanczos13}
  if egscale then lanczos := lratev(TX12(lseh), x)
  else lanczos := lratev(TX12(lsh), x);
end;


{---------------------------------------------------------------------------}
function stirf(x: extended): extended;
  {-Stirling's formula for the gamma function for x > 13.0}
  { gamma(x) = sqrt(2*Pi)*x^(x-0.5)*exp(-x)*(1 + 1/x * P(1/x))}
const
  {13 <= x <= 1024, relative peak error = 9.44E-21, relative error spread = 8.8e-4}
  SPHex: array[0..8] of THexExtW = (
           ($a1d5,$aaaa,$aaaa,$aaaa,$3ffb),  { 8.333333333333331800504E-2}
           ($c3c9,$906e,$38e3,$e38e,$3ff6),  { 3.472222222230075327854E-3}
           ($3a1c,$5ac8,$3478,$afb9,$bff6),  {-2.681327161876304418288E-3}
           ($bef3,$7023,$6a08,$f09e,$bff2),  {-2.294719747873185405699E-4}
           ($30b7,$1a21,$98b2,$cd87,$3ff4),  { 7.840334842744753003862E-4}
           ($5704,$1a39,$b11d,$9293,$3ff1),  { 6.989332260623193171870E-5}
           ($ba6f,$7c59,$5e47,$9bfb,$bff4),  {-5.950237554056330156018E-4}
           ($c395,$0295,$4443,$c64b,$bfef),  {-2.363848809501759061727E-5}
           ($6ede,$69f7,$54e3,$bb5d,$3ff4)); { 7.147391378143610789273E-4}
const
  {x > 1024, rational coefficients from the analytical expansion  }
  { exp(sum(i=1,20,x^(2*i-1)*bernfrac(2*i)/(2*i)/(2*i-1)))+O(x^7) }
  { = 1 + 1/12*x + 1/288*x^2 - 139/51840*x^3 - 571/2488320*x^4    }
  {     + 163879/209018880*x^5 + 5246819/75246796800*x^6 + O(x^7) }
  SAHex: array[0..5] of THexExtW = (
           ($aaab,$aaaa,$aaaa,$aaaa,$3ffb),  { 8.33333333333333333333E-2} {1/12}
           ($e38e,$8e38,$38e3,$e38e,$3ff6),  { 3.47222222222222222222E-3} {1/288}
           ($3df2,$d5a6,$3476,$afb9,$bff6),  {-2.68132716049382716049E-3} {-139/51840}
           ($cab1,$fd42,$7232,$f09e,$bff2),  {-2.29472093621399176955E-4} {-571/2488320}
           ($20e4,$a796,$fb43,$cd87,$3ff4),  { 7.84039221720066627474E-4} {163879/209018880}
           ($c3f2,$ce01,$0241,$923b,$3ff1)); { 6.97281375836585777429E-5} {5246819/75246796800}
const
  XMAX = 1500.0; {max. argument for direct power}
var
  SP: array[0..8] of extended absolute SPHex;
  SA: array[0..5] of extended absolute SAHex;
var
  v,w,y: extended;
begin
  {Ref: Cephes [7], function stirf in file ldouble/gammal.c}
  w := 1.0/x;
  if x > 1024.0 then w := 1.0 + w*PolEvalX(w, SA, 6)
  else w := 1.0 + w*PolEvalX(w, SP, 9);
  if x > XMAX then begin
    {Avoid overflow in power}
    y := exp(-x);
    v := power(x,0.5*x-0.25);
    y := v*(v*y);
  end
  else begin
    y := exp(x);
    y := power(x, x-0.5)/y;
  end;
  stirf := Sqrt_TwoPi*y*w;
end;


{---------------------------------------------------------------------------}
function sfc_signgamma(x: extended): extended;
  {-Return sign(gamma(x)), useless for 0 or negative integer}
begin
  if IsNanOrInf(x) or (x>0.0) or (frac(0.5*floorx(x))=0.0) then sfc_signgamma := 1.0
  else sfc_signgamma := -1.0;
end;


{---------------------------------------------------------------------------}
function sfc_gaminv_small(x: extended): extended;
  {- 1/gamma(x) for small argument values, |x| < 0.03125}
const
  {1/gamma(x) = x*P(x), 0 < x < 0.03125, peak relative error 4.2e-23}
  SPHex: array[0..8] of THexExtW = (
           ($0000,$0000,$0000,$8000,$3fff),   { 1.000000000000000000000E+0}
           ($c7a9,$7db0,$67e3,$93c4,$3ffe),   { 5.772156649015328608253E-1}
           ($7bf6,$57d1,$a013,$a7e7,$bffe),   {-6.558780715202540684668E-1}
           ($f183,$126b,$f47d,$ac0a,$bffa),   {-4.200263503403344054473E-2}
           ($6b8d,$7515,$1905,$aa89,$3ffc),   { 1.665386113720805206758E-1}
           ($10b0,$ec17,$87dc,$acd7,$bffa),   {-4.219773360705915470089E-2}
           ($9225,$dfef,$b0e9,$9da5,$bff8),   {-9.622023360406271645744E-3}
           ($fe9a,$ceb4,$c74e,$ec9a,$3ff7),   { 7.220599478036909672331E-3}
           ($baeb,$d6d3,$25e5,$9c7e,$bff5));  {-1.193945051381510095614E-3}

  {1/gamma(-x) = x*P(x), 0 < x < 0.03125, peak relative error 5.16e-23}
  SNHex: array[0..8] of THexExtW = (
           ($0000,$0000,$0000,$8000,$bfff),   {-1.000000000000000000000E+0}
           ($c7aa,$7db0,$67e3,$93c4,$3ffe),   { 5.772156649015328608727E-1}
           ($5e26,$57d1,$a013,$a7e7,$3ffe),   { 6.558780715202536547116E-1}
           ($7f64,$1234,$f47d,$ac0a,$bffa),   {-4.200263503402112910504E-2}
           ($7a5b,$d76d,$1905,$aa89,$bffc),   {-1.665386113944413519335E-1}
           ($783f,$41dd,$87d1,$acd7,$bffa),   {-4.219773343731191721664E-2}
           ($2ca1,$18f0,$386f,$9da5,$3ff8),   { 9.621911155035976733706E-3}
           ($989b,$dd68,$c5f1,$ec9c,$3ff7),   { 7.220837261893170325704E-3}
           ($5dd1,$02de,$b9f7,$948d,$3ff5));  { 1.133374167243894382010E-3}
var
  SP: array[0..8] of extended absolute SPHex;
  SN: array[0..8] of extended absolute SNHex;
var
  p: extended;
begin
  {Ref: Cephes [7], function gammal/label small in file ldouble/gammal.c}
  if x=0.0 then sfc_gaminv_small := 0.0
  else begin
    if x<0.0 then begin
      x := -x;
      p := PolEvalX(x, SN, 9);
    end
    else p := PolEvalX(x, SP, 9);
    sfc_gaminv_small := x*p;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_gamma_medium(x: extended): extended;
  {-Return gamma(x), |x| <= 13, x negative integer produces div by 0}
const
  {gamma(x+2) = P(x)/Q(x), 0 <= x <= 1, peak relative error = 1.83e-20}
  PHex: array[0..7] of THexExtW = (
           ($0000,$0000,$0000,$8000,$3fff),   { 1.000000000000000000009E+0}
           ($29cf,$19b3,$16c8,$d67a,$3ffe),   { 8.378004301573126728826E-1}
           ($8d75,$23af,$c8e4,$b9d4,$3ffd),   { 3.629515436640239168939E-1}
           ($9549,$8eb5,$8c3a,$e3f4,$3ffb),   { 1.113062816019361559013E-1}
           ($7f43,$5196,$b166,$c368,$3ff9),   { 2.385363243461108252554E-2}
           ($be6c,$3757,$c717,$861b,$3ff7),   { 4.092666828394035500949E-3}
           ($f5aa,$e82f,$335b,$ee2e,$3ff3),   { 4.542931960608009155600E-4}
           ($434a,$3f22,$2bda,$b0b2,$3ff0));  { 4.212760487471622013093E-5}

  QHex: array[0..8] of THexExtW = (
           ($0000,$0000,$0000,$8000,$3fff),   { 9.999999999999999999908E-1}
           ($e458,$2ec7,$fd57,$d47c,$3ffd),   { 4.150160950588455434583E-1}
           ($75ef,$3ab7,$4ad3,$e5bc,$bffc),   {-2.243510905670329164562E-1}
           ($3295,$3698,$d580,$bdcd,$bffa),   {-4.633887671244534213831E-2}
           ($0417,$7989,$d7bc,$e338,$3ff9),   { 2.773706565840072979165E-2}
           ($296e,$7cb1,$5dfd,$d08f,$bff4),   {-7.955933682494738320586E-4}
           ($beed,$1853,$a691,$a23d,$bff5),   {-1.237799246653152231188E-3}
           ($334b,$c2f0,$a2dd,$f60e,$3ff2),   { 2.346584059160635244282E-4}
           ($5473,$2de8,$1268,$ea67,$bfee));  {-1.397148517476170440917E-5}
var
  PX: array[0..7] of extended absolute PHex;
  QX: array[0..8] of extended absolute QHex;
var
  z: extended;
begin
  {Ref: Cephes [7], function gammal in file ldouble/gammal.c}
  {-13 <= x <= 13. Use recurrence formula to bring argument to [2,3)}
  z := 1.0;
  while x >= 3.0 do begin
    x := x-1.0;
    z := z*x;
  end;
  while x < -0.03125 do begin
    z := z/x;
    x := x+1.0;
  end;
  {Here -0.03125 <= x < 3}
  if x <= 0.03125 then begin
    {argument near a pole, use approximation of 1/gamma}
     sfc_gamma_medium := z/sfc_gaminv_small(x);
  end
  else begin
    {finish reduction to [2,3)}
    while x < 2.0 do begin
      z := z/x;
      x := x+1.0;
    end;
    if x=2.0 then  sfc_gamma_medium := z
    else begin
      x := x-2.0;
      sfc_gamma_medium := z * PolEvalX(x,PX,8) / PolEvalX(x,QX,9);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_lngcorr(x: extended): extended;
  {-Return lngamma correction term lngamma(x) - ((x-0.5)*ln(x) - x + ln(sqrt(2*Pi)), x>=8}
const
  { ln gamma(x) = (x-0.5)*ln(x) - x + ln(sqrt(2*Pi) + 1/x * A(1/x^2) }
  { x >= 8, peak relative error 1.51e-21}
  AHex: array[0..6] of THexExtW = (
          ($9fcc,$aaaa,$aaaa,$aaaa,$3ffb),   { 8.333333333333331447505E-2}
          ($4d88,$03a8,$60b6,$b60b,$bff6),   {-2.777777777750349603440E-3}
          ($f8f2,$30e5,$0092,$d00d,$3ff4),   { 7.936507795855070755671E-4}
          ($8b20,$9fce,$844e,$9c09,$bff4),   {-5.952345851765688514613E-4}
          ($3bdc,$aad1,$d492,$dc88,$3ff4),   { 8.412723297322498080632E-4}
          ($3d91,$0304,$3da1,$f685,$bff5),   {-1.880801938119376907179E-3}
          ($d984,$cc08,$91c2,$a012,$3ff7));  { 4.885026142432270781165E-3}
var
  PA: array[0..6] of extended absolute AHex;
const
  xbig = 4.294967296e9;
  xmax = 2.47860728199e4930;
begin
  {Ref: Cephes [7], part function lgaml in file ldouble/gammal.c}
  if x < XLGAE then begin
    {this is not really used but return 'correct' result}
    sfc_lngcorr := sfc_lngamma(x) - ((x-0.5)*ln(x) - x + LnSqrt2Pi);
  end
  else if x > xmax then sfc_lngcorr := 0.0
  else if x > xbig then sfc_lngcorr := 1.0/(12.0*x)
  else sfc_lngcorr := PolEvalX(1.0/(x*x),PA,7)/x;
end;


{---------------------------------------------------------------------------}
function sfc_gamma(x: extended): extended;
  {-Return gamma(x), x <= MAXGAMX; NAN/RTE if x is a non-positive integer}
var
  i: integer;
  q,z: extended;
begin
  {Ref: Cephes [7], function gammal in file ldouble/gammal.c}
  if IsNanOrInf(x) then begin
    sfc_gamma := x;
    exit;
  end;

  if (x<=0.0) and (frac(x)=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_gamma := NaN_x;
    exit;
  end;

  if x>MAXGAMX then begin
    sfc_gamma := PosInf_x;
    exit;
  end;

  q := abs(x);
  if q <= 13.0 then sfc_gamma := sfc_gamma_medium(x)
  else begin
    if x < 0.0 then begin
      if q <= MAXGAMX then begin
        i := trunc(q);
        z := q-i;
        if z > 0.5 then z := q - (i+1);
        z := abs(q * sinPi(z)) * stirf(q);
        if z <= Pi/MaxExtended then begin
          sfc_gamma := copysign(PosInf_x, sfc_signgamma(x));
          exit;
        end;
        sfc_gamma := sfc_signgamma(x)*Pi/z;
      end
      else begin
        {Use lngamma, there are only very rare 'non-zero' cases, most values}
        {underflow. gamma(-1756-0.5^24) = -2.75052000244830115684312e-4930}
        z := exp(sfc_lngammas(x,i));
        sfc_gamma := z*i;
      end;
    end
    else sfc_gamma := stirf(x);
  end
end;


{---------------------------------------------------------------------------}
function lngamma_small(x,xm1,xm2: extended; useln1p: boolean): extended;
  {-Return ln(gamma(x), 0<=x<8; xm1=x-1, xm2=x-2 are supplied for increased precision}
  { useln1p should be true if xm1 is supposed to be more accurate than x-1}
var
  r,t,z: extended;
const
  Y20: THexExtW = ($0,$0,$6000,$A2C7,$3FFC); {0.158963680267333984375 = 333371/2097152}
  P20: array[0..6] of extended = (
         -0.180355685678449379109e-1,
          0.25126649619989678683e-1,
          0.494103151567532234274e-1,
          0.172491608709613993966e-1,
         -0.259453563205438108893e-3,
         -0.541009869215204396339e-3,
         -0.324588649825948492091e-4);
  Q20: array[0..7] of extended = (
          1.0,
          0.196202987197795200688e1,
          0.148019669424231326694e1,
          0.541391432071720958364e0,
          0.988504251128010129477e-1,
          0.82130967464889339326e-2,
          0.224936291922115757597e-3,
         -0.223352763208617092964e-6);
const
  Y10: THexExtW = ($0,$0,$1000,$8735,$3FFE);  {0.52815341949462890625 = 1107618/2097152}
  P10: array[0..6] of extended = (
          0.490622454069039543534e-1,
         -0.969117530159521214579e-1,
         -0.414983358359495381969e0,
         -0.406567124211938417342e0,
         -0.158413586390692192217e0,
         -0.240149820648571559892e-1,
         -0.100346687696279557415e-2);
  Q10: array[0..6] of extended = (
          1.0,
          0.302349829846463038743e1,
          0.348739585360723852576e1,
          0.191415588274426679201e1,
          0.507137738614363510846e0,
          0.577039722690451849648e-1,
          0.195768102601107189171e-2);
const
  Y15: THexExtW = ($0,$0,$D000,$E76E,$3FFD);  {0.452017307281494140625 = 947949/2097152}
  P15: array[0..5] of extended = (
         -0.292329721830270012337e-1,
          0.144216267757192309184e0,
         -0.142440390738631274135e0,
          0.542809694055053558157e-1,
         -0.850535976868336437746e-2,
          0.431171342679297331241e-3);
  Q15: array[0..6] of extended = (
          1.0,
         -0.150169356054485044494e1,
          0.846973248876495016101e0,
         -0.220095151814995745555e0,
          0.25582797155975869989e-1,
         -0.100666795539143372762e-2,
         -0.827193521891290553639e-6);
begin
  {Based on \boost\math\special_functions\detail\lgamma_small.hpp [19]}
  {Copyright John Maddock 2006, see 3rdparty.ama for Boost license}
  if x<eps_x then begin
    lngamma_small := -ln(x);
    exit;
  end;
  if (xm1=0.0) or (xm2=0.0) then begin
    lngamma_small := 0.0;
    exit;
  end;
  {The Boost rational approximations are optimized for low absolute error.}
  {As long as their absolute error is small compared to the constants Yxx }
  {then any rounding errors in their computation will get wiped out.      }
  if x>2.0 then begin
    {Use recurrence formula to bring argument to [2,3)}
    if x>=3 then begin
      z := 1.0;
      repeat
        x := x-1.0;
        z := z*x;
      until x<3.0;
      xm2 := x-2.0;
      z := ln(z);
    end
    else z := 0.0;
    {Use the following form: lngamma(x) = (x-2)*(x+1)*(Y20 + R(x-2))}
    r := PolEvalX(xm2, P20, 7) / PolEvalX(xm2, Q20, 8);
    t := xm2*(x+1.0);
    lngamma_small := z + (t*extended(Y20) + t*r);
  end
  else begin
    if x<1.0 then begin
      {bring argument into [1,2)}
      if useln1p then z := -ln1p(xm1)
      else z := -ln(x);
      xm2 := xm1;
      xm1 := x;
      x   := x+1.0;
    end
    else z := 0.0;
    {Two approximations, one for x in [1.0,1.5] and one for x in (1.5,2]}
    if x<=1.5 then begin
      {Use the following form: lngamma(x) = (x-1)*(x-2)*(Y10 + R(x-1))}
      r := PolEvalX(xm1, P10, 7) / PolEvalX(xm1, Q10, 7);
      t := xm1*xm2;
      lngamma_small := z + (t*extended(Y10) + t*r);
    end
    else begin
      {Use the following form: lngamma(x) = (2-x)*(1-x)*(Y15 + R(2-x))}
      r := PolEvalX(-xm2, P15, 6) / PolEvalX(-xm2, Q15, 7);
      t := xm1*xm2;
      lngamma_small := z + (t*extended(Y15) + t*r);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_lngamma(x: extended): extended;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, NAN/RTE if x is a non-positive integer}
var
  t: extended;
begin
  if IsNaN(x) then sfc_lngamma := NaN_x
  else if x<=0.0 then begin
    if frac(x)=0.0 then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfc_lngamma := NaN_x;
    end
    else begin
      {use reflection formula}
      t := x*sinPi(x);
      sfc_lngamma := ln(abs(Pi/t)) - sfc_lngamma(abs(x));
    end;
    exit;
  end
  else if x <= XLGAE then begin
    {Note: will generate error for x=0}
    sfc_lngamma := lngamma_small(x,x-1.0,x-2.0,false)
  end
  else begin
    {Stirling approximation for x > 8}
    sfc_lngamma := (x-0.5)*ln(x) - x + LnSqrt2Pi + sfc_lngcorr(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_lngammas(x: extended; var s: integer): extended;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, s=-1,1 is the sign of gamma}
begin
  if IsNaN(x) then begin
    sfc_lngammas := NaN_x;
    s := 0;
    exit;
  end;
  sfc_lngammas := sfc_lngamma(x);
  if (x>0.0) or (frac(0.5*floorx(x))=0.0) then s := 1 else s := -1;
end;


{---------------------------------------------------------------------------}
function sfc_lngamma1p(x: extended): extended;
  {-Return ln(|gamma(1+x)|) with increased accuracy for x near 0}
var
  t: extended;
begin
  t := 1.0+x;
  if (t < 0.0) or (t > XLGAE) then sfc_lngamma1p := sfc_lngamma(t)
  else sfc_lngamma1p := lngamma_small(t,x,x-1.0,true)
end;


{---------------------------------------------------------------------------}
function sfc_gamma1pm1(x: extended): extended;
  {-Return gamma(1+x)-1 with increased accuracy for x near 0}
begin
  if abs(x)<=eps_x then sfc_gamma1pm1 := -EulerGamma*x
  else if (x<-0.5) or (x>2.0) then sfc_gamma1pm1 := sfc_gamma(1.0+x)-1.0
  else if x>0.0 then sfc_gamma1pm1 := expm1(lngamma_small(x+1.0,x,x-1.0,true))
  else sfc_gamma1pm1 := expm1(lngamma_small(x+2.0,x+1.0,x,false) - ln1p(x));
end;


{---------------------------------------------------------------------------}
function sfc_rgamma(x: extended): extended;
  {-Return the reciprocal gamma function rgamma = 1/gamma(x)}
begin
  if IsNanOrInf(x) then begin
    if x=PosInf_x then sfc_rgamma := 0.0
    else sfc_rgamma := NaN_x;
  end
  else if abs(x)<0.03125 then sfc_rgamma := sfc_gaminv_small(x)
  else if (x<0.0) and (frac(x)=0.0) then sfc_rgamma := 0.0
  else begin
    {
    if (x<-0.5) or (x>XLGAE) then sfc_rgamma := sfc_signgamma(x)*exp(-sfc_lngamma(x))
    else sfc_rgamma := 1.0/sfc_gamma_medium(x);
    }
    if (x<-0.5) or (x>=MAXGAMX) then sfc_rgamma := sfc_signgamma(x)*exp(-sfc_lngamma(x))
    else if x <= 13.0 then sfc_rgamma := 1.0/sfc_gamma_medium(x)
    else sfc_rgamma := 1.0/stirf(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_gstar(x: extended): extended;
  {-Return Temme's gstar(x) = gamma(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x)), x>0.}
  { For large x the asymptotic expansion is gstar(x) = 1 + 1/12x + O(1/x^2)}
var
  a,g: extended;
const
  {Chebyshev approximations calculated with Maple V R4 with Digits :=30;}
  {Coefficients for T(0,.) doubled, HexExtended with t_rcalc command xh.}
  {Definition: gstar := x-> GAMMA(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x));}

  {Range 1<=x<= 3; y := t -> t+2;}
  {chebyshev(y(t)^2*(gstar(y(t))-1-1/(12*y(t))), t=-1..1, 0.5e-20);}
  gs1n = 30;
  gs1h : array[0..gs1n-1] of THexExtW = (
           ($26A1,$14BA,$EC57,$83AB,$3FF7),  { 4.0182975773513416011E-3}
           ($8A43,$A98D,$2BD2,$B827,$3FF4),  { 7.0248799299004439952E-4}
           ($81EF,$77ED,$7E07,$A731,$BFF2),  {-1.5944798403032040299E-4}
           ($8F60,$F4D0,$DFFD,$900D,$3FF0),  { 3.4345197462305170773E-5}
           ($C34A,$42D7,$2A10,$EB6B,$BFED),  {-7.0160213288181423859E-6}
           ($7A47,$6EED,$A8E0,$B53C,$3FEB),  { 1.3503205169277659442E-6}
           ($7D24,$93CC,$D582,$8133,$BFE9),  {-2.4065836641022231323E-7}
           ($AFE6,$97BF,$5311,$A364,$3FE6),  { 3.8042639505142093235E-8}
           ($63DE,$2861,$E471,$A014,$BFE3),  {-4.6589880584690519361E-9}
           ($2FEC,$9CDC,$9ED8,$94B3,$3FDE),  { 1.3524335604397218593E-10}
           ($543F,$B1E8,$647F,$C44B,$3FDE),  { 1.7852880977454193496E-10}
           ($B937,$7641,$3927,$C00A,$BFDD),  {-8.7329651525165672802E-11}
           ($168B,$D18D,$7F0B,$84E5,$3FDC),  { 3.0217158796552201555E-11}
           ($1576,$BF4D,$8C8F,$A133,$BFDA),  {-9.1632366291192802894E-12}
           ($F8A3,$ABB3,$AEA7,$B66A,$3FD8),  { 2.5922976124555003060E-12}
           ($9848,$57E9,$3BE2,$C5E9,$BFD6),  {-7.0312136560940930257E-13}
           ($B155,$A152,$D1F5,$D0E3,$3FD4),  { 1.8553152123744796025E-13}
           ($DE6F,$8669,$0937,$D861,$BFD2),  {-4.8045799980445695572E-14}
           ($9337,$6A1D,$6404,$DD3A,$3FD0),  { 1.2280625885209937902E-14}
           ($91A1,$6ACB,$2EC2,$E020,$BFCE),  {-3.1103690941234212527E-15}
           ($4020,$F48D,$C5F0,$E19E,$3FCC),  { 7.8277734240176922378E-16}
           ($9D0A,$0C51,$259C,$E224,$BFCA),  {-1.9614622328820518705E-16}
           ($6E3B,$3B7E,$7731,$E205,$3FC8),  { 4.9010567735366476461E-17}
           ($8153,$BE84,$AC06,$E183,$BFC6),  {-1.2225157066341578107E-17}
           ($5EFA,$2FDE,$EA16,$E0CE,$3FC4),  { 3.0467200401436511414E-18}
           ($6370,$29B3,$9E32,$E013,$BFC2),  {-7.5920116218493385258E-19}
           ($6F91,$69A3,$68EE,$DF41,$3FC0),  { 1.8910477040223097493E-19}
           ($80CF,$9130,$ABE5,$DF35,$BFBE),  {-4.7266482810237095702E-20}
           ($6062,$3806,$EC04,$DCC2,$3FBC),  { 1.1687011937049585628E-20}
           ($01FE,$0907,$EE5B,$E372,$BFBA)); {-3.0102617821884302625E-21}
var
  gs1a: array[0..gs1n-1] of extended absolute gs1h;

const
  gs2n = 31; {chebyshev(gstar(t), t=3..8, 0.5e-20);}
  gs2h : array[0..gs2n-1] of THexExtW = (
           ($78D4,$FA06,$AD41,$8231,$4000),  { 2.0342820305159885945}
           ($D94B,$BBBA,$AB78,$87C7,$BFF8),  {-8.2873510863767290591E-3}
           ($D710,$9AB3,$309C,$831F,$3FF6),  { 2.0007604294772196608E-3}
           ($664B,$4102,$BD79,$FCED,$BFF3),  {-4.8242315747915459550E-4}
           ($6E6D,$ABF2,$018C,$F3A7,$3FF1),  { 1.1618250245736783082E-4}
           ($BAD0,$512E,$362E,$EA73,$BFEF),  {-2.7948623357453920514E-5}
           ($427A,$FE0E,$B94A,$E15A,$3FED),  { 6.7160841776690551181E-6}
           ($7FC4,$E615,$27FF,$D865,$BFEB),  {-1.6122694432402287168E-6}
           ($B1F8,$8BDA,$4FF3,$CF99,$3FE9),  { 3.8668303975989479479E-7}
           ($9AAE,$587F,$2775,$C6FD,$BFE7),  {-9.2661419574693789809E-8}
           ($0655,$F513,$CBA1,$BE95,$3FE5),  { 2.2187030226642767308E-8}
           ($418D,$BD3D,$83AE,$B667,$BFE3),  {-5.3086653662215709374E-9}
           ($51AD,$FD6D,$C872,$AE75,$3FE1),  { 1.2693642192007070286E-9}
           ($2E7A,$45B5,$4F47,$A6C3,$BFDF),  {-3.0333999974117580609E-10}
           ($FEA3,$E5E7,$175E,$9F52,$3FDD),  { 7.2450652218555779523E-11}
           ($7398,$B75A,$78D2,$9823,$BFDB),  {-1.7296152050716179453E-11}
           ($957B,$0418,$34AC,$9138,$3FD9),  { 4.1273879591581087674E-12}
           ($43E9,$5D5B,$8558,$8A90,$BFD7),  {-9.8456023555037793839E-13}
           ($DBD4,$B67F,$2F09,$842C,$3FD5),  { 2.3478568909854971813E-13}
           ($DA03,$7484,$1F3F,$FC15,$BFD2),  {-5.5973560909513051708E-14}
           ($8336,$2310,$5B9E,$F055,$3FD0),  { 1.3341185335663677403E-14}
           ($3600,$2A12,$8EA1,$E516,$BFCE),  {-3.1792362333608899741E-15}
           ($D32F,$A155,$BBE9,$DA55,$3FCC),  { 7.5750134834532103750E-16}
           ($7F27,$7DFA,$870F,$D00F,$BFCA),  {-1.8046385098946929066E-16}
           ($7E48,$8F57,$48FA,$C640,$3FC8),  { 4.2988857605477843878E-17}
           ($8495,$0C96,$22D2,$BCE4,$BFC6),  {-1.0239810102203542236E-17}
           ($5E48,$E817,$0E96,$B3F7,$3FC4),  { 2.4389814490930843665E-18}
           ($D507,$1A15,$ED88,$AB74,$BFC2),  {-5.8091806323722222819E-19}
           ($1C8E,$8568,$9477,$A359,$3FC0),  { 1.3836276558996084587E-19}
           ($4846,$9E37,$D63A,$9BA0,$BFBE),  {-3.2955567815269497085E-20}
           ($2038,$2762,$8C64,$9446,$3FBC)); { 7.8496438299511273955E-21}
var
  gs2a: array[0..gs2n-1] of extended absolute gs2h;
const
  gstar1: THexExtW = ($F232,$F32F,$D984,$8ACE,$3FFF);  {1.0844375514192275466}
begin
  if x>=1.79e8 then begin
    {very large x, use first terms of Stirling formula}
    a := 1.0/(12.0*x);
    sfc_gstar := 1.0 + a*(1.0 + 0.5*a);
  end
  else if x>8.0 then begin
    {medium range, recycle sfc_lngcorr}
    sfc_gstar := exp(sfc_lngcorr(x));
  end
  else if x>=3.0 then begin
    {3 <= x < 8}
    sfc_gstar := CSEvalX((2.0*x-11.0)/5.0, gs2a, gs2n);
  end
  else if x>=1.0 then begin
    {1 <= x < 3}
    a :=  CSEvalX(x-2.0, gs1a, gs1n);
    sfc_gstar := a/sqr(x) + 1.0 + 1.0/(12.0*x);
  end
  else if x>=MinExtended then begin
    {0 < x < 1: Use Temme's [26] recursion formula for gstar(x+1),}
    {but avoid recursive call if x+1 = 1 within machine precision.}
    g := 1.0 + x;
    if g=1.0 then g := extended(gstar1)
    else g := sfc_gstar(g);
    a := 1.0 + 1.0/x;
    g := g*sqrt(a);
    a := ln(a)*x - 1.0;
    sfc_gstar := g*exp(a);
  end
  else sfc_gstar := 1.0/(Sqrt_TwoPi*sqrt(x));
end;


{---------------------------------------------------------------------------}
const
  FTabH: array[0..25] of THexExtW = (        {Table of factorials}
           ($0000,$0000,$0000,$8000,$3FFF),  {1.000000000000000000000}
           ($0000,$0000,$0000,$8000,$3FFF),  {1.000000000000000000000}
           ($0000,$0000,$0000,$8000,$4000),  {2.000000000000000000000}
           ($0000,$0000,$0000,$C000,$4001),  {6.000000000000000000000}
           ($0000,$0000,$0000,$C000,$4003),  {2.400000000000000000000E+1}
           ($0000,$0000,$0000,$F000,$4005),  {1.200000000000000000000E+2}
           ($0000,$0000,$0000,$B400,$4008),  {7.200000000000000000000E+2}
           ($0000,$0000,$0000,$9D80,$400B),  {5.040000000000000000000E+3}
           ($0000,$0000,$0000,$9D80,$400E),  {4.032000000000000000000E+4}
           ($0000,$0000,$0000,$B130,$4011),  {3.628800000000000000000E+5}
           ($0000,$0000,$0000,$DD7C,$4014),  {3.628800000000000000000E+6}
           ($0000,$0000,$4000,$9845,$4018),  {3.991680000000000000000E+7}
           ($0000,$0000,$E000,$E467,$401B),  {4.790016000000000000000E+8}
           ($0000,$0000,$6600,$B994,$401F),  {6.227020800000000000000E+9}
           ($0000,$0000,$D940,$A261,$4023),  {8.717829120000000000000E+10}
           ($0000,$0000,$BBAC,$983B,$4027),  {1.307674368000000000000E+12}
           ($0000,$0000,$BBAC,$983B,$402B),  {2.092278988800000000000E+13}
           ($0000,$C000,$7766,$A1BF,$402F),  {3.556874280960000000000E+14}
           ($0000,$9800,$6653,$B5F7,$4033),  {6.402373705728000000000E+15}
           ($0000,$4480,$C983,$D815,$4037),  {1.216451004088320000000E+17}
           ($0000,$0AD0,$9DF2,$870D,$403C),  {2.432902008176640000000E+18}
           ($0000,$AE31,$DF4D,$B141,$4040),  {5.109094217170944000000E+19}
           ($6000,$CF83,$930A,$F3BA,$4044),  {1.124000727777607680000E+21}
           ($6D00,$C526,$19AF,$AF2E,$4049),  {2.585201673888497664000E+22}
           ($D1C0,$D3DC,$9343,$8362,$404E),  {6.204484017332394393600E+23}
           ($07BC,$FB09,$0619,$CD4A,$4052)); {1.551121004333098598400E+25}
var
  FacTab: array[0..25] of extended absolute FTabH;


{---------------------------------------------------------------------------}
function sfc_fac(n: integer): extended;
  {-Return the factorial n!, n<MAXGAMX-1; INF if n<0}
begin
  if n>25 then sfc_fac := stirf(n+1)
  else begin
    if n<0 then sfc_fac := PosInf_x
    else sfc_fac := FacTab[n];
  end;
end;


{---------------------------------------------------------------------------}
function sfc_lnfac(n: longint): extended;
  {-Return ln(n!), INF if n<0}
begin
  if n>25 then sfc_lnfac := sfc_lngamma(n+1.0)
  else begin
    if n<0 then sfc_lnfac := PosInf_x
    else sfc_lnfac := ln(FacTab[n]);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_dfac(n: integer): extended;
  {-Return the double factorial n!!, n<=3209; INF for even n<0}
var
  k: integer;
  t: extended;
begin
  if n<0 then begin
    if odd(n) then begin
      if n=-1 then sfc_dfac := 1
      else begin
        {(-2k-1)!! = (-1)^k/(2k-1)!!}
        t := 1.0/sfc_dfac(-n-2);
        if n and 2 = 0 then sfc_dfac := -t
        else sfc_dfac := t;
      end;
    end
    else sfc_dfac := PosInf_x;
  end
  else if n>3209 then sfc_dfac := PosInf_x
  else if odd(n) then begin
    { (2k+1)!! = (2k+1)!/(2^k * k!) }
    if n <= 25 then begin
      k := n div 2;
      sfc_dfac := FacTab[n]/ldexp(FacTab[k],k);
    end
    else begin
      {Here n >= 27 and the argument of stirf is >= 14.5}
      n := succ(n) div 2;
      t := ldexp(stirf(n+0.5), n);
      sfc_dfac := int(t/SqrtPi+0.5);
    end;
  end
  else begin
    { n!! = 2^k * k!  with k = n div 2}
    n := n div 2;
    sfc_dfac := ldexp(sfc_fac(n),n);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_binomial(n,k: integer; dbl: boolean): extended;
  {-Return the binomial coefficient 'n choose k'; set dbl=true, if used in SpecFun}
var
  t,a: extended;
  i,m: integer;
  function lnf0(m: integer): extended;
    {-Return quick and dirty approximation of ln(m!)}
  begin
    lnf0 := (0.5+m)*ln(1.0+m)-m;
  end;
begin
  (*
  Accuracy of sfc_binomial, random values 0<=k<=n, computed with t_binom
  Range (n)      Samples   RMS           Peak
      0.. 1000   20000     1.2087E-19    5.9631E-19  for (820,415)
   1000.. 1700   20000     1.4733E-19    6.5052E-19  for (1678,331)
   1700.. 8000   20000     1.4177E-18    7.4801E-18  for (7185,4779)
   8000..16000   10000     2.0377E-16    1.2166E-15  for (14557,6857)
  16000..32000   10000     1.8608E-16    1.2806E-15  for (18717,4636)
  *)

  {If n<0 and k>=0, we use the identity 1.2.6 (17) from D.E. Knuth: The art}
  {of computer programming, Vol. 1, Fundamental algorithms, 3rd ed., 1997.}
  {[*]   binomial(n,k) = (-1)^k * binomial(k-n-1,k)}

  if (k=0) or (k=n) then sfc_binomial := 1.0
  else if (k=1) or (k=pred(n)) then sfc_binomial := n
  else if k<0 then begin
    if (n>=0) or (n<k) then sfc_binomial := 0.0
    else begin
      {Here k <= n < 0:  We use binomial(n,k) = binomial(n,n-k) which is}
      {valid for all integer k, and then use [*] since n<0 and n-k >= 0:}
      {binomial(n,k) = (-1)^(n-k)*binomial(n-k-n-1,n-k)}
      i := n-k;
      t := sfc_binomial(-succ(k),i,dbl);
      if (t<>0.0) and odd(i) then sfc_binomial := -t
      else sfc_binomial := t;
    end;
  end
  else if n<0 then begin
    {Here n < 0 and k >= 0. Apply Knuth's [*] identity}
    t := sfc_binomial(k-n-1,k,dbl);
    if (t<>0.0) and odd(k) then sfc_binomial := -t
    else sfc_binomial := t;
  end
  else if k>n then sfc_binomial := 0.0
  else begin
    {here n > 2 and  1 < k < n-1}
    if n<MAXGAMX-1 then begin
      t := sfc_fac(n);
      t := t/sfc_fac(n-k);
      t := t/sfc_fac(k);
    end
    else begin
      {if k>n/2 use symmetry to reduce k}
      if k>n-k then k := n-k;
      {all (n,k) with n <= m will NOT overflow,}
      {this saves three logarithm calculations.}
      if dbl then begin
        a := 750;
        m := 1029;
      end
      else begin
        a := 11500;
        m := 16391;
      end;
      if n<=m then t:=0
      else begin
        {Calculate rough estimate t}
        t := lnf0(k);
        t := t + lnf0(n-k);
        t := lnf0(n)-t;
      end;
      if t>a then begin
        {don't waste time (especially for double), if result is INF}
        t := PosInf_x
      end
      else if (t<11300) and (k<=4000) then begin
        {Compute n*(n-1).../k! if k is 'small' and no risk of overflow}
        t := n;
        inc(n);
        for i:=2 to k do t := t*((n-i)/i);
      end
      else begin
        {Here n > 8000 and k > 4000, or a very large result > 1.7e4777.}
        {Use lnbeta function (increased inaccuracy for 'large' values)}
        t := sfc_lnbeta(k,n-k+1);
        t := k*exp(t);
        if t<MinExtended then begin
          {for t > 0.25*MinExtended inversion will not overflow}
          if 4.0*t<MinExtended then t := PosInf_x
          else t := 1.0/t;
        end
        else t := 1.0/t;
      end;
    end;
    {round to nearest integer value}
    sfc_binomial := int(t+0.5);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_lnbinomial(n,k: longint): extended;
  {-Return ln(binomial(n,k)), n >= k >= 0}
var
  t: extended;
begin
  if (n<k) or (k<0) then begin
    sfc_lnbinomial := Nan_x;
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end
  else if (k=n) or (k=0)   then sfc_lnbinomial := 0.0
  else if (k=1) or (k=n-1) then sfc_lnbinomial := ln(n)
  else if n<MAXGAMX then begin
    t := sfc_fac(n);
    t := t/sfc_fac(n-k);
    sfc_lnbinomial := ln(t/sfc_fac(k));
  end
  else begin
    {t := sfc_lngamma(n-k+1) + sfc_lngamma(k+1);
    sfc_lnbinomial := sfc_lngamma(n+1)-t;}
    sfc_lnbinomial := -(ln(k) + sfc_lnbeta(k, n-k+1));
  end;
end;


{---------------------------------------------------------------------------}
function sfc_gdr_pos(x,d: extended): extended;
  {-Return gamma(x)/gamma(x+d), x>0, x+d > 0, accurate even for |d| << |x|}
var
  xgh,res: extended;
begin
  {$ifdef debug}
    if (x <= 0.0) and (x+d <= 0.0) then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfc_gdr_pos := NaN_x;
      exit;
    end;
  {$endif}
  if x<eps_x then begin
    {x small, Gamma(x)=1/x. Result = rGamma(x+d)/x}
    {but avoid spurious underflow if x+d > MAXGAMX}
    d := x+d;
    if d<MAXGAMX then sfc_gdr_pos := sfc_rgamma(d)/x
    else begin
      {Note very accurate but better than underflow}
      res := sfc_lngamma(d);
      res := res + ln(x);
      sfc_gdr_pos := exp(-res);
    end;
    exit;
  end;
  xgh := x + lanczos_gm05;
  if abs(d) < 10.0 then begin
    res := ln1p(d/xgh);
    res := exp((0.5-x)*res);
  end
  else begin
    res := power(xgh/(xgh+d), x-0.5);
  end;
  res := res*power(exp(1.0)/(xgh + d), d);
  sfc_gdr_pos := lanczos(x,false)/lanczos(x+d,false)*res;
end;


{---------------------------------------------------------------------------}
function sfc_gamma_delta_ratio(x,d: extended): extended;
  {-Return gamma(x)/gamma(x+d), accurate even for |d| << |x|}
var
  s,t,y: extended;
begin
  if IsNanOrInf(x) or IsNanOrInf(d) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_gamma_delta_ratio := NaN_x;
    exit;
  end;
  y := x+d;
  if x=y then sfc_gamma_delta_ratio := 1.0
  else if (x>0.0) and (y>0.0) then sfc_gamma_delta_ratio := sfc_gdr_pos(x,d)
  else begin
    if (y<=0.0) and (frac(y)=0.0) then begin
      {y=x+d is a non-positive integer}
      if (x>0.0) or (frac(x)<>0.0) then begin
        {gamma(x+d) = inf, gamma(x)<>0, return zero}
        sfc_gamma_delta_ratio := 0.0;
        exit;
      end;
      {both x,y non-positive integers, no need to calculate sin for reflection}
      t := sfc_gdr_pos(1.0-x,-d);
      if frac(0.5*x)<>frac(0.5*y) then t := -t;
      sfc_gamma_delta_ratio := 1.0/t;
      exit;
    end;
    if (x<0.0) and (y<0.0) then begin
      {Here both x and y are negative, use reflection formula for gammas}
      {gamma(x)/gamma(x+d) = sinPi(x+d)/sinP(x)/ [gamma(1-x)/gamma(1-x-d]}
      t := sfc_gdr_pos(1.0-x,-d);
      t := sinPi(x)*t;
      y := sinPi(y);
      sfc_gamma_delta_ratio := y/t;
    end
    else begin
      {Remaining cases: use lngamma for both arguments}
      s := sfc_signgamma(x)*sfc_signgamma(y);
      y := sfc_lngamma(y);
      t := sfc_lngamma(x);
      t := exp(t-y);
      sfc_gamma_delta_ratio := s*t;
    end
  end;
end;


{---------------------------------------------------------------------------}
function sfc_gamma_ratio(x,y: extended): extended;
  {-Return gamma(x)/gamma(y)}
begin
  sfc_gamma_ratio := sfc_gamma_delta_ratio(x,y-x)
end;


{---------------------------------------------------------------------------}
function sfc_pochhammer(a,x: extended): extended;
  {-Return the Pochhammer symbol (a)_x = gamma(a+x)/gamma(a).}
  { Accuracy is reduced if x or x+a are near negative integers.}
var
  t: extended;
  it,ia: boolean;
  i: integer;
begin
  if IsNanOrInf(x) or IsNanOrInf(a) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_pochhammer := NaN_x;
    exit;
  end;
  if x=0.0 then sfc_pochhammer := 1.0
  else begin
    t := a+x;
    ia := (a<=0.0) and (frac(a)=0.0);
    it := (t<=0.0) and (frac(t)=0.0);
    if it or ia then begin
      if it and ia then begin
        {both a and a+x are negative integers}
        t := sfc_gamma_delta_ratio(1-a,-x);
        if frac(0.5*x)<>0 then t := -t;
        sfc_pochhammer := t;
      end
      else if ia then sfc_pochhammer := 0.0
      else begin
        {$ifopt R+}
          if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
        {$endif}
        sfc_pochhammer := NaN_x;
      end;
    end
    else begin
      if (frac(x)=0.0) and (x>0.0) and (x<=100.0) then begin
        {x is a small positive integer, use simple multiply loop}
        t := a;
        for i:=1 to trunc(x)-1 do t := t*(a+i);
        sfc_pochhammer := t;
      end
      else begin
        t := sfc_gamma_delta_ratio(a,x);
        if t<>0 then sfc_pochhammer := 1.0/t
        else sfc_pochhammer := PosInf_x;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_poch1(a,x: extended): extended;
  {-Return (pochhammer(a,x)-1)/x, psi(a) if x=0; accurate even for small |x|}
const
  sqtbig = 0.111324015423318091e2466;  {=1/sqrt(24*MinExtended)}
  lneps  = -44.3614195558364998;       {=ln(eps_x/2)}
var
  gbern : array[1..21] of extended;
var
  t,b,bp,v,v2,lv,q,poly,res,term: extended;
  nterms,incr,j,k: integer;
begin
  {Ref: W. Fullerton [14] and [20], files dpoch1.f}
  if IsNanOrInf(x) or IsNanOrInf(a) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_poch1 := NaN_x;
    exit;
  end;

  if x=0.0 then begin
    sfc_poch1 := sfc_psi(a);
    exit;
  end;

  t := abs(x);
  q := abs(a);
  if (t > 0.1*q) or (t*ln(maxx(q,2.0)) > 0.1) then begin
    res := sfc_pochhammer(a,x);
    sfc_poch1 := (res-1.0)/x;
    exit;
  end;

  if a  < -0.5 then bp := 1.0-a-x else bp := a;
  if bp < 10.0 then incr := trunc(11.0 - bp) else incr := 0;

  b  := bp + incr;
  v  := b + 0.5*(x-1.0);
  lv := ln(v);
  q  := x*lv;

  if v < sqtbig then begin
    nterms := trunc(-0.5*lneps/lv + 1.0);
    if nterms>20 then begin
      {$ifdef debug}
        sfc_write_debug_str('*** sfc_poch1: nterms too big');
      {$endif}
      nterms := 20;
    end;
    v2  := 1.0/sqr(v);
    res := 0.5*(x+1.0);
    gbern[1] := 1.0;
    gbern[2] := -0.5*res/SIXX;
    term := v2;
    poly := gbern[2]*term;
    for k:=2 to nterms do begin
      t := 0.0;
      for j:=1 to k do t := t + extended(BoFHex[k-j])*gbern[j];
      gbern[k+1] := -res*t/k;
      term := term * ((2*k-2)-x)*((2*k-1)-x)*v2;
      poly := poly + gbern[k+1]*term;
    end;
  end
  else poly := 0.0;

  poly := (x-1.0)*poly;
  res  := exprel(q)*(lv + q*poly) + poly;

  if incr > 0 then begin
    {Get poch1(bp,x) for small bp from poch1(b,x) with backwards recursion}
    for j:=pred(incr) downto 0 do begin
      t := 1.0/(bp+j);
      res := (res-t)/(x*t+1.0);
    end;
  end;

  if bp=a then sfc_poch1 := res
  else begin
    {Use reflection formula to get poch1(a,x) with a<-0.5 from poch1(bp,x)}
    {q = cot(Pi*b)}
    sincosPi(b,t,q);
    q := q/t;
    t := sinPi(x)/x;
    t := t*q - 2.0*sqr(sinPi(0.5*x))/x;
    sfc_poch1 := t + (1.0 + x*t)*res;
  end;
end;


{---------------------------------------------------------------------------}
const
  BNNHex: array[0..8] of THexExtW = (      {bernreal(2k+2)/(2k+2), k=0..7}
            ($AAAB, $AAAA, $AAAA, $AAAA, $3FFB),   { 8.33333333333333333333333E-2}
            ($8889, $8888, $8888, $8888, $BFF8),   {-8.33333333333333333333333E-3}
            ($8208, $0820, $2082, $8208, $3FF7),   { 3.96825396825396825396825E-3}
            ($8889, $8888, $8888, $8888, $BFF7),   {-4.16666666666666666666667E-3}
            ($3E10, $E0F8, $0F83, $F83E, $3FF7),   { 7.57575757575757575757576E-3}
            ($ACCB, $CACC, $CCAC, $ACCA, $BFF9),   {-2.10927960927960927960928E-2}
            ($AAAB, $AAAA, $AAAA, $AAAA, $3FFB),   { 8.33333333333333333333333E-2}
            ($F2F3, $F2F2, $F2F2, $E2F2, $BFFD),   {-4.43259803921568627450980E-1}
            ($3FCE, $FF37, $FCDC, $C373, $4000));  { 3.05395433027011974380395   }
var
  BNN: array[0..8] of extended absolute BNNHex;


{---------------------------------------------------------------------------}
function sfc_psi(x: extended): extended;
  {-Return the psi (digamma) function of x, INF if x is a non-positive integer}
const
  PPHex: array[0..5] of THexExtW = (
           ($5D8F, $71C5, $5E47, $F7B9, $3FFE),  {+0.967672245447621170427444761712   }
           ($93CE, $3DA8, $4992, $9274, $3FFF),  {+1.14417380943934132177443086511    }
           ($6D52, $5F65, $79EF, $F169, $3FFD),  {+0.471507845373433246720039154386   }
           ($9EF2, $DF82, $1865, $A71D, $3FFB),  {+0.815984636391829369377116423465e-1}
           ($D564, $B95D, $BF51, $B443, $3FF7),  {+0.550124017486102827456755050064e-2}
           ($05E4, $B416, $530F, $C9AD, $3FF1)); {+0.961671107604462646458076964138e-4}
  PQHex: array[0..6] of THexExtW = (
           ($0000, $0000, $0000, $8000, $3FFF),  {+1.00000000000000000000000000000    }
           ($A4B3, $CD52, $FAA6, $D1E9, $3FFF),  {+1.63995297569876643247016291140    }
           ($2D1C, $1958, $35D7, $F872, $3FFE),  {+0.970492711080937113643931469693   }
           ($E38A, $0D1E, $6FB5, $84F8, $3FFD),  {+0.259707918978674869994308187160   }
           ($A4B6, $ADD5, $3BC7, $81BF, $3FFA),  {+0.316765151172736832756607905517e-1}
           ($C902, $DFDF, $6A6A, $C67A, $3FF5),  {+0.151426838914381207406214536130e-2}
           ($B01C, $4DC7, $7496, $8E9D, $3FEF)); {+0.170010400090619970545999542634e-4}
const
  x0hHex: THexExtW = ($0000,$0000,$8000,$BB16,$3FFF);  {1.4616241455078125}
  x0lHex: THexExtW = ($A9C9,$F6E1,$6BE3,$8635,$3FEE);  {7.9994605498412626595423257213e-6}
var
  PP : array[0..5] of extended absolute PPHex;
  PQ : array[0..6] of extended absolute PQHex;
  x0h: extended absolute x0hHex;
  x0l: extended absolute x0lHex;
var
  p,q,nz,s,w,y,z: extended;
  i,n: integer;
  subneg: boolean;

begin

  if IsNanOrInf(x) then begin
    if x=PosInf_x then sfc_psi := PosInf_x
    else sfc_psi := NaN_x;
    exit;
  end;

  {General algorithm layout is as in Cephes [7] double/psi.c. Coefficients }
  {of the asymptotic expansion are calculated with Pari/GP's bernreal, the }
  {Pad‚ approximation polynomials with Maple V's "pade(Psi(x+x0),x,[6,6])".}
  {Conversions to THexExtW form are done with MPArith/t_rcalc's xh command.}

  subneg := false;
  nz := 0.0;

  if x<=0.0 then begin
    q := x;
    p := floorx(q);
    if p=q then begin
      sfc_psi := PosInf_x;
      exit;
    end;
    {Remove the zeros of tan(Pi*x) by subtracting the nearest integer from x}
    nz := q - p;
    if nz<>0.5 then begin
      if nz>0.5 then begin
        p  := p + 1.0;
        nz := q - p;
      end;
      {nz := Pi/tan(Pi*nz)}
      sincosPi(nz,s,z);
      nz := Pi*z/s;
    end
    else nz := 0.0;
    x := 1.0 - x;
    subneg := nz<>0.0;
  end;

  {Check for small positive integers}
  if (x <= 12.0) and (x=floorx(x)) then begin
    y := 0.0;
    n := trunc(x);
    for i:=n-1 downto 1 do y := y+1.0/i;
    y := y - EulerGamma;
  end
  else if abs(x-x0h)<0.2 then begin
    {Pad‚ approximation for x in (x0-0.2, x0+0.2), abs. error < 3e-19}
    {where x0=1.46163214496836234126.., psi(x0)=0}
    x := x-x0h;
    x := x-x0l;
    y := x*(PolEvalX(x,PP,6)/PolEvalX(x,PQ,7));
  end
  else begin
    s := x;
    w := 0.0;
    while s<12.0 do begin
      w := w + 1.0/s;
      s := s + 1.0;
    end;
    {asymptotic expansion}
    if s>1e10 then y := 0.0
    else begin
      z := 1.0/(s*s);
      y := z*PolEvalX(z,BNN,8);
    end;
    y := ln(s) - (0.5/s) - y - w;
  end;
  if subneg then y := y-nz;
  sfc_psi := y;
end;


{---------------------------------------------------------------------------}
function sfc_psistar(x: extended): extended;
  {-Return psi(x) - ln(x), x > 0}
var
  s,w,y,z: extended;
begin
  if IsNan(x) then sfc_psistar := NaN_x
  else if x < 2.0 then sfc_psistar := sfc_psi(x) - ln(x)
  else begin
    s := x;
    w := 0.0;
    {Use recurrence relation to make s >= 12}
    while s<12.0 do begin
      {psi*(1+x) = psi(x+1)-ln(1+x) = psi(x) + 1/x - ln(1+x)}
      {          = psi*(x) + 1/x - ln(1+x) + ln(x)          }
      {          = psi*(x) + 1/x - (ln(1+x) - ln(x))        }
      {          = psi*(x) + 1/x - ln(1+1/x)                }
      {          = psi*(x) - ln1pmx(1.0/x);                 }
      w := w + ln1pmx(1.0/s);
      s := s + 1.0;
    end;
    {asymptotic expansion}
    if s >= 1e20 then y := 0.0
    else begin
      z := 1.0/(s*s);
      y := z*PolEvalX(z,BNN,9);
    end;
    sfc_psistar := (w - y) - 0.5/s;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_trigamma(x: extended): extended;
  {-Return the trigamma function psi'(x), INF if x is a negative integer}
var
  s: extended;
begin
  {trigamma = psi(1,x) = zetah(2,x) for x >0}
  if IsNanOrInf(x) then begin
    if x=PosInf_x then sfc_trigamma := 0.0
    else sfc_trigamma := NaN_x;
  end
  else if x>0.0 then begin
    if x >= 1e20 then sfc_trigamma := 1.0/x
    else sfc_trigamma := sfc_zetah(2.0,x)
  end
  else if frac(x)=0.0 then sfc_trigamma := PosInf_x
  else begin
    {Reflection formula HMF [1] 6.4.7 with n=1}
    {psi(1,1-x) + psi(1,x) = Pi^2*(1+cot(Pi*x)^2) = Pi^2/sin(Pi*x)^2}
    s := sinPi(x);
    sfc_trigamma := PiSqr/sqr(s) - sfc_zetah(2.0,1.0-x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_tetragamma(x: extended): extended;
  {-Return the tetragamma function psi''(x), NAN/RTE if x is a negative integer}
var
  c,s: extended;
const
  PP3H: THexExtW = ($E0F0,$C4EB,$DAC9,$F80C,$4004);  {2*Pi^3}
begin
  {tetragamma = psi(2,x) = -2*zetah(3,x) for x >0}
  if IsNanOrInf(x) then begin
    if x=PosInf_x then sfc_tetragamma := 0.0
    else sfc_tetragamma := NaN_x;
  end
  else if x>0.0 then begin
    if x >= 1e20 then begin
      x := 1.0/x;
      sfc_tetragamma := -x*x;
    end
    else sfc_tetragamma := -2.0*sfc_zetah(3.0,x)
  end
  else if frac(x)=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_tetragamma := NaN_x;
  end
  else begin
    {Reflection formula HMF [1] 6.4.7 with n=2, 2nd derivative of cot(Pi*x) with Maple}
    {psi(2,1-x) - psi(2,x) = 2*Pi^3*cot(Pi*x)*(1+cot(Pi*x)^2)}
    sincosPi(x,s,c);
    c := c/s;
    c := c*(1.0 + c*c);
    s := -2.0*sfc_zetah(3.0,1.0-x);
    sfc_tetragamma := s - extended(PP3H)*c;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_pentagamma(x: extended): extended;
  {-Return the pentagamma function psi'''(x), INF if x is a negative integer}
var
  c,s: extended;
const
  PP4H: THexExtW = ($2C68,$4841,$7461,$C2D1,$4006);  {2*Pi^4}
begin
  {pentagamma = psi(3,x) = 6*zetah(4,x) for x >0}
  if IsNanOrInf(x) then begin
    if x=PosInf_x then sfc_pentagamma := 0.0
    else sfc_pentagamma := NaN_x;
  end
  else if x>0.0 then begin
    if x >= 1e20 then begin
      x := 1.0/x;
      sfc_pentagamma := 2.0*x*x*x;
    end
    else sfc_pentagamma := 6.0*sfc_zetah(4.0,x);
  end
  else if frac(x)=0.0 then sfc_pentagamma := PosInf_x
  else begin
    {Reflection formula HMF [1] 6.4.7 with n=3, 3rd derivative of cot(Pi*x) with Maple}
    {psi(3,1-x) + psi(3,x) = 2*Pi^4*(1+4*cot(Pi*x)^2 + 3*cot(Pi*x)^4)}
    sincosPi(x,s,c);
    c := sqr(c/s);
    c := 1.0 + c*(4.0 + 3.0*c);
    s := 6.0*sfc_zetah(4.0,1.0-x);
    sfc_pentagamma := extended(PP4H)*c - s;
  end;
end;


const
  NXPGPMAX = 100;  {max. n for polygamma(n,x) with x < 0 via cot polynomial}

{---------------------------------------------------------------------------}
function polygam_cotpoly(n: integer; x: extended): extended;
  {-Polygamma for negative x and 1 <= n <= NXPGPMAX}
type
  TCArray = array[-1..NXPGPMAX+2] of extended;
  PCArray = ^TCArray;
  { for i from 1 to 12 do sort(simplify(expand(diff(cot(x),x$i)))); od;
     -C^2 - 1             with C = cot(x)
     2*C^3 + 2*C
     -6*C^4 -8*C^2 - 2
     24*C^5 + 40*C^3 + 16*C
     -120*C^6 - 240*C^4 - 136*C^2 - 16
     720*C^7 + 1680*C^5 + 1232*C^3 + 272*C
     -5040*C^8 - 13440*C^6 - 12096*C^4 - 3968*C^2 - 272
     40320*C^9 + 120960*C^7 + 129024*C^5 + 56320*C^3 + 7936*C
     -362880*C^10 - 1209600*C^8 - 1491840*C^6 - 814080*C^4 - 176896*C^2 - 7936
     3628800*C^11 + 13305600*C^9 + 18627840*C^7 + 12207360*C^5 + 3610112*C^3 + 353792*C
     -39916800*C^12 - 159667200*C^10 - 250145280*C^8 - 191431680*C^6 - 71867136*C^4 - 11184128*C^2 - 353792
     479001600*C^13 + 2075673600*C^11 + 3597834240*C^9 + 3149752320*C^7 + 1436058624*C^5 + 309836800*C^3 + 22368256*C
  }
const
  ca: array[0..53] of single = ({Above coefficients reflected and positive,}
        { 1} 1,1,               {all are exactly representable as single!}
        { 2} 2,2,
        { 3} 2,8,6,
        { 4} 16,40,24,
        { 5} 16,136,240,120,
        { 6} 272,1232,1680,720,
        { 7} 272,3968,12096,13440,5040,
        { 8} 7936,56320,129024,120960,40320,
        { 9} 7936,176896,814080,1491840,1209600,362880,
        {10} 353792,3610112,12207360,18627840,13305600,3628800,
        {11} 353792,11184128,71867136,191431680,250145280,159667200,39916800,
        {12} 22368256,309836800,1436058624,3149752320.0,3597834240.0,2075673600,479001600);
  ko: array[1..12] of integer = (0,2,4,7,10,14,18,23,28,34,40,47); {index of first coeff for n}
var
  c,s,t: extended;
  i,k: integer;
  pa: PCArray;
  sa: word;
begin
  sincosPi(x,s,c);
  if (n<1) or (n>NXPGPMAX) or (x>=0.0) or (s=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    polygam_cotpoly := NaN_x;
    exit;
  end;

  {Reflection formula: HMF[1] (6.4.7) or  NIST[30], (5.15.6)}
  {psi(n,1-x) + (-1)^(n+1)*psi(n,x) = (-1)^n*Pi*d^n/dx^n(cot(Pi*x))}

  {psi(n,x) = - (-1)^(n+1)*psi(n,1-x) + (-1)^(n+n+1)*Pi*d^n/dx^n(cot(Pi*x))}
  {psi(n,x) = (-1)^n*psi(n,1-x) - Pi*d^n/dx^n(cot(Pi*x))}

  if (c=0.0) and (n and 1 = 0) then begin
    {early exit for cot=0: for even n cot recflection term -t*c below will be 0}
    polygam_cotpoly := sfc_polygamma(n,1.0-x);
    exit;
  end;

  {c = cot(Pi*x), s = cot(Pi*x)^2}
  c := c/s;
  s := c*c;

  {Compute polynomial in cot(Pi*x)^2 in t}
  if n<=12 then begin
    {use hard-coded coefficients}
    k := ko[n];
    t := 0.0;
    for i := (succ(n) div 2) downto 0 do t := t*s + ca[k+i];
  end
  else begin
    {Recursive computation using cot'(x) = -(1+cot^2(x)). As above the}
    {coefficient*(-1)^n are evaluated, the working array contains the }
    {interlaced coefficients for n-1 and n, indices -1...n+2 are used.}
    sa := (n+4)*sizeof(extended);
    {allocate and clear to zero}
    pa := calloc(sa);
    if pa=nil then begin
      {$ifopt R+}
        {Heap overflow error}
        RunError(203);
      {$endif}
      polygam_cotpoly := NaN_x;
      exit;
    end;
    {pa^[k] contains the positive coefficient for cot^k(x)}
    pa^[1] := 1.0;
    for i:=2 to n+1 do begin
      k := i and 1;
      while k <= i do begin
        pa^[k] := pa^[k-1]*(k-1) + pa^[k+1]*(k+1);
        inc(k,2);
      end;
    end;
    {Coefficients done, evaluate polynomial in cot(Pi*x)^2}
    t := pa^[n+1];
    i := n-1;
    repeat
      t := t*s + pa^[i];
      dec(i,2);
    until i<0;
    mfree(pa,sa);
  end;
  s := sfc_polygamma(n,1.0-x);
  if odd(n) then s := -s
  else begin
    {t is already correct for odd n, for even n multiply with -cot(Pi*x)}
    {the minus compensates the fact that the coefficients are positive. }
    t := -t*c;
  end;
  t := t*intpower(Pi,n+1);
  polygam_cotpoly := s + t;
end;


{---------------------------------------------------------------------------}
function polygam_negx(n: integer; x: extended): extended;
  {-Polygamma for finite negative x and 1 <= n <= MAXGAMX}
var
  n1,k: integer;
  z,s,t,a,b,d: extended;
const
  tmin = 13;   {should > 12}
  tmax = 70.0; {must by <= NXPGPMAX}
begin
  z  := frac(x);
  if (z=0.0) or (x>0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    polygam_negx := NaN_x;
    exit;
  end;
  if n <= 12 then begin
    {Use hard-coded polynomials}
    polygam_negx := polygam_cotpoly(n,x);
    exit;
  end;

  n1 := n+1;
  if (z=-0.5) and odd(n1) then begin
    {cot term is zero for even n}
    s := 0.0;
  end
  else begin
    if z=-0.5 then begin
      t := tmin;
      z := 0.5;
    end
    else begin
      if z<-0.5 then z := z+1.0;
      {-0.5 < z < 0.5}
      t := -5.0*log2(0.5 - abs(z));
      if t<tmin then t := tmin
      else if t>tmax then t := tmax;
    end;
    if n < t then begin
      {presumed slow convergence, use polynomial}
      polygam_negx := polygam_cotpoly(n,x);
      exit;
    end;
    s := power(z,-n1);
    k := 1;
    repeat
      a := power(z+k,-n1);
      b := power(z-k,-n1);
      d := a + b;
      s := s+d;
      inc(k);
    until abs(d)<=eps_x*abs(s);
    if odd(n) then s := -s;
  end;
  t := sfc_zetah(n1, 1.0-x);
  a := s+t;
  b := sfc_fac(n);
  if abs(a) <= MaxExtended/b then polygam_negx := -b*a
  else polygam_negx := copysign(PosInf_x, -a);
end;


{---------------------------------------------------------------------------}
function polygam_sum(n: integer; x: extended): extended;
  {-Direct summation polygamma(n,x) = n!/x^(n+1)*(1 + sum(1/(1+k/x)^(n+1), k=1...))}
const
  KMAX = 200;
var
  s,t,e: extended;
  k,m: integer;
begin
  m := -(n+1);
  e := 0.125*eps_x;
  s := 0.0;
  k := 1;
  repeat
    t := power(1.0+k/x, m);
    s := s + t;
    inc(k);
  until (t<=e) or (k>KMAX);
  if k>KMAX then begin
    {no convergence}
    polygam_sum := NaN_x;
  end
  else begin
    t := sfc_lngamma(n+1) - (n+1)*ln(x);
    if t+ln1p(s) >= ln_MaxExt then polygam_sum := PosInf_x
    else polygam_sum := (1.0+s)*exp(t);
  end;
end;


{---------------------------------------------------------------------------}
function polygam_ase(n: integer; x: extended): extended;
  {-Polygamma asymptotic expansion}
var
  z,f,t,b: extended;
  k,j: integer;
const
  K2MAX =  120; {doubled maximum index in asymptotic sum}
  NPMAX = 2000; {Max. n for calculation of (n-1)!/x^n as product}
  TFMAX = 2.6E4913;  {~ 2*MaxExtended*eps_s, used for 'no-convergence' check}
begin
  f := sfc_lngamma(n);
  z := n*ln(x);
  t := f-z;
  if t > ln_MaxExt then z := PosInf_x
  else if t < ln_MinExt then z := 0.0
  else begin
    {calculate f = (n-1)!/x^n}
    if (abs(z)<11356.0) and (n<MAXGAMX-1) then begin
      f := sfc_fac(n-1)*power(x,-n);
    end
    else begin
      if n <= NPMAX then begin
        {avoid exp/ln inaccuracies, example: sfc_polygamma(1800,980)}
        f := 1.0;
        for k:=1 to n-1 do f := f*k/x;
        f := f/x;
      end
      else f := exp(t);
    end;
    {First two terms}
    z := f + 0.5*f*(n/x);
    {Add sum terms with Bernoulli numbers}
    k := 2;
    repeat
      j := k-1;
      t := (n+j)/k*(n+j-1)/j/x/x;
      b := sfc_bernoulli(k);
      if (k>K2MAX) or (abs(f) > abs(TFMAX/t/b)) then begin
        {Convergence failed if k too large or a summand ~> MaxExtended}
        {Try direct summation, complete failure if result is NaN}
        z := polygam_sum(n,x);
        {$ifopt R+}
          if (IsNaN(z)) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
        {$endif}
        polygam_ase := z;
        exit;
      end;
      f := f*t;
      t := f*b;
      z := z + t;
      k := k+2;
    until abs(t)<abs(z)*eps_x;
  end;
  polygam_ase := z
end;


{---------------------------------------------------------------------------}
function sfc_polygamma(n: integer; x: extended): extended;
  {-Return the polygamma function: n'th derivative of psi; n>=0, x>0 for n>MAXGAMX}
  { Note: The accuracy may be reduced for n >= MAXGAMX due to ln/exp operations.}
var
  z,f,t: extended;
  k: integer;
const
  NPMAX = 2500; {Max. n > 1752 for calculation of n! as product}
begin
  if n<4 then begin
    case n of
        0: sfc_polygamma := sfc_psi(x);
        1: sfc_polygamma := sfc_trigamma(x);
        2: sfc_polygamma := sfc_tetragamma(x);
        3: sfc_polygamma := sfc_pentagamma(x);
      else begin
             {RTE or NAN if n < 0}
             {$ifopt R+}
               if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
             {$endif}
             sfc_polygamma := NaN_x;
           end;
    end;
  end
  else if IsNanOrInf(x) then begin
    if x=PosInf_x then begin
      if n>0 then sfc_polygamma := 0.0
      else sfc_polygamma := PosInf_x
    end
    else sfc_polygamma := NaN_x;
  end
  else if x > 0.0 then begin
    {polygamma(n,x) = (-1)^(n+1)*n!*sum(1/(x+k)^(n+1), k>=0), HMF[1] 6.4.10}
    {or with the Hurwitz zeta:  polygamma(n,x) = (-1)^(n+1)*n!*zetah(n+1,x)}
    z := sfc_zetah(n+1,x);
    if z=0.0 then begin
      {If zetah=0, use HMF[1] 6.4.10 or the asymptotic formula 6.4.11}
      if x>0.25*n then z := polygam_ase(n,x)
      else z := polygam_sum(n,x)
    end
    else if n<MAXGAMX-1 then begin
      {n! is OK, but n!*zetah may overflow for small x}
      f := sfc_fac(n);
      t := ln(f)+ln(z);
      if t<ln_MaxExt then z := z*f
      else z := PosInf_x;
    end
    else begin
      f := sfc_lngamma(n+1);
      t := f + ln(z);
      if t > ln_MaxExt then z := PosInf_x
      else if t < ln_MinExt then z := 0.0
      else begin
        if n>NPMAX then z := exp(t)
        else begin
          {example: sfc_polygamma(2001,32)}
          z := sfc_fac(1752)*z;
          for k:=1753 to n do z := z*k;
        end;
      end;
    end;
    if odd(n) then sfc_polygamma := z
    else sfc_polygamma := -z;
  end
  else if (x<0.0) and (n>MAXGAMX) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_polygamma := NaN_x;
  end
  else begin
    {x<0 and n<=nmax}
    sfc_polygamma := polygam_negx(n,x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_batemang(x: extended): extended;
  {-Return the Bateman function G(x) = psi((x+1)/2) - psi(x/2), x<>0,-1,-2,...}
const
  z0 = 1e5;
  z1 = 1e19;
var
  s,t,r,z: extended;
  k: integer;
begin
  if IsNanOrInf(x) or ((x<=0.0) and (frac(x)=0.0)) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_batemang := NaN_x;
    exit;
  end;
  if abs(x)<=sqrt_epsh then begin
    {G(x)  = 2/x - 2*ln(2) + Pi^2/6*x + O(x^2)}
    sfc_batemang := 2.0*(1.0/x - ln2);
    exit;
  end;

  {If x<0 use reflection formula}
  if x>0.0 then z := x
  else z := 1.0 - x;

  if z>z0 then begin
    {large z}
    if z>=z1 then s := 1.0/z
    else begin
      {G(z) = 1/z + 0.5/z^2 - 0.25/z^4 + O(1/z^6)}
      t := 0.5/sqr(z);
      s := 1.0/z + t*(1-t);
    end;
  end
  else begin
    r := 0;
    {Double step recursion Erdelyi et al.[50], 1.8(7)}
    {Use G(z) = G(z+2) + 2/z/(z+1) to make z >= 20}
    while z<20.0 do begin
      r := r + 2.0/(z*(1.0+z));
      z := z + 2.0;
    end;
    s := 1.0;
    t := 1.0;
    k := 0;
    {Erdelyi et al.[50], 1.8(6): G(z) = 2/z * 2F1(1,z,1+z,-1). Using a linear}
    {transformation gives  G(z) = 2F1(1,1,1+z,1/2)/z.  This is a very simple }
    {sum with positive terms, cf. unit sfHyperg. Max. k is about 25 for z=20.}
    repeat
      inc(k);
      t := 0.5*t/(z+k)*k;
      s := s + t;
    until t < eps_x*s;
    s := s/z + r;
  end;
  if x<0.0 then begin
    {Reflection formula Erdelyi et al.[50], 1.8(8)}
    t := TwoPi/sinPi(x);
    s := t - s;
  end;
  sfc_batemang := s;
end;


{---------------------------------------------------------------------------}
function sfc_lnbeta(x,y: extended): extended;
  {-Return the logarithm of |beta(x,y)|=|gamma(x)*gamma(y)/gamma(x+y)|}
var
  c,s: extended;
begin
  if IsNanOrInf(x) or IsNanOrInf(x) then begin
    sfc_lnbeta := NaN_x;
    exit;
  end;
  {based on SLATEC [14] routine dlbeta.f}
  {force x <= y}
  if x>y then begin
    c := x;
    x := y;
    y := c;
  end;
  s := x+y;
  if x >= XLGAE then begin
    {x and y >= XLGAE}
    c := sfc_lngcorr(x) + sfc_lngcorr(y) - sfc_lngcorr(s);
    sfc_lnbeta := -0.5*ln(y) + lnsqrt2pi + c + (x-0.5)*ln(x/s) + y*ln1p(-x/s);
  end
  else if (y >= XLGAE) and (s >= XLGAE) then begin
    {Fixed in 1.02.01: check s >= XLGAE, otherwise lnbeta(-10.1,10) fails}
    {x < XLGAE, y >= XLGAE}
    c := sfc_lngcorr(y) - sfc_lngcorr(s);
    sfc_lnbeta := sfc_lngamma(x) + c + x - x*ln(s) + (y-0.5)*ln1p(-x/s)
  end
  else begin
    {x,y < XLGAE; or s<=0}
    sfc_lnbeta := sfc_lngamma(x) + sfc_lngamma(y) - sfc_lngamma(s);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_beta(x,y: extended): extended;
  {-Return the function beta(x,y)=gamma(x)*gamma(y)/gamma(x+y)}
var
  r,t: extended;
begin
  {If one arg is 1, return the inverse of the other because}
  {gamma(1)*gamma(y)/gamma(y+1) = 1/y, analogue for y=1}
  if IsNanOrInf(x) or IsNanOrInf(x) then sfc_beta := NaN_x
  else if x=1.0 then sfc_beta := 1.0/y
  else if y=1.0 then sfc_beta := 1.0/x
  else begin
    {force x >= y}
    if x<y then begin
      t := x;
      x := y;
      y := t;
    end;
    t := x+y;
    if (t<=0.0) and (frac(t)=0.0) then begin
      if (frac(x)<>0.0) and (frac(y)<>0.0) then begin
        {gamma(x+y) = Inf, gamma(x)*gamma(y)<>Inf}
        sfc_beta := 0.0;
      end
      else begin
        t := sfc_gamma(x);
        r := sfc_pochhammer(y,x);
        sfc_beta := t/r;
      end;
      exit;
    end;
    if t=x then begin
      {gamma(x)/gamma(x+y)=1}
      sfc_beta := sfc_gamma(y);
    end
    else begin
      r := sfc_lnbeta(x,y);
      t := exp(r);
      if (x<0.0) or (y<0.0) then begin
        {get signs only if necessary}
        r := sfc_signgamma(x)*sfc_signgamma(y)/sfc_signgamma(x+y);
        sfc_beta := r*t;
      end
      else sfc_beta := t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_taylor(x: extended; n: integer): extended;
  {-Return the Taylor coefficient x^n/n!, n>=0}
var
  z,f: extended;
  k: integer;
const
  zmin = -10991.3977; { > ln(MinExtended*100!)}
begin
  if IsNanOrInf(x) or (n<0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_taylor := NaN_x;
    exit;
  end;
  if n<=0 then sfc_taylor := 1.0
  else if (n=1) or (x=0.0) then sfc_taylor := x
  else begin
    z := n*ln(abs(x));
    if z<ln_MinExt then sfc_taylor := 0.0
    else if (n<=100) and (z<ln_MaxExt) and (z>zmin) then begin
      {do not waste time for sfc_lnfac(n)}
      f := 1.0;
      for k:=1 to n do f := f*(x/k);
      sfc_taylor := f;
    end
    else begin
      f := sfc_lnfac(n);
      f := z-f;
      if f<ln_MinExt then sfc_taylor := 0.0
      else if f>ln_MaxExt then sfc_taylor := copysign(PosInf_x,x)
      else begin
        {no over/underflow}
        if n>2000 then begin
          f := exp(f);
          if (x<0.0) and odd(n) then f := -f;
        end
        else begin
          f := 1.0;
          for k:=1 to n do f := f*(x/k);
        end;
        sfc_taylor := f;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function lnbg1p(x: extended): extended;
  {-Return ln(BarnesG(1+x) via Taylor series, -0.5 <= x <= 0.5}
var
  s,t,y: extended;
  k: integer;
const
  c1 = 0.4189385332046727418;   {(ln(2*Pi)-1)/2}
  c2 = 0.7886078324507664303;   {(1+Eulergamma)/2}
  KMAX = 100;
begin
  {Ref: Formula (8) at http://mathworld.wolfram.com/BarnesG-Function.html}
  s := x*(c1 - c2*x);
  y := x*x*x;
  for k:=2 to KMAX do begin
    t := sfc_zetaint(k)/(k+1)*y;
    s := s + t;
    if abs(t) <= eps_x*abs(s) then break;
    y := -y*x;
  end;
  lnbg1p := s;
end;

(*
{---------------------------------------------------------------------------}
function lnbg1p(x: extended): extended;
  {-Return ln(BarnesG(1+x) via Taylor series, -0.5 <= x <= 0.5}
var
  s,t,y: extended;
  k: integer;
const
  c1 = 0.4189385332046727418;   {(ln(2*Pi)-1)/2}
  c2 = -0.2113921675492335697;  {(Eulergamma-1)/2}
  KMAX = 100;
begin
  {Ferreira/Lopez: An Asymptotic Expansion of the Double Gamma Function,}
  {formula (2): About 10% less iterations, sligtly less accurate}
  s := x*(c1 + c2*x) + x*sfc_lngamma1p(x);
  y := x*x*x;
  for k:=2 to KMAX do begin
    t := sfc_zetaint(k)/(k*(k+1))*y;
    s := s - t;
    if abs(t) <= eps_x*abs(s) then break;
    y := -y*x;
  end;
  lnbg1p := s;
end;
*)


{---------------------------------------------------------------------------}
function sfc_lnbg(x: extended): extended;
  {-Return ln(BarnesG(x)), real part for x < 0}
const
  c1 = 0.9189385332046727418;  {-Zeta'(0) = ln(2*Pi)/2}
  c2 = 0.1654211437004509292;  {-Zeta'(1)}
  XAE  = 10;    {Min x for asymptotic expansion}
  KMAX = 64;    {Max terms in AE}
var
  g,s,y,z,t: extended;
  k: integer;
begin
  if IsNanOrInf(x) then begin
    sfc_lnbg := NaN_x;
    exit;
  end;
  if (frac(x)=0.0) and (x<=28.0) then begin
    if x<=0.0 then sfc_lnbg := NegInf_x
    else begin
      {Logarithmic version of functional equation}
      s := 0.0;
      for k:=2 to round(x)-2 do s := s + sfc_lnfac(k);
      sfc_lnbg := s;
    end;
    exit;
  end;
  if x < -0.5 then begin
    {Use reflection formula, Adamchik[78], (35)}
    z := x-2.0;
    g := sfc_lnbg(-z);
    s := (x-1.0) * ln(abs(sinPi(z)/Pi));
    {We can use frac(z) because Cl_2 has period 2*Pi}
    t := sfc_cl2(frac(z)*TwoPi)/TwoPi;
    sfc_lnbg := g - s - t;
  end
  else if x < XAE then begin
    if x < 0.5 then begin
      {lnG(x) = lnG(x+1) - lnGamma(x)}
      sfc_lnbg := lnbg1p(x) - sfc_lngamma(x);
    end
    else if x <= 1.5 then begin
      {x is already in [0.5 , 1.5]}
      sfc_lnbg := lnbg1p(x-1.0);
    end
    else begin
      {Use functional equation to shift x downto [0.5 , 1.5]}
      g := 0.0;
      while x>1.5 do begin
        {lnG(x) = lnG(x-1) + lnGamma(x-1)}
        x := x - 1.0;
        g := g + sfc_lnGamma(x);
      end;
      sfc_lnbg := lnbg1p(x-1.0) + g;
    end;
  end
  else begin
    {Here x > XAE}
    x := x - 1;
    {Asymptotic expansion, Adamchik[78], (28)}
    z := sqr(x);
    g := (c1*x - c2) + (0.5*z - 1/12)*ln(x) - 0.75*z;
    y := z; {=x^2k}
    s := 0;
    for k:=1 to KMAX do begin
      t := sfc_bernoulli(2*k+2);
      t := t / (4*k*(k+1)*y);
      if abs(t) <= eps_x*abs(s) then break;
      y := y * z;
      s := s + t;
    end;
    sfc_lnbg := g + s;
  end;
end;



end.
