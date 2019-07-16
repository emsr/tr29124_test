unit sdGamma;

{Double precision special functions: Gamma function and related}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Double precision special functions: Gamma function and related

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in damath_info.txt/references

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
                 [30] [NIST]: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST Handbook
                      of Mathematical Functions, Cambridge, 2010. Online resource: NIST Digital
                      Library of Mathematical Functions, http://dlmf.nist.gov/
                 [32] D.E. Knuth: The Art of computer programming; Volume 1, Fundamental Algorithms,
                      3rd ed., 1997; Volume 2, Seminumerical Algorithms, 3rd ed., 1998.
                 [50] A. Erdelyi et al., Higher Transcendental Functions Vol. I-III, California
                      Institute of Technology - Bateman Manuscript Project, 1953-1955,
                      Available via http://en.wikipedia.org/wiki/Bateman_Manuscript_Project
                 [78] V.S. Adamchik, Contributions to the Theory of the Barnes function,
                      2003, https://arxiv.org/pdf/math/0308086v1.pdf

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  07.02.13  W.Ehrhardt  Initial BP7 version from AMath.sfgamma
 1.00.01  07.02.13  we          stirf, sfd_gaminv_small, sfd_gamma_medium
 1.00.02  08.02.13  we          sfd_lngcorr, sfd_gstar, lngamma_small, lanczos
 1.00.02  08.02.13  we          sfd_fac, sfd_dfac, sfd_binomial, sfd_psi, sfd_poch1
 1.00.03  09.02.13  we          sfd_trigamma/tetragamma/pentagamma/polygamma
 1.00.04  10.02.13  we          sfd_lngamma uses FacTab
 1.00.05  11.02.13  we          igam_aux, sfd_igprefix, sfd_incgamma_inv
 1.00.06  14.02.13  we          improved sfd_git
 1.00.07  01.03.13  we          Chebyshev degrees reduced in sfd_gstar

 1.01.00  30.03.13  we          improved sfd_binomial

 1.02.00  18.04.13  we          improved sfd_beta
 1.02.01  18.04.13  we          ibeta_series with parameter normalised
 1.02.02  18.04.13  we          sfd_nnbeta
 1.02.03  21.04.13  we          lanczos_gm05 = lanczos_g - 0.5 in implementation
 1.02.04  27.04.13  we          polygam_negx for n=11,12

 1.03.00  11.05.13  we          improved sfd_polygamma for large order: e.g. (3001, 871)
 1.03.01  26.05.13  we          sfd_nnbeta for a<=0 or b<=0

 1.05.00  28.07.13  we          sfd_git for x<0

 1.10.00  18.04.14  we          Improved polygam_negx (larger n for x<0}

 1.12.00  16.06.14  we          Incomplete/inverse functions moved to sdGamma2
 1.12.01  16.06.14  we          const lanczos_gm05: double

 1.13.00  06.08.14  we          First check IsNaND(x) in sfd_lngammas
 1.13.01  06.08.14  we          improved sfd_gdr_pos for small x
 1.13.02  11.08.14  we          special case x=x+y in sfd_beta
 1.13.03  20.08.14  we          sfd_gamma_delta_ratio returns 1 if x=x+d

 1.16.00  15.12.14  we          small changes in sfd_batemang

 1.21.00  11.10.15  we          polygamma for x<0 uses Euler series for Pi*cot(Pi*x)

 1.33.00  08.06.18  we          sfd_lnbinomial

 1.35.00  23.08.18  we          Check NaNs in sfd_gamma
 1.35.01  01.09.18  we          sfd_lnbg

 1.35.02  12.10.18  we          sfd_psistar

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
  XLGAE   = 8.0;                     {threshold for asymptotic expansion of sfd_lngamma}
  MAXGAMD = 171.62437695630272;      {max. argument for gamma}
  MAXLGM  = 2.559983327851638e305;   {max. argument for lngamma}
  MAXDFAC = 300;                     {max. argument for dfac}

function stirf(x: double): double;
  {-Stirling's formula for the gamma function for x > 13.0}
  { gamma(x) = sqrt(2*Pi)*x^(x-0.5)*exp(-x)*(1 + 1/x * P(1/x))}

function sfd_lngcorr(x: double): double;
  {-Return lngamma correction term lngamma(x) - ((x-0.5)*ln(x) - x + ln(sqrt(2*Pi)), x>=8}

function sfd_rgamma(x: double): double;
  {-Return the reciprocal gamma function rgamma = 1/gamma(x)}

function sfd_binomial(n,k: integer): double;
  {-Return the binomial coefficient 'n choose k'; set dbl=true, if used in SpecFun}

function sfd_lnbinomial(n,k: longint): double;
  {-Return ln(binomial(n,k)), n >= k >= 0}

function sfd_fac(n: integer): double;
  {-Return the factorial n!, n<MAXGAMD-1; INF if n<0}

function sfd_lnfac(n: longint): double;
  {-Return ln(n!), INF if n<0}

function sfd_dfac(n: integer): double;
  {-Return the double factorial n!!, n<=300; INF for even n<0}

function sfd_gamma(x: double): double;
  {-Return gamma(x), |x| <= MAXGAMD; NAN/RTE if x is a non-positive integer}

function sfd_gamma_medium(x: double): double;
  {-Return gamma(x), |x| <= 13, x negative integer produces div by 0}

function sfd_gaminv_small(x: double): double;
  {- 1/gamma(x) for small argument values, |x| < 0.03125}

function sfd_gamma1pm1(x: double): double;
  {-Return gamma(1+x)-1 with increased accuracy for x near 0}

function sfd_gamma_delta_ratio(x,d: double): double;
  {-Return gamma(x)/gamma(x+d), accurate even for |d| << |x|}

function sfd_gamma_ratio(x,y: double): double;
  {-Return gamma(x)/gamma(y)}

function sfd_pochhammer(a,x: double): double;
  {-Return the Pochhammer symbol (a)_x = gamma(a+x)/gamma(a).}

function sfd_poch1(a,x: double): double;
  {-Return (pochhammer(a,x)-1)/x, psi(a) if x=0; accurate even for small |x|}

function sfd_lngamma(x: double): double;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, NAN/RTE if x is a non-positive integer}

function sfd_lngammas(x: double; var s: integer): double;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, s=-1,1 is the sign of gamma}

function sfd_lngamma1p(x: double): double;
  {-Return ln(|gamma(1+x)|) with increased accuracy for x near 0}

function sfd_signgamma(x: double): double;
  {-Return sign(gamma(x)), useless for 0 or negative integer}

function sfd_gstar(x: double): double;
  {-Return Temme's gstar(x) = gamma(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x)), x>0.}
  { For large x the asymptotic expansion is gstar(x) = 1 + 1/12x + O(1/x^2)}

function sfd_psi(x: double): double;
  {-Return the psi (digamma) function of x, INF if x is a non-positive integer}

function sfd_psistar(x: double): double;
  {-Return psi(x) - ln(x), x > 0}

function sfd_trigamma(x: double): double;
  {-Return the trigamma function psi'(x), INF if x is a negative integer}

function sfd_tetragamma(x: double): double;
  {-Return the tetragamma function psi''(x), NAN/RTE if x is a negative integer}

function sfd_pentagamma(x: double): double;
  {-Return the pentagamma function psi'''(x), INF if x is a negative integer}

function sfd_polygamma(n: integer; x: double): double;
  {-Return the polygamma function: n'th derivative of psi; n>=0, x>0 for n>MAXGAMD.}
  { Note: The accuracy may be reduced for n>=MAXGAMD due to ln/exp operations.}

function sfd_lnbg(x: double): double;
  {-Return ln(BarnesG(x)), real part for x < 0}

function sfd_batemang(x: double): double;
  {-Return the Bateman function G(x) = psi((x+1)/2) - psi(x/2), x<>0,-1,-2,...}

function sfd_beta(x,y: double): double;
  {-Return the function beta(x,y)=gamma(x)*gamma(y)/gamma(x+y)}

function sfd_lnbeta(x,y: double): double;
  {-Return the logarithm of |beta(x,y)|=|gamma(x)*gamma(y)/gamma(x+y)|}

function sfd_taylor(x: double; n: integer): double;
  {-Return the Taylor coefficient x^n/n!, n>=0}

{#Z+}
const
  lanczos_gm05: double = 12.644565;  {=($4AE0,$6C76,$4A04,$4029)}

function lanczos(x: double; egscale: boolean): double;
  {-Return the Lanczos sum for x, exp(g) scaled if egscale=true}
{#Z-}


implementation

uses
  memh,
  DAMath,
  sdBasic,  {Basic common code}
  sdZeta;   {Zeta functions and polylogarithms}


{---------------------------------------------------------------------------}
function lanczos(x: double; egscale: boolean): double;
  {-Return the Lanczos sum for x, exp(g) scaled if egscale=true}
const
  lsh:  array[0..12] of THexDblW = (
          ($324E,$F03B,$03AF,$42C4),  {+4.4012138428004608955E+13}
          ($009A,$1968,$E9C4,$42C2),  {+4.1590453358593200516E+13}
          ($2DFF,$B34B,$622C,$42B0),  {+1.8013842787117996778E+13}
          ($CD8E,$94B4,$33FA,$4291),  {+4.7287362634753888969E+12}
          ($8CF3,$A705,$62EA,$4268),  {+8.3791008362840464705E+11}
          ($6E10,$8489,$9547,$4238),  {+1.0558370727342993449E+11}
          ($F5C2,$7D13,$11F8,$4202),  {+9.7013636184949994935E+9}
          ($2B97,$AEC6,$8499,$41C3),  {+6.5491439754820526409E+8}
          ($FAA5,$2F12,$BEAF,$417E),  {+3.2238322942133565306E+7}
          ($F61D,$3830,$3842,$4131),  {+1.1285142194970914380E+6}
          ($EAC1,$CD5D,$0A72,$40DA),  {+2.6665793784598589447E+4}
          ($03D1,$FDCF,$DE14,$4077),  {+3.8188012486329268705E+2}
          ($2706,$1FF6,$0D93,$4004)); {+2.5066282746310005025}
  lseh: array[0..12] of THexDblW = (
          ($017D,$2623,$869C,$4194),  {+8.6091529534185372176E+7}
          ($DDB4,$24B6,$657C,$4193),  {+8.1354505178580112428E+7}
          ($F0CF,$931A,$CD58,$4180),  {+3.5236626388154619108E+7}
          ($E57F,$DF9D,$A482,$4161),  {+9.2498149880244712949E+6}
          ($CF12,$3778,$0270,$4139),  {+1.6390242166871469602E+6}
          ($571B,$86AF,$3616,$4109),  {+2.0653081576412250327E+5}
          ($0C74,$EC82,$882C,$40D2),  {+1.8976701935302889156E+4}
          ($F341,$9051,$0446,$4094),  {+1.2810689099125594799E+3}
          ($90E7,$AAB1,$87CC,$404F),  {+6.3060933434202345361E+1}
          ($7133,$8223,$A8E6,$4001),  {+2.2074709097925276381}
          ($CB06,$DE25,$B4CA,$3FAA),  {+5.2160586946135054274E-2}
          ($440D,$A0D0,$7A35,$3F48),  {+7.4699038089154483163E-4}
          ($41AC,$036E,$90C0,$3ED4)); {+4.9031805734598718626E-6}
type
  TX12 = array[0..12] of double;

  function lratev(const num: TX12; z: double): double;
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
    s1,s2: double;
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
function stirf(x: double): double;
  {-Stirling's formula for the gamma function for x > 13.0}
  { gamma(x) = sqrt(2*Pi)*x^(x-0.5)*exp(-x)*(1 + 1/x * P(1/x))}
const
  {13 <= x <= 1024, relative peak error = 9.44E-21, relative error spread = 8.8e-4}
  SPHex: array[0..8] of THexDblW = (
           ($5554,$5555,$5555,$3FB5),  {+8.3333333333333314830E-2}
           ($0DD8,$1C72,$71C7,$3F6C),  {+3.4722222222300751226E-3}
           ($5907,$8F0B,$F726,$BF65),  {-2.6813271618763043040E-3}
           ($0478,$410E,$13CD,$BF2E),  {-2.2947197478731854413E-4}
           ($4426,$1643,$B0F3,$3F49),  {+7.8403348427447529072E-4}
           ($472B,$23A3,$5276,$3F12),  {+6.9893322606231933383E-5}
           ($8B37,$C8EF,$7F6B,$BF43),  {-5.9502375540563298261E-4}
           ($52B8,$8860,$C968,$BEF8),  {-2.3638488095017589101E-5}
           ($3EEE,$9C6D,$6BAA,$3F47)); {+7.1473913781436109426E-4}
const
  XMAX = 143.016; {max. argument for direct power}
var
  SP: array[0..8] of double absolute SPHex;
var
  v,w,y: double;
begin
  {Ref: Cephes [7], function stirf in file ldouble/gammal.c}
  w := 1.0/x;
  w := 1.0 + w*PolEval(w, SP, 9);
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
function sfd_signgamma(x: double): double;
  {-Return sign(gamma(x)), useless for 0 or negative integer}
begin
  if IsNanOrInfD(x) or (x>0.0) or (frac(0.5*floord(x))=0.0) then sfd_signgamma := 1.0
  else sfd_signgamma := -1.0;
end;


{---------------------------------------------------------------------------}
function sfd_gaminv_small(x: double): double;
  {- 1/gamma(x) for small argument values, |x| < 0.03125}
const
  {1/gamma(x) = x*P(x), 0 < x < 0.03125, peak relative error 4.2e-23}
  SPHex: array[0..8] of THexDblW = (
           ($0000,$0000,$0000,$3FF0),   { 1.000000000000000000000E+0}
           ($B619,$FC6F,$788C,$3FE2),   { 5.772156649015328608253E-1}
           ($FA2F,$026A,$FCF4,$BFE4),   {-6.558780715202540684668E-1}
           ($4D7E,$8FA2,$815E,$BFA5),   {-4.200263503403344054473E-2}
           ($A2AD,$20AE,$5123,$3FC5),   { 1.665386113720805206758E-1}
           ($82E2,$FB9D,$9AF0,$BFA5),   {-4.219773360705915470089E-2}
           ($FDF2,$1D3B,$B4B6,$BF83),   {-9.622023360406271645744E-3}
           ($D6A0,$E9D9,$9358,$3F7D),   { 7.220599478036909672331E-3}
           ($DA77,$BCBA,$8FC4,$BF53));  {-1.193945051381510095614E-3}

  {1/gamma(-x) = x*P(x), 0 < x < 0.03125, peak relative error 5.16e-23}
  SNHex: array[0..8] of THexDblW = (
           ($0000,$0000,$0000,$BFF0),   {-1.000000000000000000000E+0}
           ($B619,$FC6F,$788C,$3FE2),   { 5.772156649015328608727E-1}
           ($FA2C,$026A,$FCF4,$3FE4),   { 6.558780715202536547116E-1}
           ($4690,$8FA2,$815E,$BFA5),   {-4.200263503402112910504E-2}
           ($EDAF,$20BA,$5123,$BFC5),   {-1.665386113944413519335E-1}
           ($3BAF,$FA28,$9AF0,$BFA5),   {-4.219773343731191721664E-2}
           ($1E06,$0DE3,$B4A7,$3F83),   { 9.621911155035976733706E-3}
           ($AD13,$BE3B,$9398,$3F7D),   { 7.220837261893170325704E-3}
           ($5BCC,$3EE0,$91B7,$3F52));  { 1.133374167243894382010E-3}
var
  SP: array[0..8] of double absolute SPHex;
  SN: array[0..8] of double absolute SNHex;
var
  p: double;
begin
  {Ref: Cephes [7], function gammal/label small in file ldouble/gammal.c}
  if x=0.0 then sfd_gaminv_small := 0.0
  else begin
    if x<0.0 then begin
      x := -x;
      p := PolEval(x, SN, 9);
    end
    else p := PolEval(x, SP, 9);
    sfd_gaminv_small := x*p;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_gamma_medium(x: double): double;
  {-Return gamma(x), |x| <= 13, x negative integer produces div by 0}
const
  {gamma(x+2) = P(x)/Q(x), 0 <= x <= 1, peak relative error = 1.83e-20}
  PHex: array[0..7] of THexDblW = (
          ($0000,$0000,$0000,$3FF0),   { 1.000000000000000000009E+0}
          ($3665,$D903,$CF42,$3FEA),   { 8.378004301573126728826E-1}
          ($75F2,$1C84,$3A99,$3FD7),   { 3.629515436640239168939E-1}
          ($D6B3,$8751,$7E91,$3FBC),   { 1.113062816019361559013E-1}
          ($32D0,$2CCA,$6D16,$3F98),   { 2.385363243461108252554E-2}
          ($EAF8,$E2E6,$C378,$3F70),   { 4.092666828394035500949E-3}
          ($05FF,$6B7D,$C5C6,$3F3D),   { 4.542931960608009155600E-4}
          ($E448,$7B47,$1645,$3F06));  { 4.212760487471622013093E-5}

  QHex: array[0..8] of THexDblW = (
          ($0000,$0000,$0000,$3FF0),   { 9.999999999999999999908E-1}
          ($D8FD,$AAE5,$8F9F,$3FDA),   { 4.150160950588455434583E-1}
          ($56EF,$5A67,$B789,$BFCC),   {-2.243510905670329164562E-1}
          ($D306,$B006,$B9BA,$BFA7),   {-4.633887671244534213831E-2}
          ($3121,$F78F,$671A,$3F9C),   { 2.773706565840072979165E-2}
          ($9625,$BFAF,$11EB,$BF4A),   {-7.955933682494738320586E-4}
          ($0A78,$D223,$47B4,$BF54),   {-1.237799246653152231188E-3}
          ($5E06,$5BB8,$C1D4,$3F2E),   { 2.346584059160635244282E-4}
          ($BD0B,$4D05,$4CE2,$BEED));  {-1.397148517476170440917E-5}
var
  PX: array[0..7] of double absolute PHex;
  QX: array[0..8] of double absolute QHex;
var
  z: double;
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
    sfd_gamma_medium := z/sfd_gaminv_small(x);
  end
  else begin
    {finish reduction to [2,3)}
    while x < 2.0 do begin
      z := z/x;
      x := x+1.0;
    end;
    if x=2.0 then sfd_gamma_medium := z
    else begin
      x := x-2.0;
      sfd_gamma_medium := z * PolEval(x,PX,8) / PolEval(x,QX,9);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_lngcorr(x: double): double;
  {-Return lngamma correction term lngamma(x) - ((x-0.5)*ln(x) - x + ln(sqrt(2*Pi)), x>=8}
const
  { ln gamma(x) = (x-0.5)*ln(x) - x + ln(sqrt(2*Pi) + 1/x * A(1/x^2) }
  { x >= 8, peak relative error 1.51e-21}
  AHex: array[0..6] of THexDblW = (
          ($5554,$5555,$5555,$3FB5),   { 8.333333333333331447505E-2}
          ($750A,$16C0,$C16C,$BF66),   {-2.777777777750349603440E-3}
          ($1CBF,$1246,$01A0,$3F4A),   { 7.936507795855070755671E-4}
          ($F9D1,$89D3,$8130,$BF43),   {-5.952345851765688514613E-4}
          ($5A27,$9255,$911A,$3F4B),   { 8.412723297322498080632E-4}
          ($6088,$B420,$D0A7,$BF5E),   {-1.880801938119376907179E-3}
          ($811B,$3859,$0252,$3F74));  { 4.885026142432270781165E-3}
var
  PA: array[0..6] of double absolute AHex;
const
  xbig = 5.48e8;
  xmax = 1e307;
begin
  {Ref: Cephes [7], part function lgaml in file ldouble/gammal.c}
  if x < XLGAE then begin
    {this is not really used but return 'correct' result}
    sfd_lngcorr := sfd_lngamma(x) - ((x-0.5)*ln(x) - x + LnSqrt2Pi);
  end
  else if x >= xmax then sfd_lngcorr := 0.0
  else if x >= xbig then sfd_lngcorr := 1.0/(12.0*x)
  else sfd_lngcorr := PolEval(1.0/(x*x),PA,7)/x;
end;


{---------------------------------------------------------------------------}
function sfd_gamma(x: double): double;
  {-Return gamma(x), x <= MAXGAMD; NAN/RTE if x is a non-positive integer}
var
  i: integer;
  q,z: double;
begin
  {Ref: Cephes [7], function gammal in file ldouble/gammal.c}
  if IsNanorInfD(x) then begin
    sfd_gamma := x;
    exit;
  end;

  if (x<=0.0) and (frac(x)=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_gamma := NaN_d;
    exit;
  end;

  if x>MAXGAMD then begin
    sfd_gamma := PosInf_d;
    exit;
  end;

  q := abs(x);
  if q <= 13.0 then sfd_gamma := sfd_gamma_medium(x)
  else begin
    if x < 0.0 then begin
      if q <= MAXGAMD-1 then begin
        i := trunc(q);
        z := q-i;
        if z > 0.5 then z := q - (i+1);
        z := abs(q * sinPi(z)) * stirf(q);
        if z <= Pi/MaxDouble then begin
          sfd_gamma := copysignd(PosInf_d, sfd_signgamma(x));
          exit;
        end;
        sfd_gamma := sfd_signgamma(x)*Pi/z;
      end
      else begin
        {Use lngamma, there are only few non-zero cases, most values underflow:}
        {e.g. gamma(-171.624376956302712) = 1.101257447554064E-310 (denormal)}
        z := exp(sfd_lngammas(x,i));
        sfd_gamma := z*i;
      end;
    end
    else sfd_gamma := stirf(x);
  end
end;


{---------------------------------------------------------------------------}
const
  FTabH: array[0..25] of THexDblW = (  {Table of factorials}
           ($0000,$0000,$0000,$3FF0),  {1.000000000000000000000}
           ($0000,$0000,$0000,$3FF0),  {1.000000000000000000000}
           ($0000,$0000,$0000,$4000),  {2.000000000000000000000}
           ($0000,$0000,$0000,$4018),  {6.000000000000000000000}
           ($0000,$0000,$0000,$4038),  {2.400000000000000000000E+1}
           ($0000,$0000,$0000,$405E),  {1.200000000000000000000E+2}
           ($0000,$0000,$8000,$4086),  {7.200000000000000000000E+2}
           ($0000,$0000,$B000,$40B3),  {5.040000000000000000000E+3}
           ($0000,$0000,$B000,$40E3),  {4.032000000000000000000E+4}
           ($0000,$0000,$2600,$4116),  {3.628800000000000000000E+5}
           ($0000,$0000,$AF80,$414B),  {3.628800000000000000000E+6}
           ($0000,$0000,$08A8,$4183),  {3.991680000000000000000E+7}
           ($0000,$0000,$8CFC,$41BC),  {4.790016000000000000000E+8}
           ($0000,$C000,$328C,$41F7),  {6.227020800000000000000E+9}
           ($0000,$2800,$4C3B,$4234),  {8.717829120000000000000E+10}
           ($0000,$7580,$0777,$4273),  {1.307674368000000000000E+12}
           ($0000,$7580,$0777,$42B3),  {2.092278988800000000000E+13}
           ($0000,$ECD8,$37EE,$42F4),  {3.556874280960000000000E+14}
           ($0000,$CA73,$BEEC,$4336),  {6.402373705728000000000E+15}
           ($9000,$3068,$02B9,$437B),  {1.216451004088320000000E+17}
           ($5A00,$BE41,$E1B3,$43C0),  {2.432902008176640000000E+18}
           ($C620,$E9B5,$283B,$4406),  {5.109094217170944000000E+19}
           ($F06C,$6159,$7752,$444E),  {1.124000727777607680000E+21}
           ($A4CE,$35F8,$E5C3,$4495),  {2.585201673888497664000E+22}
           ($7B9A,$687A,$6C52,$44E0),  {6.204484017332394393600E+23}
           ($6121,$C33F,$A940,$4529)); {1.551121004333098598400E+25}
var
  FacTab: array[0..25] of double absolute FTabH; {0..22 exact, others correctly rounded}


{---------------------------------------------------------------------------}
function sfd_fac(n: integer): double;
  {-Return the factorial n!, n<MAXGAMD-1; INF if n<0}
begin
  if n>25 then sfd_fac := stirf(n+1)
  else begin
    if n<0 then sfd_fac := PosInf_d
    else sfd_fac := FacTab[n];
  end;
end;


{---------------------------------------------------------------------------}
function lngamma_small(x,xm1,xm2: double; useln1p: boolean): double;
  {-Return ln(gamma(x), 0<=x<8; xm1=x-1, xm2=x-2 are supplied for increased precision}
  { useln1p should be true if xm1 is supposed to be more accurate than x-1}
var
  r,t,z: double;
const
  Y20: THexDblW = ($0000,$0000,$58EC,$3FC4); {0.158963680267333984375 = 333371/2097152}
  P20: array[0..6] of double = (
         -0.180355685678449379109e-1,
          0.25126649619989678683e-1,
          0.494103151567532234274e-1,
          0.172491608709613993966e-1,
         -0.259453563205438108893e-3,
         -0.541009869215204396339e-3,
         -0.324588649825948492091e-4);
  Q20: array[0..7] of double = (
          1.0,
          0.196202987197795200688e1,
          0.148019669424231326694e1,
          0.541391432071720958364e0,
          0.988504251128010129477e-1,
          0.82130967464889339326e-2,
          0.224936291922115757597e-3,
         -0.223352763208617092964e-6);
const
  Y10: THexDblW = ($0000,$0000,$E6A2,$3FE0);  {0.52815341949462890625 = 1107618/2097152}
  P10: array[0..6] of double = (
          0.490622454069039543534e-1,
         -0.969117530159521214579e-1,
         -0.414983358359495381969e0,
         -0.406567124211938417342e0,
         -0.158413586390692192217e0,
         -0.240149820648571559892e-1,
         -0.100346687696279557415e-2);
  Q10: array[0..6] of double = (
          1.0,
          0.302349829846463038743e1,
          0.348739585360723852576e1,
          0.191415588274426679201e1,
          0.507137738614363510846e0,
          0.577039722690451849648e-1,
          0.195768102601107189171e-2);
const
  Y15: THexDblW = ($0000,$0000,$EDDA,$3FDC);  {0.452017307281494140625 = 947949/2097152}
  P15: array[0..5] of double = (
         -0.292329721830270012337e-1,
          0.144216267757192309184e0,
         -0.142440390738631274135e0,
          0.542809694055053558157e-1,
         -0.850535976868336437746e-2,
          0.431171342679297331241e-3);
  Q15: array[0..6] of double = (
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
  if x<eps_d then begin
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
    r := PolEval(xm2, P20, 7) / PolEval(xm2, Q20, 8);
    t := xm2*(x+1.0);
    lngamma_small := z + (t*double(Y20) + t*r);
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
      r := PolEval(xm1, P10, 7) / PolEval(xm1, Q10, 7);
      t := xm1*xm2;
      lngamma_small := z + (t*double(Y10) + t*r);
    end
    else begin
      {Use the following form: lngamma(x) = (2-x)*(1-x)*(Y15 + R(2-x))}
      r := PolEval(-xm2, P15, 6) / PolEval(-xm2, Q15, 7);
      t := xm1*xm2;
      lngamma_small := z + (t*double(Y15) + t*r);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_lngamma(x: double): double;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, NAN/RTE if x is a non-positive integer}
var
  t: double;
begin
  if IsNaND(x) then begin
    sfd_lngamma := NaN_d;
    exit;
  end;
  if x<=0.0 then begin
    if frac(x)=0.0 then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_lngamma := NaN_d;
    end
    else begin
      {use reflection formula}
      t := x*sinPi(x);
      sfd_lngamma := ln(abs(Pi/t)) - sfd_lngamma(abs(x));
    end;
    exit;
  end;
  if (frac(x)=0.0) and (x<=26.0) then begin
    sfd_lngamma := ln(FacTab[trunc(x)-1]);
  end
  else if x <= XLGAE then begin
    sfd_lngamma := lngamma_small(x,x-1.0,x-2.0,false)
  end
  else begin
    {Stirling approximation for x > 8}
    sfd_lngamma := (x-0.5)*ln(x) - x + LnSqrt2Pi + sfd_lngcorr(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_lngammas(x: double; var s: integer): double;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, s=-1,1 is the sign of gamma}
begin
  if IsNaND(x) then begin
    sfd_lngammas := NaN_d;
    s := 0;
    exit;
  end;
  sfd_lngammas := sfd_lngamma(x);
  if (x>0.0) or (frac(0.5*floord(x))=0.0) then s := 1 else s := -1;
end;


{---------------------------------------------------------------------------}
function sfd_lngamma1p(x: double): double;
  {-Return ln(|gamma(1+x)|) with increased accuracy for x near 0}
var
  t: double;
begin
  t := 1.0+x;
  if (t < 0.0) or (t > XLGAE) then sfd_lngamma1p := sfd_lngamma(t)
  else sfd_lngamma1p := lngamma_small(t,x,x-1.0,true)
end;


{---------------------------------------------------------------------------}
function sfd_gamma1pm1(x: double): double;
  {-Return gamma(1+x)-1 with increased accuracy for x near 0}
begin
  if abs(x)<=eps_d then sfd_gamma1pm1 := -EulerGamma*x
  else if (x<-0.5) or (x>2.0) then sfd_gamma1pm1 := sfd_gamma(1.0+x)-1.0
  else if x>0.0 then sfd_gamma1pm1 := expm1(lngamma_small(x+1.0,x,x-1.0,true))
  else sfd_gamma1pm1 := expm1(lngamma_small(x+2.0,x+1.0,x,false) - ln1p(x));
end;


{---------------------------------------------------------------------------}
function sfd_rgamma(x: double): double;
  {-Return the reciprocal gamma function rgamma = 1/gamma(x)}
begin
  if IsNanOrInfD(x) then begin
    if x=PosInf_d then sfd_rgamma := 0.0
    else sfd_rgamma := NaN_d;
  end
  else if abs(x)<0.03125 then sfd_rgamma := sfd_gaminv_small(x)
  else if (x<0.0) and (frac(x)=0.0) then sfd_rgamma := 0.0
  else begin
    {
    if (x<-0.5) or (x>XLGAE) then sfd_rgamma := sfd_signgamma(x)*exp(-sfd_lngamma(x))
    else sfd_rgamma := 1.0/sfd_gamma_medium(x);
    }
    if (x<-0.5) or (x>=MAXGAMD) then sfd_rgamma := sfd_signgamma(x)*exp(-sfd_lngamma(x))
    else if x <= 13.0 then sfd_rgamma := 1.0/sfd_gamma_medium(x)
    else sfd_rgamma := 1.0/stirf(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_gstar(x: double): double;
  {-Return Temme's gstar(x) = gamma(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x)), x>0.}
  { For large x the asymptotic expansion is gstar(x) = 1 + 1/12x + O(1/x^2)}
var
  a,g: double;
const
  {Chebyshev approximations calculated with Maple V R4 with Digits :=30;}
  {Coefficients for T(0,.) doubled, Hexdouble with t_rcalc command dh.}
  {Definition: gstar := x-> GAMMA(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x));}

  {Range 1<=x<= 3; y := t -> t+2;}
  {chebyshev(y(t)^2*(gstar(y(t))-1-1/(12*y(t))), t=-1..1, 0.5e-20);}
  gs1n = 25;
  gs1h : array[0..gs1n-1] of THexDblW = (
           ($9745,$8AE2,$757D,$3F70),  { 4.0182975773513416011E-3 }
           ($31B1,$7A55,$04E5,$3F47),  { 7.0248799299004439952E-4 }
           ($FDB0,$C0EE,$E62F,$BF24),  {-1.5944798403032040299E-4 }
           ($9A12,$FFBE,$01BB,$3F02),  { 3.4345197462305170773E-5 }
           ($5AF8,$4208,$6D65,$BEDD),  {-7.0160213288181423859E-6 }
           ($DDAF,$1C0D,$A795,$3EB6),  { 1.3503205169277659442E-6 }
           ($7990,$B052,$267A,$BE90),  {-2.4065836641022231323E-7 }
           ($F7F6,$6232,$6C8A,$3E64),  { 3.8042639505142093235E-8 }
           ($0C2C,$8E25,$029C,$BE34),  {-4.6589880584690519361E-9 }
           ($9B86,$DB13,$9673,$3DE2),  { 1.3524335604397218593E-10}
           ($3D0B,$8FF6,$896C,$3DE8),  { 1.7852880977454193496E-10}
           ($C837,$24EE,$0147,$BDD8),  {-8.7329651525165672802E-11}
           ($31A3,$E17A,$9CAF,$3DC0),  { 3.0217158796552201555E-11}
           ($E9A3,$91F7,$2671,$BDA4),  {-9.1632366291192802894E-12}
           ($767F,$D4F5,$CD55,$3D86),  { 2.5922976124555003060E-12}
           ($FD33,$7C4A,$BD27,$BD68),  {-7.0312136560940930257E-13}
           ($2A56,$3EB4,$1C7A,$3D4A),  { 1.8553152123744796025E-13}
           ($CD3C,$26F0,$0C21,$BD2B),  {-4.8045799980445695572E-14}
           ($43B2,$808D,$A74C,$3D0B),  { 1.2280625885209937902E-14}
           ($5972,$D84D,$0405,$BCEC),  {-3.1103690941234212527E-15}
           ($91A8,$BE1E,$33D8,$3CCC),  { 7.8277734240176922378E-16}
           ($8A34,$B381,$4484,$BCAC),  {-1.9614622328820518705E-16}
           ($6FCE,$E627,$40AE,$3C8C),  { 4.9010567735366476461E-17}
           ($D090,$80D7,$3075,$BC6C),  {-1.2225157066341578107E-17}
           ($FBCC,$42C5,$19DD,$3C4C)); { 3.0467200401436511414E-18}
        (* ($366C,$C645,$0273,$BC2C),  {-7.5920116218493385258E-19}
           ($346E,$1DCD,$E82D,$3C0B),  { 1.8910477040223097493E-19}
           ($2610,$7CB2,$E6B5,$BBEB),  {-4.7266482810237095702E-20}
           ($00CC,$8087,$985D,$3BCB),  { 1.1687011937049585628E-20}
           ($20E0,$CB61,$6E5D,$BBAC)); {-3.0102617821884302625E-21} *)
var
  gs1a: array[0..gs1n-1] of double absolute gs1h;

const
  gs2n = 26; {chebyshev(gstar(t), t=3..8, 0.5e-20);}
  gs2h : array[0..gs2n-1] of THexDblW = (
           ($40CF,$A83F,$4635,$4000),  { 2.0342820305159885945    }
           ($775B,$6F17,$F8F5,$BF80),  {-8.2873510863767290591E-3 }
           ($567B,$1393,$63E6,$3F60),  { 2.0007604294772196608E-3 }
           ($204D,$AF28,$9DB7,$BF3F),  {-4.8242315747915459550E-4 }
           ($7E4E,$3195,$74E0,$3F1E),  { 1.1618250245736783082E-4 }
           ($25D7,$C5CA,$4E66,$BEFD),  {-2.7948623357453920514E-5 }
           ($C1C8,$295F,$2B57,$3EDC),  { 6.7160841776690551181E-6 }
           ($C2B0,$FFFC,$0CA4,$BEBB),  {-1.6122694432402287168E-6 }
           ($7B56,$FE71,$F329,$3E99),  { 3.8668303975989479479E-7 }
           ($0FF3,$EEAB,$DFA4,$BE78),  {-9.2661419574693789809E-8 }
           ($A261,$743E,$D2B9,$3E57),  { 2.2187030226642767308E-8 }
           ($A7A8,$75D7,$CCF0,$BE36),  {-5.3086653662215709374E-9 }
           ($ADAA,$0E5F,$CEB9,$3E15),  { 1.2693642192007070286E-9 }
           ($B6A6,$E8E8,$D869,$BDF4),  {-3.0333999974117580609E-10}
           ($BD00,$EBDC,$EA42,$3DD3),  { 7.2450652218555779523E-11}
           ($EB4E,$1A56,$046F,$BDB3),  {-1.7296152050716179453E-11}
           ($8313,$9580,$2706,$3D92),  { 4.1273879591581087674E-12}
           ($AB68,$AB0B,$5210,$BD71),  {-9.8456023555037793839E-13}
           ($CFFB,$E136,$8585,$3D50),  { 2.3478568909854971813E-13}
           ($909B,$E7EE,$82A3,$BD2F),  {-5.5973560909513051708E-14}
           ($6210,$73C4,$0AAB,$3D0E),  { 1.3341185335663677403E-14}
           ($4247,$D425,$A2D1,$BCEC),  {-3.1792362333608899741E-15}
           ($2ABA,$7D34,$4AB7,$3CCB),  { 7.5750134834532103750E-16}
           ($BF50,$E1EF,$01F0,$BCAA),  {-1.8046385098946929066E-16}
           ($EAF0,$1F51,$C809,$3C88),  { 4.2988857605477843878E-17}
           ($92D1,$5A41,$9C84,$BC67)); {-1.0239810102203542236E-17}
        (* ($02EC,$D2DD,$7EE1,$3C46),  { 2.4389814490930843665E-18}
           ($42BB,$B103,$6E9D,$BC25),  {-5.8091806323722222819E-19}
           ($AD04,$8EF0,$6B32,$3C04),  { 1.3836276558996084587E-19}
           ($C6E9,$C753,$741A,$BBE3),  {-3.2955567815269497085E-20}
           ($EC44,$8C84,$88D1,$3BC2)); { 7.8496438299511273955E-21} *)
var
  gs2a: array[0..gs2n-1] of double absolute gs2h;
const
  gstar1: THexDblW = ($65FE,$309E,$59DB,$3FF1);   {1.0844375514192275466}
begin
  if x>=2.41e7 then begin
    {very large x, use first terms of Stirling formula}
    a := 1.0/(12.0*x);
    sfd_gstar := 1.0 + a*(1.0 + 0.5*a);
  end
  else if x>8.0 then begin
    {medium range, recycle sfd_lngcorr}
    sfd_gstar := exp(sfd_lngcorr(x));
  end
  else if x>=3.0 then begin
    {3 <= x < 8}
    sfd_gstar := CSEvalD((2.0*x-11.0)/5.0, gs2a, gs2n);
  end
  else if x>=1.0 then begin
    {1 <= x < 3}
    a :=  CSEvalD(x-2.0, gs1a, gs1n);
    sfd_gstar := a/sqr(x) + 1.0 + 1.0/(12.0*x);
  end
  else if x>=MinDouble then begin
    {0 < x < 1: Use Temme's [26] recursion formula for gstar(x+1),}
    {but avoid recursive call if x+1 = 1 within machine precision.}
    g := 1.0 + x;
    if g=1.0 then g := double(gstar1)
    else g := sfd_gstar(g);
    a := 1.0 + 1.0/x;
    g := g*sqrt(a);
    a := ln(a)*x - 1.0;
    sfd_gstar := g*exp(a);
  end
  else sfd_gstar := 1.0/(Sqrt_TwoPi*sqrt(x));
end;


{---------------------------------------------------------------------------}
function sfd_lnfac(n: longint): double;
  {-Return ln(n!), INF if n<0}
begin
  if n>25 then sfd_lnfac := sfd_lngamma(n+1.0)
  else begin
    if n<0 then sfd_lnfac := PosInf_d
    else sfd_lnfac := ln(FacTab[n]);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_dfac(n: integer): double;
  {-Return the double factorial n!!, n<=300; INF for even n<0}
var
  k: integer;
  t: double;
begin
  if n<0 then begin
    if odd(n) then begin
      if n=-1 then sfd_dfac := 1
      else begin
        {(-2k-1)!! = (-1)^k/(2k-1)!!}
        t := 1.0/sfd_dfac(-n-2);
        if n and 2 = 0 then sfd_dfac := -t
        else sfd_dfac := t;
      end;
    end
    else sfd_dfac := PosInf_d;
  end
  else if n>300 then sfd_dfac := PosInf_d
  else if odd(n) then begin
    { (2k+1)!! = (2k+1)!/(2^k * k!) }
    if n <= 25 then begin
      k := n div 2;
      sfd_dfac := FacTab[n]/ldexpd(FacTab[k],k);
    end
    else begin
      {Here n >= 27 and the argument of stirf is >= 14.5}
      n := succ(n) div 2;
      t := ldexpd(stirf(n+0.5), n);
      sfd_dfac := int(t/SqrtPi+0.5);
    end;
  end
  else begin
    { n!! = 2^k * k!  with k = n div 2}
    n := n div 2;
    sfd_dfac := ldexpd(sfd_fac(n),n);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_binomial(n,k: integer): double;
  {-Return the binomial coefficient 'n choose k'}
var
  t: double;
  i: integer;
const
  km = 4000;

  function lnf0(m: integer): double;
    {-Return quick and dirty approximation of ln(m!)}
  begin
    lnf0 := (0.5+m)*ln(1.0+m)-m;
  end;
begin

  (* Accuracy of sfd_binomial, random values 0<=k<=n, computed with
     t_binomd, compiled with D17(XE3)-64 and executed under Win7/64.

     Range (n)      Samples   RMS           Peak
        0..  100    20000     1.3237E-16    8.8818E-16  for (68,41)
      100..  170    20000     2.5254E-16    1.1102E-15  for (147,137)
      170..  800    20000     7.9306E-16    3.5527E-15  for (776,460)
      800.. 2000    10000     8.7523E-15    1.8874E-13  for (1275,307)
     2000..10000    10000     6.8175E-15    1.0913E-13  for (8723,138)
  *)

  {If n<0 and k>=0, we use the identity 1.2.6 (17) from D.E. Knuth: The art}
  {of computer programming, Vol. 1, Fundamental algorithms, 3rd ed., 1997.}
  {[*]   binomial(n,k) = (-1)^k * binomial(k-n-1,k)}

  if (k=0) or (k=n) then sfd_binomial := 1.0
  else if (k=1) or (k=pred(n)) then sfd_binomial := n
  else if k<0 then begin
    if (n>=0) or (n<k) then sfd_binomial := 0.0
    else begin
      {Here k <= n < 0:  We use binomial(n,k) = binomial(n,n-k) which is}
      {valid for all integer k, and then use [*] since n<0 and n-k >= 0:}
      {binomial(n,k) = (-1)^(n-k)*binomial(n-k-n-1,n-k)}
      i := n-k;
      t := sfd_binomial(-succ(k),i);
      if (t<>0.0) and odd(i) then sfd_binomial := -t
      else sfd_binomial := t;
    end;
  end
  else if n<0 then begin
    {Here n < 0 and k >= 0. Apply Knuth's [*] identity}
    t := sfd_binomial(k-n-1,k);
    if (t<>0.0) and odd(k) then sfd_binomial := -t
    else sfd_binomial := t;
  end
  else if k>n then sfd_binomial := 0.0
  else begin
    {here n > 2 and  1 < k < n-1}
    if n<MAXGAMD-1 then begin
      t := sfd_fac(n);
      t := t/sfd_fac(n-k);
      t := t/sfd_fac(k);
    end
    else begin
      {if k>n/2 use symmetry to reduce k}
      if k>n-k then k := n-k;
      {all (n,k) with n <= 1029 will NOT overflow,}
      {this saves three logarithm calculations.   }
      if n<=1029 then t:=0
      else begin
        {Calculate rough estimate t}
        t := lnf0(k);
        t := t + lnf0(n-k);
        t := lnf0(n)-t;
      end;
      if t>750.0 then begin
        {don't waste time (especially for double), if result is INF}
        t := PosInf_d
      end
      else if (t<700) and (k<=km) then begin
        {Compute n*(n-1).../k! if k is 'small' and no risk of overflow}
        t := n;
        inc(n);
        for i:=2 to k do t := t*((n-i)/i);
      end
      else begin
        {Here n > 2*km and k > km, or a very large result}
        {Use lnbeta function (increased inaccuracy for 'large' values)}
        t := sfd_lnbeta(k,n-k+1);
        t := k*exp(t);
        if t<MinDouble then begin
          {for t > 0.25*MinDouble inversion will not overflow}
          if 4.0*t<MinDouble then t := PosInf_d
          else t := 1.0/t;
        end
        else t := 1.0/t;
      end;
    end;
    {round to nearest integer value}
    sfd_binomial := int(t+0.5);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_lnbinomial(n,k: longint): double;
  {-Return ln(binomial(n,k)), n >= k >= 0}
var
  t: double;
begin
  if (n<k) or (k<0) then begin
    sfd_lnbinomial := Nan_d;
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end
  else if (k=n) or (k=0)   then sfd_lnbinomial := 0.0
  else if (k=1) or (k=n-1) then sfd_lnbinomial := ln(n)
  else if n<MAXGAMD then begin
    t := sfd_fac(n);
    t := t/sfd_fac(n-k);
    sfd_lnbinomial := ln(t/sfd_fac(k));
  end
  else begin
    {t := sfd_lngamma(n-k+1) + sfd_lngamma(k+1);
    sfd_lnbinomial := sfd_lngamma(n+1)-t;}
    sfd_lnbinomial := -(ln(k) + sfd_lnbeta(k, n-k+1));
  end;
end;


{---------------------------------------------------------------------------}
function sfd_gdr_pos(x,d: double): double;
  {-Return gamma(x)/gamma(x+d), x>0, x+d > 0, accurate even for |d| << |x|}
var
  xgh,res: double;
begin
  {$ifdef debug}
    if (x <= 0.0) and (x+d <= 0.0) then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_gdr_pos := NaN_d;
      exit;
    end;
  {$endif}
  if x<eps_d then begin
    {x small, Gamma(x)=1/x. Result = rGamma(x+d)/x}
    {but avoid spurious underflow if x+d > MAXGAMX}
    d := x+d;
    if d<MAXGAMD then sfd_gdr_pos := sfd_rgamma(d)/x
    else begin
      {Note very accurate but better than underflow}
      res := sfd_lngamma(d);
      res := res + ln(x);
      sfd_gdr_pos := exp(-res);
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
  sfd_gdr_pos := lanczos(x,false)/lanczos(x+d,false)*res;
end;


{---------------------------------------------------------------------------}
function sfd_gamma_delta_ratio(x,d: double): double;
  {-Return gamma(x)/gamma(x+d), accurate even for |d| << |x|}
var
  s,t,y: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(d) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_gamma_delta_ratio := NaN_d;
    exit;
  end;
  y := x+d;
  if x=y then sfd_gamma_delta_ratio := 1.0
  else if (x>0.0) and (y>0.0) then sfd_gamma_delta_ratio := sfd_gdr_pos(x,d)
  else begin
    if (y<=0.0) and (frac(y)=0.0) then begin
      {y=x+d is a non-positive integer}
      if (x>0.0) or (frac(x)<>0.0) then begin
        {gamma(x+d) = inf, gamma(x)<>0, return zero}
        sfd_gamma_delta_ratio := 0.0;
        exit;
      end;
      {both x,y non-positive integers, no need to calculate sin for reflection}
      t := sfd_gdr_pos(1.0-x,-d);
      if frac(0.5*x)<>frac(0.5*y) then t := -t;
      sfd_gamma_delta_ratio := 1.0/t;
      exit;
    end;
    if (x<0.0) and (y<0.0) then begin
      {Here both x and y are negative, use reflection formula for gammas}
      {gamma(x)/gamma(x+d) = sinPi(x+d)/sinP(x)/ [gamma(1-x)/gamma(1-x-d]}
      t := sfd_gdr_pos(1.0-x,-d);
      t := sinPi(x)*t;
      y := sinPi(y);
      sfd_gamma_delta_ratio := y/t;
    end
    else begin
      {Remaing cases: use lngamma for both arguments}
      s := sfd_signgamma(x)*sfd_signgamma(y);
      y := sfd_lngamma(y);
      t := sfd_lngamma(x);
      t := exp(t-y);
      sfd_gamma_delta_ratio := s*t;
    end
  end;
end;


{---------------------------------------------------------------------------}
function sfd_gamma_ratio(x,y: double): double;
  {-Return gamma(x)/gamma(y)}
begin
  sfd_gamma_ratio := sfd_gamma_delta_ratio(x,y-x)
end;


{---------------------------------------------------------------------------}
function sfd_pochhammer(a,x: double): double;
  {-Return the Pochhammer symbol (a)_x = gamma(a+x)/gamma(a).}
  { Accuracy is reduced if x or x+a are near negative integers.}
var
  t: double;
  it,ia: boolean;
  i: integer;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(a) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_pochhammer := NaN_d;
    exit;
  end;
  if x=0.0 then sfd_pochhammer := 1.0
  else begin
    t := a+x;
    ia := (a<=0.0) and (frac(a)=0.0);
    it := (t<=0.0) and (frac(t)=0.0);
    if it or ia then begin
      if it and ia then begin
        {both a and a+x are negative integers}
        t := sfd_gamma_delta_ratio(1-a,-x);
        if frac(0.5*x)<>0 then t := -t;
        sfd_pochhammer := t;
      end
      else if ia then sfd_pochhammer := 0.0
      else begin
        {$ifopt R+}
          if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
        {$endif}
        sfd_pochhammer := NaN_d;
      end;
    end
    else begin
      if (frac(x)=0.0) and (x>0.0) and (x<=100.0) then begin
        {x is a small positive integer, use simple multiply loop}
        t := a;
        for i:=1 to trunc(x)-1 do t := t*(a+i);
        sfd_pochhammer := t;
      end
      else begin
        t := sfd_gamma_delta_ratio(a,x);
        if t<>0 then sfd_pochhammer := 1.0/t
        else sfd_pochhammer := PosInf_d;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_poch1(a,x: double): double;
  {-Return (pochhammer(a,x)-1)/x, psi(a) if x=0; accurate even for small |x|}
const
  sqtbig = 0.13684286665667228e154;  {=1/sqrt(24*MinDouble)}
  lneps  = -36.736800569677101399;   {=ln(eps_d/2)}
var
  gbern : array[1..21] of double;
var
  t,b,bp,v,v2,lv,q,poly,res,term: double;
  nterms,incr,j,k: integer;
begin
  {Ref: W. Fullerton [14] and [20], files dpoch1.f}
  if IsNanOrInfD(x) or IsNanOrInfD(a) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_poch1 := NaN_d;
    exit;
  end;

  if x=0.0 then begin
    sfd_poch1 := sfd_psi(a);
    exit;
  end;

  t := abs(x);
  q := abs(a);
  if (t > 0.1*q) or (t*ln(maxd(q,2.0)) > 0.1) then begin
    res := sfd_pochhammer(a,x);
    sfd_poch1 := (res-1.0)/x;
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
        sfd_write_debug_str('*** sfd_poch1: nterms too big');
      {$endif}
      nterms := 20;
    end;
    v2  := 1.0/sqr(v);
    res := 0.5*(x+1.0);
    gbern[1] := 1.0;
    gbern[2] := -res/12.0;
    term := v2;
    poly := gbern[2]*term;
    for k:=2 to nterms do begin
      t := 0.0;
      for j:=1 to k do t := t + double(BoFHex[k-j])*gbern[j];
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

  if bp=a then sfd_poch1 := res
  else begin
    {Use reflection formula to get poch1(a,x) with a<-0.5 from poch1(bp,x)}
    {q = cot(Pi*b)}
    sincosPi(b,t,q);
    q := q/t;
    t := sinPi(x)/x;
    t := t*q - 2.0*sqr(sinPi(0.5*x))/x;
    sfd_poch1 := t + (1.0 + x*t)*res;
  end;
end;


{---------------------------------------------------------------------------}
const
  BNNHex: array[0..8] of THexDblW = (      {bernreal(2k+2)/(2k+2), k=0..8}
            ($5555,$5555,$5555,$3FB5),   { 8.33333333333333333333333E-2}
            ($1111,$1111,$1111,$BF81),   {-8.33333333333333333333333E-3}
            ($0410,$1041,$4104,$3F70),   { 3.96825396825396825396825E-3}
            ($1111,$1111,$1111,$BF71),   {-4.16666666666666666666667E-3}
            ($1F08,$F07C,$07C1,$3F7F),   { 7.57575757575757575757576E-3}
            ($5996,$9599,$9959,$BF95),   {-2.10927960927960927960928E-2}
            ($5555,$5555,$5555,$3FB5),   { 8.33333333333333333333333E-2}
            ($5E5E,$5E5E,$5E5E,$BFDC),   {-4.43259803921568627450980E-1}
            ($E6E8,$9B9F,$6E7F,$4008));  { 3.05395433027011974380395   }
var
  BNN: array[0..8] of double absolute BNNHex;


{---------------------------------------------------------------------------}
function sfd_psi(x: double): double;
  {-Return the psi (digamma) function of x, INF if x is a non-positive integer}

const
  PPHex: array[0..5] of THexDblW = (
           ($38AC,$C8EE,$F72B,$3FEE),  {+0.967672245447621170427444761712   }
           ($B512,$3247,$4E89,$3FF2),  {+1.14417380943934132177443086511    }
           ($ECAE,$3DEB,$2D2F,$3FDE),  {+0.471507845373433246720039154386   }
           ($F054,$0CBB,$E3A3,$3FB4),  {+0.815984636391829369377116423465e-1}
           ($2BBB,$EA37,$8877,$3F76),  {+0.550124017486102827456755050064e-2}
           ($82C1,$61F6,$35AA,$3F19)); {+0.961671107604462646458076964138e-4}

  PQHex: array[0..6] of THexDblW = (
           ($0000,$0000,$0000,$3FF0),  {+1.00000000000000000000000000000    }
           ($AA55,$54D9,$3D3F,$3FFA),  {+1.63995297569876643247016291140    }
           ($2B06,$BAE3,$0E46,$3FEF),  {+0.970492711080937113643931469693   }
           ($A3DC,$F6A1,$9F0D,$3FD0),  {+0.259707918978674869994308187160   }
           ($BAB5,$78F5,$37E7,$3FA0),  {+0.316765151172736832756607905517e-1}
           ($FBF9,$4D5B,$CF4D,$3F58),  {+0.151426838914381207406214536130e-2}
           ($B8F6,$92C9,$D3AE,$3EF1)); {+0.170010400090619970545999542634e-4}
const
  x0hHex: THexDblW = ($0000,$0000,$62D0,$3FF7);  {1.4616241455078125}
  x0lHex: THexDblW = ($DC35,$7C7E,$C6AD,$3EE0);  {7.9994605498412626595423257213e-6}
var
  PP : array[0..5] of double absolute PPHex;
  PQ : array[0..6] of double absolute PQHex;
  x0h: double absolute x0hHex;
  x0l: double absolute x0lHex;
var
  p,q,nz,s,w,y,z: double;
  i,n: integer;
  subneg: boolean;

begin

  if IsNanOrInfD(x) then begin
    if x=PosInf_d then sfd_psi := PosInf_d
    else sfd_psi := NaN_d;
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
    p := floord(q);
    if p=q then begin
      sfd_psi := PosInf_d;
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
  if (x <= 12.0) and (x=floord(x)) then begin
    y := 0.0;
    n := trunc(x);
    for i:=n-1 downto 1 do begin
      {FPC nonsense: 1.0/i is evaluated in single precision!!!}
      x := i;
      y := y + 1.0/x;
    end;
    y := y - EulerGamma;
  end
  else if abs(x-x0h)<0.2 then begin
    {Pad‚ approximation for x in (x0-0.2, x0+0.2), abs. error < 3e-19}
    {where x0=1.46163214496836234126.., psi(x0)=0}
    x := x-x0h;
    x := x-x0l;
    y := x*(PolEval(x,PP,6)/PolEval(x,PQ,7));
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
      y := z*PolEval(z,BNN,8);
    end;
    y := ln(s) - (0.5/s) - y - w;
  end;
  if subneg then y := y-nz;
  sfd_psi := y;
end;


{---------------------------------------------------------------------------}
function sfd_psistar(x: double): double;
  {-Return psi(x) - ln(x), x > 0}
var
  s,w,y,z: double;
begin
  if IsNanD(x) then sfd_psistar := NaN_d
  else if x < 2.0 then sfd_psistar := sfd_psi(x) - ln(x)
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
      y := z*PolEval(z,BNN,9);
    end;
    sfd_psistar := (w - y) - 0.5/s;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_trigamma(x: double): double;
  {-Return the trigamma function psi'(x), INF if x is a negative integer}
var
  s: double;
begin
  {trigamma = psi(1,x) = zetah(2,x) for x >0}
  if IsNanOrInfD(x) then begin
    if x=PosInf_d then sfd_trigamma := 0.0
    else sfd_trigamma := NaN_d;
  end
  else if x>0.0 then begin
    if x<1e17 then sfd_trigamma := sfd_zetah(2.0,x)
    else sfd_trigamma := 1.0/x;
  end
  else if frac(x)=0.0 then sfd_trigamma := PosInf_d
  else begin
    {Reflection formula HMF [1] 6.4.7 with n=1}
    {psi(1,1-x) + psi(1,x) = Pi^2*(1+cot(Pi*x)^2) = Pi^2/sin(Pi*x)^2}
    s := sinPi(x);
    x := 1.0 - x;
    if x>1e17 then x := 1.0/x else x := sfd_zetah(2.0,x);
    sfd_trigamma := PiSqr/sqr(s) - x;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_tetragamma(x: double): double;
  {-Return the tetragamma function psi''(x), NAN/RTE if x is a negative integer}
var
  c,s: double;
const
  PP3H: THexDblW = ($9D7C,$5938,$019B,$404F);  {2*Pi^3}
begin
  {tetragamma = psi(2,x) = -2*zetah(3,x) for x >0}
  if IsNanOrInfD(x) then begin
    if x=PosInf_d then sfd_tetragamma := 0.0
    else sfd_tetragamma := NaN_d;
  end
  else if x>0.0 then begin
    if x>1e17 then begin
      x := 1.0/x;
      sfd_tetragamma := -x*x;
    end
    else sfd_tetragamma := -2.0*sfd_zetah(3.0,x)
  end
  else if frac(x)=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_tetragamma := NaN_d;
  end
  else begin
    {Reflection formula HMF [1] 6.4.7 with n=2, 2nd derivative of cot(Pi*x) with Maple}
    {psi(2,1-x) - psi(2,x) = 2*Pi^3*cot(Pi*x)*(1+cot(Pi*x)^2)}
    sincosPi(x,s,c);
    c := c/s;
    c := c*(1.0 + c*c);
    s := -2.0*sfd_zetah(3.0,1.0-x);
    sfd_tetragamma := s - double(PP3H)*c;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_pentagamma(x: double): double;
  {-Return the pentagamma function psi'''(x), INF if x is a negative integer}
var
  c,s: double;
const
  PP4H: THexDblW =  ($0826,$8C29,$5A2E,$4068);  {2*Pi^4}
begin
  {pentagamma = psi(3,x) = 6*zetah(4,x) for x >0}
  if IsNanOrInfD(x) then begin
    if x=PosInf_d then sfd_pentagamma := 0.0
    else sfd_pentagamma := NaN_d;
  end
  else if x>0.0 then begin
    if x>1e17 then begin
      x := 1.0/x;
      sfd_pentagamma := 2.0*x*x*x;
    end
    else sfd_pentagamma := 6.0*sfd_zetah(4.0,x);
  end
  else if frac(x)=0.0 then sfd_pentagamma := PosInf_d
  else begin
    {Reflection formula HMF [1] 6.4.7 with n=3, 3rd derivative of cot(Pi*x) with Maple}
    {psi(3,1-x) + psi(3,x) = 2*Pi^4*(1+4*cot(Pi*x)^2 + 3*cot(Pi*x)^4)}
    sincosPi(x,s,c);
    c := sqr(c/s);
    c := 1.0 + c*(4.0 + 3.0*c);
    s := 6.0*sfd_zetah(4.0,1.0-x);
    sfd_pentagamma := double(PP4H)*c - s;
  end;
end;


const
  NXPGPMAX = 100;  {max. n for polygamma(n,x) with x < 0 via cot polynomial}

{---------------------------------------------------------------------------}
function polygam_cotpoly(n: integer; x: double): double;
  {-Polygamma for negative x and 1 <= n <= NXPGPMAX}
type
  TCArray = array[-1..NXPGPMAX+2] of double;
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
  c,s,t: double;
  i,k: integer;
  pa: PCArray;
  sa: word;
begin
  sincosPi(x,s,c);
  if (n<1) or (n>NXPGPMAX) or (x>=0.0) or (s=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    polygam_cotpoly := NaN_d;
    exit;
  end;

  {Reflection formula: HMF[1] (6.4.7) or  NIST[30], (5.15.6)}
  {psi(n,1-x) + (-1)^(n+1)*psi(n,x) = (-1)^n*Pi*d^n/dx^n(cot(Pi*x))}

  {psi(n,x) = - (-1)^(n+1)*psi(n,1-x) + (-1)^(n+n+1)*Pi*d^n/dx^n(cot(Pi*x))}
  {psi(n,x) = (-1)^n*psi(n,1-x) - Pi*d^n/dx^n(cot(Pi*x))}

  if (c=0.0) and (n and 1 = 0) then begin
    {early exit for cot=0: for even n cot recflection term -t*c below will be 0}
    polygam_cotpoly := sfd_polygamma(n,1.0-x);
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
    sa := (n+4)*sizeof(double);
    {allocate and clear to zero}
    pa := calloc(sa);
    if pa=nil then begin
      {$ifopt R+}
        {Heap overflow error}
        RunError(203);
      {$endif}
      polygam_cotpoly := NaN_d;
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
  s := sfd_polygamma(n,1.0-x);
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
function polygam_negx(n: integer; x: double): double;
  {-Polygamma for finite negative x and 1 <= n <= MAXGAMD}
var
  n1,k: integer;
  z,s,t,a,b,d: double;
const
  tmin = 13;   {should > 12}
  tmax = 70.0; {must by <= NXPGPMAX}
begin
  z  := frac(x);
  if (z=0.0) or (x>0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    polygam_negx := NaN_d;
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
    until abs(d)<=eps_d*abs(s);
    if odd(n) then s := -s;
  end;
  t := sfd_zetah(n1, 1.0-x);
  a := s+t;
  b := sfd_fac(n);
  if abs(a) <= Maxdouble/b then polygam_negx := -b*a
  else polygam_negx := copysignd(PosInf_d, -a);
end;


{---------------------------------------------------------------------------}
function polygam_sum(n: integer; x: double): double;
  {-Direct summation polygamma(n,x) = n!/x^(n+1)*(1 + sum(1/(1+k/x)^(n+1), k=1...))}
const
  KMAX = 200;
var
  s,t,e: double;
  k,m: integer;
begin
  m := -(n+1);
  e := 0.125*eps_d;
  s := 0.0;
  k := 1;
  repeat
    t := power(1.0+k/x, m);
    s := s + t;
    inc(k);
  until (t<=e) or (k>KMAX);
  if k>KMAX then begin
    {no convergence}
    polygam_sum := NaN_d;
  end
  else begin
    t := sfd_lngamma(n+1) - (n+1)*ln(x);
    if t+ln1p(s) >= ln_MaxDbl then polygam_sum := PosInf_d
    else polygam_sum := (1.0+s)*exp(t);
  end;
end;


{---------------------------------------------------------------------------}
function polygam_ase(n: integer; x: double): double;
  {-Polygamma asymptotic expansion}
var
  z,f,t,b: double;
  k,j: integer;
const
  K2MAX =  120; {doubled maximum index in asymptotic sum}
  NPMAX = 1000; {Max. n for calculation of (n-1)!/x^n as product}
  TFMAX = 0.798336124e293;  {~ 2*MaxDouble*eps_d, used for 'no-convergence' check}
begin
  f := sfd_lngamma(n);
  z := n*ln(x);
  t := f-z;
  if t > ln_MaxDbl then z := PosInf_d
  else if t < ln_MinDbl then z := 0.0
  else begin
    {calculate f = (n-1)!/x^n}
    if (abs(z)<709.0) and (n<MAXGAMD-1) then begin
      f := sfd_fac(n-1)*power(x,-n);
    end
    else begin
      if n <= NPMAX then begin
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
      b := sfd_bernoulli(k);
      if (k>K2MAX) or (abs(f) > abs(TFMAX/t/b)) then begin
        {Convergence failed if k too large or a summand ~> Maxdouble}
        {Try direct summation, complete failure if result is NaN}
        z := polygam_sum(n,x);
        {$ifopt R+}
          if (IsNaND(z)) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
        {$endif}
        polygam_ase := z;
        exit;
      end;
      f := f*t;
      t := f*b;
      z := z + t;
      k := k+2;
    until abs(t)<abs(z)*eps_d;
  end;
  polygam_ase := z
end;


{---------------------------------------------------------------------------}
function sfd_polygamma(n: integer; x: double): double;
  {-Return the polygamma function: n'th derivative of psi; n>=0, x>0 for n>MAXGAMD.}
  { Note: The accuracy may be reduced for n>=MAXGAMD due to ln/exp operations.}
var
  z,f,t: double;
  k: integer;
const
  NPMAX = 500; {Max. n > MAXGAMD for calculation of n! as product}
begin
  if n<4 then begin
    case n of
        0: sfd_polygamma := sfd_psi(x);
        1: sfd_polygamma := sfd_trigamma(x);
        2: sfd_polygamma := sfd_tetragamma(x);
        3: sfd_polygamma := sfd_pentagamma(x);
      else begin
             {RTE or NAN if n < 0}
             {$ifopt R+}
               if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
             {$endif}
             sfd_polygamma := NaN_d;
           end;
    end;
  end
  else if IsNaNOrInfD(x) then begin
    if x=PosInf_d then begin
      if n>0 then sfd_polygamma := 0.0
      else sfd_polygamma := PosInf_d
    end
    else sfd_polygamma := NaN_d;
  end
  else if x > 0.0 then begin
    {polygamma(n,x) = (-1)^(n+1)*n!*sum(1/(x+k)^(n+1), k>=0), HMF[1] 6.4.10}
    {or with the Hurwitz zeta:  polygamma(n,x) = (-1)^(n+1)*n!*zetah(n+1,x)}
    f := ln(x);
    t := n*f;
    if -t-f < ln_MinDbl then begin
      {force asymptotic formula, zetah may have }
      {some errors in the "near underflow range"}
      z := 0.0;
    end
    else z := sfd_zetah(n+1,x);
    if z=0.0 then begin
      {polygamma(n,x) = (-1)^(n+1)*n!*sum(1/(x+k)^(n+1), k>=0), HMF[1] 6.4.10}
      {or with the Hurwitz zeta:  polygamma(n,x) = (-1)^(n+1)*n!*zetah(n+1,x)}
      if x>0.25*n then z := polygam_ase(n,x)
      else z := polygam_sum(n,x)
    end
    else if n<MAXGAMD-1 then begin
      {n! is OK, but n!*zetah may overflow for small x}
      f := sfd_fac(n);
      t := ln(f)+ln(z);
      if t<ln_MaxDbl then z := z*f
      else z := PosInf_d;
    end
    else begin
      f := sfd_lngamma(n+1);
      t := f + ln(z);
      if t > ln_MaxDbl then z := PosInf_d
      else if t < ln_MinDbl then z := 0.0
      else begin
        if n>MAXGAMD+NPMAX then z := exp(t)
        else begin
          {example: sfd_polygamma(240,5)}
          z := sfd_fac(170)*z;
          for k:=171 to n do z := z*k;
        end;
      end;
    end;
    if odd(n) then sfd_polygamma := z
    else sfd_polygamma := -z;
  end
  else if (x<0.0) and (n>MAXGAMD) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_polygamma := NaN_d;
  end
  else begin
    {x<0 and n<=nmax}
    sfd_polygamma := polygam_negx(n,x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_batemang(x: double): double;
  {-Return the Bateman function G(x) = psi((x+1)/2) - psi(x/2), x<>0,-1,-2,...}
const
  z0 = 1e5;
  z1 = 1e19;
var
  s,t,r,z: double;
  k: integer;
begin
  if IsNanOrInfD(x) or ((x<=0.0) and (frac(x)=0.0)) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_batemang := NaN_d;
    exit;
  end;
  if abs(x)<=sqrt_epsh then begin
    {G(x)  = 2/x - 2*ln(2) + Pi^2/6*x + O(x^2)}
    sfd_batemang := 2.0*(1.0/x - ln2);
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
    until t < eps_d*s;
    s := s/z + r;
  end;
  if x<0.0 then begin
    {Reflection formula Erdelyi et al.[50], 1.8(8)}
    t := TwoPi/sinPi(x);
    s := t - s;
  end;
  sfd_batemang := s;
end;


{---------------------------------------------------------------------------}
function sfd_lnbeta(x,y: double): double;
  {-Return the logarithm of |beta(x,y)|=|gamma(x)*gamma(y)/gamma(x+y)|}
var
  c,s: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(x) then begin
    sfd_lnbeta := NaN_d;
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
    c := sfd_lngcorr(x) + sfd_lngcorr(y) - sfd_lngcorr(s);
    sfd_lnbeta := -0.5*ln(y) + lnsqrt2pi + c + (x-0.5)*ln(x/s) + y*ln1p(-x/s);
  end
  else if (y >= XLGAE) and (s >= XLGAE) then begin
    {Fixed in 1.02.01: check s >= XLGAE, otherwise lnbeta(-10.1,10) fails}
    {x < XLGAE, y >= XLGAE}
    c := sfd_lngcorr(y) - sfd_lngcorr(s);
    sfd_lnbeta := sfd_lngamma(x) + c + x - x*ln(s) + (y-0.5)*ln1p(-x/s)
  end
  else begin
    {x,y < XLGAE; or s<=0}
    sfd_lnbeta := sfd_lngamma(x) + sfd_lngamma(y) - sfd_lngamma(s);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_beta(x,y: double): double;
  {-Return the function beta(x,y)=gamma(x)*gamma(y)/gamma(x+y)}
var
  r,t: double;
begin
  {If one arg is 1, return the inverse of the other because}
  {gamma(1)*gamma(y)/gamma(y+1) = 1/y, analogue for y=1}
  if IsNanOrInfD(x) or IsNanOrInfD(x) then sfd_beta := NaN_d
  else if x=1.0 then sfd_beta := 1.0/y
  else if y=1.0 then sfd_beta := 1.0/x
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
        sfd_beta := 0.0;
      end
      else begin
        t := sfd_gamma(x);
        r := sfd_pochhammer(y,x);
        sfd_beta := t/r;
      end;
      exit;
    end;
    if t=x then begin
      {gamma(x)/gamma(x+y)=1}
      sfd_beta := sfd_gamma(y);
    end
    else begin
      r := sfd_lnbeta(x,y);
      t := exp(r);
      if (x<0.0) or (y<0.0) then begin
        {get signs only if necessary}
        r := sfd_signgamma(x)*sfd_signgamma(y)/sfd_signgamma(x+y);
        sfd_beta := r*t;
      end
      else sfd_beta := t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_taylor(x: double; n: integer): double;
  {-Return the Taylor coefficient x^n/n!, n>=0}
var
  z,f: double;
  k: integer;
const
  zmin = -346.04333733782; { > ln(MinDouble*100!)}
begin
  if IsNanOrInfD(x) or (n<0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_taylor := NaN_d;
    exit;
  end;
  if n<=0 then sfd_taylor := 1.0
  else if (n=1) or (x=0.0) then sfd_taylor := x
  else begin
    z := n*ln(abs(x));
    if z<ln_MinDbl then sfd_taylor := 0.0
    else if (n<=100) and (z<ln_MaxDbl) and (z>zmin) then begin
      {do not waste time for sfd_lnfac(n)}
      f := 1.0;
      for k:=1 to n do f := f*(x/k);
      sfd_taylor := f;
    end
    else begin
      f := sfd_lnfac(n);
      f := z-f;
      if f<ln_MinDbl then sfd_taylor := 0.0
      else if f>ln_MaxDbl then sfd_taylor := copysignd(PosInf_d,x)
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
        sfd_taylor := f;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function lnbg1p(x: double): double;
  {-Return ln(BarnesG(1+x) via Taylor series, -0.5 <= x <= 0.5}
var
  s,t,y: double;
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
    t := sfd_zetaint(k)/(k+1)*y;
    s := s + t;
    if abs(t) <= eps_d*abs(s) then break;
    y := -y*x;
  end;
  lnbg1p := s;
end;


(*
{---------------------------------------------------------------------------}
function lnbg1p(x: double): double;
  {-Return ln(BarnesG(1+x) via Taylor series, -0.5 <= x <= 0.5}
var
  s,t,y: double;
  k: integer;
const
  c1 = 0.4189385332046727418;   {(ln(2*Pi)-1)/2}
  c2 = -0.2113921675492335697;  {(Eulergamma-1)/2}
  KMAX = 100;
begin
  {Ferreira/Lopez: An Asymptotic Expansion of the Double Gamma Function,}
  {formula (2): About 10% less iterations, sligtly less accurate}
  s := x*(c1 + c2*x) + x*sfd_lngamma1p(x);
  y := x*x*x;
  for k:=2 to KMAX do begin
    t := sfd_zetaint(k)/(k*(k+1))*y;
    s := s - t;
    if abs(t) <= eps_d*abs(s) then break;
    y := -y*x;
  end;
  lnbg1p := s;
end;
*)


{---------------------------------------------------------------------------}
function sfd_lnbg(x: double): double;
  {-Return ln(BarnesG(x)), real part for x < 0}
const
  c1 = 0.9189385332046727418;  {-Zeta'(0) = ln(2*Pi)/2}
  c2 = 0.1654211437004509292;  {-Zeta'(1)}
  XAE  = 10;    {Min x for asymptotic expansion}
  KMAX = 64;    {Max terms in AE}
var
  g,s,y,z,t: double;
  k: integer;
begin
  if IsNanOrInfD(x) then begin
    sfd_lnbg := NaN_d;
    exit;
  end;
  if (frac(x)=0.0) and (x<=28.0) then begin
    if x<=0.0 then sfd_lnbg := NegInf_d
    else begin
      {Logarithmic version of functional equation}
      s := 0.0;
      for k:=2 to round(x)-2 do s := s + sfd_lnfac(k);
      sfd_lnbg := s;
    end;
    exit;
  end;
  if x < -0.5 then begin
    {Use reflection formula, Adamchik[78], (35)}
    z := x-2.0;
    g := sfd_lnbg(-z);
    s := (x-1.0) * ln(abs(sinPi(z)/Pi));
    {We can use frac(z) because Cl_2 has period 2*Pi}
    t := sfd_cl2(frac(z)*TwoPi)/TwoPi;
    sfd_lnbg := g - s - t;
  end
  else if x < XAE then begin
    if x < 0.5 then begin
      {lnG(x) = lnG(x+1) - lnGamma(x)}
      sfd_lnbg := lnbg1p(x) - sfd_lngamma(x);
    end
    else if x <= 1.5 then begin
      {x is already in [0.5 , 1.5]}
      sfd_lnbg := lnbg1p(x-1.0);
    end
    else begin
      {Use functional equation to shift x downto [0.5 , 1.5]}
      g := 0.0;
      while x>1.5 do begin
        {lnG(x) = lnG(x-1) + lnGamma(x-1)}
        x := x - 1.0;
        g := g + sfd_lnGamma(x);
      end;
      sfd_lnbg := lnbg1p(x-1.0) + g;
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
      t := sfd_bernoulli(2*k+2);
      t := t / (4*k*(k+1)*y);
      if abs(t) <= eps_d*abs(s) then break;
      y := y * z;
      s := s + t;
    end;
    sfd_lnbg := g + s;
  end;
end;




end.
