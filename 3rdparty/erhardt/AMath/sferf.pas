unit sfErf;

{Common code for special functions: Error function and related}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}

(*************************************************************************

 DESCRIPTION   :  Common code for special functions: Error function and related

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in amath_info.txt/references

                  [7] Cephes Mathematical Library, Version 2.8
                      http://www.moshier.net/#Cephes or http://www.netlib.org/cephes/
                 [13] W.H. Press et al, Numerical Recipes in C, 2nd ed., Cambridge, 1992
                 [14] SLATEC Common Mathematical Library, Version 4.1, July 1993
                      (general purpose mathematical and statistical routines written in Fortran 77)
                      http://www.netlib.org/slatec
                 [19] Boost C++ Libraries, Release 1.42.0, 2010.
                      http://www.boost.org/
                 [20] Special functions by Wayne Fullerton,
                      http://www.netlib.org/fn
                      Almost identical to the FNLIB subset of SLATEC [14]
                 [22] A.J. MacLeod, MISCFUN: A software package to compute uncommon special functions.
                      ACM Trans. on Math. Soft. 22 (1996), pp.288-301.
                      Fortran source: http://netlib.org/toms/757
                 [30] [NIST]: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST Handbook
                      of Mathematical Functions, Cambridge, 2010. Online resource: NIST Digital
                      Library of Mathematical Functions, http://dlmf.nist.gov/
                 [33] http://functions.wolfram.com/: Formulas and graphics about
                      mathematical functions for the mathematical and scientific
                      community and/or http://mathworld.wolfram.com/ ("/the web's
                      most extensive mathematical resource/")
                 [48] D. Dijkstra, A Continued Fraction Expansion for a Generalization of
                      Dawson's Integral, Mathematics of Computation, Vol.31, 503-510 (1977).
                      Available as http://doc.utwente.nl/75001/1/Dijkstra77continued.pdf
                 [49] W. Gautschi, Evaluation of the Repeated Integrals of the Coerror
                      Function, ACM TOMS, Vol.3, No.3, 1977, pp.240-252
                      Fortran source available from http://netlib.org/toms/521
                 [70] A.R. Didonato, Significant Digit Computation of the Elliptical
                      Coverage Function, NSWC TR 90 513, 1990. Available from
                      http://www.dtic.mil/dtic/tr/fulltext/u2/a230523.pdf
                 [72] M. Patefield, D. Tandy, Fast and Accurate Calculation of Owen's T Function,
                      Journal of Statistical Software V5 N5, 2000, pp. 1-25.
                      https://www.jstatsoft.org/article/view/v005i05/t.pdf


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  17.08.10  W.Ehrhardt  Common version number after split
 1.00.01  09.09.10  we          sfc_fresnel handles NANs

 1.04.00  12.02.11  we          Goodwin-Staton integral sfc_gsi
 1.04.01  13.02.11  we          sfc_expint3
 1.04.02  14.02.11  we          imaginary error function erfi

 1.07.00  10.06.11  we          sfc_erfce, renamed sfc_erfc_scaled

 1.08.00  15.08.11  we          Range check sfc_erf(c)_inv

 1.13.00  15.06.12  we          sfc_gendaw (generalized Dawson integral)
 1.13.01  16.06.12  we          sfc_inerfc (repeated integrals of erfc)
 1.13.02  26.06.12  we          sfc_erfg (generalized error function)

 1.15.00  15.02.13  we          Fix sfc_expint3 for 0.4e-6 .. 6e-6

 1.18.00  19.05.13  we          Prevent some wrong compiler optimizations for div by 3

 1.20.00  18.08.13  we          sfc_erf_p, sfc_erf_q, sfc_erf_z

 1.37.00  14.02.16  we          Constants for trivial ranges of erf/erfc/erfce

 1.41.00  07.07.17  we          new sfc_fresnel_fg, improved sfc_fresnel
 1.41.01  10.07.17  we          sfc_fresnel_fg for x<0
 1.41.02  18.07.17  we          sfc_erfh
 1.41.03  19.07.17  we          sfc_erf2

 1.42.00  17.08.17  we          sfc_gsi for x < 0

 1.43.00  01.12.17  we          sfc_erfi_inv
 1.43.01  02.12.17  we          sfc_expint3 with sfc_igammal

 1.44.00  17.01.18  we          minor changes in sfc_erf_p/q
 1.44.01  18.01.18  we          Owen's T function sfc_owent
 1.44.01  23.01.18  we          sfc_erfce_inv

 1.50.00  14.08.18  we          IsNanOrInf check sfc_erfi_inv
 1.50.01  23.08.18  we          Marcum Q from sfMisc

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


function sfc_dawson(x: extended): extended;
  {-Return Dawson's integral: dawson(x) = exp(-x^2)*integral(exp(t^2), t=0..x)}

function sfc_gendaw(p,x: extended): extended;
  {-Return the generalized Dawson integral F(p,x) = exp(-x^p)*integral(exp(t^p), t=0..x); x,p >= 0}

function sfc_erf(x: extended): extended;
  {-Return the error function erf(x) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..x)}

function sfc_erfc(x: extended): extended;
  {-Return the complementary error function erfc(x) = 1-erf(x)}

function sfc_erfce(x: extended): extended;
  {-Return the exponentially scaled complementary error function erfce(x) = exp(x^2)*erfc(x)}

function sfc_erfg(p,x: extended): extended;
  {-Return the generalized error function integral(exp(-t^p), t=0..x); x,p >= 0}

function sfc_erfi(x: extended): extended;
  {-Return the imaginary error function erfi(x) = erf(ix)/i}

function sfc_erfh(x,h: extended): extended;
  {-Accurately compute erf(x+h) - erf(x-h)}

function sfc_erf2(x1,x2: extended):extended;
  {-Accurately compute erf(x2) - erf(x1)}

function sfc_erf_inv(x: extended): extended;
  {-Return the inverse function of erf, erf(erf_inv(x)) = x, -1 < x < 1}

function sfc_erfc_inv(x: extended): extended;
  {-Return the inverse function of erfc, erfc(erfc_inv(x)) = x, 0 < x < 2}

function sfc_erfce_inv(x: extended): extended;
  {-Return the functional inverse of erfce, erfce(erfce_inv(x)) = x, x > 0}

function sfc_erfi_inv(y: extended): extended;
  {-Return the inverse imaginary error function: erfi(erf_inv(y))=y}

function sfc_erf_p(x: extended): extended;
  {-Return the probability function erf_p = integral(exp(-t^2/2)/sqrt(2*Pi), t=-Inf..x)}

function sfc_erf_q(x: extended): extended;
  {-Return the probability function erf_q = integral(exp(-t^2/2)/sqrt(2*Pi), t=x..Inf)}

function sfc_erf_z(x: extended): extended;
  {-Return the probability function erf_z = exp(-x^2/2)/sqrt(2*Pi)}

function sfc_inerfc(n: integer; x: extended): extended;
  {-Return the repeated integrals of erfc, n >= -1; scaled with exp(x^2) for x>0}

function sfc_expint3(x: extended): extended;
  {-Return the integral(exp(-t^3), t=0..x), x >= 0}

procedure sfc_fresnel(x: extended; var s,c: extended);
  {-Return the Fresnel integrals S(x)=integral(sin(Pi/2*t^2),t=0..x) and C(x)=integral(cos(Pi/2*t^2),t=0..x)}

procedure sfc_fresnel_fg(x: extended; var f,g: extended);
  {-Return the Fresnel auxiliary functions f,g}

function sfc_gsi(x: extended): extended;
  {-Return the Goodwin-Staton integral = integral(exp(-t*t)/(t+x), t=0..Inf), x <> 0}

function sfc_mq(m: integer; a,b: extended): extended;
  {-Return the generalized Marcum Q function Q(m,a,b), a,b >= 0}

function sfc_owent(h,a: extended): extended;
  {-Return Owen's T function T(h,a)}

{#Z+}
procedure sfc_erfc_iusc(x: extended; var p,q: extended);
  {-(internal) unevaluated scaled erfc, p/q := exp(x^2)*erfc(x), 1<=x<=128}
{#Z-}


implementation

uses
  AMath,
  sfBasic,  {Basic common code}
  sfBessel, {Bessel functions}
  sfExpInt, {Exponential integrals and related}
  sfGamma,  {Gamma function and related}
  sfGamma2; {inverse/incomplete Gamma/Beta}

{---------------------------------------------------------------------------}
function erf_small(x: extended): extended;
  {-Return erf(x) for |x| <= 1}
const
  {erf(x) = x*T(x^2)/U(x^2),  0 <= x <= 1, Peak relative error 7.6e-23}
  THex: array[0..6] of THexExtW = (
          ($c446,$6bab,$0b2a,$86d0,$4013),  {1.104385395713178565288E6}
          ($7a56,$e45a,$a4bd,$975b,$4010),  {1.549905740900882313773E5}
          ($b954,$a987,$c60c,$bc83,$400e),  {4.825977363071025440855E4}
          ($6118,$6059,$9093,$a757,$400a),  {2.677472796799053019985E3}
          ($9517,$4e93,$540e,$8f97,$4007),  {2.871822526820825849235E2}
          ($3128,$c337,$3716,$ace5,$4001),  {5.402980370004774841217E0}
          ($fd7a,$3a1a,$705b,$e0c4,$3ffb)); {1.097496774521124996496E-1}
  UHex: array[0..6] of THexExtW = (
          ($71a7,$1cad,$012e,$eef3,$4012),  {9.787360737578177599571E5}
          ($9ad5,$1aef,$45b1,$e25e,$4011),  {4.636021778692893773576E5}
          ($481d,$445b,$c807,$c232,$400f),  {9.942956272177178491525E4}
          ($ffe8,$9cac,$3b84,$c2ac,$400c),  {1.245905812306219011252E4}
          ($71ac,$b12f,$21ca,$f2e2,$4008),  {9.715333124857259246107E2}
          ($3453,$1f8e,$f688,$b507,$4004),  {4.525777638142203713736E1}
          ($0000,$0000,$0000,$8000,$3fff)); {1.000000000000000000000E0}
var
  z: extended;
  T: array[0..6] of extended absolute THex;
  U: array[0..6] of extended absolute UHex;
begin
  {Ref: Cephes [7], file ldouble\ndtrl.c}
  if x=0.0 then erf_small := x
  else begin
    z := x*x;
    erf_small := x*PolEvalx(z,T,7)/PolEvalx(z,U,7);
  end;
end;


{---------------------------------------------------------------------------}
procedure sfc_erfc_iusc(x: extended; var p,q: extended);
  {-(internal) unevaluated scaled erfc, p/q := exp(x^2)*erfc(x), 1<=x<=128}
const
  {erfc(x) = exp(-x^2)*P(1/x)/Q(1/x), 1/8<=1/x<=1, Peak relative error 5.8e-21}
  PHex: array[0..9] of THexExtW = (
          ($333b,$d9e6,$d404,$986f,$bfee),   {-9.085943037416544232472E-6}
          ($0144,$489e,$be68,$9c31,$4011),   {3.198859502299390825278E5}
          ($3b58,$3da2,$af02,$9780,$4015),   {4.964439504376477951135E6}
          ($5840,$554d,$37a3,$9239,$4018),   {3.833161455208142870198E7}
          ($4e13,$caee,$9e31,$b258,$401a),   {1.870095071120436715930E8}
          ($ada8,$356a,$4982,$94a6,$401c),   {6.234814405521647580919E8}
          ($b6d0,$c92b,$5417,$acb1,$401d),   {1.448651275892911637208E9}
          ($d025,$cfd5,$8494,$88d3,$401e),   {2.295563412811856278515E9}
          ($df23,$d843,$4032,$8881,$401e),   {2.290171954844785638925E9}
          ($4bf0,$9ad8,$7a03,$86c7,$401d));  {1.130609921802431462353E9}
  QHex: array[0..10] of THexExtW = (
          ($24fa,$96f6,$7153,$8a6c,$4012),   {5.669830829076399819566E5}
          ($3708,$33b1,$07fa,$8644,$4016),   {8.799239977351261077610E6}
          ($c4cb,$305a,$bf78,$8220,$4019),   {6.822450775590265689648E7}
          ($a75d,$436f,$30dd,$a027,$401b),   {3.358653716579278063988E8}
          ($c39d,$e415,$c43d,$87c0,$401d),   {1.138778654945478547049E9}
          ($4991,$cfda,$52f1,$a2a9,$401e),   {2.729005809811924550999E9}
          ($00e7,$7595,$cd06,$88bb,$401f),   {4.588018188918609726890E9}
          ($8eae,$8dad,$6eb4,$9aa2,$401f),   {5.188672873106859049556E9}
          ($f817,$9128,$c0f8,$d48b,$401e),   {3.565928696567031388910E9}
          ($0e43,$302d,$79ed,$86c7,$401d),   {1.130609910594093747762E9}
          ($0000,$0000,$0000,$8000,$3fff));  {1.000000000000000000000E0}
const
  {erfc(x) = exp(-x^2)*1/x*R(1/x^2)/S(1/x^2), 1/128<=1/x<1/8, Peak relative error 1.9e-21}
  RHex: array[0..4] of THexExtW = (
          ($22cf,$c711,$6c5b,$dcfb,$3ff9),   {2.697535671015506686136E-2}
          ($521d,$8527,$3435,$8dc2,$3ffe),   {5.537445669807799246891E-1}
          ($0615,$4b00,$575f,$dc7b,$4000),   {3.445028155383625172464E0}
          ($4761,$613e,$df6d,$e58e,$4001),   {7.173690522797138522298E0}
          ($260a,$ab95,$2fc7,$e7c4,$4000));  {3.621349282255624026891E0}
  SHex: array[0..5] of THexExtW = (
          ($f5af,$2fb2,$1e57,$c3d7,$3ffa),   {4.781257488046430019872E-2}
          ($3684,$3798,$b793,$80b0,$3fff),   {1.005392977603322982436E0}
          ($b611,$8f76,$f020,$d255,$4001),   {6.572990478128949439509E0}
          ($55d5,$d300,$e71e,$f564,$4002),   {1.533713447609627196926E1}
          ($5de6,$17d7,$54d6,$aba9,$4002),   {1.072884067182663823072E1}
          ($0000,$0000,$0000,$8000,$3fff));  {1.000000000000000000000E0}
var
  PP: array[0..9]  of extended absolute PHex;
  QQ: array[0..10] of extended absolute QHex;
  RR: array[0..4]  of extended absolute RHex;
  SS: array[0..5]  of extended absolute SHex;
var
  y: extended;
begin
  {Ref: Cephes [7], file ldouble\ndtrl.c}
  y := 1.0/x;
  if x<8.0 then begin
    p := PolEvalx(y,PP,10);
    q := PolEvalx(y,QQ,11);
  end
  else begin
    q := y*y;
    p := y*PolEvalx(q,RR,5);
    q := PolEvalx(q,SS,6);
  end;
end;


const
  x0e = 6.5;         {erf(x)  = 1 for |x| > x0e}
  x1e = 106.75;      {erfc(x) = 0 for  x  > x1e}

{---------------------------------------------------------------------------}
function sfc_erf(x: extended): extended;
  {-Return the error function erf(x) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..x)}
begin
  if IsInf(x) then sfc_erf := copysign(1.0,x)
  else if abs(x) <= 1.0 then sfc_erf := erf_small(x)
  else if x > +x0e then sfc_erf := +1.0
  else if x < -x0e then sfc_erf := -1.0
  else sfc_erf := 1.0 - sfc_erfc(x);
end;


{---------------------------------------------------------------------------}
function sfc_erfc(x: extended): extended;
  {-Return the complementary error function erfc(x) = 1-erf(x)}
var
  p,q,a,z: extended;
begin
  {expx2 is used to suppress error amplification in computing exp(-x^2)}
  {Ref: Cephes [7], file ldouble\ndtrl.c}
  if IsInf(x) then begin
    sfc_erfc := 1.0-copysign(1.0,x);
    exit;
  end;
  a := abs(x);
  if a < 1.0 then sfc_erfc := 1.0-erf_small(x)
  else if x < -x0e then sfc_erfc := 2.0
  else if x >  x1e then sfc_erfc := 0.0
  else begin
    {Here x is in one of the ranges -x0e..-1 or 1..x1e}
    {Compute accurate z = exp(-a^2)}
    z := expx2(-a);
    {Compute p/q := exp(a^2)*erfc(a)}
    sfc_erfc_iusc(a,p,q);
    z := (z*p)/q;
    if x>=0 then sfc_erfc := z
    else sfc_erfc := 2.0-z;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_erfce(x: extended): extended;
  {-Return the exponentially scaled complementary error function erfce(x) = exp(x^2)*erfc(x)}
var
  y,z: extended;
const
  xmin = -106.5637380121098417;
begin
  if IsNan(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
     sfc_erfce := NaN_x;
    exit;
  end;
  y := abs(x);
  if y <= 1e-10 then sfc_erfce := 1.0 - 2.0*x/SqrtPi
  else if y <= 1.0 then begin
    z := 1.0-erf_small(x);
    sfc_erfce := exp(x*x)*z;
  end
  else if x < -x0e then begin
    {here erfc(x)=2.0 accurate to extended precision}
    if x < xmin then sfc_erfce := PosInf_x
    else sfc_erfce := 2.0*expx2(y);
  end
  else if y <= 128.0 then begin
    { -x0e < x <= 128: Recycle sfc_erfc_iusc}
    sfc_erfc_iusc(y,y,z);
    y := y/z;
    if x>0.0 then sfc_erfce := y
    else begin
      z := 2.0*expx2(-x);
      sfc_erfce := z-y;
    end;
  end
  else begin
    {x > 128: use asymptotic expansion}
    if x > 0.5e10 then z := 1.0
    else begin
      y := 0.5/sqr(x);
      {erfce = (1 - y + 3y^2 - 15y^3 + 105y^4 - 945y^5 + O(y^6))/x/sqrt(pi)}
      z := (105.0 - 945.0*y)*y - 15.0;
      z := (z*y + 3.0)*y - 1.0;
      z := z*y + 1.0;
    end;
    sfc_erfce := z/x/sqrtpi;
  end;
end;


{---------------------------------------------------------------------------}
function erfgsmallp(p,z,x: extended): extended;
  {-Return erfg(p,x), p small, z=x^p}
var
  a,t,s: extended;
begin
  {Special code for p < 1 avoiding scaling/rescaling with x^p. We have }
  {erfg(p,x) = igammal(a,z)*a  with  z = x^p, a = 1/p > 1, a > z + 1/4.}
  {Use non-normalised Temme/Gautschi code for P(a,z), cf. igam_ptaylor.}
  a := 1.0/p;
  t := 1.0;
  s := 1.0;
  repeat
    a := a+1.0;
    t := t*z/a;
    s := s+t;
  until t < s*eps_x;
  {Here igammal(a,z) = s*z^a*exp(-z)/a and therefore}
  {erfg(p,x) = igammal(a,z)*a = s*x*exp(-z)}
  erfgsmallp := s*x*exp(-z);
end;


{---------------------------------------------------------------------------}
function sfc_erfg(p,x: extended): extended;
  {-Return the generalized error function integral(exp(-t^p), t=0..x); x,p >= 0}
var
  t,u,z: extended;
const
  em1 = 0.36787944117144232159552;  {exp(-1)}
begin
  {Ref; NIST[30], 7.16}
  if x=0.0 then sfc_erfg := 0.0
  else if p<MinExtended then begin
    {t^p ~ 1}
    sfc_erfg := x*em1;
  end
  else if p=1.0 then sfc_erfg := -expm1(-x)
  else begin
    {sfc_erfg = igammal(t,x^p)*t}
    t := 1.0/p;
    u := p*ln(x);
    if u >= ln_MaxExt then begin
      {igammal(t,e^u)*t ~ gamma(t)*t = gamma(t+1)}
      sfc_erfg := sfc_gamma(t+1.0)
    end
    else if u <= ln_MinExt then begin
      {t^p ~ 0, exp(-t^p) ~ 1}
      sfc_erfg := x;
    end
    else begin
      {z = x^p = exp(u)}
      z := exp(u);
      if (p>1.0) and (z>ln_Maxext) then begin
        {if z^t*exp(-z)/gamma(t) underflows then igammal(t,z) ~ gamma(t)}
        if t*u - z - sfc_lngamma(t) < ln_MinExt then begin
          sfc_erfg := sfc_gamma(t+1.0);
          exit;
        end;
      end;
      if (p<1.0) and (t > z + 0.25) then begin
        {avoid scaling/rescaling with x^p in igammal}
        sfc_erfg := erfgsmallp(p,z,x);
      end
      else sfc_erfg := sfc_igammal(t,z)*t;
    end;
  end
end;


{---------------------------------------------------------------------------}
function sfc_erfi(x: extended): extended;
  {-Return the imaginary error function erfi(x) = erf(ix)/i}
const
  tbsph: THexExtW = ($688D,$14DB,$BA82,$906E,$3FFF);  {1.1283791670955125738}
var
  z: extended;
  TwoBySqrtPi: extended absolute tbsph; {=2/sqrt(Pi)}
begin
  {Ref: http://mathworld.wolfram.com/Erfi.html}
  if IsNanOrInf(x) then sfc_erfi := x
  else begin
    z := abs(x);
    if z<=1e-5 then sfc_erfi := x*TwoBySqrtPi*(1.0+sqr(z)/THREE)
    else if z >= 106.567 then sfc_erfi := copysign(PosInf_x,x)
    else sfc_erfi := TwoBySqrtPi*sfc_dawson(x)*expx2(abs(x));
  end;
end;


{---------------------------------------------------------------------------}
function inverror(p,q: extended): extended;
  {-Return value of inverse error function: erf_inv(p) if p <= 0.5, erfc_inv(q) otherwise}
var
  r,x,z: extended;
const
  Y1= 747689.0/8388608.0; {0.0891314744949340820313f}
  P1: array[0..7] of extended = (
         -0.000508781949658280665617,
         -0.00836874819741736770379,
          0.0334806625409744615033,
         -0.0126926147662974029034,
         -0.0365637971411762664006,
          0.0219878681111168899165,
          0.00822687874676915743155,
         -0.00538772965071242932965);
  Q1: array[0..9] of extended = (
          1.0,
         -0.970005043303290640362,
         -1.56574558234175846809,
          1.56221558398423026363,
          0.662328840472002992063,
         -0.71228902341542847553,
         -0.0527396382340099713954,
          0.0795283687341571680018,
         -0.00233393759374190016776,
          0.000886216390456424707504);
const
  Y2= 73711.0/32768.0; {2.249481201171875f}
  P2: array[0..8] of extended = (
         -0.202433508355938759655,
          0.105264680699391713268,
          8.37050328343119927838,
         17.6447298408374015486,
        -18.8510648058714251895,
        -44.6382324441786960818,
         17.445385985570866523,
         21.1294655448340526258,
         -3.67192254707729348546);
  Q2: array[0..8] of extended = (
          1.0,
          6.24264124854247537712,
          3.9713437953343869095,
        -28.6608180499800029974,
        -20.1432634680485188801,
         48.5609213108739935468,
         10.8268667355460159008,
        -22.6436933413139721736,
          1.72114765761200282724);
const
  Y3= 26451.0/32768.0; {0.807220458984375f}
  P3: array[0..10] of extended = (
         -0.131102781679951906451,
         -0.163794047193317060787,
          0.117030156341995252019,
          0.387079738972604337464,
          0.337785538912035898924,
          0.142869534408157156766,
          0.0290157910005329060432,
          0.00214558995388805277169,
         -0.679465575181126350155e-6,
          0.285225331782217055858e-7,
         -0.681149956853776992068e-9);
  Q3: array[0..7] of extended = (
          1.0,
          3.46625407242567245975,
          5.38168345707006855425,
          4.77846592945843778382,
          2.59301921623620271374,
          0.848854343457902036425,
          0.152264338295331783612,
          0.01105924229346489121);
const
  Y4= 985615.0/1048576.0; {0.93995571136474609375f}
  P4: array[0..8] of extended = (
         -0.0350353787183177984712,
         -0.00222426529213447927281,
          0.0185573306514231072324,
          0.00950804701325919603619,
          0.00187123492819559223345,
          0.000157544617424960554631,
          0.460469890584317994083e-5,
         -0.230404776911882601748e-9,
          0.266339227425782031962e-11);
  Q4: array[0..6] of extended = (
          1.0,
          1.3653349817554063097,
          0.762059164553623404043,
          0.220091105764131249824,
          0.0341589143670947727934,
          0.00263861676657015992959,
          0.764675292302794483503e-4);
const
  Y5= 1031409.0/1048576.0; {0.98362827301025390625f}
  P5: array[0..8] of extended = (
         -0.0167431005076633737133,
         -0.00112951438745580278863,
          0.00105628862152492910091,
          0.000209386317487588078668,
          0.149624783758342370182e-4,
          0.449696789927706453732e-6,
          0.462596163522878599135e-8,
         -0.281128735628831791805e-13,
          0.99055709973310326855e-16);
  Q5: array[0..6] of extended = (
          1.0,
          0.591429344886417493481,
          0.138151865749083321638,
          0.0160746087093676504695,
          0.000964011807005165528527,
          0.275335474764726041141e-4,
          0.282243172016108031869e-6);
const
  Y6= 1045583.0/1048576.0;  {0.99714565277099609375f}
  P6: array[0..7] of extended = (
         -0.0024978212791898131227,
         -0.779190719229053954292e-5,
          0.254723037413027451751e-4,
          0.162397777342510920873e-5,
          0.396341011304801168516e-7,
          0.411632831190944208473e-9,
          0.145596286718675035587e-11,
         -0.116765012397184275695e-17);
  Q6: array[0..6] of extended = (
          1.0,
          0.207123112214422517181,
          0.0169410838120975906478,
          0.000690538265622684595676,
          0.145007359818232637924e-4,
          0.144437756628144157666e-6,
          0.509761276599778486139e-9);
const
  Y7= 1047961.0/1048576.0;  {0.99941349029541015625f}
  P7: array[0..7] of extended = (
         -0.000539042911019078575891,
         -0.28398759004727721098e-6,
          0.899465114892291446442e-6,
          0.229345859265920864296e-7,
          0.225561444863500149219e-9,
          0.947846627503022684216e-12,
          0.135880130108924861008e-14,
         -0.348890393399948882918e-21);
  Q7: array[0..6] of extended = (
          1.0,
          0.0845746234001899436914,
          0.00282092984726264681981,
          0.468292921940894236786e-4,
          0.399968812193862100054e-6,
          0.161809290887904476097e-8,
          0.231558608310259605225e-11);

begin
  {Ref: Boost [19], erf_inv.hpp/erf_inv_imp}
  if p <= 0.5 then begin
    z := p * (p + 10);
    r := PolEvalX(p, P1, 8) / PolEvalX(p, Q1, 10);
    inverror := z*Y1 + z*r;
  end
  else if q >= 0.25 then begin
    z := q - 0.25;
    r := PolEvalX(z, P2, 9) / PolEvalX(z, Q2, 9);
    inverror := sqrt(-2.0*ln(q)) / (Y2 + r);
   end
   else begin
     x := sqrt(-ln(q));
     if x<3.0 then begin
       z := x - 1.125;
       r := PolEvalX(z, P3, 11) / PolEvalX(z, Q3, 8);
       inverror := Y3*x + r*x;
     end
     else if x<6.0 then begin
       z := x - 3.0;
       r := PolEvalX(z, P4, 9) / PolEvalX(z, Q4, 7);
       inverror := Y4*x + r*x;
     end
     else if x<18.0 then begin
       z := x - 6.0;
       r := PolEvalX(z, P5, 9) / PolEvalX(z, Q5, 7);
       inverror := Y5*x + r*x;
     end
     else if x<44.0 then begin
       z := x - 18.0;
       r := PolEvalX(z, P6, 8) / PolEvalX(z, Q6, 7);
       inverror := Y6*x + r*x;
     end
     else begin
       z := x - 44.0;
       r := PolEvalX(z, P7, 8) / PolEvalX(z, Q7, 7);
       inverror := Y7*x + r*x;
     end;
   end;
end;


{---------------------------------------------------------------------------}
function sfc_erfc_inv(x: extended): extended;
  {-Return the inverse function of erfc, erfc(erfc_inv(x)) = x, 0 < x < 2}
begin
  {Ref: Boost [19], erf_inv.hpp/erfc_inv}
  if IsNanOrInf(x) or (x<=0) or (x>=2.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_erfc_inv := NaN_x;
    exit;
  end;
  if x>1.0 then sfc_erfc_inv := -inverror(x-1.0, 2.0-x)
  else sfc_erfc_inv := inverror(1.0-x, x)
end;


{---------------------------------------------------------------------------}
function sfc_erf_inv(x: extended): extended;
  {-Return the inverse function of erf, erf(erf_inv(x)) = x, -1 < x < 1}
begin
  {Ref: Boost [19], erf_inv.hpp/erf_inv}
  if IsNanOrInf(x) or (abs(x) >=1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_erf_inv := NaN_x;
    exit;
  end;
  if x < 0.0 then sfc_erf_inv := -inverror(-x, 1.0+x)
  else sfc_erf_inv := inverror(x, 1.0-x)
end;


{---------------------------------------------------------------------------}
function sfc_erfce_inv(x: extended): extended;
  {-Return the functional inverse of erfce, erfce(erfce_inv(x)) = x, x > 0}
var
  z,d,f: extended;
const
  cp  = 1.125 + 0.003379167095512573896;  {2/sqrt(pi)}
  eps = 1e-14;          {tolerance in Newton loop}
  zm  = -106.53857691;  {early out, if starting value <= zm}
  x1  = 2;    {lower bound for using polynomial strarting value}
  x0  = 0.5;  {upper bound}
  c1  = -0.8862269254527580136; {polynomial coefficients from Maple with  }
  c2  = +0.6960409996039634807; {solve(series(exp(x^2)*erfc(x),x,4)=y, x);}
  c3  = -0.6293113124072449408;
begin
  if IsNanOrInf(x) or (x < MinExtended) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
     sfc_erfce_inv := NaN_x;
    exit;
  end;

  {Get starting values}
  if x > x1 then begin
    z := -sqrt(ln(0.5*(x+0.5)));
    if z <= zm then begin
      sfc_erfce_inv := z;
      exit;
    end;
  end
  else if x < x0 then begin
    z := 0.5*cp/x;    { = 1/(x*SqrtPi)}
  end
  else begin
    {x near 1: solve(series(exp(x^2)*erfc(x),x,4)=y, x);}
    d := x - 1.0;
    z := ((c3*d + c2)*d + c1)*d;
  end;

  {Newton loop}
  repeat
    f := sfc_erfce(z);
    d := f - x;
    if d<>0.0 then begin
      f := 2.0*z*f - cp;
      if f<>0.0 then d := d/f;
      z := z - d;
    end;
  until abs(d) <= eps*abs(z);
  sfc_erfce_inv := z;
end;


{---------------------------------------------------------------------------}
function sfc_erfi_inv(y: extended): extended;
  {-Return the inverse imaginary error function: erfi(erf_inv(y))=y}
var
  x,z,q,d: extended;
const
  eps = 1e-13; { < (eps_x/2)^(2/3}
begin
  if IsNanOrInf(y) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
     sfc_erfi_inv := NaN_x;
    exit;
  end;
  z := abs(y);
  if z>1.0 then begin
    {http://functions.wolfram.com/06.28.06.0017.01}
    {for large x: erfi(x) ~ exp(x^2)/x/Sqrt(Pi))  }
    q := ln(z);
    x := sqrt(0.5*lnPi + q);
    if x >= 1.0 then begin
      x := sqrt(ln(SqrtPi*x) + q);
      {if x >= 1.0 then x := sqrt(ln(SqrtPi*x) + q);}
    end;
  end
  else begin
    {https://en.wikipedia.org/wiki/Error_function#Inverse_functions}
    {erfi_inv = w - w^3/3 + 7/30*w^5 + ,,,}
    x := SqrtPi*0.5*z;
    x := x*(1.0 - sqr(x)/THREE);
  end;
  {Halley iteration if z is not small, otherwise the first}
  {two terms of the Maclaurin series are used for erfi_inv}
  if z > sqrt_epsh then begin
    repeat
      d := sfc_erfi(x)-z;
      q := expx2(abs(x))/SqrtPi;
      d := d/(2.0*q-x*d);
      x := x - d;
    until abs(d)<=eps*abs(x);
  end;
  {erfi_inv is odd}
  if y >= 0.0 then sfc_erfi_inv := x
  else sfc_erfi_inv := -x;
end;


{---------------------------------------------------------------------------}
function sfc_erf_z(x: extended): extended;
  {-Return the probability function erf_z = exp(-x^2/2)/sqrt(2*Pi)}
begin
  if IsNan(x) then sfc_erf_z := x
  else if abs(x) > 151.0 then sfc_erf_z := 0.0
  else sfc_erf_Z := expmx2h(x)/Sqrt_TwoPi;
end;


{---------------------------------------------------------------------------}
function sfc_erf_p(x: extended): extended;
  {-Return the probability function erf_p = integral(exp(-t^2/2)/sqrt(2*Pi), t=-Inf..x)}
const
  SQRTH = 0.7071067811865475244008;  {sqrt(0.5)=($6484,$F9DE,$F333,$B504,$3FFE)}
var
  y,z: extended;
begin
  if IsNan(x) then begin
    sfc_erf_p := Nan_x;
    exit;
  end;
  {erf_p = (1 + erf(z))/2 = erfc(-z)/2, where z = x/sqrt(2)}
  y := x*SQRTH;
  z := abs(y);
  if z<1.0 then sfc_erf_p := 0.5 + 0.5*sfc_erf(y)
  else if y < -x1e then sfc_erf_p := 0.0
  else if y >  x0e then sfc_erf_p := 1.0
  else begin
    {expmx2h is used to damp error amplification in computing exp(-x^2/2)}
    sfc_erfc_iusc(z,y,z);
    y := 0.5*y/z;
    y := y*expmx2h(x);
    if x<0.0 then sfc_erf_p := y
    else sfc_erf_p := 1.0 - y;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_erf_q(x: extended): extended;
  {-Return the probability function erf_q = integral(exp(-t^2/2)/sqrt(2*Pi), t=x..Inf)}
begin
  sfc_erf_q := sfc_erf_p(-x);  {HMF[1], 26.2.6}
end;


{---------------------------------------------------------------------------}
function sfc_erfh(x,h: extended): extended;
  {-Accurately compute erf(x+h) - erf(x-h)}
var
  a2,ah,ax,dn2,h2,h3,hf,hg,s,st,t,u,v,x2,xmh,xph,z: extended;
  j,n,n1: integer;
const
  c = 1.0 + 0.1283791670955125739;   {2/sqrt(pi)}
  p = 2.7695895202609194698622;      {ln(9*sqrt(pi))}
  nlneps = 43.668272375276554493286; {-ln(eps_x)}
begin
  {[70] A.R. Didonato, Significant Digit Computation of the Elliptical Coverage Function}
  if IsNanOrInf(x) or IsNanOrInf(h) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_erfh := NaN_x;
    exit;
  end;

  if h<=0.0 then begin
    if h=0.0 then sfc_erfh := 0
    else sfc_erfh := -sfc_erfh(x,-h);
    exit;
  end;

  {Note: since h >= 0 we have ah=h, but using ah}
  {is closer to Fortran and uses a local variable}
  ah  := abs(h);
  ax  := abs(x);
  xph := ax + ah;
  xmh := ax - ah;

  {*** Special cases 1..10 ***}
  {***************************}

  t := sqr(maxx(ah,ax));
  if (1.6*t*t < eps_x) then begin
    {case 1} {covered}
    sfc_erfh := 2.0*c*h*(0.5 + (0.5 - (x*x + h*h/THREE)));
    exit;
  end;
  if sqr(ax*ah) < 0.5*eps_x then begin
    {the value is 2.0*exp(-x*x)*erf(h)}
    t := 2.0*expx2(-ax);
    if (h*h >= 3.0*eps_x) then begin
      {case 2} {covered}
      sfc_erfh := t*sfc_erf(h);
    end
    else begin
      {case 3} {covered}
      sfc_erfh := c*h*t;
    end;
    exit;
  end;

  if ax <= ah then begin
    if xph < x0e then begin;
      {case 4} {covered}
     sfc_erfh := sfc_erf(xph) - sfc_erf(xmh);
    end
    else if (xmh > -x0e) then begin
      {case 5} {covered}
      sfc_erfh := sfc_erfc(xmh);
    end
    else begin
      {case 6} {covered}
      sfc_erfh := 2.0;
    end;
    exit;
  end;

  if (xmh >= 9.0) then begin
    if xmh*xmh + p > -ln_minext then begin
      {case 7} {covered}
      sfc_erfh := 0.0;
      exit;
    end;
  end;

  if (4.0*ah*ax > nlneps) then begin
    {case 8} {covered}
    sfc_erfh := sfc_erfc(xmh);
    exit;
  end;

  if (ax <= 3.0*ah) then begin
    if xph < 1.0 then begin
      {case 9} {covered}
      sfc_erfh := sfc_erf(xph) - sfc_erf(xmh);
    end
    else begin
      {case 10} {covered}
      sfc_erfh := sfc_erfc(xmh) - sfc_erfc(xph);
    end;
    exit;
  end;

  {***** Main cases *****}
  {**********************}
  if (ax <= 0.4) then begin
    h2 := xph*xph;
    a2 := xmh*xmh;
    x2 := ax + ax;
    st := 1.0;
    hf := xmh;
    n  := 0;
    n1 := 1;
    dn2:= 1.0;
    s  := 0.0;
    repeat
      inc(n);
      inc(n1,2);
      dn2:= -dn2/n;
      st := h2*st + x2*hf;
      hf := a2*hf;
      t  := st*dn2/n1;
      s  := s+t;
    until (abs(t) <= eps_x*abs(s));
    s := 0.5+(0.5+s);
    sfc_erfh := 2.0*c*ah*s;
  end
  else begin
    n := 1;
    h2 := 0.0;
    z  := expmx2h(x);
    u  := 2.0*ah*c*z;
    h3 := z;
    v  := 2.0*h*h;
    hf := 2.0*ax*ah;
    s  := 0.0;
    for j:=0 to 1 do begin
      repeat
        h2 := (hf*h3 - v*h2)/n;
        inc(n);
        h3 := (hf*h2 - v*h3)/n;
        inc(n);
        hg := h3/n;
        s  := s + hg;
      until (abs(hg) <= eps_x*abs(s));
    end;
    sfc_erfh := u*(s+z);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_erf2(x1,x2: extended):extended;
  {-Accurately compute erf(x2) - erf(x1)}
begin
  if IsNanOrInf(x1) or IsNanOrInf(x2) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_erf2 := NaN_x;
    exit;
  end;
  sfc_erf2 := sfc_erfh(0.5*(x2 + x1), 0.5*(x2 - x1));
end;


{---------------------------------------------------------------------------}
function ine_as(n: integer; x: extended): extended;
  {-Scaled repeated integrals of erfc, asymptotic expansion, x > 8.5 + n/3.5}
var
  s,t,z,u: extended;
  m: integer;
begin
  {Compute exp(x^2) * i^n erfc, using HMF[1], 7.2.14}
  t := 1.0;
  s := t;
  m := 0;
  u := n;
  z := 4.0*sqr(x);
  repeat
    inc(m);
    u := u + 2.0;
    t := -t*(u/m)*((u-1.0)/z);
    s := s + t;
  until abs(t) < eps_x*abs(s);
  u := 2.0*power(2.0*x,-n-1);
  ine_as := u*s/SqrtPi;
  {For unscaled version use:}
  {t := expx2(-x)/SqrtPi; ine_as := u*t*s;}
end;


{---------------------------------------------------------------------------}
function ine_fr(n: integer; x: extended): extended;
  {-Repeated integrals of erfc, forward recursion for x<0}
var
  f0,f1,fk: extended;
  k: integer;
begin
  {HMF[1], 7.2.5}
  f0 := 2.0*expx2(-abs(x))/SqrtPi;
  if n=-1 then begin
    ine_fr := f0;
    exit;
  end;
  f1 := sfc_erfc(x);
  for k:=1 to n do begin
    if f1<1e4000 then begin
      fk := 0.5*f0 - x*f1;
      f0 := f1;
      f1 := fk/k;
    end
    else begin
      fk := 0.5*f0/k - (x/k)*f1;
      f0 := f1;
      f1 := fk;
    end;
  end;
  ine_fr := f1;
end;


{---------------------------------------------------------------------------}
function ine_frp(n: integer; x: extended): extended;
  {-Scaled repeated integrals of erfc, forward recursion for x>0}
var
  f0,f1,fk: extended;
  k: integer;
begin
  {HMF[1], 7.2.5}
  f0 := 2.0/SqrtPi;
  if n=-1 then begin
    ine_frp := f0;
    exit;
  end;
  f1 := sfc_erfce(x);
  for k:=1 to n do begin
    fk := 0.5*f0 - x*f1;
    f0 := f1;
    f1 := fk/k;
  end;
  ine_frp := f1;
end;


{---------------------------------------------------------------------------}
function ine_cf(n: integer; x: extended; var ok: boolean): extended;
  {-Return the continued fraction i(n,x)/i(n-1,x), x>0}
var
  a,a0,c,d,f,t,tiny,tol: extended;
  k: integer;
const
  MAXIT = 1000;
begin
  tol  := 2.0*eps_x;
  tiny := Sqrt_MinExt;

  {Evaluate CF NIST[30], 7.18.13 using modified Lentz method}
  f := tiny;
  c := f;
  d := 0.0;
  a0:= 0.5*n;

  for k:=0 to MAXIT do begin
    if k=0 then a := 0.5
    else a := a0 + 0.5*k;
    d := x + a*d;
    c := x + a/c;
    if c=0.0 then c := tiny;
    if d=0.0 then d := tiny;
    d := 1.0/d;
    t := c*d;
    f := f*t;
    if abs(t-1.0) < tol then begin
      ine_cf := f;
      ok := true;
      exit;
    end;
  end;
  {No convergence}
  ine_cf := f;
  ok := false;
end;


{---------------------------------------------------------------------------}
procedure ine_br(n,m: integer; x,fm,fm1: extended; var fe,fn: extended);
  {-Scaled repeated integrals of erfc, backward recursion for x>0}
var
  f0,f1,fk,sc: extended;
  k: integer;
begin
  {Recursion starts at k=m > n, fm=inerfc(m), fm1=inerfc(m-1)}
  f1 := fm1;
  f0 := fm;
  fn := 1.0;
  sc := ldexp(1,1024);
  if n=m+1 then fn := fm1
  else if n=m then fn := fm
  else fn := 1;
  {HMF[1], 7.2.5}
  {i(k,x) = 2*[(k+2)*i(k+2,x) + x*i(k+1,x)]}
  for k:=m-1 downto -1 do begin
    fk := (k+2)*f1 + x*f0;
    f1 := f0;
    f0 := 2.0*fk;
    if k=n then fn := f0;
    if f0>sc then begin
      {rescale to prevent overflow}
      f0 := f0/sc;
      f1 := f1/sc;
      fn := fn/sc;
    end;
  end;
  fe := f0;
end;


{---------------------------------------------------------------------------}
function ine_xp(n: integer; x: extended): extended;
  {-Return the scaled repeated integral for x>0, ine_xp = exp(x^2) * i^nerfc(x)}
var
  f0,f1,fe,fn: extended;
  m: integer;
  cf_ok: boolean;
begin
  if x > 8.5 + n/3.5 then begin
    {asymptotic expansion if x is 'large'}
    ine_xp := ine_as(n,x);
    exit;
  end;
  if x>0.5 then begin
    {Try continued fraction first}
    f1 := ine_cf(n+1,x,cf_ok);
    if cf_ok then begin
      {Continued fraction OK, use backward recursion starting with fn}
      f0 := Sqrt_MinExt;
      f1 := f1*Sqrt_MinExt;
      ine_br(n,n,x,f0,f1,fe,fn);
      f1 := 2.0/SqrtPi;
      ine_xp := fn*(f1/fe);
      exit;
    end;
  end;
  fe := sqrt(n) + 23.0/x;
  if fe>250.0 then begin
    {Backward iteration count to large, uses forward iteration}
    ine_xp := ine_frp(n,x);
  end
  else begin
    if fe>180.0 then fe := 180.0;
    m := round(fe*fe+10);
    f0 := 0;
    f1 := MinExtended;
    ine_br(n,m,x,f0,f1,fe,fn);
    f1 := 2.0/SqrtPi;
    ine_xp := fn*(f1/fe);
  end
end;


{fn(x) = exp(x^2)*erfc(n,x)
   fn' = exp(x^2) * (2*x*erfc(n,x) - erfc(n-1,x)) < 0 because
         erfc(n,x)/erfc(n-1,x) = 1/(2x +a(n,x)) with a(n,x) > 0, (see NIST CF)
   fn' = exp(x^2) * ( -a(n,x)*erfc(n,x))) = - a(n,x)*fn;}

{---------------------------------------------------------------------------}
function sfc_inerfc(n: integer; x: extended): extended;
  {-Return the repeated integrals of erfc, n >= -1; scaled with exp(x^2) for x>0}
begin
  if IsNanOrInf(x) or (n < -1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_inerfc := NaN_x;
    exit;
  end;
  if n<1 then begin
    if x<0.0 then begin
      if n=0 then sfc_inerfc := sfc_erfc(x)
      else sfc_inerfc := 2.0*expx2(-abs(x))/SqrtPi;
    end
    else begin
      if n=0 then sfc_inerfc := sfc_erfce(x)
      else sfc_inerfc := 2.0/SqrtPi;
    end;
    exit;
  end;
  if abs(x) < 1e-38 then sfc_inerfc := ldexp(sfc_rgamma(0.5*n + 1.0),-n)
  else if x>0.0 then begin
    if n>2965 then sfc_inerfc := 0.0  {n>279 for double}
    else sfc_inerfc := ine_xp(n,x)
  end
  else begin
    sfc_inerfc := ine_fr(n,x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_dawson(x: extended): extended;
  {-Return Dawson's integral: dawson(x) = exp(-x^2)*integral(exp(t^2), t=0..x)}
var
  y: extended;
const
  xsml = 2.85158136717879458E-0010; {sqrt(0.75*eps_x)}
  xbig = 3.03700049997604969E+0009; {sqrt(1/eps_x)}
  xmax = 1.48716436919653959E+4931; {exp(min(-ln(2.0*minextended), ln(maxextended)))}
const
  n_daw  = 15;
  n_daw2 = 33;
  n_dawa = 37;
const
  daw: array[0..n_daw-1] of extended = (
         -0.6351734375145949201065127736293e-02,
         -0.2294071479677386939899824125866e+00,
         +0.2213050093908476441683979161786e-01,
         -0.1549265453892985046743057753375e-02,
         +0.8497327715684917456777542948066e-04,
         -0.3828266270972014924994099521309e-05,
         +0.1462854806250163197757148949539e-06,
         -0.4851982381825991798846715425114e-08,
         +0.1421463577759139790347568183304e-09,
         -0.3728836087920596525335493054088e-11,
         +0.8854942961778203370194565231369e-13,
         -0.1920757131350206355421648417493e-14,
         +0.3834325867246327588241074439253e-16,
         -0.7089154168175881633584099327999e-18,
         +0.1220552135889457674416901120000e-19);
const
  daw2: array[0..n_daw2-1] of extended = (
          -0.56886544105215527114160533733674e-01,
          -0.31811346996168131279322878048822e+00,
          +0.20873845413642236789741580198858e+00,
          -0.12475409913779131214073498314784e+00,
          +0.67869305186676777092847516423676e-01,
          -0.33659144895270939503068230966587e-01,
          +0.15260781271987971743682460381640e-01,
          -0.63483709625962148230586094788535e-02,
          +0.24326740920748520596865966109343e-02,
          -0.86219541491065032038526983549637e-03,
          +0.28376573336321625302857636538295e-03,
          -0.87057549874170423699396581464335e-04,
          +0.24986849985481658331800044137276e-04,
          -0.67319286764160294344603050339520e-05,
          +0.17078578785573543710504524047844e-05,
          -0.40917551226475381271896592490038e-06,
          +0.92828292216755773260751785312273e-07,
          -0.19991403610147617829845096332198e-07,
          +0.40963490644082195241210487868917e-08,
          -0.80032409540993168075706781753561e-09,
          +0.14938503128761465059143225550110e-09,
          -0.26687999885622329284924651063339e-10,
          +0.45712216985159458151405617724103e-11,
          -0.75187305222043565872243727326771e-12,
          +0.11893100052629681879029828987302e-12,
          -0.18116907933852346973490318263084e-13,
          +0.26611733684358969193001612199626e-14,
          -0.37738863052129419795444109905930e-15,
          +0.51727953789087172679680082229329e-16,
          -0.68603684084077500979419564670102e-17,
          +0.88123751354161071806469337321745e-18,
          -0.10974248249996606292106299624652e-18,
          +0.13261199326367178513595545891635e-19);
const
   dawa: array[0..n_dawa-1] of extended = (
          +0.1690485637765703755422637438849e-01,
          +0.8683252278406957990536107850768e-02,
          +0.2424864042417715453277703459889e-03,
          +0.1261182399572690001651949240377e-04,
          +0.1066453314636176955705691125906e-05,
          +0.1358159794790727611348424505728e-06,
          +0.2171042356577298398904312744743e-07,
          +0.2867010501805295270343676804813e-08,
          -0.1901336393035820112282492378024e-09,
          -0.3097780484395201125532065774268e-09,
          -0.1029414876057509247398132286413e-09,
          -0.6260356459459576150417587283121e-11,
          +0.8563132497446451216262303166276e-11,
          +0.3033045148075659292976266276257e-11,
          -0.2523618306809291372630886938826e-12,
          -0.4210604795440664513175461934510e-12,
          -0.4431140826646238312143429452036e-13,
          +0.4911210272841205205940037065117e-13,
          +0.1235856242283903407076477954739e-13,
          -0.5788733199016569246955765071069e-14,
          -0.2282723294807358620978183957030e-14,
          +0.7637149411014126476312362917590e-15,
          +0.3851546883566811728777594002095e-15,
          -0.1199932056928290592803237283045e-15,
          -0.6313439150094572347334270285250e-16,
          +0.2239559965972975375254912790237e-16,
          +0.9987925830076495995132891200749e-17,
          -0.4681068274322495334536246507252e-17,
          -0.1436303644349721337241628751534e-17,
          +0.1020822731410541112977908032130e-17,
          +0.1538908873136092072837389822372e-18,
          -0.2189157877645793888894790926056e-18,
          +0.2156879197938651750392359152517e-20,
          +0.4370219827442449851134792557395e-19,
          -0.8234581460977207241098927905177e-20,
          -0.7498648721256466222903202835420e-20,
          +0.3282536720735671610957612930039e-20);
begin
  {Ref: W. Fullerton [14] and [20], file ddaws.f}
  y := abs(x);
  if y<xsml then sfc_dawson := x
  else if y<1.0  then sfc_dawson := x*(0.75 + CSEvalX(2.0*y*y-1.0, daw, n_daw))
  else if y<4.0  then sfc_dawson := x*(0.25 + CSEvalX(0.125*y*y-1.0, daw2, n_daw2))
  else if y<xbig then sfc_dawson := (0.5 + CSEvalX(32.0/(y*y)-1.0, dawa, n_dawa))/x
  else if y<xmax then sfc_dawson := 0.5/x
  else sfc_dawson := 0.0;
end;


{---------------------------------------------------------------------------}
function gd_cf(p,x: extended): extended;
  {-Compute continued fraction for the generalized Dawson integral F(p,x); p,x>0}
var
  q,z: extended;
var
  a,b,c,d,f,t,tiny,tol: extended;
  k: integer;
const
  MAXIT = 500;
begin
  {Compute the CF from Dijkstra [48], (2.6) using modified Lentz method}
  tol  := 2.0*eps_x;
  tiny := Sqrt_MinExt;

  q := 1.0/p;
  z := power(x,p);
  k := 0;
  a := q*x;
  d := 0.0;
  if a > 0.5*MaxExtended*tiny then begin
    c := 16.0*(a/MaxExtended);
    f := c;
  end
  else begin
    f := tiny;
    c := f;
  end;

  repeat
    if k>0 then a := -k*z;
    b := q+z+k;
    d := b + a*d;
    c := b + a/c;
    if c=0.0 then c := tiny;
    if d=0.0 then d := tiny;
    d := 1.0/d;
    t := c*d;
    f := f*t;
    k := k+1;
  until (abs(t-1.0) < tol) or (k>MAXIT);
  gd_cf := f;
  if (k>MAXIT) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
end;


{---------------------------------------------------------------------------}
function sfc_gendaw(p,x: extended): extended;
  {-Return the generalized Dawson integral F(p,x) = exp(-x^p)*integral(exp(t^p), t=0..x); x,p >= 0}
var
  a,z,s,t: extended;
  k: integer;
const
  em1 = 0.36787944117144232159552; {exp(-1)}
begin
  if IsNanOrInf(p) or IsNanOrInf(x) or (p<0.0) or (x<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_gendaw := NaN_x;
    exit;
  end;
  {Dijkstra [48]: A Continued Fraction Expansion for a Generalization of Dawson's Integral}
  if x=0.0 then sfc_gendaw := 0.0
  else if p<=Sqrt_MinExt then sfc_gendaw := x
  else if p=1.0 then sfc_gendaw := -expm1(-x)
  else begin
    a := 1.0/p;
    if x=1.0 then begin
      {Avoid problems for very large/small p with x = x^p = 1}
      if a<1.0 then begin
        {sum M(a,a+1,1)/e with Taylor series; Dijkstra[48] (2.3)}
        s := 1.0;
        z := a;
        k := 0;
        {max k<=20, max fr a~1}
        repeat
          inc(k);
          z := z/k;
          t := z/(a+k);
          s := s + t;
        until t < eps_x*s;
        sfc_gendaw := s*em1;
        exit;
      end
      else if p<=1e-10 then begin
        {F(p,1) = M(1,1+a,-1) ~ 1 - p}
        sfc_gendaw := 1.0 - p;
        exit;
      end;
      {fall through for unproblematic p values}
    end;
    z := ln(x)*p;
    t := 45.0*(1.0+a);
    if z>t then begin
      {Asymptotic expansion for large x^p, Dijkstra[48] (2.5)}
      sfc_gendaw := a*power(x,1.0-p);
    end
    else if z<-t then begin
      {One term Taylor series for small x^p, Dijkstra[48] (2.4)}
      sfc_gendaw := x;
    end
    else begin
      {Continued fraction, Dijkstra (2.6)}
      sfc_gendaw := gd_cf(p,x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure fresnel_cfrac(x: extended; var fs,fc,ff,fg: extended);
  {-Continued fraction for Fresnel functions, x>=1.5}
  { fs=S(x), fc=C(x), ff=f(x), fg=g(x)}
const
  maxit2 = 2*200;
  eps = 1e-19;
type
  TComplex = record r,i: extended; end;
var
  a,p,q: extended;
  b,c,d,h,del: TComplex;
  n: longint;
begin
  {Ref: Numerical Recipes [13], ch. 6.9, p.256, function frenel}

  {Evaluate the continued fraction for the complex erfc function using}
  {the modified Lentz method. Complex arithmetic is done inline.}

  {b = complex(1.0,-pix2)}
  b.r  := 1.0;
  b.i  := -Pi*x*x;

  {c = complex(1.0/fpmin,0.0)}
  c.r := 1/Sqrt_MinExt;
  c.i := 0.0;

  {d = h = cdiv(one,b)}
  p   := b.r/b.i;
  q   := b.i+p*b.r;
  d.r := p/q;
  d.i := -1.0/q;
  h.r := d.r;
  h.i := d.i;
  n := -1;

  repeat
    inc(n,2);
    a := -n*(n+1);

    {b = cadd(b,complex(4.0,0.0))}
    b.r := b.r + 4.0;

    {d = cdiv(one,cadd(rcmul(a,d),b))}
    d.r := d.r*a + b.r;
    d.i := d.i*a + b.i;
    if abs(d.r) >= abs(d.i) then begin
      p   := d.i/d.r;
      q   := d.r + p*d.i;
      d.r := 1.0/q;
      d.i := -p/q;
    end
    else begin
      p   := d.r/d.i;
      q   := d.i + p*d.r;
      d.r := p/q;
      d.i := -1.0/q;
    end;

    {c = cadd(b,cdiv(complex(a,0.0),c))}
    if abs(c.r) >= abs(c.i) then begin
      p   := c.i/c.r;
      q   := c.r + p*c.i;
      c.r := b.r + a/q;
      c.i := b.i - p*a/q;
    end
    else begin
      p   := c.r/c.i;
      q   := c.i + p*c.r;
      c.r := b.r + a*p/q;
      c.i := b.i - a/q;
    end;

    {del = cmul(c,d)}
    del.r := c.r*d.r - c.i*d.i;
    del.i := c.i*d.r + c.r*d.i;

    {h = cmul(h,del)}
    q   := h.r*del.r - h.i*del.i;
    h.i := h.i*del.r + h.r*del.i;
    h.r := q;
  until (abs(del.r-1.0)+abs(del.i) < eps) or (n>maxit2);

  {No convergence of CF}
  if (n>maxit2) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));

  {fg + i*ff = x*h}
  fg  := x*h.r;
  ff  := x*h.i;

  {d = cmul(complex(x,-x),h)}
  d.r := x*(h.r + h.i);
  d.i := x*(h.i - h.r);

  {h = csub(1,cmul(complex(cos(0.5*pix2),sin(0.5*pix2)),d))}
  sincosPix2(x, p, q);
  h.r := 1.0 - (q*d.r - p*d.i);
  h.i :=     - (p*d.r + q*d.i);

  {(fc,fs) = cmul(complex(0.5,0.5), h)}
  fc := 0.5*(h.r-h.i);
  fs := 0.5*(h.r+h.i);
end;


{---------------------------------------------------------------------------}
procedure fresnel_series(x: extended; var fs,fc: extended);
  {-Power series for Fresnel functions, 0<=x<=1.5}
var
  fact,sign,sum,sumc,sums,term,test: extended;
  k,n: longint;
const
  maxit = 50;
  eps   = 1e-19;
begin
  {Ref: Numerical Recipes [13], ch. 6.9, function frenel}
  {Note: term for x=1.5, k=50 is of order (1.125*Pi)^50/50! ~ 1e-37}
  {and therefore maxit=50 is more than safe. NR parameters maxit=100}
  {together with eps=6e-8 are very very pessimistic :)}
  n := 3;
  sum  := 0.0;
  sums := 0.0;
  sumc := x;
  sign := 1.0;
  fact := Pi_2*x*x;
  term := x;
  for k:=1 to maxit do begin
    term := term*fact/k;
    sum  := sum + sign*term/n;
    test := abs(sum)*eps;
    if odd(k) then begin
      sign := -sign;
      sums := sum;
      sum  := sumc;
    end
    else begin
      sumc := sum;
      sum  := sums;
    end;
    if term < test then begin
      fs := sums;
      fc := sumc;
      exit;
    end;
    inc(n,2);
  end;
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
end;


{---------------------------------------------------------------------------}
procedure fg_large(x: extended; var f,g: extended);
  {-Return Fresnel auxiliary functions f,g for x >= 50}
var
  y,z: extended;
begin
  if x >= 65000.0 then begin
    f := 1.0;
    g := 1.0/Pi/x/x;
  end
  else begin
    {HMF 7.3.27/28}
    y := 1.0/(Pi*sqr(x));
    z := sqr(y);
    f := 1.0 -  3.0*z*(1.0 - 35.0*z);
    g := 1.0 - 15.0*z*(1.0 - 63.0*z);
    g := g*y;
  end;
  f := f/Pi/x;
  g := g/Pi/x;
end;


{---------------------------------------------------------------------------}
procedure sfc_fresnel(x: extended; var s,c: extended);
  {-Return the Fresnel integrals S(x)=integral(sin(Pi/2*t^2),t=0..x) and C(x)=integral(cos(Pi/2*t^2),t=0..x)}
var
  ax,cx,sx,f,g: extended;
const
  xsmall = 0.886e-1650;  {~ cbrt(6*succ(0)/Pi)}
  xlarge = 50;
begin
  if IsNan(x) then begin
    s := NaN_x;
    c := NaN_x;
    exit;
  end;
  ax := abs(x);
  if ax<xsmall then begin
    {C(x) = x + (-1/40*Pi^2)*x^5 + O(x^9) }
    {S(x) = (1/6*Pi)*x^3 + (-1/336*Pi^3)*x^7 + O(x^11) }
    s := 0.0;
    c := x;
  end
  else begin
    if ax<1.5 then fresnel_series(ax,s,c)
    else begin
      if ax >= xlarge then begin
        fg_large(ax,f,g);
        sincosPix2(ax,sx,cx);
        c := 0.5 + f*sx - g*cx;
        s := 0.5 - f*cx - g*sx;
      end
      else fresnel_cfrac(ax,s,c,f,g);
    end;
    if x<0.0 then begin
      s := -s;
      c := -c;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure sfc_fresnel_fg(x: extended; var f,g: extended);
  {-Return the Fresnel auxiliary functions f,g}
var
  s,c: extended;
const
  x0 = 1e-10;  {double 1e-8}
  x1 = 50.0;   {double 30.0}
begin
  if IsNan(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    f := NaN_x;
    g := NaN_x;
    exit;
  end;
  if abs(x) <= x0 then begin
    f := 0.5;
    g := 0.5 - x;
  end
  else if x <= 1.5 then begin
    {HMF 7.3.5/6]}
    if x<0.0 then sfc_fresnel(x,f,g) else fresnel_series(x,f,g);
    sincosPix2(x,s,c);
    f := 0.5 - f;
    g := 0.5 - g;
    x := f*c - g*s;
    g := g*c + f*s;
    f := x;
  end
  else if x < x1 then fresnel_cfrac(x,s,c,f,g)
  else fg_large(x,f,g);
end;


{---------------------------------------------------------------------------}
function sfc_gsi(x: extended): extended;
  {-Return the Goodwin-Staton integral = integral(exp(-t*t)/(t+x), t=0..Inf), x <> 0}
const
  gsl: array[0..28] of extended = (
         0.63106560560398446247e0,
         0.25051737793216708827e0,
        -0.28466205979018940757e0,
         0.8761587523948623552e-1,
         0.682602267221252724e-2,
        -0.1081129544192254677e-1,
         0.169101244117152176e-2,
         0.50272984622615186e-3,
        -0.18576687204100084e-3,
        -0.428703674168474e-5,
         0.1009598903202905e-4,
        -0.86529913517382e-6,
        -0.34983874320734e-6,
         0.6483278683494e-7,
         0.757592498583e-8,
        -0.277935424362e-8,
        -0.4830235135e-10,
         0.8663221283e-10,
        -0.394339687e-11,
        -0.209529625e-11,
         0.21501759e-12,
         0.3959015e-13,
        -0.692279e-14,
        -0.54829e-15,
         0.17108e-15,
         0.376e-17,
        -0.349e-17,
         0.7e-19,
         0.6e-19);
  gsh: array[0..23] of extended = (
         1.81775467984718758767e0,
        -0.9921146570744097467e-1,
        -0.894058645254819243e-2,
        -0.94955331277726785e-3,
        -0.10971379966759665e-3,
        -0.1346694539578590e-4,
        -0.172749274308265e-5,
        -0.22931380199498e-6,
        -0.3127844178918e-7,
        -0.436197973671e-8,
        -0.61958464743e-9,
        -0.8937991276e-10,
        -0.1306511094e-10,
        -0.193166876e-11,
        -0.28844270e-12,
        -0.4344796e-13,
        -0.659518e-14,
        -0.100801e-14,
        -0.15502e-15,
        -0.2397e-16,
        -0.373e-17,
        -0.58e-18,
        -0.9e-19,
        -0.1e-19);
var
  t: extended;
begin
  {For x>0: Pascal translation of function GOODST from MacLeod's MISCFUN [22]}
  if IsNan(x) or (x = 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_gsi := NaN_x;
    exit;
  end;
  if x<0.0 then begin
    t := SqrtPi*sfc_dawson(x);
    x := 0.5*sfc_eisx2(x);
    sfc_gsi := t - x;
  end
  else if x<=2.0 then begin
    if x<=eps_x then sfc_gsi := (-0.5)*EulerGamma - ln(x)
    else begin
      t := CSEvalX(x-1.0, gsl, 29);
      sfc_gsi := t - exp(-sqr(x))*ln(x);
    end;
  end
  else begin
    t := 0.5*SqrtPi/x;
    if x<0.2e20 then begin
      x := (6.0-x)/(x+2.0);
      sfc_gsi := t*CSEvalX(x, gsh, 24);
    end
    else sfc_gsi := t;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_expint3(x: extended): extended;
  {-Return the integral(exp(-t^3), t=0..x), x >= 0}
const
  x0 = 0.4e-6;
  x1 = 3.6;
  e3inf = 0.875 + 0.1797951156924921122e-1;  {2/9*Pi*3^(1/2)/Gamma(2/3)}
begin
  if IsNan(x) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_expint3 := NaN_x;
  end
  else if x <= x0 then sfc_expint3 := x
  else if x >= x1 then sfc_expint3 := e3inf
  else sfc_expint3 := sfc_igammal(one_x/THREE,x*x*x)/THREE;
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function znorm1(x: extended): extended;
begin
  znorm1 := 0.5*sfc_erf(x/sqrt2);
end;


{---------------------------------------------------------------------------}
function znorm2(x: extended): extended;
begin
  znorm2 := sfc_erf_p(-x);   {= sfc_erf_q(x), avoid one call-level}
end;


{---------------------------------------------------------------------------}
function tfrr(h,a,ah: extended): extended;
  {-Return Owen's T function T(h,a) for h >= 0, 0 <= a <= 1}
  { input ah must equal a*h, no argcheck for a,h}
var
  aj,tf,dj,gj,hh,dhs,aa,dv,vi: extended;  {use aa and hh, 'as' is keyword in Delphi}
  z,y,r,eps: extended;
  i,j,jj,m,mc,maxii,ii,ia,ih: integer;
const
  arange: array[0..6] of double = (
            0.025, 0.09, 0.15, 0.36, 0.5, 0.9, 0.99999);
  hrange: array[0..13] of double = (
            0.02, 0.06, 0.09, 0.125, 0.26, 0.4, 0.6,
            1.6, 1.7, 2.33, 2.4, 3.36, 3.4, 4.8);
  select: array[0..15*8-1] of byte = (
            1, 1, 2,13,13,13,13,13,13,13,13,16,16,16, 9,
            1, 2, 2, 3, 3, 5, 5,14,14,15,15,16,16,16, 9,
            2, 2, 3, 3, 3, 5, 5,15,15,15,15,16,16,16,10,
            2, 2, 3, 5, 5, 5, 5, 7, 7,16,16,16,16,16,10,
            2, 3, 3, 5, 5, 6, 6, 8, 8,17,17,17,12,12,11,
            2, 3, 5, 5, 5, 6, 6, 8, 8,17,17,17,12,12,12,
            2, 3, 4, 4, 6, 6, 8, 8,17,17,17,17,17,12,12,
            2, 3, 4, 4, 6, 6,18,18,18,18,17,17,17,12,12);
  meth: array[1..18] of byte = (
            1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 4, 4, 4, 4, 5, 6);
  order: array[1..18] of byte = (
            2, 3, 4, 5, 7,10,12,18,10,20,30,20, 4, 7, 8,20,13, 0);

const
  (*{orginal, accurate for double}
  nc2 = 20;
  c2: array[0..nc2] of extended = (
        0.99999999999999987510,     -0.99999999999988796462,
        0.99999999998290743652,     -0.99999999896282500134,
        0.99999996660459362918,     -0.99999933986272476760,
        0.99999125611136965852,     -0.99991777624463387686,
        0.99942835555870132569,     -0.99697311720723000295,
        0.98751448037275303682,     -0.95915857980572882813,
        0.89246305511006708555,     -0.76893425990463999675,
        0.58893528468484693250,     -0.38380345160440256652,
        0.20317601701045299653,     -0.82813631607004984866E-01,
        0.24167984735759576523E-01, -0.44676566663971825242E-02,
        0.39141169402373836468E-03);
  *)

  nc2 = 26;
  c2h: array[0..nc2] of THexExtw = ( {computed from chebyshev(1/(1+x^2, x=-1..1, 3e-20)}
         ($0000,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000000}
         ($FFAC,$FFFF,$FFFF,$FFFF,$BFFE),  {-9.9999999999999999545E-1}
         ($AD29,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999885037E-1}
         ($9B0B,$FFDF,$FFFF,$FFFF,$BFFE),  {-9.9999999999988491210E-1}
         ($550C,$F942,$FFFF,$FFFF,$3FFE),  {+9.9999999999386919645E-1}
         ($6C78,$2292,$FFFF,$FFFF,$BFFE),  {-9.9999999979861237771E-1}
         ($43B9,$D8DD,$FFEC,$FFFF,$3FFE),  {+9.9999999554062407207E-1}
         ($687D,$9F8B,$FED0,$FFFF,$BFFE),  {-9.9999992936458926407E-1}
         ($1DF0,$C7F0,$F1FE,$FFFF,$3FFE),  {+9.9999916525115436962E-1}
         ($3157,$5AF3,$809D,$FFFF,$BFFE),  {-9.9999240724259832097E-1}
         ($590B,$B6F9,$6EBB,$FFFC,$3FFE),  {+9.9994556506041955339E-1}
         ($B523,$D65B,$767B,$FFEB,$BFFE),  {-9.9968662761089784242E-1}
         ($5D1E,$86C2,$A5AF,$FF9F,$3FFE),  {+9.9852977309525088763E-1}
         ($117D,$7B4B,$1F8F,$FE8B,$BFFE),  {-9.9431035283059212416E-1}
         ($D2DF,$DCCD,$B1C4,$FB4D,$3FFE),  {+9.8165427261556396742E-1}
         ($31EC,$F906,$9396,$F345,$BFFE),  {-9.5028040347917722033E-1}
         ($5F0E,$88D7,$73C1,$E2C5,$3FFE),  {+8.8582538104023232715E-1}
         ($73DE,$2195,$297F,$C6A9,$BFFE),  {-7.7601870874203432253E-1}
         ($D29A,$2CC6,$868A,$9F03,$3FFE),  {+6.2114754556094035345E-1}
         ($5572,$24C9,$0C95,$E1D3,$BFFD),  {-4.4106330223657389981E-1}
         ($88B5,$8C5A,$6E2B,$8A10,$3FFD),  {+2.6965660363764693488E-1}
         ($FD98,$363F,$8DB9,$8CE1,$BFFC),  {-1.3757916872924624006E-1}
         ($BF78,$33D5,$7681,$E789,$3FFA),  {+5.6527579220746415151E-2}
         ($F41B,$2A95,$47E9,$9260,$BFF9),  {-1.7868175936501579486E-2}
         ($A824,$9488,$9AE0,$8515,$3FF7),  {+4.0614134060528344294E-3}
         ($C8B7,$3AC1,$73F7,$9A6C,$BFF4),  {-5.8907945440967510190E-4}
         ($6523,$11DB,$2AE0,$AB62,$3FF0)); {+4.0861002618012121744E-5}

   pts: array[1..13] of extended = (
          0.35082039676451715489E-02, 0.31279042338030753740E-01,
          0.85266826283219451090E-01, 0.16245071730812277011,
          0.25851196049125434828,     0.36807553840697533536,
          0.48501092905604697475,     0.60277514152618576821,
          0.71477884217753226516,     0.81475510988760098605,
          0.89711029755948965867,     0.95723808085944261843,
          0.99178832974629703586);

   wts: array[1..13] of extended = (
          0.18831438115323502887E-01, 0.18567086243977649478E-01,
          0.18042093461223385584E-01, 0.17263829606398753364E-01,
          0.16243219975989856730E-01, 0.14994592034116704829E-01,
          0.13535474469662088392E-01, 0.11886351605820165233E-01,
          0.10070377242777431897E-01, 0.81130545742299586629E-02,
          0.60419009528470238773E-02, 0.38862217010742057883E-02,
          0.16793031084546090448E-02);

begin

  ah  := a*h;
  eps := 4*eps_x;

  {Internal error, invalid mc}
  tfrr := Nan_x;

  {Determine appropriate method from t1...t6}
  ih := 15;
  for i:=1 to 14 do begin
    if h <= hrange[i-1] then begin
      ih := i;
      break;
    end;
  end;

  ia := 8;
  for i:=1 to 7 do begin
    if a <= arange[i-1] then begin
      ia := i;
      break;
    end;
  end;

  ii := select[ih-1 + (ia-1)*15];
  m  := order[ii];
  mc := meth[ii];

  if mc=1 then begin
    { t1(h, a, m) ; m = 2, 3, 4, 5, 7, 10, 12 or 18      }
    {  jj = 2j - 1 ; gj = exp(-h*h/2) * (-h*h/2)**j / j  }
    {  aj = a**(2j-1) / (2*pi)                           }
    hh  := -0.5*h*h;
    dj  := expm1(hh);
    dhs := 1.0+dj;
    aa  := a*a;
    j   := 1;
    jj  := 1;
    aj  := a/TwoPi;
    tf  := arctan(a)/TwoPi;
    gj  := hh * dhs;
    repeat
      dv := dj*aj/jj;
      tf := tf + dv;
      if (m <= j) and (abs(dv) < eps*abs(tf)) then break;
      j  := j + 1;
      jj := jj + 2;
      aj := aj*aa;
      dj := gj - dj;
      gj := gj*hh/j;
    until false;
    tfrr := tf;
  end
  else if mc=2 then begin
    { t2(h, a, m) ; m = 10, 20 or 30                             }
    {  z = (-1)^(i-1) * zi ; ii = 2i - 1                         }
    {  vi = (-1)^(i-1) * a^(2i-1) * exp[-(a*h)^2/2] / sqrt(2*pi) }
    maxii := m + m + 1;
    ii := 1;
    tf := 0.0;
    hh := h*h;
    aa := -a*a;
    vi := a*expmx2h(ah)/Sqrt_TwoPi;
    z  := znorm1(ah)/h;
    y  := 1.0/hh;
    repeat
      tf := tf + z;
      if (maxii <= ii) and (abs(z) <= eps*abs(tf)) then break;
      z  := y*(vi - ii*z);
      vi := aa * vi;
      ii := ii + 2;
    until false;
    tfrr := tf*expmx2h(h)/Sqrt_TwoPi;
  end
  else if mc=3 then begin
    { t3(h, a, m)                                     }
    {  ii = 2i - 1                                    }
    {  vi = a**(2i-1) * exp[-(a*h)**2/2] / sqrt(2*pi) }
    ii := 1;
    tf := 0.0;
    hh := h*h;
    aa := a*a;
    vi := a/Sqrt_TwoPi*expmx2h(ah);
    z  := znorm1(ah)/h;
    y  := 1.0/hh;
    for i:= 0 to nc2 do begin
      tf := tf + z*extended(c2h[i]);
      z  := y*(ii*z - vi);
      vi := aa * vi;
      ii := ii + 2;
    end;
    tfrr := tf*expmx2h(h)/Sqrt_TwoPi;
  end
  else if mc=4 then begin
    { t4(h, a, m) ; m = 4, 7, 8 or 20;  ii = 2i + 1      }
    {  ai = a * exp[-h*h*(1+a*a)/2] * (-a*a)**i / (2*pi) }
    maxii := m + m + 1;
    ii := 1;
    hh := h*h;
    aa := -a*a;
    tf := 0.0;
    aj := a*exp(-0.5*hh*(1.0-aa))/TwoPi;
    y  := 1.0;
    repeat
      z  := aj*y;
      tf := tf + aj*y;
      if (maxii <= ii) and (abs(z) <= eps*abs(tf)) then break;
      ii := ii + 2;
      y  := (1.0-hh*y)/ii;
      aj := aj*aa;
    until false;
    tfrr := tf;
  end
  else if mc=5 then begin
    { t5(h, a, m); m = 13; 2m - point gaussian quadrature }
    tf := 0.0;
    aa := a*a;
    hh := -0.5*h*h;
    for i:=1 to 13 do begin
      r  := 1.0 + aa*pts[i];
      tf := tf + wts[i]*exp(hh*r) / r;
    end;
    tfrr := a * tf;
  end
  else if mc=6 then begin
    { t6(h, a);  approximation for a near 1, (a<=1) }
    z  := znorm2(h);
    tf := 0.5*z*(1.0 - z);
    y  := 1.0 - a;
    r  := arctan2(y, 1.0+a);
    if r <> 0.0 then tf := tf - r/TwoPi * exp(-0.5*y*h*h/r);
    tfrr := tf;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_owent(h,a: extended): extended;
  {-Return Owen's T function T(h,a)}
const
  cut = 0.67;
var
  absa, absh, ah, normah, normh, tf: extended;
begin

  if IsNanOrInf(h) or IsNanOrInf(a) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_owent := NaN_x;
    exit;
  end;

  absh := abs(h);
  absa := abs(a);
  ah   := absa * absh;

  if absh = 0.0 then tf := arctan(absa)/TwoPi
  else if absa <= 1.0 then begin
    if absa = 0.0 then tf := 0.0
    else if absa = 1.0 then begin
      {a=1 should be dispatched to t6, but there are h values}
      {which select other methods, e.g. h=1/8 -> 1, h=10 -> 3}
      normh := znorm2(absh);
      tf := 0.5*normh*(1.0-normh);
    end
    else tf := tfrr(absh,absa,ah);
  end
  else if absh <= cut then begin
    tf := 0.25 - znorm1(absh) *znorm1(ah) - tfrr(ah, 1.0/absa, absh);
  end
  else begin
    normh  := znorm2(absh);
    normah := znorm2(ah);
    tf := 0.5*(normh + normah) - normh*normah - tfrr(ah, 1.0/absa, absh);
  end;

  if a < 0.0 then tf := -tf;

  sfc_owent := tf;
end;


{---------------------------------------------------------------------------}
function sfc_mq(m: integer; a,b: extended): extended;
  {-Return the generalized Marcum Q function Q(m,a,b), a,b >= 0}
var
  ab,t,s,q,x,eps: extended;
  k,i: integer;
  inv,neg: boolean;
const
  MINV = 50;  {Use sum for 1 - Qm for m > MINV}
  IMAX = 500; {Use 500 terms in repeat loop}
begin

  {Ref: R.T. Short [71] and http://www.phaselockedsystems.com/marcum_1.0.1.zip}
  {D. Morales-Jimenez et al., Connections between the Generalized Marcum Q-Function}
  {and a class of Hypergeometric Functions, https://arxiv.org/pdf/1303.5464}

  if IsNanOrInf(a) or IsNanOrInf(b) or (a < 0.0) or (b < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_mq := NaN_x;
    exit;
  end;

  neg := m<1;
  if neg then begin
    {Q(m,a,b) = 1 - Q(1-m,b,a); see Morales-Jimenez, formula (2)}
    t := a;
    a := b;
    b := t;
    m := 1-m;
  end;

  {special case a=0 or b=0}
  if b=0.0 then begin
    if neg then sfc_mq := 0.0
    else sfc_mq := 1.0;
    exit;
  end;

  if a=0.0 then begin
    if neg then sfc_mq := sfc_igammap(m, 0.5*b*b)
    else sfc_mq := sfc_igammaq(m, 0.5*b*b);
    exit;
  end;

  ab  := a*b;
  inv := neg or (a >= b) or (m >= MINV);
  eps := 8*eps_x;

  {Bessel series summation. If a<b and m<MINV compute sum, otherwise 1 - sum}
  if inv then begin
    {Compute 1 - sum. Prepare variables for loop k=m..infinity}
    q := b/a;
    x := power(q,m);
    s := 0.0;
  end
  else begin
    q := a/b;
    x := q;
    {sum from k=-(m-1) .. m-1. k=0 term is computed separately}
    {otherwise use Ie(-k,z) = Ie(k,z) to combine two terms.   }
    s := sfc_i0e(ab);
    for k:=1 to m-1 do begin
      {term = x^k*Ie(k,ab) + x^(-k)*Ie(-k,ab)}
      t := (x + 1.0/x)*sfc_ive(k,ab);
      s := s + t;
      x := x*q;
    end;
  end;

  {sum k=m..infinity}
  k := m;
  i := 1;
  repeat
    t := x*sfc_ive(k,ab);
    s := s + t;
    x := x*q;
    k := k+1;
    i := i+1;
  until (abs(t) <= abs(s)*eps) or (i>IMAX);
  {Test no convergence}
  if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));

  {Remove exponential scaling from sum}
  s := expmx2h(a-b)*s;

  if inv then begin
    if neg then sfc_mq := s
    else sfc_mq := 1.0 - s
  end
  else sfc_mq := s;
end;



end.
