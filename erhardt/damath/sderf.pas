unit sdErf;

{Double precision special functions: Error function and related}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Double precision special functions: Error function and related

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in damath_info.txt/references

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

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  11.02.13  W.Ehrhardt  Initial BP7 version from AMath.sferf
 1.00.01  11.02.13  we          erf_small, sfd_erfc_iusc, sfd_erfi
 1.00.02  11.02.13  we          sfd_fresnel
 1.00.03  01.03.13  we          Chebyshev degrees reduced in sfd_dawson, sfd_gsi, sfd_expint3

 1.05.00  15.08.13  we          sfd_erf_z, sfd_erf_p, sfd_erf_q

 1.20.00  25.08.15  we          Use THREE & SIXX to avoid 'optimization'

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


function sfd_dawson(x: double): double;
  {-Return Dawson's integral: dawson(x) = exp(-x^2)*integral(exp(t^2), t=0..x)}

function sfd_gendaw(p,x: double): double;
  {-Return the generalized Dawson integral F(p,x) = exp(-x^p)*integral(exp(t^p), t=0..x); x,p >= 0}

function sfd_erf(x: double): double;
  {-Return the error function erf(x) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..x)}

function sfd_erfc(x: double): double;
  {-Return the complementary error function erfc(x) = 1-erf(x)}

function sfd_erfce(x: double): double;
  {-Return the exponentially scaled complementary error function erfce(x) = exp(x^2)*erfc(x)}

function sfd_erfg(p,x: double): double;
  {-Return the generalized error function integral(exp(-t^p), t=0..x); x,p >= 0}

function sfd_erfi(x: double): double;
  {-Return the imaginary error function erfi(x) = erf(ix)/i}

function sfd_erf_inv(x: double): double;
  {-Return the inverse function of erf, erf(erf_inv(x)) = x, -1 < x < 1}

function sfd_erfc_inv(x: double): double;
  {-Return the inverse function of erfc, erfc(erfc_inv(x)) = x, 0 < x < 2}

function sfd_erf_p(x: double): double;
  {-Return the probability function erf_p = integral(exp(-t^2/2)/sqrt(2*Pi), t=-Inf..x)}

function sfd_erf_q(x: double): double;
  {-Return the probability function erf_q = integral(exp(-t^2/2)/sqrt(2*Pi), t=x..Inf)}

function sfd_erf_z(x: double): double;
  {-Return the probability function erf_z = exp(-x^2/2)/sqrt(2*Pi)}

function sfd_inerfc(n: integer; x: double): double;
  {-Return the repeated integrals of erfc, n >= -1; scaled with exp(x^2) for x>0}

function sfd_expint3(x: double): double;
  {-Return the integral(exp(-t^3), t=0..x), x >= 0}

procedure sfd_fresnel(x: double; var s,c: double);
  {-Return the Fresnel integrals S(x)=integral(sin(Pi/2*t^2),t=0..x) and C(x)=integral(cos(Pi/2*t^2),t=0..x)}

function sfd_gsi(x: double): double;
  {-Return the Goodwin-Staton integral = integral(exp(-t*t)/(t+x), t=0..Inf), x > 0}


implementation

uses
  DAMath,
  sdBasic,  {Basic common code}
  sdGamma,  {Gamma function and related}
  sdGamma2; {Gamma function and related}

{---------------------------------------------------------------------------}
function erf_small(x: double): double;
  {-Return erf(x) for |x| <= 1}
const
  {erf(x) = x*T(x^2)/U(x^2),  0 <= x <= 1, Peak relative error 7.6e-23}
  THex: array[0..6] of THexDblW = (
          ($7579,$654D,$DA01,$4130),  {1.104385395713178565288E6 }
          ($8B4F,$97BC,$EB74,$4102),  {1.549905740900882313773E5 }
          ($30F7,$C195,$9078,$40E7),  {4.825977363071025440855E4 }
          ($0B2C,$126C,$EAF2,$40A4),  {2.677472796799053019985E3 }
          ($D273,$81C9,$F2EA,$4071),  {2.871822526820825849235E2 }
          ($66E6,$E2D8,$9CA6,$4015),  {5.402980370004774841217E0 }
          ($4360,$0B67,$188E,$3FBC)); {1.097496774521124996496E-1}

  UHex: array[0..6] of THexDblW = (
          ($95AE,$25C3,$DE60,$412D),  {9.787360737578177599571E5 }
          ($5DF3,$B623,$4BC8,$411C),  {4.636021778692893773576E5 }
          ($8B69,$00E8,$4659,$40F8),  {9.942956272177178491525E4 }
          ($95A0,$7093,$5587,$40C8),  {1.245905812306219011252E4 }
          ($25EE,$3956,$5C44,$408E),  {9.715333124857259246107E2 }
          ($F1C7,$D103,$A0FE,$4046),  {4.525777638142203713736E1 }
          ($0000,$0000,$0000,$3FF0)); {1.000000000000000000000E0 }
var
  z: double;
  T: array[0..6] of double absolute THex;
  U: array[0..6] of double absolute UHex;
begin
  {Ref: Cephes [7], file ldouble\ndtrl.c}
  if x=0.0 then erf_small := x
  else begin
    z := x*x;
    erf_small := x*PolEval(z,T,7)/PolEval(z,U,7);
  end;
end;


{---------------------------------------------------------------------------}
procedure sfd_erfc_iusc(x: double; var p,q: double);
  {-(internal) unevaluated scaled erfc, p/q := exp(x^2)*erfc(x), 1<=x<=128}
const
  {erfc(x) = exp(-x^2)*P(1/x)/Q(1/x), 1/8<=1/x<=1, Peak relative error 5.8e-21}
  PHex: array[0..9] of THexDblW = (
          ($3CC6,$809B,$0DFA,$BEE3),  {-9.085943037416544232472E-6}
          ($13C0,$CD09,$8637,$4113),  { 3.198859502299390825278E5 }
          ($B447,$E047,$F015,$4152),  { 4.964439504376477951135E6 }
          ($A9AB,$F46A,$4726,$4182),  { 3.833161455208142870198E7 }
          ($5DCA,$C639,$4B13,$41A6),  { 1.870095071120436715930E8 }
          ($AD56,$3046,$94C9,$41C2),  { 6.234814405521647580919E8 }
          ($2577,$82F9,$962A,$41D5),  { 1.448651275892911637208E9 }
          ($FABA,$9299,$1A70,$41E1),  { 2.295563412811856278515E9 }
          ($087C,$065B,$1028,$41E1),  { 2.290171954844785638925E9 }
          ($5B09,$4073,$D8EF,$41D0)); { 1.130609921802431462353E9 }
  QHex: array[0..10] of THexDblW = (
          ($DEC5,$2A72,$4D8E,$4121),  { 5.669830829076399819566E5 }
          ($7627,$FF46,$C880,$4160),  { 8.799239977351261077610E6 }
          ($0B59,$EF06,$4417,$4190),  { 6.822450775590265689648E7 }
          ($6DF5,$1BA8,$04E6,$41B4),  { 3.358653716579278063988E8 }
          ($82B8,$87BC,$F818,$41D0),  { 1.138778654945478547049E9 }
          ($FB49,$5E39,$552A,$41E4),  { 2.729005809811924550999E9 }
          ($B2A0,$A0CE,$1779,$41F1),  { 4.588018188918609726890E9 }
          ($B5B2,$D691,$544D,$41F3),  { 5.188672873106859049556E9 }
          ($251F,$1F12,$9178,$41EA),  { 3.565928696567031388910E9 }
          ($05A2,$3DA6,$D8EF,$41D0),  { 1.130609910594093747762E9 }
          ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0 }
const
  {erfc(x) = exp(-x^2)*1/x*R(1/x^2)/S(1/x^2), 1/128<=1/x<1/8, Peak relative error 1.9e-21}
  RHex: array[0..4] of THexDblW = (
          ($E224,$8B78,$9F6D,$3F9B),  { 2.697535671015506686136E-2}
          ($A4EA,$86B0,$B846,$3FE1),  { 5.537445669807799246891E-1}
          ($6001,$EBE9,$8F6A,$400B),  { 3.445028155383625172464E0 }
          ($27C9,$EDAC,$B1DB,$401C),  { 7.173690522797138522298E0 }
          ($72A5,$F8F5,$F885,$400C)); { 3.621349282255624026891E0 }
  SHex: array[0..5] of THexDblW = (
          ($F65F,$CAE5,$7AE3,$3FA8),  { 4.781257488046430019872E-2}
          ($F307,$F266,$1616,$3FF0),  { 1.005392977603322982436E0 }
          ($EED7,$0411,$4ABE,$401A),  { 6.572990478128949439509E0 }
          ($600B,$E3DA,$AC9C,$402E),  { 1.533713447609627196926E1 }
          ($FAEC,$9AC2,$752A,$4025),  { 1.072884067182663823072E1 }
          ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0 }
var
  PP: array[0..9]  of double absolute PHex;
  QQ: array[0..10] of double absolute QHex;
  RR: array[0..4]  of double absolute RHex;
  SS: array[0..5]  of double absolute SHex;
var
  y: double;
begin
  {Ref: Cephes [7], file ldouble\ndtrl.c}
  y := 1.0/x;
  if x<8.0 then begin
    p := PolEval(y,PP,10);
    q := PolEval(y,QQ,11);
  end
  else begin
    q := y*y;
    p := y*PolEval(q,RR,5);
    q := PolEval(q,SS,6);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_erf(x: double): double;
  {-Return the error function erf(x) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..x)}
begin
  if IsInfD(x) then sfd_erf := copysignd(1.0,x)
  else if abs(x)<=1 then sfd_erf := erf_small(x)
  else if x > +7.0 then sfd_erf := +1.0
  else if x < -7.0 then sfd_erf := -1.0
  else sfd_erf := 1.0 - sfd_erfc(x);
end;


{---------------------------------------------------------------------------}
function sfd_erfc(x: double): double;
  {-Return the complementary error function erfc(x) = 1-erf(x)}
var
  p,q,a,z: double;
begin
  {expx2 is used to suppress error amplification in computing exp(-x^2)}
  {Ref: Cephes [7], file ldouble\ndtrl.c}
  if IsInfD(x) then begin
    sfd_erfc := 1.0-copysignd(1.0,x);
    exit;
  end;
  a := abs(x);
  if a < 1.0 then sfd_erfc := 1.0-erf_small(x)
  else if x < -7.0 then sfd_erfc := 2.0
  else if x > 106.8 then sfd_erfc := 0.0
  else begin
    {Here x is in one of the ranges -7..-1 or 1..107}
    {Compute z = expl(-a^2)}
    z := expx2(-a);
    {Compute p/q := exp(a^2)*erfc(a)}
    sfd_erfc_iusc(a,p,q);
    z := (z*p)/q;
    if x>=0 then sfd_erfc := z
    else sfd_erfc := 2.0-z;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_erfce(x: double): double;
  {-Return the exponentially scaled complementary error function erfce(x) = exp(x^2)*erfc(x)}
var
  y,z: double;
const
  xmin = -26.6287357137514895486518534833;
begin
  if IsNanD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
     sfd_erfce := NaN_d;
    exit;
  end;
  y := abs(x);
  if y <= 1e-9 then sfd_erfce := 1.0 - 2.0*x/SqrtPi
  else if y <= 1.0 then begin
    z := 1.0-erf_small(x);
    sfd_erfce := exp(x*x)*z;
  end
  else if x <= -6.75 then begin
    {here erfc(x)=2.0 accurate to double precision}
    if x < xmin then sfd_erfce := PosInf_d
    else sfd_erfce := 2.0*expx2(y);
  end
  else if y <= 128.0 then begin
    { -6.75 <= x <= 128: Recycle sfd_erfc_iusc}
    sfd_erfc_iusc(y,y,z);
    y := y/z;
    if x>0.0 then sfd_erfce := y
    else begin
      z := 2.0*expx2(-x);
      sfd_erfce := z-y;
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
    sfd_erfce := z/x/sqrtpi;
  end;
end;


{---------------------------------------------------------------------------}
function erfgsmallp(p,z,x: double): double;
  {-Return erfg(p,x), p small, z=x^p}
var
  a,t,s: double;
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
  until t < s*eps_d;
  {Here igammal(a,z) = s*z^a*exp(-z)/a and therefore}
  {erfg(p,x) = igammal(a,z)*a = s*x*exp(-z)}
  erfgsmallp := s*x*exp(-z);
end;


{---------------------------------------------------------------------------}
function sfd_erfg(p,x: double): double;
  {-Return the generalized error function integral(exp(-t^p), t=0..x); x,p >= 0}
var
  t,u,z: double;
const
  em1 = 0.36787944117144232159552;  {exp(-1)}
begin
  {Ref; NIST[30], 7.16}
  if x=0.0 then sfd_erfg := 0.0
  else if p<MinDouble then begin
    {t^p ~ 1}
    sfd_erfg := x*em1;
  end
  else if p=1.0 then sfd_erfg := -expm1(-x)
  else begin
    {sfd_erfg = igammal(t,x^p)*t}
    t := 1.0/p;
    u := p*ln(x);
    if u >= ln_MaxDbl then begin
      {igammal(t,e^u)*t ~ gamma(t)*t = gamma(t+1)}
      sfd_erfg := sfd_gamma(t+1.0)
    end
    else if u <= ln_MinDbl then begin
      {t^p ~ 0, exp(-t^p) ~ 1}
      sfd_erfg := x;
    end
    else begin
      {z = x^p = exp(u)}
      z := exp(u);
      if (p>1.0) and (z>ln_MaxDbl) then begin
        {if z^t*exp(-z)/gamma(t) underflows then igammal(t,z) ~ gamma(t)}
        if t*u - z - sfd_lngamma(t) < ln_MinDbl then begin
          sfd_erfg := sfd_gamma(t+1.0);
          exit;
        end;
      end;
      if (p<1.0) and (t > z + 0.25) then begin
        {avoid scaling/rescaling with x^p in igammal}
        sfd_erfg := erfgsmallp(p,z,x);
      end
      else sfd_erfg := sfd_igammal(t,z)*t;
    end;
  end
end;


{---------------------------------------------------------------------------}
function sfd_erfi(x: double): double;
  {-Return the imaginary error function erfi(x) = erf(ix)/i}
const
  tbsph: THexDblW = ($9B6D,$5042,$0DD7,$3FF2);  { 1.12837916709551E+0000}
var
  z: double;
  TwoBySqrtPi: double absolute tbsph; {=2/sqrt(Pi)}
begin
  {Ref: http://mathworld.wolfram.com/Erfi.html}
  if IsNanOrInfD(x) then sfd_erfi := x
  else begin
    z := abs(x);
    if z<=1e-5 then sfd_erfi := x*TwoBySqrtPi*(1.0+sqr(z)/THREE)
    else if z > 26.6417388 then sfd_erfi := copysignd(PosInf_d,x)
    else sfd_erfi := TwoBySqrtPi*sfd_dawson(x)*expx2(abs(x));
  end;
end;


{---------------------------------------------------------------------------}
function inverror(p,q: double): double;
  {-Return value of inverse error function: erf_inv(p) if p <= 0.5, erfc_inv(q) otherwise}
var
  r,x,z: double;
const
  Y1= 747689.0/8388608.0; {0.0891314744949340820313f}
  P1: array[0..7] of double = (
         -0.000508781949658280665617,
         -0.00836874819741736770379,
          0.0334806625409744615033,
         -0.0126926147662974029034,
         -0.0365637971411762664006,
          0.0219878681111168899165,
          0.00822687874676915743155,
         -0.00538772965071242932965);
  Q1: array[0..9] of double = (
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
  P2: array[0..8] of double = (
         -0.202433508355938759655,
          0.105264680699391713268,
          8.37050328343119927838,
         17.6447298408374015486,
        -18.8510648058714251895,
        -44.6382324441786960818,
         17.445385985570866523,
         21.1294655448340526258,
         -3.67192254707729348546);
  Q2: array[0..8] of double = (
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
  P3: array[0..10] of double = (
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
  Q3: array[0..7] of double = (
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
  P4: array[0..8] of double = (
         -0.0350353787183177984712,
         -0.00222426529213447927281,
          0.0185573306514231072324,
          0.00950804701325919603619,
          0.00187123492819559223345,
          0.000157544617424960554631,
          0.460469890584317994083e-5,
         -0.230404776911882601748e-9,
          0.266339227425782031962e-11);
  Q4: array[0..6] of double = (
          1.0,
          1.3653349817554063097,
          0.762059164553623404043,
          0.220091105764131249824,
          0.0341589143670947727934,
          0.00263861676657015992959,
          0.764675292302794483503e-4);
const
  Y5= 1031409.0/1048576.0; {0.98362827301025390625f}
  P5: array[0..8] of double = (
         -0.0167431005076633737133,
         -0.00112951438745580278863,
          0.00105628862152492910091,
          0.000209386317487588078668,
          0.149624783758342370182e-4,
          0.449696789927706453732e-6,
          0.462596163522878599135e-8,
         -0.281128735628831791805e-13,
          0.99055709973310326855e-16);
  Q5: array[0..6] of double = (
          1.0,
          0.591429344886417493481,
          0.138151865749083321638,
          0.0160746087093676504695,
          0.000964011807005165528527,
          0.275335474764726041141e-4,
          0.282243172016108031869e-6);
const
  Y6= 1045583.0/1048576.0;  {0.99714565277099609375f}
  P6: array[0..7] of double = (
         -0.0024978212791898131227,
         -0.779190719229053954292e-5,
          0.254723037413027451751e-4,
          0.162397777342510920873e-5,
          0.396341011304801168516e-7,
          0.411632831190944208473e-9,
          0.145596286718675035587e-11,
         -0.116765012397184275695e-17);
  Q6: array[0..6] of double = (
          1.0,
          0.207123112214422517181,
          0.0169410838120975906478,
          0.000690538265622684595676,
          0.145007359818232637924e-4,
          0.144437756628144157666e-6,
          0.509761276599778486139e-9);
const
  Y7= 1047961.0/1048576.0;  {0.99941349029541015625f}
  P7: array[0..7] of double = (
         -0.000539042911019078575891,
         -0.28398759004727721098e-6,
          0.899465114892291446442e-6,
          0.229345859265920864296e-7,
          0.225561444863500149219e-9,
          0.947846627503022684216e-12,
          0.135880130108924861008e-14,
         -0.348890393399948882918e-21);
  Q7: array[0..6] of double = (
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
    r := PolEval(p, P1, 8) / PolEval(p, Q1, 10);
    inverror := z*Y1 + z*r;
  end
  else if q >= 0.25 then begin
    z := q - 0.25;
    r := PolEval(z, P2, 9) / PolEval(z, Q2, 9);
    inverror := sqrt(-2.0*ln(q)) / (Y2 + r);
   end
   else begin
     x := sqrt(-ln(q));
     if x<3.0 then begin
       z := x - 1.125;
       r := PolEval(z, P3, 11) / PolEval(z, Q3, 8);
       inverror := Y3*x + r*x;
     end
     else if x<6.0 then begin
       z := x - 3.0;
       r := PolEval(z, P4, 9) / PolEval(z, Q4, 7);
       inverror := Y4*x + r*x;
     end
     else if x<18.0 then begin
       z := x - 6.0;
       r := PolEval(z, P5, 9) / PolEval(z, Q5, 7);
       inverror := Y5*x + r*x;
     end
     else if x<44.0 then begin
       z := x - 18.0;
       r := PolEval(z, P6, 8) / PolEval(z, Q6, 7);
       inverror := Y6*x + r*x;
     end
     else begin
       z := x - 44.0;
       r := PolEval(z, P7, 8) / PolEval(z, Q7, 7);
       inverror := Y7*x + r*x;
     end;
   end;
end;


{---------------------------------------------------------------------------}
function sfd_erfc_inv(x: double): double;
  {-Return the inverse function of erfc, erfc(erfc_inv(x)) = x, 0 < x < 2}
begin
  {Ref: Boost [19], erf_inv.hpp/erfc_inv}
  if IsNanOrInfD(x) or (x<=0) or (x>=2.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_erfc_inv := NaN_d;
    exit;
  end;
  if x>1.0 then sfd_erfc_inv := -inverror(x-1.0, 2.0-x)
  else sfd_erfc_inv := inverror(1.0-x, x)
end;


{---------------------------------------------------------------------------}
function sfd_erf_inv(x: double): double;
  {-Return the inverse function of erf, erf(erf_inv(x)) = x, -1 < x < 1}
begin
  {Ref: Boost [19], erf_inv.hpp/erf_inv}
  if IsNanOrInfD(x) or (abs(x) >=1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_erf_inv := NaN_d;
    exit;
  end;
  if x < 0.0 then sfd_erf_inv := -inverror(-x, 1.0+x)
  else sfd_erf_inv := inverror(x, 1.0-x)
end;


{---------------------------------------------------------------------------}
function sfd_erf_z(x: double): double;
  {-Return the probability function erf_z = exp(-x^2/2)/sqrt(2*Pi)}
begin
  if IsNanD(x) then sfd_erf_z := x
  else if abs(x) > 151.0 then sfd_erf_z := 0.0
  else sfd_erf_Z := expmx2h(x)/Sqrt_TwoPi;
end;


{---------------------------------------------------------------------------}
function sfd_erf_p(x: double): double;
  {-Return the probability function erf_p = integral(exp(-t^2/2)/sqrt(2*Pi), t=-Inf..x)}
const
  SQRTH = 0.7071067811865475244008;  {sqrt(0.5)=($6484,$F9DE,$F333,$B504,$3FFE)}
var
  y,z: double;
begin
  if IsNanD(x) then begin
    sfd_erf_p := x;
    exit;
  end;
  {erf_p = (1 + erf(z))/2 = erfc(z)/2, where z = x/sqrt(2)}
  y := x*SQRTH;
  z := abs(y);
  if z<1.0 then sfd_erf_p := 0.5 + 0.5*sfd_erf(y)
  else if x<-150.7 then sfd_erf_p := 0.0
  else if x>=9.5 then sfd_erf_p := 1.0
  else begin
    sfd_erfc_iusc(z,y,z);
    y := 0.5*y/z;
    y := y*expmx2h(x);
    if x<0.0 then sfd_erf_p := y
    else sfd_erf_p := 1.0 - y;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_erf_q(x: double): double;
  {-Return the probability function erf_q = integral(exp(-t^2/2)/sqrt(2*Pi), t=x..Inf)}
begin
  if IsNanD(x) then sfd_erf_q := x
  else sfd_erf_q := sfd_erf_p(-x);  {HMF[1], 26.2.6}
end;


{---------------------------------------------------------------------------}
function ine_as(n: integer; x: double): double;
  {-Scaled repeated integrals of erfc, asymptotic expansion, x > 8.5 + n/3.5}
var
  s,t,z,u: double;
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
  until abs(t) < eps_d*abs(s);
  u := 2.0*power(2.0*x,-n-1);
  ine_as := u*s/SqrtPi;
  {For unscaled version use:}
  {t := expx2(-x)/SqrtPi; ine_as := u*t*s;}
end;


{---------------------------------------------------------------------------}
function ine_fr(n: integer; x: double): double;
  {-Repeated integrals of erfc, forward recursion for x<0}
var
  f0,f1,fk: double;
  k: integer;
begin
  {HMF[1], 7.2.5}
  f0 := 2.0*expx2(-abs(x))/SqrtPi;
  if n=-1 then begin
    ine_fr := f0;
    exit;
  end;
  f1 := sfd_erfc(x);
  for k:=1 to n do begin
    if f1<1e250 then begin
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
function ine_frp(n: integer; x: double): double;
  {-Scaled repeated integrals of erfc, forward recursion for x>0}
var
  f0,f1,fk: double;
  k: integer;
begin
  {HMF[1], 7.2.5}
  f0 := 2.0/SqrtPi;
  if n=-1 then begin
    ine_frp := f0;
    exit;
  end;
  f1 := sfd_erfce(x);
  for k:=1 to n do begin
    fk := 0.5*f0 - x*f1;
    f0 := f1;
    f1 := fk/k;
  end;
  ine_frp := f1;
end;


{---------------------------------------------------------------------------}
function ine_cf(n: integer; x: double; var ok: boolean): double;
  {-Return the continued fraction i(n,x)/i(n-1,x), x>0}
var
  a,a0,c,d,f,t,tiny,tol: double;
  k: integer;
const
  MAXIT = 1000;
begin
  tol  := 2.0*eps_d;
  tiny := Sqrt_MinDbl;

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
procedure ine_br(n,m: integer; x,fm,fm1: double; var fe,fn: double);
  {-Scaled repeated integrals of erfc, backward recursion for x>0}
var
  f0,f1,fk,sc: double;
  k: integer;
begin
  {Recursion starts at k=m > n, fm=inerfc(m), fm1=inerfc(m-1)}
  f1 := fm1;
  f0 := fm;
  fn := 1.0;
  sc := ldexpd(1,256);
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
function ine_xp(n: integer; x: double): double;
  {-Return the scaled repeated integral for x>0, ine_xp = exp(x^2) * i^nerfc(x)}
var
  f0,f1,fe,fn: double;
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
      f0 := Sqrt_MinDbl;
      f1 := f1*Sqrt_MinDbl;
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
    f1 := MinDouble;
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
function sfd_inerfc(n: integer; x: double): double;
  {-Return the repeated integrals of erfc, n >= -1; scaled with exp(x^2) for x>0}
begin
  if IsNanOrInfD(x) or (n < -1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_inerfc := NaN_d;
    exit;
  end;
  if n<1 then begin
    if x<0.0 then begin
      if n=0 then sfd_inerfc := sfd_erfc(x)
      else sfd_inerfc := 2.0*expx2(-abs(x))/SqrtPi;
    end
    else begin
      if n=0 then sfd_inerfc := sfd_erfce(x)
      else sfd_inerfc := 2.0/SqrtPi;
    end;
    exit;
  end;
  if abs(x) < 1e-38 then sfd_inerfc := ldexpd(sfd_rgamma(0.5*n + 1.0),-n)
  else if x>0.0 then begin
    if n>279 then sfd_inerfc := 0.0
    else sfd_inerfc := ine_xp(n,x)
  end
  else begin
    sfd_inerfc := ine_fr(n,x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_dawson(x: double): double;
  {-Return Dawson's integral: dawson(x) = exp(-x^2)*integral(exp(t^2), t=0..x)}
var
  y: double;
const
  xsml = 1.29047841397589236E-8;  {sqrt(0.75*eps_d)}
  xbig = 6.7108864e7;             {sqrt(1/eps_d)}
  xmax = 2.24711641857795676E307; {exp(min(-ln(2.0*MinDouble), ln(MaxDouble)))}
const
  n_daw  = 15;
  n_daw2 = 33;
  n_dawa = 32;
const
  daw: array[0..n_daw-1] of double = (
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
  daw2: array[0..n_daw2-1] of double = (
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
   dawa: array[0..n_dawa-1] of double = (
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
          -0.2189157877645793888894790926056e-18);
        { +0.2156879197938651750392359152517e-20,
          +0.4370219827442449851134792557395e-19,
          -0.8234581460977207241098927905177e-20,
          -0.7498648721256466222903202835420e-20,
          +0.3282536720735671610957612930039e-20);}
begin
  {Ref: W. Fullerton [14] and [20], file ddaws.f}
  y := abs(x);
  if y<xsml then sfd_dawson := x
  else if y<1.0  then sfd_dawson := x*(0.75 + CSEvalD(2.0*y*y-1.0, daw, n_daw))
  else if y<4.0  then sfd_dawson := x*(0.25 + CSEvalD(0.125*y*y-1.0, daw2, n_daw2))
  else if y<xbig then sfd_dawson := (0.5 + CSEvalD(32.0/(y*y)-1.0, dawa, n_dawa))/x
  else if y<xmax then sfd_dawson := 0.5/x
  else sfd_dawson := 0.0;
end;


{---------------------------------------------------------------------------}
function gd_cf(p,x: double): double;
  {-Compute continued fraction for the generalized Dawson integral F(p,x); p,x>0}
var
  q,z: double;
var
  a,b,c,d,f,t,tiny,tol: double;
  k: integer;
const
  MAXIT = 500;
begin
  {Compute the CF from Dijkstra [48], (2.6) using modified Lentz method}
  tol  := 2.0*eps_d;
  tiny := Sqrt_MinDbl;

  q := 1.0/p;
  z := power(x,p);
  k := 0;
  a := q*x;
  d := 0.0;
  if a > 0.5*MaxDouble*tiny then begin
    c := 16.0*(a/MaxDouble);
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
function sfd_gendaw(p,x: double): double;
  {-Return the generalized Dawson integral F(p,x) = exp(-x^p)*integral(exp(t^p), t=0..x); x,p >= 0}
var
  a,z,s,t: double;
  k: integer;
const
  em1 = 0.36787944117144232159552; {exp(-1)}
begin
  if IsNanOrInfD(p) or IsNanOrInfD(x) or (p<0.0) or (x<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_gendaw := NaN_d;
    exit;
  end;
  {Dijkstra [48]: A Continued Fraction Expansion for a Generalization of Dawson's Integral}
  if x=0.0 then sfd_gendaw := 0.0
  else if p<=Sqrt_MinDbl then sfd_gendaw := x
  else if p=1.0 then sfd_gendaw := -expm1(-x)
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
        until t < eps_d*s;
        sfd_gendaw := s*em1;
        exit;
      end
      else if p<=1e-10 then begin
        {F(p,1) = M(1,1+a,-1) ~ 1 - p}
        sfd_gendaw := 1.0 - p;
        exit;
      end;
      {fall through for unproblematic p values}
    end;
    z := ln(x)*p;
    t := 45.0*(1.0+a);
    if z>t then begin
      {Asymptotic expansion for large x^p, Dijkstra[48] (2.5)}
      sfd_gendaw := a*power(x,1.0-p);
    end
    else if z<-t then begin
      {One term Taylor series for small x^p, Dijkstra[48] (2.4)}
      sfd_gendaw := x;
    end
    else begin
      {Continued fraction, Dijkstra (2.6)}
      sfd_gendaw := gd_cf(p,x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure fresnel_cfrac(x: double; var fs,fc: double);
  {-Continued fraction for Fresnel functions, x>=1.5}
const
  maxit2 = 2*200;
  eps = 1e-16;
type
  TComplex = record r,i: double; end;
var
  a,p,q: double;
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
  c.r := 1/Sqrt_MinDbl;
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

  {d = cmul(complex(ax,-ax),h)}
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
procedure fresnel_series(x: double; var fs,fc: double);
  {-Power series for Fresnel functions, 0<=x<=1.5}
var
  fact,sign,sum,sumc,sums,term,test: double;
  k,n: longint;
const
  maxit = 50;
  eps   = 1e-16;
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
procedure fresnel_large(x: double; var fs,fc: double);
  {-Asymptotic expansion for Fresnel functions, |x| >= 1.562e6}
var
  sx,cx,y: double;
begin
  (* Maple: Asymptotic expansion of S(x) and C(x) for large x:  *)
  (*                                                            *)
  (*                             2                2             *)
  (*                 cos(1/2 Pi x )   sin(1/2 Pi x )       1    *)
  (*   S(x) =  1/2 - -------------- - --------------  + O(----) *)
  (*                      Pi x              2  3            5   *)
  (*                                      Pi  x            x    *)
  (*                                                            *)
  (*                             2                2             *)
  (*                 sin(1/2 Pi x )   cos(1/2 Pi x )       1    *)
  (*   C(x) =  1/2 + -------------- - --------------  + O(----) *)
  (*                      Pi x              2  3            5   *)
  (*                                      Pi  x            x    *)
  (*                                                            *)
  x := abs(x);
  if x > 2.8671e15 then begin
    {All terms other than 1/2 are negligible}
    fs := 0.5;
    fc := 0.5;
  end
  else begin
    sincosPix2(x,sx,cx);
    y  := 1.0/(Pi*x);
    fs := 0.5 - cx*y;
    fc := 0.5 + sx*y;
  end;
end;


{---------------------------------------------------------------------------}
procedure sfd_fresnel(x: double; var s,c: double);
  {-Return the Fresnel integrals S(x)=integral(sin(Pi/2*t^2),t=0..x) and C(x)=integral(cos(Pi/2*t^2),t=0..x)}
var
  ax: double;
const
  xsmall = 2.1458772928E-108;  {~ cbrt(6*succ(0)/Pi)}
  xlarge = 1.2221e5;           {~ cbrt(4/Pi^2/eps_d)}
begin
  if IsNanD(x) then begin
    s := NaN_d;
    c := NaN_d;
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
      if ax >= xlarge then fresnel_large(x,s,c)
      else fresnel_cfrac(ax,s,c);
    end;
    if x<0.0 then begin
      s := -s;
      c := -c;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_gsi(x: double): double;
  {-Return the Goodwin-Staton integral = integral(exp(-t*t)/(t+x), t=0..Inf), x > 0}
const
  gsl: array[0..26] of double = (
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
        -0.349e-17);
      {  0.7e-19,
         0.6e-19);}

  gsh: array[0..22] of double = (
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
        -0.9e-19);
      { -0.1e-19);}
var
  t: double;
begin
  {Pascal translation of function GOODST from A.J. MacLeod's MISCFUN [22]}
  if IsNanD(x) or (x <= 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_gsi := NaN_d;
    exit;
  end;
  if x<=2.0 then begin
    if x<=eps_d then sfd_gsi := (-0.5)*EulerGamma - ln(x)
    else begin
      t := CSEvalD(x-1.0, gsl, 27);
      sfd_gsi := t - exp(-sqr(x))*ln(x);
    end;
  end
  else begin
    t := 0.5*SqrtPi/x;
    if x<0.2e20 then begin
      x := (6.0-x)/(x+2.0);
      sfd_gsi := t*CSEvalD(x, gsh, 23);
    end
    else sfd_gsi := t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_expint3(x: double): double;
  {-Return the integral(exp(-t^3), t=0..x), x >= 0}
const
  ae3s: array[0..22] of double = (
          1.26919841422112601434,
         -0.24884644638414098226,
          0.8052622071723104125e-1,
         -0.2577273325196832934e-1,
          0.759987887307377429e-2,
         -0.203069558194040510e-2,
          0.49083458669932917e-3,
         -0.10768223914202077e-3,
          0.2155172626428984e-4,
         -0.395670513738429e-5,
          0.66992409338956e-6,
         -0.10513218080703e-6,
          0.1536258019825e-7,
         -0.209909603636e-8,
          0.26921095381e-9,
         -0.3251952422e-10,
          0.371148157e-11,
         -0.40136518e-12,
          0.4123346e-13,
         -0.403375e-14,
          0.37658e-15,
         -0.3362e-16,
          0.288e-17);
       { -0.24e-18,
          0.2e-19);}

  ae3l: array[0..21] of double = (
          1.92704649550682737293,
         -0.3492935652048138054e-1,
          0.145033837189830093e-2,
         -0.8925336718327903e-4,
          0.705423921911838e-5,
         -0.66717274547611e-6,
          0.7242675899824e-7,
         -0.878258256056e-8,
          0.116722344278e-8,
         -0.16766312812e-9,
          0.2575501577e-10,
         -0.419578881e-11,
          0.72010412e-12,
         -0.12949055e-12,
          0.2428703e-13,
         -0.473311e-14,
          0.95531e-15,
         -0.19914e-15,
          0.4277e-16,
         -0.944e-17,
          0.214e-17,
         -0.50e-18);
      {   0.12e-18,
         -0.3e-19,
          0.1e-19);}
const
  e3infh: THexDblW = ($03ED,$C4C6,$9349,$3FEC);  { 8.92979511569249E-0001}
var
  funinf: double absolute e3infh; {2/9*Pi*3^(1/2)/Gamma(2/3)}
  t,y: double;
begin
  {Pascal translation of function EXP3 from A.J. MacLeod's MISCFUN [22]}
  if IsNanD(x) or (x < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_expint3 := NaN_d;
    exit;
  end
  else if x <= 6e-6 then begin
    {expint3(x) = x - 1/4*x^4 + 1/14*x^7 + O(x^10)}
    sfd_expint3 := x;
  end
  else if x <= 2.0 then begin
    t := x*sqr(0.5*x) - 1.0;
    sfd_expint3 := x*CSEvalD(t, ae3s, 23);
  end
  else if x > 3.6 then sfd_expint3 := funinf
  else begin
    y := x*sqr(x);
    t := CSEvalD(16.0/y - 1.0, ae3l, 22);
    y := exp(-y)/(3.0*sqr(x));
    sfd_expint3 := funinf - y*t;
  end;
end;


end.
