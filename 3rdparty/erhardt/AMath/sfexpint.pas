unit sfExpInt;

{Common code for special functions: Exponential integrals and related}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Common code for special functions: Exponential integrals and related

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in amath_info.txt/references

                  [1] [HMF]: M. Abramowitz, I.A. Stegun. Handbook of Mathematical Functions, New York, 1970
                      http://www.math.sfu.ca/~cbm/aands/
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
                 [30] [NIST]: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST Handbook
                      of Mathematical Functions, Cambridge, 2010. Online resource: NIST Digital
                      Library of Mathematical Functions, http://dlmf.nist.gov/

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  17.08.10  W.Ehrhardt  Common version number after split
 1.00.01  07.09.10  we          avoid infinite loop if x is NAN in sfc_ei/sfc_e1

 1.05.00  31.03.11  we          uses sfBasic for RTE_ArgumentRange ifopt R+

 1.07.00  01.06.11  we          sfc_ein

 1.08.00  17.07.11  we          range sfc_li(x) now x>=0, x<>1

 1.13.00  20.06.12  we          corrected bound in sfc_en (11140.0 -> 11390.2)
 1.13.01  22.06.12  we          uniform asymptotic expansion ep_largep for n>=10000
 1.13.02  23.06.12  we          sfc_gei: generalized exponential integral E_p(x)

 1.14.00  23.07.12  we          more compact sfc_eiex

 1.18.00  19.05.13  we          Prevent some wrong compiler optimizations for div by 3

 1.24.00  11.03.14  we          Adjust arrays sizes in sfc_ci, si_medium

 1.32.00  25.01.15  we          sfc_ali from sfMisc
 1.32.01  25.01.15  we          sfc_ei_inv
 1.32.02  01.05.15  we          fix sfc_ei_inv for x Nan

 1.42.00  16.08.17  we          sfc_eisx2

 1.43.00  12.12.17  we          separate functions e1_small, e1ex_large
 1.43.01  13.12.17  we          sfc_e1s(x) = exp(x)*E1(x)
 1.43.02  14.12.17  we          function ei_asymp, sfc_eisx2 with ei_asymp
 1.43.03  14.12.17  we          sfc_eis(x) = exp(-x)*Ei(x)

 1.50.00  23.08.18  we          NanOrInf check in sfc_en
 1.50.01  07.09.18  we          sfc_en for n<0 and x>0
 1.50.02  07.09.18  we          sfc_eibeta

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

function sfc_ali(x: extended): extended;
  {-Return the functional inverse of li(x), li(sfc_ali(x))=x}

function sfc_chi(x: extended): extended;
  {-Return the hyperbolic cosine integral = EulerGamma + ln(|x|) + integral((cosh(t)-1)/t, t=0..|x|)}

function sfc_ci(x: extended): extended;
  {-Return the cosine integral, ci(x) = EulerGamma + ln(|x|) + integral((cos(t)-1)/t, t=0..|x|)}

function sfc_cin(x: extended): extended;
  {-Return the entire cosine integral, cin(x) = integral((1-cos(t))/t, t=0..x)}

function sfc_cinh(x: extended): extended;
  {-Return the entire hyperbolic cosine integral, cinh(x) = integral((cosh(t)-1)/t, t=0..x)}

function sfc_e1(x: extended): extended;
  {-Return the exponential integral E1(x) = integral(exp(-x*t)/t, t=1..Inf), x <> 0}

function sfc_e1s(x: extended): extended;
  {-Return E1s(x) = exp(x)*E1(x), x <> 0}

function sfc_ei(x: extended): extended;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}

function sfc_ein(x: extended): extended;
  {-Return the entire exponential integral ein(x) = integral((1-exp(-t))/t, t=0..x)}

function sfc_eis(x: extended): extended;
  {-Return exp(-x)*Ei(x), x <> 0}

function sfc_eisx2(x: extended): extended;
  {-Return exp(-x^2)*Ei(x^2), x <> 0}

function sfc_ei_inv(x: extended): extended;
  {-Return the functional inverse of Ei(x), ei_inv(ei(x))=x}

function sfc_en(n: longint; x: extended): extended;
  {-Return the exponential integral E_n(x) = integral(exp(-x*t)/t^n, t=1..Inf), x > 0}

function sfc_eibeta(n: integer; x:extended): extended;
  {-Return the exponential integral beta(n,x) = int(t^n*exp(-x*t), t=-1..1) }

function sfc_gei(p,x: extended): extended;
  {-Return the generalized exponential integral E_p(x) = integral(exp(-x*t)/t^p, t=1..Inf), x >= 0}

function sfc_li(x: extended): extended;
  {-Return the logarithmic integral li(x) = PV-integral(1/ln(t), t=0..x), x >= 0, x <> 1}

function sfc_shi(x: extended): extended;
  {-Return the hyperbolic sine integral, shi(x) = integral(sinh(t)/t, t=0..x)}

function sfc_si(x: extended): extended;
  {-Return the sine integral, si(x) = integral(sin(t)/t, t=0..x)}

function sfc_ssi(x: extended): extended;
  {-Return the shifted sine integral, ssi(x) = si(x) - Pi/2}


implementation

uses
  AMath,
  {$ifopt R+}
    sfBasic,  {for RTE_ArgumentRange}
  {$endif}
  sfGamma,
  sfGamma2;




{---------------------------------------------------------------------------}
function e1_small(x: extended): extended;
  {-Return the exponential integral E1(x) for 0 < x <= 1}
var
  y: extended;
const
  P01: array[0..5] of extended = (
         0.08651972480793979568216,
         0.0275114007037026844633,
        -0.246594388074877139824,
        -0.0237624819878732642231,
        -0.00259113319641673986276,
         0.30853660894346057053e-4);
  Q01: array[0..6] of extended = (
         1.0,
         0.317978365797784100273,
         0.0393622602554758722511,
         0.00204062029115966323229,
         0.732512107100088047854e-5,
        -0.202872781770207871975e-5,
         0.52779248094603709945e-7);
const
  y01 = 695977.0/1048576.0;
begin
  {Based on parts of boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
  y :=  PolEvalX(x,P01,6) / PolEvalX(x,Q01,7);
  e1_small := y + (x - ln(x) - y01);
end;


{---------------------------------------------------------------------------}
function e1ex_large(x: extended): extended;
  {-Return the exponential integral exp(x)*E1(x) for x > 1}
var
  y,z: extended;
const
  Pix: array[0..13] of extended = (
         -0.534401189080684443046e-23,
         -0.999999999999999999905,
         -62.1517806091379402505,
         -1568.45688271895145277,
         -21015.3431990874009619,
         -164333.011755931661949,
         -777917.270775426696103,
         -2244188.56195255112937,
         -3888702.98145335643429,
         -3909822.65621952648353,
         -2149033.9538897398457,
         -584705.537139793925189,
         -65815.2605361889477244,
         -2038.82870680427258038);
  Qix: array[0..13] of extended = (
         1.0,
         64.1517806091379399478,
         1690.76044393722763785,
         24035.9534033068949426,
         203679.998633572361706,
         1074661.58459976978285,
         3586552.65020899358773,
         7552186.84989547621411,
         9853333.79353054111434,
         7689642.74550683631258,
         3385553.35146759180739,
         763218.072732396428725,
         73930.2995984054930821,
         2063.86994219629165937);
begin
  {Based on parts of boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
  y := 1.0/x;
  z := 1.0 + PolEvalX(y,Pix,14) / PolEvalX(y,Qix,14);
  e1ex_large := z*y;
end;


{---------------------------------------------------------------------------}
function sfc_e1(x: extended): extended;
  {-Return the exponential integral E1(x) = integral(exp(-x*t)/t, t=1..Inf), x <> 0}
var
  y: extended;
begin
  {Based on boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
  if IsNaN(x) then sfc_e1 := x
  else if x <  0.0 then sfc_e1 := -sfc_ei(-x)
  else if x <= 1.0 then sfc_e1 := e1_small(x)
  else if x >= ln_MaxExt then sfc_e1 := 0.0
  else begin
    y := e1ex_large(x);
    sfc_e1 := exp(-x)*y;
  end;
end;

{---------------------------------------------------------------------------}
function ei_asymp(x: extended): extended;
  {-Return asymptotic expansion for exp(-x)*Ei(x), x -> Inf}
var
  s,t,d,n: extended;
  k: integer;
begin
  {Asymptotic expansion NIST[30], 6.12.2:  1/x * sum(k!/x^k)}
  t := 1.0;
  s := 1.0;
  n := 1.0;
  k := 1;
  repeat
    n := n*k;
    t := t/x;
    d := n*t;
    s := s+d;
    inc(k);
  until d <= eps_x*s;
  ei_asymp := s/x;
end;


{---------------------------------------------------------------------------}
function sfc_eis(x: extended): extended;
  {-Return exp(-x)*Ei(x), x <> 0}
const
  xas = 100.0; {double: 80}
begin
  if IsNaN(x) then sfc_eis := x
  else if x < -1.0 then sfc_eis := -sfc_e1s(-x)
  else if x < xas then sfc_eis := exp(-x)*sfc_ei(x)
  else sfc_eis := ei_asymp(x);
end;


{---------------------------------------------------------------------------}
function sfc_e1s(x: extended): extended;
  {-Return E1s(x) = exp(x)*E1(x), x <> 0}
begin
  {Note that E1s(x) = Re(U(1,1,x))}
  if IsNaN(x) then sfc_e1s := x
  else if x <  0.0 then sfc_e1s := -sfc_eis(-x)
  else if x <= 1.0 then sfc_e1s := exp(x)*e1_small(x)
  else sfc_e1s := e1ex_large(x);
end;


{---------------------------------------------------------------------------}
function sfc_eiex(x,expx: extended): extended;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}
  { internal version for x>6, expx must be = exp(x), used by ei and li}
var
  y,z: extended;
const
  {Based on boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
  {Copyright John Maddock 2007, see 3rdparty.ama for Boost license}
  P10: array[0..8] of extended = (
         0.00139324086199409049399,
        -0.0345238388952337563247,
        -0.0382065278072592940767,
        -0.0156117003070560727392,
        -0.00383276012430495387102,
        -0.000697070540945496497992,
        -0.877310384591205930343e-4,
        -0.623067256376494930067e-5,
        -0.377246883283337141444e-6);
  Q10: array[0..9] of extended = (
         1.0,
         1.08073635708902053767,
         0.553681133533942532909,
         0.176763647137553797451,
         0.0387891748253869928121,
         0.0060603004848394727017,
         0.000670519492939992806051,
         0.4947357050100855646e-4,
         0.204339282037446434827e-5,
         0.146951181174930425744e-7);
  P20: array[0..9] of extended = (
        -0.00893891094356946995368,
        -0.0487562980088748775943,
        -0.0670568657950041926085,
        -0.0509577352851442932713,
        -0.02551800927409034206,
        -0.00892913759760086687083,
        -0.00224469630207344379888,
        -0.000392477245911296982776,
        -0.44424044184395578775e-4,
        -0.252788029251437017959e-5);
  Q20: array[0..9] of extended = (
         1.0,
         2.00323265503572414261,
         1.94688958187256383178,
         1.19733638134417472296,
         0.513137726038353385661,
         0.159135395578007264547,
         0.0358233587351620919881,
         0.0056716655597009417875,
         0.000577048986213535829925,
         0.290976943033493216793e-4);
  P40: array[0..11] of extended = (
        -0.00356165148914447278177,
        -0.0240235006148610849678,
        -0.0516699967278057976119,
        -0.0586603078706856245674,
        -0.0409960120868776180825,
        -0.0185485073689590665153,
        -0.00537842101034123222417,
        -0.000920988084778273760609,
        -0.716742618812210980263e-4,
        -0.504623302166487346677e-9,
         0.712662196671896837736e-10,
        -0.533769629702262072175e-11);
  Q40: array[0..8] of extended = (
         1.0,
         3.13286733695729715455,
         4.49281223045653491929,
         3.84900294427622911374,
         2.15205199043580378211,
         0.802912186540269232424,
         0.194793170017818925388,
         0.0280128013584653182994,
         0.00182034930799902922549);
const
  Pix: array[0..8] of extended = (
        -0.0130653381347656250004,
         0.644487780349757303739,
         143.995670348227433964,
        -13918.9322758014173709,
         476260.975133624194484,
        -7437102.15135982802122,
         53732298.8764767916542,
        -160695051.957997452509,
         137839271.592778020028);
  Qix: array[0..8] of extended = (
         1.0,
         27.2103343964943718802,
        -8785.48528692879413676,
         397530.290000322626766,
        -7356441.34957799368252,
         63050914.5343400957524,
        -246143779.638307701369,
         384647824.678554961174,
        -166288297.874583961493);
const
  y10 = 303821.0/262144.0;
  y20 = 569887.0/524288.0;
  y40 = 136233.0/131072.0;
  yix = 265569.0/262144.0;
begin
  if x>ln_MaxExt then sfc_eiex := PosInf_x
  else begin
    if x<=10.0 then begin
      y := 0.5*x - 4.0;
      z := y10 + PolEvalX(y,P10,9) / PolEvalX(y,Q10,10);
    end
    else if x<=20.0 then begin
      y := x/5.0 - 3.0;
      z := y20 + PolEvalX(y,P20,10) / PolEvalX(y,Q20,10);
    end
    else if x<=40.0 then begin
      y := x/10.0 - 3.0;
      z := y40 + PolEvalX(y,P40,12) / PolEvalX(y,Q40,9);
    end
    else begin
      y := 1.0/x;
      z := yix + PolEvalX(y,Pix,9) / PolEvalX(y,Qix,9);
    end;
    z := z*(expx/x);
    sfc_eiex := z + x;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ei(x: extended): extended;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}
var
  y,z: extended;
{Based on boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
{Copyright John Maddock 2007, see 3rdparty.ama for Boost license}

{WE note: Do not change P6/Q6 to decimal! This change will increase}
{         the relative error around x=6 to about 15..20 eps_x}
const
  P6H: array[0..10] of THexExtW = (
         ($7A8C,$8A1E,$46C4,$BF27,$4000),  { 2.98677224343598593764}
         ($FD57,$D5A6,$A7DC,$8490,$3FFD),  { 0.25891613550886736592}
         ($EFE7,$87E0,$1C47,$CA11,$3FFE),  { 0.789323584998672832285}
         ($016D,$549D,$4BEE,$BD4D,$3FFB),  { 0.092432587824602399339}
         ($EE5D,$A73C,$A7C9,$D2A1,$3FFA),  { 0.0514236978728625906656}
         ($D5AB,$1103,$17FD,$D7C5,$3FF7),  { 0.00658477469745132977921}
         ($AD92,$84F1,$5D20,$A3BA,$3FF5),  { 0.00124914538197086254233}
         ($3827,$AC7C,$62CE,$89D0,$3FF2),  { 0.000131429679565472408551}
         ($0427,$B705,$7D19,$BD78,$3FEE),  { 0.11293331317982763165e-4}
         ($5B5D,$FA07,$DC80,$A8FA,$3FEA),  { 0.629499283139417444244e-6}
         ($98C6,$B2D3,$E672,$98C1,$3FE5)); { 0.177833045143692498221e-7}
  Q6H: array[0..8] of THexExtW = (
         ($0000,$0000,$0000,$8000,$3FFF),  { 1.0}
         ($726A,$E11E,$1134,$9A0D,$BFFF),  {-1.20352377969742325748}
         ($21BC,$21B8,$B14D,$AAC5,$3FFE),  { 0.66707904942606479811}
         ($04C7,$8509,$EBDF,$E45D,$BFFC),  {-0.223014531629140771914}
         ($788A,$488A,$7362,$CA12,$3FFA),  { 0.0493340022262908008636}
         ($403E,$5311,$F531,$F31D,$BFF7),  {-0.00741934273050807310677}
         ($05C2,$B8B2,$D5AD,$C2E9,$3FF4),  { 0.00074353567782087939294}
         ($2CD6,$4141,$C7BD,$BF33,$BFF0),  {-0.455861727069603367656e-4}
         ($628A,$9DB7,$5B80,$B084,$3FEB)); { 0.131515429329812837701e-5}
const
  rh : THexExtW = ($E1E5,$A5AF,$4A95,$BEB9,$3FFD); { = 0.37250741078136663446199186658.. ; root of Ei: Ei(r)=0}
  r1h: THexExtW = ($E000,$A5AF,$4A95,$BEB9,$3FFD); { = 1677624236387711.0/4503599627370496.0 = 3.7250741078136662132E-1}
  r2 = 0.13140183414386028200928e-16;              { = r - r1}
var
  r : extended absolute rh;
  r1: extended absolute r1h;
  P06: array[0..10]of extended absolute P6H;
  Q06: array[0..8] of extended absolute Q6H;
begin
  {avoid infinite loop if x is (negative) NAN}
  if IsNaN(x) then sfc_ei := x
  else if x<0 then sfc_ei := -sfc_e1(-x)
  else if x<=eps_x then begin
    {Boost's rational approximation is suboptimal for small x}
    {use Ei(x) = (+EulerGamma + ln(x)) + x + 1/4*x^2 + O(x^3)}
    {for x <= eps_x, i.e. |ln(x)| > 43, the x term is negligible}
    sfc_ei := EulerGamma + ln(x);
  end
  else if x<=6.0 then begin
    y := x/THREE - 1.0;
    z := PolEvalX(y,P06,11) / PolEvalX(y,Q06,9);
    y := (x - r1) - r2;
    z := y*z;
    if abs(y)<0.1 then sfc_ei := z + ln1p(y/r)
    else sfc_ei := z + ln(x/r);
  end
  else if x>ln_MaxExt then sfc_ei := PosInf_x
  else sfc_ei := sfc_eiex(x,exp(x));
end;


{---------------------------------------------------------------------------}
function sfc_eisx2(x: extended): extended;
  {-Return exp(-x^2)*Ei(x^2)}
const
  xei  = 20.0;  {double: 10}
  xmax = 1e20;  {double: 1e16}
var
  z,s: extended;
begin
  if IsNaN(x) then begin
    sfc_eisx2 := x;
    exit;
  end;
  x := abs(x);
  if x <= xei then begin
    {return exp(-x^2)*Ei(x^2)}
    s := expx2(-x);
    z := sfc_ei(sqr(x));
    sfc_eisx2 := s*z;
  end
  else begin
    if x >= xmax then begin
      {Only one term of AE, avoid overflow}
      sfc_eisx2 := (1.0/x)/x;
    end
    else begin
      {Asymptotic expansion NIST[30], 6.12.2:  1/x^2 * sum(k!/x^(2k)}
      sfc_eisx2 := ei_asymp(sqr(x));
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ein(x: extended): extended;
  {-Return the entire exponential integral ein(x) = integral((1-exp(-t))/t, t=0..x)}
const
  csh: array[0..11] of THexExtw = (
         ($B286,$BFBB,$F80D,$8038,$4000),  {+2.003477109362331949501523562122308610572     }
         ($11B1,$136B,$12FA,$8040,$BFFB),  {-0.6262221170111595696745652349876399745063e-1 }
         ($B1C3,$0C96,$8ECC,$E3FB,$3FF5),  {+0.1739369565279582463382068320644185874193e-2 }
         ($FFD0,$1FDE,$95D1,$AAF6,$BFF0),  {-0.4076080883130112580271691400428456876713e-4 }
         ($CAC2,$44E8,$49C9,$DACD,$3FEA),  {+0.8151006219468242274136388199101306590141e-6 }
         ($32F0,$38DA,$EC9F,$F314,$BFE4),  {-0.1414921923800105206557069349069463653541e-7 }
         ($D203,$586B,$F6EE,$EE17,$3FDE),  {+0.2165448782930998188223426094586382995172e-9 }
         ($5F35,$9C4E,$9E9C,$D04F,$BFD8),  {-0.2960277555152250090455326359589754482108e-11}
         ($76D7,$7327,$B417,$A493,$3FD2),  {+0.3654342755726724875204322769851602486817e-13}
         ($1B2B,$980C,$F977,$ECF8,$BFCB),  {-0.4110818570387228515941851406185706129602e-15}
         ($4067,$C135,$5831,$9CAA,$3FC5),  {+0.4246424399152657058898802616196579805661e-17}
         ($365F,$9FC8,$E7EF,$BF77,$BFBE)); {-0.4054500715002447619255009744247776959447e-19}
var
  csa: array[0..11] of extended absolute csh;
var
  t,z: extended;
begin
  {ein(x) = gamma + ln(|x|) - Ei(-x). Ref: NIST[30], 6.2(i)}
  t := abs(x);
  if t<0.25 then begin
    if t<=eps_x then sfc_ein := x
    else begin
      {Chebyshev approximation calculated with Maple from NIST[30], 6.6.4}
      t := CSEvalx(4.0*x, csa, 12);
      sfc_ein := t*x;
    end;
  end
  else if x>0 then begin
    t := ln(t);
    if x >= 45.0 then sfc_ein := t + EulerGamma
    else begin
      z := EulerGamma + sfc_e1(x);
      sfc_ein := z+t;
    end;
  end
  else begin
    z := -sfc_ei(-x);
    if x <= -50 then sfc_ein := z
    else begin
      t := ln(t) + EulerGamma;
      sfc_ein := z+t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ei_inv(x: extended): extended;
  {-Return the functional inverse of Ei(x), ei_inv(ei(x))=x}
var
  z,d: extended;
const
  x0 = -43.1;         { < ln(eps) + Eulergamma, double = -35.5}
  x1 = 0.10477e4929;  {~ max. argument without overflow ~ MaxDouble/ln_MaxDbl}
  z0 = 0.37249755859375 + 0.985218761663446199e-5; {Zero of Ei}
  b  = 3.89621573390716731;  {Taylor expansion Ei(z0+z) = b*(z-z0) + O((z-z0)^2)}
  eg = 0.5614594835668851698; {exp(-Eulergamma)}
begin
  {Basic formula: li(x) = Ei(ln(x)) -> ei_inv(x) = ln(li_inv(x))}
  if IsNan(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_ei_inv := NaN_x;
  end
  else if x<=x0 then sfc_ei_inv := eg*exp(x) {exp(x-Eulergamma)}
  else if abs(x) < sqrt_epsh then sfc_ei_inv := z0 + x/b
  else if x>0.0 then begin
    if x>x1 then sfc_ei_inv := PosInf_x
    else sfc_ei_inv := ln(sfc_ali(x))
  end
  else begin
    z := ln(sfc_ali(x));
    if x < -2.0 then begin
      {Here the basic formula is used as starting value for Newton iteration}
      {because numerically it does not always give satisfactory accuracy.   }
      repeat
        d := sfc_ei(z)-x;
        d := d*z/exp(z);
        z := z-d;
      until abs(d) <= sqrt_epsh*abs(z);
    end;
    sfc_ei_inv := z;
  end
end;


{---------------------------------------------------------------------------}
function sfc_li(x: extended): extended;
  {-Return the logarithmic integral li(x) = PV-integral(1/ln(t), t=0..x), x >= 0, x <> 1}
var
  t: extended;
begin
  {li(x) = Ei(ln(x)). Note that some authors use Li(x) = li(x)-li(2)}
  if x=0.0 then sfc_li := 0.0
  else begin
    t := ln(x);
    if t<6.0 then sfc_li := sfc_ei(t)
    else begin
      {Avoid the inaccuracies of exp(ln(x)) for large x}
      {and use the reformulated Ei code for ln(x) > 6}
      sfc_li := sfc_eiex(t,x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ali(x: extended): extended;
  {-Return the functional inverse of li(x), li(sfc_ali(x))=x}
var
  d,z: extended;
  k: integer;
begin
  k := 0;
  {Initial approximation}
  if x > 3.5 then z := x*ln(x)
  else if x > 0.75 then z := 1.0 + x
  else if x > -0.5 then z := 1.45137 + 0.37251*x  {linear near zero 1.45137 of li}
  else if x > -43.8 then z := 1.0 + exp(x-EulerGamma)
  else z := 1.0;

  {Skip iteration at the singularity of li}
  if z>1.0 then begin
    repeat
      {Halley iteration for li(z)-x=0. We have li'(x) = 1/ln(x)}
      {and li''(x) = -1/(x*ln(x)^2)}
      inc(k);
      d := sfc_li(z) - x;
      if d<>0.0 then begin
        d := d*ln(z)/(1.0 + 0.5*d/z);  {Newton if /(1+..) is omitted}
        z := z - d;
      end;
      {If |d| < z*sqrt_epsh the correction from the next iteration step}
      {can be neglected because the convergence is at least quadratic. }
    until (abs(d) < z*sqrt_epsh) or (k>4);
  end;
  sfc_ali := z;
end;


{---------------------------------------------------------------------------}
function en_series(n: longint; x: extended): extended;
  {-Calculate E_n(x), 0<x<=1, n>=1, via series expansion}
var
  m,k: longint;
  f,p,s,z: extended;
begin

  {HMF [1], formula 5.1.12}
  m := 0;
  k := 1-n;
  f := 1;
  z := -x;

  {initialize sum s with m=0 term}
  if n=1 then s := 0.0
  else s := 1.0/k;

  {sum until first term with relative contribution is <= eps_x}
  repeat
    inc(m);
    inc(k);
    f := f*z/m;
    {skip x^(n-1) term}
    if k<>0 then s := s + f/k;
  until abs(f) < eps_x*abs(s);

  if k<0 then begin
    {sum converged before reaching x^(n-1) term, skip psi calculation}
    en_series := -s;
  end
  else begin
    if n<13 then begin
      {for small n calculate psi(n) and z^(n-1)/(n-1)! via loop}
      p := 0.0;
      f := 1.0;
      for k:=n-1 downto 1 do begin
        f := f*z/k;
        p := p + 1.0/k;
      end;
      p := p - (EulerGamma + ln(x));
      en_series := p*f - s;
    end
    else begin
      {since |z| <= 1, z^(n-1) cannot overflow}
      if n > MAXGAMX then en_series := -s
      else begin
        {en_series := (psi(n)-ln(x))*z^(n-1)/gamma(n) - s;}
        p := sfc_psi(n)-ln(x);
        f := intpower(z,n-1)/sfc_gamma(n);;
        en_series := p*f - s;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function en_cfrac(n: longint; x: extended): extended;
  {-Calculate E_n(x), x>1, n>1, via continued fraction}
var
  k: longint;
  big,cf,pk,p1,p2,qk,q1,q2,r,t: extended;
  ev: boolean;
begin
  {Algorithm based on Cephes [7] function expn}

  big:= ldexp(1.0,68); {Rescaling threshold, should be power of 2}

  k  := 0;
  p2 := 1.0;
  q2 := x;
  p1 := 1.0;
  q1 := x+n;
  cf := p1/q1;
  ev := true;

  repeat
    if ev then begin
      inc(k);
      pk := p1*x + p2*k;
      qk := q1*x + q2*k;
    end
    else begin
      t  := n  + k;
      pk := p1 + p2*t;
      qk := q1 + q2*t;
    end;
    ev := not ev;

    if qk<>0.0 then begin
      r := pk/qk;
      t := abs((cf-r)/r);
      cf:= r;
    end
    else t := 1;

    p2 := p1;
    p1 := pk;
    q2 := q1;
    q1 := qk;
    if abs(pk) > big then begin
      p1 := p1/big;
      p2 := p2/big;
      q1 := q1/big;
      q2 := q2/big;
    end;
  until t<=eps_x;
  en_cfrac := cf*exp(-x);
end;

(*
{---------------------------------------------------------------------------}
function en_largen(n: longint; x: extended): extended;
  {-Value of en(n,x) for large n, n>100000, x < 11400}
var
  s,t,xn,yn: extended;
begin
  {HMF [1], 5.1.52, R(n,x) <= (1+1/(x+n+1))/n^4 < 1e-20}
  {see also Cephes [7] function expn}
  xn := x + n;
  yn := 1.0/sqr(xn);
  t  := n;
  s  := yn*t*(6.0*x*x - 8.0*t*x + t*t);
  s  := yn*(s + t*(t-2.0*x));
  s  := yn*(s + t);
  en_largen := (1.0+s)*exp(-x)/xn;
end;
*)

{---------------------------------------------------------------------------}
function ep_largep(p,x: extended): extended;
  {-Uniform asymptotic expansion E_p(x) for large p; called with p>=10000, 0 <= x <= 11390.2}
var
  s,z,a2,a3,a4: extended;
begin
  {NIST 8.20.6}
  z  := x/p;
  a2 := 1.0 - 2.0*z;
  a3 := (6.0*z - 8.0)*z + 1.0;
  a4 := ((58.0 - 24.0*z)*z - 22.0)*z + 1.0;
  z  := 1.0/(sqr(1.0 + z)*p);
  s  := (a4*z + a3)*z + a2;
  s  := (s*z + 1.0)*z + 1.0;
  ep_largep := exp(-x)*s/(x+p);
end;


{---------------------------------------------------------------------------}
function sfc_en(n: longint; x: extended): extended;
  {-Return the exponential integral E_n(x) = integral(exp(-x*t)/t^n, t=1..Inf), x > 0}
begin
  if IsNaNorInf(x)  then sfc_en := NaN_x
  else if n<0       then sfc_en := sfc_gei(n,x)
  else if n=0       then sfc_en := exp(-x)/x {HMF [1], 5.1.24}
  else if n=1       then sfc_en := sfc_e1(x)
  else if x=0.0     then sfc_en := 1.0/(n-1) {HMF [1], 5.1.23}
  else if x<0.0     then sfc_en := ln(x)     {force error}
  else if x>11390.2 then sfc_en := 0.0       {HMF [1], 5.1.19/20: 0 < en <= exp(-x)/x}
  else if n>=10000  then sfc_en := ep_largep(n,x)
  else if x<=1.0    then sfc_en := en_series(n,x)
  else                   sfc_en := en_cfrac(n,x);
end;


{---------------------------------------------------------------------------}
function eib_small(n: integer; x: extended): extended;
  {-Exponential integral beta(n,x) for 'small' x}
var
  k,q,m: longint;
  s,t,z: extended;
const
  KMAX = 100;
begin
  {Beta via igamman suffers from cancellation, use series computed with Maple:}
  {beta(n,s) = ((-1)^n+1)/(n+1) + sum(((-1)^n+(-1)^k)/k!/(k+n+1)*x^k, k=1..INF)}

  if odd(n) then q := -1 else q := 1; {q = (-1)^n}
  s := (q+1)/(n+1);
  if x<>0.0 then begin
    m := -1;  {q = (-1)^k}
    z := x;   {z = x^k   }
    for k:=1 to KMAX do begin
      t := q+m;
      if t<>0.0 then begin
        t := t/(n+k+1);
        t := t*z;
        s := s + t;
        if abs(t)<=eps_x*abs(s) then break;
      end;
      m := -m;
      z := z*x/(k+1);
      {$ifopt R+}
        if (k=KMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
      {$endif}
    end;
  end;
  eib_small := s;
end;


{---------------------------------------------------------------------------}
function sfc_eibeta(n: integer; x:extended): extended;
  {-Return the exponential integral beta(n,x) = int(t^n*exp(-x*t), t=-1..1) }
var
  f,g,y: extended;
begin
  {HMF[1], 5.1.46}
  if n<=0 then begin
    if n<>0 then sfc_eibeta := Nan_x
    else sfc_eibeta := 2.0*sinhc(x);
  end
  else if abs(x) <= 2.0 then begin
    sfc_eibeta := eib_small(n,x);
  end
  else begin
    n := n+1;
    f := sfc_igamma(n,-x);
    g := sfc_igamma(n,x);
    y := power(x, -n);
    sfc_eibeta := (f-g)*y;
  end;
end;

{---------------------------------------------------------------------------}
function sfc_gei(p,x: extended): extended;
  {-Return the generalized exponential integral E_p(x) = integral(exp(-x*t)/t^p, t=1..Inf), x >= 0}
var
  a,q,t: extended;
begin
  if (x<0.0) or IsNanOrInf(x) or IsNanOrInf(p) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_gei := NaN_x;
    exit;
  end;
  {E_p(x) = x^(p-1)*GAMMA(1-p,x), if x>0}
  if (p>=0.0) and (frac(p)=0.0) and (p<MaxLongInt) then sfc_gei := sfc_en(round(p),x)
  else begin
    if x=0.0 then begin
      if p>1.0 then sfc_gei := 1.0/(p-1.0)
      else sfc_gei := PosInf_x;
    end
    else if (x>11390.2) and (p>=0.0) then sfc_gei := 0.0
    else if p>=10000.0 then sfc_gei := ep_largep(p,x)
    else if ((x>=1.0) and (p>=-0.25)) or ((p>1.0) and (x>0.25)) then begin
      {sfc_igamma will use igam_qfraction for these (p,x) values, so}
      {compute CF directly without scaling & rescaling with x^(p-1)!}
      t := exp(-x);
      if t<>0 then igam_qfraction(1.0-p,x,t,a,t);
      sfc_gei := t;
    end
    else begin
      a := 1.0 - p;
      if a >= MAXGAMX then begin
        {Gamma(a) will overflow, compute Q(a,x)*Gamma(a)/x^a using logarithmic forms}
        sfc_incgamma_ex(a,x,p,q,t,false);
        if q=0.0 then sfc_gei := 0.0
        else begin
          t := sfc_lngamma(a) - a*ln(x);
          if q<0.75 then t := t + ln(q)
          else t := t + ln1p(-p);
          if t>ln_MaxExt then sfc_gei := PosInf_x
          else sfc_gei := exp(t);
        end;
      end
      else begin
        t := power(x,-a);
        a := sfc_igamma(a,x);
        {use PosInf_x if result would overflow}
        if (a<=1.0) or (t<MaxExtended/a) then sfc_gei := a*t
        else sfc_gei := PosInf_x;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_chi(x: extended): extended;
  {-Return the hyperbolic cosine integral = EulerGamma + ln(|x|) + integral((cosh(t)-1)/t, t=0..|x|)}
begin
  x := abs(x);
  {e1(x)/ei(x) < 1e-20 for x>23: skip calculation of e1 in this case}
  if x>23.0 then sfc_chi := 0.5*sfc_ei(x)
  else sfc_chi := 0.5*(sfc_ei(x) - sfc_e1(x));
end;


{---------------------------------------------------------------------------}
function sfc_shi(x: extended): extended;
  {-Return the hyperbolic sine integral, shi(x) = integral(sinh(t)/t, t=0..x)}
var
  z: extended;
const
  ShiCS : array[0..6] of extended = (
            7.83726856889009506952e-3,
            3.92276649342345639727e-3,
            4.13467878876172667467e-6,
            2.47074803728827421351e-9,
            9.37929559076363045716e-13,
            2.451817019520867353e-16,
            4.67416155257592e-20);
begin
  {W. Fullerton [20], file dshi.f}
  z := abs(x);
  if z<=1e-5 then begin
    {shi(x) = x + 1/18*x^3 + 1/600*x^5 + O(x^6)}
    if z<=1e-10 then sfc_shi := x
    else sfc_shi := x*(1.0 + x*x/18.0);
  end
  else begin
    if z<=0.375 then z := z*(1.0 + CSEvalX(128.0*sqr(z)/9.0-1.0, ShiCS, 7))
    else if z>23.0 then z := 0.5*sfc_ei(z)
    else z := 0.5*(sfc_ei(z)+sfc_e1(z));
    {shi(-x) = -shi(x)}
    if x>0 then sfc_shi := z
    else sfc_shi := -z;
  end;
end;


{---------------------------------------------------------------------------}
procedure auxfg(x: extended; var f,g: extended);
  {-Return the auxiliary functions f,g for ci,si, x>=4}
var
  z: extended;
const
  xbig  = 4.29496729600000000E9;     {sqrt(2/eps_x)}
  xmaxf = 2.97432873839307945E4931;  {exp(min(-ln(Minextended), ln(Maxextended) - 0.01)}
  xmaxg = 5.45374067809707965E2465;  {1/sqrt(Minextended)}
  xbnd  = 7.07106781186547524;       {sqrt(50)}
  xbndg = 14.1421356237309505;       {sqrt(200)}
const
  n_f1 = 25;
  n_f2 = 44;
  n_g1 = 27;
  n_g2 = 25;
  n_g3 = 26;
const
  f1: array[0..n_f1-1] of extended = (
        -0.1191081969051363610348201965828918,
        -0.0247823144996236247590074150823133,
        +0.0011910281453357821268120363054457,
        -0.0000927027714388561748308600360706,
        +0.0000093373141568270996868204582766,
        -0.0000011058287820557143938979426306,
        +0.0000001464772071460162169336550799,
        -0.0000000210694496287689532601227548,
        +0.0000000032293492366848236382857374,
        -0.0000000005206529617529375828014986,
        +0.0000000000874878884570278750268316,
        -0.0000000000152176187056123668294574,
        +0.0000000000027257192405419573900583,
        -0.0000000000005007053075968556290255,
        +0.0000000000000940240902726068511779,
        -0.0000000000000180014444791803678336,
        +0.0000000000000035062621432741785826,
        -0.0000000000000006935282926769149709,
        +0.0000000000000001390925136454216568,
        -0.0000000000000000282486885074170585,
        +0.0000000000000000058031305693579081,
        -0.0000000000000000012046901573375820,
        +0.0000000000000000002525052443655940,
        -0.0000000000000000000533980268805594,
        +0.0000000000000000000113855786274122);
const
  f2: array[0..n_f2-1] of extended = (
        -0.03484092538970132330836049733745577,
        -0.01668422056779596873246786312278676,
        +0.00067529012412377385045207859239727,
        -0.00005350666225447013628785577557429,
        +0.00000626934217790075267050759431626,
        -0.00000095266388019916680677790414293,
        +0.00000017456292242509880425504427666,
        -0.00000003687954030653093307097646628,
        +0.00000000872026777051395264075816938,
        -0.00000000226019703919738748530423167,
        +0.00000000063246249765250612520444877,
        -0.00000000018889118884717869240911480,
        +0.00000000005967746729997813372620472,
        -0.00000000001980443117372239011196007,
        +0.00000000000686413954772103383713264,
        -0.00000000000247310193070199106074890,
        +0.00000000000092263594549941404196042,
        -0.00000000000035523634999261784497297,
        +0.00000000000014076049625351591461820,
        -0.00000000000005726228499747652794311,
        +0.00000000000002386537545413171810106,
        -0.00000000000001017141890764597142232,
        +0.00000000000000442594531078364424968,
        -0.00000000000000196344933049189761979,
        +0.00000000000000088688748314810461024,
        -0.00000000000000040743345027311546948,
        +0.00000000000000019016837215675339859,
        -0.00000000000000009009707297478042442,
        +0.00000000000000004329211274095668667,
        -0.00000000000000002108144465322479526,
        +0.00000000000000001039637907026452274,
        -0.00000000000000000518891007948931936,
        +0.00000000000000000261955324869899371,
        -0.00000000000000000133690399951301570,
        +0.00000000000000000068941057702931664,
        -0.00000000000000000035905362610437250,
        +0.00000000000000000018878077255791706,
        -0.00000000000000000010016125265594380,
        +0.00000000000000000005360725691578228,
        -0.00000000000000000002893198974944827,
        +0.00000000000000000001574065100202625,
        -0.00000000000000000000863027106431206,
        +0.00000000000000000000476715602862288,
        -0.00000000000000000000265222739998504);
const
  g1: array[0..n_g1-1] of extended = (
        -0.3040578798253495954499726682091083,
        -0.0566890984597120587731339156118269,
        +0.0039046158173275643919984071554082,
        -0.0003746075959202260618619339867489,
        +0.0000435431556559843679552220840065,
        -0.0000057417294453025046561970723475,
        +0.0000008282552104502629741937616492,
        -0.0000001278245892594642727883913223,
        +0.0000000207978352948687884439257529,
        -0.0000000035313205921990798042032682,
        +0.0000000006210824236308951068631449,
        -0.0000000001125215474446292649336987,
        +0.0000000000209088917684421605267019,
        -0.0000000000039715831737681727689158,
        +0.0000000000007690431314272089939005,
        -0.0000000000001514696742731613519826,
        +0.0000000000000302892146552359684119,
        -0.0000000000000061399703834708825400,
        +0.0000000000000012600605829510933553,
        -0.0000000000000002615029250939483683,
        +0.0000000000000000548278844891796821,
        -0.0000000000000000116038182129526571,
        +0.0000000000000000024771654107129795,
        -0.0000000000000000005330672753223389,
        +0.0000000000000000001155666075598465,
        -0.0000000000000000000252280547744957,
        +0.0000000000000000000055429038550786);
const
  g2: array[0..n_g2-1] of extended = (
        -0.1211802894731646263541834046858267,
        -0.0316761386394950286701407923505610,
        +0.0013383199778862680163819429492182,
        -0.0000895511011392252425531905069518,
        +0.0000079155562961718213115249467924,
        -0.0000008438793322241520181418982080,
        +0.0000001029980425677530146647227274,
        -0.0000000139295750605183835795834444,
        +0.0000000020422703959875980400677594,
        -0.0000000003196534694206427035434752,
        +0.0000000000528147832657267698615312,
        -0.0000000000091339554672671033735289,
        +0.0000000000016426251238967760444819,
        -0.0000000000003055897039322660002410,
        +0.0000000000000585655825785779717892,
        -0.0000000000000115229197730940120563,
        +0.0000000000000023209469119988537310,
        -0.0000000000000004774355834177535025,
        +0.0000000000000001000996765800180573,
        -0.0000000000000000213533778082256704,
        +0.0000000000000000046277190777367671,
        -0.0000000000000000010175807410227657,
        +0.0000000000000000002267657399884672,
        -0.0000000000000000000511630776076426,
        +0.0000000000000000000116767014913108);
const
  g3: array[0..n_g3-1] of extended = (
        -0.0280574367809472928402815264335299,
        -0.0137271597162236975409100508089556,
        +0.0002894032638760296027448941273751,
        -0.0000114129239391197145908743622517,
        +0.0000006813965590726242997720207302,
        -0.0000000547952289604652363669058052,
        +0.0000000055207429918212529109406521,
        -0.0000000006641464199322920022491428,
        +0.0000000000922373663487041108564960,
        -0.0000000000144299088886682862611718,
        +0.0000000000024963904892030710248705,
        -0.0000000000004708240675875244722971,
        +0.0000000000000957217659216759988140,
        -0.0000000000000207889966095809030537,
        +0.0000000000000047875099970877431627,
        -0.0000000000000011619070583377173759,
        +0.0000000000000002956508969267836974,
        -0.0000000000000000785294988256492025,
        +0.0000000000000000216922264368256612,
        -0.0000000000000000062113515831676342,
        +0.0000000000000000018384568838450977,
        -0.0000000000000000005610887482137276,
        +0.0000000000000000001761862805280062,
        -0.0000000000000000000568111050541451,
        +0.0000000000000000000187786279582313,
        -0.0000000000000000000063531694151124);

begin
  {Ref: W. Fullerton [20], file d9sifg.f}

  {$ifopt R+}
    if abs(x)<4.0 then begin
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      exit;
    end;
  {$endif}

  if x <= xbnd then begin
    z := x*x;
    f := (1.0 + CSEvalX((1.0/z - 0.04125)/0.02125, f1, n_f1))/x;
    g := (1.0 + CSEvalX((1.0/z - 0.04125)/0.02125, g1, n_g1))/z;
  end
  else if x <= xbig then begin
    z := x*x;
    f := (1.0 + CSEvalX (100.0/z - 1.0, f2, n_f2))/x;
    if x <= xbndg then begin
      g := (1.0 + CSEvalX((10000.0/z - 125.0)/75.0, g2, n_g2))/z;
    end
    else begin
      g := (1.0 + CSEvalX(400.0/z - 1.0, g3, n_g3))/z;
    end;
  end
  else begin
    if x<xmaxf then f := 1.0/x else f := 0.0;
    if x<xmaxg then g := 1.0/sqr(x) else g := 0.0;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ci(x: extended): extended;
  {-Return the cosine integral, ci(x) = EulerGamma + ln(|x|) + integral((cos(t)-1)/t, t=0..|x|)}
var
  f,g,s,c: extended;
const
  n_ci = 14;
  ci_cs: array[0..n_ci-1] of extended = (
           -0.34004281856055363156281076633129873,
           -1.03302166401177456807159271040163751,
           +0.19388222659917082876715874606081709,
           -0.01918260436019865893946346270175301,
           +0.00110789252584784967184098099266118,
           -0.00004157234558247208803840231814601,
           +0.00000109278524300228715295578966285,
           -0.00000002123285954183465219601280329,
           +0.00000000031733482164348544865129873,
           -0.00000000000376141547987683699381798,
           +0.00000000000003622653488483964336956,
           -0.00000000000000028911528493651852433,
           +0.00000000000000000194327860676494420,
           -0.00000000000000000001115183182650184);
          {+0.00000000000000000000005527858887706,
           -0.00000000000000000000000023907013943,
           +0.00000000000000000000000000091001612,
           -0.00000000000000000000000000000307233,
           +0.00000000000000000000000000000000926);}
begin
  {Ref: W. Fullerton [20], file dci.f}
  x := abs(x);
  if x<=4.0 then begin
    if x<=1e-10 then s := -1.0
    else s := (x*x-8.0)*0.125;
    sfc_ci := ln(x) - 0.5 + CSEvalX(s, ci_cs, n_ci)
  end
  else begin
    auxfg(x,f,g);
    sincos(x,s,c);
    sfc_ci := f*s - g*c;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_cin(x: extended): extended;
  {-Return the entire cosine integral, cin(x) = integral((1-cos(t))/t, t=0..x)}
var
  f,g,s,c: extended;
const
  n_cin = 12;
  cin_cs: array[0..n_cin-1] of extended = (
          +0.37074501750909688741654801228564992,
          -0.05893574896364446831956864397363697,
          +0.00538189642113569124048745326203340,
          -0.00029860052841962135319594906563410,
          +0.00001095572575321620077031054467306,
          -0.00000028405454877346630491727187731,
          +0.00000000546973994875384912457861806,
          -0.00000000008124187461318157083277452,
          +0.00000000000095868593117706609013181,
          -0.00000000000000920266004392351031377,
          +0.00000000000000007325887999017895024,
          -0.00000000000000000049143726675842909);
begin
  {Based on W. Fullerton [20], file dcin.f}
  x := abs(x);
  if x<=4.0 then begin
   {Unlike Fullerton we use the first term of the cin series for x^2 < eps_x}
   {cin(x) = 1/4*x^2 - 1/96*x^4 + 1/4320*x^6 - 1/322560*x^8 + O(x^10)}
   s := sqr(x);
   if s<eps_x then sfc_cin := 0.25*s
   else sfc_cin := s*CSEvalX((s-8.0)*0.125, cin_cs, n_cin);
  end
  else begin
    auxfg(x,f,g);
    sincos(x,s,c);
    sfc_cin := (g*c - f*s) + ln(x) + EulerGamma;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_cinh(x: extended): extended;
  {-Return the entire hyperbolic cosine integral, cinh(x) = integral((cosh(t)-1)/t, t=0..x)}
const
  n_cinh = 11;
  cinh_cs: array[0..n_cinh-1] of extended = (
             0.1093291636520734431407425199795917,
             0.0573928847550379676445323429825108,
             0.0028095756978830353416404208940774,
             0.0000828780840721356655731765069792,
             0.0000016278596173914185577726018815,
             0.0000000227809519255856619859083591,
             0.0000000002384484842463059257284002,
             0.0000000000019360829780781957471028,
             0.0000000000000125453698328172559683,
             0.0000000000000000663637449497262300,
             0.0000000000000000002919639263594744);
begin
  {Ref: W. Fullerton [20], file dcinh.f}
  x := abs(x);
  if x<=3.0 then begin
   {Unlike Fullerton we use the first term of the cinh series for x^2 < eps_x}
   {cinh(x) = 1/4*x^2 + 1/96*x^4 + 1/4320*x^6 + 1/322560*x^8 + O(x^10)}
   x := sqr(x);
   if x<eps_x then sfc_cinh := 0.25*x
   else begin
     {Note the bug in [20]: the 9 must be replaced by 4.5!}
     x := x*(0.25 + CSEvalX(x/4.5-1.0, cinh_cs, n_cinh));
     sfc_cinh := x;
   end;
  end
  else sfc_cinh := sfc_chi(x) - (ln(x) + EulerGamma);
end;


{---------------------------------------------------------------------------}
function si_medium(x: extended): extended;
  {-Return si(x), internal function, |x| <= 4}
const
  n_si = 13;
  si_cs: array[0..n_si-1] of extended = (
           -0.1315646598184841928904275173000457,
           -0.2776578526973601892048287660157299,
           +0.0354414054866659179749135464710086,
           -0.0025631631447933977658752788361530,
           +0.0001162365390497009281264921482985,
           -0.0000035904327241606042670004347148,
           +0.0000000802342123705710162308652976,
           -0.0000000013562997692540250649931846,
           +0.0000000000179440721599736775567759,
           -0.0000000000001908387343087145490737,
           +0.0000000000000016669989586824330853,
           -0.0000000000000000121730988368503042,
           +0.0000000000000000000754181866993865);
         { -0.0000000000000000000004014178842446,
           +0.0000000000000000000000018553690716,
           -0.0000000000000000000000000075166966,
           +0.0000000000000000000000000000269113,
           -0.0000000000000000000000000000000858);}
begin
  {Ref: W. Fullerton [20], file dsi.f}
  si_medium := x*(0.75 + CSEvalX(0.125*(x*x-8.0), si_cs, n_si));
end;


{---------------------------------------------------------------------------}
function sfc_si(x: extended): extended;
  {-Return the sine integral, si(x) = integral(sin(t)/t, t=0..x)}
var
  f,g,s,c,z: extended;
begin
  {Ref: W. Fullerton [20], file dsi.f}
  z := abs(x);
  if z<1e-5 then begin
    {si(x) = x*(1 - 1/18*x^2 + 1/600*x^4) + O(x^6)}
    if z<1e-10 then sfc_si := x
    else sfc_si := x*(1.0 - x*x/18.0);
  end
  else if z<=4.0 then sfc_si := si_medium(x)
  else begin
    auxfg(z,f,g);
    sincos(z,s,c);
    z := Pi_2 - f*c - g*s;
    if x>0.0 then sfc_si := z
    else sfc_si := -z;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ssi(x: extended): extended;
  {-Return the shifted sine integral, ssi(x) = si(x) - Pi/2}
var
  f,g,s,c,z: extended;
const
  x0 = 1.926447660317370582023;  {first root of ssi, ssi(x0)=0}
const
  x0hHex: THexExtW = ($0000,$0000,$8000,$F695,$3FFF);  {1.92643737792968750000}
  x0lHex: THexExtW = ($651B,$A01D,$84D1,$AC82,$3FEE);  {1.02823876830820229999E-5}
var
  x0h: extended absolute x0hHex;
  x0l: extended absolute x0lHex;
const
  ssi0: array[0..12] of extended = (
          -0.1351285853611960776140e-1,
          +0.1215796000069655265351e+0,
          -0.6751597084794109525789e-2,
          -0.2392545546343720701774e-4,
          +0.4830474749548895782624e-5,
          -0.3755674460815287806513e-8,
          -0.1708157404013677694186e-8,
          +0.3079320056206561258058e-11,
          +0.3586916508589400202065e-12,
          -0.7339809169354881830878e-15,
          -0.4978973685818888496567e-16,
          +0.1019993223595804689554e-18,
          +0.4902801608220350017913e-20);
begin
  z := abs(x);
  if z<1e-5 then begin
    {ssi(x) = -Pi/2 + x - 1/18*x^3 + 1/600*x^5) + O(x^6)}
    if z<1e-20 then sfc_ssi := -Pi_2
    else begin
      s := 1.0;
      if z>1e-10 then s := s - x*x/18.0;
      sfc_ssi := x*s - Pi_2;
    end;
  end
  else if abs(x-x0)<0.250 then begin
    {x near first root. Since x0<4 the standard algorithm without}
    {f and g for si(x) - Pi/2 is not accurate. Use the expansion }
    {chebyshev(Ssi(x+x0), x=-0.25..0.25) calculated with Maple V }
    z := x-x0h;
    z := z-x0l;
    sfc_ssi := CSEvalX(4.0*z, ssi0, 13);
  end
  else if z<=4.0 then sfc_ssi := si_medium(x) - Pi_2
  else begin
    auxfg(z,f,g);
    sincos(z,s,c);
    z := f*c + g*s;
    if x<0.0 then sfc_ssi := z-Pi
    else sfc_ssi := -z;
  end;
end;


end.
