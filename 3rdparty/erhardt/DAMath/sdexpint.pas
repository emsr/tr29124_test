unit sdExpInt;

{Double precision special functions: Exponential integrals and related}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}

(*************************************************************************

 DESCRIPTION   :  Double precision special functions: Exponential integrals and related

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in damath_info.txt/references

                  [1] [HMF]: M. Abramowitz, I.A. Stegun. Handbook of Mathematical Functions, New York, 1970
                      http://www.math.sfu.ca/~cbm/aands/
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
                 [30] [NIST]: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST Handbook
                      of Mathematical Functions, Cambridge, 2010. Online resource: NIST Digital
                      Library of Mathematical Functions, http://dlmf.nist.gov/

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  10.02.13  W.Ehrhardt  Initial BP7 version from AMath.sfExpInt
 1.00.01  10.02.13  we          sfd_ei
 1.00.02  10.02.13  we          modified sfd_eiex to compute exp(x) if expx < 0
 1.00.03  10.02.13  we          en_series avoids FPC nonsense
 1.00.04  10.02.13  we          sfd_ssi, auxfg
 1.00.05  01.03.13  we          Chebyshev degrees reduced in auxfg, sfd_ci, sfd_cin, sfd_cinh, si_medium, sfd_ssi

 1.06.00  25.09.13  we          use const one_d

 1.17.00  25.01.15  we          sfd_ali from sdMisc
 1.17.01  25.01.15  we          sfd_ei_inv
 1.17.02  01.05.15  we          fix sfd_ei_inv for x Nan

 1.20.00  25.08.15  we          Use THREE & SIXX to avoid 'optimization'

 1.27.00  17.08.17  we          sfd_eisx2

 1.28.00  13.12.17  we          separate functions e1_small, e1ex_large
 1.28.01  13.12.17  we          sfd_e1s(x) = exp(x)*E1(x)
 1.28.01  14.12.17  we          sfd_eis(x) = exp(-x)*Ei(x)
 1.28.02  14.12.17  we          e1_small, e1ex_large with adjusted rational approximations

 1.29.00  14.02.18  we          fixed parameter/return types in ei_asymp

 1.35.00  23.08.18  we          NanOrInf check in sfd_en
 1.35.01  07.09.18  we          sfd_en for n<0 and x>0
 1.35.02  07.09.18  we          sfd_eibeta

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


function sfd_ali(x: double): double;
  {-Return the functional inverse of li(x), li(sfd_ali(x))=x}

function sfd_chi(x: double): double;
  {-Return the hyperbolic cosine integral = EulerGamma + ln(|x|) + integral((cosh(t)-1)/t, t=0..|x|)}

function sfd_ci(x: double): double;
  {-Return the cosine integral, ci(x) = EulerGamma + ln(|x|) + integral((cos(t)-1)/t, t=0..|x|)}

function sfd_cin(x: double): double;
  {-Return the entire cosine integral, cin(x) = integral((1-cos(t))/t, t=0..x)}

function sfd_cinh(x: double): double;
  {-Return the entire hyperbolic cosine integral, cinh(x) = integral((cosh(t)-1)/t, t=0..x)}

function sfd_e1(x: double): double;
  {-Return the exponential integral E1(x) = integral(exp(-x*t)/t, t=1..Inf), x <> 0}

function sfd_e1s(x: double): double;
  {-Return E1s(x) = exp(x)*E1(x), x <> 0}

function sfd_ei(x: double): double;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}

function sfd_ein(x: double): double;
  {-Return the entire exponential integral ein(x) = integral((1-exp(-t))/t, t=0..x)}

function sfd_eis(x: double): double;
  {-Return exp(-x)*Ei(x), x <> 0}

function sfd_eisx2(x: double): double;
  {-Return exp(-x^2)*Ei(x^2), x <> 0}

function sfd_ei_inv(x: double): double;
  {-Return the functional inverse of Ei(x), ei_inv(ei(x))=x}

function sfd_en(n: longint; x: double): double;
  {-Return the exponential integral E_n(x) = integral(exp(-x*t)/t^n, t=1..Inf), x > 0}

function sfd_eibeta(n: integer; x: double): double;
  {-Return the exponential integral beta(n,x) = int(t^n*exp(-x*t), t=-1..1), n >= 0}

function sfd_gei(p,x: double): double;
  {-Return the generalized exponential integral E_p(x) = integral(exp(-x*t)/t^p, t=1..Inf), x >= 0}

function sfd_li(x: double): double;
  {-Return the logarithmic integral li(x) = PV-integral(1/ln(t), t=0..x), x >= 0, x <> 1}

function sfd_shi(x: double): double;
  {-Return the hyperbolic sine integral, shi(x) = integral(sinh(t)/t, t=0..x)}

function sfd_si(x: double): double;
  {-Return the sine integral, si(x) = integral(sin(t)/t, t=0..x)}

function sfd_ssi(x: double): double;
  {-Return the shifted sine integral, ssi(x) = si(x) - Pi/2}


implementation

uses
  DAMath,
  {$ifopt R+}
    sdBasic,  {for RTE_ArgumentRange}
  {$endif}
  sdGamma,
  sdGamma2;


{---------------------------------------------------------------------------}
function e1_small(x: double): double;
  {-Return the exponential integral E1(x) for 0 < x <= 1}
var
  y: double;
const
  P01: array[0..5] of double = (
          0.0865197248079397976498,
          0.0320913665303559189999,
         -0.245088216639761496153,
         -0.0368031736257943745142,
         -0.00399167106081113256961,
         -0.000111507792921197858394);
  Q01: array[0..5] of double = (
          1.0,
          0.37091387659397013215,
          0.056770677104207528384,
          0.00427347600017103698101,
          0.000131049900798434683324,
         -0.528611029520217142048e-6);
const
  y01 = 695977.0/1048576.0;
begin
  {Based on parts of boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
  y :=  PolEval(x,P01,6) / PolEval(x,Q01,6);
  e1_small := y + (x - ln(x) - y01);
end;

{---------------------------------------------------------------------------}
function e1ex_large(x: double): double;
  {-Return the exponential integral exp(x)*E1(x) for x > 1}
var
  y,z: double;
const
  Pix: array[0..10] of double = (
         -0.121013190657725568138e-18,
         -0.999999999999998811143,
         -43.3058660811817946037,
         -724.581482791462469795,
         -6046.8250112711035463,
         -27182.6254466733970467,
         -66598.2652345418633509,
         -86273.1567711649528784,
         -54844.4587226402067411,
         -14751.4895786128450662,
         -1185.45720315201027667);
  Qix: array[0..11] of double = (
         1.0,
         45.3058660811801465927,
         809.193214954550328455,
         7417.37624454689546708,
         38129.5594484818471461,
         113057.05869159631492,
         192104.047790227984431,
         180329.498380501819718,
         86722.3403467334749201,
         18455.4124737722049515,
         1229.20784182403048905,
         -0.776491285282330997549);
begin
  {Based on parts of boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
  y := 1.0/x;
  z := 1.0 + PolEval(y,Pix,11) / PolEval(y,Qix,12);
  e1ex_large := z*y;
end;


{---------------------------------------------------------------------------}
function sfd_e1(x: double): double;
  {-Return the exponential integral E1(x) = integral(exp(-x*t)/t, t=1..Inf), x <> 0}
var
  y: double;
begin
  {Based on boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
  if IsNaND(x) then sfd_e1 := x
  else if x < 0.0 then sfd_e1 := -sfd_ei(-x)
  else if x <= 1.0 then sfd_e1 := e1_small(x)
  else if x>=ln_MaxDbl then sfd_e1 := 0.0
  else begin
    y := e1ex_large(x);
    sfd_e1 := exp(-x)*y;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_e1s(x: double): double;
  {-Return E1s(x) = exp(x)*E1(x), x <> 0}
begin
  {Note that E1s(x) = Re(U(1,1,x))}
  if IsNaND(x) then sfd_e1s := x
  else if x <  0.0 then sfd_e1s := -sfd_eis(-x)
  else if x <= 1.0 then sfd_e1s := exp(x)*e1_small(x)
  else sfd_e1s := e1ex_large(x);
end;


{---------------------------------------------------------------------------}
function ei_asymp(x: double): double;
  {-Return asymptotic expansion for exp(-x)*Ei(x), x -> Inf}
var
  s,t,d,n: double;
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
  until d <= eps_d*s;
  ei_asymp := s/x;
end;


{---------------------------------------------------------------------------}
function sfd_eis(x: double): double;
  {-Return exp(-x)*Ei(x), x <> 0}
const
  xas = 80.0;
begin
  if IsNaND(x) then sfd_eis := x
  else if x < -1.0 then sfd_eis := -sfd_e1s(-x)
  else if x < xas then sfd_eis := exp(-x)*sfd_ei(x)
  else sfd_eis := ei_asymp(x);
end;



{---------------------------------------------------------------------------}
function sfd_eiex(x,expx: double): double;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}
  { internal version for x>6, if expx > 0  it must be = exp(x), used by ei and li}
var
  y,z: double;
const
  {Based on boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
  {Copyright John Maddock 2007, see 3rdparty.ama for Boost license}
  P10: array[0..8] of double = (
         0.00139324086199409049399,
        -0.0345238388952337563247,
        -0.0382065278072592940767,
        -0.0156117003070560727392,
        -0.00383276012430495387102,
        -0.000697070540945496497992,
        -0.877310384591205930343e-4,
        -0.623067256376494930067e-5,
        -0.377246883283337141444e-6);
  Q10: array[0..9] of double = (
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
  P20: array[0..9] of double = (
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
  Q20: array[0..9] of double = (
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
  P40: array[0..11] of double = (
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
  Q40: array[0..8] of double = (
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
  Pix: array[0..8] of double = (
        -0.0130653381347656250004,
         0.644487780349757303739,
         143.995670348227433964,
        -13918.9322758014173709,
         476260.975133624194484,
        -7437102.15135982802122,
         53732298.8764767916542,
        -160695051.957997452509,
         137839271.592778020028);
  Qix: array[0..8] of double = (
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
const
  exp40h: THexDblW = ($72EB,$3979,$220D,$438A); {exp(40)}
begin
  if x<=10.0 then begin
    y := 0.5*x - 4.0;
    z := y10 + PolEval(y,P10,9) / PolEval(y,Q10,10);
  end
  else if x<=20.0 then begin
    y := x/5.0 - 3.0;
    z := y20 + PolEval(y,P20,10) / PolEval(y,Q20,10);
  end
  else if x<=40.0 then begin
    y := x/10.0 - 3.0;
    z := y40 + PolEval(y,P40,12) / PolEval(y,Q40,9);
  end
  else begin
    y := 1.0/x;
    z := yix + PolEval(y,Pix,9) / PolEval(y,Qix,9);
    if (x > 41.0) and (expx<0.0) then begin
      {exp(x) must be computed and possible overflow in Ei}
      y := x-40.0;
      if (y>ln_MaxDbl) or (x>716.35549054245) then sfd_eiex := PosInf_d
      else begin
        z := z*exp(y)/x;
        z := z*double(exp40h);
        sfd_eiex := z + x;
      end;
      exit;
    end;
  end;
  if expx < 0.0 then expx := exp(x);
  z := z*(expx/x);
  sfd_eiex := z + x;
end;


{---------------------------------------------------------------------------}
function sfd_ei(x: double): double;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}
var
  y,z: double;
{Based on boost_1_42_0\boost\math\special_functions\expint.hpp [19]}
{Copyright John Maddock 2007, see 3rdparty.ama for Boost license}
const
  P6H: array[0..9] of THexDblW = (
         ($43CF,$D891,$E4E8,$4007),  {2.98677224343598593013    }
         ($2D59,$7730,$CE55,$3FD6),  {0.356343618769377415068   }
         ($9E64,$F065,$FC9B,$3FE8),  {0.780836076283730801839   }
         ($6AAF,$E660,$5B12,$3FBD),  {0.114670926327032002811   }
         ($9C47,$03CC,$9231,$3FA9),  {0.0499434773576515260534  }
         ($E474,$4C78,$BF04,$3F7D),  {0.00726224593341228159561 }
         ($B5A6,$2177,$EB82,$3F52),  {0.00115478237227804306827 }
         ($E5B6,$2FA2,$84C8,$3F1E),  {0.000116419523609765200999}
         ($166D,$2BAC,$BDD2,$3EE0),  {0.798296365679269702435e-5}
         ($1352,$41D0,$A2F2,$3E92)); {0.2777056254402008721e-6  }
  Q6H: array[0..7] of THexDblW = (
         ($0000,$0000,$0000,$3FF0),  { 1.0                       }
         ($897E,$F65D,$BC05,$BFF2),  {-1.17090412365413911947    }
         ($D216,$6BA1,$E8A9,$3FE3),  { 0.62215109846016746276    }
         ($31F8,$6BFA,$F985,$BFC8),  {-0.195114782069495403315   }
         ($AA51,$A2B0,$0BC6,$3FA4),  { 0.0391523431392967238166  }
         ($A053,$9B98,$AD36,$BF74),  {-0.00504800158663705747345 }
         ($B83D,$56C6,$7EE8,$3F39),  { 0.000389034007436065401822}
         ($0F34,$BA45,$2508,$BEED)); {-0.138972589601781706598e-4}
const
  rh : THexDblW = ($B5FC,$52B4,$D729,$3FD7); { = 0.372507410781366634461991866580119133535689497771.. root of Ei, Ei(r)=0}
  r1h: THexDblW = ($8000,$52B4,$D729,$3FD7); { = 0.372507410780599457211792469024658203125 = 409576229586 / 1099511627776}
  r2 = 7.6717725019939755546093E-13;         { = r - r1}
var
  r : double absolute rh;
  r1: double absolute r1h;
  P06: array[0..9]of double absolute P6H;
  Q06: array[0..7] of double absolute Q6H;
begin
  {avoid infinite loop if x is (negative) NAN}
  if IsNaND(x) then sfd_ei := x
  else if x<0 then sfd_ei := -sfd_e1(-x)
  else if x<=eps_d then begin
    {Boost's rational approximation is suboptimal for small x}
    {use Ei(x) = (+EulerGamma + ln(x)) + x + 1/4*x^2 + O(x^3)}
    {for x <= eps_d, i.e. |ln(x)| > 43, the x term is negligible}
    sfd_ei := EulerGamma + ln(x);
  end
  else if x<=6.0 then begin
    y := x/THREE - 1.0;
    z := PolEval(y,P06,10) / PolEval(y,Q06,8);
    y := (x - r1) - r2;
    z := y*z;
    if abs(y)<0.1 then sfd_ei := z + ln1p(y/r)
    else sfd_ei := z + ln(x/r);
  end
  else sfd_ei := sfd_eiex(x,-1);
end;


{---------------------------------------------------------------------------}
function sfd_ei_inv(x: double): double;
  {-Return the functional inverse of Ei(x), ei_inv(ei(x))=x}
var
  z,d: double;
const
  x0 = -35.5;  { < ln(eps) + Eulergamma}
  x1 = 2.53631e305; {~ max. argument without overflow ~ MaxDouble/ln_MaxDbl}
  z0 = 0.37249755859375 + 0.985218761663446199e-5; {Zero of Ei}
  b  = 3.89621573390716731;  {Taylor expansion Ei(z0+z) = b*(z-z0) + O((z-z0)^2)}
  eg = 0.5614594835668851698; {exp(-Eulergamma)}
begin
  {Basic formula: li(x) = Ei(ln(x)) -> ei_inv(x) = ln(li_inv(x))}
  if IsNanD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ei_inv := NaN_d;
  end
  else if x<=x0 then sfd_ei_inv := eg*exp(x) {exp(x-Eulergamma)}
  else if abs(x) < sqrt_epsh then sfd_ei_inv := z0 + x/b
  else if x>0.0 then begin
    if x>x1 then sfd_ei_inv := PosInf_d
    else sfd_ei_inv := ln(sfd_ali(x))
  end
  else begin
    z := ln(sfd_ali(x));
    if x < -2.0 then begin
      {Here the basic formula is used as starting value for Newton iteration}
      {because numerically it does not always give satisfactory accuracy.   }
      repeat
        d := sfd_ei(z)-x;
        d := d*z/exp(z);
        z := z-d;
      until abs(d) <= sqrt_epsh*abs(z);
    end;
    sfd_ei_inv := z;
  end
end;


{---------------------------------------------------------------------------}
function sfd_ein(x: double): double;
  {-Return the entire exponential integral ein(x) = integral((1-exp(-t))/t, t=0..x)}
const
  csh: array[0..10] of THexDblw = (
         ($F776,$01B7,$071F,$4000),  {+2.003477109362331949501523562122308610572     }
         ($6D62,$5F42,$0802,$BFB0),  {-0.6262221170111595696745652349876399745063e-1 }
         ($92D6,$D981,$7F71,$3F5C),  {+0.1739369565279582463382068320644185874193e-2 }
         ($FBE0,$BA23,$5ED2,$BF05),  {-0.4076080883130112580271691400428456876713e-4 }
         ($9D19,$3928,$59A9,$3EAB),  {+0.8151006219468242274136388199101306590141e-6 }
         ($1B46,$93E7,$629D,$BE4E),  {-0.1414921923800105206557069349069463653541e-7 }
         ($0D7A,$DDCB,$C2FE,$3DED),  {+0.2165448782930998188223426094586382995172e-9 }
         ($89CC,$D393,$09F3,$BD8A),  {-0.2960277555152250090455326359589754482108e-11}
         ($64EF,$82EE,$9276,$3D24),  {+0.3654342755726724875204322769851602486817e-13}
         ($0183,$2EF3,$9F1F,$BCBD),  {-0.4110818570387228515941851406185706129602e-15}
         ($26A8,$0638,$954B,$3C53)); {+0.4246424399152657058898802616196579805661e-17}
      (* ($F907,$FDF3,$EEFC,$BBE7)); {-0.4054500715002447619255009744247776959447e-19}*)
var
  csa: array[0..11] of double absolute csh;
var
  t,z: double;
begin
  {ein(x) = gamma + ln(|x|) - Ei(-x). Ref: NIST[30], 6.2(i)}
  t := abs(x);
  if t<0.25 then begin
    if t<=eps_d then sfd_ein := x
    else begin
      {Chebyshev approximation calculated with Maple from NIST[30], 6.6.4}
      t := CSEvalD(4.0*x, csa, 11);
      sfd_ein := t*x;
    end;
  end
  else if x>0 then begin
    t := ln(t);
    if x >= 45.0 then sfd_ein := t + EulerGamma
    else begin
      z := EulerGamma + sfd_e1(x);
      sfd_ein := z+t;
    end;
  end
  else begin
    z := -sfd_ei(-x);
    if x <= -50 then sfd_ein := z
    else begin
      t := ln(t) + EulerGamma;
      sfd_ein := z+t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_eisx2(x: double): double;
  {-Return exp(-x^2)*Ei(x^2)}
const
  xei  = 10.0;  {double: 10}
  xmax = 1e16;  {double: 1e16}
var
  z,s: double;
begin
  if IsNaND(x) then begin
    sfd_eisx2 := x;
    exit;
  end;
  x := abs(x);
  if x <= xei then begin
    {return exp(-x^2)*Ei(x^2)}
    s := expx2(-x);
    z := sfd_ei(sqr(x));
    sfd_eisx2 := s*z;
  end
  else begin
    {Asymptotic expansion NIST[30], 6.12.2:  1/x^2 * sum(k!/x^(2k)}
    if x >= xmax then begin
      {Only one term of AE, avoid overflow}
      sfd_eisx2 := (1.0/x)/x;
    end
    else begin
      {Asymptotic expansion NIST[30], 6.12.2:  1/x^2 * sum(k!/x^(2k)}
      sfd_eisx2 := ei_asymp(sqr(x));
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_li(x: double): double;
  {-Return the logarithmic integral li(x) = PV-integral(1/ln(t), t=0..x), x >= 0, x <> 1}
var
  t: double;
begin
  {li(x) = Ei(ln(x)). Note that some authors use Li(x) = li(x)-li(2)}
  if x=0.0 then sfd_li := 0.0
  else begin
    t := ln(x);
    if t<6.0 then sfd_li := sfd_ei(t)
    else begin
      {Avoid the inaccuracies of exp(ln(x)) for large x}
      {and use the reformulated Ei code for ln(x) > 6}
      sfd_li := sfd_eiex(t,x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ali(x: double): double;
  {-Return the functional inverse of li(x), li(sfd_ali(x))=x}
var
  d,z: double;
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
      d := sfd_li(z) - x;
      if d<>0.0 then begin
        d := d*ln(z)/(1.0 + 0.5*d/z);  {Newton if /(1+..) is omitted}
        z := z - d;
      end;
      {If |d| < z*sqrt_epsh the correction from the next iteration step}
      {can be neglected because the convergence is at least quadratic. }
    until (abs(d) < z*sqrt_epsh) or (k>4);
  end;
  sfd_ali := z;
end;



{---------------------------------------------------------------------------}
function en_series(n: longint; x: double): double;
  {-Calculate E_n(x), 0<x<=1, n>=1, via series expansion}
var
  m,k: longint;
  f,p,s,z: double;
begin

  {HMF [1], formula 5.1.12}
  m := 0;
  k := 1-n;
  f := 1;
  z := -x;

  {initialize sum s with m=0 term}
  if n=1 then s := 0.0
  else s := one_d/k;

  {sum until first term with relative contribution is <= eps_d}
  repeat
    inc(m);
    inc(k);
    f := f*z/m;
    {skip x^(n-1) term}
    if k<>0 then s := s + f/k;
  until abs(f) < eps_d*abs(s);

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
        p := p + one_d/k;
      end;
      p := p - (EulerGamma + ln(x));
      en_series := p*f - s;
    end
    else begin
      {since |z| <= 1, z^(n-1) cannot overflow}
      if n > MAXGAMD then en_series := -s
      else begin
        {en_series := (psi(n)-ln(x))*z^(n-1)/gamma(n) - s;}
        p := sfd_psi(n)-ln(x);
        f := intpower(z,n-1)/sfd_gamma(n);;
        en_series := p*f - s;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function en_cfrac(n: longint; x: double): double;
  {-Calculate E_n(x), x>1, n>1, via continued fraction}
var
  k: longint;
  big,cf,pk,p1,p2,qk,q1,q2,r,t: double;
  ev: boolean;
begin
  {Algorithm based on Cephes [7] function expn}

  big:= ldexpd(1.0,68); {Rescaling threshold, should be power of 2}

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
  until t<=eps_d;
  en_cfrac := cf*exp(-x);
end;


{---------------------------------------------------------------------------}
function ep_largep(p,x: double): double;
  {-Uniform asymptotic expansion E_p(x) for large p; called with p>=10000, 0 <= x <= 740}
var
  s,z,a2,a3,a4: double;
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
function sfd_en(n: longint; x: double): double;
  {-Return the exponential integral E_n(x) = integral(exp(-x*t)/t^n, t=1..Inf), x > 0}
begin
  if IsNaNorInfD(x) then sfd_en := NaN_d
  else if n<0       then sfd_en := sfd_gei(n,x)
  else if n=0       then sfd_en := exp(-x)/x {HMF [1], 5.1.24}
  else if n=1       then sfd_en := sfd_e1(x)
  else if x=0.0     then sfd_en := one_d/(n-1) {HMF [1], 5.1.23}
  else if x<0.0     then sfd_en := ln(x)       {force error}
  else if x>740     then sfd_en := 0.0         {HMF [1], 5.1.19/20: 0 < en <= exp(-x)/x}
  else if n>=10000  then sfd_en := ep_largep(n,x)
  else if x<=1.0    then sfd_en := en_series(n,x)
  else                   sfd_en := en_cfrac(n,x);
end;


{---------------------------------------------------------------------------}
function sfd_gei(p,x: double): double;
  {-Return the generalized exponential integral E_p(x) = integral(exp(-x*t)/t^p, t=1..Inf), x >= 0}
var
  a,q,t: double;
begin
  if (x<0.0) or IsNanOrInfD(x) or IsNanOrInfD(p) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_gei := NaN_d;
    exit;
  end;
  {E_p(x) = x^(p-1)*GAMMA(1-p,x), if x>0}
  if (p>=0.0) and (frac(p)=0.0) and (p<MaxLongInt) then sfd_gei := sfd_en(round(p),x)
  else begin
    if x=0.0 then begin
      if p>1.0 then sfd_gei := 1.0/(p-1.0)
      else sfd_gei := PosInf_d;
    end
    else if (x>740) and (p>=0.0) then sfd_gei := 0.0
    else if p>=10000.0 then sfd_gei := ep_largep(p,x)
    else if ((x>=1.0) and (p>=-0.25)) or ((p>1.0) and (x>0.25)) then begin
      {sfd_igamma will use igam_qfraction for these (p,x) values, so}
      {compute CF directly without scaling & rescaling with x^(p-1)!}
      t := exp(-x);
      if t<>0 then igam_qfraction(1.0-p,x,t,a,t);
      sfd_gei := t;
    end
    else begin
      a := 1.0 - p;
      if a >= MAXGAMD then begin
        {Gamma(a) will overflow, compute Q(a,x)*Gamma(a)/x^a using logarithmic forms}
        sfd_incgamma_ex(a,x,p,q,t,false);
        if q=0.0 then sfd_gei := 0.0
        else begin
          t := sfd_lngamma(a) - a*ln(x);
          if q<0.75 then t := t + ln(q)
          else t := t + ln1p(-p);
          if t>ln_MaxDbl then sfd_gei := PosInf_d
          else sfd_gei := exp(t);
        end;
      end
      else begin
        t := power(x,-a);
        a := sfd_igamma(a,x);
        {use PosInf_d if result would overflow}
        if (a<=1.0) or (t<Maxdouble/a) then sfd_gei := a*t
        else sfd_gei := PosInf_d;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function eib_small(n: integer; x: double): double;
  {-Exponential integral beta(n,x) for 'small' x}
var
  k,q,m: longint;
  s,t,z: double;
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
        if abs(t)<=eps_d*abs(s) then break;
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
function sfd_eibeta(n: integer; x: double): double;
  {-Return the exponential integral beta(n,x) = int(t^n*exp(-x*t), t=-1..1), n >= 0}
var
  f,g,y: double;
begin
  {HMF[1], 5.1.46}
  if n<=0 then begin
    if n<>0 then sfd_eibeta := Nan_d
    else sfd_eibeta := 2.0*sinhc(x);
  end
  else if abs(x) <= 2.0 then begin
    sfd_eibeta := eib_small(n,x);
  end
  else begin
    n := n+1;
    f := sfd_igamma(n,-x);
    g := sfd_igamma(n,x);
    y := power(x, -n);
    sfd_eibeta := (f-g)*y;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_chi(x: double): double;
  {-Return the hyperbolic cosine integral = EulerGamma + ln(|x|) + integral((cosh(t)-1)/t, t=0..|x|)}
begin
  x := abs(x);
  {e1(x)/ei(x) < 0.3e-16 for x>19: skip calculation of e1 in this case}
  if x>23.0 then sfd_chi := 0.5*sfd_ei(x)
  else sfd_chi := 0.5*(sfd_ei(x) - sfd_e1(x));
end;


{---------------------------------------------------------------------------}
function sfd_shi(x: double): double;
  {-Return the hyperbolic sine integral, shi(x) = integral(sinh(t)/t, t=0..x)}
var
  z: double;
const
  ShiCS : array[0..6] of double = (
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
  if z<=1e-4 then begin
    {shi(x) = x + 1/18*x^3 + 1/600*x^5 + O(x^6)}
    if z<=1e-8 then sfd_shi := x
    else sfd_shi := x*(1.0 + x*x/18.0);
  end
  else begin
    if z<=0.375 then z := z*(1.0 + CSEvalD(128.0*sqr(z)/9.0-1.0, ShiCS, 7))
    else if z>19.0 then z := 0.5*sfd_ei(z)
    else z := 0.5*(sfd_ei(z)+sfd_e1(z));
    {shi(-x) = -shi(x)}
    if x>0 then sfd_shi := z
    else sfd_shi := -z;
  end;
end;


{---------------------------------------------------------------------------}
procedure auxfg(x: double; var f,g: double);
  {-Return the auxiliary functions f,g for ci,si, x>=4}
var
  z: double;
const
  xbig  = 94906265.6242515528891;    {sqrt(2/eps_d)}
  xmaxf = 4.49423283715566653E307;   {exp(min(-ln(Mindouble), ln(Maxdouble) - 0.01)}
  xmaxg = 6.70390396497129855E153;   {1/sqrt(Mindouble)}
  xbnd  = 7.07106781186547524;       {sqrt(50)}
  xbndg = 14.1421356237309505;       {sqrt(200)}
const
  n_f1 = 22;
  n_f2 = 34;
  n_g1 = 23;
  n_g2 = 22;
  n_g3 = 21;
const
  f1: array[0..n_f1-1] of double = (
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
        -0.0000000000000000012046901573375820);
     {  +0.0000000000000000002525052443655940,
        -0.0000000000000000000533980268805594,
        +0.0000000000000000000113855786274122);}
const
  f2: array[0..n_f2-1] of double = (
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
        -0.00000000000000000133690399951301570);
      { +0.00000000000000000068941057702931664,
        -0.00000000000000000035905362610437250,
        +0.00000000000000000018878077255791706,
        -0.00000000000000000010016125265594380,
        +0.00000000000000000005360725691578228,
        -0.00000000000000000002893198974944827,
        +0.00000000000000000001574065100202625,
        -0.00000000000000000000863027106431206,
        +0.00000000000000000000476715602862288,
        -0.00000000000000000000265222739998504);}
const
  g1: array[0..n_g1-1] of double = (
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
        +0.0000000000000000024771654107129795);
      { -0.0000000000000000005330672753223389,
        +0.0000000000000000001155666075598465,
        -0.0000000000000000000252280547744957,
        +0.0000000000000000000055429038550786);}

const
  g2: array[0..n_g2-1] of double = (
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
        -0.0000000000000000010175807410227657);
      { +0.0000000000000000002267657399884672,
        -0.0000000000000000000511630776076426,
        +0.0000000000000000000116767014913108);}
const
  g3: array[0..n_g3-1] of double = (
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
        +0.0000000000000000018384568838450977);
      { -0.0000000000000000005610887482137276,
        +0.0000000000000000001761862805280062,
        -0.0000000000000000000568111050541451,
        +0.0000000000000000000187786279582313,
        -0.0000000000000000000063531694151124);}
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
    f := (1.0 + CSEvalD((1.0/z - 0.04125)/0.02125, f1, n_f1))/x;
    g := (1.0 + CSEvalD((1.0/z - 0.04125)/0.02125, g1, n_g1))/z;
  end
  else if x <= xbig then begin
    z := x*x;
    f := (1.0 + CSEvalD (100.0/z - 1.0, f2, n_f2))/x;
    if x <= xbndg then begin
      g := (1.0 + CSEvalD((10000.0/z - 125.0)/75.0, g2, n_g2))/z;
    end
    else begin
      g := (1.0 + CSEvalD(400.0/z - 1.0, g3, n_g3))/z;
    end;
  end
  else begin
    if x<xmaxf then f := 1.0/x else f := 0.0;
    if x<xmaxg then g := 1.0/sqr(x) else g := 0.0;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ci(x: double): double;
  {-Return the cosine integral, ci(x) = EulerGamma + ln(|x|) + integral((cos(t)-1)/t, t=0..|x|)}
var
  f,g,s,c: double;
const
  n_ci = 13;
  ci_cs: array[0..n_ci-1] of double = (
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
           +0.00000000000000000194327860676494420);
         { -0.00000000000000000001115183182650184,
           +0.00000000000000000000005527858887706,
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
    sfd_ci := ln(x) - 0.5 + CSEvalD(s, ci_cs, n_ci)
  end
  else begin
    auxfg(x,f,g);
    sincos(x,s,c);
    sfd_ci := f*s - g*c;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_cin(x: double): double;
  {-Return the entire cosine integral, cin(x) = integral((1-cos(t))/t, t=0..x)}
var
  f,g,s,c: double;
const
  n_cin = 11;
  cin_cs: array[0..n_cin-1] of double = (
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
          +0.00000000000000007325887999017895024);
        { -0.00000000000000000049143726675842909);}
begin
  {Based on W. Fullerton [20], file dcin.f}
  x := abs(x);
  if x<=4.0 then begin
   {Unlike Fullerton we use the first term of the cin series for x^2 < eps_d}
   {cin(x) = 1/4*x^2 - 1/96*x^4 + 1/4320*x^6 - 1/322560*x^8 + O(x^10)}
   s := sqr(x);
   if s<eps_d then sfd_cin := 0.25*s
   else sfd_cin := s*CSEvalD((s-8.0)*0.125, cin_cs, n_cin);
  end
  else begin
    auxfg(x,f,g);
    sincos(x,s,c);
    sfd_cin := (g*c - f*s) + ln(x) + EulerGamma;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_cinh(x: double): double;
  {-Return the entire hyperbolic cosine integral, cinh(x) = integral((cosh(t)-1)/t, t=0..x)}
const
  n_cinh = 10;
  cinh_cs: array[0..n_cinh-1] of double = (
             0.1093291636520734431407425199795917,
             0.0573928847550379676445323429825108,
             0.0028095756978830353416404208940774,
             0.0000828780840721356655731765069792,
             0.0000016278596173914185577726018815,
             0.0000000227809519255856619859083591,
             0.0000000002384484842463059257284002,
             0.0000000000019360829780781957471028,
             0.0000000000000125453698328172559683,
             0.0000000000000000663637449497262300);
           { 0.0000000000000000002919639263594744}
begin
  {Ref: W. Fullerton [20], file dcinh.f}
  x := abs(x);
  if x<=3.0 then begin
   {Unlike Fullerton we use the first term of the cinh series for x^2 < eps_d}
   {cinh(x) = 1/4*x^2 + 1/96*x^4 + 1/4320*x^6 + 1/322560*x^8 + O(x^10)}
   x := sqr(x);
   if x<eps_d then sfd_cinh := 0.25*x
   else begin
     {Note the bug in [20]: the 9 must be replaced by 4.5!}
     x := x*(0.25 + CSEvalD(x/4.5-1.0, cinh_cs, n_cinh));
     sfd_cinh := x;
   end;
  end
  else sfd_cinh := sfd_chi(x) - (ln(x) + EulerGamma);
end;


{---------------------------------------------------------------------------}
function si_medium(x: double): double;
  {-Return si(x), internal function, |x| <= 4}
const
  n_si = 12;
  si_cs: array[0..n_si-1] of double = (
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
           -0.0000000000000000121730988368503042);
        {  +0.0000000000000000000754181866993865,
           -0.0000000000000000000004014178842446,
           +0.0000000000000000000000018553690716,
           -0.0000000000000000000000000075166966,
           +0.0000000000000000000000000000269113,
           -0.0000000000000000000000000000000858);}
begin
  {Ref: W. Fullerton [20], file dsi.f}
  si_medium := x*(0.75 + CSEvalD(0.125*(x*x-8.0), si_cs, n_si));
end;


{---------------------------------------------------------------------------}
function sfd_si(x: double): double;
  {-Return the sine integral, si(x) = integral(sin(t)/t, t=0..x)}
var
  f,g,s,c,z: double;
begin
  {Ref: W. Fullerton [20], file dsi.f}
  z := abs(x);
  if z<=1e-4 then begin
    {si(x) = x*(1 - 1/18*x^2 + 1/600*x^4) + O(x^6)}
    if z<=1e-8 then sfd_si := x
    else sfd_si := x*(1.0 - x*x/18.0);
  end
  else if z<=4.0 then sfd_si := si_medium(x)
  else begin
    auxfg(z,f,g);
    sincos(z,s,c);
    z := Pi_2 - f*c - g*s;
    if x>0.0 then sfd_si := z
    else sfd_si := -z;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ssi(x: double): double;
  {-Return the shifted sine integral, ssi(x) = si(x) - Pi/2}
var
  f,g,s,c,z: double;
const
  x0 = 1.926447660317370582;  {first root of ssi, ssi(x0)=0}
const
  x0hHex: THexDblW = ($4000,$C828,$D2BA,$3FFE);  {1.92644766031662584282457828521728515625}
  x0lHex: THexDblW = ($6F35,$AC93,$3403,$3D6A);  {7.44739198316135E-13}
var
  x0h: double absolute x0hHex;
  x0l: double absolute x0lHex;
const
  ssi0: array[0..11] of double = (
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
          +0.1019993223595804689554e-18);
       {  +0.4902801608220350017913e-20);}
begin
  z := abs(x);
  if z<1e-4 then begin
    {ssi(x) = -Pi/2 + x - 1/18*x^3 + 1/600*x^5) + O(x^6)}
    if z<1e-16 then sfd_ssi := -Pi_2
    else begin
      s := 1.0;
      if z>1e-8 then s := s - x*x/18.0;
      sfd_ssi := x*s - Pi_2;
    end;
  end
  else if abs(x-x0)<0.250 then begin
    {x near first root. Since x0<4 the standard algorithm without}
    {f and g for si(x) - Pi/2 is not accurate. Use the expansion }
    {chebyshev(Ssi(x+x0), x=-0.25..0.25) calculated with Maple V }
    z := x-x0h;
    z := z-x0l;
    sfd_ssi := CSEvalD(4.0*z, ssi0, 12);
  end
  else if z<=4.0 then sfd_ssi := si_medium(x) - Pi_2
  else begin
    auxfg(z,f,g);
    sincos(z,s,c);
    z := f*c + g*s;
    if x<0.0 then sfd_ssi := z-Pi
    else sfd_ssi := -z;
  end;
end;


end.
