unit sfZeta2;

{Common code for Zeta functions, part 2: less common functions}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}

{$ifdef NOBASM}
  {$undef BASM}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Common code for Zeta functions, part 2: less common functions}

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  The unit can be compiled with TP5-6 and D1, but there may be
                  some restrictions with linking large EXEs and a few functions
                  may generate invalid operations due to TP5's brain-damaged
                  usage of the FPU.

 REFERENCES    :  References used in this unit, main index in amath_info.txt/references

                 [22] A.J. MacLeod, MISCFUN: A software package to compute uncommon special functions.
                      ACM Trans. on Math. Soft. 22 (1996), pp.288-301.
                 [33] http://functions.wolfram.com/: Formulas and graphics about
                      mathematical functions for the mathematical and scientific
                      community and/or http://mathworld.wolfram.com/ ("/the web's
                      most extensive mathematical resource/")
                 [47] M. Goano, Algorithm 745: Computation of the complete and incomplete
                      Fermi-Dirac integral. ACM TOMS, Vol.21, No.3, 1995, pp.221-232.
                      Fortran source available from http://netlib.org/toms/745

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.17.00  15.06.14  W.Ehrhardt  New unit with Fermi/Dirac, Lobachewsky, harmonic functions

 1.35.00  25.08.15  we          Use THREE & SIXX to avoid 'optimization'

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

function sfc_fermi_dirac(n: integer; x: extended): extended;
  {-Return the integer order Fermi-Dirac integral F_n(x) = 1/n!*integral(t^n/(exp(t-x)+1), t=0..INF)}

function sfc_fdm05(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}

function sfc_fdp05(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}

function sfc_fdp15(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}

function sfc_fdp25(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}

function sfc_harmonic(x: extended): extended;
  {-Return the harmonic number function H(x) = psi(x+1) + EulerGamma}

function sfc_harmonic2(x,r: extended): extended;
  {-Return the generalized harmonic function H(x,r) = zeta(r)-zetah(r,x+1); x >= -1}

function sfc_llci(x: extended): extended;
  {-Return the Lobachevski function L(x) = integral(-ln(|cos(t)|), t=0..x)}

function sfc_llsi(x: extended): extended;
  {-Return the Lobachevski function Lambda(x) = integral(-ln(|2sin(t)|), t=0..x)}


implementation

uses
  AMath,
  sfBasic, sfGamma, sfZeta;

{---------------------------------------------------------------------------}
function fd_asymp_exp_int(n: integer; x: extended): extended;
  {-Fermi-Dirac asymptotic expansion for integer order, n>0, x > ln_MaxExt}
var
  p,s,t,f,z: extended;
  k: integer;
begin
  {Asymptotic expansion Goano [47], (9)}
  if n > x - 10 then begin
    {The AE is valid only if x >> n. The loop below converges e.g. for}
    {n = trunc(ln_MaxExt), x = ln_MaxExt+10, but the result is Inf. So}
    {return PosInf_x in this invalid situation for all these cases.}
    fd_asymp_exp_int := PosInf_x;
    exit;
  end;
  s := 1.0;
  p := 1.0;
  k := 0;
  z := sqr(1.0/x);
  f := 2.0*n*(n+1);
  repeat
    {loop entry invariants: p,f,s,t >= 0;  p = pochhammer(1-n,k)}
    f := f*z;
    inc(k,2);
    t := p*sfc_etaint(k)*f;
    s := s+t;
    p := p*(k-n-1)*(k-n);
  until t < eps_x*s;
  f := sfc_taylor(x,n+1);
  fd_asymp_exp_int := f*s;
end;


{---------------------------------------------------------------------------}
function fd_polylog_invf(n: integer; x,z: extended): extended;
  {-Return polylog(n,x), x<-1, n>2, z = ln(-x)}
var
  p,s,t: extended;
  j,k: integer;
begin
  {$ifdef debug}
    if (x >= -1) or (n <2) then RunError(byte(RTE_ArgumentRange));
  {$endif}
  {This is part of the standard polylog code, the only difference is,}
  {that z = ln(-x) is externally supplied to increase accuracy}
  s := 0.0;
  for k:=1 to n div 2 do begin
    j := 2*k;
    p := -sfc_etaint(j);
    j := n-j;
    if j>0 then p := sfc_taylor(z,j)*p;
    s := s + p;
  end;
  t := sfc_polylog(n, 1.0/x);
  if odd(n) then t := -t;
  fd_polylog_invf := 2.0*s - t - sfc_taylor(z,n);
end;


{---------------------------------------------------------------------------}
function sfc_fermi_dirac(n: integer; x: extended): extended;
  {-Return the integer order Fermi-Dirac integral F_n(x) = 1/n!*integral(t^n/(exp(t-x)+1), t=0..INF)}
var
  z: extended;
begin
  {Basic formula is F_n(x) = -polylog(n+1,-exp(x)}
  if x < -ln_MaxExt then begin
    sfc_fermi_dirac := 0.0;
  end
  else if x > ln_MaxExt then begin
    {asymptotic Goano [47], (9)}
    if n>0 then sfc_fermi_dirac := fd_asymp_exp_int(n, x)
    else case n of
        0: sfc_fermi_dirac := x;     {ln(1+e^x)   -> x}
       -1: sfc_fermi_dirac := 1.0;   {e^x/(1+e^x) -> 1}
      else sfc_fermi_dirac := 0.0;
    end;
  end
  else begin
    z := -exp(x);
    n := n+1;
    if z<-1.0 then begin
      {case x>0}
      if n<0 then begin
        {Li(-n,z) = (-1)^(n-1) * Li(-n,1/z}
        z := sfc_polylog(n, 1.0/z);
        if odd(n) then z := -z;
        sfc_fermi_dirac := z;
      end
      else if n>2 then sfc_fermi_dirac := -fd_polylog_invf(n, z, x)
      else begin
        {n=0,1,2: sfc_polylog can handle full range}
        sfc_fermi_dirac := -sfc_polylog(n, z);
      end;
    end
    else begin
      { -1 <= z < 0}
      sfc_fermi_dirac := -sfc_polylog(n, z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function fd_asymp_exp_real(q, x: extended): extended;
  {-Fermi-Dirac asymptotic expansion for real order q with q > -1, x >> q}
var
  p,s,t,r,f,z: extended;
  k: integer;
begin
  {Asymptotic expansion Goano [47], (9)}
  {$ifdef debug}
    if q + 10 > x then begin
      {The AE is valid only if x >> q. Normally used only}
      {for q=-1/2, 1/2, (3/2, 5/2) and x > ~40}
      fd_asymp_exp_real := PosInf_x;
      exit;
    end;
  {$endif}
  s := 1.0;
  p := 1.0;
  k := 0;
  z := sqr(1.0/x);
  f := 2.0*q*(q+1.0);
  t := f*1e20;
  repeat
    {p = pochhammer(1-n,k)}
    f := f*z;
    inc(k,2);
    r := p*sfc_etaint(k)*f;
    if abs(r) <= abs(t) then begin
      {terms of AE are still decreasing}
      t := r;
      s := s+t;
      p := p*((k-1)-q)*(k-q);
    end
    else t := 0.0;
  until abs(t) <= eps_x*abs(s);
  f := power(x,q+1.0)*sfc_rgamma(q+2.0);
  fd_asymp_exp_real := f*s;
end;


{---------------------------------------------------------------------------}
function sfc_fdp05(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}
var
  t: extended;

{Note: THexExtW are needed with 16-bit compilers for relative errors}
{less than 5 eps_x, they are not needed for the 32+ bit compilers.}
const
  { t05 := t -> fd(1/2, 12*(t+2))/(12*(t+2))^(3/2);}
  { x=12*(t+2), x=12..36, t=x/12-2                 }
  { chebyshev(t05(t),t=-1..1,1e-20);               }
  np05_2 = 35;
  fdp05x2: array[0..np05_2 - 1] of THexExtW = (
             ($0325,$93D4,$D6AE,$C136,$3FFF),  {2*0.754743020650337869444270498780,    }
             ($1332,$AD56,$7AC6,$A38D,$BFF6),  { -0.249561545813512079811039605579e-2, }
             ($31E0,$2E66,$52CB,$80BA,$3FF5),  { +0.982115370416634900398147501155e-3, }
             ($A8A0,$B79D,$9E81,$B671,$BFF3),  { -0.347983979931683442238436985807e-3, }
             ($6CAC,$2A13,$2F9A,$F3D6,$3FF1),  { +0.116270381856917279446003317185e-3, }
             ($F249,$CDCD,$408E,$9CFC,$BFF0),  { -0.374282264886563175294182249759e-4, }
             ($547E,$54BC,$A0FE,$C50C,$3FEE),  { +0.117450554111933012129375923807e-4, }
             ($11DD,$D90C,$0D37,$F2DA,$BFEC),  { -0.361877328445089611740599641640e-5, }
             ($B3E1,$433A,$A033,$93A2,$3FEB),  { +0.109996838094100301284634122410e-5, }
             ($3117,$E381,$EF27,$B1AA,$BFE9),  { -0.330931901425378656153557493753e-6, }
             ($F14C,$5DFF,$6DFD,$D41D,$3FE7),  { +0.987737251258776378533216391827e-7, }
             ($4DFD,$B405,$59FB,$FBA2,$BFE5),  { -0.292940746917047935446228969668e-7, }
             ($91D1,$55E6,$C7EC,$9475,$3FE4),  { +0.864151410267920646768938404232e-8, }
             ($3CFE,$C040,$5C3C,$AE54,$BFE2),  { -0.253682857771686895946534432944e-8, }
             ($FB73,$7874,$C1CD,$CBBB,$3FE0),  { +0.741177885877868223911943372712e-9, }
             ($6DCD,$715E,$0C92,$ECE9,$BFDE),  { -0.215468706377707924304414487698e-9, }
             ($494B,$801D,$C632,$88FF,$3FDD),  { +0.622999859844165129697643004086e-10,}
             ($38FE,$2B65,$BC9B,$9D81,$BFDB),  { -0.179064482098473796784605475808e-10,}
             ($23F5,$0871,$71CC,$B3EC,$3FD9),  { +0.511373660300697402133999999999e-11,}
             ($11E2,$C31F,$D8E2,$CC23,$BFD7),  { -0.145050214090160134093432434719e-11,}
             ($1915,$6061,$2A01,$E602,$3FD5),  { +0.408577089407887189170064356955e-12,}
             ($D795,$E7C5,$F123,$80AE,$BFD4),  { -0.114293789551114639267975259137e-12,}
             ($3BA0,$F4DF,$2AFC,$8F08,$3FD2),  { +0.317594630475703386901919669155e-13,}
             ($7CE1,$4717,$D5D7,$9DFD,$BFD0),  { -0.877029250354717665156906845619e-14,}
             ($289F,$725E,$A63D,$AD87,$3FCE),  { +0.240821085811104890886594570029e-14,}
             ($EBF1,$CB47,$2005,$BDA5,$BFCC),  { -0.657963336164063747675726341907e-15,}
             ($FE35,$40C0,$B51C,$CE5D,$3FCA),  { +0.178994011257686697681940768666e-15,}
             ($FA95,$B066,$076E,$DFC0,$BFC8),  { -0.485180718104061909634313862435e-16,}
             ($A0FE,$B3A4,$A3F0,$F1E0,$3FC6),  { +0.131122056316856541233113653119e-16,}
             ($7B28,$4D87,$4BC6,$826C,$BFC5),  { -0.353512334501432389833325039827e-17,}
             ($D3DF,$68EE,$19DC,$8C62,$3FC3),  { +0.951273613288044223998852213460e-18,}
             ($D286,$C4CB,$1F59,$96E1,$BFC1),  { -0.255599619674429801332567738120e-18,}
             ($81DB,$D217,$D5FA,$A1F9,$3FBF),  { +0.685994710108645567896702114342e-19,}
             ($30D7,$33A1,$EDF0,$ADBD,$BFBD),  { -0.183956404622969457008106663899e-19,}
             ($2063,$51A5,$A578,$BA40,$3FBB)); { +0.493006328361511158404496844422e-20)}
var
  fdp05_2: array[0..np05_2 - 1] of extended absolute fdp05x2;

const
  { t01 := t -> fd(1/2, 4*(t+2));    }
  { x = 4*(t+2), x=4..12, t = x/4-2  }
  { chebyshev(t01(t),t=-1..1,1e-20); }
  np05_1 = 29;
  fdp05x1: array[0..np05_1 - 1] of THexExtW = (
             ($9036,$2A89,$1E59,$917C,$4004),  {2*18.1856047597986575196785666713,      }
             ($CE1E,$584A,$A0C0,$C8E6,$4002),  { +12.5563056481811830249567824342,      }
             ($48DB,$74EC,$C332,$D5D7,$3FFE),  {  +0.835323524302041720793490203916,    }
             ($350E,$7608,$78A3,$A8F2,$BFFA),  {  -0.412468635980776559531897297026e-1, }
             ($8CAF,$FBFB,$F802,$9ED8,$3FF7),  {  +0.484764203972064404543231999918e-2, }
             ($9E82,$FEC5,$6318,$BEEF,$BFF4),  {  -0.728359626235761181607282375542e-3, }
             ($D277,$BBEA,$8373,$E53F,$3FF1),  {  +0.109314012312724226246017534127e-3, }
             ($D011,$0B98,$04F1,$DAEC,$BFEE),  {  -0.130487650872812128272229904414e-4, }
             ($95AE,$67EC,$7A01,$D140,$3FE9),  {  +0.389761965081220490099761304103e-6, }
             ($CB0D,$E245,$01BA,$E0BB,$3FE9),  {  +0.418593166671574474687704852229e-6, }
             ($3EEB,$A5AB,$8D3E,$BF62,$BFE8),  {  -0.178241140882781326674081586696e-6, }
             ($6752,$1E45,$A71B,$D759,$3FE6),  {  +0.501401270935653913653534971074e-7, }
             ($41C0,$19D3,$3F38,$C495,$BFE4),  {  -0.114426363665601208334168909836e-7, }
             ($D0DD,$F32C,$7F70,$96A1,$3FE2),  {  +0.219196737232486715284670447279e-8, }
             ($C520,$0180,$2F58,$B8C4,$BFDF),  {  -0.336088028068196988951049728505e-9, }
             ($63D6,$88E9,$6123,$87F0,$3FDC),  {  +0.309089460235051920170794233970e-10,}
             ($2A7F,$C55D,$0F79,$FD9B,$3FD8),  {  +0.360395382713747844307256979851e-11,}
             ($AB3A,$0FBA,$168C,$CF46,$BFD8),  {  -0.294553759575921784436627612442e-11,}
             ($AE44,$1830,$9490,$8F75,$3FD7),  {  +0.101933962170251609766666666667e-11,}
             ($05A1,$DC32,$D55A,$9820,$BFD5),  {  -0.270234067190522926130300324494e-12,}
             ($5190,$1535,$5E00,$86C3,$3FD3),  {  +0.598468621686930875020197021442e-13,}
             ($94D4,$AA55,$68D0,$C665,$BFD0),  {  -0.110131976095155602873469185472e-13,}
             ($3C2E,$FE60,$C044,$D9AC,$3FCD),  {  +0.151042240381539761772702788261e-14,}
             ($782A,$C38A,$7E28,$9173,$BFC9),  {  -0.630793784338055501351350020007e-16,}
             ($7617,$CA85,$75E8,$ED02,$BFC8),  {  -0.513932671709885284753382122899e-16,}
             ($7515,$CA08,$2EE7,$E368,$3FC7),  {  +0.246555126271970717149547742394e-16,}
             ($4387,$E419,$6D72,$8C9B,$BFC6),  {  -0.762232826802428016612044204458e-17,}
             ($57FC,$2AA3,$B920,$8C61,$3FC4),  {  +0.190252722217962867356929160937e-17,}
             ($CAB7,$C88B,$F3B5,$EB11,$BFC1)); {  -0.398224281531035871813174429148e-18)}
var
  fdp05_1: array[0..np05_1 - 1] of extended absolute fdp05x1;

const
  { t00 := t -> fd(1/2, 2*(t+1));    }
  { x = 2*(t+1), x=0..4, t = x/2 - 1 }
  { chebyshev(t00(t),t=-1..1,1e-20); }
  np05_0 = 31;
  fdp05x0: array[0..np05_0 - 1] of THexExtW = (
             ($5CAE,$FEB4,$A28E,$CEF8,$4001),  {2*3.23392547573422976421866059174,     }
             ($9337,$80EE,$7BC1,$B8A7,$4000),  { +2.88522237679328803577304079965,     }
             ($2700,$C1B3,$1A93,$D096,$3FFD),  { +0.407395201241206365841823162781,    }
             ($65B1,$69A3,$AD8A,$D352,$BFF8),  { -0.128981299145840220950575776539e-1, }
             ($A854,$62A5,$BBDD,$BDA4,$BFF6),  { -0.289373003523482132846797286683e-2, }
             ($4C64,$8FBE,$767C,$EDCC,$3FF4),  { +0.907129985662663210233804487629e-3, }
             ($9CA9,$41AF,$63E6,$A169,$BFF1),  { -0.769670870754619578482136591047e-4, }
             ($1225,$B186,$0662,$AC14,$BFEF),  { -0.205133226436761049829261853507e-4, }
             ($1C19,$1377,$110A,$FD78,$3FED),  { +0.755396515130974716466786294881e-5, }
             ($6962,$D882,$C194,$B060,$BFEA),  { -0.657059080221453744877807438359e-6, }
             ($AAAA,$E200,$EF49,$EBA7,$BFE8),  { -0.219471747991019243265527185292e-6, }
             ($8D00,$F29B,$ADEE,$AC92,$3FE7),  { +0.803605497409122797720816119775e-7, }
             ($BDA2,$B8FE,$A9FF,$E227,$BFE3),  { -0.658197496457748650644446730955e-8, }
             ($E6FA,$9CCC,$7935,$BAD5,$BFE2),  { -0.271879079459221500311994000407e-8, }
             ($D6FF,$5A41,$B672,$844F,$3FE1),  { +0.962691975768783000937454605612e-9, }
             ($65A7,$6AD4,$C992,$9D5F,$BFDD),  { -0.715654866752415373221538623120e-10,}
             ($8BC1,$0A9B,$D68A,$A082,$BFDC),  { -0.364959955979361731397687347572e-10,}
             ($8D1A,$4AE9,$39D5,$D9C1,$3FDA),  { +0.123779266638284628266366968565e-10,}
             ($3F3C,$4BDF,$E79D,$E3E2,$BFD6),  { -0.809614940976175061407222222221e-12,}
             ($DEB2,$2C38,$B7FC,$9110,$BFD6),  { -0.515375501880531288979262837777e-12,}
             ($A6C7,$AEF0,$E86C,$BBB5,$3FD4),  { +0.166720484322176377906478380684e-12,}
             ($D5BD,$0F8E,$5A42,$A6D7,$BFD0),  { -0.926154824934769235381753853068e-14,}
             ($E767,$3BCB,$E74E,$87B5,$BFD0),  { -0.753344945832780041274713386623e-14,}
             ($1876,$90E5,$FF34,$A728,$3FCE),  { +0.231981301008194528825032891928e-14,}
             ($E9BB,$329F,$4C76,$F07A,$BFC9),  { -0.104290590578891200329403020712e-15,}
             ($8C77,$1995,$9394,$8224,$BFCA),  { -0.112880951868285154135078466421e-15,}
             ($A00C,$6101,$C59A,$986F,$3FC8),  { +0.330544205143281718828522986914e-16,}
             ($FD33,$8938,$9038,$A45A,$BFC3),  { -0.111370441904753970375841475032e-17,}
             ($5A2A,$F971,$8A77,$FE41,$BFC3),  { -0.172290580181768367727739031180e-17,}
             ($C00D,$5C97,$3D12,$8D83,$3FC2),  { +0.479463510149535766603385225170e-18,}
             ($426C,$D991,$2518,$FBE7,$BFBD)); { -0.266712580144022466545676300702e-19)}
var
  fdp05_0: array[0..np05_0 - 1] of extended absolute fdp05x0;
begin
  if x < ln_MinExt then sfc_fdp05 := 0.0
  else if x<=0.0 then begin
    t := -exp(x);
    sfc_fdp05 := -sfc_polylogr(1.5,t);
  end
  else if x<=4.0 then begin
    t := 0.5*x - 1.0;
    sfc_fdp05 := CSEvalX(t,fdp05_0,np05_0);
  end
  else if x<=12.0 then begin
    t := 0.25*x - 2.0;
    sfc_fdp05 := CSEvalX(t,fdp05_1,np05_1);
  end
  else if x<=36.0 then begin
    t := 0.5*x/SIXX - 2.0;
    sfc_fdp05 := CSEvalX(t,fdp05_2,np05_2)*(x*sqrt(x));
  end
  else sfc_fdp05 := fd_asymp_exp_real(0.5, x);
end;


{---------------------------------------------------------------------------}
function sfc_fdp15(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}
var
  t: extended;
{Note: THexExtW are needed with 16-bit compilers for relative errors}
{less than 5 eps_x, they are not needed for the 32+ bit compilers.}
const
  { t05 := t -> fd(3/2, 12*(t+2))/(12*(t+2))^(1/2);}
  { x=12*(t+2), x=12..36, t=x/12-2                 }
  { chebyshev(t05(t),t=-1..1,1e-20);               }
  np15_2 = 29;
  fdp15x2: array[0..np15_2 - 1] of THexExtW = (
             ($E207,$4F91,$ACD8,$C4D6,$4007),  {+3.9367714980969664415E+2}
             ($CEB6,$2A9E,$0CEA,$AD52,$4006),  {+1.7332050956287894902E+2}
             ($5BA4,$DC89,$79BE,$AD50,$4003),  {+2.1664294711213342222E+1}
             ($05BA,$478F,$2C2C,$DCF3,$3FF2),  {+2.1071423795317094221E-4}
             ($EE34,$17C2,$FCB3,$96AD,$BFF1),  {-7.1849649972728518409E-5}
             ($B36B,$4181,$FEBC,$C6FB,$3FEF),  {+2.3720783675115948540E-5}
             ($0519,$6F14,$C84A,$80C7,$BFEE),  {-7.6759099927198066013E-6}
             ($88B8,$DBD0,$582E,$A496,$3FEC),  {+2.4525416354303170927E-6}
             ($EE2C,$E0A9,$7415,$D092,$BFEA),  {-7.7699156039780901828E-7}
             ($A360,$CFB5,$A1B1,$8344,$3FE9),  {+2.4450587530326532190E-7}
             ($E27B,$7E6F,$607A,$A40F,$BFE7),  {-7.6396421480437313314E-8}
             ($336B,$BC60,$2D36,$CB27,$3FE5),  {+2.3650125793416319129E-8}
             ($770C,$CF48,$8756,$F861,$BFE3),  {-7.2288376791631086374E-9}
             ($2356,$3897,$8029,$9549,$3FE2),  {+2.1724133960664972412E-9}
             ($EBB9,$0777,$8582,$AFA5,$BFE0),  {-6.3899849344775435431E-10}
             ($D9A7,$825E,$AC5B,$C960,$3FDE),  {+1.8315188750975328184E-10}
             ($6462,$5940,$BE83,$E005,$BFDC),  {-5.0936805166409347329E-11}
             ($B67D,$1185,$6C51,$F0D0,$3FDA),  {+1.3688699754156613657E-11}
             ($138E,$20AF,$FAEF,$F90D,$BFD8),  {-3.5392788817089458740E-12}
             ($2D8A,$5F9D,$AE78,$F689,$3FD6),  {+8.7587827996357287805E-13}
             ($FE62,$A84D,$4B0B,$E7EC,$BFD4),  {-2.0598902148806642654E-13}
             ($4EC6,$5EB9,$BC82,$CCF7,$3FD2),  {+4.5511976448198640562E-14}
             ($46CB,$7098,$B707,$A690,$BFD0),  {-9.2462311586020405740E-15}
             ($745C,$0AA6,$080C,$ED2A,$3FCD),  {+1.6456571195203601819E-15}
             ($1D65,$D535,$8794,$FE27,$BFCA),  {-2.2044381299342686014E-16}
             ($A1ED,$43FC,$DCDE,$86CE,$3FC5),  {+3.6539797251750199626E-18}
             ($61A1,$472A,$2005,$E2FD,$3FC6),  {+1.2305085869833743187E-17}
             ($B276,$0238,$0DEB,$EA41,$BFC5),  {-6.3494706080876196549E-18}
             ($59BF,$ED88,$52F5,$AF53,$3FC4)); {+2.3761033915874246600E-18}
var
  fdp15_2: array[0..np15_2 - 1] of extended absolute fdp15x2;

const
  { t01 := t -> fd(3/2, 4*(t+2));    }
  { x = 4*(t+2), x=4..12, t = x/4-2  }
  { chebyshev(t01(t),t=-1..1,1e-20); }
  np15_1 = 28;
  fdp15x1: array[0..np15_1 - 1] of THexExtW = (
             ($D0FD,$CE4A,$1F4F,$909F,$4006),  {+1.4462157152925903243E+2}
             ($DF13,$60B5,$BF4C,$8E24,$4005),  {+7.1071771990590546640E+1}
             ($D653,$FBC0,$9338,$C98F,$4002),  {+1.2597552511779260681E+1}
             ($8BD7,$9F4D,$0B81,$8DBC,$3FFE),  {+5.5365058817488071785E-1}
             ($1E94,$120D,$BB17,$A5F6,$BFF9),  {-2.0259251985920947386E-2}
             ($9BD4,$E1AC,$C321,$F86C,$3FF5),  {+1.8953312109631679276E-3}
             ($4F03,$A374,$EEB1,$FA04,$BFF2),  {-2.3843695371615998960E-4}
             ($3A38,$8F26,$2645,$8288,$3FF0),  {+3.1121214385040858788E-5}
             ($FE69,$E2AA,$DCFE,$E1F1,$BFEC),  {-3.3668395634881968254E-6}
             ($5082,$A856,$E3D5,$8787,$3FE8),  {+1.2622291243644484818E-7}
             ($4B4F,$CB97,$D712,$9E3F,$3FE7),  {+7.3690607915601816666E-8}
             ($D950,$8E95,$E3D9,$8240,$BFE6),  {-3.0327000821131128333E-8}
             ($BC2E,$5F61,$5F6D,$894A,$3FE4),  {+7.9913599535400873690E-9}
             ($0E44,$6E7F,$8710,$EAD7,$BFE1),  {-1.7086997443833728992E-9}
             ($9F34,$2A53,$D90E,$A9B8,$3FDF),  {+3.0872263232876599441E-10}
             ($54C1,$C11D,$8E5C,$C732,$BFDC),  {-4.5292264252711262319E-11}
             ($0E8A,$49E5,$C28C,$94E4,$3FD9),  {+4.2318104524080512329E-12}
             ($2149,$EA9B,$6E4C,$AB2D,$3FD5),  {+3.0407225946293674652E-13}
             ($5EB2,$BCA0,$18C7,$A757,$BFD5),  {-2.9725594761874387986E-13}
             ($6DAD,$C831,$18BC,$E36E,$3FD3),  {+1.0099923784566557989E-13}
             ($2829,$BE32,$4353,$E97C,$BFD1),  {-2.5922086958100736660E-14}
             ($1901,$6F98,$AF3B,$C82B,$3FCF),  {+5.5558514061788275762E-15}
             ($2F0C,$4096,$2FE0,$8F76,$BFCD),  {-9.9546529373470499049E-16}
             ($F65C,$93F7,$0DE1,$9C94,$3FCA),  {+1.3581005834664229560E-16}
             ($C873,$795A,$5BEC,$86DE,$BFC6),  {-7.3112409217500735597E-18}
             ($1264,$00F5,$62E9,$8130,$BFC5),  {-3.5016751122371363843E-18}
             ($E86D,$655C,$E73D,$8124,$3FC4),  {+1.7502296465398255808E-18}
             ($2BBE,$9519,$6F35,$9DF0,$BFC2)); {-5.3511881381428642648E-19}
var
  fdp15_1: array[0..np15_1 - 1] of extended absolute fdp15x1;

const
  { t00 := t -> fd(3/2, 2*(t+1));    }
  { x = 2*(t+1), x=0..4, t = x/2 - 1 }
  { chebyshev(t00(t),t=-1..1,1e-20); }
  np15_0 = 30;
  fdp15x0: array[0..np15_0 - 1] of THexExtW = (
             ($AA8A,$8D18,$5613,$B3C7,$4002),  {+1.1236166073199948765E+1}
             ($2A3E,$C299,$40E5,$C1EF,$4001),  {+6.0604557502272531624}
             ($369D,$0B58,$CE6F,$B97A,$3FFF),  {+1.4490602533539360290}
             ($4C36,$A850,$9807,$8C0B,$3FFC),  {+1.3676297709214706240E-1}
             ($4A77,$329F,$74F2,$E22F,$BFF6),  {-3.4513149750616713262E-3}
             ($55F2,$5A13,$F3CB,$93AD,$BFF4),  {-5.6335258963187187413E-4}
             ($FDF9,$C386,$0F1F,$A21E,$3FF2),  {+1.5460721805105655254E-4}
             ($350D,$0351,$4ED1,$CA93,$BFEE),  {-1.2074436032395957859E-5}
             ($FEDA,$0AC1,$0056,$A691,$BFEC),  {-2.4820329454318314047E-6}
             ($1DB7,$9078,$0EAE,$E7DA,$3FEA),  {+8.6371521103341848976E-7}
             ($AF35,$9244,$12A8,$9E5C,$BFE7),  {-7.3741962996236602468E-8}
             ($54E6,$D486,$0D12,$A63F,$BFE5),  {-1.9353615729676523342E-8}
             ($9AF5,$2A02,$CCF6,$EDE1,$3FE3),  {+6.9232783779587078981E-9}
             ($8215,$CE7F,$2238,$9F87,$BFE0),  {-5.8035899541125150058E-10}
             ($BB33,$C1AB,$67E5,$CFE7,$BFDE),  {-1.8908752199406953327E-10}
             ($F285,$93D3,$857E,$927B,$3FDD),  {+6.6612531424447944938E-11}
             ($974B,$1431,$F0CD,$B897,$BFD9),  {-5.2464633336918750095E-12}
             ($4A22,$E6FC,$CE0A,$93B7,$BFD8),  {-2.0991988621741175341E-12}
             ($3476,$7DEF,$7684,$C99E,$3FD6),  {+7.1629456476161078420E-13}
             ($0782,$7232,$1C65,$E76C,$BFD2),  {-5.1386075015702707325E-14}
             ($9EF4,$D32A,$F751,$E3EE,$BFD1),  {-2.5305697681559179825E-14}
             ($EBA8,$DD77,$E0AB,$957A,$3FD0),  {+8.2978063705001989763E-15}
             ($9F96,$3CEE,$58C5,$97BB,$BFCC),  {-5.2642551179225625266E-16}
             ($D1B4,$208E,$67F1,$BA33,$BFCB),  {-3.2300690729343082800E-16}
             ($459C,$C353,$B5E7,$E9B9,$3FC9),  {+1.0136224841459294120E-16}
             ($5530,$915C,$7A68,$CAAF,$BFC5),  {-5.4938004437287676996E-18}
             ($5440,$9452,$743F,$9E98,$BFC5),  {-4.2987402865091326124E-18}
             ($214C,$693F,$1EBB,$BE15,$3FC3),  {+1.2880491228202175925E-18}
             ($B48E,$1FC3,$3F4A,$8659,$BFBF),  {-5.6898854614185145094E-20}
             ($1422,$0469,$F0EC,$8B6D,$BFBF)); {-5.9050672129271432028E-20}
var
  fdp15_0: array[0..np15_0 - 1] of extended absolute fdp15x0;
begin
  if x < ln_MinExt then sfc_fdp15 := 0.0
  else if x<=0.0 then begin
    t := -exp(x);
    sfc_fdp15 := -sfc_polylogr(2.5,t);
  end
  else if x<=4.0 then begin
    t := 0.5*x - 1.0;
    sfc_fdp15 := CSEvalX(t,fdp15_0,np15_0);
  end
  else if x<=12.0 then begin
    t := 0.25*x - 2.0;
    sfc_fdp15 := CSEvalX(t,fdp15_1,np15_1);
  end
  else if x<=36.0 then begin
    t := 0.5*x/SIXX - 2.0;
    sfc_fdp15 := CSEvalX(t,fdp15_2,np15_2)*sqrt(x);
  end
  else sfc_fdp15 := fd_asymp_exp_real(1.5, x);
end;


{---------------------------------------------------------------------------}
function sfc_fdp25(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}
var
  t: extended;
const
  { t05 := t -> fd(5/2, 12*(t+2))/(12*(t+2))^(1/2);}
  { x=12*(t+2), x=12..36, t=x/12-2                 }
  { chebyshev(t05(t),t=-1..1,0.5e-20);             }
  np25_2 = 29;
  fdp25x2: array[0..np25_2 - 1] of THexExtW = (
             ($AD17,$6248,$CEEB,$CFFC,$400A),  {+3327.80051744835038169140925270     }
             ($1823,$A517,$6B96,$EE9E,$4009),  {+1908.95063335651571308582790355     }
             ($0A3E,$CFEA,$AE96,$DED7,$4007),  {+445.685015536804785376765540320     }
             ($2E7C,$91D2,$13E4,$948D,$4004),  {+37.1377711977006996644597888501     }
             ($0A1A,$5609,$8FDD,$8FC0,$3FF4),  {+0.548371105955101354822268394464e-3 }
             ($CBE5,$334C,$8FDA,$9C58,$BFF2),  {-0.149103112232773501712826143475e-3 }
             ($937B,$A86D,$1A65,$AAAA,$3FF0),  {+0.406895793170452009031534941051e-4 }
             ($576D,$05C7,$CD61,$BB14,$BFEE),  {-0.111509119769784098254022652378e-4 }
             ($C9FA,$F3B1,$8991,$CE0A,$3FEC),  {+0.307025256236664354012898996850e-5 }
             ($936A,$522B,$E4EB,$E406,$BFEA),  {-0.849466512091940332152581052807e-6 }
             ($BDB7,$ED86,$D3CA,$FD7B,$3FE8),  {+0.236075092526091114586850623705e-6 }
             ($84BF,$FE3B,$EFF3,$8D58,$BFE7),  {-0.658200175472105542226242341390e-7 }
             ($1C9B,$8B7D,$20DA,$9DCF,$3FE5),  {+0.183713965883408450720821603225e-7 }
             ($45CC,$9550,$4D19,$AFDB,$BFE3),  {-0.511810198663846546954610261831e-8 }
             ($4A34,$1BDE,$79FD,$C2E7,$3FE1),  {+0.141811473549707654091087714101e-8 }
             ($A935,$ACD5,$4833,$D608,$BFDF),  {-0.389322579781747804677424437470e-9 }
             ($94EE,$94A7,$6F50,$E80A,$3FDD),  {+0.105519921375420776489151857017e-9 }
             ($5258,$55A5,$CEBA,$F793,$BFDB),  {-0.281462886467752204120551871327e-10}
             ($D355,$A6B0,$EF77,$81A9,$3FDA),  {+0.737053427612575970277777777777e-11}
             ($86F3,$28EC,$0211,$851D,$BFD8),  {-0.189165394872340800854140566911e-11}
             ($DEB1,$2A6B,$B718,$85D2,$3FD6),  {+0.475435180287273717546904411718e-12}
             ($2F66,$C156,$8547,$83C1,$BFD4),  {-0.117022782527244997933790212172e-12}
             ($DB32,$AD39,$6BF7,$FE52,$3FD1),  {+0.282354095607738906708933732942e-13}
             ($BDB0,$1DD6,$F9AB,$F110,$BFCF),  {-0.669093418553307654841147190614e-14}
             ($B86D,$48C1,$0F77,$E111,$3FCD),  {+0.156171355182400009442588371477e-14}
             ($A7CF,$EC31,$02FF,$CFC7,$BFCB),  {-0.360436315380567126678385530540e-15}
             ($85D8,$54AA,$3369,$BE97,$3FC9),  {+0.826555092713945341135435212925e-16}
             ($5E4B,$7363,$759C,$AEAE,$BFC7),  {-0.189390042389614168443291988842e-16}
             ($48AC,$9D1B,$5E07,$A0EA,$3FC5)); {+0.436162329343858784694936828122e-17}
var
  fdp25_2: array[0..np25_2 - 1] of extended absolute fdp25x2;

const
  { t01 := t -> fd(5/2, 4*(t+2));     }
  { x = 4*(t+2), x=4..12, t = x/4-2   }
  { chebyshev(t01(t),t=-1..1,0.5e-20);}
  np25_1 = 26;
  fdp25x1: array[0..np25_1 - 1] of THexExtW = (
             ($BA1B,$577F,$3C0D,$DF20,$4007),  {+446.251832645153130169370289260     }
             ($C398,$3E8E,$261C,$8406,$4007),  {+264.048038034959543495722343526     }
             ($43FB,$5D77,$4735,$8D09,$4005),  {+70.5181214024156659192952476879     }
             ($9342,$AF86,$09B9,$8697,$4002),  {+8.41187450917678775219717559168     }
             ($B589,$0E5C,$D520,$8D3F,$3FFD),  {+0.275877628481958774952182034146    }
             ($C4C5,$8BD1,$5A94,$8335,$BFF8),  {-0.800832601288191495847914821399e-2 }
             ($9C02,$2377,$ABF5,$A2E5,$3FF4),  {+0.621403332192709022962371969718e-3 }
             ($1405,$CDAA,$CD47,$8CD9,$BFF1),  {-0.671628897579062265051072717362e-4 }
             ($E3E7,$B97D,$9E61,$8200,$3FEE),  {+0.774874786815110348480827930088e-5 }
             ($A43C,$A478,$5186,$CD3C,$BFEA),  {-0.764562260311955253809364235157e-6 }
             ($6BDE,$3CC9,$B0A3,$8679,$3FE6),  {+0.313099826515151953013264804033e-7 }
             ($79DF,$FF84,$1035,$CD38,$3FE4),  {+0.119453178112839507810392037265e-7 }
             ($2BD3,$642D,$DF81,$A3E2,$BFE3),  {-0.476971684612462590592051464946e-8 }
             ($9C8F,$6012,$F6F2,$A271,$3FE1),  {+0.118194420326328021143663686090e-8 }
             ($8211,$E043,$6610,$82A3,$BFDF),  {-0.237629640018665948118995031809e-9 }
             ($4FF9,$C584,$288C,$B28E,$3FDC),  {+0.405987762501810590898667771913e-10}
             ($8B03,$5AF2,$E939,$C888,$BFD9),  {-0.569954206402177488280324465353e-11}
             ($1348,$9BF0,$8B62,$95FA,$3FD6),  {+0.532831341179622954769972002079e-12}
             ($42F4,$481A,$4751,$CB3C,$3FD1),  {+0.225636690685856850370370370370e-13}
             ($6914,$3684,$1ECC,$80A1,$BFD2),  {-0.285614590169098049435176590815e-13}
             ($B9E6,$742C,$7E3A,$ABEF,$3FD0),  {+0.954433864394867464106024538843e-14}
             ($BECA,$CAC9,$E7DF,$AB0F,$BFCE),  {-0.237396396803486035394564496255e-14}
             ($00ED,$0229,$2209,$8E05,$3FCC),  {+0.492731031621107824201139948487e-15}
             ($E723,$BCB6,$0464,$C622,$BFC9),  {-0.859264393750395508067640840442e-16}
             ($5DA4,$7AA9,$6BF6,$D627,$3FC6),  {+0.116093111215732467860177907410e-16}
             ($DFD7,$0357,$3A8A,$D5F5,$BFC2));  {-0.724917645463107688269430437231e-18}
var
  fdp25_1: array[0..np25_1 - 1] of extended absolute fdp25x1;

const
  { t00 := t -> fd(5/2, 2*(t+1));     }
  { x = 2*(t+1), x=0..4, t = x/2 - 1  }
  { chebyshev(t00(t),t=-1..1,0.5e-20);}
  np25_0 = 29;
  fdp25x0: array[0..np25_0 - 1] of THexExtW = (
             ($B76E,$8778,$94F2,$833B,$4003),  {16.4040926883309016968622975910      }
             ($A3B7,$ABAD,$FC45,$9C97,$4002),  {+9.78710581984601273608560471355     }
             ($A7DD,$8556,$E425,$BD8E,$4000),  {+2.96184638656755305010270048774     }
             ($B303,$05EC,$32E2,$F7E5,$3FFD),  {+0.484170522776332566753457344729    }
             ($5F8C,$73AA,$45FB,$8C9F,$3FFA),  {+0.343315824204447335660552492644e-1 }
             ($61DF,$2579,$DE50,$BD0D,$BFF4),  {-0.721184438622545575771793462833e-3 }
             ($6C28,$1D5D,$88C0,$C0AF,$BFF1),  {-0.918796922665793193774947944856e-4 }
             ($2B60,$4439,$5F01,$BC40,$3FEF),  {+0.224413215709269119915539751533e-4 }
             ($B6E8,$EC58,$EFBB,$D910,$BFEB),  {-0.161726890542867204359418139982e-5 }
             ($4F72,$CC9C,$FFC7,$8FA8,$BFE9),  {-0.267587886937288311365673143759e-6 }
             ($0D18,$527D,$38DF,$BDA3,$3FE7),  {+0.883068826763095013112111498740e-7 }
             ($52C0,$071E,$8ADD,$FBF7,$BFE3),  {-0.733320376129048276026566510888e-8 }
             ($CBC7,$BD6D,$C556,$D703,$BFE1),  {-0.156443806118877265343058350412e-8 }
             ($ACCE,$5D93,$18BE,$9663,$3FE0),  {+0.547105069227136725477278628827e-9 }
             ($2506,$B7F9,$5EC0,$CB3E,$BFDC),  {-0.462122519168356746790325305029e-10}
             ($0F6F,$A33E,$6F54,$D79C,$BFDA),  {-0.122560705773585105508047077120e-10}
             ($D4D6,$EB0B,$43EE,$9719,$3FD9),  {+0.429448314291387890453274064781e-11}
             ($1C18,$E596,$781B,$C574,$BFD5),  {-0.350750464614910929036192564942e-12}
             ($01CB,$2BBC,$36C1,$8017,$BFD4),  {-0.113767377064356379270370370370e-12}
             ($9C28,$C72F,$5CD6,$AFC8,$3FD2),  {+0.390315927601668402160150340116e-13}
             ($01EE,$2140,$AA0D,$D708,$BFCE),  {-0.298419406931014531473839059196e-14}
             ($D6BF,$FADD,$098E,$AA0D,$BFCD),  {-0.117996534141747254896765099982e-14}
             ($8EB9,$4F40,$9CF9,$E1E3,$3FCB),  {+0.391855148990619535770411408698e-15}
             ($A5CD,$0770,$C7AA,$FBC0,$BFC7),  {-0.272951200089934340373538288281e-16}
             ($60E2,$530C,$8D5F,$F40B,$BFC6),  {-0.132297127854042464240388752538e-16}
             ($1A6D,$124D,$97CA,$9BED,$3FC5),  {+0.422643954804408384069683831290e-17}
             ($D24C,$B86D,$14FA,$99F9,$BFC1),  {-0.260840367944197757647285622463e-18}
             ($CA61,$59E0,$250C,$B97A,$BFC0),  {-0.157105238218333575135078217011e-18}
             ($E4C1,$E607,$469D,$E332,$3FBE)); {+0.481107069624679302676225882053e-19}
var
  fdp25_0: array[0..np25_0 - 1] of extended absolute fdp25x0;
begin
  if x < ln_MinExt then sfc_fdp25 := 0.0
  else if x<=0.0 then begin
    t := -exp(x);
    sfc_fdp25 := -sfc_polylogr(3.5,t);
  end
  else if x<=4.0 then begin
    t := 0.5*x - 1.0;
    sfc_fdp25 := CSEvalX(t,fdp25_0,np25_0);
  end
  else if x<=12.0 then begin
    t := 0.25*x - 2.0;
    sfc_fdp25 := CSEvalX(t,fdp25_1,np25_1);
  end
  else if x<=36.0 then begin
    t := 0.5*x/SIXX - 2.0;
    sfc_fdp25 := CSEvalX(t,fdp25_2,np25_2)*sqrt(x);
  end
  else sfc_fdp25 := fd_asymp_exp_real(2.5, x);
end;


{---------------------------------------------------------------------------}
function sfc_fdm05(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}
var
  t: extended;
{Note: THexExtW are needed with 16-bit compilers for relative errors}
{less than 5 eps_x, they are not needed for the 32+ bit compilers.}
const
  { t05 := t -> fd(-1/2, 14*(t+2)-2)/(14*(t+2)-2)^(1/2); }
  { x = 14*(t+2)-2, x=12..40, t = (x+2)/14-2             }
  { chebyshev(t05(t),t=-1..1,1e-20);                     }
  nm05_2 = 36;
  fdm05x2: array[0..nm05_2 - 1] of THexExtW = (
             ($39CB,$588B,$5F29,$9048,$4000),  {2*1.12720860976068042985126031469,      }
             ($EF98,$4015,$5936,$A6EC,$3FF5),  { +0.127352322590483468748105626469e-2,  }
             ($C2CB,$28D0,$8398,$9058,$BFF4),  { -0.550635368235688347038739063159e-3,  }
             ($CCBD,$142D,$349D,$E222,$3FF2),  { +0.215657821006449458356062287107e-3,  }
             ($187E,$C82E,$1C87,$A7E5,$BFF1),  { -0.800585587450641857019601765123e-4,  }
             ($18D4,$F791,$E75B,$F170,$3FEF),  { +0.287820136859342153053116112911e-4,  }
             ($4E41,$9762,$7838,$AA28,$BFEE),  { -0.101422121779796149201023928609e-4,  }
             ($1709,$3850,$27A3,$ECD1,$3FEC),  { +0.352884845542076071119669637285e-5,  }
             ($0C55,$0821,$754B,$A371,$BFEB),  { -0.121774670482100079163136672498e-5,  }
             ($C921,$FF00,$DD5E,$E042,$3FE9),  { +0.417719018374740984434001320844e-6,  }
             ($A682,$1C8B,$B733,$98FD,$BFE8),  { -0.142484043403418353187352871031e-6,  }
             ($3B14,$CB5D,$1979,$CF4B,$3FE6),  { +0.482642458470626719951342032248e-7,  }
             ($DCA7,$DE58,$F41C,$8B15,$BFE5),  { -0.161917130596971286058653227519e-7,  }
             ($FC0E,$2385,$3937,$B82D,$3FE3),  { +0.536024612047194258605653270535e-8,  }
             ($097A,$2926,$6F3E,$EFA9,$BFE1),  { -0.174376948911707269461907228582e-8,  }
             ($2D95,$9F93,$7ED2,$9890,$3FE0),  { +0.555026181833641050021767137440e-9,  }
             ($A662,$6548,$2EF3,$BD38,$BFDE),  { -0.172094102173839496709895590169e-9,  }
             ($B081,$C2C2,$B0CD,$E39E,$3FDC),  { +0.517547699271684434423133548903e-10, }
             ($9E8D,$458C,$0C22,$842D,$BFDB),  { -0.150266676426121966101287037037e-10, }
             ($6206,$B5EE,$3062,$9367,$3FD9),  { +0.418944756732711716944423047911e-11, }
             ($A9CA,$26C7,$AEBD,$9CB5,$BFD7),  { -0.111348937219875590656629403579e-11, }
             ($65E6,$3122,$2843,$9D0C,$3FD5),  { +0.278972381829081126481744940963e-12, }
             ($8286,$BBFF,$3CB9,$915D,$BFD3),  { -0.645546762036873772076350520773e-13, }
             ($954D,$E342,$5E99,$EDBA,$3FD0),  { +0.131965552928934103834300756587e-13, }
             ($D4D8,$E7DE,$0CA9,$9669,$BFCE),  { -0.208736291421929647980739225035e-14, }
             ($4E6F,$4C39,$B018,$C962,$3FC9),  { +0.873370384238809757090654794064e-16, }
             ($A9DA,$3F4D,$0930,$8C9A,$3FCA),  { +0.121952537222219205275296984895e-15, }
             ($301F,$34B4,$0170,$AE73,$BFC9),  { -0.756552983007583718512206829240e-16, }
             ($831A,$E68D,$DA2D,$9665,$3FC8),  { +0.326123373969832726648499559418e-16, }
             ($37C7,$E169,$55EF,$E141,$BFC6),  { -0.122111098113026111719660711879e-16, }
             ($9929,$EEEF,$5F4B,$9C4F,$3FC5),  { +0.423679233680363253635062387618e-17, }
             ($6278,$BFF2,$85BF,$CECF,$BFC3),  { -0.140140337077744622807440276665e-17, }
             ($8304,$2C42,$89BC,$846B,$3FC2),  { +0.448656650094168829272055386020e-18, }
             ($B535,$0E51,$F2CA,$A599,$BFC0),  { -0.140269808809750966318869781429e-18, }
             ($0C9F,$83C5,$0325,$CB59,$3FBE),  { +0.430605513146981698934338690424e-19, }
             ($CFAB,$2EFA,$877A,$F5FB,$BFBC)); { -0.130222070259439058246940543349e-19, }
var
  fdm05_2: array[0..nm05_2 - 1] of extended absolute fdm05x2;

const
  { t01 := t -> fd(-1/2, 4*(t+2));   }
  { x = 4*(t+2), x=4..12, t = x/4-2  }
  { chebyshev(t01(t),t=-1..1,1e-20); }
  nm05_1 = 30;
  fdm05x1: array[0..nm05_1 - 1] of THexExtW = (
             ($81A4,$8736,$836A,$C6DC,$4001),  {2*3.10720906642095693549141663063,     }
             ($45F7,$5FCA,$AF0A,$D868,$3FFE),  { +0.845347347290762745698264830651,    }
             ($1E8C,$4513,$5574,$8287,$BFFB),  { -0.637346912486776414955579558212e-1, }
             ($46F4,$B77F,$F5FA,$A43A,$3FF8),  { +0.100238229887210249047746267365e-1, }
             ($B06D,$233D,$BE6B,$F45E,$BFF5),  { -0.186439585156115756577336126739e-2, }
             ($48AE,$7077,$BEF7,$AC3F,$3FF3),  { +0.328538909279736813909986738098e-3, }
             ($4941,$98D5,$517C,$B670,$BFF0),  { -0.434967859717546117551553285448e-4, }
             ($A907,$0E55,$C147,$A038,$3FEA),  { +0.596872341564135171934135718907e-6, }
             ($CCE1,$1506,$2D65,$91E3,$3FEC),  { +0.217389183372963314012513808309e-5, }
             ($C12B,$E0C1,$195D,$8124,$BFEB),  { -0.962175518760746788464909336684e-6, }
             ($554C,$ADFC,$E35F,$9BCF,$3FE9),  { +0.290222583707548004030466363153e-6, }
             ($9D7E,$936B,$116B,$9868,$BFE7),  { -0.709698143468401550945013829553e-7, }
             ($C85D,$2593,$0FA1,$F848,$3FE4),  { +0.144518846929383515210221101027e-7, }
             ($85AC,$01A5,$4030,$9F04,$BFE2),  { -0.231399614747943009400043966811e-8, }
             ($3DA3,$9ABD,$2266,$E068,$3FDE),  { +0.204096772826715027518323729570e-9, }
             ($B411,$EAC9,$4F3F,$A9DA,$3FDC),  { +0.386200489979488286568147019637e-10,}
             ($EC72,$2FC2,$9DE0,$F3D4,$BFDB),  { -0.277203223495739126097345324431e-10,}
             ($7B46,$206C,$1E0C,$AC33,$3FDA),  { +0.978841838084900111240803825558e-11,}
             ($B41F,$F864,$2F5A,$BCD1,$BFD8),  { -0.268325278562056093253240740740e-11,}
             ($739A,$535B,$6EA0,$ACED,$3FD6),  { +0.614361785526356233322403621196e-12,}
             ($D6FF,$2E94,$1FF1,$82A3,$BFD4),  { -0.116029147310593134362697691583e-12,}
             ($B4A5,$1B19,$2407,$8F27,$3FD1),  { +0.158931638394253582586663743109e-13,}
             ($AC3A,$995A,$507D,$E126,$BFCB),  { -0.390572410679751374699284640444e-15,}
             ($9A87,$790E,$050A,$CFF4,$BFCC),  { -0.721482602544015576484549245755e-15,}
             ($6DE8,$98D3,$BA37,$C105,$3FCB),  { +0.334840441309012455905148106174e-15,}
             ($11C8,$0955,$64DF,$F191,$BFC9),  { -0.104763396492153632437572704125e-15,}
             ($2E77,$2FCC,$5928,$F5C5,$3FC7),  { +0.266465334690489637299281537629e-16,}
             ($64F7,$6AC2,$3045,$D14D,$BFC5),  { -0.567312900783804305906703267603e-17,}
             ($658A,$7DA0,$0C1E,$8E07,$3FC3),  { +0.962415969624123143993085069845e-18,}
             ($728A,$D261,$C196,$E75E,$BFBF)); { -0.979890664035765026812190968018e-19)}
var
  fdm05_1: array[0..nm05_1 - 1] of extended absolute fdm05x1;

const
  { t00 := t -> fd(-1/2, 2*(t+1));   }
  { x = 2*(t+1), x=0..4, t = x/2 - 1 }
  { chebyshev(t00(t),t=-1..1,1e-20); }
  nm05_0 = 33;
  fdm05x0: array[0..nm05_0 - 1] of THexExtW = (
             ($C466,$7A74,$6418,$B675,$4000),  {2*1.42545748896819736142428828509,     }
             ($F317,$B08B,$13D3,$CD85,$3FFE),  { +0.802811850721397311632302091077,    }
             ($B45B,$9E73,$EA41,$8C85,$BFFA),  { -0.343073988568933129244642294564e-1, }
             ($FA43,$49CC,$B004,$C441,$BFF8),  { -0.119785517610154200513442344857e-1, }
             ($BF4F,$8636,$BF31,$8FC0,$3FF7),  { +0.438699088685875336070850350529e-2, }
             ($3DF7,$E4EA,$84DC,$D39E,$BFF3),  { -0.403631620076134737472343018429e-3, }
             ($0DF4,$7406,$5B8A,$9BE1,$BFF2),  { -0.148659041454562690460518932844e-3, }
             ($683B,$ECE9,$87E3,$F3FC,$3FF0),  { +0.581709023766370096169389361999e-4, }
             ($C27D,$1617,$BE83,$A9FA,$BFED),  { -0.506578294882995558003563539688e-5, }
             ($3DDC,$68DB,$9262,$97B8,$BFEC),  { -0.226081883384096770040396740859e-5, }
             ($A089,$EBD7,$DA22,$E390,$3FEA),  { +0.847748773163128123864631556434e-6, }
             ($663E,$735A,$9A87,$8DF3,$BFE7),  { -0.661013539307752677486955427753e-7, }
             ($156A,$1BE4,$505A,$9B8D,$BFE6),  { -0.362172739869069536282661782116e-7, }
             ($3FE0,$BB25,$27C2,$DD51,$3FE4),  { +0.128823456441545703286380320048e-7, }
             ($6FB1,$1F66,$77A3,$EFF7,$BFE0),  { -0.872993657208158587706978751116e-9, }
             ($85E0,$2CDF,$7856,$A3A5,$BFE0),  { -0.595342016608391684486334590438e-9, }
             ($3BA3,$A39F,$893D,$DC70,$3FDE),  { +0.200488642920464472125344875726e-9, }
             ($87C0,$90FE,$7308,$C8A8,$BFDA),  { -0.114060870414129142501077323872e-10,}
             ($A47D,$BF8E,$434F,$AECC,$BFDA),  { -0.993611036461939592748999999999e-11,}
             ($7D8D,$9172,$4621,$DEDB,$3FD8),  { +0.316698189615823685519545765746e-11,}
             ($0253,$C2C7,$3512,$A21A,$BFD4),  { -0.143975828889301436878708068523e-12,}
             ($2EC9,$93A5,$C664,$BC81,$BFD4),  { -0.167427790285290702945465438147e-12,}
             ($D8D3,$C691,$C58F,$E381,$3FD2),  { +0.505166843470001025439981895705e-13,}
             ($4D48,$B668,$2C80,$F3D4,$BFCD),  { -0.169190220207909386709653151896e-14,}
             ($A558,$9B77,$94BC,$CC92,$BFCE),  { -0.283901488488463908424013117515e-14,}
             ($22DD,$2B0F,$8C61,$E9C6,$3FCC),  { +0.811071971814294943755263204708e-15,}
             ($CDC5,$BEFC,$1299,$9CB7,$BFC7),  { -0.169910881775102181849224823698e-16,}
             ($14B0,$2938,$4B9B,$DEF1,$BFC8),  { -0.483429615582375261924591357742e-16,}
             ($AF92,$5125,$818B,$F143,$3FC6),  { +0.130789311367733023108178561868e-16,}
             ($335A,$EAC4,$E45D,$EFE4,$BFBF),  { -0.101599107342478990354348239320e-18,}
             ($B0EF,$6AD3,$D25F,$F3A5,$BFC2),  { -0.825510657563228475743326075971e-18,}
             ($C04D,$6680,$FD3F,$F9AE,$3FC0),  { +0.211490194721918999612264245007e-18,}
             ($9607,$35A4,$9CEC,$8562,$BFBD)); { -0.141227081204025177342262136401e-19)}
var
  fdm05_0: array[0..nm05_0 - 1] of extended absolute fdm05x0;
begin
  if x < ln_MinExt then sfc_fdm05 := 0.0
  else if x<=0.0 then begin
    t := -exp(x);
    sfc_fdm05 := -sfc_polylogr(0.5,t);
  end
  else if x<=4.0 then begin
    t := 0.5*x - 1.0;
    sfc_fdm05 := CSEvalX(t,fdm05_0,nm05_0);
  end
  else if x<=12.0 then begin
    t := 0.25*x - 2.0;
    sfc_fdm05 := CSEvalX(t,fdm05_1,nm05_1);
  end
  else if x<=40.0 then begin
    t := (x+2.0)/14.0 - 2.0;
    sfc_fdm05 := CSEvalX(t,fdm05_2,nm05_2)*sqrt(x);
  end
  else sfc_fdm05 := fd_asymp_exp_real(-0.5, x);
end;

{---------------------------------------------------------------------------}
function sfc_harmonic(x: extended): extended;
  {-Return the harmonic number function H(x) = psi(x+1) + EulerGamma}
const
  nch = 22;
  chah: array[0..nch-1] of THexExtW = (    {chebyshev((Psi(x+1)+gamma)/x, x=-1/4..1/4, 0.2e-20);}
          ($F6A2,$926A,$FC75,$D714,$4000), {+3.36065589410433977112951870486     }
          ($5FDE,$33E7,$08D8,$A06B,$BFFD), {-0.313316608802735484720854737533    }
          ($8260,$0440,$04BD,$932D,$3FFA), {+0.359316048709531508723414593246e-1 }
          ($44F0,$F981,$F3AC,$8F82,$BFF7), {-0.437962434984843410853415224307e-2 }
          ($A67E,$5757,$14E3,$8F35,$3FF4), {+0.546292686372267141546977073957e-3 }
          ($062A,$3C93,$0D31,$904F,$BFF1), {-0.688117957348770917553696097164e-4 }
          ($92EC,$DBDA,$5E65,$920F,$3FEE), {+0.870585645123897195412611399047e-5 }
          ($27F1,$EF11,$C316,$9422,$BFEB), {-0.110369763769828900468418914839e-5 }
          ($291C,$25D3,$3821,$9663,$3FE8), {+0.140059343742157854915062263457e-6 }
          ($B194,$8B66,$CF73,$98BE,$BFE5), {-0.177818994219795504449052141370e-7 }
          ($F8AA,$12D6,$F085,$9B2C,$3FE2), {+0.225810137695486223694241640848e-8 }
          ($F0BF,$8F90,$797C,$9DA9,$BFDF), {-0.286785525509693105456693572379e-9 }
          ($2E1E,$F307,$74CF,$A032,$3FDC), {+0.364246022659874814089675465140e-10}
          ($C0E8,$EC99,$FEFC,$A2C6,$BFD9), {-0.462640992702404867971084536579e-11}
          ($09AE,$B22C,$BC9F,$A566,$3FD6), {+0.587623516675841187923320281923e-12}
          ($3A67,$E81B,$966E,$A811,$BFD3), {-0.746374969291138217069111742907e-13}
          ($01E8,$7468,$977C,$AAC7,$3FD0), {+0.948017526953988879040131837225e-14}
          ($9F33,$E18B,$DC66,$AD88,$BFCD), {-0.120413826871905020093317205030e-14}
          ($D7AF,$C348,$8AF5,$B055,$3FCA), {+0.152945496208458532222222222222e-15}
          ($074A,$BB31,$CDE8,$B32D,$BFC7), {-0.194266177768451363581139754366e-16}
          ($A99C,$96E0,$D2DF,$B611,$3FC4), {+0.246750352250719046709039799083e-17}
          ($C284,$4957,$C94F,$B901,$BFC1));{-0.313414011663104260448782764313e-18}
       (* ($9DF1,$9E75,$E1FD,$BBFD,$3FBE), {+0.398087971944463065142369984751e-19}
          ($58CE,$7DDF,$4EC0,$BF06,$BFBB));{-0.505638009703096728955241279549e-20} *)
var
  cha: array[0..nch-1] of extended absolute chah;
var
  y: extended;
begin
  if IsNanOrInf(x) then begin
    if x=PosInf_x then sfc_harmonic := PosInf_x
    else sfc_harmonic := NaN_x;
  end
  else if frac(x)=0.0 then begin
    if x<0.0 then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfc_harmonic := Nan_x;
    end
    else if x<=12.0 then begin
      {sfc_psi uses a similar loop, avoid (y - EulerGamma) + EulerGamma}
      y := 0.0;
      while x>=1.0 do begin
        y := y + 1.0/x;
        x := x - 1.0;
      end;
      sfc_harmonic := y;
    end
    else sfc_harmonic := sfc_psi(x+1.0) + EulerGamma;
  end
  else begin
    if abs(x)<=0.25 then begin
      y := CSEvalX(4.0*x, cha, nch);
      sfc_harmonic := y*x;
    end
    else begin
      y := sfc_psi(x+1.0);
      sfc_harmonic :=  y + EulerGamma;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function harm2core(x,r: extended): extended;
  {-Core routine for H(x,r), x >= -1, no checks, r<>0,1,-1}
  { The alternative computations use Hurwitz zeta or series for |x|<0.5}
var
  a,s,t: extended;
  k: integer;
begin
  {Get threshold value for series}
  if r > 1.0 then t := 0.125 else t := 0.5;
  if abs(x)<=t then begin
    {http://functions.wolfram.com/06.17.06.0018.01 for z0=0, note that in}
    {.../06.17.06.0002.01 there are the restriction Re(r) > 1 and |x| < 1}
    a := r*x;
    s := a*sfc_zeta1p(r);
    k := 1;
    repeat
      a := (r+k)*x*a;
      k := k+1;
      a := -a/k;
      t := s;
      s := t + a*sfc_zeta(r+k);
    until t=s;
    harm2core := s;
  end
  else begin
    s := sfc_zeta(r);
    a := sfc_zetah(r, x+1.0);
    harm2core := s - a;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_harmonic2(x,r: extended): extended;
  {-Return the generalized harmonic function H(x,r) = zeta(r)-zetah(r,x+1); x >= -1}
const
  XREC = 20.0;  {Max x for recursive call}
var
  a,s,t: extended;
  n: integer;
begin
  if IsNanOrInf(x) or IsNanOrInf(x) or (x < -1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_harmonic2 := NaN_x;
    exit;
  end;

  if (frac(x)=0.0) and (x <= 6.0)  then begin
    {Direct summation with exp<x> functions for small integer x}
    n := round(x);
    if n<2 then sfc_harmonic2 := n
    else if n=2 then sfc_harmonic2 := 1.0 + exp2(-r)
    else begin
      a := exp2(-r);
      s := exp3(-r);
      {Avoid brain-damaged Delphi warning: t IS defined for n=5!}
      if n>4 then t := exp5(-r) {$ifdef DELPHI}else t:=0.0{$endif};
      case n of
          3: sfc_harmonic2 := 1.0 + a + s;
          4: sfc_harmonic2 := 1.0 + a + s + a*a;
          5: sfc_harmonic2 := 1.0 + a + s + a*a + t;
        else sfc_harmonic2 := 1.0 + a + s + a*a + t + s*a;
      end;
    end;
  end
  else if (r<=1.0) and (r > -MaxInt) and (frac(r)=0.0)  then begin
    {r integer <= 1}
    n := 1 - round(r);
    if n=0 then begin
      {r=1, ordinary harmonic number function}
      sfc_harmonic2 := sfc_harmonic(x);
    end
    else begin
      {http://functions.wolfram.com/06.17.27.0005.01}
      a := sfc_bernoulli(n);
      if n and 1 = 0 then a := -a;
      s := sfc_bernpoly(n, x + 1.0);
      sfc_harmonic2 := (s+a)/n;
    end;
  end
  else if (-1.0 < x) and (x < -0.5) then begin
    {http://functions.wolfram.com/06.17.16.0004.01}
    {This is a critical range with possible cancellation: e.g.}
    {the relative error for H(-0.75, -3.5) is about 3.8E-18.}
    t := x + 1.0;
    s := harm2core(t,r);
    a := power(t,-r);
    sfc_harmonic2 := s - a;
  end
  else if (x > 0.5) and (x <= XREC) then begin
    {http://functions.wolfram.com/06.17.16.0005.01}
    {Make |x| <= 0.5}
    s := 0.0;
    repeat
      s := s + power(x,-r);
      x := x - 1.0;
    until x <= 0.5;
    a := harm2core(x,r);
    sfc_harmonic2 := s + a;
  end
  else begin
    {No special case or recursion}
    sfc_harmonic2 := harm2core(x,r);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_llci(x: extended): extended;
  {-Return the Lobachevski function L(x) = integral(-ln(|cos(t)|), t=0..x)}
const
  lobpi2 = -1.482846963978694993e-4; {lower bits L(Pi)=Pi*ln(2)}
  lbpb22 = -7.414234819893474966e-5; {lower bits L(Pi)/2}
  pi2    = 9.676535897932384626e-4;  {lower bits Pi}
  piby22 = 4.838267948966192313e-4;  {lower bits Pi/2}
const
  lobpi1 : single = 1115.0/512.0;    {2.177734375}  {upper bits L(Pi)}
  pi1    : single = 201.0/64.0;      {3.140625}     {upper bits Pi}
const
  nl1=16;
  arlob1 : array[0..nl1-1] of extended = (
             0.3446488495348130051,
             0.584198357190277669e-2,
             0.19175029694600330e-3,
             0.787251606456769e-5,
             0.36507477415804e-6,
             0.1830287272680e-7,
             0.96890333005e-9,
             0.5339055444e-10,
             0.303408025e-11,
             0.17667875e-12,
             0.1049393e-13,
             0.63359e-15,
             0.3878e-16,
             0.240e-17,
             0.15e-18,
             0.1e-19);
  nl2=11;
  arlob2 : array[0..nl2-1] of extended = (
             2.034594180361328511,
             0.1735185882027407681e-1,
             0.5516280426090521e-4,
             0.39781646276598e-6,
             0.369018028918e-8,
             0.3880409214e-10,
             0.44069698e-12,
             0.527674e-14,
             0.6568e-16,
             0.84e-18,
             0.1e-19);
const
  xlow1 = 0.42576519751e-19;    {double: 8.719671245e-17} {Pi/4*eps/2}
  xlow2 = 0.104125029291e-8;    {double: 4.71216091e-8}   {sqrt(10*eps)}
  xlow3 = 0.139698386192e-8;    {double: 6.32202727e-8}   {sqrt(18*eps)}
  xhigh = 0.92233720368547758e19;  {double: 4.503599627370496e15}  {1/eps}
var
  xr,fval,xcub,t,npi: extended;
  sx: integer;
  gtpi2: boolean;
begin

  {Ref: MISCFUN [22], function LOBACH}
  sx := isign(x);
  xr := abs(x);

  if xr >= xhigh then begin
    {Contribution from reduced x is negligible compared to x*ln(2)}
    {Note that Miscfun returns an error for these cases}
    sfc_llci := ln2*x;
    exit;
  end;

  {Reduce argument to [0,pi]}
  npi := int(xr/pi);
  xr  := (xr - npi*pi1) - npi*pi2;

  {Reduce argument to [0,pi/2]}
  gtpi2 := xr > pi_2;
  if gtpi2 then xr := (pi1 - xr) + pi2;

  {Code for argument in [0,pi/4]}
  if xr <= Pi_4 then begin
    xcub := xr*xr*xr;
    if xr < xlow2 then fval := xcub/SIXX
    else begin
      t := 2.0*sqr(xr/Pi_4) - 1.0;
      t := CSEvalX(t, arlob1, nl1);
      fval := t*xcub;
    end;
  end
  else begin
    {Code for argument in [pi/4,pi/2]}
    xr := (0.5*pi1 - xr) + piby22;
    if xr < xlow3 then begin
      if xr < xlow1 then begin
        {Reduction may produce negative or very small positive values,}
        {the Miscfun test xr=0 does not catch them, e.g. for x=Pi_2   }
        {extended: xr = -2.517276040e-20, double: xr = 6.123224999e-17}
        fval := 0.0;
      end
      else fval := xr*(1.0 - ln(xr))
    end
    else begin
      t := 2.0*sqr(xr/Pi_4) - 1.0;
      t := CSEvalX(t, arlob2, nl2);
      fval := (t - ln(xr))*xr;
    end;
    fval := (0.5*lobpi1 - fval) + lbpb22;
  end;

  {Compute value for argument in [pi/2,pi]}
  if gtpi2 then fval := (lobpi1 - fval) + lobpi2;

  {Add quasi-periodic contribution}
  if npi > 0.0 then fval := (fval + npi*lobpi2 ) + npi*lobpi1;

  {Adjust sign}
  sfc_llci := sx*fval;
end;


{---------------------------------------------------------------------------}
function sfc_llsi(x: extended): extended;
  {-Return the Lobachevski function Lambda(x) = integral(-ln(|2sin(t)|), t=0..x)}
begin
  {http://mathworld.wolfram.com/LobachevskysFunction.html}
  {Lambda(x) = 0.5*cl_2(2x) = Im(polylog(2, exp(2ix)))/2}
  sfc_llsi := 0.5*sfc_cl2(2.0*x);
end;

end.
