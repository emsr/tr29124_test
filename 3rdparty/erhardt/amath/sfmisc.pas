unit sfMisc;

{Common code for special functions: Other functions}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}

{$ifdef NOBASM}
  {$undef BASM}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Common code for special functions: Other functions

 REQUIREMENTS  :  TP5.5-7, D2-D7/D9-D10/D12/D17-D18, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  The unit can be compiled with TP5 but some functions generate
                  invalid operations due to TP5's brain-damaged usage of the FPU

 REFERENCES    :  References used in this unit, main index in amath_info.txt/references

                  [1] [HMF]: M. Abramowitz, I.A. Stegun. Handbook of Mathematical Functions, New York, 1970
                      http://www.math.sfu.ca/~cbm/aands/
                 [22] A.J. MacLeod, MISCFUN: A software package to compute uncommon special functions.
                      ACM Trans. on Math. Soft. 22 (1996), pp.288-301.
                      Fortran source: http://netlib.org/toms/757
                 [23] R.M. Corless, G.H. Gonnet, D.E.G. Hare, D.J. Jeffrey, D.E. Knuth,
                      On the Lambert W Function, Adv. Comput. Math., 5 (1996), pp. 329-359
                      http://www.apmaths.uwo.ca/~rcorless/frames/PAPERS/LambertW/LambertW.ps
                      or http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.33.3583
                 [24] D. Veberic, Having Fun with Lambert W(x) Function, 2010
                      http://arxiv.org/abs/1003.1628
                 [38] T. Ooura's Fortran and C source code for automatic quadrature
                      using Double Exponential transformation; available from
                      http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  17.08.10  W.Ehrhardt  Common version number after split
 1.00.01  08.09.10  we          Improved arg checking and NAN/INF handling

 1.01.00  16.10.10  we          Fix sfc_LambertW1 for WIN16

 1.02.00  05.11.10  we          Tolerance 1.0078125*eps_x in LambertW/1

 1.03.00  03.12.10  we          sfc_polylog
 1.03.01  07.12.10  we          sfc_debye

 1.04.00  17.02.11  we          sfc_ri

 1.06.00  12.05.11  we          sfc_sncndn moved to sfEllInt

 1.07.00  19.06.11  we          Move IsNanOrInf test to top in LambertW1
 1.07.01  26.06.11  we          simplified range check in sfc_LambertW/1

 1.10.00  16.12.11  we          sfc_debye with DE integration for x>4 and n>6
 1.10.01  21.12.11  we          sfc_pz
 1.10.02  22.12.11  we          sfc_pz: Cohen's speed-up $ifdef BASM
 1.10.03  02.01.12  we          improved sfc_polylog for n<-9 and -4 < x < -11.5

 1.11.00  03.02.12  we          Moved to sfZeta: sfc_dilog, sfc_polylog, sfc_pz

 1.12.00  25.02.12  we          Avoid over/underflow for large x in debye_dei

 1.21.00  24.09.13  we          sfc_cosint, sfc_sinint
 1.21.01  26.09.13  we          Fixed quasi-periodic code for sfc_cos/sinint

 1.26.00  24.05.14  we          sfc_ali functional inverse of sfc_li
 1.26.01  01.06.14  we          sfc_fpoly
 1.26.02  02.06.14  we          sfc_lpoly

 1.28.00  17.08.14  we          Catalan function sfc_catf

 1.31.00  09.01.15  we          Euler numbers sfc_euler

 1.32.00  25.01.15  we          sfc_ali moved to sfExpInt

 1.33.00  04.06.15  we          sfc_kepler

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
  MaxEuler = 1866; {Maximum index for Euler numbers for extended}

function sfc_agm(x,y: extended): extended;
  {-Return the arithmetic-geometric mean of |x| and |y|; |x|,|y| < sqrt(MaxExtended)}

function sfc_agm2(x,y: extended; var s: extended): extended;
  {-Return agm(x,y) and s = sum(2^i*c_i^2, i=1...), |x|,|y| < sqrt(MaxExtended)}

function sfc_catf(x: extended): extended;
  {-Return the Catalan function C(x) = binomial(2x,x)/(x+1)}

function sfc_cosint(n: integer; x: extended): extended;
  {-Return cosint(n, x) = integral(cos(t)^n, t=0..x), n>=0}

function sfc_debye(n: integer; x: extended): extended;
  {-Return the Debye function D(n,x) = n/x^n*integral(t^n/(exp(t)-1),t=0..x) of order n>0, x>=0}

function sfc_euler(n: integer): extended;
  {-Return the nth Euler number, 0 if n<0 or odd n}

function sfc_fpoly(n: integer; x: extended): extended;
  {-Return the Fibonacci polynomial F_n(x)}

function sfc_lpoly(n: integer; x: extended): extended;
  {-Return the Lucas polynomial L_n(x)}

function sfc_LambertW(x: extended): extended;
  {-Return the Lambert W function (principal branch), x >= -1/e}

function sfc_LambertW1(x: extended): extended;
  {-Return the Lambert W function (-1 branch), -1/e <= x < 0}

function sfc_kepler(M,e: extended): extended;
  {-Solve Kepler's equation, result x is the eccentric anomaly from the mean anomaly M and the }
  { eccentricity e >= 0; x - e*sin(x) = M, x + x^3/3 = M, or e*sinh(x) - x = M for e <1, =1, >1}

function sfc_ri(x: extended): extended;
  {-Return the Riemann prime counting function R(x), x >= 1/16}

function sfc_sinint(n: integer; x: extended): extended;
  {-Return sinint(n, x) = integral(sin(t)^n, t=0..x), n>=0}

implementation

uses
  AMath,
  sfBasic, sfGamma, sfZeta, sfExpInt;

{---------------------------------------------------------------------------}
function sfc_agm(x,y: extended): extended;
  {-Return the arithmetic-geometric mean of |x| and |y|; |x|,|y| < sqrt(MaxExtended)}
var
  a,b,t: extended;
begin
  {Ref: HMF [1], M. Abramowitz, I. Stegun; Handbook of Mathematical Functions}
  {     Chapter 17: Elliptic Integrals, Sec. 17.6: Arithmetic-Geometric Mean }
  {a := max(|x|,|y|), b := min(|x|,|y|)}
  if abs(x)<abs(y) then begin
    b := abs(x);
    a := abs(y);
  end
  else begin
    a := abs(x);
    b := abs(y);
  end;
  if (a=0.0) or (b=0.0) then sfc_agm := 0.0
  else if a=b then sfc_agm := a
  else begin
    {Invariant a >= b. Since agm is quadratically convergent, stop}
    {the iteration if the relative error is less than ~sqrt(eps_x)}
    while a-b > 2e-10*a do begin
      t := a;
      a := 0.5*(t+b);
      b := sqrt(t*b);
    end;
    {Final step to make relative error less then eps_x}
    sfc_agm := 0.5*(a+b);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_agm2(x,y: extended; var s: extended): extended;
  {-Return agm(x,y) and s = sum(2^i*c_i^2, i=1...), |x|,|y| < sqrt(MaxExtended)}
var
  a,b,c,t: extended;
  i: integer;
begin
  s := 0;
  {a := max(|x|,|y|), b := min(|x|,|y|)}
  if abs(x)<abs(y) then begin
    b := abs(x);
    a := abs(y);
  end
  else begin
    a := abs(x);
    b := abs(y);
  end;
  if (a=0.0) or (b=0.0) then sfc_agm2 := 0.0
  else if a=b then sfc_agm2 := a
  else begin
    i := 0;
    {Invariant a >= b. Since agm is quadratically convergent, stop}
    {the iteration if the relative error is less than ~sqrt(eps_x)}
    repeat
      inc(i);
      c := 0.5*(a-b);
      s := s + ldexp(c*c,i);
      t := a;
      a := 0.5*(t+b);
      b := sqrt(t*b);
      {done if the PREVIOUS relative error is < ~sqrt(eps_x)}
    until c < 2e-10*t;
    sfc_agm2 := a;
  end;
end;


{---------------------------------------------------------------------------}
function debye_dei(n: integer; x: extended): extended;
  {-Debye function D(n,x) using double exponential integration}
const
  mmax = 256;
  eps  = 1e-17;
  efs  = 0.1;
  hoff = 8.5;
var
  epsln,epsh,h0,ehp,ehm,epst,ir,iv,h,iback,
  irback,t,ep,em,xw,xa,wg,fa,fb,err,errt,errh,errd: extended;
var
  m: integer;

  function dnf(t: extended): extended;
    {-The Debye integrand}
  var
    p: extended;
  begin
    {Avoid over/underflow for large x}
    if t > ln_MaxExt then dnf := 0.0
    else begin
      p  := power(t/x,n);
      if t >= 44.0 then dnf := p*exp(-t)
      else dnf := p/expm1(t);
    end;
  end;

begin
  {This is the customised AMTools' function intdee with eps=1e-17}
  {for the integral((t/x)^n/(exp(t)-1),t=0..x), Ref [38]}
  debye_dei := 0;
  if x=0 then exit;

  epsln := 1.0 - ln(efs*eps);
  epsh  := sqrt(efs*eps);
  h0    := hoff/epsln;
  ehp   := exp(h0);
  ehm   := 1.0/ehp;
  epst  := exp(-ehm*epsln);
  ir    := dnf(0.5*x);
  ir    := 0.25*ir*x;
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
        xa := x*xw;
        wg := xa*(1.0-xw);
        fa := dnf(xa)*wg;
        fb := dnf(x-xa)*wg;
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
  debye_dei := h*iv*n;
end;


{---------------------------------------------------------------------------}
function sfc_debye(n: integer; x: extended): extended;
  {-Return the Debye function D(n,x) = n/x^n*integral(t^n/(exp(t)-1),t=0..x) of order n>0, x>=0}
var
  p,s,t,q,xk,sk: extended;
  j,k: integer;
begin
  if (x<0.0) or (n<1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_debye := NaN_x;
    exit;
  end;
  j := n+1;
  if x <= 4.0 then begin
    {HMF [1], 27.1.1; the series converges for |x| < 2Pi}
    s := 1.0 - 0.5*n*x/j;
    q := x*x;
    k := 2;
    p := n;
    repeat
      t := sfc_bernoulli(k);
      p := p*q/k/(k-1);
      t := t/(k+n)*p;
      s := s + t;
      inc(k,2);
    until abs(t)<eps_x*abs(s);
    sfc_debye := s;
  end
  else begin
    if n<=6 then begin
      {See HMF [1], 27.1.2 and 27.1.3. See also  A.J. MacLeod's MISCFUN [22]  }
      {and R.J. Mathar's http://www.strw.leidenuniv.nl/~mathar/progs/debye3.c }

      {HMF 27.1.2 gives the complementary integral; this function is computed }
      {by subtracting the HMF partial sums from the complete integral 27.1.3. }
      {This is the main source of the decreasing accuracy due to cancellation.}

      s := sfc_fac(n);
      t := sfc_zetaint(j);
      p := power(x,-j);
      s := (s*p)*t;
      if (s <> 0.0) and (x < -ln_MinExt) then begin
        q := 0.25*eps_x;  {tolerance}
        k := 1;
        repeat
          xk := x*k;
          sk := 1.0/xk;
          t  := sk;
          for j:=n downto 1 do begin
            t  := t*j/xk;
            sk := sk + t;
          end;
          t := sk*exp(-xk);
          s := s - t;
          inc(k);
        until (t<q*s) or (s<=0.0);
        {Correct rounding errors}
        if s<0.0 then s := 0.0;
      end;
      sfc_debye := s*x*n;
    end
    else begin
      {n>6 and x>4, use double exponential integration}
      sfc_debye := debye_dei(n,x)
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_LambertW(x: extended): extended;
  {-Return the Lambert W function (principal branch), x >= -1/e}
const
  em1h = 0.36785888671875;           {high part of 1/e}
  em1l = 0.205544526923215955238e-4; {low  part of 1/e}
  xlow =-0.36767578125;              {threshold for branch-point series}
const
  PCHex: array[0..11] of THexExtW = ( {Calculated with MPArith:}
           {Branch point series W(x) = sum(a[i]*y^i, i=0..), y=sqrt(x+1/e),}
           {a[i] = mu[i]*(sqrt(2e))^i, with mu[i] from [23] (4.23/4).}
           ($0000,$0000,$0000,$8000,$BFFF),   {-1.0000000000000000000}
           ($F0C8,$B7FD,$A7AD,$9539,$4000),   { 2.3316439815971242034}
           ($6379,$83A4,$C5CB,$E7F5,$BFFF),   {-1.8121878856393634902}
           ($5E29,$9492,$8742,$F7E3,$3FFF),   { 1.9366311144923597554}
           ($F892,$62C0,$9538,$96A0,$C000),   {-2.3535512018816145169}
           ($AAFC,$5D5F,$6A8E,$C447,$4000),   { 3.0668589010506319128}
           ($699C,$9DDE,$5967,$859C,$C001),   {-4.1753356002581771388}
           ($87F2,$6021,$EE2E,$BB74,$4001),   { 5.8580237298747741487}
           ($842F,$2E86,$A0C2,$866A,$C002),   {-8.4010322175239773709}
           ($B325,$780C,$161A,$C403,$4002),   { 12.250753501314460424}
           ($0BB1,$3AFE,$3A3C,$90CE,$C003),   {-18.100697012472442756}
           ($3C33,$EEF4,$7BD6,$D83B,$4003)(*, { 27.029044799010561650}
           ($526B,$4BEF,$A248,$A2DC,$C004),   {-40.715462808260627286}
           ($F1BB,$4D82,$A26E,$F721,$4004),   { 61.782846187096525741}
           ($802C,$57E8,$5D3D,$BCAC,$C005),   {-94.336648861866933963}
           ($B2E6,$5CD6,$39F1,$90D1,$4006)*));{ 144.81729038731164003}
var
  PC : array[0..11] of extended absolute PCHex;
const
  p1 = 19.0/10.0;
  p2 = 17.0/60.0;
  q1 = 29.0/10.0;
  q2 = 101.0/60.0;
const
  IMAX = 10;
  eps  = 1e-16;
var
  e,q,w,z: extended;
  i: integer;
begin
  if IsNanOrInf(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_LambertW  := NaN_x;
    exit;
  end;

  if abs(x)<=1e-10 then begin
    {LambertW = x*(1-x+1.5*x^2+O(x^3)}
    sfc_LambertW := x*(1.0-x);
  end
  else if x <= xlow then begin
    z := (x+em1h)+em1l;
    if abs(z) <= 1.0078125*eps_x then begin
      {Tolerate some rounding error, this will give W(-exp(-1)) = -1}
      sfc_LambertW := -1.0;
    end
    else begin
      if z<0.0 then begin
        {$ifopt R+}
          if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
        {$endif}
        sfc_LambertW  := NaN_x;
        exit;
      end;
      sfc_LambertW := PolEvalX(sqrt(z),PC,12);
    end;
  end
  else begin
    {Initial approximations for iteration}
    if x<3.0 then begin
      {pade(LambertW(x), x, [3,2])}
      q := (p2*x+p1)*x + 1.0;
      z := (q2*x+q1)*x + 1.0;
      w := x*(q/z);
    end
    else begin
      {Asymptotic series, see [23] (4.19) or [24] (40)}
      z := ln(x);
      q := ln(z);
      w := z - q + q/z;
      if q>3.0 then w := w + 0.5*q*(q-2.0)/sqr(z);
    end;
    {Fritsch iteration, see [24], section 2.3}
    for i:=1 to IMAX do begin
      z := ln(x/w) - w;
      e := 1.0 + w;
      q := 2.0*e*(e + z/1.5) - z;
      e := (z/e)*(q/(q-z));
      w := w*(1.0+e);
      if abs(e) <= eps then begin
        sfc_LambertW := w;
        exit;
      end;
    end;
    sfc_LambertW := w;
    if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  end;
end;


{---------------------------------------------------------------------------}
function sfc_LambertW1(x: extended): extended;
  {-Return the Lambert W function (-1 branch), -1/e <= x < 0}
const
  IMAX = 10;
const
  eh = 2.71826171875;            {high part of e}
  el = 0.2010970904523536029e-4; {low  part of e}
const
  a0 = -3.069497752294975357;
  a1 = -0.084452335662190688;
  a2 = 12.74914818472526989;
  b0 =  0.526928257471981805;
  b1 = -3.528391962619355995;
  b2 = -5.815166249896325102;
const
  PCHex: array[0..9] of THexExtW = (         {Branch-point series coefficients}
           ($0000,$0000,$0000,$8000,$BFFF),  {-1                    = -1.0000000000000000000   }
           ($0000,$0000,$0000,$8000,$3FFF),  { 1                    =  1.0000000000000000000   }
           ($AAAB,$AAAA,$AAAA,$AAAA,$BFFD),  {-1/3                  = -3.3333333333333333334E-1}
           ($1C72,$71C7,$C71C,$9C71,$3FFC),  {11/72                 =  1.5277777777777777778E-1}
           ($4DC0,$6A31,$DBF8,$A314,$BFFB),  {-43/540               = -7.9629629629629629633E-2}
           ($80F3,$9D64,$0F2B,$B648,$3FFA),  {769/17280             =  4.4502314814814814816E-2}
           ($FACB,$7F14,$E592,$D4DD,$BFF9),  {-221/8505             = -2.5984714873603762493E-2}
           ($1AD0,$1DE2,$4C4C,$8016,$3FF9),  {680863/43545600       =  1.5635632532333921222E-2}
           ($86C7,$5B95,$3D08,$FBB7,$BFF4),  {-196/204120           = -9.6021947873799725651E-4}
           ($39E2,$3187,$A549,$C515,$3FF7)); {226287557/37623398400 =  6.0145432529561178610E-3}
var
  PC: array[0..9] of extended absolute PCHex;
var
  z,e,t,w: extended;
  i: integer;
  done: boolean;
begin
  if IsNanOrInf(x) or (x>=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_LambertW1 := NaN_x;
    exit;
  end;
  {For non-BASM Pascal versions the system ln/exp functions are too}
  {inaccurate for iterations at arguments very close to -1/e and 0.}
  if x<-0.275 then begin
    z := (eh*x + el*x) + 1.0;
    if abs(z) <= 1.0078125*eps_x then begin
      {Tolerate some rounding error, this will give W_1(-exp(-1)) = -1}
      sfc_LambertW1 := -1.0;
      exit;
    end;
    {Branch-point series expansion about -1/e}
    if z<0.0 then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfc_LambertW1 := NaN_x;
      exit;
    end;
    z := -sqrt(2.0*z);
    w := PolEvalX(z,PC,8);
    {One Newton step}
    e := exp(w);
    t := w*e;
    z := (t-x)/(t+w);
    w := w - z;
  end
  else if x<-0.125 then begin
    {Pade approximation: expand(pade(LambertW(-1,x),x=-0.2,[2,2]));}
    z := sqr(x);
    w := (a0+a1*x+a2*z) / (b0+b1*x+b2*z);
  end
  else begin
    {Asymptotic series near zero, see [23] (4.19) or [24] (40)}
    z := ln(-x);
    t := ln(-z);
    w := z - t + t/z;
    w := w + 0.5*t*(t-2.0)/sqr(z);
    if x > {$ifdef BASM}-MinExtended{$else}-5.5E-4928{$endif} then begin
      {x denormal/small, use another term instead of iteration}
      sfc_LambertW1 := w + t*(6.0-9.0*t+2.0*sqr(t))/(6.*z*sqr(z));
      exit;
    end;
  end;
  done := false;
  for i:=1 to IMAX do begin
    {Halley iteration, see [24] section 2.2 or [23] formula (5.9)}
    e := exp(w);
    t := w*e-x;
    if t<>0.0 then begin
      z := w+1.0;
      {For WIN16 it is possible that t/e raises division by zero!}
      e := (e*z-0.5*(z+1.0)*t/z);
      if e<>0.0 then t := t/e;
      w := w-t;
    end
    else done := true;
    if done then begin
      sfc_LambertW1 := w;
      exit;
    end
    else done := abs(t)<1e-16*(abs(w))
  end;
  sfc_LambertW1 := w;
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
end;


{---------------------------------------------------------------------------}
function sfc_ri(x: extended): extended;
  {-Return the Riemann prime counting function R(x), x >= 1/16}
var
  f,s,t,y: extended;
  n: integer;
begin
  if IsNanOrInf(x) or (x < 0.0625) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_ri := NaN_x;
    exit;
  end;
  {R(x) = sum(æ(n)/n*li(x^(1/n))), n=1..Inf)}
  if x>=1e40 then sfc_ri := sfc_li(x) {only one significant term}
  else if x>=1e19 then begin
    {If x is large use R(x) = li(x) - 0.5*li(sqrt(x)) ...}
    s := sfc_li(x);
    t := 0.0;
    for n:=2 to 35 do begin
      if moebius[n] <> 0 then begin
        y := sfc_li(nroot(x,n))/n;
        if moebius[n]>0 then s := s + y
        else t := t + y;
        if y <= eps_x*s then begin
          sfc_ri := s-t;
          exit;
        end;
      end;
    end;
    {No fast convergence, but the maximum relative truncation error }
    {is about 4.5e-19 and this result is more accurate than the Gram}
    {series result for x >= 1e19.}
    sfc_ri := s-t;
  end
  else begin
    {For smaller x >= 1/16 use the Gram series for R(x)}
    {R(x) = 1 + sum(ln(x)^n/(n!*n*zeta(n+1)), n=1..Inf)}
    y := ln(x);
    s := 0.0;
    f := 1.0;
    n := 1;
    repeat
      f := f*(y/n);
      t := n*sfc_zetaint(n+1);
      t := f/t;
      s := s+t;
      inc(n);
    until abs(t) <= eps_x*abs(s);
    sfc_ri := 1.0+s;
  end;
end;


{---------------------------------------------------------------------------}
{------------------- Integrals of sin/cos powers ---------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function cosint1(n: integer; x: extended): extended;
  {-Compute cosint(n,x) using recurrence formula; intern n>1}
var
  s,c,c2,pk,ck: extended;
  k: integer;
begin
  {$ifdef debug}
    if n<2 then RunError(byte(RTE_ArgumentRange));
  {$endif}
  if x=0.0 then cosint1 := x
  else begin
    sincos(x,s,c);
    c2 := c*c;
    if odd(n) then begin
      {Sum terms for k=1,3, .., n}
      {Setup pk = s*c^(k+1) and ck = sfc_cosint(k,x) for k=1}
      k  := 1;
      pk := s*c2;
      ck := s;
    end
    else begin
      {Sum terms for k=0,2, .., n}
      {Setup pk = s*c^(k+1) and ck = cosint(k,x) for k=0}
      k  := 0;
      pk := s*c;
      ck := x;
    end;
    while k<n do begin
      k  := k+2;
      {compute cosint(k,x) from cosint(k-2,x)}
      ck := ((k-1)*ck + pk)/k;
      {update s*c^(k-1)}
      pk := pk*c2;
    end;
    cosint1 := ck;
  end;
end;


{---------------------------------------------------------------------------}
function cipi2(n: integer): extended;
  {-Return cosint(n,Pi/2)}
begin
  cipi2 := 0.5*SqrtPi*sfc_gamma_delta_ratio(0.5*(n+1),0.5)
end;


{---------------------------------------------------------------------------}
function sfc_cosint(n: integer; x: extended): extended;
  {-Return cosint(n, x) = integral(cos(t)^n, t=0..x), n>=0}
var
  a,y,z: extended;
begin
  if IsNanOrInf(x) or (n<0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_cosint := NaN_x;
    exit;
  end;
  if n=0 then sfc_cosint := x
  else if n=1 then sfc_cosint := sin(x)
  else if x=0.0 then sfc_cosint := 0
  else if odd(n) then begin
    z := rem_2pi_sym(x);
    sfc_cosint := cosint1(n,z);
  end
  else begin
    {cosint is quasi-periodic mod Pi}
    y := abs(x);
    {$ifdef ExtendedSyntax_on}
      rem_pio2(0.5*y,z);
    {$else}
      if 0=rem_pio2(0.5*y,z) then;
    {$endif}
    z := 2.0*z;
    a := int((y-z)/Pi + 0.5);
    z := cosint1(n,z);
    if a<>0.0 then z := z + 2.0*a*cipi2(n);
    if x<0.0 then z := -z;
    sfc_cosint := z;
  end;
end;


{---------------------------------------------------------------------------}
function sinintser(n: integer; x: extended): extended;
  {-Compute sinint(n,x) using series, |x| <~ 1.25}
var
  a,s,f,y,t,p: extended;
  k: integer;
const
  kmax = 800;
begin
  y := sin(x);
  s := 1.0/(1.0+n);
  p := power(y,n+1);
  if p<>0.0 then begin
    a := 1.0;
    f := 1.0;
    y := y*y;
    k := 2;
    repeat
      a := ((k-1)/k)*a;
      f := f*y;
      t := (a*f)/(k+n+1);
      s := s + t;
      k := k + 2;
    until (t < eps_x*s) or (k > kmax);
    {$ifdef debug}
      if k>kmax then sfc_write_debug_str('** No convergence in sinintser');
    {$endif}
  end;
  sinintser := s*p;
end;


{---------------------------------------------------------------------------}
function sinint1(n: integer; x: extended): extended;
  {-Return sinint(n, x) for abs(x) <= Pi; intern n>1}
var
  s1,s2,z: extended;
const
  z0 = 1.0;
begin
  {$ifdef debug}
    if n<2 then RunError(byte(RTE_ArgumentRange));
  {$endif}
  z := abs(x);
  if z <= z0 then begin
    if z=0.0 then s2 := 0.0
    else s2 := sinintser(n,z);
  end
  else begin
    s1 := cosint1(n,Pi_2-z);
    s2 := cipi2(n);
    s2 := s2-s1;
  end;
  if (x<0.0) and (n and 1 =0) then s2 := -s2;
  sinint1 := s2;
end;


{---------------------------------------------------------------------------}
function sfc_sinint(n: integer; x: extended): extended;
  {-Return sinint(n, x) = integral(sin(t)^n, t=0..x), n>=0}
var
  a,y,z: extended;
begin
  if IsNanOrInf(x) or (n<0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_sinint := NaN_x;
    exit;
  end;
  if n=0 then sfc_sinint := x
  else if n=1 then sfc_sinint := vers(x)
  else if x=0 then sfc_sinint := 0.0
  else if odd(n) then begin
    {sinint is even, and periodic mod 2Pi}
    z := abs(rem_2pi_sym(x));
    sfc_sinint := sinint1(n,z);
  end
  else begin
    {sinint is quasi-periodic mod Pi}
    y := abs(x);
    {$ifdef ExtendedSyntax_on}
      rem_pio2(0.5*y,z);
    {$else}
      if 0=rem_pio2(0.5*y,z) then;
    {$endif}
    z := 2.0*z;
    a := int((y-z)/Pi + 0.5);
    z := sinint1(n,z);
    if a<>0.0 then z := z + 2.0*a*cipi2(n);
    if x<0.0 then z := -z;
    sfc_sinint := z;
  end;
end;


{---------------------------------------------------------------------------}
procedure fibmatpow(n: integer; x: extended; var f1,f2: extended);
  {-Computes [x,1; 1,0]^n and returns the first column as f1,f2}
  { assumes n > 2, x > 0; i.e. f1 = F(n-1,x), f2 = F(n-2,x)}
const
  imax = maxint xor (maxint shr 1);  {2^14 or 2^30 depending on sizeof(integer)}
var
  i: integer;
  x11,x12,x22,tmp: extended;
begin
  {Use a fast O(log(n)) algorithm: compute the left-to-right binary matrix  }
  {power X^(n-1) with X = [x,1; 1,0]. f1, f2 are the (1,1) and (1,2) compo- }
  {nents. During the algorithm the powers X^k remain symmetric, and x21 is  }
  {not used. The loop uses about 4.5*log2(n) multiplications on the average.}
  {This is better then iteration for n>31, but it is slightly less accurate.}

  {get highest bit of |n|-1}
  i := imax;
  while n and i = 0 do i := i shr 1;
  i := i shr 1;

  {X = [x,1; 1,0]}
  x11 := x;
  x22 := 0.0;
  x12 := 1.0;

  while i > 0 do begin
    {X = X*X}
    tmp := sqr(x12);
    x12 := (x11 + x22)*x12;
    x11 := sqr(x11) + tmp;
    x22 := sqr(x22) + tmp;
    if n and i <> 0 then begin
      {X = X*[x,1; 1,0]}
      x22 := x12;
      x12 := x11;
      x11 := x*x12 + x22;
    end;
    i := i shr 1;
  end;
  f1 := x11;
  f2 := x12;
end;


{---------------------------------------------------------------------------}
function sfc_fpoly(n: integer; x: extended): extended;
  {-Return the Fibonacci polynomial F_n(x)}
var
  i,m: integer;
  f1,f2,t,z: extended;
  neg: boolean;
begin
  {32-bit integer note: |n|>2^15 requires |x| <~ sqrt(2)/2}
  m := abs(n);
  z := abs(x);

  {Prepare sign change flag}
  neg := (n and 1 = 0) and ((x<0.0) <> (n<0));

  {f1 will be F(n,x)}
  if z=0.0 then begin
    if odd(m) then f1 := 1.0
    else f1 := 0.0;
  end
  else if m < 32 then begin
    if m <= 0 then f1 := 0.0
    else if m=1 then f1 := 1.0
    else begin
      {Use iteration from definition F(n,x) = x*F(n-1,x) + F(n-2,x)}
      f1 := z;
      f2 := 1.0;
      {This uses m-2 multiplications}
      for i:=3 to m do begin
        t  := z*f1 + f2;
        f2 := f1;
        f1 := t;
      end;
    end;
  end
  else begin
    {Use fast matrix algorithm}
    fibmatpow(m-1,z,f1,f2);
  end;

  {here f1 = sfc_fpoly(|n|, |x|), adjust sign from signs of n and x}
  if (f1<>0.0) and neg then sfc_fpoly := -f1
  else sfc_fpoly := f1;
end;


{---------------------------------------------------------------------------}
function sfc_lpoly(n: integer; x: extended): extended;
  {-Return the Lucas polynomial L_n(x)}
var
  i,m: integer;
  l1,l2,t,z: extended;
  neg: boolean;
begin
  m := abs(n);
  z := abs(x);
  {Prepare sign change flag}
  neg := odd(n) and ((x<0.0) <> (n<0));
  if z=0.0 then begin
    if odd(m) then l1 := 0
    else l1 := 2.0;
  end
  else if m < 32 then begin
    if m <= 0 then l1 := 2.0
    else if m=1 then l1 := z
    else begin
      {Use iteration from definition L(n,x) = x*F(n-1,x) + F(n-2,x)}
      l1 := z;
      l2 := 2.0;
      for i:=2 to m do begin
        t  := z*l1 + l2;
        l2 := l1;
        l1 := t;
      end;
    end;
  end
  else begin
    {Use fast matrix algorithm and L(m,z) = z*F(m,z) + 2*F(m-1,z)}
    fibmatpow(m-1,z,l1,l2);
    l1 := z*l1 + 2.0*l2;
  end;
  {here l1 = sfc_lpoly(|n|, |x|), adjust sign from signs of n and x}
  if (l1<>0.0) and neg then sfc_lpoly := -l1
  else sfc_lpoly := l1;
end;


{---------------------------------------------------------------------------}
function sfc_catf(x: extended): extended;
  {-Return the Catalan function C(x) = binomial(2x,x)/(x+1)}
var
  f,z: extended;
const
  x2 = 8202.16431466394231; {double 519.17950298281473}
  x1 = 8192.0;              {double 512}
  x0 = 35.0;                {double 32}
begin
  {C(x) = 4^x*Gamma(x+1/2)/[sqrt(Pi)*Gamma(x+2)]  }
  {     = 4^x*gamma_delta_ratio(x+1/2, 3/2)/SqrtPi}
  f := frac(x);
  if x>x2 then sfc_catf := PosInf_x
  else if (x<0.0) and ((f=0.0) or (f=-0.5)) then begin
    if x=-1.0 then sfc_catf := -0.5
    else if f=0.0 then sfc_catf := 0.0  {Infinity in denominator}
    else sfc_catf := PosInf_x;          {Infinity in numerator}
  end
  else begin
    z := sfc_gamma_delta_ratio(x+0.5, 1.5)/SqrtPi;
    if f=0.0 then begin
      {Here x is an integer >= 0, C(x) is also an integer}
      {'round' result if there may be a fractional part. }
      z := ldexp(z, 2*trunc(x));
      if x > x0 then sfc_catf := z
      else sfc_catf := int(z+0.125);
    end
    else begin
      if x < x1 then sfc_catf := exp2(2.0*x)*z
      else begin
        {exp2(2x) would overflow}
        f := exp2(2.0*f);
        z := ldexp(z, 2*trunc(x));
        sfc_catf := f*z;
      end;
    end
  end;
end;


{---------------------------------------------------------------------------}
function sfc_euler(n: integer): extended;
  {-Return the nth Euler number, 0 if n<0 or odd n}
const
  enhex: array[0..13] of THexExtw = (
           ($0000,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000000}
           ($0000,$0000,$0000,$8000,$BFFF),  {-1.0000000000000000000}
           ($0000,$0000,$0000,$A000,$4001),  {+5.0000000000000000000}
           ($0000,$0000,$0000,$F400,$C004),  {-6.1000000000000000000E+1}
           ($0000,$0000,$0000,$AD20,$4009),  {+1.3850000000000000000E+3}
           ($0000,$0000,$0000,$C559,$C00E),  {-5.0521000000000000000E+4}
           ($0000,$0000,$B400,$A4F6,$4014),  {+2.7027650000000000000E+6}
           ($0000,$0000,$1D50,$BE20,$C01A),  {-1.9936098100000000000E+8}
           ($0000,$2000,$5FCA,$907A,$4021),  {+1.9391512145000000000E+10}
           ($0000,$0C40,$7FEC,$8BFB,$C028),  {-2.4048796754410000000E+12}
           ($8000,$126A,$E18E,$A86C,$402F),  {+3.7037118823752500000E+14}
           ($ED00,$C6B2,$6F0F,$F660,$C036),  {-6.9348874393137901000E+16}
           ($76B9,$4A7B,$B29F,$D74E,$403E),  {+1.5514534163557086905E+19}
           ($0C75,$F6C1,$8668,$DD8F,$C046)); {-4.0870725092931238925E+21}
const
  nz = 63;
  nd = 38;
var
  b,d,z: extended;
begin
  if odd(n) or (n<0) then sfc_euler := 0.0
  else if n<=26 then sfc_euler := extended(enhex[n div 2])
  else if n>MaxEuler then sfc_euler := PosInf_x
  else begin
    {E_n/B_n = -2/Pi*4^n*DirichletBeta(n+1)/Zeta(n)}
    b := -sfc_bernoulli(n);
    if n <= nz then begin
      if n > nd then d := 1.0
      else begin
        d := intpower(3.0,n+1);
        d := (d-1.0)/d;
      end;
      z := sfc_zetaint(n);
      b := (b/z)*d;
    end;
    sfc_euler := ldexp(b/Pi_2, n+n);
  end
end;


{---------------------------------------------------------------------------}
function kepler_elliptic(M,e: extended): extended;
  {-Solve x - e*sin(x) = |M|, 0 <= e < 1}
var
  x,mr,s,d,f,f1,f2: extended;
  ic: integer;
const
  imax = 30;
  eps  = 5.8e-11; {~ sqrt_epsh/4}
begin
  M  := abs(M);
  mr := rem_2pi(M);
  ic := 0;
  if mr<0.1 then x := mr + sqrt(mr)*e
  else x := mr + 0.8*e;
  if (e>0.9) or (x>Pi) then x := Pi;
  repeat
    inc(ic);
    if ic > imax then begin
      if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
      kepler_elliptic := Nan_x;
      exit;
    end;
    sincos(x,s,f1);
    f  := x - e*s - mr;
    if f=0.0 then d := 0.0
    else begin
      {Halley iteration step}
      f1 := 1.0 - e*f1;
      d  := f/f1;
      f2 := e*s;
      s  := 1.0 - 0.5*d*f2/f1;
      d  := d/s;
      x  := x - d;
    end;
  until abs(d) <= eps*abs(x);
  kepler_elliptic := x+(M-mr);
end;


{---------------------------------------------------------------------------}
function kepler_hyperbolic(M,e: extended): extended;
  {-Solve e*sinh(x) - x = |M|, e > 1}
var
  x,s,d,f,f1,f2: extended;
  ic: integer;
const
  imax = 30;
  eps  = 5.8e-11;     {~ sqrt_epsh/4 }
  m1   = 0.397e4932;  {<MaxExtended/3}
begin
  {assumes M >= 0, e > 1}
  M  := abs(M);
  if M>m1 then  x := ln(M/e) + ln2
  else x  := ln(2.0*M/e + 1.8);
  ic := 0;
  repeat
    inc(ic);
    if ic > imax then begin
      if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
      kepler_hyperbolic := Nan_x;
      exit;
    end;
    sinhcosh(x,s,f1);
    f := e*s - x - M;
    if f=0.0 then d := 0.0
    else begin
      {Halley iteration step}
      f1 := e*f1 - 1.0;
      d  := f/f1;
      f2 := e*s;
      s  := 1.0 - 0.5*d*f2/f1;
      d  := d/s;
      x  := x - d;
    end;
  until abs(d) <= eps*abs(x);
  kepler_hyperbolic := x;
end;


{---------------------------------------------------------------------------}
function kepler_parabolic(M: extended): extended;
  {-Solve Barker's equation x+x^3/3 = |M|}
var
  w,s: extended;
const
  w0 = 1.125e30;    {~40 eps(-3/2) }
  w1 = 0.397e4932;  {<MaxExtended/3}
  cbrt3 = 1.4375 + 0.4749570307408382322e-2;
begin
  {Ref: T. Fukushima, Theoretical Astrometry, Draft March 19, 2003 }
  {Ch. 8.8 Analytic Solution of Barker's Equation, formula (8.209) }
  {http://chiron.mtk.nao.ac.jp/~toshio/education/main.pdf          }
  w := abs(M);
  if w>=w0 then begin
    if w<w1 then kepler_parabolic := cbrt(3.0*M)
    else kepler_parabolic := cbrt3*cbrt(M);
  end
  else if w<=sqrt_epsh then kepler_parabolic := w
  else begin
    w := 1.5*w;
    s := hypot(1.0,w);
    s := sqr(cbrt(w+s));
    w := 2.0*s*w;
    s := 1.0 + s + s*s;
    kepler_parabolic := w/s;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_kepler(M,e: extended): extended;
  {-Solve Kepler's equation, result x is the eccentric anomaly from the mean anomaly M and the }
  { eccentricity e >= 0; x - e*sin(x) = M, x + x^3/3 = M, or e*sinh(x) - x = M for e <1, =1, >1}
var
  x: extended;
begin
  if IsNanOrInf(M) or IsNanOrInf(e) or (e<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_kepler := NaN_x;
    exit;
  end;
  if e=0.0 then sfc_kepler := M
  else begin
    if e=1.0 then x := kepler_parabolic(M)
    else if e<1.0 then x := kepler_elliptic(M, e)
    else x := kepler_hyperbolic(M, e);
    if (M<0.0) and (x<>0.0) then x := -x;
    sfc_kepler := x;
  end;
end;



end.
