unit sfMisc;

{Common code for special functions: Other functions}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Common code for special functions: Other functions

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

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
                      https://arxiv.org/abs/1003.1628
                 [30] [NIST]: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST
                      Handbook of Mathematical Functions, Cambridge, 2010. Online
                      resource: NIST Digital Library of Mathematical Functions,
                      http://dlmf.nist.gov/
                 [38] T. Ooura's Fortran and C source code for automatic quadrature
                      using Double Exponential transformation; available from
                      http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html
                 [71] R.T. Short, Computation of Rice and Noncentral Chi-Squared Probabilities, 2012.
                      Available as http://www.phaselockedsystems.com/NoncentralChiSquared.pdf
                 [74] H.H. Chan, On Ramanujan's cubic continued fraction,
                      Acta Arithmetica LXXIII.4 (1995), available from
                      http://www.math.nus.edu.sg/~chanhh/papers/4.pdf
                 [77] H.L. Krall, O. Frink, A new class of orthogonal polynomials:
                      the Bessel polynomials, Trans. Amer. Math. Soc. 65 (1949) 100-115.
                      https://www.ams.org/tran/1949-065-01/S0002-9947-1949-0028473-1/S0002-9947-1949-0028473-1.pdf
                 [79] W.J. Cody, K.E. Hillstrom, Chebyshev approximations for the
                      Coulomb phase shift, Math. Comp. 24 (1970), 671-677,
                      https://doi.org/10.1090/S0025-5718-1970-0273785-4

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

 1.38.00  17.06.16  we          Bring radical sfc_br
 1.38.01  20.06.16  we          Improved sfc_br
 1.38.02  25.06.16  we          Wright omega sfc_wo
 1.38.04  26.06.16  we          Nan/Inf handling for sfc_wo
 1.38.05  27.06.16  we          Nan/Inf handling for sfc_br

 1.39.00  28.09.16  we          General Fibonacci function sfc_fibf

 1.40.00  25.05.17  we          sfc_eulerq
 1.40.01  26.05.17  we          Improved AE for sfc_eulerq if q ~ -1
 1.40.02  20.06.17  we          sfc_detai
 1.40.03  21.06.17  we          Nan/Inf handling for sfc_detai
 1.40.04  28.06.17  we          EulerQ_aen: Anti FPC -O4 code
 1.40.05  29.06.17  we          Code removed for TP5-TP6, TPW1-D1

 1.43.00  08.12.17  we          sfc_mq: generalized Marcum Q function

 1.44.00  04.01.18  we          Improved sfc_ri
 1.44.01  10.01.18  we          sfc_ari: inverse of RiemannR
 1.44.02  29.01.18  we          sfc_einstein
 1.44.03  01.02.18  we          rewrite sfc_debye (always integrate if x>4}
 1.44.04  02.02.18  we          Transport integrals sfc_trint

 1.45.00  20.02.18  we          sfc_langevin: Langevin function from AMath
 1.45.01  21.02.18  we          sfc_llinv: Inverse Langevin function
 1.45.02  11.03.18  we          Rogers-Ramanujan continued fraction sfc_rrcf

 1.47.00  19.04.18  we          Fixes for FPC311
 1.47.01  23.06.18  we          sfc_detai moved to sfellint

 1.49.00  29.07.18  we          sfc_exprel

 1.50.00  21.08.18  we          sfc_epoly
 1.50.01  23.08.18  we          Move Marcum Q to sfErf
 1.50.02  26.08.18  we          sfc_besspoly

 1.51.00  01.10.18  we          sfc_expn

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
  MaxEuler = 1866; {Maximum index for Euler numbers for extended}

function sfc_agm(x,y: extended): extended;
  {-Return the arithmetic-geometric mean of |x| and |y|; |x|,|y| < sqrt(MaxExtended)}

function sfc_agm2(x,y: extended; var s: extended): extended;
  {-Return agm(x,y) and s = sum(2^i*c_i^2, i=1...), |x|,|y| < sqrt(MaxExtended)}

function sfc_ari(x: extended): extended;
  {-Return the functional inverse of R(x), R(sfc_ari(x))=x, x >= 1.125}

function sfc_besspoly(n: integer; x: extended): extended;
  {-Return yn(x), the nth Bessel polynomial}

function sfc_br(x: extended): extended;
  {-Return the Bring radical b := BR(x) with b^5 + b + x = 0}

function sfc_catf(x: extended): extended;
  {-Return the Catalan function C(x) = binomial(2x,x)/(x+1)}

function sfc_cosint(n: integer; x: extended): extended;
  {-Return cosint(n, x) = integral(cos(t)^n, t=0..x), n>=0}

function sfc_debye(n: integer; x: extended): extended;
  {-Return the Debye function D(n,x) = n/x^n*integral(t^n/(exp(t)-1),t=0..x) of order n>0, x>=0}

function sfc_einstein(n: integer; x: extended): extended;
  {-Return the Einstein function E_n, n=1..4, x > 0 for n=3,4}

function sfc_euler(n: integer): extended;
  {-Return the nth Euler number, 0 if n<0 or odd n}

function sfc_epoly(n: integer; x: extended): extended;
  {-Return the Euler polynomial E_n(x), 0 <= n < MaxBernoulli}

function sfc_eulerq(q: extended): extended;
  {-Return the EulerQ function product(1-q^n, n=1..Inf)}

function sfc_expn(n: integer; x: extended): extended;
  {-Return the truncated exponential sum function e_n = sum(x^k/k!, k=0..n), 0 <= n < MAXGAM-1}

function sfc_exprel(n: longint; x: extended): extended;
  {-Return the relative exponential = (e^x-sum(x^k/k!, k=0..n-1)*n!/x^n}

function sfc_fpoly(n: integer; x: extended): extended;
  {-Return the Fibonacci polynomial F_n(x)}

function sfc_fibf(v,x: extended): extended;
  {-Return the general Fibonacci function F_v(x)}

function sfc_kepler(M,e: extended): extended;
  {-Solve Kepler's equation, result x is the eccentric anomaly from the mean anomaly M and the }
  { eccentricity e >= 0; x - e*sin(x) = M, x + x^3/3 = M, or e*sinh(x) - x = M for e <1, =1, >1}

function sfc_LambertW(x: extended): extended;
  {-Return the Lambert W function (principal branch), x >= -1/e}

function sfc_LambertW1(x: extended): extended;
  {-Return the Lambert W function (-1 branch), -1/e <= x < 0}

function sfc_langevin(x: extended): extended;
  {-Return the Langevin function L(x) = coth(x) - 1/x, L(0) = 0}

function sfc_llinv(x: extended): extended;
  {-Return the functional inverse of the Langevin function, |x| < 1}

function sfc_lpoly(n: integer; x: extended): extended;
  {-Return the Lucas polynomial L_n(x)}

function sfc_ri(x: extended): extended;
  {-Return the Riemann prime counting function R(x), x >= 1/1024}

function sfc_rrcf(q: extended): extended;
  {-Return the Rogers-Ramanujan continued fraction for |q| < 1}

function sfc_sinint(n: integer; x: extended): extended;
  {-Return sinint(n, x) = integral(sin(t)^n, t=0..x), n>=0}

function sfc_trint(n: integer; x: extended): extended;
  {-Return the transport integral J_n(x) for x >= 0, n >= 2,}
  { J_n(x) = integral(t^n*exp(t)/(exp(t)-1)^2, t=0..x)      }

function sfc_wo(x: extended): extended;
  {-Return the Wright omega function, i.e. the solution w of w + ln(w) = x}


implementation


uses
  AMath, sfBasic, sfGamma, sfGamma2, sfZeta, sfExpInt;


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
type
  ti_rec  = record
              n: integer;
              x: extended;
            end;
  pti_rec = ^ti_rec;


{---------------------------------------------------------------------------}
function dnf(t: extended; param: pointer): extended; {$ifdef BIT16} far;{$endif}
  {-The Debye integrand, 1/x^n*integral(t^n/(exp(t)-1) }
var
  p: extended;
begin
  with pti_rec(param)^ do begin
    if t > ln_MaxExt then dnf := 0.0
    else begin
      p  := power(t/x,n);
      if t >= 44.0 then dnf := p*exp(-t)
      else dnf := p/expm1(t);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_debye(n: integer; x: extended): extended;
  {-Return the Debye function D(n,x) = n/x^n*integral(t^n/(exp(t)-1),t=0..x) of order n>0, x>=0}
var
  p,s,t,q: extended;
  ierr,j,k: integer;
  para: ti_rec;
  neval: longint;
const
  eps = 1e-17;
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
    {x>4: use double exponential integration}
    para.n := n;
    para.x := x;
    sfc_intde_p({$ifdef FPC_ProcVar}@{$endif}dnf, @para, 0, x, eps, s, q, neval, ierr);
    if ierr<>0 then sfc_debye := NaN_x
    else sfc_debye := n*s;
  end;
end;


{---------------------------------------------------------------------------}
function tjncore(t: extended; param: pointer): extended; {$ifdef BIT16} far;{$endif}
  {-integrand t^n*exp(t)/(exp(t)-1)^2 for transport integrals}
var
  z,q: extended;
begin
  with pti_rec(param)^ do begin
    if t > ln_MaxExt then tjncore := 0.0
    else if t <= 1.0 then begin
      {power cannot overflow, uses accurate function near 0}
      z := exp(t)/sqr(exprel(t));
      tjncore := z*power(t,n-2);
    end
    else begin
      q := exp(-t);
      z := q / sqr(1.0-q);
      q := ln(t);
      if n*q < ln_MaxExt then tjncore := z*power(t,n)
      else begin
        {Power would overflow, try logarithms}
        q := n*q + ln(z);
        tjncore := exp(q);
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function ti_int(n: integer; x: extended): extended;
  {-Transport integral J(n,x) for x>=4 with DE integration}
const
  eps = 1e-17;
var
  para: ti_rec;
  res, abserr: extended;
  neval: longint;
  ierr: integer;
begin
  para.n := n;
  sfc_intde_p({$ifdef FPC_ProcVar}@{$endif}tjncore, @para, 0, x, eps, res, abserr, neval, ierr);
  if ierr<>0 then ti_int := NaN_x
  else ti_int := res;
end;


{---------------------------------------------------------------------------}
function ti_sum(n: integer; x: extended): extended;
  {-Transport integral J(n,x) for 0 <= x <= 4, n > 1}
var
  s,q,t,f: extended;
  k,m: integer;
begin
  s := n-1.0;
  s := 1.0/s;
  f := 1;
  q := x*x;
  k := 0;
  {(1/(n-1)+sum(bernoulli(2*k)*(1-2*k)/(2*k+n-1)/(2*k)!*x^(2*k), k=1..Inf))*x^(n-1)}
  repeat
    k := k+2;
    m := k-1;
    f := f*q/(k*m);
    t := sfc_bernoulli(k)*m/(n+m);
    t := t*f;
    s := s - t;
  until abs(t) <= eps_x*s;
  q := ln(x);
  if n*q < ln_MaxExt then ti_sum := s*power(x,n-1)
  else begin
    t := (n-1)*q + ln(s);
    ti_sum := exp(t);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_trint(n: integer; x: extended): extended;
  {-Return the transport integral J_n(x) for x >= 0, n >= 2,}
  { J_n(x) = integral(t^n*exp(t)/(exp(t)-1)^2, t=0..x)      }
begin
  if IsNanorInf(x) or (x<0.0) or (n<2) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_trint := NaN_x;
  end
  else if x <= 4.0  then sfc_trint := ti_sum(n,x)
  else sfc_trint := ti_int(n,x);
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
  p1 = 19/10;  {Fix311}
  p2 = 17/60;
  q1 = 29/10;
  q2 = 101/60;
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
    if x > -MinExtended then begin
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
  {-Return the Riemann prime counting function R(x), x >= 1/1024}
var
  f,s,t,y: extended;
  n: integer;
begin
  if IsNanOrInf(x) or (x < 0.0009765625) then begin
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
    {For smaller x >= 1/1024 use the Gram series for R(x)}
    {R(x) = 1 + sum(ln(x)^n/(n!*n*zeta(n+1)), n=1..Inf)}
    y := ln(x);
    s := 1.0;
    f := 1.0;
    n := 1;
    repeat
      f := (f*y)/n;
      t := n*sfc_zetaint(n+1);
      t := f/t;
      s := s+t;
      inc(n);
    until abs(t) <= eps_x*abs(s);
    sfc_ri := s;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ari(x: extended): extended;
  {-Return the functional inverse of R(x), R(sfc_ari(x))=x, x >= 1.125}
var
  d,z,a,b,eps: extended;
  i: integer;
const
  IMAX = 200;
  XBIG = 1e36; {double 1e36}
  FAC  = 0.8;
begin
  {Based on Dana Jacobson's Math-Prime-Util function inverse_R}
  {in util.c from https://github.com/danaj/Math-Prime-Util    }

  if IsNanOrInf(x) or (x < 1.125) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_ari := NaN_x;
    exit;
  end;

  {Starting value for iterations is Inverse_Li(x)}
  z := sfc_ali(x);

  {if x >= XBIG we have Li(x) = Ri(x) for extended precision}
  if x < XBIG then begin
    b := MaxExtended;
    eps := 2*eps_x;
    for i:=0 to IMAX do begin
      d := sfc_ri(z) - x;
      if d=0.0 then break;
      {Use Quasi-Newton correction for R(x) ~ Li(x) - 0.5*Li(sqrt(x))}
      d := d*ln(z) / (1.0 - 0.5/sqrt(z));
      a := abs(d);
      {The first condition is essential for fast computation. If it}
      {is omitted the average loop count is about 9 times larger!  }
      if (a > b) or (a <= eps*abs(z)) then begin
        {The factor FAC is an experimental result which is used if }
        {the current correction term is greater than the previous. }
        {If the eps condition is true, it is (mostly) negligible.  }
        z := z - FAC*d;
        break;
      end;
      b := a;
      z := z - d;
    end;
  end;
  sfc_ari := z;
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
function sfc_fibf(v,x: extended): extended;
  {-Return the general Fibonacci function F_v(x)}
var
  y,z: extended;
begin
  if x=0.0 then begin
    {http://functions.wolfram.com/07.06.03.0001.01}
    y := sinpi(0.5*v);
    sfc_fibf := sqr(y);
  end
  else begin
    {http://functions.wolfram.com/07.06.02.0001.01}
    z := sqrt1pmx((-0.5)*x);
    z := power(z,v);
    y := cospi(v);
    z := z - y/z;
    y := hypot(2,x);
    sfc_fibf := z/y;
  end;
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
function EulerQ_sum(q: extended): extended;
  {-EulerQ from summation}
var
  k: integer;
  s,t1,t2,q3,qk: extended;
begin
  {Basic algorithm: use the Euler Pentagonal number theorem}
  {see https://en.wikipedia.org/wiki/Euler_function        }
  k  := 1;
  s  := 1.0;
  t1 := q;
  t2 := q*q;
  q3 := t2*q;
  qk := t1;
  while abs(t1) > eps_x*abs(s) do begin
    if odd(k) then s := s - (t1+t2)
    else s := s + (t1+t2);
    qk := qk*q3;
    t1 := t1*qk;
    t2 := t2*(qk*q);
    inc(k);
  end;
  EulerQ_sum := s;
end;


{---------------------------------------------------------------------------}
function EulerQ_ae(q: extended): extended;
  {-Asymptotic expression for EulerQ near 1}
const
  q1 = 0.999857;
var
  t,x: extended;
begin
  if q>q1 then EulerQ_ae := 0.0
  else begin
    {http://mathworld.wolfram.com/q-PochhammerSymbol.html, formula 8}
    t := -ln(q);
    x := 0.25*t - PiSqr/t;
    x := exp(x/SIXX);
    EulerQ_ae := sqrt(TwoPi/t)*x;
  end;
end;


{---------------------------------------------------------------------------}
function EulerQ_aen(q: extended): extended;
  {-Asymptotic expression for EulerQ near -1}
const
  q1 = -0.999964;
var
  t,x: extended;
begin
  if q < q1 then EulerQ_aen := 0.0
  else begin
    {Asymptotic expression for EulerQ near -1, based on euler_ae. Uses }
    {the identity EulerQ(-q) = EulerQ(q^2)^3/EulerQ(q)/EulerQ(q^4). See}
    {https://math.stackexchange.com/questions/2294697/asymptotics-of-q-pochhammer-euler-function-for-q-rightarrow-1}
    {Algebraically simplified with Maple}
    t := ln(abs(q));
    x := PiSqr/t - t;
    x := exp(0.25*x/SIXX);  {Anti FPC -O4 code}
    EulerQ_aen := sqrt(-Pi/t)*x;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_eulerq(q: extended): extended;
  {-Return the EulerQ function product(1-q^n, n=1..Inf)}
const
  qae = -0.85; {Use asymptotic expression for q <= qae}
begin
  if IsNanOrInf(q) or (q<-1) or (q>1) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_eulerq := NaN_x;
    exit;
  end;
  if q <= qae then sfc_eulerq := EulerQ_aen(q)
  else if q <= 0.5 then sfc_eulerq := EulerQ_sum(q)
  else sfc_eulerq := EulerQ_ae(q);
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


{---------------------------------------------------------------------------}
function sfc_br(x: extended): extended;
  {-Return the Bring radical b := BR(x) with b^5 + b + x = 0}
var
  f0,f1,f2,a,b,d: extended;
const
  eps =  0.3e-6;  {double 0.3e-5}
  c0  = -0.1348363788;
  c1  = -0.7996803229;
  c2  =  0.1866401933;
begin
  if IsNanOrInf(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_br := NaN_x;
    exit;
  end;
  if x=0.0 then sfc_br := 0.0
  else begin
    a := abs(x);
    if a <= 0.5 then b:= -a*(1.0 - sqr(sqr(a)))
    else if a <= 2 then b := (a*c2 + c1)*a + c0
    else begin
      {a > 2, only a good approximation of -(a-1)^(1/5) is needed here}
      b := -exp(0.2*ln(a-1.0));
    end;
    repeat
      {Halley iteration}
      d  := sqr(b);
      f2 := 20.0*b*d;
      d  := sqr(d);
      f1 := 5.0*d+1.0;
      f0 := ((d+1.0)*b+a)/f1;
      d  := f0/(1.0 - 0.5*f0*(f2/f1));
      b  := b - d;
    until abs(d) < eps*abs(b);
    if x>0.0 then sfc_br := b
    else sfc_br := -b;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_wo(x: extended): extended;
  {-Return the Wright omega function, i.e. the solution w of w + ln(w) = x}
var
  t,w,d: extended;
const
  xmax = 1.3e7;     {double = 9e5}
  xmin = -45.1;     {double = -37.5}
  x0   = 2.9e-3;    {double = 1e-2}
  eps  = 2.3841e-7; { ~ cbrt(eps_x)/2 }
const
  MLS: array[0..5] of extended = (  {Maclaurin series}
         +0.5671432904097838730,
         +0.3618962566348892215,
         +0.7367780517637275782e-1,
         -0.1342859654990087248e-2,
         -0.1636065147912497178e-2,
         +0.2321496555699633807e-3);
begin
  if IsNanOrInf(x) then begin
    if IsNan(x) then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfc_wo := NaN_x;
    end
    else if x<0.0 then sfc_wo := 0.0 {x=-inf}
    else sfc_wo := x;  {x=+inf}
    exit;
  end;
  {http://www.orcca.on.ca/TechReports/2000/TR-00-12.html}
  {compute starting values}
  if x>4.5 then begin
    t := ln(x);
    w := x - t + t/x;
    if x >= xmax then begin
      sfc_wo := w;
      exit;
    end;
  end
  else if x < -1.5 then begin
    w := exp(x);
    if x <= xmin then begin
      sfc_wo := w;
      exit;
    end
    else if w<>0.0 then begin
      t := (1.5 - (8/3)*w)*w - 1.0;   {Fix311}
      w := w + sqr(w)*t;
    end;
  end
  else begin
    w := PolEvalX(x,MLS,6);
    if abs(x) <= x0 then begin
      sfc_wo := w;
      exit;
    end;
  end;
  repeat
    {Halley iteration, observed max number = 2}
    t := (w + ln(w) - x)/(1.0+w);
    d := w*t/(1.0 + 0.5*t/(1.0+w));
    w := w - d;
  until abs(d) < eps*w;
  sfc_wo := w;
end;


{---------------------------------------------------------------------------}
function sfc_einstein(n: integer; x: extended): extended;
  {-Return the Einstein function E_n, n=1..4, x > 0 for n=3,4}
const
  XMAX = 45.0;  {DAMath = 38}
var
  err: boolean;
begin
  err := IsNanOrInf(x);
  sfc_einstein := NaN_x;
  if not err then begin
    if n=1 then begin
      {e1(x) = x^2 * exp(x) / (exp(x) - 1)^2}
      x := abs(x);
      if x>=XMAX then sfc_einstein := sqr(x)*exp(-x)
      else sfc_einstein := exp(x)/sqr(exprel(x));
    end
    else if n=2 then begin
      {e2(x) = x / (exp(x) - 1)}
      if x <= -XMAX then sfc_einstein := -x
      else sfc_einstein := 1.0/exprel(x);
    end
    else if (n=3) and (x>0.0) then begin
      {e3(x) = ln(1 - exp(-x))}
      sfc_einstein := ln1mexp(-x);
    end
    else if (n=4) and (x>0.0) then begin
      {e4(x) = x / (exp(x) - 1) - ln(1 - exp(-x))}
      if x>ln_MaxExt then sfc_einstein := 0
      else if x>=XMAX then sfc_einstein := (x+1.0)*exp(-x)
      else sfc_einstein := 1.0/exprel(x) - ln1mexp(-x);
    end
    else err := true;
  end;
  if err then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
  end;
end;


{---------------------------------------------------------------------------}
const
  nlv = 13;
  lvxh: array[0..nlv-1] of THexExtW = ( {chebyshev((coth(x)-1/x)/x, x=-1..1, 0.5e-20);}
          ($1668,$616C,$5EB9,$A55A,$3FFE),  {+0.645910187014507952778305035000    }
          ($6B57,$70C2,$3BA4,$A631,$BFF8),  {-0.101435739941388114059644663986e-1 }
          ($1CD9,$3F5E,$8295,$F036,$3FF2),  {+0.229084901846336481169924923413e-3 }
          ($6D31,$FD6D,$5BBB,$B6C3,$BFED),  {-0.544676537826296496965829810972e-5 }
          ($77B5,$5BE1,$0FC6,$8C96,$3FE8),  {+0.130931081441182878982618509970e-6 }
          ($0A40,$94AD,$BAF9,$D8DA,$BFE2),  {-0.315564707141221862295501357833e-8 }
          ($E95E,$7D1E,$032A,$A75C,$3FDD),  {+0.761062543956059982739459695821e-10}
          ($AEBF,$A763,$D606,$812E,$BFD8),  {-0.183580018070672141940213332479e-11}
          ($69AB,$B060,$4401,$C770,$3FD2),  {+0.442842513100509316169830548635e-13}
          ($4F3F,$8E48,$F02D,$99F3,$BFCD),  {-0.106826272532550650314814814815e-14}
          ($53E4,$15EC,$CEB3,$EDAE,$3FC7),  {+0.257696253100437571692598505481e-16}
          ($B73B,$4517,$BCB8,$B779,$BFC2),  {-0.621639295641256569625648606622e-18}
          ($3E53,$D2C7,$8ED0,$8DA1,$3FBD)); {+0.149957744659012073828688846172e-19}
var
  lvcs: array[0..nlv-1] of extended absolute lvxh;


{---------------------------------------------------------------------------}
function sfc_langevin(x: extended): extended;
  {-Return the Langevin function L(x) = coth(x) - 1/x, L(0) = 0}
var
  t,y: extended;
begin
  t := abs(x);
  if t <= 1.0 then sfc_langevin := x*CSEvalX(2.0*sqr(x)-1.0, lvcs, nlv)
  else begin
    {L(x) = coth(x)-1/x = (1+exp(-2x))/(1-exp(-2x) - 1/x}
    {L(x) = (1 + exp(-2x))/(1 - exp(-2x)) - 1 + 1 - 1/x }
    {L(x) = 2*exp(-2x)/(1-exp(-2x)) + (x-1)/x           }
    if t>23.0 then y := 1.0 - 1.0/t
    else begin
      y := exp(-2.0*t);
      y := 2.0*y/(1.0-y) + (t-1.0)/t;
    end;
    if x<0.0 then sfc_langevin := -y
    else sfc_langevin := y;
  end;
end;


{---------------------------------------------------------------------------}
procedure LLDer(x: extended; var y0,y1: extended);
  {-Compute Langevin function y0=L(x), y1=L'(x)}
var
  t,y,z: extended;
begin
  t := abs(x);
  if t <= 1.0 then begin
    CSEvalXDer(x, lvcs, nlv, y, y1);
    y1 := sqr(2.0*x)*y1 + y;
    y0 := x*y;
  end
  else begin
    {See sf_ll}
    if t>27.0 then begin
      y0 := 1.0 - 1.0/t;
      y1 := 1.0/sqr(t);
    end
    else begin
      y  := exp(-2.0*t);
      z  := 2.0/(1.0-y);
      y0 := y*z + (t-1.0)/t;
      y1 := 1.0/sqr(t) - y*sqr(z);
    end;
    if x<0.0 then y0 := -y0;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_llinv(x: extended): extended;
  {-Return the functional inverse of the Langevin function, |x| < 1}
var
  z,t,y0,y1,d: extended;
const
  t0 = 2e-5;
  t1 = 31/32;  {Fix311}
  eps = 1e-17;
begin
  if IsNanOrInf(x) or (abs(x) >= 1.0) then begin
    sfc_llinv := x;
    exit;
  end;
  t := abs(x);
  if t <= t0 then z := t*(3.0 + 1.8*sqr(t))
  else if t >= t1 then z := 1.0/(1-t)
  else begin
    {Newton iteration with Redynak approximation as start value}
    {https://en.wikipedia.org/Brillouin_and_Langevin_functions }
    z := (1.0 - t) * (1.0 + 0.1*t);
    z := (3.0 - 2.6*t + 0.7*t*t)*t/z;
    repeat
      LLDer(z,y0,y1);
      d := (y0-t)/y1;
      z := z - d;
    until abs(d) <= eps*z;
  end;
  if x<0.0 then sfc_llinv := -z
  else sfc_llinv := z;
end;


{---------------------------------------------------------------------------}
function rr_cf(q: extended): extended;
  {-Return the Rogers-Ramanujan continued fraction for |q| < 1}
var
  a,b,c,d,f,t: extended;
  k: integer;
const
  MAXIT = 50;
begin
  {Modified Lentz method with pre-calculated values for k=1}
  a := 1.0;
  b := 1.0;
  d := 1.0;
  f := 1.0;
  c := Sqrt_MaxExt;
  for k:=2 to MAXIT do begin
    a := a*q;
    d := b + a*d;
    c := b + a/c;
    if c=0.0 then c := Sqrt_MinExt;
    if d=0.0 then d := Sqrt_MinExt;
    d := 1.0/d;
    t := c*d;
    f := f*t;
    if abs(t-1.0) <= eps_x then break;
  end;
  rr_cf := f*nroot(q,5);;
end;


{---------------------------------------------------------------------------}
function sfc_rrcf(q: extended): extended;
  {-Return the Rogers-Ramanujan continued fraction for |q| < 1}
const
  g = 1.5 + 0.1180339887498948482;     { (sqrt(5)+1)/2 }
  h = 0.5 + 0.1180339887498948482;     { (sqrt(5)-1)/2 }
var
  a,f: extended;
begin
  if IsNanOrInf(q) or (abs(q) >= 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_rrcf := NaN_x;
    exit;
  end;
  if q >= 0.84375 then begin
    {limit q -> 1}
    sfc_rrcf := h;
  end
  else if q > 1/128 then begin
    {Chan [74] (1.4)}
    a := (4.0*PiSqr)/ln(q);
    a := exp(a);
    f := rr_cf(a);
    sfc_rrcf := (1.0 - g*f)/(g+f);
  end
  else if q > -1/16 then begin
    sfc_rrcf := rr_cf(q);
  end
  else if q > -0.9921875 then begin
    {Chan [74] (1.5)}
    a := PiSqr/ln(abs(q));
    a := -exp(a);
    f := -rr_cf(a);
    sfc_rrcf := (h*f-1.0)/(h+f);
  end
  else begin
    {limit q -> -1}
    sfc_rrcf := -g;
  end;
end;


{---------------------------------------------------------------------------}
function expreln_cf(n: longint; x: extended): extended;
  {-Continued fraction for expreln}
const
  bigh : THexExtw = ($0000,$0000,$0000,$8000,$5FA3);  {2^8100 = 2.2027483177801303050E+2438}
  kmax = 5000;
var
  big: extended absolute bigh;
var
  a2,b2,a1,b1,ak,bk,fk,a,b,tol: extended;
  k: longint;
begin
  {Ref HMF[1], 4.2.41 and GSL[21], function exp.c}
  k  := 2;
  ak := n+1;
  bk := ak - x;
  fk := ak/bk;
  a1 := 1.0;
  b1 := 1.0;
  tol:= 2.0*eps_x;

  while k < kmax do begin
    inc(k);
    a2 := a1;
    b2 := b1;
    a1 := ak;
    b1 := bk;
    if odd(k) then a := ((k-1) shr 1)*x
    else a := (1 - n - (k shr 1))*x;
    b := n + k - 1;
    ak := b*a1 + a*a2;
    bk := b*b1 + a*b2;

    if (abs(ak) > big) or (abs(bk) > big) then begin
      ak := ak / big;
      bk := bk / big;
      a1 := a1 / big;
      b1 := b1 / big;
      {No need to scale a2,b2; they will be set in the next loop round}
    end;

    a  := fk;
    fk := ak/bk;
    b  := a/fk;
    if (abs(b - 1.0) < tol) then begin
      expreln_cf := fk;
      exit;
    end;
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  expreln_cf := Nan_x;
end;


{---------------------------------------------------------------------------}
function sfc_exprel(n: longint; x: extended): extended;
  {-Return the relative exponential = (e^x-sum(x^k/k!, k=0..n-1)*n!/x^n}
var
  s,d,t: extended;
  k: longint;
begin
  {expreln = 1F1(1,1+n,x) = e^x * x^(-n) * [Gamma(1 + n) - n*Gamma(n, x)]}
  if (n<0) or IsNaNorInf(x) then begin
    sfc_exprel := Nan_x;
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
  end
  else if n=0 then sfc_exprel := exp(x)
  else if n=1 then sfc_exprel := exprel(x)
  else begin
    if x > 0.25*n then begin
      {Incomplete Gamma representation}
      sfc_incgamma_ex(n,x,s,t,d,true);
      sfc_exprel := s/d*n;
    end
    else if abs(x) < n*sqrt_epsh then begin
      {First to terms of 1F1(1,1+n,x)}
      sfc_exprel := 1.0 + x/(n+1);
    end
    else if x >= -10*n then begin
      {Continued fraction HMF[1], 4.2.41}
      sfc_exprel := expreln_cf(n, x);
    end
    else begin
      {The asymptotic expression for x -> -infinity is derived from the first}
      {part in HMF[1], 13.5.1 with a=n, b=n+1, the 2nd part is exponentially }
      {small: M(1,1+n,x) ~ -n/x*[1 + (n-1)/x + (n-1)(n-2)/x^2 + ...] }
      s := 1.0;
      t := 1.0;
      k := n-1;
      repeat
        t := t*k/x;
        s := s + t;
        dec(k);
      until abs(t) <= abs(s)*eps_x;
      sfc_exprel := -s*n/x;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_epoly(n: integer; x: extended): extended;
  {-Return the Euler polynomial E_n(x), 0 <= n < MaxBernoulli}
var
  s,t,w,y,z,e: extended;
  k,m: integer;
begin
  sfc_epoly := Nan_x;
  if IsNan(x) or (n<0) or (n >= MaxBernoulli) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;
  m := n+1;
  if n<3 then begin
    {NIST[30], Table 24.2.2}
    if n=0 then sfc_epoly := 1.0
    else if n=1 then sfc_epoly := x-0.5
    else sfc_epoly := x*(x-1.0)
  end
  else if x=0.0 then begin
    {NIST[30], 24.4.26}
    y := ldexp(1,m)-1.0;
    z := sfc_bernoulli(m)/m;
    sfc_epoly := -2.0*y*z;
  end
  else if x=0.5 then begin
    {NIST[30], 24.4.28}
    z := sfc_euler(n);
    sfc_epoly := ldexp(z, -n);
  end
  else if x=1.0 then begin
    {NIST[30], 24.4.26}
    y := ldexp(1,m)-1.0;
    z := sfc_bernoulli(m)/m;
    sfc_epoly := 2.0*y*z;
  end
  else begin
    if (n >= 100) and (abs(x) >= 2.0) then begin
      {NIST[30], 24.4.22}
      y := sfc_bernpoly(m,x);
      z := sfc_bernpoly(m,0.5*x);
      z := ldexp(z,m);
      sfc_epoly := ((y - z)/m)*2.0;
    end
    else begin
      {NIST[30], 24.4.14}
      z := 1.0;
      s := 0.0;
      e := 0.0;
      k := 0;
      w := ldexp(1,n+1);
      repeat
        m := n-k+1;
        if (m=1) or (m and 1 = 0) then begin
          t := (1.0-w)*sfc_bernoulli(m)*z;
          {Kahan sum s := s + t}
          t := t - e;
          y := s + t;
          e := (y - s) - t;
          s := y;
        end;
        inc(k);
        if k > n then break;
        z := z*x*(m/k);
        w := 0.5*w;
      until false;
      sfc_epoly := 2.0*s/(n+1);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_besspoly(n: integer; x: extended): extended;
  {-Return yn(x), the nth Bessel polynomial}
var
  y0,y1,yk: extended;
  k: integer;
begin
  {Ref: [77] H.L. Krall, O. Frink}
  if (x=0.0) or (n=0) then begin
    sfc_besspoly := 1.0;
    exit;
  end;

  {Map negative orders}
  if n<0 then n := -(n+1);

  {Recurrence relation :  y(0,x) = 1, y(1,x) = x+1}
  y0 := 1.0;
  y1 := x + 1.0;
  for k:=2 to n do begin
    yk := (2*k-1)*x*y1 + y0;
    y0 := y1;
    y1 := yk;
  end;
  sfc_besspoly := y1;
end;


{---------------------------------------------------------------------------}
function sfc_expn(n: integer; x: extended): extended;
  {-Return the truncated exponential sum function e_n = sum(x^k/k!, k=0..n), 0 <= n < MAXGAM-1}
var
  e,f,g: extended;
begin
  if IsNanOrInf(x) or (n<0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_expn := NaN_x;
  end
  else if x=0.0 then sfc_expn := 1.0
  else if n<=2 then begin
    case n of
        0: sfc_expn := 1.0;
        1: sfc_expn := 1.0 + x;
      else sfc_expn := 1.0 + x*(1.0 + 0.5*x);
    end;
  end
  else begin
    {http://mathworld.wolfram.com/ExponentialSumFunction.html}
    g := sfc_igamma(n+1,x);
    f := sfc_fac(n);
    e := exp(x);
    if abs(e) > abs(g) then
      sfc_expn := (e/f)*g
    else
      sfc_expn := (g/f)*e;
  end;
end;


end.
