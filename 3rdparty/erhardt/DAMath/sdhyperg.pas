unit sdHyperG;

{Double precision special functions: Hypergeometric functions and related}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Double precision special functions: Hypergeometric functions

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARKS       :  - sfd_chu is still beta
                  - In sfd_whitm the sfd_1f1 call may produce overflows

 REFERENCES    :  References used in this unit, main index in damath_info.txt/references

                  [1] [HMF]: M. Abramowitz, I.A. Stegun. Handbook of Mathematical Functions, New York, 1970
                      http://www.math.sfu.ca/~cbm/aands/
                  [7] Cephes Mathematical Library, Version 2.8
                      http://www.moshier.net/#Cephes or http://www.netlib.org/cephes/
                 [14] SLATEC Common Mathematical Library, Version 4.1, July 1993
                      (general purpose mathematical and statistical routines written in Fortran 77)
                      http://www.netlib.org/slatec
                 [20] Special functions by Wayne Fullerton,
                      http://www.netlib.org/fn
                      Almost identical to the FNLIB subset of SLATEC [14]
                 [30] [NIST]: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST Handbook
                      of Mathematical Functions, Cambridge, 2010. Online resource: NIST Digital
                      Library of Mathematical Functions, http://dlmf.nist.gov/
                 [33] http://functions.wolfram.com/: Formulas and graphics about
                      mathematical functions for the mathematical and scientific
                      community and/or http://mathworld.wolfram.com/ ("/the web's
                      most extensive mathematical resource/")
                 [42] Y.L. Luke, Algorithms for the Computation of Mathematical Functions,
                      Academic Press, 1977
                 [43] J. Pearson, Computation of Hypergeometric Functions,
                      Master's thesis, University of Oxford, 2009. Available as
                      http://people.maths.ox.ac.uk/porterm/research/pearson_final.pdf
                 [44] K.E. Muller, Computing the confluent hypergeometric function
                      M(a,b,x), Numerische Mathematik 90, 179-196, 2001
                 [45] N.N. Lebedev, Special Functions and Their Applications,
                      Dover, New York, 1972
                 [53] F.G. Tricomi, Fonctions hypergeometriques confluentes. Memorial
                      des sciences mathematiques, 140 (1960), p. 1-86. Available from
                      http://www.numdam.org/item?id=MSM_1960__140__1_0
                 [54] R.C. Forrey, Computing the hypergeometric function, J. Comp. Phys. 137, 79-100 (1997).
                      Available as http://physics.bk.psu.edu/pub/papers/hyper.pdf,
                      Fortran code from http://physics.bk.psu.edu/codes.html.
                 [55] N.M. Temme, The Numerical Computation of the Confluent Hypergeometric
                      Function U(a,b,z), Numerische Mathematik, 41, p.63-82, 1983. Available
                      from http://www.digizeitschriften.de/en/dms/toc/?PPN=PPN362160546_0041
                      or http://oai.cwi.nl/oai/asset/10717/10717D.pdf


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.03.00  25.05.13  W.Ehrhardt  Initial BP7 version from AMath.sfhyperg
 1.03.01  25.05.13  we          Fix overflows in h1f1_acf
 1.03.02  31.05.13  we          AMath changes for sfd_chu

 1.04.00  15.06.13  we          CHU: Temme's chu and uabx from AMath
 1.04.01  16.06.13  we          Whittaker functions
 1.04.02  28.06.13  we          sfd_0f1
 1.04.03  01.07.13  we          sfd_0f1r

 1.05.00  28.07.13  we          MAXIT=64 in h1f1_tricomi

 1.09.00  24.03.14  we          CHU: chu_luke only with $ifdef in dchu
 1.09.01  28.03.14  we          InfOrNan tests in sfd_1f1

 1.12.00  13.07.14  we          sfc_pcfd,sfc_pcfu,sfc_pcfhh

 1.13.00  11.08.14  we          fix: sfc_pcfd,sfc_pcfu,sfc_pcfhh use double
 1.13.01  19.08.14  we          sfc_pcfv

 1.14.00  01.09.14  we          sfd_chu for x<0
 1.14.01  04.09.14  we          fix h1f1_tricomi (division by t=0)

 1.15.00  26.10.14  we          No Kummer transformation in 1F1 if b=0,-1,-2,...
 1.15.01  29.10.14  we          1F1 = 1 + a/b*(exp(x)-1) for very small a,b

 1.16.00  05.01.15  we          Special case b=1 in gen_0f1

 1.18.00  26.06.15  we          Fix: double variables in chu_negx

 1.19.00  01.08.15  we          Avoid some overflows in h2f1_sum, improved h2f1_trans

 1.20.00  25.08.15  we          Use THREE & SIXX to avoid 'optimization'

 1.21.00  29.09.15  we          sfd_pcfv(a,x) for a < 0

 1.34.00  03.08.18  we          Fix tmax update in hyp2f0
 1.34.01  03.08.18  we          sfd_2f0

 1.36.00  05.10.18  we          Changed threshold in h1f1_acf

***************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2010-2016 Wolfgang Ehrhardt

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


function sfd_1f1(a,b,x: double): double;
  {-Return the confluent hypergeometric function 1F1(a,b,x); Kummer's function M(a,b,x)}

function sfd_1f1r(a,b,x: double): double;
  {-Return the regularized Kummer hypergeometric function 1F1(a,b,x)/Gamma(b)}

function sfd_2f1(a,b,c,x: double): double;
  {-Return the Gauss hypergeometric function 2F1(a,b;c;x)}

function sfd_2f1r(a,b,c,x: double): double;
  {-Return the regularized Gauss hypergeometric function 2F1(a,b,c,x)/Gamma(c)}

function sfd_chu(a,b,x: double): double;
  {-Return Tricomi's confluent hypergeometric function U(a,b,x). If}
  { x<0, then a must be an integer and a<0 or 1+a-b an integer < 0.}

function sfd_whitm(k,m,x: double): double;
  {-Return the Whittaker M function = exp(-x/2)*x^(0.5+m) * 1F1(m-k-0.5,2m+1,x)}

function sfd_whitw(k,m,x: double): double;
  {-Return the Whittaker W function = exp(-x/2)*x^(0.5+m) * U(m-k-0.5,2m+1,x)}

function sfd_0f1(b,x: double): double;
  {-Return the confluent hypergeometric limit function 0F1(;b;x)}

function sfd_0f1r(b,x: double): double;
  {-Return the regularized confluent hypergeometric limit function 0F1(;b;x)/Gamma(b)}

function sfd_2f0(a,b,x: double): double;
  {-Return 2F0(a,b,x), if x>0 then a or b must be a negative integer}

function sfd_pcfd(v,x: double): double;
  {-Return Whittaker's parabolic cylinder function D_v(x)}

function sfd_pcfu(a,x: double): double;
  {-Return the parabolic cylinder function U(a,x)}

function sfd_pcfv(a,x: double): double;
  {-Return the parabolic cylinder function V(a,x) with 2a integer}

function sfd_pcfhh(v,x: double): double;
  {-Return the Hermite function H_v(x) of degree v}


implementation


uses
  DAMath,
  sdBasic, sdGamma, sdPoly, sdBessel;


const
  ETHRESH = 1.1e-14; {~  50 eps_d}
  IEPS    = 3.3e-16; {~ 1.5 eps_d}


{---------------------------------------------------------------------------}
function is_negint(x: double): boolean;
  {-Return true if x is a negative integer or 0}
begin
  is_negint := (x<=0.0) and (x > -MaxLongint) and (abs(x-rint(x)) <= IEPS);
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function h1f1_b2a(a,x: double; var OK: boolean): double;
  {-Return M(a,2a,x) using bessel_iv, 2a no negative integer}
var
  s,t,v,y: double;
  gsig: integer;
begin
  OK := true;
  if (x=0.0) or (a=0.0) then h1f1_b2a := 1.0
  else begin
    {HMF[1], 13.6.3, Lebedev[45] 9.13.14. This is valid for x>=0. For x<0 use}
    {Kummer: M(a,2a,x) = exp(x)*M(a,2a,-x) and exp(x)*exp(-x/2) = exp(x/2}
    {With v = a-1/2, y = |x/2| we get:}
    {x>0: M = (y/2)^(-v) * Gamma(a+0.5) * exp(y)  * I_v(y) }
    {       = (y/2)^(-v) * Gamma(a+0.5) * exp(x)  * I_ve(y)}
    {x<0: M = (y/2)^(-v) * Gamma(a+0.5) * exp(-y) * I_v(y) }
    {       = (y/2)^(-v) * Gamma(a+0.5) *           I_ve(y)}
    if (x>=ln_MaxDbl) or (abs(a) > MAXGAMD-5) then begin
      {Calculate all factors except I_ve using logarithmic terms}
      y := abs(0.5*x);
      v := a-0.5;
      s := sfd_ive(v,y);
      if s=0.0 then begin
        {scaled i_v underflows use unscaled as last resort}
        t := sfd_lngammas(a+0.5,gsig);
        t := t + 0.5*x;
        t := t - ln(0.5*y)*v;
        s := sfd_iv(v,y);
      end
      else begin
        t := sfd_lngammas(a+0.5,gsig);
        if x>0 then t := t+x;
        t := t - ln(0.5*y)*v;
      end;
      if t>=ln_MaxDbl then begin
        OK := false;
        h1f1_b2a := PosInf_d;
      end
      else h1f1_b2a := exp(t)*s*gsig;
    end
    else begin
      v := a - 0.5;
      y := abs(0.5*x);
      s := power(0.5*y,-v);
      t := sfd_gamma(a+0.5);
      s := s*t;
      if x>0.0 then s := s*exp(x);
      t := sfd_ive(v,y);
      h1f1_b2a := s*t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function h1f1_tricomi(a,b,x: double; var OK: boolean): double;
  {-Compute M(a,b,x) with Tricomi version of Bessel expansion, b no negative integer}
var
  j: integer;
  d0,d1,d2,f,q,r,s,t,w,y: double;
const
  MAXIT = 64;
  EPST  = 0.243e-62; {~ eps_d^4}
begin

  w  := b - 2.0*a;
  if (w=0.0) or (x=0.0) then begin
    h1f1_tricomi := h1f1_b2a(a,x,OK);
    exit;
  end;

  h1f1_tricomi := 0;
  OK := false;
  y  := 2.0*x*w;
  if y<=0.0 then exit;

  {Ref: Tricomi[53], 1.6(9); HMF[1], 13.3.7, p.506; Pearson[43], 3.20/21}
  d0 := 1.0;
  d1 := 0.;
  d2 := 0.5*b;
  f  := sfd_lngammas(b,j);
  t  := f + 0.5*x + b*ln2;
  if t >= ln_MaxDbl then exit;
  f := 0.5*j*exp(t);

  y := sqrt(y);
  t := power(y,b-1);
  {Exit 'not OK' for small t}
  if abs(t) < EPST then exit;

  r := x/y;
  q := sqr(r)/t;
  s := sfd_jv(b-1,y)/t + d2*q*sfd_jv(b+1,y);
  for j:=3 to MAXIT do begin
    t  := (b+(j-2))*d1 - w*d0;
    d0 := d1;
    d1 := d2;
    d2 := t/j;
    q  := q*r;
    t  := d2*q*sfd_jv(b+(j-1),y);
    s  := s + t;
    if abs(t)<abs(s)*eps_d then begin
      OK := true;
      h1f1_tricomi := s*f;
      exit;
    end;
  end;
  {No convergence}
  h1f1_tricomi := s*f;
end;


{---------------------------------------------------------------------------}
function h1f1_luke(a,b,x: double): double;
  {-Compute M(a,b,x) with Luke rational approximation, x < 0}
const
  NMAX = 10000;
const
  sc_large = 1152921504606846976.0;  {2^60}
  sc_small = 1.0/sc_large;
var
  z,z3,t,h,e,e0: double;
  an,a2,a1,a0,bn,b2,b1,b0: double;
  na,nb,f1,f2,f3: double;
  n: integer;
begin
  {Luke [42], Ch. 15: Rational Approximations for 1F1(a,c;-z)}
  z  := -x;
  z3 := z*z*z;
  e0 := 1.0;

  e  := a/b;
  t  := 0.5*(a+1.0)/b;
  h  := 0.5*(a+2.0)/(b+1.0);

  {Starting values, [42] 15(6)}
  b0 := 1.0;
  b1 := 1.0 + t*z;
  b2 := 1.0 + h*z*(1.0 + t*z/THREE);
  a0 := 1.0;
  a1 := b1 - e*z;
  a2 := b2 - e*(1.0 + h*z)*z + e*t*(b/(b+1.0))*z*z;

  for n:=3 to NMAX do begin
    {A_n and B_n via recurrence formulas [42] 15(7)}
    na := n-1+a;
    nb := n-1+b;

    h  := nb - 1.0;
    t  := 2*n - 3;
    e  := -na*(n-1-b) / (2.0*t*h*nb);

    f1 := 1.0 + (n-2-a)/(2.0*t*nb)*z;
    f2 := z*(e + 0.25*(n+a)*na/((2*n-1)*t*h*nb)*z);
    f3 := 0.125*z3*(1.0-na)*na*(n-2-a)/(t*t*(t-2.0)*(n-3+b)*h*nb);

    an := f1*a2 + f2*a1 + f3*a0;
    bn := f1*b2 + f2*b1 + f3*b0;
    e  := an/bn;

    if abs(e-e0)<eps_d*abs(e) then begin
      h1f1_luke := e;
      exit;
    end;

    if (abs(an) > sc_large) or (abs(bn) > sc_large) then begin
      an := an * sc_small;
      bn := bn * sc_small;
      a2 := a2 * sc_small;
      b2 := b2 * sc_small;
      a1 := a1 * sc_small;
      b1 := b1 * sc_small;
    end
    else if (abs(an) < sc_small) or (abs(bn) < sc_small) then begin
      an := an * sc_large;
      bn := bn * sc_large;
      a2 := a2 * sc_large;
      b2 := b2 * sc_large;
      a1 := a1 * sc_large;
      b1 := b1 * sc_large;
    end;

    e0 := e;
    b0 := b1;
    b1 := b2;
    b2 := bn;
    a0 := a1;
    a1 := a2;
    a2 := an;
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  h1f1_luke := e0;
end;


{---------------------------------------------------------------------------}
function h1f1_acf(a,b,x: double): double;
  {-Compute M(a,b,x) with continued fraction for asymptotic series}
var
  rj, a0, a1, a2,f: double;
  j,k: integer;
const
  MAXIT = 1000;
begin
  {Continued fraction for asymptotic series in 1/x, Muller [44], Method 2.C}
  if x<0 then begin
    {Perform Kummer Transformation M(a,b,x) = exp(x)*M(b-a,b,-x)}
    x := -x;
    a := b-a;
    f := sfd_gamma_ratio(b,a)*power(x,a-b);
  end
  else begin
    {AMath uses x>=Ln_MaxExt, but DAMath with x>=Ln_MaxDbl overflows on some}
    {test cases, so use a safer path and get a little increased inaccuracy.}
    if x >= 0.5*Ln_MaxDbl then begin
      a0 := sfd_lngammas(b,j);
      a1 := sfd_lngammas(a,k);
      a2 := (a-b)*ln(x);
      f  := exp(x + a0 - a1 + a2);
      f  := (j*k)*f;
    end
    else begin
      f := exp(x)*sfd_gamma_ratio(b,a)*power(x,a-b);
    end;
  end;
  a0 := 1.0;
  a1 := 1.0 + (b-a)*(1.0-a)/x;
  for j:=2 to MAXIT do begin
    rj:= (b-a+j-1)*(j-a)/j;
    a2 := a1 + (a1-a0)*rj/x;
    if (abs(a1-a0)<eps_d*abs(a0)) and (abs(a2-a1)<eps_d*abs(a1)) then begin
      h1f1_acf := f*a2;
      exit;
    end;
    a0 := a1;
    a1 := a2;
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  h1f1_acf := f*a2;
end;


{---------------------------------------------------------------------------}
function h1f1_sum(a,b,x: double; var loss: double): double;
  {-Power series or polynomial summation of confluent hypergeometric}
  { function; compensated sum with TwoSum, loss = estimated relative error}
var
  m,s,u,umax: double;
  e,w,y: double;
  i: integer;
  done: boolean;
const
  IMAX = 32000;
  lfac = 0.03125;
  ifac = 0.01;
begin
  if abs(b)<IEPS then begin
    loss := 1.0;
    h1f1_sum := PosInf_d;
    exit;
  end;
  i := 0;
  s := 1.0;
  u := 1.0;
  e := 0.0;
  umax := 0.0;
  done := false;
  repeat
    {Next term}
    w := a+i;
    y := (b+i)*(i+1);
    if y=0.0 then begin
      if w=0.0 then u := 0.0
      else begin
        {w/0 is to be computed! Signal no convergence:}
        i := IMAX;
      end;
    end
    else u := x*(w/y)*u;
    if u<>0.0 then begin
      {remember largest term summed}
      m := abs(u);
      if m>umax then umax := m;
      {(w,y) = TwoSum(u,s)}
      w := u+s;
      m := w-u;
      y := (u-(w-m)) + (s-m);
      {sum of errors}
      e := e+y;
      s := w;
      if abs(u)<eps_d*abs(s) then begin
        {Small term, but check if b+i > 0  (if |b| is not too large)}
        {I.e. keep on summing if b+i may cross zero for larger index}
        {This is especially important if an accidental small term is}
        {encountered in the polynomial case when u=0 for a larger i.}
        done := (b < -IMAX) or (b+i > 0.0);
      end;
      inc(i);
    end
    else done := true;
  until done or (i>IMAX);
  {loss = estimated relative error}
  s := s+e;
  h1f1_sum := s;
  if s=0.0 then s:=1.0;
  if i>IMAX then loss := 1.0
  else begin
    loss := maxd(ifac*i,lfac*umax/abs(s));
    loss := eps_d*maxd(1.0,loss);
  end;
end;


{---------------------------------------------------------------------------}
function hyp2f0(a,b,x: double; var OK: boolean): double;
  {-Asymptotic sum, OK if convergence, x should be small: 0 < x < 0.02}
var
  an,bn,a0,al,t,tlst,tmax,u,sum: double;
  n: integer;
const
  nmax = 100;     {Error is O(x^n) so this should be enough}
begin
  hyp2f0 := 0;
  {Ref. HMF[1], 13.5 and Cephes[7], function hyp2f0l in hypergl.c.}
  {My version uses only the x -> Inf part and omits the converging}
  {factors that "do not actually seem to accomplish very much".   }
  OK := false;
  an := a;
  bn := b;
  a0 := 1.0;
  al := 1.0;
  sum  := 0.0;
  tlst := 0.5*MaxDouble;
  tmax := 0.0;
  for n:=1 to nmax do begin
    u := an*((bn/n)*x);
    {check for blowup}
    t := abs(u);
    if (t>1.0) and (tmax > MaxDouble/t) then exit;
    a0 := a0*u;
    t  := abs(a0);
    {terminating condition for asymptotic series}
    if t>tlst then exit else tlst := t;
    {the sum is one term behind}
    sum := sum + al;
    al  := a0;
    if t<=eps_d then begin
      OK := true;
      hyp2f0 := sum + a0;
      exit;
    end;
    an := an+1.0;
    bn := bn+1.0;
    if t>tmax then tmax := t; {fix added in V1.34.00}
  end;
end;


{---------------------------------------------------------------------------}
function h1f1_asym(a,b,x: double; var OK: boolean): double;
  {- Returnm M(a,b,x) = Gamma(b)/Gamma(a)*e^x*x^(a-b)*2F0(b-a,1-a,1/x) for x -> INF}
var
  sa,sb: integer;
  f,t,la,lb: double;
begin
  {Ref. HMF[1], 13.5 and Cephes[7], function hy1f1al in hypergl.c.}
  {My version uses only the x -> Inf part and omits the converging}
  {factors that "do not actually seem to accomplish very much".   }
  h1f1_asym := 0.0;
  OK  := false;
  if is_negint(a) or is_negint(b) then exit;
  f := hyp2f0(b-a,1.0-a,1.0/x,OK);
  if OK then begin
    lb := sfd_lngammas(b,sb);
    la := sfd_lngammas(a,sa);
    t  := (a-b)*ln(x);
    t  := lb - la + x + t;
    h1f1_asym := exp(t)*f*(sa*sb);
  end;
end;


const
  THRESH_LAGUERRE_A   = -1000;
  THRESH_LAGUERRE_Err = 1e-12;

{---------------------------------------------------------------------------}
function h1f1_sum_laguerre(a,b,x: double): double;
  {-Polynomial case, a negative integer, use compensated sum or Laguerre polynomial}
var
  r: double;
  n: integer;
begin
  if ((x<0.0) and (b>0.0)) or (b<=0) then begin
    {Laguerre: parameter b-1 <= -1 or simplified 1F1 code, so try}
    {compensated summation and switch to Laguerre if error is large}
    h1f1_sum_laguerre := h1f1_sum(a,b,x,r);
    if r<THRESH_LAGUERRE_Err then exit;
    {Very inaccurate, so try Laguerre}
  end;
  {Use generalized Laguerre polynomials (HMF [1], 13.6.9)}
  {if the error estimate for compensated sum is to large }
  n := -round(a);
  {a = -n; Laguerre M(a,b,x) = n*beta(b,n)*Ln(b-1,x)}
  r := sfd_laguerre(n,b-1.0,x);
  h1f1_sum_laguerre := n*sfd_beta(b,n)*r;
  exit;
end;


{---------------------------------------------------------------------------}
function sfd_1f1(a,b,x: double): double;
  {-Return the confluent hypergeometric function 1F1(a,b,x); Kummer's function M(a,b,x)}
var
  r,err: double;
  a_ni,OK: boolean;
begin
  if IsNanOrInfD(a) or IsNanOrInfD(b) or IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_1F1 := NaN_d;
    exit;
  end;

  if (x=0.0) or (a=0.0) then begin
    sfd_1f1 := 1.0;
    exit;
  end
  else if b=a then begin
    sfd_1f1 := exp(x);
    exit;
  end
  else if (abs(b)<=IEPS) and (b<>0.0) and (abs(a)<=IEPS) then begin
    {a and b ~ 0, b <> 0: 1F1 = 1 + a/b*(exp(x)-1), a/b cannot overflow:  }
    {M(a,b,x) = 1 + a/b*x + a/b*(a+1)/(b+1)x^2/2 + (a)_3/(b_3)x^3/3! + ...}
    {   = 1 + a/b[x + (a+1)/(b+1)x^2/2 + (a+1)/(b+1)(a+2/(b+2)x^3/3! + ...}
    {   ~ 1 + a/b[x + (1)/(1)x^2/2 + (1)/(1)(2)/(2)x^3/3! + ...           }
    {   = 1 + a/b[-1 + 1 + x + x^2/2 + x^3/3! + ... = 1 + a/b[exp(x) - 1] }

    sfd_1f1 := 1.0 + (a/b)*expm1(x);
    exit;
  end;

  a_ni := is_negint(a);

  if is_negint(b) then begin
    {b non-positive integer, defined if a negative integer b <= a <= 0}
    if (b=0.0) or (not a_ni) or (a<b) then begin
      sfd_1f1 := NaN_d;
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
    end
    else begin
      {here b < a < 0}
      if a>THRESH_LAGUERRE_A then begin
        {Compensated sum or Laguerre}
        sfd_1f1 := h1f1_sum_laguerre(a,b,x);
      end
      else begin
        sfd_1f1 := h1f1_sum(a,b,x,err);
        {$ifdef debug}
          {Kummer transformation not valid if b is non-positive integer, no}
          {alternative method available, report large errors in debug mode!}
          if err >= sqrt_epsh then begin
            sfd_write_debug_str('*** sfd_1f1: loss of accuracy');
          end;
        {$endif}
      end;
    end;
    exit;
  end;

  if a_ni then begin
    {Polynomial case}
    if a>THRESH_LAGUERRE_A then begin
      {Compensated sum or Laguerre}
      sfd_1f1 := h1f1_sum_laguerre(a,b,x);
      exit;
    end
    else begin
      sfd_1f1 := h1f1_sum(a,b,x,err);
      if err<ETHRESH then exit;
    end;
  end;

  if b=2.0*a then begin
    {b=2a, use Bessel_Iv}
    sfd_1f1 := h1f1_b2a(a, x, OK);
    if OK then exit;
  end;

  r := b-a;
  if is_negint(r) and (r>THRESH_LAGUERRE_A) then begin
    {a-b negative integer, Kummer transformation and compensated sum or Laguerre}
    r := h1f1_sum_laguerre(r,b,-x);
    sfd_1f1 := r*exp(x);
    exit;
  end;

  if abs(x*a) >= 30.0*abs(b) then begin
    if (abs(x) > 50.0) and (a>0.0) and (b>a) then begin
      if abs(x)>=100.0 then r := 2.0 else r := 1.0;
      if abs(a*(b-a)/x) <= r then begin
        {Muller [44] conclusion 6}
        sfd_1f1 := h1f1_acf(a,b,x);
        exit;
      end;
    end
  end;

  if x>=50.0 then begin
    r := maxd(1.0,abs(b-a)) * maxd(1.0, abs(1.0-a));
    if r < 0.5*x then begin
      sfd_1f1 := h1f1_asym(a,b,x,OK);
      if OK then exit;
    end;
  end;

  {No special method seems to apply}
  if isign(x)<>isign(a) then begin
    sfd_1f1 := h1f1_tricomi(a,b,x,OK);
    if OK then exit;
  end;
  if x>0.0 then sfd_1f1 := h1f1_sum(a,b,x,err)
  else begin
    if a<0.0 then begin
      {Try h1f1_sum first, because Luke has no error estimate}
      {and will fail dramatically e.g. for M(-1000.5, 1, -2)!}
      sfd_1f1 := h1f1_sum(a,b,x,err);
      if err<ETHRESH then exit;
    end;
    sfd_1f1 := h1f1_luke(a,b,x);
  end;

end;


{---------------------------------------------------------------------------}
function sfd_1f1r(a,b,x: double): double;
  {-Return the regularized Kummer hypergeometric function 1F1(a,b,x)/Gamma(b)}
var
  a_ni, b_ni: boolean;
  f,t: double;
  m, k: longint;
begin
  if IsNanOrInfD(a) or IsNanOrInfD(b) or IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_1f1r := NaN_d;
    exit;
  end;
  a_ni := is_negint(a);
  b_ni := is_negint(b);
  if b_ni then begin
    {NIST[30] 13.2.5}
    if (a_ni and (a-b > -0.1))  then begin
      {Pochhammer(a,m+1) is zero}
      sfd_1f1r := 0.0;
      exit;
    end;
    m := -round(b);
    t := a*x;
    k := 1;
    while (k<=m) and (t<>0.0) do begin
      f := (a+k)/(k+1)*x;
      t := t*f;
      inc(k);
    end;
    if t=0.0 then sfd_1f1r := 0.0
    else begin
      f := sfd_1f1(a-b+1.0, 2.0-b, x);
      sfd_1f1r := f*t;
    end;
  end
  else begin
    f := sfd_1f1(a,b,x);
    t := sfd_rgamma(b);
    sfd_1f1r := f*t;
  end;
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function h2f1_sum(a,b,c,x: double; var loss: double): double;
  {-Power series or polynomial summation of Gauss hypergeometric function;}
  { loss = estimated relative error}
var
  m,s,u,umax: double;
  e,w,y: double;
  i: integer;
  done: boolean;
const
  IMAX = 32000;
  lfac = 0.03125;
begin
  if abs(c)<IEPS then begin
    loss := 1.0;
    h2f1_sum:= PosInf_d;
    exit;
  end;
  i := 0;
  s := 1.0;
  u := 1.0;
  e := 0.0;
  umax := 0.0;
  done := false;
  repeat
    {Next term}
    w := (a+i)*(b+i);
    y := (c+i)*(i+1);
    if y=0.0 then begin
      if w=0.0 then u := 0.0
      else begin
        {w/0 is to be computed! Signal no convergence:}
        i := IMAX+1;
        u := 0.0;
      end;
    end
    else begin
      y := x*(w/y);
      w := abs(y);
      if w > 1.0 then begin
        if abs(u) >= MaxDouble/w then begin
          {Overflow; Signal no convergence:}
          i := IMAX+1;
          u := 0.0;
        end;
      end;
      u := u*y;
    end;
    if u<>0.0 then begin
      {remember largest term summed}
      m := abs(u);
      if m>umax then umax := m;
      {(w,y) = TwoSum(u,s)}
      w := u+s;
      m := w-u;
      y := (u-(w-m)) + (s-m);
      {sum of errors}
      e := e+y;
      s := w;
      if abs(u)<eps_d*abs(s) then begin
        {Small term, but check if c+i > 0  (if |c| is not too large)}
        {I.e. keep on summing if c+i may cross zero for larger index}
        {This is especially important if an accidental small term is}
        {encountered in the polynomial case when u=0 for a larger i.}
        done := (c < -IMAX) or (c+i > 0.0);
      end;
      inc(i);
    end
    else done := true;
  until done or (i>IMAX);
  {loss = estimated relative error}
  s := s+e;
  h2f1_sum := s;
  if s=0.0 then s:=1.0;
  if i>IMAX then loss := 1.0
  else begin
    loss := maxd(0.5*i,lfac*umax/abs(s));
    loss := eps_d*maxd(1.0,loss);
  end;
end;


{---------------------------------------------------------------------------}
function h2f1_luke(a,b,c,x: double; var loss: double): double;
  {-Return Luke's rational approximation for 2F1, (x<0?)}
const
  NMAX = 10000;
const
  sc_large = 1152921504606846976.0;  {2^60}
  sc_small = 1.0/sc_large;
var
  z,z3,t,h,e,e0: double;
  ab,an,a2,a1,a0,bn,b2,b1,b0: double;
  na,nb,nc,f1,f2,f3,err: double;
  n: integer;
begin
  {Ref: Luke [42], Ch. 13: Rational Approximations for 2F1(a,b;c;-z). }
  {Powerful but somewhat obscure approximation. Assumes c<>a and c<>b.}
  {Error estimation is not reliable, use only if other methods failed.}

  z  := -x;
  z3 := z*z*z;
  ab := a+b;

  e  := a*b/c;
  t  := 0.5*(a+1.0)*(b+1.0)/c;
  h  := 0.5*(a+2.0)*(b+2.0)/(c+1.0);

  {Starting values, [42] 13(6)}
  n  := 3;
  a0 := 1.0;
  b0 := 1.0;
  e0 := 1.0;
  b1 := 1.0 + t*z;
  b2 := 1.0 + h*z*(1.0 + t*z/THREE);
  a1 := b1 - e*z;
  a2 := b2 - e*(1.0 + h*z)*z + e*t*(c/(c+1.0))*z*z;

  repeat
    {A_n and B_n via recurrence formulas [42] 13(7)}
    na := n-1+a;
    nb := n-1+b;
    nc := n-1+c;
    an := (3.0*n)*n;
    h  := n-2+c;
    t  := 2*n-3;
    e  := -0.5*na*nb*(n-1-c)/(t*h*nc);

    f1 := 1.0 + 0.5*z*(an + (ab-6.0)*n + 2.0 - a*b - 2.0*ab)/(t*nc);
    f2 := z*(e-0.25*z*(an - (ab+6.0)*n + 2.0 - a*b)*na*nb / ((2*n-1)*t*h*nc));
    f3 := 0.125*z3*((na-1.0)*na*(nb-1.0)*nb*(n-2-a)*(n-2-b)) / (t*t*(2*n-5)*(n-3+c)*h*nc);

    an := f1*a2 + f2*a1 + f3*a0;
    bn := f1*b2 + f2*b1 + f3*b0;
    e  := an/bn;

    if e0=0.0 then err := abs(e)
    else err := abs((e0 - e)/e0);

    if (err < eps_d) or (n > NMAX) then begin
      h2f1_luke := e;
      loss := 2.0*(err + eps_d*n);
      exit;
    end;

    if (abs(an) > sc_large) or (abs(bn) > sc_large) then begin
      an := an * sc_small;
      bn := bn * sc_small;
      a2 := a2 * sc_small;
      b2 := b2 * sc_small;
      a1 := a1 * sc_small;
      b1 := b1 * sc_small;
    end
    else if (abs(an) < sc_small) or (abs(bn) < sc_small) then begin
      an := an * sc_large;
      bn := bn * sc_large;
      a2 := a2 * sc_large;
      b2 := b2 * sc_large;
      a1 := a1 * sc_large;
      b1 := b1 * sc_large;
    end;

    inc(n);
    e0 := e;
    b0 := b1;
    b1 := b2;
    b2 := bn;
    a0 := a1;
    a1 := a2;
    a2 := an;
  until false;

end;


{---------------------------------------------------------------------------}
function h2f1_trans(a,b,c,x: double; var loss: double): double;
  {-Conditionally apply linear transformation and/or call the power series}
var
  p,q,r,s,t,y,d,ax,d1,d2,e,y1,l2,r1,l1: double;
  i,id,aid: longint;
  s1,s2,s3,s4: integer;
const
  TOL   = 0.44E-14;   {~ 20*eps_d}
  TDSUM = 1000000.0;  {use h2f1_sum if |c-a-b| >= TDSUM}

  function reflectx: double;
    {-Compute with s=1-x via HMF[1], 15.3.6}
  begin
    q := h2f1_sum(a, b, 1.0-d, s, l2);
    e := sfd_lngammas(c,s4);
    t := e + sfd_lngammas(d,s1) - (sfd_lngammas(c-a,s2) + sfd_lngammas(c-b,s3));
    q := q*exp(t)*(s1*s2*s3*s4);
    r := power(s,d)*h2f1_sum(c-a, c-b, 1.0+d, s, loss);
    t := e + sfd_lngammas(-d,s1) - (sfd_lngammas(a,s2)+sfd_lngammas(b,s3));
    r := r*exp(t)*(s1*s2*s3*s4);
    y := q+r;
    {estimate combined error}
    q := abs(q);
    r := abs(r);
    if q>r then r := q;
    loss := loss + l2 + eps_d*r/y;
    reflectx := y;
  end;

begin
  {Based on the Cephes[7] function hyt2f1 in hyp2f1.c}
  {0 < x < 1; a,b,c no negative integer}
  d := c-a-b;
  s := 1.0 - x;
  if x <= 0.75 then begin
    if (a<0.0) or (b<0.0) then begin
      p := c-a;
      q := c-b;
      if ((p>0.0) and (p<a)) or ((q>0.0) and (q<b)) then begin
        {Use linear transformation HMF[1], 15.3.3}
        t := h2f1_sum(p, q, c, x, loss);
        r := power(s,d);
        h2f1_trans := r*t;
        exit;
      end
    end;
    h2f1_trans := h2f1_sum(a, b, c, x, loss)
  end
  else begin
    {test for integer c-a-b}
    if abs(d-rint(d)) > IEPS then begin
      {Use series in x or two series in s}
      if s<=0.0625 then begin
        {Try two series first for small s}
        r1 := reflectx;
        l1 := abs(loss);
        if abs(loss) < 0.25*ETHRESH then h2f1_trans := r1
        else begin
          h2f1_trans := h2f1_sum(a,b,c,x,loss);
          if l1 < abs(loss) then begin
            h2f1_trans := r1;
            loss := l1;
          end;
        end;
      end
      else begin
        r1 := h2f1_sum(a, b, c, x, loss);
        l1 := abs(loss);
        if abs(loss) < 0.25*ETHRESH then h2f1_trans := r1
        else begin
          h2f1_trans := reflectx;
          if l1 < abs(loss) then begin
            h2f1_trans := r1;
            loss := l1;
          end;
        end;
      end;
    end
    else if abs(d) > TDSUM then begin
      {d is (appr) integer but too large, give up and use sum}
      h2f1_trans := h2f1_sum(a, b, c, x, loss);
    end
    else begin
      {d is an (appr.) integer, use psi function: HMF[1] 15.3.10/11/12}
      id  := round(d);
      aid := abs(id);
      if id>=0.0 then begin
        e  := d;
        d1 := d;
        d2 := 0.0;
      end
      else begin
        e  := -d;
        d1 := 0.0;
        d2 := d;
      end;
      ax := ln(s);
      {sum for t = 0}
      r := sfd_psi(1.0+e) - Eulergamma - sfd_psi(a+d1) - sfd_psi(b+d1) - ax;
      y := r*sfd_rgamma(e+1.0);
      p := (a+d1)*(b+d1)*s*sfd_rgamma(e+2.0);      {Poch for t=1}
      t := 1.0;
      repeat
        {recursion for psi(z+t+1) = psi(z+t) + 1/(z+t)}
        {r := psi(1+t) + psi(1+t+e) - psi(a+t+d1) - psi(b+t+d1) - ax;}
        r := r + 1.0/t + 1.0/(t+e) - 1.0/(a+t+d1-1.0) - 1.0/(b+t+d1-1.0);
        q := p*r;
        y := y+q;
        p := p*(s*(a+t+d1)/(t+1.0));
        p := p*((b+t+d1)/(t+1.0+e));
        t := t+1.0;
      until abs(q)<=TOL*abs(y);
      if id=0.0 then begin
        t := sfd_lngammas(c,s1) - (sfd_lngammas(a,s2)+sfd_lngammas(b,s3));
        h2f1_trans := y*exp(t)*(s1*s2*s3);
        loss := 0.0;
        exit;
      end;
      y1 := 1.0;
      if aid<>1 then begin
        t := 0.0;
        p := 1.0;
        for i:=1 to aid-1 do begin
          r := 1.0-e+t;
          p := p*(s*(a+t+d2)*(b+t+d2)/r);
          t := t+1.0;
          p := p/t;
          y1:= y1+p;
        end;
      end;
      p  := sfd_lngammas(c,s4);
      t  := sfd_lngammas(e,s1) - (sfd_lngammas(a+d1,s2)+sfd_lngammas(b+d1,s3));
      y1 := y1*exp(p+t)*(s1*s2*s3*s4);
      q  := sfd_rgamma(a+d2)*sfd_rgamma(b+d2);
      y  := y*q;
      if y<>0.0 then begin
        if odd(aid) then y := -y;
        y := y*exp(p)*s4;
      end;
      q := power(s, id);
      if id>0.0 then y := y*q
      else y1 := y1*q;
      h2f1_trans := y+y1;
      loss := 0.0;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function h2f1_ptprefix(m: longint; b,c,z: double): double;
  {-Return prefix for 2F1 polynomial transformations, result = (b)_m/(c)_m*z^m}
var
  f,t: double;
  k: longint;
begin
  t := 1.0;
  if z=1.0 then begin
    for k:=0 to m-1 do t := ((b+k)/(c+k))*t;
  end
  else begin
    for k:=0 to m-1 do begin
      f := (b+k)/(c+k)*z;
      t := t*f;
    end;
  end;
  h2f1_ptprefix := t;
end;


{---------------------------------------------------------------------------}
function h2f1poly(a,b,c,x: double; var err: double): double;
  {-Return h2f1 for polynomial case a=-m}
var
  t,f,r: double;
  m: longint;
begin
  m := -round(a);
  err := 1.0;

  {$ifdef debug}
    if abs(a+m) >= IEPS then begin
      sfd_write_debug_str('a no negative integer');
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      h2f1poly := NaN_d;
      exit;
    end;
  {$endif}

  if x <= -2.0 then begin
    h2f1poly := h2f1_sum(a,b,c,x, err);
    if err>ETHRESH then begin
      {Transformation NIST[30] 15.8.6 (1st)}
      f := h2f1_sum(a, (1-m)-c, (1-m)-b, 1.0/x, r);
      r := r+m*eps_d;
      if r<err then begin
        t := h2f1_ptprefix(m,b,c,-x);
        err := r;
        h2f1poly := t*f;
      end;
    end;
  end
  else if x <= -0.5 then begin
    { -2 < x <= -0.5}
    h2f1poly := h2f1_sum(a,b,c,x, err);
    if err>ETHRESH then begin
      {Transformation NIST[30] 15.8.6 (2nd)}
      f := h2f1_sum(a, c-b, (1-m)-b, 1.0/(1.0-x), r);
      r := r+m*eps_d;
      if r<err then begin
        t := h2f1_ptprefix(m,b,c,1.0-x);
        err := r;
        h2f1poly := t*f;
      end;
    end;
  end
  else if x < 0.5 then begin
    {-0.5 < x < 0.5}
    h2f1poly := h2f1_sum(a,b,c,x, err);
  end
  else if x < 1.0 then begin
    {0.5 <= x < 1}
    h2f1poly := h2f1_sum(a,b,c,x, err);
    if err>ETHRESH then begin
      {Transformation NIST[30] 15.8.7 (1st)}
      f := h2f1_sum(a, b, (1-m)+b-c, 1.0-x, r);
      r := r+m*eps_d;
      if r<err then begin
        t := h2f1_ptprefix(m,c-b,c,1.0);
        err := r;
        h2f1poly := t*f;
      end;
    end
  end
  else if x=1.0 then begin
    {special case: do transformation first because h2f1_sum=0}
    err := m*eps_d;
    h2f1poly := h2f1_ptprefix(m,c-b,c,x);
    if err>ETHRESH then begin
      f := h2f1_sum(a,b,c,x, r);
      if r<err then begin
        err := r;
        h2f1poly := f;
      end;
    end
  end
  else begin
    {x > 1}
    h2f1poly := h2f1_sum(a,b,c,x, err);
    if err>ETHRESH then begin
      {Transformation NIST[30] 15.8.7 (2nd)}
      f := h2f1_sum(a, (1-m)-c, (1-m)+b-c, 1.0 - 1.0/x, r);
      r := r+m*eps_d;
      if r<err then begin
        t := h2f1_ptprefix(m,c-b,c,x);
        err := r;
        h2f1poly := t*f;
      end;
    end
  end;
end;


{---------------------------------------------------------------------------}
function hyp2f1e(a,b,c,x: double; var err: double): double;
  {-Return the Gauss hypergeometric function 2F1(a,b;c;x)}
var
  a_ni, b_ni, c_ni: boolean;
  d,r,s,t: double;
begin

  {Based on the Cephes[7] function hyp2f1 in hyp2f1.c}
  err := 0.0;

  {Handle easy special cases}
  if (x=0.0) or (a=0.0) or (b=0.0) then begin
    hyp2f1e := 1.0;
    exit;
  end;

  d := c-a-b;
  s := 1.0-x;

  if (x<1.0) and ((abs(a-c)<IEPS) or (abs(b-c)<IEPS)) then begin
    {Handle degenerated function for x<1, follows from HMF[1], 15.3.3}
    hyp2f1e := power(s, d);;
    exit;
  end;

  a_ni := is_negint(a);
  b_ni := is_negint(b);
  c_ni := is_negint(c);

  if a_ni and b_ni and (b>a) then begin
    {make b <= a}
    t := a;
    a := b;
    b := t;
  end;

  {Handle cases with c and a or b negative integer}
  if c_ni then begin
    if a_ni and (a>=c) then hyp2f1e := h2f1poly(a,b,c,x,err)
    else if b_ni and (b>=c) then hyp2f1e := h2f1poly(b,a,c,x,err)
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      hyp2f1e := NaN_d;
    end;
    exit;
  end;

  {Handle remaining cases with a,b negative integer}
  if a_ni then begin
    {if b_ni then b <= a}
    hyp2f1e := h2f1poly(a,b,c,x,err);
    exit;
  end
  else if b_ni then begin
    hyp2f1e := h2f1poly(b,a,c,x,err);
    exit;
  end;

  {Handle case with c-a or c-b negative integer}
  a_ni := is_negint(c-a);
  b_ni := is_negint(c-b);
  if a_ni or b_ni then begin
    {Use linear transformation HMF[1], 15.3.3}
    a := c-a;
    b := c-b;
    if a_ni and b_ni then begin
      if a>=b then t := h2f1poly(a,b,c,x,err)
      else t := h2f1poly(b,a,c,x,err);
    end
    else if a_ni then t := h2f1poly(a,b,c,x,err)
    else t := h2f1poly(b,a,c,x,err); {b_ni}
    r := power(s,d);
    hyp2f1e := r*t;
    exit;
  end;

  {Handle x~1}
  if abs(s) < IEPS then begin
    {x~1, use HMF[1], 15.1.20}
    if d<=0.0 then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      err := 1.0;
      hyp2f1e := NaN_d;
    end
    else begin
      t := sfd_rgamma(c-a)*sfd_rgamma(c-b);
      if t<>0.0 then begin
        {d and signgamma(d) >0}
        r:= sfd_lngamma(c)+sfd_lngamma(d);
        t := t*exp(r)*sfd_signgamma(c);
      end;
      hyp2f1e := t;
    end;
    exit;
  end;

  if x<0.0 then begin
    {Use linear transformations HMF[1], 15.3.4 or 15.3.5}
    if b>a then begin
      t := h2f1_trans(a, c-b, c, -x/s, err);
      r := power(s, -a);
    end
    else begin
      t := h2f1_trans(c-a, b, c, -x/s, err);
      r := power(s, -b);
    end;
    hyp2f1e := r*t;
    if err>ETHRESH then begin
      t := h2f1_luke(a,b,c,x,r);
      if r < err then begin
        hyp2f1e := t;
        err := r;
      end;
    end;
  end
  else if x>1.0 then begin
    {x > 1 is allowed only if the function is a polynomial}
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    err := 1.0;
    hyp2f1e := NaN_d;
    exit;
  end
  else begin
    {use direct sum or linear transformation HMF[1] 15.3.6}
    hyp2f1e := h2f1_trans(a,b,c,x,err);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_2f1(a,b,c,x: double): double;
  {-Return the Gauss hypergeometric function 2F1(a,b;c;x)}
var
  err: double;
begin
  if IsNanOrInfD(a) or IsNanOrInfD(b) or IsNanOrInfD(c) or IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_2f1 := NaN_d;
    exit;
  end;
  sfd_2f1 := hyp2f1e(a,b,c,x,err);
  {$ifdef debug}
    if err >= sqrt_epsh then begin
      sfd_write_debug_str('*** sfd_2f1: loss of accuracy');
    end;
  {$endif}
end;


{---------------------------------------------------------------------------}
function sfd_2f1r(a,b,c,x: double): double;
  {-Return the regularized Gauss hypergeometric function 2F1(a,b,c,x)/Gamma(c)}
var
  f,t: double;
  m,k: longint;
begin
  if IsNanOrInfD(a) or IsNanOrInfD(b) or IsNanOrInfD(c) or IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_2f1r := NaN_d;
    exit;
  end;
  if is_negint(c) then begin
    {HMF[1] 15.1.2}
    if (is_negint(a) and (a-c > -0.1)) or (is_negint(b) and (b-c > -0.1)) then begin
      {at least one of the Pochhammer symbols is zero}
      sfd_2f1r := 0.0;
      exit;
    end;
    m := -round(c);
    t := a*b*x;
    k := 1;
    while (k<=m) and (t<>0.0) do begin
      f := (a+k)*(b+k)/(k+1)*x;
      t := t*f;
      inc(k);
    end;
    if t=0.0 then sfd_2f1r := 0.0
    else begin
      f := sfd_2f1(a-c+1.0, b-c+1.0, 2.0-c, x);
      sfd_2f1r := f*t;
    end;
  end
  else begin
    f := sfd_2f1(a,b,c,x);
    t := sfd_rgamma(c);
    sfd_2f1r := f*t;
  end;
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

const
  IMAX2F0 = 10000;     { Max i in 2F0 integer version}
  MAXAN1  = 30000;     { Max n for chu_an1}
  T_XMILL = 1.0;       { min x for Miller algorithm, note: Temme uses 1.4}
  T_MAXA  = 150;       { Temmes code underflows for a >~ T_MAXA ~ MAXGAMD}
  T_EPS   = 1.33e-15;  { ~ 6*eps_d}


{---------------------------------------------------------------------------}
function chu_an1(a,x: double; n: longint): double;
  {-Return U(a,a+n+1,x), n>0}
var
  s,t: double;
  k: longint;
begin
  {NIST[30], 13.2.8}
  t := power(x,-a);
  s := t;
  for k:=0 to n-1 do begin
    t := t*(n-k)/(k+1); {update binom(n,k)    }
    t := t*(a+k)/x;     {update (a)_k/x^(a+k) }
    s := s + t;
  end;
  chu_an1 := s;
end;


{---------------------------------------------------------------------------}
function u2f0(a,b,x: double): double;
  {-Evaluate U(a,b,x) = 2F0(a,a-b+1;;-1/x); a or b negative integers > -IMAX2F0}
var
  z,t,s,q: double;
  i: integer;
begin
  {NIST[30]. 13.7.3/4}
  t := 1.0;
  s := 1.0;
  z := -x;
  i := 0;
  while (i <= IMAX2F0) and (t<>0.0) do begin
    q := ((a+i)/z)*((b+i)/(i+1));
    t := t*q;
    s := s+t;
    inc(i);
  end;
  u2f0 := s;
end;


{---------------------------------------------------------------------------}
function d9chu(a,b,z: double; var Err: integer): double;
  {-Return Luke's rational approximation for z^a*U(a,b,z)}
  { Err = 0: OK; = 1: loss of precision; = 2: No convergence}
var
  ab,anbn,bp,ct1,ct2,ct3,c2,d1z,eps,g1,g2,g3,sab,sqeps,x2i1,res: double;
  i,j: integer;
var
  aa,bb: array[1..4] of double;
const
  MAXIT = 300;
  IMAX  = -IMAX2F0+1;
begin

  {Ref: Luke[42], Ch. XX Rational Approximations for z^a*U(a;1+a-b;z).}
  {See also NIST[30], Ch. 13.3(iii) Rational Approximations.          }

  bp := 1.0 + a - b;
  {*WE: Luke assumes that a and a+1-b are no negative integers. In this case }
  {the asymptotic expansion is a polynomial in 1/z, which is computed in u2f0}
  if (is_negint(bp) and (bp>IMAX)) or (is_negint(a) and (a>IMAX)) then begin
    d9chu := u2f0(a,bp,z);
    Err := 0;
    exit;
  end;

  {Pascal translation of W. Fullerton's [14,20] Fortran function D9CHU.F}

  {Evaluate for large z  z^a*U(a,b,z) where U is the logarithmic confluent}
  {hypergeometric function. A rational approximation due to Y.  L. Luke is}
  {used. When U is not in the asymptotic region, i.e. when a or b is large}
  {compared with z, considerable significance loss occurs. A warning is   }
  {provided when the computed result is less than half precision.         }

  eps   := 4*eps_d;
  err   := 0;
  sqeps := sqrt(eps_d);

  ab    := a*bp;
  ct2   := 2.0*(z-ab);
  sab   := a + bp;

  bb[1] := 1.0;
  aa[1] := 1.0;

  ct3   := sab + 1.0 + ab;
  bb[2] := 1.0 + 2.0*z/ct3;
  aa[2] := 1.0 + ct2/ct3;

  anbn  := ct3 + sab + 3.0;
  ct1   := 1.0 + 2.0*z/anbn;
  bb[3] := 1.0 + 6.0*ct1*z/ct3;
  aa[3] := 1.0 + 6.0*ab/anbn + 3.0*ct1*ct2/ct3;

  for i:=4 to MAXIT do begin
    x2i1:= 2*i - 3;
    ct1 := x2i1/(x2i1-2.0);
    anbn:= anbn + x2i1 + sab;
    ct2 := (x2i1 - 1.0)/anbn;
    c2  := x2i1*ct2 - 1.0;
    d1z := x2i1*2.0*z/anbn;

    ct3 := sab*ct2;
    g1  := d1z + ct1*(c2+ct3);
    g2  := d1z - c2;
    g3  := ct1*(1.0 - ct3 - 2.0*ct2);

    bb[4] := g1*bb[3] + g2*bb[2] + g3*bb[1];
    aa[4] := g1*aa[3] + g2*aa[2] + g3*aa[1];

    res := aa[4]/bb[4];
    if abs(res-aa[1]/bb[1]) < eps*abs(res) then begin
      d9chu := res;
      if (res<sqeps) or (res >1.0/sqeps) then Err := 1;
      exit;
    end;
    {if overflows or underflows prove to be a problem, the statements below}
    {could be altered to incorporate a dynamically adjusted scale factor.  }
    for j:=1 to 3 do begin
      aa[j] := aa[j+1];
      bb[j] := bb[j+1];
    end;
  end;
  {No convergence}
  Err := 2;
  d9chu := res;
end;


{---------------------------------------------------------------------------}
function chu_luke(a,b,x: double; var Err: integer): double;
  {-Check if Luke can be used: If yes and result is OK then Err=0}
var
  t: double;
begin
  t := maxd(abs(a),1.0)*maxd(abs(1.0+a-b),1.0);
  if t < 0.9*abs(x) then begin
    {use Luke's rational approximation in the asymptotic region}
    t := d9chu(a,b,x,Err);
    if Err=0 then t := power(x,-a)*t;
  end
  else Err := -1;
  chu_luke := t;
end;


{---------------------------------------------------------------------------}
function dchu(a,b,x: double): double;
  {-Return Tricomi's confluent hypergeometric function U(a,b,x), x>0}
var
  a0,aintb,alnx,b0,beps,c0,eps,factor,gamri1,gamrni,
  pch1ai,pch1i,pochai,sum,t,value,xeps1,xi,xi1,xn,xtoeps: double;
  i,istrt,m,n: integer;
const
  MAXIT = 300;
begin
  {Pascal translation of W. Fullerton's [14,20] Fortran function DCHU.F}

{$ifdef standalone_dchu}
  {dchu is called after chu_luke has already been called from sfc_chu;}
  {avoid a second call, if dchu is not used as a standalone function. }
  {Try Luke and exit if OK}
  dchu := chu_luke(a,b,x,i);
  if i=0 then exit;
{$endif}

  {The ascending series will be used, because the descending rational  }
  {approximation (which is based on the asymptotic series) is unstable.}
  eps := 0.5*eps_d;
  if b >= 0.0 then aintb := trunc(b + 0.5)
  else aintb := trunc(b-0.5);

  beps := b - aintb;
  n := round(aintb);

  alnx := ln(x);
  xtoeps := exp(-beps*alnx);

  {Evaluate the finite sum.}
  if n<1 then begin
    {Consider the case b < 1.0 first}
    sum := 1.0;
    t := 1.0;
    m := - n;
    for i:=1 to m do begin
      xi1 := i-1;
      t := t*((a+xi1)*x/(i*(b+xi1)));
      sum := sum + t;
    end;
    sum := sfd_pochhammer(1.0+a-b,-a)*sum;
  end
  else begin
    {Now consider the case b >= 1}
    m := n - 2;
    if m<0 then sum := 0.0
    else begin
      t := 1.0;
      sum := 1.0;
      for i:=1 to m do begin
        t := t*((a-b+i)*x/(((1+i)-b)*i));
        sum := sum + t;
      end;
      t := sfd_lngammas(b-1,i);
      t := t - sfd_lngammas(a,m);
      t := t + (1-n)*alnx;
      t := (i*m)*exp(t);
      sum := t*xtoeps*sum;
    end;
  end;

  {Next evaluate the infinite sum}
  if n<1 then istrt := 1 - n
  else istrt := 0;
  xi := istrt;

  {factor := sfd_rgamma(1.0+a-b)*power(x,xi);}
  if (xi=0.0) or (x=1.0) then factor := sfd_rgamma(1.0+a-b)
  else begin
    t := sfd_lngammas(1.0+a-b,i);
    factor := i*exp(xi*alnx-t);
  end;

  if odd(n) then factor := -factor;
  if beps<>0.0 then factor := factor*beps*Pi/sinPi(beps);

  pochai := sfd_pochhammer(a,xi);
  gamri1 := sfd_rgamma(xi+1.0);
  gamrni := sfd_rgamma(aintb+xi);
  b0 := factor*sfd_pochhammer(a,xi-beps);
  b0 := b0*gamrni*sfd_rgamma(xi+1.0-beps);

  if abs(xtoeps-1.0) <= 0.5  then begin
    {x^(-beps) is close to 1, so we must be careful in evaluating the differences}
    pch1ai := sfd_poch1(a+xi, -beps);
    pch1i  := sfd_poch1(xi+1.0-beps, beps);
    t  := pch1ai - sfd_poch1(b+xi, -beps) - pch1i + beps*pch1ai*pch1i;
    c0 := factor*pochai*gamrni*gamri1*t;

    {xeps1 = (1 - x^(-beps))/beps = (x^(-beps) - 1)/(-beps)}
    {xeps1 := -powm1(x,-beps)/beps;}
    xeps1 := alnx*exprel(-beps*alnx);
    value := sum + c0 + xeps1*b0;
    xn := n;
    for i:=1 to 1000 do begin
      xi := istrt+i;
      xi1 := istrt+i-1;
      b0 := (a+xi1-beps)*b0*x / ((xn+xi1)*(xi-beps));
      c0 := (a+xi1)*c0*x / ((b+xi1)*xi)
            - ((a-1.0)*(xn+2.0*xi-1.0) + xi*(xi-beps))*b0
             /(xi*(b+xi1)*(a+xi1-beps));
      t := c0 + xeps1*b0;
      value := value + t;
      if abs(t) < eps*abs(value) then begin
        dchu := value;
        exit;
      end;
    end;
  end
  else begin
    {x^(-beps) is very different from 1, so the straightforward formulation is stable}
    a0 := factor*pochai*sfd_rgamma(b+xi)*gamri1/beps;
    b0 := xtoeps*b0/beps;
    value := sum + a0 - b0;

    for i:=1 to MAXIT do begin
      xi  := istrt+i;
      xi1 := istrt+i-1;
      a0  := (a+xi1)*a0*x / ((b+xi1)*xi);
      b0  := (a+xi1-beps)*b0*x / ((aintb+xi1)*(xi-beps));
      t   := a0 - b0;
      value := value + t;
      if abs(t) < eps*abs(value) then begin
        dchu := value;
        exit;
      end;
    end;
  end;
  {$ifdef debug}
    sfd_write_debug_str('No convergence in MAXIT terms of the ascending series');
  {$endif}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  dchu := value;
end;


{---------------------------------------------------------------------------}
{----------------- Functions from Temme's CHU code -------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure tchu_miller(a,b,x: double; var u,up: double);
  {-Compute U, U' with Temme's Miller algorithm; a>=0, 0<=b<=1, x>=T_XMILL}
var
  ar,br,cr,c,er,m0,m1,mr,p0,p1,p2,u1,u2,u3,v,w: double;
  n,r: longint;
  largex: boolean;

  procedure recursion;
  begin
    p2 := (br*p1 - ar*p0)/cr;
    er := er*ar/cr;
    r  := r + 1;
    if largex then mr := mr*(1+c/r);
    v  := er/p2;
    br := br + 2.0;
    cr := cr + 1.0;
    p0 := p1;
    p1 := p2;
  end;

begin
  {Ref: Miller part (x>1.4) of Temme's [55] procedure chu.}
  {Assumes: a>=0, 0 <= b <= 1, x >= T_XMILL}
{$ifdef debug}
  if (a<0) or (b<0) or (x<T_XMILL) then begin
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
  end;
{$endif}
  n := trunc(a);
  if a=n then n := n-1;
  a := a - n;

  m0 := a;
  m1 := 1.0;
  mr := 0.0;

  largex := (x>6.5) and (a<>b);
  if largex then mr := 1.0
  else if a<>b then begin
    m0 := 0.0;
    v  := 1.0;
    r  := 1;
    while v >= m1*T_EPS do begin
      {There is a typo in Temme[55], who writes v := v*v/r!}
      v  := v*x/r;
      m0 := m0 + v;
      v  := v*(a+r)/(b+r);
      m1 := m1 + v;
      inc(r);
    end;
    v  := exp(-x)*sfd_gamma(a+1.0)/sfd_gamma(b+1.0);
    m0 := v*(b+a*m0);
    m1 := v*m1;
  end;

  c  := a-b;
  cr := 2.0+c;
  br := x+a+cr;
  p0 := 0.0;
  v  := 1.0;
  p1 := 1.0;
  er := 1.0;
  r  := 0;
  while r<=n do begin
    ar := a+r;
    recursion;
  end;
  w := p0*p1/er;
  ar := a+r;
  while v*(w/p0 + mr*(2.0+a/r)) >= T_EPS do begin
    recursion;
    ar := ar + 1.0;
  end;

  c  := 1.0+c;
  v  := x+c;
  w  := 0.0;
  u2 := 1.0;
  u3 := -2.0*r/(x+sqrt(x*(x+4*r)));
  r  := r-1;
  while r>0 do begin
    if largex then begin
      w  := w + mr*u2;
      mr := mr*(r+1)/(c+r);
    end;
    u1 := (-x*u3 + (v+r)*u2)/(a+r);
    u3 := u3 - u2;
    u2 := u1;
    if r=n then begin
      u  := u2;
      up := u3;
    end;
    r := r-1;
  end;
  u1 := -x*u3 + v*u2;
  u3 := u3 - u2;
  v  := a;
  if n=0 then up := u3;
  if largex then w := power(x,-a)/(a*(w+c*u2)+u1)
  else w := power(x,-b)/(u1*m1 - u3*m0);
  {Warning: Underflow for a > ~ MAXGAMD}
  for r:=0 to n-1 do v := v/(a+r);
  v  := v*w;
  up := v*up;
  if n=0 then u := w*u1
  else u := v*u;
end;


{---------------------------------------------------------------------------}
procedure tchu_bessel(a,b,x: double; var u, up: double);
  {-Compute U, U' with Temme's BesselK code; 0<=a, 0<=b<=1, 0<x<T_XMILL}
const
  n=9;
var
  d,delta,e,f,p,q,r,s,t,t0,t1,u0,u1,u2,u3,v,w,x2,z: double;
  bb,bb1,bx: array[0..8] of double;
  fi: array[-1..8] of double;
  i: integer;
  nu,j: longint;
begin
  {Ref: BesselK part (x<=1.4) of Temme's [55] procedure chu.}
  {Assumes: 0 <= a, 0 <= b <= 1, 0 < x < T_XMILL}

{$ifdef debug}
  if (a<0.0) or (b<0.0) or (b>1.0) or (x<=0.0) or (x>T_XMILL) then begin
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
  end;
{$endif}

  v := 12.56637;   {~ 4*Pi}
  r := (x - v*(b+1))/2;
  d := (v*n + r - sqrt(r*r + 4*n*x))/(2*n);
  if d<4.7124 then begin
     d := 4.7124;  {~ 1.5*Pi}
     v := 1.0;
     w := 1.0;
  end
  else begin
    v := abs(sin(d));
    w := power(v,-1.0 - b);
  end;
  w := w*exp(0.5*x*(1.0/v + 1.0/d));
  delta := T_EPS*exp(-0.5*x + (n-1-b)*ln(d))/w;
  z := 0.5/delta;
  v := 0.5 - b;
  for i:=1 to n-1 do z := z*(i+v);

  i := 0;
  t := sqrt(x)*power(z, 0.5/n);
  e := ln(delta) + n - n*ln(x);

  repeat
    r := n + b + t;
    s := 1.0 + b + t;
    p := ln(r);
    q := ln(s);
    f := (1.0/(t+0.5) - 2*n/t + 0.5*(n-1)/(r*s) + p - q);
    f := (ln(t+0.5) - 2*n*ln(t) + (r-0.5)*p - (s-0.5)*q - e) / f;
    if f < 0.0 then begin
      t := t - f;
      t := sqrt(2.0*x + t*t);
      i := i+1;
    end;
  until (f>=0.0) or (i>=9);

  nu := 1 + floor(t*t/x -a);
  if nu<=0 then nu := 1;
  r := a + nu;
  w := sqrt(x/r);
  v := 2.0*r*w;
  sfd_bess_kv2(-b,v,t0,t1);

  v  := power(w, -b);
  u1 := v*t0;
  u0 := v*w*t1;
  x2 := x*x;
  fi[-1] := u1;
  fi[0]  := u0;
  bb1[0] := 1.0;

  bx[0] := 1.0;
  bx[1] := -x/12;
  bx[2] :=  x2/288;
  bx[3] := -x*(5*x2 - 72)/51840.0;
  bx[4] :=  x2*(5*x2 - 288)/2488320.0;
  bx[5] := -x*(x2*(7*x2 - 1008) + 6912)/209018880.0;
  bx[6] :=  x2*(x2*(35*x2 - 10080) + 279936.0)/75246796800.0;
  bx[7] := -x*(x2*(x2*(5*x2 - 2520) + 176256.0) - 746496.0)/902961561600.0;
  bx[8] :=  x2*(x2*(x2*(5*x2 - 4032) + 566784.0) - 9953280.0)/86684309913600.0;

  bb[0] := 1.0;
  bb[1] := 0.5;
  bb[2] := (3*b-1)/24;
  bb[3] := b*(b-1)/48;
  bb[4] := (b*(b*(15*b - 30) + 5) + 2)/5760;
  bb[5] := b*(b*(b*(3*b - 10) + 5) + 2)/11520;
  bb[6] := (b*(b*(b*(b*(63*b - 315) + 315) + 91) - 42) - 16)/2903040.0;
  bb[7] := b*(b*(b*(b*(b*(9*b - 63) + 105) + 7) - 42) - 16)/5806080.0;
  bb[8] := (b*(b*(b*(b*(b*(b*(135*b - 1260) + 3150) - 840) - 2345) - 540) + 404) + 144)/1393459200.0;

  for i:=1 to n-1 do begin
    t0 := bb[i];
    t1 := (b-i)*t0;
    bb1[i] := t1;
    fi[i]  := (x*fi[i-2] + (i-b)*fi[i-1])/r;
    for j:=1 to i-1 do begin
      t0 := t0 +  bb[i-j]*bx[j];
      t1 := t1 + bb1[i-j]*bx[j];
    end;
    t0 := bx[i] + b*t0;
    t1 := t1 + bx[i];
    u0 := u0 + t0*fi[i];
    u1 := u1 + t1*fi[i-1];
  end;

  w  := 2.0*exp(0.5*x)/sfd_gamma(1.0+a);
  u2 := w*u0;
  u3 := -w*u1;
  v  := a+1.0-b+x;
  for j:=nu-1 downto 1 do begin
    u1 := (-x*u3 + (v+j)*u2)/(a+j);
    u3 := u3 - u2;
    u2 := u1;
  end;
  u  := -x*u3+v*u2;
  up := a*(u3-u2);
end;


{---------------------------------------------------------------------------}
procedure tchu(a,b,x: double; var u,up: double);
  {-Compute u = U(a,b,x), up = U'(a,b,x) for a>=0, x>0 with Temme's method}
var
  c,d,v,w: double;
  i,k: longint;
begin
{$ifdef debug}
  if (a<0.0) or (x<=0.0) then begin
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
  end;
{$endif}
  {This is based on Temme's [55] code in procedure chu, but it}
  {is heavily modified: There is no array of function values, }
  {but procedures for separate Miller and Bessel branches.}
  if a=0.0 then begin
    u  := 1.0;
    up := 0.0;
  end
  else if b < 0.0 then begin
    tchu(a-b+1.0, 1.0-b, x, up,w);
    d  := a*power(x,-b);
    up := -d*up;
    u  := -x*(up + d*w)/a;
  end
  else if b > 1.0 then begin
    {Dramatically simplified compared to Temme, because}
    {his booleans r,s are always true for kmax=0.}
    c  := b - a - 1.0;
    k  := floor(c);
    if (c>=0.0) and (c=k) then begin
      v := power(x, -a);
      w := -a*v/x;
      d := a + 1.0;
    end
    else begin
      k := trunc(b);
      d := frac(b);
      tchu(a,d,x,v,w);
    end;
    for i:=0 to k-1 do begin
      c := v-w;
      w := ((i+d)*w - a*v)/x;
      v := c;
    end;
    u  := v;
    up := w;
  end
  else begin
    if x >= T_XMILL then tchu_miller(a,b,x,u,up)
    else tchu_bessel(a,b,x,u,up);
  end;
end;


{---------------------------------------------------------------------------}
function tuabx(a,b,x: double): double;
  {-Compute U(a,b,x) with Temme's method; x > 0}
var
  a1,c,p,q,r,u: double;
  j,n: longint;
begin
{$ifdef debug}
  if (x<=0.0) then begin
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
  end;
{$endif}
  {Pascal translation of Temme's [55] function uabx}
  if a<0 then n := floor(a) else n := 0;
  a1 := a-n;
  q  := a1;
  u  := 1.0;
  if (n<0) and (a=b) then begin
    tchu(a1,a1,x,u,q);
    p := u;
    r := p-q;
    for j:=1 to -n do begin
      r := x*r;
      q := (a1-j)*p;
      p := r-q;
    end;
  end
  else begin
    if a1>0 then begin
      tchu(a1,b,x,u,q);
    end;
    c  := 1.0 + a1 - b + x;
    {Deleted from Temme: a1 := a1 - 1;}
    p  := u;
    for j:=1 to -n do begin
      r := (c-j)*p - x*q;
      q := (a1-j)*(q-p);
      p := r;
    end;
  end;
  tuabx := p;
end;
{----------------------- End of Temme code ---------------------------------}


{---------------------------------------------------------------------------}
function drchu_alt0(a,b,x: double): double;
  {-Compute U(a,b,x) for a<0, frac(a)<>0 with double recursion on a- and b+}
var
  a1,b1,c,p,t,r: double;
  j,n,m: longint;
begin
{$ifdef debug}
  if (a>=0.0) or (x<=0.0) then begin
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
  end;
{$endif}
  if b<1.0 then begin
    {U(a,b,x) = x^(1-b) U(1+a-b,2-b,x)}
    c := drchu_alt0(1+a-b, 2-b,x);
    t := power(x,1-b);
    drchu_alt0 := c*t;
    exit;
  end;

  {Here b >= 1}
  n  := floor(a);
  m  := floor(b);
  a1 := a - n;
  b1 := b - (m-1);

  {Recursion on a-}
  tchu(a1,b1,x,p,c);
  {HMF[1], 13.4.26: U(a1-1,b1,x) = (a1-b1+x)*U(a1,b1,x) - x*U'(a1,b1,x)}
  c := (a1-b1+x)*p - x*c;
  {p = previous = U(a1,  b1,x)}
  {c = current  = U(a1-1,b1,x) }
  for j:=1 to -n-1 do begin
    {HMF[1], 13.4.15:  U(a-) + (b-2a-x)U(a) + a(a+1-b)U(a+1) = 0}
    r := a1-j;
    t := (x + 2.0*r - b1)*c - r*(r + 1.0 - b1)*p;
    p := c;
    c := t;
  end;

  if b>=2 then begin
    {Recursion on b+, this is empty if b < 2}
    {HMF[1], 13.4.19:  xU(a,b1+1,x) = (a+x)U(a,b1,x) + a(b1-a-1)U(a+1,b1,x)}
    t := (a+x)*c + a*(b1-a-1.0)*p;
    p := c;    {U(a,b1)}
    c := t/x;  {U(a,b1+1}
    for j:=0 to m-3 do begin
      {HMF[1], 13.4.16:  xU(b+) + (b-a-1)U(b) + (1-b-x)U(b-) = 0}
      r := b1+j;
      t := (a-r)*p + (x+r)*c;
      p := c;
      c := t/x;
    end;
  end;

  drchu_alt0 := c;
end;


{---------------------------------------------------------------------------}
function chu_negx(a,b,x: double; var Err: integer): double;
  {-Return U(a,b,x) for x < 0,  Err<>0 if arguments cannot be handled}
var
  t,f: double;
  m: longint;
begin
  if frac(a)=0.0 then begin
    Err := 0;
    t := a-b+1.0;
    if is_negint(t) and (t >= -MAXAN1) then begin
      {For U(a,a+n+1,x), n>0  apply NIST[30], 13.2.8. Note:}
      {this uses x^(-a) and therefore a must be an integer!}
      chu_negx := chu_an1(a,x,round(-t));
      exit;
    end;
    if a < 0.0 then begin
      if a >= -MAXAN1 then begin
        m := round(-a);
        if a=b then begin
          {Same as above but use Kummer transformation:}
          {U(a,b,x) = U(1,2-b,x)*x^(1-b)}
          f := chu_an1(1.0, x, m);
          t := power(x,1.0-b);
        end
        else begin
          {Tricomi [53], 2.4(6): U(-m,b,x) = (-1)^m * (b)_m * M(-m,b,x)}
          t := sfd_pochhammer(b,m);
          if odd(m) then t := -t;
          f := sfd_1f1(a,b,x);
        end;
        chu_negx := f*t;
        exit;
      end;
    end;
  end;
  chu_negx := NaN_d;
  Err := 1;
end;


{---------------------------------------------------------------------------}
function sfd_chu(a,b,x: double): double;
  {-Return Tricomi's confluent hypergeometric function U(a,b,x). If}
  { x<0, then a must be an integer and a<0 or 1+a-b an integer < 0.}
var
  s,t: double;
  Err: integer;
begin
  if IsNanOrInfD(a) or IsNanOrInfD(b) or IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_chu := NaN_d;
    exit;
  end;
  if a=0.0 then begin
    {Tricomi[53], 2.4(6); http://functions.wolfram.com/07.33.03.0013.01}
    sfd_chu := 1.0;
    exit;
  end
  else if a=-1.0 then begin
    {Tricomi[53], 2.4(6); http://functions.wolfram.com/07.33.03.0012.01}
    sfd_chu := x-b;
    exit;
  end
  else if x=0.0 then begin
    {http://functions.wolfram.com/07.33.03.0001.01}
    {U(a,b,0) = Gamma(1-b)/Gamma(a-b+1) for b < 1}
    sfd_chu := sfd_gamma_delta_ratio(1.0-b,a);
    exit;
  end;

  if x<0.0 then begin
    sfd_chu := chu_negx(a,b,x,Err);
    if Err<>0 then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_chu := NaN_d;
    end;
    exit;
  end;

  {Here x > 0. Try Luke (even for a or b negative integer) and exit if OK}
  sfd_chu := chu_luke(a,b,x,Err);
  if Err=0 then exit;

  t := a-b+1.0;

  if t=0.0 then begin
    {NIST[30], 13.6.4}
    sfd_chu := power(x,-a);
    exit;
  end
  else if b=2.0*a then begin
    {NIST[30], 13.6.10}
    s := sfd_kve(a-0.5,0.5*x);
    t := power(x,0.5-a)/sqrtpi;
    sfd_chu := s*t;
    exit;
  end;

  (*
  if is_negint(b) then begin
    {Lebedev [45] 9.10.9}
    m := 1-round(b);
    t := sfd_chu(a+m,m+1,x);
    s := power(x,m);
    sfd_chu := s*t;
    exit;
  end;
  *)
  (*
  if is_negint(a) then begin
    {Tricomi [53], 2.4(6)}
    m := round(-a);
    s := sfd_pochhammer(b,m);
    t := sfd_1f1(a,b,x);
    if odd(m) then sfd_chu := -s*t else sfd_chu := s*t;
    exit;
  end;
  *)
  if is_negint(t) and (frac(a)<>0) then begin
    if (a <= 0.0) or (t >= -MAXAN1) then begin
      sfd_chu := chu_an1(a,x,round(-t));
      exit;
    end;
  end;

  {No more implemented special cases: general function}
  if a<0.0 then begin
    {Temme code for a<0 should be avoided except in special cases}
    if t>0.0 then begin
      {Use Kummer transformation with positive t (= new a)}
      s := tuabx(t,2.0-b,x);
      t := power(x,1.0-b);
      sfd_chu := s*t;
    end
    else if (frac(a)<>0) or ((a<0.0) and (b<0.0)) then begin
      {Use double recursion on a- and b+}
      sfd_chu := drchu_alt0(a,b,x);
    end
    else begin
      {Temme's recursion is stable if b=a < 0. And empirically it }
      {is better than double recursion if a is a negative integer.}
      sfd_chu := tuabx(a,b,x);
    end;
  end
  else if (a>T_MAXA) and ((t>0.0) or (frac(t)<>0)) then begin
    {Temme's code will underflow. dchu extends the range slighty, but ...}
    sfd_chu := dchu(a,b,x);
  end
  else begin
    if (b<0.0) and (t>T_MAXA) then begin
      {Temme will underflow, use SLATEC with Kummer transformation, see above}
      s := dchu(t,2.0-b,x);
      t := power(x,1.0-b);
      sfd_chu := s*t;
    end
    else begin
      {Use plain Temme code for a > 0}
      sfd_chu := tuabx(a,b,x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_whitm(k,m,x: double): double;
  {-Return the Whittaker M function = exp(-x/2)*x^(0.5+m) * 1F1(m-k-0.5,2m+1,x)}
var
  a,b: double;
begin
  if IsNanOrInfD(k) or IsNanOrInfD(k) or IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_whitm := NaN_d;
    exit;
  end;
  b := m + 0.5;
  if (x<0.0) and (frac(b)<>0.0) then begin
    {x^b will fail}
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_whitm := NaN_d;
    exit;
  end;
  a := b - k;
  a := sfd_1f1(a,2.0*b,x);
  b := exp(-0.5*x)*power(x,b);
  sfd_whitm := a*b;
end;


{---------------------------------------------------------------------------}
function sfd_whitw(k,m,x: double): double;
  {-Return the Whittaker W function = exp(-x/2)*x^(0.5+m) * U(m-k-0.5,2m+1,x)}
var
  a,b: double;
begin
  if IsNanOrInfD(k) or IsNanOrInfD(k) or IsNanOrInfD(x) or (x<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_whitw := NaN_d;
    exit;
  end;
  b := m + 0.5;
  a := b - k;
  a := sfd_chu(a,2.0*b,x);
  b := exp(-0.5*x)*power(x,b);
  sfd_whitw := a*b;
end;


{---------------------------------------------------------------------------}
function cf_0f1(b,x: double; var OK: boolean): double;
  {-Evaluate continued fraction for 0F1}
var
  a,c,d,f,t,tiny,tol: double;
  k: integer;
const
  MAXIT = 1000;
begin
  {Ref: http://functions.wolfram.com/07.17.10.0001.01}
  {Evaluation with modified Lentz method}
  tol  := eps_d;
  tiny := Sqrt_MinDbl;
  f := 1.0;
  c := f;
  d := 0.0;
  k := 1;

  repeat
    if k=1 then begin
      a := x/b;
      t := 1.0;
    end
    else begin
      a := -x/k/(b+(k-1));
      t := 1.0 - a;
    end;
    d := t + a*d;
    c := t + a/c;
    if c=0.0 then c := tiny;
    if d=0.0 then d := tiny;
    d := 1.0/d;
    t := c*d;
    f := f*t;
    inc(k);
  until (abs(t-1.0) < tol) or (k>MAXIT);
  cf_0f1 := f;
  OK := (k<=MAXIT)
end;


{---------------------------------------------------------------------------}
function gen_0f1(b,x: double; reg: boolean): double;
  {-Return the (regularized if reg) confluent hypergeometric limit function}
var
  c,f,g,t: double;
  s: integer;
  UseCF, OK: boolean;
const
  XSMALL = 0.125;
begin
{$ifdef debug}
  if ((b<0.0) and (frac(b)=0.0)) then begin
    {Should be handled in calling functions}
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    gen_0f1 := NaN_d;
    exit;
  end;
{$endif}

  OK := true;

  {Preset result for x=0; http://functions.wolfram.com/07.17.03.0001.01}
  {This also makes some compilers happy}
  gen_0f1 := 1.0;
  if x=0.0 then begin
    if reg then gen_0f1 := sfd_rgamma(b);
    exit;
  end
  else if abs(x)<=XSMALL then begin
    f := cf_0f1(b,x,OK);
    if OK then begin
      if reg then f := f*sfd_rgamma(b);
      gen_0f1 := f;
      exit;
    end;
  end;

  c := b-1.0;
  t := 2.0*sqrt(abs(x));
  if x<0.0 then begin
    {http://functions.wolfram.com/07.17.27.0002.01}
    f := sfd_jv(c, t);
  end
  else begin
    {http://functions.wolfram.com/07.17.27.0003.01}
    f := sfd_iv(c, t);
  end;
  if b=1.0 then begin
    {For b=1 the gamma/power prefix is 1}
    gen_0f1 := f;
    exit;
  end;

  UseCF := IsNanOrInfD(f);

  if not UseCF then begin
    if (abs(b) < MAXGAMD - 10) then begin
      t := power(abs(x), -0.5*c);
      if reg then g := 1.0
      else g := sfd_gamma(b);
      if (t>1.0) and (abs(g) >= MaxDouble/t) then begin
        {Avoid overflow}
        UseCF := true;
      end
      else t := g*t;
    end
    else begin
      t := 0.5*c*ln(abs(x));
      if reg then begin
        s := 1;
        t := -t;
      end
      else begin
        g := sfd_lngammas(b,s);
        t := g-t;
      end;
      if abs(t)<ln_MaxDbl then t := s*exp(t)
      else begin
        {Avoid overflow or underflow}
        UseCF := true;
      end;
    end;
    if not UseCF then begin
      {here f is the Bessel function result and t is the gamma/power prefix}
      f := t*f;
      if not IsNanOrInfD(f) then begin
        gen_0f1 := f;
        exit;
      end;
      UseCF := true;
    end;
  end;
  if UseCF then begin
    {if OK is false the CF has been called already without success}
    if OK then begin
      f := cf_0f1(b,x,OK);
      if OK and reg then f := f*sfd_rgamma(b);
      gen_0f1 := f;
    end;
    {No convergence}
    if (not OK) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  end;
end;


{---------------------------------------------------------------------------}
function sfd_0f1r(b,x: double): double;
  {-Return the regularized confluent hypergeometric limit function 0F1(;b;x)/Gamma(b)}
var
  t,f: double;
begin
  if IsNanOrInfD(b) or IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_0f1r := NaN_d;
    exit;
  end;
  if (b<=0.0) and (frac(b)=0.0) then begin
    {http://functions.wolfram.com/07.18.17.0006.01}
    if x=0.0 then sfd_0f1r := 0.0
    else begin
      b := 1.0-b;
      t := power(x, b);
      if (t<>0.0) and not IsNanOrInfD(t) then begin
        f := gen_0f1(1.0+b,x,true);
        if (f<>0.0) and not IsNanOrInfD(f) then begin
          sfd_0f1r := t*f;
          exit;
        end;
      end;
      {Problem with power or 0F1r, use alternate form and logs}
      f := gen_0f1(1.0+b,x,false);
      t := sfd_lngamma(1.0+b);
      {power is negative, change sign of f}
      if (x<0.0) and (frac(0.5*b)<>0.0) then f := -f;
      t := b*ln(abs(x)) - t;
      if t<ln_MaxDbl then t := exp(t)
      else t := PosInf_d;
      sfd_0f1r := t*f;
    end;
  end
  else begin
    sfd_0f1r := gen_0f1(b,x,true);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_0f1(b,x: double): double;
  {-Return the confluent hypergeometric limit function 0F1(;b;x)}
begin
  if IsNanOrInfD(b) or IsNanOrInfD(x) or ((b<0.0) and (frac(b)=0.0)) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_0f1 := NaN_d;
    exit;
  end;
  sfd_0f1 := gen_0f1(b,x,false);
end;


{---------------------------------------------------------------------------}
function sfd_2f0(a,b,x: double): double;
  {-Return 2F0(a,b,x), if x>0 then a or b must be a negative integer}
var
  t,z: double;
  OK: boolean;
begin
  {HMF[1], 13.1.11: 2F0(a,b,x) = (-1/x)^a U(a,1+a-b,-1/x)}
  if (x=0.0) or (a=0.0) or (b=0.0) then sfd_2f0 := 1.0
  else begin
    z := -1.0/x;
    if z>0.0 then begin
      t := sfd_chu(a, 1.0+a-b, z);
      sfd_2f0 := power(z,a)*t;
    end
    else if (a < 0.0) and (frac(a)=0.0) then begin
      t := sfd_chu(a, 1.0+a-b, z);
      sfd_2f0 := power(z,a)*t;
    end
    else if (b < 0.0) and (frac(b)=0.0) then begin
      t := sfd_chu(b, 1.0+b-a, z);
      sfd_2f0 := power(z,b)*t;
    end
    else begin
      {x > 0, a and b no negative integer, use asymptotic}
      {series, but x must be relative small, max ~ 0.02}
      sfd_2f0 := hyp2f0(a,b,x,OK);
      if not OK then sfd_2f0 := Nan_d;
    end;
  end;
end;



{---------------------------------------------------------------------------}
{-------------------  Parabolic cylinder functions  ------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function sfd_pcfd(v,x: double): double;
  {-Return Whittaker's parabolic cylinder function D_v(x)}
var
  f,u,w,y,z: double;
begin
  if IsNanOrInfD(v) or IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_pcfd := NaN_d;
    exit;
  end;
  {Erdelyi[50], II 8.2, and Tricomi[53], Sec.4.8}
  z := abs(x);
  y := 0.5*x*x;
  w := 0.5*v;
  u := sfd_chu(-w, 0.5, y);
  {http://functions.wolfram.com/07.41.26.0007.01}
  if x<0.0 then begin
    if (frac(v)=0.0) and (v>0.0) then begin
      {Case v positive integer D_v(-x) = (-1)^v * D_v(x)}
      {http://functions.wolfram.com/07.41.16.0007.01}
      if frac(w) <> 0.0 then u := -u;
    end
    else if frac(w)<>0.0 then begin
      {http://functions.wolfram.com/07.41.16.0006.01}
      f := sfd_1f1(0.5-w,1.5, y);
      y := sqrt2*sinPi(w);
      f := y*f;
      y := v*sfd_gamma(w)/sqrtpi;
      u := u - z*y*f
    end;
  end;
  y := expx2(-0.5*z);
  f := exp2(w);
  sfd_pcfd := u*f*y;
end;


{---------------------------------------------------------------------------}
function sfd_pcfu(a,x: double): double;
  {-Return the parabolic cylinder function U(a,x)}
begin
  {HMF[1], 19.3.7}
  sfd_pcfu := sfd_pcfd(-(a+0.5), x);
end;


{---------------------------------------------------------------------------}
function pcfv_axpos(a,x: double): double;
  {-Compute parabolic cylinder V(a,x) with a,x positive, no checks}
var
  k: longint;
  v0,v1,vk,b: double;
  z,bi0,bk0,bi1,bk1: double;
const
  a0 = 0.686212627559326157189284947975; {Taylor: V(0,x) = a0 + b0*x + O(x^4)}
  b0 = 0.328001948666876466396411383466;
  a1 = 0.328001948666876466396411383466; {Taylor: V(1,x) = a1 + b0*x + O(x^2)}
  b1 = 0.343106313779663078594642473988;
const
  c10 = 0.450158158078553034777599595503; {sqrt(2)/Pi}
  c05 = 0.797884560802865355879892119869; {sqrt(2/Pi)}
const
  eps1 = 1.05e-8;   {~sqrt(eps)}
  eps2 = 1.05e-4;   {~eps^0.25}
begin
  {Evaluate three-term recurrence relation from HMF[1], 19.6.8}
  {First compute starting values from formulas in HMF[1], 19.5}
  if frac(a)<>0.0 then begin
    {Prepare v0=V(1/2,x) and v1=V(3/2,x) and exit if a=1/2 or 3/2}
    v0 := expx2(0.5*x)*c05;
    if a=0.5 then begin
      pcfv_axpos := v0;
      exit;
    end;
    v1 := v0*x;
    if a=1.5 then begin
      pcfv_axpos := v1;
      exit;
    end;
    b := 1.0;
  end
  else begin
    if a>0.0 then begin
      {Prepare v0=V(0,x) and v1=V(1,x) and exit if a=1}
      if x<=eps1 then begin
        v1 := a1 + b1*x;
        if a=1.0 then begin
          pcfv_axpos := v1;
          exit;
        end;
        v0 := a0 + b0*x;
      end
      else begin
        z := sqr(0.5*x);
        sfd_bess_ikv(0.25,z,bi0,bk0);
        sfd_bess_ikv(0.75,z,bi1,bk1);
        b  := sqrt(x);
        z  := 2.0*(bi0+bi1) + c10*(bk0+bk1);
        v1 := 0.25*b*x*z;
        if a=1.0 then begin
          pcfv_axpos := v1;
          exit;
        end;
        z  := 2.0*bi0 + c10*bk0;
        v0 := 0.5*b*z;
      end;
    end
    else begin
       {V(0,x) only}
       if x<=eps2 then pcfv_axpos := a0 + b0*x
       else begin
         z := sqr(0.5*x);
         sfd_bess_ikv(0.25,z,v0,vk);
         z := 2.0*v0 + c10*vk;
         pcfv_axpos := 0.5*sqrt(x)*z;
       end;
       exit;
    end;
    b := 0.5;
  end;
  {Step the recurrence in forward direction (increasing a)}
  for k:=2 to trunc(2.0*a) div 2 do begin
    vk := x*v1 + b*v0;
    v0 := v1;
    v1 := vk;
    b  := b+1.0;
  end;
  pcfv_axpos := v1;
end;


{---------------------------------------------------------------------------}
function pcfv_gen(a,x: double): double;
  {-General parabolic cylinder V(a,x), used for a<0}
var
  f,t,v0,v1,u1,u2,s,c: double;
begin
  f := exp2(0.5*a+0.25);
  {HMF[1], 19.3.6: v0 = V(a, 0)}
  sincospi(0.25-0.5*a,s,c);
  t  := sfd_rgamma(0.75-0.5*a);
  v0 := c*t;
  if x=0.0 then pcfv_gen := f*v0
  else begin
    t := 0.5*sqr(x);
    if v0<>0.0 then begin
      {NIST[30], 12.7.12: u1(a,x)}
      u1 := sfd_1f1(0.5*a+0.25,0.5,t);
      v0 := v0*u1;
    end;
    if s=0.0 then v1 := 0.0
    else begin
      {HMF[1], 19.3.6: v1 = V'(a, 0)}
      v1 := sfd_rgamma(0.25-0.5*a);
      if v1<>0.0 then begin
        {NIST[30], 12.7.13: u2(a,x)}
        u2 := x*sfd_1f1(0.5*a+0.75,1.5,t);
        v1 := sqrt2*v1*s*u2;
      end;
    end;
    t := expx2(-0.5*abs(x));
    {NIST[30], 12.4.2}
    pcfv_gen := (f*(v0 + v1))*t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_pcfv(a,x: double): double;
  {-Return the parabolic cylinder function V(a,x) with 2a integer}
var
  g,t: double;
const
  amax = 302.0;  {Max a for V(a,0) without overflow}
begin
  if IsNanOrInfD(a) or IsNanOrInfD(x) or (frac(2.0*a)<>0) or (abs(a)>amax) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_pcfv := NaN_d;
    exit;
  end;
  if a<0.0 then begin
    if (frac(a)=0.0) and (a > -MAXGAMD) then begin
      {HMF[1], 19.3.8, with sin(a*Pi)=0!}
      g := a+0.5;
      t := sfd_pcfd(-g, -x);
      g := sfd_gamma(g);
      sfd_pcfv := t*g/Pi;
    end
    else begin
      sfd_pcfv := pcfv_gen(a,x);
    end;
    exit;
  end;
  {Here a >= 0}
  if x>=0.0 then sfd_pcfv := pcfv_axpos(a,x)
  else begin
    {x<0, use reflection formulas}
    if frac(a)=0.0 then begin
      a := a+0.5;
      {$ifdef VER50}
        t := sfd_pcfd(-a, -x);
        t := t/Pi;
      {$else}
        t := sfd_pcfd(-a, -x)/Pi;
      {$endif}
      if t=0.0 then sfd_pcfv := 0.0
      else if a<MAXGAMD then sfd_pcfv := sfd_gamma(a)*t
      else begin
        {Use ln(Gamma(a)) because Gamma(a) overflows. As usual, }
        {this has much larger errors but is better than nothing.}
        g := sfd_lngamma(a);
        g := exp(g+ln(abs(t)));
        if t<0.0 then g := -g;
        sfd_pcfv := g;
      end;
    end
    else begin
      t := pcfv_axpos(a,-x);
      if odd(trunc(a)) then t := -t;
      sfd_pcfv := t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_pcfhh(v,x: double): double;
  {-Return the Hermite function H_v(x) of degree v}
var
  g,t,u,z: double;
begin
  if IsNanOrInfD(v) or IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_pcfhh := NaN_d;
    exit;
  end;
  {See Lebedev[45], Sec. 10.2}
  t := -0.5*v;
  z := x*x;
  u := sfd_chu(t, 0.5, z);
  if x < 0.0 then begin
    if (frac(v)=0.0) and (v>0.0) then begin
      {Case v positive integer H_v(-x) = (-1)^v * H_v(x)}
      {http://functions.wolfram.com/07.01.16.0001.01}
      if frac(t) <> 0.0 then u := -u;
    end
    else begin
      {http://functions.wolfram.com/07.01.16.0012.01}
      g := sfd_rgamma(t);
      if g<>0.0 then begin
        t := sfd_1f1(0.5+t, 1.5, z);
        u := u - 4.0*x*g*t*SqrtPi;
      end;
    end;
  end;
  sfd_pcfhh := exp2(v)*u;
end;

end.
