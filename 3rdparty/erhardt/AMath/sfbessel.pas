unit sfBessel;

{Common code for special functions: Bessel functions and related}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Common code: Bessel functions and related

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
                 [13] [NR]: W.H. Press et al, Numerical Recipes in C, 2nd ed., Cambridge, 1992,
                      http://www.nrbook.com/a/bookcpdf.html
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
                 [51] N.M. Temme, On the Numerical Evaluation of the Ordinary Bessel Function
                      of the Second Kind. J. Comput. Phys., 21(3): 343-350 (1976), section 2.
                      Available as http://oai.cwi.nl/oai/asset/10710/10710D.pdf
                 [52] N.M. Temme, On the Numerical Evaluation of the Modified Bessel
                      Function of the Third Kind, 2nd edition. Preprint, available
                      as http://oai.cwi.nl/oai/asset/7885/7885A.pdf

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  17.08.10  W.Ehrhardt  Common version number after split
 1.00.01  11.09.10  we          Improved arg checking and NAN/INF handling
 1.00.02  12.09.10  we          Extended over/underflow check in sfc_jn

 1.02.00  31.10.10  we          sfc_yn: improved pre-checks
 1.02.01  01.11.10  we          sfc_jv, sfc_yv, sfc_bess_jyv
 1.02.02  03.11.10  we          sfc_iv, sfc_kv, sfc_bess_ikv
 1.02.03  04.11.10  we          NAN/INF handling for real order functions
 1.02.04  05.11.10  we          Airy functions Ai, Ai', Bi, Bi'
 1.02.05  06.11.10  we          sfc_sph_jn, sfc_sph_yn
 1.02.06  06.11.10  we          increased MAXIT in CF1_j
 1.02.07  14.11.10  we          Fixed sfc_i0e for small x
 1.02.08  15.11.10  we          sfc_ive, sfc_kve

 1.04.00  10.02.11  we          Improved J0,J1 for |x| >= 500; Y0,Y1 for x >= 1600

 1.05.00  10.04.11  we          Zero order Kelvin functions ber,bei,ker,kei
 1.05.01  10.04.11  we          sfc_kei(0) = -Pi/4

 1.06.00  05.05.11  we          sfc_struve_h0
 1.06.01  05.05.11  we          sfc_struve_h1
 1.06.02  06.05.11  we          sfc_struve_l1
 1.06.03  07.05.11  we          sfc_struve_l0
 1.06.04  26.05.11  we          fix for very large x in kelvin_large

 1.08.00  13.08.11  we          special cases |v|=0,1 in sfc_kve
 1.08.01  13.08.11  we          v=-1 in sfc_iv/e

 1.10.00  28.12.11  we          sfc_sph_in, sfc_sph_ine
 1.10.01  29.12.11  we          sfc_sph_kn, sfc_sph_kne
 1.10.02  29.12.11  we          Fix bug for small negative x in sfc_sph_ine

 1.15.00  19.02.13  we          improved quick check in sfc_jn
 1.15.01  20.02.13  we          Yv_series
 1.15.02  20.02.13  we          handle some near overflows in bessel_jy
 1.15.03  20.02.13  we          improved check in sfc_yn

 1.18.00  09.05.13  we          Airy/Scorer functions sfc_airy_gi/hi
 1.18.01  09.05.13  we          Prevent some wrong compiler optimizations for div by 3

 1.19.00  06.06.13  we          Fix typo in h1v_large
 1.19.01  08.06.13  we          sfc_bess_kv2
 1.19.02  27.06.13  we          Check overflow / return INF in bessel_jy

 1.24.00  26.03.14  we          sfc_yn with LnPi from AMath

 1.28.00  08.08.14  we          Iv := PosInf_x if Kv,Kv1=0 in bessel_ik or if x >= IvMaxXH in sfc_iv

 1.33.00  07.06.15  we          new IJ_series replaces Jv_series
 1.33.01  07.06.15  we          rewrite of sfc_in: use IJ_series and CF1_I
 1.33.02  08.06.15  we          improved bess_m0p0/bess_m1p1 with rem_2pi_sym and Kahan summation
 1.33.03  09.06.15  we          sfc_struve_l
 1.33.04  10.06.15  we          sfc_struve_h
 1.33.05  10.06.15  we          sfc_struve_l0/1 use sfc_struve_l
 1.33.06  14.06.15  we          sfc_struve_h/l with real v >= 0
 1.33.07  19.06.15  we          scaled Airy functions sfc_airy_ais, sfc_airy_bis

 1.34.00  05.07.15  we          special cases Yv(+- 1/2, x)
 1.34.01  06.07.15  we          sfc_struve_h for v < 0
 1.34.02  07.07.15  we          sfc_struve_l for v < 0
 1.34.03  09.07.15  we          fix sfc_struve_h/l(v,0) for v<=-1

 1.35.00  26.08.15  we          sfBessel: Sqrt2 in kelvin functions (anti-'optimizing')

 1.40.00  29.06.17  we          Code removed for TP5-TP6, TPW1-D1

 1.41.00  16.07.17  we          Bessel lambda sfc_blam

 1.43.00  02.12.17  we          Suppress warnings: Local variable does not seem to be initialized0

 1.44.00  06.02.18  we          sfc_struve_h with sfc_intdei_p

 1.46.00  21.03.18  we          sfc_j0int
 1.46.01  22.03.18  we          sfc_y0int
 1.46.02  23.03.18  we          sfc_k0int
 1.46.03  24.03.18  we          sfc_i0int
 1.46.04  27.03.18  we          improved sfc_[i,j,y]0int

 1.47.00  19.04.18  we          Fixes for FPC311

 1.50.00  23.08.18  we          NanOrInf check in sfc_jn, sfc_in

 1.51.00  36.09.18  we          Moved Airy/Kelvin/Struve functions to new unit sfBessl2

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


function sfc_i0(x: extended): extended;
  {-Return I0(x), the modified Bessel function of the 1st kind, order zero}

function sfc_i0e(x: extended): extended;
  {-Return I0(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order zero}

function sfc_i1(x: extended): extended;
  {-Return I1(x), the modified Bessel function of the 1st kind, order one}

function sfc_i1e(x: extended): extended;
  {-Return I1(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order one}

function sfc_in(n: integer; x: extended): extended;
  {-Return I_n(x), the modified Bessel function of the 1st kind, order n.}

function sfc_j0(x: extended): extended;
  {-Return J0(x), the Bessel function of the 1st kind, order zero}

function sfc_j1(x: extended): extended;
  {-Return J1(x), the Bessel function of the 1st kind, order one}

function sfc_jn(n: integer; x: extended): extended;
  {-Return J_n(x), the Bessel function of the 1st kind, order n; not suitable for large n or x.}

function sfc_k0(x: extended): extended;
  {-Return K0(x), the modified Bessel function of the 2nd kind, order zero, x>0}

function sfc_k0e(x: extended): extended;
  {-Return K0(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order zero, x>0}

function sfc_k1(x: extended): extended;
  {-Return K1(x), the modified Bessel function of the 2nd kind, order one, x>0}

function sfc_k1e(x: extended): extended;
  {-Return K1(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order one, x>0}

function sfc_kn(n: integer; x: extended): extended;
  {-Return K_n(x), the modified Bessel function of the 2nd kind, order n, x>0, not suitable for large n}

function sfc_y0(x: extended): extended;
  {-Return Y0(x), the Bessel function of the 2nd kind, order zero; x>0}

function sfc_y1(x: extended): extended;
  {-Return Y1(x), the Bessel function of the 2nd kind, order one; x>0}

function sfc_yn(n: integer; x: extended): extended;
  {-Return Y_n(x), the Bessel function of the 2nd kind, order n, x>0, not suitable for large n or x}


function sfc_jv(v, x: extended): extended;
  {-Return J_v(x), the Bessel function of the 1st kind, order v; not suitable for large v.}

function sfc_yv(v, x: extended): extended;
  {-Return Y_v(x), the Bessel function of the 2nd kind, order v; x > 0; not suitable for large v.}

function sfc_iv(v, x: extended): extended;
  {-Return I_v(x), the modified Bessel function of the 1st kind, order v.}

function sfc_kv(v, x: extended): extended;
  {-Return K_v(x), the modified Bessel function of the 2nd kind, order v, x>0}

function sfc_ive(v, x: extended): extended;
  {-Return I_v(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order v.}

function sfc_kve(v, x: extended): extended;
  {-Return K_v(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order v, x>0}

procedure sfc_bess_ikv(v,x: extended; var Iv,Kv: extended);
  {-Return I_v(x) and K_v(x), no checks, x>0, |v| < MaxLongint}

procedure sfc_bess_jyv(v,x: extended; var Jv,Yv: extended);
  {-Return J_v(x) and Y_v(x), no checks, x>0, |v| < MaxLongint}

procedure sfc_bess_kv2(v,x: extended; var Kv, Kv1: extended);
  {-Return K(v,x) and K(v+1,x), no checks, x>0, |v| < MaxLongint}


function sfc_sph_jn(n: integer; x: extended): extended;
  {-Return j_n(x), the spherical Bessel function of the 1st kind, order n}

function sfc_sph_yn(n: integer; x: extended): extended;
  {-Return y_n(x), the spherical Bessel function of the 2nd kind, order n >=0 , x<>0}

function sfc_sph_in(n: integer; x: extended): extended;
  {-Return i_n(x), the modified spherical Bessel function of the 1st/2nd kind, order n}

function sfc_sph_ine(n: integer; x: extended): extended;
  {-Return i_n(x)*exp(-|x|), the exponentially scaled modified spherical Bessel function of the 1st/2nd kind, order n}

function sfc_sph_kn(n: integer; x: extended): extended;
  {-Return k_n(x), the modified spherical Bessel function of the 3rd kind, order n, x>0}

function sfc_sph_kne(n: integer; x: extended): extended;
  {-Return k_n(x)*exp(x), the exponentially scaled modified spherical Bessel function of the 3rd kind, order n, x>0}

function sfc_blam(v,x: extended): extended;
  {-Compute lambda(v,x) = Gamma(v+1)*J(v,x)/(0.5x)^v for v,x >= 0}


function sfc_i0int(u: extended): extended;
  {-Return the integral int(bessel_i0(x), x = 0..u)}

function sfc_j0int(u: extended): extended;
  {-Return the integral int(bessel_j0(x), x = 0..u)}

function sfc_y0int(u: extended): extended;
  {-Return the integral int(bessel_y0(x), x = 0..u), u >= 0}

function sfc_k0int(u: extended): extended;
  {-Return the integral int(bessel_k0(x), x = 0..u), u >= 0}


procedure bessel_ik(v,x: extended; CalcI, escale: boolean; var Iv,Kv: extended);
  {-Return I_v(x) and/or K_v(x) depending on CalcI, x>0, |v| < MaxLongint}
  { If escale=true the values are exponentially scaled.}

implementation


uses
  AMath,
  sfBasic,  {Basic common code}
  sfGamma;  {Gamma function and related}


{---------------------------------------------------------------------------}
function IJ_series(v,x: extended; CalcJv: boolean): extended;
  {-Power series for Bessel J_v(x) or I_v(x), 0 <= v < MAXGAMX-1, 0 <= x <= 2}
var
  f,s,t: extended;
  n: integer;
begin
  f := 0.5*x;
  t := power(f,v)/sfc_gamma(v+1.0);
  if CalcJv then f := -f*f else f := f*f;
  s := t;
  n := 0;
  repeat
    inc(n);
    t := t*f/n/(v+n);
    s := s + t;
  until abs(t) <= 0.5*eps_x*abs(s);
  IJ_series := s;
end;


{---------------------------------------------------------------------------}
procedure CF1_j(v,x: extended; var fv: extended; var s: integer);
  {-Return J_(v+1)(x) / J_v(x), efficient only if |x| <= |v|}
var
  c,d,f,b,t,tiny,tol: extended;
  k: longint;

const
  MAXIT = longint(32000)*100;  {see note below}

begin
  s := 1;

  {Evaluate HMF [1], 9.1.73 using modified Lentz's method. s keeps track }
  {of sign changes in the denominator. Ref: NR [13] (6.7.2) and p. 244,  }
  {function bessjy and Boost [19] file bessel_jy.hpp, function CF1_jy.   }

  {Note that CF1_j needs about O(|x|) iterations if |x| > |v|. But unless}
  {there is a better implementation below the asymptotic range,  CF1_j is}
  {is used in the last resort branch of bessel_jy. If a better algorithm }
  {(like Olver/Temme uniform Airy type asymptotic expansion) is available}
  {the factor 100 in the MAXIT declaration should be removed.}

  tol  := 2.0*eps_x;
  tiny := Sqrt_MinExt;
  c := tiny;
  f := tiny;
  d := 0.0;
  for k:=1 to MAXIT do begin
    b := 2.0*(v + k)/x;
    c := b - 1.0/c;
    d := b - d;
    if c=0.0 then c := tiny;
    if d=0.0 then d := tiny;
    d := 1.0/d;
    t := c * d;
    f := f*t;
    if d<0 then s := -s;
    if abs(t-1.0) < tol then begin
      fv := -f;
      exit;
    end;
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  fv := -f;
end;


{---------------------------------------------------------------------------}
procedure CF1_I(v,x: extended; var fv: extended);
  {-Return I_(v+1)(x) / I_v(x) using continued fraction}
var
  c,d,f,b,t,tiny,tol: extended;
  k: integer;
const
  MAXIT = 30000;
begin
  {Evaluate NIST[30], 10.33.1 using modified Lentz's method.}
  {Ref: NR [13] (6.7.21) and p.248, function bessik         }
  {and Boost [19] file bessel_ik.hpp, function CF1_Ik       }
  {If |x| <= |v|, CF1_I converges rapidly but if |x| > |v|  }
  {then CF1_I needs O(|x|) iterations to converge!          }

  tol  := 2.0*eps_x;
  tiny := Sqrt_MinExt;
  c := tiny;
  f := tiny;
  d := 0.0;
  for k:=1 to MAXIT do begin
    b := 2.0*(v + k)/x;
    c := b + 1.0/c;
    d := b + d;
    if c=0.0 then c := tiny;
    if d=0.0 then d := tiny;
    d := 1.0/d;
    t := c * d;
    f := f*t;
    if abs(t-1.0) < tol then begin
      fv := f;
      exit;
    end;
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  fv := f;
end;


{---------------------------------------------------------------------------}
procedure h1v_large(v, x: extended; var mv,tmx: extended);
  {-Return modulus and (phase - x) of the Hankel function H1_v(x), x > 0 large}
var
  s,m,m2,y: extended;
const
  c56 = 5.0/6.0;
begin
  {Modulus Mv: asymptotic expansion from HMF[1] 9.2.28}
  y  := sqr(0.5/x);
  m  := 4.0*sqr(v);
  m2 := sqr(m);
  s  := 1.0 + 0.5*y*(m-1.0)*(1.0 + 0.75*y*(m-9.0)*(1.0 + c56*y*(m-25.0)));
  mv := sqrt(2.0*s/(Pi*x));
  {Phase theta_v - x: asymptotic expansion from HMF[1] 9.2.29}
  y  := 0.25*y;
  s  := (5.0*m*m2 - 1535.0*m2 + 54703.0*m - 375733.0)/14.0;
  s  := s*y + (m2 - 114.0*m + 1073)/5.0;
  s  := s*y + (m-25.0)/6.0;
  tmx:= (m-1.0)*(s*y + 0.5)/(4.0*x) - Pi*(0.5*v+0.25)
end;


{---------------------------------------------------------------------------}
function bessj_large(v, x: extended): extended;
  {-Return J_v(x) via modulus/phase asymptotic expansion, x large}
var
  mv,tv,st,ct,sx,cx: extended;
begin
  h1v_large(v,x,mv,tv);
  sincos(tv,st,ct);
  sincos(x,sx,cx);
  {J_v := mv*cos(x+tv); cos(x+tv) = cos(x)cos(tv) - sin(x)sin(tv)}
  bessj_large := mv*(cx*ct - sx*st);
end;


{---------------------------------------------------------------------------}
function bessy_large(v, x: extended): extended;
  {-Return Y_v(x) via modulus/phase asymptotic expansion, x large}
var
  mv,tv,st,ct,sx,cx: extended;
begin
  h1v_large(v,x,mv,tv);
  sincos(tv,st,ct);
  sincos(x,sx,cx);
  {Y_v := mv*sin(x+tv); sin(x+tv) = cos(x)sin(tv) + sin(x)cos(tv)}
  bessy_large := mv*(st*cx + ct*sx);
end;


{---------------------------------------------------------------------------}
procedure bess_m0p0(x: extended; var m0,p0: extended);
  {-Modulus and phase for J0(x) and Y0(x), x >= 9.0}
var
  y,z,s: extended;
const
  m0nhex: array[0..7] of THexExtW = (
             ($cb2b,$4b73,$8075,$8aa3,$3ff8),
             ($6b78,$4cc6,$25b7,$b912,$3ffb),
             ($cfe9,$74e0,$67a1,$b75e,$3ffe),
             ($b1cd,$4e5e,$2274,$b5ec,$4000),
             ($5e4b,$e3af,$59bb,$f409,$4001),
             ($b343,$2673,$4e51,$a36b,$4002),
             ($38a3,$a663,$7b91,$dbab,$4001),
             ($8559,$f552,$3a38,$ca1d,$3ffd));
  m0dhex: array[0..7] of THexExtW = (
             ($8a83,$1b80,$003e,$adc2,$3ff8),
             ($ed5a,$31cd,$b3ac,$e7f3,$3ffb),
             ($7e3f,$b8dd,$04df,$e5fd,$3ffe),
             ($775a,$1b79,$7d9c,$e475,$4000),
             ($b837,$3075,$dbc0,$99ce,$4002),
             ($5e3d,$b5f4,$9848,$d032,$4002),
             ($4498,$3d2a,$f3fb,$91df,$4002),
             ($0000,$0000,$0000,$8000,$3fff));
  p0nhex: array[0..5] of THexExtW = (
             ($4c2f,$2dd8,$79c3,$e65d,$bfe9),
             ($dc17,$325e,$8baf,$9d35,$bff1),
             ($e514,$8866,$25a9,$8309,$bff7),
             ($8d8a,$84e7,$dbd5,$9e75,$bffb),
             ($1e30,$04da,$b769,$800a,$bffe),
             ($5106,$12a6,$4dd2,$bc55,$bffe));
  p0dhex: array[0..6] of THexExtW = (
             ($4c8c,$2dd8,$79c3,$e65d,$3fec),
             ($ac5c,$4806,$8709,$9dad,$3ff4),
             ($37bf,$fcc8,$9b9f,$844b,$3ffa),
             ($6f25,$2a95,$2dc6,$a285,$3ffe),
             ($4b69,$3f87,$131f,$891f,$4001),
             ($f3e9,$b2a5,$6652,$ec17,$4001),
             ($0000,$0000,$0000,$8000,$3fff));
var
  m0n: array[0..7] of extended absolute m0nhex;
  m0d: array[0..7] of extended absolute m0dhex;
  p0n: array[0..5] of extended absolute p0nhex;
  p0d: array[0..6] of extended absolute p0dhex;
begin
  {Ref: Cephes [7], file ldouble\j0l.c}
  {See also HMF[1], sections 9.2.17 .. 9.2.31}

  {Calculate the modulus m0(x) = sqrt(J0(x)^2 + Y0(x)^2) and the}
  {phase p0(x) = arctan(Y0(x)/J0(x)) with rational approximations}
  {For x>=9: J0(x) = m0(x)*cos(p0(x)) and Y0(x) = m0(x)*sin(p0(x))}
  z  := sqr(1.0/x);
  y  := abs(x);
  s  := rem_2pi_sym(y);
  p0 := PolEvalX(z,p0n,6)/PolEvalX(z,p0d,7)/y;
  z  := 1.0/y;
  m0 := PolEvalX(z,m0n,8)/PolEvalX(z,m0d,8)/sqrt(y);

  {Compute p0 := rem_2pi_sym(y) - Pi_4 + p0 with optimized Kahan}
  {summation. Without rem_p2pi and Kahan we get a relative error}
  {of 1.4e-18 for J0(18), with this code it is 5.4e-20.}

  {Note that this improves only J0/Y0 for arguments x near the  }
  {zeros. With the recurrence relations for Jn/Yn only absolute }
  {(NOT relative) accuracies of order eps_x can be achieved!    }
  z  := -Pi_4;
  y  := s + z;
  z  := (y - s) - z;
  p0 := y + (p0 - z);
end;


{---------------------------------------------------------------------------}
function sfc_j0(x: extended): extended;
  {-Return J0(x), the Bessel function of the 1st kind, order zero}
var
  y,z: extended;
const
  {Squares of first three roots of J0, calculated with Maple and t_rcalc/xh}
  j1h: THexExtW = ($DBB7,$315D,$DC02,$B90F,$4001);  {5.7831859629467845213}
  j2h: THexExtW = ($7778,$0EEB,$2531,$F3C5,$4003);  {3.0471262343662086400E+1}
  j3h: THexExtW = ($40F6,$0ABB,$25C1,$95C6,$4005);  {7.4887006790695183442E+1}
var
  jz1: extended absolute j1h;
  jz2: extended absolute j2h;
  jz3: extended absolute j3h;
const
  j0nhex:  array[0..7] of THexExtW = (
             ($b96c,$c486,$fb95,$9f47,$c03a),
             ($4018,$ad26,$71ba,$e643,$4034),
             ($1b0b,$6331,$7add,$8753,$c02e),
             ($943a,$69b7,$36ca,$a996,$4026),
             ($008c,$7b60,$d119,$f792,$c01d),
             ($fe10,$b608,$4829,$d503,$4014),
             ($a9a8,$e62b,$3b28,$ca73,$c00a),
             ($f759,$4208,$23d6,$a5ff,$3fff));
  j0dhex: array[0..8] of THexExtW = (
             ($00ac,$fb2b,$6f62,$804b,$4048),
             ($fdce,$a4ca,$2ed8,$88b8,$4041),
             ($3d2c,$ed55,$20e1,$9105,$4039),
             ($0841,$8cb6,$5a46,$c9e3,$4030),
             ($fed1,$086d,$3425,$cc0a,$4027),
             ($66d2,$93fe,$0762,$9b79,$401e),
             ($e1a0,$923f,$cb5c,$b1a2,$4014),
             ($bdfe,$c832,$5b9f,$8e9f,$400a),
             ($0000,$0000,$0000,$8000,$3fff));
var
  j0n: array[0..7] of extended absolute j0nhex;
  j0d: array[0..8] of extended absolute j0dhex;
begin
  {Ref: Cephes [7], file ldouble\j0l.c}
  x := abs(x);
  if x < 9.0 then begin
    {In the interval [0,9) a rational approximation of the form }
    {J0(x) = (x^2 - r^2) (x^2 - s^2) (x^2 - t^2) P7(x^2)/Q8(x^2)}
    {is used, where r, s, t are the first three zeros of J0.}
    z := sqr(x);
    y := (z - jz1)*(z - jz2)*(z - jz3);
    sfc_j0 := y*PolEvalX(z,j0n,8)/PolEvalX(z,j0d,9);
  end
  else begin
    {For x>=9 the common rational approximations to modulus}
    {and phase are used J0(x) = modulus * cos(phase).}
    if x >= 500.0 then sfc_j0 := bessj_large(0,x)
    else begin
      bess_m0p0(x,y,z);
      sfc_j0 := y*cos(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_y0(x: extended): extended;
  {-Return Y0(x), the Bessel function of the 2nd kind, order zero; x>0}
var
  y, z: extended;
const
  {The first four roots of Y0, calculated with Maple and t_rcalc/xh}
  y1h: THexExtW = ($DE10,$C2B4,$9A6C,$FD4A,$4000);  {3.9576784193148578685}
  y2h: THexExtW = ($1734,$3908,$EE27,$E2C0,$4001);  {7.0860510603017726975}
  y3h: THexExtW = ($BC69,$23EA,$B9AD,$A38E,$4002);  {10.222345043496417019}
  y4h: THexExtW = ($4D04,$0F3A,$0E25,$D5C7,$4002);  {13.361097473872763478}
var
  y0z1: extended absolute y1h;
  y0z2: extended absolute y2h;
  y0z3: extended absolute y3h;
  y0z4: extended absolute y4h;
const
  y0nhex: array[0..7] of THexExtW = (
            ($5fbd,$0171,$135a,$8340,$c035),
            ($501f,$6264,$bdf4,$9d17,$4036),
            ($23c9,$6b29,$4244,$c4c9,$c032),
            ($b219,$37ba,$5142,$9f1f,$402d),
            ($3e3c,$b343,$46c9,$e45f,$c026),
            ($2fdd,$4b27,$ca98,$a1c3,$401f),
            ($2ec0,$7b95,$297f,$df70,$c016),
            ($126c,$20be,$647f,$f344,$400c));
  y0dhex: array[0..7] of THexExtW = (
            ($241a,$8f2b,$629a,$de4b,$4038),
            ($04d3,$a629,$d61d,$b410,$4032),
            ($6732,$8c1b,$c5ab,$9384,$402b),
            ($553b,$4dc8,$8695,$a0c3,$4023),
            ($97a4,$90fa,$a7e9,$801c,$401b),
            ($d938,$b6b2,$71d8,$98be,$4012),
            ($9057,$7f25,$59b7,$8219,$4009),
            ($0000,$0000,$0000,$8000,$3fff));
  y059nh: array[0..9] of THexExtW = (
            ($f90c,$3510,$0be9,$87e7,$c012),
            ($fd54,$b2fe,$0a23,$e37e,$c00f),
            ($8c07,$29e3,$11be,$9f10,$4012),
            ($49e2,$fb52,$02af,$be8a,$4010),
            ($e8fa,$4b44,$4a39,$dc5b,$400b),
            ($62e0,$c25b,$2cb3,$8f12,$c00b),
            ($d5a3,$f673,$4e59,$9a8c,$4005),
            ($5504,$035a,$59fa,$ca14,$4003),
            ($1207,$46ea,$c3db,$bc88,$bfff),
            ($992f,$ab45,$90b6,$c20b,$3ff9));
  y059dh: array[0..9] of THexExtW = (
            ($3b3b,$ea0b,$b8d1,$8bd7,$401d),
            ($ceb6,$3463,$5ddb,$d1b5,$401e),
            ($e26b,$76b9,$250a,$a7fb,$c01c),
            ($27ff,$ca92,$3d78,$cea1,$4019),
            ($ec8a,$4697,$ddde,$a742,$c016),
            ($e8b6,$d705,$da91,$d62c,$4012),
            ($a28c,$5563,$d19f,$c75e,$c00e),
            ($ad09,$8e6a,$a502,$8b0c,$400a),
            ($debf,$a468,$8a55,$f96b,$c004),
            ($0000,$0000,$0000,$8000,$3fff));
var
  y059n: array[0..9] of extended absolute y059nh;
  y059d: array[0..9] of extended absolute y059dh;
  y0n:   array[0..7] of extended absolute y0nhex;
  y0d:   array[0..7] of extended absolute y0dhex;
begin
  {Ref: Cephes [7], file ldouble\j0l.c}
  if x<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_y0 := NaN_x;
    exit;
  end;
  if x < 9.0 then begin
    z := sqr(x);
    if z < 20.25 then begin
      {In the interval [0,4.5) a rational approximation of the}
      {form Y0(x) = P7(x)/Q7(x) + 2/Pi*ln(x)*J0(x) is used.   }
      y := ln(x)*sfc_j0(x)/Pi_2;
      sfc_y0 := y + PolEvalX(z,y0n,8)/PolEvalX(z,y0d,8);
    end
    else begin
      {In the interval [4.5,9) a rational approximation of the}
      {form Y0(x) = (x - p)(x - q)(x - r)(x - s)P9(x)/Q9(x) is}
      {is used where p, q, r, s are first four zeros of Y0(x).}
      y := (x - y0z1)*(x - y0z2)*(x - y0z3)*(x - y0z4);
      sfc_y0 := y * PolEvalX(x,y059n,10)/PolEvalX(x,y059d,10);
    end;
  end
  else begin
    {For x>=9 the common rational approximations to modulus}
    {and phase are used Y0(x) = modulus * sin(phase).}
    if x >= 1600 then sfc_y0 := bessy_large(0,x)
    else begin
      bess_m0p0(x,y,z);
      sfc_y0 := y*sin(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure bess_m1p1(x: extended; var m1,p1: extended);
  {-Modulus and phase for J1(x) and Y1(x), x >= 9.0}
var
  y,z,s: extended;
const
  m1nhex: array[0..7] of THexExtW = (
            ($a905,$05fb,$3101,$82c9,$3ff9),
            ($6de4,$8fae,$fe26,$8097,$3ffc),
            ($2cb0,$c657,$be70,$81e0,$3fff),
            ($71e6,$88a5,$0a53,$b702,$4000),
            ($e5e2,$6914,$3a08,$e582,$4001),
            ($7d55,$db8c,$e825,$a1c2,$4000),
            ($3111,$863a,$3a61,$c8a0,$3ffd),
            ($3d53,$b598,$f3bf,$a155,$c001));
  m1dhex: array[0..8] of THexExtW = (
            ($1237,$cc6c,$7356,$a3ea,$3ff9),
            ($fc82,$02c7,$17a4,$a12b,$3ffc),
            ($37ce,$79ae,$2f15,$a24c,$3fff),
            ($77b6,$34e2,$501a,$e37a,$4000),
            ($0260,$746b,$d030,$8c14,$4002),
            ($6420,$97ce,$8e44,$a208,$4000),
            ($77b5,$8f2d,$b6bf,$ebe1,$bffe),
            ($2603,$640e,$7d8d,$c775,$c001),
            ($0000,$0000,$0000,$8000,$3fff));
  p1nhex: array[0..5] of THexExtW = (
            ($540c,$c1d5,$b096,$e54f,$3feb),
            ($f74f,$be87,$7e7d,$9741,$3ff3),
            ($a830,$f4a3,$2c60,$f144,$3ff8),
            ($e907,$28b9,$7cb7,$895c,$3ffd),
            ($6050,$98aa,$3500,$cb2f,$3fff),
            ($ebc0,$5506,$512f,$80ab,$4000));
  p1dhex: array[0..6] of THexExtW = (
            ($8d72,$2be3,$cb0f,$98df,$3fed),
            ($a853,$55fb,$6c79,$ca32,$3ff4),
            ($98f8,$d610,$3c35,$a235,$3ffa),
            ($c39e,$9c8c,$5428,$bb65,$3ffe),
            ($b1f2,$e0d2,$5ab5,$9098,$4001),
            ($efe3,$292c,$0d43,$d9e6,$4001),
            ($0000,$0000,$0000,$8000,$3fff));
var
  m1n: array[0..7] of extended absolute m1nhex;
  m1d: array[0..8] of extended absolute m1dhex;
  p1n: array[0..5] of extended absolute p1nhex;
  p1d: array[0..6] of extended absolute p1dhex;
begin
  {Ref: Cephes [7], file ldouble\j1l.c}

  {Calculate the modulus m1(x) = sign(x)*sqrt(J1(x)^2 + Y1(x)^2) and }
  {the phase p1(x) = arctan(Y1(x)/J1(x)) with rational approximations}
  {For x>=9: J1(x) = m1(x)*cos(p1(x)) and Y1(x) = m1(x)*sin(p1(x))}
  z  := sqr(1.0/x);
  y  := abs(x);
  s  := rem_2pi_sym(y);
  p1 := PolEvalX(z,p1n,6)/PolEvalX(z,p1d,7)/y;
  z  := 1.0/y;
  m1 := PolEvalX(z,m1n,8)/PolEvalX(z,m1d,9)/sqrt(y);
  if x<0.0 then m1 := -m1;

  {Compute p1 := rem_2pi_sym(y) - 3Pi_4 + p1 with optimized Kahan summation}
  z  := -3.0*Pi_4;
  y  := s + z;
  z  := (y - s) - z;
  p1 := y + (p1 - z);
end;


{---------------------------------------------------------------------------}
function sfc_j1(x: extended): extended;
  {-Return J1(x), the Bessel function of the 1st kind, order one}
var
  y,z: extended;
const
  {Squares of first three roots of J1, calculated with Maple and t_rcalc/xh}
  j1h: THexExtW = ($5F8E,$4C11,$5A0C,$EAE9,$4002);  {1.4681970642123893257E+1}
  j2h: THexExtW = ($9093,$9521,$B303,$C4DF,$4004);  {4.9218456321694603672E+1}
  j3h: THexExtW = ($5EBF,$C2F1,$B86B,$CEFF,$4005);  {1.0349945389513658033E+2}
var
  jz1: extended absolute j1h;
  jz2: extended absolute j2h;
  jz3: extended absolute j3h;
const
  j1nhex:  array[0..8] of THexExtW = (
             ($d8d8,$7311,$a7d2,$97a4,$c039),
             ($d3c2,$f8f0,$f852,$c144,$4033),
             ($636c,$4d29,$9f71,$cebb,$c02c),
             ($038e,$bd23,$a7fa,$f49c,$4024),
             ($1ac8,$c825,$3c9c,$b0b6,$c01c),
             ($38f5,$f72b,$0a5c,$a122,$4013),
             ($29f3,$496b,$a54c,$b6d9,$c009),
             ($6dc3,$c850,$a096,$ee6b,$3ffe),
             ($f72f,$18cc,$50b2,$8a22,$bff3));
  j1dhex: array[0..8] of THexExtW = (
             ($dd67,$f5b3,$0522,$ad0f,$404a),
             ($665d,$b178,$242e,$9af7,$4043),
             ($e6c0,$a725,$3d56,$88f7,$403b),
             ($0122,$56c0,$f2ef,$9d6e,$4032),
             ($b498,$fdd5,$209e,$820e,$4029),
             ($6041,$c9fe,$6890,$a033,$401f),
             ($6a17,$e162,$4e86,$9218,$4015),
             ($baf9,$146e,$df50,$b88a,$400a),
             ($0000,$0000,$0000,$8000,$3fff));
var
  j1n: array[0..8] of extended absolute j1nhex;
  j1d: array[0..8] of extended absolute j1dhex;
begin
  {Ref: Cephes [7], file ldouble\j1l.c}
  z := abs(x);
  if z < 9.0 then begin
    z := sqr(x);
    {In the interval [0,9) a rational approximation of the form }
    {J1(x) = x*(x^2 - r^2)*(x^2 - s^2)*(x^2 - t^2)*P8(x^2)/Q8(x^2)}
    {is used, where r, s, t are the first three zeros of J1.}
    y := x*(z - jz1)*(z - jz2)*(z - jz3);
    sfc_j1 := y*PolEvalX(z,j1n,9)/PolEvalX(z,j1d,9);
  end
  else begin
    {For x>=9 the common rational approximations to modulus}
    {and phase are used J1(x) = modulus * cos(phase).}
    if z >= 500.0 then begin
      y := bessj_large(1,z);
      if x<0.0 then sfc_j1 := -y else sfc_j1 := y;
    end
    else begin
      bess_m1p1(x,y,z);
      sfc_j1 := y*cos(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_y1(x: extended): extended;
  {-Return Y1(x), the Bessel function of the 2nd kind, order one; x>0}
var
  y, z: extended;
const
  {The first four roots of Y1, calculated with Maple and t_rcalc/xh}
  y1h: THexExtW = ($1721,$FF92,$F6A6,$8C9D,$4000); {2.1971413260310170351}
  y2h: THexExtW = ($73C0,$3D81,$F274,$ADBF,$4001); {5.4296810407941351328}
  y3h: THexExtW = ($A148,$0B4D,$3D73,$8989,$4002); {8.5960058683311689268}
  y4h: THexExtW = ($3022,$A190,$89C6,$BBFC,$4002); {11.749154830839881243}
var
  y1z1: extended absolute y1h;
  y1z2: extended absolute y2h;
  y1z3: extended absolute y3h;
  y1z4: extended absolute y4h;
const
  y1nhex: array[0..6] of THexExtW = (
            ($3a10,$0848,$5930,$9965,$c035),
            ($7f8b,$4757,$75bd,$a196,$4033),
            ($69fd,$1242,$f62d,$de75,$c02e),
            ($5633,$aa6b,$79e5,$e62c,$4028),
            ($7607,$a687,$af0a,$d892,$c021),
            ($53e4,$194c,$befa,$bd19,$4019),
            ($5b16,$f7f8,$0d7e,$fbbd,$c00f));
  y1dhex: array[0..7] of THexExtW = (
            ($7302,$b91b,$de7e,$c399,$4037),
            ($8c6a,$397e,$0963,$ad7a,$4031),
            ($aaf0,$342b,$d098,$9ca5,$402a),
            ($57e0,$1d92,$90a9,$bd99,$4022),
            ($0e86,$117b,$36d6,$a94a,$401a),
            ($298c,$29ef,$0630,$e482,$4011),
            ($dd1a,$3b8e,$ab73,$df28,$4008),
            ($0000,$0000,$0000,$8000,$3fff));
  y159nh: array[0..9] of THexExtW = (
            ($539b,$f305,$c3d8,$97f6,$4011),
            ($f62f,$d968,$8c66,$8d15,$c013),
            ($3811,$a3da,$413f,$dc24,$c013),
            ($cd43,$2f50,$1118,$d972,$c013),
            ($a33b,$8229,$1561,$f1fc,$c00f),
            ($b2bf,$4296,$65af,$a3d1,$400d),
            ($df40,$226b,$7e37,$c0d9,$400b),
            ($e917,$8486,$0ebd,$e6c3,$c008),
            ($fdf1,$41e5,$4beb,$ac44,$4004),
            ($b5e5,$bb42,$f667,$ae3f,$bffe));
  y159dh: array[0..10] of THexExtW = (
            ($a231,$6ab0,$7952,$cdb2,$4019),
            ($b3ad,$1c6d,$0f07,$8ba8,$c01d),
            ($8e0e,$e148,$5ab3,$ff44,$401e),
            ($a46a,$0273,$bc0f,$c358,$c01c),
            ($6de5,$b797,$ea1c,$e66b,$4019),
            ($e5e5,$4172,$8863,$b4a0,$c016),
            ($3c4f,$dc46,$b802,$e107,$4012),
            ($bed4,$3ad5,$2da1,$cc1d,$c00e),
            ($d0fe,$2487,$01c0,$8be3,$400a),
            ($1a6c,$1c93,$612a,$f742,$c004),
            ($0000,$0000,$0000,$8000,$3fff));
var
  y159n: array[0..9]  of extended absolute y159nh;
  y159d: array[0..10] of extended absolute y159dh;
  y1n:   array[0..6]  of extended absolute y1nhex;
  y1d:   array[0..7]  of extended absolute y1dhex;
begin
  {Ref: Cephes [7], file ldouble\j1l.c}
  if x<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_y1 := NaN_x;
    exit;
  end;
  if x < 9.0 then begin
    z := sqr(x);
    if z < 20.25 then begin
      {In the interval [0,4.5) a rational approximation of the form}
      {Y1(x) = x*P6(x)/Q7(x) + 2/Pi*(ln(x)*J1(x) - 1/x) is used.}
      y := (ln(x)*sfc_j1(x) - 1.0/x)/Pi_2;
      sfc_y1 := y + x*PolEvalX(z,y1n,7)/PolEvalX(z,y1d,8);
    end
    else begin
      {In the interval [4.5,9) a rational approximation of the form}
      {Y1(x) = (x - p)*(x - q)*(x - r)*(x - s)*P9(x)/Q10(x) is used}
      {where p, q, r, s are first four zeros of Y1(x).}
      y := (x - y1z1)*(x - y1z2)*(x - y1z3)*(x - y1z4);
      sfc_y1 := y * PolEvalX(x,y159n,10)/PolEvalX(x,y159d,11);
    end;
  end
  else begin
    {For x>=9 the common rational approximations to modulus}
    {and phase are used Y1(x) = modulus * sin(phase).}
    if x >= 1600 then sfc_y1 := bessy_large(1,x)
    else begin
      bess_m1p1(x,y,z);
      sfc_y1 := y*sin(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_jn(n: integer; x: extended): extended;
  {-Return J_n(x), the Bessel function of the 1st kind, order n; not suitable for large n or x.}
var
  curr,prev,q,temp,init,xh: extended;
  k: integer;
  neg: boolean;
const
  small = 7.143435E-1651;  {~ cbrt(succx(0))}
  lnsml = ln_succx0;
begin
  if IsNaNorInf(x) then begin
    sfc_jn := NaN_x;
    exit;
  end;
  {Based on boost_1_42_0\boost\math\special_functions\detail\bessel_jn.hpp [19]}
  {Copyright 2006 Xiaogang Zhang, see 3rdparty.ama for Boost license}
  init := Sqrt_MinExt;
  {Flag to negate result for |n|}
  neg := (n<0) and odd(n);
  n := abs(n);

  if n=0 then curr := sfc_j0(x)
  else if n=1 then curr := sfc_j1(x)
  else if abs(x) <= small then begin
    if (x=0.0) or (n>2) then curr := 0.0
    else curr := 0.125*sqr(x);
  end
  else begin
    xh := 0.5*x;
    if abs(x) > n then begin
      {forward recurrence}
      prev := sfc_j0(x);
      curr := sfc_j1(x);
      for k:=1 to n-1 do begin
        temp := curr*k/xh - prev;
        prev := curr;
        curr := temp;
      end;
    end
    else begin
      {Quick check if |J_n(x)| < MinExtended from HMF[1] 9.1.63}
      {solution of z*exp(sqrt(1-z^2))=1 is z = 0.39989.. }
      q := abs(x/n);
      if n<=100 then temp := 1e-40
      else if n<1000 then temp := 1e-4
      else if n<5000 then temp := 0.07
      else temp := 0.3999;
      if q < temp then begin
        {Jn(x) <= [q*exp(sqrt(1-q^2))/(1+sqrt(1-q^2))]^n}
        temp := sqrt(1.0 - q*q);
        temp := ln(q/(1.0+temp)) + temp;
        if temp < lnsml/n then begin
          sfc_jn := 0.0;
          exit;
        end;
      end;
      {set overflow threshold for iteration}
      q := 0.5*MaxExtended*q;
      {backward recurrence}
      CF1_j(n,x,prev,k);
      prev := prev*init;
      curr := init;
      for k:=n downto 1 do begin
        if abs(curr) > q then begin
          {prevent overflow and set result to zero}
          sfc_jn := 0.0;
          exit;
        end;
        temp := curr*k/xh - prev;
        prev := curr;
        curr := temp;
      end;
      curr := (init/curr)*sfc_j0(x);
    end;
  end;
  if neg then sfc_jn := -curr else sfc_jn := curr;
end;


{---------------------------------------------------------------------------}
function sfc_yn(n: integer; x: extended): extended;
  {-Return Y_n(x), the Bessel function of the 2nd kind, order n, x>0, not suitable for large n or x}
var
  yn,yn1,t: extended;
  k: integer;
  neg: boolean;
begin
  {Flag to negate result for |n|}
  neg := (n<0) and odd(n);
  n := abs(n);
  if x<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_yn := NaN_x;
    exit;
  end;
  if n=0 then yn := sfc_y0(x)
  else if n=1 then yn := sfc_y1(x)
  else if (n>MAXGAMX) and (x<=2.0) then yn := NegInf_x
  else begin
    if x<1.0 then begin
      {NIST[30] 10.7.4}
      t := sfc_lnfac(n-1) - n*ln(0.5*x) - LnPi;
      if t > ln_MaxExt then begin
        if neg then sfc_yn := PosInf_x else sfc_yn := NegInf_x;
        exit;
      end;
    end;
    {forward recurrence}
    yn1 := sfc_y0(x);
    yn  := sfc_y1(x);
    x := 0.5*x;
    for k:=1 to n-1 do begin
      t  := yn*k/x - yn1;
      yn1:= yn;
      yn := t;
    end;
  end;
  if neg then sfc_yn := -yn else sfc_yn := yn;
end;


{---------------------------------------------------------------------------}
function bess_i0_small(x: extended): extended;
  {-Return Bessel function I0(x) for abs(x)<=3, x assumed >= 0}
const
  xsml = 0.65854450798271924667e-9;  {sqrt(4*eps_x)}
const
  nbi0 = 13;
  bi0h : array[0..nbi0-1] of THexExtW = (
           ($4003,$A1EE,$5479,$9CE3,$BFFB), {-0.7660547252839144951081894976243285e-1}
           ($48B3,$5F19,$0294,$F6B3,$3FFF), {+0.1927337953993808269952408750881196e+1}
           ($B8D1,$AF86,$2883,$E9BE,$3FFC), {+0.2282644586920301338937029292330415e+0}
           ($B9F3,$936B,$1D6F,$D5CB,$3FF8), {+0.1304891466707290428079334210691888e-1}
           ($5C27,$B9AE,$D127,$E3C3,$3FF3), {+0.4344270900816487451378682681026107e-3}
           ($F200,$B849,$01B0,$9E16,$3FEE), {+0.9422657686001934663923171744118766e-5}
           ($CFA2,$6F56,$AA2C,$99F9,$3FE8), {+0.1434006289510691079962091878179957e-6}
           ($189D,$34A0,$4423,$DDCE,$3FE1), {+0.1613849069661749069915419719994611e-8}
           ($DA55,$222E,$86B5,$F5B3,$3FDA), {+0.1396650044535669699495092708142522e-10}
           ($4C24,$393C,$C78C,$D7B5,$3FD3), {+0.9579451725505445344627523171893333e-13}
           ($B4AE,$B291,$D6DC,$99BD,$3FCC), {+0.5333981859862502131015107744000000e-15}
           ($1BB6,$1CA4,$D573,$B56B,$3FC4), {+0.2458716088437470774696785919999999e-17}
           ($1BB2,$06F5,$B92D,$B41F,$3FBC));{+0.9535680890248770026944341333333333e-20}
var
  bi0: array[0..nbi0-1] of extended absolute bi0h;
begin
  {Ref: W. Fullerton [14] and [20], files dbesi0.f and dbsi0e.f}
  {Hex Chebyshev values calculated with mp_arith/t_rcalc}
  if x<=xsml then bess_i0_small := 1.0
  else bess_i0_small := 2.75 + CSEvalX(x*x/4.5-1.0,bi0,nbi0);
end;


{---------------------------------------------------------------------------}
function sfc_i0e(x: extended): extended;
  {-Return I0(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order zero}
const
  xsml = 0.23283064365386962891e-9;  {sqrt(0.5*eps_x)}
const
  nai0 = 28;
   ai0: array[0..nai0-1] of extended = (
          +0.7575994494023795942729872037438e-1,
          +0.7591380810823345507292978733204e-2,
          +0.4153131338923750501863197491382e-3,
          +0.1070076463439073073582429702170e-4,
          -0.7901179979212894660750319485730e-5,
          -0.7826143501438752269788989806909e-6,
          +0.2783849942948870806381185389857e-6,
          +0.8252472600612027191966829133198e-8,
          -0.1204463945520199179054960891103e-7,
          +0.1559648598506076443612287527928e-8,
          +0.2292556367103316543477254802857e-9,
          -0.1191622884279064603677774234478e-9,
          +0.1757854916032409830218331247743e-10,
          +0.1128224463218900517144411356824e-11,
          -0.1146848625927298877729633876982e-11,
          +0.2715592054803662872643651921606e-12,
          -0.2415874666562687838442475720281e-13,
          -0.6084469888255125064606099639224e-14,
          +0.3145705077175477293708360267303e-14,
          -0.7172212924871187717962175059176e-15,
          +0.7874493403454103396083909603327e-16,
          +0.1004802753009462402345244571839e-16,
          -0.7566895365350534853428435888810e-17,
          +0.2150380106876119887812051287845e-17,
          -0.3754858341830874429151584452608e-18,
          +0.2354065842226992576900757105322e-19,
          +0.1114667612047928530226373355110e-19,
          -0.5398891884396990378696779322709e-20);
const
  nai2 = 33;
  ai02: array[0..nai2-1] of extended = (
          +0.5449041101410883160789609622680e-1,
          +0.3369116478255694089897856629799e-2,
          +0.6889758346916823984262639143011e-4,
          +0.2891370520834756482966924023232e-5,
          +0.2048918589469063741827605340931e-6,
          +0.2266668990498178064593277431361e-7,
          +0.3396232025708386345150843969523e-8,
          +0.4940602388224969589104824497835e-9,
          +0.1188914710784643834240845251963e-10,
          -0.3149916527963241364538648629619e-10,
          -0.1321581184044771311875407399267e-10,
          -0.1794178531506806117779435740269e-11,
          +0.7180124451383666233671064293469e-12,
          +0.3852778382742142701140898017776e-12,
          +0.1540086217521409826913258233397e-13,
          -0.4150569347287222086626899720156e-13,
          -0.9554846698828307648702144943125e-14,
          +0.3811680669352622420746055355118e-14,
          +0.1772560133056526383604932666758e-14,
          -0.3425485619677219134619247903282e-15,
          -0.2827623980516583484942055937594e-15,
          +0.3461222867697461093097062508134e-16,
          +0.4465621420296759999010420542843e-16,
          -0.4830504485944182071255254037954e-17,
          -0.7233180487874753954562272409245e-17,
          +0.9921475412173698598880460939810e-18,
          +0.1193650890845982085504399499242e-17,
          -0.2488709837150807235720544916602e-18,
          -0.1938426454160905928984697811326e-18,
          +0.6444656697373443868783019493949e-19,
          +0.2886051596289224326481713830734e-19,
          -0.1601954907174971807061671562007e-19,
          -0.3270815010592314720891935674859e-20);
begin
  {Ref: W. Fullerton [14] and [20], file dbsi0e.f}
  x := abs(x);
  if x<=3.0 then begin
    {Note that there is bug in dbsi0e.f from [20] for small x. We use the}
    {Taylor series for I(0,x)*exp(-x) = 1 - x + 3/4*x^2 -5/12*x^3 + O(x^4)}
    if x<=xsml then sfc_i0e := 1 - x
    else sfc_i0e := exp(-x)*bess_i0_small(x);
  end
  else if x<=8.0 then begin
    sfc_i0e := (0.375 + CSEvalX((48.0/x-11.0)/5.0, ai0, nai0))/sqrt(x);
  end
  else begin
    sfc_i0e := (0.375 + CSEvalX(16.0/x-1.0, ai02, nai2))/sqrt(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_i0(x: extended): extended;
  {-Return I0(x), the modified Bessel function of the 1st kind, order zero}
begin
  x := abs(x);
  if x<=3.0 then sfc_i0 := bess_i0_small(x)
  else if x>ln_MaxExt then sfc_i0 := PosInf_x
  else sfc_i0 := sfc_i0e(x)*exp(x);
end;


{---------------------------------------------------------------------------}
function bess_i1_small(x: extended): extended;
  {-Return Bessel function I1(x) for abs(x)<=3}
var
  y: extended;
const
  xsml = 0.23283064365386962891e-9;  {sqrt(0.5*eps_x)}
const
  nbi1 = 12;
  bi1: array[0..nbi1-1] of extended = (
         -0.19717132610998597316138503218149e-2,
         +0.40734887667546480608155393652014e+0,
         +0.34838994299959455866245037783787e-1,
         +0.15453945563001236038598401058489e-2,
         +0.41888521098377784129458832004120e-4,
         +0.76490267648362114741959703966069e-6,
         +0.10042493924741178689179808037238e-7,
         +0.99322077919238106481371298054863e-10,
         +0.76638017918447637275200171681349e-12,
         +0.47414189238167394980388091948160e-14,
         +0.24041144040745181799863172032000e-16,
         +0.10171505007093713649121100799999e-18);
begin
  {Ref: W. Fullerton [14] and [20], files dbesi1.f and dbsi1e.f}
  y := abs(x);
  if y=0.0 then bess_i1_small := 0.0
  else if y<=xsml then bess_i1_small := 0.5*x
  else bess_i1_small := x*(0.875 + CSEvalX(x*x/4.5-1.0,bi1,nbi1));
end;


{---------------------------------------------------------------------------}
function sfc_i1e(x: extended): extended;
  {-Return I1(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order one}
var
  y: extended;
const
  nai1 = 28;
   ai1: array[0..nai1-1] of extended = (
          -0.2846744181881478674100372468307e-1,
          -0.1922953231443220651044448774979e-1,
          -0.6115185857943788982256249917785e-3,
          -0.2069971253350227708882823777979e-4,
          +0.8585619145810725565536944673138e-5,
          +0.1049498246711590862517453997860e-5,
          -0.2918338918447902202093432326697e-6,
          -0.1559378146631739000160680969077e-7,
          +0.1318012367144944705525302873909e-7,
          -0.1448423418183078317639134467815e-8,
          -0.2908512243993142094825040993010e-9,
          +0.1266388917875382387311159690403e-9,
          -0.1664947772919220670624178398580e-10,
          -0.1666653644609432976095937154999e-11,
          +0.1242602414290768265232168472017e-11,
          -0.2731549379672432397251461428633e-12,
          +0.2023947881645803780700262688981e-13,
          +0.7307950018116883636198698126123e-14,
          -0.3332905634404674943813778617133e-14,
          +0.7175346558512953743542254665670e-15,
          -0.6982530324796256355850629223656e-16,
          -0.1299944201562760760060446080587e-16,
          +0.8120942864242798892054678342860e-17,
          -0.2194016207410736898156266643783e-17,
          +0.3630516170029654848279860932334e-18,
          -0.1695139772439104166306866790399e-19,
          -0.1288184829897907807116882538222e-19,
          +0.5694428604967052780109991073109e-20);
const
  nai2 = 33;
  ai12: array[0..nai2-1] of extended = (
         +0.2857623501828012047449845948469e-1,
         -0.9761097491361468407765164457302e-2,
         -0.1105889387626237162912569212775e-3,
         -0.3882564808877690393456544776274e-5,
         -0.2512236237870208925294520022121e-6,
         -0.2631468846889519506837052365232e-7,
         -0.3835380385964237022045006787968e-8,
         -0.5589743462196583806868112522229e-9,
         -0.1897495812350541234498925033238e-10,
         +0.3252603583015488238555080679949e-10,
         +0.1412580743661378133163366332846e-10,
         +0.2035628544147089507224526136840e-11,
         -0.7198551776245908512092589890446e-12,
         -0.4083551111092197318228499639691e-12,
         -0.2101541842772664313019845727462e-13,
         +0.4272440016711951354297788336997e-13,
         +0.1042027698412880276417414499948e-13,
         -0.3814403072437007804767072535396e-14,
         -0.1880354775510782448512734533963e-14,
         +0.3308202310920928282731903352405e-15,
         +0.2962628997645950139068546542052e-15,
         -0.3209525921993423958778373532887e-16,
         -0.4650305368489358325571282818979e-16,
         +0.4414348323071707949946113759641e-17,
         +0.7517296310842104805425458080295e-17,
         -0.9314178867326883375684847845157e-18,
         -0.1242193275194890956116784488697e-17,
         +0.2414276719454848469005153902176e-18,
         +0.2026944384053285178971922860692e-18,
         -0.6394267188269097787043919886811e-19,
         -0.3049812452373095896084884503571e-19,
         +0.1612841851651480225134622307691e-19,
         +0.3560913964309925054510270904620e-20);
begin
  {Ref: W. Fullerton [14] and [20], file dbsi1e.f}
  y := abs(x);
  if y<=3.0 then sfc_i1e := exp(-y)*bess_i1_small(x)
  else begin
    if y<=8.0 then begin
      y := (0.375 + CSEvalX((48.0/y-11.0)/5.0, ai1, nai1))/sqrt(y)
    end
    else begin
      y := (0.375 + CSEvalX(16.0/y-1.0, ai12, nai2))/sqrt(y);
    end;
    if x>0 then sfc_i1e := y else sfc_i1e := -y;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_i1(x: extended): extended;
  {-Return I1(x), the modified Bessel function of the 1st kind, order one}
var
  y: extended;
begin
  y := abs(x);
  if y<=3.0 then sfc_i1 := bess_i1_small(x)
  else if x>ln_MaxExt then sfc_i1 := PosInf_x
  else sfc_i1 := sfc_i1e(x)*exp(y);
end;


{---------------------------------------------------------------------------}
function bess_k0_small(x: extended): extended;
  {-Return Bessel function K0(x) for 0 < x <= 2}
var
  y: extended;
const
  xsml = 0.46566128730773925781e-9;  {sqrt(2*eps_x)}
const
  nbk0 = 12;
  bk0: array[0..nbk0-1] of extended = (
         -0.353273932339027687201140060063153e-1,
         +0.344289899924628486886344927529213e+0,
         +0.359799365153615016265721303687231e-1,
         +0.126461541144692592338479508673447e-2,
         +0.228621210311945178608269830297585e-4,
         +0.253479107902614945730790013428354e-6,
         +0.190451637722020885897214059381366e-8,
         +0.103496952576336245851008317853089e-10,
         +0.425981614279108257652445327170133e-13,
         +0.137446543588075089694238325440000e-15,
         +0.357089652850837359099688597333333e-18,
         +0.763164366011643737667498666666666e-21);
begin
  {Ref: W. Fullerton [14] and [20], files dbesk0.f and dbsk0e.f}
  if x>xsml then y := x*x
  else begin
    if x=0.0 then begin
      bess_k0_small := PosInf_x;
      exit;
    end
    else y := 0.0;
  end;
  bess_k0_small := -ln(0.5*x)*bess_i0_small(x) - 0.25 + CSEvalX(0.5*y-1.0,bk0,nbk0);
end;


{---------------------------------------------------------------------------}
function sfc_k0e(x: extended): extended;
  {-Return K0(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order zero, x>0}
const
  nak0 = 22;
   ak0: array[0..nak0-1] of extended = (
          -0.7643947903327941424082978270088e-1,
          -0.2235652605699819052023095550791e-1,
          +0.7734181154693858235300618174047e-3,
          -0.4281006688886099464452146435416e-4,
          +0.3081700173862974743650014826660e-5,
          -0.2639367222009664974067448892723e-6,
          +0.2563713036403469206294088265742e-7,
          -0.2742705549900201263857211915244e-8,
          +0.3169429658097499592080832873403e-9,
          -0.3902353286962184141601065717962e-10,
          +0.5068040698188575402050092127286e-11,
          -0.6889574741007870679541713557984e-12,
          +0.9744978497825917691388201336831e-13,
          -0.1427332841884548505389855340122e-13,
          +0.2156412571021463039558062976527e-14,
          -0.3349654255149562772188782058530e-15,
          +0.5335260216952911692145280392601e-16,
          -0.8693669980890753807639622378837e-17,
          +0.1446404347862212227887763442346e-17,
          -0.2452889825500129682404678751573e-18,
          +0.4233754526232171572821706342400e-19,
          -0.7427946526454464195695341294933e-20);
const
  nak2 = 18;
  ak02: array[0..nak2-1] of extended = (
          -0.1201869826307592239839346212452e-1,
          -0.9174852691025695310652561075713e-2,
          +0.1444550931775005821048843878057e-3,
          -0.4013614175435709728671021077879e-5,
          +0.1567831810852310672590348990333e-6,
          -0.7770110438521737710315799754460e-8,
          +0.4611182576179717882533130529586e-9,
          -0.3158592997860565770526665803309e-10,
          +0.2435018039365041127835887814329e-11,
          -0.2074331387398347897709853373506e-12,
          +0.1925787280589917084742736504693e-13,
          -0.1927554805838956103600347182218e-14,
          +0.2062198029197818278285237869644e-15,
          -0.2341685117579242402603640195071e-16,
          +0.2805902810643042246815178828458e-17,
          -0.3530507631161807945815482463573e-18,
          +0.4645295422935108267424216337066e-19,
          -0.6368625941344266473922053461333e-20);
begin
  {Ref: W. Fullerton [14] and [20], file dbsk0e.f}
  if x<=2.0 then sfc_k0e := exp(x)*bess_k0_small(x)
  else if x<=8.0 then begin
    sfc_k0e := (1.25 + CSEvalX((16.0/x-5.0)/THREE, ak0, nak0))/sqrt(x);
  end
  else begin
    sfc_k0e := (1.25 + CSEvalX(16.0/x-1.0, ak02, nak2))/sqrt(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_k0(x: extended): extended;
  {-Return K0(x), the modified Bessel function of the 2nd kind, order zero, x>0}
begin
  if x<=2.0 then sfc_k0 := bess_k0_small(x)
  else sfc_k0 := sfc_k0e(x)*exp(-x);
end;


{---------------------------------------------------------------------------}
function bess_k1_small(x: extended): extended;
  {-Return Bessel function K1(x) for 0 < x <= 2}
var
  y: extended;
const
  xsml = 0.46566128730773925781e-9;  {sqrt(2*eps_x)}
const
  nbk1 = 12;
  bk1: array[0..nbk1-1] of extended = (
         +0.25300227338947770532531120868533e-1,
         -0.35315596077654487566723831691801e+0,
         -0.12261118082265714823479067930042e+0,
         -0.69757238596398643501812920296083e-2,
         -0.17302889575130520630176507368979e-3,
         -0.24334061415659682349600735030164e-5,
         -0.22133876307347258558315252545126e-7,
         -0.14114883926335277610958330212608e-9,
         -0.66669016941993290060853751264373e-12,
         -0.24274498505193659339263196864853e-14,
         -0.70238634793862875971783797120000e-17,
         -0.16543275155100994675491029333333e-19);
begin
  {Ref: W. Fullerton [14] and [20], files dbesk1.f and dbsk1e.f}
  if x>xsml then y := x*x
  else begin
    if x=0.0 then begin
      bess_k1_small := PosInf_x;
      exit;
    end
    else y := 0.0;
  end;
  bess_k1_small := ln(0.5*x)*bess_i1_small(x) + (0.75 + CSEvalX(0.5*y-1.0,bk1,nbk1))/x;
end;


{---------------------------------------------------------------------------}
function sfc_k1e(x: extended): extended;
  {-Return K1(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order one, x>0}
const
  nak1 = 22;
   ak1: array[0..nak1-1] of extended = (
          +0.27443134069738829695257666227266e+0,
          +0.75719899531993678170892378149290e-1,
          -0.14410515564754061229853116175625e-2,
          +0.66501169551257479394251385477036e-4,
          -0.43699847095201407660580845089167e-5,
          +0.35402774997630526799417139008534e-6,
          -0.33111637792932920208982688245704e-7,
          +0.34459775819010534532311499770992e-8,
          -0.38989323474754271048981937492758e-9,
          +0.47208197504658356400947449339005e-10,
          -0.60478356628753562345373591562890e-11,
          +0.81284948748658747888193837985663e-12,
          -0.11386945747147891428923915951042e-12,
          +0.16540358408462282325972948205090e-13,
          -0.24809025677068848221516010440533e-14,
          +0.38292378907024096948429227299157e-15,
          -0.60647341040012418187768210377386e-16,
          +0.98324256232648616038194004650666e-17,
          -0.16284168738284380035666620115626e-17,
          +0.27501536496752623718284120337066e-18,
          -0.47289666463953250924281069568000e-19,
          +0.82681500028109932722392050346666e-20);
const
  nak2 = 18;
  ak12: array[0..nak2-1] of extended = (
          +0.6379308343739001036600488534102e-1,
          +0.2832887813049720935835030284708e-1,
          -0.2475370673905250345414545566732e-3,
          +0.5771972451607248820470976625763e-5,
          -0.2068939219536548302745533196552e-6,
          +0.9739983441381804180309213097887e-8,
          -0.5585336140380624984688895511129e-9,
          +0.3732996634046185240221212854731e-10,
          -0.2825051961023225445135065754928e-11,
          +0.2372019002484144173643496955486e-12,
          -0.2176677387991753979268301667938e-13,
          +0.2157914161616032453939562689706e-14,
          -0.2290196930718269275991551338154e-15,
          +0.2582885729823274961919939565226e-16,
          -0.3076752641268463187621098173440e-17,
          +0.3851487721280491597094896844799e-18,
          -0.5044794897641528977117282508800e-19,
          +0.6888673850418544237018292223999e-20);
begin
  {Ref: W. Fullerton [14] and [20], file dbsk1e.f}
  if x<=2.0 then sfc_k1e := exp(x)*bess_k1_small(x)
  else if x<=8.0 then begin
    sfc_k1e := (1.25 + CSEvalX((16.0/x-5.0)/THREE, ak1, nak1))/sqrt(x);
  end
  else begin
    sfc_k1e := (1.25 + CSEvalX(16.0/x-1.0, ak12, nak2))/sqrt(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_k1(x: extended): extended;
  {-Return K1(x), the modified Bessel function of the 2nd kind, order one, x>0}
begin
  if x<=2.0 then sfc_k1 := bess_k1_small(x)
  else sfc_k1 := sfc_k1e(x)*exp(-x);
end;


{---------------------------------------------------------------------------}
function sfc_in(n: integer; x: extended): extended;
  {-Return I_n(x), the modified Bessel function of the 1st kind, order n}
var
  curr,prev,temp,init,y: extended;
  k: integer;
const
  NMax = 256;     {double: 160}
  XMax = 1024.0;  {double: 512}
begin
  if IsNaNorInf(x)  then begin
    sfc_in := NaN_x;
    exit;
  end;
  {HMF[1], 9.6.6: I(-n,x) = I(n,x)}
  n := abs(n);
  if n=0 then sfc_in := sfc_i0(x)
  else if n=1 then sfc_in := sfc_i1(x)
  else if x=0.0 then sfc_in := 0.0
  else begin
    y := abs(x);
    {If n or x are not small then use real order function}
    if (n>NMax) or (y>XMax) then sfc_in := sfc_iv(n,x)
    else begin
      {Simplified calculation for I_n only with existing functions}
      if y <= 2.0 then begin
        curr := IJ_series(n,y,false);
      end
      else begin
        {I_(n+1)(y)/I_n(y) using continued fraction and recurse backwards}
        CF1_I(n,y,temp);
        y := 0.5*y;
        init := Sqrt_MinExt;
        curr := init;
        prev := init*temp;
        for k:=n downto 1 do begin
          temp := curr*k/y + prev;
          prev := curr;
          curr := temp;
        end;
        curr := (init/curr)*sfc_i0(x);
      end;
      {Adjust sign}
      if (x<0.0) and odd(n) then curr := -curr;
      sfc_in := curr;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_kn(n: integer; x: extended): extended;
  {-Return K_n(x), the modified Bessel function of the 2nd kind, order n, x>0, not suitable for large n}
var
  kn,knm1,knm2: extended;
  k: integer;
begin
  {HMF[1], 9.6.6: K(-n,x) = K(n,x)}
  n := abs(n);
  {Range error for x<=0 is generated in k0 or k1}
  if n=0 then kn := sfc_k0(x)
  else if n=1 then kn := sfc_k1(x)
  else begin
    {avoid false warning "Variable 'kn' might not have been initialized"}
    kn := 0.0;
    {forward recurrence, K(n+1,x) = 2n/x*K(n,x) + K(n-1,x)}
    knm2 := sfc_k0(x);
    knm1 := sfc_k1(x);
    x := 0.5*x;
    for k:=1 to n-1 do begin
      kn   := knm1*k/x + knm2;
      knm2 := knm1;
      knm1 := kn;
    end;
  end;
  sfc_kn := kn;
end;


{---------------------------------------------------------------------------}
{-------------------- Bessel functions of real order -----------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure temme_y(v,x: extended; var Y,Y1: extended);
  {-Calculate Y(v, x) and Y(v+1, x) by Temme's method for small |x|}
var
  k,g,h,p,q,f,c,d,s,s1,tol,a,e: extended;
  g1,g2,gp,gm,v2: extended;
const
  MAXIT = 30000;
begin

{$ifdef debug}
  if (abs(v) > 0.5) or (abs(x) > 2) then begin
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
  end;
{$endif}

  {N.M.Temme [51], On the Numerical Evaluation of the Ordinary}
  {Bessel Function of the Second Kind, Section 2.}
  gp := sfc_gamma1pm1(v);
  gm := sfc_gamma1pm1(-v);
  a  := ln(0.5*x);
  s  := -a*v;

  if abs(v) < eps_x then begin
    e := 0.5*v*sqr(Pi);
    d := 1.0/Pi;
  end
  else begin
    e := 2.0*sqr(sinpi(0.5*v))/v;
    d := v/sinpi(v);
  end;

  if v=0.0 then g1 := -EulerGamma
  else g1 := 0.5*(gp-gm)/((1.0+gp)*(1.0+gm)*v);
  g2 := 0.5*(2.0+gp+gm)/((1.0+gp)*(1.0+gm));

  f := 2.0*(g1*cosh(s) - g2*a*sinhc(s))*d;
  c := power(0.5*x, v);
  p := d/(c*(1.0 + gm));
  q := d*c/(1.0 + gp);

  g := f + e*q;
  c := 1.0;
  s := c*g;
  s1:= c*p;

  v2 := v*v;
  d  := -0.25*x*x;

  {series summation}
  tol := 0.5*eps_x;
  {use extended k because otherwise k*k may overflow}
  k := 1.0;
  while k <= MAXIT do begin
    c := c*d/k;
    f := (k*f + p + q) / (k*k - v2);
    p := p/(k - v);
    q := q/(k + v);
    g := f + e*q;
    h := p - k*g;
    s := s  + c*g;
    s1:= s1 + c*h;
    if abs(c*g) < abs(s)*tol then begin
      Y  := -s;
      Y1 := -2.0*s1/x;
      exit;
    end;
    k := k + 1.0;
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  Y  := -s;
  Y1 := -2.0*s1/x;
end;


{---------------------------------------------------------------------------}
procedure CF2_jy(v,x: extended; var p,q: extended);
  {-Return the continued fraction p + i*q = (J' + iY') / (J + iY)}
var
  a,br,bi,cr,ci,dr,di,er,ei,fr,fi,t: extended;
  i: integer;
const
  MAXIT = 30000;
begin
  {Ref: Numerical Recipes [13], ch. 6.7, p.244, function bessjy}
  {Evaluate the continued fraction p + iq = (J' + iY') / (J + iY)}
  {NR [13] (6.7.3) using the (complex) modified Lentz's method.}

{$ifdef debug}
  if abs(x) < 1.0 then begin
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
  end;
{$endif}

  t  := 1.0/x;
  a  := 0.25-v*v;
  fr := -0.5*t;
  fi := 1.0;
  br := 2.0*x;
  bi := 2.0;
  t  := a*t/(fr*fr+fi*fi);
  cr := br+t*fi;
  ci := bi+t*fr;
  t  := br*br+bi*bi;
  dr := +br/t;
  di := -bi/t;
  er := cr*dr-ci*di;
  ei := cr*di+ci*dr;
  t  := fr*er-fi*ei;
  fi := fr*ei+fi*er;
  fr := t;
  for i:=1 to MAXIT do begin
    a  := a + 2.0*i;
    bi := bi + 2.0;
    dr := a*dr+br;
    di := a*di+bi;
    if abs(dr)+abs(di) = 0.0 then dr := Sqrt_MinExt;
    t  := a/(cr*cr+ci*ci);
    cr := br+t*cr;
    ci := bi-t*ci;
    if abs(cr)+abs(ci) = 0.0 then cr := Sqrt_MinExt;
    t  := dr*dr+di*di;
    dr := +dr/t;
    di := -di/t;
    er := cr*dr-ci*di;
    ei := cr*di+ci*dr;
    t  := fr*er-fi*ei;
    fi := fr*ei+fi*er;
    fr := t;
    if abs(er-1.0)+abs(ei) < 8*eps_x then begin
      p := fr;
      q := fi;
      exit;
    end;
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  p := fr;
  q := fi;
end;


{---------------------------------------------------------------------------}
function Yv_series(v,x: extended): extended;
  {-Series for Yv for 'small' x < 1;  frac(v)<>0}
var
  h,k,v1,v2,t1,t2,xx,Yv: extended;
const
  hsmall: THexExtW = ($4153,$E4E4,$9836,$0A2F,$0000); {2.675476672E-4933 = ~ 1/(MaxExtended*Pi)}
begin

  {Use Yv(x) = ( Jv(v,x)cos(vx) - Jv(-v,x) )/sin(vx) }
  {and Gamma reflection to sum the two series for Jv }

  if v<MAXGAMX then h := power(0.5*x, v)/sfc_gamma(v)
  else h := exp(v*ln(0.5*x) - sfc_lngamma(v));

  if h <= extended(hsmall) then begin
    {1.0/(h*Pi) will overflow}
    Yv_series := NegInf_x;
    exit;
  end;

  sincosPi(v, v1, v2);
  t2 := 1.0/(h*Pi);
  t1 := h*(v2/(v*v1));

  Yv := t1-t2;
  xx := 0.25*x*x;
  v2 := v;
  v1 := v;
  k  := 0.0;
  repeat
    k  := k + 1.0;
    v1 := v1 + 1.0;
    v2 := v2 - 1.0;
    t1 := -t1*xx/(k*v1);
    t2 :=  t2*xx/(k*v2);
    h  := t1 - t2;
    Yv := Yv + h;
  until (abs(h)<eps_x*abs(Yv));
  Yv_series := Yv;
end;


const
  BT_J = 1; BT_Y = 2;

{---------------------------------------------------------------------------}
procedure bessel_jy(v,x: extended; BT: byte; var Jv,Yv: extended);
  {-Return J_v(x) and/or Y_v(x) depending on BT, x>0, |v| < MaxLongint, INF if overflow}
var
  n,k: longint;
  u, Ju, Yv1, Yu, Yu1, fv, fu: extended;
  w, p, q, g, curr, prev, next, init, t, x2: extended;
  reflect: boolean;
  s: integer;
const
  lnepsh = -44.3614195558364998;  {ln(0.5*eps_x)}

  {--------------------------------------------}
  function rec_overflow(a,b: extended): boolean;
    {-Test if a*b overflows, if yes set Yv and Jv = PosInf_x}
  begin
    if (abs(a) > 1.0) and (abs(b) >= MaxExtended/abs(a)) then begin
      rec_overflow := true;
      Jv := PosInf_x;
      Yv := PosInf_x;
    end
    else rec_overflow := false;
  end;

begin

  {Ref: Boost [19] file bessel_jy.hpp, function bessel_jy}
  {and NR [13], Ch.6.7, section Ordinary Bessel Functions}

  {For x < 0 the functions Jv and Yv are in general complex; and Yv}
  {is singular for x=0. |v| < MaxLongint is assumed, so we can work}
  {with longint in the recurrence iterative, but these routines are}
  {are not suitable for large v values anyway.}

  if abs(v) >= MaxLongint - 1 then begin
    Jv := NaN_x;
    Yv := NaN_x;
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    exit;
  end;
  reflect := v < 0.0;
  if reflect then begin
    v  := -v;
    BT := BT_J + BT_Y;   {J and Y needed for reflection formula}
  end;

  x2 := 0.5*x;
  n  := round(v);
  u  := v-n;
  w  := 2.0/(Pi*x); {Wronskian}

  if x<=2.0 then begin
    if v>MAXGAMX then begin
      Yv := PosInf_x;
      if reflect then Jv := PosInf_x else Jv := 0.0;
      exit;
    end
    else begin
      {Check very 'small' x and v case with (near) overflow for Yv}
      if (x<1.0) and (u<>0.0) and (lnepsh > v*ln(0.25*sqr(x)/v)) then begin
        if BT and BT_J <> 0 then Jv := IJ_series(v, x, true);
        if BT and BT_Y <> 0 then Yv := Yv_series(v, x);
      end
      else begin
        temme_y(u, x, Yu,Yu1);
        if n=0 then begin
          Yv  := Yu;
          Yv1 := Yu1;
        end
        else begin
          prev := Yu;
          curr := Yu1;
          {skip last next calculation if J is not needed: it}
          {overflows in some cases like bessel_yv(1755.45,2)}
          {bessel_jv(1755.45,2) will produce an overflow!   }
          for k:=1 to n-1 do begin
            {forward recurrence for Y}
            t := (u+k)/x2;
            if rec_overflow(t,curr) then exit;
            next := t*curr - prev;
            prev := curr;
            curr := next;
          end;
          Yv  := curr;
          if BT and BT_J = 0 then Yv1 := 0.0  {keep some compilers quiet!}
          else Yv1 := ((u+n)/x2)*curr - prev;
        end;
        if BT and BT_J <> 0 then begin
          CF1_j(v, x, fv, s);
          Jv := w/(Yv*fv - Yv1);   {Wronskian relation}
        end;
      end;
    end;
  end
  else begin
    {calculate the lower limit t=t(v) for asymptotic range}
    if BT=BT_Y then t:= 1552
    else begin
      t := maxx(3.0, v*v)*121.0;
      if BT and BT_Y <> 0 then t := maxx(t,1552);
    end;
    if x>t then begin
      {Use asymptotic expansion of Hankel function H1v(x)}
      if BT and BT_J <> 0 then Jv := bessj_large(v, x);
      if BT and BT_Y <> 0 then begin
        Yu  := bessy_large(u, x);
        Yu1 := bessy_large(u + 1, x);
      end;
    end
    else begin
      CF1_j(v, x, fv, s);
      {tiny initial value to prevent overflow}
      init := Sqrt_MinExt;
      curr := s*init;
      prev := fv*curr;
      for k:=n downto 1 do begin
        {backward recurrence for J}
        t := (u+k)/x2;
        if rec_overflow(t,curr) then exit;
        next := t*curr - prev;
        prev := curr;
        curr := next;
      end;
      CF2_jy(u, x, p, q);
      fu := prev/curr;
      t  := u/x - fu;   {t = J'/J}
      g  := (p-t)/q;
      Ju := sqrt(w/(q + g*(p-t)));
      if curr<0.0 then Ju := -Ju;
      Jv := s*Ju*(init/curr); {normalization}
      Yu := g*Ju;
      Yu1:= Yu*(u/x - p - q/g);
    end;
    if BT and BT_Y <> 0 then begin
      prev := Yu;
      curr := Yu1;
      if n=0 then Yv := prev
      else begin
        for k:=1 to n-1 do begin
          {forward recurrence for Y}
          t := (u+k)/x2;
          if rec_overflow(t,curr) then exit;
          next := t*curr - prev;
          prev := curr;
          curr := next;
        end;
        Yv := curr;
      end;
    end;
  end;

  if reflect then begin
    {For negative v use reflection formula NR [13], 6.7.19}
    {J(-v,x) = cos(Pi*v)*J(v,x) - sin(Pi*v)*Y(v,x)}
    {Y(-v,x) = sin(Pi*v)*J(v,x) + cos(Pi*v)*Y(v,x)}
    sincosPi(v,p,q);
    t  := q*Jv - p*Yv;
    Yv := p*Jv + q*Yv;
    Jv := t;
  end;
end;


{---------------------------------------------------------------------------}
procedure sfc_bess_jyv(v,x: extended; var Jv,Yv: extended);
  {-Return J_v(x) and Y_v(x), no checks, x>0, |v| < MaxLongint}
begin
  bessel_jy(v,x, BT_J + BT_Y, jv,yv);
end;


{---------------------------------------------------------------------------}
function sfc_jv(v, x: extended): extended;
  {-Return J_v(x), the Bessel function of the 1st kind, order v; not suitable for large v.}
var
  r,n: extended;
begin

  if IsNaNorInf(x) or IsNaNorInf(v) then begin
    sfc_jv := NaN_x;
    exit;
  end;

  n := int(v);
  if x<0.0 then begin
    if n=v then begin
      r := sfc_jv(v, -x);
      if frac(0.5*v)<>0 then r := -r;
      sfc_jv := r;
    end
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfc_jv := NaN_x;
    end;
  end
  else if x=0.0 then begin
    if v=0.0 then sfc_jv := 1.0
    else if v>0.0 then sfc_jv := 0.0
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfc_jv := NaN_x;
    end;
  end
  else begin
    {here x > 0}
    if n=v then begin
      {integer order}
      if abs(n)<200.0 then begin
        if x > maxx(3.0, v*v)*121.0 then begin
          r := bessj_large(abs(v),x);
          if frac(0.5*v)<0.0 then r := -r;
          sfc_jv := r;
        end
        else sfc_jv := sfc_jn(round(n),x);
        exit;
      end;
    end;
    {Here v no integer or |v| > 200}
    if ((v >= 0.0) and (v<MAXGAMX-1)) and ((x < 1.0) or (v > 0.25*x*x)) then begin
      sfc_jv := IJ_series(v, x, true);
    end
    else begin
      bessel_jy(v,x,BT_J,r,n);
      sfc_jv := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_yv(v, x: extended): extended;
  {-Return Y_v(x), the Bessel function of the 2nd kind, order v; x > 0; not suitable for large v.}
var
  r,n: extended;
begin

  if IsNaNorInf(x) or IsNaNorInf(v) or (x<=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_yv := NaN_x;
    exit;
  end;

  r := abs(v);
  if r=0.5 then begin
    r := sqrt(Pi_2*x);
    if v<0.0 then sfc_yv := sin(x)/r
    else sfc_yv := -cos(x)/r;
  end
  else if frac(v) = -0.5 then begin
    {Reflection would be used but cos(Pi*v)=0: Y(v,x) = sin(Pi*v)*J(-v,x)}
    n := sfc_jv(-v,x);
    if frac(0.5*(r-0.5))=0.0 then sfc_yv := n else sfc_yv := -n;
  end
  else begin
    n := int(v);
    if n=v then begin
      {integer order}
      if (x>1552.0) and (x>5.0*abs(v)) then begin
        r := bessy_large(abs(v),x);
        if frac(0.5*v)<0.0 then r := -r;
        sfc_yv := r;
      end
      else if abs(n)<2000 then begin
        sfc_yv := sfc_yn(round(n),x);
      end
      else begin
        {Call general routine but avoid Jv calculation for v<0}
        bessel_jy(abs(v),x,BT_Y,n,r);
        if frac(0.5*v)<0.0 then r := -r;
        sfc_yv := r;
      end;
      exit;
    end
    else begin
      bessel_jy(v,x,BT_Y,n,r);
      sfc_yv := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure temme_k(v,x: extended; var K0,K1: extended);
  {-Calculate K(v, x) and K(v+1, x) by Temme's method for small |x|}
var
  k,h,p,q,f,c,d,s,s1,tol,a: extended;
  g1,g2,gp,gm: extended;
const
  MAXIT = 30000;
begin

{$ifdef debug}
  if (abs(v) > 0.5) or (abs(x) > 2.0) then begin
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
  end;
{$endif}

  {N.M. Temme [52], On the numerical evaluation of the }
  {modified Bessel function of the third kind. See also}
  {Boost [19], file bessel_ik.hpp / function temme_ik  }
  gp := sfc_gamma1pm1(v);
  gm := sfc_gamma1pm1(-v);
  a  := ln(0.5*x);
  s  := -a*v;
  h  := exp(a*v);
  if abs(v) < eps_x then begin
    c  := 1.0;
    g1 := -EulerGamma;
  end
  else begin
    c  := sinPi(v)/(v*Pi);
    g1 := (0.5*c/v)*(gp-gm);
  end;
  g2 := 0.5*c*(2.0+gp+gm);
  {initial values}
  p := 0.5*(1.0+gp)/h;
  q := 0.5*(1.0+gm)*h;
  f := (g1*cosh(s) - a*g2*sinhc(s))/c;
  h := p;
  c := 1.0;
  s := c*f;
  s1:= c*h;

  a := v*v;
  d := 0.25*x*x;

  {series summation}
  tol := 0.5*eps_x;
  k := 1.0;
  while k <= MAXIT do begin
    f := (k*f + p + q) / (k*k - a);
    p := p/(k - v);
    q := q/(k + v);
    h := p - k*f;
    c := c*d/k;
    s := s  + c*f;
    s1:= s1 + c*h;
    if abs(c*f) < abs(s)*tol then begin
      K0 := s;
      K1 := 2.0*s1/x;
      exit;
    end;
    k := k + 1.0;
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  K0 := s;
  K1 := 2.0*s1/x;
end;


{---------------------------------------------------------------------------}
procedure CF2_K(v,x: extended; escale: boolean; var K0,K1: extended);
  {-Compute K(v,x) and K(v+1,x) via continued fraction, |v| <= 0.5, |x| > 1}
  { If escale=true the values are multiplied by exp(x)}
var
  a,a1,b,c,d,dh,ds,h,q,q1,q2,s,t: extended;
  k: integer;
const
  MAXIT = 30000;
label
  done;
begin
  {Ref: Numerical Recipes [13], Ch. 6.7, section Modified Bessel Functions}
  {and p.249, function bessik. It is based on I.J. Thompson, A.R. Barnett,}
  {1987, Computer Physics Communications, vol. 47, pp. 245-257.           }

{$ifdef debug}
  if (abs(v) > 0.5) or (abs(x) <= 1.0) then begin
    if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
  end;
{$endif}

  b  := 2.0*(1.0+x);
  d  := 1.0/b;
  h  := d;
  dh := d;
  q1 := 0.0;
  q2 := 1.0;
  a  := v*v - 0.25;
  a1 := -a;
  q  := a1;
  c  := a1;
  s  := 1.0+q*dh;
  for k:=2 to MAXIT do begin
    a  := a - 2*(k-1);
    c  := -a*c/k;
    t  := (q1-b*q2)/a;
    q1 := q2;
    q2 := t;
    q  := q + c*t;
    b  := b + 2.0;
    d  := 1.0/(b+a*d);
    dh := (b*d-1.0)*dh;
    h  := h + dh;
    ds := q*dh;
    s  := s + ds;
    if abs(ds) < abs(s)*eps_x then goto done;
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));

done:

  K0 := sqrt(Pi_2/x)/s;
  if not escale then K0 := K0*exp(-x);
  K1 := K0*(v+x+0.5-a1*h)/x;
end;


{---------------------------------------------------------------------------}
function bess_i_large(v,x: extended; var Ivs: extended): boolean;
  {-Compute I_v(x) for large x >= 100,  sqrt(x) >= 2v. Return true}
  { if 'convergence' within 50 iterations, Ivs = I_v(x)*exp(-x)*sqrt(2Pix)}
var
  s,t,u,w: extended;
  k: integer;
begin
  {Hankel expansion of Iv, NIST[30], 10.40.1 and 10.17.1}
  u := 4.0*sqr(v);
  w := 8.0*x;
  {Typical values: w >= 800,  u/w < 1/8}
  t := 1.0;
  s := 1.0;
  bess_i_large := false;
  for k:=1 to 50 do begin
    t := t*(sqr(2*k-1)-u)/k/w;
    s := s + t;
    if abs(t)<eps_x then begin
      bess_i_large := true;
      Ivs := s;
      exit;
    end;
  end;
  Ivs := s;
end;


{---------------------------------------------------------------------------}
procedure bessel_ik(v,x: extended; CalcI, escale: boolean; var Iv,Kv: extended);
  {-Return I_v(x) and/or K_v(x) depending on CalcI, x>0, |v| < MaxLongint}
  { If escale=true the values are exponentially scaled.}
var
  n,k: longint;
  u, Kv1, Ku, Ku1, fv: extended;
  w, curr, prev, next, t, x2: extended;
  reflect,OK,kzero: boolean;
begin

  {Ref: Boost [19] file bessel_ik.hpp, function bessel_ik}
  {and NR [13], Ch.6.7, section Modified Bessel Functions}

  reflect := v < 0.0;
  v  := abs(v);

  x2 := 0.5*x;
  n  := round(v);
  u  := v-n;

  if x <= 2.0 then begin
    temme_k(u, x, Ku, Ku1);
    if escale then begin
      fv := exp(x);
      Ku := Ku*fv;
      Ku1:= Ku1*fv;
    end;
  end
  else CF2_K(u, x, escale, Ku, Ku1);

  kzero := (abs(Ku1)+abs(Ku))=0.0;
  if not kzero then begin
    prev := Ku;
    curr := Ku1;
    for k:=1 to n do begin
      {forward recurrence for K}
      t := (u+k)/x2;
      if (t > 1.0) and (curr >= MaxExtended/t) then begin
        Kv := PosInf_x;
        if CalcI then begin
          if Reflect then Iv := PosInf_x else Iv := 0.0;
        end;
        exit;
      end;
      next := t*curr + prev;
      prev := curr;
      curr := next;
    end;
    Kv  := prev;
    Kv1 := curr;
    kzero := (abs(Kv1)+abs(Kv))=0.0;
  end
  else begin
    Kv  := 0.0;
    Kv1 := 0.0;
  end;

  if CalcI then begin
    if (x >= 100.0) and (2.0*v <= sqrt(x)) then begin
      {Asymptotic expansion HMF[1], 9.7.1 or NIST[30], 10.40.4}
      OK := bess_i_large(v,x,t);
      if not escale then begin
        {Even if no convergence the result is used for ovrflow check}
        if (t <= 0.0) or (x+ln(t)-0.5*ln(TwoPi*x) >= ln_MaxExt) then begin
          Iv := PosInf_x;
          exit;
        end;
      end;
      if OK then begin
        if escale then Iv := t/sqrt(TwoPi*x)
        else begin
          u  := exp(0.5*x);
          Iv := u*(u*(t/sqrt(TwoPi*x)));
        end;
      end;
    end
    else begin
      if kzero then begin
        Iv := PosInf_x;
        exit;
      end;
      CF1_I(v, x, fv);
      w  := 1.0/x;
      Iv := w/(Kv*fv + Kv1); {Wronskian relation}
    end;
    if reflect then begin
      {Kv contribution to reflection}
      t := sinPi(v)*Kv/Pi_2;
      {Note: Kv is scaled with exp(|x|), so we have to multiply with exp(-2|x|)}
      if escale then t := t*exp(-2.0*abs(x));
      Iv := Iv + t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_iv(v, x: extended): extended;
  {-Return I_v(x), the modified Bessel function of the 1st kind, order v.}
var
  r,t: extended;
const
  IvMaxXH: THexExtW = ($7D70,$8F6F,$7209,$B188,$400C);  {+1.1362111364594631013E+4}
begin

  if IsNaNorInf(x) or IsNaNorInf(v) then begin
    sfc_iv := NaN_x;
    exit;
  end;

  if x<0.0 then begin
    {if v is not an integer I(v, x) is complex}
    if frac(v)=0.0 then begin
      r := sfc_iv(v,-x);
      if frac(0.5*v)<>0 then r := -r;
      sfc_iv := r;
    end
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfc_iv := NaN_x;
    end;
  end
  else if abs(v)=0.5 then begin
    {NIST[30] 10.39.1: I(0.5,x)=sinh(x)/R, I(-0.5,x)=cosh(x)/R, R=sqrt(Pi*x/2)}
    if x >= Ln_MaxExt then begin
      if x >= extended(IvMaxXH) then sfc_iv := PosInf_x
      else begin
        {Avoid overflow for x in range 11356.52 .. 11362.11}
        r := exp(0.5*x);
        sfc_iv := r*(r/sqrt(TwoPi*x));
      end;
    end
    else begin
      r := sqrt(Pi_2*x);
      if v<0.0 then sfc_iv := cosh(x)/r
      else sfc_iv := sinh(x)/r;
    end;
  end
  else if x=0.0 then begin
    if v=0.0 then sfc_iv := 1.0
    else sfc_iv := 0.0;
  end
  else begin
    {x>0}
    if v=0.0 then sfc_iv := sfc_i0(x)
    else if abs(v)=1.0 then sfc_iv := sfc_i1(x)
    else if x >= extended(IvMaxXH) then sfc_iv := PosInf_x
    else begin
      bessel_ik(v,x,true,false,r,t);
      sfc_iv := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_kv(v, x: extended): extended;
  {-Return K_v(x), the modified Bessel function of the 2nd kind, order v, x>0}
var
  r,t: extended;
begin

  if IsNaNorInf(x) or IsNaNorInf(v) then begin
    sfc_kv := NaN_x;
    exit;
  end;

  if (frac(v)=0.0) and (v<MaxInt) then sfc_kv := sfc_kn(round(v),x)
  else if x<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_kv := NaN_x;
  end
  else begin
    if abs(v)=0.5 then begin
      {NIST[30] 10.39.2: K(0.5,x) = K(-0.5,x) = exp(-x)*sqrt(Pi/2/x)}
      sfc_kv := exp(-x)*sqrt(Pi_2/x);
    end
    else begin
      bessel_ik(v, x, false, false, t, r);
      sfc_kv := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure sfc_bess_kv2(v,x: extended; var Kv, Kv1: extended);
  {-Return K(v,x) and K(v+1,x), no checks, x>0, |v| < MaxLongint}
var
  n,k: longint;
  u, Ku, Ku1: extended;
  curr, prev, next, t, x2: extended;
begin

  {Ref: Boost [19] file bessel_ik.hpp, function bessel_ik}
  {and NR [13], Ch.6.7, section Modified Bessel Functions}
  u := v;
  if v < -0.5 then u := abs(u) - 1.0;
  x2 := 0.5*x;
  n  := round(u);
  u  := u-n;

  if x <= 2.0 then temme_k(u, x, Ku, Ku1)
  else CF2_K(u, x, false, Ku, Ku1);

  prev := Ku;
  curr := Ku1;
  for k:=1 to n do begin
    {forward recurrence for K}
    t := (u+k)/x2;
    if (t > 1.0) and (curr >= MaxExtended/t) then begin
      Kv  := PosInf_x;
      Kv1 := PosInf_x;
      exit;
    end;
    next := t*curr + prev;
    prev := curr;
    curr := next;
  end;
  if v < -0.5 then begin
    Kv1 := prev;
    Kv  := curr;
  end
  else begin
    Kv  := prev;
    Kv1 := curr;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ive(v, x: extended): extended;
  {-Return I_v(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order v.}
var
  r,t: extended;
begin

  if IsNaNorInf(x) or IsNaNorInf(v) then begin
    sfc_ive := NaN_x;
    exit;
  end;

  if x<0.0 then begin
    {if v is not an integer I(v, x) is complex}
    if frac(v)=0.0 then begin
      r := sfc_ive(v,-x);
      if frac(0.5*v)<>0 then r := -r;
      sfc_ive := r;
    end
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfc_ive := NaN_x;
    end;
  end
  else if x=0.0 then begin
    if v=0.0 then sfc_ive := 1.0
    else sfc_ive := 0.0;
  end
  else begin
    {x>0}
    if v=0.0 then sfc_ive := sfc_i0e(x)
    else if abs(v)=1.0 then sfc_ive := sfc_i1e(x)
    else begin
      bessel_ik(v,x,true,true,r,t);
      sfc_ive := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_kve(v, x: extended): extended;
  {-Return K_v(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order v, x>0}
var
  r,t: extended;
begin

  if IsNaNorInf(x) or IsNaNorInf(v) then begin
    sfc_kve := NaN_x;
    exit;
  end;

  if x<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_kve := NaN_x;
  end
  else begin
    v := abs(v);
    if v=0.0 then sfc_kve := sfc_k0e(x)
    else if v=1.0 then sfc_kve := sfc_k1e(x)
    else if abs(v)=0.5 then begin
      {NIST[30] 10.39.1: K(0.5,x) = K(-0.5,x) = exp(-x)*sqrt(Pi/2/x)}
      sfc_kve := sqrt(Pi_2/x);
    end
    else begin
      bessel_ik(v, x, false, true, t, r);
      sfc_kve := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure sfc_bess_ikv(v,x: extended; var Iv,Kv: extended);
  {-Return I_v(x) and K_v(x), no checks, x>0, |v| < MaxLongint}
begin
  bessel_ik(v,x,true,false,Iv,Kv);
end;


{---------------------------------------------------------------------------}
function sfc_sph_jn(n: integer; x: extended): extended;
  {-Return j_n(x), the spherical Bessel function of the 1st kind, order n}
var
  r,z: extended;
begin
  if IsNaNorInf(x) then begin
    sfc_sph_jn := NaN_x;
    exit;
  end;
  if n=0 then sfc_sph_jn := sinc(x)
  else begin
    if x=0.0 then sfc_sph_jn := 0.0
    else begin
      z := abs(x);
      r := sqrt(Pi_2/z)*sfc_jv(0.5+n,z);
      if (x<0.0) and odd(n) then r := -r;    {NIST 10.47.14}
      sfc_sph_jn := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_sph_yn(n: integer; x: extended): extended;
  {-Return y_n(x), the spherical Bessel function of the 2nd kind, order n >=0 , x<>0}
var
  r,z: extended;
begin
  if IsNaNorInf(x) then begin
    sfc_sph_yn := NaN_x;
    exit;
  end;
  z := abs(x);
  if x=0.0 then r := NegInf_x
  else r := sqrt(Pi_2/z)*sfc_yv(0.5+n,z);
  if (x<0.0) and odd(n+1) then r := -r;      {NIST 10.47.14}
  sfc_sph_yn := r;
end;


{---------------------------------------------------------------------------}
function sfc_sph_in(n: integer; x: extended): extended;
  {-Return i_n(x), the modified spherical Bessel function of the 1st/2nd kind, order n}
var
  r: extended;
begin
  if IsNaNorInf(x) then begin
    sfc_sph_in := NaN_x;
    exit;
  end;
  if n=0 then sfc_sph_in := sinhc(x) {i_0 = sinh(x)/x}
  else begin
    if x=0.0 then begin
      if n>0 then sfc_sph_in := 0.0
      else if odd(n) then sfc_sph_in := PosInf_x
      else sfc_sph_in := NegInf_x
    end
    else begin
      r := abs(x);
      r := sqrt(Pi_2/r)*sfc_iv(0.5+n,r);
      if (x<0.0) and odd(n) then r := -r;  {NIST 10.47.16}
      sfc_sph_in := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_sph_ine(n: integer; x: extended): extended;
  {-Return i_n(x)*exp(-|x|), the exponentially scaled modified spherical Bessel function of the 1st/2nd kind, order n}
var
  r,z: extended;
begin
  if IsNaNorInf(x) then begin
    sfc_sph_ine := NaN_x;
    exit;
  end;
  z := abs(x);
  if n=0 then begin
    {i_0e = exp(-z)*sinh(z)/z}
    if z<=1e-10 then begin
      {exp(-z)*sinh(z)/z = 1 - z + 2/3*z^2 - 1/3*z^3 + 2/15*z^4 + O(z^5)}
      sfc_sph_ine := 1.0 - z;
    end
    else begin
      {exp(-z)*sinh(z)/z = -0.5*(exp(-2z)-1)/z}
      if z>=25 then sfc_sph_ine := 0.5/z
      else sfc_sph_ine := (-0.5)*expm1(-2.0*z)/z;
    end;
  end
  else begin
    if x=0.0 then begin
      if n>0 then sfc_sph_ine := 0.0
      else if odd(n) then sfc_sph_ine := PosInf_x
      else sfc_sph_ine := NegInf_x
    end
    else begin
      r := sqrt(Pi_2/z)*sfc_ive(0.5+n,z);
      if (x<0.0) and odd(n) then r := -r;         {NIST 10.47.16}
      sfc_sph_ine := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_sph_kn(n: integer; x: extended): extended;
  {-Return k_n(x), the modified spherical Bessel function of the 3rd kind, order n, x>0}
begin
  if IsNaNorInf(x) then begin
    sfc_sph_kn := NaN_x;
    exit;
  end;
  if x<0.0 then begin
    sfc_sph_kn := NaN_x;
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;
  if n<0 then n := -n-1;     {NIST 10.47.9}
  if n=0 then sfc_sph_kn := Pi_2*exp(-x)/x
  else sfc_sph_kn := sqrt(Pi_2/x)*sfc_kv(0.5+n,x);
end;


{---------------------------------------------------------------------------}
function sfc_sph_kne(n: integer; x: extended): extended;
  {-Return k_n(x)*exp(x), the exponentially scaled modified spherical Bessel function of the 3rd kind, order n, x>0}
begin
  if IsNaNorInf(x) then begin
    sfc_sph_kne := NaN_x;
    exit;
  end;
  if x<0.0 then begin
    sfc_sph_kne := NaN_x;
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;
  if n<0 then n := -n-1;                {NIST 10.47.9}
  if n=0 then sfc_sph_kne := Pi_2/x
  else sfc_sph_kne := sqrt(Pi_2/x)*sfc_kve(0.5+n,x);
end;

{---------------------------------------------------------------------------}
{---------------- Integrals of 0 order Bessel functions --------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
{ Common Chebyshev coefficients for j0_int and y0_int                       }
{---------------------------------------------------------------------------}
const
  nj01  = 24;
  arj01 : array[0..nj01-1] of extended = (
            0.38179279321690173518,
           -0.21275636350505321870,
            0.16754213407215794187,
           -0.12853209772196398954,
            0.10114405455778847013,
           -0.9100795343201568859e-1,
            0.6401345264656873103e-1,
           -0.3066963029926754312e-1,
            0.1030836525325064201e-1,
           -0.255670650399956918e-2,
            0.48832755805798304e-3,
           -0.7424935126036077e-4,
            0.922260563730861e-5,
           -0.95522828307083e-6,
            0.8388355845986e-7,
           -0.633184488858e-8,
            0.41560504221e-9,
           -0.2395529307e-10,
            0.122286885e-11,
           -0.5569711e-13,
            0.227820e-14,
           -0.8417e-16,
            0.282e-17,
           -0.9e-19);
const
  nj0a1 = 22;
  arj0a1: array[0..nj0a1-1] of extended = (
            1.24030133037518970827,
           -0.478125353632280693e-2,
            0.6613148891706678e-4,
           -0.186042740486349e-5,
            0.8362735565080e-7,
           -0.525857036731e-8,
            0.42606363251e-9,
           -0.4211761024e-10,
            0.488946426e-11,
           -0.64834929e-12,
            0.9617234e-13,
           -0.1570367e-13,
            0.278712e-14,
           -0.53222e-15,
            0.10844e-15,
           -0.2342e-16,
            0.533e-17,
           -0.127e-17,
            0.32e-18,
           -0.8e-19,
            0.2e-19,
           -0.1e-19);

  nj0a2 = 19;
  arj0a2: array[0..nj0a2-1] of extended = (
            1.99616096301341675339,
           -0.190379819246668161e-2,
            0.1539710927044226e-4,
           -0.31145088328103e-6,
            0.1110850971321e-7,
           -0.58666787123e-9,
            0.4139926949e-10,
           -0.365398763e-11,
            0.38557568e-12,
           -0.4709800e-13,
            0.650220e-14,
           -0.99624e-15,
            0.16700e-15,
           -0.3028e-16,
            0.589e-17,
           -0.122e-17,
            0.27e-18,
           -0.6e-19,
            0.1e-19);


{---------------------------------------------------------------------------}
function sfc_j0int(u: extended): extended;
  {-Return the integral int(bessel_j0(x), x = 0..u)}
const
  xlow   = 0.806549e-9;   {double: 0.365e-7;}  {~ sqrt(6*eps)}
  xhigh  = 1.84467e19;    {double: 9.008e15;}  {~ 2/eps}
var
  x,t,c,s,f,g: extended;
begin
  {Ref: MacLeod's MISCFUN [22], function J0INT}
  x := abs(u);
  if x < xlow then t := x
  else if x <= 16.0 then begin
    t := sqr(x)/128.0 - 1.0;
    t := x*CSEvalX(t, arj01, nj01);
  end
  else if x <= xhigh then begin
    {x > 16. WE: do not compute x - Pi/4, use  }
    { cos(x-Pi/4) = (cos(x)+sin(x))*sqrt(2)/2, }
    { sin(x-Pi/4) = (sin(x)-cos(x))*sqrt(2)/2  }
    sincos(x,s,c);
    t := 512.0/sqr(x) - 1.0;
    f := CSEvalX(t,arj0a1,nj0a1) / x;
    g := CSEvalX(t,arj0a2,nj0a2);
    t := c*(f+g) + s*(f-g);
    t := 1.0 - t/sqrtpi/sqrt(x);
  end
  else begin
    sfc_j0int := Nan_x;
    exit;
  end;
  if u >= 0.0 then sfc_j0int := t
  else sfc_j0int := -t;
end;


{---------------------------------------------------------------------------}
function sfc_y0int(u: extended): extended;
  {-Return the integral int(bessel_y0(x), x = 0..u), u >= 0}
const
  ny01 = 25;
  ary01: array[0..ny01-1] of extended = (
           0.5449269630272436549,
          -0.1495732358868478216,
           0.11085634486254842337,
          -0.949533001868377711e-1,
           0.6820817786991456963e-1,
          -0.1032465338336820041,
           0.1062570328753442549,
          -0.6258367679961681990e-1,
           0.2385645760338293285e-1,
          -0.644864913015404481e-2,
           0.131287082891002331e-2,
          -0.20988088174989640e-3,
           0.2716042484138347e-4,
          -0.291199114014694e-5,
           0.26344333093795e-6,
          -0.2041172069780e-7,
           0.137124781317e-8,
          -0.8070680792e-10,
           0.419883057e-11,
          -0.19459104e-12,
           0.808782e-14,
          -0.30329e-15,
           0.1032e-16,
          -0.32e-18,
           0.1e-19);
const
  gal2m1 = -1.11593151565841244881;   {Eulergamma-ln(2)-1}
  gamln2 = -0.11593151565841244881;   {Eulergamma-ln(2)}
  twobpi = 0.63661977236758134308;    {2/Pi}
  xlow   = 0.69849e-9;                {double: 0.3161e-7;} {~ sqrt(4.5*eps)}
  xhigh  = 1.84467e19;                {double: 9.008e15;}  {~ 2/eps}
var
  x,t,c,s,f,g: extended;
begin
  {Ref: MacLeod's MISCFUN [22], function Y0INT}
  x := u;
  if (x < 0.0) or (x > xhigh) then begin
    sfc_y0int := Nan_x;
    exit;
  end;
  if x=0.0 then sfc_y0int := 0.0
  else if x < xlow then sfc_y0int := (ln(x) + gal2m1)*twobpi*x
  else if x <= 16.0  then begin
    t := sqr(x)/128.0 - 1.0;
    f := (ln(x) + gamln2) * CSEvalX(t,arj01,nj01);
    g := CSEvalX(t,ary01,ny01);
    sfc_y0int := twobpi*x*(f-g);
  end
  else begin
    {x > 16. WE: do not compute x - Pi/4, use  }
    { cos(x-Pi/4) = (cos(x)+sin(x))*sqrt(2)/2, }
    { sin(x-Pi/4) = (sin(x)-cos(x))*sqrt(2)/2  }
    sincos(x,s,c);
    t := 512.0/sqr(x) - 1.0;
    f := CSEvalX(t,arj0a1,nj0a1) / x;
    g := CSEvalX(t,arj0a2,nj0a2);
    t := c*(f-g) - s*(f+g);
    sfc_y0int := t/sqrtpi/sqrt(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_i0int(u: extended): extended;
  {-Return the integral int(bessel_i0(x), x = 0..u)}
const
  ni01 = 29;
  ari01: array[0..ni01-1] of extended = (
           0.41227906926781516801,
          -0.34336345150081519562,
           0.22667588715751242585,
          -0.12608164718742260032,
           0.6012484628777990271e-1,
          -0.2480120462913358248e-1,
           0.892773389565563897e-2,
          -0.283253729936696605e-2,
           0.79891339041712994e-3,
          -0.20053933660964890e-3,
           0.4416816783014313e-4,
          -0.822377042246068e-5,
           0.120059794219015e-5,
          -0.11350865004889e-6,
           0.69606014466e-9,
           0.180622772836e-8,
          -0.26039481370e-9,
          -0.166188103e-11,
           0.510500232e-11,
          -0.41515879e-12,
          -0.7368138e-13,
           0.1279323e-13,
           0.103247e-14,
          -0.30379e-15,
          -0.1789e-16,
           0.673e-17,
           0.44e-18,
          -0.14e-18,
          -0.1e-19);

const
  ni0a = 34;
  ari0a: array[0..ni0a-1] of extended = (
           2.03739654571143287070,
           0.1917631647503310248e-1,
           0.49923334519288147e-3,
           0.2263187103659815e-4,
           0.158682108285561e-5,
           0.16507855636318e-6,
           0.2385058373640e-7,
           0.392985182304e-8,
           0.46042714199e-9,
          -0.7072558172e-10,
          -0.6747183961e-10,
          -0.2026962001e-10,
          -0.87320338e-12,
           0.175520014e-11,
           0.60383944e-12,
          -0.3977983e-13,
          -0.8049048e-13,
          -0.1158955e-13,
           0.827318e-14,
           0.282290e-14,
          -0.77667e-15,
          -0.48731e-15,
           0.7279e-16,
           0.7873e-16,
          -0.785e-17,
          -0.1281e-16,
           0.121e-17,
           0.214e-17,
          -0.27e-18,
          -0.36e-18,
           0.7e-19,
           0.6e-19,
          -0.2e-19,
          -0.1e-19);
const
  xlow  = 0.11406325e-8;   {double: 0.5161913e-7;}  {~ sqrt(12*eps)}
  xhigh = 11361.885327278; {double: 713.758339516}
var
  x,t,z: extended;
begin
  {Ref: MacLeod's MISCFUN [22], function I0INT}
  x := abs(u);
  if x < xlow then t := x
  else if x <= 18.0 then begin
    t := (3.0*x - 18.0)/(x + 18.0);
    z := CSEvalX(t,ari01,ni01);
    t := x*exp(x)*z;
  end
  else if x <= xhigh then begin
    t := (36.0/x - 0.5) - 0.5;
    z := CSEvalX(t,ari0a,ni0a);
    t := z*(exp(x)/sqrt(TwoPi*x));
    {MISCFUN original: exp(x - 0.5*ln(x) - LnSqrt2Pi + ln(z));}
  end
  else begin
    {Overflow}
    t := PosInf_x;
  end;
  if u >= 0.0 then sfc_i0int := t
  else sfc_i0int := -t;
end;


{---------------------------------------------------------------------------}
function sfc_k0int(u: extended): extended;
  {-Return the integral int(bessel_k0(x), x = 0..u), u >= 0}
const
  nk0in1 = 16;
  ak0in1: array[0..nk0in1-1] of extended = (
            16.79702714464710959477,
             9.79134687676889407070,
             2.80501316044337939300,
             0.45615620531888502068,
             0.4716224457074760784e-1,
             0.335265148269698289e-2,
             0.17335181193874727e-3,
             0.679951889364702e-5,
             0.20900268359924e-6,
             0.516603846976e-8,
             0.10485708331e-9,
             0.177829320e-11,
             0.2556844e-13,
             0.31557e-15,
             0.338e-17,
             0.3e-19);
  nk0in2 = 16;
  ak0in2: array[0..nk0in2-1] of extended = (
            10.76266558227809174077,
             5.62333479849997511550,
             1.43543664879290867158,
             0.21250410143743896043,
             0.2036537393100009554e-1,
             0.136023584095623632e-2,
             0.6675388699209093e-4,
             0.250430035707337e-5,
             0.7406423741728e-7,
             0.176974704314e-8,
             0.3485775254e-10,
             0.57544785e-12,
             0.807481e-14,
             0.9747e-16,
             0.102e-17,
             0.1e-19);
  nk0ina = 28;
  ak0ina: array[0..nk0ina-1] of extended = (
            1.91172065445060453895,
           -0.4183064565769581085e-1,
            0.213352508068147486e-2,
           -0.15859497284504181e-3,
            0.1497624699858351e-4,
           -0.167955955322241e-5,
            0.21495472478804e-6,
           -0.3058356654790e-7,
            0.474946413343e-8,
           -0.79424660432e-9,
            0.14156555325e-9,
           -0.2667825359e-10,
            0.528149717e-11,
           -0.109263199e-11,
            0.23518838e-12,
           -0.5247991e-13,
            0.1210191e-13,
           -0.287632e-14,
            0.70297e-15,
           -0.17631e-15,
            0.4530e-16,
           -0.1190e-16,
            0.319e-17,
           -0.87e-18,
            0.24e-18,
           -0.7e-19,
            0.2e-19,
           -0.1e-19);

const
  xlow   = 0.98781676e-9;           {double: 4.47034836e-8;}  {sqrt(9*eps)}
  xhigh  = 43.6682723752;           {double: 36.0436534;   }  {-ln(eps)}
  rt2bpi = 0.7978845608028653559;   {sqrt(2/Pi)}
  gamln2 = -(0.0625+0.05343151565841244881); {Eulergamma-ln(2)}
  gal2m1 = -(1.0+0.11593151565841244881);    {Eulergamma-ln(2)-1}
var
  x,t,z: extended;
begin
  {Ref: MacLeod's MISCFUN [22], function K0INT}
  x := u;
  if x <= 0.0 then begin
    if x=0.0 then sfc_k0int := 0.0
    else sfc_k0int := Nan_x;
  end
  else if x < xlow then begin
    sfc_k0int := -x*(gal2m1 + ln(x));
  end
  else if x <= 6.0 then begin
    {NOTE: The 16 bit compilers have slightly larger errors, causes}
    {      by suboptimal parsing of the Chebyshev literal arrays !!}
    t := (sqr(x)/18.0 - 0.5) - 0.5;
    z := (gamln2 + ln(x)) * CSEvalX(t,ak0in2, nk0in2);
    sfc_k0int := x*(CSEvalX(t,ak0in1,nk0in1) - z);
  end
  else if x < xhigh then begin
    t := (12.0/x - 0.5) - 0.5;
    z := exp(-x)*CSEvalX(t,ak0ina,nk0ina);
    sfc_k0int := Pi_2 - z /(sqrt(x)*rt2bpi);
  end
  else sfc_k0int := Pi_2;
end;


{---------------------------------------------------------------------------}
function sfc_blam(v,x: extended): extended;
  {-Compute lambda(v,x) = Gamma(v+1)*J(v,x)/(0.5x)^v for v,x >= 0}
var
  z,p: extended;
begin
  if IsNaNorInf(x) or IsNaNorInf(v) or (x<0.0) or (v<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_blam := NaN_x;
    exit;
  end;
  if x<=0.002 then begin
    {use 3 terms of Maclaurin series}
    z := sqr(0.5*x);
    z := z/(v+1.0)*(1-0.5*z/(v+2.0));
    sfc_blam := 1.0 - z;
  end
  else begin
    {http://mathworld.wolfram.com/LambdaFunction.html}
    z := sfc_jv(v,x);
    z := z*sfc_gamma(v+1.0);
    p := power(0.5*x, -v);
    if p=0.0 then begin
      {try 'smaller' power}
      p := power(0.5*x, -0.25*v);
      sfc_blam := (((z*p)*p)*p)*p;
    end
    else sfc_blam := z*p;
  end;
end;



end.
