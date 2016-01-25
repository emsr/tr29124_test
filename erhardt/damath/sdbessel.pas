unit sdBessel;

{Double precision special functions: Bessel functions and related}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}

(*************************************************************************

 DESCRIPTION   :  Double precision Bessel functions and related

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in damath_info.txt/references

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
                      Available as http://oai.cwi.nl/oai/asset/10710/10710A.pdf
                 [52] N.M. Temme, On the Numerical Evaluation of the Modified Bessel
                      Function of the Third Kind, 2nd edition. Preprint, available
                      as http://oai.cwi.nl/oai/asset/7885/7885A.pdf

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  12.02.13  W.Ehrhardt  Initial BP7 version from AMath.sfbessel
 1.00.01  12.02.13  we          Constant IvMaxX in sfd_iv
 1.00.02  12.02.13  we          sfd_i0
 1.00.03  12.02.13  we          sfd_j0, sfd_y0, bess_m0p0
 1.00.04  13.02.13  we          sfd_j1, sfd_y1, bess_m1p1
 1.00.05  13.02.13  we          sfd_in
 1.00.06  13.02.13  we          fix two near overflow cases in bessel_jy
 1.00.07  14.02.13  we          Airy functions
 1.00.08  14.02.13  we          Kelvin functions
 1.00.09  14.02.13  we          Struve functions
 1.00.10  01.03.13  we          Chebyshev degrees reduced in bess_??_small, sfd_??e, ker_kei_med
 1.00.11  02.02.13  we          Fixed value of IvMaxX in sfd_iv

 1.03.00  09.05.13  we          Airy/Scorer functions sfd_airy_gi/hi

 1.04.00  15.06.13  we          Fix typo in h1v_large
 1.04.01  15.06.13  we          sfd_bess_kv2
 1.04.02  28.06.13  we          Check overflow / return INF in bessel_jy

 1.06.00  25.09.13  we          use const one_d

 1.09.00  28.03.14  we          sfd_yn with LnPi from DAMath

 1.13.00  11.08.14  we          Iv := PosInf_d if Kv,Kv1=0 in bessel_ik or if x >= IvMaxX in sfd_iv

 1.18.00  10.06.15  we          new IJ_series replaces Jv_series
 1.18.01  10.06.15  we          rewrite of sfc_in: use IJ_series and CF1_I
 1.18.02  10.06.15  we          improved bess_m0p0/bess_m1p1 with rem_2pi_sym and Kahan summation
 1.18.03  10.06.15  we          sfd_struve_l/h, sfd_struve_l0/1 use sfd_struve_l
 1.18.04  14.06.15  we          sfd_struve_h/l with real v >= 0
 1.18.05  19.06.15  we          scaled Airy functions sfd_airy_ais, sfd_airy_bis

 1.19.00  08.07.15  we          special cases Yv(+- 1/2, x)
 1.19.01  08.07.15  we          sfd_struve_h/l for v < 0
 1.19.02  09.07.15  we          fix sfd_struve_h/l(v,0) for v<=-1

 1.20.00  25.08.15  we          Use THREE & SIXX to avoid 'optimization'
 1.20.01  26.08.15  we          Sqrt2 in kelvin functions (anti-'optimizing')


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


function sfd_i0(x: double): double;
  {-Return I0(x), the modified Bessel function of the 1st kind, order zero}

function sfd_i0e(x: double): double;
  {-Return I0(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order zero}

function sfd_i1(x: double): double;
  {-Return I1(x), the modified Bessel function of the 1st kind, order one}

function sfd_i1e(x: double): double;
  {-Return I1(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order one}

function sfd_in(n: integer; x: double): double;
  {-Return I_n(x), the modified Bessel function of the 1st kind, order n.}

function sfd_j0(x: double): double;
  {-Return J0(x), the Bessel function of the 1st kind, order zero}

function sfd_j1(x: double): double;
  {-Return J1(x), the Bessel function of the 1st kind, order one}

function sfd_jn(n: integer; x: double): double;
  {-Return J_n(x), the Bessel function of the 1st kind, order n; not suitable for large n or x.}

function sfd_k0(x: double): double;
  {-Return K0(x), the modified Bessel function of the 2nd kind, order zero, x>0}

function sfd_k0e(x: double): double;
  {-Return K0(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order zero, x>0}

function sfd_k1(x: double): double;
  {-Return K1(x), the modified Bessel function of the 2nd kind, order one, x>0}

function sfd_k1e(x: double): double;
  {-Return K1(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order one, x>0}

function sfd_kn(n: integer; x: double): double;
  {-Return K_n(x), the modified Bessel function of the 2nd kind, order n, x>0, not suitable for large n}

function sfd_y0(x: double): double;
  {-Return Y0(x), the Bessel function of the 2nd kind, order zero; x>0}

function sfd_y1(x: double): double;
  {-Return Y1(x), the Bessel function of the 2nd kind, order one; x>0}

function sfd_yn(n: integer; x: double): double;
  {-Return Y_n(x), the Bessel function of the 2nd kind, order n, x>0, not suitable for large n or x}


function sfd_jv(v, x: double): double;
  {-Return J_v(x), the Bessel function of the 1st kind, order v; not suitable for large v.}

function sfd_yv(v, x: double): double;
  {-Return Y_v(x), the Bessel function of the 2nd kind, order v; x > 0; not suitable for large v.}

function sfd_iv(v, x: double): double;
  {-Return I_v(x), the modified Bessel function of the 1st kind, order v.}

function sfd_kv(v, x: double): double;
  {-Return K_v(x), the modified Bessel function of the 2nd kind, order v, x>0}

function sfd_ive(v, x: double): double;
  {-Return I_v(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order v.}

function sfd_kve(v, x: double): double;
  {-Return K_v(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order v, x>0}

procedure sfd_bess_ikv(v,x: double; var Iv,Kv: double);
  {-Return I_v(x) and K_v(x), no checks, x>0, |v| < MaxLongint}

procedure sfd_bess_jyv(v,x: double; var Jv,Yv: double);
  {-Return J_v(x) and Y_v(x), no checks, x>0, |v| < MaxLongint}

procedure sfd_bess_kv2(v,x: double; var Kv, Kv1: double);
  {-Return K(v,x) and K(v+1,x), no checks, x>0, |v| < MaxLongint}

function sfd_sph_jn(n: integer; x: double): double;
  {-Return j_n(x), the spherical Bessel function of the 1st kind, order n}

function sfd_sph_yn(n: integer; x: double): double;
  {-Return y_n(x), the spherical Bessel function of the 2nd kind, order n >=0 , x<>0}

function sfd_sph_in(n: integer; x: double): double;
  {-Return i_n(x), the modified spherical Bessel function of the 1st/2nd kind, order n}

function sfd_sph_ine(n: integer; x: double): double;
  {-Return i_n(x)*exp(-|x|), the exponentially scaled modified spherical Bessel function of the 1st/2nd kind, order n}

function sfd_sph_kn(n: integer; x: double): double;
  {-Return k_n(x), the modified spherical Bessel function of the 3rd kind, order n, x>0}

function sfd_sph_kne(n: integer; x: double): double;
  {-Return k_n(x)*exp(x), the exponentially scaled modified spherical Bessel function of the 3rd kind, order n, x>0}


function sfd_airy_ai(x: double): double;
  {-Return the Airy function Ai(x)}

function sfd_airy_aip(x: double): double;
  {-Return the Airy function Ai'(x)}

function sfd_airy_ais(x: double): double;
  {-Return the scaled Airy function Ai(x) if x <= 0, Ai(x)*exp(2/3*x^1.5) for x > 0}

function sfd_airy_bi(x: double): double;
  {-Return the Airy function Bi(x)}

function sfd_airy_bip(x: double): double;
  {-Return the Airy function Bi'(x)}

function sfd_airy_bis(x: double): double;
  {-Return the scaled Airy function Bi(x) if x <= 0, Bi(x)*exp(-2/3*x^1.5) for x > 0}

function sfd_airy_gi(x: double): double;
  {-Return the Airy/Scorer function Gi(x) = 1/Pi*integral(sin(x*t+t^3/3), t=0..INF)}

function sfd_airy_hi(x: double): double;
  {-Return the Airy/Scorer function Hi(x) = 1/Pi*integral(exp(x*t-t^3/3), t=0..INF)}


function sfd_ber(x: double): double;
  {-Return the Kelvin function ber(x)}

function sfd_bei(x: double): double;
  {-Return the Kelvin function bei(x)}

function sfd_ker(x: double): double;
  {-Return the Kelvin function ker(x), x > 0}

function sfd_kei(x: double): double;
  {-Return the Kelvin function kei(x), x >= 0}

procedure sfd_ker_kei(x: double; var kr, ki: double);
  {-Return the Kelvin functions kr=ker(x), ki=kei(x), x > 0}

procedure sfd_ber_bei(x: double; var br, bi: double);
  {-Return the Kelvin functions br=ber(x), bi=bei(x)}


function sfd_struve_h0(x: double): double;
 {-Return H0(x), the Struve function of order 0}

function sfd_struve_h1(x: double): double;
  {-Return H1(x), the Struve function of order 1}

function sfd_struve_h(v,x: double): double;
  {-Return H_v(x), the Struve function of order v, x < 0 only if v is an integer.}

function sfd_struve_l0(x: double): double;
  {-Return L0(x), the modified Struve function of order 0}

function sfd_struve_l1(x: double): double;
  {-Return L1(x), the modified Struve function of order 1}

function sfd_struve_l(v, x: double): double;
  {-Return L_v(x), the modified Struve function of order v, x < 0 only if v is an integer.}


implementation


uses
  DAMath,
  sdBasic,  {Basic common code}
  sdGamma;  {Gamma function and related}


{---------------------------------------------------------------------------}
function IJ_series(v,x: double;  CalcJv: boolean): double;
  {-Power series for Bessel J_v(x), 0 <= v < MAXGAMD-1, |x|<1 or v > x^2/4}
var
  f,s,t: double;
  n: integer;
begin
  f := 0.5*x;
  t := power(f,v)/sfd_gamma(v+1.0);
  if CalcJv then f := -f*f else f := f*f;
  s := t;
  n := 0;
  repeat
    inc(n);
    t := t*f/n/(v+n);
    s := s + t;
  until abs(t) <= 0.5*eps_d*abs(s);
  IJ_series := s;
end;


{---------------------------------------------------------------------------}
procedure CF1_j(v,x: double; var fv: double; var s: integer);
  {-Return J_(v+1)(x) / J_v(x), efficient only if |x| <= |v|}
var
  c,d,f,b,t,tiny,tol: double;
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

  tol  := 2.0*eps_d;
  tiny := Sqrt_MinDbl;
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
procedure CF1_I(v,x: double; var fv: double);
  {-Return I_(v+1)(x) / I_v(x), use only if |x| <= |v|}
var
  c,d,f,b,t,tiny,tol: double;
  k: integer;
const
  MAXIT = 30000;
begin
  {Evaluate NIST[30], 10.33.1 using modified Lentz's method.}
  {Ref: NR [13] (6.7.21) and p.248, function bessik         }
  {and Boost [19] file bessel_ik.hpp, function CF1_Ik       }
  {If |x| <= |v|, CF1_I converges rapidly but if |x| > |v|  }
  {then CF1_I needs O(|x|) iterations to converge!          }

  tol  := 2.0*eps_d;
  tiny := Sqrt_MinDbl;
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
procedure h1v_large(v, x: double; var mv,tmx: double);
  {-Return modulus and (phase - x) of the Hankel function H1_v(x), x > 0 large}
var
  s,m,m2,y: double;
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
function bessj_large(v, x: double): double;
  {-Return J_v(x) via modulus/phase asymptotic expansion, x large}
var
  mv,tv,st,ct,sx,cx: double;
begin
  h1v_large(v,x,mv,tv);
  sincos(tv,st,ct);
  sincos(x,sx,cx);
  {J_v := mv*cos(x+tv); cos(x+tv) = cos(x)cos(tv) - sin(x)sin(tv)}
  bessj_large := mv*(cx*ct - sx*st);
end;


{---------------------------------------------------------------------------}
function bessy_large(v, x: double): double;
  {-Return Y_v(x) via modulus/phase asymptotic expansion, x large}
var
  mv,tv,st,ct,sx,cx: double;
begin
  h1v_large(v,x,mv,tv);
  sincos(tv,st,ct);
  sincos(x,sx,cx);
  {Y_v := mv*sin(x+tv); sin(x+tv) = cos(x)sin(tv) + sin(x)cos(tv)}
  bessy_large := mv*(st*cx + ct*sx);
end;


{---------------------------------------------------------------------------}
procedure bess_m0p0(x: double; var m0,p0: double);
  {-Modulus and phase for J0(x) and Y0(x), x >= 9.0}
var
  y,z,s: double;
const
  m0nhex: array[0..7] of THexDblW = (
             ($6E79,$0EA9,$5470,$3F81),  { 8.461833426898867839659E-3}
             ($98CD,$B6E9,$2244,$3FB7),  { 9.036664453160200052296E-2}
             ($9C1A,$F42E,$EBCC,$3FE6),  { 7.162842530423205720962E-1}
             ($CBD6,$4E89,$BD84,$4006),  { 2.842537511425216145635E0 }
             ($75EC,$377C,$812B,$401E),  { 7.626141421290849630523E0 }
             ($CE76,$CA24,$6D69,$4024),  { 1.021369773577974343844E1 }
             ($CC67,$7234,$756F,$401B),  { 6.864682945702134624126E0 }
             ($AA51,$471E,$43A7,$3FD9)); { 3.947542376069224461532E-1}
  m0dhex: array[0..7] of THexDblW = (
             ($7011,$07C3,$B840,$3F85),  { 1.060533546154121770442E-2}
             ($39BE,$7586,$FE76,$3FBC),  { 1.132577931332212304986E-1}
             ($1BB0,$9BF7,$BFA0,$3FEC),  { 8.983920141407590632423E-1}
             ($6F2F,$B383,$8EAF,$400C),  { 3.569671060989910901903E0 }
             ($0EB7,$7806,$39DB,$4023),  { 9.613002539386213788182E0 }
             ($BE8C,$0916,$0653,$402A),  { 1.301235226061478261481E1 }
             ($A549,$7F67,$3BFE,$4022),  { 9.117176038171821115904E0 }
             ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0 }
  p0nhex: array[0..5] of THexDblW = (
             ($BB0A,$3865,$CBAF,$BE9C),  {-4.290885090773112963542E-7}
             ($4BDC,$75E6,$A6B1,$BF13),  {-7.496317036829964150970E-5}
             ($0CDD,$B531,$6124,$BF70),  {-3.998893155826990642730E-3}
             ($9CF2,$7AB0,$CEBB,$BFB3),  {-7.737323518141516881715E-2}
             ($9B44,$ED20,$0156,$BFE0),  {-5.001635199922493694706E-1}
             ($54CA,$BA42,$8AA9,$BFE7)); {-7.356766355393571519038E-1}
  p0dhex: array[0..6] of THexDblW = (
             ($BB0A,$3865,$CBAF,$3ECC),  { 3.432708072618490390095E-6}
             ($00D6,$E129,$B5B0,$3F43),  { 6.014932317342190404134E-4}
             ($9907,$73FF,$8973,$3FA0),  { 3.229866782185025048457E-2}
             ($52AE,$B8C5,$50A5,$3FE4),  { 6.348446472935245102890E-1}
             ($F0E9,$63E7,$23E2,$4011),  { 4.285043297797736118069E0 }
             ($54BE,$CA56,$82EC,$401D),  { 7.377856408614376072745E0 }
             ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0 }
var
  m0n: array[0..7] of double absolute m0nhex;
  m0d: array[0..7] of double absolute m0dhex;
  p0n: array[0..5] of double absolute p0nhex;
  p0d: array[0..6] of double absolute p0dhex;
begin
  {Ref: Cephes [7], file ldouble\j0l.c}
  {See also HMF[1], sections 9.2.17 .. 9.2.31}

  {Calculate the modulus m0(x) = sqrt(J0(x)^2 + Y0(x)^2) and the}
  {phase p0(x) = arctan(Y0(x)/J0(x)) with rational approximations}
  {For x>=9: J0(x) = m0(x)*cos(p0(x)) and Y0(x) = m0(x)*sin(p0(x))}
  z  := sqr(1.0/x);
  y  := abs(x);
  s  := rem_2pi_sym(y);
  p0 := PolEval(z,p0n,6)/PolEval(z,p0d,7)/y;
  z  := 1.0/y;
  m0 := PolEval(z,m0n,8)/PolEval(z,m0d,8)/sqrt(y);

  {Compute p0 := rem_2pi_sym(y) - Pi_4 + p0 with optimized Kahan}
  {summation. Without rem_p2pi and Kahan we get a relative error}
  {of 1.4e-18 for J0(18), with this code it is 5.4e-20.}

  {Note that this improves only J0/Y0 for arguments x near the  }
  {zeroes. With the recurrence relations for Jn/Yn only absolute}
  {(NOT relative) accuracies of order eps_x can be achieved!    }
  z  := -Pi_4;
  y  := s + z;
  z  := (y - s) - z;
  p0 := y + (p0 - z);
end;


{---------------------------------------------------------------------------}
function sfd_j0(x: double): double;
  {-Return J0(x), the Bessel function of the 1st kind, order zero}
var
  y,z: double;
const
  {Squares of first three roots of J0, calculated with Maple and t_rcalc/xh}
  j1h: THexDblW = ($2BBB,$8046,$21FB,$4017);  {5.78318596294678452117599575847}
  j2h: THexDblW = ($DD6F,$A621,$78A4,$403E);  {30.4712623436620863990781631750}
  j3h: THexDblW = ($5768,$B821,$B8C4,$4052);  {74.8870067906951834448890413101}
var
  jz1: double absolute j1h;
  jz2: double absolute j2h;
  jz3: double absolute j3h;
const
  j0nhex:  array[0..7] of THexDblW = (
             ($90D7,$72B8,$E8FF,$C3A3),  {-7.173386747526788067407E17}
             ($A4C8,$3755,$C86E,$434C),  { 1.620335009643150402368E16}
             ($6623,$5BAC,$EA6F,$C2E0),  {-1.487926133645291056388E14}
             ($36F3,$D94D,$32C6,$4265),  { 7.283696461857171054941E11}
             ($6C00,$232F,$F25A,$C1DE),  {-2.076797068740966813173E9 }
             ($C120,$0536,$A069,$414A),  { 3.490002040733471400107E6 }
             ($C575,$651C,$4E67,$C0A9),  {-3.239201943301299801018E3 }
             ($411F,$7AC8,$BFE4,$3FF4)); { 1.296848754518641770562E0 }
  j0dhex: array[0..8] of THexDblW = (
             ($6560,$EC5F,$096D,$4480),  { 9.466475654163919450528E21}
             ($9960,$DB14,$1705,$4411),  { 7.881340554308432241892E19}
             ($AAA8,$1C3D,$20A4,$4392),  { 3.265560832845194013669E17}
             ($96C1,$48D1,$3C6B,$4309),  { 8.879132373286001289461E14}
             ($0DC0,$84A1,$8146,$4279),  { 1.752689035792859338860E12}
             ($7FCD,$EC52,$6F20,$41E3),  { 2.608400226578100610991E9 }
             ($47FC,$6B92,$3459,$4146),  { 2.910386840401647706984E6 }
             ($0658,$73F9,$D3EB,$40A1),  { 2.281959869176887763845E3 }
             ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0 }
var
  j0n: array[0..7] of double absolute j0nhex;
  j0d: array[0..8] of double absolute j0dhex;
begin
  {Ref: Cephes [7], file ldouble\j0l.c}
  x := abs(x);
  if x < 9.0 then begin
    {In the interval [0,9) a rational approximation of the form }
    {J0(x) = (x^2 - r^2) (x^2 - s^2) (x^2 - t^2) P7(x^2)/Q8(x^2)}
    {is used, where r, s, t are the first three zeros of J0.}
    z := sqr(x);
    y := (z - jz1)*(z - jz2)*(z - jz3);
    sfd_j0 := y*PolEval(z,j0n,8)/PolEval(z,j0d,9);
  end
  else begin
    {For x>=9 the common rational approximations to modulus}
    {and phase are used J0(x) = modulus * cos(phase).}
    if x >= 500.0 then sfd_j0 := bessj_large(0,x)
    else begin
      bess_m0p0(x,y,z);
      sfd_j0 := y*cos(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_y0(x: double): double;
  {-Return Y0(x), the Bessel function of the 2nd kind, order zero; x>0}
var
  y, z: double;
const
  {The first four roots of Y0, calculated with Maple and t_rcalc/xh}
  y1h: THexDblW = ($569C,$4D98,$A953,$400F);  {3.957678419314857868376}
  y2h: THexDblW = ($2103,$C4E7,$581D,$401C);  {7.086051060301772697624}
  y3h: THexDblW = ($7D58,$35A4,$71D7,$4024);  {10.22234504349641701900}
  y4h: THexDblW = ($E74A,$C4A1,$B8E1,$402A);  {13.36109747387276347827}
var
  y0z1: double absolute y1h;
  y0z2: double absolute y2h;
  y0z3: double absolute y3h;
  y0z4: double absolute y4h;

const
  y0nhex: array[0..7] of THexDblW = (
            ($2E2C,$6B40,$6802,$C350),  {-1.847183690384811186958E16}
            ($4C8A,$BE8C,$A2F7,$4363),  { 4.421767595991969611983E16}
            ($6524,$488D,$9928,$C328),  {-3.461898868011666236539E15}
            ($F756,$2846,$E3EA,$42D3),  { 8.747842804834934784972E13}
            ($6868,$D936,$8BE8,$C26C),  {-9.808510181632626683952E11}
            ($64E6,$5309,$3879,$41F4),  { 5.427926320587133391307E9 }
            ($72A6,$2FEF,$EE05,$C16B),  {-1.464324149797947303151E7 }
            ($17C2,$8FE4,$688C,$40CE)); { 1.556909814120445353691E4 }
  y0dhex: array[0..7] of THexDblW = (
            ($E565,$5351,$C96C,$438B),  { 2.502813268068711844040E17}
            ($C521,$C3B4,$821A,$4326),  { 3.167750475899536301562E15}
            ($836D,$B571,$7098,$42B2),  { 2.027480766502742538763E13}
            ($B90B,$D2A9,$1870,$4234),  { 8.630939306572281881328E10}
            ($1F53,$FD32,$0394,$41B0),  { 2.686702051957904669677E8 }
            ($D65B,$3B16,$17CE,$4123),  { 6.256391154086099882302E5 }
            ($E4B2,$36EF,$432B,$4090),  { 1.040792201755841697889E3 }
            ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0 }
  y059nh: array[0..9] of THexDblW = (
            ($A21F,$7D26,$FCE1,$C120),  {-5.566567444353735925323E5 }
            ($5FE0,$4476,$6FC1,$C0FC),  {-1.164760792144532266855E5 }
            ($3C72,$37C5,$E202,$4123),  { 6.515211089266670755622E5 }
            ($6A49,$55FF,$D140,$4107),  { 1.951120419910720443331E5 }
            ($689D,$4729,$8B69,$40BB),  { 7.051411242092171161986E3 }
            ($4B6C,$9678,$E245,$C0B1),  {-4.578271827238477598563E3 }
            ($CE7B,$CB3E,$5189,$4053),  { 7.727403527387097461580E1 }
            ($6B4B,$3F40,$428B,$4039),  { 2.525993724177105060507E1 }
            ($DD42,$7B68,$9118,$BFF7),  {-1.472923738545276751402E0 }
            ($68B3,$16D5,$4172,$3F98)); { 2.368715538373384869796E-2}
  y059dh: array[0..9] of THexDblW = (
            ($4167,$1A3D,$7AF7,$41D1),  { 1.173085288957116938494E9}
            ($8C7A,$BB66,$36AB,$41EA),  { 3.518324187204647941098E9}
            ($D73C,$A14E,$FF64,$C1C4),  {-7.045635226159434678833E8}
            ($5245,$AF19,$D427,$4199),  { 1.083335477747278958468E8}
            ($D2FE,$BBC8,$E85B,$C164),  {-1.096162986826467060921E7}
            ($E0BD,$523A,$C59B,$412A),  { 8.772616606054526158657E5}
            ($AC74,$33EA,$EBDA,$C0E8),  {-5.103881883748705381186E4}
            ($CD56,$A051,$6194,$40A1),  { 2.224790285641017194158E3}
            ($8D1C,$4AB4,$2D71,$C04F),  {-6.235501989189125881723E1}
            ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0}
var
  y059n: array[0..9] of double absolute y059nh;
  y059d: array[0..9] of double absolute y059dh;
  y0n:   array[0..7] of double absolute y0nhex;
  y0d:   array[0..7] of double absolute y0dhex;

begin
  {Ref: Cephes [7], file ldouble\j0l.c}
  if x<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_y0 := NaN_d;
    exit;
  end;
  if x < 9.0 then begin
    z := sqr(x);
    if z < 20.25 then begin
      {In the interval [0,4.5) a rational approximation of the}
      {form Y0(x) = P7(x)/Q7(x) + 2/Pi*ln(x)*J0(x) is used.   }
      y := ln(x)*sfd_j0(x)/Pi_2;
      sfd_y0 := y + PolEval(z,y0n,8)/PolEval(z,y0d,8);
    end
    else begin
      {In the interval [4.5,9) a rational approximation of the}
      {form Y0(x) = (x - p)(x - q)(x - r)(x - s)P9(x)/Q9(x) is}
      {is used where p, q, r, s are first four zeros of Y0(x).}
      y := (x - y0z1)*(x - y0z2)*(x - y0z3)*(x - y0z4);
      sfd_y0 := y * PolEval(x,y059n,10)/PolEval(x,y059d,10);
    end;
  end
  else begin
    {For x>=9 the common rational approximations to modulus}
    {and phase are used Y0(x) = modulus * sin(phase).}
    if x >= 1600 then sfd_y0 := bessy_large(0,x)
    else begin
      bess_m0p0(x,y,z);
      sfd_y0 := y*sin(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure bess_m1p1(x: double; var m1,p1: double);
  {-Modulus and phase for J1(x) and Y1(x), x >= 9.0}
var
  y,z,s: double;
const
  m1nhex: array[0..7] of THexDblW = (
            ($BF75,$2020,$5926,$3F90),  { 1.596507617085714650238E-2}
            ($F5CE,$C4D1,$12FF,$3FC0),  { 1.255798064266130869132E-1}
            ($CAE6,$CE18,$3C17,$3FF0),  { 1.014671139779858141347E0 }
            ($14AE,$4A71,$E041,$4006),  { 2.859499532295180940060E0 }
            ($229D,$410D,$B047,$401C),  { 7.172146812845906480743E0 }
            ($7190,$04BB,$385D,$4004),  { 2.527521168680500659056E0 }
            ($C746,$4C30,$1407,$3FD9),  { 3.918474430130242177355E-1}
            ($B308,$77F6,$2ABE,$C014)); {-5.041742205078442098874E0 }
  m1dhex: array[0..8] of THexDblW = (
            ($8D82,$6AD9,$7D4E,$3F94),  { 2.000925566825407466160E-2}
            ($5900,$F480,$2562,$3FC4),  { 1.573909467558180942219E-1}
            ($35C7,$E2AF,$4985,$3FF4),  { 1.267949948774331531237E0 }
            ($9C4F,$0346,$6F4A,$400C),  { 3.554340386955608261463E0 }
            ($8D60,$060E,$829A,$4021),  { 8.755081357265851765640E0 }
            ($F9CD,$C892,$4111,$4004),  { 2.531772200570435289832E0 }
            ($E5AF,$D7F1,$7C36,$BFED),  {-9.214128701852838347002E-1}
            ($81C5,$B1AC,$EEAF,$C018),  {-6.233092094568239317498E0 }
            ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0 }
  p1nhex: array[0..5] of THexDblW = (
            ($3AAB,$12D8,$A9F6,$3EBC),  { 1.708502235134706284899E-6}
            ($D0FF,$CFB7,$E82F,$3F32),  { 2.884976126715926258586E-4}
            ($9475,$8C1E,$2885,$3F8E),  { 1.472572645054468815027E-2}
            ($173D,$96E5,$2B8F,$3FD1),  { 2.682837461073751055565E-1}
            ($154C,$A013,$65E6,$3FF9),  { 1.587378144541918176658E0 }
            ($A0DD,$25EA,$156A,$4000)); { 2.010456367705144783933E0 }
  p1dhex: array[0..6] of THexDblW = (
            ($7C72,$61E5,$1BF9,$3ED3),  { 4.556005960359216767984E-6}
            ($BF75,$8F2A,$464D,$3F49),  { 7.713202197319040439861E-4}
            ($C213,$86BA,$46A7,$3FA4),  { 3.960155028960712309814E-2}
            ($9198,$8513,$6CAA,$3FE7),  { 7.320149039410806471101E-1}
            ($1A56,$56BC,$130B,$4012),  { 4.518597941618813112665E0 }
            ($259E,$A865,$3CC1,$401B),  { 6.809332495854873089362E0 }
            ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0 }
var
  m1n: array[0..7] of double absolute m1nhex;
  m1d: array[0..8] of double absolute m1dhex;
  p1n: array[0..5] of double absolute p1nhex;
  p1d: array[0..6] of double absolute p1dhex;
begin
  {Ref: Cephes [7], file ldouble\j1l.c}

  {Calculate the modulus m1(x) = sign(x)*sqrt(J1(x)^2 + Y1(x)^2) and }
  {the phase p1(x) = arctan(Y1(x)/J1(x)) with rational approximations}
  {For x>=9: J1(x) = m1(x)*cos(p1(x)) and Y1(x) = m1(x)*sin(p1(x))}
  z  := sqr(1.0/x);
  y  := abs(x);
  s  := rem_2pi_sym(y);
  p1 := PolEval(z,p1n,6)/PolEval(z,p1d,7)/y;
  z  := 1.0/y;
  m1 := PolEval(z,m1n,8)/PolEval(z,m1d,9)/sqrt(y);
  if x<0.0 then m1 := -m1;

  {Compute p1 := rem_2pi_sym(y) - 3Pi_4 + p1 with optimized Kahan summation}
  z  := -3.0*Pi_4;
  y  := s + z;
  z  := (y - s) - z;
  p1 := y + (p1 - z);
end;


{---------------------------------------------------------------------------}
function sfd_j1(x: double): double;
  {-Return J1(x), the Bessel function of the 1st kind, order one}
var
  y,z: double;
const
  {Squares of first three roots of J1, calculated with Maple and t_rcalc/xh}
  j1h: THexDblW = ($822C,$4189,$5D2B,$402D);  {14.6819706421238932572197777686}
  j2h: THexDblW = ($A432,$6072,$9BF6,$4048);  {49.2184563216946036702670828464}
  j3h: THexDblW = ($5E2C,$0D78,$DFF7,$4059);  {103.499453895136580332223632536}
var
  jz1: double absolute j1h;
  jz2: double absolute j2h;
  jz3: double absolute j3h;
const
  j1nhex:  array[0..8] of THexDblW = (
             ($623B,$FA4E,$F494,$C392),  {-3.41470097444474566748E17 }
             ($1E1A,$0A5F,$289F,$4338),  { 6.80006297997263446982E15 }
             ($A52C,$EE29,$D773,$C2C9),  {-5.68263073022183470933E13 }
             ($A460,$FF57,$9394,$424E),  { 2.626500686552841932403E11}
             ($04A3,$9399,$16C7,$C1C6),  {-7.41183271195454042842E8  }
             ($E567,$4B9E,$2441,$4134),  { 1.32000129539331214495E6  }
             ($2D65,$A989,$DB34,$C096),  {-1.46280142797793933909E3  }
             ($0A0E,$12D9,$CD74,$3FED),  { 9.31329762279632791262E-1 }
             ($199F,$1643,$444A,$BF31)); {-2.63469779622127762897E-4 }
  j1dhex: array[0..8] of THexDblW = (
             ($B67C,$A45E,$A1E0,$44A5),  { 5.10779045516141578461E22}
             ($2F0D,$85D6,$5EE4,$4433),  { 3.57325874689695599524E20}
             ($E4BD,$AAD4,$1EE7,$43B1),  { 1.23367806884831151194E18}
             ($D800,$5DEA,$ADDE,$4323),  { 2.76959756375961607085E15}
             ($BAB7,$13DF,$41C4,$4290),  { 4.46866213886267829490E12}
             ($3FCC,$1219,$066D,$41F4),  { 5.37544732957807543920E9 }
             ($2C4D,$D0DC,$4309,$4152),  { 4.78723926343829674773E6 }
             ($8DD7,$EA02,$115B,$40A7),  { 2.95267951972943745733E3 }
             ($0000,$0000,$0000,$3FF0)); { 1.00000000000000000000E0 }
var
  j1n: array[0..8] of double absolute j1nhex;
  j1d: array[0..8] of double absolute j1dhex;
begin
  {Ref: Cephes [7], file ldouble\j1l.c}
  z := abs(x);
  if z < 9.0 then begin
    z := sqr(x);
    {In the interval [0,9) a rational approximation of the form }
    {J1(x) = x*(x^2 - r^2)*(x^2 - s^2)*(x^2 - t^2)*P8(x^2)/Q8(x^2)}
    {is used, where r, s, t are the first three zeros of J1.}
    y := x*(z - jz1)*(z - jz2)*(z - jz3);
    sfd_j1 := y*PolEval(z,j1n,9)/PolEval(z,j1d,9);
  end
  else begin
    {For x>=9 the common rational approximations to modulus}
    {and phase are used J1(x) = modulus * cos(phase).}
    if z >= 500.0 then begin
      y := bessj_large(1,z);
      if x<0.0 then sfd_j1 := -y else sfd_j1 := y;
    end
    else begin
      bess_m1p1(x,y,z);
      sfd_j1 := y*cos(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_y1(x: double): double;
  {-Return Y1(x), the Bessel function of the 2nd kind, order one; x>0}
var
  y, z: double;
const
  {The first four roots of Y1, calculated with Maple and t_rcalc/xh}

  y1h: THexDblW = ($F243,$D4DF,$93BE,$4001);  {2.19714132603101703515}
  y2h: THexDblW = ($B02E,$4E87,$B7FE,$4015);  {5.42968104079413513277}
  y3h: THexDblW = ($69B4,$AE61,$3127,$4021);  {8.59600586833116892643}
  y4h: THexDblW = ($3206,$38D4,$7F91,$4027);  {11.7491548308398812434}
var
  y1z1: double absolute y1h;
  y1z2: double absolute y2h;
  y1z3: double absolute y3h;
  y1z4: double absolute y4h;
const
  y1nhex: array[0..6] of THexDblW = (
            ($0907,$2601,$2CAB,$C353),  {-2.158855258453711703120E16}
            ($EAF0,$B7A8,$32CE,$4334),  { 5.685362960165615942886E15}
            ($484D,$C5A2,$CEBE,$C2EB),  {-2.445982226888344140154E14}
            ($4D6B,$3CB5,$C58F,$428C),  { 3.954354656937677136266E12}
            ($D0EF,$E154,$1255,$C21B),  {-2.906793378120403577274E10}
            ($298A,$DF43,$A337,$4197),  { 9.914315981558815369372E7 }
            ($FF0B,$AFDE,$77A1,$C0FF)); {-1.288901054372751879531E5 }
  y1dhex: array[0..7] of THexDblW = (
            ($236E,$CFD7,$733B,$4378),  { 1.101136026928555260168E17}
            ($2FD2,$2C67,$AF41,$4315),  { 1.525917240904692387994E15}
            ($8575,$1306,$94BA,$42A3),  { 1.076474894829072923244E13}
            ($B24B,$1523,$B332,$4227),  { 5.089532584184822833416E10}
            ($2F62,$DAC2,$2946,$41A5),  { 1.775133253792677466651E8 }
            ($3DE5,$C605,$9040,$411C),  { 4.679841933793707979659E5 }
            ($71DC,$6E67,$E515,$408B),  { 8.926354644853231136073E2 }
            ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0 }
  y159nh: array[0..9] of THexDblW = (
            ($60AA,$7B1E,$FED8,$4112),  { 3.112221202330688509818E5 }
            ($2D1F,$8CDB,$A2B1,$C131),  {-1.155761550219364178627E6 }
            ($7B47,$27F4,$8488,$C13B),  {-1.803400156074242435454E6 }
            ($EA1A,$2305,$2E42,$C13B),  {-1.781314136808997406109E6 }
            ($4534,$AC30,$3F82,$C0FE),  {-1.238961670382216747944E5 }
            ($52D6,$B5E8,$7A2C,$40D4),  { 2.096869860275353982829E4 }
            ($4D7C,$C6E4,$1B2F,$40B8),  { 6.171186628598134035237E3 }
            ($90DD,$D7B0,$D861,$C08C),  {-9.230477746767243316014E2 }
            ($3CC0,$7D68,$8889,$4045),  { 4.306669585790359450532E1 }
            ($6857,$CCF7,$C7FE,$BFE5)); {-6.806634906054210550896E-1}
  y159dh: array[0..10] of THexDblW = (
            ($5614,$2A4D,$B64F,$4199),  { 1.078445545755236785692E8}
            ($8DB6,$E0E3,$7501,$C1D1),  {-1.171523459555524458808E9}
            ($2912,$567C,$E88B,$41EF),  { 4.282669747880013349981E9}
            ($4E75,$81E0,$6B17,$C1C8),  {-8.193431077523942651173E8}
            ($F2EE,$4396,$CD7D,$419C),  { 1.208072488974110742912E8}
            ($2E5D,$0C68,$9411,$C166),  {-1.183757638771741974521E7}
            ($88C8,$005B,$20F7,$412C),  { 9.217235006983512475118E5}
            ($5AB8,$B427,$83A5,$C0E9),  {-5.225317824142187494326E4}
            ($90FA,$3804,$7C60,$40A1),  { 2.238187927382180589099E3}
            ($9263,$2543,$E84C,$C04E),  {-6.181482377814679766978E1}
            ($0000,$0000,$0000,$3FF0)); { 1.000000000000000000000E0}
var
  y159n: array[0..9]  of double absolute y159nh;
  y159d: array[0..10] of double absolute y159dh;
  y1n:   array[0..6]  of double absolute y1nhex;
  y1d:   array[0..7]  of double absolute y1dhex;
begin
  {Ref: Cephes [7], file ldouble\j1l.c}
  if x<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_y1 := NaN_d;
    exit;
  end;
  if x < 9.0 then begin
    z := sqr(x);
    if z < 20.25 then begin
      {In the interval [0,4.5) a rational approximation of the form}
      {Y1(x) = x*P6(x)/Q7(x) + 2/Pi*(ln(x)*J1(x) - 1/x) is used.}
      y := (ln(x)*sfd_j1(x) - 1.0/x)/Pi_2;
      sfd_y1 := y + x*PolEval(z,y1n,7)/PolEval(z,y1d,8);
    end
    else begin
      {In the interval [4.5,9) a rational approximation of the form}
      {Y1(x) = (x - p)*(x - q)*(x - r)*(x - s)*P9(x)/Q10(x) is used}
      {where p, q, r, s are first four zeros of Y1(x).}
      y := (x - y1z1)*(x - y1z2)*(x - y1z3)*(x - y1z4);
      sfd_y1 := y * PolEval(x,y159n,10)/PolEval(x,y159d,11);
    end;
  end
  else begin
    {For x>=9 the common rational approximations to modulus}
    {and phase are used Y1(x) = modulus * sin(phase).}
    if x >= 1600 then sfd_y1 := bessy_large(1,x)
    else begin
      bess_m1p1(x,y,z);
      sfd_y1 := y*sin(z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jn(n: integer; x: double): double;
  {-Return J_n(x), the Bessel function of the 1st kind, order n; not suitable for large n or x.}
var
  curr,prev,q,temp,init,xh: double;
  k: integer;
  neg: boolean;
const
  small = 1.703183936E-108;  {~ cbrt(succd(0))}
  lnsml = ln_succd0;
begin
  {Based on boost_1_42_0\boost\math\special_functions\detail\bessel_jn.hpp [19]}
  {Copyright 2006 Xiaogang Zhang, see 3rdparty.ama for Boost license}
  init := Sqrt_MinDbl;
  {Flag to negate result for |n|}
  neg := (n<0) and odd(n);
  n := abs(n);

  if n=0 then curr := sfd_j0(x)
  else if n=1 then curr := sfd_j1(x)
  else if abs(x) <= small then begin
    if (x=0.0) or (n>2) then curr := 0.0
    else curr := 0.125*sqr(x);
  end
  else begin
    xh := 0.5*x;
    if abs(x) > n then begin
      {forward recurrence}
      prev := sfd_j0(x);
      curr := sfd_j1(x);
      for k:=1 to n-1 do begin
        temp := curr*k/xh - prev;
        prev := curr;
        curr := temp;
      end;
    end
    else begin
      {Quick check if |J_n(x)| < succd(0) from HMF[1] 9.1.63}
      {solution of z*exp(sqrt(1-z^2))=1 is z = 0.39989.. }
      q := abs(x/n);
      if n<=50 then temp := 1e-6
      else if n<180 then temp := 1e-2
      else if n<1000 then temp := 0.2
      else temp := 0.3999;
      if q < temp then begin
        {Jn(x) <= [q*exp(sqrt(1-q^2))/(1+sqrt(1-q^2))]^n}
        temp := sqrt(1.0 - q*q);
        temp := ln(q/(1.0+temp)) + temp;
        if temp < lnsml/n then begin
          sfd_jn := 0.0;
          exit;
        end;
      end;
      {set overflow threshold for iteration}
      q := 0.5*MaxDouble*q;
      {backward recurrence}
      CF1_j(n,x,prev,k);
      prev := prev*init;
      curr := init;
      for k:=n downto 1 do begin
        if abs(curr) > q then begin
          {prevent overflow and set result to zero}
          sfd_jn := 0.0;
          exit;
        end;
        temp := curr*k/xh - prev;
        prev := curr;
        curr := temp;
      end;
      curr := (init/curr)*sfd_j0(x);
    end;
  end;
  if neg then sfd_jn := -curr else sfd_jn := curr;
end;


{---------------------------------------------------------------------------}
function sfd_yn(n: integer; x: double): double;
  {-Return Y_n(x), the Bessel function of the 2nd kind, order n, x>0, not suitable for large n or x}
var
  yn,yn1,t: double;
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
    sfd_yn := NaN_d;
    exit;
  end;
  if n=0 then yn := sfd_y0(x)
  else if n=1 then yn := sfd_y1(x)
  else if (n>MAXGAMD) and (x<=2.0) then yn := NegInf_d
  else begin
    if x<1.0 then begin
      {NIST[30] 10.7.4}
      t := sfd_lnfac(n-1) - n*ln(0.5*x) - lnpi;
      if t > ln_MaxDbl then begin
        if neg then sfd_yn := PosInf_d else sfd_yn := NegInf_d;
        exit;
      end;
    end;
    {forward recurrence}
    yn1 := sfd_y0(x);
    yn  := sfd_y1(x);
    x := 0.5*x;
    for k:=1 to n-1 do begin
      t  := yn*k/x - yn1;
      yn1:= yn;
      yn := t;
    end;
  end;
  if neg then sfd_yn := -yn else sfd_yn := yn;
end;


{---------------------------------------------------------------------------}
function bess_i0_small(x: double): double;
  {-Return Bessel function I0(x) for abs(x)<=3, x assumed >= 0}
const
  xsml = 0.298023223876953125e-7;  {sqrt(4*eps_d)}
const
  nbi0 = 12;
  bi0h : array[0..nbi0-1] of THexDblW = (
           ($3DC8,$8F34,$9C6A,$BFB3),  {-0.7660547252839144951081894976243285e-1 }
           ($E329,$528B,$D660,$3FFE),  {+0.1927337953993808269952408750881196e+1 }
           ($F0D7,$1075,$37C5,$3FCD),  {+0.2282644586920301338937029292330415e+0 }
           ($6D77,$ADF2,$B963,$3F8A),  {+0.1304891466707290428079334210691888e-1 }
           ($35CC,$24F7,$787A,$3F3C),  {+0.4344270900816487451378682681026107e-3 }
           ($093E,$3617,$C2C0,$3EE3),  {+0.9422657686001934663923171744118766e-5 }
           ($EADA,$458D,$3F35,$3E83),  {+0.1434006289510691079962091878179957e-6 }
           ($9403,$8466,$B9C8,$3E1B),  {+0.1613849069661749069915419719994611e-8 }
           ($45DB,$D6A4,$B670,$3DAE),  {+0.1396650044535669699495092708142522e-10}
           ($278A,$F187,$F6B8,$3D3A),  {+0.9579451725505445344627523171893333e-13}
           ($5237,$DB96,$37BA,$3CC3),  {+0.5333981859862502131015107744000000e-15}
           ($9483,$AE63,$AD7A,$3C46)); {+0.2458716088437470774696785919999999e-17}
         { ($DEA3,$25A0,$83F7,$3BC6);} {+0.9535680890248770026944341333333333e-20}
var
  bi0: array[0..nbi0-1] of double absolute bi0h;
begin
  {Ref: W. Fullerton [14] and [20], files dbesi0.f and dbsi0e.f}
  {Hex Chebyshev values calculated with mp_arith/t_rcalc}
  if x<=xsml then bess_i0_small := 1.0
  else bess_i0_small := 2.75 + CSEvalD(x*x/4.5-1.0,bi0,nbi0);
end;


{---------------------------------------------------------------------------}
function sfd_i0e(x: double): double;
  {-Return I0(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order zero}
const
  xsml = 0.23283064365386962891e-9;  {sqrt(0.5*eps_d)}
const
  nai0 = 24;
   ai0: array[0..nai0-1] of double = (
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
          +0.2150380106876119887812051287845e-17);
        { -0.3754858341830874429151584452608e-18,
          +0.2354065842226992576900757105322e-19,
          +0.1114667612047928530226373355110e-19,
          -0.5398891884396990378696779322709e-20);}
const
  nai2 = 27;
  ai02: array[0..nai2-1] of double = (
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
          +0.1193650890845982085504399499242e-17);
        { -0.2488709837150807235720544916602e-18,
          -0.1938426454160905928984697811326e-18,
          +0.6444656697373443868783019493949e-19,
          +0.2886051596289224326481713830734e-19,
          -0.1601954907174971807061671562007e-19,
          -0.3270815010592314720891935674859e-20);}
begin
  {Ref: W. Fullerton [14] and [20], file dbsi0e.f}
  x := abs(x);
  if x<=3.0 then begin
    {Note that there is bug in dbsi0e.f from [20] for small x. We use the}
    {Taylor series for I(0,x)*exp(-x) = 1 - x + 3/4*x^2 -5/12*x^3 + O(x^4)}
    if x<=xsml then sfd_i0e := 1 - x
    else sfd_i0e := exp(-x)*bess_i0_small(x);
  end
  else if x<=8.0 then begin
    sfd_i0e := (0.375 + CSEvalD((48.0/x-11.0)/5.0, ai0, nai0))/sqrt(x);
  end
  else begin
    sfd_i0e := (0.375 + CSEvalD(16.0/x-1.0, ai02, nai2))/sqrt(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_i0(x: double): double;
  {-Return I0(x), the modified Bessel function of the 1st kind, order zero}
begin
  x := abs(x);
  if x<=3.0 then sfd_i0 := bess_i0_small(x)
  else if x>ln_MaxDbl then sfd_i0 := PosInf_d
  else sfd_i0 := sfd_i0e(x)*exp(x);
end;


{---------------------------------------------------------------------------}
function bess_i1_small(x: double): double;
  {-Return Bessel function I1(x) for abs(x)<=3}
var
  y: double;
const
  xsml = 0.105367121277235e-7;  {sqrt(0.5*eps_d)}
const
  nbi1 = 11;
  bi1: array[0..nbi1-1] of double = (
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
         +0.24041144040745181799863172032000e-16);
       { +0.10171505007093713649121100799999e-18);}
begin
  {Ref: W. Fullerton [14] and [20], files dbesi1.f and dbsi1e.f}
  y := abs(x);
  if y=0.0 then bess_i1_small := 0.0
  else if y<=xsml then bess_i1_small := 0.5*x
  else bess_i1_small := x*(0.875 + CSEvalD(x*x/4.5-1.0,bi1,nbi1));
end;


{---------------------------------------------------------------------------}
function sfd_i1e(x: double): double;
  {-Return I1(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order one}
var
  y: double;
const
  nai1 = 24;
   ai1: array[0..nai1-1] of double = (
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
          -0.2194016207410736898156266643783e-17);
        { +0.3630516170029654848279860932334e-18,
          -0.1695139772439104166306866790399e-19,
          -0.1288184829897907807116882538222e-19,
          +0.5694428604967052780109991073109e-20);}
const
  nai2 = 27;
  ai12: array[0..nai2-1] of double = (
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
         -0.1242193275194890956116784488697e-17);
      {  +0.2414276719454848469005153902176e-18,
         +0.2026944384053285178971922860692e-18,
         -0.6394267188269097787043919886811e-19,
         -0.3049812452373095896084884503571e-19,
         +0.1612841851651480225134622307691e-19,
         +0.3560913964309925054510270904620e-20);}
begin
  {Ref: W. Fullerton [14] and [20], file dbsi1e.f}
  y := abs(x);
  if y<=3.0 then sfd_i1e := exp(-y)*bess_i1_small(x)
  else begin
    if y<=8.0 then begin
      y := (0.375 + CSEvalD((48.0/y-11.0)/5.0, ai1, nai1))/sqrt(y)
    end
    else begin
      y := (0.375 + CSEvalD(16.0/y-1.0, ai12, nai2))/sqrt(y);
    end;
    if x>0 then sfd_i1e := y else sfd_i1e := -y;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_i1(x: double): double;
  {-Return I1(x), the modified Bessel function of the 1st kind, order one}
var
  y: double;
begin
  y := abs(x);
  if y<=3.0 then sfd_i1 := bess_i1_small(x)
  else if x>ln_MaxDbl then sfd_i1 := PosInf_d
  else sfd_i1 := sfd_i1e(x)*exp(y);
end;


{---------------------------------------------------------------------------}
function bess_k0_small(x: double): double;
  {-Return Bessel function K0(x) for 0 < x <= 2}
var
  y: double;
const
  xsml = 0.21073424255447e-7;  {sqrt(2*eps_d)}
const
  nbk0 = 10;
  bk0: array[0..nbk0-1] of double = (
         -0.353273932339027687201140060063153e-1,
         +0.344289899924628486886344927529213e+0,
         +0.359799365153615016265721303687231e-1,
         +0.126461541144692592338479508673447e-2,
         +0.228621210311945178608269830297585e-4,
         +0.253479107902614945730790013428354e-6,
         +0.190451637722020885897214059381366e-8,
         +0.103496952576336245851008317853089e-10,
         +0.425981614279108257652445327170133e-13,
         +0.137446543588075089694238325440000e-15);
       { +0.357089652850837359099688597333333e-18,
         +0.763164366011643737667498666666666e-21);}
begin
  {Ref: W. Fullerton [14] and [20], files dbesk0.f and dbsk0e.f}
  if x>xsml then y := x*x
  else begin
    if x=0.0 then begin
      bess_k0_small := PosInf_d;
      exit;
    end
    else y := 0.0;
  end;
  bess_k0_small := -ln(0.5*x)*bess_i0_small(x) - 0.25 + CSEvalD(0.5*y-1.0,bk0,nbk0);
end;


{---------------------------------------------------------------------------}
function sfd_k0e(x: double): double;
  {-Return K0(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order zero, x>0}
const
  nak0 = 19;
   ak0: array[0..nak0-1] of double = (
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
          +0.1446404347862212227887763442346e-17);
        { -0.2452889825500129682404678751573e-18,
          +0.4233754526232171572821706342400e-19,
          -0.7427946526454464195695341294933e-20);}
const
  nak2 = 16;
  ak02: array[0..nak2-1] of double = (
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
          -0.3530507631161807945815482463573e-18);
        { +0.4645295422935108267424216337066e-19,
          -0.6368625941344266473922053461333e-20);}
begin
  {Ref: W. Fullerton [14] and [20], file dbsk0e.f}
  if x<=2.0 then sfd_k0e := exp(x)*bess_k0_small(x)
  else if x<=8.0 then begin
    sfd_k0e := (1.25 + CSEvalD((16.0/x-5.0)/THREE, ak0, nak0))/sqrt(x);
  end
  else begin
    sfd_k0e := (1.25 + CSEvalD(16.0/x-1.0, ak02, nak2))/sqrt(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_k0(x: double): double;
  {-Return K0(x), the modified Bessel function of the 2nd kind, order zero, x>0}
begin
  if x<=2.0 then sfd_k0 := bess_k0_small(x)
  else sfd_k0 := sfd_k0e(x)*exp(-x);
end;


{---------------------------------------------------------------------------}
function bess_k1_small(x: double): double;
  {-Return Bessel function K1(x) for 0 < x <= 2}
var
  y: double;
const
  xsml = 0.46566128730773925781e-9;  {sqrt(2*eps_d)}
const
  nbk1 = 11;
  bk1: array[0..nbk1-1] of double = (
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
         -0.70238634793862875971783797120000e-17);
       { -0.16543275155100994675491029333333e-19);}
begin
  {Ref: W. Fullerton [14] and [20], files dbesk1.f and dbsk1e.f}
  if x>xsml then y := x*x
  else begin
    if x=0.0 then begin
      bess_k1_small := PosInf_d;
      exit;
    end
    else y := 0.0;
  end;
  bess_k1_small := ln(0.5*x)*bess_i1_small(x) + (0.75 + CSEvalD(0.5*y-1.0,bk1,nbk1))/x;
end;


{---------------------------------------------------------------------------}
function sfd_k1e(x: double): double;
  {-Return K1(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order one, x>0}
const
  nak1 = 19;
   ak1: array[0..nak1-1] of double = (
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
          -0.16284168738284380035666620115626e-17);
        { +0.27501536496752623718284120337066e-18,
          -0.47289666463953250924281069568000e-19,
          +0.82681500028109932722392050346666e-20);}
const
  nak2 = 16;
  ak12: array[0..nak2-1] of double = (
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
          +0.3851487721280491597094896844799e-18);
        { -0.5044794897641528977117282508800e-19,
          +0.6888673850418544237018292223999e-20);}
begin
  {Ref: W. Fullerton [14] and [20], file dbsk1e.f}
  if x<=2.0 then sfd_k1e := exp(x)*bess_k1_small(x)
  else if x<=8.0 then begin
    sfd_k1e := (1.25 + CSEvalD((16.0/x-5.0)/THREE, ak1, nak1))/sqrt(x);
  end
  else begin
    sfd_k1e := (1.25 + CSEvalD(16.0/x-1.0, ak12, nak2))/sqrt(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_k1(x: double): double;
  {-Return K1(x), the modified Bessel function of the 2nd kind, order one, x>0}
begin
  if x<=2.0 then sfd_k1 := bess_k1_small(x)
  else sfd_k1 := sfd_k1e(x)*exp(-x);
end;


{---------------------------------------------------------------------------}
function sfd_in(n: integer; x: double): double;
  {-Return I_n(x), the modified Bessel function of the 1st kind, order n}
var
  curr,prev,temp,init,y: double;
  k: integer;
const
  NMax = 160;
  XMax = 512.0;
begin
  {HMF[1], 9.6.6: I(-n,x) = I(n,x)}
  n := abs(n);
  if n=0 then sfd_in := sfd_i0(x)
  else if n=1 then sfd_in := sfd_i1(x)
  else if x=0.0 then sfd_in := 0.0
  else begin
    y := abs(x);
    {If n or x are not small then use real order function}
    if (n>NMax) or (y>XMax) then sfd_in := sfd_iv(n,x)
    else begin
      {Simplified calculation for I_n only with existing functions}
      if y <= 2.0 then begin
        curr := IJ_series(n,y,false);
      end
      else begin
        {I_(n+1)(y)/I_n(y) using continued fraction and recurse backwards}
        CF1_I(n,y,temp);
        y := 0.5*y;
        init := Sqrt_MinDbl;
        curr := init;
        prev := init*temp;
        for k:=n downto 1 do begin
          temp := curr*k/y + prev;
          prev := curr;
          curr := temp;
        end;
        curr := (init/curr)*sfd_i0(x);
      end;
      {Adjust sign}
      if (x<0.0) and odd(n) then curr := -curr;
      sfd_in := curr;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_kn(n: integer; x: double): double;
  {-Return K_n(x), the modified Bessel function of the 2nd kind, order n, x>0, not suitable for large n}
var
  kn,knm1,knm2: double;
  k: integer;
begin
  {HMF[1], 9.6.6: K(-n,x) = K(n,x)}
  n := abs(n);
  {Range error for x<=0 is generated in k0 or k1}
  if n=0 then kn := sfd_k0(x)
  else if n=1 then kn := sfd_k1(x)
  else begin
    {$ifdef Delphi}
      {avoid false warning "Variable 'kn' might not have been initialized"}
      kn := 0.0;
    {$endif}
    {forward recurrence, K(n+1,x) = 2n/x*K(n,x) + K(n-1,x)}
    knm2 := sfd_k0(x);
    knm1 := sfd_k1(x);
    x := 0.5*x;
    for k:=1 to n-1 do begin
      kn   := knm1*k/x + knm2;
      knm2 := knm1;
      knm1 := kn;
    end;
  end;
  sfd_kn := kn;
end;


{---------------------------------------------------------------------------}
{-------------------- Bessel functions of real order -----------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure temme_y(v,x: double; var Y,Y1: double);
  {-Calculate Y(v, x) and Y(v+1, x) by Temme's method for small |x|}
var
  k,g,h,p,q,f,c,d,s,s1,tol,a,e: double;
  g1,g2,gp,gm,v2: double;
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
  gp := sfd_gamma1pm1(v);
  gm := sfd_gamma1pm1(-v);
  a  := ln(0.5*x);
  s  := -a*v;

  if abs(v) < eps_d then begin
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
  tol := 0.5*eps_d;
  {use double k because otherwise k*k may overflow}
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
procedure CF2_jy(v,x: double; var p,q: double);
  {-Return the continued fraction p + i*q = (J' + iY') / (J + iY)}
var
  a,br,bi,cr,ci,dr,di,er,ei,fr,fi,t: double;
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
    if abs(dr)+abs(di) = 0.0 then dr := Sqrt_MinDbl;
    t  := a/(cr*cr+ci*ci);
    cr := br+t*cr;
    ci := bi-t*ci;
    if abs(cr)+abs(ci) = 0.0 then cr := Sqrt_MinDbl;
    t  := dr*dr+di*di;
    dr := +dr/t;
    di := -di/t;
    er := cr*dr-ci*di;
    ei := cr*di+ci*dr;
    t  := fr*er-fi*ei;
    fi := fr*ei+fi*er;
    fr := t;
    if abs(er-1.0)+abs(ei) < 8*eps_d then begin
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
function Yv_series(v,x: double): double;
  {-Series for Yv for 'small' x < 1;  frac(v)<>0}
var
  h,k,v1,v2,t1,t2,xx,Yv: double;
const
  hsmall = 0.177065751663e-308;  { ~ 1/(MaxDouble*Pi)}
begin

  {Use Yv(x) = ( Jv(v,x)cos(vx) - Jv(-v,x) )/sin(vx) }
  {and Gamma reflection to sum the two series for Jv }

  if v<MAXGAMD then h := power(0.5*x, v)/sfd_gamma(v)
  else h := exp(v*ln(0.5*x) - sfd_lngamma(v));

  if h <= hsmall then begin
    {1.0/(h*Pi) will overflow}
    Yv_series := NegInf_d;
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
  until (abs(h)<eps_d*abs(Yv));
  Yv_series := Yv;
end;


{---------------------------------------------------------------------------}
const
  BT_J = 1; BT_Y = 2;

{---------------------------------------------------------------------------}
procedure bessel_jy(v,x: double; BT: byte; var Jv,Yv: double);
  {-Return J_v(x) and/or Y_v(x) depending on BT, x>0, |v| < MaxLongint, INF if overflow}
var
  n,k: longint;
  u, Ju, Yv1, Yu, Yu1, fv, fu: double;
  w, p, q, g, curr, prev, next, init, t, x2: double;
  reflect: boolean;
  s: integer;
const
  lnepsh = -36.7368005696771013991133024373;

  {--------------------------------------------}
  function rec_overflow(a,b: double): boolean;
    {-Test if a*b overflows, if yes set Yv and Jv = PosInf_d}
  begin
    if (abs(a) > 1.0) and (abs(b) >= MaxDouble/abs(a)) then begin
      rec_overflow := true;
      Jv := PosInf_d;
      Yv := PosInf_d;
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
    if v>MAXGAMD then begin
      Yv := PosInf_d;
      if reflect then Jv := PosInf_d else Jv := 0.0;
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
          {overflows in rare cases like bessel_yv(170.9, 2).}
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
      t := maxd(3.0, v*v)*121.0;
      if BT and BT_Y <> 0 then t := maxd(t,1552);
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
      init := Sqrt_MinDbl;
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
procedure sfd_bess_jyv(v,x: double; var Jv,Yv: double);
  {-Return J_v(x) and Y_v(x), no checks, x>0, |v| < MaxLongint}
begin
  bessel_jy(v,x, BT_J + BT_Y, jv,yv);
end;


{---------------------------------------------------------------------------}
function sfd_jv(v, x: double): double;
  {-Return J_v(x), the Bessel function of the 1st kind, order v; not suitable for large v.}
var
  r,n: double;
begin

  if IsNaNorInfD(x) or IsNaNorInfD(v) then begin
    sfd_jv := NaN_d;
    exit;
  end;

  n := int(v);
  if x<0.0 then begin
    if n=v then begin
      r := sfd_jv(v, -x);
      if frac(0.5*v)<>0 then r := -r;
      sfd_jv := r;
    end
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_jv := NaN_d;
    end;
  end
  else if x=0.0 then begin
    if v=0.0 then sfd_jv := 1.0
    else if v>0.0 then sfd_jv := 0.0
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_jv := NaN_d;
    end;
  end
  else begin
    {here x > 0}
    if n=v then begin
      {integer order}
      if abs(n)<200.0 then begin
        if x > maxd(3.0, v*v)*121.0 then begin
          r := bessj_large(abs(v),x);
          if frac(0.5*v)<0.0 then r := -r;
          sfd_jv := r;
        end
        else sfd_jv := sfd_jn(round(n),x);
        exit;
      end;
    end;
    {Here v no integer or |v| > 200}
    if ((v >= 0.0) and (v<MAXGAMD-1)) and ((x < 1.0) or (v > 0.25*x*x)) then begin
      sfd_jv := IJ_series(v, x, true);
    end
    else begin
      bessel_jy(v,x,BT_J,r,n);
      sfd_jv := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_yv(v, x: double): double;
  {-Return Y_v(x), the Bessel function of the 2nd kind, order v; x > 0; not suitable for large v.}
var
  r,n: double;
begin

  if IsNaNorInfD(x) or IsNaNorInfD(v) or (x<=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_yv := NaN_d;
    exit;
  end;
  r := abs(v);
  if r=0.5 then begin
    r := sqrt(Pi_2*x);
    if v<0.0 then sfd_yv := sin(x)/r
    else sfd_yv := -cos(x)/r;
  end
  else if frac(v) = -0.5 then begin
    {Reflection would be used but cos(Pi*v)=0: Y(v,x) = sin(Pi*v)*J(-v,x)}
    n := sfd_jv(-v,x);
    if frac(0.5*(r-0.5))=0.0 then sfd_yv := n else sfd_yv := -n;
  end
  else begin
    n := int(v);
    if n=v then begin
      {integer order}
      if (x>1552.0) and (x>5.0*abs(v)) then begin
        r := bessy_large(abs(v),x);
        if frac(0.5*v)<0.0 then r := -r;
        sfd_yv := r;
      end
      else if abs(n)<2000 then begin
        sfd_yv := sfd_yn(round(n),x);
      end
      else begin
        {Call general routine but avoid Jv calculation for v<0}
        bessel_jy(abs(v),x,BT_Y,n,r);
        if frac(0.5*v)<0.0 then r := -r;
        sfd_yv := r;
      end;
    end
    else begin
      bessel_jy(v,x,BT_Y,n,r);
      sfd_yv := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure temme_k(v,x: double; var K0,K1: double);
  {-Calculate K(v, x) and K(v+1, x) by Temme's method for small |x|}
var
  k,h,p,q,f,c,d,s,s1,tol,a: double;
  g1,g2,gp,gm: double;
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
  gp := sfd_gamma1pm1(v);
  gm := sfd_gamma1pm1(-v);
  a  := ln(0.5*x);
  s  := -a*v;
  h  := exp(a*v);
  if abs(v) < eps_d then begin
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
  tol := 0.5*eps_d;
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
procedure CF2_K(v,x: double; escale: boolean; var K0,K1: double);
  {-Compute K(v,x) and K(v+1,x) via continued fraction, |v| <= 0.5, |x| > 1}
  { If escale=true the values are multiplied by exp(x)}
var
  a,a1,b,c,d,dh,ds,h,q,q1,q2,s,t: double;
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
    if abs(ds) < abs(s)*eps_d then goto done
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));

done:

  K0 := sqrt(Pi_2/x)/s;
  if not escale then K0 := K0*exp(-x);
  K1 := K0*(v+x+0.5-a1*h)/x;
end;


{---------------------------------------------------------------------------}
function bess_i_large(v,x: double; var Ivs: double): boolean;
  {-Compute I_v(x) for large x >= 100,  sqrt(x) >= 2v. Return true}
  { if 'convergence' within 50 iterations, Ivs = I_v(x)*exp(-x)*sqrt(2Pix)}
var
  s,t,u,w: double;
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
    if abs(t)<eps_d then begin
      bess_i_large := true;
      Ivs := s;
      exit;
    end;
  end;
  Ivs := s;
end;


{---------------------------------------------------------------------------}
procedure bessel_ik(v,x: double; CalcI, escale: boolean; var Iv,Kv: double);
  {-Return I_v(x) and/or K_v(x) depending on CalcI, x>0, |v| < MaxLongint}
  { If escale=true the values are exponentially scaled.}
var
  n,k: longint;
  u, Kv1, Ku, Ku1, fv: double;
  w, curr, prev, next, t, x2: double;
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
      if (t > 1.0) and (curr >= MaxDouble/t) then begin
        Kv := PosInf_d;
        if CalcI then begin
          if Reflect then Iv := PosInf_d else Iv := 0.0;
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
        if (t <= 0.0) or (x+ln(t)-0.5*ln(TwoPi*x) >= ln_MaxDbl) then begin
          Iv := PosInf_d;
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
        Iv := PosInf_d;
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
function sfd_iv(v, x: double): double;
  {-Return I_v(x), the modified Bessel function of the 1st kind, order v.}
var
  r,t: double;
const
  IvMaxX = 713.987083862762618541896578221;
begin

  if IsNaNorInfD(x) or IsNaNorInfD(v) then begin
    sfd_iv := NaN_d;
    exit;
  end;

  if x<0.0 then begin
    {if v is not an integer I(v, x) is complex}
    if frac(v)=0.0 then begin
      r := sfd_iv(v,-x);
      if frac(0.5*v)<>0 then r := -r;
      sfd_iv := r;
    end
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_iv := NaN_d;
    end;
  end
  else if abs(v)=0.5 then begin
    {NIST[30] 10.39.1: I(0.5,x)=sinh(x)/R, I(-0.5,x)=cosh(x)/R, R=sqrt(Pi*x/2)}
    if x >= ln_MaxDbl then begin
      if x >= IvMaxX then sfd_iv := PosInf_d
      else begin
        {Avoid overflow for x in range 709.78 .. 713.987}
        r := exp(0.5*x);
        sfd_iv := r*(r/sqrt(TwoPi*x));
      end;
    end
    else begin
      r := sqrt(Pi_2*x);
      if v<0.0 then sfd_iv := cosh(x)/r
      else sfd_iv := sinh(x)/r;
    end;
  end
  else if x=0.0 then begin
    if v=0.0 then sfd_iv := 1.0
    else sfd_iv := 0.0;
  end
  else begin
    {x>0}
    if v=0.0 then sfd_iv := sfd_i0(x)
    else if abs(v)=1.0 then sfd_iv := sfd_i1(x)
    else if x >= IvMaxX then sfd_iv := PosInf_d
    else begin
      bessel_ik(v,x,true,false,r,t);
      sfd_iv := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_kv(v, x: double): double;
  {-Return K_v(x), the modified Bessel function of the 2nd kind, order v, x>0}
var
  r,t: double;
begin

  if IsNaNorInfD(x) or IsNaNorInfD(v) then begin
    sfd_kv := NaN_d;
    exit;
  end;

  if (frac(v)=0.0) and (v<MaxInt) then sfd_kv := sfd_kn(round(v),x)
  else if x<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_kv := NaN_d;
  end
  else begin
    if abs(v)=0.5 then begin
      {NIST[30] 10.39.2: K(0.5,x) = K(-0.5,x) = exp(-x)*sqrt(Pi/2/x)}
      sfd_kv := exp(-x)*sqrt(Pi_2/x);
    end
    else begin
      bessel_ik(v, x, false, false, t, r);
      sfd_kv := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure sfd_bess_kv2(v,x: double; var Kv, Kv1: double);
  {-Return K(v,x) and K(v+1,x), no checks, x>0, |v| < MaxLongint}
var
  n,k: longint;
  u, Ku, Ku1: double;
  curr, prev, next, t, x2: double;
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
    if (t > 1.0) and (curr >= MaxDouble/t) then begin
      Kv  := PosInf_d;
      Kv1 := PosInf_d;
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
function sfd_ive(v, x: double): double;
  {-Return I_v(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order v.}
var
  r,t: double;
begin

  if IsNaNorInfD(x) or IsNaNorInfD(v) then begin
    sfd_ive := NaN_d;
    exit;
  end;

  if x<0.0 then begin
    {if v is not an integer I(v, x) is complex}
    if frac(v)=0.0 then begin
      r := sfd_ive(v,-x);
      if frac(0.5*v)<>0 then r := -r;
      sfd_ive := r;
    end
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_ive := NaN_d;
    end;
  end
  else if x=0.0 then begin
    if v=0.0 then sfd_ive := 1.0
    else sfd_ive := 0.0;
  end
  else begin
    {x>0}
    if v=0.0 then sfd_ive := sfd_i0e(x)
    else if abs(v)=1.0 then sfd_ive := sfd_i1e(x)
    else begin
      bessel_ik(v,x,true,true,r,t);
      sfd_ive := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_kve(v, x: double): double;
  {-Return K_v(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order v, x>0}
var
  r,t: double;
begin

  if IsNaNorInfD(x) or IsNaNorInfD(v) then begin
    sfd_kve := NaN_d;
    exit;
  end;

  if x<=0.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_kve := NaN_d;
  end
  else begin
    v := abs(v);
    if v=0.0 then sfd_kve := sfd_k0e(x)
    else if v=1.0 then sfd_kve := sfd_k1e(x)
    else if abs(v)=0.5 then begin
      {NIST[30] 10.39.1: K(0.5,x) = K(-0.5,x) = exp(-x)*sqrt(Pi/2/x)}
      sfd_kve := sqrt(Pi_2/x);
    end
    else begin
      bessel_ik(v, x, false, true, t, r);
      sfd_kve := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure sfd_bess_ikv(v,x: double; var Iv,Kv: double);
  {-Return I_v(x) and K_v(x), no checks, x>0, |v| < MaxLongint}
begin
  bessel_ik(v,x,true,false,Iv,Kv);
end;


{---------------------------------------------------------------------------}
function sfd_sph_jn(n: integer; x: double): double;
  {-Return j_n(x), the spherical Bessel function of the 1st kind, order n}
var
  r,z: double;
begin
  if IsNaNorInfD(x) then begin
    sfd_sph_jn := NaN_d;
    exit;
  end;
  if n=0 then sfd_sph_jn := sinc(x)
  else begin
    if x=0.0 then sfd_sph_jn := 0.0
    else begin
      z := abs(x);
      r := sqrt(Pi_2/z)*sfd_jv(0.5+n,z);
      if (x<0.0) and odd(n) then r := -r;    {NIST 10.47.14}
      sfd_sph_jn := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_sph_yn(n: integer; x: double): double;
  {-Return y_n(x), the spherical Bessel function of the 2nd kind, order n >=0 , x<>0}
var
  r,z: double;
begin
  if IsNaNorInfD(x) then begin
    sfd_sph_yn := NaN_d;
    exit;
  end;
  z := abs(x);
  if x=0.0 then r := NegInf_d
  else r := sqrt(Pi_2/z)*sfd_yv(0.5+n,z);
  if (x<0.0) and odd(n+1) then r := -r;      {NIST 10.47.14}
  sfd_sph_yn := r;
end;


{---------------------------------------------------------------------------}
function sfd_sph_in(n: integer; x: double): double;
  {-Return i_n(x), the modified spherical Bessel function of the 1st/2nd kind, order n}
var
  r: double;
begin
  if IsNaNorInfD(x) then begin
    sfd_sph_in := NaN_d;
    exit;
  end;
  if n=0 then sfd_sph_in := sinhc(x) {i_0 = sinh(x)/x}
  else begin
    if x=0.0 then begin
      if n>0 then sfd_sph_in := 0.0
      else if odd(n) then sfd_sph_in := PosInf_d
      else sfd_sph_in := NegInf_d
    end
    else begin
      r := abs(x);
      r := sqrt(Pi_2/r)*sfd_iv(0.5+n,r);
      if (x<0.0) and odd(n) then r := -r;  {NIST 10.47.16}
      sfd_sph_in := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_sph_ine(n: integer; x: double): double;
  {-Return i_n(x)*exp(-|x|), the exponentially scaled modified spherical Bessel function of the 1st/2nd kind, order n}
var
  r,z: double;
begin
  if IsNaNorInfD(x) then begin
    sfd_sph_ine := NaN_d;
    exit;
  end;
  z := abs(x);
  if n=0 then begin
    {i_0e = exp(-z)*sinh(z)/z}
    if z<=0.5e-8 then begin
      {exp(-z)*sinh(z)/z = 1 - z + 2/3*z^2 - 1/3*z^3 + 2/15*z^4 + O(z^5)}
      sfd_sph_ine := 1.0 - z;
    end
    else begin
      {exp(-z)*sinh(z)/z = -0.5*(exp(-2z)-1)/z}
      if z>=19 then sfd_sph_ine := 0.5/z
      else sfd_sph_ine := (-0.5)*expm1(-2.0*z)/z;
    end;
  end
  else begin
    if x=0.0 then begin
      if n>0 then sfd_sph_ine := 0.0
      else if odd(n) then sfd_sph_ine := PosInf_d
      else sfd_sph_ine := NegInf_d
    end
    else begin
      r := sqrt(Pi_2/z)*sfd_ive(0.5+n,z);
      if (x<0.0) and odd(n) then r := -r;         {NIST 10.47.16}
      sfd_sph_ine := r;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_sph_kn(n: integer; x: double): double;
  {-Return k_n(x), the modified spherical Bessel function of the 3rd kind, order n, x>0}
begin
  if IsNaNorInfD(x) then begin
    sfd_sph_kn := NaN_d;
    exit;
  end;
  if x<0.0 then begin
    sfd_sph_kn := NaN_d;
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;
  if n<0 then n := -n-1;     {NIST 10.47.9}
  if n=0 then sfd_sph_kn := Pi_2*exp(-x)/x
  else sfd_sph_kn := sqrt(Pi_2/x)*sfd_kv(0.5+n,x);
end;


{---------------------------------------------------------------------------}
function sfd_sph_kne(n: integer; x: double): double;
  {-Return k_n(x)*exp(x), the exponentially scaled modified spherical Bessel function of the 3rd kind, order n, x>0}
begin
  if IsNaNorInfD(x) then begin
    sfd_sph_kne := NaN_d;
    exit;
  end;
  if x<0.0 then begin
    sfd_sph_kne := NaN_d;
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;
  if n<0 then n := -n-1;                {NIST 10.47.9}
  if n=0 then sfd_sph_kne := Pi_2/x
  else sfd_sph_kne := sqrt(Pi_2/x)*sfd_kve(0.5+n,x);
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function Airy_small(x,f0,f1: double): double;
  {-Compute Ai(x) or Bi(x) via Maclaurin series, x small, f0=f(0), f1=f'(0)}
var
  ai,da,t1,t2,x3: double;
  k: integer;
const
  MAXIT3 = 120;
begin
  {Maclaurin series NIST[30], 9.4.1/9.4.3}
  t1 := f0;
  t2 := f1*x;
  ai := t1 + t2;
  x3 := x*sqr(x);
  k  := 3;
  while k<MAXIT3 do begin
    t1 := t1*x3/(k*(k-1));
    t2 := t2*x3/(k*(k+1));
    da := t1 + t2;
    ai := ai + da;
    if abs(da) < abs(ai)*eps_d then begin
      Airy_small := ai;
      exit;
    end;
    inc(k,3);
  end;
  {No convergence}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  Airy_small := ai;
end;


const
  {  Ai(0) =  0.35502805388781723926 =  3^(-2/3)/gamma(2/3)  NIST[30], 9.2.3}
  { Ai'(0) = -0.25881940379280679840 = -3^(-1/3)/gamma(1/3)  NIST[30], 9.2.4}
  {  Bi(0) =  0.61492662744600073516 =  3^(-1/6)/gamma(2/3)  NIST[30], 9.2.5}
  { Bi'(0) =  0.44828835735382635789 =  3^(+1/6)/gamma(1/3)  NIST[30], 9.2.6}
  rt3h: THexDblW = ($4CAA,$E858,$B67A,$3FFB);  { 1.73205080756888}
  ai0h: THexDblW = ($15B8,$9627,$B8C7,$3FD6);  {+3.5502805388781721874E-1}
  ai1h: THexDblW = ($0F8B,$42B7,$907F,$BFD0);  {-2.5881940379280682363E-1}
  bi0h: THexDblW = ($3EA9,$9B4A,$AD7A,$3FE3);  {+6.1492662744600068425E-1}
  bi1h: THexDblW = ($C8A1,$A680,$B0C1,$3FDC);  {+4.4828835735382638328E-1}
var
  rt3: double absolute rt3h; { = 3^0.5  }
  ai0: double absolute ai0h; { =  Ai(0) }
  ai1: double absolute ai1h; { = Ai'(0) }
  bi0: double absolute bi0h; { =  Bi(0) }
  bi1: double absolute bi1h; { = Bi'(0) }


{---------------------------------------------------------------------------}
function sfd_airy_ai(x: double): double;
  {-Return the Airy function Ai(x)}
var
  z,Jv,Yv: double;
begin
  if IsNaNorInfD(x) then begin
    sfd_airy_ai := NaN_d;
    exit;
  end;
  if (x <= 1.0) and (x >= -2.0) then begin
    sfd_airy_ai := Airy_small(x,ai0,ai1);
  end
  else begin
    z := abs(x);
    z := 2.0*z*sqrt(z)/THREE;
    if x>0.0 then begin
      {NR[13], 6.7.41}
      z := sfd_kv(1.0/3.0, z);
      sfd_airy_ai := sqrt(x/THREE)*z/Pi;
    end
    else begin
      {NR[13], 6.7.46}
      sfd_bess_jyv(1.0/3.0, z, Jv, Yv);
      z := Jv - Yv/rt3;
      sfd_airy_ai := 0.5*sqrt(-x)*z;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_airy_bi(x: double): double;
  {-Return the Airy function Bi(x)}
var
  z,bij,bky: double;
begin
  if IsNaNorInfD(x) then begin
    sfd_airy_bi := NaN_d;
    exit;
  end;
  z := abs(x);
  if z < 1.0 then begin
    sfd_airy_bi := Airy_small(x,bi0,bi1);
  end
  else begin
    z := 2.0*z*sqrt(z)/THREE;
    if x>0.0 then begin
      {NR[13], 6.7.44}
      sfd_bess_ikv(1.0/3.0,z,bij,bky);
      z := 2.0*bij/rt3 + bky/Pi;
      sfd_airy_bi := sqrt(x)*z;
    end
    else begin
      {NR[13], 6.7.46}
      sfd_bess_jyv(1.0/3.0, z, bij, bky);
      z := bij/rt3 + bky;
      sfd_airy_bi := (-0.5)*sqrt(-x)*z;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_airy_ais(x: double): double;
  {-Return the scaled Airy function Ai(x) if x <= 0, Ai(x)*exp(2/3*x^1.5) for x > 0}
var
  z: double;
const
  xs = 1.0;            { if x<=xs use Ai(x) else BesselK}
  xl = 17179869184.0;  {> (3d_1/eps)^(2/3)}
begin
  if IsNaNorInfD(x) then begin
    sfd_airy_ais := NaN_d;
    exit;
  end;
  if x<=0.0 then sfd_airy_ais := sfd_airy_ai(x)
  else if x < xl then begin
    z := 2.0*x*sqrt(x)/THREE;
    if x<=xs then begin
      z := exp(z);
      sfd_airy_ais := sfd_airy_ai(x)*z;
    end
    else begin
      z := sfd_kve(1.0/3.0, z);
      sfd_airy_ais := sqrt(x/THREE)*z/Pi;
    end;
  end
  else begin
    {HMF[1], 10.4.59}
    z := sqrt(x);
    sfd_airy_ais := 0.5/sqrtPi/sqrt(z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_airy_bis(x: double): double;
  {-Return the scaled Airy function Bi(x) if x <= 0, Bi(x)*exp(-2/3*x^1.5) for x > 0}
var
  z,bij,bky: double;
const
  xt = 1e-20;          {for 0<=x<xt, exp(..)=1       }
  xk = 9.5;            {K(1/3,..) negligible if x>xk }
  xl = 17179869184.0;  {> (3d_1/eps)^(2/3)}
begin
  if IsNaNorInfD(x) then begin
    sfd_airy_bis := NaN_d;
    exit;
  end;
  if x<=xt then sfd_airy_bis := sfd_airy_bi(x)
  else if x<xl then begin
    z := 2.0*x*sqrt(x)/THREE;
    bessel_ik(1.0/3.0,z,true,true,bij,bky);
    bij := 2.0*bij/rt3;
    if x<xk then
      bij := bij + bky/Pi*exp(-2*z);
    sfd_airy_bis := sqrt(x)*bij;
  end
  else begin
    z := sqrt(x);
    sfd_airy_bis := 1.0/sqrtPi/sqrt(z);
  end;
end;


{---------------------------------------------------------------------------}
function AiryP_small(x,f0,f1: double): double;
  {-Compute Ai'(x) or Bi'(x) via Maclaurin series, x small, f0=f(0), f1=f'(0)}
var
  ai,da,t1,t2,x3: double;
  k: integer;
const
  MAXIT3 = 120;
begin
  {Maclaurin series NIST[30], 9.4.2/9.4.4}
  t1 := f1;
  t2 := 0.5*f0*x*x;
  ai := t1 + t2;
  x3 := x*sqr(x);
  k  := 3;
  while k<MAXIT3 do begin
    t1 := t1*x3/(k*(k-2));
    t2 := t2*x3/(k*(k+2));
    da := t1 + t2;
    ai := ai + da;
    if abs(da) < abs(ai)*eps_d then begin
      AiryP_small := ai;
      exit;
    end;
    inc(k,3);
  end;
  AiryP_small := ai;
end;


{---------------------------------------------------------------------------}
function sfd_airy_aip(x: double): double;
  {-Return the Airy function Ai'(x)}
var
  z,Jv,Yv: double;
begin
  if IsNaNorInfD(x) then begin
    sfd_airy_aip := NaN_d;
    exit;
  end;
  z := abs(x);
  if z <= 1.0 then begin
    sfd_airy_aip := AiryP_small(x,ai0,ai1);
  end
  else begin
    z := 2.0*z*sqrt(z)/THREE;
    if x>0.0 then begin
      {NR[13], 6.7.45}
      z := sfd_kv(2.0/3.0, z);
      sfd_airy_aip := -x/rt3*z/Pi;
    end
    else begin
      {NR[13], 6.7.46}
      sfd_bess_jyv(2.0/3.0, z, Jv, Yv);
      z := Jv + Yv/rt3;
      sfd_airy_aip := (-0.5)*x*z;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_airy_bip(x: double): double;
  {-Return the Airy function Bi'(x)}
var
  z,bij,bky: double;
begin
  if IsNaNorInfD(x) then begin
    sfd_airy_bip := NaN_d;
    exit;
  end;
  z := abs(x);
  if z <= 1.0 then begin
    sfd_airy_bip := AiryP_small(x,bi0,bi1);
  end
  else begin
    z := 2.0*z*sqrt(z)/THREE;
    if x>0.0 then begin
      {NR[13], 6.7.45}
      sfd_bess_ikv(2.0/3.0,z,bij,bky);
      z := 2.0*bij/rt3 + bky/Pi;
      sfd_airy_bip := x*z;
    end
    else begin
      {NR[13], 6.7.46}
      sfd_bess_jyv(2.0/3.0, z, bij, bky);
      z := bij/rt3 - bky;
      sfd_airy_bip := (-0.5)*x*z;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function gipos7(x: double): double;
  {-Return Airy/Scorer Gi(x) for x >= 7}
const
  ntip2 = 23;
  argip2: array[0..ntip2-1] of double = (
            2.00473712275801486391e0,
            0.294184139364406724e-2,
            0.71369249006340167e-3,
            0.17526563430502267e-3,
            0.4359182094029882e-4,
            0.1092626947604307e-4,
            0.272382418399029e-5,
            0.66230900947687e-6,
            0.15425323370315e-6,
            0.3418465242306e-7,
            0.728157724894e-8,
            0.151588525452e-8,
            0.30940048039e-9,
            0.6149672614e-10,
            0.1202877045e-10,
            0.233690586e-11,
            0.43778068e-12,
            0.7996447e-13,
            0.1494075e-13,
            0.246790e-14,
            0.37672e-15,
            0.7701e-16,
            0.354e-17);
         { -0.49e-18,
            0.62e-18,
           -0.40e-18,
           -0.1e-19,
            0.2e-19,
           -0.3e-19,
            0.1e-19);}
const
  xmax = 208063.8307; {> cbrt(2/eps_d)}
var
  x3, t: double;
begin
  {Ref: MacLeod's [22], MISCFUN, function airy_gi}
  gipos7 := 0.0;
  if x<7.0 then exit
  else if x <= xmax then begin
    x3 := x*x*x;
    t := (1200.0 - x3)/(514.0 + x3);
    t := CSEvalD(t, argip2, ntip2);
    gipos7 := t/x/Pi;
  end
  else gipos7 := 1.0/Pi/x;
end;


{---------------------------------------------------------------------------}
function hineg8(x: double): double;
  {-Return Airy/Scorer Hi(x) for x <= -8}
var
  x3, t: double;
const
  xmin = -208063.831;  {< cbrt(2/eps_d)}
const
  ntin2 = 16;
  arhin2: array[0..ntin2-1] of double = (
             1.99647720399779650525e0,
             -0.187563779407173213e-2,
             -0.12186470897787339e-3,
             -0.814021609659287e-5,
             -0.55050925953537e-6,
             -0.3763008043303e-7,
             -0.258858362365e-8,
             -0.17931829265e-9,
             -0.1245916873e-10,
             -0.87171247e-12,
             -0.6084943e-13,
             -0.431178e-14,
             -0.29787e-15,
             -0.2210e-16,
             -0.136e-17,
             -0.14e-18);
begin
  {Ref: MacLeod's [22], MISCFUN, function airy_hi}
  hineg8 := 0.0;
  if x > -8.0 then exit;
  if x <= xmin then  hineg8 := (-1.0)/Pi/x
  else begin
    x3 := x*x*x;
    t  := (x3 + 1200.0)/(176.0 - x3);
    t  := CSEvalD(t, arhin2, ntin2);
    hineg8 := -t/Pi/x;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_airy_hi(x: double): double;
  {-Return the Airy/Scorer function Hi(x) = 1/Pi*integral(exp(x*t-t^3/3), t=0..INF)}
const
  hizero = 0.4099510849640004901;  {Hi(0) = 2/Gamma(2/3)/3^(7/6)}
var
  t: double;
const
  ntin4 = 20;
  arhin1: array[0..ntin4-1] of double = (
            0.31481017206423404116e0,
           -0.16414499216588964341e0,
            0.6176651597730913071e-1,
           -0.1971881185935933028e-1,
            0.536902830023331343e-2,
           -0.124977068439663038e-2,
            0.24835515596994933e-3,
           -0.4187024096746630e-4,
            0.590945437979124e-5,
           -0.68063541184345e-6,
            0.6072897629164e-7,
           -0.367130349242e-8,
            0.7078017552e-10,
            0.1187894334e-10,
           -0.120898723e-11,
            0.1189656e-13,
            0.594128e-14,
           -0.32257e-15,
           -0.2290e-16,
            0.253e-17);
          { 0.9e-19,
           -0.2e-19);}
const
  ntip1 = 30;
  arhip : array[0..ntip1-1] of double = (
            1.24013562561762831114,
            0.64856341973926535804,
            0.55236252592114903246,
            0.20975122073857566794,
            0.12025669118052373568,
            0.3768224931095393785e-1,
            0.1651088671548071651e-1,
            0.455922755211570993e-2,
            0.161828480477635013e-2,
            0.40841282508126663e-3,
            0.12196479721394051e-3,
            0.2865064098657610e-4,
            0.742221556424344e-5,
            0.163536231932831e-5,
            0.37713908188749e-6,
            0.7815800336008e-7,
            0.1638447121370e-7,
            0.319857665992e-8,
            0.61933905307e-9,
            0.11411161191e-9,
            0.2064923454e-10,
            0.360018664e-11,
            0.61401849e-12,
            0.10162125e-12,
            0.1643701e-13,
            0.259084e-14,
            0.39931e-15,
            0.6014e-16,
            0.886e-17,
            0.128e-17);
          { 0.18e-18,
            0.3e-19);}
const
  xmax: double = 104.4362038448; {overflow threshold}
begin
  {Pascal port of A.J. MacLeod's [22] MISCFUN function airy_hi. }
  {AMath uses airy_bi instead of another Chebyshev approximation}
  if x < 0.0 then begin
    if x > -eps_d then sfd_airy_hi := hizero
    else if x <= -8.0 then sfd_airy_hi := hineg8(x)
    else begin
      t := (4.0*x + 12.0)/(x - 12.0);
      sfd_airy_hi := CSEvalD(t,arhin1, ntin4);
    end
  end
  else begin
    if x < eps_d then sfd_airy_hi := hizero
    else if x > 7.0 then begin
      if x>xmax then sfd_airy_hi := PosInf_d
      else begin
        t := gipos7(x);
        sfd_airy_hi := sfd_airy_bi(x) - t;
      end;
    end
    else begin
      t := 2.0*x/7.0 - 1.0;
      t := CSEvalD(t,arhip, ntip1);
      sfd_airy_hi := t*exp(1.5*x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_airy_gi(x: double): double;
  {-Return the Airy/Scorer function Gi(x) = 1/Pi*integral(sin(x*t+t^3/3), t=0..INF)}
const
  ntin3 = 40;
  argin1: array[0..ntin3-1] of double = (
           -0.20118965056732089130e0,
           -0.7244175303324530499e-1,
            0.4505018923894780120e-1,
           -0.24221371122078791099e0,
            0.2717884964361678294e-1,
           -0.5729321004818179697e-1,
           -0.18382107860337763587e0,
            0.7751546082149475511e-1,
            0.18386564733927560387e0,
            0.2921504250185567173e-1,
           -0.6142294846788018811e-1,
           -0.2999312505794616238e-1,
            0.585937118327706636e-2,
            0.822221658497402529e-2,
            0.132579817166846893e-2,
           -0.96248310766565126e-3,
           -0.45065515998211807e-3,
            0.772423474325474e-5,
            0.5481874134758052e-4,
            0.1245898039742876e-4,
           -0.246196891092083e-5,
           -0.169154183545285e-5,
           -0.16769153169442e-6,
            0.9636509337672e-7,
            0.3253314928030e-7,
            0.5091804231e-10,
           -0.209180453553e-8,
           -0.41237387870e-9,
            0.4163338253e-10,
            0.3032532117e-10,
            0.340580529e-11,
           -0.88444592e-12,
           -0.31639612e-12,
           -0.1505076e-13,
            0.1104148e-13,
            0.246508e-14,
           -0.3107e-16,
           -0.9851e-16,
           -0.1453e-16,
            0.118e-17);
          { 0.67e-18,
            0.6e-19,
           -0.1e-19);}
const
  ntip1 = 28;
  argip1: array[0..ntip1-1] of double = (
            0.26585770795022745082e0,
           -0.10500333097501922907e0,
            0.841347475328454492e-2,
            0.2021067387813439541e-1,
           -0.1559576113863552234e-1,
            0.564342939043256481e-2,
           -0.59776844826655809e-3,
           -0.42833850264867728e-3,
            0.22605662380909027e-3,
           -0.3608332945592260e-4,
           -0.785518988788901e-5,
            0.473252480746370e-5,
           -0.59743513977694e-6,
           -0.15917609165602e-6,
            0.6336129065570e-7,
           -0.276090232648e-8,
           -0.256064154085e-8,
            0.47798676856e-9,
            0.4488131863e-10,
           -0.2346508882e-10,
            0.76839085e-12,
            0.73227985e-12,
           -0.8513687e-13,
           -0.1630201e-13,
            0.356769e-14,
            0.25001e-15,
           -0.10859e-15,
           -0.158e-17);
         {  0.275e-17,
           -0.5e-19,
           -0.6e-19);}
const
  gizero = 0.204975542482000245050307456365; {Gi(0) = 1/Gamma(2/3)/3^(7/6)}
var
  t: double;
begin
  {Pascal port of A.J. MacLeod's [22] MISCFUN function airy_gi. }
  {AMath uses airy_bi instead of another Chebyshev approximation}
  if x < -8.0 then begin
    t := hineg8(x);
    sfd_airy_gi := sfd_airy_bi(x) - t;
  end
  else if x <= -eps_d then begin
    t := -(x+4.0)/4.0;
    sfd_airy_gi := CSEvalD(t,argin1, ntin3);
  end
  else if x < eps_d then begin
    {abs(x) < eps_d}
    sfd_airy_gi := gizero;
  end
  else if x < 7.0 then begin
    t := (9.0*x - 28.0)/(x + 28.0);
    sfd_airy_gi := CSEvalD(t,argip1,ntip1);
  end
  else begin
    {x > 7}
    sfd_airy_gi := gipos7(x);
  end;
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function ber_sm(x: double): double;
  {-Return the Kelvin function ber(x), 0<=x<20}
var
  f,s,t,x4: double;
  k: integer;
begin
  {Ref: HMF [1], 9.9.10}
  k := 1;
  s := 1.0;
  t := 1.0;
  x4:= -0.0625*sqr(sqr(x));
  repeat
    f := sqr(2*k);
    f := f*sqr(2*k-1);
    t := t/f*x4;
    s := s+t;
    k := k+1;
  until abs(t) < eps_d*abs(s);
  ber_sm := s;
end;


{---------------------------------------------------------------------------}
function bei_sm(x: double): double;
  {-Return the Kelvin function bei(x), 0<=x<20}
var
  f,s,t,x4: double;
  k: integer;
begin
  {Ref: HMF [1], 9.9.10}
  k := 1;
  t := 0.25*sqr(x);
  s := t;
  x4:= -sqr(t);
  repeat
    f := sqr(2*k);
    f := f*sqr(2*k+1);
    t := t/f*x4;
    s := s+t;
    k := k+1;
  until abs(t) < eps_d*abs(s);
  bei_sm := s;
end;


{---------------------------------------------------------------------------}
function ker_sm(x: double): double;
  {-Return the Kelvin function ker(x), 0<x<3}
var
  f,h,s,t,x4: double;
  k: integer;
begin
  {Ref: HMF [1], 9.9.12}
  k := 1;
  h := 0.0;
  t := 1.0;
  x := abs(x);
  x4:= -0.0625*sqr(sqr(x));
  s := Pi_4*bei_sm(x);
  f := ln(0.5*x) + EulerGamma;
  {bei/ber are accurate, but the cancellation for s is too large for x>3}
  s := s - f*ber_sm(x);
  repeat
    f := sqr(2*k);
    f := f*sqr(2*k-1);
    t := t/f*x4;
    h := h + (one_d/(2*k) + one_d/(2*k-1));
    f := t*h;
    s := s+f;
    k := k+1;
  until abs(f) < eps_d*abs(s);
  ker_sm := s;
end;


{---------------------------------------------------------------------------}
function kei_sm(x: double): double;
  {-Return the Kelvin function kei(x), 0<x<3}
var
  f,h,s,t,x4: double;
  k: integer;
begin
  {Ref: HMF [1], 9.9.12}
  k := 1;
  h := 1.0;
  t := 0.25*sqr(x);
  x := abs(x);
  x4:= -sqr(t);
  f := ln(0.5*x) + EulerGamma;
  s := t - Pi_4*ber_sm(x);
  {bei/ber are accurate, but the cancellation for s is too large for x>3}
  s := s - f*bei_sm(x);
  repeat
    f := sqr(2*k);
    f := f*sqr(2*k+1);
    t := t/f*x4;
    h := h + (one_d/(2*k) + one_d/(2*k+1));
    f := t*h;
    s := s+f;
    k := k+1;
  until abs(f) < eps_d*abs(s);
  kei_sm := s;
end;


{---------------------------------------------------------------------------}
procedure kelvin_large(x: double; cb: boolean; var br, bi, kr, ki: double);
  {-Return all four Kelvin functions for x >= 20 using asymptotic expansions}
  { ker/kei are always calculated, ber/bei only if cb=true}
var
  f0p, {f0(+x)}
  f0n, {f0(-x)}
  g0p, {g0(+x)}
  g0n, {g0(-x)}
  tk,  {term k}
  tc,  {cos() term}
  ts,  {sin() term}
  fac,xt: double;
  k: integer;
  cbp,          {calc f0p,f0p}
  ckn: boolean; {calc f0n,g0n}

const
  kmax   = 40;           {Note: for k > 2*x the terms are increasing}
  xlarge = 1003.7843389; {>~ ln_MaxDbl*sqrt(2)}
const
  skp4hx: array[0..7] of THexDblW = (   { sin(k*Pi/4), k=0..7   }
            ($0000,$0000,$0000,$0000),  {+0.0000000000000000000 }
            ($3BCD,$667F,$A09E,$3FE6),  {+0.70710678118654757274}
            ($0000,$0000,$0000,$3FF0),  {+1.0000000000000000000 }
            ($3BCD,$667F,$A09E,$3FE6),  {+0.70710678118654757274}
            ($0000,$0000,$0000,$0000),  {+0.0000000000000000000 }
            ($3BCD,$667F,$A09E,$BFE6),  {-0.70710678118654757274}
            ($0000,$0000,$0000,$BFF0),  {-1.0000000000000000000 }
            ($3BCD,$667F,$A09E,$BFE6)); {-0.70710678118654757274}
var
  sinkp4: array[0..7] of double absolute skp4hx;
begin
  if x > xlarge then begin
    kr := 0.0;
    ki := 0.0;
    if cb then begin
      {may NaN_d?}
      br := PosInf_d;
      bi := PosInf_d;
    end;
    exit;
  end;

  {The functions are calculated using formula from HMF[1], section 9.10}
  {with nu=0 and x >= 20.}

  f0p := 1.0;    {k=0 term of sum for function f0(+x): 9.10.6}
  f0n := 1.0;    {k=0 term of sum for function f0(-x): 9.10.6}
  g0p := 0.0;    {k=0 term of sum for function g0(+x): 9.10.7}
  g0n := 0.0;    {k=0 term of sum for function g0(-x): 9.10.7}
  tk  := 1.0;    {term k of sum = prod(2j-1, j=1..k)/k!/(8x)^k}
  fac := 1.0;    {factor (+-)^k}
  cbp := cb;     {calc ber/bei}
  ckn := true;   {calc ker/kei}
  xt  := 8.0*x;

  k := 0;
  while (k<kmax) and (cbp or ckn) do begin
    fac := -fac;
    inc(k);
    tk := tk*sqr(2*k-1)/xt/k;       {term k = prod(2j-1, j=1..k)/k!/(8x)^k}
    tc := tk*sinkp4[(k+2) and 7];   {tk*cos(k*Pi/4)}
    ts := tk*sinkp4[k and 7];       {tk*sin(k*Pi/4)}
    if ckn then begin
      f0n := f0n + fac*tc;          {update and conv check f0(-x), g0(-x)}
      g0n := g0n + fac*ts;
      ckn := tk >= eps_d*mind(abs(f0n), abs(g0n));
    end;
    if cbp then begin
      f0p := f0p + tc;              {update and conv check f0(+x), g0(+x)}
      g0p := g0p + ts;
      cbp := tk >= eps_d*mind(abs(f0p),abs(g0p));
    end;
  end;

  xt := x/Sqrt2;
  if xt>=ln_MaxDbl then tk := PosInf_d else tk := exp(xt);

  {get sin(beta), cos(beta), beta = x/sqrt(2) + Pi/8, HMF 9.10.5}
  sincos(xt + 0.125*Pi, ts, tc);
  fac := sqrt(Pi_2/x)/tk;
  kr  :=  fac*(f0n*tc - g0n*ts);    {HMF 9.10.3}
  ki  := -fac*(f0n*ts + g0n*tc);    {HMF 9.10.4}

  if cb then begin
    {get sin(alpha), cos(alpha), alpha = x/sqrt(2) - Pi/8, HMF 9.10.5}
    sincos(xt - 0.125*Pi, ts, tc);
    fac := tk/sqrt(TwoPi*x);
    br  := fac*(f0p*tc + g0p*ts) - ki/Pi;   {HMF 9.10.1}
    bi  := fac*(f0p*ts - g0p*tc) + kr/Pi;   {HMF 9.10.2}
  end;
end;


{---------------------------------------------------------------------------}
procedure ker_kei_med(x: double; var kr, ki: double);
  {-Return Kelvin function kr=ker(x), ki=kei(x), 3<=x<=20}
var
  fac,ts,tc,f0,g0,z: double;
const
  csfh: array[0..22] of THexDblW = (
          ($C04A,$5098,$4C28,$BFC6),   {-1.7419914183983770261E-1}
          ($B127,$B921,$57F3,$3F5A),   {+1.6078834637723389439E-3}
          ($4377,$1064,$51B7,$3F31),   {+2.6427001320521815227E-4}
          ($698F,$2A53,$0BCA,$BF09),   {-4.7771555992520623486E-5}
          ($2A55,$CEE7,$54A8,$3EDA),   {+6.2777282736162430502E-6}
          ($2C87,$2D00,$225A,$BEA8),   {-7.1925486544564774167E-7}
          ($6E0A,$2300,$FF00,$3E71),   {+6.7040681229020080687E-8}
          ($ECEF,$CF02,$C320,$BE26),   {-2.6498710934646215641E-9}
          ($9120,$8D49,$B2D2,$BE10),   {-9.7198209650154837839E-10}
          ($7CBC,$2077,$7F40,$3DFB),   {+4.0013506437635697668E-10}
          ($8509,$BEDF,$2A95,$BDDD),   {-1.0610655385424269996E-10}
          ($9DB5,$94AC,$13D6,$3DBA),   {+2.3717341712233160606E-11}
          ($79AB,$951F,$9973,$BD94),   {-4.6837658144710679937E-12}
          ($605A,$0E80,$13C8,$3D6C),   {+7.9800404853834490867E-13}
          ($C572,$4674,$8F29,$BD3C),   {-1.0146274419705259356E-13}
          ($A5B4,$0E9E,$E14A,$3CD6),   {+1.2700971536601990071E-15}
          ($5DCD,$5F41,$F7F8,$3CF7),   {+5.3221057799364269748E-15}
          ($CA08,$E4CE,$4992,$BCE7),   {-2.5854205078178475222E-15}
          ($A386,$5612,$5440,$3CD0),   {+9.0644751109753703700E-16}
          ($88D8,$6EBA,$BE83,$BCB3),   {-2.7400572090870253519E-16}
          ($CDBA,$2896,$86FC,$3C95),   {+7.4687773792482686324E-17}
          ($95E9,$EFC6,$3F6A,$BC75),   {-1.8429464094897889406E-17}
          ($540D,$56CE,$67AE,$3C52));  {+3.9909490541622123671E-18}
       (* ($1EC1,$5D56,$E272,$BC28),   {-6.7449728433925005994E-19}
          ($EA4E,$857E,$3FC1,$3BE6),   {+3.7691351120612689369E-20}
          ($4C61,$1C8C,$3954,$3BE6));  {+3.7648818270236777832E-20} *)
const
  csgh: array[0..22] of THexDblW = (
          ($9306,$B82F,$1F46,$BFC4),   {-1.5720447534035757322E-1}
          ($3083,$F78E,$D1DB,$3F82),   {+9.1893372462197369193E-3}
          ($EBAB,$C263,$91B8,$BF42),   {-5.6668778850565061864E-4}
          ($A996,$BCE4,$DA66,$3EFF),   {+3.0377512126347748295E-5}
          ($1700,$202A,$721A,$BE8A),   {-1.9703590233373688640E-7}
          ($16B1,$6792,$EC85,$BE98),   {-3.7139520931598000832E-7}
          ($BE82,$6625,$D21D,$3E79),   {+9.6189830799939401844E-8}
          ($BF8F,$06C4,$D6C9,$BE53),   {-1.8476513139939364061E-8}
          ($668E,$CDA3,$527C,$3E2A),   {+3.0643093454233132330E-9}
          ($AA49,$DA51,$FBA3,$BDFD),   {-4.3630962238885554609E-10}
          ($1584,$33AF,$1CC4,$3DC9),   {+4.5679132751061958419E-11}
          ($3B94,$3CD8,$C9BA,$3D3A),   {+9.5170086962745579158E-14}
          ($A98E,$7079,$9E74,$BD81),   {-2.0030443265088914950E-12}
          ($9816,$AF58,$9B88,$3D6C),   {+8.1307559857898340689E-13}
          ($6D25,$1572,$2BAB,$BD51),   {-2.4400860754197408382E-13}
          ($E764,$60D2,$D560,$3D31),   {+6.3357326016347255195E-14}
          ($8A2A,$82BE,$96B9,$BD10),   {-1.4733785897064301269E-14}
          ($944F,$C8D9,$4931,$3CEB),   {+3.0293452082666044179E-15}
          ($CE9D,$5193,$4D06,$BCC2),   {-5.0795139386674074074E-16}
          ($6924,$6579,$0168,$3C8A),   {+4.5112349988246225192E-17}
          ($766C,$3E3B,$231D,$3C6F),   {+1.3503592759696792310E-17}
          ($A631,$01D4,$16BC,$BC68),   {-1.0446854432502786340E-17}
          ($51DB,$0CB6,$D987,$3C54));  {+4.5210616813281396687E-18}
       (* ($954A,$0257,$C3BE,$BC3D),   {-1.6135431781668171227E-18}
          ($0540,$41D2,$FD8F,$3C22),   {+5.1473764432777121521E-19}
          ($CE28,$3698,$2AD3,$BC06),   {-1.5021136840019478519E-19}
          ($7BDB,$FB50,$9A05,$3BE7));  {+3.9982756711602287398E-20} *)
var
  csf0: array[0..22] of double absolute csfh;
  csg0: array[0..22] of double absolute csgh;

begin
  {Medium range ker/kei: 3 <= x <= 20. The series for small x become}
  {inaccurate for x >= 3 due to the difference of the ber/bei terms.}
  z := x/Sqrt2;
  {get sin(beta), cos(beta), beta = x/sqrt(2) + Pi/8}
  sincos(z + 0.125*Pi, ts, tc);
  fac:= sqrt(Pi_2/x)*exp(-z);
  {get (asymptotic functions) f0,g0 via Chebyshev expansions}
  z  := 6.0/x - 1.0;
  f0 := 1.0 + CSEvalD(z,csf0,23)/x;
  g0 := CSEvalD(z,csg0,23)/x;
  kr :=  fac*(f0*tc - g0*ts);
  ki := -fac*(f0*ts + g0*tc);
end;


{---------------------------------------------------------------------------}
procedure sfd_ker_kei(x: double; var kr, ki: double);
  {-Return the Kelvin functions kr=ker(x), ki=kei(x), x > 0}
var
  br, bi: double;
const
  Omg = 0.422784335098467139393488;  {1-gamma}
begin
  if IsNaND(x) or (x<=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    kr := NaN_d;
    ki := NaN_d;
  end
  else if x<=1e-10 then begin
    kr := -(ln(0.5*x)+EulerGamma);
    ki := -Pi_4;
  end
  else if x<=1e-5 then begin
    {ker(x) = -gamma - ln(x/2) + 1/16*(Pi-3/8*x^2)*x^2 }
    {kei(x) = 0.25*(-Pi + (1-gamma-ln(x/2))*x^2)       }
    br := ln(0.5*x);
     x := sqr(x);
    kr := -Eulergamma - br + 0.0625*(Pi-0.375*x)*x;
    ki := 0.25*((Omg - br)*x - Pi);
  end
  else if x<3.0 then begin
    kr := ker_sm(x);
    ki := kei_sm(x);
  end
  else if x<=20.0 then ker_kei_med(x,kr,ki)
  else kelvin_large(x,false,br,bi,kr,ki);
end;


{---------------------------------------------------------------------------}
procedure sfd_ber_bei(x: double; var br, bi: double);
  {-Return the Kelvin functions br=ber(x), bi=bei(x)}
var
  kr, ki: double;
begin
  if IsNaNorInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    br := NaN_d;
    bi := NaN_d;
  end
  else begin
    x := abs(x);
    if x<=1e-5 then begin
      br := 1.0;
      bi := 0.25*sqr(x);
    end
    else if x<20.0 then begin
      br := ber_sm(x);
      bi := bei_sm(x);
    end
    else kelvin_large(x,true,br,bi,kr,ki);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ber(x: double): double;
  {-Return the Kelvin function ber(x)}
var
  br,bi,kr,ki: double;
begin
  if IsNaNorInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ber := NaN_d;
  end
  else begin
    x := abs(x);
    if x<=1e-5 then sfd_ber := 1.0
    else if x<20.0 then sfd_ber := ber_sm(x)
    else begin
      kelvin_large(x,true,br,bi,kr,ki);
      sfd_ber := br;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_bei(x: double): double;
  {-Return the Kelvin function bei(x)}
var
  br,bi,kr,ki: double;
begin
  if IsNaNorInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_bei := NaN_d;
  end
  else begin
    x := abs(x);
    if x<=1e-5 then sfd_bei := 0.25*sqr(x)
    else if x<20.0 then sfd_bei := bei_sm(x)
    else begin
      kelvin_large(x,true,br,bi,kr,ki);
      sfd_bei := bi;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ker(x: double): double;
  {-Return the Kelvin function ker(x), x > 0}
var
  kr, ki: double;
begin
  if IsNaND(x) or (x<=1e-5) or (x>=3.0) then begin
    sfd_ker_kei(x,kr,ki);
    sfd_ker := kr;
  end
  else sfd_ker := ker_sm(x);
end;


{---------------------------------------------------------------------------}
function sfd_kei(x: double): double;
  {-Return the Kelvin function kei(x), x >= 0}
var
  kr, ki: double;
begin
  if IsNaND(x) or (x<=1e-5) or (x>=3.0) then begin
    if x=0.0 then sfd_kei := -Pi_4
    else begin
      sfd_ker_kei(x,kr,ki);
      sfd_kei := ki;
    end;
  end
  else sfd_kei := kei_sm(x);
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function sfd_struve_h0(x: double): double;
  {-Return H0(x), the Struve function of order 0}
var
  t,y: double;
const
  xlow = 0.3161013638317052384022e-7;  {3*sqrt(eps_d/2)}
const
  nh0 = 19;
  arrh0:  array[0..nh0-1]  of double = (
            +0.28696487399013225740e0,
            -0.25405332681618352305e0,
            +0.20774026739323894439e0,
            -0.20364029560386585140e0,
            +0.12888469086866186016e0,
            -0.4825632815622261202e-1,
            +0.1168629347569001242e-1,
            -0.198118135642418416e-2,
            +0.24899138512421286e-3,
            -0.2418827913785950e-4,
            +0.187437547993431e-5,
            -0.11873346074362e-6,
            +0.626984943346e-8,
            -0.28045546793e-9,
            +0.1076941205e-10,
            -0.35904793e-12,
            +0.1049447e-13,
            -0.27119e-15,
            +0.624e-17);
           {-0.13e-18);}
const
  nh0a = 16;
  arrh0a: array[0..nh0a-1] of double = (
            +1.99291885751992305515e0,
            -0.384232668701456887e-2,
            -0.32871993712353050e-3,
            -0.2941181203703409e-4,
            -0.267315351987066e-5,
            -0.24681031075013e-6,
            -0.2295014861143e-7,
            -0.215682231833e-8,
            -0.20303506483e-9,
            -0.1934575509e-10,
            -0.182773144e-11,
            -0.17768424e-12,
            -0.1643296e-13,
            -0.171569e-14,
            -0.13368e-15,
            -0.2077e-16);
          { +0.2e-19,
            -0.55e-18,
            +0.10e-18,
            -0.4e-19,
            +0.1e-19);}
begin
  if IsNaND(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_struve_h0 := NaN_d;
    exit;
  end;
  {See MISCFUN [22], function STRVH0 for Chebyshev coefficients and hints. }
  {As pointed out the phase of Y0(x) is very sensitive for |x| > 1/eps_d.  }
  {Therefore in high accuracy tests x must have exactly the same value in  }
  {both calculations. My implementation uses Y0, [22] has two other series.}
  t := abs(x);
  if t <= 11.0 then begin
    if t<xlow then begin
      {H0(x) = 2/Pi*x*(1-1/9*x^2+1/225*x^4+O(x^5))}
      sfd_struve_h0 := x/Pi_2;
    end
    else begin
      t := sqr(x)/60.5 - 1.0;
      t := CSEvalD(t,arrh0,nh0);
      sfd_struve_h0 := 2.0*x*t/Pi;
    end;
  end
  else begin
    y := sfd_y0(t);
    {Correct sign if x<0 since H0 is an odd function}
    if x<0.0 then y := -y;
    if t >= 1e10 then begin
      {Avoid squaring overflow and/or unnecessary calculations}
      {H0(x) = Y0(x) + 2/Pi/x*(1 - 1/x^2 + ..}
      sfd_struve_h0 := y + 2.0/Pi/x;
    end
    else begin
      t := sqr(x);
      t := (302.5 - t) / (60.5 + t);
      t := CSEvalD(t,arrh0a,nh0a);
      t := t/Pi_2/x;
      sfd_struve_h0 := y + t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_struve_h1(x: double): double;
  {-Return H1(x), the Struve function of order 1}
var
  t,y: double;
const
  xlow = 0.9017492053581906681e-9;  {sqrt(15*eps_d/2)}
const
  nh1  =  16;
  arrh1:  array[0..nh1-1] of double = (
            +0.17319061083675439319e0,
            -0.12606917591352672005e0,
             0.7908576160495357500e-1,
            -0.3196493222321870820e-1,
             0.808040581404918834e-2,
            -0.136000820693074148e-2,
             0.16227148619889471e-3,
            -0.1442352451485929e-4,
             0.99219525734072e-6,
            -0.5441628049180e-7,
             0.243631662563e-8,
            -0.9077071338e-10,
             0.285926585e-11,
            -0.7716975e-13,
             0.180489e-14,
            -0.3694e-16);
          {  0.67e-18,
            -0.1e-19);}
const
  nh1a  = 18;
  arrh1a: array[0..nh1a-1] of double = (
            +2.01083504951473379407e0,
             0.592218610036099903e-2,
             0.55274322698414130e-3,
             0.5269873856311036e-4,
             0.506374522140969e-5,
             0.49028736420678e-6,
             0.4763540023525e-7,
             0.465258652283e-8,
             0.45465166081e-9,
             0.4472462193e-10,
             0.437308292e-11,
             0.43568368e-12,
             0.4182190e-13,
             0.441044e-14,
             0.36391e-15,
             0.5558e-16,
            -0.4e-19,
             0.163e-17);
         {  -0.34e-18,
             0.13e-18,
            -0.4e-19,
             0.1e-19);}
begin
  if IsNaND(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_struve_h1 := NaN_d;
    exit;
  end;
  {See MISCFUN [22], function STRVH1 for Chebyshev coefficients and hints.}
  {As pointed out the phase of Y1(x) is very sensitive for |x| > 1/eps_d. }
  {Therefore in high accuracy tests x must have exactly the same value in }
  {both calculations; for |x| > 1/eps_d^2 the uncertainty vanishes and the}
  {value of H1 is 2/Pi to machine precision, not 0 as in [22]. My function}
  {uses Y1, whereas MISCFUN has two other Chebyshev series.}
  t := abs(x);
  if t <= 9.0 then begin
    y := sqr(t);
    if t < xlow then begin
      {H1(x) = 2/Pi*x^2*(1/3 - 1/45*x^2 + O(x^4))}
      {Note that the factor 1/3 is missing in [22], but xlow is correct}
      sfd_struve_h1 := y/THREE/Pi_2;
    end
    else begin
      t := y/40.5 - 1.0;
      t := CSEvalD(t,arrh1,nh1);
      sfd_struve_h1 := y*t/Pi_2;
    end;
  end
  else begin
    if t >= 1e40 then begin
      {Avoid squaring overflow and/or unnecessary calculations}
      {H1(x) = Y1(x) + 2/Pi(1 - 1/x^2 + ...) }
      {H1(x) = 2/Pi - (2/Pi)^(1/2)*sin(x+Pi/4)*(1/x)^(1/2) + O((1/x)^(3/2))}
      sfd_struve_h1 := 2.0/Pi;
    end
    else begin
      y := sfd_y1(t);
      t := sqr(t);
      t := (202.5 - t) / (40.5 + t);
      t := CSEvalD(t,arrh1a,nh1a);
      sfd_struve_h1 := y + t/Pi_2;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function strv_sum(c, x: double; var err: double): double;
  {-Sum from Struve power series: sum(x^k/((1.5)_k*(c)_k)}
var
  bn,cn,max,z,sum,t: double;
  e,y,w: double;
  n: integer;
const
  NMAX = 200;
label
  done;
begin
  n := 0;
  t := 1.0;
  e := 0.0;
  bn  := 1.5;
  cn  := c;
  sum := 1.0;
  max := 0.0;

  {Kahan summation}
  while n<NMAX do begin
    z := bn*cn;
    t := t*(x/z);
    w := t - e;
    y := sum + w;
    e := (y - sum) - w;
    sum := y;
    z := abs(t);
    if z>max then max := z;
    if z <= eps_d*abs(sum) then goto done;
    bn := bn+1.0;
    cn := cn+1.0;
    inc(n);
  end;

done:
  if n>=NMAX then begin
    {No convergence}
    if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  end;
  err := max*eps_d/abs(sum);
  strv_sum := sum;
end;


{---------------------------------------------------------------------------}
function strv_ae(v,x: double): double;
  {-Asymptotic expansion for H_v and L_v, x large (about 40 for H_v, 30 for L_v)}
var
  s,b,c,max,z,sum,conv,conv1: double;
  k: integer;
const
  KMAX = 1000;
label
  done;
begin
  k := 0;
  b := 0.5;
  c := 0.5-v;
  s := 1.0;
  sum := 1.0;
  max := 0.0;
  conv  := Maxdouble;
  conv1 := conv;

  while k<KMAX do begin
    s := s*b*c*x;
    sum := sum + s;
    z := abs(s);
    if z > max then max := z
    else if (z >= conv) and (z > conv1) then goto done;
    if z <= eps_d*abs(sum) then goto done;
    conv1 := conv;
    conv  := z;
    b := b+1.0;
    c := c+1.0;
    k := k+1;
  end;

done:
  if k>=KMAX then begin
    {No convergence}
    if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  end;
  strv_ae := sum;
end;


{---------------------------------------------------------------------------}
procedure strv_int(v, x: double; var result,abserr: double; var ier: integer);
  {-Compute integral(e^(-xt)(1+t^2)^(v-0.5), t=0..INF) with custom version of }
  { AMTools.fdidei: Automatic quadrature over (a,INF) using Double Exponential.}
const
  mmax = 2*64;
  efs  = 0.1;
  hoff = 11.0;
var
  eps: double;
  neval: longint;

  function f(t: double): double;
  var
    d,z: double;
  begin
    z := -t*x;
    if z >= ln_MinDbl then begin
      inc(neval);
      d := exp(z);
      z := power(1+t*t,v-0.5);
      f := d*z;
    end
    else f := 0;
  end;

var
  m: integer;
var
  epsln,epsh,h0,ehp,ehm,epst,ir,iv,h,iback,irback,
  t,ep,em,xp,xm,fp,fm,err,errt,errh,errd: double;
begin

  result := 0.0;
  abserr := 0.0;

  neval := 0;
  ier   := 0;
  eps   := 8*eps_d;
  epsln := 1.0 - ln(efs*eps);
  epsh  := sqrt(efs*eps);
  h0    := hoff/epsln;
  ehp   := exp(h0);
  ehm   := 1.0/ehp;
  epst  := exp(-ehm*epsln);
  ir    := f(1.0);
  neval := 1;
  iv    := ir*pi_2;
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
      ep := pi_4*em;
      em := pi_4/em;
      repeat
        xp  := exp(ep-em);
        xm  := 1.0/xp;
        fp  := f(xp);
        fm  := f(xm);
        fp  := fp*xp;
        fm  := fm*xm;
        ir  := ir + (fp+fm);
        iv  := iv + (fp+fm)*(ep+em);
        errt:= (abs(fp) + abs(fm))*(ep+em);
        if m=1 then err := err + errt*epst;
        ep  := ep*ehp;
        em  := em*ehm;
      until (errt <= err) and (xm <= epsh);
      t := t + h;
    until t >= h0;
    if m=1 then begin
      errh := (err/epst)*epsh*h0;
      errd := 1.0 + 2.0*errh;
    end
    else errd := h*(abs(iv - 2.0*iback) + 4.0*abs(ir - 2.0*irback));
    h := 0.5*h;
    m := m+m;
  until (errd <= errh) or (m >= mmax);

  result := h*iv;
  if errd > errh then begin
    ier := 3;
    abserr := errd*m;
  end
  else abserr := 0.5*errh*epsh*m/efs;
end;


{---------------------------------------------------------------------------}
function sfd_struve_h(v,x: double): double;
  {-Return H_v(x), the Struve function of order v, x < 0 only if v is an integer.}
var
  g,p,f,y,z,err: double;
  ier: integer;
begin
  if IsNaNorInfD(x) or IsNaNorInfD(v) or ((x<0.0) and (frac(v)<>0.0)) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_struve_h := NaN_d;
    exit;
  end;
  if x=0.0 then begin
    if (v>-1.0) or (frac(v)=-0.5) then sfd_struve_h := 0.0
    else if v=-1.0 then sfd_struve_h := 2.0/Pi  {http://functions.wolfram.com/03.09.03.0018.01}
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_struve_h := NaN_d;
    end;
  end
  else if v=0.0 then sfd_struve_h := sfd_struve_h0(x)
  else if v=1.0 then sfd_struve_h := sfd_struve_h1(x)
  else if v=0.5 then sfd_struve_h := vers(x)/sqrt(Pi_2*x)  {HMF[1],12.1.16}
  else if frac(v)=-0.5 then sfd_struve_h := sfd_yv(v,x)
  else begin
    {First compute H_v for |x|}
    z := abs(x);
    y := abs(v);
    p := power(0.5*z, v-1.0);
    if (z>=maxd(20,y)) or (z >= maxd(4,3*y)) then begin
      g := sfd_gamma(v+0.5);
      if z < maxd(40, 1.5*y) then begin
        strv_int(v,z,f,err,ier);
        f := z*f*p/(sqrtPI*g);
      end
      else begin
        f := strv_ae(v, -sqr(2.0/z));
        f := f*p/(sqrtPI*g);
      end;
      y := sfd_yv(v,z);
      f := f + y;
    end
    else begin
      y := 0.25*z*z;
      f := strv_sum(1.5+v, -y, err);
      g := sfd_gamma(v + 1.5);
      f := y*p*f/(0.5*sqrtPI*g);
    end;
    if (x<0.0) and (f<>0.0) and (frac(0.5*v)=0.0) then f := -f;
    sfd_struve_h := f;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_struve_l(v, x: double): double;
  {-Return L_v(x), the modified Struve function of order v, x < 0 only if v is an integer.}
var
  f,g,z,p,err: double;
begin
  if IsNaNorInfD(x) or IsNaNorInfD(v) or ((x<0.0) and (frac(v)<>0.0)) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_struve_l := NaN_d;
    exit;
  end;
  if x=0.0 then begin
    if (v>-1.0) or (frac(v)=-0.5) then sfd_struve_l := 0.0
    else if v=-1.0 then sfd_struve_l := 2.0/Pi  {http://functions.wolfram.com/03.10.03.0018.01}
    else begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_struve_l := NaN_d;
    end;
  end
  else if v=0.5 then sfd_struve_l := coshm1(x)/sqrt(Pi_2*x)  {NIST[30],11.4.7}
  else if frac(v)=-0.5 then sfd_struve_l := sfd_iv(-v,x)
  else begin
    z := abs(x);
    p := power(0.5*z, v-1.0);
    if z >= maxd(30,1.5*abs(v)) then begin
      g := sfd_gamma(v+0.5);
      f := strv_ae(v,sqr(2.0/z));
      f := f*p/(sqrtPI*g);
      z := sfd_iv(-v,z);
      f := z - f;
    end
    else begin
      g := sfd_gamma(v + 1.5);
      z := sqr(0.5*z);
      f := strv_sum(1.5+v, z, err);
      f := f*p*z/(0.5*SqrtPi*g);
    end;
    if (x<0.0) and (f<>0.0) and (frac(0.5*v)=0.0) then f := -f;
    sfd_struve_l := f;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_struve_l0(x: double): double;
  {-Return L0(x), the modified Struve function of order 0}
begin
  sfd_struve_l0 := sfd_struve_l(0,x);
end;


{---------------------------------------------------------------------------}
function sfd_struve_l1(x: double): double;
  {-Return L1(x), the modified Struve function of order 1}
begin
  sfd_struve_l1 := sfd_struve_l(1,x);
end;


end.
