unit DAMCmplx;

{DAMath based complex routines}

interface

{$i STD.INC}

{$ifdef BIT16}
{$N+,F+}
{$endif}

{Robust scaled complex division should be enabled for 64-bit or SSE code!  }

{$define robust_div} {Use robust complex division algorithm by Baudin/Smith}
{$define scaled_div} {Use scaling step before internal robust division     }


(*************************************************************************

 DESCRIPTION   :  DAMath based complex routines


 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  References used in this unit, main index in damath_info.txt/references

                  [1] [HMF]: M. Abramowitz, I.A. Stegun. Handbook of Mathematical Functions. Dover, 1970
                      http://www.math.sfu.ca/~cbm/aands/
                 [13] W.H. Press et al, Numerical Recipes in C, 2nd ed., Cambridge, 1992,
                      http://www.nrbook.com/a/bookcpdf.html
                 [19] Boost C++ Libraries, Release 1.42.0, 2010, (or newer)
                      http://www.boost.org/
                 [30] [NIST]: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST Handbook
                      of Mathematical Functions, Cambridge, 2010. Online resource: NIST Digital
                      Library of Mathematical Functions, http://dlmf.nist.gov/
                 [32] D.E. Knuth: The Art of computer programming; Volume 1, Fundamental
                      Algorithms, 3rd ed., 1997; Volume 2, Seminumerical Algorithms, 3rd ed., 1998.
                 [33] http://functions.wolfram.com/: Formulas and graphics about
                      mathematical functions for the mathematical and scientific
                      community and/or http://mathworld.wolfram.com/ ("/the web's
                      most extensive mathematical resource/")
                 [35] L.C. Maximon, The dilogarithm function for complex argument, 2003,
                      Proc. R. Soc. Lond. A, 459, 2807-2819, doi: 10.1098/rspa.2003.1156;
                      http://rspa.royalsocietypublishing.org/content/459/2039/2807.full.pdf
                 [61] W. Kahan, "Branch Cuts for Complex Elementary Functions, or Much Ado
                      About Nothing's Sign Bit", in The State of Art in Numerical Analysis,
                      ed. by A. Iserles and M.J.D. Powell, 1987, pp. 165-211.
                      Available as http://people.freebsd.org/~das/kahan86branch.pdf
                 [62] PARI/GP: Open Source Number Theory-oriented Computer Algebra System,
                      available from http://pari.math.u-bordeaux.fr/
                 [63] R.M. Corless, J.H. Davenport, D.J. Jeffrey, and S.M. Watt, "According to
                      Abramowitz and Stegun" or arccoth needn't be uncouth. SIGSAM BULLETIN:
                      Communications on Computer Algebra, 34(2), June 2000.
                      https://dl.acm.org/citation.cfm?doid=362001.362023
                 [67] B.C. Carlson, Numerical computation of real or complex elliptic integrals,
                      1994, https://arxiv.org/pdf/math/9409227v1.pdf
                 [68] M. Baudin, R.L. Smith, A Robust Complex Division in Scilab,
                      2013, https://arxiv.org/pdf/1210.4539.pdf
                 [81] A. Banuelos, R.A. Depine, A program for computing the Riemann Zeta function
                      for complex argument, Computer Physics Communications 20 (1980) 441-445,
                      https://doi.org/10.1016/0010-4655(80)90021-1.
                 [82] V. Pegoraro, P. Slusallek, On the Evaluation of the Complex-Valued Exponential Integral,
                      Journal of Graphics, GPU, and Game Tools, 15(3), 183-198, 2011, available from
                      http://www.cs.utah.edu/~vpegorar/research/2011_JGT/paper.pdf


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     22.12.13  W.Ehrhardt  Initial version from CAMath
 0.11     27.12.13  we          csc,csch,sec,sech,arccsc,arccsch,arcsec,arcsech
 0.12     28.12.13  we          special case for cpow
 0.13     29.12.13  we          ccbrt as cnroot(z,3,w)
 0.14     31.12.13  we          carccot/h for very small z
 0.15     03.01.14  we          carccos, carcsin for small z
 0.16     04.01.14  we          'uses damath' moved to implementation
 0.17     07.01.14  we          Fix carctanh for very small |x| and on the cut
 0.18     28.01.14  we          DAMCmplx_Version
 0.19     28.01.14  we          branch point series
 0.20     28.01.14  we          use arccsch(z) = i*arccsc(i*z), arccothc(z) = i*arccotc(i*z)
 0.21     16.02.14  we          cagm1, cagm, csurd, ccis, cnroot1
 0.22     22.02.14  we          cset
 0.23     25.02.14  we          Improved cagm
 0.24     16.03.14  we          Improved carccsc/sec/sech without branch point series
 0.25     16.03.14  we          Improved carccoth, carccot(z) = i*carccoth(i*z)
 0.26     20.03.14  we          csgn, changed agmstep1 to allow in/out overlap
 0.27     29.03.14  we          clngamma
 0.28     30.03.14  we          cgamma
 0.29     31.03.14  we          Improved clngamma near 0,1,2
 0.30     19.04.14  we          cpowx, clog10, clogbase
 0.31     24.04.14  we          cln1p, cexpm1
 0.32     26.04.14  we          cpolar
 0.33     24.05.14  we          improved cln1p, cexpm1
 0.34     07.10.14  we          cpoly, cpolyr
 0.35     09.10.14  we          cpsi
 0.36     07.11.14  we          fix value of C_1.re
 0.37     07.11.14  we          cdilog
 0.38     10.11.14  we          csqrt1mz2
 0.39     10.11.14  we          cellk/cellck
 0.40     15.12.14  we          const z in agm1sz
 0.41     26.06.15  we          Fix: double y in csqrt1mz2
 0.42     17.07.15  we          cexp2, cexp10, cellke, celle
 0.43     17.07.15  we          improved cpolyr with faster Knuth procedure
 0.44     25.07.15  we          scaled/robust cdiv
 0.45     25.08.15  we          Use THREE & SIXX to avoid 'optimization'
 0.46     14.02.16  we          cerf/cerfc
 0.47     05.05.16  we          fix missing exit in cellke for |k|=1
 0.48     07.05.16  we          improved celle and cellke
 0.49     08.05.16  we          ek_branch: improved accuracy on the branch cut
 0.50     09.05.16  we          ek_branch: fix for k < 0
 0.51     21.05.16  we          csncndn
 0.52     30.10.17  we          cLambertWk, cLambertW
 0.53     02.12.17  we          Suppress warnings: Local variable does not seem to be initialized
 0.54     08.04.18  we          crgamma, crstheta
 0.55     19.04.18  we          Fixes for FPC311
 0.56     30.04.18  we          Fix for FPC1 illegal $warn
 0.57     04.10.18  we          cpowi
 0.58     19.10.18  we          separate function csn, ccn, cdn
 0.59     20.10.18  we          LamberWK: Changed order of arguments and initial approximations
 0.60     20.10.18  we          Fix "off by 2Pi" for some arguments in lngamma
 0.61     24.10.18  we          Fix bug in Psi(-x + 0*i)
 0.62     25.10.18  we          crstheta with Taylor series for |z| < 1/64
 0.63     26.10.18  we          improved lngamma  0 <= z.re < 2

 0.64     05.11.18  we          agm(x,-x) = 0
 0.65     05.11.18  we          czeta
 0.66     06.11.18  we          csinpi
 0.67     09.11.18  we          cei, ce1, cli
 0.68     10.11.18  we          improved expint_taylor, adjust E1 branch cut
 0.69     14.11.18  we          Fix ln1p/expm1 for small for pure imaginary z
 0.70     20.11.18  we          nroot special cases n=1,2

***************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2013-2018 Wolfgang Ehrhardt

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

const
  DAMCmplx_Version = '0.70';


type
  complex = record
              re: double; {real part     }
              im: double; {imaginary part}
            end;

const
  C_I: complex = (re: 0.0; im: 1.0);  {i = (0,1)}
  C_0: complex = (re: 0.0; im: 0.0);  {complex 0}
  C_1: complex = (re: 1.0; im: 0.0);  {complex 1}


function  cabs(const z: complex): double;          {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex absolute value |z| = sqrt(z.re^2 + z.im^2)}

procedure cadd(const x,y: complex; var z: complex);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex sum z = x + y}

procedure cagm(const x,y: complex; var w: complex);
  {-Return the 'optimal' arithmetic-geometric mean w = AGM(x,y)}

procedure cagm1(const z: complex; var w: complex);
  {-Return the 'optimal' arithmetic-geometric mean w = AGM(1,z)}

procedure carccos(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular cosine w = arccos(z)}

procedure carccosh(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic cosine w = arccosh(z)}

procedure carccot(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular cotangent w = arccot(z) = arctan(1/z)}

procedure carccotc(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular cotangent w = arccotc(z) = Pi/2 - arctan(z)}

procedure carccoth(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic cotangent w = arccoth(z) = arctanh(1/z)}

procedure carccothc(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic cotangent w = arccothc(z) = arctanh(z) + i*Pi/2}

procedure carccsc(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular cosecant w = arccsc(z) = arcsin(1/z)}

procedure carccsch(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic cosecant w = arccsch(z) = arcsinh(1/z)}

procedure carcsec(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular secant w = arcsec(z) = arccos(1/z)}

procedure carcsech(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic secant w = arcsech(z) = arccosh(1/z)}

procedure carcsin(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular sine w = arcsin(z)}

procedure carcsinh(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic sine w = arcsinh(z)}

procedure carctan(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular tangent w = arctan(z)}

procedure carctanh(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic tangent w = arctanh(z)}

function  carg(const z: complex): double;          {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the principle value of the argument or phase angle arg(z) = arctan2(z.im, z.re)}

procedure ccbrt(const z: complex; var w: complex);
  {-Return the complex principal cube root w = cbrt(z) = z^(1/3)}

procedure ccis(x: double; var z: complex);           {$ifdef HAS_INLINE} inline;{$endif}
  {-Return z = exp(i*x) = cos(x) + i*sin(x)}

procedure ccn(const z: complex; k: double; var cn: complex);
  {-Return the Jacobi elliptic function cn(z,k)}

procedure cconj(const z: complex; var w: complex);   {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex conjugate w = z.re - i*z.im}

procedure ccos(const z: complex; var w: complex);
  {-Return the complex circular cosine w = cos(z)}

procedure ccosh(const z: complex; var w: complex);
  {-Return the complex hyperbolic cosine w = cosh(z)}

procedure ccot(const z: complex; var w: complex);
  {-Return the complex circular cotangent w = cot(z)}

procedure ccoth(const z: complex; var w: complex);
  {-Return the complex hyperbolic cotangent w = coth(z)}

procedure ccsc(const z: complex; var w: complex);
  {-Return the complex circular cosecant w = csc(z) = 1/sin(z)}

procedure ccsch(const z: complex; var w: complex);
  {-Return the complex hyperbolic cosecant w = csch(z) = 1/sinh(z)}

procedure cdilog(const z: complex; var w: complex);
  {-Return the principal branch of the complex dilogarithm w = -integral(ln(1-t)/t, t=0..z)}

procedure cdiv(const x,y: complex; var z: complex);
  {-Return the quotient z = x/y}

procedure cdn(const z: complex; k: double; var dn: complex);
  {-Return the Jacobi elliptic function dn(z,k)}

procedure ce1(const z: complex; var w: complex);
  {-Return the complex exponential integral E1(z), z <> 0}

procedure cei(const z: complex; var w: complex);
  {-Return the complex exponential integral Ei(z), z <> 0}

procedure cellck(const k: complex; var w: complex);
  {-Return w = K'(k), the complementary complete elliptic integral of the first kind}

procedure celle(const k: complex; var w: complex);
  {-Return w = E(k), the complete elliptic integral of the second kind}

procedure cellk(const k: complex; var w: complex);
  {-Return w = K(k), the complete elliptic integral of the first kind}

procedure cellke(const k: complex; var kk, ek: complex);
  {-Return the complete elliptic integrals kk = K(k), ek = E(k); kk=INF if k^2=1}

procedure cerf(const z: complex; var w: complex);
  {-Return the complex error function w = erf(z) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..z)}

procedure cerfc(const z: complex; var w: complex);
  {-Return the complex complementary error function w = erfc(z) = 1-erf(z)}

procedure cexp(const z: complex; var w: complex);
  {-Return the complex exponential function w = exp(z)}

procedure cexpm1(const z: complex; var w: complex);
  {-Return w = exp(z)-1, accuracy improved for z near 0}

procedure cexp2(const z: complex; var w: complex);
  {-Return w = 2^z = exp(z*ln(2))}

procedure cexp10(const z: complex; var w: complex);
  {-Return w = 10^z = exp(z*ln(10))}

procedure cgamma(const z: complex; var w: complex);
  {-Return the complex Gamma function w = Gamma(z)}

procedure cinv(const z: complex; var w: complex);
  {-Return the complex inverse w = 1/z}

procedure cLambertW(const z: complex; var w: complex);
  {-Return the principal branch w = W(z) = W(0,z) of the Lambert W function}

procedure cLambertWk(k: integer; const z: complex; var wk: complex);
  {-Return the k'th branch wk = W(k,z) of the Lambert W function}

procedure cli(const z: complex; var w: complex);
  {-Return the complex logarithmic integral w = li(z) = Ei(ln(z)), z<>1}

procedure cln(const z: complex; var w: complex);
  {-Return the complex natural logarithm w = ln(z); principal branch ln(|z|) + i*arg(z), accurate near |z|=1}

procedure cln1p(const z: complex; var w: complex);
  {-Return the principal branch of ln(1+z), accuracy improved for z near 0}

procedure clngamma(const z: complex; var w: complex);
  {-Return w = lnGamma(z), the principal branch of the log-Gamma function}

procedure clog10(const z: complex; var w: complex);
  {-Return the principal branch of the base 10 logarithm of z, w=ln(z)/ln(10)}

procedure clogbase(const b,z: complex; var w: complex);
  {-Return the principal branch of the base b logarithm of z, w=ln(z)/ln(b)}

procedure cmul(const x,y: complex; var z: complex);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex product z = x*y}

procedure cneg(const z: complex; var w: complex);    {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the negative w = -z}

procedure cnroot(const z: complex; n: integer; var w: complex);
  {-Return the complex principal n'th root w = z^(1/n)}

procedure cnroot1(n: integer; var z: complex);
  {-Return the principal nth root of unity z = exp(2*Pi*i/n)}

procedure cpolar(const z: complex; var r,theta: double);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the polar form z = r*exp(i*theta) with r = |z|, theta = arg z}

procedure cpoly(const z: complex; const a: array of complex; n: integer; var w: complex);
  {-Evaluate polynomial; return a[0] + a[1]*z + ... + a[n-1]*z^(n-1)}

procedure cpolyr(const z: complex; const a: array of double; n: integer; var w: complex);
  {-Evaluate polynomial; return a[0] + a[1]*z + ... + a[n-1]*z^(n-1)}

procedure cpow(const z,a: complex; var w: complex);
  {-Return the principal value of the complex power w = z^a = exp(a*ln(z))}

procedure cpowi(const z: complex; n: longint; var w: complex);
  {-Return the integer power w = z^n}

procedure cpowx(const z: complex; x: double; var w: complex);
  {-Return the principal value w = z^x = |z|^x * exp(i*x*arg(z))}

procedure cpsi(const z: complex; var w: complex);
  {-Return the complex digamma function w = psi(z), z <> 0,-1,-2...}

procedure crgamma(const z: complex; var w: complex);
  {-Return the reciprocal Gamma function w = 1/Gamma(z)}

procedure crstheta(const z: complex; var w: complex);
  {-Return the Riemann-Siegel function w = theta(z)}

procedure csec(const z: complex; var w: complex);
  {-Return the complex circular secant w = sec(z) = 1/cos(z)}

procedure csech(const z: complex; var w: complex);
  {-Return the complex hyperbolic secant w = sech(z) = 1/cosh(z)}

procedure cset(var z: complex; x,y: double);   {$ifdef HAS_INLINE} inline;{$endif}
  {-Set real and imaginary part of z = x+iy}

function  csgn(const z: complex): integer;
  {-Return the sign of z. Result = isign(z.re) if z.re<>0, isign(z.im) otherwise}

procedure csin(const z: complex; var w: complex);
  {-Return the complex circular sine w = sin(z)}

procedure csinh(const z: complex; var w: complex);
  {-Return the complex hyperbolic sine w = sinh(z)}

procedure csinpi(const z: complex; var w: complex);
  {-Return the complex circular sine w = sin(Pi*z)}

procedure csn(const z: complex; k: double; var sn: complex);
  {-Return the Jacobi elliptic function sn(z,k)}

procedure csncndn(const z: complex; k: double; var sn,cn,dn: complex);
  {-Return the Jacobi elliptic functions sn(z,k), cn(z,k), dn(z,k)}
  { for complex argument z and real modulus k}

procedure csqr(const z: complex; var w: complex);    {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the square w = z^2}

procedure csqrt(const z: complex; var w: complex);
  {-Return the complex principal square root w = sqrt(z)}

procedure csqrt1mz2(const z: complex; var w: complex);
  {-Return the complex principal square root w = sqrt(1-z^2)}

procedure csub(const x,y: complex; var z: complex);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex difference z = x - y}

procedure csurd(const z: complex; n: integer; var w: complex);
  {-Return the complex n'th root w = z^(1/n) with arg(w) closest to arg(z)}

procedure ctan(const z: complex; var w: complex);
  {-Return the complex circular tangent w = tan(z)}

procedure ctanh(const z: complex; var w: complex);
  {-Return the complex hyperbolic tangent w = tanh(z)}

procedure czeta(const s: complex; var w: complex);
  {-Return the complex Riemann Zeta function, w=Zeta(s), s <> 1}

{#Z+}
{-------------------------------------------------------------------}
{---------------------- Internal functions -------------------------}
{-------------------------------------------------------------------}
{#Z-}
procedure rdivc(x: double; const y: complex; var z: complex);
  {-Return the quotient z = x/y for real x}

procedure cx_sqrt(a,b: double; var u,v: double);
  {-Return u + iv := sqrt(a + bi)}

procedure coshsinhmult(y,a,b: double; var u,v: double);
  {-Return u = a*cosh(y), v = b*sinh(y) with |a|,|b| <= 1}


implementation


uses
  DAmath;

const
  MV_4  : double = 0.44942328371e308;  {< MaxDouble/4}
  MV_125: double = 0.14381545079e309;  { ~MaxDouble/1.25}
  SLB   : double = 0.29833362925e-153; {> 2*sqrt_MinDbl}
  SUB   : double = 0.670390396497e154; {< sqrt_MaxDbl/2}


{---------------------------------------------------------------------------}
function cabs(const z:complex): double;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex absolute value |z| = sqrt(z.re^2 + z.im^2)}
begin
  cabs := hypot(z.re, z.im);
end;


{---------------------------------------------------------------------------}
function carg(const z: complex): double;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the principle value of the argument or phase angle arg(z) = arctan2(z.im, z.re)}
begin
  carg := arctan2(z.im, z.re);
end;


{---------------------------------------------------------------------------}
procedure cconj(const z: complex; var w: complex);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex conjugate w = z.re - i*z.im}
begin
  w.re :=  z.re;
  w.im := -z.im;
end;


{---------------------------------------------------------------------------}
procedure cneg(const z: complex; var w: complex);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the negative w = -z}
begin
  w.re := -z.re;
  w.im := -z.im;
end;


{---------------------------------------------------------------------------}
procedure cpolar(const z: complex; var r,theta: double);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the polar form z = r*exp(i*theta) with r = |z|, theta = arg z}
begin
  r := hypot(z.re, z.im);
  theta := arctan2(z.im, z.re);
end;


{---------------------------------------------------------------------------}
procedure csqr(const z: complex; var w: complex);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the square w = z^2}
var
  t: double;
begin
  t := (z.re - z.im)*(z.re + z.im);
  w.im := 2.0*z.re*z.im;
  w.re := t;
end;


{---------------------------------------------------------------------------}
procedure cadd(const x,y: complex; var z: complex);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex sum z = x + y}
begin
  z.re := x.re + y.re;
  z.im := x.im + y.im;
end;


{---------------------------------------------------------------------------}
procedure csub(const x,y: complex; var z: complex);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex difference z = x - y}
begin
  z.re := x.re - y.re;
  z.im := x.im - y.im;
end;


{---------------------------------------------------------------------------}
procedure cmul(const x,y: complex; var z: complex);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the complex product z = x*y}
var
  t: double;
begin
  t    := x.re*y.re - x.im*y.im;
  z.im := x.re*y.im + x.im*y.re;
  z.re := t;
end;


{---------------------------------------------------------------------------}
procedure cset(var z: complex; x,y: double);   {$ifdef HAS_INLINE} inline;{$endif}
  {-Set real and imaginary part of z = x+iy}
begin
  z.re := x;
  z.im := y;
end;


{$ifdef robust_div}
{---------------------------------------------------------------------------}
procedure cdiv(const x,y: complex; var z: complex);
  {-Return the quotient z = x/y}
var
  a,b,c,d,e,f,q,r,h: double;

{$ifdef scaled_div}
const
  OVH : THexDblW = ($0000,$0000,$0000,$7FE0);  {8.98846567431158E+0307 = MaxDouble/2 = 2^1023}
  UNH : THexDblW = ($0000,$0000,$0000,$0360);  {2.00416836000897E-0292 = MinDouble*2/eps_d = 2^(-969)}
  Beh : THexDblW = ($0000,$0000,$0000,$4680);  {4.05648192073033E+0031 = 2/eps_d^2 = 2^105}
var
  OVT: double absolute OVH;   {Overflow  threshold   }
  UNT: double absolute UNH;   {Underflow threshold   }
  Be:  double absolute Beh;   {Underflow scale factor}
var
  ab,cd,s: double;
{$endif}

begin
  a := x.re;
  b := x.im;
  c := y.re;
  d := y.im;

  e := abs(c);
  f := abs(d);

{$ifdef scaled_div}
  {Robust algorithm from Baudin/Smith [68], Fig.9. Part 1: Scaling}
  s := 1.0;
  if abs(a) > abs(b) then ab := abs(a) else ab := abs(b);
  if e>f then cd := e else cd := f;

  if ab >= OVT then begin
    {Scale down a, b}
    a := 0.5*a;
    b := 0.5*b;
    s := 2.0*s;;
  end;
  if cd >= OVT then begin
    {Scale down c, d}
    c := 0.5*c;
    d := 0.5*d;
    s := 0.5*s;
  end;
  if ab <= UNT then begin
    {Scale up a, b}
    a := a*Be;
    b := b*Be;
    s := s/Be;
  end;
  if cd <= UNT then begin
    {Scale up c, d}
    c := c*Be;
    d := d*Be;
    s := s*Be;
  end;
{$endif}

  {This is the robust internal algorithm from [68], Fig.10. DAMath implements}
  {a completely expanded code without calling other functions (except abs)!  }
  if f <= e {abs(d) <= abs(c)} then begin
    r := d/c;
    if r=0.0 then begin
      e := (a + d*(b/c))/c;
      f := (b - d*(a/c))/c;
    end
    else begin
      {/q is slightly more accurate than the original *t and also gives}
      {compatible results to the non-robust code for the 'easy' cases. }
      q := (c + d*r);
      h := b*r;  if h<>0.0 then e := (a + h)/q else e := a/q + (b/q)*r;
      h := a*r;  if h<>0.0 then f := (b - h)/q else f := b/q - (a/q)*r;
    end;
    {$ifdef scaled_div}
      z.re := s*e;
      z.im := s*f;
    {$else}
      z.re := e;
      z.im := f;
    {$endif}
  end
  else begin
    r := c/d;
    if r=0.0 then begin
      e := (b + c*(a/d))/d;
      f := (a - c*(b/d))/d;
    end
    else begin
      q := d + c*r;
      h := a*r;  if h<>0.0 then e := (b + h)/q else e := b/q + (a/q)*r;
      h := b*r;  if h<>0.0 then f := (a - h)/q else f := a/q - (b/q)*r;
    end;
    {$ifdef scaled_div}
      z.re := s*e;
      z.im := -s*f;
    {$else}
      z.re := e;
      z.im := -f;
    {$endif}
  end;
end;
{$else}
{---------------------------------------------------------------------------}
procedure cdiv(const x,y: complex; var z: complex);
  {-Return the quotient z = x/y}
var
  d,q,t: double;
begin
  {Smith's method: see Knuth[32], Exercise 4.2.1.16 and NR[13], (5.4.5)}
  if abs(y.re) >= abs(y.im) then begin
    q := y.im/y.re;
    d := y.re + q*y.im;
    t := (x.re + q*x.im)/d;
    z.im := (x.im - q*x.re)/d;
    z.re := t;
  end
  else begin
    q := y.re/y.im;
    d := y.im + q*y.re;
    t := (q*x.re + x.im)/d;
    z.im := (q*x.im - x.re)/d;
    z.re := t;
  end;
end;
{$endif}


{---------------------------------------------------------------------------}
procedure rdivc(x: double; const y: complex; var z: complex);
  {-Return the quotient z = x/y for real x}
var
  d,q: double;
begin
  {Stripped-down version of cdiv}
  if abs(y.re) >= abs(y.im) then begin
    q := y.im/y.re;
    if abs(y.re) < MV_4 then d := y.re + q*y.im
    else begin
      d := 0.5*y.re + 0.5*q*y.im;
      x := 0.5*x;
    end;
    z.im := -(q*x)/d;
    z.re := x/d;
  end
  else begin
    q := y.re/y.im;
    if abs(y.re) < MV_4 then d := y.im + q*y.re
    else begin
      d := 0.5*y.im + 0.5*q*y.re;
      x := 0.5*x;
    end;
    z.re := (q*x)/d;
    z.im := -x/d;
  end;
end;


{---------------------------------------------------------------------------}
procedure cinv(const z: complex; var w: complex);
  {-Return the complex inverse w = 1/z}
begin
  rdivc(1.0,z,w);
end;


{---------------------------------------------------------------------------}
procedure cx_sqrt(a,b: double; var u,v: double);
  {-Return u + iv := sqrt(a + bi)}
var
  x,y,r,t: double;
begin
  x := abs(a);
  y := abs(b);
  if (x=0.0) and (y=0.0) then begin
    u := 0.0;
    v := 0.0;
  end
  else begin
    {Ref: NR[13], (5.4.6/7), see also HMF[1], 3.7.27}
    if x >= y then begin
      r := y/x;
      t := x;
      r := 0.5*(1.0+sqrt(1.0+r*r));
    end
    else begin
      r := x/y;
      t := y;
      r := 0.5*(r+sqrt(1.0+r*r));
    end;
    if t<=MV_125 then t := sqrt(t*r)
    else t := sqrt(t)*sqrt(r);
    if a >= 0.0 then begin
      u := t;
      v := 0.5*(b/t);
    end
    else begin
      if b < 0.0 then t := -t;
      u := 0.5*(b/t);
      v := t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure csqrt(const z: complex; var w: complex);
  {-Return the complex principal square root w = sqrt(z)}
begin
  cx_sqrt(z.re,z.im,w.re,w.im);
end;


{---------------------------------------------------------------------------}
procedure csqrt1mz2(const z: complex; var w: complex);
  {-Return the complex principal square root w = sqrt(1-z^2)}
var
  u: complex;
  y: double;
begin
  y := z.im;
  if (abs(z.re) < SUB) and (abs(y) < SUB) then begin
    u.re := (1.0-z.re)*(1.0+z.re) + sqr(y);
    u.im := -2.0*z.re*y;
    csqrt(u,w);
  end
  else begin
    {w = +- Iz with w.re >= 0}
    if y >= 0.0 then begin
      w.im := -z.re;
      w.re :=  y;
    end
    else begin
      w.im :=  z.re;
      w.re := -y;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure cpolyr(const z: complex; const a: array of double; n: integer; var w: complex);
  {-Evaluate polynomial; return a[0] + a[1]*z + ... + a[n-1]*z^(n-1)}
var
  uj,vj,r,s,t: double;
  j: integer;
begin
  if n<=1 then begin
    w.im := 0.0;
    if n=1 then w.re := a[0]
    else w.re := 0.0;
    exit;
  end;
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('cpolyr:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
      if n>high(a)+1 then n := high(a)+1;
    end;
  {$endif}

  {Use the procedure from  Knuth [32], 4.6.4 (3). This saves about}
  {2n multiplications and n additions compared to the old version.}
  uj := a[n-1];
  vj := a[n-2];
  if n>2 then begin
    r := 2.0*z.re;
    s := sqr(z.re) + sqr(z.im);
    for j:=n-3 downto 0 do begin
      t  := vj   + r*uj;
      vj := a[j] - s*uj;
      uj := t;
    end;
  end;
  w.re := z.re*uj + vj;
  w.im := z.im*uj;
end;


{---------------------------------------------------------------------------}
procedure cpoly(const z: complex; const a: array of complex; n: integer; var w: complex);
  {-Evaluate polynomial; return a[0] + a[1]*z + ... + a[n-1]*z^(n-1)}
var
  u: complex;
  x: double;
  i: integer;
begin
  if n<=1 then begin
    if n=1 then w := a[0]
    else w := C_0;
    exit;
  end;
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('cpolyr:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
      if n>high(a)+1 then n := high(a)+1;
    end;
  {$endif}
  u := a[n-1];
  for i:=n-2 downto 0 do begin
    {u = u*z + a[i]}
    x := u.re*z.re - u.im*z.im;
    u.im := u.re*z.im + u.im*z.re + a[i].im;
    u.re := x + a[i].re;
  end;
  w.re := u.re;
  w.im := u.im;
end;


{---------------------------------------------------------------------------}
procedure cln(const z: complex; var w: complex);
  {-Return the complex natural logarithm w = ln(z); principal branch ln(|z|) + i*arg(z), accurate near |z|=1}
var
  a,x,y: double;
begin
  {Ref: HMF[1], 4.1.2/3, accuracy improved for |z| near 1}
  x := abs(z.re);
  y := abs(z.im);
  if x<y then begin
    a := x;
    x := y;
    y := a;
  end;
  if (x<0.5) or (x>1.5) or (x+y>2.0) then a := ln(hypot(x, y))
  else begin
    {avoid inaccuracies for |z|~1, eg ln(1e-20 + i) ~ 0.5e-40 + Pi/2}
    {ln(|z|) = 0.5*ln(x^2 + y^2) = 0.5*ln1p[(x^2-1) + y^2]}
    a := ln1p((x-1.0)*(x+1.0) + y*y)*0.5;
  end;
  w.im := arctan2(z.im, z.re);
  w.re := a;
end;


{---------------------------------------------------------------------------}
procedure cln1p(const z: complex; var w: complex);
  {-Return the principal branch of ln(1+z), accuracy improved for z near 0}
var
  a: double;
begin
  {w.re = ln(|1+z|), w.im = arg(1+z), z = x + i*y}
  a := maxd(abs(z.re), abs(z.im));
  if a <= 0.5*sqrt_epsh then begin
    {Use two terms ln(1+z) = z - 0.5*z^2, otherwise the result will be }
    {inaccurate for pure imaginary z, (with one term w.re would be = 0)}
    a := (z.im + z.re)*(z.im - z.re);
    w.im := z.im*(1.0 - z.re);
    w.re := z.re + 0.5*a;
    exit;
  end;
  if (a > 0.75) or (z.re < -0.5) then begin
    {This is simply the 'standard' real part; cf. cln with x<0.5 or x>1.75}
    a := ln(hypot(1.0+z.re, z.im));
  end
  else begin
    {|1+z|^2 = (1+x)^2 + y^2 = 1 + 2x + x^2 + y^2 = 1 + 2x + |z|^2}
    a := sqr(z.re) + sqr(z.im);
    {Note: This becomes inaccurate if x < 0 and 2x+a ~ 0, here the relative}
    {error of w.re can be of order eps_x/x, eg for y=1e-5, x=-0.5e-10 where}
    {it is 1.442e-9. Even when using the correctly rounded extended values }
    {the relative error of 2x+a (and of w.re) is 2.984e-10 in this case.   }
    {But also note that in many of these cases the complete relative error }
    {cabs((f-w)/f) with f=ln(1+z) remains small because |w.im| >> |w.re| !!}
    {**** If x>0 there is no cancellation and the error of w.re is small.  }
    {As for cexpm1 in the Taylor series for the case x=-y^2/2 there is no  }
    {y^2 term: i*y + i*y^3/6 + y^4/4 - i*y^5/20 - i*y^7/56 - y^8/64 ...    }
    a := 0.5*ln1p(2.0*z.re + a);
  end;
  w.im := arctan2(z.im, 1.0 + z.re);
  w.re := a;
end;


{---------------------------------------------------------------------------}
procedure clog10(const z: complex; var w: complex);
  {-Return the principal branch of the base 10 logarithm of z, w=ln(z)/ln(10)}
begin
  cln(z,w);
  w.re := w.re*log10e;
  w.im := w.im*log10e;
end;


{---------------------------------------------------------------------------}
procedure clogbase(const b,z: complex; var w: complex);
  {-Return the principal branch of the base b logarithm of z, w=ln(z)/ln(b)}
var
  u: complex;
begin
  cln(b,u);
  cln(z,w);
  cdiv(w,u,w);
end;


{---------------------------------------------------------------------------}
procedure cexp(const z: complex; var w: complex);
  {-Return the complex exponential function w = exp(z)}
var
  s,c,x: double;

  function expmul(const y: double): double;
    {-Return exp(x)*y}
  var
    t: double;
  begin
    if y=0.0 then expmul := 0.0
    else begin
      t := x + ln(abs(y));
      if t<=ln_MaxDbl then t := exp(t)
      else t := PosInf_d;
      if y<0.0 then expmul := -t
      else expmul := t;
    end;
  end;
begin
  {HMF[1], 4.3.47: exp(x + iy) = cos(y)*exp(x) + i*sin(y)*exp(x)}
  sincos(z.im,s,c);
  x := z.re;
  if x<=ln_MaxDbl then begin
    {No overflow}
    x := exp(x);
    w.re := c*x;
    w.im := s*x;
  end
  else begin
    {exp(x) will overflow, but product(s) may be finite}
    w.re := expmul(c);
    w.im := expmul(s);
  end;
end;


{---------------------------------------------------------------------------}
procedure cexpm1(const z: complex; var w: complex);
  {-Return w = exp(z)-1, accuracy improved for z near 0}
var
  a,x,y: double;
  u: complex;
  k: integer;
const
  a1 = 6.1e-5;    {~ (eps/8)^0.25}
begin
  {Note that like cln1p there are z = x+iy with relative errors for w.re of }
  {order eps_x/x if x = y^2/2. For these x, y the Taylor series starts with }
  {exp(y^2/2 + i*y)-1 = i*y + i*y^3/3 - y^4/12 + i*y^5/20 - y^6/45 + O(y^7) }
  {i.e. the (real) 2nd order term vanishes and there will be 'large' errors.}

  {BUT: Also note that in many of these cases the 'complete' relative error }
  {cabs((f-w)/f) with f = exp(z)-1 remains small because |w.im| >> |w.re| !!}

  a := maxd(abs(z.re), abs(z.im));
  if a >= a1 then begin
    if a >= 1.0 then begin
      {Use definition if |x| and |y| are large enough}
      cexp(z,w);
      w.re := w.re - 1.0;
    end
    else begin
      {With z = x + iy and a = exp(x)-1 the value for exp(z)-1 is}
      {exp(z)-1 = exp(x)*exp(iy)-1 = exp(x)*[cos(y)+i*sin(y)] - 1}
      {= exp(x)cos(y) - 1         + i*exp(x)sin(y)}
      {= exp(x)(1-vers(y)) - 1    + i*exp(x)sin(y)}
      {= exp(x)-1 - exp(x)vers(y) + i*(1+a)sin(y) }
      {= a - (1+a)vers(y)         + i*(1+a)sin(y) }
      a := expm1(z.re);
      x := 1.0 + a;
      y := vers(z.im);
      w.re := a - x*y;
      w.im := x*sin(z.im);
    end
  end
  else begin
    {Use the truncated Taylor series with backwards summation: }
    {exp(z)-1 = ((((z/k + 1)z/(k-1) + 1)z/(k-2) + 1)...)z/2 + z}
    {Number of terms from a heuristic based on Stirling formula}
    {Note that k is >= 2}
    k := 1 - floor(68/(log2(a)-3));
    u.re := 1.0;
    u.im := 0.0;
    while k>1 do begin
      x := u.re*z.re - u.im*z.im;
      y := u.re*z.im + u.im*z.re;
      u.re := 1.0 + x/k;
      u.im := y/k;
      dec(k);
    end;
    cmul(z,u,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure cexp2(const z: complex; var w: complex);
  {-Return w = 2^z = exp(z*ln(2))}
var
  r: double;
begin
  r := exp2(z.re);
  sincos(z.im*ln2, w.im, w.re);
  w.re := r*w.re;
  w.im := r*w.im;
end;


{---------------------------------------------------------------------------}
procedure cexp10(const z: complex; var w: complex);
  {-Return w = 10^z = exp(z*ln(10))}
var
  r: double;
begin
  r := exp10(z.re);
  sincos(z.im/log10e, w.im, w.re);
  w.re := r*w.re;
  w.im := r*w.im;
end;


{---------------------------------------------------------------------------}
procedure coshsinhmult(y,a,b: double; var u,v: double);
  {-Return u = a*cosh(y), v = b*sinh(y) with |a|,|b| <= 1}
var
  t: double;
begin
  {u = a*cosh(y), v = b*sinh(y)}
  if abs(y)<=ln_MaxDbl then begin
    if y=0.0 then begin
      u := a;
      v := 0.0;
    end
    else begin
      {v=sinh(y), u=cosh(y)}
      sinhcosh(y,v,u);
      u := a*u;
      v := b*v;
    end;
  end
  else begin
    {extreme case: exp(|y|) will overflow}
    t := abs(y);
    if a=0.0 then u := 0.0
    else begin
      {compute a*cosh(y)}
      u := t + ln(abs(a)) - ln2;
      if u<=ln_MaxDbl then u := exp(u)
      else u := PosInf_d;
      if a<0.0 then u := -u;
    end;
    if b=0.0 then v := 0.0
    else begin
      {compute b*sinh(y)}
      v := t + ln(abs(b)) - ln2;
      if v<=ln_MaxDbl then v := exp(v)
      else v := PosInf_d;
      if (b<0.0) <> (y<0.0) then v := -v;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure ccos(const z: complex; var w: complex);
  {-Return the complex circular cosine w = cos(z)}
var
  c,s: double;
begin
  {HMF[1], 4.3.56: cos(x + iy) = cos(x)*cosh(y) - i*sin(x)*sinh(y)}
  sincos(z.re, s,c);
  coshsinhmult(z.im, c, -s, w.re, w.im);
end;


{---------------------------------------------------------------------------}
procedure csin(const z: complex; var w: complex);
  {-Return the complex circular sine w = sin(z)}
var
  c,s: double;
begin
  {HMF[1], 4.3.55: sin(x + iy) = sin(x)*cosh(y) + i*cos(x)*sinh(y)}
  sincos(z.re, s,c);
  coshsinhmult(z.im, s, c, w.re, w.im);
end;


{---------------------------------------------------------------------------}
procedure csinpi(const z: complex; var w: complex);
  {-Return the complex circular sine w = sin(Pi*z)}
var
  c,s: double;
begin
  {HMF[1], 4.3.55: sin(x + iy) = sin(x)*cosh(y) + i*cos(x)*sinh(y)}
  sincospi(z.re, s,c);
  coshsinhmult(z.im*Pi, s, c, w.re, w.im);
end;



{---------------------------------------------------------------------------}
procedure ccosh(const z: complex; var w: complex);
  {-Return the complex hyperbolic cosine w = cosh(z)}
var
  c,s: double;
begin
  {HMF[1], 4.5.50: cosh(x + iy) = cos(y)*cosh(x) + i*sin(y)*sinh(x)}
  sincos(z.im, s,c);
  coshsinhmult(z.re, c, s, w.re, w.im);
end;


{---------------------------------------------------------------------------}
procedure csinh(const z: complex; var w: complex);
  {-Return the complex hyperbolic sine w = sinh(z)}
var
  c,s: double;
begin
  {HMF[1], 4.5.49: sinh(x + iy) = cos(y)*sinh(x) + i*sin(y)*cosh(x)}
  sincos(z.im, s,c);
  coshsinhmult(z.re, s, c, w.im, w.re);
end;


{---------------------------------------------------------------------------}
procedure ctanh(const z: complex; var w: complex);
  {-Return the complex hyperbolic tangent w = tanh(z)}
var
  x,y,c,s,h,t: double;
const
  t0 = 19.0;
begin
  {HMF[1], 4.5.51: tanh(x + iy) = (sinh(2x) + i*sin(2y))/(cosh(2x) + cos(2y))}
  {See AMath reference manual for implementation notes}
  x := z.re;
  y := z.im;
  t := abs(x);
  if y=0.0 then begin
    w.re := tanh(x);
    w.im := 0.0;
  end
  else if t=0.0 then begin
    w.re := 0.0;
    w.im := tan(y);
  end
  else if t >= t0 then begin
    w.re := isign(x);
    w.im := 2.0*sin(2.0*y)*exp(-2.0*t);
  end
  else begin
    {Note: The argument y for sincos is correct! NOT 2y,}
    {the argument doubling is implicit in the formulas. }
    sincos(y,s,c);
    {accurately compute h=exp(2x)-1 and t=4(h+1)c}
    if t<=1.0 then begin
      h := expm1(2.0*x);
      t := 4*(1.0 + h)*c;
    end
    else begin
      t := exp(2.0*x);
      h := t - 1.0;
      t := 4.0*t*c;
    end;
    x := h*h;
    {Note: theoretically y = h^2 + 4(h+1)c^2 > 0 because h+1>0}
    {but underflow may occur, therefore multiply by 1/y or Inf}
    y := x + t*c;
    if abs(y) < Mindouble then y := PosInf_d else y := 1.0/y;
    w.re := (x+2.0*h)*y;
    w.im := (t*s)*y;
  end;
end;


{---------------------------------------------------------------------------}
procedure ccoth(const z: complex; var w: complex);
  {-Return the complex hyperbolic cotangent w = coth(z)}
var
  x,y,c,s,t: double;
const
  t0 = 19.0;
begin
  {HMF[1], 4.5.52: coth(x + iy) = (sinh(2x) - i*sin(2y))/(cosh(2x) - cos(2y))}
  x := z.re;
  y := z.im;
  t := abs(x);
  if y=0.0 then begin
    w.re := coth(x);
    w.im := 0.0;
  end
  else if t=0.0 then begin
    w.re := 0.0;
    w.im := -cot(y);
  end
  else if t >= t0 then begin
    w.re := isign(x);
    w.im := -2.0*sin(2.0*y)*exp(-2.0*t);
  end
  else begin
    x := 2.0*x;
    y := 2.0*y;
    if t<0.5 then begin
      {accurately compute cosh(2x) - cos(2y) for small x}
      {cosh(2x) = coshm1(2x)+1 and -cos(2y) = vers(2y)-1}
      t := coshm1(x) + vers(y);
      {Note: theoretically  t > 0 because  cosh(2x) > 1, but}
      {underflow may occur, therefore multiply by 1/y or Inf}
      if abs(t) < Mindouble then t := PosInf_d else t := 1.0/t;
      w.re := sinh(x)*t;
      w.im := -sin(y)*t;
    end
    else begin
      sincos(y,s,c);
      sinhcosh(x,x,y);  {x=sinh(2x), y=cosh(2x)}
      {here cosh >= 1.5 and t will be > 0.5}
      t := y - c;
      w.re :=  x/t;
      w.im := -s/t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure ctan(const z: complex; var w: complex);
  {-Return the complex circular tangent w = tan(z)}
var
  u: complex;
begin
  {tan(z) = - i*tanh(iz)}
  u.re := -z.im;
  u.im :=  z.re;
  ctanh(u,u);
  w.re :=  u.im;
  w.im := -u.re;
end;


{---------------------------------------------------------------------------}
procedure ccot(const z: complex; var w: complex);
  {-Return the complex circular cotangent w = cot(z)}
var
  u: complex;
begin
  {cot(z) = i*coth(iz)}
  u.re := -z.im;
  u.im :=  z.re;
  ccoth(u,u);
  w.re := -u.im;
  w.im :=  u.re;
end;


{---------------------------------------------------------------------------}
procedure cpow(const z,a: complex; var w: complex);
  {-Return the principal value of the complex power w = z^a = exp(a*ln(z))}
var
  u: complex;
  t: double;
begin
  {Ref: NIST[30], 4.2(iv) and Kahan[61], Table 1]}
  if (a.re=0.0) and (a.im=0.0) then begin
    {z^0 = 1}
    w.re := 1.0;
    w.im := 0.0;
  end
  else if (z.re=0.0) and (z.im=0.0) and (a.re > 0.0) then begin
    {0^a = 0 if re(a) > 0}
    w.re := 0.0;
    w.im := 0.0;
  end
  else begin
    {w = exp(a*ln(z))}
    cln(z,u);
    {u := u*a}
    t    := u.re*a.re - u.im*a.im;
    u.im := u.re*a.im + u.im*a.re;
    u.re := t;
    cexp(u,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure cpowi(const z: complex; n: longint; var w: complex);
  {-Return the integer power w = z^n}
var
  x: complex;
begin
  w := c_1;
  if n=0 then exit;
  if n < 0 then begin
    n := -n;
    cinv(z,x);
  end
  else x := z;
  repeat
    if odd(n) then cmul(w,x,w);
    n := n shr 1;
    if n<>0 then csqr(x,x);
  until n=0;
end;


{$ifdef FPC}
  {$ifndef DEBUG}
    {$ifndef VER1}
      {$warn 5036 OFF}  {Local variable "xx" does not seem to be initialized}
    {$endif}
  {$endif}
{$endif}


{---------------------------------------------------------------------------}
procedure cpowx(const z: complex; x: double; var w: complex);
  {-Return the complex principal value w = z^x = |z|^x * exp(i*x*arg(z))}
var
  r,s,c,t: double;
begin
  if x=0.0 then w := C_1
  else begin
    cpolar(z,r,t);
    r := power(r,x);
    sincos(t*x,s,c);
    w.re := r*c;
    w.im := r*s;
  end;
end;


{---------------------------------------------------------------------------}
procedure ccbrt(const z: complex; var w: complex);
  {-Return the complex principal cube root w = cbrt(z) = z^(1/3)}
begin
  cnroot(z,3,w);
end;


{---------------------------------------------------------------------------}
procedure cnroot(const z: complex; n: integer; var w: complex);
  {-Return the complex principal n'th root w = z^(1/n)}
var
  r,s,c,t: double;
begin
  if n=1 then w := z
  else if n=2 then csqrt(z,w)
  else begin
    cpolar(z,r,t);
    r := nroot(r,n);
    sincos(t/n,s,c);
    w.re := r*c;
    w.im := r*s;
  end;
end;


{---------------------------------------------------------------------------}
procedure cnroot1(n: integer; var z: complex);
  {-Return the principal nth root of unity z = exp(2*Pi*i/n)}
begin
  sincos(TwoPi/n, z.im, z.re);
end;


{---------------------------------------------------------------------------}
procedure ccis(x: double; var z: complex); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return z = exp(i*x) = cos(x) + i*sin(x)}
begin
  sincos(x, z.im, z.re);
end;


{---------------------------------------------------------------------------}
procedure carctanh(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic tangent w = arctanh(z)}
var
  x,x1,y,u,v: double;
begin
  {Ref: Boost[19], function math\complex\atanh.hpp}
  x  := abs(z.re);
  if x=0.0 then begin
    {arctanh(iy)=i*arctan(y);}
    w.re := 0.0;
    w.im := arctan(z.im);
    exit;
  end;
  {here x > 0}
  x1 := 1.0-x;
  y  := abs(z.im);
  if (x>SLB) and (x<SUB) and (y>SLB) and (y<SUB) then begin
    {x and y in standard safe range}
    v := y*y;
    u := ln1p(4.0*x/(sqr(x1) + v));
    v := arctan2(2.0*y, x1*(1.0+x) - v);
    w.im := copysignd(0.5 ,z.im)*v;
    w.re := copysignd(0.25,z.re)*u;
  end
  else begin
    {special cases where standard formulas may over/underflow}
    {safe-compute (w.re) = ln1p(4*x/((1-x)^2 + y^2)}
    if x>=SUB then begin
      if (x=PosInf_d) or (y=PosInf_d) then u := 0.0
      else if y >= SUB then u := ln1p((4.0/y)/(x/y + y/x))
      else if y > 1.0 then u := ln1p(4.0/(x + y*y/x))
      else u := ln1p(4.0/x);
    end
    else if y >= SUB then begin
      if x > 1.0 then u := ln1p((4.0*x/y)/(y + x1*x1/y))
      else u := 4.0*x/y/y;
    end
    else if x <> 1.0 then begin
      u := x1*x1;
      if y > SLB then u := u + y*y;
      u := ln1p(4.0*x/u);
    end
    else u := 2.0*(ln2 - ln(y));
    {safe-compute (w.im) = arctan2(2y, (1-x^2) - y^2)}
    if (x >= SUB) or (y >= SUB) then v := Pi
    else if (x <= SLB) then begin
       if y <= SLB then v := arctan2(2.0*y, 1.0)
       else begin
         if (x=0.0) and (y=0.0) then v := 0.0
         else v := arctan2(2.0*y, 1.0 - y*y);
       end;
    end
    else begin
      {The next statement adjusts the sign on the cut: if z.im=0 then}
      {w.im=Pi/2 for z.re < -1, 0 for |z.re| < 1, -Pi/2 for z.re > 1 }
      if (y=0.0) and (z.re>1.0) then v := -Pi
      else v := arctan2(2.0*y, x1*(1.0+x));
    end;
    if (z.im < 0.0) and (v<>0.0) then v := -v;
    w.im := 0.5*v;
    w.re := copysignd(0.25,z.re)*u;
  end;
end;


{---------------------------------------------------------------------------}
procedure carctan(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular tangent w = arctan(z)}
var
  u: complex;
begin
  {Ref HMF[1], 4.4.22: arctan(z) = -i*arctanh(iz)}
  u.re := -z.im;
  u.im :=  z.re;
  carctanh(u,u);
  w.re :=  u.im;
  w.im := -u.re;
end;


{---------------------------------------------------------------------------}
procedure carccosh(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic cosine w = arccosh(z)}
var
  x,y,xm,ym,xp,yp: double;
begin
  {Ref: Kahan[61], procedure CACOSH                  }
  {u = arcsinh(re(sqrt(cconj(z)-1.0)*sqrt(z+1.0)))   }
  {v = 2.0*arctan2(im(sqrt(z-1.0), re(sqrt(z+1.0)))) }
  x := z.re;
  y := z.im;
  if (abs(x) > MV_4) or (abs(y) > MV_4) then begin
    w.re := ln(hypot(0.5*x, 0.5*y)) + 2.0*ln2;
    w.im := arctan2(y,x);
  end
  else begin
    {xp + i*yp = sqrt(z+1)}
    cx_sqrt(x+1.0, y, xp, yp);
    {xm + i*ym = sqrt(z-1)}
    cx_sqrt(x-1.0, y, xm, ym);
    {use im(sqrt(cconj(z)-1)) = -im(sqrt(z-1))}
    y := xp*xm + yp*ym;
    w.re := arcsinh(y);
    w.im := 2.0*arctan2(ym, xp);
  end;
end;


{---------------------------------------------------------------------------}
procedure carccos(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular cosine w = arccos(z)}
var
  x,y,xm,ym,xp,yp: double;
begin
  {Ref: Kahan[61], procedure CACOS               }
  {u = 2.0*arctan2(re(sqrt(1-z), re(sqrt(1+z)))) }
  {v = arcsinh(im(sqrt(1+cconj(z))*sqrt(1-z)))   }
  x  := z.re;
  xp := abs(x);
  y  := z.im;
  yp := abs(y);

  if (xp > MV_4) or (xp > MV_4) then begin
    w.re := arctan2(y,x);
    w.im := ln(hypot(0.5*x, 0.5*y)) + 2.0*ln2;
  end
  else if (xp <= sqrt_epsh) and (yp <= sqrt_epsh) then begin
    {arccos(z) = Pi/2 - z - 1/6*z^3 - 3/40*z^5 +O(z^6)}
    w.re := Pi_2-x;
    w.im := -y;
  end
  else begin
    {xp + i*yp = sqrt(1+z)}
    cx_sqrt(1.0+x, y, xp, yp);
    {xm + i*ym = sqrt(1-z)}
    cx_sqrt(1.0-x, -y, xm, ym);
    {use im(sqrt(1+cconj(z))) = -im(sqrt(1+z))}
    y := xp*ym - yp*xm;
    w.re:= 2.0*arctan2(xm, xp);
    w.im := arcsinh(y);
  end;
end;


{---------------------------------------------------------------------------}
procedure carcsin(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular sine w = arcsin(z)}
var
  x,y,xm,ym,xp,yp: double;
begin
  {Ref: Kahan[61], procedure CASIN             }
  {u = arctan2(re(z), re(sqrt(1-z)*sqrt(1+z))) }
  {v = arcsinh(im(sqrt(1-cconj(z))*sqrt(1+z))) }
  x  := z.re;
  xp := abs(x);
  y  := z.im;
  yp := abs(y);
  if (abs(x) > MV_4) or (abs(y) > MV_4) then begin
    if y >= 0.0 then w.re := arctan2(x,y)
    else w.re := arctan2(-x,-y);
    w.im := ln(hypot(0.5*x, 0.5*y)) + 2.0*ln2;
    w.re := copysignd(w.re,x);
    w.im := copysignd(w.im,y);
  end
  else if (xp <= sqrt_epsh) and (yp <= sqrt_epsh) then begin
    {arcsin(z) = z + 1/6*z^3 + 3/40*z^5 +O(z^6) }
    w.re := x;
    w.im := y;
  end
  else begin
    {xp + i*yp = sqrt(1+z)}
    cx_sqrt(1.0+x, y, xp, yp);
    {xm + i*ym = sqrt(1-z)}
    cx_sqrt(1.0-x, -y, xm, ym);
    {use im(sqrt(1-cconj(z))) = -im(sqrt(1-z)) = -ym}
    y := yp*xm - xp*ym;
    w.re := arctan2(x, xp*xm - yp*ym);
    w.im := arcsinh(y);
  end;
end;


{---------------------------------------------------------------------------}
procedure carcsinh(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic sine w = arcsinh(z)}
var
  u: complex;
begin
  {HMF[1], 4.6.14: arcsinh(z) = - i*arcsin(iz)}
  u.re := -z.im;
  u.im :=  z.re;
  carcsin(u,u);
  w.re :=  u.im;
  w.im := -u.re;
end;


{---------------------------------------------------------------------------}
procedure carccot(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular cotangent w = arccot(z) = arctan(1/z)}
var
  u: complex;
begin
  {arccot(z) = i*arccoth(i*z))}
  u.re := -z.im;
  u.im :=  z.re;
  carccoth(u,u);
  w.re := -u.im;
  w.im :=  u.re;
end;


{---------------------------------------------------------------------------}
procedure carccotc(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular cotangent w = arccotc(z) = Pi/2 - arctan(z)}
begin
  if z.re > 0.0 then carccot(z, w)
  else if (abs(z.re) < eps_d) and (abs(z.im) < eps_d) then begin
    {arccotc(z) = Pi_2 - z + z^3/3 + O(z^5)}
    w.re := Pi_2;
    w.im := -z.im;
  end
  else begin
    carctan(z,w);
    w.re := Pi_2 - w.re;
    w.im := -w.im;
  end;
end;


{---------------------------------------------------------------------------}
procedure carccoth(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic cotangent w = arccoth(z) = arctanh(1/z)}
var
  u: complex;
  x: double;
begin
  x := abs(z.re);
  if z.im=0.0 then begin
    {z is real}
    if x >= 1.0 then begin
      w.re := arctanh(1.0/z.re);
      w.im := 0.0;
    end
    else begin
      w.re := arctanh(z.re);
      if z.re>0.0 then w.im := -Pi_2 else w.im := Pi_2;
    end;
  end
  else if (abs(1.0-x) < 0.5) and (abs(z.im) < 0.5) then begin
    {See [63] 5.4: Definition of arccoth, and appendix Lemma 3}
    {arccoth(z) = 0.5*[ln(-1-z) - ln(1-z)]}
    u.re := -1.0 - z.re;
    u.im := -z.im;
    w.re := 1.0 - z.re;
    w.im := -z.im;
    cln(u,u);
    cln(w,w);
    w.re := 0.5*(u.re - w.re);
    w.im := 0.5*(u.im - w.im);
  end
  else begin
    {arccoth(z) = arctanh(1/z)}
    rdivc(1.0, z, w);
    carctanh(w,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure carccothc(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic cotangent w = arccothc(z) = arctanh(z) + i*Pi/2}
begin
  if z.im < 0.0 then carccoth(z,w)
  else begin
    carctanh(z,w);
    w.re := w.re;
    w.im := w.im + Pi_2;
  end;
end;


{---------------------------------------------------------------------------}
procedure csec(const z: complex; var w: complex);
  {-Return the complex circular secant w = sec(z) = 1/cos(z)}
begin
  ccos(z,w);
  rdivc(1.0, w, w);
end;


{---------------------------------------------------------------------------}
procedure ccsc(const z: complex; var w: complex);
  {-Return the complex circular cosecant w = csc(z) = 1/sin(z)}
begin
  csin(z,w);
  rdivc(1.0, w, w);
end;


{---------------------------------------------------------------------------}
procedure csech(const z: complex; var w: complex);
  {-Return the complex hyperbolic secant w = sech(z) = 1/cosh(z)}
begin
  ccosh(z,w);
  rdivc(1.0, w, w);
end;


{---------------------------------------------------------------------------}
procedure ccsch(const z: complex; var w: complex);
  {-Return the complex hyperbolic cosecant w = csch(z) = 1/sinh(z)}
begin
  csinh(z,w);
  rdivc(1.0, w, w);
end;


{---------------------------------------------------------------------------}
procedure carccsc(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular cosecant w = arccsc(z) = arcsin(1/z)}
var
  u: complex;
  xm,ym,xp,yp: double;
begin
  xp := abs(z.re);
  if (abs(xp-1.0) < 0.5) and (abs(z.im) < 0.5) then begin
    {This is the code for arcsin(1/z) with 1 +- 1/z = (z +- 1)/z}
    {xp + i*yp = sqrt(1+1/z) = sqrt((z+1)/z)}
    u.re := z.re + 1.0;
    u.im := z.im;
    cdiv(u,z,u);
    cx_sqrt(u.re, u.im, xp, yp);
    {xm + i*ym = sqrt(1-1/z) = sqrt((z-1)/z)}
    u.re := z.re - 1.0;
    u.im := z.im;
    cdiv(u,z,u);
    cx_sqrt(u.re, u.im, xm, ym);
    rdivc(1.0,z,u);
    w.re := arctan2(u.re, xp*xm - yp*ym);
    w.im := arcsinh(yp*xm - xp*ym);
  end
  else begin
    rdivc(1.0, z, w);
    carcsin(w,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure carcsec(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse circular secant w = arcsec(z) = arccos(1/z)}
var
  u: complex;
  xm,ym,xp,yp: double;
begin
  xp := abs(z.re);
  if (abs(xp-1.0) < 0.5) and (abs(z.im) < 0.5) then begin
    {This is the code for arccos(1/z) with 1 +- 1/z = (z +- 1)/z}
    {xp + i*yp = sqrt(1+1/z) = sqrt((1+z)/z}
    u.re := z.re + 1.0;
    u.im := z.im;
    cdiv(u,z,u);
    cx_sqrt(u.re, u.im, xp, yp);
    {xm + i*ym = sqrt(1-1/z)=sqrt(z-1)/z}
    u.re := z.re - 1.0;
    u.im := z.im;
    cdiv(u,z,u);
    cx_sqrt(u.re, u.im, xm, ym);
    w.re := 2.0*arctan2(xm, xp);
    w.im := arcsinh(xp*ym - yp*xm);
  end
  else begin
    rdivc(1.0, z, w);
    carccos(w,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure carcsech(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic secant w = arcsech(z) = arccosh(1/z)}
var
  u: complex;
  xm,ym,xp,yp: double;
begin
  xp := abs(z.re);
  if (abs(xp-1.0) < 0.5) and (abs(z.im) < 0.5) then begin
    {This is the code for arccosh(1/z) with 1/z +- 1 = (1 +- z)/z}
    {xp + i*yp = sqrt(1/z+1) = sqrt((1+z)/z}
    u.re := z.re + 1.0;
    u.im := z.im;
    cdiv(u,z,u);
    cx_sqrt(u.re, u.im, xp, yp);
    {xm + i*ym = sqrt(1/z-1) = sqrt((1-z)/z)}
    u.re := 1.0 - z.re;
    u.im := -z.im;
    cdiv(u,z,u);
    cx_sqrt(u.re, u.im, xm, ym);
    w.re := arcsinh(xp*xm + yp*ym);
    w.im := 2.0*arctan2(ym, xp);
  end
  else begin
    rdivc(1.0, z, w);
    carccosh(w,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure carccsch(const z: complex; var w: complex);
  {-Return the principal value of the complex inverse hyperbolic cosecant w = arccsch(z) = arcsinh(1/z)}
var
  u: complex;
begin
  {Use arccsch(z) = i*arccsc(i*z)}
  u.re := -z.im;
  u.im :=  z.re;
  carccsc(u,u);
  w.re := -u.im;
  w.im :=  u.re;
end;

(*
Wolfram Alpha:
AGM(1,-2)
-0.42296620840880168736459740606  +0.66126618346180476446723986556*I
AGM(1,-1/2)
 0.21148310420440084368229870303  +0.33063309173090238223361993278*I

Maple
AGM(1,-2)
-0.422966208408801687364597406061 -0.661266183461804764467239865563*I
AGM(1,-1/2)
 0.211483104204400843682298703030 +0.330633091730902382233619932782*I

Pari 2.6.2
agm(1,-2)
-0.4229662084088016873645974061   -0.6612661834618047644672398656*I
agm(1,-0.5)
 0.2114831042044008436822987030   -0.3306330917309023822336199328*I
*)

{---------------------------------------------------------------------------}
procedure agm1sz(const z: complex; var w: complex);
  {-Compute 'optimal' w = AGM(1,z), 'small' z:  |z| <= 2*Sqrt_MaxExt}
var
  a,b: complex;
  r: integer;
  done: boolean;
begin

  if (z.im=0.0) and ((z.re=0.0) or (z.re+1.0=0.0)) then begin
    w.re := 0.0;
    w.im := 0.0;
    exit;
  end;

  {Compute suitable starting values a, b for optimal AGM, see Pari/GP }
  {source code and the discussion in the pari-dev thread 'Complex AGM'}
  {http://pari.math.u-bordeaux.fr/archives/pari-dev-1202/msg00045.html}
  {or http://comments.gmane.org/gmane.comp.mathematics.pari.devel/3543}
  if z.re < 0.0 then begin
    {If the second condition is included, then the values on the cut are}
    {conjugated if re.z<1 like Maple, without it we have Wolfram values.}
    if (z.im < 0.0) {or ((z.im=0.0) and (z.re<-1.0))} then begin
      {a := +I*a}
      r := -1;
      a.re := -0.5*z.im;
      a.im := +0.5*(z.re + 1.0);
    end
    else begin
      {a := -I*a}
      r := 1;
      a.re := +0.5*z.im;
      a.im := -0.5*(z.re + 1.0);
    end;
    cx_sqrt(-z.re, -z.im, b.re, b.im);
  end
  else begin
    {no rotation}
    r := 0;
    a.re := 0.5*(z.re + 1.0);
    a.im := 0.5*z.im;
    cx_sqrt(z.re, z.im, b.re, b.im);
  end;
  {Here a.re >= 0 and b.re >= 0}

  {standard AGM iteration}
  repeat
    {w = a*b}
    w.re := a.re*b.re - a.im*b.im;
    w.im := a.re*b.im + a.im*b.re;
    {a = (a+b)/2}
    a.re := 0.5*(a.re+b.re);
    a.im := 0.5*(a.im+b.im);
    {b = sqrt(a*b)}
    cx_sqrt(w.re, w.im, b.re, b.im);
    {loop until a-b is small of order sqrt(eps)}
    done := (abs(a.re-b.re)<=sqrt_epsh*abs(a.re)) and
            (abs(a.im-b.im)<=sqrt_epsh*abs(a.im))
  until done;

  {final iteration makes |a-b| ~ eps*|a|}
  a.re := 0.5*(a.re+b.re);
  a.im := 0.5*(a.im+b.im);

  {undo rotation}
  case r of
      1: begin
           {w=a*I}
           w.re := -a.im;
           w.im :=  a.re;
         end;
     -1: begin
           {w=-a*I}
           w.re :=  a.im;
           w.im := -a.re;
         end;
    else begin
           w.re := a.re;
           w.im := a.im;
         end;
  end; {case}
end;


{---------------------------------------------------------------------------}
procedure cagm1(const z: complex; var w: complex);
  {-Return the 'optimal' arithmetic-geometric mean w = AGM(1,z)}
var
  u: complex;
begin

  {'Optimal' means: |w| is maximal over all possible AGM sequences}
  {which can be obtained by choosing one of the two sqrt branches.}

  if (abs(z.re) > Sqrt_MaxDbl) or (abs(z.im) > Sqrt_MaxDbl) then begin
    {Avoid possible overflow in AGM iteration loop and}
    {use  AGM(1,z) = AGM(1/z,1)*z = AGM(1/z,1) / (1/z)}
    cinv(z,u);
    agm1sz(u,w);
    cdiv(w,u,w);
  end
  else agm1sz(z,w);
end;


{---------------------------------------------------------------------------}
procedure agmstep1(const x,y: complex; var a,b: complex);
  {-compute one step of optimal AGM}
var
  u,v: complex;
begin
  cmul(x,y,u);
  a.re := 0.5*(x.re + y.re);
  a.im := 0.5*(x.im + y.im);
  csqrt(u,b);
  csub(a,b,u);
  cadd(a,b,v);
  if cabs(u)>cabs(v) then begin
    {if |a-b| > |a+b| then select other branch of sqrt: b -> -b}
    b.re := -b.re;
    b.im := -b.im;
  end;
end;


{---------------------------------------------------------------------------}
procedure cagm(const x,y: complex; var w: complex);
  {-Return the 'optimal' arithmetic-geometric mean w = AGM(x,y)}
var
  u,v: complex;
  ax,ay: double;
  zerostep: boolean;
begin
  ax := cabs(x);
  ay := cabs(y);
  if (ax=0.0) or (ay=0.0) or ((x.re = -y.re) and (x.im = -y.im)) then w := C_0
  else begin
    zerostep := false;
    if ax>ay then begin
      if ay/ax>succd0 then begin
        cdiv(y,x,v);
        u := x;
      end
      else zerostep := true;
    end
    else begin
      if ax/ay>succd0 then begin
        u := y;
        cdiv(x,y,v);
      end
      else zerostep := true;
    end;
    if zerostep then begin
      {quotient might underflow, take one manual step and a recursive call}
      agmstep1(x,y,u,v);
      cagm(u,v,w);
    end
    else begin
      cagm1(v,w);
      cmul(w,u,w);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure cellk(const k: complex; var w: complex);
  {-Return w = K(k), the complete elliptic integral of the first kind}
begin
  {K(k) = Pi/2 / AGM(1,k'), see NIST[30], 19.8.5 or}
  {or http://functions.wolfram.com/08.02.26.0133.01}
  csqrt1mz2(k,w);
  cagm1(w,w);
  rdivc(Pi_2,w,w);
end;


{---------------------------------------------------------------------------}
procedure cellck(const k: complex; var w: complex);
  {-Return w = K'(k), the complementary complete elliptic integral of the first kind}
begin
  {K'(k) = K(k') = Pi/2 / AGM(1,k)}
  cagm1(k,w);
  rdivc(Pi_2,w,w);
end;


{---------------------------------------------------------------------------}
procedure ek_large(const k: complex; var ek: complex);
  {-Return E(k) for |k| >= 100 using asymptotic series}
var
  z, r, lr, s1, s2: complex;
  a: double;
  n: integer;
const
  p1: array[0..5] of double = (0.943147180559945309417232121458,
                               0.397683975699931636771540151822e-1,
                               0.90537740887474363789327556933e-2,
                               0.33930789524726231140274769236e-2,
                               0.162431188378010665734182961567e-2,
                               0.90014434723516574370573499061e-3);
  p2: array[0..5] of double = (1/4, 1/32, 3/256, 25/4096, 245/65536, 1323/524288);
const
  a1 = 1e19; {~1.0/eps_x}
  a2 = 1e5;
  a3 = 1000.0;
  a4 = 300;
  a5 = 80;
begin
  {Use series http://functions.wolfram.com/08.01.06.0010.02}
  {Because the different definition, Wolfram's z is k^2.   }
  {Therefore sqrt(-z) is computed inline using three cases }
  {and ln(-z) is based on ln(k) to avoid overflows.        }

  {lr = ln(-k^2)}
  if k.im=0.0 then begin
    lr.re := 2.0*ln(abs(k.re));
    lr.im := Pi;
  end
  else begin
    cln(k,lr);
    lr.re := 2.0*lr.re;
    lr.im := 2.0*lr.im;
    if lr.im > 0.0 then lr.im := lr.im - Pi
    else lr.im := lr.im + Pi;
  end;

  {r = sqrt(-k^2)}
  if k.im = 0.0 then begin
    r.re := 0.0;
    r.im := abs(k.re);
  end
  else if k.im > 0.0 then begin
    r.re :=  k.im;
    r.im := -k.re;
  end
  else begin
    r.re := -k.im;
    r.im :=  k.re;
  end;

  {get number of terms in series}
  a := cabs(k);
  if a >= a1 then n := 1
  else if a >= a2 then n := 2
  else if a >= a3 then n := 3
  else if a >= a4 then n := 4
  else if a >= a5 then n := 5
  else n := 6;

  if n=1 then begin
    {only constant terms, don't compute 1/k^2 and avoid overflow}
    s1.re := p1[0];
    s1.im := 0.0;
    s2.re := p2[0];
    s2.im := 0.0;
  end
  else begin
    csqr(k,z);
    cinv(z,z);
    cpolyr(z, p1, n, s1);
    cpolyr(z, p2, n, s2);
  end;

  cdiv(s1,r,s1);
  cmul(s2,lr,s2);
  cdiv(s2,r,s2);
  ek.re := r.re + s1.re + s2.re;
  ek.im := r.im + s1.im + s2.im;
end;


{---------------------------------------------------------------------------}
function cel2_10(kc: double): double;
  {-Return cel2(kc,1,0), kc <>0}
var
  a,b,c,m,h: double;
  done: boolean;
begin
  kc := abs(kc);
  {Stripped down version of cel2, see sdellint}
  m := 1.0;
  c := 1.0;
  a := 1.0;
  b := 0.0;
  done := false;
  repeat
    b := 2.0*(c*kc + b);
    c := a;
    h := m;
    m := kc+m;
    a := b/m + a;
    if abs(h-kc) <= sqrt_epsh*h then done := true
    else begin
      kc:= 2.0*sqrt(h*kc);
    end;
  until done;
  cel2_10 := Pi_4*a/m;
end;


{---------------------------------------------------------------------------}
procedure ek_branch(k: double; var ek: complex);
  {-Return E(k), k on the branch cut, i.e. k real and |k|>1}
var
  z,t: double;
begin
  {Use Legendre relation NIST[30], 19.7.3 for E(1/k). Note that the}
  {direct evaluation of E(k)-kc^2*K(k) and E(kc)-k^2*K(kc) suffers }
  {from cancellation and is replaced by single calls to cel2(,1,0),}
  {analogous to the case |k|>1 for the real function sfd_ellint_2. }
  k := abs(k);
  z := 1.0/k;
  t := sqrt((1.0-z)*(1.0+z));
  ek.re := cel2_10(t)/k;
  ek.im := cel2_10(z)*((k*k-1.0)/k);
end;


{---------------------------------------------------------------------------}
procedure cke_small(const k: complex; var kk, ek: complex);
  {-Return the complete elliptic integrals kk = K(k), ek = E(k); kk=INF if k^2=1}
var
  a,b,w,s: complex;
  f: double;
  done: boolean;
begin

  {Ref: Carlson [67], p.11: Algorithm for R_F(x,y,0) and R_G(x,y,0)}
  {and formula (42): K(k) = R_F(1-k^2,1,0), E(k) = 2*R_G(1-k^2,1,0)}

  {Starting values a_0 = sqrt(1-k^2), b_0 = 1}
  b := C_1;
  csqrt1mz2(k,a);
  if cabs(a)=0.0 then begin
    ek := C_1;
    kk.re := PosInf_d;
    kk.im := PosInf_d;
    exit;
  end;

  {Initialize sum factor for E(k)}
  s.re := 0.5*(a.re + b.re);
  s.im := 0.5*(a.im + b.im);
  csqr(s,s);
  f := 0.25;

  {AGM iteration}
  repeat
    {w = a*b}
    w.re := a.re*b.re - a.im*b.im;
    w.im := a.re*b.im + a.im*b.re;
    {a = (a+b)/2}
    a.re := 0.5*(a.re + b.re);
    a.im := 0.5*(a.im + b.im);
    {b = sqrt(a*b)}
    cx_sqrt(w.re, w.im, b.re, b.im);
    {loop until a-b is small of order sqrt(eps)}
    csub(a,b,w);
    done := (abs(w.re)<=sqrt_epsh*abs(a.re)) and (abs(w.im)<=sqrt_epsh*abs(a.im));
    {update sum factor for E(k)}
    {w = (a-b)^2}
    csqr(w,w);
    f := 2.0*f;
    s.re := s.re - f*w.re;
    s.im := s.im - f*w.im;
  until done;
  {K(k) = Pi/(a+b),  E(k) = K(k)*s}
  cadd(a,b,w);
  rdivc(Pi,w,kk);
  cmul(kk,s,ek);
end;


{---------------------------------------------------------------------------}
procedure cellke(const k: complex; var kk, ek: complex);
  {-Return the complete elliptic integrals kk = K(k), ek = E(k); kk=INF if k^2=1}
var
  a: double;
begin
  a := cabs(k);
  if a >= 40.0 then begin
    cellk(k, kk);
    ek_large(k, ek)
  end
  else if ((a >= 4.0) and (k.im=0)) then begin
    cellk(k, kk);
    ek_branch(k.re, ek)
  end
  else cke_small(k, kk, ek);
end;


{---------------------------------------------------------------------------}
procedure celle(const k: complex; var w: complex);
  {-Return w = E(k), the complete elliptic integral of the second kind}
var
  t: complex;
  a: double;
begin
  a := cabs(k);
  if a >= 40.0 then ek_large(k, w)
  else if (a >= 4.0) and (k.im=0) then ek_branch(k.re, w)
  else cke_small(k,t,w);
end;


{---------------------------------------------------------------------------}
procedure csurd(const z: complex; n: integer; var w: complex);
  {-Return the complex n'th root w = z^(1/n) with arg(w) closest to arg(z)}
var
  r,s,c,a: double;
  k: integer;
begin
  if n=1 then w:=z
  else if n=2 then csqrt(z,w)
  else begin
    r := nroot(cabs(z),n);
    if (z.re=0.0) and (n and 3 = 1) then begin
      {I^(4k+1) = I,  (-I)^(4k+1) = -I}
      if z.im < 0 then r := -r;
      w.re := 0.0;
      w.im := r;
    end
    else if (z.im=0.0) and odd(n) then begin
       if z.re < 0 then r := -r;
       w.re := r;
       w.im := 0.0;
    end
    else begin
      a := carg(z);
      k := round(a*(n-1)/TwoPi);
      sincos((k*TwoPi+a)/n,s,c);
      w.re := r*c;
      w.im := r*s;
    end;
  end
end;


{---------------------------------------------------------------------------}
function csgn(const z: complex): integer;
  {-Return the sign of z. Result = isign(z.re) if z.re<>0, isign(z.im) otherwise}
var
  s: integer;
begin
  if IsNaNd(z.re) then csgn := 0
  else begin
    s := isign(z.re);
    if s=0 then s := isign(z.im);
    csgn := s;
  end;
end;


{---------------------------------------------------------------------------}
procedure clngam_lanczos(const z: complex; var w: complex);
  {-Return lnGamma(z) using Lanczos sum, assumes z.re >= 1}
var
  k: integer;
  xr,tr: double;
  s,t: complex;
{Coefficients for g=9 from Paul Godfrey's  http://my.fit.edu/~gabdo/gammacoeff.txt}
const
  lgam = 9;
  lcoeffh: array[0..lgam+1] of THexDblW = (
             ($0001,$0000,$0000,$3FF0),  {+1.0000000000000002220}
             ($1E8A,$72BD,$5466,$40B6),  {+5.7164001882743414171E+3}
             ($544E,$F23E,$EFA6,$C0CC),  {-1.4815304267684139631E+4}
             ($8478,$134D,$E9BF,$40CB),  {+1.4291492776574785239E+4}
             ($FA4F,$0405,$CC29,$C0B8),  {-6.3481602176414589849E+3}
             ($5C6C,$E28A,$566E,$4094),  {+1.3016082860583219372E+3}
             ($5CCC,$23F6,$0B4F,$C05B),  {-1.0817670535143696497E+2}
             ($3B56,$68D7,$D877,$4004),  {+2.6056965056117560309}
             ($59AC,$DC13,$680D,$BF7E),  {-7.4234525102014163600E-3}
             ($8CA3,$989F,$E7E6,$3E6C),  {+5.3841364325095641522E-8}
             ($1056,$E5B1,$47EB,$BE31)); {-4.0235331412682361993E-9}
var
  lcoeff: array[0..lgam+1] of double absolute lcoeffh;
begin
  {Note for some arguments 1 <= z.re < 2 the w.im is off by 2Pi}
  {compute s = ln(Lanczos sum) = ln(c[0]+sum(c[k]/(z+k)))}
  s.re := lcoeff[0];
  s.im := 0.0;
  for k:=1 to lgam+1 do begin
    t.re := z.re + (k-1);
    t.im := z.im;
    rdivc(lcoeff[k],t,t);
    s.re := s.re + t.re;
    s.im := s.im + t.im;
  end;
  cln(s,s);

  {t = ln(z+lgam-0.5)}
  t.re := z.re + (lgam-0.5);
  t.im := z.im;
  cln(t,t);

  {w = -(z+lgam-0.5) + (z-0.5)*t + ln(sqrt(2Pi)) + s}
  xr := z.re - 0.5;
  tr := t.re - 1.0;
  w.re := tr*xr   - t.im*z.im + s.re + (LnSqrt2Pi-lgam);
  w.im := tr*z.im + t.im*xr   + s.im;

end;


{---------------------------------------------------------------------------}
procedure clngam1z(const z: complex; var w: complex);
  {-Return lnGamma(1+z) with power series, |z| < 2}
var
  s,p,d: complex;
  x,eps: double;
  n: integer;
const
  NMAX = 128;
const
  znhex: array[2..28] of THexDblW = ( {for n from 2 to 28 do evalf((-1)^n*(Zeta(n)-1)/n) od;}
           ($0FA6,$C4A6,$A34C,$3FD4),  {+3.2246703342411320303E-1}
           ($7607,$1A55,$3E00,$BFB1),  {-6.7352301053198102010E-2}
           ($8483,$AC7D,$1322,$3F95),  {+2.0580808427784546416E-2}
           ($F5F2,$C218,$404F,$BF7E),  {-7.3855510286739856768E-3}
           ($6C30,$EADB,$ADD6,$3F67),  {+2.8905103307415233593E-3}
           ($8E08,$C2BF,$8AC5,$BF53),  {-1.1927539117032610189E-3}
           ($96E9,$F863,$B36A,$3F40),  {+5.0966952474304245014E-4}
           ($2FC8,$C76D,$3FD4,$BF2D),  {-2.2315475845357938579E-4}
           ($D65A,$0F17,$127B,$3F1A),  {+9.9457512781808530980E-5}
           ($81EF,$BD7C,$8DE5,$BF07),  {-4.4926236738133142046E-5}
           ($EB02,$EE66,$80DC,$3EF5),  {+2.0507212775670691067E-5}
           ($2243,$63CE,$CBC9,$BEE3),  {-9.4394882752683967152E-6}
           ($4AAC,$39F3,$597A,$3ED2),  {+4.3748667899074881744E-6}
           ($9541,$B767,$1B2E,$BEC1),  {-2.0392157538013661897E-6}
           ($2F0F,$DEB2,$064C,$3EB0),  {+9.5514121304074193530E-7}
           ($FD2F,$D93C,$2600,$BE9E),  {-4.4924691987645661855E-7}
           ($7A4D,$B3F0,$76BB,$3E8C),  {+2.1207184805554664645E-7}
           ($8A97,$CBBF,$F5A6,$BE7A),  {-1.0043224823968099084E-7}
           ($0B0F,$C207,$9B93,$3E69),  {+4.7698101693639803983E-8}
           ($3EAC,$34DF,$62C7,$BE58),  {-2.2711094608943163504E-8}
           ($ADCD,$ACCF,$469D,$3E47),  {+1.0838659214896954593E-8}
           ($AEAD,$8447,$434A,$BE36),  {-5.1834750419700466442E-9}
           ($D2C3,$77FF,$55A8,$3E25),  {+2.4836745438024784752E-9}
           ($8D0E,$7925,$7B16,$BE14),  {-1.1921401405860911547E-9}
           ($C10C,$2B2F,$B15D,$3E03),  {+5.7313672416788622514E-10}
           ($E3E0,$9FAB,$F69A,$BDF2),  {-2.7595228851242333559E-10}
           ($434C,$A337,$4932,$3DE2)); {+1.3304764374244488820E-10}
begin
  {HMF[1], 6.1.33}
  x := 0.4227843350984671394; {1-EulerGamma}
  n := 2;
  s.re := z.re + 1.0;
  s.im := z.im;
  cln(s,s);
  s.re := z.re*x - s.re;
  s.im := z.im*x - s.im;
  eps  := 0.5*eps_d;
  {p = z^n}
  p.re := z.re;
  p.im := z.im;

  repeat
    cmul(p,z,p);
    {compute x = (-1)^n*(Zeta(n)-1)/n}
    if n<29 then x := double(znhex[n])
    else begin
      x := ldexpd(1,-n);
      x := x + x*x + exp3(-n); {1/2^n + 1/3^n + 1/4^n}
      if odd(n) then x := -x/n
      else x := x/n
    end;
    d.re := p.re*x;
    d.im := p.im*x;
    s.re := s.re + d.re;
    s.im := s.im + d.im;
    inc(n);
  until ((abs(d.re) <= eps*abs(s.re)) and (abs(d.im) <= eps*abs(s.im))) or (n>NMAX);
  {$ifdef debug}
    if n>NMAX then begin
      writeln('No convergence in clngam1z');
    end;
  {$endif}
  w.re := s.re;
  w.im := s.im;
end;


{---------------------------------------------------------------------------}
procedure clngamma(const z: complex; var w: complex);
  {-Return w = lnGamma(z), the principal branch of the log-Gamma function}
var
  t: complex;
  afix,s,c,y,sh,ch: double;
  si: integer;
begin
  {Note that lnGamma(z) is normally <> ln(Gamma(z)), real parts are }
  {equal but Im(Ln(Gamma(z)) is in [-Pi, Pi]. This function contains}
  {some guesswork about the multiples of 2*Pi which are to be added }
  {to w.im if z.re < 0, the reference was Maple's lnGAMMA function. }
  y := abs(z.im);
  if (y <= 0.75) and (abs(z.re-1)<=1.25) then begin
    {Use power series if z near pole at 0 or zeros at 1,2}
    {Near the zeros Lanczos gives only absolute accuracy.}
    if abs(z.re-1) <= 0.25 then begin
      {Code for z near 1}
      w.re := z.re - 1.0;
      w.im := z.im;
      clngam1z(w,w);
      exit;
    end;
    if abs(z.re) <= 0.25 then begin
      {Code for z near 0: lnGamma(z) = lnGamma(z+1) - ln(z)}
      {ln(z)}
      cln(z,t);
      w.re := z.re;
      w.im := z.im;
      {lnGamma(z+1)}
      clngam1z(w,w);
      w.re := w.re - t.re;
      w.im := w.im - t.im;
      exit;
    end;
    if abs(z.re-2.0) <= 0.25 then begin
      {Code for z near 2: lnGamma(z) = ln(z-1) + lnGamma(z-1)}
      {ln(z-1)}
      t.re := z.re - 1.0;
      t.im := z.im;
      cln(t,t);
      {compute lnGamma(z-1) with z-1 near 1}
      w.re := z.re - 2.0;
      w.im := z.im;
      clngam1z(w,w);
      w.re := w.re + t.re;
      w.im := w.im + t.im;
      exit;
    end;
  end;

  if (z.re >= 2.0) or (y>SUB) then begin
    {Here Lanczos seems OK, for z.re < 2 sometimes w.im is off by 2*Pi}
    clngam_lanczos(z,w);
  end
  else if (z.re > 1.0) then begin
    {Use lnGamma(z) = lnGamma(z+1) - ln(z) to make z.re > 2}
    cln(z,t);
    w.re := z.re + 1.0;
    w.im := z.im;
    clngam_lanczos(w,w);
    w.re := w.re - t.re;
    w.im := w.im - t.im;
  end
  else if z.re >= 0.0 then begin
    {Use lnGamma(z) = lnGamma(z+2) - ln(z+1) - ln(z) to make z.re > 2}
    t.re := z.re + 1.0;
    t.im := z.im;
    cmul(z,t,t);
    cln(t,t);
    w.re := z.re + 2.0;
    w.im := z.im;
    clngam_lanczos(w,w);
    w.re := w.re - t.re;
    w.im := w.im - t.im;
  end
  else begin
    {Use reflection formula HMF[1], 6.1.17: Gamma(z)Gamma(1-z)=Pi/sin(Pi*z}
    {and conjugation [1], 6.1.23: ln(Gamma(z)) = conj(ln(Gamma(conj(z))). }
    {Note: The real case could be simplified by using sfGamma, but then a }
    {large overhead from other special functions units would be included. }
    if y=0.0 then begin
      {Real case for z.re < 0}
      t.re := -z.re;
      t.im := 0.0;
      clngamma(t,t);
      s := z.re*sinPi(z.re);
      w.re := ln(abs(Pi/s)) - t.re;
      w.im := floord(z.re)*Pi;
    end
    else begin
      si := isign(z.im);
      {t := lngam(1-z)}
      t.re := 1.0 - z.re;
      t.im := -y;
      clngamma(t,t);
      {This is the 'magic' arg fix to make w.im compatible to Maple and}
      {Wolfram Alpha, see http://functions.wolfram.com/06.11.16.0002.01}
      afix := floord(0.5*(z.re+0.5))*TwoPi;
      {w := ln(sin(Pi*z))}
      sincosPi(z.re,s,c);
      if y >= 8.0 then begin
        {Here sinh(Pi*y)^2 > 1/eps_x, tanh(Pi*y) = 1}
        {HMF[1], 4.3.59: w.re ~ ln(sinh(y*Pi) = y*Pi - ln(2)}
        w.re := y*Pi - ln2;
        {HMF[1], 4.3.60: w.im ~ arctan(cot(z.re*Pi)tanh(Pi*y)), so}
        {w.im = arctan(cot(z.re*Pi)) = arctan(c/s) = arctan2(c,s).}
        w.im := arctan2(c,s);
      end
      else begin
        {Use full formulas from HMF[1] 4.3.59/60. An option would be:}
        {coshsinhmult(y*Pi, s,c, w.re,w.im); cln(w,w);}
        sinhcosh(y*Pi,sh,ch);
        w.re := ln(hypot(s,sh));
        w.im := arctan2(c*sh, s*ch);
      end;
      w.re :=  LnPi - w.re - t.re;
      w.im := (afix - w.im - t.im)*si;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure cgamma(const z: complex; var w: complex);
  {-Return the complex Gamma function w = Gamma(z)}
var
  rez: boolean;
begin
  rez := z.im=0.0;
  clngamma(z,w);
  cexp(w,w);
  if rez then w.im := 0.0;
end;


{---------------------------------------------------------------------------}
procedure cpsi(const z: complex; var w: complex);
  {-Return the complex digamma function w = psi(z), z <> 0,-1,-2...}
const
  NB = 6;                     {Bernoulli(2(k+1))/(2(k+1)), k=0..NB}
  B: array[0..NB] of double = (1/12,  -1/120,     1/252, -1/240,      {FIX311}
                               1/132, -691/32760, 1/12);
var
  u,v,x: complex;
  a: double;
const
  x0 = 12.0;
begin
  if z.re < 0.0 then begin
    {HMF [1], 6.3.7: psi(z) = psi(1-z) - Pi*cot(Pi*z)}
    x.re := 1.0-z.re;
    x.im := -z.im;
    cpsi(x,u);
    {compute v = Pi*cot(Pi*z)}
    if z.im=0.0 then begin
      sincosPi(z.re,v.re,v.im);
      {cot(x) = cos(x)/sin(x), will crash if frac(x)=0!}
      v.re := Pi*v.im/v.re;
      v.im := 0.0;
    end
    else if frac(z.re)=0.0 then begin
      {cot((-m + y*I)*Pi)= -I*coth(Pi*y))}
      v.re := 0.0;
      v.im := -Pi*coth(Pi*z.im);
    end
    else begin
      {Because cot has period Pi, use only the fractional part of z.re.  }
      {Note: There can be a large relative error, if z.re is not exactly }
      {representable as float, e.g. if z.re = -10 - 1e-8 then frac(z.re) }
      { = -1.0000000827403710e-8 with a relative error of 8.2740371E-8!  }
      {And then v.im, w.im will have relative errors of the same order.  }
      a := frac(z.re);
      {a < 0: if abs(a) > 0.5 take the next multiple of Pi. Note that }
      {a+1 generates no additional error according to Sterbenz' lemma.}
      if a < -0.5 then a := a + 1.0;
      x.re := Pi*a;
      x.im := Pi*z.im;
      ccot(x,v);
      v.re := Pi*v.re;
      v.im := Pi*v.im;
    end;
    w.re := u.re - v.re;
    w.im := u.im - v.im;
    exit;
  end;

  {This next code is analogous to the real digamma sfd_psi:}
  {Step 1: Make Re(x) >= x0 using HMF [1], 6.3.5}
  x.re := z.re;
  x.im := z.im;
  w.re := 0.0;
  w.im := 0.0;
  while x.re < x0 do begin
    cinv(x,v);
    cadd(w,v,w);
    x.re := x.re + 1.0;
  end;

  {Step 2: If re(x) >= x0, use asymptotic expansion from HMF [1], 6.3.18}
  {v = -1/(2x) - sum(B_2k/(2k*x^2k)}
  a := maxd(abs(x.re), abs(x.im));
  if a >= 1e10 then begin
    {Avoid overflow on square and useless polynomial evaluation}
    rdivc(-0.5, x, v);
  end
  else begin
    {TODO: optimize and use fewer terms depending on a?}
    csqr(x,u);
    cinv(u,u);
    cpolyr(u,B,NB+1,v);
    cmul(v,u,v);
    rdivc(-0.5, x, u);
    csub(u,v,v);
  end;

  {Step 3: psi(x) = ln(x) + v - w, with w from Step 1}
  cln(x,u);
  cadd(u,v,u);
  csub(u,w,w);
end;


{---------------------------------------------------------------------------}
procedure dilog_taylor(const z: complex; var w: complex);
  {-Compute w = dilog(z), |z| <= 0.5 with Taylor series}
var
  u,v,s,t: complex;
  f: double;
  n: integer;
const
  NMAX = 150;
begin
  {dilog(z) = sum(z^n/n^2, n=1..INF}
  u := z;
  t := u;
  s := t;
  for n:=2 to NMAX do begin
    cmul(t,u,t);
    f := sqr(n);
    v.re := s.re + t.re/f;
    v.im := s.im + t.im/f;
    if (v.re=s.re) and (v.im=s.im) then break;
    s.re := v.re;
    s.im := v.im;
  end;
  w.re := s.re;
  w.im := s.im;
end;


{---------------------------------------------------------------------------}
procedure dilog_bernsum(const z: complex; var w: complex);
  {-Compute w = dilog(z) with Debye function / Bernoulli sum}
  { Note: this procedure is only called if 0.5 < |z| < 2.0 !}
var
  s,y,t,v,f: complex;
  n: integer;
  b: double;
begin
  {Maximon[35], (4.3)}
  y.re := 1.0 - z.re;
  y.im := -z.im;
  cln(y,y);
  {The series converges for |ln(1-z)| < 2*Pi. Because z arguments near 1 are}
  {handled via transformation 1-z, the observed maximum of |ln(1-z)| < 3.22!}
  y.re := -y.re;
  y.im := -y.im;
  csqr(y,t);
  f.re := 1.0;
  f.im := 0.0;
  s.re := 1.0;
  s.im := 0.0;
  n := 0;
  while n < MaxB2nSmall do begin
    cmul(f,t,f);
    inc(n);
    b := n+n;
    b := sqr(b)+b; {=2n(2n+1)}
    f.re := f.re/b;
    f.im := f.im/b;
    {here f = (-ln(1-z))^(2n)/(2n+1)!}
    b := double(DAMath.B2nHex[n]);
    v.re := s.re + b*f.re;
    v.im := s.im + b*f.im;
    if (v.re=s.re) and (v.im=s.im) then break;
    s.re := v.re;
    s.im := v.im;
  end;
  {$ifdef debug}
    {Should not happen, observed maximum n was 34, MaxB2nSmall is 60!}
    if n>=MaxB2nSmall then writeln('dilog_bernsum: no convergence');
  {$endif}
  cmul(s,y,s);
  w.re := s.re - 0.25*t.re;
  w.im := s.im - 0.25*t.im;
end;


{---------------------------------------------------------------------------}
procedure cdilog(const z: complex; var w: complex);
  {-Return the principal branch of the complex dilogarithm w = -integral(ln(1-t)/t, t=0..z)}
var
  a: double;
  u,v,t: complex;
  isre,negi: boolean;
begin
  if (z.im=0.0) and (z.re=1.0) then begin
    {avoid evaluation of ln(0)}
    w.re := PiSqr/SIXX;
    w.im := 0.0;
    exit;
  end;
  isre := (z.im=0.0) and (z.re<=1.0);  {result is real}
  negi := (z.im=0.0) and (z.re>1.0);   {z on branch cut, w.im is negative}
  a := hypot(z.re, z.im);
  if a<=0.5 then begin
    {Use Taylor series: Maximon[35], (3.1)}
    dilog_taylor(z, w);
  end
  else if a>=2.0 then begin
    {Maximon[35], (3.2)}
    {compute t=ln^2(-z)}
    t.re := -z.re;
    t.im := -z.im;
    cln(t,t);
    csqr(t,t);
    {compute u = dilog(1/z)}
    cinv(z,u);
    dilog_taylor(u,u);
    {Complete transformation}
    w.re := (-0.5)*t.re - u.re - PiSqr/SIXX;
    w.im := (-0.5)*t.im - u.im;
  end
  else if hypot(1.0 - z.re,z.im) <= 0.5 then begin
    {Maximon[35], (3.3)}
    {t=1-z}
    t.re := 1.0 - z.re;
    t.im := -z.im;
    {compute v = ln(z)*ln(1-z)}
    cln(t,u);
    cln(z,v);
    cmul(u,v,v);
    {compute u = dilog(1-z)}
    dilog_taylor(t,u);
    {Complete transformation}
    w.re := PiSqr/SIXX - v.re - u.re;
    w.im := -v.im - u.im;
  end
  else begin
    {Maximon[35], (4.3)}
    dilog_bernsum(z,w);
    if negi then w.im := -w.im;
  end;
  if isre then begin
    {force real result}
    w.im := 0.0;
  end
  else if negi and (w.im > 0.0) then begin
    {Fix the sign of the imaginary part of the principal branch on}
    {the cut: Here we have Im(dilog(x+0*i)) = -Pi*ln(x) for x > 1.}
    w.im := -w.im;
  end;
end;


const
  x0e = 6.0;
  x1e = 27.25;

{---------------------------------------------------------------------------}
procedure cerfc(const z: complex; var w: complex);
  {-Return the complex complementary error function w = erfc(z) = 1-erf(z)}
const
  ph  = 9.4599403715183963963159646240033837464;
  pvh = 1.8419880743036792792631929248006767493e+01;

const
  npq = 18;
  nrs = 18;

  p: array[0..npq] of double = (
       2.1093083061644187538279122968913808152e-01,
       1.6713797949733065528971052035163045322e-01,
       1.0494102880451803704489103456267864462e-01,
       5.2209624806229062497556308453704817574e-02,
       2.0582158194044619069754225289969978299e-02,
       6.4293391618431334949721030322694803569e-03,
       1.5913908100149480106505036886507485785e-03,
       3.1212060500464898607481297780989148689e-04,
       4.8506855193831619356742051276480109919e-05,
       5.9733626677651815061875193953000538326e-06,
       5.8286735523223186734841193923199386322e-07,
       4.5066690471880700341073630648957059945e-08,
       2.7610653454261808871589140980422088532e-09,
       1.3403949680961254958324655659147790094e-10,
       5.1561314110869289752886735305964573946e-12,
       1.5716297853674992841622413214543994197e-13,
       3.7958724379391814374198930599398267061e-15,
       7.2645388829894728683842661597496303959e-17,
       1.1016397065175311063529497475707242525e-18);

  q: array[0..npq] of double = (
       2.9088820866572159615394846141476878557e-02,
       2.6179938779914943653855361527329190702e-01,
       7.2722052166430399038487115353692196393e-01,
       1.4253522224620358211543474609323670493e+00,
       2.3561944901923449288469825374596271631e+00,
       3.5197473248552313134627763831187023054e+00,
       4.9160107264506949750017289979095924762e+00,
       6.5449846949787359134638403818322976754e+00,
       8.4066692304393541288491105348868179031e+00,
       1.0501064332832549621157539457073153159e+01,
       1.2828170002158322390389127148391303444e+01,
       1.5387986238416672436543873608841268757e+01,
       1.8180513041607599759621778838423049098e+01,
       2.1205750411731104359622842837136644468e+01,
       2.4463698348787186236547065604982054867e+01,
       2.7954356852775845390394447141959280294e+01,
       3.1677725923697081821164987448068320749e+01,
       3.5633805561550895528858686523309176233e+01,
       3.9822595766337286513475544367681846745e+01);

  r: array[0..nrs] of double = (
       1.0857833597842664924141880507518822324e-01,
       1.9330394605384376865851659517702988226e-01,
       1.3634629684679999083131393208557974008e-01,
       7.6204577450604403788088248675471849578e-02,
       3.3748451951221171387030636574487677141e-02,
       1.1843000251292400843913857684704052043e-02,
       3.2930983812757205929921213829383580821e-03,
       7.2557574646222760967368809292892366817e-04,
       1.2667645289367856821531300675005701424e-04,
       1.7524438664090954799894562703988992670e-05,
       1.9210002539217758384351411001702958461e-06,
       1.6685753127194156512431242433234203465e-07,
       1.1484161468186100454021914187739990566e-08,
       6.2630784459909945721312231747915945961e-10,
       2.7065216004464856039123247342635784082e-11,
       9.2676629018209927646524369276088075961e-13,
       2.5145718215575387941885792873856026053e-14,
       5.4062104214337911880413049229210178595e-16,
       9.2099427616964263777001122276945386232e-18);

  s: array[0..nrs] of double = (
       0.0,
       1.1635528346628863846157938456590751423e-01,
       4.6542113386515455384631753826363005692e-01,
       1.0471975511965977461542144610931676281e+00,
       1.8616845354606182153852701530545202277e+00,
       2.9088820866572159615394846141476878557e+00,
       4.1887902047863909846168578443726705123e+00,
       5.7014088898481432846173898437294681972e+00,
       7.4467381418424728615410806122180809107e+00,
       9.4247779607693797153879301498385086526e+00,
       1.1635528346628863846157938456590751423e+01,
       1.4078989299420925253851105532474809222e+01,
       1.6755160819145563938467431377490682049e+01,
       1.9664042905802779900006915991638369905e+01,
       2.2805635559392573138469559374917872789e+01,
       2.6179938779914943653855361527329190702e+01,
       2.9786952567369891446164322448872323643e+01,
       3.3626676921757416515396442139547271612e+01,
       3.7699111843077518861551720599354034610e+01);

var
  x,y,t,u,v,a: complex;
  i: integer;
  re0,im0: boolean;
begin
  {Ref: archive gamerf.zip, file gamerf2a.f, function cqerfc(x)}
  {Gamma / Error Functions in Fortran, (C) 1996 Takuya OOURA   }
  {http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html}

  {Re z=0 -> Re w = 1, Im z=0 -> Im w = 0}
  re0 := z.re=0.0;
  im0 := z.im=0.0;
  if im0 then begin
    w.im := 0.0;
    if re0 then begin
      w.re := 1.0;
      exit;
    end
    else if z.re < -x0e  then begin
      w.re := 2.0;
      exit;
    end
    else if z.re > x1e then begin
      w.re := 0.0;
      exit;
    end;
  end;

  csqr(z,y);
  v.re := 0.0;
  v.im := 0.0;
  u.im := y.im;
  t.re := -y.re;
  t.im := -y.im;
  cexp(t, x);
  cmul(x,z,x);

  if abs(z.re)+abs(z.im) <= ph then begin
    t.re := pvh*z.re;
    t.im := pvh*z.im;
    cexp(t,a);
    if a.re >= 0.0 then begin
      for i:=npq downto 0 do begin
        u.re := y.re + q[i];
        rdivc(p[i], u, t);
        v.re := v.re + t.re;
        v.im := v.im + t.im;
      end;
      a.re := a.re + 1.0;
      rdivc(2.0,a,t);
    end
    else begin
      for i:=nrs downto 0 do begin
        u.re := y.re + s[i];
        rdivc(r[i], u, t);
        v.re := v.re + t.re;
        v.im := v.im + t.im;
      end;
      a.re := a.re - 1.0;
      rdivc(-2.0,a,t);
    end;
    cmul(v,x,u);
    w.re := u.re + t.re;
    w.im := u.im + t.im;
  end
  else begin
    for i:=npq downto 0 do begin
      u.re := y.re + q[i];
      rdivc(p[i], u, t);
      v.re := v.re + t.re;
      v.im := v.im + t.im;
    end;
    if z.re >= 0.0 then cmul(x,v,w)
    else begin
      cmul(x,v,w);
      w.re := w.re + 2.0;
    end;
  end;
  if re0 then w.re := 1.0;
  if im0 then w.im := 0.0;
end;


{---------------------------------------------------------------------------}
procedure cerf(const z: complex; var w: complex);
  {-Return the complex error function w = erf(z) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..z)}
var
  u: complex;
const
  nt = 8;
  cts: array[0..nt-1] of double = (
         1.1283791670955125738961589031215451717e+00,
        -3.7612638903183752463205296770718172390e-01,
         1.1283791670955125738961589031215451717e-01,
        -2.6866170645131251759432354836227265993e-02,
         5.2239776254421878421118467737108572763e-03,
        -8.5483270234508528325466583569814028158e-04,
         1.2055332981789664251027338708563516792e-04,
        -1.4925650358406250977462419353459592218e-05);
begin
  {Ref: archive gamerf.zip, file gamerf2a.f, function cqerf(x)}
  {Gamma / Error Functions in Fortran, (C) 1996 Takuya OOURA  }
  {http://www.kurims.kyoto-u.ac.jp/~ooura/gamerf.html}
  if abs(z.re)+abs(z.im) <= 0.125 then begin
    {Use Taylor series}
    csqr(z,u);
    cpolyr(u,cts,nt,u);
    cmul(z,u,w);
  end
  else begin
    {Use erf(z) = 1-erfc(z)}
    if z.re >= 0.0 then begin
      cerfc(z,u);
      w.re := 1.0 - u.re;
      w.im := -u.im;
    end
    else begin
      u.re := -z.re;
      u.im := -z.im;
      cerfc(u,w);
      w.re := w.re - 1.0;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure rsncndn(x,mc: double; var sn,cn,dn: double);
  {-Return the Jacobi elliptic functions sn,cn,dn for argument x and complementary parameter mc}
const
  NA = 17;
var
  a,b,c,d: double;
  i,l: integer;
  bo: boolean;
  m,n: array[0..NA] of double;
begin
  {This is a copy of sfc_sncndn from sfellint.pas}
  {Handle degenerated cases x=0, mc=0 or mc=1 using HMF[1], Table 16.6}
  if x=0.0 then begin
    sn := 0.0;
    cn := 1.0;
    dn := 1.0;
  end
  else if mc=0.0 then begin
    {mc=0, ie m=1}
    sn := tanh(x);
    cn := sech(x);
    dn := cn;
  end
  else if mc=1.0 then begin
    {mc=1, ie m=0}
    sincos(x,sn,cn);
    dn := 1.0;
  end
  else begin
    bo := mc<0.0;
    d  := 1.0-mc; {Avoid compiler warning, move in front of 'if'}
    if bo then begin
      {Here m=1-mc>1, apply Jacobi real transformation, HMF[1] 16.11}
      mc := -mc/d;
      d  := sqrt(d);
      x  := d*x;
    end;
    {AGM process is done with mc and 1. Maximum index l for mc=succ(0)}
    {and ca=eps_x is 16. So NA=12 is too small for extended precision.}
    {Some l values: mc=1e-99 -> l=10,  mc=1e-10 -> l=6,  mc=0.1 -> l=4}
    a  := 1.0;
    dn := 1.0;
    for i:=0 to NA do begin
      l   := i;
      m[i]:= a;
      mc  := sqrt(mc);
      n[i]:= mc;
      c   := 0.5*(a+mc);
      if abs(a-mc) <= sqrt_epsh*a then break;
      mc := a*mc;
      a  := c;
    end;
    sincos(c*x,sn,cn);
    if sn<>0.0 then begin
      a := cn/sn;
      c := a*c;
      for i:=l downto 0 do begin
        b := m[i];
        a := c*a;
        c := dn*c;
        dn:= (n[i]+a)/(b+a);
        a := c/b;
      end;
      a := 1.0/sqrt(sqr(c)+1.0);
      if sn<0.0 then sn := -a else sn := a;
      cn := c*sn;
    end;
    if bo then begin
      a  := dn;
      dn := cn;
      cn := a;
      sn := sn/d;
    end;
  end
end;


{---------------------------------------------------------------------------}
procedure csncndn(const z: complex; k: double; var sn,cn,dn: complex);
  {-Return the Jacobi elliptic functions sn(z,k), cn(z,k), dn(z,k)}
  { for complex argument z and real modulus k}
var
  s1,s2,c1,c2,d1,d2,q,m,mc: double;
begin
  m  := sqr(k);
  mc := (1.0-k)*(1.0+k);
  if z.im=0.0 then begin
    rsncndn(z.re, mc, sn.re, cn.re, dn.re);
    sn.im := 0.0;
    cn.im := 0.0;
    dn.im := 0.0;
  end
  else if z.re=0.0 then begin
    rsncndn(z.im, m, s2, c2, d2);
    sn.re := 0.0;
    sn.im := s2/c2;
    cn.re := 1.0/c2;
    cn.im := 0.0;
    dn.re := d2/c2;
    dn.im := 0.0;
  end
  else begin
    rsncndn(z.re, mc, s1, c1, d1);
    rsncndn(z.im,  m, s2, c2, d2);
    q := sqr(c2) + sqr(k*s1*s2);
    sn.re := s1*d2/q;         {http://functions.wolfram.com/09.36.19.0001.01}
    sn.im := s2*c1*c2*d1/q;   {http://functions.wolfram.com/09.36.19.0002.01}
    cn.re := c1*c2/q;         {http://functions.wolfram.com/09.26.19.0001.01}
    cn.im := -s1*s2*d1*d2/q;  {http://functions.wolfram.com/09.26.19.0002.01}
    dn.re := c2*d1*d2/q;      {http://functions.wolfram.com/09.29.19.0001.01}
    dn.im := -m*s1*s2*c1/q;   {http://functions.wolfram.com/09.29.19.0002.01}
  end;
end;


{---------------------------------------------------------------------------}
procedure csn(const z: complex; k: double; var sn: complex);
  {-Return the Jacobi elliptic function sn(z,k)}
var
  cn,dn: complex;
begin
  csncndn(z, k, sn,cn,dn);
end;


{---------------------------------------------------------------------------}
procedure ccn(const z: complex; k: double; var cn: complex);
  {-Return the Jacobi elliptic function cn(z,k)}
var
  sn,dn: complex;
begin
  csncndn(z, k, sn,cn,dn);
end;


{---------------------------------------------------------------------------}
procedure cdn(const z: complex; k: double; var dn: complex);
  {-Return the Jacobi elliptic function dn(z,k)}
var
  sn,cn: complex;
begin
  csncndn(z, k, sn,cn,dn);
end;


{---------------------------------------------------------------------------}
procedure lwk_init_val(const z: complex; k: integer; var w: complex);
  {-Compute the initial value w for LambertW(k,z) iteration}
var
  u,v,p: complex;
  d: double;
const
  eh  = 2.7177734375;
  el  = 0.5083909590452353603e-3; {e=eh+el}
  e   = eh+el;
  em1 = 0.3678794411714423216; {1/e}
const
  c1 = 0.3333333333333333333;  {= 1/3}
  c2 = 0.1527777777777777778;  {= 11/72}
  c3 : complex = (re: 23.2837992410100535; im : 202.0203166004481079690);
  c4 : complex = (re: 67.01105156624095; im: 33.1655151845347645);
  c5 : complex = (re: -34.46206; im: 21.259442);
  c6 : complex = (re: -15.23103; im: 10.629721);
  c7 = 2.48369826539555787;
  c8 = 1.654368;
begin
  {Basic numerical method from Corless et.al. [23], Pade approximation}
  {from "C++ code for the complex Lambert W function" by Istvan Mezo  }
  {available from https://sites.google.com/site/istvanmezo81/others   }
  if (abs(k)<2) and (hypot(z.re+em1,z.im) <= 0.36) then begin
    {z near -1/e: Use branch point series, see [23], (4.22)}
    {p = sqrt(2*(e*z + 1))}
    p.re := 2*e*z.re + 2.0;
    p.im := 2*e*z.im;
    csqrt(p,p);
    {u = p^2,  v = p^3}
    csqr(p,u);
    cmul(p,u,v);
    u.re := c1*u.re;
    u.im := c1*u.im;
    v.re := c2*v.re;
    v.im := c2*v.im;
    if k=0 then begin
      {w = -1 + p - 1/3*p^2 + 11/72*p^3}
      w.re := -1.0 + p.re - u.re + v.re;
      w.im :=        p.im - u.im + v.im;
      exit;
    end
    else if (k = 1) and (z.im < 0.0) then begin
      {w = -1 - p - 1/3*p^2 - 11/72*p^3}
      w.re := -1.0 - p.re - u.re - v.re;
      w.im :=      - p.im - u.im - v.im;
      exit;
    end
    else if (k = -1) and (z.im > 0.0) then begin
      {w = -1 - p - 1/3*p^2 - 11/72*p^3}
      w.re := -1.0 - p.re - u.re - v.re;
      w.im :=      - p.im - u.im - v.im;
      exit;
    end;
    {else fall through and use other approximation}
  end;

  d := hypot(z.re-0.5, z.im);
  if (k=0) and (d<1.2) then begin
      {(1,1) Pade approximant for W(0,z)}
      {w = (0.35173371 * (0.1237166 + 7.061302897*z)) / (2 + 0.827184*(1 + 2*z))}
      u.re := c7*z.re + 0.043515298706586;
      u.im := c7*z.im;
      v.re := c8*z.re + 2.827184;
      v.im := c8*z.im;
      cdiv(u,v,w);
      exit;
  end;
  if (k=-1) and (d<0.5) then begin
      {(1,1) Pade approximant for W(-1,z)}
      {w = -(((2.2591588985 + 4.22096*I) * ((-14.073271 - 33.767687754*I)*z -
             (12.7127 - 9.071643*I) * (1 + 2*z))) / (2 - (17.23103 - 10.629721*I)*(1 + 2*z)))}
      cmul(c3,z,u);
      cadd(u,c4,u);
      cmul(c5,z,v);
      cadd(v,c6,v);
      cdiv(u,v,w);
      exit;
  end;

  {Starting value from the asymptotic approximation [23], (4.20)}
  {w = ln(z) + 2*Pi*k*i - log(log(z) + 2*Pi*k*i))}
  cln(z,u);
  u.im := u.im + TwoPi*k;
  cln(u,v);
  csub(u,v,w);
end;


{---------------------------------------------------------------------------}
procedure cLambertWk(k: integer; const z: complex; var wk: complex);
  {-Return the k'th branch wk = W(k,z) of the Lambert W function}
var
  w, u, v, x, f0, f1, f2: complex;
  i: integer;
const
  eps  = 1e-10;
  IMAX = 30;
const
  eh = 2.7177734375;
  el = 0.5083909590452353603e-3; {e=eh+el}
  e  = eh+el;
begin

  if z.im=0.0 then begin
    {Special cases}
    if z.re=0.0 then begin
      if k=0 then wk.re := 0.0 else wk.re := NegInf_d;
      wk.im := 0.0;
      exit;
    end;
    if (abs(e*z.re + 1.0) <= 4*eps_d) and ((k=0) or (k=-1)) then begin
      {close to -1/e}
      wk.re := 1.0;
      wk.im := 0.0;
      exit;
    end;
  end;

  {get initial value for iteration}
  lwk_init_val(z,k,w);

  {Halley loop}
  for i:=1 to IMAX do begin
    {x = exp(w); f0 = f(w) = w*exp(w)}
    cexp(w,x);
    cmul(x,w,f0);

    {f1 = f'(w) = exp(w) + w*exp(w) = x + f0; }
    cadd(x,f0,f1);

    {f2 = f''(w) = 2*exp(w) + w*exp(w) = x + f1;}
    cadd(x,f1,f2);

    {u = [f(w) - z]*f'(w)}
    csub(f0,z,x);
    cmul(x,f1,u);

    {v = f'(w)^2 - (f(w)-z)*f''(w)/2 }
    cmul(x,f2,v);
    csqr(f1,x);
    v.re := x.re - 0.5*v.re;
    v.im := x.im - 0.5*v.im;

    {x = u/v}
    cdiv(u,v,x);

    {w = w - x}
    w.re := w.re - x.re;
    w.im := w.im - x.im;

    if abs(x.re) + abs(x.im) <= eps*(abs(w.re) + abs(w.im)) then break;
  end;

  wk.re := w.re;
  wk.im := w.im;
end;


{---------------------------------------------------------------------------}
procedure cLambertW(const z: complex; var w: complex);
  {-Return the principal branch w = W(z) = W(0,z) of the Lambert W function}
begin
  cLambertWk(0, z, w);
end;


{---------------------------------------------------------------------------}
procedure crgamma(const z: complex; var w: complex);
  {-Return the reciprocal Gamma function w = 1/Gamma(z)}
const
  tsa: array[0..5] of double = (
       +0.0,
       +1.0,                         {z}
       +0.5772156649015328606,       {z^2}
       -0.6558780715202538811,       {z^3}
       -0.4200263503409523553e-1,    {z^4}
       +0.1665386113822914895);      {z^5}
begin
  if (abs(z.re)<=1e-5) and (abs(z.im)<=1e-5) then begin
    {Maclaurin series}
    cpolyr(z, tsa, 6, w);
  end
  else begin
    clngamma(z,w);
    w.re := -w.re;
    w.im := -w.im;
    cexp(w,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure crstheta(const z: complex; var w: complex);
  {-Return the Riemann-Siegel function w = theta(z)}
const
  TS: array[0..5] of double = (
         -2.68609170961283279111647874873,
         +2.69432791536535250694454098289,
         -6.40218109079112857644745498764,
        +18.2859525169638076731561890870,
        -56.8889181687769276420032504557,
       +186.181822000864703629552237858);
var
  u,v: complex;
begin
  if hypot(z.re,z.im) <= 0.015625 then begin
    csqr(z,u);
    cpolyr(u,TS,6,v);
    cmul(v,z,w);
  end
  else begin
    {http://functions.wolfram.com/10.03.02.0001.01}
    {w = -ln(Pi)*z/2 -i/2*(lngamma(1/4 +iz/2) -lngamma(1/4 -iz/2))}
    {u = z/2}
    u.re := 0.5*z.re;
    u.im := 0.5*z.im;
    {w = lngamma(1/4 +iz/2)}
    v.re := 0.25 - u.im;
    v.im := u.re;
    clngamma(v, w);
    {w = lngamma(1/4 -iz/2)}
    v.re := 0.25 + u.im;
    v.im := -u.re;
    clngamma(v, v);
    {v = lngamma(1/4 +iz/2) -lngamma(1/4 -iz/2)}
    v.re := w.re - v.re;
    v.im := w.im - v.im;
    {w = -ln(Pi)*z/2 -i/2*(lngamma(1/4 +iz/2) -lngamma(1/4 -iz/2))}
    w.re :=  0.5*v.im - lnPi*u.re;
    w.im := -0.5*v.re - lnPi*u.im;
  end;
end;


{---------------------------------------------------------------------------}
procedure czeta_pos(const s: complex; var w: complex);
  {-Return w=Zeta(s), s.re >= 0}
var
  a,c,z,zn,t,u: complex;
  k,n,m: longint;
  x,y,eps: double;
begin

  {estimate initial n}
  if s.re >= 1.0 then begin
    x := 10 + 0.25*abs(s.im)/s.re;
    if x>100 then n := 100
    else n := round(x);
  end
  else n := 32;

  {Ref: A. Banuelos, R.A. Depine [81]}
  m := 2;
  eps := 2*eps_d;
  c := c_1;
  zn := c;
  repeat
    {zn = sum(k^(-s), k=1..N-1}
    for k:=m to n do begin
      x := ln1p(-1/k);
      t.re := x*s.re;
      t.im := x*s.im;
      cexp(t,t);
      cmul(c,t,c);
      cadd(zn,c,zn);
    end;
    n := n+1;
    x := ln(n);

    {N^(-s)/2}
    t.re := -x*s.re;
    t.im := -x*s.im;
    cexp(t,t);
    z.re := zn.re + 0.5*t.re;
    z.im := zn.im + 0.5*t.im;

    {N^(1-s)/(s-1)}
    u.re := s.re - 1.0;
    u.im := s.im;
    t.re := -x*u.re;
    t.im := -x*u.im;
    cexp(t,t);
    cdiv(t,u,u);
    cadd(z,u,z);

    {0.5*s*B_1/N^(s+1)}
    t.re := x*(s.re+1);
    t.im := x*s.im;
    cexp(t,t);
    cdiv(s,t,a);
    x := 0.5*Bern2n[1];
    a.re := x*a.re;
    a.im := x*a.im;
    cadd(z,a,z);

    y := sqr(int(n));
    for k:=1 to 30 do begin
      t.re := s.re + 2*k;
      t.im := s.im;
      u.re := t.re - 1.0;
      u.im := t.im;
      cmul(u,t,t);
      x := Bern2n[k+1]/Bern2n[k]/(2*k+2)/(2*k+1)/y;
      t.re := x*t.re;
      t.im := x*t.im;
      cmul(a,t,a);
      cadd(z,a,z);
      x := hypot(s.re + 2*k+1, s.im)*cabs(a);
      if x < eps then begin
        w := z;
        {writeln('n=',n);}
        exit;
      end;
    end;
    m := n;
    n := n+64;
  until false;
end;


{---------------------------------------------------------------------------}
procedure czeta_neg(const s: complex; var w: complex);
  {-Return w=Zeta(s), s.re < 0.5}
var
  a,t,u: complex;
const
  smin = 2e-3;  {double = 2e-3}
  TCA: array[0..5] of double = (
         -0.500000000000000000000000000000,
         -0.918938533204672741780329736405,
         -1.00317822795429242560505001337,
         -1.00078519447704240796017680223,
         -0.999879299500571164957800813655,
         -1.00000194089632045603779988198);
begin
  if cabs(s) < smin then begin
    {Maclaurin series}
    cpolyr(s, TCA, 6, w);
  end
  else begin
    {Zeta(s) = 2^s*Pi^(s-1)*sin(Pi/2*s)*gamma(1-s)*zeta(1-s);}
    u.re := 1.0-s.re;
    u.im := -s.im;
    cexp2(s,a);
    t.re := lnPi*u.re;
    t.im := lnPi*u.im;
    cexp(t,t);
    cdiv(a,t,a);
    t.re := 0.5 * s.re;
    t.im := 0.5 * s.im;
    csinpi(t,t);
    cmul(a,t,a);
    cgamma(u,t);
    cmul(a,t,a);
    czeta_pos(u,t);
    cmul(a,t,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure czeta(const s: complex; var w: complex);
  {-Return the complex Riemann Zeta function, w=Zeta(s), s <> 1}
var
  x: double;
begin
  x := s.re;
  if s.im=0.0 then begin
    if x=0.0 then begin
      w.re := -0.5;
      w.im := 0.0;
      exit;
    end
    else if x=1 then begin
      w.re := NaN;
      w.im := NaN;
      exit;
    end
    else if (x<0.0) and (frac(0.5*x)=0.0) then begin
      {s.re = -2n}
      w.re := 0.0;
      w.im := 0.0;
      exit;
    end;
  end;
  if x >= 0.375 then begin
    czeta_pos(s,w);
  end
  else begin
    czeta_neg(s,w);
  end;
end;


{---------------------------------------------------------------------------}
procedure expint_taylor(const z: complex; ei: boolean; var w: complex);
  {-Common Taylor series for E1/Ei}
var
  k: integer;
  term,s,t: complex;
const
  KMAX = 200;
begin
  if ei then begin
    cln(z,s);
    s.im := isign(z.im)*abs(s.im);
  end
  else begin
    cneg(z,t);
    cln(t,s);
    if (t.im=0.0) and (t.re<0) then s.im := Pi;
  end;
  s.re := s.re + Eulergamma;
  term := c_1;
  for k:=1 to KMAX do begin
    t.re := z.re/k;
    t.im := z.im/k;
    cmul(term,t,term);
    t.re := term.re/k;
    t.im := term.im/k;
    cadd(s,t,s);
    if cabs(t) <= eps_d*cabs(s) then break;
  end;
  w := s;
end;


{---------------------------------------------------------------------------}
procedure expint_asymp(const z: complex; var w: complex; var OK: boolean);
  {-Common asymptotic expression for E1/Ei}
var
  k: integer;
  term,s,t: complex;
const
  KMAX = 100;
begin
  OK := false;
  s.re := 0.0;
  s.im := 0.0;
  cexp(z,term);
  cdiv(term,z,term);
  for k:=1 to KMAX do begin
    cadd(s,term,s);
    if cabs(term) <= eps_d*cabs(s) then begin
      OK := true;
      break;
    end;
    rdivc(k,z,t);
    cmul(term,t,term);
  end;
  if OK then begin
    w.re := s.re;
    w.im := s.im;
  end;
end;


{---------------------------------------------------------------------------}
procedure expint_cfrac(const z: complex; var w: complex; var OK: boolean);
  {-Common continued fraction for E1/Ei with modified Lentz method, z<>1}
var
  a,b,c,d,f,t,h: complex;
  tiny,tol: double;
  k: longint;
const
  MAXIT = 1000;
begin
  OK   := false;
  tol  := eps_d;
  tiny := MinDouble;

  {Note the AMath does not include sign(z.im)*i*Pi and exp(z) in the continued}
  {fraction calculation (the are applied afterwards), because this would give}
  {rather inaccurate results, for large negative z.re}

  f.re := 1.0 - z.re;
  f.im :=     - z.im;

  c := f;
  d := c_0;

  {Imaginary parts of a and b are fix}
  a.im := 0.0;
  b.im := -z.im;

  for k:=1 to MAXIT do begin
    a.re := -sqr(k);
    b.re := (2*k+1) - z.re;
    {d := b + a*d}
    cmul(a,d,h);
    cadd(b,h,d);
    {c := b + a/c}
    cdiv(a,c,h);
    cadd(b,h,c);
    if (c.re=0) and (c.im=0) then c.re := tiny;
    if (d.re=0) and (d.im=0) then d.re := tiny;
    cinv(d,d);
    cmul(c,d,t);
    cmul(f,t,f);
    if hypot(t.re-1.0,t.im) <= tol then begin
      OK := true;
      break;
    end;
  end;
  if OK then begin
    cexp(z,h);
    cdiv(h,f,h);
    w.re := -h.re;
    w.im := -h.im;
  end;
end;


{---------------------------------------------------------------------------}
procedure ei_asymp(const z: complex; var w: complex; var OK: boolean);
var
  siz: integer;
begin
  siz := isign(z.im);
  expint_asymp(z, w, OK);
  if OK then begin
    w.im := w.im + siz*Pi;
  end;
end;

{---------------------------------------------------------------------------}
procedure ei_cfrac(const z: complex; var w: complex; var OK: boolean);
  {-Ei(z) via continued fraction with modified Lentz method, z<>1}
var
  siz: integer;
begin
  siz := isign(z.im);
  expint_cfrac(z,w,OK);
  if OK then begin
    w.im := w.im + siz*Pi;
  end;
end;


{---------------------------------------------------------------------------}
procedure cei(const z: complex; var w: complex);
  {-Return the complex exponential integral Ei(z), z <> 0}
var
  a: double;
  OK: boolean;
const
  ca = 41.0;
begin
  a := cabs(z);
  if a=0.0 then begin
    w.re := NegInfinity;
    w.im := NegInfinity;
    exit;
  end;
  if a > ca then begin
    Ei_Asymp(z,w,OK);
    if OK then exit;
  end;
  if (a > 5.0) and ((z.re < 0.0) or (abs(z.im) >= 2.0)) then begin
    ei_cfrac(z,w,OK);
    if OK then exit;
  end;
  expint_taylor(z,true,w);
end;


{---------------------------------------------------------------------------}
procedure ce1(const z: complex; var w: complex);
  {-Return the complex exponential integral E1(z), z <> 0}
var
  a: double;
  t: complex;
  OK: boolean;
label
  adjust;
const
  ca = 41.0;
begin
  a := cabs(z);
  if a=0.0 then begin
    w.re := Infinity;
    w.im := Infinity;
    exit;
  end;
  cneg(z,t);
  if a > ca then begin
    expint_asymp(t,w,OK);
    if OK then goto adjust;
  end;
  if (a > 4.0) and ((z.re > 0.0) or (abs(z.im) >= 2.0)) then begin
    expint_cfrac(t,w,OK);
    if OK then goto adjust;
  end;
  expint_taylor(t, false, w);

adjust:
  w.re := -w.re;
  {Adjust w.im to be -Pi if z is real and < 0, note t = -z}
  if (t.im=0.0) and (t.re > 0.0) then w.im := -Pi
  else w.im := -w.im
end;


{---------------------------------------------------------------------------}
procedure cli(const z: complex; var w: complex);
  {-Return the complex logarithmic integral w = li(z) = Ei(ln(z)), z<>1}
begin
  if cabs(z)=0.0 then w := c_0
  else begin
    cln(z,w);
    cei(w,w);
  end;
end;


end.
