unit SpecFunD;

{Special functions (double precision)}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


uses
  sdBasic;   {Basic common code}

(*************************************************************************

 DESCRIPTION   :  Special functions (double precision)

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  See damath_info.txt/references


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------

 1.00.00  05.02.13  W.Ehrhardt  Initial BP7 version from AMath.SpecFun
 1.00.01  20.02.13  we          SpecFunD_Version

 1.01.00  15.03.13  we          Remaining Jacobi elliptic functions
 1.01.01  25.03.13  we          Remaining inverse Jacobi elliptic functions
 1.01.02  28.03.13  we          trilogarithm function trilog(x)

 1.02.00  11.04.13  we          fermi_dirac_p15
 1.02.01  13.04.13  we          k changed to longint in poisson_pmf
 1.02.02  14.04.13  we          rayleigh_pdf/cdf/inv
 1.02.03  18.04.13  we          beta3
 1.02.04  19.04.13  we          maxwell_pdf/cdf/inv
 1.02.05  20.04.13  we          evt1_pdf/cdf/inv
 1.02.06  01.05.13  we          FresnelC/FresnelS

 1.03.00  09.05.13  we          Airy/Scorer functions sfd_airy_gi/hi
 1.03.01  27.05.13  we          Hypergeometric functions

 1.04.00  16.06.13  we          Whittaker functions
 1.04.01  28.06.13  we          hyperg_0F1
 1.04.02  01.07.13  we          hyperg_0F1r

 1.05.00  15.08.13  we          kumaraswamy_pdf/cdf/inv
 1.05.01  18.08.13  we          erf_p, erf_q, erf_z

 1.06.00  12.09.13  we          Bernoulli polynomials bernpoly
 1.06.01  18.09.13  we          chebyshev_v/w
 1.06.02  24.09.13  we          cosint, sinint
 1.06.03  25.09.13  we          BatemanG

 1.07.00  19.10.13  we          moyal_pdf/cdf/inv
 1.07.01  23.10.13  we          EllipticKim, EllipticECim

 1.08.00  26.12.13  we          fermi_dirac_p25

 1.10.00  02.05.14  we          harmonic
 1.10.01  03.05.14  we          harmonic2
 1.10.02  04.05.14  we          zipf_pmf/cdf

 1.11.00  24.05.14  we          li_inv
 1.11.01  30.05.14  we          Fix zipf_pmf/cdf function type
 1.11.02  31.05.14  we          lobachevsky_c/s
 1.11.03  04.06.14  we          fibpoly, lucpoly
 1.11.04  05.06.14  we          levy_pdf/cdf/inv

 1.12.00  22.06.14  we          lngamma_inv
 1.12.01  25.06.14  we          psi_inv
 1.12.02  03.07.14  we          ibeta_inv
 1.12.03  12.07.14  we          CylinderD, CylinderU, HermiteH

 1.13.00  17.08.14  we          Catalan function
 1.13.01  19.08.14  we          CylinderV

 1.16.00  23.12.14  we          invgamma_pdf/cdf/inv
 1.16.01  26.12.14  we          logseries_cdf/pmf: logarithmic (series) distribution
 1.16.02  03.01.15  we          wald_cdf/pdf/inv
 1.16.03  09.01.15  we          Euler numbers

 1.17.00  25.01.15  we          ei_inv
 1.17.01  28.03.15  we          comp_ellint_d, ellint_d

 1.18.00  04.06.15  we          kepler
 1.18.01  10.06.15  we          struve_l/h
 1.18.02  20.06.15  we          scaled Airy functions airy_ais, airy_bis
 1.18.03  22.06.15  we          Carlson ell_rg

 1.21.00  26.09.15  we          lemniscate functions sin_lemn, cos_lemn

***************************************************************************)

(*-------------------------------------------------------------------------
 (C) Copyright 2010-2015 Wolfgang Ehrhardt

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
  SpecFunD_Version = SD_BaseVersion;


{#Z+}
{---------------------------------------------------------------------------}
{--------------------------- Bessel functions ------------------------------}
{---------------------------------------------------------------------------}
{#Z-}

function airy_ai(x: double): double;
  {-Return the Airy function Ai(x)}

function airy_aip(x: double): double;
  {-Return the Airy function Ai'(x)}

function airy_ais(x: double): double;
  {-Return the scaled Airy function Ai(x) if x <= 0, Ai(x)*exp(2/3*x^1.5) for x > 0}

function airy_bi(x: double): double;
  {-Return the Airy function Bi(x)}

function airy_bip(x: double): double;
  {-Return the Airy function Bi'(x)}

function airy_bis(x: double): double;
  {-Return the scaled Airy function Bi(x) if x <= 0, Bi(x)*exp(-2/3*x^1.5) for x > 0}

function airy_gi(x: double): double;
  {-Return the Airy/Scorer function Gi(x) = 1/Pi*integral(sin(x*t+t^3/3), t=0..INF)}

function airy_hi(x: double): double;
  {-Return the Airy/Scorer function Hi(x) = 1/Pi*integral(exp(x*t-t^3/3), t=0..INF)}

function bessel_i0(x: double): double;
  {-Return I0(x), the modified Bessel function of the 1st kind, order zero}

function bessel_i0e(x: double): double;
  {-Return I0(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order zero}

function bessel_i1(x: double): double;
  {-Return I1(x), the modified Bessel function of the 1st kind, order one}

function bessel_i1e(x: double): double;
  {-Return I1(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order one}

function bessel_in(n: integer; x: double): double;
  {-Return I_n(x), the modified Bessel function of the 1st kind, order n.}

function bessel_iv(v, x: double): double;
  {-Return I_v(x), the modified Bessel function of the 1st kind, order v.}

function bessel_ive(v, x: double): double;
  {-Return I_v(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order v.}

function bessel_j0(x: double): double;
  {-Return J0(x), the Bessel function of the 1st kind, order zero}

function bessel_j1(x: double): double;
  {-Return J1(x), the Bessel function of the 1st kind, order one}

function bessel_jn(n: integer; x: double): double;
  {-Return J_n(x), the Bessel function of the 1st kind, order n; not suitable for large n or x.}

function bessel_jv(v, x: double): double;
  {-Return J_v(x), the Bessel function of the 1st kind, order v; not suitable for large v.}

function bessel_k0(x: double): double;
  {-Return K0(x), the modified Bessel function of the 2nd kind, order zero, x>0}

function bessel_k0e(x: double): double;
  {-Return K0(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order zero, x>0}

function bessel_k1(x: double): double;
  {-Return K1(x), the modified Bessel function of the 2nd kind, order one, x>0}

function bessel_k1e(x: double): double;
  {-Return K1(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order one, x>0}

function bessel_kn(n: integer; x: double): double;
  {-Return K_n(x), the modified Bessel function of the 2nd kind, order n, x>0, not suitable for large n}

function bessel_kv(v, x: double): double;
  {-Return K_v(x), the modified Bessel function of the 2nd kind, order v, x>0}

function bessel_kve(v, x: double): double;
  {-Return K_v(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order v, x>0}

function bessel_y0(x: double): double;
  {-Return Y0(x), the Bessel function of the 2nd kind, order zero; x>0}

function bessel_y1(x: double): double;
  {-Return Y1(x), the Bessel function of the 2nd kind, order one; x>0}

function bessel_yn(n: integer; x: double): double;
  {-Return Y_n(x), the Bessel function of the 2nd kind, order n, x>0, not suitable for large n or x}

function bessel_yv(v, x: double): double;
  {-Return Y_v(x), the Bessel function of the 2nd kind, order v; x > 0; not suitable for large v.}

function kelvin_ber(x: double): double;
  {-Return the Kelvin function ber(x)}

function kelvin_bei(x: double): double;
  {-Return the Kelvin function bei(x)}

function kelvin_ker(x: double): double;
  {-Return the Kelvin function ker(x), x > 0}

function kelvin_kei(x: double): double;
  {-Return the Kelvin function kei(x), x >= 0}

procedure kelvin_kerkei(x: double; var kr, ki: double);
  {-Return the Kelvin functions kr=ker(x), ki=kei(x), x > 0}

procedure kelvin_berbei(x: double; var br, bi: double);
  {-Return the Kelvin functions br=ber(x), bi=bei(x)}

function sph_bessel_jn(n: integer; x: double): double;
  {-Return j_n(x), the spherical Bessel function of the 1st kind, order n}

function sph_bessel_yn(n: integer; x: double): double;
  {-Return y_n(x), the spherical Bessel function of the 2nd kind, order n >=0 , x<>0}

function sph_bessel_in(n: integer; x: double): double;
  {-Return i_n(x), the modified spherical Bessel function of the 1st/2nd kind, order n}

function sph_bessel_ine(n: integer; x: double): double;
  {-Return i_n(x)*exp(-|x|), the exponentially scaled modified spherical Bessel function of the 1st/2nd kind, order n}

function sph_bessel_kn(n: integer; x: double): double;
  {-Return k_n(x), the modified spherical Bessel function of the 3rd kind, order n, x>0}

function sph_bessel_kne(n: integer; x: double): double;
  {-Return k_n(x)*exp(x), the exponentially scaled modified spherical Bessel function of the 3rd kind, order n, x>0}

function struve_h0(x: double): double;
  {-Return H0(x), the Struve function of order 0}

function struve_h1(x: double): double;
  {-Return H1(x), the Struve function of order 1}

function struve_h(v, x: double): double;
  {-Return H_v(x), the Struve function of order v, x < 0 only if v is an integer.}

function struve_l0(x: double): double;
  {-Return L0(x), the modified Struve function of order 0}

function struve_l1(x: double): double;
  {-Return L1(x), the modified Struve function of order 1}

function struve_l(v, x: double): double;
  {-Return L_v(x), the modified Struve function of order v, x < 0 only if v is an integer.}

{#Z+}
{---------------------------------------------------------------------------}
{---------------------- Gamma function and related -------------------------}
{---------------------------------------------------------------------------}
{#Z-}
const
  MAXGAM : double = 171.62437695630272;        {max. argument for gamma}
  MAXLGM : double = 2.559983327851638e305;     {max. argument for lngamma}

const
  MAXDFAC= 300;                                {max. argument for dfac}


function gamma(x: double): double;
  {-Return gamma(x), x <= MAXGAM; invalid if x is a non-positive integer}

function gamma1pm1(x: double): double;
  {-Return gamma(1+x)-1 with increased accuracy for x near 0}

function gammastar(x: double): double;
  {-Return Temme's gammastar(x) = gamma(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x)), x>0.}
  { For large x the asymptotic expansion is gammastar(x) = 1 + 1/12x + O(1/x^2)}

function gamma_delta_ratio(x,d: double): double;
  {-Return gamma(x)/gamma(x+d), accurate even for |d| << |x|}

function gamma_ratio(x,y: double): double;
  {-Return gamma(x)/gamma(y)}

function pochhammer(a,x: double): double;
  {-Return the Pochhammer symbol gamma(a+x)/gamma(a)}

function poch1(a,x: double): double;
  {-Return (pochhammer(a,x)-1)/x, psi(a) if x=0; accurate even for small |x|}

function binomial(n,k: integer): double;
  {-Return the binomial coefficient 'n choose k'}

function fac(n: integer): double;
  {-Return the factorial n!, n<MAXGAM-1; INF if n<0}

function dfac(n: integer): double;
  {-Return the double factorial n!!, n<=MAXDFAC; INF for even n<0}

function lnfac(n: longint): double;
  {-Return ln(n!), INF if n<0}

function lngamma(x: double): double;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, invalid if x is a non-positive integer}
  { Function signgamma can be used if the sign of gamma(x) is needed.}

function lngammas(x: double; var s: integer): double;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, s=-1,1 is the sign of gamma}

function lngamma1p(x: double): double;
  {-Return ln(|gamma(1+x)|) with increased accuracy for x near 0}

function lngamma_inv(y: double): double;
  {-Inverse of lngamma: return x with lngamma(x) = y, y >= -0.12142, x > 1.4616}

function signgamma(x: double): double;
  {-Return sign(gamma(x)), useless for 0 or negative integer}

function rgamma(x: double): double;
  {-Return the reciprocal gamma function rgamma = 1/gamma(x)}

procedure incgamma(a,x: double; var p,q: double);
  {-Return the normalised incomplete gamma functions P and Q, a>=0, x>=0}
  { P(a,x) = integral(exp(-t)*t^(a-1), t=0..x  )/gamma(a)}
  { Q(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf)/gamma(a)}

procedure incgamma_inv(a,p,q: double; var x: double; var ierr: integer);
  {-Return the inverse normalised incomplete gamma function, i.e. calculate}
  { x with P(a,x)=p and Q(a,x)=q. Input parameter a>0, p>=0, q>0 and p+q=1.}
  { ierr is >= 0 for success, < 0 for input errors or iterations failures. }

function igamma_inv(a,p,q: double): double;
  {-Return the inverse normalised incomplete gamma function, i.e. calculate}
  { x with P(a,x)=p and Q(a,x)=q. Input parameter a>0, p>=0, q>0 and p+q=1.}

function igammap_inv(a,p: double): double;
  {-Inverse incomplete gamma: return x with P(a,x)=p, a>=0, 0<=p<1}

function igammaq_inv(a,q: double): double;
  {-Inverse complemented incomplete gamma: return x with Q(a,x)=q, a>=0, 0<q<=1}

function igammap(a,x: double): double;
  {-Return the normalised lower incomplete gamma function P(a,x), a>=0, x>=0}
  { P(a,x) = integral(exp(-t)*t^(a-1), t=0..x)/gamma(a)}

function igammaq(a,x: double): double;
  {-Return the normalised upper incomplete gamma function Q(a,x), a>=0, x>=0}
  { Q(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf)/gamma(a)}

function igamma(a,x: double): double;
  {-Return the non-normalised upper incomplete gamma function}
  { GAMMA(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf), x>=0}

function igammal(a,x: double): double;
  {-Return the non-normalised lower incomplete gamma function}
  { gamma(a,x) = integral(exp(-t)*t^(a-1), t=0..x); x>=0, a<>0,-1,-2,..}

function igammat(a,x: double): double;
  {-Return Tricomi's entire incomplete gamma function igammat(a,x)}
  { = igammal(a,x)/gamma(a)/x^a = P(a,x)/x^a }

function psi(x: double): double;
  {-Return the psi (digamma) function of x, INF if x is a non-positive integer}

function psi_inv(y: double): double;
  {-Inverse of psi, return x with psi(x)=y, y <= ln_MaxDbl}

function trigamma(x: double): double;
  {-Return the trigamma function of x, INF if x is a negative integer}

function tetragamma(x: double): double;
  {-Return the tetragamma function psi''(x), NAN/RTE if x is a negative integer}

function pentagamma(x: double): double;
  {-Return the pentagamma function psi'''(x), INF if x is a negative integer}

function polygamma(n: integer; x: double): double;
  {-Return the polygamma function: n'th derivative of psi; n>=0, x>0 for n>MAXGAMD.}
  { Note: The accuracy may be reduced for n>=MAXGAMD due to ln/exp operations.}

function BatemanG(x: double): double;
  {-Return the Bateman function G(x) = psi((x+1)/2) - psi(x/2), x<>0,-1,-2,...}

function lnbeta(x,y: double): double;
  {-Return the logarithm of |beta(x,y)|=|gamma(x)*gamma(y)/gamma(x+y)|}

function beta(x,y: double): double;
  {-Return the function beta(x,y)=gamma(x)*gamma(y)/gamma(x+y)}

function beta3(a, b, x: double): double;
  {-Return the non-normalised incomplete beta function B_x(a,b) for}
  { a>0, b>0, 0<=x<=1, B_x = integral(t^(a-1)*(1-t)^(b-1), t=0..x).}

function ibeta(a, b, x: double): double;
  {-Return the normalised incomplete beta function, a>0, b>0, 0 <= x <= 1}
  { ibeta = integral(t^(a-1)*(1-t)^(b-1) / betax(a,b), t=0..x)}

function ibeta_inv(a, b, y: double): double;
  {-Return the functional inverse of the normalised incomplete beta function}
  { with a > 0, b > 0, and 0 <= y <= 1.}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Zeta functions and polylogarithms   -------------------}
{---------------------------------------------------------------------------}
{#Z-}

function zeta(s: double): double;
  {-Return the Riemann zeta function at s, s<>1}

function zeta1p(x: double): double;
  {-Return the Riemann zeta function at 1+x, x<>0}

function zetaint(n: integer): double;
  {-Return zeta(n) for integer arguments, n<>1}

function zetam1(s: double): double;
  {-Return Riemann zeta(s)-1, s<>1}

function zetah(s,a: double): double;
  {-Return the Hurwitz zeta function zetah(s,a) = sum(1/(i+a)^s, i=0..INF), s<>1, a>0}

function LerchPhi(z,s,a: double): double;
  {-Return the Lerch transcendent Phi(z,s,a) = sum(z^n/(n+a)^s, n=0..INF), |z|<=1, s >= -1, a>0; s<>1 if z=1}

function DirichletBeta(s: double): double;
  {-Return the Dirichlet beta function sum((-1)^n/(2n+1)^s, n=0..INF)}

function DirichletLambda(s: double): double;
  {-Return the Dirichlet lambda function sum(1/(2n+1)^s, n=0..INF), s<>1}

function LegendreChi(s, x: double): double;
  {-Return Legendre's Chi-function chi(s,x); s>=0, |x|<=1, x<>1 if s<=1}

function eta(s: double): double;
  {-Return the Dirichlet eta function}

function etaint(n: integer): double;
  {-Return the Dirichlet function eta(n) for integer arguments}

function etam1(s: double): double;
  {-Return Dirichlet eta(s)-1}

function primezeta(x: double): double;
  {-Return the prime zeta function P(x) = sum(1/p^x, p prime), x > 1}

function harmonic(x: double): double;
  {-Return the harmonic function H(x) = psi(x+1) + EulerGamma}

function harmonic2(x,r: double): double;
  {-Return the generalized harmonic function H(x,r) = zeta(r)-zetah(r,x+1); x >= -1}

function cl2(x: double): double;
  {-Return the Clausen function: integral(-ln(2*|sin(t/2)|),t=0..x) = Im(Li_2(exp(ix)))}

function ti2(x: double): double;
  {-Return the inverse tangent integral, ti2(x) = integral(arctan(t)/t, t=0..x)}

function lobachevsky_c(x: double): double;
  {-Return the Lobachevski function L(x) = integral(-ln(|cos(t)|), t=0..x)}

function lobachevsky_s(x: double): double;
  {-Return the Lobachevski function Lambda(x) = integral(-ln(|2sin(t)|), t=0..x)}

function dilog(x: double): double;
  {-Return dilog(x) = Re(Li_2(x)), Li_2(x) = -integral(ln(1-t)/t, t=0..x)}

function trilog(x: double): double;
  {-Return the trilogarithm function trilog(x) = Re(Li_3(x))}

function polylog(n: integer; x: double): double;
  {-Return the polylogarithm Li_n(x) of integer order; x<1 for n >= 0}

function polylogr(s, x: double): double;
  {-Return the polylogarithm Li_s(x) of real order; s >= -1, |x|<=1, x<>1 if s=1}

function fermi_dirac(n: integer; x: double): double;
  {-Return the integer order Fermi-Dirac integral F_n(x) = 1/n!*integral(t^n/(exp(t-x)+1), t=0..INF)}

function fermi_dirac_m05(x: double): double;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}

function fermi_dirac_p05(x: double): double;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}

function fermi_dirac_p15(x: double): double;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}

function fermi_dirac_p25(x: double): double;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Legendre style) -------------------}
{---------------------------------------------------------------------------}
{#Z-}
function comp_ellint_1(k: double): double;
  {-Return the complete elliptic integral of the 1st kind, |k| < 1}

function comp_ellint_2(k: double): double;
  {-Return the complete elliptic integral of the 2nd kind, |k| <= 1}

function comp_ellint_3(nu,k: double): double;
  {-Return the complete elliptic integral of the 3rd kind, |k|<1, nu<>1}

function comp_ellint_d(k: double): double;
  {-Return the complete elliptic integral D(k) = (K(k) - E(k))/k^2, |k| < 1}

function ellint_1(phi,k: double): double;
  {-Return the Legendre elliptic integral F(phi,k) of the 1st kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}

function ellint_2(phi,k: double): double;
  {-Return the Legendre elliptic integral E(phi,k) of the 2nd kind}
  { = integral(sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}

function ellint_3(phi,nu,k: double): double;
  {-Return the Legendre elliptic integral PI(phi,nu,k) of the 3rd kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2)/(1-nu*sin(x)^2),x=0..phi) with  }
  { |k*sin(phi)|<=1, returns Cauchy principal value if nu*sin(phi)^2>1}

function ellint_d(phi,k: double): double;
  {-Return the Legendre elliptic integral D(phi,k) = (F(phi,k) - E(phi,k))/k^2 }
  { = integral(sin(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }

function heuman_lambda(phi,k: double): double;
  {-Return Heuman's function Lambda_0(phi,k) = F(phi,k')/K(k') + 2/Pi*K(k)*Z(phi,k'), |k|<=1}

function jacobi_zeta(phi,k: double): double;
  {-Return the Jacobi Zeta function Z(phi,k) = E(phi,k) - E(k)/K(k)*F(phi,k), |k|<=1}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Carlson style) --------------------}
{---------------------------------------------------------------------------}
{#Z-}
function ell_rc(x,y: double): double;
  {-Return Carlson's degenerate elliptic integral RC; x>=0, y<>0}

function ell_rf(x,y,z: double): double;
  {-Return Carlson's elliptic integral of the 1st kind; x,y,z >=0, at most one =0}

function ell_rd(x,y,z: double): double;
  {-Return Carlson's elliptic integral of the 2nd kind; z>0; x,y >=0, at most one =0}

function ell_rg(x,y,z: double): double;
  {-Return Carlson's completely symmetric elliptic integral of the 2nd kind; x,y,z >= 0}

function ell_rj(x,y,z,r: double): double;
  {-Return Carlson's elliptic integral of the 3rd kind; r<>0; x,y,z >=0, at most one =0}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Bulirsch style) -------------------}
{---------------------------------------------------------------------------}
{#Z-}
function cel1(kc: double): double;
  {-Return Bulirsch's complete elliptic integral of the 1st kind, kc<>0}

function cel2(kc, a, b: double): double;
  {-Return Bulirsch's complete elliptic integral of the 2nd kind, kc<>0}

function cel(kc, p, a, b: double): double;
  {-Return Bulirsch's general complete elliptic integral, kc<>0, Cauchy principle value if p<0}

function el1(x,kc: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 1st kind}

function el2(x,kc,a,b: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 2nd kind, kc<>0}

function el3(x,kc,p: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 3rd kind, 1+p*x^2<>0}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Maple V style) --------------------}
{---------------------------------------------------------------------------}
{#Z-}
function EllipticF(z,k: double): double;
  {-Return the incomplete elliptic integral of the 1st kind; |z|<=1, |k*z|<1}

function EllipticK(k: double): double;
  {-Return the complete elliptic integral of the 1st kind, |k| < 1}

function EllipticKim(k: double): double;
  {-Return K(i*k), the complete elliptic integral of the 1st kind with}
  { imaginary modulus = integral(1/sqrt(1-x^2)/sqrt(1+k^2*x^2),x=0..1)}

function EllipticCK(k: double): double;
  {-Return the complementary complete elliptic integral of the 1st kind, k<>0}

function EllipticE(z,k: double): double;
  {Return the incomplete elliptic integrals of the 2nd kind, |z|<=1, |k*z| <= 1}

function EllipticEC(k: double): double;
  {-Return the complete elliptic integral of the 2nd kind, |k| < 1}

function EllipticECim(k: double): double;
  {-Return E(i*k), the complete elliptic integral of the 2nd kind with}
  { imaginary modulus = integral(sqrt(1+k^2*x^2)/sqrt(1-x^2),x=0..1)  }

function EllipticCE(k: double): double;
  {-Return the complementary complete elliptic integral of the 2nd kind}

function EllipticPi(z,nu,k: double): double;
  {-Return the incomplete elliptic integral of the 3rd kind, |z|<=1, |k*z|<1}

function EllipticPiC(nu,k: double): double;
  {-Return the complete elliptic integral of the 3rd kind, |k|<1, nu<>1}

function EllipticCPi(nu,k: double): double;
  {-Return the complementary complete elliptic integral of the 3rd kind, k<>0, nu<>1}

{#Z+}
{-----------------------------------------------------------------------------}
{--------------------- Jacobi elliptic and theta functions -------------------}
{-----------------------------------------------------------------------------}
{#Z-}
function EllipticModulus(q: double): double;
  {-Return the elliptic modulus k(q) = theta_2(q)^2/theta_3(q)^2, 0 <= q <= 1}

function EllipticNome(k: double): double;
  {-Return the elliptic nome q(k) = exp(-Pi*EllipticCK(k)/EllipticK(k)), |k| < 1}

function jacobi_am(x,k: double): double;
  {-Return the Jacobi amplitude am(x,k)}

function jacobi_arccn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccn(x,k), |x| <= 1, x >= sqrt(1 - 1/k^2) if k >= 1}

function jacobi_arccd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccd(x,k); |x| <= 1 if |k| < 1; |x| >= 1 if |k| > 1 }

function jacobi_arccs(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccs(x,k), |x| >= sqrt(k^2-1) for |k|>1}

function jacobi_arcdc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcdc(x,k); |x| >= 1 if |k| < 1; |x| <= 1 if |k| > 1 }

function jacobi_arcdn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcdn(x,k), 0 <= x <= 1, x^2 + k^2 > 1 if |k| < 1;  |x| <= 1 if |k| > 1}

function jacobi_arcds(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcds(x,k), x^2 + k^2 >= 1}

function jacobi_arcnc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcnc(x,k), x >= 1, x^2 <= k^2/(k^2-1) for |k|>1}

function jacobi_arcnd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcnd(x,k), x >= 1, x^2 <= k^2/(1-k^2) if k < 1}

function jacobi_arcns(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcns(x,k), |x| >= 1, |x| >= k if k>=1}

function jacobi_arcsc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsc(x,k), |x| <= 1/sqrt(k^2-1) for |k|>1}

function jacobi_arcsd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsd(x,k), x^2*(1-k^2) <= 1}

function jacobi_arcsn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsn(x,k), |x| <= 1 and |x*k| <= 1}

function jacobi_sn(x,k: double): double;
  {-Return the Jacobi elliptic function sn(x,k)}

function jacobi_cn(x,k: double): double;
  {-Return the Jacobi elliptic function cn(x,k)}

function jacobi_dn(x,k: double): double;
  {-Return the Jacobi elliptic function dn(x,k)}

function jacobi_nc(x,k: double): double;
  {-Return the Jacobi elliptic function nc(x,k)}

function jacobi_sc(x,k: double): double;
  {-Return the Jacobi elliptic function sc(x,k)}

function jacobi_dc(x,k: double): double;
  {-Return the Jacobi elliptic function dc(x,k)}

function jacobi_nd(x,k: double): double;
  {-Return the Jacobi elliptic function nd(x,k)}

function jacobi_sd(x,k: double): double;
  {-Return the Jacobi elliptic function sd(x,k)}

function jacobi_cd(x,k: double): double;
  {-Return the Jacobi elliptic function cd(x,k)}

function jacobi_ns(x,k: double): double;
  {-Return the Jacobi elliptic function ns(x,k)}

function jacobi_cs(x,k: double): double;
  {-Return the Jacobi elliptic function cs(x,k)}

function jacobi_ds(x,k: double): double;
  {-Return the Jacobi elliptic function ds(x,k)}

function jacobi_theta(n: integer; x,q: double): double;
  {-Return the Jacobi theta function theta_n(x,q), n=1..4, 0 <= q < 1}

procedure sncndn(x,mc: double; var sn,cn,dn: double);
  {-Return the Jacobi elliptic functions sn,cn,dn for argument x and}
  { complementary parameter mc.}

function theta1p(q: double): double;
  {-Return the derivative  theta1p(q) := d/dx(theta_1(x,q)) at x=0,}
  { = 2*q^(1/4)*sum((-1)^n*(2n+1)*q^(n*(n+1)),n=0..Inf), 0 <= q < 1}

function theta2(q: double): double;
  {-Return Jacobi theta_2(q) = 2*q^(1/4)*sum(q^(n*(n+1)),n=0..Inf) 0 <= q < 1}

function theta3(q: double): double;
  {-Return Jacobi theta_3(q) = 1 + 2*sum(q^(n*n)),n=1..Inf); |q| < 1}

function theta4(q: double): double;
  {-Return Jacobi theta_4(q) = 1 + 2*sum((-1)^n*q^(n*(n+1)),n=1..Inf); |q| < 1}

procedure sincos_lemn(x: double; var sl,cl: double);
  {-Return the lemniscate functions sl = sin_lemn(x), cl = cos_lemn(x)}

function sin_lemn(x: double): double;
  {-Return the lemniscate sine functions sl = sin_lemn(x)}

function cos_lemn(x: double): double;
  {-Return the lemniscate cosine functions cl = cos_lemn(x)}

{#Z+}
{---------------------------------------------------------------------------}
{----------------------- Error function and related ------------------------}
{---------------------------------------------------------------------------}
{#Z-}

function dawson(x: double): double;
  {-Return Dawson's integral: dawson(x) = exp(-x^2)*integral(exp(t^2), t=0..x)}

function dawson2(p,x: double): double;
  {-Return the generalized Dawson integral F(p,x) = exp(-x^p)*integral(exp(t^p), t=0..x); x,p >= 0}

function erf(x: double): double;
  {-Return the error function erf(x) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..x)}

function erfc(x: double): double;
  {-Return the complementary error function erfc(x) = 1-erf(x)}

function erfce(x: double): double;
  {-Return the exponentially scaled complementary error function erfce(x) = exp(x^2)*erfc(x)}

function inerfc(n: integer; x: double): double;
  {-Return the repeated integrals of erfc, n >= -1; scaled with exp(x^2) for x>0}

function erfg(p,x: double): double;
  {-Return the generalized error function integral(exp(-t^p), t=0..x); x,p >= 0}

function erfi(x: double): double;
  {-Return the imaginary error function erfi(x) = erf(ix)/i}

function erf_inv(x: double): double;
  {-Return the inverse function of erf, erf(erf_inv(x)) = x, -1 < x < 1}

function erfc_inv(x: double): double;
  {-Return the inverse function of erfc, erfc(erfc_inv(x)) = x, 0 < x < 2}

function erf_p(x: double): double;
  {-Return the probability function erf_p = integral(exp(-t^2/2)/sqrt(2*Pi), t=-Inf..x)}

function erf_q(x: double): double;
  {-Return the probability function erf_q = integral(exp(-t^2/2)/sqrt(2*Pi), t=x..Inf)}

function erf_z(x: double): double;
  {-Return the probability function erf_z = exp(-x^2/2)/sqrt(2*Pi)}

function expint3(x: double): double;
  {-Return the integral(exp(-t^3), t=0..x), x >= 0}

procedure Fresnel(x: double; var s,c: double);
  {-Return the Fresnel integrals S(x)=integral(sin(Pi/2*t^2),t=0..x) and C(x)=integral(cos(Pi/2*t^2),t=0..x)}

function FresnelC(x: double): double;
  {-Return the Fresnel integral C(x)=integral(cos(Pi/2*t^2),t=0..x)}

function FresnelS(x: double): double;
  {-Return the Fresnel integral S(x)=integral(sin(Pi/2*t^2),t=0..x)}

function gsi(x: double): double;
  {-Return the Goodwin-Staton integral = integral(exp(-t*t)/(t+x), t=0..Inf), x > 0}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Exponential integrals and related ---------------------}
{---------------------------------------------------------------------------}
{#Z-}

function chi(x: double): double;
  {-Return the hyperbolic cosine integral = EulerGamma + ln(|x|) + integral((cosh(t)-1)/t, t=0..|x|)}

function ci(x: double): double;
  {-Return the cosine integral, ci(x) = EulerGamma + ln(|x|) + integral((cos(t)-1)/t, t=0..|x|)}

function e1(x: double): double;
  {-Return the exponential integral E1(x) = integral(exp(-x*t)/t, t=1..Inf), x <> 0}

function cin(x: double): double;
  {-Return the entire cosine integral, cin(x) = integral((1-cos(t))/t, t=0..x)}

function cinh(x: double): double;
  {-Return the entire hyperbolic cosine integral, cinh(x) = integral((cosh(t)-1)/t, t=0..x)}

function ei(x: double): double;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}

function ein(x: double): double;
  {-Return the entire exponential integral ein(x) = integral((1-exp(-t))/t, t=0..x)}

function ei_inv(x: double): double;
  {-Return the functional inverse of Ei(x), ei_inv(ei(x))=x}

function en(n: longint; x: double): double;
  {-Return the exponential integral E_n(x) = integral(exp(-x*t)/t^n, t=1..Inf), x > 0}

function gei(p,x: double): double;
  {-Return the generalized exponential integral E_p(x) = integral(exp(-x*t)/t^p, t=1..Inf), x >= 0}

function li(x: double): double;
  {-Return the logarithmic integral li(x) = PV-integral(1/ln(t), t=0..x), x >= 0, x <> 1}

function shi(x: double): double;
  {-Return the hyperbolic sine integral, integral(sinh(t)/t, t=0..x)}

function si(x: double): double;
  {-Return the sine integral, si(x) = integral(sin(t)/t, t=0..x)}

function ssi(x: double): double;
  {-Return the shifted sine integral, ssi(x) = si(x) - Pi/2}

{#Z+}
{---------------------------------------------------------------------------}
{------------------ Orthogonal polynomials and related ---------------------}
{---------------------------------------------------------------------------}
{#Z-}

function chebyshev_t(n: integer; x: double): double;
  {-Return Tn(x), the Chebyshev polynomial of the first kind, degree n}

function chebyshev_u(n: integer; x: double): double;
  {-Return Un(x), the Chebyshev polynomial of the second kind, degree n}

function chebyshev_v(n: integer; x: double): double;
  {-Return V_n(x), the Chebyshev polynomial of the third kind, degree n>=0}

function chebyshev_w(n: integer; x: double): double;
  {-Return W_n(x), the Chebyshev polynomial of the fourth kind, degree n>=0}

function gegenbauer_c(n: integer; a,x: double): double;
  {-Return Cn(a,x), the nth Gegenbauer (ultraspherical) polynomial with}
  { parameter a. The degree n must be non-negative; a should be > -0.5 }
  { When a = 0,   C0(0,x) = 1,  and   Cn(0,x) = 2/n*Tn(x)   for n <> 0.}

function hermite_h(n: integer; x: double): double;
  {-Return Hn(x), the nth Hermite polynomial, degree n >= 0}

function jacobi_p(n: integer; a,b,x: double): double;
  {-Return Pn(a,b,x), the nth Jacobi polynomial with parameters a,b. Degree n}
  { must be >= 0; a,b should be > -1 (a+b must not be an integer < -1).}

function laguerre(n: integer; a,x: double): double;
  {-Return Ln(a,x), the nth generalized Laguerre polynomial with parameter a;}
  { degree n must be >= 0. x >=0 and a > -1 are the standard ranges.}

function laguerre_ass(n,m: integer; x: double): double;
  {-Return the associated Laguerre polynomial Ln(m,x); n,m >= 0}

function laguerre_l(n: integer; x: double): double;
  {-Return the nth Laguerre polynomial Ln(0,x); n >= 0}

function legendre_p(l: integer; x: double): double;
  {-Return P_l(x), the Legendre polynomial/function P_l, degree l}

function legendre_q(l: integer; x: double): double;
  {-Return Q_l(x), the Legendre function of the 2nd kind, degree l >=0, |x| <> 1}

function legendre_plm(l,m: integer; x: double): double;
  {-Return the associated Legendre polynomial P_lm(x)}

function legendre_qlm(l,m: integer; x: double): double;
  {-Return Q(l,m,x), the associated Legendre function of the second kind; l >= 0, l+m >= 0, |x|<>1}

procedure spherical_harmonic(l, m: integer; theta, phi: double; var yr,yi: double);
  {-Return Re and Im of the spherical harmonic function Y_lm(theta,phi)}

function toroidal_plm(l,m: integer; x: double): double;
  {-Return the toroidal harmonic function P(l-0.5,m,x); l,m=0,1; x >= 1}

function toroidal_qlm(l,m: integer; x: double): double;
  {-Return the toroidal harmonic function Q(l-0.5,m,x); l=0,1; x > 1}

function zernike_r(n,m: integer; r: double): double;
  {-Return the Zernike radial polynomial Rnm(r), r >= 0, n >= m >= 0, n-m even}

{#Z+}
{---------------------------------------------------------------------------}
{---------------------- Statistical distributions --------------------------}
{---------------------------------------------------------------------------}
{#Z-}

function beta_pdf(a, b, x: double): double;
  {-Return the probability density function of the beta distribution with}
  { parameters a and b: beta_pdf = x^(a-1)*(1-x)^(b-1) / beta(a,b)}

function beta_cdf(a, b, x: double): double;
  {-Return the cumulative beta distribution function, a>0, b>0}

function beta_inv(a, b, y: double): double;
  {-Return the functional inverse of the beta distribution function. a>0, b>0;}
  { 0 <= y <= 1. Given y the function finds x such that beta_cdf(a, b, x) = y}

function binomial_cdf(p: double; n, k: longint): double;
  {-Return the cumulative binomial distribution function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}

function binomial_pmf(p: double; n, k: longint): double;
  {-Return the binomial distribution probability mass function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}

function cauchy_pdf(a, b, x: double): double;
  {-Return the Cauchy probability density function with location a }
  { and scale b > 0, 1/(Pi*b*(1+((x-a)/b)^2))}

function cauchy_cdf(a, b, x: double): double;
  {-Return the cumulative Cauchy distribution function with location a}
  { and scale b > 0, = 1/2 + arctan((x-a)/b)/Pi}

function cauchy_inv(a, b, y: double): double;
  {-Return the functional inverse of Cauchy distribution function}
  { with location a and scale b > 0}

function chi2_pdf(nu: longint; x: double): double;
  {-Return the probability density function of the chi-square distribution, nu>0}

function chi2_cdf(nu: longint; x: double): double;
  {-Return the cumulative chi-square distribution with nu>0 degrees of freedom, x >= 0}

function chi2_inv(nu: longint; p: double): double;
  {-Return the functional inverse of the chi-square distribution, nu>0, 0 <= p < 1}

function evt1_pdf(a, b, x: double): double;
  {-Return the probability density function of the Extreme Value Type I distribution}
  { with location a and scale b > 0, result = exp(-(x-a)/b)/b * exp(-exp(-(x-a)/b)) }

function evt1_cdf(a, b, x: double): double;
  {-Return the cumulative Extreme Value Type I distribution function}
  { with location a and scale b > 0;  result = exp(-exp(-(x-a)/b)). }

function evt1_inv(a, b, y: double): double;
  {-Return the functional inverse of the Extreme Value Type I distribution}
  { function with location a and scale b > 0;  result = a - b*ln(ln(-y)). }

function exp_pdf(a, alpha, x: double): double;
  {-Return the exponential probability density function with location a }
  { and rate alpha > 0, = alpha*exp(-alpha*(x-a)) if x >= a, 0 if x < a.}

function exp_cdf(a, alpha, x: double): double;
  {-Return the cumulative exponential distribution function with location a}
  { and rate alpha > 0, = 1 - exp(-alpha*(x-a)) if x >= a, 0 if x < a.}

function exp_inv(a, alpha, y: double): double;
  {-Return the functional inverse of exponential distribution function with}
  { location a and rate alpha > 0}

function f_pdf(nu1, nu2: longint; x: double): double;
  {-Return the probability density function of the F distribution; x >= 0, nu1, nu2 > 0}

function f_cdf(nu1, nu2: longint; x: double): double;
  {-Return the cumulative F distribution function; x >= 0, nu1, nu2 > 0}

function f_inv(nu1, nu2: longint; y: double): double;
  {-Return the functional inverse of the F distribution, nu1, nu2 > 0, 0 <= y <= 1}

function gamma_pdf(a, b, x: double): double;
  {-Return the probability density function of a gamma distribution with shape}
  { a>0, scale b>0: gamma_pdf = x^(a-1)*exp(-x/b)/gamma(a)/b^a, x>0}

function gamma_cdf(a, b, x: double): double;
  {-Return the cumulative gamma distribution function, shape a>0, scale b>0}

function gamma_inv(a, b, p: double): double;
  {-Return the functional inverse of the gamma distribution function, shape a>0,}
  { scale b>0, 0 <= p <= 1, i.e. finds x such that gamma_cdf(a, b, x) = p}

function invgamma_pdf(a, b, x: double): double;
  {-Return the probability density function of an inverse gamma distribution}
  { with shape a>0, scale b>0: result = (b/x)^a/x*exp(-b/x)/Gamma(a), x >= 0}

function invgamma_cdf(a, b, x: double): double;
  {-Return the cumulative inverse gamma distribution function, shape a>0, scale}
  { b>0: result = Gamma(a,b/x)/Gamma(a) = Q(a,b/x) = igammaq(a,b/x), x >= 0}

function invgamma_inv(a, b, y: double): double;
  {-Return the functional inverse of the inverse gamma distribution function, shape}
  { a>0, scale b>0, 0 <= y <= 1, i.e. find x such that invgamma_cdf(a, b, x) = y  }

function hypergeo_pmf(n1,n2,n,k: longint): double;
  {-Return the hypergeometric distribution probability mass function; n,n1,n2 >= 0, n <= n1+n2;}
  { i.e. the probability that among n randomly chosen samples from a container}
  { with n1 type1 objects and n2 type2 objects are exactly k type1 objects:}

function hypergeo_cdf(n1,n2,n,k: longint): double;
  {-Return the cumulative hypergeometric distribution function; n,n1,n2 >= 0, n <= n1+n2}

function kumaraswamy_pdf(a, b, x: double): double;
  {-Return the Kumaraswamy probability density function with shape}
  { parameters a,b>0, 0<=x<=1; result = a*b*x^(a-1)*(1-x^a)^(b-1) }

function kumaraswamy_cdf(a, b, x: double): double;
  {-Return the cumulative Kumaraswamy distribution function with}
  { shape parameters a,b > 0, 0 <= x <= 1; result = 1-(1-x^a)^b}

function kumaraswamy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Kumaraswamy distribution}
  { with shape parameters a,b > 0; result = [1-(1-y)^(1/b)]^(1/a)}

function laplace_cdf(a, b, x: double): double;
  {-Return the cumulative Laplace distribution function with location a and scale b > 0}

function laplace_inv(a, b, y: double): double;
  {-Return the functional inverse of the Laplace distribution with location a and scale b > 0}

function laplace_pdf(a, b, x: double): double;
  {-Return the Laplace probability density function with location a}
  { and scale b > 0, result = exp(-abs(x-a)/b) / (2*b)}

function levy_pdf(a, b, x: double): double;
  {-Return the Levy probability density function with}
  { location a and scale parameter b > 0}

function levy_cdf(a, b, x: double): double;
  {-Return the cumulative Levy distribution function with}
  { location a and scale parameter b > 0}

function levy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Levy distribution}
  { with location a and scale parameter b > 0}

function logistic_pdf(a, b, x: double): double;
  {-Return the logistic probability density function with location a}
  { and scale parameter b > 0, exp(-(x-a)/b)/b/(1+exp(-(x-a)/b))^2}

function logistic_cdf(a, b, x: double): double;
  {-Return the cumulative logistic distribution function with}
  { location a and scale parameter b > 0}

function logistic_inv(a, b, y: double): double;
  {-Return the functional inverse of the logistic distribution}
  { with location a and scale parameter b > 0}

function lognormal_pdf(a, b, x: double): double;
  {-Return the log-normal probability density function with}
  { location a and scale parameter b > 0, zero for x <= 0.}

function lognormal_cdf(a, b, x: double): double;
  {-Return the cumulative log-normal distribution function with}
  { location a and scale parameter b > 0, zero for x <= 0.}

function lognormal_inv(a, b, y: double): double;
  {-Return the functional inverse of the log-normal distribution}
  { with location a and scale parameter b > 0, 0 < y < 1.}

function logseries_pmf(a: double; k: longint): double;
  {-Return the logarithmic (series) probability mass function}
  { with shape 0 < a < 1, k > 0; result = -a^k/(k*ln(1-a))   }

function logseries_cdf(a: double; k: longint): double;
  {-Return the cumulative logarithmic (series) distribution function with shape 0 < a < 1, k > 0}

function maxwell_pdf(b, x: double): double;
  {-Return the Maxwell probability density function with scale b > 0, x >= 0}

function maxwell_cdf(b, x: double): double;
  {-Return the cumulative Maxwell distribution function with scale b > 0, x >= 0}

function maxwell_inv(b, y: double): double;
  {-Return the functional inverse of the Maxwell distribution with scale b > 0}

function moyal_pdf(a, b, x: double): double;
  {-Return the Moyal probability density function with}
  { location a and scale parameter b > 0}

function moyal_cdf(a, b, x: double): double;
  {-Return the cumulative Moyal distribution function with}
  { location a and scale parameter b > 0}

function moyal_inv(a, b, y: double): double;
  {-Return the functional inverse of the Moyal distribution}
  { with location a and scale parameter b > 0}

function negbinom_cdf(p,r: double; k: longint): double;
  {-Return the cumulative negative binomial distribution function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}

function negbinom_pmf(p,r: double; k: longint): double;
  {-Return the negative binomial distribution probability mass function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}

function normal_pdf(mu, sd, x: double): double;
  {-Return the normal (Gaussian) probability density function with mean mu}
  { and standard deviation sd>0, exp(-0.5*(x-mu)^2/sd^2) / sqrt(2*Pi*sd^2)}

function normal_cdf(mu, sd, x: double): double;
  {-Return the normal (Gaussian) distribution function}
  { with mean mu and standard deviation sd > 0        }

function normal_inv(mu, sd, y: double): double;
  {-Return the functional inverse of the normal (Gaussian) distribution}
  { with mean mu and standard deviation sd > 0, 0 < y < 1.}

function normstd_pdf(x: double): double;
  {-Return the std. normal probability density function exp(-x^2/2)/sqrt(2*Pi)}

function normstd_cdf(x: double): double;
  {-Return the standard normal distribution function}

function normstd_inv(y: double): double;
  {-Return the inverse standard normal distribution function, 0 < y < 1.}
  { For x=normstd_inv(y) and y from (0,1), normstd_cdf(x) = y}

function pareto_pdf(k, a, x: double): double;
  {-Return the Pareto probability density function with minimum value k > 0}
  { and shape a, x >= a > 0, result = (a/x)*(k/x)^a}

function pareto_cdf(k, a, x: double): double;
  {-Return the cumulative Pareto distribution function minimum value k > 0}
  { and shape a, x >= a > 0, result = 1-(k/x)^a}

function pareto_inv(k, a, y: double): double;
  {-Return the functional inverse of the Pareto distribution with minimum}
  { value k > 0 and shape a, x >= a > 0, result = k/(1-x)^(1/a)}

function poisson_cdf(mu: double; k: longint): double;
  {-Return the cumulative Poisson distribution function with mean mu >= 0}

function poisson_pmf(mu: double; k: longint): double;
  {-Return the Poisson distribution probability mass function with mean mu >= 0}

function rayleigh_pdf(b, x: double): double;
  {-Return the Rayleigh probability density function with}
  { scale b > 0, x >= 0; result = x*exp(-0.5*(x/b)^2)/b^2}

function rayleigh_cdf(b, x: double): double;
  {-Return the cumulative Rayleigh distribution function with scale b > 0, x >= 0}

function rayleigh_inv(b, y: double): double;
  {-Return the functional inverse of the Rayleigh distribution with scale b > 0}

function triangular_pdf(a, b, c, x: double): double;
  {-Return the triangular probability density function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}

function triangular_cdf(a, b, c, x: double): double;
  {-Return the cumulative triangular distribution function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}

function triangular_inv(a, b, c, y: double): double;
  {-Return the functional inverse of the triangular distribution with}
  { lower limit a, upper limit b, mode c; a<b, a <= c <= b, 0 <= y <= 1}

function t_pdf(nu: longint; x: double): double;
  {-Return the probability density function of Student's t distribution, nu>0}

function t_cdf(nu: longint; t: double): double;
  {-Return the cumulative Student t distribution with nu>0 degrees of freedom}

function t_inv(nu: longint; p: double): double;
  {-Return the functional inverse of Student's t distribution, nu>0, 0 <= p <= 1}

function uniform_pdf(a, b, x: double): double;
  {-Return the uniform probability density function on [a,b], a<b}

function uniform_cdf(a, b, x: double): double;
  {-Return the cumulative uniform distribution function on [a,b], a<b}

function uniform_inv(a, b, y: double): double;
  {-Return the functional inverse of the uniform distribution on [a,b], a<b}

function wald_pdf(mu, b, x: double): double;
  {-Return the Wald (inverse Gaussian) probability density}
  { function with mean mu > 0, scale b > 0 for x >= 0     }

function wald_cdf(mu, b, x: double): double;
  {-Return the Wald (inverse Gaussian) distribution  }
  { function with mean mu > 0, scale b > 0 for x >= 0}

function wald_inv(mu, b, y: double): double;
  {-Return the functional inverse of the Wald (inverse Gaussian)}
  { distribution with mean mu > 0, scale b > 0, 0 <= y < 1.     }

function weibull_pdf(a, b, x: double): double;
  {-Return the Weibull probability density function with shape a > 0}
  { and scale b > 0, result = a*x^(a-1)*exp(-(x/b)^a)/ b^a, x > 0}

function weibull_cdf(a, b, x: double): double;
  {-Return the cumulative Weibull distribution function with}
  { shape parameter a > 0 and scale parameter b > 0}

function weibull_inv(a, b, y: double): double;
  {-Return the functional inverse of the Weibull distribution}
  { shape parameter a > 0 and scale parameter b > 0}

function zipf_pmf(r: double; k: longint): double;
  {-Return the Zipf distribution probability mass function k^(-(r+1))/zeta(r+1), r>0, k>0}

function zipf_cdf(r: double; k: longint): double;
  {-Return the cumulative Zipf distribution function H(k,r+1)/zeta(r+1), r>0, k>0}

{#Z+}
{---------------------------------------------------------------------------}
{-------------------- Hypergeometric functions -----------------------------}
{---------------------------------------------------------------------------}
{#Z-}

function hyperg_1F1(a,b,x: double): double;
  {-Return the confluent hypergeometric function 1F1(a,b,x); Kummer's function M(a,b,x)}

function hyperg_1F1r(a,b,x: double): double;
  {-Return the regularized Kummer hypergeometric function 1F1(a,b,x)/Gamma(b)}

function hyperg_u(a,b,x: double): double;
  {-Return Tricomi's confluent hypergeometric function U(a,b,x). If}
  { x<0, then a must be an integer and a<0 or 1+a-b an integer < 0.}

function hyperg_2F1(a,b,c,x: double): double;
  {-Return the Gauss hypergeometric function 2F1(a,b;c;x)}

function hyperg_2F1r(a,b,c,x: double): double;
  {-Return the regularized Gauss hypergeometric function 2F1(a,b,c,x)/Gamma(c)}

function WhittakerM(k,m,x: double): double;
  {-Return the Whittaker M function = exp(-x/2)*x^(0.5+m) * 1F1(m-k-0.5,2m+1,x)}

function WhittakerW(k,m,x: double): double;
  {-Return the Whittaker W function = exp(-x/2)*x^(0.5+m) * U(m-k-0.5,2m+1,x)}

function hyperg_0F1(b,x: double): double;
  {-Return the confluent hypergeometric limit function 0F1(;b;x)}

function hyperg_0F1r(b,x: double): double;
  {-Return the regularized confluent hypergeometric limit function 0F1(;b;x)/Gamma(b)}

function CylinderD(v,x: double): double;
  {-Return Whittaker's parabolic cylinder function D_v(x)}

function CylinderU(a,x: double): double;
  {-Return the parabolic cylinder function U(a,x)}

function CylinderV(a,x: double): double;
  {-Return the parabolic cylinder function V(a,x) with 2a integer}

function HermiteH(v,x: double): double;
  {-Return the Hermite function H_v(x) of degree v}

{#Z+}
{---------------------------------------------------------------------------}
{--------------------------- Other functions -------------------------------}
{---------------------------------------------------------------------------}
{#Z-}

function agm(x,y: double): double;
  {-Return the arithmetic-geometric mean of |x| and |y|}

function bernoulli(n: integer): double;
  {-Return the nth Bernoulli number, 0 if n<0 or odd n >= 3}

function bernpoly(n: integer; x: double): double;
  {-Return the Bernoulli polynomial B_n(x), 0 <= n <= MaxBernoulli}

function catalan(x: double): double;
  {-Return the Catalan function C(x) = binomial(2x,x)/(x+1)}

function cosint(n: integer; x: double): double;
  {-Return cosint(n, x) = integral(cos(t)^n, t=0..x), n>=0}

function debye(n: integer; x: double): double;
  {-Return the Debye function D(n,x) = n/x^n*integral(t^n/(exp(t)-1),t=0..x) of order n>0, x>=0}

function euler(n: integer): double;
  {-Return the nth Euler number, 0 if n<0 or odd n}

function fibpoly(n: integer; x: double): double;
  {-Return the Fibonacci polynomial F_n(x)}

function kepler(M,e: double): double;
  {-Solve Kepler's equation, result x is the eccentric anomaly from the mean anomaly M and the }
  { eccentricity e >= 0; x - e*sin(x) = M, x + x^3/3 = M, or e*sinh(x) - x = M for e <1, =1, >1}

function lucpoly(n: integer; x: double): double;
  {-Return the Lucas polynomial L_n(x)}

function LambertW(x: double): double;
  {-Return the Lambert W function (principal branch), x >= -1/e}

function LambertW1(x: double): double;
  {-Return the Lambert W function (-1 branch), -1/e <= x < 0}

function li_inv(x: double): double;
  {-Return the functional inverse of li(x), li(li_inv(x))=x}

function RiemannR(x: double): double;
  {-Return the Riemann prime counting function R(x), x >= 1/16}

function sinint(n: integer; x: double): double;
  {-Return sinint(n, x) = integral(sin(t)^n, t=0..x), n>=0}


implementation

uses
  DAMath,
  sdBessel,  {Bessel functions}
  sdEllInt,  {Elliptic integrals}
  sdExpInt,  {Exponential integrals and related}
  sdGamma,   {Gamma function and related}
  sdGamma2,  {Gamma function and related}
  sdZeta,    {Zeta functions and polylogarithms}
  sdZeta2,   {Zeta functions and polylogarithms}
  sdPoly,    {Orthogonal polynomials and related}
  sdHyperG,  {Hypergeometric functions}
  sdErf,     {Error function and related}
  sdSDist,   {Statistical distributions}
  sdMisc;    {Other functions}


{---------------------------------------------------------------------------}
{--------------------------- Bessel functions ------------------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function bessel_i0(x: double): double;
  {-Return I0(x), the modified Bessel function of the 1st kind, order zero}
begin
  bessel_i0 := sfd_i0(x);
end;


{---------------------------------------------------------------------------}
function bessel_i0e(x: double): double;
  {-Return I0(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order zero}
begin
  bessel_i0e := sfd_i0e(x);
end;


{---------------------------------------------------------------------------}
function bessel_i1(x: double): double;
  {-Return I1(x), the modified Bessel function of the 1st kind, order one}
begin
  bessel_i1 := sfd_i1(x);
end;


{---------------------------------------------------------------------------}
function bessel_i1e(x: double): double;
  {-Return I1(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order one}
begin
  bessel_i1e := sfd_i1e(x);
end;


{---------------------------------------------------------------------------}
function bessel_in(n: integer; x: double): double;
  {-Return I_n(x), the modified Bessel function of the 1st kind, order n.}
begin
  bessel_in := sfd_in(n,x);
end;


{---------------------------------------------------------------------------}
function bessel_iv(v, x: double): double;
  {-Return I_v(x), the modified Bessel function of the 1st kind, order v.}
begin
  bessel_iv := sfd_iv(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_ive(v, x: double): double;
  {-Return I_v(x)*exp(-|x|), the exponentially scaled modified Bessel function of the 1st kind, order v.}
begin
  bessel_ive := sfd_ive(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_j0(x: double): double;
  {-Return J0(x), the Bessel function of the 1st kind, order zero}
begin
  bessel_j0 := sfd_j0(x);
end;


{---------------------------------------------------------------------------}
function bessel_j1(x: double): double;
  {-Return J1(x), the Bessel function of the 1st kind, order one}
begin
  bessel_j1 := sfd_j1(x);
end;


{---------------------------------------------------------------------------}
function bessel_jn(n: integer; x: double): double;
  {-Return J_n(x), the Bessel function of the 1st kind, order n; not suitable for large n or x.}
begin
  bessel_jn := sfd_jn(n,x);
end;


{---------------------------------------------------------------------------}
function bessel_jv(v, x: double): double;
  {-Return J_v(x), the Bessel function of the 1st kind, order v; not suitable for large v.}
begin
  bessel_jv := sfd_jv(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_k0(x: double): double;
  {-Return K0(x), the modified Bessel function of the 2nd kind, order zero, x>0}
begin
  bessel_k0 := sfd_k0(x);
end;


{---------------------------------------------------------------------------}
function bessel_k0e(x: double): double;
  {-Return K0(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order zero, x>0}
begin
  bessel_k0e := sfd_k0e(x);
end;


{---------------------------------------------------------------------------}
function bessel_k1(x: double): double;
  {-Return K1(x), the modified Bessel function of the 2nd kind, order one, x>0}
begin
  bessel_k1 := sfd_k1(x);
end;


{---------------------------------------------------------------------------}
function bessel_k1e(x: double): double;
  {-Return K1(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order one, x>0}
begin
  bessel_k1e := sfd_k1e(x);
end;


{---------------------------------------------------------------------------}
function bessel_kn(n: integer; x: double): double;
  {-Return K_n(x), the modified Bessel function of the 2nd kind, order n, x>0, not suitable for large n}
begin
  bessel_kn := sfd_kn(n,x);
end;


{---------------------------------------------------------------------------}
function bessel_kv(v, x: double): double;
  {-Return K_v(x), the modified Bessel function of the 2nd kind, order v, x>0}
begin
  bessel_kv := sfd_kv(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_kve(v, x: double): double;
  {-Return K_v(x)*exp(x), the exponentially scaled modified Bessel function of the 2nd kind, order v, x>0}
begin
  bessel_kve := sfd_kve(v,x);
end;


{---------------------------------------------------------------------------}
function bessel_y0(x: double): double;
  {-Return Y0(x), the Bessel function of the 2nd kind, order zero; x>0}
begin
  bessel_y0 := sfd_y0(x);
end;


{---------------------------------------------------------------------------}
function bessel_y1(x: double): double;
  {-Return Y1(x), the Bessel function of the 2nd kind, order one; x>0}
begin
  bessel_y1 := sfd_y1(x);
end;


{---------------------------------------------------------------------------}
function bessel_yn(n: integer; x: double): double;
  {-Return Y_n(x), the Bessel function of the 2nd kind, order n, x>0, not suitable for large n or x}
begin
  bessel_yn := sfd_yn(n,x);
end;


{---------------------------------------------------------------------------}
function bessel_yv(v, x: double): double;
  {-Return Y_v(x), the Bessel function of the 2nd kind, order v; x > 0; not suitable for large v.}
begin
  bessel_yv := sfd_yv(v,x);
end;


{---------------------------------------------------------------------------}
function airy_ai(x: double): double;
  {-Return the Airy function Ai(x)}
begin
  airy_ai := sfd_airy_ai(x);
end;


{---------------------------------------------------------------------------}
function airy_aip(x: double): double;
  {-Return the Airy function Ai'(x)}
begin
  airy_aip := sfd_airy_aip(x);
end;


{---------------------------------------------------------------------------}
function airy_ais(x: double): double;
  {-Return the scaled Airy function Ai(x) if x <= 0, Ai(x)*exp(2/3*x^1.5) for x > 0}
begin
  airy_ais := sfd_airy_ais(x);
end;


{---------------------------------------------------------------------------}
function airy_bi(x: double): double;
  {-Return the Airy function Bi(x)}
begin
  airy_bi := sfd_airy_bi(x);
end;


{---------------------------------------------------------------------------}
function airy_bip(x: double): double;
  {-Return the Airy function Bi'(x)}
begin
  airy_bip := sfd_airy_bip(x);
end;


{---------------------------------------------------------------------------}
function airy_bis(x: double): double;
  {-Return the scaled Airy function Bi(x) if x <= 0, Bi(x)*exp(-2/3*x^1.5) for x > 0}
begin
  airy_bis := sfd_airy_bis(x);
end;


{---------------------------------------------------------------------------}
function airy_gi(x: double): double;
  {-Return the Airy/Scorer function Gi(x) = 1/Pi*integral(sin(x*t+t^3/3), t=0..INF)}
begin
  airy_gi := sfd_airy_gi(x);
end;


{---------------------------------------------------------------------------}
function airy_hi(x: double): double;
  {-Return the Airy/Scorer function Hi(x) = 1/Pi*integral(exp(x*t-t^3/3), t=0..INF)}
begin
  airy_hi := sfd_airy_hi(x);
end;


{---------------------------------------------------------------------------}
function kelvin_ber(x: double): double;
  {-Return the Kelvin function ber(x)}
begin
  kelvin_ber := sfd_ber(x);
end;


{---------------------------------------------------------------------------}
function kelvin_bei(x: double): double;
  {-Return the Kelvin function bei(x)}
begin
  kelvin_bei := sfd_bei(x);
end;


{---------------------------------------------------------------------------}
function kelvin_ker(x: double): double;
  {-Return the Kelvin function ker(x), x > 0}
begin
  kelvin_ker := sfd_ker(x);
end;


{---------------------------------------------------------------------------}
function kelvin_kei(x: double): double;
  {-Return the Kelvin function kei(x), x >= 0}
begin
  kelvin_kei := sfd_kei(x);
end;


{---------------------------------------------------------------------------}
procedure kelvin_kerkei(x: double; var kr, ki: double);
  {-Return the Kelvin functions kr=ker(x), ki=kei(x), x > 0}
begin
  sfd_ker_kei(x,kr,ki);
end;


{---------------------------------------------------------------------------}
procedure kelvin_berbei(x: double; var br, bi: double);
  {-Return the Kelvin functions br=ber(x), bi=bei(x)}
begin
  sfd_ber_bei(x,br,bi);
end;


{---------------------------------------------------------------------------}
function sph_bessel_jn(n: integer; x: double): double;
  {-Return j_n(x), the spherical Bessel function of the 1st kind, order n}
begin
  sph_bessel_jn := sfd_sph_jn(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_yn(n: integer; x: double): double;
  {-Return y_n(x), the spherical Bessel function of the 2nd kind, order n >=0 , x<>0}
begin
  sph_bessel_yn := sfd_sph_yn(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_in(n: integer; x: double): double;
  {-Return i_n(x), the modified spherical Bessel function of the 1st/2nd kind, order n}
begin
  sph_bessel_in := sfd_sph_in(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_ine(n: integer; x: double): double;
  {-Return i_n(x)*exp(-|x|), the exponentially scaled modified spherical Bessel function of the 1st/2nd kind, order n}
begin
  sph_bessel_ine := sfd_sph_ine(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_kn(n: integer; x: double): double;
  {-Return k_n(x), the modified spherical Bessel function of the 3rd kind, order n, x>0}
begin
  sph_bessel_kn := sfd_sph_kn(n,x);
end;


{---------------------------------------------------------------------------}
function sph_bessel_kne(n: integer; x: double): double;
  {-Return k_n(x)*exp(x), the exponentially scaled modified spherical Bessel function of the 3rd kind, order n, x>0}
begin
  sph_bessel_kne := sfd_sph_kne(n,x);
end;



{---------------------------------------------------------------------------}
function struve_h0(x: double): double;
 {-Return H0(x), the Struve function of order 0}
begin
  struve_h0 := sfd_struve_h0(x);
end;


{---------------------------------------------------------------------------}
function struve_h1(x: double): double;
  {-Return H1(x), the Struve function of order 1}
begin
  struve_h1 := sfd_struve_h1(x);
end;


{---------------------------------------------------------------------------}
function struve_h(v, x: double): double;
  {-Return H_v(x), the Struve function of order v, x < 0 only if v is an integer.}
begin
  struve_h := sfd_struve_h(v,x);
end;


{---------------------------------------------------------------------------}
function struve_l0(x: double): double;
  {-Return L0(x), the modified Struve function of order 0}
begin
  struve_l0 := sfd_struve_l0(x);
end;


{---------------------------------------------------------------------------}
function struve_l1(x: double): double;
  {-Return L1(x), the modified Struve function of order 1}
begin
  struve_l1 := sfd_struve_l1(x);
end;


{---------------------------------------------------------------------------}
function struve_l(v, x: double): double;
  {-Return L_v(x), the modified Struve function of order v, x < 0 only if v is an integer.}
begin
  struve_l := sfd_struve_l(v, x);
end;



{---------------------------------------------------------------------------}
{---------------------- Gamma function and related -------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function signgamma(x: double): double;
  {-Return sign(gamma(x)), useless for 0 or negative integer}
begin
  signgamma := sfd_signgamma(x);
end;


{---------------------------------------------------------------------------}
function gamma(x: double): double;
  {-Return gamma(x), x <= MAXGAM; invalid if x is a non-positive integer}
begin
  if x>MAXGAM then gamma := copysignd(PosInf_d, signgamma(x))
  else gamma := sfd_gamma(x);
end;


{---------------------------------------------------------------------------}
function gamma1pm1(x: double): double;
  {-Return gamma(1+x)-1 with increased accuracy for x near 0}
begin
  if x+1.0 > MAXGAM then gamma1pm1 := copysignd(PosInf_d, signgamma(x+1.0))
  else gamma1pm1 := sfd_gamma1pm1(x);
end;


{---------------------------------------------------------------------------}
function gammastar(x: double): double;
  {-Return Temme's gammastar(x) = gamma(x)/(sqrt(2*Pi)*x^(x-0.5)*exp(-x)), x>0.}
  { For large x the asymptotic expansion is gammastar(x) = 1 + 1/12x + O(1/x^2)}
begin
  gammastar := sfd_gstar(x);
end;


{---------------------------------------------------------------------------}
function gamma_delta_ratio(x,d: double): double;
  {-Return gamma(x)/gamma(x+d), accurate even for |d| << |x|}
begin
  gamma_delta_ratio := sfd_gamma_delta_ratio(x,d);
end;


{---------------------------------------------------------------------------}
function gamma_ratio(x,y: double): double;
  {-Return gamma(x)/gamma(y)}
begin
  gamma_ratio := sfd_gamma_ratio(x,y);
end;


{---------------------------------------------------------------------------}
function pochhammer(a,x: double): double;
  {-Return the Pochhammer symbol gamma(a+x)/gamma(a)}
begin
  pochhammer := sfd_pochhammer(a,x);
end;


{---------------------------------------------------------------------------}
function poch1(a,x: double): double;
  {-Return (pochhammer(a,x)-1)/x, psi(a) if x=0; accurate even for small |x|}
begin
  poch1 := sfd_poch1(a,x);
end;


{---------------------------------------------------------------------------}
function binomial(n,k: integer): double;
  {-Return the binomial coefficient 'n choose k'}
begin
  binomial := sfd_binomial(n,k);
end;


{---------------------------------------------------------------------------}
function fac(n: integer): double;
  {-Return the factorial n!, n<MAXGAM-1; INF if n<0}
begin
  if n > MAXGAM-1 then fac := PosInf_d
  else fac := sfd_fac(n);
end;


{---------------------------------------------------------------------------}
function dfac(n: integer): double;
  {-Return the double factorial n!!, n<=MAXDFAC; INF for even n<0}
begin
  if n > MAXDFAC then dfac := PosInf_d
  else dfac := sfd_dfac(n);
end;


{---------------------------------------------------------------------------}
function lnfac(n: longint): double;
  {-Return ln(n!), INF if n<0}
begin
  lnfac := sfd_lnfac(n);
end;


{---------------------------------------------------------------------------}
function lngamma(x: double): double;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, invalid if x is a non-positive integer}
  { Function signgamma can be used if the sign of gamma(x) is needed.}
begin
  if x=PosInf_d then lngamma := x
  else if abs(x) > MAXLGM then lngamma := PosInf_d
  else lngamma := sfd_lngamma(x);
end;


{---------------------------------------------------------------------------}
function lngammas(x: double; var s: integer): double;
  {-Return ln(|gamma(x)|), |x| <= MAXLGM, s=-1,1 is the sign of gamma}
begin
  lngammas := lngamma(x);
  if (x>0.0) or (frac(0.5*floord(x))=0.0) then s := 1 else s := -1;
end;


{---------------------------------------------------------------------------}
function lngamma1p(x: double): double;
  {-Return ln(|gamma(1+x)|) with increased accuracy for x near 0}
begin
  if IsNanOrInfD(x) then lngamma1p := x
  else if abs(x) > MAXLGM then lngamma1p := PosInf_d
  else lngamma1p := sfd_lngamma1p(x);
end;


{---------------------------------------------------------------------------}
function rgamma(x: double): double;
  {-Return the reciprocal gamma function rgamma = 1/gamma(x)}
begin
  if x>=179.0 then rgamma := 0.0
  else rgamma := sfd_rgamma(x);
end;


{---------------------------------------------------------------------------}
procedure incgamma(a,x: double; var p,q: double);
  {-Return the normalised incomplete gamma functions P and Q, a>=0, x>=0}
  { P(a,x) = integral(exp(-t)*t^(a-1), t=0..x  )/gamma(a)}
  { Q(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf)/gamma(a)}
begin
  sfd_incgamma(a,x,p,q);
end;


{---------------------------------------------------------------------------}
function igammap(a,x: double): double;
  {-Return the normalised lower incomplete gamma function P(a,x), a>=0, x>=0}
  { P(a,x) = integral(exp(-t)*t^(a-1), t=0..x)/gamma(a)}
begin
  igammap := sfd_igammap(a,x);
end;


{---------------------------------------------------------------------------}
function igammaq(a,x: double): double;
  {-Return the normalised lower incomplete gamma function Q(a,x), a>=0, x>=0}
  { Q(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf)/gamma(a)}
begin
  igammaq := sfd_igammaq(a,x);
end;


{---------------------------------------------------------------------------}
function igamma(a,x: double): double;
  {-Return the non-normalised upper incomplete gamma function}
  { GAMMA(a,x) = integral(exp(-t)*t^(a-1), t=x..Inf), x>=0}
begin
  igamma := sfd_igamma(a,x);
end;


{---------------------------------------------------------------------------}
function igammal(a,x: double): double;
  {-Return the non-normalised lower incomplete gamma function}
  { gamma(a,x) = integral(exp(-t)*t^(a-1), t=0..x); x>=0, a<>0,-1,-2,..}
begin
  igammal := sfd_igammal(a,x);
end;


{---------------------------------------------------------------------------}
function igammat(a,x: double): double;
  {-Return Tricomi's entire incomplete gamma function igammat(a,x)}
  { = igammal(a,x)/gamma(a)/x^a = P(a,x)/x^a }
begin
  igammat := sfd_git(a,x);
end;



{---------------------------------------------------------------------------}
procedure incgamma_inv(a,p,q: double; var x: double; var ierr: integer);
  {-Return the inverse normalised incomplete gamma function, i.e. calculate}
  { x with P(a,x)=p and Q(a,x)=q. Input parameter a>0, p>=0, q>0 and p+q=1.}
  { ierr is >= 0 for success, < 0 for input errors or iterations failures. }
begin
  sfd_incgamma_inv(a,p,q,x,ierr);
end;


{---------------------------------------------------------------------------}
function igamma_inv(a,p,q: double): double;
  {-Return the inverse normalised incomplete gamma function, i.e. calculate}
  { x with P(a,x)=p and Q(a,x)=q. Input parameter a>0, p>=0, q>0 and p+q=1.}
begin
  igamma_inv := sfd_igamma_inv(a,p,q);
end;


{---------------------------------------------------------------------------}
function igammap_inv(a,p: double): double;
  {-Inverse incomplete gamma: return x with P(a,x)=p, a>=0, 0<=p<1}
begin
  igammap_inv := sfd_igammap_inv(a,p);
end;


{---------------------------------------------------------------------------}
function igammaq_inv(a,q: double): double;
  {-Inverse complemented incomplete gamma: return x with Q(a,x)=q, a>=0, 0<q<=1}
begin
  igammaq_inv := sfd_igammaq_inv(a,q);
end;


{---------------------------------------------------------------------------}
function lngamma_inv(y: double): double;
  {-Inverse of lngamma: return x with lngamma(x) = y, y >= -0.12142, x > 1.4616}
begin
  lngamma_inv := sfd_ilng(y);
end;


{---------------------------------------------------------------------------}
function psi(x: double): double;
  {-Return the psi (digamma) function of x, INF if x is a non-positive integer}
begin
  psi := sfd_psi(x);
end;


{---------------------------------------------------------------------------}
function psi_inv(y: double): double;
  {-Inverse of psi, return x with psi(x)=y, y <= ln_MaxDbl}
begin
  psi_inv := sfd_ipsi(y);
end;


{---------------------------------------------------------------------------}
function trigamma(x: double): double;
  {-Return the trigamma function of x, INF if x is a negative integer}
begin
  trigamma := sfd_trigamma(x);
end;


{---------------------------------------------------------------------------}
function tetragamma(x: double): double;
  {-Return the tetragamma function psi''(x), NAN/RTE if x is a negative integer}
begin
  tetragamma := sfd_tetragamma(x);
end;


{---------------------------------------------------------------------------}
function pentagamma(x: double): double;
  {-Return the pentagamma function psi'''(x), INF if x is a negative integer}
begin
  pentagamma := sfd_pentagamma(x);
end;


{---------------------------------------------------------------------------}
function polygamma(n: integer; x: double): double;
  {-Return the polygamma function: n'th derivative of psi; n>=0, x>0 for n>MAXGAMD.}
  { Note: The accuracy may be reduced for n>=MAXGAMD due to ln/exp operations.}
begin
  polygamma := sfd_polygamma(n, x);
end;


{---------------------------------------------------------------------------}
function BatemanG(x: double): double;
  {-Return the Bateman function G(x) = psi((x+1)/2) - psi(x/2), x<>0,-1,-2,...}
begin
  BatemanG := sfd_batemang(x);
end;


{---------------------------------------------------------------------------}
function lnbeta(x,y: double): double;
  {-Return the logarithm of |beta(x,y)|=|gamma(x)*gamma(y)/gamma(x+y)|}
begin
  lnbeta := sfd_lnbeta(x,y);
end;


{---------------------------------------------------------------------------}
function beta(x,y: double): double;
  {-Return the function beta(x,y)=gamma(x)*gamma(y)/gamma(x+y)}
begin
  beta := sfd_beta(x,y);
end;


{---------------------------------------------------------------------------}
function ibeta(a, b, x: double): double;
  {-Return the normalised incomplete beta function, a>0, b>0, 0 <= x <= 1}
  { ibetax = integral(t^(a-1)*(1-t)^(b-1) / betax(a,b), t=0..x)}
begin
  ibeta := sfd_ibeta(a,b,x);
end;


{---------------------------------------------------------------------------}
function ibeta_inv(a, b, y: double): double;
  {-Return the functional inverse of the normalised incomplete beta function}
  { with a > 0, b > 0, and 0 <= y <= 1.}
begin
  ibeta_inv := sfd_ibeta_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function beta3(a, b, x: double): double;
  {-Return the non-normalised incomplete beta function B_x(a,b) for}
  { a>0, b>0, 0<=x<=1, B_x = integral(t^(a-1)*(1-t)^(b-1), t=0..x).}
begin
  beta3 := sfd_nnbeta(a,b,x);
end;


{---------------------------------------------------------------------------}
{------------------- Zeta functions and polylogarithms   -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function zeta(s: double): double;
  {-Return the Riemann zeta function at s, s<>1}
begin
  zeta := sfd_zeta(s);
end;


{---------------------------------------------------------------------------}
function zeta1p(x: double): double;
  {-Return the Riemann zeta function at 1+x, x<>0}
begin
  zeta1p := sfd_zeta1p(x);
end;


{---------------------------------------------------------------------------}
function eta(s: double): double;
  {-Return the Dirichlet eta function}
begin
  eta := sfd_eta(s);
end;


{---------------------------------------------------------------------------}
function etaint(n: integer): double;
  {-Return the Dirichlet function eta(n) for integer arguments}
begin
  etaint := sfd_etaint(n);
end;


{---------------------------------------------------------------------------}
function etam1(s: double): double;
  {-Return Dirichlet eta(s)-1}
begin
  etam1 := sfd_etam1(s);
end;


{---------------------------------------------------------------------------}
function zetam1(s: double): double;
  {-Return Riemann zeta(s)-1, s<>1}
begin
  zetam1 := sfd_zetam1(s);
end;


{---------------------------------------------------------------------------}
function zetah(s,a: double): double;
  {-Return the Hurwitz zeta function zetah(s,a) = sum(1/(i+a)^s, i=0..INF), s<>1, a>0}
begin
  zetah := sfd_zetah(s,a);
end;


{---------------------------------------------------------------------------}
function zetaint(n: integer): double;
  {-Return zeta(n) for integer arguments, n<>1}
begin
  zetaint := sfd_zetaint(n);
end;


{---------------------------------------------------------------------------}
function LerchPhi(z,s,a: double): double;
  {-Return the Lerch transcendent Phi(z,s,a) = sum(z^n/(n+a)^s, n=0..INF), |z|<=1, s >= -1, a>0; s<>1 if z=1}
begin
  LerchPhi := sfd_lerch(z,s,a);
end;


{---------------------------------------------------------------------------}
function DirichletBeta(s: double): double;
  {-Return the Dirichlet beta function sum((-1)^n/(2n+1)^s, n=0..INF)}
begin
  DirichletBeta := sfd_dbeta(s);
end;


{---------------------------------------------------------------------------}
function DirichletLambda(s: double): double;
  {-Return the Dirichlet lambda function sum(1/(2n+1)^s, n=0..INF), s<>1}
begin
  DirichletLambda := sfd_dlambda(s);
end;


{---------------------------------------------------------------------------}
function LegendreChi(s, x: double): double;
  {-Return Legendre's Chi-function chi(s,x); s>=0, |x|<=1, x<>1 if s<=1}
begin
  LegendreChi := sfd_lchi(s, x);
end;


{---------------------------------------------------------------------------}
function primezeta(x: double): double;
  {-Return the prime zeta function P(x) = sum(1/p^x, p prime), x > 1}
begin
  primezeta := sfd_pz(x);
end;


{---------------------------------------------------------------------------}
function harmonic(x: double): double;
  {-Return the harmonic function H(x) = psi(x+1) + EulerGamma}
begin
  harmonic := sfd_harmonic(x);
end;


{---------------------------------------------------------------------------}
function harmonic2(x,r: double): double;
  {-Return the generalized harmonic function H(x,r) = zeta(r)-zetah(r,x+1); x >= -1}
begin
  harmonic2 := sfd_harmonic2(x,r);
end;


{---------------------------------------------------------------------------}
function dilog(x: double): double;
  {-Return dilog(x) = Re(Li_2(x)), Li_2(x) = -integral(ln(1-t)/t, t=0..x)}
begin
  dilog := sfd_dilog(x);
end;


{---------------------------------------------------------------------------}
function trilog(x: double): double;
  {-Return the trilogarithm function trilog(x) = Re(Li_3(x))}
begin
  trilog := sfd_trilog(x);
end;


{---------------------------------------------------------------------------}
function polylog(n: integer; x: double): double;
  {-Return the polylogarithm Li_n(x) of integer order; x<1 for n >= 0}
begin
  polylog := sfd_polylog(n,x);
end;

{---------------------------------------------------------------------------}
function polylogr(s, x: double): double;
  {-Return the polylogarithm Li_s(x) of real order; s >= -1, |x|<=1, x<>1 if s=1}
begin
  polylogr := sfd_polylogr(s,x);
end;


{---------------------------------------------------------------------------}
function cl2(x: double): double;
  {-Return the Clausen function: integral(-ln(2*|sin(t/2)|),t=0..x) = Im(Li_2(exp(ix)))}
begin
  cl2 := sfd_cl2(x);
end;


{---------------------------------------------------------------------------}
function ti2(x: double): double;
  {-Return the inverse tangent integral, ti2(x) = integral(arctan(t)/t, t=0..x)}
begin
  ti2 := sfd_ti2(x);
end;


{---------------------------------------------------------------------------}
function lobachevsky_c(x: double): double;
  {-Return the Lobachevski function L(x) = integral(-ln(|cos(t)|), t=0..x)}
begin
  lobachevsky_c := sfd_llci(x);
end;


{---------------------------------------------------------------------------}
function lobachevsky_s(x: double): double;
  {-Return the Lobachevski function Lambda(x) = integral(-ln(|2sin(t)|), t=0..x)}
begin
  lobachevsky_s := sfd_llsi(x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac(n: integer; x: double): double;
  {-Return the integer order Fermi-Dirac integral F_n(x) = 1/n!*integral(t^n/(exp(t-x)+1), t=0..INF)}
begin
  fermi_dirac := sfd_fermi_dirac(n,x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac_m05(x: double): double;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}
begin
  fermi_dirac_m05 := sfd_fdm05(x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac_p05(x: double): double;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}
begin
  fermi_dirac_p05 := sfd_fdp05(x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac_p15(x: double): double;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}
begin
  fermi_dirac_p15 := sfd_fdp15(x);
end;


{---------------------------------------------------------------------------}
function fermi_dirac_p25(x: double): double;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}
begin
  fermi_dirac_p25 := sfd_fdp25(x);
end;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Legendre style) -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function comp_ellint_1(k: double): double;
  {-Return the complete elliptic integral of the 1st kind, |k| < 1}
begin
  comp_ellint_1 := sfd_EllipticK(k);
end;


{---------------------------------------------------------------------------}
function comp_ellint_2(k: double): double;
  {-Return the complete elliptic integral of the 2nd kind, |k| <= 1}
begin
  comp_ellint_2 := sfd_EllipticEC(k);
end;


{---------------------------------------------------------------------------}
function comp_ellint_3(nu,k: double): double;
  {-Return the complete elliptic integral of the 3rd kind, |k|<1, nu<>1}
begin
  comp_ellint_3 := sfd_EllipticPiC(nu,k);
end;


{---------------------------------------------------------------------------}
function comp_ellint_d(k: double): double;
  {-Return the complete elliptic integral D(k) = (K(k) - E(k))/k^2, |k| < 1}
begin
  comp_ellint_d := sfd_cel_d(k);
end;


{---------------------------------------------------------------------------}
function ellint_1(phi,k: double): double;
  {-Return the Legendre elliptic integral F(phi,k) of the 1st kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}
begin
  ellint_1 := sfd_ellint_1(phi,k);
end;


{---------------------------------------------------------------------------}
function ellint_2(phi,k: double): double;
  {-Return the Legendre elliptic integral E(phi,k) of the 2nd kind}
  { = integral(sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}
begin
  ellint_2 := sfd_ellint_2(phi,k);
end;


{---------------------------------------------------------------------------}
function ellint_3(phi,nu,k: double): double;
  {-Return the Legendre elliptic integral PI(phi,nu,k) of the 3rd kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2)/(1-nu*sin(x)^2),x=0..phi) with  }
  { |k*sin(phi)|<=1, returns Cauchy principal value if nu*sin(phi)^2>1}
begin
  ellint_3 := sfd_ellint_3(phi,nu,k);
end;


{---------------------------------------------------------------------------}
function ellint_d(phi,k: double): double;
  {-Return the Legendre elliptic integral D(phi,k) = (F(phi,k) - E(phi,k))/k^2 }
  { = integral(sin(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }
begin
  ellint_d := sfd_ellint_d(phi,k);
end;


{---------------------------------------------------------------------------}
function heuman_lambda(phi,k: double): double;
  {-Return Heuman's function Lambda_0(phi,k) = F(phi,k')/K(k') + 2/Pi*K(k)*Z(phi,k'), |k|<=1}
begin
  heuman_lambda := sfd_hlambda(phi,k);
end;


{---------------------------------------------------------------------------}
function jacobi_zeta(phi,k: double): double;
  {-Return the Jacobi Zeta function Z(phi,k) = E(phi,k) - E(k)/K(k)*F(phi,k), |k|<=1}
begin
  jacobi_zeta := sfd_jzeta(phi,k);
end;



{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Carlson style) --------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function ell_rc(x,y: double): double;
  {-Return Carlson's degenerate elliptic integral RC; x>=0, y<>0}
begin
  ell_rc := sfd_ell_rc(x,y);
end;


{---------------------------------------------------------------------------}
function ell_rf(x,y,z: double): double;
  {-Return Carlson's elliptic integral of the 1st kind; x,y,z >=0, at most one =0}
begin
  ell_rf := sfd_ell_rf(x,y,z);
end;


{---------------------------------------------------------------------------}
function ell_rd(x,y,z: double): double;
  {-Return Carlson's elliptic integral of the 2nd kind; z>0; x,y >=0, at most one =0}
begin
  ell_rd := sfd_ell_rd(x,y,z);
end;


{---------------------------------------------------------------------------}
function ell_rg(x,y,z: double): double;
  {-Return Carlson's completely symmetric elliptic integral of the 2nd kind; x,y,z >= 0}
begin
  ell_rg := sfd_ell_rg(x,y,z);
end;


{---------------------------------------------------------------------------}
function ell_rj(x,y,z,r: double): double;
  {-Return Carlson's elliptic integral of the 3rd kind; r<>0; x,y,z >=0, at most one =0}
begin
  ell_rj := sfd_ell_rj(x,y,z,r);
end;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Bulirsch style) -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function cel1(kc: double): double;
  {-Return Bulirsch's complete elliptic integral of the 1st kind, kc<>0}
begin
  cel1 := sfd_cel1(kc);
end;


{---------------------------------------------------------------------------}
function cel2(kc, a, b: double): double;
  {-Return Bulirsch's complete elliptic integral of the 2nd kind, kc<>0}
begin
  cel2 := sfd_cel2(kc,a,b);
end;


{---------------------------------------------------------------------------}
function cel(kc, p, a, b: double): double;
  {-Return Bulirsch's general complete elliptic integral, kc<>0, Cauchy principle value if p<0}
begin
  cel := sfd_cel(kc,p,a,b);
end;


{---------------------------------------------------------------------------}
function el1(x,kc: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 1st kind}
begin
  el1 := sfd_el1(x,kc);
end;


{---------------------------------------------------------------------------}
function el2(x,kc,a,b: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 2nd kind, kc<>0}
begin
  el2 := sfd_el2(x,kc,a,b);
end;


{---------------------------------------------------------------------------}
function el3(x,kc,p: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 3rd kind, 1+p*x^2<>0}
begin
  el3 := sfd_el3(x,kc,p);
end;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Maple V style) --------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function EllipticF(z,k: double): double;
  {-Return the incomplete elliptic integral of the 1st kind; |z|<=1, |k*z|<1}
begin
  EllipticF := sfd_EllipticF(z,k);
end;


{---------------------------------------------------------------------------}
function EllipticK(k: double): double;
  {-Return the complete elliptic integral of the 1st kind, |k| < 1}
begin
  EllipticK := sfd_EllipticK(k);
end;


{---------------------------------------------------------------------------}
function EllipticKim(k: double): double;
  {-Return K(i*k), the complete elliptic integral of the 1st kind with}
  { imaginary modulus = integral(1/sqrt(1-x^2)/sqrt(1+k^2*x^2),x=0..1)}
begin
  EllipticKim := sfd_EllipticKim(k);
end;


{---------------------------------------------------------------------------}
function EllipticCK(k: double): double;
  {-Return the complementary complete elliptic integral of the 1st kind, k<>0}
begin
  EllipticCK := sfd_EllipticCK(k);
end;


{---------------------------------------------------------------------------}
function EllipticE(z,k: double): double;
  {-Return the incomplete elliptic integrals of the 2nd kind, |z|<=1, |k*z| <= 1}
begin
  EllipticE := sfd_EllipticE(z,k);
end;


{---------------------------------------------------------------------------}
function EllipticEC(k: double): double;
  {-Return the complete elliptic integral of the 2nd kind, |k| < 1}
begin
  EllipticEC := sfd_EllipticEC(k);
end;


{---------------------------------------------------------------------------}
function EllipticECim(k: double): double;
  {-Return E(i*k), the complete elliptic integral of the 2nd kind with}
  { imaginary modulus = integral(sqrt(1+k^2*x^2)/sqrt(1-x^2),x=0..1)  }
begin
  EllipticECim := sfd_EllipticECim(k);
end;


{---------------------------------------------------------------------------}
function EllipticCE(k: double): double;
  {-Return the complementary complete elliptic integral of the 2nd kind}
begin
  EllipticCE := sfd_EllipticCE(k);
end;


{---------------------------------------------------------------------------}
function EllipticPi(z,nu,k: double): double;
  {-Return the incomplete elliptic integral of the 3rd kind, |z|<=1, |k*z|<1}
begin
  EllipticPi := sfd_EllipticPi(z,nu,k);
end;


{---------------------------------------------------------------------------}
function EllipticPiC(nu,k: double): double;
  {-Return the complete elliptic integral of the 3rd kind, |k|<1}
begin
  EllipticPiC := sfd_EllipticPiC(nu,k);
end;


{---------------------------------------------------------------------------}
function EllipticCPi(nu,k: double): double;
  {-Return the complementary complete elliptic integral of the 3rd kind, k<>0, nu<>1}
begin
  EllipticCPi := sfd_EllipticCPi(nu,k);
end;


{---------------------------------------------------------------------------}
{------------------- Jacobi elliptic and theta functions -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function EllipticModulus(q: double): double;
  {-Return the elliptic modulus k(q) = theta_2(q)^2/theta_3(q)^2, 0 <= q <= 1}
begin
  EllipticModulus := sfd_ellmod(q);
end;


{---------------------------------------------------------------------------}
function EllipticNome(k: double): double;
  {-Return the elliptic nome q(k) = exp(-Pi*EllipticCK(k)/EllipticK(k)), |k| < 1}
begin
  EllipticNome := sfd_nome(k);
end;


{---------------------------------------------------------------------------}
function jacobi_am(x,k: double): double;
  {-Return the Jacobi amplitude am(x,k)}
begin
  jacobi_am := sfd_jam(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_theta(n: integer; x,q: double): double;
  {-Return the Jacobi theta function theta_n(x,q), n=1..4, 0 <= q < 1}
begin
  jacobi_theta := sfd_jtheta(n,x,q);
end;


{---------------------------------------------------------------------------}
procedure sncndn(x,mc: double; var sn,cn,dn: double);
  {-Return the Jacobi elliptic functions sn,cn,dn for argument x and}
  { complementary parameter mc.}
begin
  sfd_sncndn(x,mc,sn,cn,dn);
end;


{---------------------------------------------------------------------------}
function jacobi_sn(x,k: double): double;
  {-Return the Jacobi elliptic function sn(x,k)}
begin
  jacobi_sn := sfd_jacobi_sn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_cn(x,k: double): double;
  {-Return the Jacobi elliptic function cn(x,k)}
begin
  jacobi_cn := sfd_jacobi_cn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_dn(x,k: double): double;
  {-Return the Jacobi elliptic function dn(x,k)}
begin
  jacobi_dn := sfd_jacobi_dn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_nc(x,k: double): double;
  {-Return the Jacobi elliptic function nc(x,k)}
begin
  jacobi_nc := sfd_jacobi_nc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_sc(x,k: double): double;
  {-Return the Jacobi elliptic function sc(x,k)}
begin
  jacobi_sc := sfd_jacobi_sc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_dc(x,k: double): double;
  {-Return the Jacobi elliptic function dc(x,k)}
begin
  jacobi_dc := sfd_jacobi_dc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_nd(x,k: double): double;
  {-Return the Jacobi elliptic function nd(x,k)}
begin
  jacobi_nd := sfd_jacobi_nd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_sd(x,k: double): double;
  {-Return the Jacobi elliptic function sd(x,k)}
begin
  jacobi_sd := sfd_jacobi_sd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_cd(x,k: double): double;
  {-Return the Jacobi elliptic function cd(x,k)}
begin
  jacobi_cd := sfd_jacobi_cd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_ns(x,k: double): double;
  {-Return the Jacobi elliptic function ns(x,k)}
begin
  jacobi_ns := sfd_jacobi_ns(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_cs(x,k: double): double;
  {-Return the Jacobi elliptic function cs(x,k)}
begin
  jacobi_cs := sfd_jacobi_cs(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_ds(x,k: double): double;
  {-Return the Jacobi elliptic function ds(x,k)}
begin
  jacobi_ds := sfd_jacobi_ds(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arccn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccn(x,k), |x| <= 1, x >= sqrt(1 - 1/k^2) if k >= 1}
begin
  jacobi_arccn := sfd_jacobi_arccn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arccd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccd(x,k); |x| <= 1 if |k| < 1; |x| >= 1 if |k| > 1 }
begin
  jacobi_arccd := sfd_jacobi_arccd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arccs(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccs(x,k), |x| >= sqrt(k^2-1) for |k|>1}
begin
  jacobi_arccs := sfd_jacobi_arccs(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcdc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcdc(x,k); |x| >= 1 if |k| < 1; |x| <= 1 if |k| > 1 }
begin
  jacobi_arcdc := sfd_jacobi_arcdc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcdn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcdn(x,k), 0 <= x <= 1, x^2 + k^2 > 1 if |k| < 1;  |x| <= 1 if |k| > 1}
begin
  jacobi_arcdn := sfd_jacobi_arcdn(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcds(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcds(x,k), x^2 + k^2 >= 1}
begin
  jacobi_arcds := sfd_jacobi_arcds(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcnc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcnc(x,k), x >= 1, x^2 <= k^2/(k^2-1) for |k|>1}
begin
  jacobi_arcnc := sfd_jacobi_arcnc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcnd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcnd(x,k), x >= 1, x^2 <= k^2/(1-k^2) if k < 1}
begin
  jacobi_arcnd := sfd_jacobi_arcnd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcns(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcns(x,k), |x| >= 1, |x| >= k if k>=1}
begin
  jacobi_arcns := sfd_jacobi_arcns(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcsc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsc(x,k), |x| <= 1/sqrt(k^2-1) for |k|>1}
begin
  jacobi_arcsc := sfd_jacobi_arcsc(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcsd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsd(x,k), x^2*(1-k^2) <= 1}
begin
  jacobi_arcsd := sfd_jacobi_arcsd(x,k);
end;


{---------------------------------------------------------------------------}
function jacobi_arcsn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsn(x,k), |x| <= 1 and |x*k| <= 1}
begin
  jacobi_arcsn := sfd_jacobi_arcsn(x,k);
end;


{---------------------------------------------------------------------------}
function theta1p(q: double): double;
  {-Return the derivative  theta1p(q) := d/dx(theta_1(x,q)) at x=0,}
  { = 2*q^(1/4)*sum((-1)^n*(2n+1)*q^(n*(n+1)),n=0..Inf), 0 <= q < 1}
begin
  theta1p := sfd_theta1p(q);
end;


{---------------------------------------------------------------------------}
function theta2(q: double): double;
  {-Return Jacobi theta_2(q) = 2*q^(1/4)*sum(q^(n*(n+1)),n=0..Inf) 0 <= q < 1}
begin
  theta2 := sfd_theta2(q);
end;


{---------------------------------------------------------------------------}
function theta3(q: double): double;
  {-Return Jacobi theta_3(q) = 1 + 2*sum(q^(n*n)),n=1..Inf); |q| < 1}
begin
  theta3 := sfd_theta3(q);
end;


{---------------------------------------------------------------------------}
function theta4(q: double): double;
  {-Return Jacobi theta_4(q) = 1 + 2*sum((-1)^n*q^(n*(n+1)),n=1..Inf); |q| < 1}
begin
  theta4 := sfd_theta4(q);
end;


{---------------------------------------------------------------------------}
procedure sincos_lemn(x: double; var sl,cl: double);
  {-Return the lemniscate functions sl = sin_lemn(x), cl = cos_lemn(x)}
begin
  sfd_lemn(x,sl,cl);
end;


{---------------------------------------------------------------------------}
function sin_lemn(x: double): double;
  {-Return the lemniscate sine functions sl = sin_lemn(x)}
var
  sl,cl: double;
begin
  sfd_lemn(x,sl,cl);
  sin_lemn := sl;
end;


{---------------------------------------------------------------------------}
function cos_lemn(x: double): double;
  {-Return the lemniscate cosine functions cl = cos_lemn(x)}
var
  sl,cl: double;
begin
  sfd_lemn(x,sl,cl);
  cos_lemn := cl;
end;



{---------------------------------------------------------------------------}
{----------------------- Error function and related ------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function dawson(x: double): double;
  {-Return Dawson's integral: dawson(x) = exp(-x^2)*integral(exp(t^2), t=0..x)}
begin
  dawson := sfd_dawson(x);
end;


{---------------------------------------------------------------------------}
function dawson2(p,x: double): double;
  {-Return the generalized Dawson integral F(p,x) = exp(-x^p)*integral(exp(t^p), t=0..x); x,p >= 0}
begin
  dawson2 := sfd_gendaw(p,x);
end;


{---------------------------------------------------------------------------}
function erf(x: double): double;
  {-Return the error function erf(x) = 2/sqrt(Pi)*integral((exp(-t^2), t=0..x)}
begin
  erf := sfd_erf(x);
end;


{---------------------------------------------------------------------------}
function erfc(x: double): double;
  {-Return the complementary error function erfc(x) = 1-erf(x)}
begin
  erfc := sfd_erfc(x);
end;


{---------------------------------------------------------------------------}
function erfce(x: double): double;
  {-Return the exponentially scaled complementary error function erfce(x) = exp(x^2)*erfc(x)}
begin
  erfce := sfd_erfce(x);
end;


{---------------------------------------------------------------------------}
function inerfc(n: integer; x: double): double;
  {-Return the repeated integrals of erfc, n >= -1; scaled with exp(x^2) for x>0}
begin
  inerfc := sfd_inerfc(n,x);
end;


{---------------------------------------------------------------------------}
function erfg(p,x: double): double;
  {-Return the generalized error function integral(exp(-t^p), t=0..x); x,p >= 0}
begin
  erfg := sfd_erfg(p,x);
end;


{---------------------------------------------------------------------------}
function erfi(x: double): double;
  {-Return the imaginary error function erfi(x) = erf(ix)/i}
begin
  erfi := sfd_erfi(x);
end;


{---------------------------------------------------------------------------}
function erf_inv(x: double): double;
  {-Return the inverse function of erf, erf(erf_inv(x)) = x, -1 < x < 1}
begin
  erf_inv := sfd_erf_inv(x);
end;


{---------------------------------------------------------------------------}
function erfc_inv(x: double): double;
  {-Return the inverse function of erfc, erfc(erfc_inv(x)) = x, 0 < x < 2}
begin
  erfc_inv := sfd_erfc_inv(x);
end;


{---------------------------------------------------------------------------}
function erf_p(x: double): double;
  {-Return the probability function erf_p = integral(exp(-t^2/2)/sqrt(2*Pi), t=-Inf..x)}
begin
  erf_p := sfd_erf_p(x);
end;


{---------------------------------------------------------------------------}
function erf_q(x: double): double;
  {-Return the probability function erf_q = integral(exp(-t^2/2)/sqrt(2*Pi), t=x..Inf)}
begin
  erf_q := sfd_erf_q(x);
end;


{---------------------------------------------------------------------------}
function erf_z(x: double): double;
  {-Return the probability function erf_z = exp(-x^2/2)/sqrt(2*Pi)}
begin
  erf_z := sfd_erf_z(x);
end;



{---------------------------------------------------------------------------}
function expint3(x: double): double;
  {-Return the integral(exp(-t^3), t=0..x), x >= 0}
begin
  expint3 := sfd_expint3(x);
end;


{---------------------------------------------------------------------------}
procedure Fresnel(x: double; var s,c: double);
  {-Return the Fresnel integrals S(x)=integral(sin(Pi/2*t^2),t=0..x) and C(x)=integral(cos(Pi/2*t^2),t=0..x)}
begin
  sfd_fresnel(x,s,c);
end;


{---------------------------------------------------------------------------}
function FresnelC(x: double): double;
  {-Return the Fresnel integral C(x)=integral(cos(Pi/2*t^2),t=0..x)}
var
  sx,cx: double;
begin
  sfd_fresnel(x,sx,cx);
  FresnelC := cx;
end;


{---------------------------------------------------------------------------}
function FresnelS(x: double): double;
  {-Return the Fresnel integral S(x)=integral(sin(Pi/2*t^2),t=0..x)}
var
  sx,cx: double;
begin
  sfd_fresnel(x,sx,cx);
  FresnelS := sx;
end;


{---------------------------------------------------------------------------}
function gsi(x: double): double;
  {-Return the Goodwin-Staton integral = integral(exp(-t*t)/(t+x), t=0..Inf), x > 0}
begin
  gsi := sfd_gsi(x);
end;


{---------------------------------------------------------------------------}
{------------------- Exponential integrals and related ---------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function chi(x: double): double;
  {-Return the hyperbolic cosine integral = EulerGamma + ln(|x|) + integral((cosh(t)-1)/t, t=0..|x|)}
begin
  chi := sfd_chi(x);
end;


{---------------------------------------------------------------------------}
function ci(x: double): double;
  {-Return the cosine integral, ci(x) = EulerGamma + ln(|x|) + integral((cos(t)-1)/t, t=0..|x|)}
begin
  ci := sfd_ci(x);
end;


{---------------------------------------------------------------------------}
function cin(x: double): double;
  {-Return the entire cosine integral, cin(x) = integral((1-cos(t))/t, t=0..x)}
begin
  cin := sfd_cin(x);
end;


{---------------------------------------------------------------------------}
function cinh(x: double): double;
  {-Return the entire hyperbolic cosine integral, cinh(x) = integral((cosh(t)-1)/t, t=0..x)}
begin
  cinh := sfd_cinh(x);
end;


{---------------------------------------------------------------------------}
function e1(x: double): double;
  {-Return the exponential integral E1(x) = integral(exp(-x*t)/t, t=1..Inf), x <> 0}
begin
  if x>=740.0 then e1 := 0.0
  else e1 := sfd_e1(x);
end;


{---------------------------------------------------------------------------}
function ei(x: double): double;
  {-Return the exponential integral Ei(x) = PV-integral(exp(t)/t, t=-Inf..x)}
begin
  ei := sfd_ei(x);
end;


{---------------------------------------------------------------------------}
function ein(x: double): double;
  {-Return the entire exponential integral ein(x) = integral((1-exp(-t))/t, t=0..x)}
begin
  ein := sfd_ein(x);
end;


{---------------------------------------------------------------------------}
function ei_inv(x: double): double;
  {-Return the functional inverse of Ei(x), ei_inv(ei(x))=x}
begin
  ei_inv := sfd_ei_inv(x);
end;


{---------------------------------------------------------------------------}
function en(n: longint; x: double): double;
  {-Return the exponential integral E_n(x) = integral(exp(-x*t)/t^n, t=1..Inf), x > 0}
begin
  en := sfd_en(n,x);
end;


{---------------------------------------------------------------------------}
function gei(p,x: double): double;
  {-Return the generalized exponential integral E_p(x) = integral(exp(-x*t)/t^p, t=1..Inf), x >= 0}
begin
  gei := sfd_gei(p,x);
end;


{---------------------------------------------------------------------------}
function li(x: double): double;
  {-Return the logarithmic integral li(x) = PV-integral(1/ln(t), t=0..x), x >= 0, x <> 1}
begin
  li := sfd_li(x);
end;


{---------------------------------------------------------------------------}
function shi(x: double): double;
  {-Return the hyperbolic sine integral, integral(sinh(t)/t, t=0..x)}
begin
  shi := sfd_shi(x);
end;


{---------------------------------------------------------------------------}
function si(x: double): double;
  {-Return the sine integral, si(x) = integral(sin(t)/t, t=0..x)}
begin
  si := sfd_si(x);
end;


{---------------------------------------------------------------------------}
function ssi(x: double): double;
  {-Return the shifted sine integral, ssi(x) = si(x) - Pi/2}
begin
  ssi := sfd_ssi(x);
end;


{---------------------------------------------------------------------------}
{---------------------- Statistical distributions --------------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function beta_pdf(a, b, x: double): double;
  {-Return the probability density function of the beta distribution with}
  { parameters a and b: beta_pdf = x^(a-1)*(1-x)^(b-1) / beta(a,b)}
begin
  beta_pdf := sfd_beta_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function beta_cdf(a, b, x: double): double;
  {-Return the cumulative beta distribution function, a>0, b>0}
begin
  beta_cdf := sfd_beta_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function beta_inv(a, b, y: double): double;
  {-Return the functional inverse of the beta distribution function. a>0, b>0;}
  { 0 <= y <= 1. Given y the function finds x such that beta_cdf(a, b, x) = y}
begin
  beta_inv := sfd_beta_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function binomial_cdf(p: double; n, k: longint): double;
  {-Return the cumulative binomial distribution function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}
begin
  binomial_cdf := sfd_binomial_cdf(p,n,k);
end;


{---------------------------------------------------------------------------}
function binomial_pmf(p: double; n, k: longint): double;
  {-Return the binomial distribution probability mass function with number}
  { of trials n >= 0 and success probability 0 <= p <= 1}
begin
  binomial_pmf := sfd_binomial_pmf(p,n,k);
end;


{---------------------------------------------------------------------------}
function cauchy_pdf(a, b, x: double): double;
  {-Return the Cauchy probability density function with location a }
  { and scale b > 0, 1/(Pi*b*(1+((x-a)/b)^2))}
begin
  cauchy_pdf := sfd_cauchy_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function cauchy_cdf(a, b, x: double): double;
  {-Return the cumulative Cauchy distribution function with location a}
  { and scale b > 0, = 1/2 + arctan((x-a)/b)/Pi}
begin
  cauchy_cdf := sfd_cauchy_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function cauchy_inv(a, b, y: double): double;
  {-Return the functional inverse of Cauchy distribution function}
  { with location a and scale b > 0}
begin
  cauchy_inv := sfd_cauchy_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function chi2_pdf(nu: longint; x: double): double;
  {-Return the probability density function of the chi-square distribution, nu>0}
begin
  chi2_pdf := sfd_chi2_pdf(nu,x);
end;


{---------------------------------------------------------------------------}
function chi2_cdf(nu: longint; x: double): double;
  {-Return the cumulative chi-square distribution with nu>0 degrees of freedom, x >= 0}
begin
  chi2_cdf := sfd_chi2_cdf(nu,x);
end;


{---------------------------------------------------------------------------}
function chi2_inv(nu: longint; p: double): double;
  {-Return the functional inverse of the chi-square distribution, nu>0, 0 <= p < 1}
begin
  chi2_inv := sfd_chi2_inv(nu,p);
end;


{---------------------------------------------------------------------------}
function evt1_pdf(a, b, x: double): double;
  {-Return the probability density function of the Extreme Value Type I distribution}
  { with location a and scale b > 0, result = exp(-(x-a)/b)/b * exp(-exp(-(x-a)/b)) }
begin
  evt1_pdf := sfd_evt1_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function evt1_cdf(a, b, x: double): double;
  {-Return the cumulative Extreme Value Type I distribution function}
  { with location a and scale b > 0;  result = exp(-exp(-(x-a)/b)). }
begin
  evt1_cdf := sfd_evt1_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function evt1_inv(a, b, y: double): double;
  {-Return the functional inverse of the Extreme Value Type I distribution}
  { function with location a and scale b > 0;  result = a - b*ln(ln(-y)). }
begin
  evt1_inv := sfd_evt1_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function exp_pdf(a, alpha, x: double): double;
  {-Return the exponential probability density function with location a }
  { and rate alpha > 0, = alpha*exp(-alpha*(x-a)) if x >= a, 0 if x < a.}
begin
  exp_pdf := sfd_exp_pdf(a, alpha, x);
end;


{---------------------------------------------------------------------------}
function exp_cdf(a, alpha, x: double): double;
  {-Return the cumulative exponential distribution function with location a}
  { and rate alpha > 0, = 1 - exp(-alpha*(x-a)) if x >= a, 0 if x < a.}
begin
  exp_cdf := sfd_exp_cdf(a, alpha, x);
end;


{---------------------------------------------------------------------------}
function exp_inv(a, alpha, y: double): double;
  {-Return the functional inverse of exponential distribution function with}
  { location a and rate alpha > 0}
begin
  exp_inv := sfd_exp_inv(a, alpha, y);
end;


{---------------------------------------------------------------------------}
function f_pdf(nu1, nu2: longint; x: double): double;
  {-Return the probability density function of the F distribution; x >= 0, nu1, nu2 > 0}
begin
  f_pdf := sfd_f_pdf(nu1, nu2, x);
end;


{---------------------------------------------------------------------------}
function f_cdf(nu1, nu2: longint; x: double): double;
  {-Return the cumulative F distribution function; x >= 0, nu1, nu2 > 0}
begin
  f_cdf := sfd_f_cdf(nu1, nu2, x);
end;


{---------------------------------------------------------------------------}
function f_inv(nu1, nu2: longint; y: double): double;
  {-Return the functional inverse of the F distribution, nu1, nu2 > 0, 0 <= y <= 1}
begin
  f_inv := sfd_f_inv(nu1, nu2, y);
end;


{---------------------------------------------------------------------------}
function gamma_pdf(a, b, x: double): double;
  {-Return the probability density function of a gamma distribution with shape}
  { a>0, scale b>0: gamma_pdf = x^(a-1)*exp(-x/b)/gamma(a)/b^a, x>0}
begin
  gamma_pdf := sfd_gamma_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function gamma_cdf(a, b, x: double): double;
  {-Return the cumulative gamma distribution function, shape a>0, scale b>0}
begin
  gamma_cdf := sfd_gamma_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function hypergeo_pmf(n1,n2,n,k: longint): double;
  {-Return the hypergeometric distribution probability mass function; n,n1,n2 >= 0, n <= n1+n2;}
  { i.e. the probability that among n randomly chosen samples from a container}
  { with n1 type1 objects and n2 type2 objects are exactly k type1 objects:}
begin
  hypergeo_pmf := sfd_hypergeo_pmf(n1,n2,n,k);
end;


{---------------------------------------------------------------------------}
function hypergeo_cdf(n1,n2,n,k: longint): double;
  {-Return the cumulative hypergeometric distribution function; n,n1,n2 >= 0, n <= n1+n2}
begin
  hypergeo_cdf := sfd_hypergeo_cdf(n1,n2,n,k);
end;


{---------------------------------------------------------------------------}
function invgamma_pdf(a, b, x: double): double;
  {-Return the probability density function of an inverse gamma distribution}
  { with shape a>0, scale b>0: result = (b/x)^a/x*exp(-b/x)/Gamma(a), x >= 0}
begin
  invgamma_pdf := sfd_invgamma_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function invgamma_cdf(a, b, x: double): double;
  {-Return the cumulative inverse gamma distribution function, shape a>0, scale}
  { b>0: result = Gamma(a,b/x)/Gamma(a) = Q(a,b/x) = igammaq(a,b/x), x >= 0}
begin
  invgamma_cdf := sfd_invgamma_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function invgamma_inv(a, b, y: double): double;
  {-Return the functional inverse of the inverse gamma distribution function, shape}
  { a>0, scale b>0, 0 <= y <= 1, i.e. find x such that invgamma_cdf(a, b, x) = y  }
begin
  invgamma_inv := sfd_invgamma_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function kumaraswamy_pdf(a, b, x: double): double;
  {-Return the Kumaraswamy probability density function with shape}
  { parameters a,b>0, 0<=x<=1; result = a*b*x^(a-1)*(1-x^a)^(b-1) }
begin
  kumaraswamy_pdf := sfd_kumaraswamy_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function kumaraswamy_cdf(a, b, x: double): double;
  {-Return the cumulative Kumaraswamy distribution function with}
  { shape parameters a,b > 0, 0 <= x <= 1; result = 1-(1-x^a)^b}
begin
  kumaraswamy_cdf := sfd_kumaraswamy_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function kumaraswamy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Kumaraswamy distribution}
  { with shape parameters a,b > 0; result = [1-(1-y)^(1/b)]^(1/a)}
begin
  kumaraswamy_inv := sfd_kumaraswamy_inv(a,b,y)
end;


{---------------------------------------------------------------------------}
function laplace_cdf(a, b, x: double): double;
  {-Return the cumulative Laplace distribution function with location a and scale b > 0}
begin
  laplace_cdf := sfd_laplace_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function laplace_inv(a, b, y: double): double;
  {-Return the functional inverse of the Laplace distribution with location a and scale b > 0}
begin
  laplace_inv := sfd_laplace_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function laplace_pdf(a, b, x: double): double;
  {-Return the Laplace probability density function with location a}
  { and scale b > 0, result = exp(-abs(x-a)/b) / (2*b)}
begin
  laplace_pdf := sfd_laplace_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function gamma_inv(a, b, p: double): double;
  {-Return the functional inverse of the gamma distribution function, shape a>0,}
  { scale b>0, 0 <= p <= 1, i.e. finds x such that gamma_cdf(a, b, x) = p}
begin
  gamma_inv := sfd_gamma_inv(a, b, p);
end;


{---------------------------------------------------------------------------}
function logistic_pdf(a, b, x: double): double;
  {-Return the logistic probability density function with location a}
  { and scale parameter b > 0, exp(-(x-a)/b)/b/(1+exp(-(x-a)/b))^2}
begin
  logistic_pdf := sfd_logistic_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function logistic_cdf(a, b, x: double): double;
  {-Return the cumulative logistic distribution function with}
  { location a and scale parameter b > 0}
begin
  logistic_cdf := sfd_logistic_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function logistic_inv(a, b, y: double): double;
  {-Return the functional inverse of the logistic distribution}
  { with location a and scale parameter b > 0}
begin
  logistic_inv := sfd_logistic_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function levy_pdf(a, b, x: double): double;
  {-Return the Levy probability density function with}
  { location a and scale parameter b > 0}
begin
  levy_pdf := sfd_levy_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function levy_cdf(a, b, x: double): double;
  {-Return the cumulative Levy distribution function with}
  { location a and scale parameter b > 0}
begin
  levy_cdf := sfd_levy_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function levy_inv(a, b, y: double): double;
  {-Return the functional inverse of the Levy distribution}
  { with location a and scale parameter b > 0}
begin
  levy_inv := sfd_levy_inv(a, b, y);
end;


{---------------------------------------------------------------------------}
function lognormal_pdf(a, b, x: double): double;
  {-Return the log-normal probability density function with}
  { location a and scale parameter b > 0, zero for x <= 0.}
begin
  lognormal_pdf := sfd_lognormal_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function lognormal_cdf(a, b, x: double): double;
  {-Return the cumulative log-normal distribution function with}
  { location a and scale parameter b > 0, zero for x <= 0.}
begin
  lognormal_cdf := sfd_lognormal_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function lognormal_inv(a, b, y: double): double;
  {-Return the functional inverse of the log-normal distribution}
  { with location a and scale parameter b > 0, 0 < y < 1.}
begin
  lognormal_inv := sfd_lognormal_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function logseries_pmf(a: double; k: longint): double;
  {-Return the logarithmic (series) probability mass function}
  { with shape 0 < a < 1, k > 0; result = -a^k/(k*ln(1-a))   }
begin
  logseries_pmf := sfd_ls_pmf(a,k);
end;


{---------------------------------------------------------------------------}
function logseries_cdf(a: double; k: longint): double;
  {-Return the cumulative logarithmic (series) distribution function with shape 0 < a < 1, k > 0}
begin
  logseries_cdf := sfd_ls_cdf(a,k);
end;


{---------------------------------------------------------------------------}
function maxwell_pdf(b, x: double): double;
  {-Return the Maxwell probability density function with scale b > 0, x >= 0}
begin
  maxwell_pdf := sfd_maxwell_pdf(b,x);
end;


{---------------------------------------------------------------------------}
function maxwell_cdf(b, x: double): double;
  {-Return the cumulative Maxwell distribution function with scale b > 0, x >= 0}
begin
  maxwell_cdf := sfd_maxwell_cdf(b,x);
end;


{---------------------------------------------------------------------------}
function maxwell_inv(b, y: double): double;
  {-Return the functional inverse of the Maxwell distribution with scale b > 0}
begin
  maxwell_inv := sfd_maxwell_inv(b,y);
end;


{---------------------------------------------------------------------------}
function moyal_pdf(a, b, x: double): double;
  {-Return the Moyal probability density function with}
  { location a and scale parameter b > 0}
begin
  moyal_pdf := sfd_moyal_pdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function moyal_cdf(a, b, x: double): double;
  {-Return the cumulative Moyal distribution function with}
  { location a and scale parameter b > 0}
begin
  moyal_cdf := sfd_moyal_cdf(a, b, x);
end;


{---------------------------------------------------------------------------}
function moyal_inv(a, b, y: double): double;
  {-Return the functional inverse of the Moyal distribution}
  { with location a and scale parameter b > 0}
begin
  moyal_inv := sfd_moyal_inv(a, b, y)
end;


{---------------------------------------------------------------------------}
function negbinom_cdf(p,r: double; k: longint): double;
  {-Return the cumulative negative binomial distribution function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}
begin
  negbinom_cdf := sfd_negbinom_cdf(p,r,k);
end;


{---------------------------------------------------------------------------}
function negbinom_pmf(p,r: double; k: longint): double;
  {-Return the negative binomial distribution probability mass function with target}
  { for number of successful trials r > 0 and success probability 0 <= p <= 1}
begin
  negbinom_pmf := sfd_negbinom_pmf(p,r,k);
end;


{---------------------------------------------------------------------------}
function normal_pdf(mu, sd, x: double): double;
  {-Return the normal (Gaussian) probability density function with mean mu}
  { and standard deviation sd>0, exp(-0.5*(x-mu)^2/sd^2) / sqrt(2*Pi*sd^2)}
begin
  normal_pdf := sfd_normal_pdf(mu,sd,x);
end;


{---------------------------------------------------------------------------}
function normal_cdf(mu, sd, x: double): double;
  {-Return the normal (Gaussian) distribution function}
  { with mean mu and standard deviation sd > 0        }
begin
  normal_cdf := sfd_normal_cdf(mu,sd,x);
end;


{---------------------------------------------------------------------------}
function normal_inv(mu, sd, y: double): double;
  {-Return the functional inverse of the normal (Gaussian) distribution}
  { with mean mu and standard deviation sd > 0, 0 < y < 1.}
begin
  normal_inv := sfd_normal_inv(mu, sd, y);
end;


{---------------------------------------------------------------------------}
function normstd_pdf(x: double): double;
  {-Return the std. normal probability density function exp(-x^2/2)/sqrt(2*Pi)}
begin
  if abs(x) >= 38.581 then normstd_pdf := 0.0
  else normstd_pdf := sfd_normstd_pdf(x);
end;


{---------------------------------------------------------------------------}
function normstd_cdf(x: double): double;
  {-Return the standard normal distribution function}
begin
  normstd_cdf := sfd_normstd_cdf(x);
end;


{---------------------------------------------------------------------------}
function normstd_inv(y: double): double;
  {-Return the inverse standard normal distribution function, 0 <= y <= 1.}
  { For x=normstd_inv(y) and y from [0,1], normstd_cdf(x) = y}
begin
  if y <= 0.0 then normstd_inv := -MaxDouble
  else if y >= 1.0 then normstd_inv := MaxDouble
  else normstd_inv := sfd_normstd_inv(y);
end;


{---------------------------------------------------------------------------}
function pareto_pdf(k, a, x: double): double;
  {-Return the Pareto probability density function with minimum value k > 0}
  { and shape a, x >= a > 0, result = (a/x)*(k/x)^a}
begin
  pareto_pdf := sfd_pareto_pdf(k,a,x);
end;


{---------------------------------------------------------------------------}
function pareto_cdf(k, a, x: double): double;
  {-Return the cumulative Pareto distribution function minimum value k > 0}
  { and shape a, x >= a > 0, result = 1-(k/x)^a}
begin
  pareto_cdf := sfd_pareto_cdf(k,a,x);
end;


{---------------------------------------------------------------------------}
function pareto_inv(k, a, y: double): double;
  {-Return the functional inverse of the Pareto distribution with minimum}
  { value k > 0 and shape a, x >= a > 0, result = k/(1-x)^(1/a)}
begin
  pareto_inv := sfd_pareto_inv(k,a,y);
end;


{---------------------------------------------------------------------------}
function poisson_cdf(mu: double; k: longint): double;
  {-Return the cumulative Poisson distribution function with mean mu >= 0}
begin
  poisson_cdf := sfd_poisson_cdf(mu,k);
end;


{---------------------------------------------------------------------------}
function poisson_pmf(mu: double; k: longint): double;
  {-Return the Poisson distribution probability mass function with mean mu >= 0}
begin
  poisson_pmf := sfd_poisson_pmf(mu,k);
end;


{---------------------------------------------------------------------------}
function rayleigh_pdf(b, x: double): double;
  {-Return the Rayleigh probability density function with}
  { scale b > 0, x >= 0; result = x*exp(-0.5*(x/b)^2)/b^2}
begin
  rayleigh_pdf := sfd_rayleigh_pdf(b,x);
end;


{---------------------------------------------------------------------------}
function rayleigh_cdf(b, x: double): double;
  {-Return the cumulative Rayleigh distribution function with scale b > 0, x >= 0}
begin
  rayleigh_cdf := sfd_rayleigh_cdf(b,x);
end;


{---------------------------------------------------------------------------}
function rayleigh_inv(b, y: double): double;
  {-Return the functional inverse of the Rayleigh distribution with scale b > 0}
begin
  rayleigh_inv := sfd_rayleigh_inv(b,y);
end;


{---------------------------------------------------------------------------}
function t_pdf(nu: longint; x: double): double;
  {-Return the probability density function of Student's t distribution, nu>0}
begin
  t_pdf := sfd_t_pdf(nu,x);
end;


{---------------------------------------------------------------------------}
function t_cdf(nu: longint; t: double): double;
  {-Return the cumulative Student t distribution with nu>0 degrees of freedom}
begin
  t_cdf := sfd_t_cdf(nu,t);
end;


{---------------------------------------------------------------------------}
function t_inv(nu: longint; p: double): double;
  {-Return the functional inverse of Student's t distribution, nu>0, 0 <= p <= 1}
begin
  t_inv := sfd_t_inv(nu,p);
end;


{---------------------------------------------------------------------------}
function triangular_pdf(a, b, c, x: double): double;
  {-Return the triangular probability density function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}
begin
  triangular_pdf := sfd_triangular_pdf(a,b,c,x);
end;


{---------------------------------------------------------------------------}
function triangular_cdf(a, b, c, x: double): double;
  {-Return the cumulative triangular distribution function with}
  { lower limit a, upper limit b, mode c;  a<b, a <= c <= b}
begin
  triangular_cdf := sfd_triangular_cdf(a,b,c,x);
end;


{---------------------------------------------------------------------------}
function triangular_inv(a, b, c, y: double): double;
  {-Return the functional inverse of the triangular distribution with}
  { lower limit a, upper limit b, mode c; a<b, a <= c <= b, 0 <= y <= 1}
begin
  triangular_inv := sfd_triangular_inv(a,b,c,y);
end;



{---------------------------------------------------------------------------}
function uniform_pdf(a, b, x: double): double;
  {-Return the uniform probability density function on [a,b], a<b}
begin
  uniform_pdf := sfd_uniform_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function uniform_cdf(a, b, x: double): double;
  {-Return the cumulative uniform distribution function on [a,b], a<b}
begin
  uniform_cdf := sfd_uniform_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function uniform_inv(a, b, y: double): double;
  {-Return the functional inverse of the uniform distribution on [a,b], a<b}
begin
  uniform_inv := sfd_uniform_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function wald_pdf(mu, b, x: double): double;
  {-Return the Wald (inverse Gaussian) probability density}
  { function with mean mu > 0, scale b > 0 for x >= 0     }
begin
  wald_pdf := sfd_wald_pdf(mu, b, x);
end;


{---------------------------------------------------------------------------}
function wald_cdf(mu, b, x: double): double;
  {-Return the Wald (inverse Gaussian) distribution  }
  { function with mean mu > 0, scale b > 0 for x >= 0}
begin
  wald_cdf := sfd_wald_cdf(mu, b, x);
end;


{---------------------------------------------------------------------------}
function wald_inv(mu, b, y: double): double;
  {-Return the functional inverse of the Wald (inverse Gaussian)}
  { distribution with mean mu > 0, scale b > 0, 0 <= y < 1.     }
begin
  wald_inv := sfd_wald_inv(mu, b, y);
end;


{---------------------------------------------------------------------------}
function weibull_pdf(a, b, x: double): double;
  {-Return the Weibull probability density function with shape a > 0}
  { and scale b > 0, result = a*x^(a-1)*exp(-(x/b)^a)/ b^a, x > 0}
begin
  weibull_pdf := sfd_weibull_pdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function weibull_cdf(a, b, x: double): double;
  {-Return the cumulative Weibull distribution function with}
  { shape parameter a > 0 and scale parameter b > 0}
begin
  weibull_cdf := sfd_weibull_cdf(a,b,x);
end;


{---------------------------------------------------------------------------}
function weibull_inv(a, b, y: double): double;
  {-Return the functional inverse of the Weibull distribution}
  { shape parameter a > 0 and scale parameter b > 0}
begin
  weibull_inv := sfd_weibull_inv(a,b,y);
end;


{---------------------------------------------------------------------------}
function zipf_pmf(r: double; k: longint): double;
  {-Return the Zipf distribution probability mass function k^(-(r+1))/zeta(r+1), r>0, k>0}
begin
  zipf_pmf := sfd_zipf_pmf(r,k);
end;


{---------------------------------------------------------------------------}
function zipf_cdf(r: double; k: longint): double;
  {-Return the cumulative Zipf distribution function H(k,r+1)/zeta(r+1), r>0, k>0}
begin
  zipf_cdf := sfd_zipf_cdf(r,k);
end;


{---------------------------------------------------------------------------}
{------------------ Orthogonal polynomials and related ---------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function chebyshev_t(n: integer; x: double): double;
  {-Return Tn(x), the Chebyshev polynomial of the first kind, degree n}
begin
  chebyshev_t := sfd_chebyshev_t(n,x);
end;


{---------------------------------------------------------------------------}
function chebyshev_u(n: integer; x: double): double;
  {-Return Un(x), the Chebyshev polynomial of the second kind, degree n}
begin
  chebyshev_u := sfd_chebyshev_u(n,x);
end;


{---------------------------------------------------------------------------}
function chebyshev_v(n: integer; x: double): double;
  {-Return V_n(x), the Chebyshev polynomial of the third kind, degree n>=0}
begin
  chebyshev_v := sfd_chebyshev_v(n,x);
end;


{---------------------------------------------------------------------------}
function chebyshev_w(n: integer; x: double): double;
  {-Return W_n(x), the Chebyshev polynomial of the fourth kind, degree n>=0}
begin
  chebyshev_w := sfd_chebyshev_w(n,x);
end;


{---------------------------------------------------------------------------}
function gegenbauer_c(n: integer; a,x: double): double;
  {-Return Cn(a,x), the nth Gegenbauer (ultraspherical) polynomial with}
  { parameter a. The degree n must be non-negative; a should be > -0.5 }
  { When a = 0,   C0(0,x) = 1,  and   Cn(0,x) = 2/n*Tn(x)   for n <> 0.}
begin
  gegenbauer_c := sfd_gegenbauer_c(n,a,x);
end;


{---------------------------------------------------------------------------}
function hermite_h(n: integer; x: double): double;
  {-Return Hn(x), the nth Hermite polynomial, degree n >= 0}
begin
  hermite_h := sfd_hermite_h(n,x);
end;


{---------------------------------------------------------------------------}
function jacobi_p(n: integer; a,b,x: double): double;
  {-Return Pn(a,b,x), the nth Jacobi polynomial with parameters a,b. Degree n}
  { must be >= 0; a,b should be > -1 (a+b must not be an integer < -1).}
begin
  jacobi_p := sfd_jacobi_p(n,a,b,x);
end;


{---------------------------------------------------------------------------}
function laguerre(n: integer; a,x: double): double;
  {-Return Ln(a,x), the nth generalized Laguerre polynomial with parameter a;}
  { degree n must be >= 0. x >=0 and a > -1 are the standard ranges.}
begin
  laguerre := sfd_laguerre(n,a,x);
end;


{---------------------------------------------------------------------------}
function laguerre_ass(n,m: integer; x: double): double;
  {-Return the associated Laguerre polynomial Ln(m,x); n,m >= 0}
begin
  laguerre_ass := sfd_laguerre_ass(n,m,x);
end;


{---------------------------------------------------------------------------}
function laguerre_l(n: integer; x: double): double;
  {-Return the nth Laguerre polynomial Ln(0,x); n >= 0}
begin
  laguerre_l := sfd_laguerre_l(n,x);
end;


{---------------------------------------------------------------------------}
function legendre_p(l: integer; x: double): double;
  {-Return P_l(x), the Legendre polynomial/function P_l, degree l}
begin
  legendre_p := sfd_legendre_p(l,x);
end;


{---------------------------------------------------------------------------}
function legendre_q(l: integer; x: double): double;
  {-Return Q_l(x), the Legendre function of the 2nd kind, degree l >=0, |x| <> 1}
begin
 legendre_q := sfd_legendre_q(l,x);
end;


{---------------------------------------------------------------------------}
function legendre_plm(l,m: integer; x: double): double;
  {-Return the associated Legendre polynomial P_lm(x)}
begin
  legendre_plm := sfd_legendre_plm(l,m,x);
end;


{---------------------------------------------------------------------------}
function legendre_qlm(l,m: integer; x: double): double;
  {-Return Q(l,m,x), the associated Legendre function of the second kind; l >= 0, l+m >= 0, |x|<>1}
begin
  legendre_qlm := sfd_legendre_qlm(l,m,x);
end;


{---------------------------------------------------------------------------}
procedure spherical_harmonic(l, m: integer; theta, phi: double; var yr,yi: double);
  {-Return Re and Im of the spherical harmonic function Y_lm(theta,phi)}
begin
  sfd_spherical_harmonic(l,m,theta,phi,yr,yi);
end;


{---------------------------------------------------------------------------}
function toroidal_plm(l,m: integer; x: double): double;
  {-Return the toroidal harmonic function P(l-0.5,m,x); l,m=0,1; x >= 1}
begin
  toroidal_plm := sfd_thp(l,m,x);
end;


{---------------------------------------------------------------------------}
function toroidal_qlm(l,m: integer; x: double): double;
  {-Return the toroidal harmonic function Q(l-0.5,m,x); l=0,1; x > 1}
begin
  toroidal_qlm := sfd_thq(l,m,x);
end;


{---------------------------------------------------------------------------}
function zernike_r(n,m: integer; r: double): double;
  {-Return the Zernike radial polynomial Rnm(r), r >= 0, n >= m >= 0, n-m even}
begin
  zernike_r := sfd_zernike_r(n,m,r);
end;


{---------------------------------------------------------------------------}
{-------------------- Hypergeometric functions -----------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function hyperg_1F1(a,b,x: double): double;
  {-Return the confluent hypergeometric function 1F1(a,b,x); Kummer's function M(a,b,x)}
begin
  hyperg_1F1 := sfd_1f1(a,b,x);
end;


{---------------------------------------------------------------------------}
function hyperg_1F1r(a,b,x: double): double;
  {-Return the regularized Kummer hypergeometric function 1F1(a,b,x)/Gamma(b)}
begin
  hyperg_1F1r := sfd_1f1r(a,b,x);
end;


{---------------------------------------------------------------------------}
function hyperg_u(a,b,x: double): double;
  {-Return Tricomi's confluent hypergeometric function U(a,b,x). If}
  { x<0, then a must be an integer and a<0 or 1+a-b an integer < 0.}
begin
  hyperg_u := sfd_chu(a,b,x);
end;


{---------------------------------------------------------------------------}
function hyperg_2F1(a,b,c,x: double): double;
  {-Return the Gauss hypergeometric function 2F1(a,b;c;x)}
begin
  hyperg_2F1 := sfd_2f1(a,b,c,x);
end;


{---------------------------------------------------------------------------}
function hyperg_2F1r(a,b,c,x: double): double;
  {-Return the regularized Gauss hypergeometric function 2F1(a,b,c,x)/Gamma(c)}
begin
  hyperg_2F1r := sfd_2f1r(a,b,c,x);
end;


{---------------------------------------------------------------------------}
function WhittakerM(k,m,x: double): double;
  {-Return the Whittaker M function = exp(-x/2)*x^(0.5+m) * 1F1(m-k-0.5,2m+1,x)}
begin
  WhittakerM := sfd_whitm(k,m,x);
end;


{---------------------------------------------------------------------------}
function WhittakerW(k,m,x: double): double;
  {-Return the Whittaker W function = exp(-x/2)*x^(0.5+m) * U(m-k-0.5,2m+1,x)}
begin
  WhittakerW := sfd_whitw(k,m,x);
end;


{---------------------------------------------------------------------------}
function hyperg_0F1(b,x: double): double;
  {-Return the confluent hypergeometric limit function 0F1(;b;x)}
begin
  hyperg_0F1 := sfd_0f1(b,x);
end;


{---------------------------------------------------------------------------}
function hyperg_0F1r(b,x: double): double;
  {-Return the regularized confluent hypergeometric limit function 0F1(;b;x)/Gamma(b)}
begin
  hyperg_0F1r := sfd_0f1r(b,x);
end;


{---------------------------------------------------------------------------}
function CylinderD(v,x: double): double;
  {-Return Whittaker's parabolic cylinder function D_v(x)}
begin
  CylinderD := sfd_pcfd(v,x);
end;


{---------------------------------------------------------------------------}
function CylinderU(a,x: double): double;
  {-Return the parabolic cylinder function U(a,x)}
begin
  CylinderU := sfd_pcfu(a,x);
end;


{---------------------------------------------------------------------------}
function CylinderV(a,x: double): double;
  {-Return the parabolic cylinder function V(a,x) with 2a integer}
begin
  CylinderV := sfd_pcfv(a,x);
end;


{---------------------------------------------------------------------------}
function HermiteH(v,x: double): double;
  {-Return the Hermite function H_v(x) of degree v}
begin
  HermiteH := sfd_pcfhh(v,x);
end;


{---------------------------------------------------------------------------}
{--------------------------- Other functions -------------------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function agm(x,y: double): double;
  {-Return the arithmetic-geometric mean of |x| and |y|}
begin
  agm := sfd_agm(x,y);
end;


{---------------------------------------------------------------------------}
function bernoulli(n: integer): double;
  {-Return the nth Bernoulli number, 0 if n<0 or odd n >= 3}
begin
  bernoulli := sfd_bernoulli(n);
end;


{---------------------------------------------------------------------------}
function bernpoly(n: integer; x: double): double;
  {-Return the Bernoulli polynomial B_n(x), 0 <= n <= MaxBernoulli}
begin
  bernpoly := sfd_bernpoly(n,x);
end;


{---------------------------------------------------------------------------}
function catalan(x: double): double;
  {-Return the Catalan function C(x) = binomial(2x,x)/(x+1)}
begin
  catalan := sfd_catf(x);
end;


{---------------------------------------------------------------------------}
function debye(n: integer; x: double): double;
  {-Return the Debye function D(n,x) = n/x^n*integral(t^n/(exp(t)-1),t=0..x) of order n>0, x>=0}
begin
  debye := sfd_debye(n,x);
end;


{---------------------------------------------------------------------------}
function euler(n: integer): double;
  {-Return the nth Euler number, 0 if n<0 or odd n}
begin
  euler := sfd_euler(n);
end;


{---------------------------------------------------------------------------}
function fibpoly(n: integer; x: double): double;
  {-Return the Fibonacci polynomial F_n(x)}
begin
  fibpoly := sfd_fpoly(n,x);
end;


{---------------------------------------------------------------------------}
function kepler(M,e: double): double;
  {-Solve Kepler's equation, result x is the eccentric anomaly from the mean anomaly M and the }
  { eccentricity e >= 0; x - e*sin(x) = M, x + x^3/3 = M, or e*sinh(x) - x = M for e <1, =1, >1}
begin
  kepler := sfd_kepler(M, e);
end;


{---------------------------------------------------------------------------}
function lucpoly(n: integer; x: double): double;
  {-Return the Lucas polynomial L_n(x)}
begin
  lucpoly := sfd_lpoly(n,x);
end;


{---------------------------------------------------------------------------}
function LambertW(x: double): double;
  {-Return the Lambert W function (principal branch), x >= -1/e}
begin
  LambertW := sfd_LambertW(x);
end;


{---------------------------------------------------------------------------}
function LambertW1(x: double): double;
  {-Return the Lambert W function (-1 branch), -1/e <= x < 0}
begin
  LambertW1 := sfd_LambertW1(x);
end;


{---------------------------------------------------------------------------}
function li_inv(x: double): double;
  {-Return the functional inverse of li(x), li(li_inv(x))=x}
begin
  li_inv := sfd_ali(x);
end;


{---------------------------------------------------------------------------}
function RiemannR(x: double): double;
  {-Return the Riemann prime counting function R(x), x >= 1/16}
begin
  RiemannR := sfd_ri(x);
end;


{---------------------------------------------------------------------------}
function cosint(n: integer; x: double): double;
  {-Return cosint(n, x) = integral(cos(t)^n, t=0..x), n>=0}
begin
  cosint := sfd_cosint(n,x);
end;


{---------------------------------------------------------------------------}
function sinint(n: integer; x: double): double;
  {-Return sinint(n, x) = integral(sin(t)^n, t=0..x), n>=0}
begin
  sinint := sfd_sinint(n,x);
end;


end.
