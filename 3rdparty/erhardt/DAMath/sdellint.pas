unit sdEllInt;

{Double precision Elliptic integrals, Jacobi elliptic and theta functions}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Double precision Elliptic integrals, Jacobi elliptic and theta functions

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
                 [10] R. Bulirsch, Numerical Calculation of Elliptic Integrals
                      and Elliptic Functions
                      Numerische Mathematik 7, 78-90, 1965
                 [11] R. Bulirsch, Numerical Calculation of Elliptic Integrals
                      and Elliptic Functions, part III.
                      Numerische Mathematik 13, 305-315, 1969
                 [12] B.C. Carlson, Computing Elliptic Integrals by Duplication
                      Numerische Mathematik 33, 1-16, 1979
                 [13] W.H. Press et al, Numerical Recipes in C, 2nd ed., Cambridge, 1992
                 [19] Boost C++ Libraries, Release 1.42.0, 2010.
                      http://www.boost.org/
                 [30] [NIST]: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST Handbook
                      of Mathematical Functions, Cambridge, 2010. Online resource: NIST Digital
                      Library of Mathematical Functions, http://dlmf.nist.gov/
                 [33] http://functions.wolfram.com/: Formulas and graphics about
                      mathematical functions for the mathematical and scientific
                      community and/or http://mathworld.wolfram.com/ ("/the web's
                      most extensive mathematical resource/")
                 [67] B.C. Carlson, Numerical computation of real or complex elliptic integrals,
                      1994, https://arxiv.org/pdf/math/9409227v1.pdf


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  05.02.13  W.Ehrhardt  Initial BP7 version from AMath.sfellint
 1.00.01  06.02.13  we          THexDblW constants for sfd_nome
 1.00.02  06.02.13  we          adjusted some constants for double precision
 1.00.03  06.02.13  we          changed test 0<x to x<-eps_d  in sfd_ellint_?
 1.00.04  06.02.13  we          changed test abs(t-Pi_2)
 1.00.05  07.02.13  we          adjusted Carlson eps constants
 1.00.06  07.02.13  we          sfd_EllipticPiC/CPi always with cel

 1.01.00  15.03.13  we          Remaining Jacobi elliptic functions
 1.01.01  25.03.13  we          Remaining inverse Jacobi elliptic functions

 1.06.00  28.09.13  we          special_reduce_modpi, used in sfd_ellint_1/2/3 and sfd_hlambda

 1.07.00  23.10.13  we          sfd_EllipticKim, sfd_EllipticECim

 1.17.00  28.03.15  we          sfd_cel_d, sfd_ellint_d

 1.18.00  22.06.15  we          avoid overflow in sfd_ell_rf, better arg checks in Carlson functions
 1.18.01  22.06.15  we          new sfd_ell_rg

 1.19.00  19.07.15  we          sfd_el1(x, kc) = arctan(x) for kc^2=1

 1.20.00  25.08.15  we          Use THREE & SIXX to avoid 'optimization'

 1.21.00  26.09.15  we          lemniscate functions sfd_lemn

 1.23.00  12.05.16  we          sfd_cel_b, sfd_ellint_b

 1.25.00  21.05.17  we          sfd_EllipticK, sfd_EllipticEC for k>1

 1.31.00  30.03.18  we          sfd_cel_d and sfd_cel_b for |k| >= 1
 1.31.01  31.03.18  we          sfd_EllipticPiC for |k| > 1
 1.31.02  02.04.18  we          sfd_EllipticPiCim
 1.31.03  05.04.18  we          Fix sfd_EllipticPiC for k^2 >= nu > 1
 1.31.04  05.04.18  we          Mathematica style complete elliptic integrals
 1.31.05  12.04.18  we          Neville theta functions
 1.31.06  15.04.18  we          sfd_arccl/sl

 1.32.00  19.04.18  we          Fixes for FPC311
 1.32.01  22.04.18  we          sfd_el2 computed with Carlson
 1.32.02  24.04.18  we          sfd_emlambda(i*y): elliptic modular function
 1.32.03  24.04.18  we          sfd_detai from sdmisc
 1.32.04  28.04.18  we          Weierstrass elliptic functions
 1.32.05  01.05.18  we          sfd_wpg with JacobiCN, CompAuxWP, equireduce, wpdup
 1.32.06  03.05.18  we          sfd_wpg_inv
 1.32.07  04.05.18  we          changed argument check in sfd_ell_rd
 1.32.08  04.05.18  we          sfd_M_EllipticF/sfd_M_EllipticE
 1.32.09  06.05.18  we          sfd_KleinJ

 1.33.00  14.05.18  we          Basic lemniscatic case: separate function for argument
 1.33.01  14.05.18  we          sfd_wpe_der / sfd_wpg_der
 1.33.02  23.05.18  we          sfd_M_EllipticF := phi for m=0
 1.33.03  23.05.18  we          sfd_M_EllipticPi

 1.35.00  20.08.18  we          IsNanOrInf check sfd_EllipticPiCim
 1.35.01  23.08.18  we          Implemented max iteration check in some functions

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


{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Legendre style) -------------------}
{---------------------------------------------------------------------------}
{#Z-}
function sfd_ellint_1(phi,k: double): double;
  {-Return the Legendre elliptic integral F(phi,k) of the 1st kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}

function sfd_ellint_2(phi,k: double): double;
  {-Return the Legendre elliptic integral E(phi,k) of the 2nd kind}
  { = integral(sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}

function sfd_ellint_3(phi,nu,k: double): double;
  {-Return the Legendre elliptic integral PI(phi,nu,k) of the 3rd kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2)/(1-nu*sin(x)^2),x=0..phi) with  }
  { |k*sin(phi)|<=1, returns Cauchy principal value if nu*sin(phi)^2>1}

function sfd_ellint_b(phi,k: double): double;
  {-Return the Legendre elliptic integral B(phi,k) = (E(phi,k) - kc^2*F(phi,k))/k^2}
  { = integral(cos(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }

function sfd_ellint_d(phi,k: double): double;
  {-Return the Legendre elliptic integral D(phi,k) = (F(phi,k) - E(phi,k))/k^2 }
  { = integral(sin(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }

function sfd_cel_b(k: double): double;
  {-Return the complete elliptic integral B(k) = (E(k) - kc^2*K(k))/k^2, real part for |k|>1}

function sfd_cel_d(k: double): double;
  {-Return the complete elliptic integral D(k) = (K(k) - E(k))/k^2, real part for |k|>1}

function sfd_hlambda(phi,k: double): double;
  {-Return Heuman's function Lambda_0(phi,k) = F(phi,k')/K(k') + 2/Pi*K(k)*Z(phi,k'), |k|<=1}

function sfd_jzeta(phi,k: double): double;
  {-Return the Jacobi Zeta function Z(phi,k) = E(phi,k) - E(k)/K(k)*F(phi,k), |k|<=1}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Carlson style) --------------------}
{---------------------------------------------------------------------------}
{#Z-}
function sfd_ell_rc(x,y: double): double;
  {-Return Carlson's degenerate elliptic integral RC; x>=0, y<>0}

function sfd_ell_rf(x,y,z: double): double;
  {-Return Carlson's elliptic integral of the 1st kind; x,y,z >=0, at most one =0}

function sfd_ell_rd(x,y,z: double): double;
  {-Return Carlson's elliptic integral of the 2nd kind; z>0; x,y >=0, at most one =0}

function sfd_ell_rg(x,y,z: double): double;
  {-Return Carlson's completely symmetric elliptic integral of the 2nd kind; x,y,z >= 0}

function sfd_ell_rj(x,y,z,r: double): double;
  {-Return Carlson's elliptic integral of the 3rd kind; r<>0; x,y,z >=0, at most one =0}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Bulirsch style) -------------------}
{---------------------------------------------------------------------------}
{#Z-}
function sfd_cel1(kc: double): double;
  {-Return Bulirsch's complete elliptic integral of the 1st kind, kc<>0}

function sfd_cel2(kc, a, b: double): double;
  {-Return Bulirsch's complete elliptic integral of the 2nd kind, kc<>0}

function sfd_cel(kc, p, a, b: double): double;
  {-Return Bulirsch's general complete elliptic integral, kc<>0, Cauchy principle value if p<0}

function sfd_el1(x,kc: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 1st kind}

function sfd_el2(x,kc,a,b: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 2nd kind, kc<>0}

function sfd_el3(x,kc,p: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 3rd kind, 1+p*x^2<>0}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Maple V style) --------------------}
{---------------------------------------------------------------------------}
{#Z-}
function sfd_EllipticF(z,k: double): double;
  {-Return the incomplete elliptic integral of the 1st kind; |z|<=1, |k*z|<1}

function sfd_EllipticK(k: double): double;
  {-Return the complete elliptic integral of the 1st kind}
  { = integral(1/sqrt(1-x^2)/sqrt(1-k^2*x^2),x=0..1),    }
  { calculated with agm, real part if |k| > 1            }

function sfd_EllipticKim(k: double): double;
  {-Return K(i*k), the complete elliptic integral of the 1st kind with}
  { imaginary modulus = integral(1/sqrt(1-x^2)/sqrt(1+k^2*x^2),x=0..1)}

function sfd_EllipticCK(k: double): double;
  {-Return the complementary complete elliptic integral of the 1st kind, k<>0}

function sfd_EllipticE(z,k: double): double;
  {-Return the incomplete elliptic integrals of the 2nd kind, |z|<=1, |k*z| <= 1}

function sfd_EllipticEC(k: double): double;
  {-Return the complete elliptic integral of the 2nd kind,}
  { = integral(sqrt(1-k^2*x^2)/sqrt(1-x^2),x=0..1), real part for k>1}

function sfd_EllipticECim(k: double): double;
  {-Return E(i*k), the complete elliptic integral of the 2nd kind with}
  { imaginary modulus = integral(sqrt(1+k^2*x^2)/sqrt(1-x^2),x=0..1)  }

function sfd_EllipticCE(k: double): double;
  {-Return the complementary complete elliptic integral of the 2nd kind}

function sfd_EllipticPi(z,nu,k: double): double;
  {-Return the incomplete elliptic integral of the 3rd kind, |z|<=1, |k*z|<1}

function sfd_EllipticPiC(nu,k: double): double;
  {-Return the complete elliptic integral of the 3rd kind, |k|<>1, nu<>1; real part for |k|>1}

function sfd_EllipticPiCim(nu,k: double): double;
  {-Return Pi(nu, k*i), the complete elliptic integral of the 3rd kind}
  { with imaginary modulus, nu <> 1, real part if nu > 1}

function sfd_EllipticCPi(nu,k: double): double;
  {-Return the complementary complete elliptic integral of the 3rd kind, k<>0, nu<>1}

{#Z+}
{---------------------------------------------------------------------------}
{---------------- Elliptic integrals (Mathematica style) -------------------}
{---------------------------------------------------------------------------}
{#Z-}
function sfd_M_EllipticK(m: double): double;
  {-Return the complete elliptic integral of the 1st kind,}
  { K(m) = integral(dx/sqrt(1-m*sin(x)^2),x=0..Pi/2), real part for m>1}

function sfd_M_EllipticEC(m: double): double;
  {-Return the complete elliptic integral of the 2nd kind,}
  { E(m) = integral(sqrt(1-m*sin(x)^2),x=0..Pi/2), real part for m>1}

function sfd_M_EllipticPiC(n,m: double): double;
  {-Return the complete elliptic integral of the 3rd kind, n<>1, m<>1, real part}
  { for m>1; Pi(n|m) = integral(1/(1-n*sin(x)^2/sqrt(1-m*sin(x)^2)), x=0..Pi/2)}

function sfd_M_EllipticF(phi,m: double): double;
  {-Return the incomplete elliptic integral of the 1st kind,}
  { F(phi,m) = integral(dx/sqrt(1-m*sin(x)^2),x=0..phi), m*sin(phi)^2 <= 1}

function sfd_M_EllipticE(phi,m: double): double;
  {-Return the incomplete elliptic integral of the 2nd kind,}
  { E(phi,m) = integral(sqrt(1-m*sin(x)^2),x=0..phi), m*sin(phi)^2 <= 1}

function sfd_M_EllipticPi(n,phi,m: double): double;
  {-Return the incomplete elliptic integral Pi(n,phi,m) of the 3rd kind}
  { = integral(1/sqrt(1-m*sin(x)^2)/(1-n*sin(x)^2),x=0..phi), with n<>1}
  { m*sin(phi)^2 <= 1, real part if complex.}

{#Z+}
{--------------------------------------------------------------------------}
{------------------- Jacobi elliptic and theta functions ------------------}
{--------------------------------------------------------------------------}
{#Z-}
function sfd_ellmod(q: double): double;
  {-Return the elliptic modulus k(q) = theta_2(q)^2/theta_3(q)^2, 0 <= q <= 1}

function sfd_nome(k: double): double;
  {-Return the elliptic nome q(k) = exp(-Pi*EllipticCK(k)/EllipticK(k)), |k| < 1}

function sfd_jam(x,k: double): double;
  {-Return the Jacobi amplitude am(x,k)}

procedure sfd_sncndn(x,mc: double; var sn,cn,dn: double);
  {-Return the Jacobi elliptic functions sn,cn,dn for argument x and}
  { complementary parameter mc.}

function sfd_jacobi_sn(x,k: double): double;
  {-Return the Jacobi elliptic function sn(x,k)}

function sfd_jacobi_cn(x,k: double): double;
  {-Return the Jacobi elliptic function cn(x,k)}

function sfd_jacobi_dn(x,k: double): double;
  {-Return the Jacobi elliptic function dn(x,k)}

function sfd_jacobi_nc(x,k: double): double;
  {-Return the Jacobi elliptic function nc(x,k)}

function sfd_jacobi_sc(x,k: double): double;
  {-Return the Jacobi elliptic function sc(x,k)}

function sfd_jacobi_dc(x,k: double): double;
  {-Return the Jacobi elliptic function dc(x,k)}

function sfd_jacobi_nd(x,k: double): double;
  {-Return the Jacobi elliptic function nd(x,k)}

function sfd_jacobi_sd(x,k: double): double;
  {-Return the Jacobi elliptic function sd(x,k)}

function sfd_jacobi_cd(x,k: double): double;
  {-Return the Jacobi elliptic function cd(x,k)}

function sfd_jacobi_ns(x,k: double): double;
  {-Return the Jacobi elliptic function ns(x,k)}

function sfd_jacobi_cs(x,k: double): double;
  {-Return the Jacobi elliptic function cs(x,k)}

function sfd_jacobi_ds(x,k: double): double;
  {-Return the Jacobi elliptic function ds(x,k)}

function sfd_jacobi_arccn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccn(x,k), |x| <= 1, x >= sqrt(1 - 1/k^2) if k >= 1}

function sfd_jacobi_arccd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccd(x,k); |x| <= 1 if |k| < 1; |x| >= 1 if |k| > 1 }

function sfd_jacobi_arccs(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccs(x,k), |x| >= sqrt(k^2-1) for |k|>1}

function sfd_jacobi_arcdc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcdc(x,k); |x| >= 1 if |k| < 1; |x| <= 1 if |k| > 1 }

function sfd_jacobi_arcdn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcdn(x,k), 0 <= x <= 1, x^2 + k^2 > 1 if |k| < 1;  |x| <= 1 if |k| > 1}

function sfd_jacobi_arcds(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcds(x,k), x^2 + k^2 >= 1}

function sfd_jacobi_arcnc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcnc(x,k), x >= 1, x^2 <= k^2/(k^2-1) for |k|>1}

function sfd_jacobi_arcnd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcnd(x,k), x >= 1, x^2 <= k^2/(1-k^2) if k < 1}

function sfd_jacobi_arcns(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcns(x,k), |x| >= 1, |x| >= k if k>=1}

function sfd_jacobi_arcsc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsc(x,k), |x| <= 1/sqrt(k^2-1) for |k|>1}

function sfd_jacobi_arcsd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsd(x,k), x^2*(1-k^2) <= 1}

function sfd_jacobi_arcsn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsn(x,k), |x| <= 1 and |x*k| <= 1}


function sfd_theta2(q: double): double;
  {-Return Jacobi theta_2(q) = 2*q^(1/4)*sum(q^(n*(n+1)),n=0..Inf) 0 <= q < 1}

function sfd_theta3(q: double): double;
  {-Return Jacobi theta_3(q) = 1 + 2*sum(q^(n*n)),n=1..Inf); |q| < 1}

function sfd_theta4(q: double): double;
  {-Return Jacobi theta_4(q) = 1 + 2*sum((-1)^n*q^(n*(n+1)),n=1..Inf); |q| < 1}

function sfd_theta1p(q: double): double;
  {-Return the derivative  theta1p(q) := d/dx(theta_1(x,q)) at x=0,}
  { = 2*q^(1/4)*sum((-1)^n*(2n+1)*q^(n*(n+1)),n=0..Inf), 0 <= q < 1}

function sfd_jtheta(n: integer; x,q: double): double;
  {-Return Jacobi theta_n(x,q), n=1..4, 0 <= q < 1}


function sfd_ntheta_s(x, k: double): double;
  {-Return the Neville theta_s function, |k| <= 1}

function sfd_ntheta_c(x, k: double): double;
  {-Return the Neville theta_c function, |k| <= 1}

function sfd_ntheta_d(x, k: double): double;
  {-Return the Neville theta_d function, |k| <= 1}

function sfd_ntheta_n(x, k: double): double;
  {-Return the Neville theta_n function, |k| <= 1}


procedure sfd_lemn(x: double; var sl,cl: double);
  {-Return the lemniscate functions sl = sin_lemn(x), cl = cos_lemn(x)}

function sfd_sin_lemn(x: double): double;
  {-Return the lemniscate sine function sin_lemn(x)}

function sfd_cos_lemn(x: double): double;
  {-Return the lemniscate cosine function cos_lemn(x)}

function sfd_arccl(x: double): double;
  {-Return the inverse lemniscate cosine function, |x| <= 1}

function sfd_arcsl(x: double): double;
  {-Return the inverse lemniscate sine function, |x| <= 1}


{#Z+}
{---------------------------------------------------------------------------}
{--------------- Weierstrass elliptic and modular functions ----------------}
{---------------------------------------------------------------------------}
{#Z-}
function sfd_detai(x: double): double;
  {-Return Dedekind eta(i*x), x >= 0}

function sfd_emlambda(y: double): double;
  {-Return the elliptic modular function lambda(iy), y >= 0}

function sfd_KleinJ(y: double): double;
  {-Return Klein's complete invariant J(iy), y>0}

function sfd_wpl(x: double): double;
  {-Return the Weierstrass function wp(x,1,0)=wpe(x,1/2,0), basic lemniscatic case}

function sfd_wpe(x,e1,e2: double): double;
  {-Return the Weierstrass function P(x,e1,e2) from the lattice roots e1 < e2}

function sfd_wpe_der(x,e1,e2: double): double;
  {-Return Weierstrass P'(x,e1,e2) from the lattice roots e1 < e2}

function sfd_wpe_im(y,e1,e2: double): double;
  {-Return the Weierstrass function P(iy,e1,e2) from the lattice roots e1 < e2}

function sfd_wpe_inv(y,e1,e2: double): double;
  {-Return the smallest positive x with wpe(x)=y, y >= e1}

function sfd_wpg(x,g2,g3: double): double;
  {-Return the Weierstrass function P(x,e1,e2) from lattice invariants g2, g3}

function sfd_wpg_der(x,g2,g3: double): double;
  {-Return Weierstrass P'(x,e1,e2) from lattice invariants g2, g3}

function sfd_wpg_im(y,g2,g3: double): double;
  {-Return the Weierstrass function P(iy, g2, g3)}

function sfd_wpg_inv(y,g2,g3: double): double;
  {-Return the smallest positive x with wpg(x,g2,g3)=y, y >= e2}


implementation


uses
  DAMath,
  sdBasic,  {for RTE_ArgumentRange}
  sdMisc;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Carlson style) --------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function sfd_ell_rc(x,y: double): double;
  {-Return Carlson's degenerate elliptic integral RC; x>=0, y<>0}
var
  l,m,s: double;
const
  eps = 0.0006; {if |s| < eps then |error| < 16|s|^6/(1-2|s|) < 0.75e-18}
  c1  = 3/10;
  c2  = 1/7;    {Fix311}
  c3  = 3/8;
  c4  = 9/22;
const
  IMAX = 50;
var
  i: integer;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(y) or (x<0.0) or (y=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ell_rc := NaN_d;
    exit;
  end;
  if y<0.0 then begin
    {Return Cauchy principal value if y < 0, see [12], (2.12)}
    if x=0 then sfd_ell_rc := 0
    else begin
      y := -y;
      sfd_ell_rc := sqrt(x/(x+y))*sfd_ell_rc(x+y,y);
    end;
  end
  else begin
    {This is Carlson's [12] Algorithm 2, see also [13], p.267 function rc}
    i := 1;
    repeat
      l := y + 2.0*sqrt(x*y);
      x := 0.25*(x+l);
      y := 0.25*(y+l);
      s := x + 2.0*y;
      m := s/THREE;
      s := (y-x)/s;
      inc(i);
    until (abs(s)<eps) or (i>IMAX);
    {Test no convergence}
    if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
    s := (((c4*s + c3)*s + c2)*s + c1)*s*s;
    sfd_ell_rc := (1.0+s)/sqrt(m);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ell_rf(x,y,z: double): double;
  {-Return Carlson's elliptic integral of the 1st kind; x,y,z >=0, at most one =0}
const
  eps = 0.0001;   {[12]: If s < eps then |error| < s^4.25/(1-s)/4 < 0.25e-17}
                  {Note: Carlson's original s^6 scaling does not work for   }
                  {double precision, see e.g. [19] (Boost uses s^4.25 too)  }
  c1 :double = 1/24;  {Fix311}
  c2 :double = 1/10;
  c3 :double = 3/44;
  c4 :double = 1/14;
var
  l,m,dx,dy,dz,e2,e3,rx,ry,rz,s: double;
const
  IMAX = 50;
var
  i: integer;
begin

  if IsNanOrInfD(x) or IsNanOrInfD(y) or IsNanOrInfD(z) or (x<0.0) or (y<0.0) or (z<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ell_rf := NaN_d;
    exit;
  end;

  {Force at most one of x,y,z=0}
  if isign(x)+isign(y)+isign(z) < 2 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ell_rf := NaN_d;
    exit;
  end;

  {This is Carlson's [12] Algorithm 1, see also [13], p.264 function rf}
  i := 1;
  repeat
    rx := sqrt(x);
    ry := sqrt(y);
    rz := sqrt(z);
    l  := rx*(ry+rz) + ry*rz;
    x  := 0.25*(l+x);
    y  := 0.25*(l+y);
    z  := 0.25*(l+z);
    m  := (x+y+z)/THREE;
    dx := (m-x)/m;
    dy := (m-y)/m;
    dz := (m-z)/m;
    inc(i);
    if i>IMAX then break;
  until (abs(dx)<eps) and (abs(dy)<eps) and (abs(dz)<eps);
  {Test no convergence}
  if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  e2 := dx*dy - dz*dz;
  e3 := dx*dy*dz;
  s  := e2*(c1*e2 - c2 - c3*e3) + c4*e3;
  sfd_ell_rf := (1.0+s)/sqrt(m);
end;


{---------------------------------------------------------------------------}
function sfd_ell_rd(x,y,z: double): double;
  {-Return Carlson's elliptic integral of the 2nd kind; z>0; x,y >=0, at most one =0}
const
  eps = 0.001;   {if s < eps then |error| < s^6/(1-s)^1.5) < 0.12e-16}
  c1 :double = -3/14;  {Fix311}
  c2 :double =  1/6;
  c3 :double =  9/22;
  c4 :double =  3/26;
  c5 :double =  9/88;
  c6 :double =  9/52;
var
  l,m,dx,dy,dz,ea,eb,ec,ed,ef,p4,rx,ry,rz,s,sn: double;
const
  IMAX = 50;
var
  i: integer;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(y) or IsNanOrInfD(z) or (x<0.0) or (y<0.0) or (z<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ell_rd := NaN_d;
    exit;
  end;

  {This is Carlson's [12] Algorithm 4, see also [13], p.265 function rd}
  sn := 0.0;
  p4 := 1.0;
  i := 1;
  repeat
    rx := sqrt(x);
    ry := sqrt(y);
    rz := sqrt(z);
    l  := rx*(ry+rz) + ry*rz;
    sn := sn + p4/(rz*(z+l));
    p4 := 0.25*p4;
    x  := 0.25*(l+x);
    y  := 0.25*(l+y);
    z  := 0.25*(l+z);
    m  := (x + y + 3.0*z)/5.0;
    dx := (m-x)/m;
    dy := (m-y)/m;
    dz := (m-z)/m;
    inc(i);
    if i>IMAX then break;
  until (abs(dx)<eps) and (abs(dy)<eps) and (abs(dz)<eps);
  {Test no convergence}
  if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  ea := dx*dy;
  eb := dz*dz;
  ec := ea - eb;
  ed := ea - 6.0*eb;
  ef := ed + ec + ec;
  s  := ed*(c1 + c5*ed - c6*dz*ef);
  s  := s + dz*(c2*ef + dz*(dz*c4*ea - c3*ec));
  sfd_ell_rd := 3.0*sn + p4*(1.0 + s)/(m*sqrt(m));
end;


{---------------------------------------------------------------------------}
function sfd_ell_rg(x,y,z: double): double;
  {-Return Carlson's completely symmetric elliptic integral of the 2nd kind; x,y,z >= 0}
var
  f,d,s: double;
const
  eps = 1e-22; {eps*ln(eps) small}
begin
  if IsNanOrInfD(x) or IsNanOrInfD(y) or IsNanOrInfD(z) or (x<0.0) or (y<0.0) or (z<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ell_rg := NaN_d;
    exit;
  end;

  {Order arguments, integral does not change because of symmetry. This}
  {is done to make (x-z)*(y-z) <= 0 in the general case, which avoids }
  {cancellations, see the comment after NIST[30], 19.21.10.           }
  if x<y then begin f:=x; x:=y; y:=f; end;
  if x<z then begin f:=x; x:=z; z:=f; end;
  if y>z then begin f:=y; y:=z; z:=f; end;

  {here x >= z >= y}
  if x=y then begin
    {x=z=y}
    sfd_ell_rg := sqrt(x);          {NIST[30] 19.20.4,1}
  end
  else if z=0.0 then begin
    {x>0, z=y=0}
    sfd_ell_rg := 0.5*sqrt(x);      {NIST[30] 19.20.4,4}
  end
  else if (y=0.0) then begin
    if x=z then begin
      {x=z>0, y=0}
      sfd_ell_rg := Pi_4*sqrt(x);   {NIST[30] 19.20.4,3}
    end
    else begin
      {x>z>0, y=0}
      {From Carlson[67],(56) and E(k) expressed with cel2}
      s := z/x;
      f := 0.5*sfd_cel2(sqrt(s),1.0,s);
      sfd_ell_rg := sqrt(x)*f;
    end;
  end
  else if z <= eps*x then begin
    {Asymptotic case:  NIST[30], 19.27.4}
    sfd_ell_rg := 0.5*sqrt(x);
  end
  else begin
    {general case NIST[30], 19.21.10: here x >=z >= y > 0}
    s := sqrt((x/z)*y);
    f := z*sfd_ell_rf(x,y,z);
    if (z=x) or (z=y) then d := 0.0
    else begin
      {split products to avoid overflow}
      d := sfd_ell_rd(x,y,z)/THREE;
      d := (x-z)*d;
      d := (z-y)*d;
    end;
    sfd_ell_rg := 0.5*(f+d+s);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ell_rj(x,y,z,r: double): double;
  {-Return Carlson's elliptic integral of the 3rd kind; r<>0; x,y,z >=0, at most one =0}
const
  eps = 0.001;   {if s < eps then |error| < s^6/(1-s)^1.5) < 0.12e-16}
  c1  = -3/14;   {Fix311}
  c2  = 1/3;
  c3  = 3/22;
  c4  = 3/26;
  c5  = 0.75*c3;
  c6  = 1.5*c4;
  c7  = 0.5*c2;
  c8  = -2.0*c3;
var
  l,m,dr,dx,dy,dz,ea,eb,ec,ed,ee,p4,rx,ry,rz,s,sn: double;
const
  IMAX = 50;
var
  i: integer;
begin

  if IsNanOrInfD(x) or IsNanOrInfD(y) or IsNanOrInfD(z) or IsNanOrInfD(r)
     or (x<0.0) or (y<0.0) or (z<0.0) or (r=0.0) then
  begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ell_rj := NaN_d;
    exit;
  end;

  {Force at most one of x,y,z=0}
  if isign(x)+isign(y)+isign(z) < 2 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ell_rj := NaN_d;
    exit;
  end;

  if r<0.0 then begin
    {Return Cauchy principal value if r < 0, see [12], (2.22)}
    {Order 0 <= x <= y <= z}
    if x>y then begin s:=x; x:=y; y:=s; end;
    if x>z then begin s:=x; x:=z; z:=s; end;
    if y>z then begin s:=y; y:=z; z:=s; end;
    {m = gamma from 2.22}
    m := y + (z-y)*(y-x)/(y-r);
    s := sfd_ell_rc(x*z/y, r*m/y) - sfd_ell_rf(x,y,z);
    s := (m-y)*sfd_ell_rj(x,y,z,m) + 3.0*s;
    sfd_ell_rj := s/(y-r);
    exit;
  end;

  {This is Carlson's [12] Algorithm 3, see also [13], p.266 function rj}
  sn := 0.0;
  p4 := 1.0;
  i := 1;
  repeat
    rx := sqrt(x);
    ry := sqrt(y);
    rz := sqrt(z);
    l  := rx*(ry+rz) + ry*rz;
    ea := sqr(r*(rx+ry+rz) + rx*ry*rz); { = alpha}
    eb := r*sqr(r+l);                   { = beta }
    sn := sn + p4*sfd_ell_rc(ea,eb);
    p4 := 0.25*p4;
    x  := 0.25*(l+x);
    y  := 0.25*(l+y);
    z  := 0.25*(l+z);
    r  := 0.25*(l+r);
    m  := (x+y+z+r+r)/5.0;
    dx := (m-x)/m;
    dy := (m-y)/m;
    dz := (m-z)/m;
    dr := (m-r)/m;
    inc(i);
    if i>IMAX then break;
  until (abs(dx)<eps) and (abs(dy)<eps) and (abs(dz)<eps) and (abs(dr)<eps);
  {Test no convergence}
  if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));

  ea := dx*(dy+dz) + dy*dz;
  eb := dx*dy*dz;
  ec := dr*dr;
  ed := ea - 3.0*ec;
  ee := eb + 2.0*dr*(ea-ec);
  s  := ed*(c1 + c5*ed - c6*ee);
  s  := s + eb*(c7 + dr * (c8 + dr*c4));
  s  := s + dr*ea*(c2 - dr*c3) - c2*dr*ec;
  sfd_ell_rj := 3.0*sn + p4*(1.0 + s)/(m*sqrt(m));
end;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Bulirsch style) -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function sfd_cel1(kc: double): double;
  {-Return Bulirsch's complete elliptic integral of the 1st kind, kc<>0}
  { = integral(1/sqrt((1+t^2)*(1+kc^2*t^2)), t=0..Inf)}
var
  h,m: double;
const
  IMAX = 50;
var
  i: integer;
begin
  if IsNanOrInfD(kc) or (kc=0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_cel1 := NaN_d;
    exit;
  end;
  kc := abs(kc);

  {This is my Pascal translation of the Algol procedure cel1 from [10]}
  m := 1.0;
  i := 1;
  repeat
    h := m;
    m := kc + m;
    if abs(h-kc) <= 0.5e-8*h then break;
    kc:= sqrt(h*kc);
    m := 0.5*m;
    inc(i);
  until i>IMAX;
  {Test no convergence}
  if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  sfd_cel1 := Pi/m;
end;


{---------------------------------------------------------------------------}
function sfd_cel2(kc, a, b: double): double;
  {-Return Bulirsch's complete elliptic integral of the 2nd kind, kc<>0}
  { = integral((a+b*t^2)/[(1+t^2)*sqrt((1+t^2)*(1+kc^2*t^2))], t=0..Inf)}
var
  c,m,h: double;
const
  IMAX = 50;
var
  i: integer;
begin
  {Because the following algorithm is accurate to double precision}
  {only if a*b >= 0, sfd_cel(kc,1,a,b) is called if a*b<0}

  if a*b<0 then begin
    sfd_cel2 := sfd_cel(kc, 1.0, a, b);
    exit;
  end;

  if IsNanOrInfD(kc) or (kc=0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_cel2 := NaN_d;
    exit;
  end;
  kc := abs(kc);

  {This is my Pascal translation of the Algol procedure cel2 from [10]}
  m := 1.0;
  i := 1;
  c := a;
  a := a+b;
  repeat
    b := 2.0*(c*kc + b);
    c := a;
    h := m;
    m := kc+m;
    a := b/m + a;
    if abs(h-kc) <= 0.5e-8*h then break;
    kc:= 2*sqrt(h*kc);
    inc(i);
  until i>IMAX;
  {Test no convergence}
  if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  sfd_cel2 := Pi_4*a/m;
end;


{---------------------------------------------------------------------------}
function sfd_cel(kc, p, a, b: double): double;
  {-Return Bulirsch's general complete elliptic integral, kc<>0, Cauchy principle value if p<0}
  { = integral((a+b*t^2)/[(1+p*t^2)*sqrt((1+t^2)*(1+kc^2*t^2))], t=0..Inf)}
var
  e,f,g,m: double;
const
  IMAX = 50;
var
  i: integer;
begin

  if IsNanOrInfD(kc) or ((kc=0.0) and (b<>0.0)) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_cel := NaN_d;
    exit;
  end;
  kc := abs(kc);
  if (kc=0.0) and (b=0.0) then kc := 1e-17;

  {This is my Pascal translation of the Algol procedure cel from [11]}
  e  := kc;
  if p>0.0 then begin
    p := sqrt(p);
    b := b/p;
  end
  else begin
    {note q from [11] is replaced by m}
    f := kc*kc;
    m := 1.0 - f;
    g := 1.0 - p;
    f := f - p;
    m := (b-a*p)*m;
    p := sqrt(f/g);
    a := (a-b)/g;
    b := a*p - m/(g*g*p);
  end;

  m := 1.0;
  i := 1;
  repeat
    f := a;
    a := b/p + a;
    g := e/p;
    b := 2.0*(f*g + b);
    p := g + p;
    g := m;
    m := kc+m;
    if abs(g-kc) <= 0.5e-8*g then break;
    kc:= 2.0*sqrt(e);
    e := kc*m;
    inc(i);
  until i > IMAX;
  {Test no convergence}
  if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  sfd_cel := Pi_2*(a*m+b)/(m*(m+p));
end;


{---------------------------------------------------------------------------}
function sfd_el1(x,kc: double): double;
  {-Return the incomplete elliptic integral of the 1st kind}
  { = integral(1/sqrt((1+t^2)*(1+kc^2*t^2)), t=0..x)}
const
  ca=0.5e-8;
  cb=1e-19;
var
  e,g,m,y: double;
  l: integer;
const
  IMAX = 50;
var
  i: integer;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(kc) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_el1 := NaN_d;
    exit;
  end;
  kc:= abs(kc);
  if x=0.0 then sfd_el1 := 0.0
  else if kc=0.0 then sfd_el1 := arcsinh(x)
  else if kc=1.0 then sfd_el1 := arctan(x)
  else begin
    {This is my Pascal translation of the Algol procedure el1 from [10]}
    y := abs(1.0/x);
    l := 0;
    m := 1.0;
    i := 1;
    repeat
      e := m*kc;
      g := m;
      m := kc+m;
      y := y-e/y;
      if y=0.0 then y := sqrt(e)*cb;
      if abs(g-kc) <= ca*g then break;
      kc := 2.0*sqrt(e);
      l  := l+l;
      if y<0.0 then inc(l);
      inc(i);
    until i>IMAX;
    {Test no convergence}
    if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
    if y<0.0 then l := l+1;
    e := (arctan(m/y)+Pi*l)/m;
    if x<0.0 then e := -e;
    sfd_el1 := e;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_el2(x,kc,a,b: double): double;
  {-Return the general elliptic integral of the 2nd kind, kc<>0}
  { = integral((a+b*t^2)/[(1+t^2)*sqrt((1+t^2)*(1+kc^2*t^2))], t=0..x)}
var
  y,z: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(kc) or IsNanOrInfD(a) or IsNanOrInfD(b) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_el2 := NaN_d;
    exit;
  end;
  if x=0.0 then sfd_el2 := 0.0
  else if kc=0.0 then begin
    {sfd_el2 := x*(a-b)/sqrt(1+x*x) + b*arcsinh(x)}
    y := x*(a-b)/sqrt(1.0+x*x);
    z := b*arcsinh(x);
    sfd_el2 := y+z;
  end
  else begin
    {Since V1.32.01 the Carlson form from [12] (4.17) is used because}
    {the Bulirsch algorithm el2 from [10] suffers from cancellation  }
    {for small x with certain parameters.}
    y := 1.0 + sqr(kc*x);
    z := 1.0 + sqr(x);
    sfd_el2 := a*x*sfd_ell_rf(1.0,y,z) + (b-a)/THREE*x*x*x*sfd_ell_rd(1.0,y,z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_el3(x,kc,p: double): double;
  {-Return Bulirsch's incomplete elliptic integral of the 3rd kind, 1+p*x^2<>0}
  { = integral((1+t^2)/[(1+p*t^2)*sqrt((1+t^2)*(1+kc^2*t^2))], t=0..x)}
var
  a,b,x2: double;
begin
  {Note: the el3 function from [11] replaces a former inadequate version,}
  {but it is still rather complicated and clumsy. I use the Carlson form }
  {of el3, see Carlson [12], (4.18)}
  x2:= sqr(x);
  a := 1.0+sqr(kc*x);
  b := 1.0+x2;
  sfd_el3 := x*sfd_ell_rf(1.0,a,b) + sfd_ell_rj(1.0,a,b,1.0+p*x2)*x*x2*(1.0-p)/THREE;
end;


{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Maple V style) --------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function sfd_EllipticCK(k: double): double;
  {-Return the complementary complete elliptic integral of the 1st kind, k<>0}
begin
  if k=0.0 then sfd_EllipticCK := PosInf_d
  else begin
    k := abs(k);
    if k<=1e-9 then sfd_EllipticCK := ln(4.0/k)
    else if k>=1e9 then sfd_EllipticCK := ln(4.0*k)/k
    else sfd_EllipticCK := Pi_2/sfd_agm(1.0,k);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticK(k: double): double;
  {-Return the complete elliptic integral of the 1st kind}
  { = integral(1/sqrt(1-x^2)/sqrt(1-k^2*x^2),x=0..1),    }
  { calculated with agm, real part if |k| > 1            }
var
  a,b,t: double;
const
  c = 1.039720770839917964;  {= 0.5*ln(8)}
begin
  if IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_EllipticK := NaN_d;
    exit;
  end;
  k := abs(k);
  if k>1.0 then begin
    {NIST[30], 19.7.3, first equation}
    a := k-1.0;
    if a > sqrt_epsh then sfd_EllipticK := sfd_EllipticK(1.0/k)/k
    else begin
      b := ln(a);
      {from Maple series(EllipticK(1+x), x=0, 3);}
      sfd_EllipticK := (c-0.5*b) + 0.25*a*(1-2.0*c+b);
    end;
  end
  else if k<=1.25e-4 then begin
    {K(k) = (1 + 1/4*k^2 + 9/64*k^4 + O(k^6))*Pi/2}
    sfd_EllipticK := Pi_2*(1.0 + sqr(0.5*k));
  end
  else if k=1.0 then sfd_EllipticK := PosInf_d
  else begin
    a := 1.0;
    b := sqrt((1.0-k)*(1.0+k));
    while a-b > 5e-8*a do begin
      t := a;
      a := 0.5*(t+b);
      b := sqrt(t*b);
    end;
    sfd_EllipticK := Pi/(a+b);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticKim(k: double): double;
  {-Return K(i*k), the complete elliptic integral of the 1st kind with}
  { imaginary modulus = integral(1/sqrt(1-x^2)/sqrt(1+k^2*x^2),x=0..1)}
var
  z: double;
const
  k0 = 2.01e-8;
  k1 = 2e4;
  ln4 = 1.3862943611198906188;
begin
  k := abs(k);
  {K(ik) = K(k/sqrt(1+k^2)/sqrt(1+k^2)}
  if k<=k0 then begin
    {K(ik) = Pi/2*(1 - k^2/4 + O(k^4))}
    sfd_EllipticKim := Pi_2;
  end
  else if k>=k1 then begin
    {K(ik) = [ln(4k) + 1/4(1 - ln(4k))/k^2 + O(1/k^4)]/k}
    if k>=1e12 then begin
      {Avoid overflow on 4*k or k^2}
      sfd_EllipticKim := (ln(k) + ln4)/k;
    end
    else begin
      z := ln(4.0*k);
      z := z + 0.25*(1.0 - z)/sqr(k);
      sfd_EllipticKim := z/k;
    end;
  end
  else if k>=1.0 then begin
    z := 1.0/hypot(1.0, k);
    {Avoid inaccuracies for k/sqrt(1+k^2) ~ 1 and use the complementary}
    {function with argument sqrt(1-(k/sqrt(1+k^2))^2) = 1/sqrt(1+k^2)) }
    sfd_EllipticKim := sfd_EllipticCK(z)*z
  end
  else begin
    z := hypot(1.0, k);
    sfd_EllipticKim := sfd_EllipticK(k/z)/z;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticEC(k: double): double;
  {-Return the complete elliptic integral of the 2nd kind,}
  { = integral(sqrt(1-k^2*x^2)/sqrt(1-x^2),x=0..1), real part for k>1}
var
  t: double;
begin
  if IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_EllipticEC := NaN_d;
    exit;
  end;
  k := abs(k);

  if k>1.0 then begin
    {NIST[30], 19.7.3, third equation can be written as}
    {Re(E(1/k))=k*B(k) and therefore for Re(E(k))=B(1/k)/k}
    sfd_EllipticEC := sfd_cel_b(1.0/k)/k;
    exit;
  end;

  {See Bulirsch [10], p.80}
  if abs(k)<=0.5e-8 then sfd_EllipticEC := Pi_2
  else begin
    t := (1.0-k)*(1.0+k);
    if t=0.0 then sfd_EllipticEC := 1.0
    else sfd_EllipticEC := sfd_cel2(sqrt(t),1.0,t)
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticECim(k: double): double;
  {-Return E(i*k), the complete elliptic integral of the 2nd kind with}
  { imaginary modulus = integral(sqrt(1+k^2*x^2)/sqrt(1-x^2),x=0..1)  }
var
  z: double;
const
  k0 = 2.9e-8;  { <= 2*sqrt(eps),  double = 2.9e-8}
  k1 = 4.5e8;   { double: 4.5e8}
begin
  {E(ik) = sqrt(1+k^2)*E(k/sqrt(1+k^2)}
  k := abs(k);
  if k<=k0 then begin
    {E(ik) = Pi/2*(1 + k^2/4 + O(k^4))}
    sfd_EllipticECim := Pi_2;
  end
  else if k>=k1 then begin
    {E(ik) = k*[1 + (1/4 + 1/2*ln(4) + 1/2*ln(k))/k^2 + O(1/k^4))]}
    sfd_EllipticECim := k;
  end
  else begin
    z := 1.0/hypot(1.0, k);
    sfd_EllipticECim := sfd_cel2(z,1.0,z*z)/z;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticCE(k: double): double;
  {-Return the complementary complete elliptic integral of the 2nd kind}
begin
  k := abs(k);
  if k<=0.5e-8 then sfd_EllipticCE := 1.0
  else if k>=2e8 then sfd_EllipticCE := k
  else sfd_EllipticCE := sfd_cel2(k,1.0,k*k)
end;


{---------------------------------------------------------------------------}
function sfd_EllipticF(z,k: double): double;
  {-Return the incomplete elliptic integral of the 1st kind; |z|<=1, |k*z|<1}
  { = integral(1/sqrt(1-t^2)/sqrt(1-k^2*t^2),t=0..z)}
var
  x,kc: double;
begin
  if z=0.0 then sfd_EllipticF := 0.0
  else if k=0.0 then sfd_EllipticF := arcsin(z)
  else if abs(k) > 1.0 then begin
    {change of integration variable: x = k*t}
    kc := 1.0/k;
    sfd_EllipticF := kc*sfd_EllipticF(z*k,kc);
  end
  else begin
    kc := sqrt((1.0-k)*(1.0+k));
    if z=1.0 then sfd_EllipticF := sfd_cel1(kc)
    else if z=-1.0 then sfd_EllipticF := -sfd_cel1(kc)
    else begin
      {See Bulirsch [10], p.80}
      x  := z/sqrt((1.0-z)*(1.0+z));
      sfd_EllipticF := sfd_el1(x,kc);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticE(z,k: double): double;
  {-Return the incomplete elliptic integrals of the 2nd kind, |z|<=1, |k*z| <= 1}
  { = integral(sqrt(1-k^2*x^2)/sqrt(1-x^2),x=0..z)}
var
  x,t: double;
begin
  k := abs(k);
  if z=0.0 then sfd_EllipticE := 0.0
  else if k=0.0 then sfd_EllipticE := arcsin(z)
  else if k=1.0 then sfd_EllipticE := z
  else if k<1.0 then begin
    t:= (1.0-k)*(1.0+k);
    if abs(z)=1.0 then begin
      sfd_EllipticE := copysignd(sfd_cel2(sqrt(t),1.0,t),z);
    end
    else begin
      {See Bulirsch [10], p.80}
      x  := z/sqrt((1.0-z)*(1.0+z));
      sfd_EllipticE := sfd_el2(x,sqrt(t),1.0,t);
    end
  end
  else begin
    {Use Reciprocal-Modulus Transformation, see sfd_ellint_2}
    z := z*k;
    k := 1.0/k;
    x := z/sqrt((1.0-z)*(1.0+z));
    t := sqrt((1.0-k)*(1.0+k));
    sfd_EllipticE := sfd_el2(x,t,1.0,0.0)*k;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticPi(z,nu,k: double): double;
  {-Return the incomplete elliptic integral of the 3rd kind, |z|<=1, |k*z|<1}
  { = integral(1/(1-nu*t^2)/sqrt(1-t^2)/sqrt(1-k^2*t^2),t=0..z)}
var
  x,x2,a,b,c: double;
begin
  if nu=0.0 then begin
    sfd_EllipticPi := sfd_EllipticF(z,k);
    exit;
  end;
  c := (1.0-k)*(1.0+k);
  if abs(z)=1 then begin
    {Complete integral, see Carlson [12], (4.3) with sin()=1, cos()=0}
    x := sfd_ell_rf(0.0,c,1.0) + sfd_ell_rj(0.0,c,1.0,1.0-nu)*nu/THREE;
    sfd_EllipticPi := copysignd(x,z);
  end
  else begin
    {See Carlson [12], (4.18) with Bulirsch [11], (1.1.3)}
    x := z/sqrt((1.0-z)*(1.0+z));
    x2:= x*x;
    a := 1.0 + c*x2;
    b := 1.0 + x2;
    c := 1.0 + (1.0-nu)*x2;
    sfd_EllipticPi := x*sfd_ell_rf(1.0,a,b) + sfd_ell_rj(1.0,a,b,c)*x*x2*nu/THREE;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticCPi(nu,k: double): double;
  {-Return the complementary complete elliptic integral of the 3rd kind, k<>0, nu<>1}
begin
  if nu=0.0 then sfd_EllipticCPi := sfd_EllipticCK(k)
  else if nu=1.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_EllipticCPi := NaN_d;
    exit;
  end
  else begin
    sfd_EllipticCPi := sfd_cel(k,1.0-nu,1.0,1.0);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticPiC(nu,k: double): double;
  {-Return the complete elliptic integral of the 3rd kind, |k|<>1, nu<>1; real part for |k|>1}
var
  a,q,s: double;
begin
  k := abs(k);
  if nu=0.0 then sfd_EllipticPiC := sfd_EllipticK(k)
  else if nu=1.0 then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_EllipticPiC := NaN_d;
    exit;
  end
  else if k > 1.0 then begin
    q := sqr(k);
    s := 1.0/k;
    if abs(nu-q) <= 4*eps_d*q then begin
      {Special case n = k^2, from Mathematica}
      sfd_EllipticPiC := sfd_EllipticEC(k)/((1.0+k)*(1.0-k));
    end
    else if (nu>1.0) and (q>nu) then begin
      {Include additional term a for k^2 > n > 1; (bug in Maple 7 / V)}
      {Pi(n,k) = Pi(n | k^2) = Pi(n/k^2 | 1/k^2) - a = Pi(n/k^2, 1/k) - a}
      q := nu/q;
      a := q/((nu - 1.0)*(1.0 - q));
      a := Pi_2 * sqrt(a);
      s := s * sfd_EllipticPiC(q, s);
      sfd_EllipticPiC := s - a;
    end
    else begin
      sfd_EllipticPiC := s * sfd_EllipticPiC(nu/q,s);
    end;
  end
  else begin
    a := sqrt((1.0+k)*(1.0-k));
    {See Bulirsch [11], (1.2.2)}
    sfd_EllipticPiC := sfd_cel(a,1.0-nu,1.0,1.0);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_EllipticPiCim(nu,k: double): double;
  {-Return Pi(nu, k*i), the complete elliptic integral of the 3rd kind}
  { with imaginary modulus, nu <> 1, real part if nu > 1}
var
  a,b,p,e,f,g,m,kc: double;
const
  kc_small = 1e-17;
  tol      = 0.5e-8;
const
  IMAX = 50;
var
  i: integer;
begin
  if IsNanOrInfD(nu) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_EllipticPiCim := NaN_d;
    exit;
  end;
  if nu=0.0 then begin
    sfd_EllipticPiCim := sfd_EllipticKim(k);
    exit;
  end;

  {This code is a modification of Bulirsch's [11] general}
  {complete elliptic integral cel with k replaced by i*k }
  kc := hypot(1.0,k);
  if kc=0.0 then kc := kc_small;

  e := kc;
  a := 1.0;
  b := 1.0;
  p := 1.0 - nu;

  if p>0.0 then begin
    p := sqrt(p);
    b := b/p;
  end
  else begin
    {note q from [11] is replaced by m}
    f := 1.0 + k*k;
    m := -k*k;
    g := nu;
    f := f - p;
    m := (b-a*p)*m;
    p := sqrt(f/g);
    a := (a-b)/g;
    b := a*p - m/(g*g*p);
  end;

  m := 1.0;
  i := 1;
  repeat
    f := a;
    a := b/p + a;
    g := e/p;
    b := 2.0*(f*g + b);
    p := g + p;
    g := m;
    m := kc+m;
    if abs(g-kc) <= tol*g then break;
    kc:= 2.0*sqrt(e);
    e := kc*m;
    inc(i);
  until i>IMAX;
  {Test no convergence}
  if (i>IMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  sfd_EllipticPiCim := Pi_2*(a*m+b)/(m*(m+p));
end;



{---------------------------------------------------------------------------}
{------------------- Elliptic integrals (Legendre style) -------------------}
{---------------------------------------------------------------------------}

(*
Equivalent Maple V definitions for ellint_1/2/3, jacobi_zeta, heuman_lambda.
Maple V will give wrong answers, if the normal integrals e1/2/3 are evaluated
with x < 0. Therefore signum(x) and abs(x) are used.

e1 := (x,k)   -> signum(x)*int(1/sqrt(1-k^2*sin(t)^2),t=0..abs(x));
e2 := (x,k)   -> signum(x)*int(sqrt(1-k^2*sin(t)^2),t=0..abs(x));
e3 := (x,n,k) -> signum(x)*int(1/(1-n*sin(t)^2)/sqrt(1-k^2*sin(t)^2),t=0..abs(x));

jz := (x,k) -> e2(x,k) - EllipticE(k)/EllipticK(k)*e1(x,k);
hl := (x,k) -> 2/Pi*(e1(x,sqrt(1-k^2))*(EllipticE(k)-EllipticK(k)) + EllipticK(k)*e2(x,sqrt(1-k^2)));
*)


{---------------------------------------------------------------------------}
procedure special_reduce_modpi(z: double; var a,x: double; var near_pi2: boolean);
  {-Reduce |z| mod Pi, |z| = 0.5*a*Pi + x; nearpi2 if x is near Pi_2.}
  { Internal special code for incomplete integrals and sfd_hlambda}
var
  t: double;
begin
  t := abs(z);
  {$ifdef ExtendedSyntax_on}
    rem_pio2(0.5*t,x);
  {$else}
    if 0=rem_pio2(0.5*t,x) then;
  {$endif}
  x := 2.0*x;
  if abs(abs(x)-Pi_2) <= 2*eps_d  then begin
    near_pi2 := true;
    {|z| is close to a multiple of Pi/2. Avoid rounding errors in tan(x) and}
    {assume |z| = a*Pi/2. Note that is normally applies to the odd multiples}
    {only, but the even multiples are harmless because x is near zero.}
    a := int(t/Pi_2 + 0.5);
  end
  else begin
    near_pi2 := false;
    a := 2.0*int((t-x)/Pi + 0.5);
    {If due to rounding |x| >= Pi/2 set x=0 and adjust a because}
    {tan(+-Pi_2) would have the wrong sign. Note that many cases}
    {are already captured by the near-Pi/2 test.}
    if x>=Pi_2 then begin
      a := a+1.0;
      x := 0.0;
    end
    else if x<=-Pi_2 then begin
      a := a-1.0;
      x := 0.0;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ellint_1(phi,k: double): double;
  {-Return the Legendre elliptic integral F(phi,k) of the 1st kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}
var
  a,t,x: double;
  s: integer;
  near_pi2: boolean;
begin
  if IsNanOrInfD(phi) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ellint_1 := NaN_d;
    exit;
  end;
  k := abs(k);
  if phi=0.0 then sfd_ellint_1 := 0.0
  else if k=0.0 then sfd_ellint_1 := phi
  else if k=1.0 then begin
    {If k=1 then abs(phi) must be < Pi/2, otherwise result is infinite.}
    {Basic integral F(x,1) = arcgd(x), NIST[30] 19.6.8}
    if abs(phi) < Pi_2 then sfd_ellint_1 := arcgd(phi)
    else sfd_ellint_1 := copysignd(PosInf_d,phi);
  end
  else if k>1.0 then begin
    {Use Reciprocal-Modulus Transformation - NIST[30] 19.7.4:}
    {F(phi,k) = F(arcsin(k*sin(phi)), 1/k)/k if |k*sin(phi)| <= 1}
    a := k*sin(phi);
    if (abs(a) > 1.0) or (abs(phi) > Pi_2) then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_ellint_1 := NaN_d;
      exit;
    end;
    {sfd_ellint_1 := F(arcsin(a),1/k)/k = el1(tan(arcsin(a)), sqrt(1-1/k^2))/k}
    {tan(arcsin(a)) = a/sqrt((1+a)*(1-a))}
    x := sqrt((1.0-a)*(1.0+a));
    if x=0.0 then x := sfd_EllipticK(1.0/k)
    else begin
      t := sqrt((k-1.0)*(k+1.0));
      x := sfd_el1(a/x, t/k);
    end;
    sfd_ellint_1 := x/k;
  end
  else begin
    {Reduce argument modulo Pi, x = |phi| mod Pi will be <= Pi/2.}
    {F(n*Pi+x,k) = 2*n*K(k) + F(x,k), Ref: NIST[30] 19.2.10      }
    s := isign(phi);
    special_reduce_modpi(phi, a, x, near_pi2);

    {Be careful near multiples of Pi/2. This does no really help if}
    {sfd_ellint_1 is called with double arguments. In this case the}
    {relative error will be larger than normal. If really needed,  }
    {the range reduction for double phi should be done in Specfun. }

    {Example from the test cases: ellint_1(Pi/2,.99). If calculated}
    {with the real double values the double relative error is zero:}
    {phi = 1.5707963267948965579989817342720925807952880859375}
    {k   = 0.9899999999999999911182158029987476766109466552734}
    {res = 3.3566005233611915084158314894008396005010485965489}

    if near_pi2 then begin
      {phi is close to a (odd) multiple of Pi/2. Avoid rounding errors in}
      {tan(x) and assume phi = a*Pi/2; in this case F(a*Pi/2,k) = a*K(k)}
      t := sfd_EllipticK(k);
      sfd_ellint_1 := s*a*t;
    end
    else begin
      if a<>0.0 then a := a*sfd_EllipticK(k);
      if x<>0.0 then begin
        x := tan(x);
        t := sqrt((1.0-k)*(1.0+k));
        x := sfd_el1(x, t);
      end;
      sfd_ellint_1 := s*(a+x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ellint_2(phi,k: double): double;
  {-Return the Legendre elliptic integral E(phi,k) of the 2nd kind}
  { = integral(sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1}
var
  a,t,x,w: double;
  s: integer;
  near_pi2: boolean;
begin
  if IsNanOrInfD(phi) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ellint_2 := NaN_d;
    exit;
  end;
  k := abs(k);
  if phi=0.0 then sfd_ellint_2 := 0.0
  else if k=0.0 then sfd_ellint_2 := phi
  else if k>1.0 then begin
    {Apply the Reciprocal-Modulus Transformation from NIST[30] 19.7.4:}
    {The RHS term R = (E(á,k)-k'^2 F(á,k))/k can be rewritten with the}
    {function B(á,k) = integral(cos(x)^2/sqrt(1-k^2*sin(x)^2),x=0..á))}
    {as R = B(á,k)*k. The B integral can be evaluated by a single call}
    {to the functions el2(.,k',1,0) or cel2(k',1,0), see Bulirsch [10]}
    a := k*sin(phi);
    if (abs(a) > 1.0) or (abs(phi) > Pi_2) then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_ellint_2 := NaN_d;
      exit;
    end;
    {The upper integration limit for el2 is tan(arcsin(a)) = a/w with }
    {w = sqrt((1.0-a)*(1.0+a)), see sfd_ellint_1 code. If w=0, i.e. if}
    {|k*sin(phi)| = 1 then the complete integral cel2(k',a,b) is used.}
    x := sqrt((k-1.0)*(k+1.0))/k;
    w := sqrt((1.0-a)*(1.0+a));
    if w=0.0 then x := sfd_cel2(x,1.0,0.0) else x := sfd_el2(a/w,x,1.0,0.0);
    sfd_ellint_2 := x/k;
  end
  else begin
    {Reduce argument modulo Pi, x = |phi| mod Pi will be <= Pi/2.}
    {E(n*Pi+x,k) = 2*n*E(k) + E(x,k), Ref: NIST[30] 19.2.10      }
    s := isign(phi);
    special_reduce_modpi(phi, a, x, near_pi2);
    {Be careful near multiples of Pi/2, see e.g. the problems of Cephes'}
    {qcalc for ellie(Pi/2-eps, 0.9801), which is off by E(k).           }
    {For problems with double arguments see the comment in sfd_ellint_1.}
    if near_pi2  then begin
      {phi is relative close to a multiple of Pi/2. Avoid rounding errors}
      {in tan(x) and use phi = a*Pi/2; in this case E(a*Pi/2,k) = a*E(k).}
      {Note that this test normally captures odd multiples only, but the }
      {even multiples are harmless because x is near zero.}
      t := sfd_EllipticEC(k);
      sfd_ellint_2 := s*a*t;
    end
    else begin
      t :=(1.0-k)*(1.0+k);  {kc^2}
      if t=0.0 then begin
        {k=1, basic integral E(x,1) = sin(x), E(1)=1;  NIST[30] 19.6.9}
        x := sin(x);
      end
      else begin
        if a<>0.0 then begin
          w := sfd_EllipticEC(k);
          a := a*w;
        end;
        if x<>0.0 then begin
          w := sqrt(t);
          x := tan(x);
          x := sfd_el2(x, w, 1.0, t);
        end;
      end;
      sfd_ellint_2 := s*(a+x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ellint_3(phi,nu,k: double): double;
  {-Return the Legendre elliptic integral PI(phi,nu,k) of the 3rd kind}
  { = integral(1/sqrt(1-k^2*sin(x)^2)/(1-nu*sin(x)^2),x=0..phi) with  }
  { |k*sin(phi)|<=1, returns Cauchy principal value if nu*sin(phi)^2>1}
var
  a,t,x,y,s,q: double;
  near_pi2: boolean;
begin
  if IsNanOrInfD(phi) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ellint_3 := NaN_d;
    exit;
  end;
  k := abs(k);
  if phi=0.0 then sfd_ellint_3 := 0.0
  else if nu=0.0 then begin
    {HMF[1], 17.7.18/19}
    if k=0.0 then sfd_ellint_3 := phi
    else sfd_ellint_3 := sfd_ellint_1(phi,k)
  end
  else if (k=0) and (nu=1.0) then sfd_ellint_3 := tan(phi) {HMF[1], 17.7.20}
  else if k>1.0 then begin
    {Apply the Reciprocal-Modulus Transformation from NIST[30] 19.7.4:}
    a := k*sin(phi);
    k := 1.0/k;
    if (abs(a) > 1.0) or (abs(phi) > Pi_2) then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_ellint_3 := NaN_d;
    end
    else begin
      {t := sfd_ellint_3(arcsin(a),k*k*nu,k);}
      {Because here |a|<=1 and |k|<=1 we can use sfd_EllipticPi}
      t := sfd_EllipticPi(a,k*k*nu,k);
      sfd_ellint_3 := t*k;
    end
  end
  else begin
    {Reduce argument modulo Pi, x = |phi| mod Pi will be <= Pi/2.}
    {PI(n*Pi+x,nu,k) = 2*n*PI(nu,k) + PI(x,nu,k)}
    special_reduce_modpi(phi, a, x, near_pi2);
    if near_pi2 then begin
      {phi is relative close to a multiple of Pi/2. Avoid rounding errors}
      t := sfd_EllipticPiC(nu,k);;
      if phi<0.0 then t := -t;
      sfd_ellint_3 := a*t;
    end
    else begin
      if a<>0.0 then a :=a*sfd_EllipticPiC(nu,k);
      {For the incomplete integral with |x|<=Pi/2 see Carlson [12], (4.3)}
      sincos(x,s,x);
      q := sqr(s);
      t := 1.0 - nu*q;
      if t=0.0 then begin
        {handle case nu*sin(x)^2 = 1 here, sfd_ell_rj returns NAN/RTE for t=0}
        t := PosInf_d;
      end
      else begin
        {Cauchy principal value for t < 0}
        x := sqr(x);
        y := 1.0 - sqr(k)*q;
        t := sfd_ell_rj(x, y, 1.0, t);
        t := nu*t*q/THREE;
        y := sfd_ell_rf(x, y, 1.0);
        t := s*(y+t)+a;
      end;
      if phi>0.0 then sfd_ellint_3 := t
      else sfd_ellint_3 := -t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_cel_d(k: double): double;
  {-Return the complete elliptic integral D(k) = (K(k) - E(k))/k^2, real part for |k|>1}
var
  t: double;
begin
  if IsNanOrInfd(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_cel_d := NaN_d;
    exit;
  end;
  k := abs(k);
  if k > 1.0 then begin
    t := 1.0/k;
    sfd_cel_d := sfd_cel_d(t)*t*t*t;
  end
  else if k<=2e-5 then begin
    {D(k) = (1 + 3/8*k^2 + 15/64*k^4 + O(k^6))*Pi/4}
    sfd_cel_d := Pi_4*(1.0 + 0.375*sqr(k));
  end
  else begin
    t := (1.0-k)*(1.0+k);
    if t=0.0 then sfd_cel_d := PosInf_d
    else sfd_cel_d := sfd_cel2(sqrt(t),0.0,1.0)
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ellint_d(phi,k: double): double;
  {-Return the Legendre elliptic integral D(phi,k) = (F(phi,k) - E(phi,k))/k^2 }
  { = integral(sin(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }
var
  a,t,x,w: double;
  s: integer;
  near_pi2: boolean;
begin
  k := abs(k);
  if phi=0.0 then sfd_ellint_d := 0.0
  else if k>1.0 then begin
    {Use Reciprocal-Modulus Transformation: derived from NIST[30] 19.7.4:}
    {D(phi,k) = D(arcsin(k*sin(phi)), 1/k)/k^3 if |k*sin(phi)| <= 1}
    a := k*sin(phi); { = sin(beta)}
    if (abs(a) > 1.0) or (abs(phi) > Pi_2) then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_ellint_d := NaN_d;
      exit;
    end;
    t := (1.0-a)*(1.0+a); {cos(beta)^2}
    x := a/k;
    t := sfd_ell_rd(t, (1.0-x)*(1.0+x), 1.0);
    sfd_ellint_d  := x*x*x*t/THREE;
  end
  else begin
    {Reduce argument modulo Pi, x = |phi| mod Pi will be <= Pi/2.}
    {D(n*Pi+x,k) = 2*n*D(k) + D(x,k), Ref: NIST[30] 19.2.10      }
    s := isign(phi);
    special_reduce_modpi(phi, a, x, near_pi2);
    if near_pi2  then begin
      {phi is relative close to a multiple of Pi/2; avoid rounding errors}
      t := sfd_cel_d(k);
      sfd_ellint_d := s*a*t;
    end
    else begin
      if a<>0.0 then begin
        w := sfd_cel_d(k);
        a := a*w;
      end;
      if x<>0.0 then begin
        {Carlson[12], 4.4: D = sin(phi)^3/3*R_D(cos(phi)^2, 1-(k*sin(phi))^2,1)}
        sincos(x,w,t);
        x := w*k;
        x := sfd_ell_rd(t*t, (1.0-x)*(1.0+x), 1.0);
        x := w*w*w*x/THREE;
      end;
      sfd_ellint_d := s*(a+x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_cel_b(k: double): double;
  {-Return the complete elliptic integral B(k) = (E(k) - kc^2*K(k))/k^2, real part for |k|>1}
var
  t: double;
begin
  if IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_cel_b := NaN_d;
    exit;
  end;
  k := abs(k);
  if k > 1.0 then begin
    t := 1.0/k;
    sfd_cel_b := sfd_EllipticEC(t)*t;
  end
  else if k<=2e-5 then begin
    {B(k) = (1 + 1/8*k^2 + 3/64*k^4 + O(k^6))*Pi/4}
    sfd_cel_b := Pi_4*(1.0 + 0.125*sqr(k));
  end
  else begin
    t := (1.0-k)*(1.0+k);
    if t=0.0 then sfd_cel_b := 1.0
    else sfd_cel_b := sfd_cel2(sqrt(t),1.0,0.0)
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ellint_b(phi,k: double): double;
  {-Return the Legendre elliptic integral B(phi,k) = (E(phi,k) - kc^2*F(phi,k))/k^2}
  { = integral(cos(x)^2/sqrt(1-k^2*sin(x)^2),x=0..phi), |k*sin(phi)| <= 1      }
var
  a,t,x,w: double;
  s: integer;
  near_pi2: boolean;
begin
  k := abs(k);
  if phi=0.0 then sfd_ellint_b := 0.0
  else if k>1.0 then begin
    {Use Reciprocal-Modulus Transformation: derived from NIST[30] 19.7.4:}
    {B(phi,k) = E(arcsin(k*sin(phi)), 1/k)/k if |k*sin(phi)| <= 1}
    a := k*sin(phi);
    if (abs(a) > 1.0) or (abs(phi) > Pi_2) then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_ellint_b := NaN_d;
      exit;
    end;
    x := arcsin(a);
    t := sfd_ellint_2(x,1.0/k);
    sfd_ellint_b := t/k;
  end
  else begin
    {Reduce argument modulo Pi, x = |phi| mod Pi will be <= Pi/2.}
    {B(n*Pi+x,k) = 2*n*B(k) + B(x,k)}
    s := isign(phi);
    special_reduce_modpi(phi, a, x, near_pi2);
    if near_pi2  then begin
      {phi is relative close to a multiple of Pi/2; avoid rounding errors}
      t := sfd_cel_b(k);
      sfd_ellint_b := s*a*t;
    end
    else begin
      if a<>0.0 then begin
        w := sfd_cel_b(k);
        a := a*w;
      end;
      if x<>0.0 then begin
        x := tan(x);
        x := sfd_el2(x, sqrt((1.0-k)*(1.0+k)), 1.0, 0.0);
      end;
      sfd_ellint_b := s*(a+x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_hlambda(phi,k: double): double;
  {-Return Heuman's function Lambda_0(phi,k) = F(phi,k')/K(k') + 2/Pi*K(k)*Z(phi,k'), |k|<=1}
var
  a,t,x,p: double;
  s: integer;
  near_pi2: boolean;
begin
  if IsNanOrInfD(phi) or IsNanOrInfD(k) or (abs(k) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_hlambda := NaN_d;
    exit;
  end;

  {In HMF[1] notation Lambda_0(phi\alpha) =  Lambda_0(phi,sin(alpha))}
  {Lambda_0(phi,k) = F(phi,k')/K(k') + 2/Pi*K(k)*Z(phi,k')          HMF[1] 17.4.39}
  {Lambda_0(phi,k) = 2/Pi*[K(k)*E(phi,k') - (K(k)-E(k))*F(phi,k')]  HMF[1] 17.4.40}

  k := abs(k);
  if phi=0.0 then sfd_hlambda := 0.0
  else if k=1.0 then sfd_hlambda := phi/Pi_2
  else begin
    s := isign(phi);
    {Reduce argument modulo Pi, x = |phi| mod Pi will be <= Pi/2. Since }
    {Z(x,k) is periodic and F(n*Pi+x,k') = 2*n*K(k') + F(x,k)', Lambda_0}
    {is quasi-periodic with lambda(n*Pi+x,k) = 2*n + lambda(x,k)}
    special_reduce_modpi(phi, a, x, near_pi2);
    if near_pi2 then begin
      {phi is very close to a multiple of Pi/2. Avoid rounding errors}
      {and assume phi = a*Pi/2; in this case sfd_hlambda(phi,k) = a  }
      sfd_hlambda := s*a;
    end
    else begin
      {This is the base case with |phi| < Pi/2, use a single call to cel}
      t := (1.0-k)*(1.0+k);  {=k'^2}
      {Lambda_0(phi,k) = 2/Pi*sqrt(p)*sin(phi)*cel(k',p,1,k'^2) }
      {with p = 1 + (k*tan(phi))^2), see Bulirsch[11], (1.2.3)  }
      p := 1.0 + sqr(k*tan(x));
      x := sin(x)*sqrt(p)/Pi_2;
      t := sfd_cel(sqrt(t), p, 1.0, t);
      sfd_hlambda := s*(a + t*x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jzeta(phi,k: double): double;
  {-Return the Jacobi Zeta function Z(phi,k) = E(phi,k) - E(k)/K(k)*F(phi,k), |k|<=1}
var
  a,b,p,e,f,g,m,kc,t: double;
  i: integer;
  done: boolean;
begin
  {Info: Maple_JacobiZeta(z,k) = sfd_jzeta(am(z,k),k)) }
  {Mathematica_JacobiZeta(z,m) = sfd_jzeta(z,sqrt(m))  }

  {Z is periodic for phi: Z(phi + n*Pi,k) = Z(phi,k)}
  {Ref: http://functions.wolfram.com/08.07.04.0004.01}

  if IsNanOrInfD(phi) or IsNanOrInfD(k) or (abs(k) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jzeta := NaN_d;
    exit;
  end;
  kc := sqrt((1.0-k)*(1.0+k));
  if (phi=0.0) or (k=0.0) then sfd_jzeta := 0.0
  else if kc=0.0 then begin
    {Z(phi + n*Pi, k) = Z(phi, k) is automatically handled except for k=1}
    {Z(x,1)=sin(x) if abs(x) <= Pi/2; Ref: http://functions.wolfram.com/08.07.03.0002.01}
    {$ifdef ExtendedSyntax_on}
      {avoid warning}
      rem_pio2(0.5*phi,t);
    {$else}
      i := rem_pio2(0.5*phi,t);
    {$endif}
    sfd_jzeta := sin(2.0*t);
  end
  else begin
    {The general case for kc<>0 is based on the Algol procedure}
    {cel and the formulas (1.2.3), (1.2.4a) from Bulirsch [11].}
    sincos(phi,a,b);
    p := b*b + sqr(a*kc);
    t := k*k*a*b;
    e := kc;
    p := sqrt(p);
    a := 0.0;
    b := 1.0;
    m := 1.0;
    i := 0;
    done := false;
    repeat
      inc(i);
      f := a;
      a := b/p + a;
      g := e/p;
      b := 2.0*(f*g + b);
      p := g + p;
      g := m;
      m := kc+m;
      if abs(g-kc) <= 1e-8*g then done := true
      else begin
        kc:= 2.0*sqrt(e);
        e := kc*m;
      end;
    until done;
    sfd_jzeta := t*ldexpd(a,-i);
  end;
end;

{---------------------------------------------------------------------------}
{---------------- Elliptic integrals (Mathematica style) -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function sfd_M_EllipticK(m: double): double;
  {-Return the complete elliptic integral of the 1st kind,}
  { K(m) = integral(dx/sqrt(1-m*sin(x)^2),x=0..Pi/2), real part for m>1}
begin
  if m>=0.0 then sfd_M_EllipticK := sfd_EllipticK(sqrt(m))
  else sfd_M_EllipticK := sfd_EllipticKim(sqrt(-m))
end;


{---------------------------------------------------------------------------}
function sfd_M_EllipticEC(m: double): double;
  {-Return the complete elliptic integral of the 2nd kind,}
  { E(m) = integral(sqrt(1-m*sin(x)^2),x=0..Pi/2), real part for m>1}
begin
  if m>=0.0 then sfd_M_EllipticEC := sfd_EllipticEC(sqrt(m))
  else sfd_M_EllipticEC := sfd_EllipticECim(sqrt(-m))
end;


{---------------------------------------------------------------------------}
function sfd_M_EllipticPiC(n,m: double): double;
  {-Return the complete elliptic integral of the 3rd kind, n<>1, m<>1, real part}
  { for m>1; Pi(n|m) = integral(1/(1-n*sin(x)^2/sqrt(1-m*sin(x)^2)), x=0..Pi/2)}
begin
  if m>=0.0 then sfd_M_EllipticPiC := sfd_EllipticPiC(n, sqrt(m))
  else sfd_M_EllipticPiC := sfd_EllipticPiCim(n, sqrt(-m))
end;


{---------------------------------------------------------------------------}
function sfd_M_EllipticF(phi,m: double): double;
  {-Return the incomplete elliptic integral of the 1st kind,}
  { F(phi,m) = integral(dx/sqrt(1-m*sin(x)^2),x=0..phi), m*sin(phi)^2 <= 1}
var
  x,a,c,s,y: double;
  near_pi2: boolean;
  sig: integer;
begin
  if m=0.0 then begin
    sfd_M_EllipticF := phi;
    exit;
  end;
  sig := isign(phi);
  special_reduce_modpi(phi, a, x, near_pi2);
  if near_pi2 then begin
    x := sfd_M_EllipticK(m);
    sfd_M_EllipticF := sig*a*x;
  end
  else begin
    if a<>0.0 then a := a*sfd_M_EllipticK(m);
    if x<>0.0 then begin
      {Ref: Carlson [12], (4.1)}
      sincos(x,s,c);
      y := 1.0 - m*s*s;
      {Note: there may be large relative errors,}
      {if y is close to 0 e.g. for m=2, phi=pi/4}
      if abs(y) <= eps_d then y := 0.0;
      {if y<0 RF returns RTE or Nan}
      x := s*sfd_ell_rf(c*c, y, 1.0);
    end;
    sfd_M_EllipticF := sig*(a+x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_M_EllipticE(phi,m: double): double;
  {-Return the incomplete elliptic integral of the 2nd kind,}
  { E(phi,m) = integral(sqrt(1-m*sin(x)^2),x=0..phi), m*sin(phi)^2 <= 1}
var
  x,y,a,c,s: double;
  near_pi2: boolean;
  sig: integer;
begin
  sig := isign(phi);
  special_reduce_modpi(phi, a, x, near_pi2);
  if near_pi2 then begin
    x := sfd_M_EllipticEC(m);
    sfd_M_EllipticE := sig*a*x;
  end
  else begin
    if a<>0.0 then a := a*sfd_M_EllipticEC(m);
    if x<>0.0 then begin
      {Ref: Carlson [12], (4.2)}
      sincos(x,s,c);
      x := c*c;
      y := 1.0 - m*s*s;
      {Note: there may be large relative errors,}
      {if y is close to 0 e.g. for m=2, phi=pi/4}
      if abs(y) <= eps_d then y := 0.0;
      {if y<0 RF or RD returns RTE or Nan}
      x := sfd_ell_rf(x, y, 1.0) - sfd_ell_rd(x, y, 1.0)*sqr(s)*m/THREE;
    end;
    sfd_M_EllipticE := sig*(s*x + a);
  end;
end;


{$define UseMathematicaMode}

{---------------------------------------------------------------------------}
function sfd_M_EllipticPi(n,phi,m: double): double;
  {-Return the incomplete elliptic integral Pi(n,phi,m) of the 3rd kind}
  { = integral(1/sqrt(1-m*sin(x)^2)/(1-n*sin(x)^2),x=0..phi), with n<>1}
  { m*sin(phi)^2 <= 1, real part if complex.}
var
  a,t,x,y,s,q,c: double;
  near_pi2: boolean;
begin
  sfd_M_EllipticPi := NaN_d;
  if IsNanOrInfD(phi) or IsNanOrInfD(m) or IsNanOrInfD(n) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;
  if phi=0.0 then sfd_M_EllipticPi := 0.0
  else if n=0.0 then begin
    {http://functions.wolfram.com/08.06.03.0007.01}
    sfd_M_EllipticPi := sfd_M_EllipticF(phi, m);
  end
  else if (n=1.0) and (m<>1.0) then begin
    {http://functions.wolfram.com/08.06.03.0008.01}
    sincos(phi,s,c);
    x := 1.0 - m*sqr(s);
    if x < 0.0 then exit;
    t := sqrt(x)*s/c;
    t := t - sfd_M_EllipticE(phi, m);
    sfd_M_EllipticPi := t/(1.0 - m) + sfd_M_EllipticF(phi, m);
  end
  else begin
    {Reduce argument modulo Pi, x = |phi| mod Pi will be <= Pi/2.}
    {PI(n,j*Pi+x,m) = 2*j*PI(n,k) + PI(n,x,m)}
    special_reduce_modpi(phi, a, x, near_pi2);
    if near_pi2 then begin
      {phi is relative close to a multiple of Pi/2. Avoid rounding errors}
      t := sfd_M_EllipticPiC(n,m);
      if phi<0.0 then t := -t;
      sfd_M_EllipticPi := a*t;
    end
    else begin
      {Incomplete integral with |x|<=Pi/2, see Carlson [12], (4.3)}
      sincos(x,s,c);
      q := sqr(s);
      t := 1.0 - n*q;
      y := 1.0 - m*q;
      {$ifdef UseMathematicaMode}
        {Mathematica (real part) results are the same as the (real part)}
        {results from Maple, MPMath, Maxima if n < 1 or |phi| <= Pi/2.}
        if n > 1.0 then a := 0.0;
      {$endif}
      if a<>0.0 then begin
        if n=1.0 then exit; {and return NaN}
        a := a*sfd_M_EllipticPiC(n,m);
      end;
      if (abs(t)<=eps_s) or (y<0.0) then begin
        {n*sin^2(phi)=1 or m*sin^2(phi) > 1, return NaN}
        exit;
      end;
      if y<=eps_s then begin
        {m*sin^2(phi) ~ 0, e.g. for phi=pi/4, m=2}
        {Avoid rounding artefacts and set zero.  }
        y := 0.0;
      end;
      x := sqr(c);
      {Cauchy principal value for t < 0}
      t := sfd_ell_rj(x, y, 1.0, t);
      t := n*t*q/THREE;
      y := sfd_ell_rf(x, y, 1.0);
      t := s*(y+t) + a;
      if phi>0.0 then sfd_M_EllipticPi := t
      else sfd_M_EllipticPi := -t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
{------------------- Jacobi elliptic and theta functions -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function sfd_jam(x,k: double): double;
  {-Return the Jacobi amplitude am(x,k)}
var
  a,t,z,sn,cn,dn: double;
  s: integer;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jam := NaN_d;
    exit;
  end;
  k := abs(k);
  if k=0.0 then begin
    {NIST[30] 22.16.4}
    sfd_jam := x;
  end
  else if k=1.0 then begin
    {NIST[30] 22.16.5}
    sfd_jam := gd(x);
  end
  else if k>1.0 then begin
    {For |k| > 1, am(x,k) is periodic, see NIST[30] 22.19(i)}
    {Note, that x should not be much larger than the period }
    {4*K(1/k)/k in order to obtain a small relative error!  }
    sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
    sfd_jam := arctan2(sn,cn);
  end
  else begin
    {Use am(-x,k) = -am(x,k) for negative x}
    s := isign(x);
    z := abs(x);
    {Skip range reduction if z < Pi_2 < K}
    if z < Pi_2 then a := 0.0
    else begin
      {Use quasi-periodicity am(x + 2K,k) = am(x,k) + Pi (NIST[30] 22.16.2)}
      {to reduce the range to |z| <= K (where standard sfd_sncndn applies).}
      a := 2.0*sfd_EllipticK(k);
      t := floord(z/a+0.5);
      z := z - t*a;
      a := t*pi;
    end;
    {sfd_sncndn uses the complementary parameter mc = 1 - k^2 = (1+k)*(1-k)}
    t := (1.0+k)*(1.0-k);
    sfd_sncndn(z,t,sn,cn,dn);
    {am(z,k) = arctan2(sn,cn) is slightly more accurate than arcsin(sn)}
    t := arctan2(sn,cn);
    sfd_jam := s*(a + t);
  end;
end;


{---------------------------------------------------------------------------}
procedure sfd_sncndn(x,mc: double; var sn,cn,dn: double);
  {-Return the Jacobi elliptic functions sn,cn,dn for argument x and}
  { complementary parameter mc.}
const
  ca = 0.5e-8;
  NA = 17;
var
  a,b,c,d: double;
  i,l: integer;
  bo: boolean;
  m,n: array[0..NA] of double;
begin
  (* The other Jacobi elliptic functions can be derived from sn,cn,dn:
   *
   * nc(x,mc) = 1/cn(x,mc)
   * sc(x,mc) = sn(x,mc)/cn(x,mc)
   * dc(x,mc) = dn(x,mc)/cn(x,mc)
   *
   * nd(x,mc) = 1/dn(x,mc)
   * sd(x,mc) = sn(x,mc)/dn(x,mc)
   * cd(x,mc) = cn(x,mc)/dn(x,mc)
   *
   * ns(x,mc) = 1/sn(x,mc)
   * cs(x,mc) = cn(x,mc)/sn(x,mc)
   * ds(x,mc) = dn(x,mc)/sn(x,mc)
   *
   * and the amplitude am(x,mc) = arctan2(sn,cn) for the basic
   * interval |x| < K(mc); for larger x values use sfd_jam
   *)

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
    {This is my Pascal translation of the Algol procedure sncndn from [10]}
    bo := mc<0.0;
    d  := 1.0-mc; {Avoid compiler warning, move in front of 'if'}
    if bo then begin
      {Here m=1-mc>1, apply Jacobi real transformation, HMF[1] 16.11}
      mc := -mc/d;
      d  := sqrt(d);
      x  := d*x;
    end;
    a  := 1.0;
    dn := 1.0;
    for i:=0 to NA do begin
      l   := i;
      m[i]:= a;
      mc  := sqrt(mc);
      n[i]:= mc;
      c   := 0.5*(a+mc);
      if abs(a-mc) <= ca*a then break;
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
      {Here c = cot(am(x,mc)), ie am = arccot(c)}
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
function sfd_jacobi_sn(x,k: double): double;
  {-Return the Jacobi elliptic function sn(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_sn := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_sn := sn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_cn(x,k: double): double;
  {-Return the Jacobi elliptic function cn(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_cn := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_cn := cn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_dn(x,k: double): double;
  {-Return the Jacobi elliptic function dn(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_dn := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_dn := dn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_nc(x,k: double): double;
  {-Return the Jacobi elliptic function nc(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_nc := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_nc := 1.0/cn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_sc(x,k: double): double;
  {-Return the Jacobi elliptic function sc(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_sc := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_sc := sn/cn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_dc(x,k: double): double;
  {-Return the Jacobi elliptic function dc(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_dc := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_dc := dn/cn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_nd(x,k: double): double;
  {-Return the Jacobi elliptic function nd(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_nd := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_nd := 1.0/dn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_sd(x,k: double): double;
  {-Return the Jacobi elliptic function sd(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_sd := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_sd := sn/dn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_cd(x,k: double): double;
  {-Return the Jacobi elliptic function cd(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_cd := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_cd := cn/dn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_ns(x,k: double): double;
  {-Return the Jacobi elliptic function ns(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_ns := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_ns := 1.0/sn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_cs(x,k: double): double;
  {-Return the Jacobi elliptic function cs(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_cs := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_cs := cn/sn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_ds(x,k: double): double;
  {-Return the Jacobi elliptic function ds(x,k)}
var
  sn,cn,dn: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_ds := NaN_d;
    exit;
  end;
  sfd_sncndn(x,(1.0-k)*(1.0+k),sn,cn,dn);
  sfd_jacobi_ds := dn/sn;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arcsn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsn(x,k), |x| <= 1 and |x*k| <= 1}
var
  a: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (abs(x) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcsn := NaN_d;
    exit;
  end;
  k := abs(k);
  {arcsn(x,k) = F(arcsin(x),k); arcsn(x,0)=arcsin(x); arcsn(x,1)=arctanh(x)}
  if k=1.0 then begin
    sfd_jacobi_arcsn := arctanh(x);
  end
  else begin
    a := arcsin(x);
    sfd_jacobi_arcsn := sfd_ellint_1(a,k);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arccn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccn(x,k), |x| <= 1, x >= sqrt(1 - 1/k^2) if k >= 1}
var
  a: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (abs(x) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arccn := NaN_d;
    exit;
  end;
  k := abs(k);
  {arccn(x,k) = F(acos(x),k); arccn(x,0)=arccos(x); arccn(x,1)=arcsech(x)}
  if k=1.0 then begin
    {http://functions.wolfram.com/09.38.03.0003.01}
    sfd_jacobi_arccn := arcsech(x);
  end
  else if x=0.0 then begin
    {http://functions.wolfram.com/09.38.03.0006.01}
    sfd_jacobi_arccn := sfd_EllipticK(k);
  end
  else begin
    a := arccos(x);
    sfd_jacobi_arccn := sfd_ellint_1(a,k);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arcdn(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcdn(x,k), 0 <= x <= 1, x^2 + k^2 > 1 if |k| < 1;  |x| <= 1 if |k| > 1}
var
  a,b: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (k=0.0) or (abs(x) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcdn := NaN_d;
    exit;
  end;
  k := abs(k);
  {arcdn(x,k) = F(arcsin(sqrt(1 - x*x)/k),k); arcdn(x,1)=arcsech(x)}
  if k=1.0 then begin
    {http://functions.wolfram.com/09.41.03.0003.01}
    sfd_jacobi_arcdn := arcsech(x);
  end
  else begin
    a := sqrt((1.0-x)*(1.0+x))/k;
    a := arcsin(a);
    a := sfd_ellint_1(a,k);
    if (x<0.0) and (k>1.0) then begin
      b := 2.0*sfd_EllipticK(1.0/k)/k;
      sfd_jacobi_arcdn := b-a;
    end
    else sfd_jacobi_arcdn := a;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arcsc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsc(x,k), |x| <= 1/sqrt(k^2-1) for |k|>1}
var
  a: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcsc := NaN_d;
    exit;
  end;
  k := abs(k);
  {arcsc(x,k) = F(arctan(x),k); arcsc(x,0)=arctan(x); arcsc(x,1)=arcsinh(x)}
  if k=1.0 then begin
    {http://functions.wolfram.com/09.46.03.0003.01}
    sfd_jacobi_arcsc := arcsinh(x);
  end
  else begin
    a := arctan(x);
    sfd_jacobi_arcsc := sfd_ellint_1(a,k);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arccs(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccs(x,k), |x| >= sqrt(k^2-1) for |k|>1}
var
  a: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arccs := NaN_d;
    exit;
  end;
  k := abs(k);
  {arccs(x,k) = F(arccot(x),k); arccs(x,1)=arccsch(x)}
  if k=1.0 then begin
    {http://functions.wolfram.com/09.39.03.0003.01}
    sfd_jacobi_arccs := arccsch(x);
  end
  else begin
    a := arccot(x);
    sfd_jacobi_arccs := sfd_ellint_1(a,k);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arcnc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcnc(x,k), x >= 1, x^2 <= k^2/(k^2-1) for |k|>1}
var
  a: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (x<1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcnc := NaN_d;
    exit;
  end;
  k := abs(k);
  {arcnc(x,k) = F(arcsec(x),k); arcnc(x,1)=arccosh(x)}
  if k=1.0 then begin
    {http://functions.wolfram.com/09.43.03.0003.01}
    sfd_jacobi_arcnc := arccosh(x);
  end
  else begin
    a := arcsec(x);
    sfd_jacobi_arcnc := sfd_ellint_1(a,k);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arcnd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcnd(x,k), x >= 1, x^2 <= k^2/(1-k^2) if k < 1}
var
  a,b: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (x<1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcnd := NaN_d;
    exit;
  end;
  k := abs(k);
  {arcnd(x,k) = F(arcsin(sqrt((x^2 - 1)/|kx|),k); arcnd(x,1)=arccosh(x)}
  if k=1.0 then begin
    {http://functions.wolfram.com/09.44.03.0003.01}
    sfd_jacobi_arcnd := arccosh(x);
  end
  else if (k>1.0) and (x>=2.0) then begin
    {Use arcnd(x,k) = arcnc(x, 1/k)/k; this is better in practice,}
    {especially for large x}
    a := sfd_ellint_1(arcsec(x),1.0/k);
    sfd_jacobi_arcnd := a/k;
  end
  else begin
    a := (x-1.0)*(x+1.0);
    b := sqrt(a)/abs(k*x);
    sfd_jacobi_arcnd := sfd_ellint_1(arcsin(b),k);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arcns(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcns(x,k), |x|>=1, |x|>=k if k>=1}
var
  a: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (abs(x) < 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcns := NaN_d;
    exit;
  end;
  k := abs(k);
  {arcns(x,k) = F(arctan(x),k); arcsc(x,1)=arccoth(x)}
  if k=1.0 then begin
    {http://functions.wolfram.com/09.45.03.0003.01}
    sfd_jacobi_arcns := arccoth(x);
  end
  else begin
    {http://functions.wolfram.com/09.45.27.0014.01}
    a := arccsc(x);
    sfd_jacobi_arcns := sfd_ellint_1(a,k);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arccd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arccd(x,k); |x| <= 1 if |k| < 1; |x| >= 1 if |k| > 1 }
var
  a,b,f: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arccd := NaN_d;
    exit;
  end;
  k := abs(k);
  a := abs(x);
  if (k=1.0) or ((k<1.0) and (a>1.0)) or ((k>1.0) and (a<1.0)) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arccd := NaN_d;
    exit;
  end;

  if x=1.0 then begin
    {http://functions.wolfram.com/09.37.03.0008.01}
    sfd_jacobi_arccd := 0;
    exit;
  end;

  if (k>1.0) and (x <= -1.0) then begin
    a := sfd_EllipticK(1.0/k)/k;
    b := sfd_jacobi_arccd(-x,k);
    sfd_jacobi_arccd := 2.0*a - b;
    exit;
  end;

  if k=0.0 then begin
    {http://functions.wolfram.com/09.37.03.0001.01}
    sfd_jacobi_arccd := arccos(x);
    exit;
  end;

  f := 1.0;
  if x>=16.0 then begin
    {Here k>1 and 'large' x. Use argument AND modulus transformation:}
    {arccd(x,k) = arcdc(1/x,k) = arccd(1/x,1/k)/k.}
    x := 1.0/x;
    k := 1.0/k;
    f := k;
  end
  else if x>0.0625 then begin
    a := k*x;
    b := (1.0-a)*(1.0+a);
    a := (1.0-x)*(1.0+x)/b;
    b := arcsin(sqrt(a));
    sfd_jacobi_arccd := sfd_ellint_1(b,k);
    exit;
  end;

  {Here k<1: arccd(x,k) = K(k) - arcsn(x,k)}
  {This also gives the correct value for x<0}
  b := sfd_EllipticK(k);
  if x<>0.0 then begin
    {http://functions.wolfram.com/09.37.27.0012.01}
    a := sfd_jacobi_arcsn(x,k);
    b := b - a;
  end;
  sfd_jacobi_arccd := b*f;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arcdc(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcdc(x,k); |x| >= 1 if |k| < 1; |x| <= 1 if |k| > 1 }
var
  a,b: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcdc := NaN_d;
    exit;
  end;
  k := abs(k);
  a := abs(x);

  if (k=1.0) or ((k<1.0) and (a<1.0)) or ((k>1.0) and (a>1.0)) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcdc := NaN_d;
    exit;
  end;

  if x=1.0 then begin
    {http://functions.wolfram.com/09.40.03.0008.01}
    sfd_jacobi_arcdc := 0;
    exit;
  end;

  {arcdc(x,k) = F(arcsin(sqrt(x^2 - 1)/(x^2 - k^2),k); arcdc(x,0)=arcsec(x)}
  if k=0.0 then begin
    {http://functions.wolfram.com/09.40.03.0001.01}
    sfd_jacobi_arcdc := arcsec(x);
  end
  else if a >= 16.0 then begin
    sfd_jacobi_arcdc := sfd_jacobi_arccd(1.0/x,k);
  end
  else if a <= 0.0625 then begin
    {Here k > 1, use arcdc(x,k) = (K(1/k)-arcsn(x,1/k))/k}
    {http://functions.wolfram.com/09.40.27.0012.01}
    k := 1.0/k;
    if a<>0.0 then a := sfd_jacobi_arcsn(x,k);
    b := sfd_EllipticK(k);
    sfd_jacobi_arcdc := (b-a)*k;
  end
  else begin
    a := (x-1.0)*(x+1.0);
    b := (x-k)*(x+k);
    a := a/b;
    a := sqrt(a);
    b := sfd_ellint_1(arcsin(a),k);
    if x<0.0 then begin
      if k<1.0 then a := sfd_EllipticK(k)
      else a := sfd_EllipticK(1.0/k)/k;
      sfd_jacobi_arcdc := 2.0*a-b;
    end
    else sfd_jacobi_arcdc := b;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arcsd(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcsd(x,k), x^2*(1-k^2) <= 1}
var
  a: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcsd := NaN_d;
    exit;
  end;
  {arcsd(x,k) = F(arcsin(|x|/sqrt(1+(kx)^2),k); arcsd(x,1)=arcsinh(x)}
  k := abs(k);
  if k=1.0 then begin
    {http://functions.wolfram.com/09.47.03.0003.01}
    sfd_jacobi_arcsd := arcsinh(x);
  end
  else begin
    if k>1.0 then begin
      {use arcsd(x,k) = arcsc(k*x, 1/k)/k}
      a := arctan(k*x);
      sfd_jacobi_arcsd := sfd_ellint_1(a,1.0/k)/k;
    end
    else begin
      a := x/hypot(1.0, k*x);
      sfd_jacobi_arcsd := sfd_ellint_1(arcsin(a),k);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_jacobi_arcds(x,k: double): double;
  {-Return the inverse Jacobi elliptic function arcds(x,k), x^2 + k^2 >= 1}
var
  a: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jacobi_arcds := NaN_d;
    exit;
  end;
  k := abs(k);
  {arcds(x,k) = F(arcsin(1/sqrt(k^2+x^2)),k); arcds(x,1)=arccsch(x)}
  if k=1.0 then begin
    {http://functions.wolfram.com/09.42.03.0003.01}
    sfd_jacobi_arcds := arccsch(x);
  end
  else if k > 1.0 then begin
    {use arcds(x,k) = arccs(x/k, 1/k)/k}
    a := arccot(x/k);
    a := sfd_ellint_1(a,1.0/k);
    sfd_jacobi_arcds := a/k;
  end
  else begin
    a := 1.0/hypot(k, x);
    if x<0.0 then a := -a;
    sfd_jacobi_arcds := sfd_ellint_1(arcsin(a),k);
  end;
end;


(*
Theta functions for argument 0 and small nome q use the Fourier series
given in NIST[30], 20.2.2/3/4 and HMF[1] 16.27.2/3/4:

theta2(q) = 2*q^(1/4)*sum(q^(n*(n+1)),n=0..Inf)
theta3(q) = 1 + 2*sum(q^(n*n),n=1..Inf)
theta4(q) = 1 + 2*sum((-1)^n*q^(n*(n+1)),n=1..Inf) = theta3(-q)
*)

{---------------------------------------------------------------------------}
function theta2small(q: double): double;
  {-Return theta2(q) = 2*q^(1/4)*sum(q^(n*(n+1)),n=0..Inf), small q >= 0}
var
  p,t,s,f: double;
begin
  s := 0.0;
  if q<>0.0 then begin
    t := sqr(q);
    p := 1.0;
    f := 2.0*sqrt(sqrt(q));
    while f >= eps_d*s do begin
      s := s + f;
      p := p*t;
      f := f*p;
    end;
  end;
  theta2small := s;
end;


{---------------------------------------------------------------------------}
function theta3small(q: double): double;
  {-Return theta3(q) = 1 + 2*sum(q^(n*n),n=1..Inf) for small q}
var
  p,t,s,f: double;
begin
  t := sqr(q);
  p := q;
  f := 2.0*q;
  s := 1.0;
  while abs(f) >= eps_d*s do begin
    s := s + f;
    p := p*t;
    f := f*p;
  end;
  theta3small := s;
end;


(*
As described in NIST[30] 20.7(viii) "Transformations of Lattice Parameter"
for z=0 and q = exp(i*Pi*tau) the theta q series can be restricted to very
small q, theoretically to q = 0..exp(-Pi)=0.043213918263772249773 for z=0.

See the example in NIST[20] 20.14 "Methods of Computation" for theta3(0.9).

In practice I use the following transformations for q >= 0.1:

theta2(q) = sqrt(-Pi/ln(q)) * theta4(exp(Pi^2/ln(q))   [NIST 20.7.31]
theta3(q) = sqrt(-Pi/ln(q)) * theta3(exp(Pi^2/ln(q))   [NIST 20.7.32]
theta4(q) = sqrt(-Pi/ln(q)) * theta2(exp(Pi^2/ln(q))   [NIST 20.7.33]

Additionally for q < 0 the relation theta4(q) = theta3(-q) is applied.

[NIST 20.7.30] implies that theta1p(q) = d/dz(theta1(z,q)) at z=0 has the
functional equation theta1p(q) = (-Pi/ln(q))^(3/2) * theta1p(exp(Pi^2/ln(q)).
*)

{---------------------------------------------------------------------------}
function sfd_theta2(q: double): double;
  {-Return Jacobi theta_2(q) = 2*q^(1/4)*sum(q^(n*(n+1)),n=0..Inf), 0 <= q < 1}
var
  s,z: double;
begin
  if IsNanOrInfD(q) or (q < 0.0) or (q >= 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_theta2 := NaN_d;
    exit;
  end;
  if q<0.1 then sfd_theta2 := theta2small(q)
  else begin
    {theta2(q) = sqrt(-Pi/ln(q)) * theta3(-exp(Pi^2/ln(q))}
    z := ln(q);
    s := exp(PiSqr/z);
    s := theta3small(-s);
    sfd_theta2 := s*Sqrt(-Pi/z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_theta3(q: double): double;
  {-Return Jacobi theta_3(q) = 1 + 2*sum(q^(n*n)),n=1..Inf); |q| < 1}
var
  s,z: double;
begin
  if IsNanOrInfD(q) or (abs(q) >= 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_theta3 := NaN_d;
    exit;
  end;
  if q<0.0 then sfd_theta3 := sfd_theta4(-q)
  else if q<0.1 then sfd_theta3 := theta3small(q)
  else begin
    z := ln(q);
    s := exp(PiSqr/z);
    s := theta3small(s);
    sfd_theta3 := s*Sqrt(-Pi/z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_theta4(q: double): double;
  {-Return Jacobi theta_4(q) = 1 + 2*sum((-1)^n*q^(n*(n+1)),n=1..Inf); |q| < 1}
var
  s,z: double;
begin
  if IsNanOrInfD(q) or (abs(q) >= 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_theta4 := NaN_d;
    exit;
  end;
  if q<0.0 then sfd_theta4 := sfd_theta3(-q)
  else if q<0.1 then sfd_theta4 := theta3small(-q)
  else begin
    z := ln(q);
    if q>=0.7 then begin
      {Use only one term in theta2small and avoid underflow:}
      {s = exp(PiSqr/z) may be zero but s^0.25 is nonzero!}
      s := 2.0*exp(0.25*PiSqr/z);
    end
    else begin
      s := exp(PiSqr/z);
      s := theta2small(s);
    end;
    sfd_theta4 := s*Sqrt(-Pi/z);
  end;
end;


{---------------------------------------------------------------------------}
function theta1psmall(q: double): double;
  {-Return Jacobi theta1p(q) = 2*q^(1/4)*sum((-1)^n*(2n+1)*q^(n*(n+1)),n=0..Inf), q small}
var
  p,t,s,f,r: double;
  n: integer;
begin
  s := 0.0;
  if q<>0.0 then begin
    t := sqr(q);
    p := -1.0;
    f := 2.0*sqrt(sqrt(q));
    r := f;
    n := 1;
    while abs(r) >= eps_d*abs(s) do begin
      inc(n,2);
      s := s + r;
      p := p*t;
      f := f*p;
      r := f*n;
    end;
  end;
  theta1psmall := s;
end;


{---------------------------------------------------------------------------}
function sfd_theta1p(q: double): double;
  {-Return the derivative  theta1p(q) := d/dx(theta_1(x,q)) at x=0,}
  { = 2*q^(1/4)*sum((-1)^n*(2n+1)*q^(n*(n+1)),n=0..Inf), 0 <= q < 1}
var
  s,z: double;
begin
  if IsNanOrInfD(q) or (q < 0.0) or (q >= 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_theta1p := NaN_d;
    exit;
  end;
  if q<=0.1 then sfd_theta1p := theta1psmall(q)
  else begin
    {theta1p(q) = (-Pi/ln(q))^(3/2) * theta1p(exp(Pi^2/ln(q))}
    z := ln(q);
    if q>=0.7 then begin
      {Use only one term in theta1psmall and avoid underflow:}
      {s = exp(PiSqr/z) may be zero but s^0.25 is nonzero!}
      s := 2.0*exp(0.25*PiSqr/z);
    end
    else begin
      s := exp(PiSqr/z);
      s := theta1psmall(s);
    end;
    z := abs(Pi/z);
    sfd_theta1p := s*z*sqrt(z);
  end;
end;


(*
Theta functions for argument x and small q use the Fourier series
given in NIST[30], 20.2.2/3/4 and HMF[1] 16.27.2/3/4:

theta1(x,q) = 2*q^(1/4)*sum((-1)^k*q^(k*(k+1))*sin((2k+1)*x),k=0..inf)
theta2(x,q) = 2*q^(1/4)*sum(q^(k*(k+1))*cos((2k+1)*x),k=0..inf)
theta3(x,q) = 1 + 2*sum(q^(k^2)*cos(2*k*z),k=1..inf)
theta4(x,q) = 1 + 2*sum((-1)^k*q^(k^2)*cos(2*k*z),k=1..inf)
*)


{---------------------------------------------------------------------------}
function theta_fs(n: integer; x,q: double): double;
  {-Return theta_n(x,q) via Fourier series, |q| should be < about 0.25}
var
  ck,sk,c2,s2,s,fk,t,qk,q2: double;
begin
  sincos(2.0*x, s2, c2);
  if (n=1) or (n=2) then begin
    s  := 0.0;
    fk := 2.0*sqrt(sqrt(q));
    if n=1 then qk := -1.0 else qk := 1.0;
    sincos(x,sk,ck);
  end
  else begin
    {n=3 or n=4}
    s  := 1.0;
    fk := 2.0*q;
    qk := q;
    if n<>3 then begin
      fk := -fk;
      qk := -qk;
    end;
    ck := c2;
    sk := s2;
  end;
  q2 := sqr(q);
  while abs(fk) > eps_d*abs(s) do begin
    if n=1 then s := s + fk*sk else s := s + fk*ck;
    {Use addition theorems for sin/cos to get sin/cos((2k+1)*x)}
    t  := ck*c2 - sk*s2;
    sk := sk*c2 + ck*s2;
    ck := t;
    qk := qk*q2;
    fk := fk*qk;
  end;
  theta_fs := s;
end;


(*
With non-zero arguments the "Transformations of Lattice Parameter" from
NIST[30] 20.7 become more complicated, explicit formulas can be found at

http://functions.wolfram.com/EllipticFunctions/EllipticTheta1/06/02/0001/
Unique reference URL: http://functions.wolfram.com/09.01.06.0042.01

EllipticTheta[1, z, q] == (-((2 Sqrt[Pi])/Sqrt[-Log[q]])) E^((Pi^2 + 4 z^2)/(4 Log[q]))
  Sum[(-1)^k E^((k (k + 1) Pi^2)/Log[q]) Sinh[((2 k + 1) Pi z)/Log[q]], {k, 0, Infinity}]

and similar for theta2, theta3, theta4. With the definitions p=ln(q),
y=Pi^2/ln(q), w = Pi*abs(x/ln(q)) I get for t1 = theta1(x,q) etc

t1 = 2*sqrt(-Pi/p) * exp((x^2+Pi^2/4)/p) * sum((-1)^k*exp(k(k+1)y) * sinh((2k+1)w), k=0..inf)
t4 = 2*sqrt(-Pi/p) * exp((x^2+Pi^2/4)/p) * sum(       exp(k(k+1)y) * cosh((2k+1)w), k=0..inf)

t2 =   sqrt(-Pi/p) * exp(x^2/p) * [1 + 2*sum((-1)^k*exp(k^2*y) * cosh(2kw), k=1..inf)]
t3 =   sqrt(-Pi/p) * exp(x^2/p) * [1 + 2*sum(       exp(k^2*y) * cosh(2kw), k=1..inf)]
*)


{---------------------------------------------------------------------------}
function theta_ts(n: integer; x,q: double): double;
  {-Return theta_n(x,q) via transformed series, q should be > ~ 0.25, }
  { otherwise the 3-term approximations of the sums may be inaccurate!}
var
  p,w,s,t1,t2,y: double;
  j: integer;
begin
  {Reduce modulo Pi, |result| <= Pi/2, j is needed only for theta1/theta2}
  j := rem_pio2(0.5*x,w);
  x := 2.0*w;
  if (n=1) and (x=0.0) then begin
    theta_ts := 0.0;
    exit;
  end;
  p := ln(q);
  w := Pi*abs(x/p);
  if w>=ln_MaxDbl then begin
    {Hack: exp(w) will overflow, q is very close to 1. Set result=0 is OK for}
    {most (x,q) but not for theta_1((2k+1)*Pi/2,1) or theta_2((k*Pi,1), where}
    {the results are undefined. n=2, x=k*Pi should not arrive here for q < 1.}
    if (n<>1) or (abs(x)<>Pi_2) then begin
      theta_ts := 0.0;
      exit;
    end;
    {exp(w) will overflow}
  end;
  y := PiSqr/p;
  n := n and 3;
  if n=0 then begin
    {n=0: theta_4}
    if w > 19.0 then begin
      {Here cosh(k*z) = 0.5*exp(k*z). This is very handy because separately}
      {the cosh term(s) may overflow and/or the exp term(s) may underflow! }
      s  := exp(w);
      t1 := exp(2.0*y + 3.0*w);
      t2 := exp(6.0*y + 5.0*w);
      s  := 0.5*((t2+t1)+s);
    end
    else begin
      s  := cosh(w);
      t1 := exp(2.0*y)*cosh(3.0*w);
      t2 := exp(6.0*y)*cosh(5.0*w);
      s  := (t2+t1)+s;
    end;
    y := sqr(x)+sqr(Pi_2);
    y := 2.0*exp(y/p);
  end
  else if n=1 then begin
     {theta_1(x+j*Pi,q) = (-1)^j*theta_1(x,q), theta_1(-x,q) = -theta_1(x,q)}
     {Combine multiples of Pi from reduction and sign change for odd function}
     if odd(j) then j := -isign(x) else j := isign(x);
     if w > 19.0 then begin
       {here sinh(k*z) = 0.5*exp(k*z)}
       s  :=  exp(w);
       t1 := -exp(2.0*y + 3.0*w);
       t2 :=  exp(6.0*y + 5.0*w);
       s  :=  0.5*((t2+t1)+s);
     end
     else begin
       s  :=  sinh(w);
       t1 := -exp(2.0*y)*sinh(3.0*w);
       t2 :=  exp(6.0*y)*sinh(5.0*w);
       s  :=  (t2+t1)+s;
     end;
     y := sqr(x)+sqr(Pi_2);
     y := 2.0*j*exp(y/p);
   end
   else begin
     {n=2 or n=3}
     if w>19.0 then begin
       {here cosh(k*z) = 0.5*exp(k*z)}
       t1 := 0.5*exp(y + 2.0*w);
       t2 := 0.5*exp(4.0*y + 4.0*w);
     end
     else begin
       t1 := exp(y)*cosh(2.0*w);
       t2 := exp(4.0*y)*cosh(4.0*w);
     end;
     y := exp(sqr(x)/p);
     if n=2 then begin
       s := 1.0 + 2.0*(t2-t1);
       {theta_2(x+j*Pi,q) = (-1)^j*theta_2(x,q)}
       if odd(j) then s := -s;
     end
     else s := 1.0 + 2.0*(t2+t1);
  end;
  theta_ts := sqrt(-Pi/p)*s*y;
end;


{---------------------------------------------------------------------------}
function sfd_jtheta(n: integer; x,q: double): double;
  {-Return Jacobi theta_n(x,q), n=1..4, 0 <= q < 1}
begin
  if (n<1) or (n>4) or IsNanOrInfD(x) or IsNanOrInfD(q) or (q<0.0) or (q>=1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_jtheta := NaN_d;
    exit;
  end;
  if q<=0.25 then sfd_jtheta := theta_fs(n,x,q)
  else sfd_jtheta := theta_ts(n,x,q);
end;


{---------------------------------------------------------------------------}
function sfd_ntheta_s(x, k: double): double;
  {-Return the Neville theta_s function, |k| <= 1}
var
  f,t,q,kc: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (abs(k) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ntheta_s := NaN_d;
    exit;
  end;
  k := abs(k);
  if x=0.0 then sfd_ntheta_s := 0.0
  else if k=0.0 then sfd_ntheta_s := sin(x)
  else if k=1.0 then sfd_ntheta_s := sinh(x)
  else begin
    kc:= sqrt((1.0+k)*(1.0-k));
    f := Pi_2/sfd_elliptick(k);
    q := sfd_nome(k);
    t := sfd_jtheta(1,f*x, q);
    {http://functions.wolfram.com/09.12.27.0008.01}
    sfd_ntheta_s := sqrt(f/(k*kc)) * t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ntheta_c(x, k: double): double;
  {-Return the Neville theta_c function, |k| <= 1}
var
  f,t,q: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (abs(k) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ntheta_c := NaN_d;
    exit;
  end;
  k := abs(k);
  if (x=0.0) or (k=1.0) then sfd_ntheta_c := 1.0
  else if k=0.0 then sfd_ntheta_c := cos(x)
  else begin
    f := Pi_2/sfd_elliptick(k);
    q := sfd_nome(k);
    t := sfd_jtheta(2, f*x, q);
    {http://functions.wolfram.com/09.09.27.0008.01}
    sfd_ntheta_c := sqrt(f/k) * t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ntheta_d(x, k: double): double;
  {-Return the Neville theta_d function, |k| <= 1}
var
  f,t,q: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (abs(k) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ntheta_d := NaN_d;
    exit;
  end;
  k := abs(k);
  if (x=0.0) or (k=0.0) or (k=1.0) then sfd_ntheta_d := 1.0
  else begin
    f := Pi_2/sfd_elliptick(k);
    q := sfd_nome(k);
    t := sfd_jtheta(3, f*x, q);
    {http://functions.wolfram.com/09.10.27.0008.02}
    sfd_ntheta_d := sqrt(f) * t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ntheta_n(x, k: double): double;
  {-Return the Neville theta_n function, |k| <= 1}
var
  f,t,q,kc: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(k) or (abs(k) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ntheta_n := NaN_d;
    exit;
  end;
  k := abs(k);
  if (x=0.0) or (k=0.0) then sfd_ntheta_n := 1.0
  else if k=1.0 then sfd_ntheta_n := cosh(x)
  else begin
    kc:= sqrt((1.0+k)*(1.0-k));
    f := Pi_2/sfd_elliptick(k);
    q := sfd_nome(k);
    t := sfd_jtheta(4, f*x, q);
    {http://functions.wolfram.com/09.11.27.0008.01}
    sfd_ntheta_n := sqrt(f/kc) * t;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ellmod(q: double): double;
  {-Return the elliptic modulus k(q) = theta_2(q)^2/theta_3(q)^2, 0 <= q <= 1}
var
  s,t: double;
const
  ns = 15;
  csem: array[0..ns-1] of double = (
          +1.625891465772601185014117987322,
          -0.168800247494414500563963109226,
          +0.167783989069592230112335970242e-1,
          -0.137029033311059238794404650790e-2,
          +0.983250139150171749424944050175e-4,
          -0.656505092885482142098957147599e-5,
          +0.414016166743303940386659406466e-6,
          -0.247869953553996473283396161694e-7,
          +0.142708842202339217388709463169e-8,
          -0.795758836092878516233596831460e-10,
          +0.430534892144426740027174229179e-11,
          -0.227115370335370231275267112278e-12,
          +0.117265335584543762807690016967e-13,
          -0.593363805868622919218264061794e-15,
          +0.294920599498337447452493071359e-16);
         {-0.144293466348869681410375881412e-17,
          +0.695591732809755059231957434888e-19);
          -0.330828010159576563873746101955e-20}
begin
  if IsNanOrInfD(q) or (q < 0.0) or (q > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ellmod := NaN_d;
    exit;
  end;
  if q >= 0.785 then begin
    {Result is 1.0 accurate to double precision}
    sfd_ellmod := 1.0;
  end
  else if q <= 0.125 then begin
    {For small q use a Chebyshev approximation of k(q) calculated with Maple}
    {using the truncated series from theta_2(q)^2/theta_3(q)^2/(4*sqrt(q)).}
    t := 4.0*sqrt(q);  {Error if q < 0!}
    if q <= 0.8e-10 then s := 1.0 - 4.0*q
    else s := CSEvalD(16.0*q - 1.0, csem, ns);
    sfd_ellmod := t*s;
  end
  else begin
    {Just the definition, see e.g. NIST[30], 22.2.2}
    s := sfd_theta2(q);
    t := sfd_theta3(q);
    sfd_ellmod := sqr(s/t);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_nome(k: double): double;
  {-Return the elliptic nome q(k) = exp(-Pi*EllipticCK(k)/EllipticK(k)), |k| < 1}
var
  a,b,t,r: double;
const
  k0 = 5e-5;
  k1 = 0.125;
  k2 = 0.99999999;
const
  nsn = 12;
  nshx: array[0..nsn-1] of THexDblW = (
          ($6312,$308D,$0C17,$4000),    {2.0059036057613265065256485799798524633812340147058     }
          ($0C9B,$FC51,$2520,$3F70),    {0.39416588306863359080191352399431774016597800766483e-2 }
          ($24DA,$CA58,$4A72,$3F50),    {0.99431238175134174582632234399239863694397835877573e-3 }
          ($04EE,$0F7F,$6C37,$3ED5),    {0.51075733082394804367489466126719473327615646197487e-5 }
          ($A007,$437E,$0546,$3EA6),    {0.65626512134754001464418435164938712614061717607757e-6 }
          ($3948,$23BA,$1AB5,$3E38),    {0.56122258469326131324759612202936891012813215632088e-8 }
          ($2E14,$2429,$FACA,$3E00),    {0.49417276548762464038941094464580963037820551299246e-9 }
          ($1338,$9349,$E267,$3D99),    {0.58854295580591350466762622102835403773353282065283e-11}
          ($5C6A,$5BDA,$5974,$3D5C),    {0.40286965577813384289291737568191299745169616595575e-12}
          ($3D6B,$A5F8,$527C,$3CFB),    {0.60667503204232521105790175507714700526663813002074e-14}
          ($D9A9,$A594,$FE2E,$3CB8),    {0.34684615282276839064836958110435873701557771216851e-15}
          ($EBDC,$E457,$A7B6,$3C5C));   {0.62135918068660567534275568130140096659280495938637e-17}
(*        ($4537,$CCC6,$F093,$3C16),    {0.31089165594535819628922137562248497069837738126753e-18}
          ($F736,$B234,$00BB,$3BBE),    {0.63533535832825954506732958344711833425042943145442e-20} *)
var
  ensa: array[0..nsn-1] of double absolute nshx;
begin
  if IsNanOrInfD(k) or (abs(k) >= 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_nome := NaN_d;
    exit;
  end;
  k := abs(k);
  {For accuracy tests near k=1 it is essential that the multi-precision}
  {calculation routine uses the stored (double) Pascal k values, e.g.}
  {k2 is stored as 0.999999990000000000005177941408973651959968265146!}
  if k<=k1 then begin
    {Use NIST[30], 19.5.5: With r = k^2/16 and q = exp(-Pi*K'(k)/K(k))}
    {q = r + 8*r^2 + 84*r^3 + 992*r^4. More terms using m=k^2 are given at}
    {http://functions.wolfram.com/EllipticFunctions/EllipticNomeQ/06/01/02/0001/}
    {Unique reference URL: http://functions.wolfram.com/09.53.06.0001.01}
    r := sqr(0.25*k);
    if k<=k0 then begin
      {Use only the first two terms, i.e. 1 + 84r^2 = 1}
      sfd_nome := r*(1.0+8.0*r)
    end
    else begin
      {Writing q(r)=r*f(r), I use a Chebyshev expansion calculated with Maple}
      {for f(r) with k<=1/8. Using the Wolfram m values, f is the polynomial:}
      {f(r) = 1 + 8*r + 84*r^2 + 992*r^3 + 12514*r^4 + 164688*r^5            }
      {         + 2232200*r^6 + 30920128*r^7 + 435506703*r^8                 }
      {         + 6215660600*r^9 + 89668182220*r^10 + 1305109502496*r^11     }
      {         + 19138260194422*r^12 + 282441672732656*r^13                 }
      {         + 4191287776164504*r^14 + 62496081197436736*r^15             }
      {         + 935823746406530603*r^16 + 14065763582458332888*r^17        }
      {         + 212122153814497767004*r^18 + 3208590886304243284640*r^19   }
      {The maximum relative error for r*f(r) is less than 4e-39 for k <= 1/8.}
      t := CSEvalD(16.0*k-1,ensa,nsn);
      sfd_nome := r*t;
    end;
  end
  else if k>=k2 then begin
    {Use parts of the following formula with m=k^2 and t=1-m = (1-k)*(1+k) from
     http://functions.wolfram.com/EllipticFunctions/EllipticNomeQ/06/01/03/0001/
     Unique reference URL: http://functions.wolfram.com/09.53.06.0003.01
     EllipticNomeQ[m] == Exp[Pi^2/Log[(1 - m)/16]]  *
                           (1 + (1/2) (m - 1) Pi^2 (1/Log[(1 - m)/16]^2)
                              - (13/64) (m - 1)^2 Pi^2 (1/Log[(1 - m)/16]^2)
                              + (23/192) (m - 1)^3 Pi^2 (1/Log[(1 - m)/16]^2)
                              + (1/4) (m - 1)^2 Pi^2 (1/Log[(1 - m)/16]^3)
                              - (13/64) (m - 1)^3 Pi^2 (1/Log[(1 - m)/16]^3) + ...)
    }
    t := (1.0-k)*(1.0+k);  {Do NOT change to 1.0-k*k !!}
    r := ln(0.0625*t);
    a := exp(PiSqr/r);
    b := PiSqr*t/sqr(r);
    sfd_nome := a*(1.0 - 0.5*b*(1.0 + (13/32)*t));  {Fix311}
  end
  else begin
    {get K'(k)/K(k) via two inline AGM iterations}
    {compute r = 2*agm(1,k)}
    a := 1.0;
    b := k;
    while a-b > 0.5e-8*a do begin
      t := a;
      a := 0.5*(t+b);
      b := sqrt(t*b);
    end;
    r := a+b;
    {compute t = 2*agm(1,k')}
    a := 1.0;
    b := sqrt((1.0-k)*(1.0+k));
    while a-b > 0.5e-8*a do begin
      t := a;
      a := 0.5*(t+b);
      b := sqrt(t*b);
    end;
    t := (a+b);
    {q = exp(-Pi*K'(k)/K(k)) = exp(-Pi*agm(1,k')/agm(1,k))}
    sfd_nome := exp(-(t/r)*Pi);
  end;
end;


{---------------------------------------------------------------------------}
procedure sfd_lemn(x: double; var sl,cl: double);
  {-Return the lemniscate functions sl = sin_lemn(x), cl = cos_lemn(x)}
const
  {P1+P2+P3 = 2*omega = 5.244115108584239621: period of lemniscate functions}
  HP1 : THexDblW = ($0000,$0000,$F9F9,$4014); {2749426/2^19}
  HP2 : THexDblW = ($0000,$0000,$E7D4,$3EB3); { 326133/2^38}
  P3  = 1.567026880591160005e-13;
  TOM = 5.244115108584239621;
  x0  = 1.01e-4;
var
  P1 : double absolute HP1; {= 5.244113922119140625  }
  P2 : double absolute HP2; {= 1.18646494229324162e-6}
var
  y: double;
begin
  if IsNanOrInfD(x) then begin
    sl := Nan_d;
    cl := Nan_d;
    exit;
  end;
  y := abs(x);
  if y<=x0 then begin
    sl := x;
    cl := 1.0-x*x;
  end
  else begin
    if y > TOM then begin
      {Cody/Waite style reduction mod 2*omega}
      y := rint(x/TOM);
      x := ((x - y * P1) - y * P2) - y * P3;
    end;
    sfd_sncndn(x*sqrt2,0.5,sl,cl,y);
    sl := sl/y/sqrt2;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_sin_lemn(x: double): double;
  {-Return the lemniscate sine function sin_lemn(x)}
var
  sl,cl: double;
begin
  sfd_lemn(x,sl,cl);
  sfd_sin_lemn := sl;
end;


{---------------------------------------------------------------------------}
function sfd_cos_lemn(x: double): double;
  {-Return the lemniscate cosine function cos_lemn(x)}
var
  sl,cl: double;
begin
  sfd_lemn(x,sl,cl);
  sfd_cos_lemn := cl;
end;


{---------------------------------------------------------------------------}
function sfd_arccl(x: double): double;
  {-Return the inverse lemniscate cosine function, |x| <= 1}
const
  Sqrt2h = 0.7071067811865475244;      {sqrt(2)/2}
  Ki     = 1.0+0.3110287771460599052;  {K(i)}
begin
  if IsNanOrInfD(x) or (abs(x) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_arccl := NaN_s;
  end
  else if x=0.0 then sfd_arccl := Ki
  else if x=-1.0 then sfd_arccl := 2*Ki
  else sfd_arccl := sfd_jacobi_arccn(x, Sqrt2h) / Sqrt2;
end;


{---------------------------------------------------------------------------}
function sfd_arcsl(x: double): double;
  {-Return the inverse lemniscate sine function, |x| <= 1}
const
  Sqrt2h = 0.7071067811865475244;      {sqrt(2)/2}
  Ki     = 1.0+0.3110287771460599052;  {K(i)}
begin
  if IsNanOrInfD(x) or (abs(x) > 1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_arcsl := NaN_d;
  end
  else if abs(x)=1.0 then sfd_arcsl := copysignd(ki, x)
  else sfd_arcsl := sfd_jacobi_arcsd(x*Sqrt2, Sqrt2h) / Sqrt2;
end;


{---------------------------------------------------------------------------}
{--------------- Weierstrass elliptic and modular functions ----------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function sfd_detai(x: double): double;
  {-Return Dedekind eta(i*x), x >= 0}
var
  q,t: double;
const
  x0 = 2e-5;
  x1 = 0.125;
  x2 = 8;
  x3 = 3000;
  pi12 = -0.25 -0.1179938779914943654e-1; {-Pi/12}
begin
  if IsNanOrInfD(x) or (x<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_detai := NaN_d;
  end
  else if (x <= x0) or (x >= x3) then begin
    sfd_detai := 0;
  end
  else if x <= x1 then begin
    t := exp(-Pi/(12.0*x));
    sfd_detai := t/sqrt(x);
  end
  else if x >= x2 then begin
    sfd_detai := exp(pi12*x);
  end
  else begin
    q := exp(-TwoPi*x);
    t := sfd_eulerq(q);
    sfd_detai := t*exp(pi12*x);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_emlambda(y: double): double;
  {-Return the elliptic modular function lambda(iy), y >= 0}
var
  q: double;
const
  y1 = 1/15;
  y2 = 12.5;
  yl = 240;
begin
  {tau = iy, q = exp(i*Pi*tau) = exp(-Pi*y)}
  if y <= y1 then sfd_emlambda := 1.0
  else if y >= yl then sfd_emlambda := 0.0
  else begin
    q := exp(-Pi*y);
    if y > y2 then begin
      {OEIS A029845: 16/Lambda = 1/q + 8 - 20*q + O(q^3)}
      sfd_emlambda := 16.0*q;
    end
    else begin
      {NIST[30], 23.15.2/6}
      sfd_emlambda := sqr(sfd_ellmod(q));
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_KleinJ(y: double): double;
  {-Return Klein's complete invariant J(iy), y>0}
var
  q,t1,t2,t3,t4: double;
begin
  {NIST 23.15.7}
  if y=1.0 then sfd_kleinj := 1.0
  else if y>0.0 then begin
    if y>1.0 then q := exp(-Pi*y)
    else q := exp(-Pi/y);
    t1 := 54.0*power(sfd_theta1p(q),8);
    t2 := power(sfd_theta2(q),8);
    t3 := power(sfd_theta3(q),8);
    t4 := power(sfd_theta4(q),8);
    sfd_kleinj := power(t2+t3+t4,3)/t1;
  end
  else sfd_kleinj := PosInf_d;
end;


(*  Internal development history

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.00.01  10.04.10  W.Ehrhardt  First version of wpe
 0.00.02  11.04.10  we          wp
 0.00.03  14.04.10  we          inline sorting and root polishing in wpe
 0.00.04  14.04.10  we          rename wp to wp0, new wp with homogenity
 0.00.05  14.04.10  we          wpg (calls wpe or wp)
 0.00.06  15.04.10  we          wpe_inv
 0.00.07  20.04.10  we          wpeim, wpgim
 0.00.08  20.04.10  we          wpl (basic lemniscate case) with Cosy/Waite
 0.00.09  21.04.10  we          wpe: special case e1=e3
 0.00.10  21.04.10  we          wpg: root polishing via conditional define
 0.00.11  25.04.10  we          wpg/wpe: use wpl if g3=0, g2>0
*)


{---------------------------------------------------------------------------}
function lemnreduce(x: double): double;
  {-Argument reduction for basic lemniscatic case}
const
  HP1 : THexDblW = ($0000,$0000,$AA4A,$400D); {1944138/2^19}
  HP2 : THexDblW = ($0000,$0000,$BAC0,$3E9A); { 109484/2^38}
  P3  = 0.3443996238599894830e-11;
  TOM = 3.0+0.7081493546027438367;            {Gamma(1/4)^2/2/sqrt(Pi)}
var
  P1 : double absolute HP1; {= 3.70814895629883}
  P2 : double absolute HP2; {= 3.983004717156297-7}
var
  t: double;
begin
  if abs(x) > TOM/2 then begin
    {Cody/Waite style reduction}
    t := rint(x/TOM);
    x := ((x - t * P1) - t * P2) - t * P3;
  end;
  lemnreduce := x;
end;



{---------------------------------------------------------------------------}
function sfd_wpl(x: double): double;
  {-Return the Weierstrass function wp(x,1,0)=wpe(x,1/2,0), basic lemniscatic case}
var
  s,c,d: double;
begin
  x := lemnreduce(x);
  sfd_sncndn(abs(x),0.5,s,c,d);
  sfd_wpl := sqr(d/s);
end;


{---------------------------------------------------------------------------}
function sfd_wpe(x,e1,e2: double): double;
  {-Return the Weierstrass function P(x,e1,e2) from the lattice roots e1 < e2}
var
  e3,mc,s,t: double;
var
  sn,cn,dn: double;
begin
  e3 := -e1-e2;

  {sort(e1,e2,e3);}
  if e1 < e2 then begin t:=e1; e1:=e2; e2:=t; end;
  if e2 < e3 then begin t:=e2; e2:=e3; e3:=t; end;
  if e1 < e2 then begin t:=e1; e1:=e2; e2:=t; end;

  t  := e1-e3;
  if t<=0.0 then begin
    {e1=e2=e3=0}
    sfd_wpe := sqr(1.0/x);
    exit;
  end;
  if (e1=0.0) or (e2=0.0) or (e3=0) then begin
    {here g3=0, g2 >0, use inline lemniscating case}
    t := sqrt(2.0*(e1*e1 + e2*e2 + e3*e3));
    sfd_wpe := t*sfd_wpl(sqrt(t)*x);
    exit;
  end;
  s  := sqrt(t);
  mc := (e1-e2)/t;
  sfd_sncndn(s*x,mc,sn,cn,dn);
  {NIST 23.6 17,21,22,23, use smallest e_k > 0}
  if e3 >= 0.0 then begin
    t := s/sn;
    sfd_wpe := e3 + sqr(t);
  end
  else if e2 >= 0.0 then begin
    t := s*(dn/sn);
    sfd_wpe := e2 + sqr(t);
  end
  else begin
    t := s*(cn/sn);
    sfd_wpe := e1 + sqr(t);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_wpe_im(y,e1,e2: double): double;
  {-Return the Weierstrass function P(iy,e1,e2) from the lattice roots e1 < e2}
begin
  sfd_wpe_im := -sfd_wpe(y,-e1,-e2);
end;


{---------------------------------------------------------------------------}
function sfd_wpe_inv(y,e1,e2: double): double;
  {-Return the smallest positive x with wpe(x)=y, y >= e1}
begin
  {NIST 23.6.30}
  sfd_wpe_inv := sfd_ell_rf(y - e1, y - e2, y + (e1+e2));
end;


{---------------------------------------------------------------------------}
procedure CompAuxWP(g2,g3: double; var delta,e1,e2,h2: double);
  {-Compute auxiliary values from g2,g3}
var
  q,r,t,a,b: double;
begin
  delta := g2*g2*g2 - 27.0*g3*g3;
  if delta > 0 then begin
    {solve f(e) = 4e^3 - g2*e - g3 = 0 with positive discriminant}
    {compute e1,e2,e3 from g2,g3, NR 5.6.9 .. 5.6.12}
    {a=0, b=-g2/4, c=-g3/4}
    r  := -g3/8.0;
    q  := g2/12.0;
    t  := arccos(r/sqrt(q*q*q));
    r  := -2.0*sqrt(q);
    e1 := r*cos(t/THREE);
    e2 := r*cos((t+TwoPi)/THREE);
  end
  else begin
    {negative discriminant}
    r := -g3/8.0;
    q := g2/12.0;
    t := sqrt(r*r - q*q*q);
    {Compute A,B from NR 5.6.15/16}
    a := cbrt(abs(r)+t);
    if r>0.0 then a := -a;
    if a=0.0 then b := 0.0 else b := q/a;
    {Real root e2}
    e2 := (a+b);
    {Complex roots e1,e3 = (A+B)/2 +- i/2*sqrt(3)(A-B)}
    {H2 = sqrt((e2-e1)*(e2-e3)), HMF[1],18.3.5/6}
    h2 := sqrt(3.0*(a*a + a*b + b*b));
  end;
end;


{---------------------------------------------------------------------------}
function equireduce(x: double): double;
  {-Argument reduction for equianharmonic case}
const
  HP1 : THexDblW = ($0000,$0000,$7AB1,$4008); {1604273/2^19}
  HP2 : THexDblW = ($0000,$0000,$98A0,$3E85); {  44229/2^38}
  P3  = 0.2213278467016205641e-11;
  TOM = 3.0 + 0.05990807411438574983;         {Gamma(1/3)^3/2/Pi}
var
  P1 : double absolute HP1; {= 3.0599079132080078125}
  P2 : double absolute HP2; {= 0.160904164658859372138977050781e-6}
var
  t: double;
begin
  if abs(x) > TOM/2 then begin
    {Cody/Waite style reduction}
    t := rint(x/TOM);
    x := ((x - t * P1) - t * P2) - t * P3;
  end;
  equireduce := x;
end;


{---------------------------------------------------------------------------}
function wpdup(x,g2,g3: double): double;
  {-Return the Weierstrass function P(x,g2,g3) via duplication}
var
  y,z,t,v,p: double;
  k,n: integer;
const
  xlow = 1/256;
begin
  z := abs(x);
  {Ref: R. Coquereaux et al [75]}
  n := 0;
  while z >= xlow do begin
    z := 0.5*z;
    inc(n);
  end;
  y := z*z;
  p := (g3/28.0*y + g2/20.0)*y + 1.0/y;
  {p := (((3.0/6160.0*g2*g3*y + g2*g2/1200)*y + g3/28.0)*y + g2/20.0)*y + 1.0/y;}
  for k:=1 to n do begin
    y := p*p;
    t := sqr(6.0*y - 0.5*g2);
    v := 4.0*((4.0*y -g2)*p -g3);
    p := t/v - 2.0*p;
  end;
  wpdup := p;
end;


{---------------------------------------------------------------------------}
function sfd_wpg(x,g2,g3: double): double;
  {-Return the Weierstrass function P(x,e1,e2) from lattice invariants g2, g3}
var
  delta,e1,e2,h2,z,mc,sn,cn,dn: double;
begin
  if (g2=0.0) and (g3=0.0) then sfd_wpg := sqr(1.0/x)
  else if (g2=1.0) and (g3=0.0) then sfd_wpg := sfd_wpl(x)
  else begin
    CompAuxWP(g2,g3, delta, e1,e2,h2);
    if delta > 0.0 then sfd_wpg:= sfd_wpe(x,e1,e2)
    else begin
      if (g2=0.0) and (g3=1.0) then x := equireduce(x);
      if abs(x) < 0.25 then sfd_wpg := wpdup(x,g2,g3)
      else begin
        {HMF[1], 18.9.9-11}
        z  := 2.0*x*sqrt(h2);
        mc := 0.75*e2/h2 + 0.5;
        sfd_sncndn(z,mc, sn,cn,dn);
        sfd_wpg := e2 + h2*(1.0+cn)/(1.0-cn);
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_wpg_im(y,g2,g3: double): double;
  {-Return the Weierstrass function P(iy, g2, g3)}
begin
  sfd_wpg_im := -sfd_wpg(y,g2,-g3);
end;


{---------------------------------------------------------------------------}
function sfd_wpg_inv(y,g2,g3: double): double;
  {-Return the smallest positive x with wpg(x,g2,g3)=y, y >= e2}
var
  delta,z,h2,e1,e2,k,t: double;
begin
  sfd_wpg_inv := Nan_d;
  if (g2=0.0) and (g3=0.0) then sfd_wpg_inv := 1.0/sqrt(y)
  else if (g2=1.0) and (g3=0.0) then sfd_wpg_inv := sfd_wpe_inv(y,0.5,0)
  else begin
    CompAuxWP(g2,g3, delta, e1,e2,h2);
    if delta > 0.0 then begin
      sfd_wpg_inv := sfd_wpe_inv(y,e1,e2)
    end
    else begin
      t  := 0.5 - 0.75*e2/h2;
      if t < 0.0 then exit;  {Should not happen}
      k := sqrt(t);
      t := (y-e2)/h2;
      z := t+1.0;
      if z=0.0 then exit;
      z := (t-1.0)/z;
      if abs(z) <= 1.0 then begin
        t := sfd_jacobi_arccn(z, k);
        sfd_wpg_inv := 0.5*t/sqrt(h2);
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_wpe_der(x,e1,e2: double): double;
  {-Return Weierstrass P'(x,e1,e2) from the lattice roots e1 < e2}
var
  e3,mc,s,t: double;
var
  sn,cn,dn: double;
begin
  e3 := -e1-e2;

  {sort(e1,e2,e3);}
  if e1 < e2 then begin t:=e1; e1:=e2; e2:=t; end;
  if e2 < e3 then begin t:=e2; e2:=e3; e3:=t; end;
  if e1 < e2 then begin t:=e1; e1:=e2; e2:=t; end;

  if (e1=0.5) and (e2=0.0) then x := lemnreduce(x);
  t  := e1-e3;
  if t<=0.0 then begin
    {e1=e2=e3=0}
    sfd_wpe_der := -2.0/(x*x*x);
    exit;
  end;
  s  := sqrt(t);
  mc := (e1-e2)/t;
  sfd_sncndn(s*x,mc,sn,cn,dn);
  {HMF[1], 18.9.12}
  sfd_wpe_der := -2.0*s*t*cn*dn/(sn*sn*sn);
end;


{---------------------------------------------------------------------------}
function wpgd_small(x,g2,g3: double): double;
  {-Weierstrass P' via Laurent series for small x}
var
  z,t: double;
const
  c1 = 5/5096;
  c2 = 3/770;
begin
  { P'(x,g2,g3) = -2*x^(-3) + 1/10*g2*x + 1/7*g3*x^3 + 1/200*g2^2*x^5 }
  {     + 3/770*g2*g3*x^7 + (1/15600*g2^3 + 5/5096*g3^2)*x^9 +O(x^11) }
  z := x*x;
  t := g2*g2*g3*z/15400.0 + g2*g2*g2/15600.0 + c1*sqr(g3);
  t := (((t*z + c2*g2*g3)*z + g2*g2/200.0)*z + g3/7.0)*z + g2/10.0;
  wpgd_small := -2.0/(z*x) + t*x;
end;


{---------------------------------------------------------------------------}
function sfd_wpg_der(x,g2,g3: double): double;
  {-Return Weierstrass P'(x,e1,e2) from lattice invariants g2, g3}
var
  delta,e1,e2,h2,z,mc,sn,cn,dn: double;
begin
  if (g2=0.0) and (g3=0.0) then sfd_wpg_der := -2.0/(x*x*x)
  else if (g2=1.0) and (g3=0.0) then sfd_wpg_der := sfd_wpe_der(x,0.5,0)
  else begin
    CompAuxWP(g2,g3, delta, e1,e2,h2);
    if delta > 0.0 then sfd_wpg_der := sfd_wpe_der(x,e1,e2)
    else begin
      if (g2=0.0) and (g3=1.0) then x := equireduce(x);
      if abs(x) <= 1/16 then sfd_wpg_der := wpgd_small(x,g2,g3)
      else begin
        {HMF[1], 18.9.9-12}
        z  := 2.0*x*sqrt(h2);
        mc := 0.75*e2/h2 + 0.5;
        sfd_sncndn(z,mc, sn,cn,dn);
        sfd_wpg_der := -4.0*h2*sqrt(h2)*sn*dn/sqr(1.0-cn);
      end;
    end;
  end;
end;


end.
