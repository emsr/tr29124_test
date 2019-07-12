unit sdZeta;

{Double precision Zeta functions and polylogarithms}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Double precision Zeta functions and polylogarithms

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in damath_info.txt/references

                  [7] Cephes Mathematical Library, Version 2.8
                      http://www.moshier.net/#Cephes or http://www.netlib.org/cephes/
                 [14] SLATEC Common Mathematical Library, Version 4.1, July 1993
                      (general purpose mathematical and statistical routines written in Fortran 77)
                      http://www.netlib.org/slatec
                 [19] Boost C++ Libraries, Release 1.42.0, 2010.
                      http://www.boost.org/
                 [20] Special functions by Wayne Fullerton,
                      http://www.netlib.org/fn
                      Almost identical to the FNLIB subset of SLATEC [14]
                 [21] GNU Scientific Library, GSL-1.14 (March 2010),
                      http://www.gnu.org/software/gsl
                 [22] A.J. MacLeod, MISCFUN: A software package to compute uncommon special functions.
                      ACM Trans. on Math. Soft. 22 (1996), pp.288-301.
                      Fortran source: http://netlib.org/toms/757
                 [30] [NIST]: F.W.J. Olver, D.W. Lozier, R.F. Boisvert, C.W. Clark, NIST Handbook
                      of Mathematical Functions, Cambridge, 2010. Online resource: NIST Digital
                      Library of Mathematical Functions, http://dlmf.nist.gov/
                 [31] R.E. Crandall, Note on fast polylogarithm computation, 2006.
                      http://www.reed.edu/~crandall/papers/Polylog.pdf
                 [33] http://functions.wolfram.com/: Formulas and graphics about
                      mathematical functions for the mathematical and scientific
                      community and/or http://mathworld.wolfram.com/ ("/the web's
                      most extensive mathematical resource/")
                 [36] P. Borwein, An Efficient Algorithm for the Riemann Zeta Function,
                      CMS Conference Proc. 27 (2000), pp. 29-34. Available as
                      http://www.cecm.sfu.ca/personal/pborwein/PAPERS/P155.pdf
                 [40] S.V. Aksenov et al., Application of the combined nonlinear-condensation
                      transformation to problems in statistical analysis and theoretical physics.
                      Computer Physics Communications, 150, 1-20, 2003.
                      Available from http://dx.doi.org/10.1016/S0010-4655(02)00627-6
                      or as e-print: https://arxiv.org/pdf/math/0207086v1
                 [41] C. Ferreira, J.L. Lopez, Asymptotic expansions of the Hurwitz-Lerch
                      zeta function, J. Math. Anal. Appl. 298 (2004), 210-224. Available from
                      http://dx.doi.org/10.1016/j.jmaa.2004.05.040
                 [50] A. Erdelyi et al., Higher Transcendental Functions Vol. I-III, California
                      Institute of Technology - Bateman Manuscript Project, 1953-1955,
                      Available via http://en.wikipedia.org/wiki/Bateman_Manuscript_Project
                 [73] D.C. Wood, The Computation of Polylogarithms, 1992, Technical
                      Report 15-92, http://www.cs.kent.ac.uk/pubs/1992/110


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.00.00  08.02.13  W.Ehrhardt  Initial BP7 version from AMath.sfzeta
 1.00.01  08.02.13  we          sfd_dilog, sfd_etaint, sfd_zetaint
 1.00.02  09.02.13  we          zetap with 53-bit logic from [19]
 1.00.03  09.02.13  we          etam1pos
 1.00.04  09.02.13  we          sfd_pz
 1.00.04  09.02.13  we          sfd_cl2
 1.00.05  12.02.13  we          sfd_fdm05, sfd_fdp05
 1.00.06  16.02.13  we          improved sfd_zetah
 1.00.07  01.03.13  we          Chebyshev degrees reduced in sfd_ti2

 1.01.00  27.03.13  we          sfd_zetah for s<0
 1.01.01  28.03.13  we          sfd_trilog
 1.01.02  30.03.13  we          sfd_zetah for a near 1 and s<0

 1.02.00  11.04.13  we          sfd_fdp15

 1.05.00  16.08.13  we          sfd_zeta for small s
 1.05.01  16.08.13  we          sfd_polylog with sfd_zetaint

 1.06.00  07.09.13  we          Improved sfd_zetah/hurwitz_formula
 1.06.01  12.09.13  we          Bernoulli polynomials sfd_bernpoly
 1.06.02  25.09.13  we          use const one_d

 1.08.00  26.12.13  we          sfd_fdp25

 1.10.00  02.05.14  we          sfd_harmonic
 1.10.01  03.05.14  we          sfd_harmonic2

 1.11.00  29.05.14  we          Removed redundancy in sfc_harmonic2
 1.11.01  29.05.14  we          polylogneg with array of single
 1.11.02  31.05.14  we          sfd_llci, sfd_llsi, improved harm2core

 1.12.00  16.06.14  we          Fermi/Dirac, Lobachewsky, harmonic functions move to sdZeta2
 1.12.01  24.06.14  we          Avoid Delphi-64 hint in lphi_aj

 1.13.00  30.07.14  we          Improved and expanded primezeta sfd_pz

 1.18.00  18.05.15  we          Improved sfd_etam1 (adjusted Borwein constants)
 1.18.01  18.05.15  we          sfd_polylogr and sfd_lerch for s >= -1
 1.18.02  25.05.15  we          sfd_lerch(1,s,a) for s<>1 and special case z=0

 1.19.00  19.07.15  we          eta near s=1

 1.20.00  25.08.15  we          Use THREE & SIXX to avoid 'optimization'
 1.20.01  25.08.15  we          sfd_bernpoly: avoid integer overflow and FPC optimization error

 1.28.00  02.12.17  we          Suppress warnings: Local variable does not seem to be initialized

 1.29.00  04.02.18  we          Rewrite polylog, returns real part for x>1,n>1
 1.29.01  04.02.18  we          sfd_zetaint(1) = NaN
 1.29.02  05.02.18  we          sfd_polylog(1,x) for x > 1
 1.29.03  07.02.18  we          sfd_polylogr for x < -1
 1.29.04  09.02.18  we          Merged with (old) unit sdZeta2
 1.29.05  09.02.18  we          sfd_fdr (Fermi-Dirac real order)
 1.29.06  12.02.18  we          sfd_polylogr(s,x)=x for large s

 1.30.00  24.02.18  we          sfd_lerch for z < -1
 1.30.01  25.02.18  we          fix NaN check for r in sfd_harmonic2
 1.30.02  26.02.18  we          fix overflow for LerchPhi integrand
 1.30.03  26.02.18  we          sfd_ti (inverse tangent integral of real order)
 1.30.04  01.03.18  we          improved lerch integration
 1.30.05  04.03.18  we          sfd_polylogr for 1 < x <= 256
 1.30.06  05.03.18  we          sfd_lchi for |x| > 1
 1.30.07  05.03.18  we          Bose-Einstein integral sfd_beint
 1.30.08  06.03.18  we          sfd_polylogr x > 256

 1.35.00  23.08.18  we          max iteration check in liae_pos
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


function sfd_cl2(x: double): double;
  {-Return the Clausen function: integral(-ln(2*|sin(t/2)|),t=0..x) = Im(Li_2(exp(ix)))}

function sfd_ti2(x: double): double;
  {-Return the inverse tangent integral, ti2(x) = integral(arctan(t)/t, t=0..x)}

function sfd_ti(s,x: double): double;
  {-Return the inverse tangent integral of order s >= 0}

function sfd_dilog(x: double): double;
  {-Return dilog(x) = Re(Li_2(x)), Li_2(x) = -integral(ln(1-t)/t, t=0..x)}

function sfd_trilog(x: double): double;
  {-Return the trilogarithm function trilog(x) = Re(Li_3(x))}

function sfd_polylog(n: integer; x: double): double;
  {-Return the polylogarithm Li_n(x) of integer order; real part for n>0,x>1}

function sfd_polylogr(s, x: double): double;
  {-Return the polylogarithm Li_s(x) of real order s >= -1, x <= 256;}
  { s > 0 if x > 256; real part if x > 1}

function sfd_pz(x: double): double;
  {-Return the prime zeta function P(x) = sum(1/p^x, p prime), x > 1/5; }
  { for x<1 the real part of P(x) is returned.}

function sfd_eta(s: double): double;
  {-Return the Dirichlet eta function}

function sfd_etaint(n: integer): double;
  {-Return the Dirichlet function eta(n) for integer arguments}

function sfd_etam1(s: double): double;
  {-Return Dirichlet eta(s)-1}

function sfd_lerch(z,s,a: double): double;
  {-Return the Lerch transcendent Phi(z,s,a) = sum(z^n/(n+a)^s, n=0..INF), z <= 1, s >= -1, a >= 0}

function sfd_zeta(s: double): double;
  {-Return the Riemann zeta function at s, s<>1}

function sfd_zeta1p(x: double): double;
  {-Return the Riemann zeta function at 1+x, x<>0}

function sfd_zetam1(s: double): double;
  {-Return Riemann zeta(s)-1, s<>1}

function sfd_zetah(s,a: double): double;
  {-Return the Hurwitz zeta function zetah(s,a) = sum(1/(i+a)^s, i=0..INF), s<>1, a>0}

function sfd_zetaint(n: integer): double;
  {-Return zeta(n) for integer arguments, n<>1}

function sfd_dbeta(s: double): double;
  {-Return the Dirichlet beta function sum((-1)^n/(2n+1)^s, n=0..INF)}

function sfd_dlambda(s: double): double;
  {-Return the Dirichlet lambda function sum(1/(2n+1)^s, n=0..INF), s<>1}

function sfd_lchi(s, x: double): double;
  {-Return Legendre's Chi-function chi(s,x); s>=0; real part if |x| > 1}

function sfd_bernpoly(n: integer; x: double): double;
  {-Return the Bernoulli polynomial B_n(x), 0 <= n <= MaxBernoulli}

function sfd_beint(s, x: double): double;
  {-Return the Bose-Einstein integral of real order s >= -1, real part if x > 1}

function sfd_fermi_dirac(n: integer; x: double): double;
  {-Return the integer order Fermi-Dirac integral F_n(x) = 1/n!*integral(t^n/(exp(t-x)+1), t=0..INF)}

function sfd_fdr(s,x: double): double;
  {-Return the Fermi-Dirac functions of real order s >= -1}

function sfd_harmonic(x: double): double;
  {-Return the harmonic function H(x) = psi(x+1) + EulerGamma}

function sfd_harmonic2(x,r: double): double;
  {-Return the generalized harmonic function H(x,r) = zeta(r)-zetah(r,x+1); x >= -1}

function sfd_llsi(x: double): double;
  {-Return the Lobachevski function Lambda(x) = integral(-ln(|2sin(t)|), t=0..x)}

function sfd_llci(x: double): double;
  {-Return the Lobachevski function L(x) = integral(-ln(|cos(t)|), t=0..x)}


{#Z+}
{---------------------------------------------------------------------------}
{Obsolete functions, use sfc_fdr with s = -1/2, 1/2, 3/2, 5/2 }
{#Z-}
function sfd_fdm05(x: double): double;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}

function sfd_fdp05(x: double): double;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}

function sfd_fdp15(x: double): double;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}

function sfd_fdp25(x: double): double;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}



implementation


uses
  DAMath,
  sdBasic, sdGamma;


{---------------------------------------------------------------------------}
function sfd_cl2(x: double): double;
  {-Return the Clausen function: integral(-ln(2*|sin(t/2)|),t=0..x) = Im(Li_2(exp(ix)))}
var
  t,z: double;
  n: integer;
const
  {h:= x-> 1+sum((-1)^(k-1)*bernoulli(2*k)/(2*k+1)!/(2*k)*x^(2*k),k=1..201);}
  {chebyshev(h(x),x=-Pi/2..Pi/2,0.1e-20); (Note: C1H[0] is doubled)}
  C1H: array[0..9] of THexDblW = (
         ($4549,$504A,$46D9,$4000),  {2.03459418036132851087309706802     }
         ($ACF2,$88A1,$C4AF,$3F91),  {0.173518588202740768068747978018e-1 }
         ($0FEA,$85A5,$EBD3,$3F0C),  {0.551628042609052141552962667866e-4 }
         ($617B,$4E5E,$B26F,$3E9A),  {0.397816462765976331046376760966e-6 }
         ($A82F,$D268,$B2CA,$3E2F),  {0.369018028917878121042765628603e-8 }
         ($7C16,$C253,$5530,$3DC5),  {0.388040921363679206236694833253e-10}
         ($4C25,$1633,$02E4,$3D5F),  {0.440696976897243799778208200770e-12}
         ($8C85,$92D7,$C3AA,$3CF7),  {0.527673937802511966876329709918e-14}
         ($60A1,$5B48,$EEA0,$3C92),  {0.656840358046069499876295250472e-16}
         ($ECDB,$1C8B,$1413,$3C2F)); {0.842382170379712962962962962962e-18}
        {($E0BA,$9810,$1E8C,$3BCA)}  {0.110619677185377920055584807691e-19}
const
  {g:=x->-I*(polylog(2,exp(I*(x+Pi)))-polylog(2,exp(-I*(x+Pi))))/2;}
  {chebyshev(g(x)/x,x=-Pi/2..Pi/2, 0.1e-20);(Note: C2H[0] is doubled)}
  C2H: array[0..14] of THexDblW = (
         ($8B5A,$BCF3,$737B,$BFF4),  {-1.27819417771453068261437382719     }
         ($42A4,$C622,$2669,$3FAC),  { 0.549805693018517156397035696498e-1 }
         ($A4B6,$4F74,$7FA6,$3F4F),  { 0.961261945950606429385907687403e-3 }
         ($C215,$A716,$CE4E,$3F00),  { 0.320546868225504765586825318140e-4 }
         ($F307,$97C1,$4DFE,$3EB6),  { 0.132946169542554501413438286774e-5 }
         ($D7C5,$3209,$AB0A,$3E70),  { 0.620936018243975194590942777798e-7 }
         ($B877,$8B32,$E210,$3E2A),  { 0.312960065639111267232623583835e-8 }
         ($7F89,$31C8,$DCFD,$3DE6),  { 0.166351953819266977593393018006e-9 }
         ($3082,$411F,$392F,$3DA4),  { 0.919652725071942544960123750608e-11}
         ($C079,$7054,$6FCE,$3D62),  { 0.524003773875845009365555555555e-12}
         ($AEE9,$0AA3,$3719,$3D21),  { 0.305803841873659454138326548177e-13}
         ($832D,$75DD,$63EC,$3CE0),  { 0.181969182494879509888450831696e-14}
         ($130D,$8CE5,$B781,$3C9F),  { 0.110039826319626150935008681510e-15}
         ($6F36,$D4BE,$1B4C,$3C5F),  { 0.674517757154247023242486416442e-17}
         ($C820,$8B96,$DD0E,$3C1E)); { 0.418278465157247206583308946518e-18}
      (* ($E123,$A193,$EE13,$3BDE),  { 0.261987180876106877468620455859e-19}
         ($EF56,$5827,$44A4,$3B9F),  { 0.165532116203486192256574495540e-20} *)
var
  C1: array[0..9] of double absolute C1H;
  C2: array[0..14] of double absolute C2H;
begin
  {See MISCFUN [22], function clausn for formulas and hints. As pointed}
  {out, only absolute accuracy can be guaranteed close to the zeros.   }
  {Observed relative errors are e.g. 1700 eps_d for abs(x-Pi) ~ 6.7E-4,}
  {but even for abs(x-Pi) ~ 0.03 they are still about 128 eps_d.}

  {Therefore two separate Chebyshev expansions are calculated, one for }
  {z in (-Pi/2, Pi/2) the other for z in (Pi/2, 3*Pi/2), where z is the}
  {reduced argument. For z=0 or z=Pi cl2 is zero. Calculations are done}
  {using Maple with Digits:=50; and t_rcalc to convert to Hex/double.}
  {Note that both approximations are done for even functions, and there}
  {are only the even Chebyshev polynomials, so the argument for CSEvalD}
  {is 2(x/(Pi/2))^2 - 1 with the calculated coefficients.}

  {Argument reduction x mod Pi, |z| <= Pi/2}
  n := rem_pio2(0.5*x,z);
  z := 2.0*z;
  t := abs(z);

  if t=0.0 then begin
    {if x is an exact multiple of Pi then cl2(x)=0}
    sfd_cl2 := 0.0
  end
  else begin
    if odd(n) then begin
      {Use approximation for (Pi/2, 3*Pi/2) centered at Pi}
      if t<=0.5e-8 then begin
        {cl2(x+Pi) = x*(-ln(2) + 1/24*x^2 + 1/960*x^4+ O(x^6))}
        sfd_cl2 := -ln2*z;
      end
      else begin
        {Use Chebyshev expansion for cl2(x+Pi)/x calculated from}
        {cl2(x) = Im(Li_2(exp(ix) = -i*(Li_2(exp(ix))-Li_2(exp(-ix)))/2}
        sfd_cl2 := CSEvalD(2.0*(sqr(2.0*z/Pi)-0.5),C2,15) * z;
      end;
    end
    else begin
      {Use approximation for (-Pi/2, Pi/2)}
      if t<=0.9e-8 then begin
        {cl2(x) = x*(1-ln(x) + 1/72*x^3 + 1/14400*x^5 + O(x^6) }
        sfd_cl2 := z*(1.0-ln(t));
      end
      else begin
        {Use Chebyshev expansion calculated from [22], (3) or HMF[1], 27.8.2}
        {cl2(x) = -x*ln|x| + x - sum((-1)^k*B_2k/(2k+1)!/(2k)*x^(2k+1),k=1..Inf)}
        sfd_cl2 := (CSEvalD(2.0*(sqr(2.0*z/Pi)-0.5),C1,10)-ln(t))*z;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_dilog(x: double): double;
  {-Return dilog(x) = Re(Li_2(x)), Li_2(x) = -integral(ln(1-t)/t, t=0..x)}
var
  t: double;
const
  xbig = 9007199254740992.0;
const
  nsp = 21;
  csp: array[0..nsp-1] of double = (
         0.1527365598892405872946684910028e+00,
         0.8169658058051014403501838185271e-01,
         0.5814157140778730872977350641182e-02,
         0.5371619814541527542247889005319e-03,
         0.5724704675185826233210603054782e-04,
         0.6674546121649336343607835438589e-05,
         0.8276467339715676981584391689011e-06,
         0.1073315673030678951270005873354e-06,
         0.1440077294303239402334590331513e-07,
         0.1984442029965906367898877139608e-08,
         0.2794005822163638720201994821615e-09,
         0.4003991310883311823072580445908e-10,
         0.5823462892044638471368135835757e-11,
         0.8576708692638689278097914771224e-12,
         0.1276862586280193045989483033433e-12,
         0.1918826209042517081162380416062e-13,
         0.2907319206977138177795799719673e-14,
         0.4437112685276780462557473641745e-15,
         0.6815727787414599527867359135607e-16,
         0.1053017386015574429547019416644e-16,
         0.1635389806752377100051821734570e-17);
      {  0.2551852874940463932310901642581e-18,
         0.3999020621999360112770470379519e-19,
         0.6291501645216811876514149171199e-20);}
const
  Pi26H: THexDblW = ($07D3,$6253,$51A6,$3FFA);
var
  Pi26 : double absolute Pi26H;
begin
  {Ref: W. Fullerton [14] and [20], file dspenc.f}

  {Note that there is some confusion about the naming: some authors}
  {and/or computer algebra systems use dilog(x)= Li_2(1-x) and then}
  {call Li_2(x) Spence function/integral or similar.}

  {The imaginary part Im(Li_2(x) is 0 for x<=1 and -Pi*ln(x) for x>1}

  if x>2.0 then begin
    t := 2.0*Pi26 - 0.5*sqr(ln(x));
    if x>=xbig then sfd_dilog := t
    else sfd_dilog := t - (1.0 + CSEvalD(4.0/x-1.0, csp, nsp))/x;
  end
  else if x>1.0 then begin
    {1 < x <= 2}
    sfd_dilog := Pi26 - 0.50*ln(x)*ln(sqr((x-1.0))/x)
                 + (x-1.0)*(1.0 + CSEvalD(4.0*(x-1.0)/x-1.0, csp, nsp))/x;
  end
  else if x>0.5 then begin
    {0.5 < x <= 1}
    if x=1.0 then sfd_dilog := Pi26
    else sfd_dilog := Pi26 - ln(x)*ln1p(-x) - (1.0-x)*(1.0+CSEvalD(4.0*(1.0-x)-1.0, csp, nsp));
  end
  else if x>=0.0 then begin
    {0 <= x <= 0.5}
    sfd_dilog := x*(1.0 + CSEvalD(4.0*x-1.0, csp, nsp));
  end
  else if x>-1.0 then begin
    {-1 < x < 0}
    sfd_dilog := -0.5*sqr(ln(1.0-x)) - x*(1.0+CSEvalD(4.0*x/(x-1.0)-1.0, csp, nsp))/(x-1.0);
  end
  else begin
    {x <= -1.0}
    t := ln1p(-x);
    t := -Pi26 - 0.5*t*(2.0*ln(-x)-t);
    if x <= -xbig then sfd_dilog := t
    else begin
      sfd_dilog := t + (1.0 + CSEvalD(4.0/(1.0-x)-1.0, csp, nsp))/(1.0-x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_trilog(x: double): double;
  {-Return the trilogarithm function trilog(x) = Re(Li_3(x))}
var
  a,b: double;
begin
  if IsNanD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_trilog := Nan_d;
    exit;
  end;
  if abs(x)>1.0 then begin
    b := 1.0/x;
    if abs(b) > sqrt_epsh then a := sfd_polylog(3,b)
    else begin
      {large |x|, use polylog(3,b) = b*(1 + b/8 + b^2/27 + O(b^3))}
      a := b*(1.0 + 0.125*b);
    end;
    {transform to intervall (-1,1) via inversion formulas}
    if x>0.0 then begin
      {This is the 'new' part for x>1. The real part of the RHS of}
      {Crandall's [31] formula (1.3) is Pi^2/3 * ln(x) - ln^3(x)/6}
      b := ln(x);
      b := (PiSqr - 0.5*sqr(b))*b/THREE;
      sfd_trilog := a + b;
    end
    else begin
      {Could be handled by sfd_polylog, but this branch}
      {avoids one recursive call with sfd_polylog(1/x):}
      {Re(RHS) becomes -Pi^2/6 * ln(-x) - ln^3(-x)/6 }
      b := ln(-x);
      b := (PiSqr + sqr(b))*b/SIXX;
      sfd_trilog := a - b;
    end;
  end
  else begin
    {-1 <= x <= 1}
    sfd_trilog := sfd_polylog(3,x);
  end;
end;


{---------------------------------------------------------------------------}
function polylog_series(n: integer; x: double): double;
  {-polylog via series, |n|>2, |x| < 1}
var
  e,k,p,s,t,w,y,z: double;
begin
  {Standard series sum(x^k/k^n, k=1..) calculated with sum2 algorithm.}
  {Actually sum2 is needed for x < 0, especially for large negative n.}
  p := x;
  s := p;
  e := 0.0;
  k := 1.0;
  repeat
    p := p*x;
    k := k+1.0;
    if n>0 then t := p/power(k,n)
    else t := p*power(k, -n);
    {(w,y) = TwoSum(t,s)}
    w := t+s;
    z := w-t;
    y := (t-(w-z)) + (s-z);
    {sum of errors}
    e := e+y;
    s := w;
  until abs(t) <= eps_d*abs(s);
  polylog_series := s+e;
end;


{---------------------------------------------------------------------------}
function cr15_complex(n: integer; x: double): double;
  {-Compute polylog(n,x) using Crandall[31] (1.5); n<-1, |log(x)| < 2*Pi}
type
  TComplex = record r,i: double; end;
var
  lx,s,p,t: TComplex;
  bk: double;
  k: integer;
  done: boolean;
begin
  {This function uses Crandall[31] (1.5) with (inline) complex arithmetic}
  done := false;
  k    := 0;
  {lx = ln(x)}
  lx.r := ln(abs(x));
  lx.i := Pi;
  p.r  := 1.0;
  p.i  := 0.0;
  s.r  := sfd_bernoulli(1-n)/(1-n);
  s.i  := 0.0;
  repeat
    inc(k);
    t.r := p.r*lx.r - p.i*lx.i;
    t.i := p.r*lx.i + p.i*lx.r;
    p.r := t.r/k;
    p.i := t.i/k;
    bk := sfd_bernoulli(k-n+1)/(k-n+1);
    if bk<>0.0 then begin
      t.r := p.r*bk;
      t.i := p.i*bk;
      s.r := s.r + t.r;
      s.i := s.i + t.i;
      done := (abs(t.r)+abs(t.i)) <= eps_d*(abs(s.r)+abs(s.i));
    end;
  until done;
  {rt.r = Re(power(-lx, n-1))}
  p.r := hypot(lx.r,lx.i);
  p.i := arctan2(-lx.i,-lx.r);
  t.r := power(p.r,n-1)*cos((n-1)*p.i);
  cr15_complex := sfd_fac(-n)*t.r - s.r;
end;


{---------------------------------------------------------------------------}
function polylog_neg(n: integer; x: double): double;
  {-Polylog for n<0 and x<1}
var
  s,p,lx,t: double;
  k: integer;
  done: boolean;
const
  a4: array[0..3] of single = (1, 11, 11, 1);
  a5: array[0..4] of single = (1, 26, 66, 26, 1);
  a6: array[0..5] of single = (1, 57, 302, 302, 57, 1);
  a7: array[0..6] of single = (1, 120, 1191, 2416, 1191, 120, 1);
  a8: array[0..7] of single = (1, 247, 4293, 15619, 15619, 4293, 247, 1);
  a9: array[0..8] of single = (1, 502, 14608, 88234, 156190, 88234, 14608, 502, 1);
begin
  if x = 1.0 then polylog_neg := PosInf_d
  else if n >= -9 then begin
    p := power(1.0-x,1-n);
    case n of
       -9: s := PolEvalS(x, a9, 9);
       -8: s := PolEvalS(x, a8, 8);
       -7: s := PolEvalS(x, a7, 7);
       -6: s := PolEvalS(x, a6, 6);
       -5: s := PolEvalS(x, a5, 5);
       -4: s := PolEvalS(x, a4, 4);
       -3: s := 1.0 + x*(4.0 + x);
       -2: s := 1.0 + x;
      else s := 1.0
    end;
    polylog_neg := x/p*s;
  end
  else if x=-1.0 then polylog_neg := -sfd_etaint(n)
  else if abs(x) <= 0.25 then polylog_neg := polylog_series(n,x)
  else if (x > 4.0) or (x < -11.5)  then begin
    s := polylog_neg(n,1.0/x);
    if odd(n) then polylog_neg := s
    else polylog_neg := -s;
  end
  else if x>0 then begin
    {Use Crandall[31] (1.5) with real arithmetic}
    lx := ln(x);
    k  := 0;
    p  := 1.0;
    s  := sfd_bernoulli(1-n)/(1-n);
    done := false;
    repeat
      inc(k);
      p := p*lx/k;
      t := sfd_bernoulli(k-n+1);
      if t<>0.0 then begin
        t := t/(k-n+1)*p;
        s := s + t;
        done := abs(t) <= abs(s)*eps_d;
      end;
    until done;
    polylog_neg := sfd_fac(-n)*power(-lx, n-1) - s;
  end
  else begin
    {Same formula as above but with complex arithmetic}
    polylog_neg := cr15_complex(n,x);
  end;
end;


{---------------------------------------------------------------------------}
function polylog_inv(n: integer; x: double): double;
  {-Return Li_(n,x) for x < -1 with inversion formula}
var
  s,t,z: double;
  j,k: integer;
begin
  {Inversion formula for x < -1, see Cephes[7] and use Li_k(-1) = -eta(k)}
  {Li_n(-z) + (-1)^n*Li_n(-1/z) = -ln(z)^n/n! - 2*sum(ln(z)^(n-2k)/(n-2k)!*eta(2k), k=1..n div 2)}
  {See also: L. Lewin, Polylogarithms and associated functions, 1981, (7.20)}
  {or http://functions.wolfram.com/10.08.17.0060.01}
  z := ln(-x);
  s := 0.0;
  for k:=1 to n div 2 do begin
    j := 2*k;
    t := -sfd_etaint(j);
    j := n-j;
    if j>0 then t := sfd_taylor(z,j)*t;
    s := s + t;
  end;
  t := sfd_polylog(n, 1.0/x);
  if odd(n) then t := -t;
  polylog_inv := 2.0*s - t - sfd_taylor(z,n);
end;


{---------------------------------------------------------------------------}
function polylog_xg05(n: integer; x: double): double;
  {-Return Re(polylog(n,x)), n>2, x >= 0.5}
var
  s,t,z,f,h: double;
  k,i: integer;
const
  xmax = 256.0;
  hinv = -123456789.0; {flag h not yet valid}
begin
  if x > xmax then begin
    {Crandall [31], Algorithm 3.1 for positive x > 1, with r2=xmax}
    {Recursive Li_n(x) = (Li_n(sqrt(x)) + Li_n(-sqrt(x)))/2^(n-1) }
    z := sqrt(x);
    h := sfd_polylog(n, -z);
    f := polylog_xg05(n,z);
    polylog_xg05 := ldexpd(f+h, n-1);
    exit;
  end;
  {Crandall [31], formula 1.4; valid for |x| < exp(2*pi) }
  z := ln(x);
  f := 1;
  s := sfd_zetaint(n);
  k := 0;
  h := hinv;
  repeat
    inc(k);
    f := f*z/k;
    i := n-k;
    if i=1 then begin
      {store final h value}
      h := f;
    end
    else if (i >= 0) or odd(i) then begin
      {sum only non-zero Bernoulli terms}
      t := sfd_zetaint(i)*f;
      s := s + t;
      if abs(t)<=eps_d*abs(s) then break;
    end;
  until false;
  f := sfd_harmonic(n-1) - ln(abs(z));
  if h=hinv then begin
    {h = z^(n-1)/(n-1) not yet calculated}
    h := sfd_taylor(z, n-1);
  end;
  polylog_xg05 := s + h*f;
end;


{---------------------------------------------------------------------------}
function sfd_polylog(n: integer; x: double): double;
  {-Return the polylogarithm Li_n(x) of integer order; real part for n>0,x>1}
var
  s,t: double;
begin

  sfd_polylog := Nan_d;
  {Ref: Cephes [7] function polylog in file polylog.c and}
  {Crandall [31], Note on fast polylogarithm computation.}

  if IsNanOrInfD(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;

  if x=0.0 then begin
    sfd_polylog := 0.0;
    exit;
  end
  else if abs(x)=1.0 then begin
    if x>0.0 then  sfd_polylog := sfd_zetaint(n)
    else sfd_polylog := -sfd_etaint(n);
    exit;
  end
  else if n<0 then begin
    sfd_polylog := polylog_neg(n,x);
    exit;
  end
  else if n<=2 then begin
    case n of
        0: sfd_polylog := x/(1.0-x);
        1: begin
             if x<1.0 then sfd_polylog := -ln1p(-x)
             else sfd_polylog := -ln(x-1.0);
           end;
      else sfd_polylog := sfd_dilog(x);
    end;
    exit;
  end
  else if (n>40) and (x > -1.0) then begin
    if n > 54 then sfd_polylog := x
    else sfd_polylog := x*(1.0 + x*ldexpd(1,-n));
    exit;
  end;

  {Here n>2}
  if x < -1.0 then begin
    {Inversion formula for x < -1}
    sfd_polylog := polylog_inv(n,x);
  end
  else if (x < 0.5) and (x > -0.875) then begin
    {Use standard series sum(x^k/k^n, k=1..)}
    sfd_polylog := polylog_series(n,x);
  end
  else if x < 0.0 then begin
    {For -1 < x <= -0.875 use recursive calls: Crandall[31], formula (1.2)}
    s := sfd_polylog(n,x*x);
    t := sfd_polylog(n,-x);
    sfd_polylog := ldexpd(s,1-n) - t;
  end
  else begin
    {Here x>=0.5, use Crandall[31], formula (1.4)}
    sfd_polylog := polylog_xg05(n,x);
  end;
end;


{---------------------------------------------------------------------------}
function liae_neg(s,x: double): double;
  {-Asymptotic formula for Li(s,x), x << -1}
var
  w,sum,r,t,f: double;
  k: integer;
begin
  {Ref: Wood [73], (11.1)}
  k := 0;
  w := ln(abs(x));
  f := power(w, s)*sfd_rgamma(s+1.0);
  t := MaxDouble;
  sum := 0.5*f;
  w  := w*w;
  repeat
    f := f/w*(s-k)*(s-k-1);
    k := k + 2;
    r := t;
    t := sfd_etaint(k)*f;
    sum := sum + t;
    t := abs(t);
  until (t >= r) or (t <= eps_d*abs(sum));
  liae_neg := -2.0*sum;
end;


{---------------------------------------------------------------------------}
type
  tli_rec  = record
               s,x,g: double;  {g=1/Gamma(s)}
             end;
  pli_rec = ^tli_rec;



{---------------------------------------------------------------------------}
function lif(t: double; param: pointer): double; {$ifdef BIT16} far;{$endif}
  {-Polylog integrand t^(s-1)/(exp(t)/x-1)/Gamma(s) }
var
  z,p: double;
begin
  with pli_rec(param)^ do begin
    z := exp(-t);
    if z=0.0 then lif := 0.0
    else begin
      p := power(t, s - 1.0)*g;
      lif := x*p*z/(1.0-x*z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function liae_pos(s,lnx: double): double;
  {-Asymptotic formula for Li(s,ln(x)), x >> 1}
var
  w,sum,r,t,f: double;
  k: integer;
const
  KMAX = 500;
begin
  {Ref: Wood [73], (11.2)}
  k := 0;
  w := lnx;
  f := power(w, s)*sfd_rgamma(s+1.0);
  t := Maxdouble;
  sum := -0.5*f;
  w  := w*w;
  repeat
    f := f/w*(s-k)*(s-k-1);
    k := k + 2;
    r := t;
    t := sfd_zetaint(k)*f;
    sum := sum + t;
    t := abs(t);
  until (k>KMAX) or (t >= r) or (t <= eps_d*abs(sum));
  {Test no convergence}
  if (k>KMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  liae_pos := 2.0*sum;
end;


{---------------------------------------------------------------------------}
function li_xgt0(s,x: double): double;
  {-Return Re(polylog(s,x)), s > -1, 1 < x, frac(s)<>0}
var
  sum,t,z,f,h: double;
  k: integer;
const
  xmax = 256;   {< exp(2*Pi), 256 ~ exp(5.5)}
begin
  if (x > xmax) and (s > 0.0) then begin
    if x >= 1e20 then begin
      li_xgt0 := liae_pos(s,ln(x));
    end
    else begin
      {Recursive computation with Wood [73] (15.1), x > 1 }
      {Li_s(x) = (Li_s(sqrt(x)) + Li_s(-sqrt(x)))/2^(s-1) }
      z := sqrt(x);
      h := sfd_polylogr(s, -z);
      f := li_xgt0(s,z);
      li_xgt0 := 0.5*(f+h)*exp2(s);
    end;
    exit;
  end;
  {NIST [30], 25.12.12; valid for |x| < exp(2*pi) }
  z := ln(x);
  f := z;
  sum := sfd_zeta(s-1)*f;
  for k:=2 to maxint do begin
    f := f*z/k;
    t := sfd_zeta(s-k)*f;
    sum := sum + t;
    if abs(t)<=eps_d*abs(sum) then break;
  end;
  sum := sum + sfd_zeta(s);
  f := sfd_gamma(1.0-s);
  if f<>0.0 then begin
    {compute h = Re((-z)^s-1)}
    {z=arg(-ln x) = 0 or Pi}
    if z < 0.0 then t := 1.0 else t := cospi(s-1.0);
    h := power(abs(z),s-1.0)*t;
    sum := sum + h*f;
  end;
  li_xgt0 := sum;
end;


{---------------------------------------------------------------------------}
function sfd_polylogr(s,x: double): double;
  {-Return the polylogarithm Li_s(x) of real order s >= -1, x <= 256;}
  { s > 0 if x > 256; real part if x > 1}
var
  z, aerr: double;
  ierr: integer;
  para: tli_rec;
  neval: longint;
const
  eps = 4.0e-16;
begin
  sfd_polylogr := NaN_d;
  if IsNanOrInfD(s) or IsNanOrInfD(x) or (s < -1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;
  z := abs(x);
  if z=0.0 then sfd_polylogr := 0.0
  else if z=1.0 then begin
    if x=1.0 then sfd_polylogr := sfd_zeta(s)
    else sfd_polylogr := -sfd_eta(s)
  end
  else if frac(s)=0.0 then begin
    {use integer version}
    sfd_polylogr := sfd_polylog(round(s),x)
  end
  else if z<1.0 then begin
    if s>=34.5 then begin
      {from series and z<1}
      if s>=54.0 then sfd_polylogr := x
      else sfd_polylogr := x*(1.0 + x*exp2(-s))
    end
    else begin
      {Li_s(x) = x*Phi(x,s,1), see e.g. NIST[30], 25.14.3}
      {or http://functions.wolfram.com/10.08.26.0010.01}
      sfd_polylogr := x*sfd_lerch(x,s,1.0);
    end;
  end
  else if x > 1.0 then begin
    {will return NaN if x > 256}
    sfd_polylogr := li_xgt0(s,x);
  end
  else if s>0.0 then begin
    if x <= -1e20 then begin
      {Use asymptotic formula Wood [73], (11.1)}
      sfd_polylogr := liae_neg(s,x);
    end
    else if s > log2(z) + 55.0 then begin
      sfd_polylogr := x;
    end
    else begin
      {Use integral formula, see Wood [73], (1.1) or }
      {http://functions.wolfram.com/10.08.07.0001.01}
      para.s := s;
      para.x := x;
      para.g := sfd_rgamma(s);
      sfd_intdei_p({$ifdef FPC_ProcVar}@{$endif}lif, @para, 0.0, eps, z, aerr, neval, ierr);
      if ierr=0 then sfd_polylogr := z;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ti(s,x: double): double;
  {-Return the inverse tangent integral of order s >= 0}
var
  t,z: double;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(s) or (s < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_ti := NaN_d;
    exit;
  end;
  t := abs(x);
  if t=0.0 then sfd_ti := 0.0
  else if s=0.0 then begin
    {ti(0,x) = x/(1+x*2)}
    if t > 1e10 then sfd_ti := 1.0/x
    else sfd_ti := x/(1.0 + x*x);
  end
  else if s=1.0 then begin
    {ti(1,x) = arctan(x)}
    sfd_ti := arctan(x);
  end
  else if s=2.0 then begin
    sfd_ti := sfd_ti2(x);
  end
  else begin
    z := 1.85*ln(t)+33.5;
    if s > z then begin
      {shortcut for x^2/3^s < eps_x/2}
      sfd_ti := x;
    end
    else begin
      t := sfd_lerch(-x*x,s,0.5);
      if IsNanOrInfD(t) then sfd_ti := t
      else sfd_ti := exp2(-s)*x*t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_ti2(x: double): double;
  {-Return the inverse tangent integral, ti2(x) = integral(arctan(t)/t, t=0..x)}
var
  y: double;
const
  {chebyshev(int(arctan(t)/t, t=0..x)/x, x=-1..1, 0.5e-20); aics[0] doubled}
  nai = 21;
  aics: array[0..nai-1] of double = (
          +0.191040361296235937512100085363e+1,
          -0.417635143765674693973184505426e-1,
          +0.275392550786367433987311337253e-2,
          -0.250518095262488814159890349742e-3,
          +0.266698128512117108101673214400e-4,
          -0.311890514107001322989208298185e-5,
          +0.388338531322492994856981431022e-6,
          -0.505727458496376372130362845503e-7,
          +0.681225282949264960525940397947e-8,
          -0.942125616543636974084406666666e-9,
          +0.133078788164081272284621690761e-9,
          -0.191267807507294407420492895005e-10,
          +0.278912620074783714512526059665e-11,
          -0.411748196101712508143369500279e-12,
          +0.614298719454016651500191474755e-13,
          -0.924928654021041585596995204818e-14,
          +0.140386740366392221748521257509e-14,
          -0.214598969089468659553281910650e-15,
          +0.330121857094999796296296296296e-16,
          -0.510717131686183481708845972024e-17,
          +0.794150680705224718491128609488e-18);
{         -0.124060349865156550047951820204e-18,
          +0.194621210567853358154874940741e-19,
          -0.306490766819307325157212222222e-20}
begin
  {ti2(x) = integral(arctan(t)/t,  t=0..x)}
  {ti2(x) = i/2*[dilog(-i*x) - dilog(i*x)]}
  {ti2(x) = ti2(1/x) + Pi/2*ln(x),  x<>0}

  {See references [22] (function atnint) and [21] (specfunc/atanint.c)}
  {Note that the precision of both refs are not suitable for double!}

  y := abs(x);
  if y<1.5e-4 then begin
    {ti2(x) := x*(1 - 1/9x^2 + 1/25x^4 +O(x^6)}
    if y<1e-8 then sfd_ti2 := x
    else sfd_ti2 := x*(1.0-sqr(x)/9.0);
  end
  else begin
    {ti2(x)/x is even and it's Chebyshev approximation on [-1,1] contains }
    {only the even Chebyshev polynomials. Since T(2n,x) = T(n,2x^2-1), the}
    {coefficients aics can be used with CSEvalD and the argument 2x^2-1.  }
    if y<=1.0 then sfd_ti2 := x*CSEvalD(2.0*sqr(x)-1.0, aics, nai)
    else begin
      {ti2(x) = Pi/2*ln(x) + 1/x + O(1/x^2), x->Inf}
      if y>=1.9e8 then y := Pi_2*ln(y) + 1.0/y
      else begin
        {Use ti2(x) = ti2(1/x) + sign(x)*Pi/2*ln(|x|) for |x|>1}
        y := Pi_2*ln(y) + CSEvalD(2.0/sqr(x)-1.0, aics, nai)/y;
      end;
      if x>0.0 then sfd_ti2 := y else sfd_ti2 := -y;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_pz(x: double): double;
  {-Return the prime zeta function P(x) = sum(1/p^x, p prime), x > 1/5; }
  { for x<1 the real part of P(x) is returned.}
var
  s,t,y: double;
  k,m: integer;
begin
  {Im(P(x)) = Pi   for 1/2 < x < 1
            = Pi/2 for 1/3 < x < 1/2
            = Pi/6 for 1/5 < x < 1/3
  }
  if IsNanD(x) or (x < 0.2) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_pz := NaN_d;
    exit;
  end;
  if x>=27.0 then begin
    {use only primes up to 7}
    if x>=16450.0 then sfd_pz := 0.0
    else if x>120.0 then sfd_pz := exp2(-x)
    else begin
      t := exp3(-x);
      if x<50.0 then begin
        s := exp5(-x);
        t := t + s;
        if x<36.0 then begin
          s := exp7(-x);
          t := t + s;
        end;
      end;
      s := exp2(-x);
      sfd_pz := s+t;
    end;
  end
  else if x=1.0 then sfd_pz := PosInf_d
  else begin
    {P(x) = sum(mu(k)/k*ln(zeta(kx)), k>0) by Moebius inversion, see e.g.}
    {H. Cohen, High Precision Computation of Hardy-Littlewood Constants, }
    {Section 2.1, from  http://www.math.u-bordeaux.fr/~cohen/hardylw.dvi }
    {Cohen's speed-up is used with A=7.}
    y := -x;
    s := exp2(y) + exp3(y) + exp5(y) + exp7(y);
    {double: maximum observed k is 77 for x near 1/5}
    for k:=1 to n_moeb do begin
      m := moebius[k];
      if m<>0 then begin
        y := k*x;
        if y=1.0 then begin
          sfd_pz := NegInf_d;
          exit;
        end;
        {Compute the next term ln(zeta(y)*prod(1-p^y))*mu(k)/k}
        if y<1.0 then t := ln(abs(sfd_zeta(y)))
        else t := ln1p(sfd_zetam1(y));
        t := t + ln1p(-exp2(-y));
        t := t + ln1p(-exp3(-y));
        t := t + ln1p(-exp5(-y));
        t := t + ln1p(-exp7(-y));
        s := s + m*t/k;
        if abs(t) < eps_d*abs(s) then begin
          sfd_pz := s;
          exit;
        end;
      end;
    end;
    {No convergence}
    sfd_pz := s;
    if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  end;
end;


const
  ETAEPS = 1e-9;

{---------------------------------------------------------------------------}
function etam1pos(s: double): double;
  {-Return the eta(s)-1 = -sum((-1)^k/k^s, k=2..) for s >= -ETAEPS; internal use}
const
  dh: array[2..21] of THexDblW = (
        ($FFF7,$FFFF,$FFFF,$BFEF),  {-0.99999999999999902275696874687457448602590099318551   }
        ($E475,$FFFF,$FFFF,$3FEF),  { 0.99999999999921722833196624653416330674669554159345   }
        ($9A09,$FFF1,$FFFF,$BFEF),  {-0.99999999989523856980663370125947646261237047984975   }
        ($0D42,$FCFF,$FFFF,$3FEF),  { 0.99999999440516539966907531075601109232000721978239   }
        ($8CFB,$AAAD,$FFFF,$BFEF),  {-0.99999984107526471939869454740922825058329045647251   }
        ($21D6,$2DC7,$FFFA,$3FEF),  { 0.99999722424495977611752951962413441827665769598385   }
        ($1A7E,$D101,$FFBB,$BFEF),  {-0.99996748753694905701338147661170450570128541770365   }
        ($DFB9,$EAD0,$FDC8,$3FEF),  { 0.99972959387286330418019713251226520509830719146205   }
        ($EE53,$5B69,$F262,$BFEF),  {-0.99833791593796165010606871953054529657088456794866   }
        ($65AE,$9541,$C04F,$3FEF),  { 0.99222544814545242436872039584612844343083226075262   }
        ($3F37,$7D35,$182B,$BFEF),  {-0.97170042471586881373488760368477079972897240816801   }
        ($DEA2,$A5CF,$6370,$3FED),  { 0.91838867554811917572493229936955614076310266119499   }
        ($15CD,$77E7,$F07C,$BFE9),  {-0.81060622614375577713958787977575259111471382492346   }
        ($89DC,$8C86,$817E,$3FE4),  { 0.64080741354365097690667605567720669136094433516647   }
        ($684D,$4B76,$BAC9,$BFDB),  {-0.43327553147685622106645049289009503610633718101903   }
        ($FE47,$11BC,$8B43,$3FCE),  { 0.23862493864179355351975617193114893186753322954281   }
        ($0B79,$3C28,$ECE0,$BFB9),  {-0.10127068966543884859769364706092486637644173152532   }
        ($9511,$03D1,$7ED6,$3F9F),  { 0.30757278426240711311501120924339249867646096179440e-1}
        ($CAE1,$C8CD,$34C3,$BF78),  {-0.59096954181423200773189926666852707169276342004165e-2}
        ($4DBB,$EF21,$9ABC,$3F41)); { 0.53724503801293818884718115151684279244796674549241e-3}

const
  ln2m1h: THexDblW = ($8C22,$020B,$A37A,$BFD3);  {-3.0685281944005469057E-1}
  lnpi2h: THexDblW = ($1316,$25AA,$E6BB,$3FCC);  {0.2257913526447274323630976} {ln(pi/2)/2}
var
  d: array[2..21] of double absolute dh;
  p: array[2..21] of double;
  x,sum: double;
  k: integer;
begin

  if s=1.0 then begin
    etam1pos := double(ln2m1h);
    exit;
  end
  else if abs(s) <= ETAEPS then begin
    etam1pos := double(lnpi2h)*s - 0.5;
    exit;
  end;

  x := -s;

  {Calculate p[k] := 1/k^s but only if necessary. Prime powers are}
  {evaluated with power/exp?, otherwise products of p[k] are used.}
  p[2] := exp2(x);
  if s >= 102.0 then begin
    etam1pos := -p[2];
    exit;
  end;

  p[3] := exp3(x);
  p[4] := p[2]*p[2];
  if s >= 45.3 then begin
    etam1pos := (p[3] - p[4]) - p[2];
    exit;
  end;

  p[10] := exp10(x);
  p[ 5] := p[10]/p[2];
  p[ 6] := p[2]*p[3];
  p[ 7] := exp7(x);
  p[ 8] := p[2]*p[4];
  p[ 9] := p[3]*p[3];
  if s >= 24.4 then begin
    sum := 0.0;
    for k:=10 downto 3 do begin
      if odd(k) then sum := sum + p[k]
      else sum := sum - p[k];
    end;
    etam1pos := sum - p[2];
    exit;
  end;

  p[11] := power(11,x);
  p[12] := p[2]*p[6];
  p[13] := power(13,x);
  p[14] := p[2]*p[7];
  p[15] := p[3]*p[5];
  p[16] := p[2]*p[8];
  p[17] := power(17,x);
  p[18] := p[2]*p[9];
  p[19] := power(19,x);
  p[20] := p[2]*p[10];
  p[21] := p[3]*p[7];
  if s>=17.0 then begin
    sum := p[2]*p[11];
    for k:=21 downto 3 do begin
      if odd(k) then sum := sum + p[k]
      else sum := sum - p[k];
    end;
    etam1pos := sum - p[2];
  end
  else begin
    {Convergence acceleration: see P. Borwein[36]}

    {The d[] are from P. Borwein's Algorithm 2, but scaled and shifted.}
    {Calculated with Maple VR4 and T_RCalc/dh using n:=20; Digits:=50; }

    { d(n,k) = n*sum((n+j-1)!*4^j/((n-j)!*(2*j)!),j=0..k); }
    { d[2]   := 1/d(n,n)-1.0;                              }
    { d[i]   := (-1)^(i-2)*(d(n,i-2)/d(n,n)-1), i=3..n+1;  }
    sum := d[21]*p[21];
    for k:=20 downto 3 do sum := sum + d[k]*p[k];
    etam1pos := sum + d[2]*p[2];
  end;
end;


{---------------------------------------------------------------------------}
function zetap(s,sc: double): double;
  {-Return the Riemann zeta function at s>0, s<>1, sc=1-s}
var
  y: double;
{Based on boost_1_42_0\boost\math\special_functions\zeta.hpp [19]}
{Copyright John Maddock 2007, see 3rdparty.ama for Boost license}
const
  P1: array[0..5] of double = (
        0.24339294433593750202,
       -0.49092470516353571651,
        0.0557616214776046784287,
       -0.00320912498879085894856,
        0.000451534528645796438704,
       -0.933241270357061460782e-5);
  Q1: array[0..5] of double = (
        1.0,
       -0.279960334310344432495,
        0.0419676223309986037706,
       -0.00413421406552171059003,
        0.00024978985622317935355,
       -0.101855788418564031874e-4);
  P2: array[0..5] of double = (
        0.577215664901532860516,
        0.243210646940107164097,
        0.0417364673988216497593,
        0.00390252087072843288378,
        0.000249606367151877175456,
        0.110108440976732897969e-4);
  Q2: array[0..5] of double = (
        1.0,
        0.295201277126631761737,
        0.043460910607305495864,
        0.00434930582085826330659,
        0.000255784226140488490982,
        0.10991819782396112081e-4);
  P4: array[0..5] of double = (
       -0.0537258300023595030676,
        0.0445163473292365591906,
        0.0128677673534519952905,
        0.00097541770457391752726,
        0.769875101573654070925e-4,
        0.328032510000383084155e-5);
  Q4: array[0..6] of double = (
        1.0,
        0.33383194553034051422,
        0.0487798431291407621462,
        0.00479039708573558490716,
        0.000270776703956336357707,
        0.106951867532057341359e-4,
        0.236276623974978646399e-7);
  P7: array[0..5] of double = (
       -2.49710190602259410021,
       -2.60013301809475665334,
       -0.939260435377109939261,
       -0.138448617995741530935,
       -0.00701721240549802377623,
       -0.229257310594893932383e-4);
  Q7: array[0..8] of double = (
        1.0,
        0.706039025937745133628,
        0.15739599649558626358,
        0.0106117950976845084417,
       -0.36910273311764618902e-4,
        0.493409563927590008943e-5,
       -0.234055487025287216506e-6,
        0.718833729365459760664e-8,
       -0.1129200113474947419e-9);
 P15: array[0..6] of double = (
       -4.78558028495135619286,
       -1.89197364881972536382,
       -0.211407134874412820099,
       -0.000189204758260076688518,
        0.00115140923889178742086,
        0.639949204213164496988e-4,
        0.139348932445324888343e-5);
 Q15: array[0..8] of double = (
        1.0,
        0.244345337378188557777,
        0.00873370754492288653669,
       -0.00117592765334434471562,
       -0.743743682899933180415e-4,
       -0.21750464515767984778e-5,
        0.471001264003076486547e-8,
       -0.833378440625385520576e-10,
        0.699841545204845636531e-12);
 P36: array[0..7] of double = (
      -10.3948950573308896825,
       -2.85827219671106697179,
       -0.347728266539245787271,
       -0.0251156064655346341766,
       -0.00119459173416968685689,
       -0.382529323507967522614e-4,
       -0.785523633796723466968e-6,
       -0.821465709095465524192e-8);
 Q36: array[0..7] of double = (
        1.0,
        0.208196333572671890965,
        0.0195687657317205033485,
        0.00111079638102485921877,
        0.408507746266039256231e-4,
        0.955561123065693483991e-6,
        0.118507153474022900583e-7,
        0.222609483627352615142e-14);
const
   x1: double = 81487.0/65536.0;
   x4: double = 366299.0/524288.0;
begin
  {$ifdef debug}
    if abs((s+sc)-1.0)>1.1e-15*maxd(1.0,abs(s)) then begin
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      zetap := Nan_d;
      exit;
    end;
  {$endif}
  if s<1.0 then begin
    y := PolEval(sc,P1,6)/PolEval(sc,Q1,6) - x1;
    zetap := (y+sc)/sc;
  end
  else if s<=2.0 then begin
    s := -sc;
    y := PolEval(s,P2,6)/PolEval(s,Q2,6);
    zetap := y - 1.0/sc;
  end
  else if s<=4.0 then begin
    s := s-2.0;
    y := PolEval(s,P4,6)/PolEval(s,Q4,7) + x4;
    zetap := y - 1.0/sc;
  end
  else if s<=7.0 then begin
    s := s-4.0;
    y := PolEval(s,P7,6)/PolEval(s,Q7,9);
    zetap := 1.0 + exp(y);
  end
  else if s<15.0 then begin
    s := s-7.0;
    y := PolEval(s,P15,7)/PolEval(s,Q15,9);
    zetap := 1.0 + exp(y);
  end
  else if s<36.0 then begin
    s := s-15.0;
    y := PolEval(s,P36,8)/PolEval(s,Q36,8);
    zetap := 1.0 + exp(y);
  end
  else if s<56.0 then zetap := 1.0 + exp2(-s)
  else zetap := 1.0;
end;


{---------------------------------------------------------------------------}
function sfd_etaint(n: integer): double;
  {-Return the Dirichlet function eta(n) for integer arguments}
const
  etahex: array[0..53] of THexDblW = (
            ($0000,$0000,$0000,$3FE0),  {0.5                             }
            ($39EF,$FEFA,$2E42,$3FE6),  {0.693147180559945309417232121458}
            ($07D3,$6253,$51A6,$3FEA),  {0.822467033424113218236207583323}
            ($0932,$0768,$D970,$3FEC),  {0.901542677369695714049803621133}
            ($BA7E,$CADD,$4E17,$3FEE),  {0.947032829497245917576503234474}
            ($AA02,$EBBB,$1B9A,$3FEF),  {0.972119770446909305935655143556}
            ($1B65,$7135,$89A2,$3FEF),  {0.985551091297435104098439244483}
            ($F7F7,$1D58,$C354,$3FEF),  {0.992593819922830282670425713134}
            ($E59F,$0844,$E124,$3FEF),  {0.996233001852647899227289260081}
            ($7D3D,$73C8,$F063,$3FEF),  {0.998094297541605330767783031850}
            ($7D17,$B391,$F821,$3FEF),  {0.999039507598271565639221845698}
            ($B242,$6063,$FC0B,$3FEF),  {0.999517143498060754144094174832}
            ($C2DB,$D433,$FE03,$3FEF),  {0.999757685143858190853179678718}
            ($7DD7,$4924,$FF01,$3FEF),  {0.999878542763265115492174992820}
            ($02EA,$6E54,$FF80,$3FEF),  {0.999939170345979718170954192255}
            ($9808,$24EE,$FFC0,$3FEF),  {0.999969551213099238082632932628}
            ($A431,$0C59,$FFE0,$3FEF),  {0.999984764214906106441682774963}
            ($6F38,$0420,$FFF0,$3FEF),  {0.999992378292041011976937872238}
            ($C9E9,$0160,$FFF8,$3FEF),  {0.999996187869610113479689226412}
            ($C221,$0075,$FFFC,$3FEF),  {0.999998093508171675106856492967}
            ($4B24,$0027,$FFFE,$3FEF),  {0.999999046611581522115050842561}
            ($1BAA,$000D,$FFFF,$3FEF),  {0.999999523258215542816316664328}
            ($5F36,$8004,$FFFF,$3FEF),  {0.999999761613230822547897204945}
            ($753C,$C001,$FFFF,$3FEF),  {0.999999880801318439503223824847}
            ($7C74,$E000,$FFFF,$3FEF),  {0.999999940398892394628361403148}
            ($297F,$F000,$FFFF,$3FEF),  {0.999999970198856962834415133105}
            ($0DD6,$F800,$FFFF,$3FEF),  {0.999999985099231996568787661807}
            ($049D,$FC00,$FFFF,$3FEF),  {0.999999992549550484963515852741}
            ($018A,$FE00,$FFFF,$3FEF),  {0.999999996274753400108727527673}
            ($0083,$FF00,$FFFF,$3FEF),  {0.999999998137369418112186746561}
            ($002C,$FF80,$FFFF,$3FEF),  {0.999999999068682281453978627282}
            ($000F,$FFC0,$FFFF,$3FEF),  {0.999999999534340331454217514687}
            ($0005,$FFE0,$FFFF,$3FEF),  {0.999999999767169895951490822818}
            ($0002,$FFF0,$FFFF,$3FEF),  {0.999999999883584858046030472656}
            ($0001,$FFF8,$FFFF,$3FEF),  {0.999999999941792399045315923883}
            ($0000,$FFFC,$FFFF,$3FEF),  {0.999999999970896189529809522585}
            ($0000,$FFFE,$FFFF,$3FEF),  {0.999999999985448091433884763958}
            ($0000,$FFFF,$FFFF,$3FEF),  {0.999999999992724044606584750052}
            ($8000,$FFFF,$FFFF,$3FEF),  {0.999999999996362021933168755507}
            ($C000,$FFFF,$FFFF,$3FEF),  {0.999999999998181010843208735552}
            ($E000,$FFFF,$FFFF,$3FEF),  {0.999999999999090505380478878099}
            ($F000,$FFFF,$FFFF,$3FEF),  {0.999999999999545252676530873572}
            ($F800,$FFFF,$FFFF,$3FEF),  {0.999999999999772626333695897739}
            ($FC00,$FFFF,$FFFF,$3FEF),  {0.999999999999886313165324764872}
            ($FE00,$FFFF,$FFFF,$3FEF),  {0.999999999999943156582154653371}
            ($FF00,$FFFF,$FFFF,$3FEF),  {0.999999999999971578290908083389}
            ($FF80,$FFFF,$FFFF,$3FEF),  {0.999999999999985789145397627202}
            ($FFC0,$FFFF,$FFFF,$3FEF),  {0.999999999999992894572680008745}
            ($FFE0,$FFFF,$FFFF,$3FEF),  {0.999999999999996447286333736097}
            ($FFF0,$FFFF,$FFFF,$3FEF),  {0.999999999999998223643164778613}
            ($FFF8,$FFFF,$FFFF,$3FEF),  {0.999999999999999111821581692828}
            ($FFFC,$FFFF,$FFFF,$3FEF),  {0.999999999999999555910790614255}
            ($FFFE,$FFFF,$FFFF,$3FEF),  {0.999999999999999777955395229737}
            ($FFFF,$FFFF,$FFFF,$3FEF)); {0.999999999999999888977697589079}
var
  m: integer;
begin
  if n>53 then sfd_etaint := 1.0
  else if n>=0 then sfd_etaint := double(etahex[n])
  else begin
    m := 1-n;
    if odd(m) then sfd_etaint := 0
    else sfd_etaint := sfd_bernoulli(m)/m*exp2m1(m);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_eta(s: double): double;
  {-Return the Dirichlet eta function}
var
  z,t: double;
const
  c2 = 0.3268629627944929957e-1;
  c1 = 0.1598689037424309718;
  z0 = 2.5e-5;
begin
  if s>=55.0 then sfd_eta := 1.0
  else if (frac(s)=0.0) and (abs(s)<=MaxInt) then sfd_eta := sfd_etaint(round(s))
  else begin
    z := 1.0-s;
    if abs(z) <= z0 then sfd_eta := ln2 - (c2*z + c1)*z
    else if s >= -ETAEPS then sfd_eta := etam1pos(s)+1.0
    else if s <= -8.0 then begin
      {use eta(s) = (1-2^z)*zeta(s)}
      sfd_eta := sfd_zeta(s)*(1.0-exp2(z));
    end
    else begin
      {Use reflection formula for eta with s < 0}
      t := -exp2m1(z)/exp2m1(-s);
      t := t*cospi(0.5*z);
      t := t*sfd_gamma(z);
      t := t/power(pi,z);
      sfd_eta := t*etam1pos(z)+t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_etam1(s: double): double;
  {-Return Dirichlet eta(s)-1}
begin
  if s < -ETAEPS then sfd_etam1 := sfd_eta(s)-1.0
  else sfd_etam1 := etam1pos(s);
end;


{---------------------------------------------------------------------------}
function sfd_zeta(s: double): double;
  {-Return the Riemann zeta function at s, s<>1}
var
  a,b,sc: double;
  ig: integer;
begin
  {Ref: Boost [19], file zeta.hpp}
  if IsNanOrInfD(s) then begin
    if s=PosInf_d then sfd_zeta := 1.0
    else sfd_zeta := Nan_d;
    exit;
  end;
  if frac(s)=0.0 then begin
    if (s > -MaxBernoulli) and  (s <= MaxInt) then begin
      sfd_zeta := sfd_zetaint(round(s));
      exit;
    end;
  end
  else if abs(s) <= 5.2e-9 then begin
    {Small s, especially useful for s<0 to avoid reflection machinery}
    sfd_zeta := (-0.5) - LnSqrt2Pi*s;
    exit;
  end;
  sc := 1.0-s;
  if s<0.0 then begin
    if frac(0.5*s)=0.0 then sfd_zeta := 0.0
    else begin
      {compute Gamma(sc)/(2Pi)^sc}
      if sc>MaxGAMD-1 then begin
        {Not very accurate but does not overflow!  For s=-1753.5 this branch}
        {gives a rel. error of 2.0e-16 vs 1.8e-17 for the non-lngamma branch}
        a := sfd_lngammas(sc,ig);
        b := 2.0*sc*LnSqrt2Pi;
        {here Zeta(sc)=1}
        a := a-b;
        if a<Ln_MaxDbl then a := exp(a)
        else a := PosInf_d;
        a := ig*a;
      end
      else begin
        a := sfd_gamma(sc);
        b := power(TwoPi,-sc);
        a := a*b*zetap(sc,s);
      end;
      {here a = zeta(sc)*Gamma(sc)/(2Pi)^sc}
      sfd_zeta := 2.0*sinPi(0.5*s)*a;
    end;
  end
  else sfd_zeta := zetap(s,sc);
end;


{---------------------------------------------------------------------------}
function sfd_zeta1p(x: double): double;
  {-Return the Riemann zeta function at 1+x, x<>0}
begin
  {Ref: Boost [19], file zeta.hpp}
  if abs(x)<1.0 then sfd_zeta1p := zetap(1.0+x,-x)
  else sfd_zeta1p := sfd_zeta(1.0+x);
end;


{---------------------------------------------------------------------------}
function sfd_zetam1(s: double): double;
  {-Return Riemann zeta(s)-1, s<>1}
var
  t: double;
begin
  if s <= 2.0 then sfd_zetam1 := sfd_zeta(s) - 1.0
  else begin
    if s >= 100.0 then sfd_zetam1 := exp2(-s)
    else begin
      t := 0.5*exp2(s);
      sfd_zetam1 := (1.0 + etam1pos(s)*t) / (t - 1.0);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_zetaint(n: integer): double;
  {-Return zeta(n) for integer arguments, n<>1}
const
  znhex: array[2..53] of THexDblW = (
           ($07D3,$6253,$51A6,$3FFA),  {1.64493406684822643647241516665}
           ($0621,$04F0,$3BA0,$3FF3),  {1.20205690315959428539973816151}
           ($D848,$2AC7,$5132,$3FF1),  {1.08232323371113819151600369654}
           ($7CCE,$8ECA,$9741,$3FF0),  {1.03692775514336992633136548646}
           ($9245,$84C0,$4709,$3FF0),  {1.01734306198444913971451792979}
           ($CF39,$DA14,$2232,$3FF0),  {1.00834927738192282683979754985}
           ($6397,$6AF8,$10B3,$3FF0),  {1.00407735619794433937868523851}
           ($16B5,$F3D8,$0839,$3FF0),  {1.00200839282608221441785276923}
           ($5BB9,$E33A,$0412,$3FF0),  {1.00099457512781808533714595891}
           ($48B3,$31BE,$0206,$3FF0),  {1.00049418860411946455870228253}
           ($2CD3,$0A5B,$0102,$3FF0),  {1.00024608655330804829863799805}
           ($08BC,$AC9D,$0080,$3FF0),  {1.00012271334757848914675183653}
           ($CAD4,$392B,$0040,$3FF0),  {1.00006124813505870482925854510}
           ($97E2,$12F7,$0020,$3FF0),  {1.00003058823630702049355172851}
           ($DEB2,$064C,$0010,$3FF0),  {1.00001528225940865187173257149}
           ($39B4,$0218,$0008,$3FF0),  {1.00000763719763789976227360029}
           ($654E,$00B2,$0004,$3FF0),  {1.00000381729326499983985646165}
           ($611F,$003B,$0002,$3FF0),  {1.00000190821271655393892565696}
           ($C594,$0013,$0001,$3FF0),  {1.00000095396203387279611315205}
           ($95D6,$8006,$0000,$3FF0),  {1.00000047693298678780646311672}
           ($319B,$4002,$0000,$3FF0),  {1.00000023845050272773299000365}
           ($BB1E,$2000,$0000,$3FF0),  {1.00000011921992596531107306779}
           ($3E5A,$1000,$0000,$3FF0),  {1.00000005960818905125947961244}
           ($14C7,$0800,$0000,$3FF0),  {1.00000002980350351465228018606}
           ($06ED,$0400,$0000,$3FF0),  {1.00000001490155482836504123467}
           ($024F,$0200,$0000,$3FF0),  {1.00000000745071178983542949198}
           ($00C5,$0100,$0000,$3FF0),  {1.00000000372533402478845705482}
           ($0042,$0080,$0000,$3FF0),  {1.00000000186265972351304900640}
           ($0016,$0040,$0000,$3FF0),  {1.00000000093132743241966818287}
           ($0007,$0020,$0000,$3FF0),  {1.00000000046566290650337840730}
           ($0002,$0010,$0000,$3FF0),  {1.00000000023283118336765054921}
           ($0001,$0008,$0000,$3FF0),  {1.00000000011641550172700519776}
           ($0000,$0004,$0000,$3FF0),  {1.00000000005820772087902700890}
           ($0000,$0002,$0000,$3FF0),  {1.00000000002910385044497099687}
           ($0000,$0001,$0000,$3FF0),  {1.00000000001455192189104198425}
           ($8000,$0000,$0000,$3FF0),  {1.00000000000727595983505748101}
           ($4000,$0000,$0000,$3FF0),  {1.00000000000363797954737865119}
           ($2000,$0000,$0000,$3FF0),  {1.00000000000181898965030706595}
           ($1000,$0000,$0000,$3FF0),  {1.00000000000090949478402638893}
           ($0800,$0000,$0000,$3FF0),  {1.00000000000045474737830421540}
           ($0400,$0000,$0000,$3FF0),  {1.00000000000022737368458246526}
           ($0200,$0000,$0000,$3FF0),  {1.00000000000011368684076802278}
           ($0100,$0000,$0000,$3FF0),  {1.00000000000005684341987627587}
           ($0080,$0000,$0000,$3FF0),  {1.00000000000002842170976889302}
           ($0040,$0000,$0000,$3FF0),  {1.00000000000001421085482803161}
           ($0020,$0000,$0000,$3FF0),  {1.00000000000000710542739521085}
           ($0010,$0000,$0000,$3FF0),  {1.00000000000000355271369133712}
           ($0008,$0000,$0000,$3FF0),  {1.00000000000000177635684357912}
           ($0004,$0000,$0000,$3FF0),  {1.00000000000000088817842109308}
           ($0002,$0000,$0000,$3FF0),  {1.00000000000000044408921031438}
           ($0001,$0000,$0000,$3FF0),  {1.00000000000000022204460507980}
           ($0001,$0000,$0000,$3FF0)); {1.00000000000000011102230251411}
begin
  if n>53 then sfd_zetaint := 1.0
  else if n>1 then sfd_zetaint := double(znhex[n])
  else if n<0 then sfd_zetaint := sfd_bernoulli(1-n)/(n-1)+0.0  {avoid -0}
  else if n=0 then sfd_zetaint := -0.5
  else sfd_zetaint := Nan_d;
end;


{---------------------------------------------------------------------------}
function hz_a1(s,a: double; var ok: boolean): double;
  {-Return zetah(s,a) for a close to 1, s<>1 (normally s < 0)}
var
  x,t,f,z: double;
  n: integer;
const
  MaxIter = 100;
begin
  hz_a1 := 0.0;
  ok := false;
  {compute sum from NIST[30], 25.11.10}
  z := sfd_zeta(s);
  if IsNanOrInfD(z) then exit;
  x := 1.0-a;
  if x<>0.0 then begin
    n := 0;
    f := 1;
    repeat
      f := f*(s+n);
      inc(n);
      t := s+n;
      if (n>MaxIter) or (t=1.0) then exit;
      f := f*x/n;
      t := sfd_zeta(t);
      t := t*f;
      z := z + t;
      {do not use t for exit test because zeta(s+n) may be very small}
    until (abs(f)<=eps_d*abs(z)) and (n>2);
  end;
  ok := true;
  hz_a1 := z;
end;


{---------------------------------------------------------------------------}
function hurwitz_formula(s,a: double): double;
  {-Compute zetah with range reduction and Hurwitz formula, a > 0, s < 0}
var
  n,k: longint;
  a2,f,h,r,t,tol,sh,z: double;
  done: boolean;
const
  MaxIter = 10000;
begin
  r := 0.0;
  z := -s;
  f := frac(a);
  k := trunc(a);
  if f=0.0 then begin
    {Use strict Hurwitz formula with reduced a in (0,1].}
    {Empirically it is valid also for reduced frac(a)=0.}
    f := 1.0;
    k := k-1;
  end;
  a2  := 2.0*f;
  sh  := 2.0*frac(0.25*s);  {sh = s/2 mod 2}
  tol := 0.5*eps_d;
  {NIST [30],25.11.4}
  for n:=0 to k-1 do begin
    t := power(f+n,z);
    r := r + t;
  end;
  h := 0.0;
  {Perform summation only for non-trivial arguments}
  if (frac(a2)<>0.0) or (frac(sh)<>0.0) then begin
    {A. Erdelyi [50], Higher Transcendental Functions I, 1.10 (6),}
    {see also NIST [30], 25.11.9  with s <-> 1-s and cos <-> sin. }
    n := 1;
    z := s-1.0;
    done := false;
    repeat
      t := n*a2 + sh;
      if frac(t)<>0.0 then begin
        t := sinpi(t);
        f := power(n, z);
        t := t*f;
        h := h + t;
        done := (f < tol*(abs(h) + 1e-10));  {1e-10 addition to avoid exit if h ~ 0}
      end;
      inc(n);
    until done or (n>MaxIter);
    if (n>MaxIter) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  end;
  if h=0.0 then hurwitz_formula := -r
  else begin
    {compute = 2*h*Gamma(1-s)/(2Pi)^(1-s)}
    z := 1.0 - s;
    if (z>0.0) and (z<MAXGAMD) then begin
      f := sfd_gamma(z);
      t := power(TwoPi, -z);
      t := (t*f)*(2.0*h);
    end
    else begin
      f := sfd_lngamma(z);
      t := ln(abs(2.0*h)) + f - z*ln(TwoPi);
      if t>= ln_MaxDBL then t := PosInf_d
      else t := exp(t);
      if h<0.0 then t := -t;
    end;
    hurwitz_formula := t - r;
  end
end;


{---------------------------------------------------------------------------}
function bernpoly_intern(n: integer; x: double): double;
  {-Return the Bernoulli polynomial B_n(x), internal use: for small n only}
var
  i: integer;
  s,f: double;
  ref: boolean;
begin
  if n<0 then bernpoly_intern := 0
  else if n=0 then bernpoly_intern := 1.0
  else if n=1 then bernpoly_intern := x - 0.5
  else begin
    ref := false;
    if abs(x-1.0) <= 0.125 then begin
      {if x is near 1 use NIST[30], 24.4.3: B_n(1-x) = (-1)^n B_n(x)}
      x := 1.0 - x;
      ref := odd(n);
    end;
    {Compute sum NIST[30], 24.2.5}
    f := n;
    s := x - 0.5*n;
    for i:=2 to n do begin
      f := (n+1-i)*f/i;
      s := s*x;
      if i and 1 = 0 then s := s + sfd_bernoulli(i)*f;
    end;
    if ref then s := -s;
    bernpoly_intern := s;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_zetah(s,a: double): double;
  {-Return the Hurwitz zeta function zetah(s,a) = sum(1/(i+a)^s, i=0..INF), s<>1, a>0}
var
  eps,q,r,t,u,w,z: double;
  j,k,n: integer;
  ok: boolean;
const
  MaxIter = 40;
  eps1 = 5e-4;   {Threshold for hz_a1}
  S_HF = -8;     {Threshold for Hurzwitz formula}
label
  noconv;
begin
  {Initial reference: Cephes [7], file double\zeta.c}
  if IsNanOrInfD(s) or IsNanOrInfD(a) or (s=1.0) or (a<=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_zetah := NaN_d;
    exit;
  end
  else if s=0.0 then begin
    sfd_zetah := 0.5 - a;       {NIST [30], 25.11.13}
    exit;
  end
  else if a=0.5 then begin
    if s < 1024.0 then begin
      sfd_zetah := exp2m1(s)*sfd_zeta(s);   {NIST [30], 25.11.11}
      exit;
    end;
  end
  else if a=1.0 then begin
    sfd_zetah := sfd_zeta(s);   {NIST [30], 25.11.2}
    exit;
  end
  else if a=2.0 then begin
    sfd_zetah := sfd_zetam1(s); {NIST [30], 25.11.2/3}
    exit;
  end;

  if s>0.0 then begin
    {Quick checks: avoid overflow/underflow}
    t := ln(a);
    t := -s*t;
    if t>ln_MaxDbl then begin
      sfd_zetah := PosInf_d;
      exit;
    end
    else if t<ln_MinDbl then begin
      {here s > 1, 1/a^s ~ 0, avoid underflow in some cases:}
      {zetah ~ a^(1-s)/(s-1) + 0.5/a^s = (a/(s-1) + 0.5)/a^s}
      {      = a^(1-s)/(s-1) + 0.5*a^(1-s)/a}
      u := power(a, 1-s);
      sfd_zetah := u/(s-1.0) + 0.5*u/a;
      exit;
    end;
  end
  else begin
    if (a>1.0) then begin
      if -s*ln(a-1.0) >ln_MaxDbl then begin
        sfd_zetah := PosInf_d;
        exit;
      end;
    end;
    {Test if Bernoulli polynomials or Hurwitz formula can be used}
    if (s >= S_HF) and (frac(s)=0.0) then begin
      {s small negative integer: use 'naive' Bernoulli polynomial code}
      k := 1 - round(s);
      sfd_zetah := -bernpoly_intern(k,a)/k;  {NIST [30], 25.11.14}
      exit;
    end;
    if abs(a-1.0)<eps1 then begin
      {try expansion near a=1}
      sfd_zetah := hz_a1(s,a,ok);
      if ok then exit;
    end;
    if s < S_HF then begin
      sfd_zetah := hurwitz_formula(s,a);
      exit;
    end
  end;

  {Calculate z = zetah with the Euler-Maclaurin summation formula:}
  {See e.g. http://functions.wolfram.com/10.02.06.0020.01         }
  {z = sum1(s,a,n) + (a+n)^(1-s)/(s-1) - 1/2/(a+n)^s + sum2(s,a,n)}
  {with sum1(a,n) = sum(1/(a+k)^s, k=0..n) and                    }
  {sum2 = sum(B_2k/(2k)! * s(s+1)..(s+2k)/(a+n)^(s+2k+1),k=0..INF)}

  eps := 0.5*eps_d;
  if abs(s)<0.5 then n := 8
  else n := 9;

  {Avoid silly paranoid warning: r,w,t might not have been initialized}
  r := 0; w := 0; t := 0;

  {compute z = s1(a,n) = sum(1/(a+k)^s, k=0..n)}
  z := 0.0;
  for k:=0 to n do begin
    w := a+k;
    r := power(w,-s);
    z := z+r;
    if r < eps*z then begin
      sfd_zetah := z;
      exit;
    end
  end;

  {Add the two single terms z = z + (a+n)^(1-s)/(s-1) - 1/2/(a+n)^s}
  {Here w=(a+n), r=1/(a+n)^s}
  z := z + r*w/(s-1.0);
  z := z - 0.5*r;

  {Add the terms of sum2 until a term is < e*(partial result)}
  q := 1.0;
  k := 0;
  for j:=0 to MaxIter do begin
    q := q*(s+k);
    r := r/w;
    if j>NBoF then begin
      {Some problematic (s,a) values: Check decreasing terms in}
      {asymptotic expansion and increase convergence tolerance.}
      if j=NBOF+1 then eps := 1.5*eps_d;
      u := t;
      t := sfd_bernoulli(k+2);  {TP5 fix}
      t := t/sfd_fac(k+2);
      t := (q*t)*r;
      if abs(u)<=abs(t) then goto noconv;
    end
    else t := (q*double(BoFHex[j]))*r;
    z := z + t;
    if abs(t) < eps*abs(z) then begin
      sfd_zetah := z;
      exit;
    end;
    q := q*(s+k+1);
    r := r/w;
    inc(k,2);
  end;

noconv:
  {No convergence}
  sfd_zetah := z;
  {$ifdef debug}
    writeln('sfd_zetah - Euler-Maclaurin issue: ',s:21,' ',a:21);
  {$endif}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
end;


{---------------------------------------------------------------------------}
function bernpoly_rev(n: integer; x: double): double;
  {-Return the Bernoulli polynomial B_n(x), internal use: n>1, small x}
var
  i: integer;
  s,t,u,f: double;
begin
  {This code runs NIST[30] 24.2.5 backwards}
  s := sfd_bernoulli(n);
  f := 1.0;
  t := PosInf_d;
  for i:=n-1 downto 0 do begin
    f := f*x*(i+1)/(n-i);
    if (i and 1 = 0) or (i=1) then begin
      {B_i <> 0}
      u := t;
      t := s;
      s := t + sfd_bernoulli(i)*f;
      if (s=t) and (s=u) then begin
        {return if 3 partial sums have the same FP value}
        bernpoly_rev := s;
        exit;
      end;
    end;
  end;
  bernpoly_rev := s;
end;


{---------------------------------------------------------------------------}
function sfd_bernpoly(n: integer; x: double): double;
  {-Return the Bernoulli polynomial B_n(x), 0 <= n <= MaxBernoulli}
var
  a: double;
  neg: boolean;
begin
  if IsNanD(x) or (n<0) or (n>MaxBernoulli) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_bernpoly := Nan_d;
    exit;
  end;
  if n=0 then sfd_bernpoly := 1.0
  else if n=1 then sfd_bernpoly := x-0.5
  else if (x=0.0) or (x=1.0) then sfd_bernpoly := sfd_bernoulli(n)  {n<>1 !}
  else if x=0.5 then begin
    {B_n(1/2) = (2^(1-n)-1) B_n)}
    if n>65 then a := -1
    else a := exp2m1(1-n);
    sfd_bernpoly := a*sfd_bernoulli(n)
  end
  else if (x=-1.0) then begin
    {B_n(-1) := B_n + (-1)^n*n}
    if odd(n) then sfd_bernpoly := -n
    else sfd_bernpoly := sfd_bernoulli(n) + n;
  end
  else begin
    a := 1-x;
    if abs(a) > 0.5 then neg := false
    else begin
      {B_n(1-x) = (-1)^n B_n(x)}
      x := a;
      neg := odd(n);
    end;
    a := abs(x);
    if a > 12000.0*n then begin
      {Large x: First terms from NIST[30], 24.4.12 with x=0.5 and h = x-0.5}
      {B_n(x) = (x-0.5)^n - n(n-1)/24*(x-0.5)^(n-2) + 7/240*binomial(n,4)*(x-0.5)^(n-4) + ...}
      a := n; {This avoids integer overflow and FPC optimization error}
      a := sqr(x-0.5) - 0.25*a*(a-1.0)/SIXX;
      a := a*power(x-0.5,n-2);
    end
    else begin
      if n<=10 then begin
        {small degree}
        a := bernpoly_intern(n,x);
      end
      else if a < 0.125/n then begin
        {small argument}
        a := bernpoly_rev(n,x);
      end
      else begin
        {Loop count in Hurwitz formula is ~ trunc(abs(x)),}
        if x>=0.0 then a := sfd_zetah(1-n, x)
        else begin
          a := sfd_zetah(1-n, 1.0-x);
          if odd(n) then neg := not neg;
        end;
        a := -a*n;
      end;
    end;
    if neg then a := -a;
    sfd_bernpoly := a;
  end;
end;


{---------------------------------------------------------------------------}
type
  tphi_rec = record
               z,s,a,g: double;  {g=1/Gamma(s)}
             end;
  pphi_rec = ^tphi_rec;


{---------------------------------------------------------------------------}
function phif(x: double; param: pointer): double; {$ifdef BIT16} far;{$endif}
  {-LerchPhi integrand: x^(s-1)*exp(-ax)/(1-z*exp(-x)/Gamma(s)}
var
  d,p,f: double;
begin
  {NIST[30], 25.14.5}
  with pphi_rec(param)^ do begin
    p := (s-1.0)*ln(x);
    if p >= ln_MaxDbl then begin
      {power would overflow, use logarithms}
      p := exp(p - a*x - sfd_lngamma(s));
      d := 1.0 - z*exp(-x);
      phif := g*p/d;
    end
    else begin
      f := exp(-a*x);
      if f=0.0 then phif := 0.0
      else begin
        p := power(x,s-1);
        d := 1.0 - z*exp(-x);
        phif := g*p*f/d;
      end
    end;
  end;
end;


{---------------------------------------------------------------------------}
function phi_int(z,s,a: double): double;
  {-LerchPhi from integral representation, z < 0, 0 < s < MAX_GAMX, a > 0}
var
  res,aerr: double;
  ierr: integer;
  para: tphi_rec;
  neval: longint;
const
  eps  = 2e-15;
  eps2 = 5e-12;
begin
  phi_int := Nan_d;
  para.z := z;
  para.s := s;
  para.a := a;
  para.g := sfd_rgamma(s);
  if para.g=0.0 then begin
    {No reliable computation, Gamma(s) > Maxdouble}
    exit;
  end;
  {Use integral formula, see Wood[73], (1.1) or }
  {http://functions.wolfram.com/10.08.07.0001.01}
  sfd_intdei_p({$ifdef FPC_ProcVar}@{$endif}phif, @para, 0.0, eps, res, aerr, neval, ierr);
  if ierr=0 then phi_int := res
  else begin
    if (ierr=3) and (aerr <= eps2*abs(res)) then begin
      {Max. iterations with reduced accuracy}
      phi_int := res;
    end
  end;
end;


{---------------------------------------------------------------------------}
function lphi_aj(z, s, a: double): double;
  {-Lerch Phi via Aksenov/Jentschura convergence acceleration}
const
  IMAX  = 50;  {Max. iteration count for CNCT case }
  IMAX2 = 90;  {Max. iteration for direct summation}
var
  num,den,SAj : array[0..IMAX] of double;

  {---------------------------------------------------------------------------}
  function vwaj(j: integer): double;
     {-Compute the van Wijngaarden quantities A_j from b^j_k}
  var
    sum, bjk, z2ind: double;
    ind, two2k: double;
  begin
    sum  := 0.0;
    two2k:= 1.0;
    {Sum b^j_k's over k}
    repeat
      {Index for the term of the original series}
      ind   := two2k*(j+1)-1.0;
      z2ind := power(z,ind);
      bjk   := two2k*z2ind/power(a+ind,s);
      sum   := sum + bjk;
      two2k := 2.0*two2k;
    until abs(bjk) <= eps_d*abs(sum);
    vwaj := sum;
  end;

var
  j,i,k,sign: integer;
  sn,eps0,eps,skn,skn0,omega,fk,est,t,x: double;
begin

  {Based on the C code lerchphi.c by S.V. Aksenov and U.D. Jentschura.}
  {The original specs are heavily modified: s and a must be positive, }
  {eps=0 test is moved to avoid division by zero, num[0]=den[0]=0 is  }
  {allowed on the 1st iteration etc. The code for |z|<=0.5 is not used}
  {by the driver function, but is left here for testing and comparing.}

  {The CNCT is described in [40], C source and lerchphi user guide is }
  {available from http://aksenov.freeshell.org/lerchphi.html          }

  lphi_aj := 0.0;
  x := abs(z);

  {$ifdef debug}
    if (z >= 1.0) or (z < -1.0) or (a <= 0.0) or ( s < -1.0)  then begin
      {Return NaN, this should not happen if called by sfd_lerch}
      lphi_aj := Nan_d;
      exit;
    end;
  {$endif}

  if x <= MinDouble then begin
    {return first term of series}
    lphi_aj := power(a, -s);
    exit;
  end;

  {sn denotes current partial sum of defining series:
    z >  0.5: sn is partial sum S_n of the van Wijngaarden transformed series.
    z <= 0.5: sn is the partial sum of the power series defining LerchPhi.
  skn0 and skn denote successive partial sums S^k_n that are same as sn
  in case of direct summation and delta-transformed in case of CNCT.
  eps0 and eps denote successive differences between partial sums S^k_n.}

  eps0 := 0.0;
  skn  := 0.0;
  skn0 := skn;  {Avoid hint/warning: Keep Delphi-64 happy and use skn}
  sn   := 0.0;

  {omega is the next term of a partial sum of defining power series for direct}
  {summation, of van Wijngaarden transformed series for CNCT, and also becomes}
  {a remainder estimate in the delta transformation in CNCT.}

  {For z<=0.5 van Wijngaarden transformation is not used hence no calls to aj()}
  if z <= 0.5 then omega := power(a, -s)
  else omega := vwaj(0);

  if IsInfD(omega) then begin
    lphi_aj := omega;
    exit;
  end;
  i    := -1;
  sign := -1;

  {Main loop: iterations for S^k_n}
  repeat
    {i points to current iterate}
    inc(i);
    sign := -sign;
    sn   := sn + omega;

    {Next term: omega}
    if z <= 0.5 then omega := z*power((a+i)/(a+i+1), s)*omega
    else begin
      {CNCT (z > 0.5) case}
      SAj[i] := sign*omega;
      if odd(i) then omega := vwaj(i+1)
      else begin
        j := i div 2;
        omega := 0.5*(SAj[j] - power(z,j)/power(a+j,s));
      end;
      omega := -sign*omega;
    end;

    if x <= 0.5 then begin
      {Direct summation case: store current}
      skn := sn;
    end
    else begin
      {CNCT case: Make sure omega is representable machine number}
      if abs(omega) <= MinDouble then begin
        num[i] := 0;
        den[i] := 0;
        if i=0 then begin
          {WE: Allow zeros for i=0}
          continue;
        end
        else begin
          if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
          exit;
        end;
      end
      else begin
        t := 1.0/omega;
        {Last terms in sums of numerator and denominator of i-th partial sum}
        num[i] := sn*t;
        den[i] := t;
      end;
      {Recurrence computation of numerator and denominator of a S_k^n}
      k  := i;
      fk := (k+1)*k;
      for j:=i-1 downto 0 do begin
        t := fk/(k*(k+1));
        num[j] := num[j+1] - t*num[j];
        den[j] := den[j+1] - t*den[j];
        inc(k);
      end;
      {Current approximation of the sum S_k^n}
      skn := num[0]/den[0];
    end; {CNCT case}
    lphi_aj := skn;
    eps := abs(skn - skn0);
    {Check the three termination criteria}

    {Successive iterates skn are the same. WE: must be   }
    {checked first, otherwise division by zero may occur!}
    if eps=0.0 then exit;

    {|est/skn| is less than the requested accuracy
    (est is a remainder estimate)}
    if (i>0) and (eps < eps0) then begin
      if x > 0.5 then begin
        est := eps/eps0;
        est := 2.0/est/(1.0-est)*eps;
      end
      else est := 2.0*power(x, i+1)/power(a+i+1, s);
      if abs(est) < 8.0*eps_d*abs(skn) then exit;
    end;

    {Maximum number of iterations is exceeded}
    if i>=imax then begin
      if (x>0.5) or (i>=IMAX2) then begin
        if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
        exit;
      end;
    end;

    {Go on to the next iteration}
    skn0 := skn;
    eps0 := eps;

  until false;

end;


{---------------------------------------------------------------------------}
function lphi_aea(z,s,a: double): double;
  {-Lerch Phi, asymptotic expansion for large a>0}
var
  n: integer;
  e,r,t,x,y: double;
const
  NMAX=1000;
begin
  {This function uses the asymptotic expansion given in [41], Theorem 1.}
  {Note: The simple convergence test is not accurate for |z| close to 1,}
  {in this case a test using the error term in [41] should be used.     }
  n := 1;
  r := 1.0/(1.0-z);
  x := power(a,-s);
  if x<>0 then begin
    y := 1.0;
    e := 0.5*eps_d;
    repeat
      {y holds the terms (s)_n/(n!*a^n)}
      y := y*(s/n)/a;
      t := sfd_polylog(-n,z);
      t := y*t;
      if odd(n) then r := r - t
      else r := r + t;
      inc(n);
      s := s+1.0;
    until (n>NMAX) or (abs(t) < e*abs(r));
  end;
  lphi_aea := x*r;
  if (n>NMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
end;


{---------------------------------------------------------------------------}
function lphi_sum(z,s,a: double): double;
  {-Lerch Phi, direct summation sum(z^n/(n+a)^s, n=0...)}
var
  n: integer;
  r,t,x,e: double;
const
  NMAX=1000;
begin
  n := 0;
  x := 1.0; {x = z^n}
  r := 0.0;
  e := 0.5*eps_d;
  s := -s;
  repeat
    t := x*power(a+n,s);
    r := r + t;
    n := n+1;
    x := x*z;
  until (n>NMAX) or (abs(t) <= e*abs(r));
  lphi_sum := r;
  if (n>NMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
end;


{---------------------------------------------------------------------------}
function sfd_lerch(z,s,a: double): double;
  {-Return the Lerch transcendent Phi(z,s,a) = sum(z^n/(n+a)^s, n=0..INF), z <= 1, s >= -1, a >= 0}
var
  x,t: double;
begin
  {Phi(z,s,a) = sum(z^n/(n+a)^s, n=0...)}
  if IsNanOrInfD(z) or IsNanOrInfD(s) or IsNanOrInfD(a)
     or (z>1.0) or (s<-1.0) or (a<0.0) then
  begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_lerch := Nan_d;
    exit;
  end;

  if s=0.0 then begin
    {geometric series Phi = sum(z^n)}
    sfd_lerch := 1.0/(1.0-z);
    exit;
  end
  else if z=0.0 then begin
    sfd_lerch := power(a,-s);
    exit;
  end
  else if z=1.0 then begin
    {continuation from s>1}
    sfd_lerch := sfd_zetah(s,a);
    exit;
  end
  else if a=0.0 then begin
    {http://functions.wolfram.com/10.06.03.0049.01}
    sfd_lerch := sfd_polylogr(s,z);
    exit;
  end;

  x := abs(z);

  if s<0.0 then begin
    {Note: For s<0 no asymptotic expansion for large a}
    if (x>0.5) or (z >= 0.875) then begin
      {accelerated summation}
      sfd_lerch := lphi_aj(z,s,a);
    end
    else begin
      {direct summation}
      sfd_lerch := lphi_sum(z,s,a);
    end;
    exit;
  end;

  {Here s > 0}
  if z <= -1.0 then begin
    if z=-1.0 then begin
      {Although Aksenov/Jentschura do not allow z=-1, this is a classical}
      {series transformation scenario and their code gives good results. }
      {The alternative Phi = (zetah(s,a/2)-zetah(s,(a+1)/2))/2^s is very }
      {susceptible to truncation errors and gives 0 for a>2/eps_x!       }
      sfd_lerch := lphi_aj(z,s,a);
    end
    else begin
      {Use integral NIST[30], 25.14.5 for x<-1, a>0, s>0}
      sfd_lerch := phi_int(z,s,a)
    end;
    exit;
  end;

  t := s*ln(a);
  if t > ln_MaxDbl then begin
    {1/a^s = 0 accurate to double precision}
    sfd_lerch := 0;
  end
  else if (x <= 0.3) and (a >= 100.0) and (s <= 100.0) then begin
    {Use asymptotic expansion for large a, here only if s and x are not too large}
    sfd_lerch := lphi_aea(z,s,a);
  end
  else if (s >= 30.0) or (t < -43.7) then begin
    if (x>0.5) and (s < 0.1*a) then begin
      {accelerated summation because direct summation needs many iterations}
      sfd_lerch := lphi_aj(z,s,a);
    end
    else begin
      {direct summation if s is large enough or 1/a^s > 1/eps_d}
      sfd_lerch := lphi_sum(z,s,a);
    end;
  end
  else if (1.0-x)*a >= 50.0  then begin
    {Remaining case for asymptotic expansion for large a}
    sfd_lerch := lphi_aea(z,s,a);
    exit;
  end
  else begin
    if x<=0.5 then begin
      {direct summation, lphi_aj would do the same but a lot more complicated}
      sfd_lerch := lphi_sum(z,s,a);
    end
    else begin
      {accelerated summation}
      sfd_lerch := lphi_aj(z,s,a);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_dbeta(s: double): double;
  {-Return the Dirichlet beta function sum((-1)^n/(2n+1)^s, n=0..INF)}
var
  t: double;
const
  d0: THexDblW = ($4317,$ED7E,$0FE1,$3FD9); {0.39159439270683677647194534689911102809021011577}
begin
  {http://en.wikipedia.org/wiki/Dirichlet_beta_function   }
  {http://mathworld.wolfram.com/DirichletBetaFunction.html}
  {This is also known as the Catalan beta function        }

  if abs(s) <= 1e-10 then sfd_dbeta := 0.5 + double(d0)*s
  else if s<0.0 then begin
    if frac(0.5*s)=-0.5 then sfd_dbeta := 0.0
    else begin
      t := cosPi(0.5*s);
      s := 1.0 - s;
      t := sfd_dbeta(s)*t;
      sfd_dbeta := sfd_gamma(s)*t*power(Pi_2,-s);
    end;
  end
  else if s>=37.8 then sfd_dbeta := 1.0
  else if s>=25.8 then sfd_dbeta := 1.0 - exp3(-s)
  else if s>=21.3 then sfd_dbeta := 1.0 - exp3(-s) + exp5(-s)
  else begin
    sfd_dbeta := sfd_lerch(-1.0,s,0.5)*exp2(-s);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_dlambda(s: double): double;
  {-Return the Dirichlet lambda function sum(1/(2n+1)^s, n=0..INF), s<>1}
begin
  {http://mathworld.wolfram.com/DirichletLambdaFunction.html}
  sfd_dlambda := -exp2m1(-s)*sfd_zeta(s);
end;


{---------------------------------------------------------------------------}
function sfd_lchi(s, x: double): double;
  {-Return Legendre's Chi-function chi(s,x); s>=0; real part if |x| > 1}
var
  a,b,z: double;
begin
  if IsNanOrInfD(s) or IsNanOrInfD(x) or (s<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_lchi := Nan_d;
    exit;
  end;
  z := abs(x);
  if z=0.0 then sfd_lchi := 0.0
  else if s=0.0 then sfd_lchi := x/((1.0+x)*(1.0-x))
  else if s=1.0 then begin
    if z<1.0 then sfd_lchi := arctanh(x)
    else begin
      {Compute real part of arctanh(z), see manual for carctanh}
      a := 2.0/(1.0-z);
      b := 0.25*ln1p(z*sqr(a));
      sfd_lchi := b*isign(x);
    end;
  end
  else if z <= 1.0 then begin
    {chi(s,x) = sum(x^(2n+1)/(2n+1)^s, n=0..INF)}
    {chi(s,x) = 2^(-s)*x*Phi(x^2,s,1/2) = 0.5*(Li_s(x)-Li_s(-x))}
    {see e.g. http://en.wikipedia.org/wiki/Legendre_chi_function}
    {or  http://mathworld.wolfram.com/LegendresChi-Function.html}
    if s>=21.3 then begin
      if s>=37.8 then sfd_lchi := x
      else begin
        z := sqr(x);
        if s>=25.8 then sfd_lchi := x*(1.0 + z*exp3(-s))
        else sfd_lchi := x*(1.0 + z*(exp3(-s) + z*exp5(-s)));
      end;
    end
    else begin
      sfd_lchi := sfd_lerch(x*x,s,0.5)*x*exp2(-s);
    end;
  end
  else begin
    {|x| > 1}
    a := sfd_polylogr(s, z);
    b := sfd_polylogr(s,-z);
    sfd_lchi := 0.5*(a-b)*isign(x);
  end;

end;


{---------------------------------------------------------------------------}
function sfd_beint(s, x: double): double;
  {-Return the Bose-Einstein integral of real order s >= -1, real part if x > 1}
const
  XL = 50;
begin
  if IsNanOrInfD(s) or IsNanOrInfD(x) or (s < -1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_beint := NaN_d;
    exit;
  end;
  if x > XL then begin
    {Use asymptotic expansion of polylog without computing exp(x))}
    sfd_beint := liae_pos(s+1.0,x)
  end
  else sfd_beint := sfd_polylogr(s+1.0,exp(x));
end;


{---------------------------------------------------------------------------}
function fd_asymp_exp_real(q, x: double): double;
  {-Fermi-Dirac asymptotic expansion for real order q with q > -1, x >> q}
var
  p,s,t,r,f,z: double;
  k: integer;
begin
  {Asymptotic expansion Goano [47], (9)}
  {$ifdef debug}
    if q + 10 > x then begin
      {The AE is valid only if x >> q. Normally used only}
      {for q=-1/2, 1/2, (3/2, 5/2) and x > ~40}
      fd_asymp_exp_real := PosInf_d;
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
    r := p*sfd_etaint(k)*f;
    if abs(r) <= abs(t) then begin
      {terms of AE are still decreasing}
      t := r;
      s := s+t;
      p := p*((k-1)-q)*(k-q);
    end
    else t := 0.0;
  until abs(t) <= eps_d*abs(s);
  f := power(x,q+1.0)*sfd_rgamma(q+2.0);
  fd_asymp_exp_real := f*s;
end;


{---------------------------------------------------------------------------}
function sfd_fermi_dirac(n: integer; x: double): double;
  {-Return the integer order Fermi-Dirac integral F_n(x) = 1/n!*integral(t^n/(exp(t-x)+1), t=0..INF)}
var
  z: double;
begin
  {Basic formula is F_n(x) = -polylog(n+1,-exp(x)}
  if x < -ln_MaxDbl then begin
    sfd_fermi_dirac := 0.0;
  end
  else if x > ln_MaxDbl then begin
    {asymptotic Goano [47], (9)}
    if n>0 then sfd_fermi_dirac := fd_asymp_exp_real(n, x)
    else case n of
        0: sfd_fermi_dirac := x;     {ln(1+e^x)   -> x}
       -1: sfd_fermi_dirac := 1.0;   {e^x/(1+e^x) -> 1}
      else sfd_fermi_dirac := 0.0;
    end;
  end
  else begin
    z := -exp(x);
    sfd_fermi_dirac := -sfd_polylog(n+1, z);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_fdr(s,x: double): double;
  {-Return the Fermi-Dirac functions of real order s >= -1}
var
  ae,cs,z,y: double;
begin
  if IsNanOrInfD(s) or IsNanOrInfD(x) or (s < -1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_fdr := NaN_d;
    exit;
  end;
  if frac(s)=0.0 then begin
    sfd_fdr := sfd_fermi_dirac(round(s),x);
  end
  else begin
    if x < ln_MinDbl then sfd_fdr := 0.0
    else if x >= 30.0 + abs(s) then begin
      ae := fd_asymp_exp_real(s,x);
      cs := cospi(s);
      if (cs <> 0.0) and (x >= -ln_MaxDbl) then begin
        y  := exp(-x);
        z  := -sfd_polylogr(s + 1.0, -y);
        ae := ae + z*cs;
      end;
      sfd_fdr := ae;
    end
    else begin
      y := exp(x);
      z := s+1.0;
      if x < z*ln2 - 37.0 then sfd_fdr := y
      else sfd_fdr := -sfd_polylogr(z, -y);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_fdp05(x: double): double;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}
begin
  sfd_fdp05 := sfd_fdr(0.5, x);
end;


{---------------------------------------------------------------------------}
function sfd_fdp15(x: double): double;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}
begin
  sfd_fdp15 := sfd_fdr(1.5, x);
end;


{---------------------------------------------------------------------------}
function sfd_fdp25(x: double): double;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}
begin
  sfd_fdp25 := sfd_fdr(2.5, x);
end;


{---------------------------------------------------------------------------}
function sfd_fdm05(x: double): double;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}
begin
  sfd_fdm05 := sfd_fdr(-0.5, x);
end;


{---------------------------------------------------------------------------}
function sfd_harmonic(x: double): double;
  {-Return the harmonic function H(x) = psi(x+1) + EulerGamma}
const
  nch = 20;
  chah: array[0..nch-1] of THexDblW = (    {chebyshev((Psi(x+1)+gamma)/x, x=-1/4..1/4, 0.2e-20);}
          ($4D5F,$8EB2,$E29F,$400A),  {+3.36065589410433977112951870486     }
          ($7CEC,$1B06,$0D61,$BFD4),  {-0.313316608802735484720854737533    }
          ($8810,$97A0,$65A0,$3FA2),  {+0.359316048709531508723414593246e-1 }
          ($3029,$759F,$F05E,$BF71),  {-0.437962434984843410853415224307e-2 }
          ($EAF5,$9C6A,$E6A2,$3F41),  {+0.546292686372267141546977073957e-3 }
          ($9261,$A627,$09E1,$BF12),  {-0.688117957348770917553696097164e-4 }
          ($7B52,$CCBB,$41EB,$3EE2),  {+0.870585645123897195412611399047e-5 }
          ($E225,$62DD,$8458,$BEB2),  {-0.110369763769828900468418914839e-5 }
          ($BA65,$0424,$CC67,$3E82),  {+0.140059343742157854915062263457e-6 }
          ($6CD6,$EE71,$17D9,$BE53),  {-0.177818994219795504449052141370e-7 }
          ($5ADF,$10A2,$659E,$3E23),  {+0.225810137695486223694241640848e-8 }
          ($F21E,$2F91,$B52F,$BDF3),  {-0.286785525509693105456693572379e-9 }
          ($60E6,$99FE,$064E,$3DC4),  {+0.364246022659874814089675465140e-10}
          ($9338,$DF9D,$58DF,$BD94),  {-0.462640992702404867971084536579e-11}
          ($4581,$93F6,$ACD7,$3D64),  {+0.587623516675841187923320281923e-12}
          ($0367,$CDDD,$0232,$BD35),  {-0.746374969291138217069111742907e-13}
          ($8D00,$EF8E,$58F2,$3D05),  {+0.948017526953988879040131837225e-14}
          ($3174,$8CDC,$B11B,$BCD5),  {-0.120413826871905020093317205030e-14}
          ($691B,$5EB8,$0AB1,$3CA6),  {+0.152945496208458532222222222222e-15}
          ($6621,$BD17,$65B9,$BC76)); {-0.194266177768451363581139754366e-16}
       (* ($DC15,$5BF2,$C23A,$3C46),  {+0.246750352250719046709039799083e-17}
          ($2AF8,$29E9,$2039,$BC17),  {-0.313414011663104260448782764313e-18}
          ($CEB4,$3FB3,$7FBC,$3BE7),  {+0.398087971944463065142369984751e-19}
          ($BBEB,$D80F,$E0C9,$BBB7)); {-0.505638009703096728955241279549e-20} *)
var
  cha: array[0..nch-1] of double absolute chah;
var
  y: double;
begin
  if IsNanOrInfD(x) then begin
    if x=PosInf_d then sfd_harmonic := PosInf_d
    else sfd_harmonic := NaN_d;
  end
  else if frac(x)=0.0 then begin
    if x<0.0 then begin
      {$ifopt R+}
        if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      {$endif}
      sfd_harmonic := Nan_d;
    end
    else if x<=12.0 then begin
      {sfd_psi uses a similar loop, avoid (y - EulerGamma) + EulerGamma}
      y := 0.0;
      while x>=1.0 do begin
        y := y + 1.0/x;
        x := x - 1.0;
      end;
      sfd_harmonic := y;
    end
    else sfd_harmonic := sfd_psi(x+1.0) + EulerGamma;
  end
  else begin
    if abs(x)<=0.25 then begin
      y := CSEvalD(4.0*x, cha, nch);
      sfd_harmonic := y*x;
    end
    else begin
      y := sfd_psi(x+1.0);
      sfd_harmonic :=  y + EulerGamma;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function harm2core(x,r: double): double;
  {-Core routine for H(x,r), x >= -1, no checks, r<>0,1,-1}
  { The alternative computations use Hurwitz zeta or series for |x|<0.5}
var
  a,s,t: double;
  k: integer;
begin
  {Get threshold value for series}
  if r > 1.0 then t := 0.125 else t := 0.5;
  if abs(x)<=t then begin
    {http://functions.wolfram.com/06.17.06.0018.01 for z0=0, note that in}
    {.../06.17.06.0002.01 there are the restriction Re(r) > 1 and |x| < 1}
    a := r*x;
    s := a*sfd_zeta1p(r);
    k := 1;
    repeat
      a := (r+k)*x*a;
      k := k+1;
      a := -a/k;
      t := s;
      s := t + a*sfd_zeta(r+k);
    until t=s;
    harm2core := s;
  end
  else begin
    s := sfd_zeta(r);
    a := sfd_zetah(r, x+1.0);
    harm2core := s - a;
  end;
end;


{---------------------------------------------------------------------------}
function sfd_harmonic2(x,r: double): double;
  {-Return the generalized harmonic function H(x,r) = zeta(r)-zetah(r,x+1); x >= -1}
const
  XREC = 20.0;  {Max x for recursive call}
var
  a,s,t: double;
  n: integer;
begin
  if IsNanOrInfD(x) or IsNanOrInfD(r) or (x < -1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfd_harmonic2 := NaN_d;
    exit;
  end;

  if (frac(x)=0.0) and (x <= 6.0)  then begin
    {Direct summation with exp<x> functions for small integer x}
    n := round(x);
    if n<2 then sfd_harmonic2 := n
    else if n=2 then sfd_harmonic2 := 1.0 + exp2(-r)
    else begin
      a := exp2(-r);
      s := exp3(-r);
      case n of
          3: sfd_harmonic2 := 1.0 + a + s;
          4: sfd_harmonic2 := 1.0 + a + s + a*a;
          5: sfd_harmonic2 := 1.0 + a + s + a*a + exp5(-r);
        else sfd_harmonic2 := 1.0 + a + s + a*a + exp5(-r) + s*a;
      end;
    end;
  end
  else if (r<=1.0) and (r > -MaxInt) and (frac(r)=0.0)  then begin
    {r integer <= 1}
    n := 1 - round(r);
    if n=0 then begin
      {r=1, ordinary harmonic number function}
      sfd_harmonic2 := sfd_harmonic(x);
    end
    else begin
      {http://functions.wolfram.com/06.17.27.0005.01}
      a := sfd_bernoulli(n);
      if n and 1 = 0 then a := -a;
      s := sfd_bernpoly(n, x + 1.0);
      sfd_harmonic2 := (s+a)/n;
    end;
  end
  else if (-1.0 < x) and (x < -0.5) then begin
    {http://functions.wolfram.com/06.17.16.0004.01}
    {This is a critical range with possible cancellation: e.g.}
    {the relative error for H(-0.75, -3.5) is about 3.8E-18.}
    t := x + 1.0;
    s := harm2core(t,r);
    a := power(t,-r);
    sfd_harmonic2 := s - a;
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
    sfd_harmonic2 := s + a;
  end
  else begin
    {No special case or recursion}
    sfd_harmonic2 := harm2core(x,r);
  end;
end;


{---------------------------------------------------------------------------}
function sfd_llci(x: double): double;
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
  nl1=13;
  arlob1 : array[0..nl1-1] of double = (
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
             0.3878e-16);
           { 0.240e-17,
             0.15e-18,
             0.1e-19);}
  nl2=9;
  arlob2 : array[0..nl2-1] of double = (
             2.034594180361328511,
             0.1735185882027407681e-1,
             0.5516280426090521e-4,
             0.39781646276598e-6,
             0.369018028918e-8,
             0.3880409214e-10,
             0.44069698e-12,
             0.527674e-14,
             0.6568e-16);
           { 0.84e-18,
             0.1e-19);}
const
  xlow1 = 8.719671245e-17;       {Pi/4*eps/2}
  xlow2 = 4.71216091e-8;         {sqrt(10*eps)}
  xlow3 = 6.32202727e-8;         {sqrt(18*eps)}
  xhigh = 4.503599627370496e15;  {1/eps}
var
  xr,fval,xcub,t,npi: double;
  sx: integer;
  gtpi2: boolean;
begin

  {Ref: MISCFUN [22], function LOBACH}
  sx := isign(x);
  xr := abs(x);

  if xr >= xhigh then begin
    {Contribution from reduced x is negligible compared to x*ln(2)}
    {Note that Miscfun returns an error for these cases}
    sfd_llci := ln2*x;
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
      t := CSEvalD(t, arlob1, nl1);
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
      t := CSEvalD(t, arlob2, nl2);
      fval := (t - ln(xr))*xr;
    end;
    fval := (0.5*lobpi1 - fval) + lbpb22;
  end;

  {Compute value for argument in [pi/2,pi]}
  if gtpi2 then fval := (lobpi1 - fval) + lobpi2;

  {Add quasi-periodic contribution}
  if npi > 0.0 then fval := (fval + npi*lobpi2 ) + npi*lobpi1;

  {Adjust sign}
  sfd_llci := sx*fval;
end;


{---------------------------------------------------------------------------}
function sfd_llsi(x: double): double;
  {-Return the Lobachevski function Lambda(x) = integral(-ln(|2sin(t)|), t=0..x)}
begin
  {http://mathworld.wolfram.com/LobachevskysFunction.html}
  {Lambda(x) = 0.5*cl_2(2x) = Im(polylog(2, exp(2ix)))/2}
  sfd_llsi := 0.5*sfd_cl2(2.0*x);
end;

end.
