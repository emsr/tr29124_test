unit sfZeta;

{Common code for special functions: Zeta functions and polylogarithms}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}

(*************************************************************************

 DESCRIPTION   :  Common code for Zeta functions and polylogarithms

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  ---

 REFERENCES    :  References used in this unit, main index in amath_info.txt/references

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
 1.11.00  03.02.12  W.Ehrhardt  New unit with parts of sfgamma, sfmisc, sfbasic
 1.11.01  04.02.12  we          sfc_lerch
 1.11.02  05.02.12  we          sfc_lerch for z=-1: use lphi_aj instead of zetah
 1.11.03  06.02.12  we          sfc_dbeta
 1.11.04  06.02.12  we          sfc_dlambda
 1.11.05  07.02.12  we          sfc_lchi
 1.11.06  08.02.12  we          sfc_polylogr

 1.12.00  13.03.12  we          sfc_lchi(1,x) = arctanh(x)
 1.12.01  20.03.12  we          sfc_zetah with BoFHex

 1.13.00  09.06.12  we          sfc_polylog uses sfc_taylor
 1.13.01  09.06.12  we          sfc_fermi_dirac
 1.13.02  10.06.12  we          fix n>x in fd_asymp_exp_int
 1.13.03  11.06.12  we          sfc_zeta uses sfc_zetaint
 1.13.04  12.06.12  we          sfc_etaint
 1.13.05  13.06.12  we          sfc_fdm05, sfc_fdp05

 1.15.00  16.02.13  we          improved sfc_zetah

 1.16.00  20.03.13  we          fix sfc_zetah for j>NBoF
 1.16.01  27.03.13  we          sfc_zetah for s<0
 1.16.02  28.03.13  we          sfc_trilog
 1.16.03  29.03.13  we          TP5 fix for sfc_zeta
 1.16.04  29.03.13  we          sfc_zetah for a near 1 and s<0

 1.17.00  10.04.13  we          sfc_fdp15

 1.18.00  19.05.13  we          Prevent some wrong compiler optimizations for div by 3

 1.20.00  16.08.13  we          sfc_zeta for small s
 1.20.01  16.08.13  we          sfc_polylog with sfc_zetaint

 1.21.00  07.09.13  we          Improved sfc_zetah/hurwitz_formula
 1.21.01  11.09.13  we          Bernoulli polynomials sfc_bernpoly

 1.23.00  26.12.13  we          sfc_fdp25

 1.25.00  01.05.14  we          sfc_harmonic
 1.25.01  03.05.14  we          sfc_harmonic2

 1.26.00  27.05.14  we          Removed redundancy in sfc_harmonic2
 1.26.01  28.05.14  we          sfc_harmonic with nch = 22
 1.26.02  29.05.14  we          polylogneg with array of single
 1.26.03  30.05.14  we          sfc_llci, sfc_llsi
 1.26.04  31.05.14  we          Improved harm2core

 1.27.00  15.06.14  we          Fermi/Dirac, Lobachewsky, harmonic functions move to sfZeta2

 1.28.00  29.07.14  we          Improved and expanded primezeta sfc_pz
 1.28.01  05.08.14  we          sfc_zeta for s < -MaxGAMX

 1.33.00  17.05.15  we          sfc_lerch for s >= -1
 1.33.01  18.05.15  we          sfc_polylogr for s >= -1
 1.33.02  18.05.15  we          Improved sfc_etam1 (adjusted Borwein constants)
 1.33.03  24.05.15  we          sfc_lerch(1,s,a) for s<>1 and special case z=0

 1.34.00  18.07.15  we          eta near s=1

 1.35.00  25.08.15  we          Use THREE & SIXX to avoid 'optimization'
 1.35.01  25.08.15  we          sfc_bernpoly: avoid integer overflow and FPC optimization error

 1.40.00  29.06.17  we          Code removed for TP5-TP6, TPW1-D1

 1.43.00  02.12.17  we          Suppress warnings: Local variable does not seem to be initialized

 1.44.00  04.02.18  we          Rewrite polylog, returns real part for x>1,n>1
 1.44.01  04.02.18  we          sfc_zetaint(1) = NaN
 1.44.02  05.02.18  we          sfc_polylog(n,x) for x > 1
 1.44.03  07.02.18  we          sfc_polylogr for x < -1
 1.44.04  09.02.18  we          Merged with (old) unit sfZeta2
 1.44.05  09.02.18  we          sfc_fdr (Fermi-Dirac real order)
 1.44.06  12.02.18  we          sfc_polylogr(s,x)=x for large s

 1.45.00  24.02.18  we          sfc_lerch for z < -1
 1.45.01  25.02.18  we          fix NaN check for r in sfc_harmonic2
 1.45.02  26.02.18  we          fix overflow for LerchPhi integrand
 1.45.03  26.02.18  we          sfc_ti (inverse tangent integral of real order)
 1.45.04  01.03.18  we          improved lerch integration
 1.45.05  03.03.18  we          sfc_polylogr for 1 < x <= 256
 1.45.06  05.03.18  we          sfc_lchi for |x| > 1
 1.45.07  05.03.18  we          Bose-Einstein integral sfc_beint
 1.45.08  06.03.18  we          sfc_polylogr for x > 256

 1.49.00  26.07.18  we          Avoid D25 warning

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


function sfc_cl2(x: extended): extended;
  {-Return the Clausen function: integral(-ln(2*|sin(t/2)|),t=0..x) = Im(Li_2(exp(ix)))}

function sfc_ti(s,x: extended): extended;
  {-Return the inverse tangent integral of order s >= 0}

function sfc_ti2(x: extended): extended;
  {-Return the inverse tangent integral, ti2(x) = integral(arctan(t)/t, t=0..x)}

function sfc_dilog(x: extended): extended;
  {-Return dilog(x) = Re(Li_2(x)), Li_2(x) = -integral(ln(1-t)/t, t=0..x)}

function sfc_trilog(x: extended): extended;
  {-Return the trilogarithm function trilog(x) = Re(Li_3(x))}

function sfc_polylog(n: integer; x: extended): extended;
  {-Return the polylogarithm Li_n(x) of integer order; real part for n>0,x>1}

function sfc_polylogr(s, x: extended): extended;
  {-Return the polylogarithm Li_s(x) of real order s >= -1, x <= 256;}
  { s > 0 if x > 256; real part if x > 1}

function sfc_pz(x: extended): extended;
  {-Return the prime zeta function P(x) = sum(1/p^x, p prime), x > 1/5; }
  { for x<1 the real part of P(x) is returned.}

function sfc_eta(s: extended): extended;
  {-Return the Dirichlet eta function}

function sfc_etaint(n: integer): extended;
  {-Return the Dirichlet function eta(n) for integer arguments}

function sfc_etam1(s: extended): extended;
  {-Return Dirichlet eta(s)-1}

function sfc_lerch(z,s,a: extended): extended;
  {-Return the Lerch transcendent Phi(z,s,a) = sum(z^n/(n+a)^s, n=0..INF), z <= 1, s >= -1, a >= 0}

function sfc_zeta(s: extended): extended;
  {-Return the Riemann zeta function at s, s<>1}

function sfc_zeta1p(x: extended): extended;
  {-Return the Riemann zeta function at 1+x, x<>0}

function sfc_zetam1(s: extended): extended;
  {-Return Riemann zeta(s)-1, s<>1}

function sfc_zetah(s,a: extended): extended;
  {-Return the Hurwitz zeta function zetah(s,a) = sum(1/(i+a)^s, i=0..INF), s<>1, a>0}

function sfc_zetaint(n: integer): extended;
  {-Return zeta(n) for integer arguments, n<>1}

function sfc_dbeta(s: extended): extended;
  {-Return the Dirichlet beta function sum((-1)^n/(2n+1)^s, n=0..INF)}

function sfc_dlambda(s: extended): extended;
  {-Return the Dirichlet lambda function sum(1/(2n+1)^s, n=0..INF), s<>1}

function sfc_lchi(s, x: extended): extended;
  {-Return Legendre's Chi-function chi(s,x); s>=0; real part if |x| > 1}

function sfc_bernpoly(n: integer; x: extended): extended;
  {-Return the Bernoulli polynomial B_n(x), 0 <= n <= MaxBernoulli}

function sfc_beint(s, x: extended): extended;
  {-Return the Bose-Einstein integral of real order s >= -1, real part if x > 1}

function sfc_fermi_dirac(n: integer; x: extended): extended;
  {-Return the integer order Fermi-Dirac integral F_n(x) = 1/n!*integral(t^n/(exp(t-x)+1), t=0..INF)}

function sfc_fdr(s,x: extended): extended;
  {-Return the Fermi-Dirac integral of real order s >= -1}

function sfc_harmonic(x: extended): extended;
  {-Return the harmonic number function H(x) = psi(x+1) + EulerGamma}

function sfc_harmonic2(x,r: extended): extended;
  {-Return the generalized harmonic function H(x,r) = zeta(r)-zetah(r,x+1); x >= -1}

function sfc_llci(x: extended): extended;
  {-Return the Lobachevski function L(x) = integral(-ln(|cos(t)|), t=0..x)}

function sfc_llsi(x: extended): extended;
  {-Return the Lobachevski function Lambda(x) = integral(-ln(|2sin(t)|), t=0..x)}


{#Z+}
{---------------------------------------------------------------------------}
{Obsolete functions, use sfc_fdr with s = -1/2, 1/2, 3/2, 5/2 }
{#Z-}

function sfc_fdm05(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}

function sfc_fdp05(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}

function sfc_fdp15(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}

function sfc_fdp25(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}


implementation

uses
  AMath,
  sfBasic, sfGamma;


{---------------------------------------------------------------------------}
function sfc_cl2(x: extended): extended;
  {-Return the Clausen function: integral(-ln(2*|sin(t/2)|),t=0..x) = Im(Li_2(exp(ix)))}
var
  t,z: extended;
  n: integer;
const
  {h:= x-> 1+sum((-1)^(k-1)*bernoulli(2*k)/(2*k+1)!/(2*k)*x^(2*k),k=1..201);}
  {chebyshev(h(x),x=-Pi/2..Pi/2,0.1e-20); (Note: C1H[0] is doubled)}
  C1H: array[0..10] of THexExtW = (
         ($4B40,$522A,$CA82,$8236,$4000),  {2.0345941803613285109}
         ($92C1,$0D67,$7C45,$8E25,$3FF9),  {1.7351858820274076806E-2}
         ($50C0,$287F,$9C2D,$E75E,$3FF0),  {5.5162804260905214154E-5}
         ($DA31,$F30B,$7A72,$D593,$3FE9),  {3.9781646276597633104E-7}
         ($760C,$4541,$5693,$FD96,$3FE2),  {3.6901802891787812104E-9}
         ($B375,$9BE0,$8612,$AAA9,$3FDC),  {3.8804092136367920624E-11}
         ($280D,$9A61,$20B1,$F817,$3FD5),  {4.4069697689724379977E-13}
         ($26AB,$BC64,$5496,$BE1D,$3FCF),  {5.2767393780251196689E-15}
         ($05E0,$4305,$02DA,$9775,$3FC9),  {6.5684035804606950287E-17}
         ($EC0B,$5F67,$98E4,$F8A0,$3FC2),  {8.4238217037971627742E-19}
         ($7F28,$86A9,$64C0,$D0F4,$3FBC)); {1.1061967718533347427E-20}
const
  {g:=x->-I*(polylog(2,exp(I*(x+Pi)))-polylog(2,exp(-I*(x+Pi))))/2;}
  {chebyshev(g(x)/x,x=-Pi/2..Pi/2, 0.1e-20);(Note: C2H[0] is doubled)}
  C2H: array[0..15] of THexExtW = (
         ($CF66,$9C5A,$DDE7,$A39B,$BFFF),  {-1.2781941777145306826}
         ($1E33,$1215,$4E31,$E133,$3FFA),  { 5.4980569301851715639E-2}
         ($B1E0,$A525,$327B,$FBFD,$3FF4),  { 9.6126194595060642936E-4}
         ($AA4E,$B610,$7538,$8672,$3FF0),  { 3.2054686822550476557E-5}
         ($3403,$0F98,$F4BE,$B26F,$3FEB),  { 1.3294616954255450141E-6}
         ($2652,$4EBE,$5190,$8558,$3FE7),  { 6.2093601824397519462E-8}
         ($B572,$95C3,$8459,$D710,$3FE2),  { 3.1296006563911126723E-9}
         ($45E0,$43FC,$E98E,$B6E7,$3FDE),  { 1.6635195381926697759E-10}
         ($0F2D,$F984,$7A08,$A1C9,$3FDA),  { 9.1965272507194254495E-12}
         ($C5F5,$A603,$7382,$937E,$3FD6),  { 5.2400377387584500936E-13}
         ($47E7,$1D77,$C855,$89B8,$3FD2),  { 3.0580384187365945412E-14}
         ($640A,$EC19,$63AE,$831F,$3FCE),  { 1.8196918249487950988E-15}
         ($6A2A,$2898,$0C67,$FDBC,$3FC9),  { 1.1003982631962615223E-16}
         ($A44C,$F379,$66A5,$F8DA,$3FC5),  { 6.7451775715424687730E-18}
         ($5366,$B641,$745C,$F6E8,$3FC1),  { 4.1827846515724770346E-19}
         ($4D9A,$9F08,$9D0C,$F770,$3FBD)); { 2.6198718087610612748E-20}
        {($FA33,$3EBE,$22C1,$FA25,$3FB9)}  { 1.6553211620337321920E-21}
var
  C1: array[0..10] of extended absolute C1H;
  C2: array[0..15] of extended absolute C2H;
begin
  {See MISCFUN [22], function clausn for formulas and hints. As pointed}
  {out, only absolute accuracy can be guaranteed close to the zeros.   }
  {Observed relative errors are e.g. 1700 eps_x for abs(x-Pi) ~ 6.7E-4,}
  {but even for abs(x-Pi) ~ 0.03 they are still about 128 eps_x.}

  {Therefore two separate Chebyshev expansions are calculated, one for }
  {z in (-Pi/2, Pi/2) the other for z in (Pi/2, 3*Pi/2), where z is the}
  {reduced argument. For z=0 or z=Pi cl2 is zero. Calculations are done}
  {using Maple with Digits:=50; and t_rcalc to convert to Hex/Extended.}
  {Note that both approximations are done for even functions, and there}
  {are only the even Chebyshev polynomials, so the argument for CSEvalX}
  {is 2(x/(Pi/2))^2 - 1 with the calculated coefficients.}

  {Argument reduction x mod Pi, |z| <= Pi/2}
  n := rem_pio2(0.5*x,z);
  z := 2.0*z;
  t := abs(z);

  if t=0.0 then begin
    {if x is an exact multiple of Pi then cl2(x)=0}
    sfc_cl2 := 0.0
  end
  else begin
    if odd(n) then begin
      {Use approximation for (Pi/2, 3*Pi/2) centered at Pi}
      if t<=1e-10 then begin
        {cl2(x+Pi) = x*(-ln(2) + 1/24*x^2 + 1/960*x^4+ O(x^6))}
        sfc_cl2 := -ln2*z;
      end
      else begin
        {Use Chebyshev expansion for cl2(x+Pi)/x calculated from}
        {cl2(x) = Im(Li_2(exp(ix) = -i*(Li_2(exp(ix))-Li_2(exp(-ix)))/2}
        sfc_cl2 := CSEvalX(2.0*(sqr(2.0*z/Pi)-0.5),C2,16) * z;
      end;
    end
    else begin
      {Use approximation for (-Pi/2, Pi/2)}
      if t<=1e-9 then begin
        {cl2(x) = x*(1-ln(x) + 1/72*x^3 + 1/14400*x^5 + O(x^6) }
        sfc_cl2 := z*(1.0-ln(t));
      end
      else begin
        {Use Chebyshev expansion calculated from [22], (3) or HMF[1], 27.8.2}
        {cl2(x) = -x*ln|x| + x - sum((-1)^k*B_2k/(2k+1)!/(2k)*x^(2k+1),k=1..Inf)}
        sfc_cl2 := (CSEvalX(2.0*(sqr(2.0*z/Pi)-0.5),C1,11)-ln(t))*z;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_dilog(x: extended): extended;
  {-Return dilog(x) = Re(Li_2(x)), Li_2(x) = -integral(ln(1-t)/t, t=0..x)}
var
  t: extended;
const
  xbig = 1.8446744073709551616e19;
const
  nsp = 24;
  csp: array[0..nsp-1] of extended = (
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
         0.1635389806752377100051821734570e-17,
         0.2551852874940463932310901642581e-18,
         0.3999020621999360112770470379519e-19,
         0.6291501645216811876514149171199e-20);
const
  Pi26 = 1.625+0.01993406684822643647;  {Pi^2/6}
begin
  {Ref: W. Fullerton [14] and [20], file dspenc.f}

  {Note that there is some confusion about the naming: some authors}
  {and/or computer algebra systems use dilog(x)= Li_2(1-x) and then}
  {call Li_2(x) Spence function/integral or similar.}

  {The imaginary part Im(Li_2(x) is 0 for x<=1 and -Pi*ln(x) for x>1}

  if x>2.0 then begin
    t := 2.0*Pi26 - 0.5*sqr(ln(x));
    if x>=xbig then sfc_dilog := t
    else sfc_dilog := t - (1.0 + CSEvalX(4.0/x-1.0, csp, nsp))/x;
  end
  else if x>1.0 then begin
    {1 < x <= 2}
    sfc_dilog := Pi26 - 0.50*ln(x)*ln(sqr((x-1.0))/x)
                 + (x-1.0)*(1.0 + CSEvalX(4.0*(x-1.0)/x-1.0, csp, nsp))/x;
  end
  else if x>0.5 then begin
    {0.5 < x <= 1}
    if x=1.0 then sfc_dilog := Pi26
    else sfc_dilog := Pi26 - ln(x)*ln1p(-x) - (1.0-x)*(1.0+CSEvalX(4.0*(1.0-x)-1.0, csp, nsp));
  end
  else if x>=0.0 then begin
    {0 <= x <= 0.5}
    sfc_dilog := x*(1.0 + CSEvalX(4.0*x-1.0, csp, nsp));
  end
  else if x>-1.0 then begin
    {-1 < x < 0}
    sfc_dilog := -0.5*sqr(ln(1.0-x)) - x*(1.0+CSEvalX(4.0*x/(x-1.0)-1.0, csp, nsp))/(x-1.0);
  end
  else begin
    {x <= -1.0}
    t := ln1p(-x);
    t := -Pi26 - 0.5*t*(2.0*ln(-x)-t);
    if x <= -xbig then sfc_dilog := t
    else begin
      sfc_dilog := t + (1.0 + CSEvalx(4.0/(1.0-x)-1.0, csp, nsp))/(1.0-x);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_trilog(x: extended): extended;
  {-Return the trilogarithm function trilog(x) = Re(Li_3(x))}
var
  a,b: extended;
begin
  if IsNan(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_trilog := Nan_x;
    exit;
  end;
  if abs(x)>1.0 then begin
    b := 1.0/x;
    if abs(b) > sqrt_epsh then a := sfc_polylog(3,b)
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
      sfc_trilog := a + b;
    end
    else begin
      {Could be handled by sfc_polylog, but this branch}
      {avoids one recursive call with sfc_polylog(1/x):}
      {Re(RHS) becomes -Pi^2/6 * ln(-x) - ln^3(-x)/6 }
      b := ln(-x);
      b := (PiSqr + sqr(b))*b/SIXX;
      sfc_trilog := a - b;
    end;
  end
  else begin
    {-1 <= x <= 1}
    sfc_trilog := sfc_polylog(3,x);
  end;
end;


{---------------------------------------------------------------------------}
function polylog_series(n: integer; x: extended): extended;
  {-polylog via series, |n|>2, |x| < 1}
var
  e,k,p,s,t,w,y,z: extended;
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
  until abs(t) <= eps_x*abs(s);
  polylog_series := s+e;
end;


{---------------------------------------------------------------------------}
function cr15_complex(n: integer; x: extended): extended;
  {-Compute polylog(n,x) using Crandall[31] (1.5); n<-1, |log(x)| < 2*Pi}
type
  TComplex = record r,i: extended; end;
var
  lx,s,p,t: TComplex;
  bk: extended;
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
  s.r  := sfc_bernoulli(1-n)/(1-n);
  s.i  := 0.0;
  repeat
    inc(k);
    t.r := p.r*lx.r - p.i*lx.i;
    t.i := p.r*lx.i + p.i*lx.r;
    p.r := t.r/k;
    p.i := t.i/k;
    bk := sfc_bernoulli(k-n+1)/(k-n+1);
    if bk<>0.0 then begin
      t.r := p.r*bk;
      t.i := p.i*bk;
      s.r := s.r + t.r;
      s.i := s.i + t.i;
      done := (abs(t.r)+abs(t.i)) <= eps_x*(abs(s.r)+abs(s.i));
    end;
  until done;
  {rt.r = Re(power(-lx, n-1))}
  p.r := hypot(lx.r,lx.i);
  p.i := arctan2(-lx.i,-lx.r);
  t.r := power(p.r,n-1)*cos((n-1)*p.i);
  cr15_complex := sfc_fac(-n)*t.r - s.r;
end;


{---------------------------------------------------------------------------}
function polylog_neg(n: integer; x: extended): extended;
  {-Polylog for n<0 and x<1}
var
  s,p,lx,t: extended;
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
  if x = 1.0 then polylog_neg := PosInf_x
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
  else if x=-1.0 then polylog_neg := -sfc_etaint(n)
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
    s  := sfc_bernoulli(1-n)/(1-n);
    done := false;
    repeat
      inc(k);
      p := p*lx/k;
      t := sfc_bernoulli(k-n+1);
      if t<>0.0 then begin
        t := t/(k-n+1)*p;
        s := s + t;
        done := abs(t) <= abs(s)*eps_x;
      end;
    until done;
    polylog_neg := sfc_fac(-n)*power(-lx, n-1) - s;
  end
  else begin
    {Same formula as above but with complex arithmetic}
    polylog_neg := cr15_complex(n,x);
  end;
end;


{---------------------------------------------------------------------------}
function polylog_inv(n: integer; x: extended): extended;
  {-Return Li_(n,x) for x < -1 with inversion formula}
var
  s,t,z: extended;
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
    t := -sfc_etaint(j);
    j := n-j;
    if j>0 then t := sfc_taylor(z,j)*t;
    s := s + t;
  end;
  t := sfc_polylog(n, 1.0/x);
  if odd(n) then t := -t;
  polylog_inv := 2.0*s - t - sfc_taylor(z,n);
end;


{---------------------------------------------------------------------------}
function polylog_xg05(n: integer; x: extended): extended;
  {-Return Re(polylog(n,x)), n>2, x >= 0.5}
var
  s,t,z,f,h: extended;
  k,i: integer;
const
  xmax = 256.0;
  hinv = -123456789.0; {flag h not yet valid}
begin
  if x > xmax then begin
    {Crandall [31], Algorithm 3.1 for positive x > 1, with r2=xmax}
    {Recursive Li_n(x) = (Li_n(sqrt(x)) + Li_n(-sqrt(x)))/2^(n-1) }
    z := sqrt(x);
    h := sfc_polylog(n, -z);
    f := polylog_xg05(n,z);
    polylog_xg05 := ldexp(f+h, n-1);
    exit;
  end;
  {Crandall [31], formula 1.4; valid for |x| < exp(2*pi) }
  z := ln(x);
  f := 1;
  s := sfc_zetaint(n);
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
      t := sfc_zetaint(i)*f;
      s := s + t;
      if abs(t)<=eps_x*abs(s) then break;
    end;
  until false;
  f := sfc_harmonic(n-1) - ln(abs(z));
  if h=hinv then begin
    {h = z^(n-1)/(n-1) not yet calculated}
    h := sfc_taylor(z, n-1);
  end;
  polylog_xg05 := s + h*f;
end;


{---------------------------------------------------------------------------}
function sfc_polylog(n: integer; x: extended): extended;
  {-Return the polylogarithm Li_n(x) of integer order; real part for n>0,x>1}
var
  s,t: extended;
begin
  sfc_polylog := Nan_x;
  {Ref: Cephes [7] function polylog in file polylog.c and}
  {Crandall [31], Note on fast polylogarithm computation.}

  if IsNanOrInf(x) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;

  if x=0.0 then begin
    sfc_polylog := 0.0;
    exit;
  end
  else if abs(x)=1.0 then begin
    if x>0.0 then sfc_polylog := sfc_zetaint(n)
    else sfc_polylog := -sfc_etaint(n);
    exit;
  end
  else if n<0 then begin
    sfc_polylog := polylog_neg(n,x);
    exit;
  end
  else if n<=2 then begin
    case n of
        0: sfc_polylog := x/(1.0-x);
        1: begin
             if x<1.0 then sfc_polylog := -ln1p(-x)
             else sfc_polylog := -ln(x-1.0);
           end;
      else sfc_polylog := sfc_dilog(x);
    end;
    exit;
  end
  else if (n>40) and (abs(x) < 1.0) then begin
    if n > 63 then sfc_polylog := x
    else sfc_polylog := x*(1.0 + x*ldexp(1,-n));
    exit;
  end;

  {Here n>2}
  if x < -1.0 then begin
    {Inversion formula for x < -1}
    sfc_polylog := polylog_inv(n,x);
  end
  else if (x < 0.5) and (x > -0.875) then begin
    {Use standard series sum(x^k/k^n, k=1..)}
    sfc_polylog := polylog_series(n,x);
  end
  else if x<0 then begin
    {For -1 < x <= -0.875 use recursive calls: Crandall[31], formula (1.2)}
    s := sfc_polylog(n,x*x);
    t := sfc_polylog(n,-x);
    sfc_polylog := ldexp(s,1-n) - t;
  end
  else begin
    {Here x >= 0.5, use Crandall[31], formula (1.4)}
    sfc_polylog := polylog_xg05(n,x);
  end;
end;


{---------------------------------------------------------------------------}
function liae_neg(s,x: extended): extended;
  {-Asymptotic formula for Li(s,x), x << -1}
var
  w,sum,r,t,f: extended;
  k: integer;
begin
  {Ref: Wood [73], (11.1)}
  k := 0;
  w := ln(abs(x));
  f := power(w, s)*sfc_rgamma(s+1.0);
  t := MaxExtended;
  sum := 0.5*f;
  w  := w*w;
  repeat
    f := f/w*(s-k)*(s-k-1);
    k := k + 2;
    r := t;
    t := sfc_etaint(k)*f;
    sum := sum + t;
    t := abs(t);
  until (t >= r) or (t <= eps_x*abs(sum));
  liae_neg := -2.0*sum;
end;


{---------------------------------------------------------------------------}
type
  tli_rec  = record
               s,x,g: extended;  {g=1/Gamma(s)}
             end;
  pli_rec = ^tli_rec;


{---------------------------------------------------------------------------}
function lif(t: extended; param: pointer): extended; {$ifdef BIT16} far;{$endif}
  {-Polylog integrand t^(s-1)/(exp(t)/x-1)/Gamma(s) }
var
  z,p: extended;
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
function liae_pos(s,lnx: extended): extended;
  {-Asymptotic formula for Li(s,ln(x)), x >> 1}
var
  w,sum,r,t,f: extended;
  k: integer;
const
  KMAX = 500;
begin
  {Ref: Wood [73], (11.2)}
  k := 0;
  w := lnx;
  f := power(w, s)*sfc_rgamma(s+1.0);
  t := MaxExtended;
  sum := -0.5*f;
  w  := w*w;
  repeat
    f := f/w*(s-k)*(s-k-1);
    k := k + 2;
    r := t;
    t := sfc_zetaint(k)*f;
    sum := sum + t;
    t := abs(t);
  until (k>KMAX) or (t >= r) or (t <= eps_x*abs(sum));
  {Test no convergence}
  if (k>KMAX) and (RTE_NoConvergence>0) then RunError(byte(RTE_NoConvergence));
  liae_pos := 2.0*sum;
end;


{---------------------------------------------------------------------------}
function li_xgt0(s,x: extended): extended;
  {-Return Re(polylog(s,x)), s > -1, 1 < x, frac(s)<>0}
var
  sum,t,z,f,h: extended;
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
      h := sfc_polylogr(s, -z);
      f := li_xgt0(s,z);
      li_xgt0 := 0.5*(f+h)*exp2(s);
    end;
    exit;
  end;
  if x > xmax then begin
    li_xgt0 := Nan_x;
    exit;
  end;
  {NIST [30], 25.12.12; valid for |x| < exp(2*pi) }
  z := ln(x);
  f := z;
  sum := sfc_zeta(s-1)*f;
  for k:=2 to maxint do begin
    f := f*z/k;
    t := sfc_zeta(s-k)*f;
    sum := sum + t;
    if abs(t)<=eps_x*abs(sum) then break;
  end;
  sum := sum + sfc_zeta(s);
  f := sfc_gamma(1.0-s);
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
function sfc_polylogr(s,x: extended): extended;
  {-Return the polylogarithm Li_s(x) of real order s >= -1, x <= 256;}
  { s > 0 if x > 256; real part if x > 1}
var
  z, aerr: extended;
  ierr: integer;
  para: tli_rec;
  neval: longint;
const
  eps = 4e-19;
begin
  sfc_polylogr := NaN_x;
  if IsNanOrInf(s) or IsNanOrInf(x) or (s < -1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    exit;
  end;
  z := abs(x);
  if z=0.0 then sfc_polylogr := 0.0
  else if z=1.0 then begin
    if x=1.0 then sfc_polylogr := sfc_zeta(s)
    else sfc_polylogr := -sfc_eta(s)
  end
  else if frac(s)=0.0 then begin
    {use integer version}
    sfc_polylogr := sfc_polylog(round(s),x)
  end
  else if z<1.0 then begin
    if s>=40.5 then begin
      {from series and z<1}
      if s>=64.0 then sfc_polylogr := x
      else sfc_polylogr := x*(1.0 + x*exp2(-s))
    end
    else begin
      {Li_s(x) = x*Phi(x,s,1), see e.g. NIST[30], 25.14.3}
      {or http://functions.wolfram.com/10.08.26.0010.01}
      sfc_polylogr := x*sfc_lerch(x,s,1.0);
    end;
  end
  else if x > 1.0 then begin
    {will return NaN if x > 256}
    sfc_polylogr := li_xgt0(s,x);
  end
  else if s>0.0 then begin
    if x <= -1e20 then begin
      {Use asymptotic formula Wood [73], (11.1)}
      sfc_polylogr := liae_neg(s,x);
    end
    else if s > log2(z) + 65.0 then begin
      sfc_polylogr := x;
    end
    else begin
      {Use integral formula, see Wood [73], (1.1) or }
      {http://functions.wolfram.com/10.08.07.0001.01}
      para.s := s;
      para.x := x;
      para.g := sfc_rgamma(s);
      sfc_intdei_p({$ifdef FPC_ProcVar}@{$endif}lif, @para, 0.0, eps, z, aerr, neval, ierr);
      if ierr=0 then sfc_polylogr := z;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ti(s,x: extended): extended;
  {-Return the inverse tangent integral of order s >= 0}
var
  t,z: extended;
begin
  if IsNanOrInf(x) or IsNanOrInf(s) or (s < 0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_ti := NaN_x;
    exit;
  end;
  t := abs(x);
  if t=0.0 then sfc_ti := 0.0
  else if s=0.0 then begin
    {ti(0,x) = x/(1+x*2)}
    if t > 1e10 then sfc_ti := 1.0/x
    else sfc_ti := x/(1.0 + x*x);
  end
  else if s=1.0 then begin
    {ti(1,x) = arctan(x)}
    sfc_ti := arctan(x);
  end
  else if s=2.0 then begin
    sfc_ti := sfc_ti2(x);
  end
  else begin
    z := 1.85*ln(t)+40.5;
    if s > z then begin
      {shortcut for x^2/3^s < eps_x/2}
      sfc_ti := x;
    end
    else begin
      t := sfc_lerch(-x*x,s,0.5);
      if IsNanOrInf(x) then sfc_ti := t
      else sfc_ti := exp2(-s)*x*t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_ti2(x: extended): extended;
  {-Return the inverse tangent integral, ti2(x) = integral(arctan(t)/t, t=0..x)}
var
  y: extended;
const
  {chebyshev(int(arctan(t)/t, t=0..x)/x, x=-1..1, 0.5e-20); aics[0] doubled}
  nai = 23;
  aics: array[0..nai-1] of extended = (
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
          +0.794150680705224718491128609488e-18,
          -0.124060349865156550047951820204e-18,
          +0.194621210567853358154874940741e-19{,
          -0.306490766819307325157212222222e-20});
begin
  {ti2(x) = integral(arctan(t)/t,  t=0..x)}
  {ti2(x) = i/2*[dilog(-i*x) - dilog(i*x)]}
  {ti2(x) = ti2(1/x) + Pi/2*ln(x),  x<>0}

  {See references [22] (function atnint) and [21] (specfunc/atanint.c)}
  {Note that the precision of both refs are not suitable for extended!}

  y := abs(x);
  if y<3.41e-5 then begin
    {ti2(x) := x*(1 - 1/9x^2 + 1/25x^4 +O(x^6)}
    if y<6.98e-10 then sfc_ti2 := x
    else sfc_ti2 := x*(1.0-sqr(x)/9.0);
  end
  else begin
    {ti2(x)/x is even and it's Chebyshev approximation on [-1,1] contains }
    {only the even Chebyshev polynomials. Since T(2n,x) = T(n,2x^2-1), the}
    {coefficients aics can be used with CSEvalX and the argument 2x^2-1.  }
    if y<=1.0 then sfc_ti2 := x*CSEvalX(2.0*sqr(x)-1.0, aics, nai)
    else begin
      {ti2(x) = Pi/2*ln(x) + 1/x + O(1/x^2), x->Inf}
      if y>=3.04e9 then y := Pi_2*ln(y) + 1.0/y
      else begin
        {Use ti2(x) = ti2(1/x) + sign(x)*Pi/2*ln(|x|) for |x|>1}
        y := Pi_2*ln(y) + CSEvalX(2.0/sqr(x)-1.0, aics, nai)/y;
      end;
      if x>0.0 then sfc_ti2 := y else sfc_ti2 := -y;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_pz(x: extended): extended;
  {-Return the prime zeta function P(x) = sum(1/p^x, p prime), x > 1/5; }
  { for x<1 the real part of P(x) is returned.}
var
  s,t,y: extended;
  k,m: integer;
begin
  {Im(P(x)) = Pi   for 1/2 < x < 1
            = Pi/2 for 1/3 < x < 1/2
            = Pi/6 for 1/5 < x < 1/3
  }
  if IsNan(x) or (x < 0.2) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_pz := NaN_x;
    exit;
  end;
  if x>=27.0 then begin
    {use only primes up to 7}
    if x>=16450.0 then sfc_pz := 0.0
    else if x>120.0 then sfc_pz := exp2(-x)
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
      sfc_pz := s+t;
    end;
  end
  else if x=1.0 then sfc_pz := PosInf_x
  else begin
    {P(x) = sum(mu(k)/k*ln(zeta(kx)), k>0) by Moebius inversion, see e.g.}
    {H. Cohen, High Precision Computation of Hardy-Littlewood Constants, }
    {Section 2.1, from  http://www.math.u-bordeaux.fr/~cohen/hardylw.dvi }
    {Cohen's speed-up is used with A=7.}
    y := -x;
    s := exp2(y) + exp3(y) + exp5(y) + exp7(y);
    {maximum observed k is 91 for x near 1/5}
    for k:=1 to n_moeb do begin
      m := moebius[k];
      if m<>0 then begin
        y := k*x;
        if y=1.0 then begin
          sfc_pz := NegInf_x;
          exit;
        end;
        {Compute the next term ln(zeta(y)*prod(1-p^y))*mu(k)/k}
        if y<1.0 then t := ln(abs(sfc_zeta(y)))
        else t := ln1p(sfc_zetam1(y));
        t := t + ln1p(-exp2(-y));
        t := t + ln1p(-exp3(-y));
        t := t + ln1p(-exp5(-y));
        t := t + ln1p(-exp7(-y));
        s := s + m*t/k;
        if abs(t) < eps_x*abs(s) then begin
          sfc_pz := s;
          exit;
        end;
      end;
    end;
    {No convergence}
    sfc_pz := s;
    if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
  end;
end;


{---------------------------------------------------------------------------}
const
  ETAEPS = 1e-9;

{---------------------------------------------------------------------------}
function etam1pos(s: extended): extended;
  {-Return the eta(s)-1 = -sum((-1)^k/k^s, k=2..) for s >= -ETAEPS; internal use}

const
  dh: array[2..24] of THexExtW = (
        ($FFA5,$FFFF,$FFFF,$FFFF,$BFFE),  {-0.99999999999999999506430323408864009302707409760650   }
        ($875D,$FFFE,$FFFF,$FFFF,$3FFE),  {+0.99999999999999477309712489986985851567146936528101   }
        ($D5A9,$FEFB,$FFFF,$FFFF,$BFFE),  {-0.99999999999907570687373807736430090108503647599525   }
        ($3E9B,$B83F,$FFFF,$FFFF,$3FFE),  {+0.99999999993474107123666050197526788003473422599175   }
        ($787A,$74E2,$FFF5,$FFFF,$BFFE),  {-0.99999999754516889043092198752546995530922208300482   }
        ($9FF0,$7505,$FF0B,$FFFF,$3FFE),  {+0.99999994306292316806008385807007727156754522290285   }
        ($2930,$A5AA,$F115,$FFFF,$BFFE),  {-0.99999911097044304457819242638771628351284408679994   }
        ($6065,$D3A9,$59D2,$FFFF,$3FFE),  {+0.99999009511126280553264350727993678656849990880578   }
        ($1A12,$43A3,$9FBC,$FFFA,$BFFE),  {-0.99991796823782089316825215441770081101374648485252   }
        ($0A33,$3791,$E44C,$FFDD,$3FFE),  {+0.99947954998748770036508902917665860666132371180327   }
        ($E514,$245E,$659A,$FF56,$BFFE),  {-0.99741206181749538061964608067153326403095105573943   }
        ($5B36,$DD12,$214D,$FD5F,$3FFE),  {+0.98973282004323819299371512908106770568956690464520   }
        ($544C,$0A46,$35FD,$F78F,$BFFE),  {-0.96702897479760824696922362003099562015851811010573   }
        ($E565,$E3A5,$FCA7,$E9C9,$3FFE),  {+0.91323832606180806715735142935851714059264865842761   }
        ($983E,$D706,$464E,$CF8F,$BFFE),  {-0.81077994751742677227759487569665336999099255999308   }
        ($6535,$017E,$D18E,$A766,$3FFE),  {+0.65391263691844299666858829009021366603535356790712   }
        ($6615,$FFBF,$F03A,$EC59,$BFFD),  {-0.46162367553904352979303183031457789989618319051142   }
        ($ED27,$D227,$22A5,$8C88,$3FFD),  {+0.27447613023930715026708917962385218097998528309956   }
        ($838F,$6B07,$85EE,$870D,$BFFC),  {-0.13188752429665086110446620766901353799621544888100   }
        ($9514,$E77C,$BD13,$C796,$3FFA),  {+0.48727739891972228733661913925793276796577280986824e-1}
        ($5449,$F08D,$C274,$D36F,$BFF8),  {-0.12905063533033740943161602774867625818271608663177e-1}
        ($A6BB,$3D4D,$0397,$8E43,$3FF6),  {+0.21707423941183752637782342766808453857479577229903e-2}
        ($1D19,$2582,$5684,$B618,$BFF2)); {-0.17365939152947002110225874213446763085983661783922e-3}
const
  ln2m1h: THexExtW = ($0CA8,$5C61,$D010,$9D1B,$BFFD);  {-3.0685281944005469057E-1}
  lnpi2h: THexExtW = ($AE23,$5098,$D92D,$E735,$3FFC);  {0.2257913526447274323630976} {ln(pi/2)/2}
var
  d: array[2..24] of extended absolute dh;
  p: array[2..24] of extended;
  x,sum: extended;
  k: integer;

begin

  if s=1.0 then begin
    etam1pos := extended(ln2m1h); {ln(2)-1}
    exit;
  end
  else if abs(s) <= ETAEPS then begin
    etam1pos := extended(lnpi2h)*s - 0.5;
    exit;
  end;

  x := -s;

  {Calculate p[k] := 1/k^s but only if necessary. Prime powers are}
  {evaluated with power/exp?, otherwise products of p[k] are used.}
  p[2] := exp2(x);
  if s >= 120.0 then begin
    etam1pos := -p[2];
    exit;
  end;

  p[3] := exp3(x);
  p[4] := p[2]*p[2];
  if s >= 50.0 then begin
    etam1pos := (p[3] - p[4]) - p[2];
    exit;
  end;

  p[10] := exp10(x);
  p[ 5] := p[10]/p[2];
  p[ 6] := p[2]*p[3];
  p[ 7] := exp7(x);
  p[ 8] := p[2]*p[4];
  p[ 9] := p[3]*p[3];
  if s >= 28.0 then begin
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
  p[22] := p[2]*p[11];
  if s>=19.5 then begin
    sum := 0.0;
    for k:=22 downto 3 do begin
      if odd(k) then sum := sum + p[k]
      else sum := sum - p[k];
    end;
    etam1pos := sum - p[2];
  end
  else begin
    p[23] := power(23,x);
    p[24] := p[3]*p[8];
    {Convergence acceleration: see P. Borwein[36]}

    {The d[] are from P. Borwein's Algorithm 2, but scaled and shifted.}
    {Calculated with Maple VR4 and T_RCalc/xh using n:=23; Digits:=50; }

    { d(n,k) = n*sum((n+j-1)!*4^j/((n-j)!*(2*j)!),j=0..k); }
    { d[2]   := 1/d(n,n)-1.0;                              }
    { d[i]   := (-1)^(i-2)*(d(n,i-2)/d(n,n)-1), i=3..n+1;  }
    sum := p[24]*d[24];
    for k:=23 downto 3 do sum := sum + d[k]*p[k];
    etam1pos := sum + d[2]*p[2];
  end;
end;


{---------------------------------------------------------------------------}
function zetap(s,sc: extended): extended;
  {-Return the Riemann zeta function at s>0, s<>1, sc=1-s}
var
  y: extended;
{Based on boost_1_42_0\boost\math\special_functions\zeta.hpp [19]}
{Copyright John Maddock 2007, see 3rdparty.ama for Boost license}
const
  P1: array[0..5] of extended = (
        0.243392944335937499969,
       -0.496837806864865688082,
        0.0680008039723709987107,
       -0.00511620413006619942112,
        0.000455369899250053003335,
       -0.279496685273033761927e-4);
  Q1: array[0..6] of extended = (
        1.0,
       -0.30425480068225790522,
        0.050052748580371598736,
       -0.00519355671064700627862,
        0.000360623385771198350257,
       -0.159600883054550987633e-4,
        0.339770279812410586032e-6);
  P2: array[0..5] of extended = (
        0.577215664901532860605,
        0.222537368917162139445,
        0.0356286324033215682729,
        0.00304465292366350081446,
        0.000178102511649069421904,
        0.700867470265983665042e-5);
  Q2: array[0..6] of extended = (
        1.0,
        0.259385759149531030085,
        0.0373974962106091316854,
        0.00332735159183332820617,
        0.000188690420706998606469,
        0.635994377921861930071e-5,
        0.226583954978371199405e-7);
  P4: array[0..6] of extended = (
       -0.053725830002359501027,
        0.0470551187571475844778,
        0.0101339410415759517471,
        0.00100240326666092854528,
        0.685027119098122814867e-4,
        0.390972820219765942117e-5,
        0.540319769113543934483e-7);
  Q4: array[0..7] of extended = (
        1.0,
        0.286577739726542730421,
        0.0447355811517733225843,
        0.00430125107610252363302,
        0.000284956969089786662045,
        0.116188101609848411329e-4,
        0.278090318191657278204e-6,
       -0.19683620233222028478e-8);
  P7: array[0..7] of extended = (
       -2.49710190602259407065,
       -3.36664913245960625334,
       -1.77180020623777595452,
       -0.464717885249654313933,
       -0.0643694921293579472583,
       -0.00464265386202805715487,
       -0.000165556579779704340166,
       -0.252884970740994069582e-5);
  Q7: array[0..8] of extended = (
        1.0,
        1.01300131390690459085,
        0.387898115758643503827,
        0.0695071490045701135188,
        0.00586908595251442839291,
        0.000217752974064612188616,
        0.397626583349419011731e-5,
       -0.927884739284359700764e-8,
        0.119810501805618894381e-9);
  P15: array[0..8] of extended = (
        -4.78558028495135548083,
        -3.23873322238609358947,
        -0.892338582881021799922,
        -0.131326296217965913809,
        -0.0115651591773783712996,
        -0.000657728968362695775205,
        -0.252051328129449973047e-4,
        -0.626503445372641798925e-6,
        -0.815696314790853893484e-8);
  Q15: array[0..8] of extended = (
         1.0,
         0.525765665400123515036,
         0.10852641753657122787,
         0.0115669945375362045249,
         0.000732896513858274091966,
         0.30683952282420248448e-4,
         0.819649214609633126119e-6,
         0.117957556472335968146e-7,
        -0.193432300973017671137e-12);
  P42: array[0..8] of extended = (
        -10.3948950573308861781,
         -2.82646012777913950108,
         -0.342144362739570333665,
         -0.0249285145498722647472,
         -0.00122493108848097114118,
         -0.423055371192592850196e-4,
         -0.1025215577185967488e-5,
         -0.165096762663509467061e-7,
         -0.145392555873022044329e-9);
  Q42: array[0..9] of extended = (
         1.0,
         0.205135978585281988052,
         0.0192359357875879453602,
         0.00111496452029715514119,
         0.434928449016693986857e-4,
         0.116911068726610725891e-5,
         0.206704342290235237475e-7,
         0.209772836100827647474e-9,
        -0.939798249922234703384e-16,
         0.264584017421245080294e-18);
const
   x1: extended = 81487.0/65536.0;
   x4: extended = 366299.0/524288.0;
begin
  {$ifdef debug}
    if abs((s+sc)-1.0)>5e-19*maxx(1.0,abs(s)) then begin
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
      zetap := NaN_x;
      exit;
    end;
  {$endif}
  if s<1.0 then begin
    y := PolEvalX(sc,P1,6)/PolEvalX(sc,Q1,7) - x1;
    zetap := (y+sc)/sc;
  end
  else if s<=2.0 then begin
    s := -sc;
    y := PolEvalX(s,P2,6)/PolEvalX(s,Q2,7);
    zetap := y - 1.0/sc;
  end
  else if s<=4.0 then begin
    s := s-2.0;
    y := PolEvalX(s,P4,7)/PolEvalX(s,Q4,8) + x4;
    zetap := y - 1.0/sc;
  end
  else if s<=7.0 then begin
    s := s-4.0;
    y := PolEvalX(s,P7,8)/PolEvalX(s,Q7,9);
    zetap := 1.0 + exp(y);
  end
  else if s<15.0 then begin
    s := s-7.0;
    y := PolEvalX(s,P15,9)/PolEvalX(s,Q15,9);
    zetap := 1.0 + exp(y);
  end
  else if s<42.0 then begin
    s := s-15.0;
    y := PolEvalX(s,P42,9)/PolEvalX(s,Q42,10);
    zetap := 1.0 + exp(y);
  end
  else if s<64.0 then zetap := 1.0 + exp2(-s)
  else zetap := 1.0;
end;


{---------------------------------------------------------------------------}
function sfc_etaint(n: integer): extended;
  {-Return the Dirichlet function eta(n) for integer arguments}
const
  etahex: array[0..64] of THexExtW = (
            ($0000,$0000,$0000,$8000,$3FFE),  {+5.0000000000000000000E-1}
            ($79AC,$D1CF,$17F7,$B172,$3FFE),  {+6.9314718055994530943E-1} {ln(2)}
            ($9918,$983E,$3312,$D28D,$3FFE),  {+8.2246703342411321821E-1}
            ($8EA3,$4049,$803B,$E6CB,$3FFE),  {+9.0154267736969571407E-1}
            ($F2FE,$EDD3,$BE56,$F270,$3FFE),  {+9.4703282949724591755E-1}
            ($0DC3,$DD50,$D75D,$F8DC,$3FFE),  {+9.7211977044690930596E-1}
            ($26B0,$A8DB,$1389,$FC4D,$3FFE),  {+9.8555109129743510409E-1}
            ($BAB7,$C7BF,$A0EA,$FE1A,$3FFE),  {+9.9259381992283028266E-1}
            ($FBAF,$272C,$2042,$FF09,$3FFE),  {+9.9623300185264789924E-1}
            ($EBAC,$43E9,$1B9E,$FF83,$3FFE),  {+9.9809429754160533077E-1}
            ($BA16,$8BE8,$0D9C,$FFC1,$3FFE),  {+9.9903950759827156566E-1}
            ($0FDC,$1D92,$5B03,$FFE0,$3FFE),  {+9.9951714349806075415E-1}
            ($D69C,$9E16,$1EA1,$FFF0,$3FFE),  {+9.9975768514385819088E-1}
            ($B52F,$23EE,$0A49,$FFF8,$3FFE),  {+9.9987854276326511548E-1}
            ($52EC,$A017,$0372,$FFFC,$3FFE),  {+9.9993917034597971818E-1}
            ($3D85,$74C0,$0127,$FFFE,$3FFE),  {+9.9996955121309923808E-1}
            ($871B,$CD21,$0062,$FFFF,$3FFE),  {+9.9998476421490610646E-1}
            ($BF99,$0379,$8021,$FFFF,$3FFE),  {+9.9999237829204101199E-1}
            ($49E9,$064F,$C00B,$FFFF,$3FFE),  {+9.9999618786961011349E-1}
            ($05F0,$AE11,$E003,$FFFF,$3FFE),  {+9.9999809350817167510E-1}
            ($1CDC,$3A59,$F001,$FFFF,$3FFE),  {+9.9999904661158152213E-1}
            ($52AB,$68DD,$F800,$FFFF,$3FFE),  {+9.9999952325821554283E-1}
            ($B299,$22F9,$FC00,$FFFF,$3FFE),  {+9.9999976161323082254E-1}
            ($E242,$0BA9,$FE00,$FFFF,$3FFE),  {+9.9999988080131843951E-1}
            ($9FF6,$03E3,$FF00,$FFFF,$3FFE),  {+9.9999994039889239462E-1}
            ($F529,$014B,$FF80,$FFFF,$3FFE),  {+9.9999997019885696281E-1}
            ($AC5B,$006E,$FFC0,$FFFF,$3FFE),  {+9.9999998509923199657E-1}
            ($E572,$0024,$FFE0,$FFFF,$3FFE),  {+9.9999999254955048496E-1}
            ($4CD0,$000C,$FFF0,$FFFF,$3FFE),  {+9.9999999627475340009E-1}
            ($19B0,$0004,$FFF8,$FFFF,$3FFE),  {+9.9999999813736941811E-1}
            ($5DEB,$0001,$FFFC,$FFFF,$3FFE),  {+9.9999999906868228147E-1}
            ($74A5,$0000,$FFFE,$FFFF,$3FFE),  {+9.9999999953434033146E-1}
            ($26E2,$0000,$FFFF,$FFFF,$3FFE),  {+9.9999999976716989595E-1}
            ($0CF6,$8000,$FFFF,$FFFF,$3FFE),  {+9.9999999988358485804E-1}
            ($0452,$C000,$FFFF,$FFFF,$3FFE),  {+9.9999999994179239904E-1}
            ($0171,$E000,$FFFF,$FFFF,$3FFE),  {+9.9999999997089618955E-1}
            ($007B,$F000,$FFFF,$FFFF,$3FFE),  {+9.9999999998544809144E-1}
            ($0029,$F800,$FFFF,$FFFF,$3FFE),  {+9.9999999999272404461E-1}
            ($000E,$FC00,$FFFF,$FFFF,$3FFE),  {+9.9999999999636202195E-1}
            ($0005,$FE00,$FFFF,$FFFF,$3FFE),  {+9.9999999999818101087E-1}
            ($0002,$FF00,$FFFF,$FFFF,$3FFE),  {+9.9999999999909050541E-1}
            ($0001,$FF80,$FFFF,$FFFF,$3FFE),  {+9.9999999999954525270E-1}
            ($0000,$FFC0,$FFFF,$FFFF,$3FFE),  {+9.9999999999977262632E-1}
            ($0000,$FFE0,$FFFF,$FFFF,$3FFE),  {+9.9999999999988631316E-1}
            ($0000,$FFF0,$FFFF,$FFFF,$3FFE),  {+9.9999999999994315658E-1}
            ($0000,$FFF8,$FFFF,$FFFF,$3FFE),  {+9.9999999999997157829E-1}
            ($0000,$FFFC,$FFFF,$FFFF,$3FFE),  {+9.9999999999998578914E-1}
            ($0000,$FFFE,$FFFF,$FFFF,$3FFE),  {+9.9999999999999289457E-1}
            ($0000,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999644729E-1}
            ($8000,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999822364E-1}
            ($C000,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999911182E-1}
            ($E000,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999955591E-1}
            ($F000,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999977796E-1}
            ($F800,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999988898E-1}
            ($FC00,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999994449E-1}
            ($FE00,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999997224E-1}
            ($FF00,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999998612E-1}
            ($FF80,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999999306E-1}
            ($FFC0,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999999653E-1}
            ($FFE0,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999999826E-1}
            ($FFF0,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999999913E-1}
            ($FFF8,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999999957E-1}
            ($FFFC,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999999978E-1}
            ($FFFE,$FFFF,$FFFF,$FFFF,$3FFE),  {+9.9999999999999999989E-1}
            ($FFFF,$FFFF,$FFFF,$FFFF,$3FFE)); {+9.9999999999999999994E-1}
var
  m: integer;
begin
  if n>64 then sfc_etaint := 1.0
  else if n>=0 then sfc_etaint := extended(etahex[n])
  else begin
    m := 1-n;
    if odd(m) then sfc_etaint := 0
    else sfc_etaint := sfc_bernoulli(m)/m*exp2m1(m);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_eta(s: extended): extended;
  {-Return the Dirichlet eta function}
var
  z,t: extended;
const
  c2 = 0.3268629627944929957e-1;
  c1 = 0.15625 + 0.003618903742430971757; {sum correctly rounded for 16 bit}
  z0 = 2e-6;
begin
  if s>=65.0 then sfc_eta := 1.0
  else if (frac(s)=0.0) and (abs(s)<=MaxInt) then sfc_eta := sfc_etaint(round(s))
  else begin
    z := 1.0-s;
    if abs(z) <= z0 then sfc_eta := ln2 - (c2*z + c1)*z
    else if s >= -ETAEPS then sfc_eta := etam1pos(s)+1.0
    else if s <= -8.0 then begin
      {use eta(s) = (1-2^z)*zeta(s)}
      sfc_eta := sfc_zeta(s)*(1.0-exp2(z));
    end
    else begin
      {Use reflection formula for eta with s < 0}
      t := -exp2m1(z)/exp2m1(-s);
      t := t*cospi(0.5*z);
      t := t*sfc_gamma(z);
      t := t/power(pi,z);
      sfc_eta := t*etam1pos(z)+t;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_etam1(s: extended): extended;
  {-Return Dirichlet eta(s)-1}
begin
  if s < -ETAEPS then sfc_etam1 := sfc_eta(s)-1.0
  else sfc_etam1 := etam1pos(s);
end;


{---------------------------------------------------------------------------}
function sfc_zeta(s: extended): extended;
  {-Return the Riemann zeta function at s, s<>1}
var
  sc: extended;
  a,b: extended;
  ig: integer;
begin
  {Ref: Boost [19], file zeta.hpp}
  if IsNanOrInf(s) then begin
    if s=PosInf_x then sfc_zeta := 1.0
    else sfc_zeta := NaN_x;
    exit;
  end;
  if frac(s)=0.0 then begin
    if (s > -MaxBernoulli) and  (s <= MaxInt) then begin
      sfc_zeta := sfc_zetaint(round(s));
      exit;
    end;
  end
  else if abs(s) <= 1.6e-10 then begin
    {Small s, especially useful for s<0 to avoid reflection machinery}
    sfc_zeta := (-0.5) - LnSqrt2Pi*s;
    exit;
  end;
  sc := 1.0-s;
  if s<0.0 then begin
    if frac(0.5*s)=0.0 then sfc_zeta := 0.0
    else begin
      {compute Gamma(sc)/(2Pi)^sc}
      if sc>MaxGAMX-1 then begin
        {Not very accurate but does not overflow!  For s=-1753.5 this branch}
        {gives a rel. error of 2.0e-16 vs 1.8e-17 for the non-lngamma branch}
        a := sfc_lngammas(sc,ig);
        b := 2.0*sc*LnSqrt2Pi;
        {here Zeta(sc)=1}
        a := a-b;
        if a<Ln_MaxExt then a := exp(a)
        else a := PosInf_d;
        a := ig*a;
      end
      else begin
        a := sfc_gamma(sc);
        b := power(TwoPi,-sc);
        a := a*b;
        b := zetap(sc,s);
        a := a*b;
      end;
      {here a = zeta(sc)*Gamma(sc)/(2Pi)^sc}
      b := sinPi(0.5*s);
      sfc_zeta := 2.0*a*b;
    end;
  end
  else sfc_zeta := zetap(s,sc);
end;


{---------------------------------------------------------------------------}
function sfc_zeta1p(x: extended): extended;
  {-Return the Riemann zeta function at 1+x, x<>0}
begin
  {Ref: Boost [19], file zeta.hpp}
  if abs(x)<1.0 then sfc_zeta1p := zetap(1.0+x,-x)
  else sfc_zeta1p := sfc_zeta(1.0+x);
end;


{---------------------------------------------------------------------------}
function sfc_zetam1(s: extended): extended;
  {-Return Riemann zeta(s)-1, s<>1}
var
  t: extended;
begin
  if s <= 2.0 then sfc_zetam1 := sfc_zeta(s) - 1.0
  else begin
    if s >= 120.0 then sfc_zetam1 := exp2(-s)
    else begin
      t := 0.5*exp2(s);
      sfc_zetam1 := (1.0 + etam1pos(s)*t) / (t - 1.0);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_zetaint(n: integer): extended;
  {-Return zeta(n) for integer arguments, n<>1}
const
  znhex: array[2..63] of THexExtW = (
           ($9918,$983E,$3312,$D28D,$3FFF),  {+1.6449340668482264364}
           ($09C2,$8031,$0027,$99DD,$3FFF),  {+1.2020569031595942854}
           ($41B6,$3EC2,$9156,$8A89,$3FFF),  {+1.0823232337111381915}
           ($6DBD,$53E6,$0C76,$84BA,$3FFF),  {+1.0369277551433699263}
           ($247C,$0492,$4C26,$8238,$3FFF),  {+1.0173430619844491397}
           ($C46D,$A679,$96D0,$8111,$3FFF),  {+1.0083492773819228268}
           ($B746,$C31C,$9B57,$8085,$3FFF),  {+1.0040773561979443394}
           ($AB81,$C0B5,$CF9E,$8041,$3FFF),  {+1.0020083928260822144}
           ($CBF1,$D2DD,$9719,$8020,$3FFF),  {+1.0009945751278180854}
           ($9954,$F245,$318D,$8010,$3FFF),  {+1.0004941886041194645}
           ($9821,$D966,$1052,$8008,$3FFF),  {+1.0002460865533080483}
           ($DEF6,$E845,$0564,$8004,$3FFF),  {+1.0001227133475784892}
           ($9C2B,$5E56,$01C9,$8002,$3FFF),  {+1.0000612481350587048}
           ($11BF,$BCBF,$0097,$8001,$3FFF),  {+1.0000305882363070205}
           ($9178,$66F5,$8032,$8000,$3FFF),  {+1.0000152822594086518}
           ($A19A,$C1CD,$4010,$8000,$3FFF),  {+1.0000076371976378998}
           ($6E8A,$932A,$2005,$8000,$3FFF),  {+1.0000038172932649999}
           ($F9BA,$DB08,$1001,$8000,$3FFF),  {+1.0000019082127165539}
           ($A233,$9E2C,$0800,$8000,$3FFF),  {+1.0000009539620338727}
           ($ACA0,$34AE,$0400,$8000,$3FFF),  {+1.0000004769329867878}
           ($D9D9,$118C,$0200,$8000,$3FFF),  {+1.0000002384505027277}
           ($F138,$05D8,$0100,$8000,$3FFF),  {+1.0000001192199259653}
           ($CFFF,$01F2,$0080,$8000,$3FFF),  {+1.0000000596081890513}
           ($3A95,$00A6,$0040,$8000,$3FFF),  {+1.0000000298035035146}
           ($662E,$0037,$0020,$8000,$3FFF),  {+1.0000000149015548284}
           ($76B9,$0012,$0010,$8000,$3FFF),  {+1.0000000074507117898}
           ($2768,$0006,$0008,$8000,$3FFF),  {+1.0000000037253340248}
           ($0D18,$0002,$0004,$8000,$3FFF),  {+1.0000000018626597235}
           ($AF05,$0000,$0002,$8000,$3FFF),  {+1.0000000009313274324}
           ($3A56,$0000,$0001,$8000,$3FFF),  {+1.0000000004656629064}
           ($1372,$8000,$0000,$8000,$3FFF),  {+1.0000000002328311834}
           ($067B,$4000,$0000,$8000,$3FFF),  {+1.0000000001164155017}
           ($0229,$2000,$0000,$8000,$3FFF),  {+1.0000000000582077209}
           ($00B8,$1000,$0000,$8000,$3FFF),  {+1.0000000000291038504}
           ($003D,$0800,$0000,$8000,$3FFF),  {+1.0000000000145519218}
           ($0014,$0400,$0000,$8000,$3FFF),  {+1.0000000000072759598}
           ($0007,$0200,$0000,$8000,$3FFF),  {+1.0000000000036379796}
           ($0002,$0100,$0000,$8000,$3FFF),  {+1.0000000000018189896}
           ($0001,$0080,$0000,$8000,$3FFF),  {+1.0000000000009094948}
           ($0000,$0040,$0000,$8000,$3FFF),  {+1.0000000000004547474}
           ($0000,$0020,$0000,$8000,$3FFF),  {+1.0000000000002273737}
           ($0000,$0010,$0000,$8000,$3FFF),  {+1.0000000000001136868}
           ($0000,$0008,$0000,$8000,$3FFF),  {+1.0000000000000568434}
           ($0000,$0004,$0000,$8000,$3FFF),  {+1.0000000000000284217}
           ($0000,$0002,$0000,$8000,$3FFF),  {+1.0000000000000142108}
           ($0000,$0001,$0000,$8000,$3FFF),  {+1.0000000000000071054}
           ($8000,$0000,$0000,$8000,$3FFF),  {+1.0000000000000035527}
           ($4000,$0000,$0000,$8000,$3FFF),  {+1.0000000000000017764}
           ($2000,$0000,$0000,$8000,$3FFF),  {+1.0000000000000008882}
           ($1000,$0000,$0000,$8000,$3FFF),  {+1.0000000000000004441}
           ($0800,$0000,$0000,$8000,$3FFF),  {+1.0000000000000002220}
           ($0400,$0000,$0000,$8000,$3FFF),  {+1.0000000000000001110}
           ($0200,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000555}
           ($0100,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000278}
           ($0080,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000139}
           ($0040,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000069}
           ($0020,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000035}
           ($0010,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000017}
           ($0008,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000009}
           ($0004,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000004}
           ($0002,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000002}
           ($0001,$0000,$0000,$8000,$3FFF)); {+1.0000000000000000001}
begin
  if n>63 then sfc_zetaint := 1.0
  else if n>1 then sfc_zetaint := extended(znhex[n])
  else if n<0 then sfc_zetaint := sfc_bernoulli(1-n)/(n-1)+0.0  {avoid -0}
  else if n=0 then sfc_zetaint := -0.5
  else sfc_zetaint := Nan_x;
end;


{---------------------------------------------------------------------------}
function hz_a1(s,a: extended; var ok: boolean): extended;
  {-Return zetah(s,a) for a close to 1, s<>1 (normally s < 0)}
var
  x,t,f,z: extended;
  n: integer;
const
  MaxIter = 100;
begin
  hz_a1 := 0.0;
  ok := false;
  {compute sum from NIST[30], 25.11.10}
  z := sfc_zeta(s);
  if IsNanOrInf(z) then exit;
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
      t := sfc_zeta(t);
      t := t*f;
      z := z + t;
      {do not use t for exit test because zeta(s+n) may be very small}
    until (abs(f)<=eps_x*abs(z)) and (n>2);
  end;
  ok := true;
  hz_a1 := z;
end;


{---------------------------------------------------------------------------}
function hurwitz_formula(s,a: extended): extended;
  {-Compute zetah with range reduction and Hurwitz formula, a > 0, s < 0}
var
  n,k: longint;
  a2,f,h,r,t,tol,sh,z: extended;
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
  tol := 0.5*eps_x;
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
    if (z>0.0) and (z<MAXGAMX) then begin
      f := sfc_gamma(z);
      t := power(TwoPi, -z);
      t := (t*f)*(2.0*h);
    end
    else begin
      f := sfc_lngamma(z);
      t := ln(abs(2.0*h)) + f - z*ln(TwoPi);
      if t>= ln_MaxExt then t := PosInf_x
      else t := exp(t);
      if h<0.0 then t := -t;
    end;
    hurwitz_formula := t - r;
  end
end;


{---------------------------------------------------------------------------}
function bernpoly_intern(n: integer; x: extended): extended;
  {-Return the Bernoulli polynomial B_n(x), internal use: for small n only}
var
  i: integer;
  s,f: extended;
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
      if i and 1 = 0 then s := s + sfc_bernoulli(i)*f;
    end;
    if ref then s := -s;
    bernpoly_intern := s;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_zetah(s,a: extended): extended;
  {-Return the Hurwitz zeta function zetah(s,a) = sum(1/(i+a)^s, i=0..INF), s<>1, a>0}
var
  eps,q,r,t,u,w,z: extended;
  j,k,n: integer;
  ok: boolean;
const
  MaxIter = 40;
  eps1 = 1e-4;   {Threshold for hz_a1}
  S_HF = -8;     {Threshold for Hurzwitz formula}
begin
  {Initial reference: Cephes [7], file double\zeta.c}
  if IsNanOrInf(s) or IsNanOrInf(a) or (s=1.0) or (a<=0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_zetah := NaN_x;
    exit;
  end
  else if s=0.0 then begin
    sfc_zetah := 0.5 - a;       {NIST [30], 25.11.13}
    exit;
  end
  else if a=0.5 then begin
    if s < 16384.0 then begin
      sfc_zetah := exp2m1(s)*sfc_zeta(s);   {NIST [30], 25.11.11}
      exit;
    end;
  end
  else if a=1.0 then begin
    sfc_zetah := sfc_zeta(s);   {NIST [30], 25.11.2}
    exit;
  end
  else if a=2.0 then begin
    sfc_zetah := sfc_zetam1(s); {NIST [30], 25.11.2/3}
    exit;
  end;

  if s>0.0 then begin
    {Quick checks: avoid overflow/underflow}
    t := ln(a);
    t := -s*t;
    if t>ln_MaxExt then begin
      sfc_zetah := PosInf_x;
      exit;
    end
    else if t<ln_MinExt then begin
      {here s > 1, 1/a^s ~ 0, avoid underflow in some cases:}
      {zetah ~ a^(1-s)/(s-1) + 0.5/a^s = (a/(s-1) + 0.5)/a^s}
      {      = a^(1-s)/(s-1) + 0.5*a^(1-s)/a}
      u := power(a, 1-s);
      sfc_zetah := u/(s-1.0) + 0.5*u/a;
      exit;
    end;
  end
  else begin
    if (a>1.0) then begin
      if -s*ln(a-1.0) >ln_MaxExt then begin
        sfc_zetah := PosInf_x;
        exit;
      end;
    end;
    {Test if Bernoulli polynomials or Hurwitz formula can be used}
    if (s >= S_HF) and (frac(s)=0.0) then begin
      {s small negative integer: use 'naive' Bernoulli polynomial code}
      k := 1 - round(s);
      sfc_zetah := -bernpoly_intern(k,a)/k;  {NIST [30], 25.11.14}
      exit;
    end;
    if abs(a-1.0)<eps1 then begin
      {try expansion near a=1}
      sfc_zetah := hz_a1(s,a,ok);
      if ok then exit;
    end;
    if s < S_HF then begin
      sfc_zetah := hurwitz_formula(s,a);
      exit;
    end
  end;

  {Calculate z = zetah with the Euler-Maclaurin summation formula:}
  {See e.g. http://functions.wolfram.com/10.02.06.0020.01         }
  {z = sum1(s,a,n) + (a+n)^(1-s)/(s-1) - 1/2/(a+n)^s + sum2(s,a,n)}
  {with sum1(a,n) = sum(1/(a+k)^s, k=0..n) and                    }
  {sum2 = sum(B_2k/(2k)! * s(s+1)..(s+2k)/(a+n)^(s+2k+1),k=0..INF)}

  eps := 0.5*eps_x;
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
      sfc_zetah := z;
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
      if j=NBOF+1 then eps := 1.5*eps_x;
      u := t;
      t := sfc_bernoulli(k+2);  {TP5 fix}
      t := t/sfc_fac(k+2);
      t := (q*t)*r;
      if abs(u)<=abs(t) then break;
    end
    else t := (q*extended(BoFHex[j]))*r;
    z := z + t;
    if abs(t) < eps*abs(z) then begin
      sfc_zetah := z;
      exit;
    end;
    q := q*(s+k+1);
    r := r/w;
    inc(k,2);
  end;

  {No convergence}
  sfc_zetah := z;
  {$ifdef debug}
    writeln('sfc_zetah - Euler-Maclaurin issue: ',s:21,' ',a:21);
  {$endif}
  if RTE_NoConvergence>0 then RunError(byte(RTE_NoConvergence));
end;


{---------------------------------------------------------------------------}
function bernpoly_rev(n: integer; x: extended): extended;
  {-Return the Bernoulli polynomial B_n(x), internal use: n>1, small x}
var
  i: integer;
  s,t,u,f: extended;
begin
  {This code runs NIST[30] 24.2.5 backwards}
  s := sfc_bernoulli(n);
  f := 1.0;
  t := PosInf_x;
  for i:=n-1 downto 0 do begin
    f := f*x*(i+1)/(n-i);
    if (i and 1 = 0) or (i=1) then begin
      {B_i <> 0}
      u := t;
      t := s;
      s := t + sfc_bernoulli(i)*f;
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
function sfc_bernpoly(n: integer; x: extended): extended;
  {-Return the Bernoulli polynomial B_n(x), 0 <= n <= MaxBernoulli}
var
  a: extended;
  neg: boolean;
begin
  if IsNan(x) or (n<0) or (n>MaxBernoulli) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_bernpoly := Nan_x;
    exit;
  end;
  if n=0 then sfc_bernpoly := 1.0
  else if n=1 then sfc_bernpoly := x-0.5
  else if (x=0.0) or (x=1.0) then sfc_bernpoly := sfc_bernoulli(n)  {n<>1 !}
  else if x=0.5 then begin
    {B_n(1/2) = (2^(1-n)-1) B_n)}
    if n>65 then a := -1
    else a := exp2m1(1-n);
    sfc_bernpoly := a*sfc_bernoulli(n)
  end
  else if (x=-1.0) then begin
    {B_n(-1) := B_n + (-1)^n*n}
    if odd(n) then sfc_bernpoly := -n
    else sfc_bernpoly := sfc_bernoulli(n) + n;
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
        if x>=0.0 then a := sfc_zetah(1-n, x)
        else begin
          a := sfc_zetah(1-n, 1.0-x);
          if odd(n) then neg := not neg;
        end;
        a := -a*n;
      end;
    end;
    if neg then a := -a;
    sfc_bernpoly := a;
  end;
end;


{---------------------------------------------------------------------------}
type
  tphi_rec = record
               z,s,a,g: extended;  {g=1/Gamma(s)}
             end;
  pphi_rec = ^tphi_rec;


{---------------------------------------------------------------------------}
function phif(x: extended; param: pointer): extended; {$ifdef BIT16} far;{$endif}
  {-LerchPhi integrand: x^(s-1)*exp(-ax)/(1-z*exp(-x)/Gamma(s)}
var
  d,p,f: extended;
begin
  {NIST[30], 25.14.5}
  with pphi_rec(param)^ do begin
    p := (s-1.0)*ln(x);
    if p >= ln_MaxExt then begin
      {power would overflow, use logarithms}
      p := exp(p - a*x - sfc_lngamma(s));
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
function phi_int(z,s,a: extended): extended;
  {-LerchPhi from integral representation, z < 0, 0 < s < MAX_GAMX, a > 0}
var
  res,aerr: extended;
  ierr: integer;
  para: tphi_rec;
  neval: longint;
const
  eps  = 4e-19;
  eps2 = 1e-15;
begin
  phi_int := Nan_x;
  para.z := z;
  para.s := s;
  para.a := a;
  para.g := sfc_rgamma(s);
  if para.g=0.0 then begin
    {No reliable computation, Gamma(s) > MaxExtended}
    exit;
  end;
  {Use integral formula, see Wood[73], (1.1) or }
  {http://functions.wolfram.com/10.08.07.0001.01}
  sfc_intdei_p({$ifdef FPC_ProcVar}@{$endif}phif, @para, 0.0, eps, res, aerr, neval, ierr);
  if ierr=0 then phi_int := res
  else begin
    if (ierr=3) and (aerr <= eps2*abs(res)) then begin
      {Max. iterations with reduced accuracy}
      phi_int := res;
    end
  end;
end;


{---------------------------------------------------------------------------}
function lphi_aj(z, s, a: extended): extended;
  {-Lerch Phi via Aksenov/Jentschura convergence acceleration}
const
  IMAX  = 50;  {Max. iteration count for CNCT case }
  IMAX2 = 90;  {Max. iteration for direct summation}
var
  num,den,SAj : array[0..IMAX] of extended;

  {---------------------------------------------------------------------------}
  function vwaj(j: integer): extended;
     {-Compute the van Wijngaarden quantities A_j from b^j_k}
  var
    sum, bjk, z2ind: extended;
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
    until abs(bjk) <= eps_x*abs(sum);
    vwaj := sum;
  end;

var
  j,i,k,sign: integer;
  sn,eps0,eps,skn,skn0,omega,fk,est,t,x: extended;

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
    if (z >= 1.0) or (z < -1.0) or (a <= 0.0) or (s < -1.0)  then begin
      {Return NaN, this should not happen if called by sfc_lerch}
      lphi_aj := NaN_x;
      exit;
    end;
  {$endif}

  if x <= MinExtended then begin
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
  skn0 := skn;  {Avoid hint/warning: Keep Delphi happy and use skn}
  sn   := 0.0;

  {omega is the next term of a partial sum of defining power series for direct}
  {summation, of van Wijngaarden transformed series for CNCT, and also becomes}
  {a remainder estimate in the delta transformation in CNCT.}

  {For z<=0.5 van Wijngaarden transformation is not used hence no calls to aj()}
  if z <= 0.5 then omega := power(a, -s)
  else omega := vwaj(0);

  if IsInf(omega) then begin
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
      if abs(omega) <= MinExtended then begin
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
      if abs(est) < 8.0*eps_x*abs(skn) then exit;
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
function lphi_aea(z,s,a: extended): extended;
  {-Lerch Phi, asymptotic expansion for large a>0}
var
  n: integer;
  e,r,t,x,y: extended;
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
    e := 0.5*eps_x;
    repeat
      {y holds the terms (s)_n/(n!*a^n)}
      y := y*(s/n)/a;
      t := sfc_polylog(-n,z);
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
function lphi_sum(z,s,a: extended): extended;
  {-Lerch Phi, direct summation sum(z^n/(n+a)^s, n=0...)}
var
  n: integer;
  r,t,x,e: extended;
const
  NMAX=1000;
begin
  n := 0;
  x := 1.0; {x = z^n}
  r := 0.0;
  e := 0.5*eps_x;
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
function sfc_lerch(z,s,a: extended): extended;
  {-Return the Lerch transcendent Phi(z,s,a) = sum(z^n/(n+a)^s, n=0..INF), z <= 1, s >= -1, a >= 0}
var
  x,t: extended;
begin
  {Phi(z,s,a) = sum(z^n/(n+a)^s, n=0...)}
  if IsNanOrInf(z) or IsNanOrInf(s) or IsNanOrInf(a)
     or (z>1.0) or (s<-1.0) or (a<0.0) then
  begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_lerch := NaN_x;
    exit;
  end;

  if s=0.0 then begin
    {geometric series Phi = sum(z^n)}
    sfc_lerch := 1.0/(1.0-z);
    exit;
  end
  else if z=0.0 then begin
    sfc_lerch := power(a,-s);
    exit;
  end
  else if z=1.0 then begin
    {continuation from s>1}
    sfc_lerch := sfc_zetah(s,a);
    exit;
  end
  else if a=0.0 then begin
    {http://functions.wolfram.com/10.06.03.0049.01}
    sfc_lerch := sfc_polylogr(s,z);
    exit;
  end;

  x := abs(z);

  if s<0.0 then begin
    {Note: For s<0 no asymptotic expansion for large a}
    if (x>0.5) or (z >= 0.875) then begin
      {accelerated summation}
      sfc_lerch := lphi_aj(z,s,a);
    end
    else begin
      {direct summation}
      sfc_lerch := lphi_sum(z,s,a);
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
      sfc_lerch := lphi_aj(z,s,a);
    end
    else begin
      {Use integral NIST[30], 25.14.5 for x<-1, a>0, s>0}
      sfc_lerch := phi_int(z,s,a)
    end;
    exit;
  end;

  t := s*ln(a);
  if t > ln_MaxExt then begin
    {1/a^s = 0 accurate to extended precision}
    sfc_lerch := 0;
  end
  else if (x <= 0.3) and (a >= 100.0) and (s <= 100.0) then begin
    {Use asymptotic expansion for large a, here only if s and x are not too large}
    sfc_lerch := lphi_aea(z,s,a);
  end
  else if (s >= 30.0) or (t < -43.7) then begin
    if (x>0.5) and (s < 0.1*a) then begin
      {accelerated summation because direct summation needs many iterations}
      sfc_lerch := lphi_aj(z,s,a);
    end
    else begin
      {direct summation if s is large enough or 1/a^s > 1/eps_x}
      sfc_lerch := lphi_sum(z,s,a);
    end;
  end
  else if (1.0-x)*a >= 50.0  then begin
    {Remaining case for asymptotic expansion for large a}
    sfc_lerch := lphi_aea(z,s,a);
    exit;
  end
  else begin
    if x<=0.5 then begin
      {direct summation, lphi_aj would do the same but a lot more complicated}
      sfc_lerch := lphi_sum(z,s,a);
    end
    else begin
      {accelerated summation}
      sfc_lerch := lphi_aj(z,s,a);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_dbeta(s: extended): extended;
  {-Return the Dirichlet beta function sum((-1)^n/(2n+1)^s, n=0..INF)}
var
  t: extended;
const
  d0: THexExtW = ($B90E,$F218,$0F6B,$C87F,$3FFD); {0.39159439270683677647194534689911102809021011577}
begin
  {http://en.wikipedia.org/wiki/Dirichlet_beta_function   }
  {http://mathworld.wolfram.com/DirichletBetaFunction.html}
  {This is also known as the Catalan beta function        }

  if abs(s) <= 1e-10 then sfc_dbeta := 0.5 + extended(d0)*s
  else if s<0.0 then begin
    if frac(0.5*s)=-0.5 then sfc_dbeta := 0.0
    else begin
      t := cosPi(0.5*s);
      s := 1.0 - s;
      t := sfc_dbeta(s)*t;
      sfc_dbeta := sfc_gamma(s)*t*power(Pi_2,-s);
    end;
  end
  else if s>=40.5 then sfc_dbeta := 1.0
  else if s>=27.6 then sfc_dbeta := 1.0 - exp3(-s)
  else if s>=22.8 then sfc_dbeta := 1.0 - exp3(-s) + exp5(-s)
  else begin
    sfc_dbeta := sfc_lerch(-1.0,s,0.5)*exp2(-s);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_dlambda(s: extended): extended;
  {-Return the Dirichlet lambda function sum(1/(2n+1)^s, n=0..INF), s<>1}
begin
  {http://mathworld.wolfram.com/DirichletLambdaFunction.html}
  sfc_dlambda := -exp2m1(-s)*sfc_zeta(s);
end;


{---------------------------------------------------------------------------}
function sfc_lchi(s, x: extended): extended;
  {-Return Legendre's Chi-function chi(s,x); s>=0; real part if |x| > 1}
var
  a,b,z: extended;
begin
  if IsNanOrInf(s) or IsNanOrInf(x) or (s<0.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_lchi := NaN_x;
    exit;
  end;
  z := abs(x);
  if z=0.0 then sfc_lchi := 0.0
  else if s=0.0 then sfc_lchi := x/((1.0+x)*(1.0-x))
  else if s=1.0 then begin
    if z <= 1.0 then sfc_lchi := arctanh(x)
    else begin
      {Compute real part of arctanh(z), see manual for carctanh}
      a := 2.0/(1.0-z);
      b := 0.25*ln1p(z*sqr(a));
      sfc_lchi := b*isign(x);
    end;
  end
  else if z <= 1.0 then begin
    {chi(s,x) = sum(x^(2n+1)/(2n+1)^s, n=0..INF)}
    {chi(s,x) = 2^(-s)*x*Phi(x^2,s,1/2) = 0.5*(Li_s(x)-Li_s(-x))}
    {see e.g. http://en.wikipedia.org/wiki/Legendre_chi_function}
    {or  http://mathworld.wolfram.com/LegendresChi-Function.html}
    if s>=22.8 then begin
      if s>=40.5 then sfc_lchi := x
      else begin
        z := sqr(x);
        if s>=27.6 then sfc_lchi := x*(1.0 + z*exp3(-s))
        else sfc_lchi := x*(1.0 + z*(exp3(-s) + z*exp5(-s)));
      end;
    end
    else begin
      sfc_lchi := sfc_lerch(x*x,s,0.5)*x*exp2(-s);
    end;
  end
  else begin
    {|x| > 1}
    a := sfc_polylogr(s, z);
    b := sfc_polylogr(s,-z);
    sfc_lchi := 0.5*(a-b)*isign(x);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_beint(s, x: extended): extended;
  {-Return the Bose-Einstein integral of real order s >= -1, real part if x > 1}
const
  XL = 50;
begin
  if IsNanOrInf(s) or IsNanOrInf(x) or (s < -1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_beint := NaN_x;
    exit;
  end;
  if x > XL then begin
    {Use asymptotic expansion of polylog without computing exp(x))}
    sfc_beint := liae_pos(s+1.0,x)
  end
  else sfc_beint := sfc_polylogr(s+1.0,exp(x));
end;


{---------------------------------------------------------------------------}
function fd_asymp_exp_real(q, x: extended): extended;
  {-Fermi-Dirac asymptotic expansion for real order q with q > -1, x >> q}
var
  p,s,t,r,f,z: extended;
  k: integer;
begin
  {Asymptotic expansion Goano [47], (9)}
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
    if n>0 then sfc_fermi_dirac := fd_asymp_exp_real(n, x)
    else case n of
        0: sfc_fermi_dirac := x;     {ln(1+e^x)   -> x}
       -1: sfc_fermi_dirac := 1.0;   {e^x/(1+e^x) -> 1}
      else sfc_fermi_dirac := 0.0;
    end;
  end
  else begin
    z := -exp(x);
    sfc_fermi_dirac := -sfc_polylog(n+1, z);
  end;
end;


{---------------------------------------------------------------------------}
function sfc_fdr(s,x: extended): extended;
  {-Return the Fermi-Dirac integral of real order s >= -1}
var
  ae,cs,z,y: extended;
begin
  if IsNanOrInf(s) or IsNanOrInf(x) or (s < -1.0) then begin
    {$ifopt R+}
      if RTE_ArgumentRange>0 then RunError(byte(RTE_ArgumentRange));
    {$endif}
    sfc_fdr := NaN_x;
    exit;
  end;
  if frac(s)=0.0 then begin
    sfc_fdr := sfc_fermi_dirac(round(s),x);
  end
  else begin
    if x < ln_MinExt then sfc_fdr := 0.0
    else if x >= 30.0 + abs(s) then begin
      ae := fd_asymp_exp_real(s,x);
      cs := cospi(s);
      if (cs <> 0.0) and (x >= -ln_MaxExt) then begin
        y  := exp(-x);
        z  := -sfc_polylogr(s + 1.0, -y);
        ae := ae + z*cs;
      end;
      sfc_fdr := ae;
    end
    else begin
      y := exp(x);
      z := s+1.0;
      if x < z*ln2 - 45.0 then sfc_fdr := y
      else sfc_fdr := -sfc_polylogr(z, -y);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sfc_fdp05(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(1/2,x)}
begin
  sfc_fdp05 := sfc_fdr(0.5,x);
end;


{---------------------------------------------------------------------------}
function sfc_fdp15(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(3/2,x)}
begin
  sfc_fdp15 := sfc_fdr(1.5,x);
end;


{---------------------------------------------------------------------------}
function sfc_fdp25(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(5/2,x)}
begin
  sfc_fdp25 := sfc_fdr(2.5,x);
end;


{---------------------------------------------------------------------------}
function sfc_fdm05(x: extended): extended;
  {-Return the complete Fermi-Dirac integral F(-1/2,x)}
begin
  sfc_fdm05 := sfc_fdr(-0.5,x);
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
  if IsNanOrInf(x) or IsNanOrInf(r) or (x < -1.0) then begin
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
      t := 1.0 + a + s;
      if n>4 then t := exp5(-r) + t;
      case n of
          3: sfc_harmonic2 := t;
          4: sfc_harmonic2 := t + a*a;
          5: sfc_harmonic2 := t + a*a;
        else sfc_harmonic2 := t + a*(a + s);
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

