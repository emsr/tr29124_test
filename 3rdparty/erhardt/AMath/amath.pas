unit AMath;

{Accurate floating point math unit}

interface

{$i std.inc}

{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Accurate floating point math unit

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  $define use_fast_exp if the roundmode-safe exp code should
                  NOT be used. This assumes that the exp/exp3/5/7/10 routines
                  should be called with rounding to nearest, otherwise in
                  rare case there may be wrong results.


 REFERENCES    :  References used in this unit, main index in amath_info.txt/references

                  [1] [HMF]: M. Abramowitz, I.A. Stegun. Handbook of Mathematical Functions. Dover, 1970
                      http://www.math.sfu.ca/~cbm/aands/
                  [2] Intel, IA-32 Architecture Software Developer's Manual
                      Volume 2A: Instruction Set Reference, A-M
                  [3] D. Goldberg, What Every Computer Scientist Should Know About
                      Floating-Point Arithmetic, 1991; extended and edited reprint via
                      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.22.6768
                  [4] ISO/IEC 10967-2, Information technology: Language independent
                      arithmetic, Part 2: Elementary numerical functions
                      http://standards.iso.org/ittf/PubliclyAvailableStandards/c024427_ISO_IEC_10967-2_2001(E).zip
                  [5] FDLIBM 5.3 (Freely Distributable LIBM), developed at
                      Sun Microsystems,  see http://www.netlib.org/fdlibm/ or
                      http://www.validlab.com/software/fdlibm53.tar.gz
                  [6] K.C. Ng, "Argument Reduction for Huge Arguments: Good to theLast Bit",
                      Technical report, SunPro, 1992. Available from
                      http://www.validlab.com/arg.pdf
                  [7] Cephes Mathematical Library, Version 2.8
                      http://www.moshier.net/#Cephes or http://www.netlib.org/cephes/
                  [8] T. Ogita, S.M. Rump, and S. Oishi, Accurate sum and dot product,
                      SIAM J. Sci. Comput., 26 (2005), pp. 1955-1988. Available as
                      http://www.ti3.tu-harburg.de/paper/rump/OgRuOi05.pdf
                  [9] N.J. Higham, Accuracy and Stability of Numerical Algorithms,
                      2nd ed., Philadelphia, 2002
                      http://www.maths.manchester.ac.uk/~higham/asna/
                 [25] I. Smith, Examples.xls/txt Version 3.3.4, Personal communication, 2010
                 [32] D.E. Knuth: The Art of computer programming; Volume 1, Fundamental
                      Algorithms, 3rd ed., 1997; Volume 2, Seminumerical Algorithms, 3rd ed., 1998.
                 [46] S. Graillat, P. Langlois, N Louvet, Compensated Horner Scheme,
                      Research Report No RR2005-04, 2005, Université de Perpignan.
                      Available from http://www-pequan.lip6.fr/~graillat/papers/rr2005-04.pdf or
                      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.81.2979
                 [64] T.J. Dekker, A Floating-point technique for extending the available
                      precision. Numerische Mathematik, 18, 224-242, 1971. Available from
                      http://www.digizeitschriften.de/en/dms/toc/?PPN=PPN362160546_0018
                 [65] S. Linnainmaa, Software for Doubled-Precision Floating-Point Computations.
                      ACM TOMS 7 (1981), pp. 272-283.



 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     04.10.09  W.Ehrhardt  Initial BP7 BASM version for ln1p
 0.11     04.10.09  we          Pure Pascal for ln1p
 0.12     04.10.09  we          BASM version for exmp1
 0.13     04.10.09  we          Pure Pascal for exmp1
 0.14     04.10.09  we          Machine epsilons, case y < eps_x in pp expm1
 0.15     04.10.09  we          Simplified BASM expm1
 0.16     17.10.09  we          exp, sinh, cosh
 0.17     17.10.09  we          Non-BASM code, needs std.inc
 0.18     18.10.09  we          tanh, fix pp expm1, some consts
 0.19     20.10.09  we          Pure Pascal: accurate exp, simplified expm1
 0.20     21.10.09  we          RTR 205 for pp exp overflow, more consts
 0.21     22.10.09  we          log10, log2, logN
 0.22     24.10.09  we          arctan for TP5, arctan2
 0.23     25.10.09  we          intpower and power
 0.24     25.10.09  we          arccos and arcsin
 0.25     26.10.09  we          exp2
 0.26     26.10.09  we          tan, cot, sec, csc
 0.27     27.10.09  we          coth, csch, sech
 0.28     28.10.09  we          arccosh, arccosh1p, arcsinh, arctanh
 0.29     28.10.09  we          arccot, arccsc, arcsec
 0.30     29.10.09  we          arccoth, arccsch, arcsech
 0.31     02.11.09  we          arccotc
 0.32     03.11.09  we          hypot, frexp, ldexp
 0.33     04.11.09  we          ilogb, modf, floor/x, ceil/x
 0.34     04.11.09  we          IsInf, IsNaN, IsNaNorInf
 0.35     04.11.09  we          longint($80000000) fix for FPC/D4+
 0.36     06.11.09  we          bugfix cot, always use Pascal ln1p
 0.37     06.11.09  we          improved arccosh, arcsech
 0.38     07.11.09  we          improved arccoth
 0.39     07.11.09  we          improved arcsec
 0.40     07.11.09  we          improved arccsc
 0.41     07.11.09  we          improved tanh
 0.42     08.11.09  we          improved arccsch
 0.43     09.11.09  we          improved sinh, fixed arcsec for x<0
 0.44     13.11.09  we          improved arccosh1p, arcsinh, arctanh
 0.45     15.11.09  we          exp10, bugfix intpower
 0.46     18.11.09  we          improved arccotc
 0.47     18.11.09  we          fix csch,sech for abs(x) > ln(MaxExtended)
 0.48     24.11.09  we          fix D2 internal error with fstp [@result]
 0.49                           intentionally left blank :)
 0.50     28.11.09  we          _tan, _cot (basic versions, OK for |arg| < Pi/4}
 0.51     28.11.09  we          sin,cos,tan,cot with rempio2
 0.52     29.11.09  we          Remove temp integer in sin,cos
 0.53     29.11.09  we          Type THexDblW
 0.54     03.12.09  we          More accurate power (based on Cephes)
 0.55     03.12.09  we          frexpd, scalbn
 0.56     03.12.09  we          fix BASM ln1p, improve arccosh/1p (ac_help)
 0.57     04.12.09  we          arcsech with ac_help
 0.58     05.12.09  we          AMath_Version, improved arccos and arcsin
 0.59     05.12.09  we          sum2, sum2x, dot2, dot2x
 0.60     06.12.09  we          integrated rempio2, tan uses _cot
 0.61     06.12.09  we          sum/dot ranges 0..n-1
 0.62     06.12.09  we          sincos, _sincos
 0.63     08.12.09  we          norm2, ssq, sumsqr
 0.64     09.12.09  we          improved non-basm ldexp
 0.65     09.12.09  we          copysign, fmod
 0.66     10.12.09  we          ldexpd
 0.67     10.12.09  we          PolEval, PolEvalX, non-basm exp10 from Cephes
 0.68     11.12.09  we          Mean, MeanAndStdDev
 0.69     11.12.09  we          improved sinh
 0.70     12.12.09  we          IsInfD, IsNaND, IsNaNorInfD
 0.71     12.12.09  we          nextd/x, predd/x, succd/x
 0.72     15.12.09  we          Get/Set8087CW, Get/SetRoundMode, Get/SetPrecisionMode, Get/SetExceptionMask
 0.73     17.12.09  we          Ext2Hex, Dbl2Hex, MinExtended via MinExtHex
 0.74     18.12.09  we          fisEQd/x, fisNEd/x
 0.75     19.12.09  we          Min/MaxDouble, Min/MaxSingle via Hex constants
 0.76     19.12.09  we          complete rewrite of predd/x, succd/x using bit manipulation
 0.77     19.12.09  we          improved ln1p, ac_help, arctanh
 0.78     21.12.09  we          improved arcsinh, use improved ln1p for non-basm
 0.79     25.12.09  we          fix fisNEd, frexpd
 0.80     28.12.09  we          rem_int2, simplified rem_pio2_cw, changed < Pi/4 to <= Pi/4
 0.81     28.12.09  we          cosPi, sinPi, sincosPi
 0.82     06.01.10  we          expx2, Sqrt_TwoPi
 0.83     26.01.10  we          early outs for expx2
 0.84     31.01.10  we          sqrt1pm1, improved arctanh
 0.85     07.02.10  we          allow denormal input for logarithms
 0.86     07.02.10  we          cbrt, nroot
 0.87     09.02.10  we          NoBASM routines: exp2, cbrt
 0.88     16.02.10  we          Euler's Gamma constant
 0.89     02.03.10  we          CSEvalX
 0.90     13.03.10  we          maxd, maxx, mind, minx
 0.91     03.04.10  we          Sqrt(Min/MaxExtended)
 0.92     03.04.10  we          Hex values in separate const declarations
 0.93     03.04.10  we          exp3
 0.94     09.04.10  we          sinc, sincPi
 0.95     14.04.10  we          ln1pmx, typed const with THexExtW
 0.96     15.04.10  we          ln_MaxExt, ln_MinExt
 0.97     20.04.10  we          rem_2pi/rem_2pi_sym, PiSqr/LnSqrt2Pi
 0.98     20.04.10  we          fix rem_pio2_ph for x=0.0
 0.99     07.06.10  we          meanx, MeanAndStdDevX
 1.00     08.06.10  we          rms, rmsx
 1.01     18.06.10  we          exprel
 1.02     01.07.10  we          cosf1
 1.03     02.07.10  we          coshm1
 1.04     03.07.10  we          mssqd, mssqx
 1.05     03.07.10  we          ssdev/x, psdev/x, svar/x, pvar/x
 1.06     18.07.10  we          powm1, pow1pm1
 1.07     04.08.10  we          More accurate ln1pmx
 1.08     27.08.10  we          More accurate power function (table with up to 512 entries)
 1.09     07.09.10  we          Ext2Dbl
 1.10     08.09.10  we          Bugfix frac/int/trunc for Delphi 2..5, FPC1
 1.11     12.09.10  we          Bugfix round for FPC1
 1.12     15.09.10  we          exp2m1
 1.13     29.09.10  we          types TFuncX and TFuncD
 1.14     29.09.10  we          TFuncX/D compatible int,frac,arctan,sqrt
 1.15     03.10.10  we          Interfaced mssqd, mssqx ifndef CONST
 1.16     03.10.10  we          moment/momentx
 1.17     03.10.10  we          cbrt ifdef changed from BASM16 to BASM
 1.18     06.10.10  we          SqrtPi = sqrt(Pi)
 1.19     23.10.10  we          power: handle 0^(-|y|) and avoid intpower overflow
 1.20     09.01.11  we          simplified succd/predd
 1.21     11.01.11  we          Single functions: IsInfS,IsNaNS,IsNaNorInfS,fisEQs,fisNEs,succs,preds,Sgl2Hex
 1.22     16.01.11  we          roundmode-safe exp/exp10/exp3 code
 1.23     17.01.11  we          roundmode-safe arccos
 1.24     22.01.11  we          frexps, ldexps, mins, maxs, copysigns
 1.25     03.02.11  we          renamed cosf1 to vers
 1.26     03.02.11  we          haversine hav(x)
 1.27     04.02.11  we          coversine covers(x) = 1 - sin(x)
 1.28     05.02.11  we          inverse haversine archav(x)
 1.29     05.02.11  we          Gudermannian gd(x)
 1.30     05.02.11  we          inverse Gudermannian function arcgd(x)
 1.31     13.02.11  we          functions ulpd/s/x
 1.32     19.03.11  we          isign
 1.33     25.03.11  we          RandG, RandG01
 1.34     29.03.11  we          Inline Get8087CW, Set8087CW for VER5X
 1.35     21.04.11  we          angle2
 1.36     22.04.11  we          orient2d, area_triangle, in_triangle/_ex
 1.37     28.04.11  we          improved rem_pio2
 1.38     30.04.11  we          improved rem_2pi
 1.39     24.05.11  we          case x=0 in rem_pio2_cw
 1.40     12.06.11  we          Hex2Ext, Hex2Dbl, Hex2Sgl
 1.41     13.06.11  we          hypot3
 1.42     13.06.11  we          improved vers, covers
 1.43     13.06.11  we          simplified sin,cos,tan,cot
 1.44     13.06.11  we          improved coshm1
 1.45     14.06.11  we          improved Hex2Float
 1.46     23.06.11  we          DegToRad, RadToDeg, Pi_180
 1.47     04.09.11  we          arccsch(x) for very small x
 1.48     05.09.11  we          arcsech(x) for very small x
 1.49     06.09.11  we          arccsc(x) for very large x
 1.50     22.12.11  we          exp5, exp7
 1.51     28.02.12  we          rint
 1.52     29.02.12  we          fix VP rmNearest to round to even for half-integers
 1.53     10.04.12  we          fix NOBASM exp
 1.54     11.04.12  we          Langevin function
 1.55     26.04.12  we          PolEvalEE/X
 1.56     27.04.12  we          PolEvalCHE/X
 1.57     03.06.12  we          Power: use intpower only for |y| < 8192
 1.58     07.06.12  we          Fix exp/3/5/7/10 for very large x
 1.59     29.06.12  we          const sqrt_epsh = sqrt(eps_x/2)
 1.60     20.02.13  we          Special cases |y|=1/2 in power
 1.61     22.02.13  we          Constants succx0Hex, succx0, ln_succx0
 1.62     05.04.13  we          frexp/d/s returns non-zero exponents for denormal
 1.63     05.04.13  we          ilogb returns floor(log2()) for denormal
 1.64     06.04.13  we          denormal for non-Basm: exp,exp2,exp10,ln,cbrt
 1.65     06.04.13  we          power: denormal support and intpower for y<2048
 1.66     07.04.13  we          exp10m1
 1.67     07.04.13  we          log2p1, log10p1
 1.68     08.04.13  we          denormal for non-Basm fmod
 1.69     25.04.13  we          arcsec for large x
 1.70     19.05.13  we          const THREE: double=3.0;  Used for circumventing some 'optimizations'
 1.71     29.06.13  we          Power returns +-INF for overflow
 1.72     17.08.13  we          expmx2h = exp(-0.5*x^2)
 1.73     03.10.13  we          Degrees versions of trig / invtrig functions
 1.74     04.10.13  we          trig_deg: more accurate degree reduction
 1.75     09.10.13  we          ln1p
 1.76     10.10.13  we          coshm1
 1.77     10.10.13  we          tanh
 1.78     10.10.13  we          expm1
 1.79     12.10.13  we          arcsinh
 1.80     12.10.13  we          arccsch, sinh, ln1pmx
 1.81     12.10.13  we          arcgd
 1.82     15.10.13  we          Sqrt2 as (hex)extended constant
 1.83     24.10.13  we          Improved degrees versions of trig functions
 1.84     25.10.13  we          IsInf/Nan function with {$ifdef HAS_INLINE} inline;{$endif}
 1.85     27.10.13  we          remainder
 1.86     20.11.13  we          sinhcosh, improved sinh/cosh
 1.87     26.03.14  we          LnPi = ln(Pi) as (hex)extended constant
 1.88     21.04.14  we          Improved Langevin
 1.89     29.05.14  we          PolEvalS
 1.90     07.06.14  we          ln1mexp
 1.91     20.06.14  we          logit
 1.92     05.07.14  we          basic ext2 routines
 1.93     06.07.14  we          pow1p
 1.94     06.07.14  we          compound
 1.95     06.10.14  we          logistic
 1.96     06.11.14  we          Small Bernoulli numbers as Hex constants from sfBasic
 1.97     26.12.14  we          Typed const one_x=1 to avoid some FPC problems
 1.98     05.01.15  we          Minor changes (editorial, some FP constants, ifdef BIT32 -> ifndef BIT16)
 1.99     09.01.15  we          xxto2x uses TwoSum instead of FastTwoSum
 2.00     10.01.15  we          pi2x
 2.01     26.02.15  we          powpi2k, powpi
 2.02     03.04.15  we          special case y=-1 in power
 2.03     31.05.15  we          versint
 2.04     19.06.15  we          sinhc
 2.05     25.06.15  we          sinhmx
 2.06     24.07.15  we          avoid ilogb in power if y too large for intpower
 2.07     25.08.15  we          Constant SIXX = 6.0 to avoid some 'optimizations'
 2.08     25.04.16  we          arccos1m = arccos(1-x)
 2.09     27.09.16  we          sqrt1pmx = sqrt(1+x*x)-x
 2.10     30.10.16  we          lncosh
 2.11     07.06.17  we          lnsinh
 2.12     08.06.17  we          tanPi
 2.13     21.06.17  we          pow1pf (without ext2)
 2.14     29.06.17  we          Code removed for TP5-TP6, TPW1-D1
 2.15     29.06.17  we          PolEvalC, PolEvalDeriv
 2.16     01.07.17  we          PolEvalCHEDer
 2.17     24.07.17  we          ln1pexp
 2.18     27.07.17  we          logaddexp, logsubexp
 2.19     27.07.17  we          interfaced xmul12
 2.20     27.07.17  we          xadd12
 2.21     28.07.17  we          xdivrem
 2.22     29.07.17  we          isRMNearest
 2.23     16.08.17  we          Removed old non-BASM code
 2.24     30.10.17  we          fma_x
 2.25     02.12.17  we          Suppress warnings: Local variable does not seem to be initialized
 2.26     01.01.18  we          TFuncX/DP: Functions with pointer to parameter(s)
 2.27     19.02.18  we          CSEvalXDer
 2.28     20.02.18  we          Langevin function moved to sfMisc
 2.29     19.04.18  we          Fixes for FPC311
 2.30     17.05.18  we          comprel
 2.31     07.06.18  we          eps_x as THexExtW (needed for FPC2.6.x)
 2.32     22.06.18  we          FPU control $ifdef moved to implementation
 2.33     27.06.18  we          Anti-FPC311 quotient in versint
 2.34     05.07.18  we          SafeDivX
 2.35     06.07.18  we          CoeffVar/X (coefficient of variation)
 2.36     10.07.18  we          interfaced sqr12x
 2.37     23.07.18  we          improved intpower, old ASM function renamed to intpower0
 2.38     26.07.18  we          undo 2.32 changes (some VPC special functions failed)
 2.39     16.08.18  we          remove assembler declaration from 32-bit _sincos
 2.40     17.08.18  we          signaling NaNs
 2.41     23.08.18  we          Bitmask names of single exceptions compatible to DAMath
 2.42     25.08.18  we          Total variance tvar/x
 2.43     21.09.18  we          Moved ext2 routines to new unit AMath2
 2.44     03.10.18  we          Suppress silly warnings for FPC32+ with stc.inc, remove old work-arounds
(2.45     22.09.18  we          ceil2x, floor2x, trunc2x, ldexp2x)
(2.46     23.09.18  we          abs2x, chs2x, cmp2x, exp1_2x, exp2x, ln2x, ln2_2x, nroot2x, pow2x)
 2.47     14.10.18  we          Reintegrated amath2/ext2

 2.48     05.11.18  we          Bern2n: array[0..MaxB2nSmall] absolute B2nHex
 2.49     06.11.18  we          pcov/x (population covariance)
 2.50     09.11.18  we          Generic aliases Infinity, NegInfinity, NaN

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
  AMath_Version = '2.50';

{#Z+}
{---------------------------------------------------------------------------}
{---------------------- Types, constants, variables ------------------------}
{---------------------------------------------------------------------------}
{#Z-}



type         {Functions with one real argument}
  TFuncX   = function(x: extended): extended;
  TFuncD   = function(x: double): double;

type         {Functions with one real argument and pointer to parameter(s)}
  TFuncXP  = function(x: extended; p: pointer): extended;
  TFuncDP  = function(x: double; p: pointer): double;

type
  THexExtA = packed array[0..9] of byte;  {Extended as array of bytes}
  THexDblA = packed array[0..7] of byte;  {Double   as array of bytes}
  THexSglA = packed array[0..3] of byte;  {Single   as array of bytes}
  THexExtW = packed array[0..4] of word;  {Extended as array of word}
  THexDblW = packed array[0..3] of word;  {Double   as array of word}
  THexSglW = packed array[0..1] of word;  {Single   as array of word}

type
  TExtRec  = packed record     {Extended as sign, exponent, significand}
               lm: longint;    {low  32 bit of significand}
               hm: longint;    {high 32 bit of significand}
               xp: word;       {biased exponent and sign  }
             end;

type
  TDblRec  = packed record     {Double as sign, exponent, significand}
               lm: longint;    {low  32 bit of significand}
               hm: longint;    {high bits of significand, biased exponent and sign}
             end;
type
  ext2 = record                {Double-extended as unevaluated sum of}
           l: extended;        {low  part and}
           h: extended;        {high part}
         end;


{#Z+}
{Machine epsilons: smallest (negative) powers of 2 with 1 + eps_? <> 1}
{These constants should evaluate to the hex values if viewed with ?,mh}
const eps_d: double   = 2.2204460492503131E-16;   {Hex: 000000000000B03C}
const eps_s: single   = 1.1920929E-7;             {Hex: 00000034}

const eps_xHex    : THexExtW = ($0000,$0000,$0000,$8000,$3fc0); {Machine epsilon for extended}
const PosInfXHex  : THexExtW = ($0000,$0000,$0000,$8000,$7fff); {extended +INF   as hex}
const NegInfXHex  : THexExtW = ($0000,$0000,$0000,$8000,$ffff); {extended -INF   as hex}
const NaNXHex     : THexExtW = ($ffff,$ffff,$ffff,$ffff,$7fff); {an extended quiet NaN as hex}
const sNaNXHex    : THexExtW = ($0001,$0000,$0000,$8000,$7fff); {an extended signaling NaN as hex}
const MinExtHex   : THexExtW = ($0000,$0000,$0000,$8000,$0001); {MinExtended as Hex}
const MaxExtHex   : THexExtW = ($ffff,$ffff,$ffff,$ffff,$7ffe); {MaxExtended as Hex}
const succx0Hex   : THexExtW = ($0001,$0000,$0000,$0000,$0000); {succx(0) as Hex   }
const Sqrt_MinXH  : THexExtW = ($0000,$0000,$0000,$8000,$2000); {sqrt(MinExtended) as Hex}
const Sqrt_MaxXH  : THexExtW = ($ffff,$ffff,$ffff,$ffff,$5ffe); {sqrt(MaxExtended) as Hex}
const ln_MaxXH    : THexExtW = ($79ab,$d1cf,$17f7,$b172,$400c); {predx(ln(MaxExtended))}
const ln_MinXH    : THexExtW = ($eb2f,$1210,$8c67,$b16c,$c00c); {succx(ln(MinExtended))}
const ln2hex      : THexExtW = ($79ac,$d1cf,$17f7,$b172,$3ffe); {ln(2)}
const log2ehex    : THexExtW = ($f0bc,$5c17,$3b29,$b8aa,$3fff); {log2(e)}
const log10ehex   : THexExtW = ($7195,$3728,$d8a9,$de5b,$3ffd); {log10(e)}
const TwoPihex    : THexExtW = ($c235,$2168,$daa2,$c90f,$4001); {2*Pi}
const Pi_2hex     : THexExtW = ($c235,$2168,$daa2,$c90f,$3fff); {Pi/2}
const Pi_4hex     : THexExtW = ($c235,$2168,$daa2,$c90f,$3ffe); {Pi/4}
const Pi_180hex   : THexExtW = ($C8AE,$94E9,$3512,$8EFA,$3FF9); {Pi/180}
const PiSqrHex    : THexExtW = ($F2D2,$F22E,$E64D,$9DE9,$4002); {9.8696044010893586185}
const SqrtPihex   : THexExtW = ($553D,$A77B,$C48D,$E2DF,$3FFF); {sqrt(Pi) = +1.7724538509055160273}
const Sqrt_2Pihex : THexExtW = ($2cb3,$b138,$98ff,$a06c,$4000); {sqrt(2*Pi)}
const LnPihex     : THexExtW = ($E85F,$3D0D,$8247,$9286,$3FFF); {ln(pi)}
const LnSqrt2Pihex: THexExtW = ($a535,$25f5,$8e43,$eb3f,$3ffe); {ln(sqrt(2*Pi)=}
const EulerGamHex : THexExtW = ($c7a5,$7db0,$67e3,$93c4,$3ffe); {Euler's constant}
const Sqrt2hex    : THexExtW = ($6484,$F9DE,$F333,$B504,$3FFF); {sqrt(2)}

const PosInfDHex  : THexDblW = ($0000,$0000,$0000,$7ff0); {double +INF  as hex}
const NegInfDHex  : THexDblW = ($0000,$0000,$0000,$fff0); {double -INF  as hex}
const NaNDHex     : THexDblW = ($ffff,$ffff,$ffff,$7fff); {a double quiet NaN as hex}
const sNaNDHex    : THexDblW = ($0001,$0000,$0000,$7ff0); {a double signaling NaN as hex}
const MaxDblHex   : THexDblW = ($ffff,$ffff,$ffff,$7fef); {1.797693134862315E+308}
const MinDblHex   : THexDblW = ($0000,$0000,$0000,$0010); {2.225073858507201E-308}

const PosInfSHex  : THexSglA = ($00,$00,$80,$7f);  {single +INF  as hex}
const NegInfSHex  : THexSglA = ($00,$00,$80,$ff);  {single -INF  as hex}
const NaNSHex     : THexSglA = ($ff,$ff,$ff,$7f);  {a single quiet NaN as hex}
const sNaNSHex    : THexSglA = ($01,$00,$80,$7f);  {a single signaling NaN as hex}
const MaxSglHex   : THexSglA = ($ff,$ff,$7f,$7f);  {3.4028234E+38}
const MinSglHex   : THexSglA = ($00,$00,$80,$00);  {1.1754944E-38}

{Small Bernoulli numbers as Hex constants. Formerly part }
{of unit sfBasic, they are now also used by unit AMCmplx.}

const
  MaxB2nSmall  = 60;

const
  B2nHex: array[0..MaxB2nSmall] of THexExtw = ( {Bernoulli(2n), n=0..}
            ($0000,$0000,$0000,$8000,$3FFF),  {+1.0000000000000000000}
            ($AAAB,$AAAA,$AAAA,$AAAA,$3FFC),  {+1.6666666666666666667E-1}
            ($8889,$8888,$8888,$8888,$BFFA),  {-3.3333333333333333335E-2}
            ($C30C,$0C30,$30C3,$C30C,$3FF9),  {+2.3809523809523809523E-2}
            ($8889,$8888,$8888,$8888,$BFFA),  {-3.3333333333333333335E-2}
            ($26CA,$6C9B,$C9B2,$9B26,$3FFB),  {+7.5757575757575757578E-2}
            ($8198,$9819,$1981,$8198,$BFFD),  {-2.5311355311355311355E-1}
            ($5555,$5555,$5555,$9555,$3FFF),  {+1.1666666666666666666}
            ($F2F3,$F2F2,$F2F2,$E2F2,$C001),  {-7.0921568627450980392}
            ($27C8,$9F1E,$7C78,$DBE2,$4004),  {+5.4971177944862155390E+1}
            ($67F4,$7F39,$F396,$8447,$C008),  {-5.2912424242424242426E+2}
            ($28D0,$33F1,$FC4A,$C180,$400B),  {+6.1921231884057971016E+3}
            ($6606,$0660,$2066,$A91A,$C00F),  {-8.6580253113553113550E+4}
            ($5555,$5555,$6955,$AE03,$4013),  {+1.4255171666666666666E+6}
            ($9C8B,$AE32,$DB88,$D044,$C017),  {-2.7298231067816091954E+7}
            ($FE36,$9A41,$9527,$8F6D,$401C),  {+6.0158087390064236836E+8}
            ($E5E6,$C5E5,$2B1D,$E140,$C020),  {-1.5116315767092156863E+10}
            ($5555,$EA55,$0E6E,$C80E,$4025),  {+4.2961464306116666666E+11}
            ($5309,$8E05,$E567,$C787,$C02A),  {-1.3711655205088332772E+13}
            ($9555,$CF4C,$5D33,$DE11,$402F),  {+4.8833231897359316666E+14}
            ($C84C,$1CEA,$45FA,$891C,$C035),  {-1.9296579341940068148E+16}
            ($9B70,$68B9,$B5E0,$BAE4,$403A),  {+8.4169304757368261500E+17}
            ($60ED,$6259,$67A8,$8BF3,$C040),  {-4.0338071854059455412E+19}
            ($D4E6,$F92A,$1ECC,$E551,$4045),  {+2.1150748638081991606E+21}
            ($8EB9,$1AF2,$630E,$CCC1,$C04B),  {-1.2086626522296525934E+23}
            ($C97F,$74B7,$D9A2,$C68B,$4051),  {+7.5008667460769643668E+24}
            ($50C0,$269C,$2369,$D066,$C057),  {-5.0387781014810689143E+26}
            ($D9ED,$4EC6,$CA9A,$EC0F,$405D),  {+3.6528776484818123334E+28}
            ($E7CF,$899B,$CBC7,$8FE1,$C064),  {-2.8498769302450882226E+30}
            ($5070,$15A2,$D8D6,$BC43,$406A),  {+2.3865427499683627645E+32}
            ($65EE,$E640,$2AAB,$83E3,$C071),  {-2.1399949257225333666E+34}
            ($707D,$1C11,$CE9D,$C56A,$4077),  {+2.0500975723478097570E+36}
            ($1E1C,$7AC8,$21EF,$9D85,$C07E),  {-2.0938005911346378408E+38}
            ($E927,$8376,$73E5,$85BA,$4085),  {+2.2752696488463515559E+40}
            ($8500,$E7B3,$9488,$F123,$C08B),  {-2.6257710286239576047E+42}
            ($B991,$385E,$7503,$E67C,$4092),  {+3.2125082102718032518E+44}
            ($0E52,$1EBF,$9A0F,$E92A,$C099),  {-4.1598278166794710914E+46}
            ($7C1E,$0021,$4E79,$F942,$40A0),  {+5.6920695482035280023E+48}
            ($970D,$F26E,$AF7E,$8C94,$C0A8),  {-8.2183629419784575694E+50}
            ($2C49,$524C,$2BE5,$A716,$40AF),  {+1.2502904327166993017E+53}
            ($C0A8,$EA20,$F1CB,$D0F8,$C0B6),  {-2.0015583233248370275E+55}
            ($8E7B,$4C4E,$5298,$8956,$40BE),  {+3.3674982915364374232E+57}
            ($7DC5,$B6DB,$4610,$BD7C,$C0C5),  {-5.9470970503135447718E+59}
            ($12CC,$59DB,$F874,$890D,$40CD),  {+1.1011910323627977560E+62}
            ($4590,$25D4,$A47F,$CFA5,$C0D4),  {-2.1355259545253501188E+64}
            ($A679,$4EF1,$AF84,$A492,$40DC),  {+4.3328896986641192418E+66}
            ($5E81,$3FBE,$35D8,$8854,$C0E4),  {-9.1885528241669328228E+68}
            ($8A90,$C146,$AAF0,$EBD8,$40EB),  {+2.0346896776329074493E+71}
            ($C12C,$4593,$68CD,$D4D3,$C0F3),  {-4.7003833958035731077E+73}
            ($FF4E,$F062,$48DC,$C82E,$40FB),  {+1.1318043445484249271E+76}
            ($3F3B,$AFB3,$532E,$C417,$C103),  {-2.8382249570693706958E+78}
            ($7E40,$DF17,$81CD,$C7E2,$410B),  {+7.4064248979678850632E+80}
            ($45B1,$D38C,$5DAE,$D3DC,$C113),  {-2.0096454802756604484E+83}
            ($5F41,$4EA9,$1C33,$E951,$411B),  {+5.6657170050805941445E+85}
            ($6A4B,$0BDB,$E3F8,$8563,$C124),  {-1.6584511154136216916E+88}
            ($4724,$DD5F,$F84A,$9E3F,$412C),  {+5.0368859950492377418E+90}
            ($328E,$4648,$E010,$C2A9,$C134),  {-1.5861468237658186369E+93}
            ($826C,$BA97,$AC47,$F81F,$413C),  {+5.1756743617545626984E+95}
            ($93F0,$781D,$43FD,$A3C1,$C145),  {-1.7488921840217117340E+98}
            ($0AFB,$0494,$C087,$DFB2,$414D),  {+6.1160519994952185254E+100}
            ($8384,$4095,$B028,$9E09,$C156)); {-2.2122776912707834942E+103}
{#Z-}

var
  {Absolute vars, i.e. constants with hex patterns}
  MinExtended : extended absolute MinExtHex;   {= 3.362103143112093507E-4932}  {= 2^(-16382)}
  MaxExtended : extended absolute MaxExtHex;   {= 1.189731495357231764E+4932}  {= 2^16384-2^16320 = (2^64-1)*2^16320}
  Sqrt_MinExt : extended absolute Sqrt_MinXH;  {= 1.833603867554847166E-2466}  {= 0.5^8191}
  Sqrt_MaxExt : extended absolute Sqrt_MaxXH;  {= 1.090748135619415929E+2466}  {= 2.0^8192}
  succx0      : extended absolute succx0Hex;   {= 0.364519953188247460E-4950}  {= succx(0) = 2^(-16445)}
  ln_MaxExt   : extended absolute ln_MaxXH;    {= 11356.52340629414394}
  ln_MinExt   : extended absolute ln_MinXH;    {=-11355.13711193302405}
  ln2         : extended absolute ln2hex;      {= 0.69314718055994530942}
  log2e       : extended absolute log2ehex;    {= 1.4426950408889634079 }
  log10e      : extended absolute log10ehex;   {= 0.43429448190325182765}
  TwoPi       : extended absolute TwoPihex;    {= 6.2831853071795864769 }
  Pi_2        : extended absolute Pi_2hex;     {= 1.5707963267948966192 }
  Pi_4        : extended absolute Pi_4hex;     {= 0.78539816339744830962}
  Pi_180      : extended absolute Pi_180hex;   {= 0.17453292519943295769e-1}
  PiSqr       : extended absolute PiSqrHex;    {= 9.8696044010893586185}
  SqrtPi      : extended absolute SqrtPihex;   {= 1.7724538509055160273}
  Sqrt_TwoPi  : extended absolute Sqrt_2Pihex; {= 2.5066282746310005024}
  LnPi        : extended absolute LnPihex;     {= 1.1447298858494001742}
  LnSqrt2Pi   : extended absolute LnSqrt2Pihex;{= 0.91893853320467274178}
  EulerGamma  : extended absolute EulerGamHex; {= 0.57721566490153286061}
  Sqrt2       : extended absolute Sqrt2hex;    {= 1.41421356237309504880168872421}
  eps_x       : extended absolute eps_xHex;    {1.084202172485504434E-19}

var
  PosInf_x    : extended absolute PosInfXHex;  {extended +INF  }
  NegInf_x    : extended absolute NegInfXHex;  {extended -INF  }
  NaN_x       : extended absolute NaNXHex;     {an extended quiet NaN}
  sNaN_x      : extended absolute sNaNXHex;    {an extended signaling  NaN}
  Infinity    : extended absolute PosInfXHex;  {extended +INF  }
  NaN	        : extended absolute NaNXHex;     {an extended quiet NaN}
  NegInfinity : extended absolute NegInfXHex;  {extended -INF  }

var
  MaxDouble   : double absolute MaxDblHex;     {1.797693134862315E+308} {= 2^1024 - 2^971}
  MinDouble   : double absolute MinDblHex;     {2.225073858507201E-308} {= 2^(-1022)}
  PosInf_d    : double absolute PosInfDHex;    {double +INF }
  NegInf_d    : double absolute NegInfDHex;    {double -INF }
  NaN_d       : double absolute NaNDHex;       {a double quiet NaN}
  sNaN_d      : double absolute sNaNDHex;      {a double signaling NaN}

var
  MaxSingle   : single absolute MaxSglHex;     {3.4028234E+38}
  MinSingle   : single absolute MinSglHex;     {1.1754944E-38}
  PosInf_s    : single absolute PosInfSHex;    {single +INF }
  NegInf_s    : single absolute NegInfSHex;    {single -INF }
  NaN_s       : single absolute NaNSHex;       {a single quiet NaN}
  sNaN_s      : single absolute sNaNSHex;      {a single signaling NaN}

var
  Bern2n      : array[0..MaxB2nSmall] of extended absolute B2nHex;

const
  sqrt_epsh   : extended = 2.3283064365386962890625e-10;  {sqrt(eps_x/2)}
  ln_succx0   = -11398.80538430830061;

const
  THREE       : extended = 3.0; {Used for circumventing some 'optimizations'. D2/D3}
                                {use *0.333.. instead of /3! This gives incorrectly}
                                {rounded results for 5/3, 7/3, and many others!    }
                                {Also used by FPC 3.1.1 with -O4 optimization!     }
  SIXX        : extended = 6.0; {Same reason (SIX would be sine integral)          }
  one_x       : extended = 1.0; {Avoid FPC nonsense using single for e.g. 1.0/n    }

const
  ph_cutoff   : extended = 549755813888.0;     {2^39, threshold for Payne/Hanek}
                                               {do not change or adjust tol in rem_pio2}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elementary transcendental functions -------------------}
{---------------------------------------------------------------------------}
{#Z-}

function arccos(x: extended): extended;
  {-Return the inverse circular cosine of x, |x| <= 1}

function arccos1m(x: extended): extended;
  {-Return arccos(1-x), 0 <= x <= 2, accurate even for x near 0}

function arccosd(x: extended): extended;
  {-Return the inverse circular cosine of x, |x| <= 1, result in degrees}

function arccosh(x: extended): extended;
  {-Return the inverse hyperbolic cosine, x >= 1. Note: for x near 1 the }
  { function arccosh1p(x-1) should be used to reduce cancellation errors!}

function arccosh1p(x: extended): extended;
  {-Return arccosh(1+x), x>=0, accurate even for x near 0}

function arccot(x: extended): extended;
  {-Return the sign symmetric inverse circular cotangent; arccot(x) = arctan(1/x), x <> 0}

function arccotc(x: extended): extended;
  {-Return the continuous inverse circular cotangent; arccotc(x) = Pi/2 - arctan(x)}

function arccotcd(x: extended): extended;
  {-Return the continuous inverse circular cotangent;}
  { arccotcd(x) = 90 - arctand(x), result in degrees }

function arccotd(x: extended): extended;
  {-Return the sign symmetric inverse circular cotangent,}
  { arccotd(x) = arctand(1/x), x <> 0, result in degrees }

function arccoth(x: extended): extended;
  {-Return the inverse hyperbolic cotangent of x, |x| > 1}

function arccsc(x: extended): extended;
  {-Return the inverse cosecant of x, |x| >= 1}

function arccsch(x: extended): extended;
  {-Return the inverse hyperbolic cosecant of x, x <> 0}

function arcgd(x: extended): extended;
  {-Return the inverse Gudermannian function arcgd(x), |x| < Pi/2}

function archav(x: extended): extended;
  {-Return the inverse haversine archav(x), 0 <= x <= 1}

function arcsec(x: extended): extended;
  {-Return the inverse secant of x, |x| >= 1}

function arcsech(x: extended): extended;
  {-Return the inverse hyperbolic secant of x, 0 < x <= 1}

function arcsin(x: extended): extended;
  {-Return the inverse circular sine of x, |x| <= 1}

function arcsind(x: extended): extended;
  {-Return the inverse circular sine of x, |x| <= 1, result in degrees}

function arcsinh(x: extended): extended;
  {-Return the inverse hyperbolic sine of x}

function arctan2(y, x: extended): extended;
  {-Return arctan(y/x); result in [-Pi..Pi] with correct quadrant}

function arctand(x: extended): extended;
  {-Return the inverse circular tangent of x, result in degrees}

function arctanh(x: extended): extended;
  {-Return the inverse hyperbolic tangent of x, |x| < 1}

function compound(x: extended; n: longint): extended;
  {-Return (1+x)^n; accurate version of Delphi/VP internal function}

function comprel(x,n: extended): extended;
  {-Return ((1+x)^n-1)/x; accurate version of Delphi/VP internal function}

function cos(x: extended): extended;
  {-Accurate version of circular cosine, uses system.cos for |x| <= Pi/4}

function cosd(x: extended): extended;
  {-Return cos(x), x in degrees}

function cosh(x: extended): extended;
  {-Return the hyperbolic cosine of x}

function coshm1(x: extended): extended;
  {-Return cosh(x)-1, accurate even for x near 0}

function cosPi(x: extended): extended;
  {-Return cos(Pi*x), result will be 1 for abs(x) >= 2^64}

function cot(x: extended): extended;
  {-Return the circular cotangent of x, x mod Pi <> 0}

function cotd(x: extended): extended;
  {-Return cot(x), x in degrees}

function coth(x: extended): extended;
  {-Return the hyperbolic cotangent of x, x<>0}

function covers(x: extended): extended;
  {-Return the coversine covers(x) = 1 - sin(x)}

function csc(x: extended): extended;
  {-Return the circular cosecant of x, x mod Pi <> 0}

function csch(x: extended): extended;
  {-Return the hyperbolic cosecant of x, x<>0}

function exp(x: extended): extended;
  {-Accurate exp, result good to extended precision}

function exp10(x: extended): extended;
  {-Return 10^x}

function exp10m1(x: extended): extended;
  {-Return 10^x - 1; special code for small x}

function exp2(x: extended): extended;
  {-Return 2^x}

function exp2m1(x: extended): extended;
  {-Return 2^x-1, accurate even for x near 0}

function exp3(x: extended): extended;
  {-Return 3^x}

function exp5(x: extended): extended;
  {-Return 5^x}

function exp7(x: extended): extended;
  {-Return 7^x}

function expm1(x: extended): extended;
  {-Return exp(x)-1, accurate even for x near 0}

function exprel(x: extended): extended;
  {-Return exprel(x) = (exp(x) - 1)/x,  1 for x=0}

function expx2(x: extended): extended;
  {-Return exp(x*|x|) with damped error amplification in computing exp of the product.}
  { Used for exp(x^2) = expx2(abs(x)) and exp(-x^2) = expx2(-abs(x))}

function expmx2h(x: extended): extended;
  {-Return exp(-0.5*x^2) with damped error amplification}

function gd(x: extended): extended;
  {-Return the Gudermannian function gd(x)}

function hav(x: extended): extended;
  {-Return the haversine hav(x) = 0.5*(1 - cos(x))}

function ln(x: extended): extended;
  {-Return natural logarithm of x, x may be denormal}

function ln1mexp(x: extended): extended;
  {-Return ln(1-exp(x)), x<0}

function ln1pexp(x: extended): extended;
  {-Accurately compute ln(1+exp(x)) without overflow}

function ln1p(x: extended): extended;
  {-Return ln(1+x), accurate even for x near 0}

function ln1pmx(x: extended): extended;
  {-Return ln(1+x)-x, x>-1, accurate even for -0.5 <= x <= 1.0}

function lncosh(x: extended): extended;
  {-Return ln(cosh(x)), accurate for x ~ 0 and without overflow for large x}

function lnsinh(x: extended): extended;
  {-Return ln(sinh(x)), x > 0, accurate for x ~ 0 and without overflow for large x}

function log10(x: extended): extended;
  {-Return base 10 logarithm of x}

function log10p1(x: extended): extended;
  {-Return log10(1+x), accurate even for x near 0}

function log2(x: extended): extended;
  {-Return base 2 logarithm of x}

function log2p1(x: extended): extended;
  {-Return log2(1+x), accurate even for x near 0}

function logbase(b, x: extended): extended;
  {-Return base b logarithm of x}

function logaddexp(x,y: extended): extended;
  {-Accurately compute ln[exp(x) + exp(y)]}

function logsubexp(x,y: extended): extended;
  {-Accurately compute ln[exp(x) - exp(y)], x > y}

function logistic(x: extended): extended;
  {-Return logistic(x) = 1/(1+exp(-x))}

function logit(x: extended): extended;
  {-Return logit(x) = ln(x/(1.0-x)), accurate near x=0.5}

function power(x, y : extended): extended;
  {-Return x^y; if frac(y)<>0 then x must be > 0}

function powm1(x,y: extended): extended;
  {-Return x^y - 1; special code for small x,y}

function pow1p(x,y: extended): extended;
  {-Return (1+x)^y, x > -1, with ext2 arithmetic for critical values}

function pow1pf(x,y: extended): extended;
  {-Return (1+x)^y, x > -1, without ext2, less accurate than pow1p}

function pow1pm1(x,y: extended): extended;
  {-Return (1+x)^y - 1; special code for small x,y}

function powpi2k(k,n: longint): extended;
  {-Return accurate scaled powers of Pi, result = (Pi*2^k)^n}

function powpi(n: longint): extended;
  {-Return accurate powers of Pi, result = Pi^n}

function sec(x: extended): extended;
  {-Return the circular secant of x, x mod Pi <> Pi/2}

function sech(x: extended): extended;
  {-Return the hyperbolic secant of x}

function sin(x: extended): extended;
  {-Accurate version of circular sine, uses system.sin for |x| <= Pi/4}

procedure sincos(x: extended; var s,c: extended);
  {-Return accurate values s=sin(x), c=cos(x)}

procedure sincosd(x: extended; var s,c: extended);
  {-Return sin(x) and cos(x), x in degrees}

procedure sincosPi(x: extended; var s,c: extended);
  {-Return s=sin(Pi*x), c=cos(Pi*x); (s,c)=(0,1) for abs(x) >= 2^64}

procedure sinhcosh(x: extended; var s,c: extended);
  {-Return s=sinh(x) and c=cosh(x)}

function sinc(x: extended): extended;
  {-Return the cardinal sine sinc(x) = sin(x)/x}

function sincPi(x: extended): extended;
  {-Return the normalised cardinal sine sincPi(x) = sin(Pi*x)/(Pi*x)}

function sind(x: extended): extended;
  {-Return sin(x), x in degrees}

function sinh(x: extended): extended;
  {-Return the hyperbolic sine of x, accurate even for x near 0}

function sinhc(x: extended): extended;
  {-Return sinh(x)/x, accurate even for x near 0}

function sinhmx(x: extended): extended;
  {-Return sinh(x)-x, accurate even for x near 0}

function sinPi(x: extended): extended;
  {-Return sin(Pi*x), result will be 0 for abs(x) >= 2^64}

function tan(x: extended): extended;
  {-Return the circular tangent of x, x mod Pi <> Pi/2}

function tand(x: extended): extended;
  {-Return tan(x), x in degrees}

function tanh(x: extended): extended;
  {-Return the hyperbolic tangent of x, accurate even for x near 0}

function tanPi(x: extended): extended;
  {-Return tan(Pi*x), result will be 0 for abs(x) >= 2^64}

function vers(x: extended): extended;
  {-Return the versine vers(x) = 1 - cos(x)}

function versint(x: extended): extended;
  {-Return versint(x) = integral(vers(t),t=0..x) = x - sin(x), accurate near 0}


function logN(N, x: extended): extended; {-Delphi alias for logbase}
function lnxp1(x: extended): extended;   {-Delphi alias for ln1p}


{#Z+}
{---------------------------------------------------------------------------}
{---------------------- Elementary numerical functions ---------------------}
{---------------------------------------------------------------------------}
{#Z-}
function cbrt(x: extended): extended;
  {-Return the cube root of x}

function ceil(x: extended): longint;
  {-Return the smallest integer >= x; |x|<=MaxLongint}

function ceilx(x: extended): extended;
  {-Return the smallest integer >= x}

function floor(x: extended): longint;
  {-Return the largest integer <= x; |x|<=MaxLongint}

function floorx(x: extended): extended;
  {-Return the largest integer <= x}

function fmod(x,y: extended): extended;
  {-Return x mod y, y<>0, sign(result) = sign(x)}

function hypot(x,y: extended): extended;
  {-Return sqrt(x*x + y*y)}

function hypot3(x,y,z: extended): extended;
  {-Return sqrt(x*x + y*y + z*z)}

function intpower(x: extended; n: longint): extended;
  {-Return x^n; via binary exponentiation (no overflow detection)}

function modf(x: extended; var ip: longint): extended;
  {-Return frac(x) and trunc(x) in ip, |x|<=MaxLongint}

function nroot(x: extended; n: integer): extended;
  {-Return the nth root of x; n<>0, x >= 0 if n is even}

function remainder(x,y: extended): extended;
  {-Return the IEEE754 remainder x REM y = x - rmNearest(x/y)*y}

function sqrt1pm1(x: extended): extended;
  {-Return sqrt(1+x)-1, accurate even for x near 0, x>=-1}

function sqrt1pmx(x: extended): extended;
  {-Return sqrt(1+x^2)-x}

{#Z+}
{---------------------------------------------------------------------------}
{----------------------- Floating point functions --------------------------}
{---------------------------------------------------------------------------}
{#Z-}
function  copysign(x,y: extended): extended; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return abs(x)*sign(y)}

function  copysignd(x,y: double): double;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return abs(x)*sign(y)}

function  copysigns(x,y: single): single;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return abs(x)*sign(y)}

procedure frexp(x: extended; var m: extended; var e: longint);
  {-Return the mantissa m and exponent e  of x with x = m*2^e, 0.5 < m < 1;}
  { if x is 0, +-INF, NaN, return m=x, e=0}

procedure frexpd(d: double; var m: double; var e: longint);
  {-Return the mantissa m and exponent e  of d with d = m*2^e, 0.5 < m < 1;}
  { if d is 0, +-INF, NaN, return m=d, e=0}

procedure frexps(s: single; var m: single; var e: longint);
  {-Return the mantissa m and exponent e of s with s = m*2^e, 0.5 <= abs(m) < 1;}
  { if s is 0, +-INF, NaN, return m=s, e=0}

function  ilogb(x: extended): longint;
  {-Return base 2 exponent of x. For finite x ilogb = floor(log2(|x|)), }
  { otherwise -MaxLongint for x = 0 or MaxLongint if x = +-INF or Nan.  }

function  ldexp(x: extended; e: longint): extended;
  {-Return x*2^e}

function  ldexpd(d: double; e: longint): double;
  {-Return d*2^e}

function  ldexps(s: single; e: longint): single;
  {-Return s*2^e}

function  IsInf(x: extended): boolean;          {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x is +INF or -INF}

function  IsInfD(d: double): boolean;           {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is +INF or -INF}

function  IsInfS(s: single): boolean;           {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is +INF or -INF}

function  IsNaN(x: extended): boolean;          {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x is a NaN}

function  IsNaND(d: double): boolean;           {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is a NaN}

function  IsNaNS(s: single): boolean;           {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is a NaN}

function  IsNaNorInf(x: extended): boolean;     {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x is a NaN or infinite}

function  IsNaNorInfD(d: double): boolean;      {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is a NaN or infinite}

function  IsNaNorInfS(s: single): boolean;      {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is a NaN or infinite}

function isRMNearest: boolean;
  {-Check if round_to_nearest without FPU control}

function  predd(d: double): double;
  {-Return next representable double after d in the direction -Inf}

function  predx(x: extended): extended;
  {-Return next representable extended after x in the direction -Inf}

function  preds(s: single): single;
  {-Return next representable single after s in the direction -Inf}

function  rint(x: extended): extended;
  {-Return the integral value nearest x for the current rounding mode}

function  succd(d: double): double;
  {-Return next representable double after d in the direction +Inf}

function  succx(x: extended): extended;
  {-Return next representable extended after x in the direction +Inf}

function  succs(s: single): single;
  {-Return next representable single after s in the direction +Inf}

function  scalbn(x: extended; e: longint): extended;
  {-Return x*2^e}

function  ulpd(d: double): double;
  {-Return the 'unit in the last place': ulpd(d)=|d|-predd(|d|) for finite d}

function  ulps(s: single): single;
  {-Return the 'unit in the last place': ulps(s)=|s|-preds(|s|) for finite s}

function  ulpx(x: extended): extended;
  {-Return the 'unit in the last place': ulpx(x)=|x|-predx(|x|) for finite x}

{#Z+}
{---------------------------------------------------------------------------}
{------------------------- FPU control functions ---------------------------}
{---------------------------------------------------------------------------}

{These functions try to be as compatible to Delphi/FreePascal as reasonable.}
{16 bit Pascal cannot return sets as function results, Delphi sets are not  }
{byte-compatible in FPC, not all Delphi versions support default values etc.}
{D3+ and FPC2+ have Set8087CW in system.pas. It is used here if available,  }
{because it also sets the system variable Default8087CW. D6+ and FPC+ have  }
{Get8087CW in system.pas. The following series of ifdefs finds out which of }
{the functions have to be defined here.}
{#Z-}

{$ifdef BIT16}
  {$define need_get87}
  {$define need_set87}
{$else}
  {$ifdef VirtualPascal}
    {$define need_get87}
    {$define need_set87}
    {$define buggy_round}
  {$endif}

  {$ifdef FPC}
    {$ifdef VER1}
      {$define need_get87}
      {$define need_set87}
      {$define buggy_trunc}
      {$define buggy_round}
    {$endif}
  {$endif}

  {$ifdef Delphi}
    {$ifndef CONDITIONALEXPRESSIONS}
      {$ifdef VER90 }
        {$define need_set87}
      {$endif}
      {Delphi 3+ has Set8087CW}
      {Delphi 6+ has Get8087CW}
      {$define need_get87}
      {$define buggy_trunc}
    {$endif}
  {$endif}
{$endif}

type
  TFPURoundingMode  = (rmNearest, rmDown, rmUp, rmTruncate);

type
  TFPUPrecisionMode = (pmSingle, pmReserved, pmDouble, pmExtended);

const
  {Bitmasks of single exceptions}
  bexInvalidOp    = $01;
  bexDenormalized = $02;
  bexZeroDivide   = $04;
  bexOverflow     = $08;
  bexUnderflow    = $10;
  bexPrecision    = $20;

  bexAllArithmeticExceptions = $3F;


{$ifdef need_get87}
function Get8087CW: word;
  {-Return the FPU control word}
{$endif}

{$ifdef need_set87}
procedure Set8087CW(cw: word);
  {-Set new FPU control word}
{$endif}

function  GetRoundMode: TFPURoundingMode;
  {-Return the current rounding mode}

function  SetRoundMode(NewRoundMode: TFPURoundingMode): TFPURoundingMode;
  {-Set new rounding mode and return the old mode}

function  GetPrecisionMode: TFPUPrecisionMode;
  {-Return the current precision control mode}

function  SetPrecisionMode(NewPrecision: TFPUPrecisionMode): TFPUPrecisionMode;
  {-Set new precision control mode and return the old precision}

procedure GetExceptionMask(var Mask: byte);
  {-Return the current exception mask}

procedure SetExceptionMask(NewMask: byte);
  {-Set new exception mask}


{#Z+}
{---------------------------------------------------------------------------}
{-------------- Polynomial, Vector, Statistic Operations -------------------}
{---------------------------------------------------------------------------}
{#Z-}

function  CoeffVar(const a: array of double; n: integer): extended;
  {-Return the coefficient of variation (sdev/mean) of a double vector}

function  CoeffVarX(const a: array of extended; n: integer): extended;
  {-Return the coefficient of variation (sdev/mean) of an extended vector}

function  CSEvalX(x: extended; const a: array of extended; n: integer): extended;
  {-Evaluate Chebyshev sum a[0]/2 + a[1]*T_1(x) +..+ a[n-1]*T_(n-1)(x) using Clenshaw algorithm}

procedure CSEvalXDer(x: extended; const a: array of extended; n: integer; var px,dp: extended);
  {-Evaluate Chebyshev sum p(x) = a[0]/2 + a[1]*T_1(x) +..+ a[n-1]*T_(n-1)(x),}
  { using Clenshaw algorithm, return pd = p(x) and dp = p'(x) }

function  dot2(const x,y: array of double; n: integer): extended;
  {-Accurate dot product sum(x[i]*y[i], i=0..n-1) of two double vectors}

function  dot2x(const x,y: array of extended; n: integer): extended;
  {-Accurate dot product sum(x[i]*y[i], i=0..n-1) of two extended vectors}

function  mean(const a: array of double; n: integer): extended;
  {-Compute accurate mean = sum(a[i], i=0..n-1)/n of a double vector}

function  meanx(const a: array of extended; n: integer): extended;
  {-Compute accurate mean = sum(a[i], i=0..n-1)/n of an extended vector}

procedure MeanAndStdDev(const a: array of double; n: integer; var mval, sdev: extended);
  {-Accurate mean and sample standard deviation of a double vector}

procedure MeanAndStdDevX(const a: array of extended; n: integer; var mval, sdev: extended);
  {-Accurate mean and sample standard deviation of an extended vector}

procedure moment(const a: array of double; n: integer; var m1, m2, m3, m4, skew, kurt: extended);
  {-Return the first 4 moments, skewness, and kurtosis of a double vector}

procedure momentx(const a: array of extended; n: integer; var m1, m2, m3, m4, skew, kurt: extended);
  {-Return the first 4 moments, skewness, and kurtosis of an extended vector}

procedure mssqd(const a: array of double; n: integer; var mval, scale, sumsq: extended);
  {-Calculate mean mval and ssqd sum((a[i]-mval)^2) of a double vector}

procedure mssqx(const a: array of extended; n: integer; var mval, scale, sumsq: extended);
  {-Calculate mean mval and ssqx sum((a[i]-mval)^2) of an extended vector}

function  norm2(const a: array of double; n: integer): extended;
  {-Calculate the 2-norm = sqrt(sum(a[i]^2, i=0..n-1)) of a double vector}

function  norm2x(const a: array of extended; n: integer): extended;
  {-Calculate the 2-norm = sqrt(sum(a[i]^2, i=0..n-1)) of an extended vector}

function  PolEval(x: extended; const a: array of double; n: integer): extended;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}

procedure PolEvalC(const a: array of double; n: integer; x,y: double; var u,v: double);
  {-Evaluate polynomial a[0] + a[1]*z + ... + a[n-1]*z^(n-1)}
  { for complex z = x + i*y, result is u + i*v}

procedure PolEvalDeriv(x: double; const a: array of double; n: integer; var px,dp: extended);
  {-Evaluate polynomial p(x) a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}
  { Return px = p(x) and dp = p'(x)}

function  PolEvalS(x: extended; const a: array of single; n: integer): extended;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}

function  PolEvalX(x: extended; const a: array of extended; n: integer): extended;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}

function  PolEvalEE(x: double; const a: array of double; n: integer; var e: double): double;
  {-Evaluate polynomial; return p(x) = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.   }

function  PolEvalEEX(x: extended; const a: array of extended; n: integer; var e: extended): extended;
  {-Evaluate polynomial; return p(x) = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.   }

function  PolEvalCHE(x: double; const a: array of double; n: integer; var e: double): double;
  {-Evaluate polynomial; return p(x) = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { accurate double precision version using compensated Horner scheme,    }
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.   }

procedure PolEvalCHEDer(x: double; const a: array of double; n: integer; var px,dp,e: double);
  {-Evaluate polynomial; return p = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { accurate double precision version using compensated Horner scheme, }
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.}
  { dp is the non-compensated derivative p'(x).                        }

function  PolEvalCHEX(x: extended; const a: array of extended; n: integer; var e: extended): extended;
  {-Evaluate polynomial; return p(x) = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { accurate extended precision version using compensated Horner scheme,  }
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.   }

function  rms(const a: array of double; n: integer): extended;
  {-Calculate the RMS value sqrt(sum(a[i]^2, i=0..n-1)/n) of a double vector}

function  rmsx(const a: array of extended; n: integer): extended;
  {-Calculate the RMS value sqrt(sum(a[i]^2, i=0..n-1)/n) of an extended vector}

procedure ssqx(const a: array of extended; n: integer; var scale, sumsq: extended);
  {-Calculate sum(a[i]^2, i=0..n-1) = scale^2*sumsq, scale>=0, sumsq>0}

procedure ssqd(const a: array of double; n: integer; var scale, sumsq: extended);
  {-Calculate sum(a[i]^2, i=0..n-1) = scale^2*sumsq, scale>=0, sumsq>0}

function  sum2(const a: array of double; n: integer): extended;
  {-Compute accurate sum(a[i], i=0..n-1) of a double vector}

function  sum2x(const a: array of extended; n: integer): extended;
  {-Compute accurate sum(a[i], i=0..n-1) of extended vector}

function  sumsqr(const a: array of double; n: integer): extended;
  {-Calculate sum(a[i]^2, i=0..n-1) of a double vector}

function  sumsqrx(const a: array of extended; n: integer): extended;
  {-Calculate sum(a[i]^2, i=0..n-1) of an extended vector}

function  ssdev(const a: array of double; n: integer): extended;
  {-Return the sample standard deviation of a double vector}

function  ssdevx(const a: array of extended; n: integer): extended;
  {-Return the sample standard deviation of an extended vector}

function  psdev(const a: array of double; n: integer): extended;
  {-Return the population standard deviation of a double vector}

function  psdevx(const a: array of extended; n: integer): extended;
  {-Return the population standard deviation of an extended vector}

function  svar(const a: array of double; n: integer): extended;
  {-Return the sample variance of a double vector}

function  svarx(const a: array of extended; n: integer): extended;
  {-Return the sample variance of an extended vector}

function  pvar(const a: array of double; n: integer): extended;
  {-Return the population variance of a double vector}

function  pvarx(const a: array of extended; n: integer): extended;
  {-Return the population variance of an extended vector}

function  pcov(const x,y: array of double; n: integer): extended;
  {-Calculate the population covariance = sum((x[i]-mx)(y[i]-my), i=0..n-1)/n of two double vectors}

function  pcovx(const x,y: array of extended; n: integer): extended;
  {-Calculate the population covariance = sum((x[i]-mx)(y[i]-my), i=0..n-1)/n of two extended vectors}

function  tvar(const a: array of double; n: integer): extended;
  {-Return the total variance of a double vector}

function  tvarx(const a: array of extended; n: integer): extended;
  {-Return the total variance of an extended vector}

{#Z+}
{--------------------------------------------------------------------}
{---------- Argument reduction for trigonometric functions ----------}
{--------------------------------------------------------------------}
{#Z-}
function rem_pio2_cw(x: extended; var z: extended): integer;
  {-Cody/Waite reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8}

function rem_pio2_ph(x: extended; var z: extended): integer;
  {-Payne/Hanek reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8}

function rem_pio2(x: extended; var z: extended): integer;
  {-Argument reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8.}
  { Uses Payne/Hanek if |x| > ph_cutoff, Cody/Waite otherwise}

function rem_2pi(x: extended): extended;
  {-Return x mod 2*Pi}

function rem_2pi_sym(x: extended): extended;
  {-Return x mod 2*Pi, -Pi <= result <= Pi}

function rem_int2(x: extended; var z: extended): integer;
  {-Argument reduction of x: z*Pi = x*Pi - n*Pi/2, |z|<=1/4, result = n mod 8.}
  { Used for argument reduction in sin(Pi*x) and cos(Pi*x)}

{#Z+}
{--------------------------------------------------------------------}
{------------------------- Other function ---------------------------}
{--------------------------------------------------------------------}
{#Z-}

function DegToRad(x: extended): extended;
  {-Convert angle x from degrees to radians}

function RadToDeg(x: extended): extended;
  {-Convert angle x from radians to degrees}

function angle2(x1,x2,y1,y2: extended): extended;
  {-Return the accurate angle between the vectors (x1,x2) and (y1,y2)}

function area_triangle(x1,y1,x2,y2,x3,y3: extended): extended;
  {-Return the area of the triangle defined by the points (xi,yi)}

function Dbl2Hex(d: double): string;
  {-Return d as a big-endian hex string}

function Ext2Dbl(x: extended): double;
  {-Return x as double, or +-Inf if too large}

function Ext2Hex(x: extended): string;
  {-Return x as a big-endian hex string}

function Sgl2Hex(s: single): string;
  {-Return s as a big-endian hex string}

function fisEQd(x,y: double): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are bit-identical}

function fisEQx(x,y: extended): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are bit-identical}

function fisNEd(x,y: double): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are not bit-identical}

function fisNEx(x,y: extended): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are not bit-identical}

function fisEQs(x,y: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are bit-identical}

function fisNEs(x,y: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are not bit-identical}

procedure Hex2Ext(const hex: string; var x: extended; var code: integer);
  {-Convert big-endian hex string to extended, leading $ is skipped, OK if code=0;}
  { hex must have 20 hex characters (21 if leading $), inverse of Ext2Hex.}

procedure Hex2Dbl(const hex: string; var d: double; var code: integer);
  {-Convert big-endian hex string to double, leading $ is skipped, OK if code=0;}
  { hex must have 16 hex characters (17 if leading $), inverse of Dbl2Hex.}

procedure Hex2Sgl(const hex: string; var s: single; var code: integer);
  {-Convert big-endian hex string to single, leading $ is skipped, OK if code=0;}
  { hex must have 8 hex characters (9 if leading $), inverse of Sgl2Hex.}

function in_triangle(x,y,x1,y1,x2,y2,x3,y3: extended): boolean;
  {-Return true if the point (x,y) lies strictly inside the triangle defined}
  { by the three points (xi,yi), false if it lies on a side or outside.}

function in_triangle_ex(x,y,x1,y1,x2,y2,x3,y3: extended): integer;
  {-Return +1 if the point (x,y) lies strictly inside the triangle defined by}
  { the three points (xi,yi), -1 if it lies strictly outside, 0 otherwise.}

function isign(x: extended): integer;
  {-Return the sign of x, 0 if x=0 or NAN}

function maxd(x, y: double): double;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the maximum of two doubles; x,y <> NAN}

function maxx(x, y: extended): extended; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the maximum of two extendeds; x,y <> NAN}

function maxs(x, y: single): single; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the maximum of two singles; x,y <> NAN}

function mind(x, y: double): double; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the minimum of two doubles; x,y <> NAN}

function minx(x, y: extended): extended; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the minimum of two extendeds; x,y <> NAN}

function mins(x, y: single): single; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the minimum of two singles; x,y <> NAN}

function orient2d(x1,y1,x2,y2,x3,y3: extended): extended;
  {-Return the mathematical orientation of the three points (xi,yi): >0 if}
  { the order is counterclockwise, <0 if clockwise, =0 if they are collinear.}
  { Result is twice the signed area of the triangle defined by the points.}

function RandG(Mean, StdDev: extended): extended;
  {-Random number from Gaussian (normal) distribution with given mean}
  { and standard deviation |StdDev|}

function RandG01: extended;
  {-Random number from standard normal distribution (Mean=0, StdDev=1)}

function SafeDivX(x,y: extended): extended;
  {-Safe quotient x/y, +-Infinity if overflow, NaN if x=y=0}

{#Z+}
{--------------------------------------------------------------------}
{----------------- ext2 (double-extended) functions  ----------------}
{--------------------------------------------------------------------}
{#Z-}

procedure hhto2x(const a,b: THexExtW; var x: ext2);
  {-Return x = a + b using xxto2x}

procedure xxto2x(a,b: extended; var x: ext2);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a + b using TwoSum algorithm}

procedure xto2x(a: extended; var x: ext2);     {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a}

procedure add2x(const a,b: ext2; var x: ext2); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a+b}

procedure add21x(const a: ext2; b: extended; var x: ext2);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a+b}

procedure sub2x(const a,b: ext2; var x: ext2); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a-b}

procedure mul2x(const a,b: ext2; var x: ext2);
  {-Return x = a*b}

procedure mul21x(const a: ext2; b: extended; var x: ext2);
  {-Return x = a*b}

procedure sqr2x(const a: ext2; var x: ext2);
  {-Return x = a^2}

procedure pow2xi(const a: ext2; n: extended; var y: ext2);
  {-Return y = a^n, frac(n)=0, a<>0 if n<0}

procedure div2x(const a,b: ext2; var x: ext2);
  {-Return x = a/b,  b<>0}

procedure div21x(const a: ext2; b: extended; var x: ext2);
  {-Return x = a/b,  b<>0}

procedure inv2x(const b: ext2; var x: ext2);
  {-Return x = 1/b,  b<>0}

procedure sqrt2x(const a: ext2; var x: ext2);
  {-Return x = sqrt(a),  a >= 0}

procedure pi2x(var x: ext2);
  {-Return x = Pi with double-extended precision}

procedure xadd12(const a, b: extended; var xh, xl: extended); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return [xh,xl] = a + b}

procedure xmul12(a,b: extended; var xh,xl: extended); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return [xh,xl] = a * b}

procedure sqr12x(a: extended; var xh,xl: extended); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return [xh,xl] = a^2}

procedure xdivrem(a,b: extended; var q,r: extended);
  {-Compute q,r with a = q*b + r,  assumes round to nearest}

function  fma_x(a,b,c: extended): extended;
  {-Accurately compute a*b + c}

procedure floor2x(const a: ext2; var x: ext2);
  {-Return x = floor(a)}

procedure ceil2x(const a: ext2; var x: ext2);
  {-Return x = ceil(a)}

procedure trunc2x(const a: ext2; var x: ext2);
  {-Return x = trunc(a)}

procedure ldexp2x(const a: ext2; n: longint; var x: ext2);
  {-Return x = a*2^n}

procedure abs2x(const a: ext2; var x: ext2);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = abs(a)}

procedure chs2x(const a: ext2; var x: ext2);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = -a}

function  cmp2x(const a,b: ext2): integer;
  {-Return sign(a-b)}

procedure ln2_2x(var x: ext2);
  {-Return x = ln(2) with double-extended precision}

procedure exp1_2x(var x: ext2);
  {-Return x = exp(1) with double-extended precision}

procedure exp2x(const a: ext2; var x: ext2);
  {-Return x = exp(a)}

procedure ln2x(const a: ext2; var x: ext2);
  {-Return x = ln(a)}

procedure pow2x(const a,b: ext2; var x: ext2);
  {-Return x = a^b, a > 0}

procedure nroot2x(const a: ext2; n: longint; var x: ext2);
  {-Return x = a^(1/n), a > 0 if n is even}



{#Z+}
{---------------------------------------------------------------------------}
{------------------  Internal and bugfix functions -------------------------}
{---------------------------------------------------------------------------}
{#Z-}
{$ifdef BIT32}
{Use extended argument (not CONST extended), i.e. compatible to type TFuncX}
function int(x: extended): extended;
  {-Return the integer part x; the value of x rounded toward zero}

function frac(x: extended): extended;
  {-Return the fractional part x = x - int(x)}
{$endif}

function arctan(x: extended): extended;
  {-Return the inverse circular tangent of x}

function sqrt(x: extended): extended;
  {-Return the square root of x >= 0}

{$ifdef buggy_trunc}
{Bugfix version for Delphi2..5, FPC1}
function trunc(x: extended): longint;
{$endif}

{$ifdef buggy_round}
function round(x: extended): longint;
{$endif}

function intpower0(x: extended; n: longint): extended;
  {-Return x^n; via binary exponentiation (no overflow detection)}


implementation



const
  x80000000 = longint($80000000);


{FPC and Delphi 2 do not like fstp [@result]}
{$ifdef FPC}
  {$define NO_FPU_RESULT}
{$else}
  {$ifdef VER90}
    {$define NO_FPU_RESULT}
  {$endif}
{$endif}


{Bugfix for Delphi2-5: These versions simply do a fldcw cwChop (= $1F32)}
{instead of or'ing the current control word with $0F00. Sometimes there }
{will be crashes after fldcw if operations with NAN's or INF's have been}
{performed before these functions are called. FPC1 has similar bugs.}

{For other 32-bit compilers make frac and int compatible to type TFuncX}


{$ifdef BIT32}
{---------------------------------------------------------------------------}
function int(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return the integer part x; the value of x rounded toward zero}
asm
  fld     [x]
  sub     esp,4
  fnstcw  word ptr [esp]
  fnstcw  word ptr [esp+2]
  fwait
  or      word ptr [esp+2],$0f00  {set extended precision, round toward zero}
  fldcw   word ptr [esp+2]
  frndint
  fwait
  fldcw   word ptr [esp]
  add     esp,4
end;

{---------------------------------------------------------------------------}
function frac(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return the fractional part x = x - int(x)}
asm
  fld     [x]
  fld     st(0)
  sub     esp,4
  fnstcw  word ptr [esp]
  fnstcw  word ptr [esp+2]
  fwait
  or      word ptr [esp+2],$0f00  {set extended precision, round toward zero}
  fldcw   word ptr [esp+2]
  frndint
  fwait
  fldcw   word ptr [esp]
  add     esp,4
 {$ifdef FPC}
  fsubp   st(1),st            {!!!??? brain-damage: fsub gives two warnings!}
 {$else}
  fsub
 {$endif}
end;
{$endif}

{---------------------------------------------------------------------------}
function arctan(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return the inverse circular tangent of x}
asm
  fld    [x]
  fld1
  fpatan
  fwait
end;

{---------------------------------------------------------------------------}
function sqrt(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return the square root of x >= 0}
asm
  fld     [x]
  fsqrt
  fwait
end;


{$ifdef buggy_trunc}
{---------------------------------------------------------------------------}
function trunc(x: extended): longint;
var
  cws,cw: word;
  iarr: packed array[0..1] of longint;
begin
  asm
    fld     [x]
    fstcw   [cws]
    fstcw   [cw]
    fwait
    or      [cw],$0F00
    fldcw   [cw]
    fistp   qword ptr [iarr]
    fwait
    fldcw   [cws]
  end;
  trunc := iarr[0];
end;
{$endif}


{$ifdef buggy_round}
{$ifdef FPC}
{---------------------------------------------------------------------------}
function round(x: extended): longint;
var
  iarr: packed array[0..1] of longint;
begin
  {This is a fix for the buggy FPC1 code}
  asm
    fld     [x]
    fistp   qword ptr [iarr]
    fwait
  end;
  round := iarr[0];
end;
{$endif}
{$ifdef VirtualPascal}
{---------------------------------------------------------------------------}
function round(x: extended): longint; assembler; {&Frame-} {&Uses none}
var
  TempLong: Longint;
  {VP does not round to even e.g. for x=4.5 by design(!?), see system.pas}
asm
  fld     [x]
  fistp   TempLong
  fwait
  mov     eax,TempLong
end;
{$endif}
{$endif}

const
  ln2_hi: THexDblA = ($00,$00,$E0,$FE,$42,$2E,$E6,$3F);
  ln2_lo: THexDblA = ($76,$3C,$79,$35,$EF,$39,$EA,$3D);


{$ifndef use_fast_exp}
const
  half: single = 0.5;     {used for roundmode-safe exp routines}
  two : single = 2.0;     {used for roundmode-safe exp routines}
  ebig: single = 24576.0; {used for roundmode-safe exp routines}
{$endif}

{---------------------------------------------------------------------------}
function exp(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Accurate exp, result good to extended precision}
asm
  {This version of Norbert Juffa's exp is from the VirtualPascal RTL source,}
  {discussed and explained in the VP Bugtracker system. Quote:              }
  {                                                                         }
  { ... "since the 387, F2XM1 can accecpt arguments in [-1, 1].             }
  {                                                                         }
  { So, we can split the argument into an integer and a fraction part using }
  { FRNDINT and the fraction part will always be -1 <= f <= 1 no matter what}
  { rounding control. This means we don't have to load/restore the FPU      }
  { control word (CW) which is slow on modern OOO FPUs (since FLDCW is a    }
  { serializing instruction).                                               }
  {                                                                         }
  { Note that precision is lost in doing exponentation when the fraction is }
  { subtracted from the integer part of the argument. The "naive" code can  }
  { loose up to 11 (or 15) bits of the extended precision format for large  }
  { DP or EP arguments, yielding a result good to double precision. To get a}
  { function accurate to full extended precision, we need to simulate higher}
  { precision intermediate arithmetic."                                     }
  { Ref: [Virtual Pascal 0000056]: More accurate Exp() function. URL (Oct.2009):}
  { https://admin.topica.com/lists/virtualpascal@topica.com/read/message.html?sort=a&mid=908867704&start=7}

  fld     [x]                { x                                                 }
  fldl2e                     { log2(e) | x                                       }
  fmul    st,st(1)           { z = x * log2(e) x                                 }
  frndint                    { int(z) | x                                        }
  fld     qword ptr [ln2_hi] { ln2_hi | int(z) | x                               }
  fmul    st,st(1)           { int(z)*ln2_hi | int(z) | x                        }
  fsubp   st(2),st           { int(z) | x-int(z)*ln2_hi                          }
  fld     qword ptr [ln2_lo] { ln2_lo | int(z) | x-int(z)*ln2_hi                 }
  fmul    st, st(1)          { int(z)*ln2_lo | int(z) | x-int(z)*ln2_hi          }
  fsubp   st(2),st           { int(z) | (x-int(z)*ln2_hi)-int(z)*ln2_lo          }
  fxch    st(1)              { (x-int(z)*ln2_hi)-int(z)*ln2_lo | int(z)          }
  fldl2e                     { log2(e) | (x-int(z)*ln2_hi)-int(z)*ln2_lo | int(z)}
  fmulp   st(1),st           { frac(z) | int(z)                                  }

{$ifndef use_fast_exp}
  {It may happen (especially for rounding modes other than "round to nearest")   }
  {that |frac(z)| > 1. In this case the result of f2xm1 is undefined. The next   }
  {lines will test frac(z) and use a safe algorithm if necessary.                }
  {Another problem pops up if x is very large e.g. for x=1e3000. AMath checks    }
  {int(z) and returns 2^int(z) if int(z) > 1.5*16384, result is 0 or overflow!   }
  fld     st
  fabs                       { abs(frac(z)) | frac(z) | int(z)                   }
  fld1                       { 1 | abs(frac(z)) | frac(z) | int(z)               }
  fcompp
  fstsw   ax
  sahf
  jae     @@1                { frac(z) <= 1, no special action needed            }
  fld     st(1)              { int(z) | frac(z) | int(z)                         }
  fabs                       { abs(int(z)) | frac(z) | int(z)                    }
  fcomp   [ebig]
  fstsw   ax
  sahf
  jb      @@0
  fsub    st,st              { set frac=0 and scale with too large int(z)}
  jmp     @@1
@@0:
  {Safely calculate 2^frac(z)-1 as (2^(frac(z)/2)-1)*(2^(frac(z)/2)+1) and use   }
  {2^(frac(z)/2)+1 = (2^(frac(z)/2)-1) + 2 (suggested by N. Juffa, 16.Jan.2011)  }
  fmul    dword ptr [half]   { frac(z)/2  | int(z)                               }
  f2xm1                      { 2^(frac(z)/2)-1 | int(z)                          }
  fld     st                 { 2^(frac(z)/2)-1 | 2^(frac(z)/2)-1 | int(z)        }
  fadd    dword ptr [two]    { 2^(frac(z)/2)+1 | 2^(frac(z)/2)-1 | int(z)        }
  fmulp   st(1),st           { 2^frac(z)-1 | int(z)                              }
  jmp     @@2
{$endif}

@@1:
  f2xm1                      { 2^frac(z)-1 | int(z)                              }

@@2:
  fld1                       { 1 | 2^frac(z)-1 | int(z)                          }
  faddp   st(1),st           { 2^frac(z) | int(z)                                }
  fscale                     { 2^z | int(z)                                      }
  fstp    st(1)              { 2^z = e^x                                         }
  fwait
end;


{---------------------------------------------------------------------------}
function exp10(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return 10^x}
const
  lg2_hi: THexDblA = ($00,$00,$80,$50,$13,$44,$D3,$3F);
  lg2_lo: THexDblA = ($2B,$F1,$11,$F3,$FE,$79,$DF,$3D);
asm
  fld     [x]                { x                                                  }
  fldl2t                     { log2(10) | x                                       }
  fmul    st,st(1)           { z = x * log2(10) | x                               }
  frndint                    { int(z) | x                                         }
  fld     qword ptr [lg2_hi] { lg2_hi | int(z) | x                                }
  fmul    st,st(1)           { int(z)*lg2_hi | int(z) | x                         }
  fsubp   st(2),st           { int(z) | x-int(z)*lg2_hi                           }
  fld     qword ptr [lg2_lo] { lg2_lo | int(z) | x-int(z)*lg2_hi                  }
  fmul    st, st(1)          { int(z)*lg2_lo | int(z) | x-int(z)*lg2_hi           }
  fsubp   st(2),st           { int(z) | (x-int(z)*lg2_hi)-int(z)*lg2_lo           }
  fxch    st(1)              { (x-int(z)*lg2_hi)-int(z)*lg2_lo | int(z)           }
  fldl2t                     { log2(10) | (x-int(z)*lg2_hi)-int(z)*lg2_lo | int(z)}
  fmulp   st(1),st           { frac(z) | int(z)                                   }

{$ifndef use_fast_exp}
  {See the exp code for a description of these conditional lines}
  fld     st
  fabs                       { abs(frac(z)) | frac(z) | int(z)                   }
  fld1                       { 1 | abs(frac(z)) | frac(z) | int(z)               }
  fcompp
  fstsw   ax
  sahf
  jae     @@1                { frac(z) <= 1, no special action needed            }
  fld     st(1)              { int(z) | frac(z) | int(z)                         }
  fabs                       { abs(int(z)) | frac(z) | int(z)                    }
  fcomp   [ebig]
  fstsw   ax
  sahf
  jb      @@0
  fsub    st,st              { set frac=0 and scale with too large int(z)}
  jmp     @@1
@@0:
  fmul    dword ptr [half]   { frac(z)/2  | int(z)                               }
  f2xm1                      { 2^(frac(z)/2)-1 | int(z)                          }
  fld     st                 { 2^(frac(z)/2)-1 | 2^(frac(z)/2)-1 | int(z)        }
  fadd    dword ptr [two]    { 2^(frac(z)/2)+1 | 2^(frac(z)/2)-1 | int(z)        }
  fmulp   st(1),st           { 2^frac(z)-1 | int(z)                              }
  jmp     @@2
{$endif}

@@1:
  f2xm1                      { 2^frac(z)-1 | int(z)                              }

@@2:
  fld1                       { 1 | 2^frac(z)-1 | int(z)                          }
  faddp   st(1),st           { 2^frac(z) | int(z)                                }
  fscale                     { 2^z | int(z)                                      }
  fstp    st(1)              { 2^z = 10^x                                         }
  fwait
end;


{---------------------------------------------------------------------------}
function exp2(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return 2^x}
asm
  fld     [x]                { x                       }
  fld     st(0)              { x | x                   }
  frndint                    { int(x) | x              }
  fxch    st(1)              { x | int(x)              }
  fsub    st(0),st(1)        { frac(x) | int(x)        }
  f2xm1                      { 2^frac(x)-1 | int(x)    }
  fld1                       { 1 | 2^frac(x)-1 | int(x)}
  faddp   st(1),st           { 2^frac(x) | int(x)      }
  fscale                     { 2^z | int(x)            }
  fstp    st(1)
  fwait
end;


{---------------------------------------------------------------------------}
function exp3(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return 3^x}
const
  l32_hi: THexDblA = ($00,$00,$30,$98,$93,$30,$E4,$3F);
  l32_lo: THexDblA = ($58,$1F,$02,$A7,$F4,$D4,$C4,$3D);
  log23h: THexExtW = ($43D0,$FDEB,$0D1C,$CAE0,$3FFF);  {1.5849625007211561815}
asm
  fld     [x]                { x                                                  }
  fld     tbyte ptr [log23h] { log2(3) | x                                        }
  fmul    st,st(1)           { z = x * log2(3) | x                                }
  frndint                    { int(z) | x                                         }
  fld     qword ptr [l32_hi] { l32_hi | int(z) | x                                }
  fmul    st,st(1)           { int(z)*l32_hi | int(z) | x                         }
  fsubp   st(2),st           { int(z) | x-int(z)*l32_hi                           }
  fld     qword ptr [l32_lo] { l32_lo | int(z) | x-int(z)*l32_hi                  }
  fmul    st, st(1)          { int(z)*l32_lo | int(z) | x-int(z)*l32_hi           }
  fsubp   st(2),st           { int(z) | (x-int(z)*l32_hi)-int(z)*l32_lo           }
  fxch    st(1)              { (x-int(z)*l32_hi)-int(z)*l32_lo | int(z)           }
  fld     tbyte ptr [log23h] { log2(3) | (x-int(z)*l32_hi)-int(z)*l32_lo | int(z) }
  fmulp   st(1),st           { frac(z) | int(z)                                   }

{$ifndef use_fast_exp}
  {See the exp code for a description of these conditional lines}
  fld     st
  fabs                       { abs(frac(z)) | frac(z) | int(z)                   }
  fld1                       { 1 | abs(frac(z)) | frac(z) | int(z)               }
  fcompp
  fstsw   ax
  sahf
  jae     @@1                { frac(z) <= 1, no special action needed            }
  fld     st(1)              { int(z) | frac(z) | int(z)                         }
  fabs                       { abs(int(z)) | frac(z) | int(z)                    }
  fcomp   [ebig]
  fstsw   ax
  sahf
  jb      @@0
  fsub    st,st              { set frac=0 and scale with too large int(z)}
  jmp     @@1
@@0:
  fmul    dword ptr [half]   { frac(z)/2  | int(z)                               }
  f2xm1                      { 2^(frac(z)/2)-1 | int(z)                          }
  fld     st                 { 2^(frac(z)/2)-1 | 2^(frac(z)/2)-1 | int(z)        }
  fadd    dword ptr [two]    { 2^(frac(z)/2)+1 | 2^(frac(z)/2)-1 | int(z)        }
  fmulp   st(1),st           { 2^frac(z)-1 | int(z)                              }
  jmp     @@2
{$endif}

@@1:
  f2xm1                      { 2^frac(z)-1 | int(z)                              }

@@2:
  fld1                       { 1 | 2^frac(z)-1 | int(z)                          }
  faddp   st(1),st           { 2^frac(z) | int(z)                                }
  fscale                     { 2^z | int(z)                                      }
  fstp    st(1)              { 2^z = 3^x                                          }
  fwait
end;


{---------------------------------------------------------------------------}
function exp5(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return 5^x}
const
  l52_hi: THexDblA = ($00,$00,$00,$69,$34,$90,$DB,$3F);
  l52_lo: THexDblA = ($22,$D5,$33,$94,$CB,$3D,$B4,$3D);
  log25h: THexExtW = ($8AFE,$CD1B,$784B,$949A,$4000);  {2.3219280948873623478}
asm
  fld     [x]                { x                                                  }
  fld     tbyte ptr [log25h] { log2(5) | x                                        }
  fmul    st,st(1)           { z = x * log2(5) | x                                }
  frndint                    { int(z) | x                                         }
  fld     qword ptr [l52_hi] { l52_hi | int(z) | x                                }
  fmul    st,st(1)           { int(z)*l52_hi | int(z) | x                         }
  fsubp   st(2),st           { int(z) | x-int(z)*l52_hi                           }
  fld     qword ptr [l52_lo] { l52_lo | int(z) | x-int(z)*l52_hi                  }
  fmul    st, st(1)          { int(z)*l52_lo | int(z) | x-int(z)*l52_hi           }
  fsubp   st(2),st           { int(z) | (x-int(z)*l52_hi)-int(z)*l52_lo           }
  fxch    st(1)              { (x-int(z)*l52_hi)-int(z)*l52_lo | int(z)           }
  fld     tbyte ptr [log25h] { log2(5) | (x-int(z)*l52_hi)-int(z)*l52_lo | int(z) }
  fmulp   st(1),st           { frac(z) | int(z)                                   }

{$ifndef use_fast_exp}
  {See the exp code for a description of these conditional lines}
  fld     st
  fabs                       { abs(frac(z)) | frac(z) | int(z)                   }
  fld1                       { 1 | abs(frac(z)) | frac(z) | int(z)               }
  fcompp
  fstsw   ax
  sahf
  jae     @@1                { frac(z) <= 1, no special action needed            }
  fld     st(1)              { int(z) | frac(z) | int(z)                         }
  fabs                       { abs(int(z)) | frac(z) | int(z)                    }
  fcomp   [ebig]
  fstsw   ax
  sahf
  jb      @@0
  fsub    st,st              { set frac=0 and scale with too large int(z)}
  jmp     @@1
@@0:
  fmul    dword ptr [half]   { frac(z)/2  | int(z)                               }
  f2xm1                      { 2^(frac(z)/2)-1 | int(z)                          }
  fld     st                 { 2^(frac(z)/2)-1 | 2^(frac(z)/2)-1 | int(z)        }
  fadd    dword ptr [two]    { 2^(frac(z)/2)+1 | 2^(frac(z)/2)-1 | int(z)        }
  fmulp   st(1),st           { 2^frac(z)-1 | int(z)                              }
  jmp     @@2
{$endif}

@@1:
  f2xm1                      { 2^frac(z)-1 | int(z)                              }

@@2:
  fld1                       { 1 | 2^frac(z)-1 | int(z)                          }
  faddp   st(1),st           { 2^frac(z) | int(z)                                }
  fscale                     { 2^z | int(z)                                      }
  fstp    st(1)              { 2^z = 5^x                                          }
  fwait
end;


{---------------------------------------------------------------------------}
function exp7(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return 7^x}
const
  l72_hi: THexDblA = ($00,$00,$00,$3A,$19,$CC,$D6,$3F);
  l72_lo: THexDblA = ($C1,$DE,$44,$9C,$36,$D5,$09,$3E);
  log27h: THexExtW = ($66CD,$A021,$B3FA,$B3AB,$4000);  {2.8073549220576041075}
asm
  fld     [x]                { x                                                  }
  fld     tbyte ptr [log27h] { log2(7) | x                                        }
  fmul    st,st(1)           { z = x * log2(7) | x                                }
  frndint                    { int(z) | x                                         }
  fld     qword ptr [l72_hi] { l72_hi | int(z) | x                                }
  fmul    st,st(1)           { int(z)*l72_hi | int(z) | x                         }
  fsubp   st(2),st           { int(z) | x-int(z)*l72_hi                           }
  fld     qword ptr [l72_lo] { l72_lo | int(z) | x-int(z)*l72_hi                  }
  fmul    st, st(1)          { int(z)*l72_lo | int(z) | x-int(z)*l72_hi           }
  fsubp   st(2),st           { int(z) | (x-int(z)*l72_hi)-int(z)*l72_lo           }
  fxch    st(1)              { (x-int(z)*l72_hi)-int(z)*l72_lo | int(z)           }
  fld     tbyte ptr [log27h] { log2(7) | (x-int(z)*l72_hi)-int(z)*l72_lo | int(z) }
  fmulp   st(1),st           { frac(z) | int(z)                                   }

{$ifndef use_fast_exp}
  {See the exp code for a description of these conditional lines}
  fld     st
  fabs                       { abs(frac(z)) | frac(z) | int(z)                   }
  fld1                       { 1 | abs(frac(z)) | frac(z) | int(z)               }
  fcompp
  fstsw   ax
  sahf
  jae     @@1                { frac(z) <= 1, no special action needed            }
  fld     st(1)              { int(z) | frac(z) | int(z)                         }
  fabs                       { abs(int(z)) | frac(z) | int(z)                    }
  fcomp   [ebig]
  fstsw   ax
  sahf
  jb      @@0
  fsub    st,st              { set frac=0 and scale with too large int(z)}
  jmp     @@1
@@0:
  fmul    dword ptr [half]   { frac(z)/2  | int(z)                               }
  f2xm1                      { 2^(frac(z)/2)-1 | int(z)                          }
  fld     st                 { 2^(frac(z)/2)-1 | 2^(frac(z)/2)-1 | int(z)        }
  fadd    dword ptr [two]    { 2^(frac(z)/2)+1 | 2^(frac(z)/2)-1 | int(z)        }
  fmulp   st(1),st           { 2^frac(z)-1 | int(z)                              }
  jmp     @@2
{$endif}

@@1:
  f2xm1                      { 2^frac(z)-1 | int(z)                              }

@@2:
  fld1                       { 1 | 2^frac(z)-1 | int(z)                          }
  faddp   st(1),st           { 2^frac(z) | int(z)                                }
  fscale                     { 2^z | int(z)                                      }
  fstp    st(1)              { 2^z = 7^x                                          }
  fwait
end;


{---------------------------------------------------------------------------}
function fmod(x,y: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return x mod y, y<>0, sign(result) = sign(x)}
asm
       fld    [y]
       fld    [x]
  @@1: fprem
       fstsw  ax
       sahf
       jp     @@1
       fstp   st(1)
       fwait
end;


{---------------------------------------------------------------------------}
function ln(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return natural logarithm of x}
asm
  fldln2
  fld     [x]
  fyl2x
  fwait
end;


{---------------------------------------------------------------------------}
function log10(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return base 10 logarithm of x}
asm
  fldlg2
  fld     [x]
  fyl2x
  fwait
end;


{---------------------------------------------------------------------------}
function log2(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return base 2 logarithm of x}
asm
  fld1
  fld     [x]
  fyl2x
  fwait
end;


{---------------------------------------------------------------------------}
function logbase(b, x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return base b logarithm of x}
asm
  fld1
  fld     [x]
  fyl2x
  fld1
  fld     [b]
  fyl2x
  fdivp   st(1),st
  fwait
end;


{---------------------------------------------------------------------------}
function logN(N, x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Delphi alias for logbase}
asm
  fld1
  fld     [x]
  fyl2x
  fld1
  fld     [N]
  fyl2x
  fdivp   st(1),st
  fwait
end;


{---------------------------------------------------------------------------}
function arctan2(y, x: extended): extended;  assembler; {&Frame-} {&Uses none}
  {-Return arctan(y/x); result in [-Pi..Pi] with correct quadrant}
asm
  fld     [y]
  fld     [x]
  fpatan
  fwait
end;

{$ifdef BASM16}
{---------------------------------------------------------------------------}
function intpower0(x: extended; n: longint): extended; assembler;
  {-Return x^n; via binary exponentiation (no overflow detection)}
asm
     db $66; mov   ax,word ptr [n]
             fld1            {r := 1 }
     db $66; cwd             {cdq!, edx=-1 if n<0, 0 otherwise}
     db $66; xor   ax,dx
     db $66; sub   ax,dx     {eax == i := abs(n)}
             jz    @@3       {if n=0 done}
             fld   [x]
             jmp   @@2
@@1:         fmul  st,st     {x := x*x}
@@2: db $66; shr   ax,1      {i := i shr 1}
             jnc   @@1       {i had not been odd, repeat squaring}
             fmul  st(1),st  {r := r*x}
             jnz   @@1       {if i=0 exit loop}
             fstp  st        {pop x}
             or    dx,dx     {n<0 iff dx<>0}
             jz    @@3
             fld1
             fdivr           {intpower := 1/r}
@@3:         fwait
end;


{---------------------------------------------------------------------------}
function ldexp(x: extended; e: longint): extended; assembler;
  {-Return x*2^e}
asm
  fild    [e]
  fld     [x]
  fscale
  fstp    st(1)
end;


{---------------------------------------------------------------------------}
function ldexpd(d: double; e: longint): double; assembler;
  {-Return d*2^e}
asm
  fild    [e]
  fld     [d]
  fscale
  fstp    st(1)
end;


{---------------------------------------------------------------------------}
function ldexps(s: single; e: longint): single; assembler;
  {-Return s*2^e}
asm
  fild    [e]
  fld     [s]
  fscale
  fstp    st(1)
end;


{---------------------------------------------------------------------------}
function scalbn(x: extended; e: longint): extended; assembler;
  {-Return x*2^e}
asm
  fild    [e]
  fld     [x]
  fscale
  fstp    st(1)
end;

{$else}

{---------------------------------------------------------------------------}
function intpower0(x: extended; n: longint): extended; assembler; {&Frame-} {&Uses none}
  {-Return x^n; via binary exponentiation (no overflow detection)}
asm
       mov    eax,n     {Note: this may be a mov eax,eax in Delphi/FPC}
                        {but it is needed in VirtualPascal}
       fld1             {r := 1 }
       cdq              {edx=-1 if n<0, 0 otherwise}
       xor    eax,edx
       sub    eax,edx   {eax == i := abs(n)}
       jz     @@3       {if n=0 done}
       fld    [x]
       jmp    @@2
  @@1: fmul   st,st     {x := x*x }
  @@2: shr    eax,1     {i := i shr 1}
       jnc    @@1       {i had not been odd, repeat squaring}
       fmul   st(1),st  {r := r*x }
       jnz    @@1       {if i=0 exit loop}
       fstp   st        {pop x}
       or     edx,edx   {n<0 iff dx<>0}
       jz     @@3
       fld1
       fdivrp st(1),st  {intpower := 1/r; FPC does not like fdivr}
  @@3: fwait
end;


{$ifdef NO_FPU_RESULT}
 {---------------------------------------------------------------------------}
 function ldexp(x: extended; e: longint): extended;
   {-Return x*2^e}
 begin
   asm
     fild    [e]
     fld     [x]
     fscale
     fstp    st(1)
     fstp    [x]
     fwait
   end;
   ldexp := x;
 end;

 {---------------------------------------------------------------------------}
 function ldexpd(d: double; e: longint): double;
   {-Return d*2^e}
 begin
   asm
     fild    [e]
     fld     [d]
     fscale
     fstp    st(1)
     fstp    [d]
     fwait
   end;
   ldexpd := d;
 end;


 {---------------------------------------------------------------------------}
 function ldexps(s: single; e: longint): single;
   {-Return s*2^e}
 begin
   asm
     fild    [e]
     fld     [s]
     fscale
     fstp    st(1)
     fstp    [s]
     fwait
   end;
   ldexps := s;
 end;


 {---------------------------------------------------------------------------}
 function scalbn(x: extended; e: longint): extended;
   {-Return x*2^e}
 begin
   asm
     fild    [e]
     fld     [x]
     fscale
     fstp    st(1)
     fstp    [x]
     fwait
   end;
   scalbn := x;
 end;

{$else}

 {---------------------------------------------------------------------------}
 function ldexp(x: extended; e: longint): extended;
   {-Return x*2^e}
 begin
   asm
     fild    [e]
     fld     [x]
     fscale
     fstp    st(1)
     fstp    [@result]
     fwait
   end;
 end;


 {---------------------------------------------------------------------------}
 function ldexpd(d: double; e: longint): double;
   {-Return d*2^e}
 begin
   asm
     fild    [e]
     fld     [d]
     fscale
     fstp    st(1)
     fstp    [@result]
     fwait
   end;
 end;


 {---------------------------------------------------------------------------}
 function ldexps(s: single; e: longint): single;
   {-Return s*2^e}
 begin
   asm
     fild    [e]
     fld     [s]
     fscale
     fstp    st(1)
     fstp    [@result]
     fwait
   end;
 end;


 {---------------------------------------------------------------------------}
 function scalbn(x: extended; e: longint): extended;
   {-Return x*2^e}
 begin
   asm
     fild    [e]
     fld     [x]
     fscale
     fstp    st(1)
     fstp    [@result]
     fwait
   end;
 end;
{$endif}


{$endif}


{---------------------------------------------------------------------------}
function _tan(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return the circular tangent of x, x mod Pi <> Pi/2}
asm
  {tan := sin(x)/cos(x)}
  fld    [x]
  fptan
  fstp   st(0)
  fwait
end;


{---------------------------------------------------------------------------}
function _cot(x: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return the circular cotangent of x, x mod Pi <> 0}
asm
  {cot := cos(x)/sin(x) = 1/tan(x)}
  fld    [x]
  fptan
  fdivrp st(1),st
  fwait
end;


{$ifdef BIT16}
{---------------------------------------------------------------------------}
procedure _sincos(x: extended; var s,c: extended); assembler;
asm
  fld     [x]
  db      $D9,$FB  {fsincos}
  les     di,[c]
  fstp    es:tbyte ptr[di]
  les     di,[s]
  fstp    es:tbyte ptr[di]
  fwait
end;
{$else}
{---------------------------------------------------------------------------}
procedure _sincos(x: extended; var s,c: extended);
begin
  {Delphi uses stack frame (for x). So avoid problems with}
  {register convention in VPC and always use stack frame}
  asm
    fld     [x]
    fsincos
    mov     ecx, c
    fstp    tbyte ptr [ecx]
    mov     ecx, s
    fstp    tbyte ptr [ecx]
    fwait
  end;
end;
{$endif}


{---------------------------------------------------------------------------}
function expm1(x: extended): extended;
  {-Return exmp1(x)-1, accurate even for x near 0}
const
  nce=17;
  ceh: array[0..nce-1] of THexExtW = (     {chebyshev(((exp(x)-1)/x - 1)/x, x=-1..1, 0.1e-20);}
         ($6A89,$722B,$FAC5,$8577,$3FFF),  {+1.04272398606220957726049179176     }
         ($D006,$0AFA,$F911,$B131,$3FFC),  {+0.173042194047179631675883846985    }
         ($D2C9,$D68A,$A866,$B073,$3FF9),  {+0.215395249458237651324867090935e-1 }
         ($B98A,$635F,$18B1,$8CA8,$3FF6),  {+0.214624979834132263391667992385e-2 }
         ($6731,$7DD9,$FF11,$BAFB,$3FF2),  {+0.178322182449141937355350007488e-3 }
         ($4297,$88F5,$C8BC,$D52A,$3FEE),  {+0.127057507929428665687742064639e-4 }
         ($2C99,$3975,$4676,$D4B9,$3FEA),  {+0.792457652881593907264616203751e-6 }
         ($196B,$F18F,$09DF,$BCC1,$3FE6),  {+0.439477285666343943629044429651e-7 }
         ($A60E,$D37D,$5A0B,$96C6,$3FE2),  {+0.219406227546144888207254174079e-8 }
         ($F7F8,$EF82,$8D42,$DB05,$3FDD),  {+0.995995318266659554765505887521e-10}
         ($CB48,$CD67,$F8D6,$91D8,$3FD9),  {+0.414523660148535566445986193898e-11}
         ($9AFA,$57D1,$F72D,$B352,$3FD4),  {+0.159271781651030822298380021406e-12}
         ($560F,$44FD,$5097,$CCC2,$3FCF),  {+0.568320507930678710973086211864e-14}
         ($33CD,$BF09,$610D,$DA3C,$3FCA),  {+0.189289431283787337567205644620e-15}
         ($7E97,$7C7D,$7B0A,$DA14,$3FC5),  {+0.591107031096332196225712012488e-17}
         ($1F08,$8285,$97C9,$CD1E,$3FC0),  {+0.173742977663542892171117958472e-18}
         ($0AA6,$23C8,$CF63,$B638,$3FBB)); {+0.482337391484586252954369194368e-20}
var
  ces: array[0..nce-1] of extended absolute ceh;
  t: extended;
begin
  if abs(x)<=1.0 then begin
    t := CSEvalx(x, ces, nce);
    expm1 := x + x*x*t;
  end
  else expm1 := exp(x)-1.0;
end;


{---------------------------------------------------------------------------}
function ln1p(x: extended): extended;
  {-Return ln(1+x), accurate even for x near 0}
const
  PH: array[0..8] of THexExtW = (
        ($0000,$0000,$0000,$8000,$3FFE),  {+0.5                                }
        ($71C7,$C71C,$1C71,$F1C7,$3FFF),  {+1.88888888888888888888888888889    }
        ($7E0C,$619A,$EFD3,$B8B6,$4000),  {+2.88616557734204793028322440087    }
        ($D227,$277C,$7CD2,$9227,$4000),  {+2.28366013071895424836601307190    }
        ($FE54,$53A8,$A8FE,$FE53,$3FFE),  {+0.993464052287581699346405228758   }
        ($42C5,$1A22,$7798,$ED6F,$3FFC),  {+0.231870525988173046996576408341   }
        ($D42F,$2F24,$24D4,$D42F,$3FF9),  {+0.259013861955038425626660920779e-1}
        ($83CC,$CC00,$0083,$83CC,$3FF5),  {+0.100553041729512317747611865259e-2}
        ($E03E,$D40A,$BFE2,$8852,$3FE9)); {+0.253921822549273529665686528432e-6}
const
  QH: array[0..8] of THexExtW = (
        ($0000,$0000,$0000,$8000,$3FFF),  {+1.0                                 }
        ($8E39,$38E3,$E38E,$8E38,$4001),  {+4.44444444444444444444444444444     }
        ($C3C4,$C3C3,$C3C3,$83C3,$4002),  {+8.23529411764705882352941176471     }
        ($C3C4,$C3C3,$C3C3,$83C3,$4002),  {+8.23529411764705882352941176471     }
        ($B9BA,$B9B9,$B9B9,$99B9,$4001),  {+4.80392156862745098039215686275     }
        ($D2D3,$D2D2,$D2D2,$D2D2,$3FFF),  {+1.64705882352941176470588235294     }
        ($A22C,$2C04,$04A2,$A22C,$3FFD),  {+0.316742081447963800904977375566    }
        ($F71F,$1E80,$80F7,$F71E,$3FF9),  {+0.301659125188536953242835595777e-1 }
        ($40F9,$27E9,$D1FB,$86CA,$3FF5)); {+0.102838338132455779514603044015e-2 }
var
  P: array[0..8] of extended absolute PH;
  Q: array[0..8] of extended absolute QH;
var
  y: extended;
begin
  if x >= 4.0 then ln1p := ln(1.0+x)
  else if (x>=-0.25) and (x<=0.375) then begin
    {Pade approximation of degree (8,8) for y(x) = -(ln(1+x)/x-1)/x}
    {Errors: -0.21676209e-19 for x=0.25, 0.70343473e-19 for x=0.375}
    y := PolEvalX(x,P,9)/PolEvalX(x,Q,9);
    ln1p := x - y*x*x;
  end
  else begin
    y := 1.0 + x;
    {The following formula is more accurate than Goldberg [3], Theorem 4. The}
    {Taylor series f(x) = f(x0) + (x-x0)*f'(x0) + .. for f(x) = ln(1+x) gives}
    {ln1p(x) = ln(1+x0) + (x-x0)/(1+x0) = ln(y) + (x-(y-1))/y, with y = 1+x0.}
    if y=1.0 then ln1p := x
    else ln1p := ln(y) + (x-(y-1.0))/y;
  end;
end;


{---------------------------------------------------------------------------}
function exp10m1(x: extended): extended;
  {-Return 10^x - 1; special code for small x}
var
  z: extended;
begin
  if abs(x) > 0.3 then exp10m1 := exp10(x) - 1.0
  else if x=0.0 then exp10m1 := 0.0
  else begin
    {z = x*log2(10), 10^x-1 = exp2m1(x*log2(10))}
    z := 3.0*x + 0.32192809488736234787*x;
    exp10m1 := exp2m1(z);
  end
end;


{---------------------------------------------------------------------------}
function exp2m1(x: extended): extended;
  {-Return 2^x-1, accurate even for x near 0}
begin
  if abs(x)<=1 then begin
    asm
      fld   [x]
      f2xm1
      {$ifdef NO_FPU_RESULT}
        fstp   [x]
      {$else}
        fstp   [@result]
      {$endif}
    end;
    {$ifdef NO_FPU_RESULT}
      exp2m1 := x;
    {$endif}
  end
  else exp2m1 := exp2(x)-1.0;
end;


{---------------------------------------------------------------------------}
function exprel(x: extended): extended;
  {-Return exprel(x) = (exp(x) - 1)/x,  1 for x=0}
const
  xsmall = 1e-6;
var
  z: extended;
begin
  z := abs(x);
  if z < ln2 then begin
    if z < xsmall then  begin
      {exprel = 1 + 1/2*x + 1/6*x^2 + 1/24*x^3 + O(x^4) }
      exprel := 1.0 + x*(0.5 + x/SIXX);
    end
    else begin
      {See Higham [9], 1.14.1 Computing (e^x-1)/x, Algorithm 2}
      z := exp(x);
      x := ln(z);
      if x=0.0 then exprel := 1.0
      else exprel := (z-1.0)/x;
    end;
  end
  else begin
    if x<-45.0 then exprel := -1.0/x
    else exprel := (exp(x) - 1.0)/x;
  end;
end;


{---------------------------------------------------------------------------}
function expx2(x: extended): extended;
  {-Return exp(x*|x|) with damped error amplification in computing exp of the product.}
  { Used for exp(x^2) = expx2(abs(x)) and exp(-x^2) = expx2(-abs(x))}
const
  MFAC = 32768.0;
  MINV = 3.0517578125e-5;
var
  u,u1,m,f: extended;
  neg: boolean;
begin
  if x >= 106.5669902282 then begin
    expx2 := PosInf_x;
    exit;
  end
  else if x <= -106.77 then begin
    expx2 := 0.0;
    exit;
  end;

  {Ref: Cephes [7], file ldouble\expx2l.c}
  neg := x<0;
  x := abs(x);
  if x <= 1.0 then begin
    if neg then u := -x*x else u := x*x;
    expx2 := exp(u);
  end
  else begin
    {Represent x as an exact multiple of MFAC plus a residual.}
    {MFAC is a power of 2 chosen so that exp(m * m) does not  }
    {overflow or underflow and so that |x - m| is small.      }
    m := MINV*floorx(MFAC*x + 0.5);
    f := x - m;

    {x^2 = m^2 + 2mf + f^2}
    u  := m*m;
    u1 := 2.0*m*f + f*f;

    if neg then begin
      u  := -u;
      u1 := -u1;
    end;

    if u+u1 > ln_MaxExt then expx2 := PosInf_x
    else begin
      {u is exact, u1 is small}
      expx2 := exp(u)*exp(u1);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function expmx2h(x: extended): extended;
  {-Return exp(-0.5*x^2) with damped error amplification}
const
  MFAC = 32768.0;
  MINV = 3.0517578125e-5;
var
  u,u1,m,f: extended;
begin
  {Ref: Cephes [7], file ldouble\expx2l.c}
  x := abs(x);
  if x >= 150.994 then expmx2h := 0.0
  else if x <= 2.0 then begin
    expmx2h := exp(-0.5*x*x);
  end
  else begin
    {Represent x as an exact multiple of MFAC plus a residual.}
    {MFAC is a power of 2 chosen so that exp(m * m) does not  }
    {overflow or underflow and so that |x - m| is small.      }
    m := MINV*floorx(MFAC*x + 0.5);
    f := x - m;
    {0.5*x^2 = 0.5*(m^2 + mf + f^2)}
    u  := (-0.5)*m*m;
    u1 := (-0.5)*f*f - m*f;
    {u is exact, u1 is small}
    expmx2h := exp(u)*exp(u1);
  end;
end;


(*
Transformation/reduction z = x+k*pi/2: cf. HMF [1], 4.3.44

  k   sin(z)   cos(z)   tan(z)   cot(z)
 ---------------------------------------
  0      s        c        t     cot(x)
  1      c       -s     -1/t       -t
  2     -s       -c        t     cot(x)
  3     -c        s     -1/t       -t
 ---------------------------------------

 with s=sin(x), c=cos(x), t=tan(x)
*)

{---------------------------------------------------------------------------}
function cos(x: extended): extended;
  {-Accurate version of circular cosine, uses system.cos for |x| <= Pi/4}
var
  t: extended;
begin
  {reduction mod Pi/2, |t| <= Pi/4}
  case rem_pio2(x,t) and 3 of
      0:  cos :=  system.cos(t);
      1:  cos := -system.sin(t);
      2:  cos := -system.cos(t);
    else  cos :=  system.sin(t);
  end;
end;


const
  npcv=7; {chebyshev(((1-cos(x))-x^2/2)/x^4,x=-Pi/4..Pi/4,1e-20) converted to poly}
  pcvh: array[0..npcv-1] of THexExtW = (
          ($AAAB,$AAAA,$AAAA,$AAAA,$BFFA),  {-0.416666666666666666661679544980e-1 }
          ($B26B,$0B60,$60B6,$B60B,$3FF5),  {+0.138888888888888879063613392395e-2 }
          ($B143,$0CE8,$00D0,$D00D,$BFEF),  {-0.248015873015846864625912712986e-4 }
          ($5CC9,$89D7,$7DBB,$93F2,$3FE9),  {+0.275573192214211550883212170218e-6 }
          ($6C85,$4AA7,$C6F6,$8F76,$BFE2),  {-0.208767557953773544469035900043e-8 }
          ($BDBB,$E308,$5D99,$C9CA,$3FDA),  {+0.114704613872478358447153370873e-10}
          ($234A,$7E37,$70F1,$D5BC,$BFD2)); {-0.474589475226674827431568738364e-13}
var
  pcv: array[0..npcv-1] of extended absolute pcvh;


{---------------------------------------------------------------------------}
function vers(x: extended): extended;
  {-Return the versine vers(x) = 1 - cos(x)}
var
  t: extended;
begin
  {reduction mod Pi/2, |t| <= Pi/4}
  case rem_pio2(x,t) and 3 of
      0:  begin
            {vers(x) = 1 - cos(t)}
            x := t*t;
            t := PolEvalX(x,pcv,npcv);
            vers := 0.5*x + sqr(x)*t;
          end;
      1:  vers := 1.0 + system.sin(t);
      2:  vers := 1.0 + system.cos(t);
    else  vers := 1.0 - system.sin(t)
  end;
end;


{---------------------------------------------------------------------------}
function covers(x: extended): extended;
  {-Return the coversine covers(x) = 1 - sin(x)}
var
  t: extended;
begin
  {reduction mod Pi/2, |t| <= Pi/4}
  case rem_pio2(x,t) and 3 of
      0:  covers := 1.0 - system.sin(t);
      1:  begin
            {covers(x) = 1 - cos(t)}
            x := t*t;
            t := PolEvalX(x,pcv,npcv);
            covers := 0.5*x + sqr(x)*t;
          end;
      2:  covers := 1.0 + system.sin(t);
    else  covers := 1.0 + system.cos(t);
  end;
end;


{---------------------------------------------------------------------------}
function hav(x: extended): extended;
  {-Return the haversine hav(x) = 0.5*(1 - cos(x))}
begin
  hav := 0.5*vers(x);
end;


{---------------------------------------------------------------------------}
function cosPi(x: extended): extended;
  {-Return cos(Pi*x), result will be 1 for abs(x) >= 2^64}
var
  t: extended;
  i: integer;
begin
  i := rem_int2(x,t) and 3;
  t := Pi*t;
  case i of
      0:  cosPi :=  system.cos(t);
      1:  cosPi := -system.sin(t);
      2:  cosPi := -system.cos(t);
    else  cosPi :=  system.sin(t);
  end;
end;


{---------------------------------------------------------------------------}
function sin(x: extended): extended;
  {-Accurate version of circular sine, uses system.sin for |x| <= Pi/4}
var
  t: extended;
begin
  {reduction mod Pi/2, |t| <= Pi/4}
  case rem_pio2(x,t) and 3 of
      0:  sin :=  system.sin(t);
      1:  sin :=  system.cos(t);
      2:  sin := -system.sin(t);
    else  sin := -system.cos(t);
  end;
end;


{---------------------------------------------------------------------------}
function sinPi(x: extended): extended;
  {-Return sin(Pi*x), result will be 0 for abs(x) >= 2^64}
var
  t: extended;
  i: integer;
begin
  i := rem_int2(x,t) and 3;
  t := Pi*t;
  case i of
      0:  sinPi :=  system.sin(t);
      1:  sinPi :=  system.cos(t);
      2:  sinPi := -system.sin(t);
    else  sinPi := -system.cos(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure sincos(x: extended; var s,c: extended);
  {-Return accurate values s=sin(x), c=cos(x)}
var
  t,ss,cc: extended;
  n: integer;
begin
  {reduction mod Pi/2, |t| <= Pi/4}
  n := rem_pio2(x,t) and 3;
  _sincos(t,ss,cc);
  case n of
      0: begin s:= ss; c:= cc; end;
      1: begin s:= cc; c:=-ss; end;
      2: begin s:=-ss; c:=-cc; end;
    else begin s:=-cc; c:= ss; end;
  end;
end;


{---------------------------------------------------------------------------}
procedure sincosPi(x: extended; var s,c: extended);
  {-Return s=sin(Pi*x), c=cos(Pi*x); (s,c)=(0,1) for abs(x) >= 2^64}
var
  t,ss,cc: extended;
  n: integer;
begin
  n := rem_int2(x,t) and 3;
  t := Pi*t;
  _sincos(t,ss,cc);
  case n of
      0: begin s:= ss; c:= cc; end;
      1: begin s:= cc; c:=-ss; end;
      2: begin s:=-ss; c:=-cc; end;
    else begin s:=-cc; c:= ss; end;
  end;
end;


{---------------------------------------------------------------------------}
function tan(x: extended): extended;
  {-Return the circular tangent of x, x mod Pi <> Pi/2}
var
  t: extended;
begin
  {reduction mod Pi/2, |t| <= Pi/4}
  if odd(rem_pio2(x,t)) then tan := -_cot(t)
  else tan := _tan(t);
end;


{---------------------------------------------------------------------------}
function cot(x: extended): extended;
  {-Return the circular cotangent of x, x mod Pi <> 0}
var
  t: extended;
begin
  {reduction mod Pi/2, |t| <= Pi/4}
  if odd(rem_pio2(x,t)) then cot := -_tan(t)
  else cot := _cot(t);
end;


{---------------------------------------------------------------------------}
function gd(x: extended): extended;
  {-Return the Gudermannian function gd(x)}
var
  z: extended;
begin
  {gd(x) = arctan(sinh(x)) = 2*arctan(tanh(x/2))}
  z := abs(x);
  if z < 46.0 then begin
    if z < 7.5e-3 then begin
      {gd(x) = x - 1/6*x^3 + 1/24*x^5 - 61/5040*x^7 + 277/72576*x^9 + O(x^10)}
      if z < 3e-5 then gd := x*(1.0-sqr(x)/SIXX)
      else begin
        z  := x*x;
        gd := x*(1.0 - z/SIXX*(1.0 - 0.25*z*(1.0 - 61.0*z/210.0)));
      end;
    end
    else gd := arctan(sinh(x));
  end
  else gd := copysign(Pi_2,x);
end;


{---------------------------------------------------------------------------}
function versint(x: extended): extended;
  {-Return versint(x) = integral(vers(t),t=0..x) = x - sin(x), accurate near 0}
var
  t: extended;
const
  CSN = 9;
  CSVH: array[0..CSN-1] of THexExtW = (     {chebyshev((x-sin(x))/(x^3), x=-4/3..4/3, 0.5e-20)}
          ($6B6F,$3E6D,$40D7,$A351,$3FFD),  {+0.318979288363346959407456566900    } {*2}
          ($BC29,$6F15,$A397,$E8AF,$BFF7),  {-0.710101592891792352000298121574e-2 }
          ($48A4,$7E23,$2DE8,$9E69,$3FF1),  {+0.755361827626963385962654365195e-4 }
          ($0FA8,$3E4D,$B9B5,$FB81,$BFE9),  {-0.468467809127500007796911607370e-6 }
          ($7040,$DD23,$3DD3,$8292,$3FE2),  {+0.190006184732484199713774118361e-8 }
          ($BA5C,$0D4A,$8FE6,$BF0D,$BFD9),  {-0.543005219722655848794099416504e-11}
          ($1655,$9E01,$2079,$CF8C,$3FD0),  {+0.115211934731518499595222170191e-13}
          ($7229,$34BC,$4BD3,$ADFF,$BFC7),  {-0.188648197267528207115909378480e-16}
          ($2A3B,$EA7E,$624F,$E879,$3FBD)); {+0.246141587293140761985659464726e-19}
  fs = 1.125;  {2*(3/4)^2}
  t0 = 4/3;
var
  CSV: array[0..CSN-1] of extended absolute CSVH;
begin
  t := abs(x);
  if t < t0 then begin
    if t <= sqrt_epsh then begin
      {versint(x) = x^3/6*(1 - 1/29*x^2 + O(x^4))}
      versint := x*x*x/SIXX;
    end
    else begin
      t := CSEvalX(fs*sqr(t) - 1.0, CSV, CSN);
      versint := x*x*x*t;
    end
  end
  else versint := x - sin(x);
end;


{---------------------------------------------------------------------------}
function arcgd(x: extended): extended;
  {-Return the inverse Gudermannian function arcgd(x), |x| < Pi/2}
var
  z,s,c: extended;
begin
  {arcgd(x) := 2*arctanh(tan(x/2)) = arcsinh(tan(x)) = ln(sec(x)+tan(x))}
  z := abs(x);
  if z < sqrt_epsh then arcgd := x
  else if z<=0.85 then arcgd := arcsinh(tan(x))
  else begin
    {arcgd := ln(sec(x)+tan(x)) = ln(1/cos(x)+sin(x)/cos(x))}
    sincos(z,s,c);
    {cos(x) is > 0 for |x| < Pi/2, if c <= 0 then x >= Pi/2}
    if c<=0.0 then z := PosInf_x
    else z := ln((1.0+s)/c);
    if x>0.0 then arcgd := z
    else arcgd := -z;
  end;
end;


{---------------------------------------------------------------------------}
function logcf(x,i,d: extended): extended;
  {-Calculate 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3d) .. via continued fractions}
var
  c1,c2,c3,c4,a1,b1,a2,b2: extended;
const
  hugeh: THexExtW = ($0000,$0000,$0000,$8000,$43FF); {2^1024}
  tinyh: THexExtW = ($0000,$0000,$0000,$8000,$3BFF); {1/2^1024}
var
  huge: extended absolute hugeh;
  tiny: extended absolute tinyh;
begin
  {Pascal translation of function logcf by Ian Smith [25]}
  c1 := 2.0*d;
  c2 := i + d;
  c4 := c2 + d;
  a1 := c2;
  b1 := i*(c2 - i*x);
  b2 := d*d*x;
  a2 := c4*c2 - b2;
  b2 := c4*b1 - i*b2;

  while abs(a2*b1-a1*b2) > abs(eps_x*b1*b2)  do begin
    c3 := c2*c2*x;
    c2 := c2 + d;
    c4 := c4 + d;
    a1 := c4*a2 - c3*a1;
    b1 := c4*b2 - c3*b1;

    c3 := c1*c1*x;
    c1 := c1 + d;
    c4 := c4 + d;
    a2 := c4*a1 - c3*a2;
    b2 := c4*b1 - c3*b2;

    c3 := abs(b2);
    if c3 > huge then begin
      a1 := a1*tiny;
      b1 := b1*tiny;
      a2 := a2*tiny;
      b2 := b2*tiny;
    end
    else if c3 < tiny then begin
      a1 := a1*huge;
      b1 := b1*huge;
      a2 := a2*huge;
      b2 := b2*huge;
    end;
  end;

  logcf := a2 / b2;
end;


{---------------------------------------------------------------------------}
function ln1pmx(x: extended): extended;
  {-Return ln(1+x)-x, x>-1, accurate even for -0.5 <= x <= 1.0}
const
  c9 = 2/9;  {FIX311}
  c7 = 2/7;
  c5 = 2/5;
  c3 = 2/3;
  xm = -0.525;
  xp =  1.05;
  x0 =  0.016;
var
  y,z: extended;
begin
  {Based on function log1 by Ian Smith [25]}
  if (x<xm) or (x>xp) then ln1pmx := ln(1.0+x)-x
  else begin
    z := x/(2.0 + x);
    y := z*z;
    if abs(x)>x0 then y := 2.0*y*logcf(y, 3.0, 2.0)
    else y := (((c9*y + c7)*y + c5)*y + c3)*y;
    ln1pmx := z*(y-x);
  end;
end;


{---------------------------------------------------------------------------}
function lnxp1(x: extended): extended;
  {-Delphi alias for ln1p}
begin
  lnxp1 := ln1p(x);
end;


{---------------------------------------------------------------------------}
function ln1mexp(x: extended): extended;
  {-Return ln(1-exp(x)), x<0}
const
  x0 = -44.5;  { < ln(0.5_eps), double: -36.8}
  x1 = -0.693147180559945309; {-ln(2)}
begin
  if x <= x0 then begin
    if x < ln_succx0 then ln1mexp := 0.0
    else ln1mexp := -exp(x);
  end
  else if x < x1 then ln1mexp := ln1p(-exp(x))
  else if x > -eps_x then ln1mexp := ln(-x)
  else ln1mexp := ln(-expm1(x));
end;


{---------------------------------------------------------------------------}
function ln1pexp(x: extended): extended;
  {-Accurately compute ln(1+exp(x)) without overflow}
const
  x0 = -44.5;  { < ln(0.5_eps), double = -36.8}
  x1 = 21.9;   { x + exp(-x) - ln(1+exp(x) < eps_x/2, double = 18.1}
  x2 = 40.7;   { exp(-x)/x < eps_x/2, double = 33.3}
begin
  if x <= x0 then begin
    if x < ln_succx0 then ln1pexp := 0.0
    else ln1pexp := exp(x);
  end
  else if x <= x1 then ln1pexp := ln1p(exp(x))
  else if x <= x2 then ln1pexp := x + exp(-x)
  else ln1pexp := x;
end;


{---------------------------------------------------------------------------}
function log2p1(x: extended): extended;
  {-Return log2(1+x), accurate even for x near 0}
begin
  if abs(x) > 0.5 then log2p1 := log2(1.0 + x)
  else log2p1 := log2e*ln1p(x);
end;


{---------------------------------------------------------------------------}
function log10p1(x: extended): extended;
  {-Return log10(1+x), accurate even for x near 0}
begin
  if abs(x) > 0.5 then log10p1 := log10(1.0 + x)
  else log10p1 := log10e*ln1p(x);
end;


{---------------------------------------------------------------------------}
function logaddexp(x,y: extended): extended;
  {-Accurately compute ln[exp(x) + exp(y)]}
begin
  logaddexp := x + ln1pexp(y - x)
end;


{---------------------------------------------------------------------------}
function logsubexp(x,y: extended): extended;
  {-Accurately compute ln[exp(x) - exp(y)], x > y}
begin
  logsubexp := x + ln1mexp(y - x)
end;


{---------------------------------------------------------------------------}
function logistic(x: extended): extended;
  {-Return logistic(x) = 1/(1+exp(-x))}
const
  x0 = 46.0;
  x1 = -23.0;
begin
  if abs(x) > x0 then begin
    if x>0.0 then logistic := 1.0
    else logistic := exp(x);  {x negative and 1+exp(-x) ~ exp(-x)}
  end
  else if x < x1 then begin
    {exp(x) is small: logistic = exp(x) + exp(x)^2 - exp(x)^3 + ...}
    {This is slightly more accurate than the default for x near -x0}
    x := exp(x);
    logistic := x - sqr(x);
  end
  else logistic := 1.0/(1.0+exp(-x));
end;


{---------------------------------------------------------------------------}
function logit(x: extended): extended;
  {-Return logit(x) = ln(x/(1.0-x)), accurate near x=0.5}
const
  xl = 1/3;  {logit(1/3) = -ln2}  {FIX311}
  xh = 2/3;  {logit(2/3) =  ln2}
begin
  if (x < xl) or (x > xh) then logit := ln(x/(1.0-x))
  else logit := ln1p((2.0*x-1.0)/(1.0-x));
end;


{---------------------------------------------------------------------------}
function sinc(x: extended): extended;
  {-Return the cardinal sine sinc(x) = sin(x)/x}
begin
  {sin(x)/x = 1 - 1/6*x^2 + 1/120*x^4 - 1/5040*x^6 + 1/362880*x^8 + O(x^10)}
  x := abs(x);
  if x < 6.15e-4 then begin
    if x < 2.3e-10 then sinc := 1.0
    else begin
      x := sqr(x);
      if x < 2.3e-10 then sinc := 1.0 - x/SIXX
      else sinc := 1.0 - (x/SIXX - sqr(x)/120.0)
    end;
  end
  else sinc := sin(x)/x;
end;


{---------------------------------------------------------------------------}
function sincPi(x: extended): extended;
  {-Return the normalised cardinal sine sincPi(x) = sin(Pi*x)/(Pi*x)}
var
  y: extended;
begin
  x := abs(x);
  y := Pi*x;
  if y < 6.15e-4 then sincPi := sinc(y)
  else sincPi := sinPi(x)/y;
end;


{---------------------------------------------------------------------------}
function sinh_small(x: extended; rel: boolean): extended;
  {-Internal: return sinh(x) for |x| <= 1}
const
  ncs=9;
  csh: array[0..ncs-1] of THexExtW = (     {chebyshev(sinh(x)/x-1, x=-1..1, 0.1e-20)}
         ($D006,$0AFA,$F911,$B131,$3FFC),  {+0.173042194047179631675883846985    }
         ($4EEC,$D088,$9973,$B364,$3FFB),  {+0.875942219227604771549002634544e-1 }
         ($AECD,$1FE8,$437A,$8D7D,$3FF5),  {+0.107947777456713275024272706516e-2 }
         ($D1B0,$68E6,$89C6,$D5E7,$3FED),  {+0.637484926075475048156855545659e-5 }
         ($DAE7,$9306,$8CA6,$BD2E,$3FE5),  {+0.220236640492305301591904955350e-7 }
         ($E0C5,$862E,$36BE,$DB5F,$3FDC),  {+0.498794018041584931494262802066e-10}
         ($5D46,$9B41,$8645,$B389,$3FD3),  {+0.797305355411573048133738827571e-13}
         ($D7DD,$B169,$A8B3,$DA6F,$3FC9),  {+0.947315871307254445124531745434e-16}
         ($7636,$F264,$EF5C,$CD44,$3FBF)); {+0.869349205044812416919840906342e-19}
var
  css: array[0..ncs-1] of extended absolute csh;
var
  t: extended;
const
  t1 = 2.3E-10;  {~ sqrt(2^-64)}
begin
  if abs(x)<=t1 then begin
    {sinh(x) = x*(1 + 1/6*x^2 + 1/120*x^4 + O(x^6))}
    if rel then sinh_small := 1.0
    else sinh_small := x;
  end
  else begin
    t := CSEvalX(2.0*sqr(x)-1.0, css, ncs);
    if rel then sinh_small := 1.0 + t
    else sinh_small := x + x*t;
  end;
end;


{---------------------------------------------------------------------------}
function sinh(x: extended): extended;
  {-Return the hyperbolic sine of x, accurate even for x near 0}
var
  t: extended;
const
  t0 = 23.0;     {ceil(-0.5*ln(2^-64)) = ceil(32*ln(2))}
begin
  t := abs(x);
  if t<1.0 then sinh := sinh_small(x, false)
  else begin
    {sinh(t) = 0.5*(exp(t) - exp(-t))}
    if t<=t0 then begin
      {calculate inverse only if it is not too small}
      t := exp(t);
      t := 0.5*(t - 1.0/t);
    end
    else if t<ln_MaxExt then begin
      {t>t0: exp(t)-exp(-t) ~ exp(t)}
      t := 0.5*exp(t)
    end
    else begin
      {exp(t) would overflow, use small margin below overflow}
      t := exp(t-ln2);
    end;
    if x>0.0 then sinh := t
    else sinh := -t;
  end;
end;


{---------------------------------------------------------------------------}
function lnsinh(x: extended): extended;
  {-Return ln(sinh(x)), x > 0, accurate for x ~ 0 and without overflow for large x}
const
  xl = 22.0;
var
  t: extended;
begin
  if x <= sqrt_epsh then begin
    {ln(sinh(x)) = ln(x) + x^2/6 - x^4/180 + O(x^6)}
    {Additionally this generates an error if x <= 0}
    lnsinh := ln(x) + sqr(x)/SIXX;
  end
  else if x<1.0 then begin
    t := sinh_small(x, false);
    lnsinh := ln(t);
  end
  else if x<xl  then begin
    t := 1.0 - exp(-2.0*x);
    lnsinh := x + ln(0.5*t);
  end
  else lnsinh := x - ln2;
end;


{---------------------------------------------------------------------------}
function sinhc(x: extended): extended;
  {-Return sinh(x)/x, accurate even for x near 0}
var
  t: extended;
begin
  t := abs(x);
  if abs(t) < 1.0 then sinhc := sinh_small(t, true)
  else if t < ln_MaxExt then sinhc := sinh(t)/t
  else begin
    t := t - ln(2.0*t);
    sinhc := exp(t);
  end;
end;


{---------------------------------------------------------------------------}
function sinhmx(x: extended): extended;
  {-Return sinh(x)-x, accurate even for x near 0}
var
  t: extended;
const
  ncs=10;
  csh: array[0..ncs-1] of THexExtW = (  {chebyshev((sinh(x)-x)/x^3, x=-2..2, 0.1e-20)}
         ($2BC6,$38FA,$BA5E,$BD02,$3FFD),  {+0.369161437989977107646200496798    } {*2}
         ($3E29,$DF18,$BFFE,$963C,$3FF9),  {+0.183395147241496573111025588024e-1 }
         ($0773,$C2B0,$FE21,$E224,$3FF3),  {+0.431336408128174581497124393209e-3 }
         ($D47F,$FFCC,$571C,$C6E1,$3FED),  {+0.592709289470722070924036619944e-5 }
         ($47C2,$D923,$DC79,$E56E,$3FE6),  {+0.534190451019238924481986729972e-7 }
         ($4D3E,$B7C5,$7F8F,$BAF2,$3FDF),  {+0.340055083020092921361501610500e-9 }
         ($CBBC,$75F0,$15DD,$E29C,$3FD7),  {+0.161015882323068379168404355934e-11}
         ($89E2,$978C,$9983,$D448,$3FCF),  {+0.589205330187637306105897428969e-14}
         ($51A1,$B57D,$C8A7,$9E47,$3FC7),  {+0.171607959509371641280224879040e-16}
         ($0A5E,$F081,$7613,$C04F,$3FBE)); {+0.407233102668331851851851851852e-19}
var
  csv: array[0..ncs-1] of extended absolute csh;
begin
  t := abs(x);
  if t < 2.0 then begin
    if t <= sqrt_epsh then begin
      {sinh(x)-x = x^3/6*(1 + 1/20*x^2 + O(x^4))}
      sinhmx := x*x*x/SIXX;
    end
    else begin
      t := CSEvalX(0.5*sqr(t) - 1.0, csv, ncs);
      sinhmx := x*x*x*t;
    end
  end
  else sinhmx := sinh(x)-x;
end;


{---------------------------------------------------------------------------}
function cosh(x: extended): extended;
  {-Return the hyperbolic cosine of x}
const
  x0 = 23.0;     {ceil(-0.5*ln(2^-64)) = ceil(32*ln(2))}
begin
  x := abs(x);
  if x<=1.0 then cosh := 1.0 + coshm1(x)
  else begin
    if x>x0 then begin
      {x>x0: exp(x) + exp(-x) ~ exp(x)}
      if x<ln_MaxExt then cosh := 0.5*exp(x)
      else cosh := exp(x-ln2);
    end
    else begin
      {calculate inverse only if it is not too small}
      x := exp(x);
      cosh := 0.5*(x + 1.0/x);
    end
  end;
end;


{---------------------------------------------------------------------------}
function lncosh(x: extended): extended;
  {-Return ln(cosh(x)), accurate for x ~ 0 and without overflow for large x}
const
  xl = 23.0;
begin
  x := abs(x);
  if x <= sqrt_epsh then begin
    {ln(cosh(x)) = x^2/2 - x^4/12 + O(x^6)}
    lncosh := 0.5*sqr(x);
  end
  else if x < 1.0 then lncosh := ln1p(coshm1(x))
  else if x >= xl then begin
    {ln(cosh(x)) = x - ln(2) + 1/(exp(x)^2) + O(1/(exp(x)^4))}
    lncosh := x - ln2;
  end
  else begin
    lncosh := x + ln(0.5*(1.0 + exp(-2.0*x)));
  end;
end;


{---------------------------------------------------------------------------}
function coshm1(x: extended): extended;
  {-Return cosh(x)-1, accurate even for x near 0}
const
  ncc=8;
  cch: array[0..ncc-1] of THexExtW = (     {chebyshev(((cosh(x)-1)/(x^2)-1/2)/x^2, x=-1..1, 0.1e-20);}
         ($7F5B,$DC92,$B00E,$AD8C,$3FFB),  {+0.847409967933085972383472682702e-1 }
         ($E233,$6D68,$4FB4,$B954,$3FF4),  {+0.706975331110557282636312722683e-3 }
         ($C401,$57F2,$9421,$D38C,$3FEC),  {+0.315232776534872632694769980081e-5 }
         ($60F4,$9F76,$CDC0,$9634,$3FE4),  {+0.874315531301413118696907944019e-8 }
         ($D666,$F7AB,$C086,$9172,$3FDB),  {+0.165355516210397415961496829471e-10}
         ($5E0D,$C9ED,$687C,$CC55,$3FD1),  {+0.226855895848535933261822520050e-13}
         ($0DEF,$8B25,$7570,$D9B9,$3FC7),  {+0.236057319781153495473072617928e-16}
         ($7B04,$4A7F,$2A3E,$B5FC,$3FBD)); {+0.192684134365813759102742156829e-19}
var
  ccs: array[0..ncc-1] of extended absolute cch;
  z: extended;
begin
  x := abs(x);
  if x<=1.0 then begin
    x := sqr(x);
    if x<eps_x then begin
      {cosh(x)-1 = 0.5*x^2*[1 + 1/12*x^2 + O(x^4)]}
      coshm1 := 0.5*x;
    end
    else begin
      z := CSEvalX(2.0*x-1.0, ccs, ncc);
      coshm1 := 0.5*x + x*x*z;
    end;
  end
  else coshm1 := cosh(x)-1.0;
end;

(*
function sinh(x: extended): extended;
var
  s,c: extended;
begin
  sinhcosh(x,s,c);
  sinh := s;
end;


function cosh(x: extended): extended;
var
  s,c: extended;
begin
  sinhcosh(x,s,c);
  cosh := c;
end;
*)


{---------------------------------------------------------------------------}
procedure sinhcosh(x: extended; var s,c: extended);
  {-Return s=sinh(x) and c=cosh(x)}
var
  t: extended;
const
  t0 = 23.0;     {ceil(-0.5*ln(2^-64)) = ceil(32*ln(2))}
begin
  t := abs(x);
  if t<=1.0 then begin
    s := sinh_small(x, false);
    c := 1.0 + coshm1(x)
  end
  else begin
    if t<=t0 then begin
      {calculate inverse only if it is not too small}
      t := exp(t);
      c := 1.0/t;
      s := 0.5*(t - c);
      c := 0.5*(t + c);
    end
    else if t<ln_MaxExt then begin
      {t>t0: exp(t) + exp(-t) ~ exp(t) - exp(-t) ~ exp(t)}
      s := 0.5*exp(t);
      c := s;
    end
    else begin
      {exp(t) would overflow, use small margin below overflow}
      s := exp(t-ln2);
      c := s;
    end;
    if x<0.0 then s := -s;
  end;
end;


{---------------------------------------------------------------------------}
function tanh(x: extended): extended;
  {-Return the hyperbolic tangent of x, accurate even for x near 0}
var
  t,z: extended;
const
  t0 = 23.0;       {ceil(-0.5*ln(2^-64)) = ceil(32*ln(2))}
  t1 = 1.5e-5;     {~ 2^(-64/4)}
  t2 = 0.625;
  PH: array[0..3] of THexExtW = (
        ($e3be,$bfbd,$5cbc,$a381,$c009),  {-1.3080425704712825945553e+3}
        ($b576,$ef5e,$6d57,$a81b,$c005),  {-8.4053568599672284488465e+1}
        ($5959,$9111,$9cc7,$f4e2,$bffe),  {-9.5658283111794641589011e-1}
        ($d2a4,$1b0c,$8f15,$8f99,$bff1)); {-6.8473739392677100872869e-5}
  QH: array[0..3] of THexExtW = (
        ($d5a2,$1f9c,$0b1b,$f542,$400a),  { 3.9241277114138477845780e+3}
        ($3793,$c95f,$fa2f,$e3b9,$4009),  { 1.8218117903645559060232e+3}
        ($687f,$ce24,$dd6c,$c084,$4005),  { 9.6259501838840336946872e+1}
        ($0000,$0000,$0000,$8000,$3fff)); { 1.0000000000000000000000e+0}
var
  P: array[0..3] of extended absolute PH;
  Q: array[0..3] of extended absolute QH;
begin
  t := abs(x);
  if t>t0 then begin
    if x<0.0 then tanh := -1.0 else tanh := +1.0;
  end
  else if t<t2 then begin
    z := x*x;
    if t <= t1 then begin
      {tanh(x) = x*(1 - 1/3*x^2 + 2/15*x^4 + O(x^6))}
      t := -0.33333333333333333333;
    end
    else begin
      {Ref: Cephes[7], tanhl.c}
      t := PolEvalX(z,P,4)/PolEvalX(z,Q,4);
    end;
    z := x*t*z;
    tanh := x + z;
  end
  else begin
    t := exp(2.0*t)+1.0;
    if x>0.0 then tanh := 1.0 - 2.0/t
    else tanh := -1.0 + 2.0/t;
  end;
end;


{---------------------------------------------------------------------------}
function intpower(x: extended; n: longint): extended;
  {-Return x^n; via binary exponentiation (no overflow detection)}
var
  i: longint;
  r: extended;
begin
  i := abs(n);
  r := abs(x);
  if (i=0) or (r=1.0) then begin
    if odd(i) and (x < 0) then intpower := -1
    else intpower := 1.0;
    exit;
  end;
  if (n<0) and (r>1.0) then begin
    if (1.0+ilogb(r))*i >= 16383 then begin
      {avoid spurious overflow, use x^(-n) = (1/x)^n}
      {so that e.g. 2^(-16400) is computed as denormal}
      x := 1.0/x;
      n := i;
    end;
  end;
  {Use internal assembler function}
  intpower := intpower0(x,n);
end;


{The default power function uses the version from pow_512.inc}
{.$i pow_032.inc}
{.$i pow_064.inc}
{.$i pow_128.inc}
{.$i pow_256.inc}
{.$i pow_512.inc}
{.$i pow_1024.inc}

{Power functions space / accuracy figures, eps1/2: peak relative errors}
{measured in eps_x (about 1e-19) for |x|,|y| < 1000 and |x|,|y| < 2000.}

{TabSize    Bytes    eps1 / eps2 }
{     32      500    25.6 / 42.9 }  {used in old AMath versions}
{     64      980    13.0 / 22.3 }
{    128     1940     7.1 / 12.1 }
{    256     3860     4.6 /  6.5 }
{    512     7700     1.9 /  3.1 }  {this version}
{   1024    15380     1.4 /  1.7 }


{---------------------------------------------------------------------------}
function power(x, y : extended): extended;
  {-Return x^y; if frac(y)<>0 then x must be > 0}

  {This is my Pascal translation of the ldouble/powl.c file from the   }
  {freely distributable Cephes Mathematical Library V2.8, (June 2000)  }
  {Copyright 1984, 1991, 1998 by Stephen L. Moshier                    }

  {Following Cody and Waite, the Cephes function uses a lookup table of}
  {2**-i/NXT and pseudo extended precision arithmetic to obtain several}
  {extra bits of accuracy in both the logarithm and the exponential.   }

  {My routine uses a table size of 512 entries resulting in about four }
  {additional bits of accuracy; the tables are calculated with MPArith.}

const
  NXT = 512;  {table size, must be power of 2}

const
  MEXP   =  NXT*16384.0;
  MNEXP  = -NXT*(16384.0+64.0);    {+64 for denormal support}

const
  LOG2EA : THexExtW = ($c2ef,$705f,$eca5,$e2a8,$3ffd); { = 0.44269504088896340736 = log2(e)-1}

const
  {Coefficients for log(1+x) =  x - .5x^2 + x^3 *  P(z)/Q(z) }
  {on the domain  2^(-1/32) - 1  <=  x  <=  2^(1/32) - 1     }
  PH: array[0..3] of THexExtW = (
        ($b804,$a8b7,$c6f4,$da6a,$3ff4),
        ($7de9,$cf02,$58c0,$fae1,$3ffd),
        ($405a,$3722,$67c9,$e000,$3fff),
        ($cd99,$6b43,$87ca,$b333,$3fff));
  QH: array[0..2] of THexExtW = (
        ($6307,$a469,$3b33,$a800,$4001),
        ($fec2,$62d7,$a51c,$8666,$4002),
        ($da32,$d072,$a5d7,$8666,$4001));
  {Coefficients for 2^x = 1 + x P(x),  on the interval -1/32 <= x <= 0}
  RH: array[0..6] of THexExtW = (
        ($a69b,$530e,$ee1d,$fd2a,$3fee),
        ($c746,$8e7e,$5960,$a182,$3ff2),
        ($63b6,$adda,$fd6a,$aec3,$3ff5),
        ($c104,$fd99,$5b7c,$9d95,$3ff8),
        ($e05e,$249d,$46b8,$e358,$3ffa),
        ($5d1d,$162c,$effc,$f5fd,$3ffc),
        ($79aa,$d1cf,$17f7,$b172,$3ffe));
var
  P: array[0..3] of extended absolute PH;
  Q: array[0..2] of extended absolute QH;
  R: array[0..6] of extended absolute RH;

const
  {A[i] = 2^(-i/NXT), calculated with MPArith/t_powtab.pas}
  {If i is even, A[i] + B[i/2] gives additional accuracy. }
  A: array[0..NXT] of THexExtW = (
       ($0000,$0000,$0000,$8000,$3fff), ($aed2,$1c8d,$5652,$ffa7,$3ffe), ($c8a5,$511e,$cb59,$ff4e,$3ffe),
       ($a3ca,$fb18,$5f0a,$fef6,$3ffe), ($884c,$7b8f,$115c,$fe9e,$3ffe), ($695d,$3745,$e243,$fe45,$3ffe),
       ($9f35,$96a8,$d1b4,$fded,$3ffe), ($a15e,$05d2,$dfa6,$fd95,$3ffe), ($c175,$f486,$0c0c,$fd3e,$3ffe),
       ($e654,$d630,$56de,$fce6,$3ffe), ($47bb,$21e4,$c011,$fc8e,$3ffe), ($2a57,$525a,$4799,$fc37,$3ffe),
       ($9c49,$e5f0,$ed6c,$fbdf,$3ffe), ($3211,$5ea9,$b181,$fb88,$3ffe), ($c3f4,$4227,$93cc,$fb31,$3ffe),
       ($2bca,$19b1,$9443,$fada,$3ffe), ($033a,$722a,$b2db,$fa83,$3ffe), ($6270,$dc15,$ef8a,$fa2c,$3ffe),
       ($9f35,$eb93,$4a46,$f9d6,$3ffe), ($0c81,$3861,$c305,$f97f,$3ffe), ($ba74,$5dd4,$59bb,$f929,$3ffe),
       ($36c6,$fadf,$0e5e,$f8d3,$3ffe), ($4d9c,$b209,$e0e5,$f87c,$3ffe), ($cad3,$2972,$d145,$f826,$3ffe),
       ($3bb9,$0ad1,$df73,$f7d0,$3ffe), ($b12e,$036e,$0b65,$f77b,$3ffe), ($8239,$c428,$5510,$f725,$3ffe),
       ($0f0b,$016e,$bc6c,$f6cf,$3ffe), ($846e,$733f,$416c,$f67a,$3ffe), ($9f9c,$d52c,$e407,$f624,$3ffe),
       ($7290,$e653,$a433,$f5cf,$3ffe), ($28bd,$695f,$81e6,$f57a,$3ffe), ($cc2c,$2486,$7d15,$f525,$3ffe),
       ($0b17,$e18c,$95b5,$f4d0,$3ffe), ($fddf,$6db9,$cbbe,$f47b,$3ffe), ($ed7e,$99e3,$1f24,$f427,$3ffe),
       ($1a5b,$3a64,$8fde,$f3d2,$3ffe), ($8390,$271a,$1de1,$f37e,$3ffe), ($ae9c,$3b6b,$c923,$f329,$3ffe),
       ($6f7d,$563f,$919a,$f2d5,$3ffe), ($b13a,$59ff,$773c,$f281,$3ffe), ($3eda,$2c97,$79ff,$f22d,$3ffe),
       ($8cc1,$b770,$99d8,$f1d9,$3ffe), ($8280,$e774,$d6be,$f185,$3ffe), ($4509,$ad09,$30a7,$f132,$3ffe),
       ($0154,$fc11,$a788,$f0de,$3ffe), ($b76a,$cbe8,$3b58,$f08b,$3ffe), ($05e5,$1767,$ec0d,$f037,$3ffe),
       ($f5cb,$dcda,$b99b,$efe4,$3ffe), ($c6e3,$1e0a,$a3fb,$ef91,$3ffe), ($bc6b,$e032,$ab20,$ef3e,$3ffe),
       ($ea3d,$2c03,$cf03,$eeeb,$3ffe), ($025b,$0da3,$0f98,$ee99,$3ffe), ($22ea,$94a7,$6cd5,$ee46,$3ffe),
       ($a491,$d418,$e6b1,$edf3,$3ffe), ($e946,$e26f,$7d22,$eda1,$3ffe), ($2b84,$d994,$301e,$ed4f,$3ffe),
       ($4dea,$d6da,$ff9b,$ecfc,$3ffe), ($ab41,$fb03,$eb8f,$ecaa,$3ffe), ($e6f1,$6a3c,$f3f1,$ec58,$3ffe),
       ($bddc,$4c1c,$18b6,$ec07,$3ffe), ($d79f,$cba2,$59d4,$ebb5,$3ffe), ($9840,$1736,$b743,$eb63,$3ffe),
       ($f244,$60a5,$30f7,$eb12,$3ffe), ($392f,$dd24,$c6e7,$eac0,$3ffe), ($f464,$c548,$790a,$ea6f,$3ffe),
       ($b27b,$550e,$4756,$ea1e,$3ffe), ($dcf5,$cbd1,$31c0,$e9cd,$3ffe), ($8c57,$6c4f,$3840,$e97c,$3ffe),
       ($5cb8,$7ca4,$5acb,$e92b,$3ffe), ($42ab,$464b,$9958,$e8da,$3ffe), ($6094,$161c,$f3dd,$e889,$3ffe),
       ($dc68,$3c4b,$6a50,$e839,$3ffe), ($b5d3,$0c68,$fca8,$e7e8,$3ffe), ($9cbf,$dd5b,$aada,$e798,$3ffe),
       ($c84b,$0965,$74df,$e748,$3ffe), ($ce22,$ee1f,$5aaa,$e6f8,$3ffe), ($7a41,$ec78,$5c34,$e6a8,$3ffe),
       ($a717,$68b3,$7973,$e658,$3ffe), ($1616,$ca69,$b25c,$e608,$3ffe), ($48a8,$7c83,$06e7,$e5b9,$3ffe),
       ($5988,$ed3e,$7709,$e569,$3ffe), ($d681,$8e26,$02ba,$e51a,$3ffe), ($9a96,$d418,$a9ef,$e4ca,$3ffe),
       ($a88d,$373d,$6ca0,$e47b,$3ffe), ($05e2,$330d,$4ac2,$e42c,$3ffe), ($9619,$4649,$444c,$e3dd,$3ffe),
       ($f681,$f300,$5934,$e38e,$3ffe), ($5a51,$be8a,$8972,$e33f,$3ffe), ($6732,$3185,$d4fc,$e2f0,$3ffe),
       ($1226,$d7d9,$3bc7,$e2a2,$3ffe), ($7cdd,$40b2,$bdcc,$e253,$3ffe), ($d369,$fe83,$5aff,$e205,$3ffe),
       ($2a54,$a703,$1359,$e1b7,$3ffe), ($5d23,$d329,$e6cf,$e168,$3ffe), ($ed37,$1f30,$d559,$e11a,$3ffe),
       ($e111,$2a94,$deec,$e0cc,$3ffe), ($a3fe,$980f,$037f,$e07f,$3ffe), ($e627,$0d99,$430a,$e031,$3ffe),
       ($7d00,$3469,$9d82,$dfe3,$3ffe), ($4420,$b8f0,$12de,$df96,$3ffe), ($fe78,$4ada,$a316,$df48,$3ffe),
       ($37f2,$9d10,$4e1f,$defb,$3ffe), ($276d,$65af,$13f1,$deae,$3ffe), ($9124,$5e0e,$f482,$de60,$3ffe),
       ($a96f,$42bb,$efc9,$de13,$3ffe), ($f7f0,$d378,$05bc,$ddc7,$3ffe), ($3b1a,$d33d,$3653,$dd7a,$3ffe),
       ($4c20,$0832,$8185,$dd2d,$3ffe), ($0346,$3bb4,$e747,$dce0,$3ffe), ($1c92,$3a4f,$6791,$dc94,$3ffe),
       ($1cde,$d3c0,$0259,$dc48,$3ffe), ($3755,$daf2,$b797,$dbfb,$3ffe), ($3343,$25fe,$8742,$dbaf,$3ffe),
       ($5255,$8e29,$714f,$db63,$3ffe), ($3730,$efe4,$75b6,$db17,$3ffe), ($cc72,$2ac9,$946f,$dacb,$3ffe),
       ($2c0c,$219e,$cd6f,$da7f,$3ffe), ($8704,$ba4d,$20ad,$da34,$3ffe), ($0d97,$ddeb,$8e21,$d9e8,$3ffe),
       ($d7b6,$78af,$15c2,$d99d,$3ffe), ($cdeb,$79f9,$b786,$d951,$3ffe), ($929c,$d44a,$7364,$d906,$3ffe),
       ($6bac,$7d46,$4954,$d8bb,$3ffe), ($2c84,$6db3,$394c,$d870,$3ffe), ($2070,$a177,$4343,$d825,$3ffe),
       ($f56a,$1797,$6731,$d7da,$3ffe), ($a737,$d239,$a50b,$d78f,$3ffe), ($6af4,$d69d,$fcca,$d744,$3ffe),
       ($9af2,$2d20,$6e65,$d6fa,$3ffe), ($a2fe,$e13b,$f9d1,$d6af,$3ffe), ($ed01,$0180,$9f08,$d665,$3ffe),
       ($ce07,$9f9b,$5dfe,$d61b,$3ffe), ($739a,$d04f,$36ac,$d5d1,$3ffe), ($d18a,$ab75,$2909,$d587,$3ffe),
       ($9006,$4bfe,$350c,$d53d,$3ffe), ($fa1f,$cfed,$5aab,$d4f3,$3ffe), ($eca5,$585b,$99df,$d4a9,$3ffe),
       ($c561,$0972,$f29e,$d45f,$3ffe), ($52af,$0a6e,$64df,$d416,$3ffe), ($c379,$859a,$f099,$d3cc,$3ffe),
       ($978e,$a853,$95c4,$d383,$3ffe), ($9054,$a302,$5457,$d33a,$3ffe), ($a1de,$a91e,$2c49,$d2f1,$3ffe),
       ($e45a,$f12a,$1d91,$d2a8,$3ffe), ($85e3,$b4b5,$2827,$d25f,$3ffe), ($bcab,$3056,$4c02,$d216,$3ffe),
       ($b983,$a3af,$8918,$d1cd,$3ffe), ($9ac6,$5169,$df62,$d184,$3ffe), ($5f98,$7f34,$4ed6,$d13c,$3ffe),
       ($db8d,$75c5,$d76c,$d0f3,$3ffe), ($aaa0,$80d8,$791b,$d0ab,$3ffe), ($2595,$ef2b,$33da,$d063,$3ffe),
       ($56ab,$127e,$07a2,$d01b,$3ffe), ($eeb5,$3f94,$f468,$cfd2,$3ffe), ($3a86,$ce32,$fa24,$cf8a,$3ffe),
       ($18c1,$1919,$18cf,$cf43,$3ffe), ($f001,$7e0a,$505e,$cefb,$3ffe), ($a55d,$5dc6,$a0ca,$ceb3,$3ffe),
       ($934d,$1c07,$0a0a,$ce6c,$3ffe), ($80e4,$1f84,$8c15,$ce24,$3ffe), ($9967,$d1ee,$26e2,$cddd,$3ffe),
       ($6445,$9ff0,$da6a,$cd95,$3ffe), ($bd65,$f92c,$a6a3,$cd4e,$3ffe), ($cdd2,$503d,$8b86,$cd07,$3ffe),
       ($04bd,$1ab4,$8909,$ccc0,$3ffe), ($10e5,$d115,$9f23,$cc79,$3ffe), ($da4f,$eeda,$cdcd,$cc32,$3ffe),
       ($7c5d,$f272,$14fe,$cbec,$3ffe), ($4041,$5d3b,$74ae,$cba5,$3ffe), ($97c9,$b385,$ecd3,$cb5e,$3ffe),
       ($1886,$7c92,$7d66,$cb18,$3ffe), ($774e,$4290,$265e,$cad2,$3ffe), ($8414,$929e,$e7b2,$ca8b,$3ffe),
       ($2624,$fcc7,$c15a,$ca45,$3ffe), ($58aa,$1401,$b34f,$c9ff,$3ffe), ($27a3,$6e2f,$bd86,$c9b9,$3ffe),
       ($ad1a,$a41c,$dff8,$c973,$3ffe), ($0ecc,$517f,$1a9d,$c92e,$3ffe), ($7c14,$14f3,$6d6c,$c8e8,$3ffe),
       ($2c45,$8ffe,$d85c,$c8a2,$3ffe), ($5d4c,$6709,$5b66,$c85d,$3ffe), ($52b2,$4164,$f681,$c817,$3ffe),
       ($54fd,$c942,$a9a4,$c7d2,$3ffe), ($b15d,$abb9,$74c8,$c78d,$3ffe), ($b9ba,$98c2,$57e4,$c748,$3ffe),
       ($c51e,$4336,$52f0,$c703,$3ffe), ($306b,$60cf,$65e3,$c6be,$3ffe), ($5f79,$aa24,$90b5,$c679,$3ffe),
       ($be7f,$daac,$d35e,$c634,$3ffe), ($c3d9,$b0bb,$2dd6,$c5f0,$3ffe), ($f22e,$ed80,$a014,$c5ab,$3ffe),
       ($dadd,$5506,$2a11,$c567,$3ffe), ($20d3,$ae32,$cbc3,$c522,$3ffe), ($7baa,$c2c0,$8523,$c4de,$3ffe),
       ($bb30,$5f47,$5629,$c49a,$3ffe), ($cb33,$5334,$3ecc,$c456,$3ffe), ($b7b1,$70ca,$3f04,$c412,$3ffe),
       ($b15d,$8d21,$56c9,$c3ce,$3ffe), ($1279,$8026,$8613,$c38a,$3ffe), ($6407,$2497,$ccda,$c346,$3ffe),
       ($6351,$5807,$2b15,$c303,$3ffe), ($07c9,$fad9,$a0bc,$c2bf,$3ffe), ($8940,$f03f,$2dc8,$c27c,$3ffe),
       ($6673,$1e3d,$d231,$c238,$3ffe), ($6be9,$6da3,$8ded,$c1f5,$3ffe), ($bb33,$ca0f,$60f5,$c1b2,$3ffe),
       ($d278,$21ec,$4b42,$c16f,$3ffe), ($9456,$6670,$4cca,$c12c,$3ffe), ($5028,$8b9b,$6586,$c0e9,$3ffe),
       ($ca8d,$8836,$956e,$c0a6,$3ffe), ($4652,$55d5,$dc7a,$c063,$3ffe), ($8db0,$f0d0,$3aa1,$c021,$3ffe),
       ($fbdf,$5848,$afdd,$bfde,$3ffe), ($86f8,$8e24,$3c24,$bf9c,$3ffe), ($ca35,$970d,$df6f,$bf59,$3ffe),
       ($1083,$7a73,$99b6,$bf17,$3ffe), ($5f66,$4285,$6af1,$bed5,$3ffe), ($8238,$fc37,$5317,$be93,$3ffe),
       ($15b4,$b73d,$5222,$be51,$3ffe), ($93e2,$8609,$6809,$be0f,$3ffe), ($6049,$7dcf,$94c4,$bdcd,$3ffe),
       ($d483,$b67e,$d84b,$bd8b,$3ffe), ($4d18,$4ac5,$3297,$bd4a,$3ffe), ($36bf,$580c,$a39f,$bd08,$3ffe),
       ($1bdc,$fe78,$2b5b,$bcc7,$3ffe), ($b269,$60e7,$c9c5,$bc85,$3ffe), ($ea22,$a4f2,$7ed3,$bc44,$3ffe),
       ($fb0d,$f2e9,$4a7e,$bc03,$3ffe), ($7453,$75d4,$2cbf,$bbc2,$3ffe), ($4b6f,$5b70,$258d,$bb81,$3ffe),
       ($ebac,$d430,$34e0,$bb40,$3ffe), ($45fb,$133e,$5ab2,$baff,$3ffe), ($e119,$4e73,$96f9,$babe,$3ffe),
       ($ea09,$be5f,$e9ae,$ba7d,$3ffe), ($44e1,$9e42,$52ca,$ba3d,$3ffe), ($9deb,$2c0b,$d245,$b9fc,$3ffe),
       ($7b17,$a85c,$6816,$b9bc,$3ffe), ($4dbf,$5684,$1437,$b97c,$3ffe), ($84c1,$7c80,$d69f,$b93b,$3ffe),
       ($9ee9,$62fb,$af47,$b8fb,$3ffe), ($3dab,$554c,$9e27,$b8bb,$3ffe), ($3834,$a174,$a337,$b87b,$3ffe),
       ($aec9,$981f,$be70,$b83b,$3ffe), ($1e7c,$8ca4,$efca,$b7fb,$3ffe), ($752f,$d4ff,$373d,$b7bc,$3ffe),
       ($25e9,$c9d7,$94c2,$b77c,$3ffe), ($3d80,$c677,$0851,$b73d,$3ffe), ($7791,$28d1,$91e3,$b6fd,$3ffe),
       ($53cc,$517c,$316f,$b6be,$3ffe), ($2b8f,$a3b2,$e6ee,$b67e,$3ffe), ($47d3,$8550,$b259,$b63f,$3ffe),
       ($f76c,$5ed5,$93a8,$b600,$3ffe), ($a594,$9b63,$8ad3,$b5c1,$3ffe), ($f0d2,$a8b9,$97d3,$b582,$3ffe),
       ($c224,$f738,$baa0,$b543,$3ffe), ($6484,$f9de,$f333,$b504,$3ffe), ($9cbb,$2646,$4185,$b4c6,$3ffe),
       ($c180,$f4a9,$a58c,$b487,$3ffe), ($d3ed,$dfdb,$1f43,$b449,$3ffe), ($9841,$654b,$aea2,$b40a,$3ffe),
       ($aef3,$0501,$53a1,$b3cc,$3ffe), ($ae18,$419f,$0e38,$b38e,$3ffe), ($3b0e,$a05f,$de60,$b34f,$3ffe),
       ($2489,$a911,$c412,$b311,$3ffe), ($7cdd,$e61c,$bf46,$b2d3,$3ffe), ($b4a4,$e47d,$cff5,$b295,$3ffe),
       ($b5ac,$33c5,$f618,$b257,$3ffe), ($fe3b,$6618,$31a6,$b21a,$3ffe), ($bc9f,$102e,$8299,$b1dc,$3ffe),
       ($eb09,$c94f,$e8e8,$b19e,$3ffe), ($6bbd,$2b56,$648e,$b161,$3ffe), ($2590,$d2ac,$f581,$b123,$3ffe),
       ($20b2,$5e4a,$9bbc,$b0e6,$3ffe), ($a3c9,$6fb7,$5736,$b0a9,$3ffe), ($515c,$ab09,$27e8,$b06c,$3ffe),
       ($4584,$b6e0,$0dcb,$b02f,$3ffe), ($33f9,$3c69,$08d8,$aff2,$3ffe), ($8661,$e75b,$1906,$afb5,$3ffe),
       ($7af7,$65f8,$3e50,$af78,$3ffe), ($4375,$690a,$78ad,$af3b,$3ffe), ($2457,$a3e3,$c816,$aefe,$3ffe),
       ($9465,$cc5c,$2c84,$aec2,$3ffe), ($5c8e,$9ad6,$a5f0,$ae85,$3ffe), ($b80e,$ca35,$3452,$ae49,$3ffe),
       ($74e3,$17e4,$d7a4,$ae0c,$3ffe), ($1491,$43d0,$8fdd,$add0,$3ffe), ($ed2f,$1068,$5cf7,$ad94,$3ffe),
       ($4ac6,$42a1,$3eea,$ad58,$3ffe), ($90fb,$a1ec,$35af,$ad1c,$3ffe), ($5d04,$f83e,$413f,$ace0,$3ffe),
       ($a7ed,$1209,$6194,$aca4,$3ffe), ($e929,$be3f,$96a4,$ac68,$3ffe), ($3971,$ce50,$e06a,$ac2c,$3ffe),
       ($75e9,$1626,$3edf,$abf1,$3ffe), ($639b,$6c2a,$b1fa,$abb5,$3ffe), ($d337,$a93e,$39b5,$ab7a,$3ffe),
       ($c526,$a8c0,$d609,$ab3e,$3ffe), ($8de1,$4886,$86ef,$ab03,$3ffe), ($fa98,$68de,$4c5f,$aac8,$3ffe),
       ($7629,$ec90,$2652,$aa8d,$3ffe), ($2e5f,$b8d8,$14c2,$aa52,$3ffe), ($3979,$b569,$17a7,$aa17,$3ffe),
       ($bc03,$cc6b,$2efa,$a9dc,$3ffe), ($0ef8,$ea7c,$5ab4,$a9a1,$3ffe), ($e631,$fea9,$9ace,$a966,$3ffe),
       ($771b,$fa77,$ef41,$a92b,$3ffe), ($9fbe,$d1d8,$5806,$a8f1,$3ffe), ($0e09,$7b32,$d516,$a8b6,$3ffe),
       ($6770,$ef58,$6669,$a87c,$3ffe), ($70d1,$298f,$0bfa,$a842,$3ffe), ($36a1,$2789,$c5c0,$a807,$3ffe),
       ($356a,$e965,$93b4,$a7cd,$3ffe), ($828e,$71af,$75d1,$a793,$3ffe), ($f55b,$c55f,$6c0e,$a759,$3ffe),
       ($5062,$ebd9,$7665,$a71f,$3ffe), ($6b1e,$eee8,$94cf,$a6e5,$3ffe), ($5bdf,$dac3,$c745,$a6ab,$3ffe),
       ($a20c,$be08,$0dc0,$a672,$3ffe), ($509c,$a9be,$6839,$a638,$3ffe), ($38ea,$b151,$d6a9,$a5fe,$3ffe),
       ($15cb,$ea94,$5909,$a5c5,$3ffe), ($b6ee,$6dbe,$ef53,$a58b,$3ffe), ($2c84,$556d,$997f,$a552,$3ffe),
       ($f339,$be9e,$5786,$a519,$3ffe), ($2070,$c8b6,$2962,$a4e0,$3ffe), ($8ec5,$9576,$0f0c,$a4a7,$3ffe),
       ($0ae4,$4905,$087d,$a46e,$3ffe), ($809e,$09e6,$15ae,$a435,$3ffe), ($284c,$00ff,$3698,$a3fc,$3ffe),
       ($b47c,$5991,$6b34,$a3c3,$3ffe), ($7fe0,$413e,$b37c,$a38a,$3ffe), ($bb93,$e802,$0f68,$a352,$3ffe),
       ($9d94,$8037,$7ef3,$a319,$3ffe), ($8f9e,$3e91,$0215,$a2e1,$3ffe), ($5e39,$5a1f,$98c7,$a2a8,$3ffe),
       ($6819,$0c49,$4303,$a270,$3ffe), ($cdc5,$90d0,$00c1,$a238,$3ffe), ($a188,$25ce,$d1fc,$a1ff,$3ffe),
       ($17a8,$0bb3,$b6ac,$a1c7,$3ffe), ($b6e4,$8544,$aeca,$a18f,$3ffe), ($8939,$d79f,$ba50,$a157,$3ffe),
       ($4cf7,$4a34,$d938,$a11f,$3ffe), ($a612,$26c7,$0b7a,$a0e8,$3ffe), ($4fc2,$b971,$510f,$a0b0,$3ffe),
       ($4e69,$509b,$a9f2,$a078,$3ffe), ($21be,$3d01,$161b,$a041,$3ffe), ($f745,$d1ae,$9583,$a009,$3ffe),
       ($dd06,$6400,$2825,$9fd2,$3ffe), ($f493,$4ba1,$cdf9,$9f9a,$3ffe), ($a651,$e28b,$86f8,$9f63,$3ffe),
       ($d504,$8504,$531d,$9f2c,$3ffe), ($11ae,$91a1,$3260,$9ef5,$3ffe), ($cfa5,$693f,$24bb,$9ebe,$3ffe),
       ($98ff,$6f0b,$2a27,$9e87,$3ffe), ($4337,$0879,$429e,$9e50,$3ffe), ($2420,$9d47,$6e18,$9e19,$3ffe),
       ($4720,$977c,$ac90,$9de2,$3ffe), ($a2aa,$6367,$fdff,$9dab,$3ffe), ($4e03,$6f9f,$625e,$9d75,$3ffe),
       ($b751,$2cff,$d9a7,$9d3e,$3ffe), ($d9e2,$0eaa,$63d3,$9d08,$3ffe), ($74cb,$8a07,$00db,$9cd2,$3ffe),
       ($41be,$16c0,$b0ba,$9c9b,$3ffe), ($2c2d,$2ec3,$7368,$9c65,$3ffe), ($88b4,$4e40,$48df,$9c2f,$3ffe),
       ($4cc1,$f3aa,$3118,$9bf9,$3ffe), ($4688,$9fb3,$2c0e,$9bc3,$3ffe), ($5539,$d54e,$39b9,$9b8d,$3ffe),
       ($a17e,$19ad,$5a14,$9b57,$3ffe), ($d63d,$f441,$8d16,$9b21,$3ffe), ($599b,$eeb9,$d2bb,$9aeb,$3ffe),
       ($864a,$94ff,$2afc,$9ab6,$3ffe), ($e51a,$753b,$95d2,$9a80,$3ffe), ($66ca,$1fd1,$1337,$9a4b,$3ffe),
       ($9e27,$275d,$a324,$9a15,$3ffe), ($fa65,$20b7,$4593,$99e0,$3ffe), ($01c4,$a2f1,$fa7d,$99aa,$3ffe),
       ($8c77,$4751,$c1dd,$9975,$3ffe), ($ffcf,$a959,$9bab,$9940,$3ffe), ($89aa,$66c1,$87e2,$990b,$3ffe),
       ($5c25,$1f75,$867b,$98d6,$3ffe), ($e996,$7597,$976f,$98a1,$3ffe), ($20c6,$0d80,$bab9,$986c,$3ffe),
       ($a96f,$8db8,$f051,$9837,$3ffe), ($2103,$9eff,$3832,$9803,$3ffe), ($57ab,$ec43,$9255,$97ce,$3ffe),
       ($8d99,$22a6,$feb5,$9799,$3ffe), ($b08e,$f17a,$7d49,$9765,$3ffe), ($99b3,$0a41,$0e0e,$9731,$3ffe),
       ($4ba3,$20ac,$b0fb,$96fc,$3ffe), ($30cb,$ea9a,$660a,$96c8,$3ffe), ($5a00,$2018,$2d37,$9694,$3ffe),
       ($bd5c,$7b60,$0679,$9660,$3ffe), ($7560,$b8d9,$f1cb,$962b,$3ffe), ($0054,$9714,$ef27,$95f7,$3ffe),
       ($7fef,$d6cc,$fe86,$95c3,$3ffe), ($f940,$3ae8,$1fe3,$9590,$3ffe), ($94d5,$8878,$5336,$955c,$3ffe),
       ($df2c,$86b2,$987a,$9528,$3ffe), ($0961,$fef7,$efa8,$94f4,$3ffe), ($2a1f,$bccb,$58bb,$94c1,$3ffe),
       ($7ed3,$8ddb,$d3ac,$948d,$3ffe), ($ad28,$41f9,$6075,$945a,$3ffe), ($04b6,$ab1c,$ff0f,$9426,$3ffe),
       ($c105,$9d5c,$af75,$93f3,$3ffe), ($4bc1,$eef9,$71a0,$93c0,$3ffe), ($7f3c,$7851,$458b,$938d,$3ffe),
       ($e92c,$13e6,$2b2f,$935a,$3ffe), ($0dac,$9e5c,$2285,$9327,$3ffe), ($aa7c,$f673,$2b88,$92f4,$3ffe),
       ($fa89,$fd0f,$4632,$92c1,$3ffe), ($f9ac,$9531,$727d,$928e,$3ffe), ($a8b4,$a3f8,$b062,$925b,$3ffe),
       ($51ad,$10a0,$ffdc,$9228,$3ffe), ($cc64,$c481,$60e3,$91f6,$3ffe), ($c336,$ab11,$d373,$91c3,$3ffe),
       ($f815,$b1df,$5785,$9191,$3ffe), ($89d3,$c896,$ed13,$915e,$3ffe), ($39b3,$e0f9,$9417,$912c,$3ffe),
       ($b12b,$eee4,$4c8b,$90fa,$3ffe), ($c7f9,$e84d,$1669,$90c8,$3ffe), ($ca6b,$c540,$f1ab,$9095,$3ffe),
       ($bfef,$7fe0,$de4b,$9063,$3ffe), ($b1dc,$1466,$dc43,$9031,$3ffe), ($f285,$8120,$eb8c,$8fff,$3ffe),
       ($6481,$c672,$0c21,$8fce,$3ffe), ($c23d,$e6d1,$3dfc,$8f9c,$3ffe), ($e5c4,$e6c8,$8117,$8f6a,$3ffe),
       ($10d3,$ccf4,$d56c,$8f38,$3ffe), ($3520,$a201,$3af5,$8f07,$3ffe), ($3ce9,$70af,$b1ac,$8ed5,$3ffe),
       ($53c0,$45cd,$398b,$8ea4,$3ffe), ($2f95,$303a,$d28c,$8e72,$3ffe), ($5a01,$40e3,$7ca9,$8e41,$3ffe),
       ($79d3,$8ac4,$37dc,$8e10,$3ffe), ($9cd6,$22e6,$0420,$8ddf,$3ffe), ($81d8,$205f,$e16e,$8dad,$3ffe),
       ($e2f8,$9c50,$cfc0,$8d7c,$3ffe), ($c024,$b1e7,$cf11,$8d4b,$3ffe), ($a9e6,$7e5b,$df5b,$8d1a,$3ffe),
       ($0c60,$20ee,$0098,$8cea,$3ffe), ($7a95,$bae9,$32c1,$8cb9,$3ffe), ($f9e9,$6fa0,$75d2,$8c88,$3ffe),
       ($4dde,$646f,$c9c4,$8c57,$3ffe), ($4414,$c0b6,$2e91,$8c27,$3ffe), ($0085,$adde,$a434,$8bf6,$3ffe),
       ($4a01,$5754,$2aa7,$8bc6,$3ffe), ($d6e7,$ea8b,$c1e3,$8b95,$3ffe), ($9a18,$96fb,$69e4,$8b65,$3ffe),
       ($1032,$8e1e,$22a3,$8b35,$3ffe), ($8cfe,$0370,$ec1b,$8b04,$3ffe), ($8924,$2c72,$c645,$8ad4,$3ffe),
       ($f018,$40a4,$b11c,$8aa4,$3ffe), ($6e47,$7989,$ac9a,$8a74,$3ffe), ($bf7f,$12a1,$b8ba,$8a44,$3ffe),
       ($fd9a,$496e,$d575,$8a14,$3ffe), ($ef5e,$5d70,$02c6,$89e5,$3ffe), ($57a4,$9025,$40a7,$89b5,$3ffe),
       ($44b2,$2507,$8f13,$8985,$3ffe), ($5fdd,$618e,$ee03,$8955,$3ffe), ($3d5e,$8d2e,$5d72,$8926,$3ffe),
       ($ac6b,$f155,$dd5a,$88f6,$3ffe), ($078c,$d96e,$6db6,$88c7,$3ffe), ($8527,$92da,$0e80,$8898,$3ffe),
       ($8851,$6cf7,$bfb2,$8868,$3ffe), ($f1d4,$b919,$8146,$8839,$3ffe), ($717c,$ca8e,$5337,$880a,$3ffe),
       ($d792,$f698,$357f,$87db,$3ffe), ($66a0,$9473,$2819,$87ac,$3ffe), ($256c,$fd4e,$2afe,$877d,$3ffe),
       ($3131,$8c4e,$3e2a,$874e,$3ffe), ($1010,$9e8d,$6196,$871f,$3ffe), ($03c3,$9318,$953d,$86f0,$3ffe),
       ($5c88,$caef,$d919,$86c1,$3ffe), ($cc4a,$a905,$2d25,$8693,$3ffe), ($ba04,$923f,$915b,$8664,$3ffe),
       ($9563,$ed72,$05b5,$8636,$3ffe), ($2a9f,$2364,$8a2f,$8607,$3ffe), ($f694,$9ec9,$1ec1,$85d9,$3ffe),
       ($7b15,$cc48,$c367,$85aa,$3ffe), ($9376,$1a72,$781c,$857c,$3ffe), ($c95d,$f9c8,$3cd8,$854e,$3ffe),
       ($a9c0,$dcb8,$1198,$8520,$3ffe), ($1a29,$379c,$f656,$84f1,$3ffe), ($ae30,$80b8,$eb0b,$84c3,$3ffe),
       ($fd30,$303e,$efb3,$8495,$3ffe), ($f83c,$c049,$0447,$8468,$3ffe), ($4046,$acde,$28c3,$843a,$3ffe),
       ($7c8b,$73e9,$5d21,$840c,$3ffe), ($b132,$9541,$a15b,$83de,$3ffe), ($9628,$92a4,$f56c,$83b0,$3ffe),
       ($ee37,$efb6,$594e,$8383,$3ffe), ($de5a,$3203,$ccfd,$8355,$3ffe), ($4547,$e0fc,$5071,$8328,$3ffe),
       ($1333,$85f6,$e3a7,$82fa,$3ffe), ($a1d7,$ac2b,$8698,$82cd,$3ffe), ($0ca8,$e0bb,$393f,$82a0,$3ffe),
       ($894c,$b2a5,$fb97,$8272,$3ffe), ($c048,$b2ce,$cd9a,$8245,$3ffe), ($25ec,$73fc,$af43,$8218,$3ffe),
       ($5370,$8ad4,$a08c,$81eb,$3ffe), ($6056,$8dde,$a170,$81be,$3ffe), ($3bfd,$1581,$b1ea,$8191,$3ffe),
       ($0773,$bc03,$d1f3,$8164,$3ffe), ($6f7c,$1d88,$0188,$8138,$3ffe), ($06d4,$d814,$40a1,$810b,$3ffe),
       ($a0af,$8b85,$8f3b,$80de,$3ffe), ($ab6c,$d999,$ed4f,$80b1,$3ffe), ($8b84,$65e8,$5ad9,$8085,$3ffe),
       ($f6b1,$d5e5,$d7d2,$8058,$3ffe), ($4f51,$d0e0,$6436,$802c,$3ffe), ($0000,$0000,$0000,$8000,$3ffe));

  B: array[0..NXT div 2] of THexExtW = (
       ($0000,$0000,$0000,$0000,$0000), ($75f1,$bc63,$885f,$c06e,$3fbc), ($a1ed,$30c5,$4cd4,$a45b,$bfbd),
       ($7996,$d6b2,$a131,$ee62,$bfbc), ($23fa,$9c3e,$8b4d,$f581,$bfbd), ($ed9b,$4bb4,$c430,$8aba,$3fbd),
       ($3a32,$426a,$018d,$c4b4,$bfbd), ($9352,$0f8b,$4e4d,$decd,$3fbd), ($ff99,$62ba,$7628,$f84b,$3fbd),
       ($8143,$0442,$24f9,$b4b8,$3fbc), ($8d5b,$21a9,$86c7,$d2df,$3fbc), ($f811,$5153,$4607,$8019,$bfbd),
       ($8a8b,$f824,$b494,$b795,$bfb7), ($71b9,$b976,$bc54,$b92d,$bfbc), ($b9e0,$6310,$046b,$fced,$bfbd),
       ($f03c,$3c66,$a4f9,$cbc8,$3fbd), ($1f87,$db30,$18f5,$f73a,$3fbd), ($f2a3,$b0d9,$b054,$9480,$bfbb),
       ($a758,$8798,$7cdc,$b74d,$bfbd), ($7f23,$646f,$a919,$f15b,$3fb5), ($2da7,$b85c,$ab19,$bb3f,$bfbb),
       ($b86e,$5d58,$808e,$d6f3,$3fbd), ($cfcf,$b4c8,$42f5,$ec3f,$3fbc), ($e18c,$0bc4,$2a38,$ad64,$3fbd),
       ($7226,$291b,$39ed,$8cac,$3fbd), ($7cff,$5c4a,$6191,$ab5c,$3fbd), ($9321,$30a3,$3c06,$95de,$3fbd),
       ($f8e2,$0dde,$ca3a,$8736,$3fbc), ($f624,$4c97,$5b6d,$c01a,$3fbd), ($f23b,$506c,$1942,$b491,$bfbd),
       ($97d9,$0bf0,$0910,$9f3a,$3fbc), ($2748,$d2ff,$abf4,$bc2d,$bfbc), ($ac15,$3e46,$2932,$bf4a,$bfbc),
       ($20bd,$8514,$d685,$d4ef,$3fbd), ($259d,$efed,$3b0a,$f6e3,$bfba), ($d9ce,$8932,$32e4,$e4ef,$bfbc),
       ($8fbc,$58e1,$21a1,$f22f,$3fbd), ($43fe,$5717,$6dff,$e9b8,$bfbb), ($5160,$51be,$8fac,$f895,$3fbd),
       ($fbdb,$24a1,$8024,$83e1,$bfbb), ($8765,$76dd,$7a52,$f2f4,$3fbb), ($5d6e,$3d3f,$4342,$afb4,$bfbc),
       ($ddb5,$c442,$8805,$cbc4,$3fbd), ($8e61,$ab24,$c4aa,$d777,$bfbd), ($eecf,$5980,$9079,$9bfe,$3fba),
       ($3963,$50d1,$6c4a,$f8db,$bfbb), ($ac99,$7d75,$df0e,$b207,$bfbd), ($943c,$8a3a,$0dd6,$ba7a,$3fbd),
       ($7944,$ba66,$a091,$cb12,$3fb9), ($40a9,$b97a,$c086,$b4ab,$3fbd), ($9fa0,$e344,$2518,$8d70,$3fbd),
       ($b2ac,$6d49,$83a9,$988d,$bfbb), ($5e9c,$5ee6,$7498,$8be1,$bfbc), ($dc97,$9dfa,$3696,$ad2e,$3fbd),
       ($2972,$b47f,$6af5,$cb3c,$3fbd), ($7771,$7175,$32d9,$8595,$bfbd), ($a991,$78a6,$356a,$f610,$3fbc),
       ($43dd,$6e2f,$883f,$f352,$3fbb), ($1d1c,$f187,$dd36,$efdd,$bfbc), ($8358,$36fd,$6208,$9c21,$3fbd),
       ($2a20,$0f6a,$09ae,$bc61,$bfb7), ($83db,$4d18,$7e2a,$bf6d,$bfbd), ($1cbc,$ed94,$bf8d,$8559,$3fbc),
       ($ad2d,$b32d,$cc57,$bf03,$bfbd), ($ff78,$40b4,$2ee6,$e69a,$3fbc), ($6957,$1b4c,$4e47,$c449,$bfbb),
       ($1837,$3407,$b9db,$8d32,$bfbc), ($7a1e,$aecd,$45ae,$f389,$bfbc), ($58b5,$4c4c,$bdff,$b243,$3fbd),
       ($0795,$517f,$742b,$cf65,$bfba), ($8797,$f1a9,$b158,$dfb2,$3fbd), ($70d3,$3663,$0670,$f563,$3fbc),
       ($0143,$1ef2,$72be,$9124,$3fbb), ($8b0e,$bf14,$c09d,$ff4e,$3fba), ($3fd1,$56aa,$b86d,$b8fb,$3fba),
       ($01a1,$78c9,$ba0c,$f388,$bfbc), ($abf8,$996c,$8e6a,$a4ae,$bfbc), ($7873,$a64a,$a1c0,$c62a,$3fbd),
       ($21f2,$8dc0,$1cc9,$994f,$3fbc), ($7323,$e675,$113e,$a0ae,$3fbc), ($c895,$5069,$e383,$ee53,$bfbb),
       ($3036,$2fc1,$edaa,$9c9b,$bfba), ($14e6,$8c84,$73b9,$ef64,$bfbd), ($344f,$8ec8,$24c7,$bdaa,$3fbd),
       ($9231,$a140,$370b,$b6f8,$bfba), ($3b7d,$2256,$fc67,$9659,$bfbd), ($8383,$0390,$6a5f,$b7c9,$bfbd),
       ($a683,$6654,$5f4b,$f59e,$bfbc), ($2d04,$f5dd,$0dab,$fe3c,$bfbd), ($024d,$6fed,$c73b,$abf4,$bfbd),
       ($0385,$6e1c,$d3ed,$c368,$3fbc), ($aaed,$f499,$7f8f,$b2a1,$3fbd), ($6207,$24ff,$471a,$fb17,$bfbc),
       ($5d3b,$72dc,$e4d9,$a537,$bfbc), ($84aa,$c55d,$d161,$aa1c,$3fbd), ($4ef2,$1b9e,$10d1,$d7bf,$3fbd),
       ($7cde,$9376,$4325,$f8ab,$3fbc), ($0271,$0aa0,$1d48,$e551,$3fbd), ($9ed8,$b162,$20d2,$cf43,$bfbd),
       ($34c3,$ef34,$b6d1,$b5f4,$3fbc), ($e90a,$a2e0,$1584,$83b2,$3fbc), ($5a72,$5b3e,$e2c8,$9d23,$bfbd),
       ($c8f1,$bf6a,$6839,$d094,$bfbd), ($2880,$f54f,$77ca,$e68c,$3fbd), ($0f6b,$4a01,$fab3,$f88a,$3fbd),
       ($96b3,$aca3,$bab8,$f23c,$bfbd), ($a707,$354e,$649a,$de67,$3fbd), ($7c2d,$5246,$6cee,$ee30,$3fba),
       ($3d7b,$a07a,$7aa1,$bf51,$bfbb), ($da5f,$d319,$2351,$8901,$bfbd), ($b2dd,$7562,$4593,$9334,$3fbd),
       ($81bb,$2209,$1a37,$ed61,$bfbd), ($a10c,$25e0,$c093,$aefd,$bfbd), ($a77b,$9ef7,$556e,$d431,$3fbc),
       ($9a84,$e9a8,$fef4,$a3fa,$bfbc), ($a31d,$eff0,$228f,$ee2d,$3fba), ($0718,$8b27,$33a4,$e9aa,$3fbd),
       ($fcc0,$e4bd,$6478,$a79d,$bfb5), ($5478,$3b9f,$65d5,$d96c,$bfbb), ($46fe,$1c8c,$358b,$a80b,$bfbd),
       ($2d0d,$b35b,$bbc2,$dc3c,$3fbb), ($ece5,$6684,$ae95,$a6ec,$bfbd), ($ff37,$4277,$9e7c,$fc36,$3fbc),
       ($46da,$590d,$0d64,$ba4f,$bfbc), ($4b3f,$aa83,$e1bb,$e2cb,$3fb9), ($efe7,$c23c,$bc5c,$aa6d,$3fbd),
       ($0214,$b0cc,$7ff0,$9566,$bfbd), ($a511,$5e39,$6cac,$e15a,$bfbc), ($7d3e,$ea95,$1366,$b2fb,$3fbd),
       ($838b,$8b93,$5d68,$97b3,$3fbd), ($e500,$3363,$6118,$ea37,$bfbb), ($530c,$2670,$0d8c,$e64b,$bfbd),
       ($4670,$e629,$5371,$fb3c,$3fbc), ($0259,$822c,$69cf,$f572,$bfbd), ($e410,$d9a4,$4c4e,$f871,$3fbd),
       ($477b,$1aa3,$469b,$ef41,$bfbb), ($44e3,$25bd,$902d,$f05f,$bfbd), ($4b87,$09bd,$ae13,$cf92,$3fbd),
       ($cb15,$91c4,$d5b7,$90a6,$bfbd), ($70d3,$0dc8,$8646,$a443,$3fbd), ($bf35,$d132,$bf8c,$8367,$bfbc),
       ($dba0,$dea4,$a5d4,$ac2b,$3fbc), ($07b1,$0d44,$02d3,$9637,$3fbc), ($b6fe,$cca3,$b983,$bd67,$3fba),
       ($5d89,$eb34,$5191,$9301,$3fbd), ($79c4,$83b1,$19f9,$b3b9,$bfbd), ($ec93,$bcf2,$7343,$bc2b,$3fbd),
       ($7312,$e23a,$2bdc,$c645,$bfbc), ($00f3,$eb3c,$4764,$cb00,$3fbd), ($eed6,$fbe6,$a3ba,$db85,$bfbd),
       ($7264,$4ca6,$0243,$ec62,$3fbd), ($6cc4,$df77,$4b57,$9b46,$3fbc), ($f4e6,$6a63,$49d8,$a83c,$3fbd),
       ($b45b,$11d7,$8202,$ce57,$3fbc), ($6f19,$7dde,$37b2,$d0ad,$bfbd), ($2f65,$cbd8,$36e8,$9369,$3fbc),
       ($1aa2,$f8c1,$9655,$c274,$bfbd), ($c479,$6b70,$68a1,$a0b4,$3fbd), ($7b4e,$f33b,$e561,$c910,$bfbd),
       ($9774,$519b,$7aa1,$eb4b,$bfbc), ($80d9,$b883,$fb10,$e5eb,$3fbb), ($f6f0,$fdf9,$e512,$b70f,$bfbd),
       ($461e,$9586,$f46f,$d8bc,$3fbd), ($a3fa,$a7e0,$9267,$a257,$bfbc), ($1eec,$781e,$4831,$d1db,$3fba),
       ($b6b0,$331c,$4457,$f0ed,$bfbb), ($20c0,$6268,$a6dd,$ed0b,$bfbd), ($7a5c,$038e,$e827,$df33,$3fbd),
       ($42b1,$fe60,$f620,$c90b,$bfbd), ($da65,$029a,$04c0,$be97,$3fba), ($7e1e,$3313,$6bef,$fbbc,$bfbd),
       ($6a4d,$915c,$a675,$e91f,$3fbc), ($cd4c,$d87e,$3cf6,$c96e,$3fbb), ($74e4,$eded,$42cb,$81d2,$bfbc),
       ($73e3,$cd04,$5be7,$8fe5,$bfb9), ($379c,$7943,$7193,$cc6f,$bfbd), ($045d,$288c,$c1ec,$bedd,$bfbd),
       ($2374,$c710,$c284,$8d08,$3fbd), ($f035,$9bba,$5ac7,$f914,$3fb6), ($9707,$39e2,$6400,$e647,$bfbd),
       ($93d4,$bc59,$cc3e,$86da,$bfbc), ($1ff1,$9a23,$4714,$d2fe,$bfbc), ($53bd,$d0c8,$d9be,$9cb0,$3fbd),
       ($ae8d,$7102,$4c53,$8d58,$3fbd), ($81c2,$b867,$d0ba,$baaf,$bfbd), ($e89f,$2ec0,$365c,$a250,$bfbb),
       ($1adb,$75e7,$ec6e,$c468,$3fbc), ($714a,$cb7d,$9eef,$b14d,$3fbd), ($1a8c,$5a50,$c9a6,$de7b,$bfbb),
       ($8739,$39e4,$d2c3,$85f0,$3fbd), ($443a,$e0f2,$79fe,$c61c,$bfbc), ($8e76,$65ec,$8180,$cb81,$bfbd),
       ($eded,$5c85,$4630,$8d5a,$3fbd), ($8609,$5b52,$3509,$eaab,$3fba), ($af61,$3833,$5d52,$a0f4,$3fbd),
       ($523a,$33e5,$218f,$9802,$bfbc), ($458a,$7538,$36d0,$91d5,$3fbd), ($48b9,$923c,$75a8,$a2d3,$bfbc),
       ($dde8,$aa3b,$6381,$adcd,$bfba), ($4325,$58c4,$e96f,$bce6,$3fbb), ($7951,$9547,$eb44,$ba2b,$3fbc),
       ($bd1d,$d0c0,$eff8,$e375,$3fbc), ($7fdb,$a096,$f03c,$f15c,$3fbd), ($a970,$4a37,$030c,$c1fd,$bfbc),
       ($c463,$89b6,$afc5,$b319,$bfbc), ($74dd,$962b,$618d,$d125,$3fbc), ($6a72,$5204,$7a7a,$aadf,$3fbb),
       ($ce25,$a73d,$8d59,$cdbb,$bfb8), ($9d82,$e5ac,$8e0a,$fd6d,$3fba), ($1e38,$d3d9,$aefa,$9f6a,$3fbd),
       ($49e1,$424b,$6d6b,$d02d,$bfbd), ($cda6,$6474,$7e55,$9696,$3fbc), ($0fc5,$929d,$2950,$eeb0,$3fbd),
       ($3c94,$28f4,$da14,$b6d7,$3fbd), ($826a,$e8f2,$0890,$cffb,$3fba), ($ba60,$486f,$dcd2,$89fd,$3fbd),
       ($5317,$1e0f,$5132,$b700,$3fbc), ($c524,$6473,$087c,$bbc0,$3fbc), ($0653,$81ca,$4bbd,$e18d,$bfbd),
       ($8588,$002b,$92db,$f4b6,$bfba), ($a594,$e37c,$96d2,$9670,$bfbd), ($14b4,$820b,$7c5f,$b761,$3fbd),
       ($6840,$afa3,$5c3d,$8f46,$bfb9), ($06b3,$cdef,$a9d2,$e00e,$bfba), ($6dfd,$eb58,$af14,$8373,$bfb9),
       ($c4e4,$5cd8,$d15e,$dc47,$bfbb), ($92a3,$13bf,$dd56,$d573,$3fb9), ($bbfe,$55af,$dfa2,$f41c,$bfbb),
       ($0c22,$c368,$1d92,$80ca,$3fba), ($6cb6,$e77c,$d246,$847d,$bfbd), ($08d9,$971d,$cdb2,$d452,$bfbd),
       ($3455,$2381,$3f01,$eac2,$3fbd), ($eac6,$318c,$aed9,$bbf1,$3fbd), ($d5c0,$baa7,$d25f,$dd62,$3fbd),
       ($aeb2,$ebd9,$2275,$df6e,$bfbd), ($b52b,$7922,$eaff,$9190,$3fbd), ($a559,$ed4a,$9f15,$e85c,$3fbc),
       ($4b5f,$b50b,$1f98,$a6f2,$bfbc), ($eaa7,$4742,$8ed8,$91f4,$bfbc), ($e75c,$7784,$690e,$c483,$3fbc),
       ($f938,$7aac,$91cf,$e8da,$bfbc), ($0e4d,$9bc6,$9e8e,$b643,$3fbb), ($a5d1,$6b5b,$62c2,$f030,$3fba),
       ($e760,$a1eb,$147f,$c700,$bfba), ($ac10,$9fe8,$7650,$d7c9,$3fbb), ($b51d,$bff9,$1161,$cd15,$3fbc),
       ($4be1,$1360,$589e,$eff1,$bfbb), ($1f52,$e5ac,$e670,$ded5,$bfbd), ($2ffd,$1948,$1d6d,$f8a9,$3fbc),
       ($63ed,$2320,$a834,$de4e,$3fbc), ($9130,$5b81,$5df2,$c706,$bfbd), ($d9d8,$f72d,$8a90,$bf70,$bfbd),
       ($bde9,$7a29,$ca4f,$f7ca,$3fbd), ($54fb,$4a66,$3ab1,$cef0,$3fba), ($f480,$db9b,$5c66,$94cd,$3fbc),
       ($6abc,$ee23,$ec13,$d654,$bfbd), ($0000,$0000,$0000,$0000,$0000));
var
  w, z, ya, yb, u: extended;
  F, Fa, Fb, G, Ga, Gb, H, Ha, Hb: extended;
  e,i : longint;
  k: integer;
  nflg : boolean;
begin

  if IsNanOrInf(x) or IsNanOrInf(y) then begin
    power := NaN_x;
    exit;
  end;

  {Handle easy cases}
  w := abs(y);
  if w=0.0 then begin
    power := 1.0;
    exit;
  end
  else if (x=0.0) then begin
    if y>0.0 then power := 0.0
    else begin
      {here y<0: if y is odd and x=-0 return -INF}
      if frac(0.5*y)=0.0 then power := PosInf_x
      else power := copysign(PosInf_x,x);
    end;
    exit;
  end
  else if x=1.0 then begin
    power := 1.0;
    exit;
  end
  else if w=1.0 then begin
    if y>0.0 then power := x
    else begin
      if abs(x) > 0.25*MinExtended then power := 1.0/x
      else power := copysign(PosInf_x,x);
    end;
    exit;
  end
  else if w=0.5 then begin
    if y>0.0 then power := sqrt(x)
    else power := 1.0/sqrt(x);
    exit;
  end;

  nflg := false;         {true if x<0 raised to integer power }
  w := floorx(y);
  if w=y then begin
    z := abs(w);
    if z <= 2048.0 then begin
      u := ilogb(x);
      {intpower0 does not catch overflows, therefore call it only if safe.}
      if z*(abs(u) + 1.0) < 16382.0 then begin
        power := intpower0(x,trunc(w));
        exit;
      end;
    end;
  end;
  if x<0.0 then begin
    if w<>y then begin
      {noninteger power of negative number: exception or RTE}
      power := exp(y*ln(x));
      exit;
    end;
    {Find out if the integer exponent is odd or even.}
    w := 2.0*floorx(0.5*y);
    nflg := w<>y;
    x := abs(x);
  end;

  {separate significand from exponent}
  frexp(x, x, e);

  {find significand in antilog table A[]}
  i := 1;
  if x <= extended(A[257])  then i := 257;
  if x <= extended(A[i+128])then inc(i,128);
  if x <= extended(A[i+64]) then inc(i,64);
  if x <= extended(A[i+32]) then inc(i,32);
  if x <= extended(A[i+16]) then inc(i,16);
  if x <= extended(A[i+8])  then inc(i, 8);
  if x <= extended(A[i+4])  then inc(i, 4);
  if x <= extended(A[i+2])  then inc(i, 2);
  if x >= extended(A[1])    then i := -1;
  inc(i);

  {Find (x - A[i])/A[i] in order to compute log(x/A[i]):}
  {log(x) = log( a x/a ) = log(a) + log(x/a)            }
  {log(x/a) = log(1+v),  v = x/a - 1 = (x-a)/a          }
  x := x - extended(A[i]);
  x := x - extended(B[i div 2]);
  x := x / extended(A[i]);

  {rational approximation for log(1+v) = v - v**2/2 + v**3 P(v) / Q(v)}
  z := x*x;

  {w:= polevl(x,P,3)}
  w := ((P[0]*x + P[1])*x + P[2])*x + P[3];
  {u = p1evl(x,Q,3)}
  u := ((x + Q[0])*x + Q[1])*x + Q[2];

  w := x*(z*w/u) - 0.5*z;

  {Convert to base 2 logarithm: multiply by log2(e) = 1+LOG2EA}
  z := extended(LOG2EA)*w;
  z := z + w;
  z := z + extended(LOG2EA)*x;
  z := z + x;

  {Compute exponent term of the base 2 logarithm.}
  w := -i;
  w := w/NXT + e;

  {Now base 2 log of x is w + z.  Multiply base 2 log by y, in pseudo multi}
  {precision. Separate y into large part ya and small part yb less than 1/NXT}

  {WE: mul/div with NXT is faster than original ldexp(,LNXT)}

  ya := floorx(y*NXT)/NXT;    {reduc(y)}
  yb := y - ya;

  F  := z*y + w*yb;
  Fa := floorx(F*NXT)/NXT;    {reduc(F)}
  Fb := F - Fa;

  G  := Fa + w*ya;
  Ga := floorx(G*NXT)/NXT;    {reduc(G)}
  Gb := G - Ga;

  H  := Fb + Gb;
  Ha := floorx(H*NXT)/NXT;    {reduc(H)}
  w  := NXT*(Ga+Ha);

  {Test the power of 2 for over/underflow}
  if w>MEXP then begin
    if nflg then power := NegInf_x
    else power := PosInf_x;
    exit;
  end
  else if w<MNEXP then begin
    power := 0.0;
    exit;
  end;

  e := trunc(w);
  Hb := H - Ha;
  if Hb>0.0 then begin
    inc(e);
    Hb := Hb - 1.0/NXT;
  end;

  {Now the product y * log2(x) = Hb + e/NXT. Compute base 2 exponential}
  {of Hb, where -1/NXT <= Hb <= 0.}

  {z := Hb*polevl(Hb,R,6);    z=2**Hb-1}
  z := R[0];
  for k:=1 to 6 do z := z*Hb + R[k];
  z := Hb*z;

  {Express e/NXT as an integer plus a negative number of (1/NXT)ths.}
  {Find lookup table entry for the fractional power of 2.}
  if e<0 then i := 0 else i := 1;
  inc(i, e div NXT);
  e := NXT*i - e;

  w := extended(A[e]);
  z := w + w*z;     {2**-e * ( 1 + (2**Hb-1) )     }
  z := ldexp(z,i);  {multiply by integer power of 2}

  if nflg then power := -z {odd integer exponent and x<0}
  else power := z;

end;


{---------------------------------------------------------------------------}
function powm1(x,y: extended): extended;
  {-Return x^y - 1; special code for small x,y}
var
  p: extended;
begin
  if y=0.0 then begin
    powm1 := 0.0;
    exit;
  end;
  if (x>0.0) and ((x<2.0) or (abs(y)<2.0)) then begin
    p := y*ln(x);
    if abs(p) < 4.0 then begin
      powm1 := expm1(p);
      exit;
    end;
  end;
  powm1 := power(x,y)-1.0;
end;


{---------------------------------------------------------------------------}
function pow1pm1(x,y: extended): extended;
  {-Return (1+x)^y - 1; special code for small x,y}
begin
  if (x=0.0) or (y=0.0) then pow1pm1 := 0.0
  else if y=1.0 then pow1pm1 := x
  else if y=2.0 then pow1pm1 := x*(2.0+x)
  else if abs(x)<1.0 then pow1pm1 := expm1(y*ln1p(x))
  else pow1pm1 := powm1(1.0+x,y);
end;


{---------------------------------------------------------------------------}
function pow1p_ex(x,y: extended; usexx: boolean): extended;
  {-Return (1+x)^y, x > -1}
var
  z,w,v: extended;
  xx: ext2;
const
  vmax = 11333.0; {<ln(MaxExtended/2^33): to allow ext2 splitting}
begin
  if (y=0.0) or (x=0.0) then begin
    pow1p_ex := 1.0;
    exit;
  end;
  z := abs(x);
  w := abs(y);
  v := y*ln1p(x);
  if abs(v) <= 2.5 then begin
    if (z<=0.0625) or ((z<=0.5) and (w<25.0)) then begin
      pow1p_ex := exp(v);
      exit;
    end;
  end;
  if usexx and (w>1.0) and (v<vmax) then begin
    w := frac(y);
    xxto2x(1.0,x,xx);
    pow2xi(xx,int(y),xx);
    if w<>0.0 then begin
      w := pow1p(x,w);
      mul21x(xx,w,xx);
    end;
    pow1p_ex := xx.h+xx.l;
  end
  else pow1p_ex := power(1.0+x, y);
end;


{---------------------------------------------------------------------------}
function pow1p(x,y: extended): extended;
  {-Return (1+x)^y, x > -1, with ext2 arithmetic for critical values}
begin
  pow1p := pow1p_ex(x,y,true);
end;


{---------------------------------------------------------------------------}
function pow1pf(x,y: extended): extended;
  {-Return (1+x)^y, x > -1, without ext2, less accurate than pow1p}
begin
  pow1pf := pow1p_ex(x,y,false);
end;


{---------------------------------------------------------------------------}
function powpi2k(k,n: longint): extended;
  {-Return accurate scaled powers of Pi, result = (Pi*2^k)^n}
var
  m: longint;
  x: extended;
const
  lph  : extended = 13529/8192;      {high part log2(Pi) = 1.6514892578125}
  lpl  = 0.6871659818798043279e-5;   {low  part log2(Pi)}
  lmax = +16384;
  lmin = -16447;
begin
  x := (k + (lph+lpl))*n;
  if x >= lmax then powpi2k := PosInf_x
  else if x < lmin then powpi2k := 0.0
  else begin
    x := n*lph;
    m := round(x);
    x := (x - m) + n*lpl;
    m := m + k*n;
    x := exp2(x);
    powpi2k := ldexp(x,m);
  end;
end;


{---------------------------------------------------------------------------}
function powpi(n: longint): extended;
  {-Return accurate powers of Pi, result = Pi^n}
begin
  powpi := powpi2k(0,n);
end;


{---------------------------------------------------------------------------}
function compound(x: extended; n: longint): extended;
  {-Return (1+x)^n; accurate version of Delphi/VP internal function}
var
  xx: ext2;
begin
  {Note: IEEE Std 754-2008, Table 9.1 requires x >= -1}
  {This is OK for pow1p but IMO nonsense for compound!}
  xxto2x(1.0,x,xx);
  pow2xi(xx,n,xx);
  compound := xx.h+xx.l;
end;


{---------------------------------------------------------------------------}
function comprel(x,n: extended): extended;
  {-Return ((1+x)^n-1)/x; accurate version of Delphi/VP internal function}
var
  xx,tt: ext2;
begin
  {Relative compound function: (compound(x,n)-1)/x, see math.annuity2}
  if x=0.0 then comprel := n
  else begin
    xxto2x(1.0,x,xx);
    pow2xi(xx,round(n),xx);
    tt.h := 1.0;
    tt.l := 0.0;
    sub2x(xx,tt,xx);
    tt.h := x;
    tt.l := 0.0;
    div2x(xx,tt,xx);
    comprel := xx.h+xx.l;
  end;
end;


(*
{---------------------------------------------------------------------------}
function comprel_fast(x: extended; n: extended): extended;
  {-Return the ((1+x)^n-1)/x}
begin
  {Faster comprel without ext2, but there are large errors up to 3000 eps. }
  {It is still much better than the Delphi implementation, see t_amathx.cmp}
  if n=0 then comprel_fast := n
  else if abs(x) < sqrt_epsh/n then begin
    comprel_fast := (1.0 + 0.5*(n-1)*x)*n;
  end
  else comprel_fast := pow1pm1(x,n)/x;
end;
*)


{---------------------------------------------------------------------------}
function arccos(x: extended): extended;
  {-Return the inverse circular cosine of x, |x| <= 1}
begin
  {basic formula arccos(x) = arctan(sqrt(1-x^2)/x))}
  if abs(x)=1.0 then begin
    if x<0.0 then arccos := Pi else arccos := 0.0;
  end
  else arccos := arctan2(sqrt((1.0-x)*(1.0+x)),x)
end;


{---------------------------------------------------------------------------}
function arccos1m(x: extended): extended;
  {-Return arccos(1-x), 0 <= x <= 2, accurate even for x near 0}
begin
  if x>=0.5 then arccos1m := arccos(1.0-x)
  else begin
    {use arccos(z) = 2arcsin(sqrt((1-z)/2)), ref: }
    {http://functions.wolfram.com/01.13.06.0043.01}
    arccos1m := 2.0*arcsin(sqrt(0.5*x));
  end;
end;


{---------------------------------------------------------------------------}
function archav(x: extended): extended;
  {-Return the inverse haversine archav(x), 0 <= x <= 1}
begin
  {archav(x) = 2*arcsin(sqrt(x)) = 2*arctan2(sqrt(x), sqrt(1-x))}
  archav := 2.0*arctan2(sqrt(x), sqrt(1.0-x));
end;


{---------------------------------------------------------------------------}
function arcsin(x: extended): extended;
  {-Return the inverse circular sine of x, |x| <= 1}
begin
  {basic formula arcsin(x) = arctan(x/sqrt(1-x^2))}
  arcsin := arctan2(x, sqrt((1.0-x)*(1.0+x)))
end;


{---------------------------------------------------------------------------}
function sec(x: extended): extended;
  {-Return the circular secant of x, x mod Pi <> Pi/2}
begin
  sec := 1.0/cos(x);
end;


{---------------------------------------------------------------------------}
function csc(x: extended): extended;
  {-Return the circular cosecant of x, x mod Pi <> 0}
begin
  csc := 1.0/sin(x);
end;


{---------------------------------------------------------------------------}
function coth(x: extended): extended;
  {-Return the hyperbolic cotangent of x, x<>0}
begin
  coth := 1.0/tanh(x);
end;


{---------------------------------------------------------------------------}
function sech(x: extended): extended;
  {-Return the hyperbolic secant of x}
begin
  if abs(x) > ln_MaxExt then sech := 0.0
  else sech := 1.0/cosh(x);
end;


{---------------------------------------------------------------------------}
function csch(x: extended): extended;
  {-Return the hyperbolic cosecant of x, x<>0}
begin
  if abs(x) > ln_MaxExt then csch := 0.0
  else csch := 1.0/sinh(x);
end;


{---------------------------------------------------------------------------}
function arcsinh(x: extended): extended;
  {-Return the inverse hyperbolic sine of x}
var
  t,z: extended;
const
  t0 = 0.43e10;    {sqrt(t0^2 + 1) = t0}
  t1 = 0.59e4932;  {t1 < Maxextended/2 }
const
  CSN = 24;
  CSAH: array[0..CSN-1] of THexExtW = (     {chebyshev((arcsinh(x)/x-1), x=-1..1, 1e-20);}
          ($74C4,$2C57,$F726,$8346,$BFFC),  {-0.128200399117381863433721273592    }
          ($60FB,$E565,$99EE,$F0E4,$BFFA),  {-0.588117611899517675652117571383e-1 }
          ($4327,$5038,$DAB6,$9AE8,$3FF7),  {+0.472746543221248156407252497561e-2 }
          ($E4B7,$CEAC,$CB4F,$8174,$BFF4),  {-0.493836316265361721013601747903e-3 }
          ($F880,$4F1E,$8FBD,$F564,$3FF0),  {+0.585062070585574122874948352586e-4 }
          ($AE3D,$780B,$06F9,$FA8D,$BFED),  {-0.746699832893136813547550692112e-5 }
          ($8655,$9CEE,$EACE,$865F,$3FEB),  {+0.100116935835581992659661920174e-5 }
          ($D271,$A5F8,$C535,$9549,$BFE8),  {-0.139035438587083336086164721962e-6 }
          ($2CEA,$53B7,$9C56,$AA47,$3FE5),  {+0.198231694831727935473173603337e-7 }
          ($051A,$3D90,$00CD,$C63D,$BFE2),  {-0.288474684178488436127472674074e-8 }
          ($79C8,$31D9,$DC1C,$EA98,$3FDF),  {+0.426729654671599379534575435248e-9 }
          ($7F9A,$EA05,$5578,$8CAF,$BFDD),  {-0.639760846543663578687526705064e-10}
          ($CC8C,$E66B,$2C10,$AAA1,$3FDA),  {+0.969916860890647041478747328871e-11}
          ($A487,$C2BA,$24E9,$D0EA,$BFD7),  {-0.148442769720437708302434249871e-11}
          ($77D3,$838C,$C3D7,$80EF,$3FD5),  {+0.229037379390274479880187613939e-12}
          ($BF67,$A9BA,$A045,$A046,$BFD2),  {-0.355883951327326451644472656709e-13}
          ($89C0,$7126,$8F52,$C876,$3FCF),  {+0.556396940800567901824410323292e-14}
          ($CA8D,$2192,$F0F4,$FC17,$BFCC),  {-0.874625095996246917448307610618e-15}
          ($3213,$353B,$6AE5,$9F47,$3FCA),  {+0.138152488445267754259259259259e-15}
          ($7B9A,$3E40,$512C,$CA25,$BFC7),  {-0.219166882829054218109516486403e-16}
          ($E036,$4ADC,$8495,$80C6,$3FC5),  {+0.349046585251352023737824855452e-17}
          ($BAAB,$4F80,$8CCC,$A4A6,$BFC2),  {-0.557857884196310306067662486759e-18}
          ($D7E8,$6819,$4637,$D332,$3FBF),  {+0.894451477596418179416166864814e-19}
          ($A859,$BFDE,$0201,$87D9,$BFBD)); {-0.143834333235934314272184000000e-19}
         {($5CEB,$BE05,$504B,$AF3C,$3FBA)}  {+0.231922385806700541867783333333e-20}
var
  CSA: array[0..CSN-1] of extended absolute CSAH;
begin
  t := abs(x);
  if t<=1.0 then begin
    if t <= sqrt_epsh then begin
      {arcsinh(x) = x*(1 - 1/6*x^2 + 3/40*x^4 + O(x^6))}
      arcsinh := x;
    end
    else begin
      z := 2.0*x*x - 1.0;
      t := CSEvalX(z, CSA, CSN);
      arcsinh := x + x*t;
    end;
  end
  else begin
    if t >= t0 then begin
      {skip sqrt() because sqrt(t^2+1) = t}
      if t <= t1 then z := ln(2.0*t)
      else z := ln(t)+ln2
    end
    else z := ln(t + sqrt(1.0+t*t));
    if x>0.0 then arcsinh := z
    else arcsinh := -z;
  end;
end;


{---------------------------------------------------------------------------}
function ac_help(x: extended): extended;
  {-Calculate ln1p(x+sqrt(2x+x^2)}
var
  y: extended;
begin
  x := x+sqrt(2.0*x+x*x);
  {see ln1p for the formula}
  y := 1.0 + x;
  if y=1.0 then ac_help := x
  else ac_help := ln(y) + (x-(y-1.0))/y;
end;


{---------------------------------------------------------------------------}
function arccosh(x: extended): extended;
  {-Return the inverse hyperbolic cosine, x >= 1. Note: for x near 1 the}
  { function arccosh1p(x-1) should be used to reduce cancellation errors!}
begin
  if x=1.0 then arccosh := 0
  else begin
    if x>1E10 then begin
      {skip sqrt() calculation because sqrt(x^2-1) = x}
      arccosh := ln(x)+ln2;
    end
    else if x>2.0 then begin
      {arccosh := ln(x+sqrt((x-1)*(x+1)))}
      arccosh := ln(2.0*x - 1.0/(x+sqrt(x*x-1.0)))
    end
    else begin
      {arccosh = ln1p(y+sqrt(2*y + y*y)), y=x-1}
      arccosh := ac_help(x-1.0);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function arccosh1p(x: extended): extended;
  {-Return arccosh(1+x), x>=0, accurate even for x near 0}
begin
  if x=0.0 then arccosh1p := 0
  else if x<1.0 then begin
    {arccosh(X) = ln(X + sqrt(X^2 - 1)), substituting X = 1+x gives}
    {arccosh1p(x) = ln1p(x+sqrt((X-1)*(X+1))) = ln1p(x+sqrt(x*(x+2)))}
    {             = ln1p(x+sqrt(2*x + x*x))}
    arccosh1p := ac_help(x);
  end
  else arccosh1p := arccosh(1.0+x);
end;


{---------------------------------------------------------------------------}
function arctanh(x: extended): extended;
  {-Return the inverse hyperbolic tangent of x, |x| < 1}
var
  t: extended;
const
  t0 = 2.3E-10;
begin
  t := abs(x);
  if t<t0 then begin
    {arctanh(x) = x + 1/3*x^3 + 1/5*x^5 + 1/7*x^7 + O(x^9)}
    arctanh := x;
  end
  else begin
    {arctanh(x) = 0.5*ln((1+x)/(1-x))     = 0.5*ln((1-x+2x)/(1-x)) }
    {           = 0.5*ln(1+2x/(1-x))      = 0.5*ln1p(2x/(1-x))     }
    {  or       = 0.5*(ln(1+x)-ln(1-x))   = 0.5*(ln1p(x)-ln1p(-x)) }
    {  or       = 0.5*ln(1+2x+2x^2/(1-x)) = 0.5*ln1p(2x+2x^2/(1-x))}
    if t<0.75 then t := 0.5*(ln1p(t)-ln1p(-t))
    else t := 0.5*ln((1.0+t)/(1.0-t));
    if x>0.0 then arctanh := t
    else arctanh := -t;
  end;
end;


{---------------------------------------------------------------------------}
function arccot(x: extended): extended;
  {-Return the sign symmetric inverse circular cotangent; arccot(x) = arctan(1/x), x <> 0}
begin
  if abs(x) > 1E-20 then arccot := arctan(1.0/x)
  else begin
    if x>=0.0 then arccot := Pi_2
    else arccot := -Pi_2;
  end;
end;


{---------------------------------------------------------------------------}
function arccotc(x: extended): extended;
  {-Return the continuous inverse circular cotangent; arccotc(x) = Pi/2 - arctan(x)}
begin
  if x > 1E-20 then arccotc := arctan(1.0/x)
  else arccotc := Pi_2 - arctan(x);
end;


{---------------------------------------------------------------------------}
function arccsc(x: extended): extended;
  {-Return the inverse cosecant of x, |x| >= 1}
var
  y,z: extended;
begin
  {arccsc = arcsin(1.0/x) = arctan(1/sqrt(x^2-1))}
  y := abs(x);
  if y >= 1E10 then begin
    {arccsc(x) ~ 1/x + 1/6/x^3 + 3/40/x^5 + O(1/x^6) for large x}
    z := 1.0/y;
  end
  else if y > 2.0 then z := arctan2(1.0, sqrt(y*y-1.0))
  else begin
    z := y-1.0;
    z := arctan2(1.0, sqrt(z*y+z));
  end;
  if x>0.0 then arccsc := z
  else arccsc := -z;
end;


{---------------------------------------------------------------------------}
function arcsec(x: extended): extended;
  {-Return the inverse secant of x, |x| >= 1}
var
  t: extended;
begin
  if abs(x) >= 0.5e10 then begin
    {avoid x^2 overflow and/or arctan evaluation}
    arcsec := Pi_2 - 1.0/x;
  end
  else begin
    t := arctan(sqrt((x-1.0)*(x+1.0)));
    if x>0.0 then arcsec := t
    else arcsec := Pi-t;
  end;
end;


{---------------------------------------------------------------------------}
function arccoth(x: extended): extended;
  {-Return the inverse hyperbolic cotangent of x, |x| > 1}
var
  t: extended;
begin
  t := abs(x);
  {compute t = arccoth(|x|)}
  if t<32.0 then t := 0.5*ln1p(2.0/(t-1.0))
  else t := -0.5*ln1p(-2.0/(t+1.0));
  {adjust sign}
  if x>0.0 then arccoth := t
  else arccoth := -t;
end;


{---------------------------------------------------------------------------}
function arcsech(x: extended): extended;
  {-Return the inverse hyperbolic secant of x, 0 < x <= 1}
var
  t: extended;
begin
  t := 1.0-x;
  if t=0 then arcsech := 0.0
  else if x<=1E-10 then begin
    {avoid overflow in ac_help}
    {arcsech(x) = (ln(2)-ln(x)) -1/4*x^2 -3/32*x^4 + O(x^6)}
    arcsech := -ln(0.5*x);
  end
  else begin
    {arcsech(x) = arccosh(1/x), see arccosh branch for x<=2}
    arcsech := ac_help(t/x);
  end;
end;


{---------------------------------------------------------------------------}
function arccsch(x: extended): extended;
  {-Return the inverse hyperbolic cosecant of x, x <> 0}
var
  t: extended;
begin
  t := abs(x);
  if t<=1.0 then begin
    if t<=5e-10 then begin
      {avoid overflow for 1/t^2}
      {arccsch(x) = (ln(2)-ln(x)) + 1/4*x^2 -3/32*x^4 + O(x^6)}
      t := -ln(0.5*t)
    end
    else begin
      {ln(1/t + sqrt(1+1/t^2)) = ln((1+sqrt(1+t^2)/t) = }
      t := ln1p(sqrt(1.0+t*t)) - ln(t);
    end;
    if x>0.0 then arccsch := t
    else arccsch := -t;
  end
  else arccsch := arcsinh(1.0/x);
end;


{---------------------------------------------------------------------------}
{------------- Degrees versions of trig / invtrig functions ----------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure trig_deg(x: extended; var y,z: extended; var n: integer; var m45: boolean);
  {-Reduce x in degrees mod 90; y=x mod 90, |y|<45. z=x/45, m45 if x is multiple of 45}
const
  XMAX = 1e17; {~2^64/180}
const
  c45: single = 45.0;  {Anti-optimize: avoid multiplication with suboptimal inverse}
begin
  {Basic internal reduction routine mod 90. Use Cody/Waite logic, but no}
  {pseudo-multiprecision because 45.0 has only 6 non-zero mantissa bits.}
  if x=0.0 then begin
    y := 0.0;
    z := 0.0;
    n := 0;
    m45 := true;
  end
  else begin
    if abs(x) > XMAX then begin
      {Standard method not reliable, use reduction mod 360 and continue}
      m45 := false;
      x := fmod(x,360.0);
      z := x/c45;
    end
    else begin
      z := x/c45;
      m45 := (frac(z)=0.0) and (frac(x)=0.0);
    end;
    y := floorx(z);
    n := trunc(y - 16.0*floorx(y/16.0));
    if odd(n) then begin
      inc(n);
      y := y + 1.0;
    end;
    n := (n shr 1) and 7;
    y := x - y*c45;
  end;
end;


{---------------------------------------------------------------------------}
procedure sincosd(x: extended; var s,c: extended);
  {-Return sin(x) and cos(x), x in degrees}
var
  y,ss,cc: extended;
  n: integer;
  m45: boolean;
begin
  trig_deg(x,y,ss,n,m45);
  _sincos(y*Pi_180,ss,cc);
  case n and 3 of
      0: begin s:= ss; c:= cc; end;
      1: begin s:= cc; c:=-ss; end;
      2: begin s:=-ss; c:=-cc; end;
    else begin s:=-cc; c:= ss; end;
  end;
  {Avoid -0}
  s := s+0.0;
  c := c+0.0;
end;


{---------------------------------------------------------------------------}
function sind(x: extended): extended;
  {-Return sin(x), x in degrees}
var
  y,z: extended;
  n  : integer;
  m45: boolean;
begin
  trig_deg(x,y,z,n,m45);
  z := y*Pi_180;
  case n and 3 of
      0: sind := system.sin(z);
      1: sind := system.cos(z);
      2: sind := 0.0 - system.sin(z);
    else sind := 0.0 - system.cos(z);
  end;
end;


{---------------------------------------------------------------------------}
function cosd(x: extended): extended;
  {-Return cos(x), x in degrees}
var
  y,z: extended;
  n  : integer;
  m45: boolean;
begin
  trig_deg(x,y,z,n,m45);
  z := y*Pi_180;
  case n and 3 of
      0: cosd := system.cos(z);
      1: cosd := 0.0 - system.sin(z);
      2: cosd := 0.0 - system.cos(z);
    else cosd := system.sin(z);
  end;
end;


{---------------------------------------------------------------------------}
function tand(x: extended): extended;
  {-Return tan(x), x in degrees}
var
  y,z: extended;
  n  : integer;
  m45: boolean;
begin
  trig_deg(x,y,z,n,m45);
  if m45 then begin
    z := abs(z);
    y := isign(x);
    case round(4.0*frac(0.25*z)) of
        0: tand := 0.0;
        1: tand := y;
        2: tand := PosInf_x;
      else tand := -y;          {case 3, but keep some stupid compilers happy}
    end;
  end
  else if odd(n) then tand := 0.0 - _cot(y*Pi_180)
  else tand := _tan(y*Pi_180);
end;


{---------------------------------------------------------------------------}
function tanPi(x: extended): extended;
  {-Return tan(Pi*x), result will be 0 for abs(x) >= 2^64}
var
  c,s: extended;
begin
  sincospi(x,s,c);
  if c=0.0 then tanPi := PosInf_x
  else tanPi := s/c;
end;


{---------------------------------------------------------------------------}
function cotd(x: extended): extended;
  {-Return cot(x), x in degrees}
var
  y,z: extended;
  n  : integer;
  m45: boolean;
begin
  trig_deg(x,y,z,n,m45);
  if m45 then begin
    z := abs(z);
    y := isign(x);
    case round(4.0*frac(0.25*z)) of
        0: cotd := PosInf_x;
        1: cotd := y;
        2: cotd := 0.0;
      else cotd := -y;          {case 3, but keep some stupid compilers happy}
    end;
  end
  else if odd(n) then cotd := 0.0 - _tan(y*Pi_180)
  else cotd := _cot(y*Pi_180);
end;


{---------------------------------------------------------------------------}
function arctand(x: extended): extended;
  {-Return the inverse circular tangent of x, result in degrees}
begin
  {exact for extended +-1, +-INF; 90 for x>2^65, -90 for x < -2^65}
  arctand := arctan(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
function arccosd(x: extended): extended;
  {-Return the inverse circular cosine of x, |x| <= 1, result in degrees}
begin
  arccosd := arccos(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
function arccotd(x: extended): extended;
  {-Return the sign symmetric inverse circular cotangent,}
  { arccotd(x) = arctand(1/x), x <> 0, result in degrees}
begin
  arccotd := arccot(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
function arccotcd(x: extended): extended;
  {-Return the continuous inverse circular cotangent;}
  { arccotcd(x) = 90 - arctand(x), result in degrees}
begin
  arccotcd := arccotc(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
function arcsind(x: extended): extended;
  {-Return the inverse circular sine of x, |x| <= 1, result in degrees}
begin
  arcsind := arcsin(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
{---------------------- Elementary numerical functions ---------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function cbrt(x: extended): extended;
  {-Return the cube root of x}
var
  y: extended;
begin
  if (TExtRec(x).xp and $7FFF=$7FFF) or (x=0.0) or (abs(x)=1.0) then begin
    {x is 0, +1, -1, Inf, or NaN}
    cbrt := x;
  end
  else begin
    {calculate initial approximation}
    y := copysign(exp2(log2(abs(x))/THREE),x);
    {perform one Newton step}
    cbrt := y - (y - x/sqr(y))/THREE;
  end;
end;


{---------------------------------------------------------------------------}
function ceil(x: extended): longint;
  {-Return the smallest integer >= x; |x|<=MaxLongint}
var
  i: longint;
begin
  x := modf(x,i);
  if x>0.0 then ceil := i+1
  else ceil := i;
end;


{---------------------------------------------------------------------------}
function ceilx(x: extended): extended;
  {-Return the smallest integer >= x}
var
  t: extended;
begin
  t := int(x);
  if (x<=0.0) or (x=t) then ceilx := t
  else ceilx := t + 1.0
end;


{---------------------------------------------------------------------------}
function floor(x: extended): longint;
  {-Return the largest integer <= x; |x|<=MaxLongint}
var
  i: longint;
begin
  x := modf(x,i);
  if x<0.0 then floor := i-1
  else floor := i;
end;


{---------------------------------------------------------------------------}
function floorx(x: extended): extended;
  {-Return the largest integer <= x}
var
  t: extended;
begin
  t := int(x);
  if (x>=0.0) or (x=t) then floorx := t
  else floorx := t - 1.0;
end;


{---------------------------------------------------------------------------}
function hypot(x,y: extended): extended;
  {-Return sqrt(x*x + y*y)}
begin
  x := abs(x);
  y := abs(y);
  if x>y then hypot := x*sqrt(1.0+sqr(y/x))
  else if x>0.0 then hypot := y*sqrt(1.0+sqr(x/y)) {here y >= x > 0}
  else hypot := y;                                 {here x=0}
end;


{---------------------------------------------------------------------------}
function hypot3(x,y,z: extended): extended;
  {-Return sqrt(x*x + y*y + z*z)}
var
  ax,ay,az,r,s,t: extended;
begin
  ax := abs(x);
  ay := abs(y);
  az := abs(z);
  {Find maximum of absolute values and divide the other two by the maximum}
  if ax > ay then begin
    if ax > az then begin
      r := ay/ax;
      s := az/ax;
      t := ax;
    end
    else begin
      r := ax/az;
      s := ay/az;
      t := az;
    end
  end
  else if ay > az then begin
    r := ax/ay;
    s := az/ay;
    t := ay;
  end
  else if az > 0.0 then begin
    r := ax/az;
    s := ay/az;
    t := az;
  end
  else begin
    hypot3 := 0.0;
    exit;
  end;
  hypot3 := t*sqrt((sqr(r) + sqr(s)) + 1.0);
end;


{---------------------------------------------------------------------------}
function modf(x: extended; var ip: longint): extended;
  {-Return frac(x) and trunc(x) in ip, |x|<=MaxLongint}
begin
  ip := trunc(x);
  modf := x-ip
end;


{---------------------------------------------------------------------------}
function nroot(x: extended; n: integer): extended;
  {-Return the nth root of x; n<>0, x >= 0 if n is even}
var
  y: extended;
begin
  if n<0 then nroot := 1.0/nroot(x,-n)
  else if n<4 then begin
    case n of
        3:  nroot := cbrt(x);
        2:  nroot := sqrt(x);
        1:  nroot := x;
      else  nroot := 1.0/n;
    end
  end
  else if x=0.0 then nroot := 0.0
  else begin
    {if x<0 and n even, log(y) will generate RTE}
    if odd(n) then y := abs(x) else y := x;
    {calculate initial approximation}
    y := copysign(exp2(log2(y)/n),x);
    {perform one Newton step}
    nroot := y - (y - x/intpower(y,n-1))/n
  end;
end;


{$ifdef BIT32}
{---------------------------------------------------------------------------}
function remainder(x,y: extended): extended; assembler; {&Frame-} {&Uses none}
  {-Return the IEEE754 remainder x REM y = x - rmNearest(x/y)*y}
asm
       fld    [y]
       fld    [x]
  @@1: fprem1
       fstsw  ax
       sahf
       jp     @@1
       fstp   st(1)
       fwait
end;
{$else}
{---------------------------------------------------------------------------}
function remainder(x,y: extended): extended;
  {-Return the IEEE754 remainder x REM y = x - rmNearest(x/y)*y}
var
  yh: extended;
  ey,sx: word;
begin
  {Ref: FDLIBM 5.3 [5], file e_remainder.c}

  if TExtRec(x).xp and $7FFF = $7FFF then begin
    {x is INF or NaN}
    remainder := NaN_x;
    exit;
  end;

  sx := TExtRec(x).xp and $8000;
  ey := TExtRec(y).xp and $7FFF;
  if ey=$7FFF then with TExtRec(y) do begin
    {x is INF or NaN}
    if (hm=longint($80000000)) and (lm=0) then begin
      {y is INF}
      remainder := x;
    end
    else remainder := NaN_x;
    exit;
  end;

  y := abs(y);
  if ey <$7FFE then begin
    {|y| < 0.5*MaxExtended}
    x := fmod(x,y+y);
  end;
  x := abs(x);
  if ey<2 then begin
    {|y| < 2*MinExtended}
    if x+x > y then begin
      if x=y then x := 0.0
      else x := x-y;
      if x+x >= y then x := x-y;
    end;
  end
  else begin
    yh := 0.5*y;
    if x > yh then begin
      if x=y then x := 0.0
      else x := x-y;
      if x >= yh then x := x-y;
    end;
  end;
  if sx=0 then remainder := x
  else remainder := -x;
end;
{$endif}


{---------------------------------------------------------------------------}
function sqrt1pm1(x: extended): extended;
  {-Return sqrt(1+x)-1, accurate even for x near 0, x>=-1}
begin
  if abs(x)>0.75 then sqrt1pm1 := sqrt(1.0+x)-1.0
  else begin
    {sqrt(1+x)-1 = (sqrt(1+x)-1)*(sqrt(1+x)+1)/(sqrt(1+x)+1) = x/(sqrt(1+x)-1)}
    sqrt1pm1 := x/(1.0+sqrt(1.0+x));
  end;
end;


{---------------------------------------------------------------------------}
function sqrt1pmx(x: extended): extended;
  {-Return sqrt(1+x^2)-x}
begin
  if x <= 0.0 then sqrt1pmx := hypot(1.0,x)-x
  else sqrt1pmx := 1.0/(hypot(1.0,x)+x);
end;


{---------------------------------------------------------------------}
{---------------------- Floating point functions ---------------------}
{---------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function copysign(x,y: extended): extended; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return abs(x)*sign(y)}
begin
  THexExtW(x)[4] := (THexExtW(x)[4] and $7FFF) or (THexExtW(y)[4] and $8000);
  copysign := x;
end;


{---------------------------------------------------------------------------}
function copysignd(x,y: double): double;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return abs(x)*sign(y)}
begin
  THexDblW(x)[3] := (THexDblW(x)[3] and $7FFF) or (THexDblW(y)[3] and $8000);
  copysignd := x;
end;


{---------------------------------------------------------------------------}
function copysigns(x,y: single): single;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return abs(x)*sign(y)}
begin
  THexSglA(x)[3] := (THexSglA(x)[3] and $7F) or (THexSglA(y)[3] and $80);
  copysigns := x;
end;


{---------------------------------------------------------------------------}
procedure frexp(x: extended; var m: extended; var e: longint);
  {-Return the mantissa m and exponent e of x with x = m*2^e, 0.5 <= abs(m) < 1;}
  { if x is 0, +-INF, NaN, return m=x, e=0}
var
  xh: THexExtW absolute x;  {x as array of word}
const
  H2_64: THexExtW = ($0000,$0000,$0000,$8000,$403f);  {2^64}
begin
  e := xh[4] and $7FFF;
  {First check is INF or NAN}
  if (e=$7FFF) or (x=0.0) then e := 0
  else begin
    if e=0 then begin
      {denormal}
      x := x*extended(H2_64);
      e := xh[4] and $7FFF;
      dec(e,64+$3FFE);
    end
    else dec(e,$3FFE);
    xh[4] := (xh[4] and $8000) or $3FFE;
  end;
  m := x;
end;


{---------------------------------------------------------------------------}
procedure frexpd(d: double; var m: double; var e: longint);
  {-Return the mantissa m and exponent e of d with d = m*2^e, 0.5 <= abs(m) < 1;}
  { if d is 0, +-INF, NaN, return m=d, e=0}
var
  w: integer;
const
  H2_54: THexDblW = ($0000,$0000,$0000,$4350);  {2^54}
begin
  w := THexDblW(d)[3] and $7FF0;
  e := 0;
  if (w=$7FF0) or (d=0) then begin
    {+-INF, NaN, 0 or denormal}
    m := d;
  end
  else begin
    if w < $0010 then begin
      d := d*double(H2_54);
      e := -54;
      w := THexDblW(d)[3] and $7FF0;
    end;
    inc(e, (w shr 4) - 1022);
    THexDblW(d)[3] := (THexDblW(d)[3] and $800F) or $3FE0;
    m := d;
  end;
end;


{---------------------------------------------------------------------------}
procedure frexps(s: single; var m: single; var e: longint);
  {-Return the mantissa m and exponent e of s with s = m*2^e, 0.5 <= abs(m) < 1;}
  { if s is 0, +-INF, NaN, return m=s, e=0}
var
  L: longint absolute s;
  x: double;
begin
  e := L and $7F800000;
  if (e = $7F800000) or (s=0) then begin
    {+-INF, NaN, 0}
    m := s;
    e := 0;
  end
  else begin
    frexpd(s,x,e);
    m := x;
  end;
end;


{---------------------------------------------------------------------------}
function ilogb(x: extended): longint;
  {-Return base 2 exponent of x. For finite x ilogb = floor(log2(|x|))}
  { otherwise -MaxLongint for x=0 and MaxLongint if x = +-INF or Nan. }
var
  e: integer;
  f: longint;
  m: extended;
begin
  e := THexExtW(x)[4] and $7FFF;
  if e=$7FFF then ilogb := MaxLongint
  else if x=0.0 then ilogb := -MaxLongint
  else if e<>0 then  ilogb := e-$3FFF
  else begin
    {e=0, x<>0: denormal use frexp exponent}
    frexp(x,m,f);
    ilogb := f-1;
  end;
end;


{---------------------------------------------------------------------------}
function IsInf(x: extended): boolean;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x is +INF or -INF}
begin
  with TExtRec(x) do begin
    IsInf := (xp and $7FFF=$7FFF) and (hm=longint($80000000)) and (lm=0);
  end;
end;


{---------------------------------------------------------------------------}
function IsInfD(d: double): boolean;   {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is +INF or -INF}
begin
  with TDblRec(d) do begin
    IsInfD := (hm and $7FFFFFFF=$7FF00000) and (lm=0);
  end;
end;


{---------------------------------------------------------------------------}
function IsNaN(x: extended): boolean;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x is a NaN}
begin
  with TExtRec(x) do begin
    IsNaN := (xp and $7FFF=$7FFF) and ((hm<>longint($80000000)) or (lm<>0));
  end;
end;


{---------------------------------------------------------------------------}
function IsNaND(d: double): boolean;   {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is a NaN}
begin
  with TDblRec(d) do begin
    IsNaND := (hm and $7FF00000=$7FF00000) and ((hm and $000FFFFF<>0) or (lm<>0));
  end;
end;


{---------------------------------------------------------------------------}
function IsNaNorInf(x: extended): boolean;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x is a NaN or infinite}
begin
  IsNaNorInf := TExtRec(x).xp and $7FFF=$7FFF;
end;


{---------------------------------------------------------------------------}
function IsNaNorInfD(d: double): boolean;   {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is a NaN or infinite}
begin
  IsNaNorInfD := THexDblW(d)[3] and $7FF0=$7FF0;
end;


{---------------------------------------------------------------------------}
function IsInfS(s: single): boolean;        {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is +INF or -INF}
var
  L: longint absolute s;
begin
  IsInfS := L and $7FFFFFFF = $7F800000;
end;


{---------------------------------------------------------------------------}
function IsNaNS(s: single): boolean;        {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is a NaN}
var
  L: longint absolute s;
begin
  IsNaNS := (L and $7F800000 = $7F800000) and  (L and $7FFFFF <> 0);
end;


{---------------------------------------------------------------------------}
function IsNaNorInfS(s: single): boolean;   {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is a NaN or infinite}
var
  L: longint absolute s;
begin
  IsNaNorInfS := L and $7F800000 = $7F800000;
end;


{---------------------------------------------------------------------------}
function isRMNearest: boolean;
  {-Check if round_to_nearest without FPU control}
const
  small: extended = 1e-30;
var
  x,y: extended;
begin
  x := 1.0 + small;
  y := 1.0 - small;
  isRMNearest := x=y;
end;


{---------------------------------------------------------------------------}
function predd(d: double): double;
  {-Return next representable double after d in the direction -Inf}
begin
  with TDblRec(d) do begin
    if THexDblW(d)[3] and $7FF0=$7FF0 then begin
      {Inf or Nan}
      if (hm and $7FFFFFFF=$7FF00000) and (lm=0) then begin
        {d is +- Inf}
        if d>0.0 then d := MaxDouble;
      end;
    end
    else begin
      {finite number}
      if d=0.0 then begin
        hm := x80000000;
        lm := 1;
      end
      else if d<0.0 then begin
        {d<0: increment significand}
        inc(lm);
        if lm=0 then inc(hm);
      end
      else begin
        {d>0: decrement significand}
        if lm=0 then dec(hm);
        dec(lm);
      end;
    end;
    predd := d;
  end;
end;


{---------------------------------------------------------------------------}
function predx(x: extended): extended;
  {-Return next representable extended after x in the direction -Inf}
begin
  with TExtRec(x) do begin
    if xp and $7FFF=$7FFF then begin
      {Inf or Nan}
      if (hm=x80000000) and (lm=0) then begin
        {x is +- Inf}
        if x>0.0 then x := MaxExtended;
      end;
      predx := x;
      exit;
    end;

    if xp and $7FFF = 0 then begin
      {Handle pseudo-denormal: Set exponent to +/- 1, significand is unchanged}
      if hm<0 then xp := xp or 1;
    end
    else if hm>=0 then begin
      {don't touch unnormals}
      predx := x;
      exit;
    end;

    {finite number}
    if x=0 then begin
      xp := $8000;
      lm := 1;
    end
    else if xp and $8000 <> 0 then begin
      {x<0: increment significand}
      inc(lm);
      if lm=0 then begin
        inc(hm);
        if (hm=0) or ((xp=$8000) and (hm=x80000000))  then begin
          inc(xp);
          hm := hm or x80000000;
          if xp=$FFFF then x := NegInf_x;
        end;
      end;
    end
    else begin
      {x>0: decrement significand}
      if hm<0 then begin
        if lm=0 then begin
          if hm=x80000000 then begin
            dec(xp);
            dec(hm);
            if xp>0 then hm := hm or x80000000;
          end
          else dec(hm);
        end;
        dec(lm);
      end
      else begin
        {denormal}
        if lm=0 then dec(hm);
        dec(lm);
      end;
    end;
  end;
  predx := x;
end;


{---------------------------------------------------------------------------}
function preds(s: single): single;
  {-Return next representable single after s in the direction -Inf}
var
  L: longint absolute s;
begin
  if L and $7F800000 = $7F800000 then begin
    {Inf or Nan, don't change Nan or -Inf}
    if L and $7FFFFF = 0 then begin
      {s is +- Inf}
      if s>0.0 then s := MaxSingle;
    end;
  end
  else begin
    {finite number}
    if s=0.0 then L := longint($80000001)
    else if s<0.0 then inc(L)
    else dec(L);
  end;
  preds := s;
end;


{---------------------------------------------------------------------------}
function succd(d: double): double;
  {-Return next representable double after d in the direction +Inf}
begin
  with TDblRec(d) do begin
    if THexDblW(d)[3] and $7FF0=$7FF0 then begin
      {Inf or Nan}
      if (hm and $7FFFFFFF=$7FF00000) and (lm=0) then begin
        {d is +- Inf}
        if d<0.0 then d := -MaxDouble;
      end;
    end
    else begin
      {finite number}
      if d=0.0 then begin
        hm := 0;
        lm := 1;
      end
      else if d>0.0 then begin
        {d>0: increment significand}
        inc(lm);
        {hm < $7FF00000, so inc(hm) cannot overflow and will give}
        {the correct result succd(predd(PosInf_d)) = PosInf_d}
        if lm=0 then inc(hm);
      end
      else begin
        {d<0: decrement significand}
        if lm=0 then dec(hm);
        dec(lm);
      end;
    end;
    succd := d;
  end;
end;


{---------------------------------------------------------------------------}
function succx(x: extended): extended;
  {-Return next representable extended after x in the direction +Inf}
begin
  with TExtRec(x) do begin

    if xp and $7FFF=$7FFF then begin
      {Inf or Nan}
      if (hm=x80000000) and (lm=0) then begin
        {x is +- Inf}
        if x<0.0 then x := -MaxExtended;
      end;
      succx := x;
      exit;
    end;
    if xp and $7FFF = 0 then begin
      {Handle pseudo-denormal: Set exponent to +/- 1, significand is unchanged}
      if hm<0 then xp := xp or 1;
    end
    else if hm>=0 then begin
      {don't touch unnormals}
      succx := x;
      exit;
    end;

    {finite number}
    if x=0.0 then begin
      xp := 0;
      lm := 1;
    end
    else if xp and $8000 = 0 then begin
      {x>0: increment significand}
      inc(lm);
      if lm=0 then begin
        inc(hm);
        if (hm=0) or ((xp=0) and (hm=x80000000)) then begin
          inc(xp);
          hm := hm or x80000000;
          if xp=$7FFF then x := PosInf_x;
        end;
      end;
    end
    else begin
      {x<0: decrement significand}
      if lm=0 then begin
        if (hm>=0) or (hm=x80000000)  then begin
          dec(hm);
          dec(xp);
          if xp and $7FFF > 0 then hm := hm or x80000000;
        end
        else dec(hm);
      end;
      dec(lm);
    end;
  end;
  succx := x;
end;


{---------------------------------------------------------------------------}
function succs(s: single): single;
  {-Return next representable single after s in the direction +Inf}
var
  L: longint absolute s;
begin
  if L and $7F800000 = $7F800000 then begin
    {Inf or Nan, don't change Nan or +Inf}
    if L and $7FFFFF = 0 then begin
      {s is +- Inf}
      if s<0.0 then s := -MaxSingle;
    end;
  end
  else begin
    {finite number}
    if s=0.0 then L := 1
    else if s>0.0 then inc(L)
    else dec(L);
  end;
  succs := s;
end;


{---------------------------------------------------------------------------}
function ulpd(d: double): double;
  {-Return the 'unit in the last place': ulpd(d)=|d|-predd(|d|) for finite d}
begin
  if THexDblW(d)[3] and $7FF0=$7FF0 then with TDblRec(d) do begin
    {Inf or Nan}
    if (hm and $7FFFFFFF=$7FF00000) and (lm=0) then ulpd := PosInf_d
    else ulpd := d;
  end
  else begin
    d := abs(d);
    {Note if d=0 then ulpd(0) = 0 - predd(d) = succd(0) !}
    ulpd := d - predd(d);
  end;
end;


{---------------------------------------------------------------------------}
function ulps(s: single): single;
  {-Return the 'unit in the last place': ulps(s)=|s|-preds(|s|) for finite s}
var
  L: longint absolute s;
begin
  if L and $7F800000 = $7F800000 then begin
    {Inf or Nan}
    if L and $7FFFFF = 0 then ulps := PosInf_s
    else ulps := s;
  end
  else begin
    s := abs(s);
    {Note if s=0 then ulps(0) = 0 - preds(s) = succs(0) !}
    ulps := s - preds(s);
  end;
end;


{---------------------------------------------------------------------------}
function ulpx(x: extended): extended;
  {-Return the 'unit in the last place': ulpx(x)=|x|-predx(|x|) for finite x}
begin
  if TExtRec(x).xp and $7FFF=$7FFF then with TExtRec(x) do begin
    {Inf or Nan}
    if (hm=longint($80000000)) and (lm=0) then ulpx := PosInf_x
    else ulpx := x;
  end
  else begin
    x := abs(x);
    {Note if x=0 then ulpx(0) = 0 - predx(x) = succx(0) !}
    ulpx := x - predx(x);
  end;
end;


{---------------------------------------------------------------------------}
function rint(x: extended): extended;
  {-Return the integral value nearest x for the current rounding mode}
const
  twopower63h: THexExtW = ($0000,$0000,$0000,$8000,$403E); {2^63}
var
  c: extended absolute twopower63h;
begin
  if TExtRec(x).xp and $7FFF >= $4040 then begin
    {x is either INF or NAN or has no fractional part}
    rint := x  {Inf or NAN}
  end
  else begin
    if x>=0.0 then begin
      x := x + c;
      rint := x - c;
    end
    else begin
      x := x - c;
      rint := x + c;
    end;
  end;
end;


{---------------------------------------------------------------------------}
{------------------------- FPU control functions ---------------------------}
{---------------------------------------------------------------------------}

{These functions try to be as compatible to Delphi/FreePascal as reasonable.}
{16 bit Pascal cannot return sets as function results, Delphi sets are not  }
{byte-compatible in FPC, not all Delphi versions support default values etc.}
{D3+ and FPC2+ have Set8087CW in system.pas. It is used here if available,  }
{because it also sets the system variable Default8087CW. D6+ and FPC+ have  }
{Get8087CW in system.pas.}

{$ifdef need_get87}
  {----------------------------------------------------------------}
  function Get8087CW: word;
    {-Return the FPU control word}
  var
    cw: word;
  begin
    asm
      fstcw [cw]
      fwait
    end;
    Get8087CW := cw;
  end;
{$endif}

{$ifdef need_set87}
  {----------------------------------------------------------------}
  procedure Set8087CW(cw: word);
    {-Set new FPU control word}
  begin
    asm
      fnclex
      fldcw  [cw]
      fwait
    end;
  end;
{$endif}


{---------------------------------------------------------------------------}
function GetRoundMode: TFPURoundingMode;
  {-Return the current rounding mode}
begin
  GetRoundMode := TFPURoundingMode((Get8087CW shr 10) and 3);
end;


{---------------------------------------------------------------------------}
function SetRoundMode(NewRoundMode: TFPURoundingMode): TFPURoundingMode;
  {-Set new rounding mode and return the old mode}
var
  CW: word;
begin
  CW := Get8087CW;
  Set8087CW((CW and $F3FF) or (ord(NewRoundMode) shl 10));
  SetRoundMode := TFPURoundingMode((CW shr 10) and 3);
end;


{---------------------------------------------------------------------------}
function GetPrecisionMode: TFPUPrecisionMode;
  {-Return the current precision control mode}
begin
  GetPrecisionMode := TFPUPrecisionMode((Get8087CW shr 8) and 3);
end;


{---------------------------------------------------------------------------}
function SetPrecisionMode(NewPrecision: TFPUPrecisionMode): TFPUPrecisionMode;
  {-Set new precision control mode and return the old precision}
var
  CW: word;
begin
  CW := Get8087CW;
  Set8087CW((CW and $FCFF) or (ord(NewPrecision) shl 8));
  SetPrecisionMode := TFPUPrecisionMode((CW shr 8) and 3);
end;


{---------------------------------------------------------------------------}
procedure GetExceptionMask(var Mask: byte);
  {-Return the current exception mask}
begin
  Mask := Get8087CW and $3F;
end;


{---------------------------------------------------------------------------}
procedure SetExceptionMask(NewMask: byte);
  {-Set new exception mask}
var
  CW: word;
begin
  CW := Get8087CW;
  Set8087CW((CW and $FFC0) or (NewMask and $3F));
end;


(*************************************************************************
 DESCRIPTION     :  Accurate argument reduction mod Pi/2
 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     17.11.09  W.Ehrhardt  Initial BP7 k_rem_pio2
 0.11     28.11.09  we          rem_pio2_ph, rem_pio2_cw, rem_pio2
 0.12     28.11.09  we          NOBASM version
 0.13     28.11.09  we          k_rem_pio2: make i,j integer, new t: longint
 0.14     29.11.09  we          use abs(x) in rem_pio2_ph; comments/references
 0.15     29.11.09  we          PIo2: array[0..7] of THexDblW
 0.16     03.12.09  we          used ldexp for scalbnd
***************************************************************************)

{--------------------------------------------------------------------------}
{------------------  Cody and Waite argument reduction ---------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function rem_pio2_cw(x: extended; var z: extended): integer;
  {-Cody/Waite reduction of x: z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8}
var
  i: longint;
  y: extended;
const
  HP1 : THexExtW = ($0000,$0000,$da80,$c90f,$3ffe);
  HP2 : THexExtW = ($0000,$0000,$a300,$8885,$3fe4);
  HP3 : THexExtW = ($3707,$a2e0,$3198,$8d31,$3fc8);
var
  DP1 : extended absolute HP1; {= 7.853981554508209228515625E-1;}
  DP2 : extended absolute HP2; {= 7.94662735614792836713604629E-9;}
  DP3 : extended absolute HP3; {= 3.06161699786838294306516483E-17;}
begin
  {This is my Pascal translation of the CW reduction given in}
  {sinl.c from the Cephes Math Library Release 2.7: May, 1998}
  {Copyright 1985, 1990, 1998 by Stephen L. Moshier          }
  if x=0.0 then begin
    z := 0.0;
    rem_pio2_cw := 0;
  end
  else begin
    y := floorx(x/Pi_4);
    i := trunc(y - 16.0*floorx(y/16.0));
    if odd(i) then begin
      inc(i);
      y := y + 1.0;
    end;
    rem_pio2_cw := (i shr 1) and 7;
    {Extended precision modular arithmetic}
    z := ((x - y * DP1) - y * DP2) - y * DP3;
  end;
end;


{---------------------------------------------------------------------------}
{-----------------  Payne and Hanek argument reduction ---------------------}
{---------------------------------------------------------------------------}

{k_rem_pio2 is my Pascal translation of the C function __kernel_rem_pio2}

(*
 * @(#)k_rem_pio2.c 5.1 93/09/24 */
 *
 * ====================================================
 * Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
 *
 * Developed at SunPro, a Sun Microsystems, Inc. business.
 * Permission to use, copy, modify, and distribute this
 * software is freely granted, provided that this notice
 * is preserved.
 * ====================================================
 *)

(* k_rem_pio2 return the last three bits of N with y = x - N*pi/2
 * so that |y| < pi/2.
 *
 * The method is to compute the integer (mod 8) and fraction parts of
 * (2/pi)*x without doing the full multiplication. In general we
 * skip the part of the product that are known to be a huge integer
 * (more accurately, = 0 mod 8 ). Thus the number of operations are
 * independent of the exponent of the input.
 *
 * (2/pi) is represented by an array of 24-bit integers in ipio2[].
 *
 * Input parameters:
 *      x[]     The input value (must be positive) is broken into nx
 *              pieces of 24-bit integers in double precision format.
 *              x[i] will be the i-th 24 bit of x. The scaled exponent
 *              of x[0] is given in input parameter e0 (i.e., x[0]*2^e0
 *              match x's up to 24 bits.
 *
 *              Example of breaking a double positive z into x[0]+x[1]+x[2]:
 *                      e0 = ilogb(z)-23
 *                      z  = scalbn(z,-e0)
 *              for i = 0,1,2
 *                      x[i] = floor(z)
 *                      z    = (z-x[i])*2**24
 *
 *
 *      y[]     output result in an array of double precision numbers.
 *              The dimension of y[] is:
 *                      24-bit  precision       1
 *                      53-bit  precision       2
 *                      64-bit  precision       2
 *                      113-bit precision       3
 *              The actual value is the sum of them. Thus for 113-bit
 *              precison, one may have to do something like:
 *
 *              long double t,w,r_head, r_tail;
 *              t = (long double)y[2] + (long double)y[1];
 *              w = (long double)y[0];
 *              r_head = t+w;
 *              r_tail = w - (r_head - t);
 *
 *      e0      The exponent of x[0]. Must be <= 16360 or you need to
 *              expand the ipio2 table.
 *
 *      nx      dimension of x[]
 *
 *      prec    an integer indicating the precision:
 *                      0       24  bits (single)
 *                      1       53  bits (double)
 *                      2       64  bits (extended)
 *                      3       113 bits (quad)
 *
 * Here is the description of some local variables:
 *
 *      jk      jk+1 is the initial number of terms of ipio2[] needed
 *              in the computation. The recommended value is 2,3,4,
 *              6 for single, double, extended,and quad.
 *
 *      jz      local integer variable indicating the number of
 *              terms of ipio2[] used.
 *
 *      jx      nx - 1
 *
 *      jv      index for pointing to the suitable ipio2[] for the
 *              computation. In general, we want
 *                      ( 2^e0*x[0] * ipio2[jv-1]*2^(-24jv) )/8
 *              is an integer. Thus
 *                      e0-3-24*jv >= 0 or (e0-3)/24 >= jv
 *              Hence jv = max(0,(e0-3)/24).
 *
 *      jp      jp+1 is the number of terms in PIo2[] needed, jp = jk.
 *
 *      q[]     double array with integral value, representing the
 *              24-bits chunk of the product of x and 2/pi.
 *
 *      q0      the corresponding exponent of q[0]. Note that the
 *              exponent for q[i] would be q0-24*i.
 *
 *      PIo2[]  double precision array, obtained by cutting pi/2
 *              into 24 bits chunks.
 *
 *      f[]     ipio2[] in floating point
 *
 *      iq[]    integer array by breaking up q[] in 24-bits chunk.
 *
 *      fq[]    final product of x*(2/pi) in fq[0],..,fq[jk]
 *
 *      ih      integer. If >0 it indicates q[] is >= 0.5, hence
 *              it also indicates the *sign* of the result.
 *
 *)

{PIo2[] double array, obtained by cutting pi/2 into 24 bits chunks.}
const
  PIo2: array[0..7] of THexDblW = (
          ($0000, $4000, $21FB, $3FF9),  {1.5707962512969971    }
          ($0000, $0000, $442D, $3E74),  {7.5497894158615964e-08}
          ($0000, $8000, $4698, $3CF8),  {5.3903025299577648e-15}
          ($0000, $6000, $CC51, $3B78),  {3.2820034158079130e-22}
          ($0000, $8000, $1B83, $39F0),  {1.2706557530806761e-29}
          ($0000, $4000, $2520, $387A),  {1.2293330898111133e-36}
          ($0000, $8000, $8222, $36E3),  {2.7337005381646456e-44}
          ($0000, $0000, $F31D, $3569)); {2.1674168387780482e-51}


{ipio2: 16560 bits of 2/pi, i.e. the big-endian integer trunc((2/pi)*2^16560}
{This amount of bits is needed because extended exponents have a maximum of }
{about 2^14. For double-ranges only about 1200 bits are needed.}
const
  ipio2: array[0..689] of longint = (
           $A2F983, $6E4E44, $1529FC, $2757D1, $F534DD, $C0DB62,
           $95993C, $439041, $FE5163, $ABDEBB, $C561B7, $246E3A,
           $424DD2, $E00649, $2EEA09, $D1921C, $FE1DEB, $1CB129,
           $A73EE8, $8235F5, $2EBB44, $84E99C, $7026B4, $5F7E41,
           $3991D6, $398353, $39F49C, $845F8B, $BDF928, $3B1FF8,
           $97FFDE, $05980F, $EF2F11, $8B5A0A, $6D1F6D, $367ECF,
           $27CB09, $B74F46, $3F669E, $5FEA2D, $7527BA, $C7EBE5,
           $F17B3D, $0739F7, $8A5292, $EA6BFB, $5FB11F, $8D5D08,
           $560330, $46FC7B, $6BABF0, $CFBC20, $9AF436, $1DA9E3,
           $91615E, $E61B08, $659985, $5F14A0, $68408D, $FFD880,
           $4D7327, $310606, $1556CA, $73A8C9, $60E27B, $C08C6B,
           {the following bits are needed for extended exponents only}
           $47C419, $C367CD, $DCE809, $2A8359, $C4768B, $961CA6,
           $DDAF44, $D15719, $053EA5, $FF0705, $3F7E33, $E832C2,
           $DE4F98, $327DBB, $C33D26, $EF6B1E, $5EF89F, $3A1F35,
           $CAF27F, $1D87F1, $21907C, $7C246A, $FA6ED5, $772D30,
           $433B15, $C614B5, $9D19C3, $C2C4AD, $414D2C, $5D000C,
           $467D86, $2D71E3, $9AC69B, $006233, $7CD2B4, $97A7B4,
           $D55537, $F63ED7, $1810A3, $FC764D, $2A9D64, $ABD770,
           $F87C63, $57B07A, $E71517, $5649C0, $D9D63B, $3884A7,
           $CB2324, $778AD6, $23545A, $B91F00, $1B0AF1, $DFCE19,
           $FF319F, $6A1E66, $615799, $47FBAC, $D87F7E, $B76522,
           $89E832, $60BFE6, $CDC4EF, $09366C, $D43F5D, $D7DE16,
           $DE3B58, $929BDE, $2822D2, $E88628, $4D58E2, $32CAC6,
           $16E308, $CB7DE0, $50C017, $A71DF3, $5BE018, $34132E,
           $621283, $014883, $5B8EF5, $7FB0AD, $F2E91E, $434A48,
           $D36710, $D8DDAA, $425FAE, $CE616A, $A4280A, $B499D3,
           $F2A606, $7F775C, $83C2A3, $883C61, $78738A, $5A8CAF,
           $BDD76F, $63A62D, $CBBFF4, $EF818D, $67C126, $45CA55,
           $36D9CA, $D2A828, $8D61C2, $77C912, $142604, $9B4612,
           $C459C4, $44C5C8, $91B24D, $F31700, $AD43D4, $E54929,
           $10D5FD, $FCBE00, $CC941E, $EECE70, $F53E13, $80F1EC,
           $C3E7B3, $28F8C7, $940593, $3E71C1, $B3092E, $F3450B,
           $9C1288, $7B20AB, $9FB52E, $C29247, $2F327B, $6D550C,
           $90A772, $1FE76B, $96CB31, $4A1679, $E27941, $89DFF4,
           $9794E8, $84E6E2, $973199, $6BED88, $365F5F, $0EFDBB,
           $B49A48, $6CA467, $427271, $325D8D, $B8159F, $09E5BC,
           $25318D, $3974F7, $1C0530, $010C0D, $68084B, $58EE2C,
           $90AA47, $02E774, $24D6BD, $A67DF7, $72486E, $EF169F,
           $A6948E, $F691B4, $5153D1, $F20ACF, $339820, $7E4BF5,
           $6863B2, $5F3EDD, $035D40, $7F8985, $295255, $C06437,
           $10D86D, $324832, $754C5B, $D4714E, $6E5445, $C1090B,
           $69F52A, $D56614, $9D0727, $50045D, $DB3BB4, $C576EA,
           $17F987, $7D6B49, $BA271D, $296996, $ACCCC6, $5414AD,
           $6AE290, $89D988, $50722C, $BEA404, $940777, $7030F3,
           $27FC00, $A871EA, $49C266, $3DE064, $83DD97, $973FA3,
           $FD9443, $8C860D, $DE4131, $9D3992, $8C70DD, $E7B717,
           $3BDF08, $2B3715, $A0805C, $93805A, $921110, $D8E80F,
           $AF806C, $4BFFDB, $0F9038, $761859, $15A562, $BBCB61,
           $B989C7, $BD4010, $04F2D2, $277549, $F6B6EB, $BB22DB,
           $AA140A, $2F2689, $768364, $333B09, $1A940E, $AA3A51,
           $C2A31D, $AEEDAF, $12265C, $4DC26D, $9C7A2D, $9756C0,
           $833F03, $F6F009, $8C402B, $99316D, $07B439, $15200C,
           $5BC3D8, $C492F5, $4BADC6, $A5CA4E, $CD37A7, $36A9E6,
           $9492AB, $6842DD, $DE6319, $EF8C76, $528B68, $37DBFC,
           $ABA1AE, $3115DF, $A1AE00, $DAFB0C, $664D64, $B705ED,
           $306529, $BF5657, $3AFF47, $B9F96A, $F3BE75, $DF9328,
           $3080AB, $F68C66, $15CB04, $0622FA, $1DE4D9, $A4B33D,
           $8F1B57, $09CD36, $E9424E, $A4BE13, $B52333, $1AAAF0,
           $A8654F, $A5C1D2, $0F3F0B, $CD785B, $76F923, $048B7B,
           $721789, $53A6C6, $E26E6F, $00EBEF, $584A9B, $B7DAC4,
           $BA66AA, $CFCF76, $1D02D1, $2DF1B1, $C1998C, $77ADC3,
           $DA4886, $A05DF7, $F480C6, $2FF0AC, $9AECDD, $BC5C3F,
           $6DDED0, $1FC790, $B6DB2A, $3A25A3, $9AAF00, $9353AD,
           $0457B6, $B42D29, $7E804B, $A707DA, $0EAA76, $A1597B,
           $2A1216, $2DB7DC, $FDE5FA, $FEDB89, $FDBE89, $6C76E4,
           $FCA906, $70803E, $156E85, $FF87FD, $073E28, $336761,
           $86182A, $EABD4D, $AFE7B3, $6E6D8F, $396795, $5BBF31,
           $48D784, $16DF30, $432DC7, $356125, $CE70C9, $B8CB30,
           $FD6CBF, $A200A4, $E46C05, $A0DD5A, $476F21, $D21262,
           $845CB9, $496170, $E0566B, $015299, $375550, $B7D51E,
           $C4F133, $5F6E13, $E4305D, $A92E85, $C3B21D, $3632A1,
           $A4B708, $D4B1EA, $21F716, $E4698F, $77FF27, $80030C,
           $2D408D, $A0CD4F, $99A520, $D3A2B3, $0A5D2F, $42F9B4,
           $CBDA11, $D0BE7D, $C1DB9B, $BD17AB, $81A2CA, $5C6A08,
           $17552E, $550027, $F0147F, $8607E1, $640B14, $8D4196,
           $DEBE87, $2AFDDA, $B6256B, $34897B, $FEF305, $9EBFB9,
           $4F6A68, $A82A4A, $5AC44F, $BCF82D, $985AD7, $95C7F4,
           $8D4D0D, $A63A20, $5F57A4, $B13F14, $953880, $0120CC,
           $86DD71, $B6DEC9, $F560BF, $11654D, $6B0701, $ACB08C,
           $D0C0B2, $485551, $0EFB1E, $C37295, $3B06A3, $3540C0,
           $7BDC06, $CC45E0, $FA294E, $C8CAD6, $41F3E8, $DE647C,
           $D8649B, $31BED9, $C397A4, $D45877, $C5E369, $13DAF0,
           $3C3ABA, $461846, $5F7555, $F5BDD2, $C6926E, $5D2EAC,
           $ED440E, $423E1C, $87C461, $E9FD29, $F3D6E7, $CA7C22,
           $35916F, $C5E008, $8DD7FF, $E26A6E, $C6FDB0, $C10893,
           $745D7C, $B2AD6B, $9D6ECD, $7B723E, $6A11C6, $A9CFF7,
           $DF7329, $BAC9B5, $5100B7, $0DB2E2, $24BA74, $607DE5,
           $8AD874, $2C150D, $0C1881, $94667E, $162901, $767A9F,
           $BEFDFD, $EF4556, $367ED9, $13D9EC, $B9BA8B, $FC97C4,
           $27A831, $C36EF1, $36C594, $56A8D8, $B5A8B4, $0ECCCF,
           $2D8912, $34576F, $89562C, $E3CE99, $B920D6, $AA5E6B,
           $9C2A3E, $CC5F11, $4A0BFD, $FBF4E1, $6D3B8E, $2C86E2,
           $84D4E9, $A9B4FC, $D1EEEF, $C9352E, $61392F, $442138,
           $C8D91B, $0AFC81, $6A4AFB, $D81C2F, $84B453, $8C994E,
           $CC2254, $DC552A, $D6C6C0, $96190B, $B8701A, $649569,
           $605A26, $EE523F, $0F117F, $11B5F4, $F5CBFC, $2DBC34,
           $EEBC34, $CC5DE8, $605EDD, $9B8E67, $EF3392, $B817C9,
           $9B5861, $BC57E1, $C68351, $103ED8, $4871DD, $DD1C2D,
           $A118AF, $462C21, $D7F359, $987AD9, $C0549E, $FA864F,
           $FC0656, $AE79E5, $362289, $22AD38, $DC9367, $AAE855,
           $382682, $9BE7CA, $A40D51, $B13399, $0ED7A9, $480569,
           $F0B265, $A7887F, $974C88, $36D1F9, $B39221, $4A827B,
           $21CF98, $DC9F40, $5547DC, $3A74E1, $42EB67, $DF9DFE,
           $5FD45E, $A4677B, $7AACBA, $A2F655, $23882B, $55BA41,
           $086E59, $862A21, $834739, $E6E389, $D49EE5, $40FB49,
           $E956FF, $CA0F1C, $8A59C5, $2BFA94, $C5C1D3, $CFC50F,
           $AE5ADB, $86C547, $624385, $3B8621, $94792C, $876110,
           $7B4C2A, $1A2C80, $12BF43, $902688, $893C78, $E4C4A8,
           $7BDBE5, $C23AC4, $EAF426, $8A67F7, $BF920D, $2BA365,
           $B1933D, $0B7CBD, $DC51A4, $63DD27, $DDE169, $19949A,
           $9529A8, $28CE68, $B4ED09, $209F44, $CA984E, $638270,
           $237C7E, $32B90F, $8EF5A7, $E75614, $08F121, $2A9DB5,
           $4D7E6F, $5119A5, $ABF9B5, $D6DF82, $61DD96, $023616,
           $9F3AC4, $A1A283, $6DED72, $7A8D39, $A9B882, $5C326B,
           $5B2746, $ED3400, $7700D2, $55F4FC, $4D5901, $8071E0);

const
  init_jk: array[0..3] of integer = (2,3,4,6); {initial value for jk}

const
  two24:  double = 16777216.0;             {2^24}
  twon24: double = 5.9604644775390625e-08; {1/2^24}

type
  TDA02 = array[0..2] of double; {Type definition needed for Pascal versions}
                                 {without open arrays, 2 is OK for extended!}

{---------------------------------------------------------------------------}
function k_rem_pio2(const x: TDA02; var y: TDA02; e0, nx, prec: integer): integer;
 {-Calculate y with y = x - n*Pi/2, |y| <= Pi/4 and return n mod 8}
label
  recompute;
var
  i,ih,j,jz,jx,jv,jp,jk,carry,k,m,n,q0: integer;
  t: longint; {WE: used for longint calculations, all other C-ints can be integers}
  iq: array[0..19] of longint;
  f,fq,q: array[0..19] of double;
  z,fw: double;
begin

{$ifdef FPC}
  {Suppress warnings: Local variable does not seem to be initialized}
  f[0]  := 0.0;
  fq[0] := 0.0;
  q[0]  := 0.0;
{$endif}

  {initialize jk}
  jk := init_jk[prec];
  jp := jk;

  {determine jx,jv,q0, note that 3>q0}
  jx := nx-1;
  jv := (e0-3) div 24; if jv<0 then jv := 0;
  q0 := e0-24*(jv+1);

  {set up f[0] to f[jx+jk] where f[jx+jk] = ipio2[jv+jk]}
  j := jv-jx;
  m := jx+jk;
  for i:=0 to m do begin
    if j<0 then f[i] := 0.0 else f[i] := ipio2[j];
    inc(j);
  end;

  {compute q[0],q[1],...q[jk]}
  for i:=0 to jk do begin
    fw := 0.0;
    for j:=0 to jx do begin
      fw   := fw + x[j]*f[jx+i-j];
      q[i] := fw;
    end;
  end;
  jz := jk;

recompute:

  {distill q[] into iq[] reversingly}
  i := 0;
  j := jz;
  z := q[jz];
  while j>0 do begin
    fw    := trunc(twon24*z);
    iq[i] := trunc(z-two24*fw);
    z     := q[j-1]+fw;
    inc(i);
    dec(j);
  end;

  {compute n}
  z  := ldexp(z,q0);              {actual value of z}
  z  := z - 8.0*floorx(z*0.125);  {trim off integer >= 8}
  n  := trunc(z);
  z  := z - n;
  ih := 0;
  if q0>0 then begin
    {need iq[jz-1] to determine n}
    t  := (iq[jz-1] shr (24-q0));
    inc(n,t);
    dec(iq[jz-1], t shl (24-q0));
    ih := iq[jz-1] shr (23-q0);
  end
  else if q0=0 then ih := iq[jz-1] shr 23
  else if z>=0.5 then ih := 2;

  if ih>0 then begin
    {q > 0.5}
    inc(n);
    carry := 0;
    for i:=0 to jz-1 do begin
      {compute 1-q}
      t := iq[i];
      if carry=0 then begin
        if t<>0 then begin
          carry := 1;
          iq[i] := $1000000 - t;
        end
      end
      else iq[i] := $ffffff - t;
    end;
    if q0>0 then begin
      {rare case: chance is 1 in 12}
      case q0 of
        1: iq[jz-1] := iq[jz-1] and $7fffff;
        2: iq[jz-1] := iq[jz-1] and $3fffff;
      end;
    end;
    if ih=2 then begin
      z := 1.0 - z;
      if carry<>0 then  z := z - ldexp(1.0,q0);
    end;
  end;

  {check if recomputation is needed}
  if z=0.0 then begin
    t := 0;
    for i:=jz-1 downto jk do t := t or iq[i];
    if t=0 then begin
      {need recomputation}
      k := 1;
      while iq[jk-k]=0 do inc(k);   {k = no. of terms needed}
      for i:=jz+1 to jz+k do begin
        {add q[jz+1] to q[jz+k]}
        f[jx+i] := ipio2[jv+i];
        fw := 0.0;
        for j:=0 to jx do fw := fw + x[j]*f[jx+i-j];
        q[i] := fw;
      end;
      inc(jz,k);
      goto recompute;
    end;
  end;

  {chop off zero terms}
  if z=0.0 then begin
    dec(jz);
    dec(q0,24);
    while iq[jz]=0 do begin
      dec(jz);
      dec(q0,24);
    end;
  end
  else begin
    {break z into 24-bit if necessary}
    z := ldexp(z,-q0);
    if z>=two24 then begin
      fw := trunc(twon24*z);
      iq[jz] := trunc(z-two24*fw);
      inc(jz);
      inc(q0,24);
      iq[jz] := trunc(fw);
    end
    else iq[jz] := trunc(z);
  end;

  {convert integer "bit" chunk to floating-point value}
  fw := ldexp(1.0,q0);
  for i:=jz downto 0 do begin
    q[i] := fw*iq[i];
    fw := fw*twon24;
  end;

  {compute PIo2[0,...,jp]*q[jz,...,0]}
  for i:=jz downto 0 do begin
    fw :=0.0;
    k := 0;
    while (k<=jp) and (k<=jz-i) do begin
      fw := fw + double(PIo2[k])*(q[i+k]);
      fq[jz-i] := fw;
      inc(k);
    end;
  end;

  {compress fq[] into y[]}
  case prec of
     0:  begin
           fw := 0.0;
           for i:=jz downto 0 do fw := fw + fq[i];
           if ih=0 then y[0] := fw else y[0] := -fw;
         end;
     1,
     2:  begin
           fw := 0.0;
           for i:=jz downto 0 do fw := fw + fq[i];
           if ih=0 then y[0] := fw else y[0] := -fw;
           fw := fq[0]-fw;
           for i:=1 to jz do fw := fw + fq[i];
           if ih=0 then y[1] := fw else y[1] := -fw;
         end;
     3:  begin
           {painful}
           for i:=jz downto 1 do begin
             fw     := fq[i-1]+fq[i];
             fq[i]  := fq[i] + (fq[i-1]-fw);
             fq[i-1]:= fw;
           end;
           for i:=jz downto 2 do begin
             fw      := fq[i-1]+fq[i];
             fq[i]  := fq[i] + (fq[i-1]-fw);
             fq[i-1]:= fw;
           end;
           fw := 0.0;
           for i:=jz downto 2 do fw := fw + fq[i];
           if ih=0 then begin
             y[0] :=  fq[0];
             y[1] :=  fq[1];
             y[2] :=  fw;
           end
           else begin
             y[0] := -fq[0];
             y[1] := -fq[1];
             y[2] := -fw;
           end;
         end;
  end;
  k_rem_pio2 := n and 7;
end;


{---------------------------------------------------------------------------}
function rem_pio2_ph(x: extended; var z: extended): integer;
  {-Payne/Hanek reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8}
var
  ax, ay: TDA02;
  e0: integer;
begin
  z := abs(x);
  if IsNanOrInf(x) or (x=0.0) then rem_pio2_ph := 0
  else begin
    e0 := ilogb(z)-23;
    z  := ldexp(z,-e0);
    ax[0] := trunc(z); z := ldexp(z-ax[0],24);
    ax[1] := trunc(z);
    ax[2] := trunc(ldexp(z-ax[1],24));
    e0 := k_rem_pio2(ax,ay, e0, 3, 2);
    if x>=0 then begin
      z := ay[0]+ay[1];
      rem_pio2_ph := e0;
    end
    else begin
      z := -ay[0]-ay[1];
      rem_pio2_ph := (-e0) and 7;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function rem_pio2(x: extended; var z: extended): integer;
  {-Argument reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8.}
  { Uses Payne/Hanek if |x| >= ph_cutoff, Cody/Waite otherwise}
const
  tol: extended = 5.9604644775390625e-8;  {ph_cutoff*eps_x}
begin
  z := abs(x);
  if z<Pi_4 then begin
    z := x;
    rem_pio2 := 0;
  end
  else if z > ph_cutoff then rem_pio2 := rem_pio2_ph(x,z)
  else begin
    rem_pio2 := rem_pio2_cw(x,z);
    if abs(z) <= tol then begin
      {If x is close to a multiple of Pi/2, the C/W relative error may be}
      {large, e.g. 1.0389e-6 for x = 491299012647073665.0/16777216.0 with}
      {z_cw = 5.0672257144385922E-20 and z_ph = 5.0672204500144216E-20.  }
      {In this case redo the calculation with the Payne/Hanek algorithm. }
      rem_pio2 := rem_pio2_ph(x,z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function rem_2pi(x: extended): extended;
  {-Return x mod 2*Pi}
var
  {$ifndef ExtendedSyntax_on} n: integer;{$endif}
  z: extended;
const
  {Reduction constants: 2*Pi = Pi2H + Pi2M + Pi2L}
  H2PH : THexExtW = ($0000,$0000,$da80,$c90f,$4001);
  H2PM : THexExtW = ($0000,$0000,$a300,$8885,$3fe7);
  H2PL : THexExtW = ($3707,$a2e0,$3198,$8d31,$3fcb);
var
  Pi2H: extended absolute H2PH; {= 6.28318524360656738}
  Pi2M: extended absolute H2PM; {= 6.35730188491834269E-8}
  Pi2L: extended absolute H2PL; {= 2.44929359829470635E-16}
begin
  if abs(x) <= ph_cutoff then begin
    {Direct Cody/Waite style reduction. This is more efficient}
    {than calling rem_pio2_cw with additional adjustment.}
    z := floorx(x/TwoPi);
    z := ((x - z*Pi2H) - z*Pi2M) - z*Pi2L;
    if z>TwoPi then begin
      {May be reached due to rounding, but was never observed}
      z := ((z - Pi2H) - Pi2M) - Pi2L;
    end;
  end
  else begin
    {Use Payne/Hanek code with adjustments for mod 2Pi}
    {$ifndef ExtendedSyntax_on} n:={$endif} rem_pio2_ph(0.25*x,z);
    z := 4.0*z;
    {Here |z| <= Pi, adjustment for z<0 is done below}
  end;
  {Note that in rare cases result=TwoPi, e.g. for x = -TwoPi.}
  {All z with -2eps_x <= z < 0 will have z + TwoPi = TwoPi}
  if z<0.0 then z := ((z + Pi2H) + Pi2M) + Pi2L;
  rem_2pi := z;
end;


{---------------------------------------------------------------------------}
function rem_2pi_sym(x: extended): extended;
  {-Return x mod 2*Pi, -Pi <= result <= Pi}
var
  {$ifndef ExtendedSyntax_on} n: integer;{$endif}
  z: extended;
begin
  {$ifndef ExtendedSyntax_on} n:={$endif} rem_pio2(0.25*x,z);
  rem_2pi_sym := 4.0*z;
end;


{---------------------------------------------------------------------------}
function rem_int2(x: extended; var z: extended): integer;
  {-Argument reduction of x: z*Pi = x*Pi - n*Pi/2, |z|<=1/4, result = n mod 8.}
  { Used for argument reduction in sin(Pi*x) and cos(Pi*x)}
var
  y: extended;
  i: integer;
begin
  if IsNanOrInf(x) or (abs(x)<=0.25) then begin
    rem_int2 := 0;
    z := x;
    exit;
  end;
  if frac(x)=0.0 then begin
    {Here x is an integer or abs(x) >= 2^64}
    z := 0.0;
    i := 0;
    {set i=2, if x is a odd}
    if (TExtRec(x).xp and $7FFF < $403F) and (frac(0.5*x)<>0.0) then i:=2;
  end
  else begin
    {Here x is not an integer. First calculate x mod 2,}
    {this leaves Pi*x = Pi*(x mod 2) mod 2*Pi invariant}
    x := 0.5*x;
    x := 2.0*(x-floorx(x));
    {then apply the Cody/Waite style range reduction}
    y := floorx(4.0*x);
    i := trunc(y - 16.0*floorx(y/16.0));
    if odd(i) then begin
      inc(i);
      y := y + 1.0;
    end;
    i := (i shr 1) and 7;
    z := x-0.25*y;
  end;
  rem_int2 := i;
end;


{---------------------------------------------------------------------------}
{-------------- Polynomial, Vector, Statistic Operations -------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function sum2(const a: array of double; n: integer): extended;
  {-Compute accurate sum(a[i], i=0..n-1) of a double vector}
var
  e,s,t,x,y,z: double;
  i: integer;
begin

  {empty sum}
  if n<=0 then begin
    sum2 := 0.0;
    exit;
  end;
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('sum2:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;

  {Double version of [8] Algorithm 4.4}
  s := a[0];
  e := 0.0;
  for i:=1 to n-1 do begin
    t := a[i];
    {(x,y) = TwoSum(a[i],s)}
    x := t+s;
    z := x-t;
    y := (t-(x-z)) + (s-z);
    {sum of errors}
    e := e+y;
    s := x;
  end;
  sum2 := s+e;
end;


{---------------------------------------------------------------------------}
function sum2x(const a: array of extended; n: integer): extended;
  {-Compute accurate sum(a[i], i=0..n-1) of extended vector}
var
  e,s,t,x,y,z: extended;
  i: integer;
begin

  {empty sum}
  if n<=0 then begin
    sum2x := 0.0;
    exit;
  end;

  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('sum2x:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;

  {Extended version of [8] Algorithm 4.4}
  s := a[0];
  e := 0.0;
  for i:=1 to n-1 do begin
    t := a[i];
    {(x,y) = TwoSum(t,s)}
    x := t+s;
    z := x-t;
    y := (t-(x-z)) + (s-z);
    {sum of errors}
    e := e+y;
    s := x;
  end;
  sum2x := s+e;
end;


{---------------------------------------------------------------------------}
function dot2(const x,y: array of double; n: integer): extended;
  {-Accurate dot product sum(x[i]*y[i], i=0..n-1) of two double vectors}
const
  csd = 134217729.0; {2^27+1}  {Split constant for TwoProduct}
var
  i: integer;
  p,s,h,q,r: double;
  x1,x2,y1,y2: double;
begin

  {empty sum}
  if n<=0 then begin
    dot2 := 0.0;
    exit;
  end;

  if n>high(x)+1 then n:=high(x)+1;
  if n>high(y)+1 then n:=high(y)+1;

  {Double version of [8], Algorithm 5.3. Note that the original version}
  {in [8] initializes (p,s) = TwoProduct(x[0], x[0]). This saves a few }
  {flops compared with my code, but TwoProduct source either has to be }
  {doubled or called as a subroutine (which is normally slower)}
  s := 0.0;
  p := 0.0;
  for i:=0 to n-1 do begin
    {TwoProduct(x[i],y[i],h,r);}
      q  := x[i];
      {split x[i] into x1,x2}
      r  := csd*q;
      x2 := r-q;
      x1 := r-x2;
      x2 := q-x1;
      r  := y[i];
      {h = x[i]*y[i]}
      h  := q*r;
      {split y into y1,y2}
      q  := csd*r;
      y2 := q-r;
      y1 := q-y2;
      y2 := r-y1;
      {r = x2*y2 - (((h-x1*y1)-x2*y1)-x1*y2}
      q  := x1*y1;
      q  := h-q;
      y1 := y1*x2;
      q  := q-y1;
      x1 := x1*y2;
      q  := q-x1;
      x2 := x2*y2;
      r  := x2-q;
    {(p,q) = TwoSum(p,h)}
      x1 := p+h;
      x2 := x1-p;
      y1 := x1-x2;
      y2 := h-x2;
      q  := p-y1;
      q  := q+y2;
      p  := x1;
    {s = s + (q+r))}
      q := q+r;
      s := s+q;
  end;
  dot2 := p+s;
end;


{---------------------------------------------------------------------------}
function dot2x(const x,y: array of extended; n: integer): extended;
  {-Accurate dot product sum(x[i]*y[i], i=0..n-1) of two extended vectors}
const
  csx = 4294967297.0; {2^32+1}  {Split constant for TwoProduct}
var
  i: integer;
  p,s,h,q,r: extended;
  x1,x2,y1,y2: extended;
begin

  {empty sum}
  if n<=0 then begin
    dot2x := 0.0;
    exit;
  end;

  if n>high(x)+1 then n:=high(x)+1;
  if n>high(y)+1 then n:=high(y)+1;

  {Extended version of [8], Algorithm 5.3. Note that the original version}
  {in [8] initializes (p,s) = TwoProduct(x[0], x[0]). This saves a few }
  {flops compared with my code, but TwoProduct source either has to be }
  {doubled or called as a subroutine (which is normally slower)}
  s := 0.0;
  p := 0.0;
  for i:=0 to n-1 do begin
    {TwoProdD(x[i],y[i],h,r);}
      q  := x[i];
      {split x[i] into x1,x2}
      r  := csx*q;
      x2 := r-q;
      x1 := r-x2;
      x2 := q-x1;
      r  := y[i];
      {h = x[i]*y[i]}
      h  := q*r;
      {split y into y1,y2}
      q  := csx*r;
      y2 := q-r;
      y1 := q-y2;
      y2 := r-y1;
      {r = x2*y2 - (((h-x1*y1)-x2*y1)-x1*y2}
      q  := x1*y1;
      q  := h-q;
      y1 := y1*x2;
      q  := q-y1;
      x1 := x1*y2;
      q  := q-x1;
      x2 := x2*y2;
      r  := x2-q;
    {(p,q) = TwoSum(p,h)}
      x1 := p+h;
      x2 := x1-p;
      y1 := x1-x2;
      y2 := h-x2;
      q  := p-y1;
      q  := q+y2;
      p  := x1;
    {s = s + (q+r))}
      q := q+r;
      s := s+q;
  end;
  dot2x := p+s;
end;


{---------------------------------------------------------------------------}
function rmsx(const a: array of extended; n: integer): extended;
  {-Calculate the RMS value sqrt(sum(a[i]^2, i=0..n-1)/n) of an extended vector}
var
  scale, sumsq: extended;
begin
  ssqx(a,n,scale,sumsq);
  rmsx := scale*sqrt(sumsq/n);
end;


{---------------------------------------------------------------------------}
function rms(const a: array of double; n: integer): extended;
  {-Calculate the RMS value sqrt(sum(a[i]^2, i=0..n-1)/n) of a double vector}
var
  scale, sumsq: extended;
begin
  ssqd(a,n,scale,sumsq);
  rms := scale*sqrt(sumsq/n);
end;


{---------------------------------------------------------------------------}
procedure ssqx(const a: array of extended; n: integer; var scale, sumsq: extended);
  {-Calculate sum(a[i]^2, i=0..n-1) = scale^2*sumsq, scale>=0, sumsq>0}
var
  i: integer;
  t: extended;
begin
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('ssqx:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;
  {based on http://www.netlib.org/lapack/explore-html/dlassq.f.html}
  {see also Higham [9], Problem 27.5 for a proof of correctness}
  scale := 0.0;
  sumsq := 1.0;
  if n<=1 then begin
    if n=1 then scale := abs(a[0]);
    exit;
  end;
  for i:=0 to n-1 do begin
    t := abs(a[i]);
    if t<>0 then begin
      if scale < t then begin
        sumsq := 1.0 + sumsq*sqr(scale/t);
        scale := t;
      end
      else sumsq := sumsq + sqr(t/scale);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure ssqd(const a: array of double; n: integer; var scale, sumsq: extended);
  {-Calculate sum(a[i]^2, i=0..n-1) = scale^2*sumsq, scale>=0, sumsq>0}
var
  i: integer;
  t: extended;
begin
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('ssqd:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;
  {based on http://www.netlib.org/lapack/explore-html/dlassq.f.html}
  {see also Higham [9], Problem 27.5 for a proof of correctness}
  scale := 0.0;
  sumsq := 1.0;
  if n<=1 then begin
    if n=1 then scale := abs(a[0]);
    exit;
  end;
  for i:=0 to n-1 do begin
    t := abs(a[i]);
    if t<>0 then begin
      if scale < t then begin
        sumsq := 1.0 + sumsq*sqr(scale/t);
        scale := t;
      end
      else sumsq := sumsq + sqr(t/scale);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function sumsqrx(const a: array of extended; n: integer): extended;
  {-Calculate sum(a[i]^2, i=0..n-1) of an extended vector}
var
  scale, sumsq: extended;
begin
  ssqx(a,n,scale,sumsq);
  sumsqrx := sqr(scale)*sumsq;
end;


{---------------------------------------------------------------------------}
function norm2x(const a: array of extended; n: integer): extended;
  {-Calculate the 2-norm = sqrt(sum(a[i]^2, i=0..n-1)) of an extended vector}
var
  scale, sumsq: extended;
begin
  ssqx(a,n,scale,sumsq);
  norm2x := scale*sqrt(sumsq);
end;


{---------------------------------------------------------------------------}
function sumsqr(const a: array of double; n: integer): extended;
  {-Calculate sum(a[i]^2, i=0..n-1) of a double vector}
var
  scale, sumsq: extended;
begin
  ssqd(a,n,scale,sumsq);
  sumsqr := sqr(scale)*sumsq;
end;


{---------------------------------------------------------------------------}
function norm2(const a: array of double; n: integer): extended;
  {-Calculate the 2-norm = sqrt(sum(a[i]^2, i=0..n-1)) of a double vector}
var
  scale, sumsq: extended;
begin
  ssqd(a,n,scale,sumsq);
  norm2 := scale*sqrt(sumsq);
end;


{---------------------------------------------------------------------------}
function PolEval(x: extended; const a: array of double; n: integer): extended;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}
var
  i: integer;
  s: extended;
begin
  if n<=0 then begin
    PolEval := 0.0;
    exit;
  end;
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('PolEval:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;
  s := a[n-1];
  for i:=n-2 downto 0 do s := s*x + a[i];
  PolEval := s;
end;


{---------------------------------------------------------------------------}
procedure PolEvalC(const a: array of double; n: integer; x,y: double; var u,v: double);
  {-Evaluate polynomial a[0] + a[1]*z + ... + a[n-1]*z^(n-1)}
  { for complex z = x + i*y, result is u + i*v}
var
  uj,vj,r,s,t: extended;
  j: integer;
begin
  if y=0.0 then begin
    {use real mode function}
    v := 0.0;
    u := PolEval(x,a,n);
    exit;
  end;
  if n<=1 then begin
    v := 0.0;
    if n=1 then u := a[0]
    else u := 0.0;
    exit;
  end;
  if n>high(a)+1 then n := high(a)+1;
  {Use the procedure from  Knuth [32], 4.6.4 (3). This saves about}
  {2n multiplications and n additions compared to the old version.}
  uj := a[n-1];
  vj := a[n-2];
  if n>2 then begin
    r := 2.0*x;
    s := sqr(x) + sqr(y);
    for j:=n-3 downto 0 do begin
      t  := vj   + r*uj;
      vj := a[j] - s*uj;
      uj := t;
    end;
  end;
  u := x*uj + vj;
  v := y*uj;
end;


{---------------------------------------------------------------------------}
procedure PolEvalDeriv(x: double; const a: array of double; n: integer; var px,dp: extended);
  {-Evaluate polynomial p(x) a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}
  { Return px = p(x) and dp = p'(x)}
var
  i: integer;
begin
  px := 0.0;
  dp := 0.0;
  if n<=0 then exit;
  if n>high(a)+1 then n := high(a)+1;
  {compute px = p(x) and dp = p'(x)}
  px := a[n-1];
  for i:= n-2 downto 0 do begin
    dp := dp*x + px;
    px := px*x + a[i];
  end;
end;


{---------------------------------------------------------------------------}
function PolEvalX(x: extended; const a: array of extended; n: integer): extended;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}
var
  i: integer;
  s: extended;
begin
  if n<=0 then begin
    PolEvalX := 0.0;
    exit;
  end;
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('PolEvalX:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;
  s := a[n-1];
  for i:=n-2 downto 0 do s := s*x + a[i];
  PolEvalX := s;
end;


{---------------------------------------------------------------------------}
function PolEvalS(x: extended; const a: array of single; n: integer): extended;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}
var
  i: integer;
  s: extended;
begin
  if n<=0 then begin
    PolEvalS := 0.0;
    exit;
  end;
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('PolEvalS:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;
  s := a[n-1];
  for i:=n-2 downto 0 do s := s*x + a[i];
  PolEvalS := s;
end;


{---------------------------------------------------------------------------}
function CSEvalX(x: extended; const a: array of extended; n: integer): extended;
  {-Evaluate Chebyshev sum a[0]/2 + a[1]*T_1(x) +..+ a[n-1]*T_(n-1)(x) using Clenshaw algorithm}
var
  b0,b1,b2: extended;
  i: integer;
begin
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('CSEvalX:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;
  b2 := 0.0;
  b1 := 0.0;
  b0 := 0.0;
  x  := 2.0*x;
  for i:=n-1 downto 0 do begin
    b2 := b1;
    b1 := b0;
    b0 := x*b1 - b2 + a[i];
  end;
  CSEvalX := 0.5*(b0-b2);
end;


{---------------------------------------------------------------------------}
procedure CSEvalXDer(x: extended; const a: array of extended; n: integer; var px,dp: extended);
  {-Evaluate Chebyshev sum p(x) = a[0]/2 + a[1]*T_1(x) +..+ a[n-1]*T_(n-1)(x),}
  { using Clenshaw algorithm, return pd = p(x) and dp = p'(x) }
var
  x2,b0,b1,b2,c0,c1,c2: extended;
  i: integer;
begin
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('CSEvalXDer:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;
  b0 := 0.0;
  b1 := 0.0;
  b2 := 0.0;
  c0 := 0.0;
  c1 := 0.0;
  c2 := 0.0;
  x2 := 4.0*x*x - 2.0;
  for i:=n-1 downto 0 do begin
    b2 := b1;
    b1 := b0;
    b0 := a[i] - b2 + x2 * b1;
    if 0 < i  then begin
      c2 := c1;
      c1 := c0;
      c0 := b0 - c2 + x2 * c1;
    end;
  end;
  dp := c0 - c2;
  px := 0.5*(b0 - b2);
end;


{---------------------------------------------------------------------------}
function PolEvalEE(x: double; const a: array of double; n: integer; var e: double): double;
  {-Evaluate polynomial; return p(x) = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.   }
var
  i: integer;
  s,z: double;
begin
  {empty sum}
  if n<=0 then begin
    PolEvalEE := 0.0;
    e := 0.0;
    exit;
  end;
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('PolEvalEE:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;
  {Double version of Algorithm 5.1 from Higham [9]}
  s := a[n-1];
  z := abs(x);
  e := 0.5*abs(s);
  for i:=n-2 downto 0 do begin
    s := s*x + a[i];
    e := e*z + abs(s);
  end;
  PolEvalEE := s;
  {Note: With round to nearest, u=eps_d/2 can be used}
  e := eps_d*abs(2.0*e-abs(s));
end;


{---------------------------------------------------------------------------}
function PolEvalEEX(x: extended; const a: array of extended; n: integer; var e: extended): extended;
  {-Evaluate polynomial; return p(x) = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.   }
var
  i: integer;
  s,z: extended;
begin
  {empty sum}
  if n<=0 then begin
    PolEvalEEX := 0.0;
    e := 0.0;
    exit;
  end;
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('PolEvalEEX:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;
  {Extended version of Algorithm 5.1 from Higham [9]}
  s := a[n-1];
  z := abs(x);
  e := 0.5*abs(s);
  for i:=n-2 downto 0 do begin
    s := s*x + a[i];
    e := e*z + abs(s);
  end;
  PolEvalEEX := s;
  {Note: With round to nearest, u=eps_x/2 can be used}
  e := eps_x*abs(2.0*e-abs(s));
end;


{---------------------------------------------------------------------------}
procedure PolEvalCHEDer(x: double; const a: array of double; n: integer; var px,dp,e: double);
  {-Evaluate polynomial; return p = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { accurate double precision version using compensated Horner scheme, }
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.}
  { dp is the non-compensated derivative p'(x).                        }
const
  csd = 134217729.0; {2^27+1}  {Split constant for TwoProduct}
var
  i: integer;
  p,q,s,t,h,u,w,z: double;
  sh,sl,xh,xl: double;
begin
  {empty sum}
  dp := 0.0;
  if n<=0 then begin
    px := 0.0;
    e  := 0.0;
    exit;
  end;

  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('PolEvalCHE:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;

  {Combined double version of Graillat et al [46] algorithms 7,8, and 9}
  s := a[n-1];
  {[xh,xl] = Split(x)}
  z  := x*csd;
  xh := z - x;
  xh := z - xh;
  xl := x - xh;
  {Initialize h = Horner sum for q+t}
  h := 0.0;
  {Initialize variables for error bound}
  u := 0.0;
  w := abs(x);
  for i:=n-2 downto 0 do begin
    dp := dp*x + s;
    {[sh,sl] = Split(s)}
    z  := s*csd;
    sh := z - s;
    sh := z - sh;
    sl := s - sh;
    {[p,q] = TwoProduct(s,x)}
    p  := s*x;
    q  := sh*xh;
    q  := p - q;
    z  := sl*xh;
    q  := q - z;
    z  := sh*xl;
    q  := q - z;
    z  := sl*xl;
    q  := z - q;
    {[s,t] = TwoSum(p,a[i])};
    s := p + a[i];
    z := s - p;
    t := s - z;
    t := p - t;
    z := a[i] - z;
    t := t + z;
    {Horner sum for q+t}
    h := h*x;
    z := q + t;
    h := h + z;
    {Horner sum for error bound}
    u := u*w + (abs(q)+abs(t));
  end;
  z := h + s;
  px := z;
  {Compute e with formula (15) of Graillat et al [46]}
  {Note: With round to nearest, u=eps_d/2 can be used}
  e := abs(z)*(1.0 + 2.0*eps_d);
  e := e + 4.0*(n-0.5)*u;
  e := eps_d*e;
end;


{---------------------------------------------------------------------------}
function PolEvalCHE(x: double; const a: array of double; n: integer; var e: double): double;
  {-Evaluate polynomial; return p(x) = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { accurate double precision version using compensated Horner scheme,    }
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.   }
var
  px,dp: double;
begin
  PolEvalCHEDer(x,a,n,px,dp,e);
  PolEvalCHE := px;
end;


{---------------------------------------------------------------------------}
function PolEvalCHEX(x: extended; const a: array of extended; n: integer; var e: extended): extended;
  {-Evaluate polynomial; return p(x) = a[0] + a[1]*x +...+ a[n-1]*x^(n-1);}
  { accurate extended precision version using compensated Horner scheme,  }
  { e is the dynamic absolute error estimate with |p(x) - result| <= e.   }
const
  csx = 4294967297.0; {2^32+1}  {Split constant for TwoProduct}
var
  i: integer;
  p,q,s,t,h,u,w,z: extended;
  sh,sl,xh,xl: extended;
begin
  {empty sum}
  if n<=0 then begin
    PolEvalCHEX := 0.0;
    e := 0.0;
    exit;
  end;
  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('PolEvalCHEX:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;

  {Combined extended version of Graillat et al [46] algorithms 7,8, and 9}
  s := a[n-1];
  {[xh,xl] = Split(x)}
  z  := x*csx;
  xh := z - x;
  xh := z - xh;
  xl := x - xh;
  {Initialize h = Horner sum for q+t}
  h := 0.0;
  {Initialize variables for error bound}
  u := 0.0;
  w := abs(x);
  for i:=n-2 downto 0 do begin
    {[sh,sl] = Split(s)}
    z  := s*csx;
    sh := z - s;
    sh := z - sh;
    sl := s - sh;
    {[p,q] = TwoProduct(s,x)}
    p  := s*x;
    q  := sh*xh;
    q  := p - q;
    z  := sl*xh;
    q  := q - z;
    z  := sh*xl;
    q  := q - z;
    z  := sl*xl;
    q  := z - q;
    {[s,t] = TwoSum(p,a[i])};
    s := p + a[i];
    z := s - p;
    t := s - z;
    t := p - t;
    z := a[i] - z;
    t := t + z;
    {Horner sum for q+t}
    h := h*x;
    z := q + t;
    h := h + z;
    {Horner sum for error bound}
    u := u*w + (abs(q)+abs(t));
  end;
  z := h + s;
  PolEvalCHEX := z;
  {Compute e with formula (15) of Graillat et al [46]}
  {Note: With round to nearest, u=eps_x/2 can be used}
  e := abs(z)*(1.0 + 2.0*eps_x);
  e := e + 4.0*(n-0.5)*u;
  e := eps_x*e;
end;


{---------------------------------------------------------------------------}
function mean(const a: array of double; n: integer): extended;
  {-Compute accurate mean = sum(a[i], i=0..n-1)/n of a double vector}
begin
  if n<=0 then mean := 0.0
  else mean := sum2(a,n)/n;
end;


{---------------------------------------------------------------------------}
function meanx(const a: array of extended; n: integer): extended;
  {-Compute accurate mean = sum(a[i], i=0..n-1)/n of an extended vector}
begin
  if n<=0 then meanx := 0.0
  else meanx := sum2x(a,n)/n;
end;


{---------------------------------------------------------------------------}
procedure mssqx(const a: array of extended; n: integer; var mval, scale, sumsq: extended);
  {-Calculate mean mval and ssqx sum((a[i]-mval)^2) of an extended vector}
var
  i: integer;
  t: extended;
begin

  {empty sum}
  scale := 0.0;
  sumsq := 1.0;
  if n<=0 then begin
    mval := 0.0;
    exit;
  end;

  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('mssqx:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;

  if n=1 then mval := a[0]
  else begin
    mval  := sum2x(a,n)/n;
    {calculate sum((a[i]-mval)^2), cf. ssqx}
    for i:=0 to n-1 do begin
      t := abs(mval-a[i]);
      if t<>0 then begin
        if scale < t then begin
          sumsq := 1.0 + sumsq*sqr(scale/t);
          scale := t;
        end
        else sumsq := sumsq + sqr(t/scale);
      end;
    end;
  end
end;


{---------------------------------------------------------------------------}
procedure mssqd(const a: array of double; n: integer; var mval, scale, sumsq: extended);
  {-Calculate mean mval and ssqd sum((a[i]-mval)^2) of a double vector}
var
  i: integer;
  t: extended;
begin

  {empty sum}
  scale := 0.0;
  sumsq := 1.0;
  if n<=0 then begin
    mval := 0.0;
    exit;
  end;

  {$ifdef debug}
    if n>high(a)+1 then begin
      writeln('mssqd:  n > high(a)+1, n = ',n, ' vs. ', high(a)+1);
      readln;
    end;
  {$endif}
  if n>high(a)+1 then n := high(a)+1;

  if n=1 then mval := a[0]
  else begin
    mval  := sum2(a,n)/n;
    {calculate sum((a[i]-mval)^2), cf. ssqd}
    for i:=0 to n-1 do begin
      t := abs(mval-a[i]);
      if t<>0 then begin
        if scale < t then begin
          sumsq := 1.0 + sumsq*sqr(scale/t);
          scale := t;
        end
        else sumsq := sumsq + sqr(t/scale);
      end;
    end;
  end
end;


{---------------------------------------------------------------------------}
procedure MeanAndStdDevX(const a: array of extended; n: integer; var mval, sdev: extended);
  {-Accurate mean and sample standard deviation of an extended vector}
var
  scale, sumsq: extended;
begin
  mssqx(a, n, mval, scale, sumsq);
  if n<2 then sdev := 0.0
  else sdev := scale*sqrt(sumsq/(n-1));
end;


{---------------------------------------------------------------------------}
procedure MeanAndStdDev(const a: array of double; n: integer; var mval, sdev: extended);
  {-Accurate mean and sample standard deviation of a double vector}
var
  scale, sumsq: extended;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then sdev := 0.0
  else sdev := scale*sqrt(sumsq/(n-1));
end;


{---------------------------------------------------------------------------}
function CoeffVar(const a: array of double; n: integer): extended;
  {-Return the coefficient of variation (sdev/mean) of a double vector}
var
  mval, sdev: extended;
begin
  MeanAndStdDev(a,n,mval,sdev);
  CoeffVar := SafeDivX(sdev,mval);
end;


{---------------------------------------------------------------------------}
function CoeffVarX(const a: array of extended; n: integer): extended;
  {-Return the coefficient of variation (sdev/mean) of an extended vector}
var
  mval, sdev: extended;
begin
  MeanAndStdDevX(a,n,mval,sdev);
  CoeffVarX := SafeDivX(sdev,mval);
end;


{---------------------------------------------------------------------------}
procedure moment(const a: array of double; n: integer; var m1, m2, m3, m4, skew, kurt: extended);
  {-Return the first 4 moments, skewness, and kurtosis of a double vector}
var
  q,d,scale,sumsq: extended;
  i: integer;
begin
  {Note: There many different notations for skewness and kurtosis, some scale}
  {with the sample and not the population standard deviation, some subtract 3}
  {from the kurtosis etc. These procedures are accurate versions of Delphi's }
  {MomentSkewKurtosis and virtually use the following Maple V R4 definitions:}

  { m1   := describe[moment[1,0]](a);    }
  { m2   := describe[moment[2,mean]](a); }
  { m3   := describe[moment[3,mean]](a); }
  { m4   := describe[moment[4,mean]](a); }
  { skew := describe[skewness](a);       }
  { kurt := describe[kurtosis](a);       }

  mssqd(a, n, m1, scale, sumsq);
  m2   := sqr(scale)*(sumsq/n);
  m3   := 0.0;
  m4   := 0.0;
  skew := 0.0;
  kurt := 0.0;
  if n<2 then exit;
  q := scale*sqrt(sumsq/n);

  {Compute the higher moments using recursion formulas. These are obviously}
  {more accurate than the Delphi formulas: try Delphi with (a,a+1,a+2,a+2) }
  {for a=10^k and see the disaster for k>=5! My routines are accurate up to}
  {k=15 (for k=16 the values are not exact doubles (and a+1=a!), but even  }
  {in this case moment gives correct results for the actual input vector!) }
  {m_k(i) = sum(d(j)^k, j=1..i)/i      = sum(d(j)^k,j=1..i-1)/i + d(i)^k/i }
  {       = (i-1)m_k(i-1)/i + d(i)^k/i = m_k(i-1) + (d(i)^k - m_k(i-1))/i  }
  for i:=1 to n do begin
    d  := a[i-1] - m1;
    m3 := m3 + (sqr(d)*d - m3)/i;
    m4 := m4 + (sqr(d*d) - m4)/i;
    {skew/kurt are defined only if m2<>0}
    if q<>0 then begin
      d := d/q;
      skew := skew + (sqr(d)*d - skew)/i;
      kurt := kurt + (sqr(d*d) - kurt)/i;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure momentx(const a: array of extended; n: integer; var m1, m2, m3, m4, skew, kurt: extended);
  {-Return the first 4 moments, skewness, and kurtosis of an extended vector}
var
  q,d,scale,sumsq: extended;
  i: integer;
begin
  {See procedure moment for explanations and notes}
  mssqx(a, n, m1, scale, sumsq);
  m2   := sqr(scale)*(sumsq/n);
  m3   := 0.0;
  m4   := 0.0;
  skew := 0.0;
  kurt := 0.0;
  if n<2 then exit;
  q := scale*sqrt(sumsq/n);
  {compute the higher moments using recursion formulas}
  for i:=1 to n do begin
    d  := a[i-1] - m1;
    m3 := m3 + (sqr(d)*d - m3)/i;
    m4 := m4 + (sqr(d*d) - m4)/i;
    {skew/kurt are defined only if m2<>0}
    if q<>0 then begin
      d := d/q;
      skew := skew + (sqr(d)*d - skew)/i;
      kurt := kurt + (sqr(d*d) - kurt)/i;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function ssdev(const a: array of double; n: integer): extended;
  {-Return the sample standard deviation of a double vector}
var
  mval,scale,sumsq: extended;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then ssdev := 0.0
  else ssdev := scale*sqrt(sumsq/(n-1));
end;


{---------------------------------------------------------------------------}
function ssdevx(const a: array of extended; n: integer): extended;
  {-Return the sample standard deviation of an extended vector}
var
  mval,scale,sumsq: extended;
begin
  mssqx(a, n, mval, scale, sumsq);
  if n<2 then ssdevx := 0.0
  else ssdevx := scale*sqrt(sumsq/(n-1));
end;


{---------------------------------------------------------------------------}
function psdev(const a: array of double; n: integer): extended;
  {-Return the population standard deviation of a double vector}
var
  mval,scale,sumsq: extended;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then psdev := 0.0
  else psdev := scale*sqrt(sumsq/n);
end;


{---------------------------------------------------------------------------}
function psdevx(const a: array of extended; n: integer): extended;
  {-Return the population standard deviation of an extended vector}
var
  mval,scale,sumsq: extended;
begin
  mssqx(a, n, mval, scale, sumsq);
  if n<2 then psdevx := 0.0
  else psdevx := scale*sqrt(sumsq/n);
end;


{---------------------------------------------------------------------------}
function svar(const a: array of double; n: integer): extended;
  {-Return the sample variance of a double vector}
var
  mval,scale,sumsq: extended;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then svar := 0.0
  else svar := sqr(scale)*(sumsq/(n-1));
end;


{---------------------------------------------------------------------------}
function svarx(const a: array of extended; n: integer): extended;
  {-Return the sample variance of an extended vector}
var
  mval,scale,sumsq: extended;
begin
  mssqx(a, n, mval, scale, sumsq);
  if n<2 then svarx := 0.0
  else svarx := sqr(scale)*(sumsq/(n-1));
end;


{---------------------------------------------------------------------------}
function pvar(const a: array of double; n: integer): extended;
  {-Return the population variance of a double vector}
var
  mval,scale,sumsq: extended;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then pvar := 0.0
  else pvar := sqr(scale)*(sumsq/n);
end;


{---------------------------------------------------------------------------}
function pvarx(const a: array of extended; n: integer): extended;
  {-Return the population variance of an extended vector}
var
  mval,scale,sumsq: extended;
begin
  mssqx(a, n, mval, scale, sumsq);
  if n<2 then pvarx := 0.0
  else pvarx := sqr(scale)*(sumsq/n);
end;


{---------------------------------------------------------------------------}
function tvar(const a: array of double; n: integer): extended;
  {-Return the total variance of a double vector}
var
  mval,scale,sumsq: extended;
begin
  mssqd(a, n, mval, scale, sumsq);
  tvar := sqr(scale)*(sumsq);
end;


{---------------------------------------------------------------------------}
function tvarx(const a: array of extended; n: integer): extended;
  {-Return the total variance of an extended vector}
var
  mval,scale,sumsq: extended;
begin
  mssqx(a, n, mval, scale, sumsq);
  tvarx := sqr(scale)*(sumsq);
end;


{---------------------------------------------------------------------------}
function pcov(const x,y: array of double; n: integer): extended;
  {-Calculate the population covariance = sum((x[i]-mx)(y[i]-my), i=0..n-1)/n of two double vectors}
var
  mx,my,s: extended;
  i: integer;
begin
  if n<=0 then pcov := 0.0
  else begin
    mx := mean(x,n);
    my := mean(y,n);
    s := 0.0;
    for i:=0 to n-1 do begin
      s := s + (x[i] - mx) * (y[i] - my);
    end;
    pcov := s/n;
  end;
end;


{---------------------------------------------------------------------------}
function pcovx(const x,y: array of extended; n: integer): extended;
  {-Calculate the population covariance = sum((x[i]-mx)(y[i]-my), i=0..n-1)/n of two extended vectors}
var
  mx,my,s: extended;
  i: integer;
begin
  if n<=0 then pcovx := 0.0
  else begin
    mx := meanx(x,n);
    my := meanx(y,n);
    s := 0.0;
    for i:=0 to n-1 do begin
      s := s + (x[i] - mx) * (y[i] - my);
    end;
    pcovx := s/n;
  end;
end;


{--------------------------------------------------------------------}
{------------------------- Other function ---------------------------}
{--------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function DegToRad(x: extended): extended;
  {-Convert angle x from degrees to radians}
begin
  DegToRad := x*Pi_180;
end;


{---------------------------------------------------------------------------}
function RadToDeg(x: extended): extended;
  {-Convert angle x from radians to degrees}
begin
  RadToDeg := x/Pi_180;
end;


{---------------------------------------------------------------------------}
function angle2(x1,x2,y1,y2: extended): extended;
  {-Return the accurate angle between the vectors (x1,x2) and (y1,y2)}
var
  d1,d2: extended;
begin
  {Ref: W. Kahan, Cross-Products and Rotations in Euclidean 2- and 3-Space,}
  {13, p.15. Online: http://www.cs.berkeley.edu/~wkahan/MathH110/Cross.pdf}
  {A uniformly accurate arctan() formula for the angle between two vectors:}
  d1 := hypot(x1,x2);
  d2 := hypot(y1,y2);
  if (d1=0.0) or (d2=0.0) then angle2 := 0.0
  else begin
    {x =  x/|x|}
    x1 := x1/d1;
    x2 := x2/d1;
    {y =  y/|y|}
    y1 := y1/d2;
    y2 := y2/d2;
    {d =  x/|x| - y/|y|}
    d1 := x1-y1;
    d2 := x2-y2;
    {x =  x/|x| + y/|y|}
    x1 := x1+y1;
    x2 := x2+y2;
    {y1 = |x/|x| - y/|y||}
    y1 := hypot(d1,d2);
    {y2 = |x/|x| + y/|y||}
    y2 := hypot(x1,x2);
    angle2 := 2.0*arctan2(y1,y2);
  end;
end;


{---------------------------------------------------------------------------}
function orient2d(x1,y1,x2,y2,x3,y3: extended): extended;
  {-Return the mathematical orientation of the three points (xi,yi): >0 if}
  { the order is counterclockwise, <0 if clockwise, =0 if they are collinear.}
  { Result is twice the signed area of the triangle defined by the points.}
var
  det,detl,detr: extended;
const
  eps = 3.2526065175E-19;  {= (3 + 16*eps_x)*eps_x}
var
  v1,v2: array[0..5] of extended;
begin
  {Ref: J.R. Shewchuk, Adaptive Precision Floating-Point Arithmetic and Fast}
  {Robust Geometric Predicates; cf. his public domain code (predicates.c)}
  {http://www-2.cs.cmu.edu/~quake/robust.html}

  detl := (x1 - x3)*(y2 - y3);
  detr := (y1 - y3)*(x2 - x3);
  det  := detl - detr;

  {Set default result and test if fast calculation can be used.}
  orient2d := det;
  if (detl=0.0) or (detr=0.0) then exit;

  if ((detl > 0.0) and (detr < 0.0)) or ((detl < 0.0) and (detr > 0.0)) then exit;

  if abs(det) >= eps*(abs(detl) + abs(detr)) then exit;

  {Fast calculation is not reliable and an accurate algorithm is needed}
  {to compute orient2 = x1*y2 - x1*y3 - y1*x2 + y1*x3 - x3*y2 + y3*x2. }
  {Shewchuk gives several robust formulas based on TwoSum & TwoProduct,}
  {but here I calculate orient2 as an accurate dot product using dot2x.}

  v1[0] := x1;    v2[0] :=  y2;
  v1[1] := x1;    v2[1] := -y3;
  v1[2] := y1;    v2[2] := -x2;
  v1[3] := y1;    v2[3] :=  x3;
  v1[4] := x3;    v2[4] := -y2;
  v1[5] := y3;    v2[5] :=  x2;

  orient2d := dot2x(v1,v2,6);
end;


{---------------------------------------------------------------------------}
function area_triangle(x1,y1,x2,y2,x3,y3: extended): extended;
  {-Return the area of the triangle defined by the points (xi,yi)}
begin
  area_triangle := 0.5*abs(orient2d(x1,y1,x2,y2,x3,y3));
end;


{---------------------------------------------------------------------------}
function in_triangle(x,y,x1,y1,x2,y2,x3,y3: extended): boolean;
  {-Return true if the point (x,y) lies strictly inside the triangle defined}
  { by the three points (xi,yi), false if it lies on a side or outside.}
var
  ori: integer;
begin
  {See e.g. http://www.mochima.com/articles/cuj_geometry_article/cuj_geometry_article.html}
  ori := isign(orient2d(x1,y1,x2,y2,x,y));
  if ori=0 then in_triangle := false
  else if ori <> isign(orient2d(x2,y2,x3,y3,x,y)) then in_triangle := false
  else in_triangle := ori=isign(orient2d(x3,y3,x1,y1,x,y));
end;


{---------------------------------------------------------------------------}
function in_triangle_ex(x,y,x1,y1,x2,y2,x3,y3: extended): integer;
  {-Return +1 if the point (x,y) lies strictly inside the triangle defined by}
  { the three points (xi,yi), -1 if it lies strictly outside, 0 otherwise.}
var
  ori1,ori2,ori3: integer;
begin
  in_triangle_ex := -1;

  ori1 := isign(orient2d(x1,y1,x2,y2,x,y));
  ori2 := isign(orient2d(x2,y2,x3,y3,x,y));
  if ori1*ori2 = -1 then exit;

  ori3 := isign(orient2d(x3,y3,x1,y1,x,y));
  if (ori1=ori2) and (ori2=ori3) then begin
    {all three orientations have the same value, inside if <> 0}
    if ori1<>0 then in_triangle_ex := 1
    else in_triangle_ex := 0;
  end
  else if (ori1=0) and ((ori2*ori3) >= 0) then in_triangle_ex := 0
  else if (ori2=0) and ((ori1*ori3) >= 0) then in_triangle_ex := 0
  else if (ori3=0) and ((ori1*ori2) >= 0) then in_triangle_ex := 0;
end;


{---------------------------------------------------------------------------}
function Ext2Dbl(x: extended): double;
  {-Return x as double, or +-Inf if too large}
begin
  if IsNanOrInf(x) then Ext2Dbl := x
  else if x>MaxDouble then Ext2Dbl := PosInf_d
  else if x < -MaxDouble then Ext2Dbl := NegInf_d
  else Ext2Dbl := x;
end;


{$ifndef BIT16}
  type char8=ansichar;
{$else}
  type char8=char;
{$endif}

const
  hnib: array[0..15] of char8 = '0123456789ABCDEF';

{---------------------------------------------------------------------------}
function Ext2Hex(x: extended): string;
  {-Return x as a big-endian hex string}
var
  ax : THexExtA absolute x;
  i,j: integer;
  s  : string[20];
begin
  {$ifndef BIT16}
    setlength(s,20);
  {$else}
    s[0] := #20;
  {$endif}
  j := 1;
  for i:=9 downto 0 do begin
    s[j] := hnib[ax[i] shr 4 ]; inc(j);
    s[j] := hnib[ax[i] and $f]; inc(j);
  end;
  Ext2Hex := {$ifdef D12Plus} string {$endif} (s);
end;


{---------------------------------------------------------------------------}
function Dbl2Hex(d: double): string;
  {-Return d as a big-endian hex string}
var
  ad : THexDblA absolute d;
  i,j: integer;
  s  : string[16];
begin
  {$ifndef BIT16}
    setlength(s,16);
  {$else}
    s[0] := #16;
  {$endif}
  j := 1;
  for i:=7 downto 0 do begin
    s[j] := hnib[ad[i] shr 4 ]; inc(j);
    s[j] := hnib[ad[i] and $f]; inc(j);
  end;
  Dbl2Hex := {$ifdef D12Plus} string {$endif} (s);
end;


{---------------------------------------------------------------------------}
function Sgl2Hex(s: single): string;
  {-Return s as a big-endian hex string}
var
  ad : THexSglA absolute s;
  i,j: integer;
  t  : string[8];
begin
  {$ifndef BIT16}
    setlength(t,8);
  {$else}
    t[0] := #8;
  {$endif}
  j := 1;
  for i:=3 downto 0 do begin
    t[j] := hnib[ad[i] shr 4 ]; inc(j);
    t[j] := hnib[ad[i] and $f]; inc(j);
  end;
  Sgl2Hex := {$ifdef D12Plus} string {$endif} (t);
end;


{---------------------------------------------------------------------------}
procedure Hex2Float(const hex: string; n: integer; var a: THexExtA; var code: integer);
  {-Common code for hex to float conversion, internal use only}
var
  i,j,m: integer;
  b: byte;
  c: char;
const
  c0  = ord('0');
  cal = ord('a') - 10;
  cau = ord('A') - 10;
begin
  if hex='' then code := -1
  else if (n=9) or (n=7) or (n=3) then begin
    if hex[1]='$' then m := 1 else m := 0;
    if length(hex)<>2*n+2+m then code := -2
    else begin
      b := 0;
      j := n;
      for i:=0 to 2*n+1 do begin
        inc(m);
        c := hex[m];
        if      (c>='0') and (c<='9') then b := (b shl 4) or ((ord(c)-c0 ) and $0F)
        else if (c>='A') and (c<='F') then b := (b shl 4) or ((ord(c)-cau) and $0F)
        else if (c>='a') and (c<='f') then b := (b shl 4) or ((ord(c)-cal) and $0F)
        else begin
          code := -4;
          exit;
        end;
        if odd(i) then begin
          a[j] := b;
          b := 0;
          dec(j);
        end;
      end;
      code := 0;
    end;
  end
  else code := -3;
end;


{---------------------------------------------------------------------------}
procedure Hex2Ext(const hex: string; var x: extended; var code: integer);
  {-Convert big-endian hex string to extended, leading $ is skipped, OK if code=0;}
  { hex must have 20 hex characters (21 if leading $), inverse of Ext2Hex.}
var
  a: THexExtA;
begin
  Hex2Float(hex, 9, a, code);
  if code=0 then x := extended(a);
end;


{---------------------------------------------------------------------------}
procedure Hex2Dbl(const hex: string; var d: double; var code: integer);
  {-Convert big-endian hex string to double, leading $ is skipped, OK if code=0;}
  { hex must have 16 hex characters (17 if leading $), inverse of Dbl2Hex.}
var
  a: THexExtA;
  t: double absolute a;
begin
  Hex2Float(hex, 7, a, code);
  if code=0 then d := t;
end;


{---------------------------------------------------------------------------}
procedure Hex2Sgl(const hex: string; var s: single; var code: integer);
  {-Convert big-endian hex string to single, leading $ is skipped, OK if code=0;}
  { hex must have 8 hex characters (9 if leading $), inverse of Sgl2Hex.}
var
  a: THexExtA;
  t: single absolute a;
begin
  Hex2Float(hex, 3, a, code);
  if code=0 then s := t;
end;


{---------------------------------------------------------------------------}
function fisEQd(x,y: double): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are bit-identical}
var
  xr: TDblRec absolute x;
  yr: TDblRec absolute y;
begin
  fisEQd := (xr.hm=yr.hm) and (xr.lm=yr.lm);
end;


{---------------------------------------------------------------------------}
function fisEQx(x,y: extended): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are bit-identical}
var
  xr: TExtRec absolute x;
  yr: TExtRec absolute y;
begin
  fisEQx := (xr.xp=yr.xp) and (xr.hm=yr.hm) and (xr.lm=yr.lm);
end;


{---------------------------------------------------------------------------}
function fisNEd(x,y: double): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are not bit-identical}
var
  xr: TDblRec absolute x;
  yr: TDblRec absolute y;
begin
  fisNEd := (xr.hm<>yr.hm) or (xr.lm<>yr.lm);
end;


{---------------------------------------------------------------------------}
function fisNEx(x,y: extended): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are not bit-identical}
var
  xr: TExtRec absolute x;
  yr: TExtRec absolute y;
begin
  fisNEx := (xr.xp<>yr.xp) or (xr.hm<>yr.hm) or (xr.lm<>yr.lm);
end;


{---------------------------------------------------------------------------}
function fisEQs(x,y: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are bit-identical}
var
  xr: longint absolute x;
  yr: longint absolute y;
begin
  fisEQs := xr=yr;
end;


{---------------------------------------------------------------------------}
function fisNEs(x,y: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are not bit-identical}
var
  xr: longint absolute x;
  yr: longint absolute y;
begin
  fisNEs := xr<>yr;
end;


{---------------------------------------------------------------------------}
function isign(x: extended): integer;
  {-Return the sign of x, 0 if x=0 or NAN}
begin
  if IsNaN(x) or (x=0.0) then isign := 0
  else if x>0.0 then isign := 1
  else isign := -1;
end;


{---------------------------------------------------------------------------}
function maxd(x, y: double): double; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the maximum of two doubles; x,y <> NAN}
begin
  if x>y then maxd := x
  else maxd:= y;
end;


{---------------------------------------------------------------------------}
function maxx(x, y: extended): extended; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the maximum of two extendeds; x,y <> NAN}
begin
  if x>y then maxx := x
  else maxx := y;
end;


{---------------------------------------------------------------------------}
function maxs(x, y: single): single; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the maximum of two singles; x,y <> NAN}
begin
  if x>y then maxs := x
  else maxs := y;
end;


{---------------------------------------------------------------------------}
function mind(x, y: double): double; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the minimum of two doubles; x,y <> NAN}
begin
  if x<y then mind := x
  else mind := y;
end;


{---------------------------------------------------------------------------}
function minx(x, y: extended): extended; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the minimum of two extendeds; x,y <> NAN}
begin
  if x<y then minx := x
  else minx := y;
end;


{---------------------------------------------------------------------------}
function mins(x, y: single): single; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the minimum of two singles; x,y <> NAN}
begin
  if x<y then mins := x
  else mins := y;
end;


{---------------------------------------------------------------------------}
function RandG01: extended;
  {-Random number from standard normal distribution (Mean=0, StdDev=1)}
var
  s,v1,v2: extended;
{$ifdef J_OPT}
  {$J+}
{$endif}
const
  RGNxt01: extended = 0.0;
  RGValid: boolean  = false;
begin
  {Generate two random numbers from the standard normal distribution with}
  {the polar rejection method, see Knuth [32], Alg. P, Section 3.4.1. One}
  {number is returned immediately, the second is saved for the next call.}

  {Note: there may be issues with multi-threading and static typed const.}
  {But here we are using random which has the same problems via randseed }
  {and/or the Mersenne twister static variables in newer FPC versions.}

  if RGValid then RandG01 := RGNxt01
  else begin
    repeat
      v1 := 2.0*random - 1.0;
      v2 := 2.0*random - 1.0;
      s  := sqr(v1) + sqr(v2);
    until (s<1.0) and (s>0.0);
    s := sqrt(-2.0*ln(s)/s);
    RandG01 := v1*s;
    RGNxt01 := v2*s;
  end;
  RGValid := not RGValid;
end;


{---------------------------------------------------------------------------}
function RandG(Mean, StdDev: extended): extended;
  {-Random number from Gaussian (normal) distribution with given mean}
  { and standard deviation |StdDev|}
begin
  RandG := Mean + StdDev*RandG01;
end;


{---------------------------------------------------------------------------}
function SafeDivX(x,y: extended): extended;
  {-Safe quotient x/y, +-Infinity if overflow, NaN if x=y=0}
var
  a: extended;
begin
  a := abs(y);
  if a >= 1.0 then SafeDivX := x/y
  else begin
    if (a=0.0) and (x=0.0) then SafeDivX := Nan_x
    else begin
      if abs(x) > a*MaxExtended then begin
        if (TExtRec(x).xp xor TExtRec(y).xp) and $8000 = 0 then begin
          SafeDivX := PosInf_x;
        end
        else begin
          SafeDivX := NegInf_x;
        end;
      end
      else SafeDivX := Nan_x;
    end;
  end;
end;


{--------------------------------------------------------------------}
{----------------- ext2 (double-extended) functions  ----------------}
{--------------------------------------------------------------------}

{$ifdef FPC2Plus}
  {$ifndef DEBUG}
    {$warn 5036 OFF}  {Local variable "xx" does not seem to be initialized}
  {$endif}
{$endif}

const
  CSX = 4294967297.0; {2^32+1}  {Split constant for extended}


{---------------------------------------------------------------------------}
procedure xxto2x(a,b: extended; var x: ext2); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a + b using TwoSum algorithm}
var
  z: extended;
begin
  {This the TwoSum algorithm, see Knuth [32], sect. 4.2.2, Theorem A/B}
  {or Ogita et al. [8], Algorithm 3.1, or Graillat [46], Algorithm 2.1}
  {FastTwosum from older AMATH versions is correct only if |a| >= |b| }
  x.h := a + b;
  z   := x.h - a;
  x.l := (a - (x.h - z)) + (b - z);
end;


{---------------------------------------------------------------------------}
procedure hhto2x(const a,b: THexExtW; var x: ext2);
  {-Return x = a + b using xxto2x}
begin
  xxto2x(extended(a), extended(b), x);
end;


{---------------------------------------------------------------------------}
procedure xadd12(const a, b: extended; var xh, xl: extended); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return [xh,xl] = a + b}
var
  z: extended;
begin
  xh := a + b;
  z  := xh - a;
  xl := (a - (xh - z)) + (b - z);
end;


{---------------------------------------------------------------------------}
procedure xto2x(a: extended; var x: ext2);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a}
begin
  x.h := a;
  x.l := 0.0;
end;


{---------------------------------------------------------------------------}
procedure xmul12(a,b: extended; var xh,xl: extended); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a * b}
var
  t,a1,a2,b1,b2: extended;
begin
  {G.W. Veltkamp's multiplication routine, see Dekker[64] mul12}
  {See also Linnainmaa[65]: exactmul2; and [8],[46]: TwoProduct}
  t  := a*CSX; a1 := (a - t) + t; a2 := a - a1;
  t  := b*CSX; b1 := (b - t) + t; b2 := b - b1;
  xh := a*b;
  xl := (((a1*b1 - xh) + a1*b2) + a2*b1) + a2*b2;
end;


{---------------------------------------------------------------------------}
procedure sqr12x(a: extended; var xh,xl: extended); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return [xh,xl] = a^2}
var
  t,a1,a2: extended;
begin
  {mul12 with _b=_a, avoids second split}
  t  := a*CSX; a1 := (a - t) + t; a2 := a - a1;
  xh := sqr(a);
  xl := ((sqr(a1) - xh) + 2.0*a1*a2) + sqr(a2);
end;


{---------------------------------------------------------------------------}
procedure mul2x(const a,b: ext2; var x: ext2);
  {-Return x = a*b}
var
  zh,zl: extended;
begin
  {Linnainmaa[65]: longmul}
  xmul12(a.h, b.h, zh, zl);
  zl  := ((a.h + a.l)*b.l + a.l*b.h) + zl;
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure mul21x(const a: ext2; b: extended; var x: ext2);
  {-Return x = a*b}
var
  zh,zl: extended;
begin
  {mul2 with b.h=b, b.l=0}
  xmul12(a.h, b, zh, zl);
  zl  := a.l*b + zl;
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure sqr2x(const a: ext2; var x: ext2);
  {-Return x = a^2}
var
  zh,zl: extended;
begin
  {Simplified mul2}
  sqr12x(a.h, zh, zl);
  zl  := ((a.h + a.l)*a.l + a.l*a.h) + zl;
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure div2x(const a,b: ext2; var x: ext2);
  {-Return x = a/b,  b<>0}
var
  qh,ql,zh,zl: extended;
begin
  {Linnainmaa[65]: longdiv}
  zh := a.h/b.h;
  xmul12(b.h, zh, qh, ql);
  zl := ((((a.h - qh) - ql) + a.l) - zh*b.l) / (b.h+b.l);
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure div21x(const a: ext2; b: extended; var x: ext2);
  {-Return x = a/b,  b<>0}
var
  qh,ql,zh,zl: extended;
begin
  {div2 with b.l=0}
  zh := a.h/b;
  xmul12(b, zh, qh, ql);
  zl := (((a.h - qh) - ql) + a.l) / b;
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure inv2x(const b: ext2; var x: ext2);
  {-Return x = 1/b,  b<>0}
var
  qh,ql,zh,zl: extended;
begin
  {div2 with a,h=1, a.l=0}
  zh := 1.0/b.h;
  xmul12(b.h, zh, qh, ql);
  zl := (((1.0 - qh) - ql) - zh*b.l) / (b.h+b.l);
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure xdivrem(a,b: extended; var q,r: extended);
  {-Compute q,r with a = q*b + r,  assumes round to nearest}
var
  x,y: extended;
begin
  q := a/b;
  xmul12(q,b,x,y);
  r := (a-x)-y;
end;


{---------------------------------------------------------------------------}
procedure add2x(const a,b: ext2; var x: ext2); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a+b}
var
  zh,zl: extended;
begin
  {Linnainmaa[65]: longadd}
  zh := a.h + b.h;
  zl := a.h - zh;
  zl := (((zl + b.h) + (a.h-(zl + zh))) + a.l) + b.l;
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure add21x(const a: ext2; b: extended; var x: ext2); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a+b}
var
  zh,zl: extended;
begin
  {add2 with b.l=0}
  zh := a.h + b;
  zl := a.h - zh;
  zl := ((zl + b) + (a.h-(zl + zh))) + a.l;
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure sub2x(const a,b: ext2; var x: ext2); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a - b}
var
  zh,zl: extended;
begin
  {add2 with a + (-b)}
  zh := a.h - b.h;
  zl := a.h - zh;
  zl := (((zl - b.h) + (a.h-(zl + zh))) + a.l) - b.l;
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure sqrt2x(const a: ext2; var x: ext2);
  {-Return x = sqrt(a),  a >= 0}
var
  qh,ql,zh,zl: extended;
begin
  {Dekker[64] sqrt2, with mul12 replaced by sqr12}
  zh := sqrt(a.h);
  sqr12x(zh, qh, ql);
  zl := 0.5*(((a.h - qh) - ql) + a.l)/zh;
  x.h := zh + zl;
  x.l := (zh - x.h) + zl;
end;


{---------------------------------------------------------------------------}
procedure pow2xi(const a: ext2; n: extended; var y: ext2);
  {-Return y = a^n, frac(n)=0, a<>0 if n<0}
var
  x: extended;
  z: ext2;
begin
  {$ifdef debug}
    if frac(n) <> 0.0 then begin
      writeln('pow2xi:  frac(n) <>0');
    end;
  {$endif}
  x := int(n);
  if x<0.0 then begin
    {1/a^x instead of (1/a)^x maybe a bit more accurate but can overflow}
    inv2x(a,y);
    pow2xi(y,-x,y);
  end
  else if x<=2.0 then begin
    if x=1.0 then y := a
    else if x=2.0 then sqr2x(a,y)
    else begin
      {a^0 = 1}
      y.h := 1.0;
      y.l := 0.0;
    end;
  end
  else if ((a.h=0.0) and (a.l=0.0)) then y := a {=0}
  else begin
    {here a<>0, x>2}
    z.h := a.h;
    z.l := a.l;
    y.h := 1.0;
    y.l := 0.0;
    while true do begin
      x := 0.5*x;
      if frac(x)<>0.0 then mul2x(y,z,y);
      x := int(x);
      if x=0.0 then exit
      else sqr2x(z,z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure pi2x(var x: ext2);
  {-Return x = Pi with double-extended precision}
const
  pih: THexExtW = ($C235,$2168,$DAA2,$C90F,$4000);  {3.1415926535897932385}
  pil: THexExtW = ($8CBB,$FC8F,$75D1,$ECE6,$BFBE);  {-5.0165576126683320234E-20}
begin
  x.l := extended(pil);
  x.h := extended(pih);
end;


{---------------------------------------------------------------------------}
function fma_x(a,b,c: extended): extended;
  {-Accurately compute a*b + c}
var
  x,y,z: extended;
  ah,al,bh,bl: extended;
begin

  {(x,y) = TwoProdct(a,b) [46] Algorithm 5}

  {Split a, [46] Algorithm 4}
  x  := a * b;
  z  := a*CSX;
  ah := z - a;
  ah := z - ah;
  al := a - ah;

  {Split b}
  z  := b*CSX;
  bh := z - b;
  bh := z - bh;
  bl := b - bh;

  {Finish TwoProduct}
  y := ah*bh;
  y := x - y;
  z := al*bh;
  y := y - z;
  z := ah*bl;
  y := y - z;
  z := al*bl;
  y := z - y;

  {add_dd_d(x,y,c) from [46] Algorithm 11}
  {(ah, al) = TwoSum(x,y) from [43] Algorithm 2}
  ah := x  + c;
  z  := ah - x;
  al := ah - z;
  al := x  - al;
  z  := c  - z;
  al := (al + z) + y; {merged step: al+y}

  {(x,y) = FastTwoSum(ah, al+y) from [46] Algorithm 3}
  x  := ah + al;
  z  := x  - ah;
  y  := al - z;

  fma_x := x + y;

end;


{---------------------------------------------------------------------------}
procedure floor2x(const a: ext2; var x: ext2);
  {-Return x = floor(a)}
var
  t: extended;
begin
  t := floorx(a.h);
  if t<>a.h then begin
    x.h := t;
    x.l := 0.0;
  end
  else xxto2x(t,floorx(a.l),x);
end;


{---------------------------------------------------------------------------}
procedure ceil2x(const a: ext2; var x: ext2);
  {-Return x = ceil(a)}
var
  t: extended;
begin
  t := ceilx(a.h);
  if t<>a.h then begin
    x.h := t;
    x.l := 0.0;
  end
  else xxto2x(t,ceilx(a.l),x);
end;


{---------------------------------------------------------------------------}
procedure trunc2x(const a: ext2; var x: ext2);
  {-Return x = trunc(a)}
begin
  if a.h >= 0.0 then floor2x(a,x)
  else ceil2x(a,x);
end;


{---------------------------------------------------------------------------}
procedure ldexp2x(const a: ext2; n: longint; var x: ext2);
  {-Return x = a*2^n}
begin
  x.h := ldexp(a.h, n);
  x.l := ldexp(a.l, n);
end;


{---------------------------------------------------------------------------}
procedure abs2x(const a: ext2; var x: ext2);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = abs(a)}
begin
  if a.h < 0.0 then begin
    x.h := -a.h;
    x.l := -a.l;
  end
  else begin
    x.h := a.h;
    x.l := a.l;
  end;
end;


{---------------------------------------------------------------------------}
procedure chs2x(const a: ext2; var x: ext2);   {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = -a}
begin
  x.h := -a.h;
  x.l := -a.l;
end;


{---------------------------------------------------------------------------}
function cmp2x(const a,b: ext2): integer;
  {-Return sign(a-b)}
var
  t: ext2;
begin
  sub2x(a,b,t);
  cmp2x := isign(t.h + t.l);
end;


{---------------------------------------------------------------------------}
procedure ln2_2x(var x: ext2);
  {-Return x = ln(2) with double-extended precision}
const
  ln2h: THexExtW = ($79AC,$D1CF,$17F7,$B172,$3FFE);  {6.9314718055994530943E-1}
  ln2l: THexExtW = ($2543,$F034,$319F,$D871,$BFBC);  {-1.1458352726798732811E-20}
begin
  x.l := extended(ln2l);
  x.h := extended(ln2h);
end;


{---------------------------------------------------------------------------}
procedure exp1_2x(var x: ext2);
  {-Return x = exp(1) with double-extended precision}
const
  e2h: THexExtW = ($4A9B,$A2BB,$5458,$ADF8,$4000);  {2.7182818284590452354}
  e2l: THexExtW = ($861C,$B185,$53BF,$A047,$BFBF);  {-6.7880636641277841170E-20}
begin
  x.l := extended(e2l);
  x.h := extended(e2h);
end;


{---------------------------------------------------------------------------}
procedure exp2x(const a: ext2; var x: ext2);
  {-Return x = exp(a)}
var
  m: longint;
  r,s,t,p: ext2;
  f: extended;
  i: integer;
const
  IMAX = 25; {Max iteration index, IMAX! must be exact} {double = 22}
  eps  = 0.2295887403949780289e-40;  {eps_x^2/512}
  shift = 10;
  FAC   = 1 shl shift;
begin

  {based on  D.H. Bailey's QD package, function exp in dd_real.cpp}
  {available from https://www.davidhbailey.com/dhbsoftware/       }

  if a.h < ln_MinExt then begin
    x.h := 0.0;
    x.l := 0.0;
    exit;
  end
  else if a.h > ln_MaxExt then begin
    x.h := PosInf_x;
    x.l := 0.0;
    exit;
  end;

  {Range reduction}
  m := round(a.h/ln2);
  ln2_2x(r);
  mul21x(r,m,r);
  sub2x(a,r,r);
  r.h := r.h/FAC;
  r.l := r.l/FAC;

  {Taylor series}
  sqr2x(r,p);
  s.h := 0.5*p.h;
  s.l := 0.5*p.l;
  add2x(r,s,s);
  mul2x(p,r,p);
  f := 6.0;
  div21x(p,f,t);

  for i:=4 to IMAX do begin
    add2x(s,t,s);
    mul2x(p,r,p);
    f := f*i;
    div21x(p,f,t);
    if abs(t.h) < eps then break;
    {$ifdef debug}
      if i=IMAX then writeln('exp2x: i=IMAX');
    {$endif}
  end;
  add2x(s,t,s);
  {undo range reduction}
  for i:=1 to shift do begin
    sqr2x(s,t);
    s.h := 2*s.h;
    s.l := 2*s.l;
    add2x(s,t,s);
  end;
  add21x(s,1,s);
  ldexp2x(s, m, x);
end;


{---------------------------------------------------------------------------}
procedure ln2x(const a: ext2; var x: ext2);
  {-Return x = ln(a)}
var
  y,z: ext2;
begin
  {x_1 = x - f(x)/f'(x) = x - (e^x-a)/e^x = x - 1 + ae^(-x)}
  {    = ae^(-x) - 1 -(-x) = ae^y - 1 - y with y=-ln(x)    }
  {Initial approximation)}
  y.h := -ln(a.h);
  y.l := 0;
  {One Newton step}
  exp2x(y,z);
  mul2x(z,a,z);
  add21x(z,-1.0,z);
  sub2x(z,y,x);
end;


{---------------------------------------------------------------------------}
procedure pow2x(const a,b: ext2; var x: ext2);
  {-Return x = a^b, a > 0}
var
  t: ext2;
begin
  {x = exp(b*ln(a))}
  ln2x(a,t);
  mul2x(t,b,t);
  exp2x(t,x)
end;


{---------------------------------------------------------------------------}
procedure nroot2x(const a: ext2; n: longint; var x: ext2);
  {-Return x = a^(1/n), a > 0 if n is even}
var
  k: longint;
  z: extended;
begin
  k := abs(n);
  z := a.h;
  if (k and 1 = 0) and (z<0.0) then begin
    x.h := nroot(z,n);
    x.l := 0;
    exit;
  end
  else abs2x(a,x);
  if k=2 then sqrt2x(x,x)
  else begin
    ln2x(x,x);
    div21x(x,k,x);
    exp2x(x,x);
  end;
  if z < 0 then chs2x(x,x);
  if k < 0 then inv2x(x,x);
end;

end.
