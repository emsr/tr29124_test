unit DAMath;

{Double precision accurate math unit}

interface

{$i std.inc}


{$ifdef BIT16}
{$N+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Double precision accurate math unit

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REMARK        :  Assumes IEEE-754 53 bit double precision (binary64)
                  and rounding to nearest

 REFERENCES    :  References used in this unit, main index in damath_info.txt/references

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
                 [19] Boost C++ Libraries, Release 1.42.0, 2010 or later version,
                      http://www.boost.org/
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
 0.10     12.01.13  W.Ehrhardt  Initial BP7 version from AMath V1.59
 0.11     12.01.13  we          BP7 compiles without errors
 0.12     12.01.13  we          D2+/FPC compile without errors
 0.13     12.01.13  we          64-bit D17/FPC260 compile without errors
 0.14     12.01.13  we          frexpd/ilogb
 0.15     13.01.13  we          PurePascal: exp
 0.16     13.01.13  we          PurePascal: expm1
 0.17     14.01.13  we          PurePascal: exp2
 0.18     14.01.13  we          PurePascal: exp10; constant ln10
 0.19     14.01.13  we          PurePascal: exp3, exp5, exp7; constants ln3, ln5, ln7
 0.20     14.01.13  we          double: exprel,cbrt,sinh,cosh,tanh,coth,sech,csch
 0.21     15.01.13  we          renamed: floorx to floord, ceilx to ceild
 0.22     15.01.13  we          ph_cutoff, rem_pio2_xx routines, rem_int2
 0.23     15.01.13  we          PurePascal: _tan and _cot
 0.24     16.01.13  we          double: ln1pmx, logcf, vers, covers, rint
 0.25     16.01.13  we          double: arccosh, arccosh1p; removed ac_help
 0.26     17.01.13  we          double: arctan2, orient2d
 0.27     17.01.13  we          NegZero_d; some functions as inline
 0.28     17.01.13  we          double: power with table size 32
 0.29     18.01.13  we          power table size 128, poly coefficients as THexDblW
 0.30     18.01.13  we          double: sincos, langevin, coshm1, exp2m1
 0.31     19.01.13  we          double: rem_2pi, remaining constants, nroot without Newton
 0.32     20.01.13  we          removed FPU control functions, improved ln1p
 0.33     22.01.13  we          bugfix sincos/Pi, correctly rounded LnSqrt2Pihex
 0.34     23.01.13  we          some more HexDbl values corrected
 0.35     28.01.13  we          Nan/Inf handling for most elementary transcendental functions
 0.36     28.01.13  we          Remove remaining extendeds and AMath compiler fixes
 0.37     29.01.13  we          Improved coshm1
 0.38     30.01.13  we          Changed ldexp/exp? to generate overflow instead of returning INF
 0.39     08.02.13  we          Special cases |y|=1/2 in power
 0.40     09.02.13  we          cosh parameter adjusted to double precision
 0.41     22.02.13  we          Constants succd0Hex, succd0, ln_succd0
 0.42     26.02.13  we          gd/arcgd parameter adjusted to double precision
 0.43     05.04.13  we          frexpd/s return non-zero exponents for denormal
 0.44     07.04.13  we          ilogb returns floor(log2()) for denormal
 0.45     07.04.13  we          exp10m1, log2p1, log10p1
 0.46     08.04.13  we          fmod support for denormal
 0.47     25.04.13  we          arcsec for large x
 0.48     15.06.13  we          Fix for VPC round
 0.49     29.06.13  we          Power returns +-INF for overflow
 0.50     17.08.13  we          expmx2h = exp(-0.5*x^2)
 0.51     25.09.13  we          const one_d = 1.0 to avoid FPC nonsense using single for e.g. 1.0/n
 0.52     28.09.13  we          const THREE: double=3.0;  used for circumventing some 'optimizations'
 0.53     03.10.13  we          adjust arccot/arccotc for x=+-Inf
 0.54     04.10.13  we          Degree versions of trig / invtrig functions
 0.55     13.10.13  we          coshm1, arcsinh, arcgd
 0.56     13.10.13  we          ln1p
 0.57     13.10.13  we          tanh, sinh
 0.58     19.10.13  we          Sqrt2 as (hex)double constant
 0.59     24.10.13  we          Improved degrees versions of trig functions
 0.60     28.10.13  we          remainder
 0.61     21.11.13  we          sinhcosh, improved sinh/cosh
 0.62     28.03.14  we          LnPi = ln(Pi) as (hex)double constant
 0.63     21.04.14  we          Improved Langevin
 0.64     30.05.14  we          PolEvalS
 0.65     07.06.14  we          ln1mexp
 0.66     20.06.14  we          logit
 0.67     05.07.14  we          basic dbl2 routines
 0.68     06.07.14  we          pow1p
 0.69     06.07.14  we          compound
 0.70     06.07.14  we          Fix dbl2 routines for Delphi 64-bit
 0.71     06.07.14  we          Fix fisEQd/fisNEd for FPC64 with -O2/3
 0.72     06.10.14  we          logistic
 0.73     07.11.14  we          Small Bernoulli numbers as Hex constants from sdBasic
 0.74     05.01.15  we          Minor changes (editorial, some FP constants)
 0.75     10.01.15  we          ddto2d uses TwoSum instead of FastTwoSum
 0.76     10.01.15  we          pi2d
 0.77     26.02.15  we          powpi2k, powpi
 0.78     03.04.15  we          special case y=-1 in power
 0.79     17.06.15  we          versint
 0.80     20.06.15  we          sinhc
 0.81     25.06.15  we          sinhmx
 0.82     25.07.15  we          avoid ilogb in power if y to large for intpower
 0.83     25.08.15  we          Constant SIXX = 6.0 to avoid some 'optimizations'
 0.84     25.04.16  we          arccos1m = arccos(1-x)
 0.85     28.09.16  we          sqrt1pmx
 0.86     09.06.17  we          lncosh, lnsinh, tanPi
 0.87     11.06.17  we          improved lncosh
 0.88     22.06.17  we          pow1pf (without dbl2)
 0.89     29.06.17  we          PolEvalC, PolEvalDeriv
 0.90     01.07.17  we          PolEvalCHEDer
 0.91     01.07.17  we          ln1pexp
 0.92     28.07.17  we          logaddexp, logsubexp
 0.93     28.07.17  we          dadd12, ddivrem, interfaced dmul12
 0.94     29.07.17  we          isRMNearest
 0.95     31.10.17  we          fma_d
 0.96     02.12.17  we          Suppress warnings: Local variable does not seem to be initialized
 0.97     09.01.18  we          TFuncDP: Functions with pointer to parameter(s)
 0.98     21.02.18  we          CSEvalDer
 0.99     20.02.18  we          Langevin function moved to sdMisc
 1.00     19.04.18  we          Fixes for FPC311

 1.01     17.05.18  we          comprel
 1.02     26.06.18  we          Anti-FPC311 quotient in versint
 1.03     06.07.18  we          SafeDiv, CoeffVar
 1.04     12.07.18  we          interfaced sqr12d
 1.05     23.07.18  we          improved intpower
 1.06     02.08.18  we          improved ldexpd
 1.07     21.08.18  we          signaling NaNs
 1.08     23.08.18  we          fixes for buggy trunc/round/frac for some 32-it compilers
 1.09     23.08.18  we          Fix ln for some compilers
 1.10     26.08.18  we          Some more Inf/Nan tests
 1.11     27.08.18  we          tvar (total variance of a double vector)
 1.12     24.09.18  we          Moved dbl2 routines to new unit DAMath2
(1.13     24.09.18  we          ceil2d, floor2d, trunc2d, ldexp2d, abs2d, chs2d,
                                chs2d, exp1_2d, exp2d, ln2d, ln2_2d, nroot2d, pow2d)
 1.14     14.10.18  we          Reintegrated damath2/dbl2

 1.15     05.11.18  we          Bern2n: array[0..MaxB2nSmall] absolute B2nHex
 1.16     06.11.18  we          pcov (population covariance)
 1.17     09.11.18  we          Generic aliases Infinity, NegInfinity, NaN


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

(*-------------------------------------------------------------------------
  This Pascal code uses material and ideas from open source and public
  domain libraries, see the file '3rdparty.ama' for the licenses.
---------------------------------------------------------------------------*)

const
  DAMath_Version = '1.17';

{#Z+}
{---------------------------------------------------------------------------}
{---------------------- Types, constants, variables ------------------------}
{---------------------------------------------------------------------------}
{#Z-}


type
  TFuncD   = function(x: double): double;
               {Functions with one real argument}
type
  TFuncDP  = function(x: double; p: pointer): double;
               {Functions with one real argument and pointer to parameter(s)}

type
  THexDblA = packed array[0..7] of byte;  {Double   as array of bytes}
  THexSglA = packed array[0..3] of byte;  {Single   as array of bytes}
  THexDblW = packed array[0..3] of word;  {Double   as array of word}
  THexSglW = packed array[0..1] of word;  {Single   as array of word}


type
  TDblRec  = packed record     {Double as sign, exponent, significand}
               lm: longint;    {low  32 bit of significand}
               hm: longint;    {high bits of significand, biased exponent and sign}
             end;
type
  dbl2 = record                {Double-double as unevaluated sum of}
           l: double;          {low  part and}
           h: double;          {high part}
         end;


{#Z+}
{Machine epsilons: smallest (negative) powers of 2 with 1 + eps_? <> 1}
{These constants should evaluate to the hex values if viewed with ?,mh}
const eps_d: double   = 2.2204460492503131E-16;   {Hex: 000000000000B03C}
const eps_s: single   = 1.1920929E-7;             {Hex: 00000034}

const PosInfSHex  : THexSglA = ($00,$00,$80,$7f); {single +INF  as hex}
const NegInfSHex  : THexSglA = ($00,$00,$80,$ff); {single -INF  as hex}
const NaNSHex     : THexSglA = ($ff,$ff,$ff,$7f); {a single NaN as hex}
const sNaNSHex    : THexSglA = ($01,$00,$80,$7f);  {a single signaling NaN as hex}
const MaxSglHex   : THexSglA = ($ff,$ff,$7f,$7f); {3.4028234E+38}
const MinSglHex   : THexSglA = ($00,$00,$80,$00); {1.1754944E-38}

const PosInfDHex  : THexDblW = ($0000,$0000,$0000,$7ff0); {double +INF  as hex}
const NegInfDHex  : THexDblW = ($0000,$0000,$0000,$fff0); {double -INF  as hex}
const NaNDHex     : THexDblW = ($ffff,$ffff,$ffff,$7fff); {a double NaN as hex}
const sNaNDHex    : THexDblW = ($0001,$0000,$0000,$7ff0); {a double signaling NaN as hex}
const MaxDblHex   : THexDblW = ($ffff,$ffff,$ffff,$7fef); {1.797693134862315E+308}
const MinDblHex   : THexDblW = ($0000,$0000,$0000,$0010); {2.225073858507201E-308}
const succd0Hex   : THexDblW = ($0001,$0000,$0000,$0000); {succd(0) as Hex}
const Neg0DblHex  : THexDblW = ($0000,$0000,$0000,$8000); {-0}

const ln2hex      : THexDblW = ($39EF,$FEFA,$2E42,$3FE6); {ln(2)}
const ln3hex      : THexDblW = ($030A,$7AAD,$93EA,$3FF1); {ln(3)}
const ln5hex      : THexDblW = ($8D33,$F7ED,$C041,$3FF9); {ln(5)}
const ln7hex      : THexDblW = ($5A57,$AE32,$2272,$3FFF); {ln(7)}
const ln10hex     : THexDblW = ($5516,$BBB5,$6BB1,$4002); {ln(10)}
const Sqrt_MinDH  : THexDblW = ($0000,$0000,$0000,$2000); {sqrt(MinDouble) as Hex}
const Sqrt_MaxDH  : THexDblW = ($FFFF,$FFFF,$FFFF,$5FEF); {sqrt(MaxDouble) as Hex}
const ln_MaxDH    : THexDblW = ($39EE,$FEFA,$2E42,$4086); {predd(ln(MaxDouble))}
const ln_MinDH    : THexDblW = ($BCD1,$DD7A,$232B,$C086); {succd(ln(MinDouble))}
const Pi_2hex     : THexDblW = ($2D18,$5444,$21FB,$3FF9); {Pi/2}
const Pi_4hex     : THexDblW = ($2D18,$5444,$21FB,$3FE9); {Pi/4}
const TwoPihex    : THexDblW = ($2D18,$5444,$21FB,$4019); {2*Pi}
const Pi_180hex   : THexDblW = ($9D39,$A252,$DF46,$3F91); {Pi/180}
const Sqrt_2Pihex : THexDblW = ($2706,$1FF6,$0D93,$4004); {sqrt(2*Pi)}
const log2ehex    : THexDblW = ($82FE,$652B,$1547,$3FF7); {log2(e)}
const log10ehex   : THexDblW = ($E50E,$1526,$CB7B,$3FDB); {log10(e)}
const LnPihex     : THexDblW = ($A1BD,$48E7,$50D0,$3FF2); {ln(pi)}
const LnSqrt2Pihex: THexDblW = ($BEB5,$C864,$67F1,$3FED); {ln(sqrt(2*Pi)=}
const PiSqrHex    : THexDblW = ($45DE,$C9BE,$BD3C,$4023); {Pi^2}
const SqrtPihex   : THexDblW = ($EF6B,$91B4,$5BF8,$3FFC); {sqrt(Pi)}
const EulerGamHex : THexDblW = ($B619,$FC6F,$788C,$3FE2); {Euler's constant}
const Sqrt2hex    : THexDblW = ($3BCD,$667F,$A09E,$3FF6); {sqrt(2)}


{Small Bernoulli numbers as Hex constants. Formerly part  }
{of unit sdBasic, they are now also used by unit DAMCmplx.}

const
  MaxB2nSmall  = 60;

const
  B2nHex: array[0..MaxB2nSmall] of THexDblW = ( {Bernoulli(2n), n=0..}
            ($0000,$0000,$0000,$0000),  {+0.0000000000000000000}
            ($5555,$5555,$5555,$3FC5),  {+1.6666666666666665741E-1}
            ($1111,$1111,$1111,$BFA1),  {-3.3333333333333332871E-2}
            ($8618,$1861,$6186,$3F98),  {+2.3809523809523808202E-2}
            ($1111,$1111,$1111,$BFA1),  {-3.3333333333333332871E-2}
            ($9365,$364D,$64D9,$3FB3),  {+7.5757575757575759678E-2}
            ($0330,$3033,$3303,$BFD0),  {-2.5311355311355310249E-1}
            ($AAAB,$AAAA,$AAAA,$3FF2),  {+1.1666666666666667407}
            ($5E5E,$5E5E,$5E5E,$C01C),  {-7.0921568627450977118}
            ($E3C5,$8F13,$7C4F,$404B),  {+5.4971177944862155584E+1}
            ($E72D,$72CF,$88FE,$C080),  {-5.2912424242424242493E+2}
            ($7E25,$8946,$301F,$40B8),  {+6.1921231884057970092E+3}
            ($CC0D,$0CC0,$2344,$C0F5),  {-8.6580253113553117146E+4}
            ($AAAB,$2AAA,$C06D,$4135),  {+1.4255171666666667443E+6}
            ($C654,$7115,$089B,$C17A),  {-2.7298231067816093564E+7}
            ($4840,$A4F3,$EDB2,$41C1),  {+6.0158087390064239502E+8}
            ($BCBD,$63B8,$2805,$C20C),  {-1.5116315767092157364E+10}
            ($4AAB,$CDDD,$01C1,$4259),  {+4.2961464306116668701E+11}
            ($C0AA,$ACF1,$F0FC,$C2A8),  {-1.3711655205088332031E+13}
            ($E993,$A679,$C22B,$42FB),  {+4.8833231897359318750E+14}
            ($9D59,$BF43,$2388,$C351),  {-1.9296579341940068000E+16}
            ($1733,$BC0D,$5C96,$43A7),  {+8.4169304757368256000E+17}
            ($4B2C,$F50C,$7E6C,$C401),  {-4.0338071854059454464E+19}
            ($255B,$D99F,$AA23,$445C),  {+2.1150748638081992622E+21}
            ($5E52,$61C3,$982C,$C4B9),  {-1.2086626522296526202E+23}
            ($96F9,$344E,$D17B,$4518),  {+7.5008667460769641660E+24}
            ($D38A,$6D24,$0CC4,$C57A),  {-5.0387781014810688499E+26}
            ($D8DB,$5349,$81F9,$45DD),  {+3.6528776484818122276E+28}
            ($337D,$78F1,$FC39,$C641),  {-2.8498769302450882361E+30}
            ($B44A,$1AC2,$887B,$46A7),  {+2.3865427499683627448E+32}
            ($C80D,$557C,$7C65,$C710),  {-2.1399949257225334859E+34}
            ($822E,$D3A3,$AD59,$4778),  {+2.0500975723478097390E+36}
            ($5904,$3DEF,$B0A4,$C7E3),  {-2.0938005911346379301E+38}
            ($6EDD,$7CB0,$B74E,$4850),  {+2.2752696488463514863E+40}
            ($F671,$911C,$2472,$C8BE),  {-2.6257710286239577207E+42}
            ($0BD7,$A067,$CF8E,$492C),  {+3.2125082102718031743E+44}
            ($D7E2,$41E3,$2553,$C99D),  {-4.1598278166794711978E+46}
            ($0430,$CF20,$2849,$4A0F),  {+5.6920695482035283174E+48}
            ($4DD3,$EFDE,$9295,$CA81),  {-8.2183629419784577665E+50}
            ($4986,$7CAA,$E2C5,$4AF4),  {+1.2502904327166994004E+53}
            ($4418,$397D,$1F1E,$CB6A),  {-2.0015583233248370052E+55}
            ($89D2,$5309,$2ACA,$4BE1),  {+3.3674982915364375556E+57}
            ($DB70,$C216,$AF88,$CC57),  {-5.9470970503135450205E+59}
            ($3B62,$0E8B,$21BF,$4CD1),  {+1.1011910323627976762E+62}
            ($BA89,$8FE4,$F4B4,$CD49),  {-2.1355259545253502079E+64}
            ($DE35,$F089,$9255,$4DC4),  {+4.3328896986641193847E+66}
            ($F7CC,$BB07,$0A86,$CE41),  {-9.1885528241669331811E+68}
            ($28D1,$5E18,$7B15,$4EBD),  {+2.0346896776329073708E+71}
            ($B278,$19A8,$9A6D,$CF3A),  {-4.7003833958035730158E+73}
            ($0C60,$1B9E,$05C9,$4FB9),  {+1.1318043445484249411E+76}
            ($F668,$65D5,$82EA,$D038),  {-2.8382249570693707354E+78}
            ($E2F0,$39BB,$FC50,$50B8),  {+7.4064248979678852935E+80}
            ($7189,$B5DA,$7B8B,$D13A),  {-2.0096454802756605262E+83}
            ($D52C,$8669,$2A23,$51BD),  {+5.6657170050805942089E+85}
            ($7B6D,$7F01,$AC7C,$D240),  {-1.6584511154136215904E+88}
            ($ABE9,$095B,$C7FF,$52C3),  {+5.0368859950492378390E+90}
            ($C906,$0208,$553C,$D348),  {-1.5861468237658185630E+93}
            ($52F0,$88F7,$03F5,$53CF),  {+5.1756743617545625189E+95}
            ($03B2,$7FAF,$7828,$D454),  {-1.7488921840217115846E+98}
            ($9281,$10E0,$F658,$54DB),  {+6.1160519994952182359E+100}
            ($12B0,$0508,$C136,$D563)); {-2.2122776912707833194E+103}
{#Z-}

var
  {Absolute vars, i.e. constants with hex patterns}
  Sqrt_MinDbl : double absolute Sqrt_MinDH;  {= 1.49166814624004E-0154} {=0.5^511}
  Sqrt_MaxDbl : double absolute Sqrt_MaxDH;  {= 1.34078079299426E+0154} {=2.0^512}
  ln_MaxDbl   : double absolute ln_MaxDH;    {= 709.782712893384}
  ln_MinDbl   : double absolute ln_MinDH;    {=-708.396418532264}

  ln2         : double absolute ln2hex;      {= 0.69314718055994530942}
  ln3         : double absolute ln3hex;      {= 1.09861228866810969140}
  ln5         : double absolute ln5hex;      {= 1.60943791243410037460}
  ln7         : double absolute ln7hex;      {= 1.94591014905531330511}
  ln10        : double absolute ln10hex;     {= 2.30258509299404568402}
  TwoPi       : double absolute TwoPihex;    {= 6.2831853071795864769 }
  Pi_2        : double absolute Pi_2hex;     {= 1.5707963267948966192 }
  Pi_4        : double absolute Pi_4hex;     {= 0.78539816339744830962}
  Pi_180      : double absolute Pi_180hex;   {= 0.17453292519943295769e-1}
  log2e       : double absolute log2ehex;    {= 1.4426950408889634079 }
  log10e      : double absolute log10ehex;   {= 0.43429448190325182765}
  Sqrt_TwoPi  : double absolute Sqrt_2Pihex; {= 2.5066282746310005024 }
  LnPi        : double absolute LnPihex;     {= 1.1447298858494001742}
  LnSqrt2Pi   : double absolute LnSqrt2Pihex;{= 0.91893853320467274178}

  PiSqr       : double absolute PiSqrHex;    {= 9.8696044010893586185}
  SqrtPi      : double absolute SqrtPihex;   {= 1.7724538509055160273}
  EulerGamma  : double absolute EulerGamHex; {= 0.57721566490153286061}
  Sqrt2       : double absolute Sqrt2hex;    {= 1.41421356237309504880168872421}

var
  MaxDouble   : double absolute MaxDblHex;   {1.797693134862315E+308} {=2^1024 - 2^971}
  MinDouble   : double absolute MinDblHex;   {2.225073858507201E-308} {=2^(-1022)}
  PosInf_d    : double absolute PosInfDHex;  {double +INF }
  NegInf_d    : double absolute NegInfDHex;  {double -INF }
  NaN_d       : double absolute NaNDHex;     {a double NaN}
  sNaN_d      : double absolute sNaNDHex;    {a double signaling NaN}
  succd0      : double absolute succd0Hex;   {= 4.94065645841247E-324} {= succd(0) = 2^(-1074)}
  NegZero_d   : double absolute Neg0DblHex;  {-0}
  Infinity    : double absolute PosInfDHex;  {double +INF  }
  NaN	      : double absolute NaNDHex;     {a double quiet NaN}
  NegInfinity : double absolute NegInfDHex;  {double -INF  }

var
  MaxSingle   : single absolute MaxSglHex;   {3.4028234E+38}
  MinSingle   : single absolute MinSglHex;   {1.1754944E-38}
  PosInf_s    : single absolute PosInfSHex;  {single +INF }
  NegInf_s    : single absolute NegInfSHex;  {single -INF }
  NaN_s       : single absolute NaNSHex;     {a single NaN}
  sNaN_s      : single absolute sNaNSHex;    {a single signaling NaN}

var
  Bern2n      : array[0..MaxB2nSmall] of double absolute B2nHex;

const
  sqrt_epsh   : double = 1.0536712127723507946742242E-8;  {sqrt(eps_d/2)}
  ln_succd0   = -744.440071921381;

const
  one_d       : double = 1.0;  {Avoid FPC nonsense using single for e.g. 1.0/n}

  THREE       : double = 3.0;  {Used for circumventing some 'optimizations'. D2/D3}
                               {use *0.333.. instead of /3! This gives incorrectly}
                               {rounded results for 5/3, 7/3, and many others!    }
                               {Also used by FPC 3.1.1 with -O4 optimization!     }
  SIXX        : double = 6.0;  {Same reason (SIX would be AMath sine integral)    }

const
  ph_cutoff   : double = 1073741824.0;   {2^30, threshold for Payne/Hanek}
                                         {do not change or adjust tol in rem_pio2}

{#Z+}
{---------------------------------------------------------------------------}
{------------------- Elementary transcendental functions -------------------}
{---------------------------------------------------------------------------}
{#Z-}

function arccos(x: double): double;
  {-Return the inverse circular cosine of x, |x| <= 1}

function arccos1m(x: double): double;
  {-Return arccos(1-x), 0 <= x <= 2, accurate even for x near 0}

function arccosd(x: double): double;
  {-Return the inverse circular cosine of x, |x| <= 1, result in degrees}

function arccosh(x: double): double;
  {-Return the inverse hyperbolic cosine, x >= 1. Note: for x near 1 the }
  { function arccosh1p(x-1) should be used to reduce cancellation errors!}

function arccosh1p(x: double): double;
  {-Return arccosh(1+x), x>=0, accurate even for x near 0}

function arccot(x: double): double;
  {-Return the sign symmetric inverse circular cotangent; arccot(x) = arctan(1/x), x <> 0}

function arccotc(x: double): double;
  {-Return the continuous inverse circular cotangent; arccotc(x) = Pi/2 - arctan(x)}

function arccotcd(x: double): double;
  {-Return the continuous inverse circular cotangent;}
  { arccotcd(x) = 90 - arctand(x), result in degrees }

function arccotd(x: double): double;
  {-Return the sign symmetric inverse circular cotangent,}
  { arccotd(x) = arctand(1/x), x <> 0, result in degrees }

function arccoth(x: double): double;
  {-Return the inverse hyperbolic cotangent of x, |x| > 1}

function arccsc(x: double): double;
  {-Return the inverse cosecant of x, |x| >= 1}

function arccsch(x: double): double;
  {-Return the inverse hyperbolic cosecant of x, x <> 0}

function arcgd(x: double): double;
  {-Return the inverse Gudermannian function arcgd(x), |x| < Pi/2}

function archav(x: double): double;
  {-Return the inverse haversine archav(x), 0 <= x <= 1}

function arcsec(x: double): double;
  {-Return the inverse secant of x, |x| >= 1}

function arcsech(x: double): double;
  {-Return the inverse hyperbolic secant of x, 0 < x <= 1}

function arcsin(x: double): double;
  {-Return the inverse circular sine of x, |x| <= 1}

function arcsind(x: double): double;
  {-Return the inverse circular sine of x, |x| <= 1, result in degrees}

function arcsinh(x: double): double;
  {-Return the inverse hyperbolic sine of x}

function arctan2(y, x: double): double;
  {-Return arctan(y/x); result in [-Pi..Pi] with correct quadrant}

function arctan(x: double): double;
  {-Return the inverse circular tangent of x}

function arctand(x: double): double;
  {-Return the inverse circular tangent of x, result in degrees}

function arctanh(x: double): double;
  {-Return the inverse hyperbolic tangent of x, |x| < 1}

function compound(x: double; n: longint): double;
  {-Return (1+x)^n; accurate version of Delphi/VP internal function}

function comprel(x,n: double): double;
  {-Return ((1+x)^n-1)/x; accurate version of Delphi/VP internal function}

function cos(x: double): double;
  {-Accurate version of circular cosine, uses system.cos for |x| <= Pi/4}

function cosd(x: double): double;
  {-Return cos(x), x in degrees}

function cosh(x: double): double;
  {-Return the hyperbolic cosine of x}

function coshm1(x: double): double;
  {-Return cosh(x)-1, accurate even for x near 0}

function cosPi(x: double): double;
  {-Return cos(Pi*x), result will be 1 for abs(x) >= 2^52}

function cot(x: double): double;
  {-Return the circular cotangent of x, x mod Pi <> 0}

function cotd(x: double): double;
  {-Return cot(x), x in degrees}

function coth(x: double): double;
  {-Return the hyperbolic cotangent of x, x<>0}

function covers(x: double): double;
  {-Return the coversine covers(x) = 1 - sin(x)}

function csc(x: double): double;
  {-Return the circular cosecant of x, x mod Pi <> 0}

function csch(x: double): double;
  {-Return the hyperbolic cosecant of x, x<>0}

function exp(x: double): double;
  {-Return exp(x), overflow if x>ln_MaxDbl}

function exp10(x: double): double;
  {-Return 10^x}

function exp10m1(x: double): double;
  {-Return 10^x - 1; special code for small x}

function exp2(x: double): double;
  {-Return 2^x}

function exp2m1(x: double): double;
  {-Return 2^x-1, accurate even for x near 0}

function exp3(x: double): double;
  {-Return 3^x}

function exp5(x: double): double;
  {-Return 5^x}

function exp7(x: double): double;
  {-Return 7^x}

function expm1(x: double): double;
  {-Return exp(x)-1, accurate even for x near 0}

function exprel(x: double): double;
  {-Return exprel(x) = (exp(x) - 1)/x,  1 for x=0}

function expx2(x: double): double;
  {-Return exp(x*|x|) with damped error amplification in computing exp of the product.}
  { Used for exp(x^2) = expx2(abs(x)) and exp(-x^2) = expx2(-abs(x))}

function expmx2h(x: double): double;
  {-Return exp(-0.5*x^2) with damped error amplification}

function gd(x: double): double;
  {-Return the Gudermannian function gd(x)}

function hav(x: double): double;
  {-Return the haversine hav(x) = 0.5*(1 - cos(x))}

function ln(x: double): double;
  {-Return natural logarithm of x, x may be denormal}

function lncosh(x: double): double;
  {-Return ln(cosh(x)), accurate for x ~ 0 and without overflow for large x}

function lnsinh(x: double): double;
  {-Return ln(sinh(x)), x > 0, accurate for x ~ 0 and without overflow for large x}

function ln1mexp(x: double): double;
  {-Return ln(1-exp(x)), x<0}

function ln1pexp(x: double): double;
  {-Accurately compute ln(1+exp(x)) without overflow}

function ln1p(x: double): double;
  {-Return ln(1+x), accurate even for x near 0}

function ln1pmx(x: double): double;
  {-Return ln(1+x)-x, x>-1, accurate even for -0.5 <= x <= 1.0}

function log10(x: double): double;
  {-Return base 10 logarithm of x}

function log10p1(x: double): double;
  {-Return log10(1+x), accurate even for x near 0}

function log2(x: double): double;
  {-Return base 2 logarithm of x}

function log2p1(x: double): double;
  {-Return log2(1+x), accurate even for x near 0}

function logaddexp(x,y: double): double;
  {-Accurately compute ln[exp(x) + exp(y)]}

function logsubexp(x,y: double): double;
  {-Accurately compute ln[exp(x) - exp(y)], x > y}

function logistic(x: double): double;
  {-Return logistic(x) = 1/(1+exp(-x))}

function logit(x: double): double;
  {-Return logit(x) = ln(x/(1.0-x)), accurate near x=0.5}

function logbase(b, x: double): double;
  {-Return base b logarithm of x}

function power(x, y : double): double;
  {-Return x^y; if frac(y)<>0 then x must be > 0}

function powm1(x,y: double): double;
  {-Return x^y - 1; special code for small x,y}

function pow1p(x,y: double): double;
  {-Return (1+x)^y, x > -1, with dbl2 arithmetic for critical values}

function pow1pf(x,y: double): double;
  {-Return (1+x)^y, x > -1, without dbl2, less accurate than pow1p}

function pow1pm1(x,y: double): double;
  {-Return (1+x)^y - 1; special code for small x,y}

function powpi2k(k,n: longint): double;
  {-Return accurate scaled powers of Pi, result = (Pi*2^k)^n}

function powpi(n: longint): double;
  {-Return accurate powers of Pi, result = Pi^n}

function sec(x: double): double;
  {-Return the circular secant of x, x mod Pi <> Pi/2}

function sech(x: double): double;
  {-Return the hyperbolic secant of x}

function sin(x: double): double;
  {-Accurate version of circular sine, uses system.sin for |x| <= Pi/4}

procedure sincos(x: double; var s,c: double);
  {-Return accurate values s=sin(x), c=cos(x)}

procedure sincosd(x: double; var s,c: double);
  {-Return sin(x) and cos(x), x in degrees}

procedure sincosPi(x: double; var s,c: double);
  {-Return s=sin(Pi*x), c=cos(Pi*x); (s,c)=(0,1) for abs(x) >= 2^52}

procedure sinhcosh(x: double; var s,c: double);
  {-Return s=sinh(x) and c=cosh(x)}

function sinc(x: double): double;
  {-Return the cardinal sine sinc(x) = sin(x)/x}

function sincPi(x: double): double;
  {-Return the normalised cardinal sine sincPi(x) = sin(Pi*x)/(Pi*x)}

function sind(x: double): double;
  {-Return sin(x), x in degrees}

function sinh(x: double): double;
  {-Return the hyperbolic sine of x, accurate even for x near 0}

function sinhc(x: double): double;
  {-Return sinh(x)/x, accurate even for x near 0}

function sinhmx(x: double): double;
  {-Return sinh(x)-x, accurate even for x near 0}

function sinPi(x: double): double;
  {-Return sin(Pi*x), result will be 0 for abs(x) >= 2^52}

function tan(x: double): double;
  {-Return the circular tangent of x, x mod Pi <> Pi/2}

function tand(x: double): double;
  {-Return tan(x), x in degrees}

function tanh(x: double): double;
  {-Return the hyperbolic tangent of x, accurate even for x near 0}

function tanPi(x: double): double;
  {-Return tan(Pi*x), result will be 0 for abs(x) >= 2^52}

function vers(x: double): double;
  {-Return the versine vers(x) = 1 - cos(x)}

function versint(x: double): double;
  {-Return versint(x) = integral(vers(t),t=0..x) = x - sin(x), accurate near 0}


function logN(N, x: double): double; {-Delphi alias for logbase}
function lnxp1(x: double): double;   {-Delphi alias for ln1p}


{#Z+}
{---------------------------------------------------------------------------}
{---------------------- Elementary numerical functions ---------------------}
{---------------------------------------------------------------------------}
{#Z-}
function cbrt(x: double): double;
  {-Return the cube root of x}

function ceil(x: double): longint;
  {-Return the smallest integer >= x; |x|<=MaxLongint}

function ceild(x: double): double;
  {-Return the smallest integer >= x}

function floor(x: double): longint;
  {-Return the largest integer <= x; |x|<=MaxLongint}

function floord(x: double): double;
  {-Return the largest integer <= x}

function fmod(x,y: double): double;
  {-Return x mod y, y<>0, sign(result) = sign(x)}

function hypot(x,y: double): double;
  {-Return sqrt(x*x + y*y)}

function hypot3(x,y,z: double): double;
  {-Return sqrt(x*x + y*y + z*z)}

function intpower(x: double; n: longint): double;
  {-Return x^n; via binary exponentiation (no overflow detection)}

function modf(x: double; var ip: longint): double;
  {-Return frac(x) and trunc(x) in ip, |x|<=MaxLongint}

function nroot(x: double; n: integer): double;
  {-Return the nth root of x; n<>0, x >= 0 if n is even}

function remainder(x,y: double): double;
  {-Return the IEEE754 remainder x REM y = x - rmNearest(x/y)*y}

function sqrt1pm1(x: double): double;
  {-Return sqrt(1+x)-1, accurate even for x near 0, x>=-1}

function sqrt1pmx(x: double): double;
  {-Return sqrt(1+x^2)-x}

{#Z+}
{---------------------------------------------------------------------------}
{----------------------- Floating point functions --------------------------}
{---------------------------------------------------------------------------}
{#Z-}
function  copysignd(x,y: double): double;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return abs(x)*sign(y)}

function  copysigns(x,y: single): single;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return abs(x)*sign(y)}

procedure frexpd(d: double; var m: double; var e: longint);
  {-Return the mantissa m and exponent e  of d with d = m*2^e, 0.5 <= abs(m) < 1;}
  { if d is 0, +-INF, NaN, return m=d, e=0}

procedure frexps(s: single; var m: single; var e: longint);
  {-Return the mantissa m and exponent e of s with s = m*2^e, 0.5 <= abs(m) < 1;}
  { if s is 0, +-INF, NaN, return m=s, e=0}

function  ilogb(x: double): longint;
  {-Return base 2 exponent of x. For finite x ilogb = floor(log2(|x|))}
  { otherwise -MaxLongint for x=0 and MaxLongint if x = +-INF or Nan. }

function  ldexpd(x: double; e: longint): double;
  {-Return x*2^e}

function  ldexps(s: single; e: longint): single;
  {-Return s*2^e}

function  IsInfD(d: double): boolean;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is +INF or -INF}

function  IsInfS(s: single): boolean;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is +INF or -INF}

function  IsNaND(d: double): boolean;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is a NaN}

function  IsNaNS(s: single): boolean;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is a NaN}

function  IsNaNorInfD(d: double): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is a NaN or infinite}

function  IsNaNorInfS(s: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is a NaN or infinite}

function  isRMNearest: boolean;
  {-Check if round_to_nearest without FPU control}

function  predd(d: double): double;
  {-Return next representable double after d in the direction -Inf}

function  preds(s: single): single;
  {-Return next representable single after s in the direction -Inf}

function  rint(x: double): double;
  {-Return the integral value nearest x for the current rounding mode}

function  succd(d: double): double;
  {-Return next representable double after d in the direction +Inf}

function  succs(s: single): single;
  {-Return next representable single after s in the direction +Inf}

function  scalbn(x: double; e: longint): double;
  {-Return x*2^e}

function  ulpd(d: double): double;
  {-Return the 'unit in the last place': ulpd(d)=|d|-predd(|d|) for finite d}

function  ulps(s: single): single;
  {-Return the 'unit in the last place': ulps(s)=|s|-preds(|s|) for finite s}


{#Z+}
{---------------------------------------------------------------------------}
{-------------- Polynomial, Vector, Statistic Operations -------------------}
{---------------------------------------------------------------------------}
{#Z-}

function  CoeffVar(const a: array of double; n: integer): double;
  {-Return the coefficient of variation (sdev/mean) of a double vector}

function  CSEvalD(x: double; const a: array of double; n: integer): double;
  {-Evaluate Chebyshev sum a[0]/2 + a[1]*T_1(x) +..+ a[n-1]*T_(n-1)(x) using Clenshaw algorithm}

procedure CSEvalDDer(x: double; const a: array of double; n: integer; var px,dp: double);
  {-Evaluate Chebyshev sum p(x) = a[0]/2 + a[1]*T_1(x) +..+ a[n-1]*T_(n-1)(x),}
  { using Clenshaw algorithm, return pd = p(x) and dp = p'(x) }

function  dot2(const x,y: array of double; n: integer): double;
  {-Accurate dot product sum(x[i]*y[i], i=0..n-1) of two double vectors}

function  mean(const a: array of double; n: integer): double;
  {-Compute accurate mean = sum(a[i], i=0..n-1)/n of a double vector}

procedure MeanAndStdDev(const a: array of double; n: integer; var mval, sdev: double);
  {-Accurate mean and sample standard deviation of a double vector}

procedure moment(const a: array of double; n: integer; var m1, m2, m3, m4, skew, kurt: double);
  {-Return the first 4 moments, skewness, and kurtosis of a double vector}

procedure mssqd(const a: array of double; n: integer; var mval, scale, sumsq: double);
  {-Calculate mean mval and ssqd sum((a[i]-mval)^2) of a double vector}

function  norm2(const a: array of double; n: integer): double;
  {-Calculate the 2-norm = sqrt(sum(a[i]^2, i=0..n-1)) of a double vector}

function  PolEval(x: double; const a: array of double; n: integer): double;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}

procedure PolEvalC(const a: array of double; n: integer; x,y: double; var u,v: double);
  {-Evaluate polynomial a[0] + a[1]*z + ... + a[n-1]*z^(n-1)}
  { for complex z = x + i*y, result is u + i*v}

procedure PolEvalDeriv(x: double; const a: array of double; n: integer; var px,dp: double);
  {-Evaluate polynomial p(x) a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}
  { Return px = p(x) and dp = p'(x)}

function  PolEvalS(x: double; const a: array of single; n: integer): double;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}

function  PolEvalEE(x: double; const a: array of double; n: integer; var e: double): double;
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

function  rms(const a: array of double; n: integer): double;
  {-Calculate the RMS value sqrt(sum(a[i]^2, i=0..n-1)/n) of a double vector}

procedure ssqd(const a: array of double; n: integer; var scale, sumsq: double);
  {-Calculate sum(a[i]^2, i=0..n-1) = scale^2*sumsq, scale>=0, sumsq>0}

function  sum2(const a: array of double; n: integer): double;
  {-Compute accurate sum(a[i], i=0..n-1) of a double vector}

function  sumsqr(const a: array of double; n: integer): double;
  {-Calculate sum(a[i]^2, i=0..n-1) of a double vector}

function  ssdev(const a: array of double; n: integer): double;
  {-Return the sample standard deviation of a double vector}

function  psdev(const a: array of double; n: integer): double;
  {-Return the population standard deviation of a double vector}

function  svar(const a: array of double; n: integer): double;
  {-Return the sample variance of a double vector}

function  pvar(const a: array of double; n: integer): double;
  {-Return the population variance of a double vector}

function  pcov(const x,y: array of double; n: integer): double;
  {-Calculate the population covariance = sum((x[i]-mx)(y[i]-my), i=0..n-1)/n of two double vectors}

function  tvar(const a: array of double; n: integer): double;
  {-Return the total variance of a double vector}

{#Z+}
{--------------------------------------------------------------------}
{---------- Argument reduction for trigonometric functions ----------}
{--------------------------------------------------------------------}
{#Z-}
function rem_pio2_cw(x: double; var z: double): integer;
  {-Cody/Waite reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8}

function rem_pio2_ph(x: double; var z: double): integer;
  {-Payne/Hanek reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8}

function rem_pio2(x: double; var z: double): integer;
  {-Argument reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8.}
  { Uses Payne/Hanek if |x| > ph_cutoff, Cody/Waite otherwise}

function rem_2pi(x: double): double;
  {-Return x mod 2*Pi}

function rem_2pi_sym(x: double): double;
  {-Return x mod 2*Pi, -Pi <= result <= Pi}

function rem_int2(x: double; var z: double): integer;
  {-Argument reduction of x: z*Pi = x*Pi - n*Pi/2, |z|<=1/4, result = n mod 8.}
  { Used for argument reduction in sin(Pi*x) and cos(Pi*x)}

{#Z+}
{--------------------------------------------------------------------}
{------------------------- Other function ---------------------------}
{--------------------------------------------------------------------}
{#Z-}

function DegToRad(x: double): double;
  {-Convert angle x from degrees to radians}

function RadToDeg(x: double): double;
  {-Convert angle x from radians to degrees}

function angle2(x1,x2,y1,y2: double): double;
  {-Return the accurate angle between the vectors (x1,x2) and (y1,y2)}

function area_triangle(x1,y1,x2,y2,x3,y3: double): double;
  {-Return the area of the triangle defined by the points (xi,yi)}

function Dbl2Hex(d: double): string;
  {-Return d as a big-endian hex string}

function Sgl2Hex(s: single): string;
  {-Return s as a big-endian hex string}

function fisEQd(x,y: double): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are bit-identical}

function fisNEd(x,y: double): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are not bit-identical}

function fisEQs(x,y: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are bit-identical}

function fisNEs(x,y: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are not bit-identical}

procedure Hex2Dbl(const hex: string; var d: double; var code: integer);
  {-Convert big-endian hex string to double, leading $ is skipped, OK if code=0;}
  { hex must have 16 hex characters (17 if leading $), inverse of Dbl2Hex.}

procedure Hex2Sgl(const hex: string; var s: single; var code: integer);
  {-Convert big-endian hex string to single, leading $ is skipped, OK if code=0;}
  { hex must have 8 hex characters (9 if leading $), inverse of Sgl2Hex.}

function in_triangle(x,y,x1,y1,x2,y2,x3,y3: double): boolean;
  {-Return true if the point (x,y) lies strictly inside the triangle defined}
  { by the three points (xi,yi), false if it lies on a side or outside.}

function in_triangle_ex(x,y,x1,y1,x2,y2,x3,y3: double): integer;
  {-Return +1 if the point (x,y) lies strictly inside the triangle defined by}
  { the three points (xi,yi), -1 if it lies strictly outside, 0 otherwise.}

function isign(x: double): integer;
  {-Return the sign of x, 0 if x=0 or NAN}

function maxd(x, y: double): double;  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the maximum of two doubles; x,y <> NAN}

function maxs(x, y: single): single; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the maximum of two singles; x,y <> NAN}

function mind(x, y: double): double; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the minimum of two doubles; x,y <> NAN}

function mins(x, y: single): single; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the minimum of two singles; x,y <> NAN}

function orient2d(x1,y1,x2,y2,x3,y3: double): double;
  {-Return the mathematical orientation of the three points (xi,yi): >0 if}
  { the order is counterclockwise, <0 if clockwise, =0 if they are collinear.}
  { Result is twice the signed area of the triangle defined by the points.}

function RandG(Mean, StdDev: double): double;
  {-Random number from Gaussian (normal) distribution with given mean}
  { and standard deviation |StdDev|}

function RandG01: double;
  {-Random number from standard normal distribution (Mean=0, StdDev=1)}

function SafeDiv(x,y: double): double;
  {-Safe quotient x/y, +-Infinity if overflow, NaN if x=y=0}


{#Z+}
{--------------------------------------------------------------------}
{-----------------  dbl2 (double-double) functions  -----------------}
{--------------------------------------------------------------------}
{#Z-}

procedure ddto2d(a,b: double; var x: dbl2);
  {-Return x = a + b using TwoSum algorithm}

procedure hhto2d(const a,b: THexDblW; var x: dbl2);
  {-Return x = a + b using xxto2x}

procedure dto2d(a: double; var x: dbl2);
  {-Return x = a}

procedure add2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a+b}

procedure add21d(const a: dbl2; b: double; var x: dbl2);
  {-Return x = a+b}

procedure sub2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a-b}

procedure mul2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a*b}

procedure mul21d(const a: dbl2; b: double; var x: dbl2);
  {-Return x = a*b}

procedure sqr2d(const a: dbl2; var x: dbl2);
  {-Return x = a^2}

procedure pi2d(var dpi: dbl2);
  {-Return x = Pi with double-double precision}

procedure pow2di(const a: dbl2; n: double; var y: dbl2);
  {-Return y = a^n, frac(n)=0, a<>0 if n<0}

procedure div2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a/b,  b<>0}

procedure div21d(const a: dbl2; b: double; var x: dbl2);
  {-Return x = a/b,  b<>0}

procedure inv2d(const b: dbl2; var x: dbl2);
  {-Return x = 1/b,  b<>0}

procedure sqrt2d(const a: dbl2; var x: dbl2);
  {-Return x = sqrt(a),  a >= 0}

procedure sqr12d(a: double; var xh,xl: double);
  {-Return [xh,xl] = a^2}

procedure dadd12(const a, b: double; var xh, xl: double); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return [xh,xl] = a + b}

procedure dmul12(a,b: double; var xh,xl: double); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return [xh,xl] = a * b}

procedure ddivrem(a,b: double; var q,r: double);
  {-Compute q,r with a = q*b + r,  assumes round to nearest}

function  fma_d(a,b,c: double): double;
  {-Accurately compute a*b + c}

procedure floor2d(const a: dbl2; var x: dbl2);
  {-Return x = floor(a)}

procedure ceil2d(const a: dbl2; var x: dbl2);
  {-Return x = ceil(a)}

procedure trunc2d(const a: dbl2; var x: dbl2);
  {-Return x = trunc(a)}

procedure ldexp2d(const a: dbl2; n: longint; var x: dbl2);
  {-Return x = a*2^n}

procedure abs2d(const a: dbl2; var x: dbl2);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = abs(a)}

procedure chs2d(const a: dbl2; var x: dbl2);  {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = -a}

function  cmp2d(const a,b: dbl2): integer;
  {-Return sign(a-b)}

procedure ln2_2d(var x: dbl2);
  {-Return x = ln(2) with double-double precision}

procedure exp1_2d(var x: dbl2);
  {-Return x = exp(1) with double-double precision}

procedure exp2d(const a: dbl2; var x: dbl2);
  {-Return x = exp(a)}

procedure ln2d(const a: dbl2; var x: dbl2);
  {-Return x = ln(a)}

procedure pow2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a^b, a > 0}

procedure nroot2d(const a: dbl2; n: longint; var x: dbl2);
  {-Return x = a^(1/n), a > 0 if n is even}

{#Z+}
{---------------------------------------------------------------------------}
{------------------  Internal and bugfix functions -------------------------}
{---------------------------------------------------------------------------}


{$ifndef BIT64}
  {$ifdef VirtualPascal}
    {$define buggy_round}
  {$endif}

  {$ifdef FPC}
    {$ifdef VER1}
      {$define buggy_trunc_int_frac}
      {$define buggy_round}
    {$endif}
  {$endif}

  {$ifdef Delphi}
    {$ifndef CONDITIONALEXPRESSIONS}
      {$define buggy_trunc_int_frac}
    {$endif}
  {$endif}

  {$ifdef CPUARM}
    {$define buggy_trunc_int_frac}
  {$endif}

  {$ifdef buggy_round}
  function round(x: double): longint;
  {$endif}

  {$ifdef buggy_trunc_int_frac}
  function trunc(x: double): longint;
  function int(x: double): double;
  function frac(x: double): double;
  {$endif}
{$endif}
{#Z-}

{$ifdef VER3_1_1}
  {$ifdef CPU64}
    function frac(d: double): double;
      {-FPC311 crashes for argument > 2^63}
  {$endif}
{$endif}

implementation


{---------------------------------------------------------------------------}
{------------------------ BUG FIXES ----------------------------------------}
{---------------------------------------------------------------------------}


{$ifdef buggy_trunc_int_frac}
{$ifdef CPUARM}
function trunc(x: double): longint;
begin
  if IsNanOrInfD(x) then trunc := 0
  else trunc := system.trunc(x);
end;
function int(x: double): double;
begin
  if IsNanOrInfD(x) then int := x
  else int := system.int(x);
end;
function frac(x: double): double;
begin
  if IsNanOrInfD(x) then frac := x
  else frac := system.frac(x);
end;
{$else}
{---------------------------------------------------------------------------}
function trunc(x: double): longint;
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
{---------------------------------------------------------------------------}
function int(x: double): double; assembler; {&Frame-} {&Uses none}
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
function frac(x: double): double; assembler; {&Frame-} {&Uses none}
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
{$endif}


{$ifdef buggy_round}
{$ifdef FPC}
{---------------------------------------------------------------------------}
function round(x: double): longint;
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
function round(x: double): longint; assembler; {&Frame-} {&Uses none}
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

{$ifdef VER3_1_1}
  {$ifdef CPU64}
    function frac(d: double): double;
      {-FPC311 crashes for argument > 2^63}
    begin
      frac := d - int(d);
    end;
  {$endif}
{$endif}


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}



{---------------------------------------------------------------------------}
function arctan(x: double): double;
  {-Return the inverse circular tangent of x}
var
  t: double;
const
  t0 = 2E16;
  t1 = 1E-8;
begin
  if IsNand(x) then arctan := Nan_d
  else begin
    t := abs(x);
    if t>t0 then begin
      if x<0 then arctan := -Pi_2 else arctan := Pi_2;
    end
    else if t<t1 then begin
      {arctan(x) = x*(1 - 1/3*x^2 + 1/5*x^4 + O(x^6))}
      arctan := x;
    end
    else arctan := system.arctan(x);
  end;
end;


{---------------------------------------------------------------------------}
function exp(x: double): double;
  {-Return exp(x), overflow if x>ln_MaxDbl}
const
  halF  : array[0..1] of double = (0.5,-0.5);
  invln2: THexDblW = ($82FE,$652B,$1547,$3FF7);  {1.44269504088896E+0000}
  twom1k: THexDblW = ($0000,$0000,$0000,$0170);  {2^(-1000)}
  ln2HI : array[0..1] of THexDblW = (($0000,$FEE0,$2E42,$3FE6),  { 6.93147180369124E-0001}
                                     ($0000,$FEE0,$2E42,$BFE6)); {-6.93147180369124E-0001}
  ln2LO : array[0..1] of THexDblW = (($3C76,$3579,$39EF,$3DEA),  { 1.90821492927059E-0010}
                                     ($3C76,$3579,$39EF,$BDEA)); {-1.90821492927059E-0010}
  PHex  : array[1..5] of THexDblW = (($553E,$5555,$5555,$3FC5),  { 1.66666666666666E-0001}
                                     ($BD93,$16BE,$C16C,$BF66),  {-2.77777777770156E-0003}
                                     ($DE2C,$AF25,$566A,$3F11),  { 6.61375632143793E-0005}
                                     ($6BF1,$C5D2,$BD41,$BEBB),  {-1.65339022054653E-0006}
                                     ($A4D0,$72BE,$3769,$3E66)); { 4.13813679705724E-0008}
var
  P: array[1..5] of double absolute PHex;
var
  y,hi,lo,c,t: double;
  xsb: integer;
  k,hx: longint;
begin
  {Algorithm is from FDLIBM [5], file e_exp.c: }
  {* <quote from [5]>
   * Method
   *   1. Argument reduction:
   *      Reduce x to an r so that |r| <= 0.5*ln2 ~ 0.34658.
   *      Given x, find r and integer k such that
   *
   *               x = k*ln2 + r,  |r| <= 0.5*ln2.
   *
   *      Here r will be represented as r = hi-lo for better
   *      accuracy.
   *
   *   2. Approximation of exp(r) by a special rational function on
   *      the interval [0,0.34658]:
   *      Write
   *          R(r**2) = r*(exp(r)+1)/(exp(r)-1) = 2 + r*r/6 - r**4/360 + ...
   *      We use a special Remes algorithm on [0,0.34658] to generate
   *      a polynomial of degree 5 to approximate R. The maximum error
   *      of this polynomial approximation is bounded by 2**-59. In
   *      other words,
   *          R(z) ~ 2.0 + P1*z + P2*z**2 + P3*z**3 + P4*z**4 + P5*z**5
   *      (where z=r*r, and the values of P1 to P5 are listed below)
   *      and
   *          |                  5          |     -59
   *          | 2.0+P1*z+...+P5*z   -  R(z) | <= 2
   *          |                             |
   *      The computation of exp(r) thus becomes
   *                             2*r
   *              exp(r) = 1 + -------
   *                            R - r
   *                                 r*R1(r)
   *                     = 1 + r + ----------- (for better accuracy)
   *                                2 - R1(r)
   *      where
   *                               2       4             10
   *              R1(r) = r - (P1*r  + P2*r  + ... + P5*r   ).
   *
   *   3. Scale back to obtain exp(x):
   *      From step 1, we have
   *         exp(x) = 2^k * exp(r)
   *
   * Special cases:
   *      exp(INF) is INF, exp(NaN) is NaN;
   *      exp(-INF) is 0, and
   *      for finite argument, only exp(0)=1 is exact.
   *
   * Accuracy:
   *      according to an error analysis, the error is always less than
   *      1 ulp (unit in the last place).
   *
   * Misc. info.
   *      For IEEE double
   *          if x >  7.09782712893383973096e+02 then exp(x) overflow
   *          if x < -7.45133219101941108420e+02 then exp(x) underflow
   *
   * Constants:
   * The hexadecimal values are the intended ones for the following
   * constants. The decimal values may be used, provided that the
   * compiler will convert from decimal to binary accurately enough
   * to produce the hexadecimal values shown.
   *}

  xsb := (THexDblA(x)[7] shr 7) and 1;   {sign bit of x   }
  hx  := TDblRec(x).hm and $7fffffff;    {high word of |x|}
  {filter out non-finite argument}
  if hx >= $40862E42 then begin
    {|x|>=709.78, or NaN}
    if hx >= $7ff00000 then begin
      if THexDblW(x)[3] and $7FF0=$7FF0 then begin
        {Inf or NaN}
        if IsNaND(x) or (x>0.0) then exp := x
        else exp := 0.0;
        exit;
      end;
    end;
    if x>ln_MaxDbl then begin
      exp := x*MaxDouble;
      exit;
    end
    else if x < -745.133219101941108420 then begin
      exp := 0.0;
      exit;
    end;
  end;
  {argument reduction}
  if hx > $3fd62e42 then begin      {if  |x| > 0.5 ln2}
    if hx < $3FF0A2B2 then begin    {and |x| < 1.5 ln2}
      hi := x - double(ln2HI[xsb]);
      lo := double(ln2LO[xsb]);
      k  := 1 - xsb -xsb;
    end
    else begin
      k  := trunc(double(invln2)*x+halF[xsb]);
      t  := k;
      hi := x - t*double(ln2HI[0]);    {t*ln2HI is exact here}
      lo := t*double(ln2LO[0]);
    end;
    x  := hi - lo;
  end
  else if hx < $3e300000 then begin
    {|x|<2**-28}
    exp := 1.0 + x;
    exit;
  end
  else begin
    k  := 0;
    {keep some compilers quiet: if k=0, then hi & lo are not used}
    hi := 0.0;
    lo := 0.0;
  end;
  {x is now in primary range}
  t := x*x;
  c := x - t*(P[1] + t*(P[2] + t*(P[3] + t*(P[4] + t*P[5]))));
  if k=0 then begin
    exp := 1.0 - ((x*c)/(c-2.0)-x);
    exit;
  end
  else y := 1.0 - ((lo-(x*c)/(2.0-c))-hi);
  {scale back with k'th power of 2}
  if (k >= -1021) then begin
    inc(TDblRec(y).hm, k shl 20); {add k to y's exponent}
    exp := y;
  end
  else begin
    {special treatment for denormals}
    inc(TDblRec(y).hm, (k+1000) shl 20);
    exp := y*double(twom1k);
  end;
end;


{---------------------------------------------------------------------------}
function exp2(x: double): double;
  {-Return 2^x}
var
  i: integer;
  a: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then exp2 := x
    else exp2 := 0.0;
    exit;
  end;
  a := abs(x);
  if a < 0.5 then begin
    if a < 5e-17 then exp2 := 1.0
    else exp2 := exp(ln2*x);
  end
  else if a=1.0 then begin
    if x>0.0 then exp2 := 2.0
    else exp2 := 0.5;
  end
  else if x >= 1024.0 then exp2 := x*MaxDouble
  else if x <= -1075.0 then exp2 := 0.0
  else begin
    i := trunc(x);
    if x <> i then begin
      if x<0.0 then i := trunc(x - 0.5)
      else i := trunc(x+0.5);
      a := exp(ln2*(x-i));
      exp2 := ldexpd(a,i);
    end
    else exp2 := ldexpd(1,i);
  end;
end;


{---------------------------------------------------------------------------}
function exp3(x: double): double;
  {-Return 3^x}
var
  i: integer;
  y: double;
const
  ld3:    THexDblW = ($BD68,$A39F,$5C01,$3FF9);  { 1.58496250072116E+0000} {log2(3)}
  l32_hi: THexDblW = ($0000,$9834,$3093,$3FE4);  { 6.30929753562668E-0001} {log3(2) high part}
  l32_lo: THexDblW = ($7D5E,$9C08,$53D2,$3DA3);  { 8.78909868152022E-0012} {log3(2) low  part}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then exp3 := x
    else exp3 := 0.0;
    exit;
  end;
  if abs(x) >= 679.0 then begin
    if x < 0.0 then exp3 := 0.0
    else exp3 := x*MaxDouble;
  end
  else begin
    i := trunc(x);
    if (i=0) and (abs(x)<1e-17) then exp3 := 1.0
    else begin
      y := double(ld3)*x;
      if y<0.0 then i := trunc(y - 0.5)
      else i := trunc(y+0.5);
      y := ln3*((x - i*double(l32_hi)) - i*double(l32_lo));
      exp3 := ldexpd(exp(y),i);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function exp5(x: double): double;
  {-Return 5^x}
var
  i: integer;
  y: double;
const
  ld5:    THexDblW = ($A371,$0979,$934F,$4002);  { 2.32192809488736E+0000} {log2(5)}
  l52_hi: THexDblW = ($0000,$6904,$9034,$3FDB);  { 4.30676558069536E-0001} {log5(2) high part}
  l52_lo: THexDblW = ($5488,$50CF,$F72E,$3D90);  { 3.85751956649558E-0012} {log5(2) low  part}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then exp5 := x
    else exp5 := 0.0;
    exit;
  end;
  if abs(x) >= 463.0 then begin
    if x < 0.0 then exp5 := 0.0
    else exp5 := x*MaxDouble;
  end
  else begin
    i := trunc(x);
    if (i=0) and (abs(x)<1e-17) then exp5 := 1.0
    else begin
      y := double(ld5)*x;
      if y<0.0 then i := trunc(y - 0.5)
      else i := trunc(y+0.5);
      y := ln5*((x - i*double(l52_hi)) - i*double(l52_lo));
      exp5 := ldexpd(exp(y),i);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function exp7(x: double): double;
  {-Return 7^x}
var
  i: integer;
  y: double;
const
  ld7:    THexDblW = ($042D,$7F54,$7576,$4006);  { 2.80735492205760E+0000} {log2(7)}
  l72_hi: THexDblW = ($0000,$3AC0,$CC19,$3FD6);  { 3.56207187054679E-0001} {log7(2) high part}
  l72_lo: THexDblW = ($EC16,$C44D,$5369,$3DCD);  { 5.33433787923143E-0011} {log7(2) low  part}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then exp7 := x
    else exp7 := 0.0;
    exit;
  end;
  if abs(x) >= 383.0 then begin
    if x < 0.0 then exp7 := 0.0
    else exp7 := x*MaxDouble;
  end
  else begin
    i := trunc(x);
    if (i=0) and (abs(x)<1e-17) then exp7 := 1.0
    else begin
      y := double(ld7)*x;
      if y<0.0 then i := trunc(y - 0.5)
      else i := trunc(y+0.5);
      y := ln7*((x - i*double(l72_hi)) - i*double(l72_lo));
      exp7 := ldexpd(exp(y),i);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function exp10(x: double): double;
  {-Return 10^x}
var
  i: integer;
  y: double;
const
  ld10:   THexDblW = ($A371,$0979,$934F,$400A);  { 3.32192809488736E+0000} {log2(10)}
  lg2_hi: THexDblW = ($0000,$509E,$4413,$3FD3);  { 3.01029995658610E-0001} {log10(2) high part}
  lg2_lo: THexDblW = ($12B3,$311F,$9FEF,$3D97);  { 5.37164476746700E-0012} {log10(2) low  part}
const
  p10tab: array[-22..22] of THexDblW = (
            ($5EE6,$1017,$3920,$3B5E),  {1.0E-22}
            ($6495,$E179,$FD7F,$3DA5),  {1.0E-21}
            ($4223,$0C92,$9CA1,$3BC7),  {1.0E-20}
            ($D2AC,$4FB6,$83C9,$3BFD),  {1.0E-19}
            ($43AC,$D1D2,$725D,$3C32),  {1.0E-18}
            ($D496,$4646,$0EF5,$3C67),  {1.0E-17}
            ($89BC,$97D8,$D2B2,$3C9C),  {1.0E-16}
            ($5616,$9EE7,$03AF,$3CD2),  {1.0E-15}
            ($2B9B,$86A1,$849B,$3D06),  {1.0E-14}
            ($7682,$6849,$25C2,$3D3C),  {1.0E-13}
            ($EA11,$812D,$9799,$3D71),  {1.0E-12}
            ($6495,$E179,$FD7F,$3DA5),  {1.0E-11}
            ($BDBA,$D9D7,$7CDF,$3DDB),  {1.0E-10}
            ($D694,$E826,$2E0B,$3E11),  {1.0E-09}
            ($8C3A,$E230,$798E,$3E45),  {1.0E-08}
            ($AF48,$9ABC,$D7F2,$3E7A),  {1.0E-07}
            ($ED8D,$A0B5,$C6F7,$3EB0),  {1.0E-06}
            ($68F0,$88E3,$F8B5,$3EE4),  {1.0E-05}
            ($432C,$EB1C,$36E2,$3F1A),  {1.0E-04}
            ($A9FC,$D2F1,$624D,$3F50),  {1.0E-03}
            ($147A,$47AE,$7AE1,$3F84),  {1.0E-02}
            ($999A,$9999,$9999,$3FB9),  {1.0E-01}
            ($0000,$0000,$0000,$3FF0),  {1.0E+00}
            ($0000,$0000,$0000,$4024),  {1.0E+01}
            ($0000,$0000,$0000,$4059),  {1.0E+02}
            ($0000,$0000,$4000,$408F),  {1.0E+03}
            ($0000,$0000,$8800,$40C3),  {1.0E+04}
            ($0000,$0000,$6A00,$40F8),  {1.0E+05}
            ($0000,$0000,$8480,$412E),  {1.0E+06}
            ($0000,$0000,$12D0,$4163),  {1.0E+07}
            ($0000,$0000,$D784,$4197),  {1.0E+08}
            ($0000,$0000,$CD65,$41CD),  {1.0E+09}
            ($0000,$2000,$A05F,$4202),  {1.0E+10}
            ($0000,$E800,$4876,$4237),  {1.0E+11}
            ($0000,$A200,$1A94,$426D),  {1.0E+12}
            ($0000,$E540,$309C,$42A2),  {1.0E+13}
            ($0000,$1E90,$BCC4,$42D6),  {1.0E+14}
            ($0000,$2634,$6BF5,$430C),  {1.0E+15}
            ($8000,$37E0,$C379,$4341),  {1.0E+16}
            ($A000,$85D8,$3457,$4376),  {1.0E+17}
            ($C800,$674E,$C16D,$43AB),  {1.0E+18}
            ($3D00,$6091,$58E4,$43E1),  {1.0E+19}
            ($8C40,$78B5,$AF1D,$4415),  {1.0E+20}
            ($EF50,$D6E2,$1AE4,$444B),  {1.0E+21}
            ($D592,$064D,$F0CF,$4480)); {1.0E+22}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then exp10 := x
    else exp10 := 0.0;
    exit;
  end;
  if abs(x) >= 324.0 then begin
    if x < 0.0 then exp10 := 0.0
    else exp10 := x*MaxDouble;
  end
  else begin
    i := trunc(x);
    if (i=x) and (abs(i)<=22) then exp10 := double(p10tab[i])
    else begin
      if (i=0) and (abs(x)<5e-18) then exp10 := 1.0
      else begin
        y := double(ld10)*x;
        if y<0.0 then i := trunc(y - 0.5)
        else i := trunc(y+0.5);
        y := ln10*((x - i*double(lg2_hi)) - i*double(lg2_lo));
        exp10 := ldexpd(exp(y),i);
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function expm1(x: double): double;
  {-Return exp(x)-1, accurate even for x near 0}
const
  PHex: array[0..5] of THexDblW = (
          ($0000,$0000,$CD80,$BF9C),  {-2.8127670288085937500E-2}
          ($DB56,$8270,$68B5,$3FE0),  {+5.1278186299064532072E-1}
          ($8D18,$34E5,$2757,$BFB0),  {-6.3100290693501981387E-2}
          ($92D1,$630D,$D5E7,$3F87),  {+1.1638457975729295593E-2}
          ($EC4D,$C9DC,$161A,$BF41),  {-5.2143390687520998431E-4}
          ($3689,$BF3B,$890D,$3EF6)); {+2.1491399776965686808E-5}

  QHex: array[0..5] of THexDblW = (
          ($0000,$0000,$0000,$3FF0),  {+1.0000000000000000000}
          ($D381,$9B03,$1544,$BFDD),  {-4.5442309511354755935E-1}
          ($3266,$9C09,$41F8,$3FB7),  {+9.0850389570911710413E-2}
          ($87CA,$C6B9,$A985,$BF84),  {-1.0088963629815501238E-2}
          ($C8CD,$DF8F,$A51B,$3F44),  {+6.3003407478692265934E-4}
          ($DB8C,$7BF1,$D98C,$BEF2)); {-1.7976570003654402936E-5}
const
  Y: single = 0.10281276702880859e1;
var
  P: array[0..5] of double absolute PHex;
  Q: array[0..5] of double absolute QHex;
var
  a: double;
begin
  {Ref: Boost[19], file expm1.hpp}
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then expm1 := x
    else expm1 := 0.0;
    exit;
  end;
  a := abs(x);
  if a < eps_d then expm1 := x
  else if a > 0.5 then expm1 := exp(x) - 1.0
  else expm1 := x*Y + x*(PolEval(x,P,6)/PolEval(x,Q,6));
end;


{---------------------------------------------------------------------------}
function ln1p(x: double): double;
  {-Return ln(1+x), accurate even for x near 0}
const
  PH: array[0..8] of THexDblW = (
        ($0000,$0000,$0000,$3FE0),  {+0.5                                }
        ($E38E,$8E38,$38E3,$3FFE),  {+1.88888888888888888888888888889    }
        ($3350,$FA6C,$16DD,$4007),  {+2.88616557734204793028322440087    }
        ($EF9A,$9A44,$44EF,$4002),  {+2.28366013071895424836601307190    }
        ($7520,$1FCA,$CA75,$3FEF),  {+0.993464052287581699346405228758   }
        ($4448,$F303,$ADEE,$3FCD),  {+0.231870525988173046996576408341   }
        ($E49B,$9A85,$85E4,$3F9A),  {+0.259013861955038425626660920779e-1}
        ($8010,$1079,$7980,$3F50),  {+0.100553041729512317747611865259e-2}
        ($815C,$FC5A,$0A57,$3E91)); {+0.253921822549273529665686528432e-6}
const
  QH: array[0..8] of THexDblW = (
        ($0000,$0000,$0000,$3FF0),  {+1.0                                }
        ($1C72,$71C7,$C71C,$4011),  {+4.44444444444444444444444444444    }
        ($7878,$7878,$7878,$4020),  {+8.23529411764705882352941176471    }
        ($7878,$7878,$7878,$4020),  {+8.23529411764705882352941176471    }
        ($3737,$3737,$3737,$4013),  {+4.80392156862745098039215686275    }
        ($5A5A,$5A5A,$5A5A,$3FFA),  {+1.64705882352941176470588235294    }
        ($8094,$9445,$4580,$3FD4),  {+0.316742081447963800904977375566   }
        ($D01F,$1EE3,$E3D0,$3F9E),  {+0.301659125188536953242835595777e-1}
        ($FD28,$3F64,$D95A,$3F50)); {+0.102838338132455779514603044015e-2}
var
  P: array[0..8] of double absolute PH;
  Q: array[0..8] of double absolute QH;
var
  y: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then ln1p := x
    else ln1p := Nan_d;
    exit;
  end;
  if abs(x) >= 1.25 then ln1p := ln(1.0+x)
  else if (x>=0.0625) and (x<=0.625) then begin
    {Pade approximation of degree (8,8) for y(x) = -(ln(1+x)/x-1)/x}
    {Errors: 0.22493623e-16 for x=-0.35, 0.74058664e-16 for x=0.625}
    y := PolEval(x,P,9)/PolEval(x,Q,9);
    ln1p := x - y*x*x;
  end
  else begin
    {Goldberg [3], Theorem 4}
    y := 1.0 + x;
    if y=1.0 then ln1p := x
    else ln1p := x*ln(y)/(y-1.0);
  end;
end;

(*

{---------------------------------------------------------------------------}
function ln1p(x: double): double;
  {-Return ln(1+x), accurate even for x near 0}
var
  y: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then ln1p := x
    else ln1p := Nan_d;
    exit;
  end;
  if abs(x) >= 0.6 then ln1p := ln(1.0 + x)
  else begin
    {Goldberg [3], Theorem 4}
    y := 1.0 + x;
    if y=1.0 then ln1p := x
    else ln1p := x*ln(y)/(y-1.0);
  end;
end;
*)


{$ifdef BIT16}
  {$define buggy_ln_inf}
{$else}
  {$undef buggy_ln_inf}
  {$ifdef FPC}
    {$ifdef VER1}
      {$define buggy_ln_inf}
    {$else}
      {$ifdef CPUARM}
        {$ifdef VER2}
          {$define buggy_ln_inf}
        {$endif}
      {$endif}
    {$endif}
  {$endif}
{$endif}


{---------------------------------------------------------------------------}
function ln(x: double): double;
  {-Return natural logarithm of x, x may be denormal}
begin
  {$ifdef buggy_ln_inf}
    if (THexDblW(x)[3] and $7FF0=$7FF0) or (x<=0.0) then begin
      {x is Nan, +-INF, or <= 0}
      if IsInfd(x) and (x>0) then ln := x   {ln(Inf)=inf}
      else ln := system.ln(x) + Nan_d;      {force NaN or exception}
      exit;
    end;
  {$endif}
  ln := system.ln(x);
end;


{---------------------------------------------------------------------------}
function log10(x: double): double;
  {-Return base 10 logarithm of x}
begin
  log10 := ln(x)*log10e;
end;


{---------------------------------------------------------------------------}
function log2(x: double): double;
  {-Return base 2 logarithm of x}
begin
  log2 := ln(x)*log2e;
end;


{---------------------------------------------------------------------------}
function logbase(b, x: double): double;
  {-Return base b logarithm of x}
begin
  logbase := ln(x)/ln(b);
end;


{---------------------------------------------------------------------------}
function logN(N, x: double): double;
  {-Delphi alias for logbase}
begin
  logN := ln(x)/ln(N);
end;


{---------------------------------------------------------------------------}
function arctan2(y, x: double): double;
  {-Return arctan(y/x); result in [-Pi..Pi] with correct quadrant}
var
  z: double;
  hx,hy,ix,iy,lx,ly: longint;
  k,m: integer;
const
  pi_lo = 1.2246467991473531772E-16;
begin
  {Pascal translation of FDLIBM[5] routine __ieee754_atan2 in e_atan2.c}
  with TDblRec(x) do begin
    hx := hm;
    lx := lm;
    ix := hx and $7fffffff;
  end;
  with TDblRec(y) do begin
    hy := hm;
    ly := lm;
    iy := hy and $7fffffff;
  end;
  if (hx and $7FF00000=$7FF00000) and ((hx and $000FFFFF<>0) or (lx<>0)) then begin
    {x is NaN}
    arctan2 := x;
    exit;
  end;
  if (hy and $7FF00000=$7FF00000) and ((hy and $000FFFFF<>0) or (ly<>0)) then begin
    {y is NaN}
    arctan2 := y;
    exit;
  end;

  if x=1.0 then begin
    arctan2 := arctan(y);
    exit;
  end;

  m := ((hy shr 31) and 1) or ((hx shr 30) and 2);  {2*sign(x)+sign(y)}
  if (iy or ly)=0 then begin
    {y=0}
    case m of
        1: arctan2 := y;       { arctan2( -0,+anything) =  -0 }
        2: arctan2 := Pi;      { arctan2( +0,-anything) =  Pi }
        3: arctan2 := -Pi;     { arctan2( -0,-anything) = -Pi }
      else arctan2 := y;       { arctan2( +0,+anything) =   0 }
    end;
    exit;
  end;
  if (ix or lx)=0 then begin
    {x=0}
    if hy<0 then arctan2 := -Pi_2 else arctan2 := Pi_2;
    exit;
  end;
  if ix=$7ff00000 then begin
    {x is INF}
    arctan2 := 0.0; {keep some compilers happy}
    if iy=$7ff00000 then begin
      case m of
        0: arctan2 :=  Pi_4;      { arctan2(+INF,+INF) }
        1: arctan2 := -Pi_4;      { arctan2(-INF,+INF) }
        2: arctan2 :=  3.0*Pi_4;  { arctan2(+INF,-INF) }
        3: arctan2 := -3.0*Pi_4;  { arctan2(-INF,-INF) }
       end
    end
    else begin
      case m of
       0: arctan2 := 0.0;         { arctan2(+...,+INF) }
       1: arctan2 := NegZero_d;   { arctan2(-...,+INF) }
       2: arctan2 := Pi;          { arctan2(+...,-INF) }
       3: arctan2 := -Pi;         { arctan2(-...,-INF) }
      end;
    end;
    exit;
  end;
  if iy=$7ff00000 then begin
    {y is INF}
    if hy<0 then arctan2 := -Pi_2 else arctan2 := Pi_2;
    exit;
  end;

  {compute y/x}
  k := (iy shr 20) - (ix shr 20);
  if k > 60 then z := Pi_2 + 0.5*pi_lo        { |y/x| >  2**60 }
  else if (hx<0) and (k<-60) then z := 0.0    { |y|/x < -2**60 }
  else z := arctan(abs(y/x));                 { safe to do y/x }
  case m of
      0: arctan2 := z;              { arctan2(+,+) }
      1: arctan2 := -z;             { arctan2(-,+) }
      2: arctan2 := Pi - (z-pi_lo); { arctan2(+,-) }
    else arctan2 := (z-pi_lo) - Pi; { arctan2(-,-) }
  end;
end;


{---------------------------------------------------------------------------}
function intpower0(x: double; n: longint): double;
  {-Return x^n; via binary exponentiation (no overflow detection)}
var
  i: longint;
  r: double;
begin
  i := abs(n);
  r := 1.0;
  while i > 0 do begin
    while integer(i) and 1 = 0 do begin
      i := i shr 1;
      x := x*x
    end;
    dec(i);
    r := r*x;
  end;
  if n>=0 then intpower0 := r
  else intpower0 := 1.0/r;
end;


{---------------------------------------------------------------------------}
function intpower(x: double; n: longint): double;
  {-Return x^n; via binary exponentiation (no overflow detection)}
var
  i: longint;
  r: double;
begin
  i := abs(n);
  r := abs(x);
  if (i=0) or (r=1.0) then begin
    if odd(i) and (x < 0) then intpower := -1
    else intpower := 1.0;
    exit;
  end;
  if (n<0) and (r>1.0) then begin
    if (1.0+ilogb(r))*i >= 1023 then begin
      {avoid spurious overflow, use x^(-n) = (1/x)^n}
      {so that e.g. 2^(-1030) is computed as denormal}
      x := 1.0/x;
      n := i;
    end;
  end;
  intpower := intpower0(x,n);
end;


{---------------------------------------------------------------------------}
function ldexpd(x: double; e: longint): double;
  {-Return d*2^e}
var
  i: integer;
const
  H2_54: THexDblW = ($0000,$0000,$0000,$4350);  {2^54}
begin
  {if +-INF, NaN, 0 or if e=0 return d}
  i := (THexDblW(x)[3] and $7FF0) shr 4;
  if (i=$7FF) or (e=0) or (x=0.0) then ldexpd := x
  else if i=0 then begin
    {Denormal: result = d*2^54*2^(e-54)}
    ldexpd := ldexpd(x*double(H2_54), e-54);
  end
  else begin
    e := e+i;
    if e>$7FE then begin
      {overflow}
      ldexpd := MaxDouble*copysignd(MaxDouble,x);
    end
    else if e<1 then begin
      {underflow or denormal}
      if e<-53 then ldexpd := succd0*copysignd(succd0,x) {will set correct result for rounding mode}
      else begin
        {Denormal: result = d*2^(e+54)/2^54}
        inc(e,54);
        THexDblW(x)[3] := (THexDblW(x)[3] and $800F) or (e shl 4 and $7FF0);
        ldexpd := x/double(H2_54);
      end;
    end
    else begin
      THexDblW(x)[3] := (THexDblW(x)[3] and $800F) or (e shl 4 and $7FF0);
      ldexpd := x;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function  ldexps(s: single; e: longint): single;
  {-Return s*2^e}
var
  L: longint absolute s;
begin
  if (L and $7F800000 = $7F800000) or (s=0.0) then ldexps := s
  else begin
    (*
    {use this if infinities and no overflow should be returned}
    x := ldexpd(s,e);
    if x>MaxSingle then ldexps := PosInf_s
    else if x<-MaxSingle then ldexps := NegInf_s
    else ldexps := x;
    *)
    ldexps := ldexpd(s,e);
  end;
end;


{---------------------------------------------------------------------------}
function scalbn(x: double; e: longint): double;
  {-Return x*2^e}
begin
  scalbn := ldexpd(x,e);
end;


{---------------------------------------------------------------------------}
{internal tan/cot routines for |arg| <= Pi/4. Based on Cephes 2.8 functions}
{from double/tan.c, Copyyright 1984, 1995, 2000 by Stephen L. Moshier}
const
  PHex: array[0..2] of THexDblW = (
          ($9176,$d329,$1fea,$c171),
          ($9ddd,$a5fc,$99ec,$4131),
          ($3f38,$d24f,$92d8,$c0c9));

  QHex: array[0..4] of THexDblW = (
          ($5a31,$3cbe,$afe0,$c189),
          ($d8ef,$c2ea,$d98f,$4177),
          ($bc96,$582a,$27bc,$c134),
          ($6572,$eeb3,$b8a5,$40ca),
          ($0000,$0000,$0000,$3ff0));


{---------------------------------------------------------------------------}
function _tan(x: double): double;
  {-Return the circular tangent of x, |x| <= Pi/4}
var
  z, zz: double;
  P: array[0..2] of double absolute PHex;
  Q: array[0..4] of double absolute QHex;
begin
  z  := abs(x);
  zz := z*z;
  if zz > 1e-14 then begin
    z := z + z*(zz*PolEval(zz, P, 3)/PolEval(zz, Q, 5));
  end;
  if x >= 0.0 then _tan := z
  else _tan := -z;
end;


{---------------------------------------------------------------------------}
function _cot(x: double): double;
  {-Return the circular cotangent of x, 0 < |x| <= Pi/4}
var
  z, zz: double;
  P: array[0..2] of double absolute PHex;
  Q: array[0..4] of double absolute QHex;
begin
  z  := abs(x);
  zz := z*z;
  if zz > 1e-14 then begin
    z := z + z*(zz*PolEval(zz, P, 3)/PolEval(zz, Q, 5));
  end;
  if x >= 0.0 then _cot := 1.0/z
  else _cot := (-1.0)/z;
end;


{---------------------------------------------------------------------------}
function fmod(x,y: double): double;
  {-Return x mod y, y<>0, sign(result) = sign(x)}
var
  tx,ty,mx,my: double;
  ex,ey: longint;
begin
  {Based on  mod.go  from the math package of the Go programming }
  {language, available from https://code.google.com/p/go/source/.}
  {The Go source files are distributed under a BSD-style license.}
  if (THexDblW(x)[3] and $7FF0=$7FF0) or (THexDblW(y)[3] and $7FF0=$7FF0) then fmod := NaN_d
  else if x=0.0 then fmod := 0.0
  else if y=0.0 then fmod := x/y
  else begin
    tx := abs(x);
    ty := abs(y);
    frexpd(ty,my,ey);
    {Warning: loop count depends on ilogb(x) & ilogb(y), the worst case}
    {is about 1000 loops, e.g. 1049 for x=MaxDouble, y=18*succd(0)!    }
    while ty <= tx do begin
      frexpd(tx,mx,ex);
      if mx<my then dec(ex);
      mx := ldexpd(ty, ex-ey);
      tx := tx - mx;
    end;
    if x>=0.0 then fmod := tx
    else fmod := -tx;
  end;
end;


{---------------------------------------------------------------------------}
function remainder(x,y: double): double;
  {-Return the IEEE754 remainder x REM y = x - rmNearest(x/y)*y}
var
  yh: double;
  ey,sx: word;
begin
  {Ref: FDLIBM 5.3 [5], file e_remainder.c}

  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {x is INF or NaN}
    remainder := NaN_d;
    exit;
  end;

  sx := THexDblW(x)[3] and $8000;
  ey := THexDblW(y)[3] and $7FF0;

  if ey=$7FF0 then with TDblRec(y) do begin
    {y is INF or NaN}
    if (hm and $7FFFFFFF=$7FF00000) and (lm=0) then begin
      {y is INF}
      remainder := x;
    end
    else remainder := NaN_d;
    exit;
  end;

  y := abs(y);
  if ey <$7FE0 then begin
    {|y| < 0.5*Maxdouble}
    x := fmod(x,y+y);
  end;
  x := abs(x);
  if ey<$0020 then begin
    {|y| < 2*Mindouble}
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


{---------------------------------------------------------------------------}
function exp10m1(x: double): double;
  {-Return 10^x - 1; special code for small x}
var
  z: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then exp10m1 := x
    else exp10m1 := 0.0;
    exit;
  end;
  if abs(x) > 0.3 then exp10m1 := exp10(x) - 1.0
  else if x=0.0 then exp10m1 := 0.0
  else begin
    {z = x*ln(10), 10^x-1 = expm1(x*ln(10))}
    z := 2.0*x + 0.30258509299404568402*x;
    exp10m1 := expm1(z);
  end
end;


{---------------------------------------------------------------------------}
function exp2m1(x: double): double;
  {-Return 2^x-1, accurate even for x near 0}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then exp2m1 := x
    else exp2m1 := 0.0;
    exit;
  end;
  if abs(x)<1 then exp2m1 := expm1(x*ln2)
  else exp2m1 := exp2(x)-1.0;
end;


{---------------------------------------------------------------------------}
function exprel(x: double): double;
  {-Return exprel(x) = (exp(x) - 1)/x,  1 for x=0}
const
  xsmall = 1e-5;
var
  z: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then exprel := Nan_D
    else exprel := 0.0;
    exit;
  end;
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
    if x<-38.5 then exprel := -1.0/x
    else exprel := (exp(x) - 1.0)/x;
  end;
end;


{---------------------------------------------------------------------------}
function expx2(x: double): double;
  {-Return exp(x*|x|) with damped error amplification in computing exp of the product.}
  { Used for exp(x^2) = expx2(abs(x)) and exp(-x^2) = expx2(-abs(x))}
const
  MFAC = 32768.0;
  MINV = 3.0517578125e-5;
var
  u,u1,m,f: double;
  neg: boolean;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNaND(x) or (x>0.0) then expx2 := x
    else expx2 := 0.0;
    exit;
  end;
  if x >= 26.6417475570463282 then begin
    expx2 := {PosInf_d;}x*MaxDouble;
    exit;
  end
  else if x <= -27.2844291111502141 then begin
    expx2 := 0.0;
    exit;
  end;

  {Ref: Cephes [7], file ldouble\expx2l.c}
  neg := x<0.0;
  x := abs(x);
  if x <= 1.0 then begin
    if neg then u := -x*x else u := x*x;
    expx2 := exp(u);
  end
  else begin
    {Represent x as an exact multiple of MFAC plus a residual.}
    {MFAC is a power of 2 chosen so that exp(m * m) does not  }
    {overflow or underflow and so that |x - m| is small.      }
    m := MINV*floord(MFAC*x + 0.5);
    f := x - m;

    {x^2 = m^2 + 2mf + f^2}
    u  := m*m;
    u1 := 2.0*m*f + f*f;

    if neg then begin
      u  := -u;
      u1 := -u1;
    end;

    if u+u1 > ln_MaxDbl then expx2 := {PosInf_d}x*MaxDouble
    else begin
      {u is exact, u1 is small}
      expx2 := exp(u)*exp(u1);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function expmx2h(x: double): double;
  {-Return exp(-0.5*x^2) with damped error amplification}
const
  MFAC = 32768.0;
  MINV = 3.0517578125e-5;
var
  u,u1,m,f: double;
begin
  {Ref: Cephes [7], file ldouble\expx2l.c}
  if IsNaND(x) then begin
    expmx2h := x;
    exit;
  end;
  x := abs(x);
  if x >= 38.604 then expmx2h := 0.0
  else if x <= 2.0 then begin
    expmx2h := exp(-0.5*x*x);
  end
  else begin
    {Represent x as an exact multiple of MFAC plus a residual.}
    {MFAC is a power of 2 chosen so that exp(m * m) does not  }
    {overflow or underflow and so that |x - m| is small.      }
    m := MINV*floord(MFAC*x + 0.5);
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
function cos(x: double): double;
  {-Accurate version of circular cosine, uses system.cos for |x| <= Pi/4}
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    cos := Nan_D;
    exit;
  end;
  {reduction mod Pi/2, |t| <= Pi/4}
  case rem_pio2(x,t) and 3 of
      0:  cos :=  system.cos(t);
      1:  cos := -system.sin(t);
      2:  cos := -system.cos(t);
    else  cos :=  system.sin(t);
  end;
end;


{---------------------------------------------------------------------------}
const
  npcv=7; {chebyshev(((1-cos(x))-x^2/2)/x^4,x=-Pi/4..Pi/4,1e-20) converted to poly}
  pcvh: array[0..npcv-1] of THexDblW = (
          ($5555,$5555,$5555,$BFA5),  {-0.416666666666666666661679544980e-1 }
          ($6C16,$16C1,$C16C,$3F56),  {+0.138888888888888879063613392395e-2 }
          ($9D16,$1A01,$01A0,$BEFA),  {-0.248015873015846864625912712986e-4 }
          ($3AEC,$B771,$7E4F,$3E92),  {+0.275573192214211550883212170218e-6 }
          ($54EE,$DEC9,$EED8,$BE21),  {-0.208767557953773544469035900043e-8 }
          ($6118,$B33C,$394B,$3DA9),  {+0.114704613872478358447153370873e-10}
          ($C6E4,$1E2F,$B78E,$BD2A)); {-0.474589475226674827431568738364e-13}
var
  pcv: array[0..npcv-1] of double absolute pcvh;

{---------------------------------------------------------------------------}
function vers(x: double): double;
  {-Return the versine vers(x) = 1 - cos(x)}
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    vers := Nan_D;
    exit;
  end;
  {reduction mod Pi/2, |t| <= Pi/4}
  case rem_pio2(x,t) and 3 of
      0:  begin
            {vers(x) = 1 - cos(t)}
            x := t*t;
            t := PolEval(x,pcv,npcv);
            vers := 0.5*x + sqr(x)*t;
          end;
      1:  vers := 1.0 + system.sin(t);
      2:  vers := 1.0 + system.cos(t);
    else  vers := 1.0 - system.sin(t)
  end;
end;


{---------------------------------------------------------------------------}
function covers(x: double): double;
  {-Return the coversine covers(x) = 1 - sin(x)}
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    covers := Nan_D;
    exit;
  end;
  {reduction mod Pi/2, |t| <= Pi/4}
  case rem_pio2(x,t) and 3 of
      0:  covers := 1.0 - system.sin(t);
      1:  begin
            {covers(x) = 1 - cos(t)}
            x := t*t;
            t := PolEval(x,pcv,npcv);
            covers := 0.5*x + sqr(x)*t;
          end;
      2:  covers := 1.0 + system.sin(t);
    else  covers := 1.0 + system.cos(t);
  end;
end;


{---------------------------------------------------------------------------}
function hav(x: double): double;
  {-Return the haversine hav(x) = 0.5*(1 - cos(x))}
begin
  hav := 0.5*vers(x);
end;


{---------------------------------------------------------------------------}
function cosPi(x: double): double;
  {-Return cos(Pi*x), result will be 1 for abs(x) >= 2^52}
var
 t: double;
 i: integer;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    cosPi := Nan_D;
    exit;
  end;
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
function sin(x: double): double;
  {-Accurate version of circular sine, uses system.sin for |x| <= Pi/4}
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    sin := Nan_D;
    exit;
  end;
  {reduction mod Pi/2, |t| <= Pi/4}
  case rem_pio2(x,t) and 3 of
      0:  sin :=  system.sin(t);
      1:  sin :=  system.cos(t);
      2:  sin := -system.sin(t);
    else  sin := -system.cos(t);
  end;
end;


{---------------------------------------------------------------------------}
function sinPi(x: double): double;
  {-Return sin(Pi*x), result will be 0 for abs(x) >= 2^52}
var
 t: double;
 i: integer;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    sinPi := Nan_D;
    exit;
  end;
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
procedure sincos(x: double; var s,c: double);
  {-Return accurate values s=sin(x), c=cos(x)}
var
  t,ss,cc: double;
  n: integer;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    s := Nan_D;
    c := Nan_D;
    exit;
  end;
  {reduction mod Pi/2, |t| <= Pi/4}
  n := rem_pio2(x,t) and 3;
  ss := system.sin(t);
  cc := system.cos(t);
  case n of
      0: begin s:= ss; c:= cc; end;
      1: begin s:= cc; c:=-ss; end;
      2: begin s:=-ss; c:=-cc; end;
    else begin s:=-cc; c:= ss; end;
  end;
end;


{---------------------------------------------------------------------------}
procedure sincosPi(x: double; var s,c: double);
  {-Return s=sin(Pi*x), c=cos(Pi*x); (s,c)=(0,1) for abs(x) >= 2^52}
var
  t,ss,cc: double;
  n: integer;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    s := Nan_D;
    c := Nan_D;
    exit;
  end;
  n := rem_int2(x,t) and 3;
  t := Pi*t;
  ss := system.sin(t);
  cc := system.cos(t);
  case n of
      0: begin s:= ss; c:= cc; end;
      1: begin s:= cc; c:=-ss; end;
      2: begin s:=-ss; c:=-cc; end;
    else begin s:=-cc; c:= ss; end;
  end;
end;


{---------------------------------------------------------------------------}
function tan(x: double): double;
  {-Return the circular tangent of x, x mod Pi <> Pi/2}
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    tan := Nan_D;
    exit;
  end;
  {reduction mod Pi/2, |t| <= Pi/4}
  if odd(rem_pio2(x,t)) then tan := -_cot(t)
  else tan := _tan(t);
end;


{---------------------------------------------------------------------------}
function cot(x: double): double;
  {-Return the circular cotangent of x, x mod Pi <> 0}
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    cot := Nan_D;
    exit;
  end;
  {reduction mod Pi/2, |t| <= Pi/4}
  if odd(rem_pio2(x,t)) then cot := -_tan(t)
  else cot := _cot(t);
end;


{---------------------------------------------------------------------------}
function gd(x: double): double;
  {-Return the Gudermannian function gd(x)}
var
  z: double;
begin
  if IsNanD(x) then begin
    gd := x;
    exit;
  end;
  {gd(x) = arctan(sinh(x)) = 2*arctan(tanh(x/2))}
  z := abs(x);
  if z < 38.0 then begin
    if z < 1.25e-2 then begin
      {gd(x) = x - 1/6*x^3 + 1/24*x^5 - 61/5040*x^7 + 277/72576*x^9 + O(x^10)}
      if z < 1e-4 then gd := x*(1.0-sqr(x)/SIXX)
      else begin
        z  := x*x;
        gd := x*(1.0 - z/SIXX*(1.0 - 0.25*z*(1.0 - 61.0*z/210.0)));
      end;
    end
    else gd := arctan(sinh(x));
  end
  else gd := copysignd(Pi_2,x);
end;


{---------------------------------------------------------------------------}
function versint(x: double): double;
  {-Return versint(x) = integral(vers(t),t=0..x) = x - sin(x), accurate near 0}
var
  t: double;
const
  CSN = 8;
  CSVH: array[0..CSN-1] of THexDblW = ( {chebyshev((x-sin(x))/(x^3), x=-4/3..4/3, 0.5e-20)}
          ($CDAD,$1AE7,$6A28,$3FD4),    {+0.318979288363346959407456566900    } {*2}
          ($E2B8,$72ED,$15F4,$BF7D),    {-0.710101592891792352000298121574e-2 }
          ($C469,$BD0F,$CD25,$3F13),    {+0.755361827626963385962654365195e-4 }
          ($C9A2,$36A7,$7037,$BE9F),    {-0.468467809127500007796911607370e-6 }
          ($A46E,$BA7B,$5247,$3E20),    {+0.190006184732484199713774118361e-8 }
          ($A957,$FCC1,$E1B1,$BD97),    {-0.543005219722655848794099416504e-11}
          ($C023,$0F33,$F184,$3D09),    {+0.115211934731518499595222170191e-13}
          ($978E,$7A66,$BFE9,$BC75));   {-0.188648197267528207115909378480e-16}
         {($4FC5,$49FD,$0F2C,$3BDD))}   {+0.246141587293140761985659464726e-19}
  fs = 1.125;  {2*(3/4)^2}
  t0 = 4/3;
var
  CSV: array[0..CSN-1] of double absolute CSVH;
begin
  if IsNanD(x) then begin
    versint := x;
    exit;
  end;
  t := abs(x);
  if t < t0 then begin
    if t <= sqrt_epsh then begin
      {versint(x) = x^3/6*(1 - 1/20*x^2 + O(x^4))}
      versint := x*x*x/SIXX;
    end
    else begin
      t := CSEvalD(fs*sqr(t) - 1.0, CSV, CSN);
      versint := x*x*x*t;
    end
  end
  else versint := x - sin(x);
end;


{---------------------------------------------------------------------------}
function arcgd(x: double): double;
  {-Return the inverse Gudermannian function arcgd(x), |x| < Pi/2}
var
  z,s,c: double;
begin
  if IsNanD(x) then begin
    arcgd := x;
    exit;
  end;
  {arcgd(x) := 2*arctanh(tan(x/2)) = arcsinh(tan(x)) = ln(sec(x)+tan(x))}
  z := abs(x);
  if z < sqrt_epsh then arcgd := x
  else if z<=0.85 then arcgd := arcsinh(tan(x))
  else begin
    {arcgd := ln(sec(x)+tan(x)) = ln(1/cos(x)+sin(x)/cos(x))}
    sincos(z,s,c);
    {cos(x) is > 0 for |x| < Pi/2, if c <= 0 then x >= Pi/2}
    if c<=0.0 then z := PosInf_d
    else z := ln((1.0+s)/c);
    if x>0.0 then arcgd := z
    else arcgd := -z;
  end;
end;


{---------------------------------------------------------------------------}
function logcf(x,i,d: double): double;
  {-Calculate 1/i + x/(i+d) + x^2/(i+2*d) + x^3/(i+3d) .. via continued fractions}
var
  c1,c2,c3,c4,a1,b1,a2,b2: double;
const
  hugeh: THexDblW = ($0000,$0000,$0000,$4FF0); {2^256}
  tinyh: THexDblW = ($0000,$0000,$0000,$2FF0); {1/2^256}
var
  huge: double absolute hugeh;
  tiny: double absolute tinyh;
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

  while abs(a2*b1-a1*b2) > abs(eps_d*b1*b2)  do begin
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
function ln1pmx(x: double): double;
  {-Return ln(1+x)-x, x>-1, accurate even for -0.5 <= x <= 1.0}
const
  c9 = 2/9;  {FIX311}
  c7 = 2/7;
  c5 = 2/5;
  c3 = 2/3;
  x0 = -0.525;
  x1 = 1.5;
var
  y,z: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsInfD(x) and (x>0) then ln1pmx := -x else ln1pmx := x;
    exit;
  end;
  {Based on function log1 by Ian Smith [25]}
  if (x<x0) or (x>x1) then ln1pmx := ln(1.0+x) - x
  else begin
    z := x/(2.0 + x);
    y := z*z;
    if abs(x)>1.0e-2 then y := 2.0*y*logcf(y, 3.0, 2.0)
    else y := (((c9*y + c7)*y + c5)*y + c3)*y;
    ln1pmx := z*(y-x);
  end;
end;


{---------------------------------------------------------------------------}
function ln1mexp(x: double): double;
  {-Return ln(1-exp(x)), x<0}
const
  x0 = -36.75;  { < ln(eps/2)}
  x1 = -0.693147180559945309; {-ln(2)}
begin
  if IsNanD(x) then begin
    ln1mexp := x;
    exit;
  end;
  if x <= x0 then begin
    if x < ln_succd0 then ln1mexp := 0.0
    else ln1mexp := -exp(x);
  end
  else if x < x1 then ln1mexp := ln1p(-exp(x))
  else if x > -eps_d then ln1mexp := ln(-x)
  else ln1mexp := ln(-expm1(x));
end;


{---------------------------------------------------------------------------}
function ln1pexp(x: double): double;
  {-Accurately compute ln(1+exp(x)) without overflow}
const
  x0 = -36.8;  { < ln(0.5_eps), double = -36.8}
  x1 = 18.1;   { x + exp(-x) - ln(1+exp(x) < eps_x/2, double = 18.1}
  x2 = 33.3;   { exp(-x)/x < eps_x/2, double = 33.3}
begin
  if IsNanD(x) then begin
    ln1pexp := x;
    exit;
  end;
  if x <= x0 then begin
    if x < ln_succd0 then ln1pexp := 0.0
    else ln1pexp := exp(x);
  end
  else if x <= x1 then ln1pexp := ln1p(exp(x))
  else if x <= x2 then ln1pexp := x + exp(-x)
  else ln1pexp := x;
end;


{---------------------------------------------------------------------------}
function lnxp1(x: double): double;
  {-Delphi alias for ln1p}
begin
  lnxp1 := ln1p(x);
end;


{---------------------------------------------------------------------------}
function log2p1(x: double): double;
  {-Return log2(1+x), accurate even for x near 0}
begin
  if IsNanD(x) then begin
    log2p1 := x;
    exit;
  end;
  if abs(x) > 0.5 then  log2p1 := log2(1.0 + x)
  else log2p1 := log2e*ln1p(x);
end;


{---------------------------------------------------------------------------}
function log10p1(x: double): double;
  {-Return log10(1+x), accurate even for x near 0}
begin
  if IsNanD(x) then begin
    log10p1 := x;
    exit;
  end;
  if abs(x) > 0.5 then  log10p1 := log10(1.0 + x)
  else log10p1 := log10e*ln1p(x);
end;


{---------------------------------------------------------------------------}
function logaddexp(x,y: double): double;
  {-Accurately compute ln[exp(x) + exp(y)]}
begin
  logaddexp := x + ln1pexp(y - x)
end;


{---------------------------------------------------------------------------}
function logsubexp(x,y: double): double;
  {-Accurately compute ln[exp(x) - exp(y)], x > y}
begin
  logsubexp := x + ln1mexp(y - x)
end;


{---------------------------------------------------------------------------}
function logistic(x: double): double;
  {-Return logistic(x) = 1/(1+exp(-x))}
const
  x0 = 38.0;
  x1 = -19.0;
begin
  if IsNanD(x) then begin
    logistic := x;
    exit;
  end;
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
function logit(x: double): double;
  {-Return logit(x) = ln(x/(1.0-x)), accurate near x=0.5}
const
  xl = 1/3;  {logit(1/3) = -ln2}  {FIX311}
  xh = 2/3;  {logit(2/3) =  ln2}
begin
  if (x < xl) or (x > xh) then logit := ln(x/(1.0-x))
  else logit := ln1p((2.0*x-1.0)/(1.0-x));
end;


{---------------------------------------------------------------------------}
function sinc(x: double): double;
  {-Return the cardinal sine sinc(x) = sin(x)/x}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNanD(x) then sinc := x else sinc := 0.0;
    exit;
  end;
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
function sincPi(x: double): double;
  {-Return the normalised cardinal sine sincPi(x) = sin(Pi*x)/(Pi*x)}
var
  y: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNanD(x) then sincPi := x else sincPi := 0.0;
    exit;
  end;
  x := abs(x);
  y := Pi*x;
  if y < 6.15e-4 then sincPi := sinc(y)
  else sincPi := sinPi(x)/y;
end;


{---------------------------------------------------------------------------}
function sinh_small(x: double; rel: boolean): double;
  {-Internal: return sinh(x) for |x| <= 1}
const
  ncs=8;
  csh: array[0..ncs-1] of THexDblW = ({chebyshev(sinh(x)/x-1, x=-1..1, 0.1e-20)}
         ($5F5A,$2221,$263F,$3FC6),   {+0.173042194047179631675883846985    }
         ($110A,$2E7A,$6C93,$3FB6),   {+0.875942219227604771549002634544e-1 }
         ($FD16,$6F43,$AFA8,$3F51),   {+0.107947777456713275024272706516e-2 }
         ($1CDA,$38CD,$BCF1,$3EDA),   {+0.637484926075475048156855545659e-5 }
         ($60DB,$94D2,$A5D1,$3E57),   {+0.220236640492305301591904955350e-7 }
         ($C5DC,$D7D0,$6BE6,$3DCB),   {+0.498794018041584931494262802066e-10}
         ($682C,$C8B3,$7130,$3D36),   {+0.797305355411573048133738827571e-13}
         ($2D3B,$1676,$4DF5,$3C9B));  {+0.947315871307254445124531745434e-16}
        {($4C8F,$EB9E,$A89D,$3BF9)}   {+0.869349205044812416919840906342e-19}
var
  css: array[0..ncs-1] of double absolute csh;
var
  t: double;
const
  t1 = 1E-8;
begin
  if abs(x)<=t1 then begin
    {sinh(x) = x*(1 + 1/6*x^2 + 1/120*x^4 + O(x^6))}
    if rel then sinh_small := 1.0
    else sinh_small := x;
  end
  else begin
    t := CSEvalD(2.0*sqr(x)-1.0, css, ncs);
    if rel then sinh_small := 1.0 + t
    else sinh_small := x + x*t;
  end;
end;


{---------------------------------------------------------------------------}
function sinh(x: double): double;
  {-Return the hyperbolic sine of x, accurate even for x near 0}
var
  t: double;
const
  t0 = 19.0;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    sinh := x;
    exit;
  end;
  t := abs(x);
  if t<1.0 then sinh := sinh_small(x, false)
  else begin
    {sinh(t) = 0.5*(exp(t) - exp(-t))}
    if t<=t0 then begin
      {calculate inverse only if it is not too small}
      t := exp(t);
      t := 0.5*(t - 1.0/t);
    end
    else if t<ln_MaxDbl then begin
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
function lnsinh(x: double): double;
  {-Return ln(sinh(x)), x > 0, accurate for x ~ 0 and without overflow for large x}
const
  xl = 19.0;
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNanD(x) or (x<=0) then lnsinh := Nan_d else lnsinh := PosInf_D;
    exit;
  end;
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
function sinhc(x: double): double;
  {-Return sinh(x)/x, accurate even for x near 0}
var
  t: double;
begin
  if IsNanD(x) then begin
    sinhc := Nan_d;
    exit;
  end;
  t := abs(x);
  if abs(t) < 1.0 then sinhc := sinh_small(t, true)
  else if t < ln_MaxDbl then sinhc := sinh(t)/t
  else begin
    t := t - ln(2.0*t);
    sinhc := exp(t);
  end;
end;


{---------------------------------------------------------------------------}
function sinhmx(x: double): double;
  {-Return sinh(x)-x, accurate even for x near 0}
var
  t: double;
const
  ncs=9;
  csh: array[0..ncs-1] of THexDblW = (  {chebyshev((sinh(x)-x)/x^3, x=-2..2, 0.1e-20)}
         ($1F45,$4BC7,$A057,$3FD7),  {+0.369161437989977107646200496798    } {*2}
         ($E308,$FFDB,$C797,$3F92),  {+0.183395147241496573111025588024e-1 }
         ($5601,$C438,$449F,$3F3C),  {+0.431336408128174581497124393209e-3 }
         ($F99B,$E39F,$DC2A,$3ED8),  {+0.592709289470722070924036619944e-5 }
         ($2469,$8F3B,$ADDB,$3E6C),  {+0.534190451019238924481986729972e-7 }
         ($F8AA,$F1F6,$5E4F,$3DF7),  {+0.340055083020092921361501610500e-9 }
         ($BE19,$BBAE,$5382,$3D7C),  {+0.161015882323068379168404355934e-11}
         ($F191,$3072,$8913,$3CFA),  {+0.589205330187637306105897428969e-14}
         ($AFAA,$14F6,$C8F9,$3C73)); {+0.171607959509371641280224879040e-16}
        {($1021,$C27E,$09EE,$3BE8)}  {+0.407233102668331851851851851852e-19}

var
  csv: array[0..ncs-1] of double absolute csh;
begin
  if IsNanOrInfD(x) then begin
    sinhmx := Nan_d;
    exit;
  end;
  t := abs(x);
  if t < 2.0 then begin
    if t <= sqrt_epsh then begin
      {sinh(x)-x = x^3/6*(1 + 1/20*x^2 + O(x^4))}
      sinhmx := x*x*x/SIXX;
    end
    else begin
      t := CSEvalD(0.5*sqr(t) - 1.0, csv, ncs);
      sinhmx := x*x*x*t;
    end
  end
  else sinhmx := sinh(x)-x;
end;


{---------------------------------------------------------------------------}
function cosh(x: double): double;
  {-Return the hyperbolic cosine of x}
const
  x0 = 19.0;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNanD(x) then cosh := x else cosh := PosInf_D;
    exit;
  end;
  if x<=1.0 then cosh := 1.0 + coshm1(x)
  else begin
    if x>x0 then begin
      {x>x0: exp(x) + exp(-x) ~ exp(x)}
      if x<ln_MaxDbl then cosh := 0.5*exp(x)
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
function lncosh(x: double): double;
  {-Return ln(cosh(x)), accurate for x ~ 0 and without overflow for large x}
const
  xl = 19.0;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNanD(x) then lncosh := x else lncosh := PosInf_D;
    exit;
  end;
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
function coshm1(x: double): double;
  {-Return cosh(x)-1, accurate even for x near 0}
const
  ncc=7;
  cch: array[0..ncc-1] of THexDblW = (     {chebyshev(((cosh(x)-1)/(x^2)-1/2)/x^2, x=-1..1, 0.1e-20);}
        ($9250,$01DB,$B196,$3FB5),   {+0.847409967933085972383472682702e-1 }
        ($AD1C,$F68D,$2A89,$3F47),   {+0.706975331110557282636312722683e-3 }
        ($FE59,$842A,$7192,$3ECA),   {+0.315232776534872632694769980081e-5 }
        ($EECC,$B813,$C699,$3E42),   {+0.874315531301413118696907944019e-8 }
        ($F57B,$10DE,$2E58,$3DB2),   {+0.165355516210397415961496829471e-10}
        ($3DAC,$0F99,$8AAD,$3D19),   {+0.226855895848535933261822520050e-13}
        ($64A2,$AE11,$372E,$3C7B));  {+0.236057319781153495473072617928e-16}
       {($4FEF,$47C9,$BF85,$3BD6)}   {+0.192684134365813759102742156829e-19}
var
  ccs: array[0..ncc-1] of double absolute cch;
  z: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNanD(x) then coshm1 := x else coshm1 := PosInf_D;
    exit;
  end;
  x := abs(x);
  if x<=1.0 then begin
    x := sqr(x);
    if x<eps_d then begin
      {cosh(x)-1 = 0.5*x^2*[1 + 1/12*x^2 + O(x^4)]}
      coshm1 := 0.5*x;
    end
    else begin
      z := CSEvalD(2.0*x-1.0, ccs, ncc);
      coshm1 := 0.5*x + x*x*z;
    end;
  end
  else coshm1 := cosh(x)-1.0;
end;

(*
function sinh(x: double): double;
var
  s,c: double;
begin
  sinhcosh(x,s,c);
  sinh := s;
end;


function cosh(x: double): double;
var
  s,c: double;
begin
  sinhcosh(x,s,c);
  cosh := c;
end;
*)


{---------------------------------------------------------------------------}
procedure sinhcosh(x: double; var s,c: double);
  {-Return s=sinh(x) and c=cosh(x)}
var
  t: double;
const
  t0 = 19.0;       {~ 0.5*ln(2^-54)}
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
    else if t<ln_MaxDbl then begin
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
function tanh(x: double): double;
  {-Return the hyperbolic tangent of x, accurate even for x near 0}
var
  t,z: double;
const
  t0 = 19.0;       {~ 0.5*ln(2^-54)}
  t1 = 0.86e-4;    {~ 2^(-54/4)}
  t2 = 0.7;
  PH: array[0..3] of THexDblW = (
        ($F7BC,$9797,$702B,$C094),   {-1.3080425704712825945553e+3}
        ($EBD7,$AAFD,$036D,$C055),   {-8.4053568599672284488465e+1}
        ($222B,$98F2,$9C53,$BFEE),   {-9.5658283111794641589011e-1}
        ($619A,$E2A3,$F331,$BF11));  {-6.8473739392677100872869e-5}
  QH: array[0..3] of THexDblW = (
        ($F39B,$6363,$A841,$40AE),   { 3.9241277114138477845780e+3}
        ($2BE7,$45F9,$773F,$409C),   { 1.8218117903645559060232e+3}
        ($C48D,$AD99,$109B,$4058),   { 9.6259501838840336946872e+1}
        ($0000,$0000,$0000,$3FF0));  { 1.0000000000000000000000e+0}
var
  P: array[0..3] of double absolute PH;
  Q: array[0..3] of double absolute QH;
begin
  if IsNanD(x) then begin
    tanh := x;
    exit;
  end;
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
      t := PolEval(z,P,4)/PolEval(z,Q,4);
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
function power(x, y : double): double;
  {-Return x^y; if frac(y)<>0 then x must be > 0}

  {This is my Pascal translation of the ldouble/powl.c file from the  }
  {freely distributable Cephes Mathematical Library V2.8, (June 2000) }
  {Copyright 1984, 1991, 1998 by Stephen L. Moshier                   }

  {Following Cody and Waite, this function uses a lookup table of     }
  {2**-i/32 and pseudo extended precision arithmetic to obtain several}
  {extra bits of accuracy in both the logarithm and the exponential.  }

const
  NXT = 128; {table size}

const
  MEXP   =  NXT*1024.0;
  MNEXP  = -NXT*1074.0;

const
  LOG2EA : THexDblW = ($0BF8,$94AE,$551D,$3FDC); { = 0.44269504088896340736 = log2(e)-1}

const
  {Coefficients for log(1+x) =  x - .5x^2 + x^3 *  P(z)/Q(z) }
  {on the domain  2^(-1/32) - 1  <=  x  <=  2^(1/32) - 1     }
  P: array[0..3] of THexDblW = (
       ($16F7,$DE95,$4D58,$3F4B),  {+8.3319510773868690322E-4}
       ($E050,$1819,$5C2B,$3FDF),  {+4.9000050881978030048E-1}
       ($E448,$F926,$000C,$3FFC),  {+1.7500123722550302574}
       ($687A,$F94D,$6670,$3FF6)); {+1.4000100839971580946}

  Q: array[0..2] of THexDblW = (
       ($8D2C,$6674,$0007,$4015),  {+5.2500282295834885815}
       ($5B00,$A38C,$CCD4,$4020),  {+8.4000598057587012590}
       ($0E5B,$BAFA,$CCD4,$4010)); {+4.2000302519914738397}

const
  {Coefficients for 2^x = 1 + x P(x),  on the interval -1/32 <= x <= 0}
  R: array[0..6] of THexDblW = (
       ($61D5,$C3AA,$A55D,$3EEF),  {+1.5089970579127660196E-5}
       ($CFD9,$2C11,$304B,$3F24),  {+1.5402715328927013322E-4}
       ($BB4C,$AD55,$D87F,$3F55),  {+1.3333556028915670086E-3}
       ($B338,$6F9F,$B2AB,$3F83),  {+9.6181291046036759829E-3}
       ($93BC,$D704,$6B08,$3FAC),  {+5.5504108664798462724E-2}
       ($C58C,$FF82,$BFBD,$3FCE),  {+2.4022650695910063856E-1}
       ($39EF,$FEFA,$2E42,$3FE6)); {+6.9314718055994528623E-1}

const
  {A[i] = 2^(-i/32), rounded to IEEE double precision. }
  {If i is even, A[i]+B[i/2] gives additional accuracy.}
  A: array[0..NXT] of THexDblW = (
       ($0000,$0000,$0000,$3ff0),  ($71f1,$2b8f,$d3c2,$3fef),  ($90d8,$819e,$a7c1,$3fef),
       ($be14,$ad9c,$7bfd,$3fef),  ($4540,$5b6e,$5076,$3fef),  ($ba97,$376b,$252b,$3fef),
       ($5a27,$ee61,$fa1b,$3fee),  ($67f1,$2d8e,$cf48,$3fee),  ($90da,$a2a4,$a4af,$3fee),
       ($4c83,$fbc7,$7a51,$3fee),  ($3ff6,$e78b,$502e,$3fee),  ($a129,$14f5,$2646,$3fee),
       ($9b5f,$337b,$fc97,$3fed),  ($b460,$f301,$d321,$3fed),  ($3285,$03db,$a9e6,$3fed),
       ($8398,$16c9,$80e3,$3fed),  ($a487,$dcfb,$5818,$3fed),  ($89f2,$080d,$2f87,$3fed),
       ($897c,$4a07,$072d,$3fed),  ($c3fa,$555d,$df0b,$3fec),  ($9069,$dcef,$b720,$3fec),
       ($e7b5,$9406,$8f6d,$3fec),  ($d14b,$2e57,$67f1,$3fec),  ($d07a,$5fff,$40ab,$3fec),
       ($529c,$dd85,$199b,$3fec),  ($1e09,$5bd7,$f2c2,$3feb),  ($c1d2,$904b,$cc1e,$3feb),
       ($064a,$30a1,$a5b0,$3feb),  ($5e47,$f2fb,$7f76,$3feb),  ($593a,$8de5,$5972,$3feb),
       ($15fb,$b84f,$33a2,$3feb),  ($b666,$298d,$0e07,$3feb),  ($d3ad,$995a,$e89f,$3fea),
       ($f37a,$bfd3,$c36b,$3fea),  ($fdbf,$5579,$9e6b,$3fea),  ($b358,$1330,$799e,$3fea),
       ($255d,$b23e,$5503,$3fea),  ($2d33,$ec4a,$309b,$3fea),  ($e565,$7b5d,$0c66,$3fea),
       ($2323,$19e3,$e863,$3fe9),  ($f090,$82a3,$c491,$3fe9),  ($07ba,$70ca,$a0f1,$3fe9),
       ($4e50,$9fde,$7d82,$3fe9),  ($520f,$cbc8,$5a44,$3fe9),  ($c5e5,$b0cd,$3737,$3fe9),
       ($ffc6,$0b91,$145b,$3fe9),  ($7736,$9915,$f1ae,$3fe8),  ($448c,$16b5,$cf32,$3fe8),
       ($a0db,$422a,$ace5,$3fe8),  ($6699,$d98a,$8ac7,$3fe8),  ($92ed,$9b44,$68d9,$3fe8),
       ($c7ad,$4623,$471a,$3fe8),  ($ce13,$994c,$2589,$3fe8),  ($1a12,$543e,$0427,$3fe8),
       ($4e62,$36cf,$e2f3,$3fe7),  ($c132,$0130,$c1ed,$3fe7),  ($0187,$73eb,$a114,$3fe7),
       ($5d3f,$4fde,$8069,$3fe7),  ($67c9,$5642,$5feb,$3fe7),  ($8174,$48a5,$3f9a,$3fe7),
       ($5f74,$e8ec,$1f75,$3fe7),  ($9484,$f951,$ff7d,$3fe6),  ($1a2f,$3c65,$dfb2,$3fe6),
       ($dabf,$750b,$c012,$3fe6),  ($3bcd,$667f,$a09e,$3fe6),  ($a973,$d44c,$8155,$3fe6),
       ($2225,$8255,$6238,$3fe6),  ($c320,$34cc,$4346,$3fe6),  ($5585,$b03a,$247e,$3fe6),
       ($dc09,$b976,$05e1,$3fe6),  ($2148,$15ad,$e76f,$3fe5),  ($46b7,$8a59,$c926,$3fe5),
       ($5429,$dd48,$ab07,$3fe5),  ($c7fd,$d497,$8d12,$3fe5),  ($27da,$36b5,$6f47,$3fe5),
       ($920f,$ca5d,$51a4,$3fe5),  ($4f82,$569d,$342b,$3fe5),  ($6642,$a2cf,$16da,$3fe5),
       ($2ca7,$769d,$f9b2,$3fe4),  ($dd0d,$99fd,$dcb2,$3fe4),  ($2a27,$d536,$bfda,$3fe4),
       ($d3de,$f0d7,$a32a,$3fe4),  ($3cd0,$b5c1,$86a2,$3fe4),  ($0057,$ed1d,$6a41,$3fe4),
       ($892d,$6061,$4e08,$3fe4),  ($a897,$d950,$31f5,$3fe4),  ($2e2a,$21f7,$160a,$3fe4),
       ($801c,$04ac,$fa45,$3fe3),  ($3422,$4c12,$dea6,$3fe3),  ($a8e5,$c313,$c32d,$3fe3),
       ($9ff7,$34e5,$a7db,$3fe3),  ($d866,$6d05,$8cae,$3fe3),  ($a9cb,$373a,$71a7,$3fe3),
       ($9ff1,$5f92,$56c5,$3fe3),  ($16ff,$b264,$3c08,$3fe3),  ($d831,$fc4c,$2170,$3fe3),
       ($b715,$0a31,$06fe,$3fe3),  ($2f56,$a93e,$ecaf,$3fe2),  ($030b,$a6e4,$d285,$3fe2),
       ($d990,$d0da,$b87f,$3fe2),  ($dee1,$f51f,$9e9d,$3fe2),  ($6381,$e1f5,$84df,$3fe2),
       ($7cdd,$65e2,$6b45,$3fe2),  ($a63f,$4fb2,$51ce,$3fe2),  ($6238,$6e75,$387a,$3fe2),
       ($dc96,$917d,$1f49,$3fe2),  ($8cd6,$8862,$063b,$3fe2),  ($d91d,$22fc,$ed50,$3fe1),
       ($b9aa,$3168,$d487,$3fe1),  ($5cd4,$8404,$bbe0,$3fe1),  ($cb75,$eb6f,$a35b,$3fe1),
       ($8dea,$388c,$8af9,$3fe1),  ($517b,$3c7d,$72b8,$3fe1),  ($8e51,$c8a5,$5a98,$3fe1),
       ($2de0,$aea9,$429a,$3fe1),  ($31cc,$c06c,$2abd,$3fe1),  ($5b51,$d012,$1301,$3fe1),
       ($d31b,$affe,$fb66,$3fe0),  ($d1a2,$32d3,$e3ec,$3fe0),  ($47f7,$2b72,$cc92,$3fe0),
       ($890f,$6cf9,$b558,$3fe0),  ($f383,$cac6,$9e3e,$3fe0),  ($9bc8,$1875,$8745,$3fe0),
       ($f6de,$29dd,$706b,$3fe0),  ($8574,$d315,$59b0,$3fe0),  ($7f85,$e86e,$4315,$3fe0),
       ($8061,$3e77,$2c9a,$3fe0),  ($3335,$a9fb,$163d,$3fe0),  ($0000,$0000,$0000,$3fe0));

  B: array[0..NXT div 2] of THexDblW = (
       ($0000,$0000,$0000,$0000),  ($931e,$f3a5,$4853,$3c77),  ($a18b,$2dd8,$d3e1,$3c89),
       ($b6b0,$86a4,$c7f4,$3c8d),  ($2893,$179c,$9c23,$bc8e),  ($cc8f,$80a9,$9e89,$3c73),
       ($4b5c,$4f18,$a5cd,$bc81),  ($b532,$696d,$2300,$3c8c),  ($3707,$d75b,$ed02,$3c72),
       ($7a9c,$4379,$bc37,$bc8c),  ($49db,$d1e9,$03cb,$3c65),  ($3cad,$ff48,$884d,$3c82),
       ($48dd,$8950,$1065,$3c71),  ($9e84,$7a2d,$3dd0,$3c72),  ($ac3b,$7e54,$584f,$bc65),
       ($d708,$3084,$805e,$bc52),  ($cc81,$345d,$a1cd,$3c87),  ($fd31,$0ef7,$fac9,$3c80),
       ($41e1,$db8d,$2f6e,$bc8d),  ($d533,$5d1c,$5949,$bc83),  ($f2be,$b071,$7c46,$3c6c),
       ($85d1,$7c1b,$185b,$bc8d),  ($7ebc,$81b5,$5fc7,$bc57),  ($976c,$a2e3,$cc13,$3c75),
       ($4b27,$5686,$e9f1,$3c86),  ($f6ba,$9bd4,$c6f8,$bc8f),  ($32d8,$d415,$4c1d,$bc8d),
       ($797e,$ba15,$5d02,$3c60),  ($992f,$ee04,$1577,$bc74),  ($6dd3,$5731,$2459,$bc80),
       ($7a99,$8688,$6e47,$bc71),  ($88ab,$683c,$be3a,$bc5b),  ($6456,$13b2,$dd34,$bc8b),
       ($1c34,$8759,$b609,$bc8b),  ($b497,$7e40,$83c1,$bc83),  ($e65e,$3080,$a6f9,$3c8b),
       ($47ad,$0546,$324c,$3c86),  ($93ad,$011d,$bb2c,$3c89),  ($3cad,$1db1,$7abe,$bc70),
       ($57e3,$d259,$b309,$bc84),  ($42e2,$afec,$4397,$3c6d),  ($62f0,$b690,$c1a3,$3c63),
       ($80d0,$04ef,$9b7a,$3c38),  ($9278,$1c30,$f369,$bc4e),  ($9ebc,$11f0,$da09,$3c7a),
       ($f5e3,$d661,$e436,$bc65),  ($eae2,$bf42,$3aea,$bc86),  ($59a6,$8436,$2721,$3c83),
       ($82e4,$d231,$f46a,$3c76),  ($41d5,$54db,$0247,$3c80),  ($1255,$afad,$12e8,$3c76),
       ($e9d9,$9940,$bd33,$3c72),  ($0573,$b6c7,$b07e,$3c89),  ($8495,$814a,$c775,$3c7d),
       ($643c,$00a2,$016e,$3c8e),  ($68e4,$7b49,$5b4c,$3c7e),  ($8a76,$b9d7,$9041,$bc71),
       ($369e,$9af1,$2fbf,$bc83),  ($9b3a,$3944,$c510,$bc86),  ($7b53,$27c5,$3a17,$3c30),
       ($610b,$4adc,$a62e,$3c88),  ($84ff,$4bb2,$86be,$3c51),  ($b465,$a475,$73e2,$3c7d),
       ($085d,$535b,$9083,$bc61),  ($0000,$0000,$0000,$0000));

var
  w, z, ya, yb, u: double;
  F, Fa, Fb, G, Ga, Gb, H, Ha, Hb: double;
  e,i : longint;
  k: integer;
  nflg : boolean;
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) or (THexDblW(y)[3] and $7FF0=$7FF0) then begin
    power := NaN_d;
    exit;
  end;

  {Handle easy cases}
  w := abs(y);
  if y=0.0 then begin
    power := 1.0;
    exit;
  end
  else if x=0.0 then begin
    if y>0.0 then power := 0.0
    else begin
      {here y<0: if y is odd and x=-0 return -INF}
      {if frac(0.5*y)=0.0 then power := PosInf_d
      else power := copysignd(PosInf_d,x);
      }
      {generate error}
      power := 1.0/x;
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
      if abs(x) > 0.25*MinDouble then power := 1.0/x
      else power := copysignd(PosInf_d,x);
    end;
    exit;
  end
  else if abs(y)=0.5 then begin
    if y>0.0 then power := sqrt(x)
    else power := 1.0/sqrt(x);
    exit;
  end;

  nflg := false;         {true if x<0 raised to integer power }
  w := floord(y);
  if w=y then begin
    z := abs(w);
    if z<512.0 then begin
      u := ilogb(x);
      {intpower0 does not catch overflows, therefore call it only if safe.}
      if z*(abs(u) + 1.0) < 1022.0 then begin
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
    w := 2.0*floord(0.5*y);
    nflg := w<>y;
    x := abs(x);
  end;

  {separate significand from exponent}
  frexpd(x, x, e);

  {find significand in antilog table A[]}
  i := 1;
  if x <= double(A[65])   then i := 65;
  if x <= double(A[i+32]) then inc(i,32);
  if x <= double(A[i+16]) then inc(i,16);
  if x <= double(A[i+8])  then inc(i, 8);
  if x <= double(A[i+4])  then inc(i, 4);
  if x <= double(A[i+2])  then inc(i, 2);
  if x >= double(A[1])    then i := -1;
  inc(i);

  {Find (x - A[i])/A[i] in order to compute log(x/A[i]):}
  {log(x) = log( a x/a ) = log(a) + log(x/a)            }
  {log(x/a) = log(1+v),  v = x/a - 1 = (x-a)/a          }
  x := x - double(A[i]);
  x := x - double(B[i div 2]);
  x := x / double(A[i]);

  {rational approximation for log(1+v) = v - v**2/2 + v**3 P(v) / Q(v)}
  z := x*x;

  {w:= polevl(x,P,3)}
  w := double(P[0]);  for k:=1 to 3 do w := w*x + double(P[k]);

  {u = p1evl(x,Q,3)}
  u := ((x + double(Q[0]))*x + double(Q[1]))*x + double(Q[2]);

  w := x*(z*w/u) - 0.5*z;

  {Convert to base 2 logarithm: multiply by log2(e) = 1+LOG2EA}
  z := double(LOG2EA)*w;
  z := z + w;
  z := z + double(LOG2EA)*x;
  z := z + x;

  {Compute exponent term of the base 2 logarithm.}
  w := -i;
  w := w/NXT + e;

  {Now base 2 log of x is w + z.  Multiply base 2 log by y, in extended}
  {precision. Separate y into large part ya and small part yb less than 1/NXT}

  ya := floord(y*NXT)/NXT;    {reduc(y)}
  yb := y - ya;

  F  := z*y + w*yb;
  Fa := floord(F*NXT)/NXT;    {reduc(F)}
  Fb := F - Fa;

  G  := Fa + w*ya;
  Ga := floord(G*NXT)/NXT;    {reduc(G)}
  Gb := G - Ga;

  H  := Fb + Gb;
  Ha := floord(H*NXT)/NXT;    {reduc(H)}
  w  := NXT*(Ga+Ha);

  {Test the power of 2 for over/underflow}
  if w>MEXP then begin
    if nflg then power := NegInf_d
    else power := PosInf_d;
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
  z := double(R[0]);
  for k:=1 to 6 do z := z*Hb + double(R[k]);
  z := Hb*z;

  {Express e/NXT as an integer plus a negative number of (1/NXT)ths.}
  {Find lookup table entry for the fractional power of 2.    }
  if e<0 then i := 0 else i := 1;
  inc(i, e div NXT);
  e := NXT*i - e;

  w := double(A[e]);
  z := w + w*z;     {2**-e * ( 1 + (2**Hb-1) )     }
  z := ldexpd(z,i);  {multiply by integer power of 2}

  if nflg then power := -z {odd integer exponent and x<0}
  else power := z;

end;


{---------------------------------------------------------------------------}
function powm1(x,y: double): double;
  {-Return x^y - 1; special code for small x,y}
var
  p: double;
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) or (THexDblW(y)[3] and $7FF0=$7FF0) then begin
    powm1 := NaN_d;
    exit;
  end;
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
function pow1pm1(x,y: double): double;
  {-Return (1+x)^y - 1; special code for small x,y}
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) or (THexDblW(y)[3] and $7FF0=$7FF0) then begin
    pow1pm1 := NaN_d;
    exit;
  end;
  if (x=0.0) or (y=0.0) then pow1pm1 := 0.0
  else if y=1.0 then pow1pm1 := x
  else if y=2.0 then pow1pm1 := x*(2.0+x)
  else if abs(x)<1.0 then pow1pm1 := expm1(y*ln1p(x))
  else pow1pm1 := powm1(1.0+x,y);
end;


{---------------------------------------------------------------------------}
function pow1p_ex(x,y: double; usexx: boolean): double;
  {-Return (1+x)^y, x > -1}
var
  z,w,v: double;
  xx: dbl2;
const
  vmax = 690; {<ln(MaxDouble/2^28): to allow dbl2 splitting}
begin
  if IsNanD(x) or IsNanD(y) then begin
    pow1p_ex := Nan_d;
    exit;
  end;
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
    {$ifdef FPC}
      xx.h := 0.0;   {Avoid silly warning}
    {$endif}
    ddto2d(1.0,x,xx);
    pow2di(xx,int(y),xx);
    if w<>0.0 then begin
      w := pow1p(x,w);
      mul21d(xx,w,xx);
    end;
    pow1p_ex := xx.h+xx.l;
  end
  else pow1p_ex := power(1.0+x, y);
end;


{---------------------------------------------------------------------------}
function pow1p(x,y: double): double;
  {-Return (1+x)^y, x > -1, with dbl2 arithmetic for critical values}
begin
  pow1p := pow1p_ex(x,y,true);
end;


{---------------------------------------------------------------------------}
function pow1pf(x,y: double): double;
  {-Return (1+x)^y, x > -1, without dbl2, less accurate than pow1p}
begin
  pow1pf := pow1p_ex(x,y,false);
end;


{---------------------------------------------------------------------------}
function powpi2k(k,n: longint): double;
  {-Return accurate scaled powers of Pi, result = (Pi*2^k)^n}
var
  m: longint;
  x: double;
const
  lph  : double = 13529/8192;       {high part log2(Pi) = 1.6514892578125}
  lpl  = 0.6871659818798043279e-5;  {low  part log2(Pi)}
  lmax = +1024;
  lmin = -1077;
begin
  x := (k + (lph+lpl))*n;
  if x >= lmax then powpi2k := PosInf_d
  else if x < lmin then powpi2k := 0.0
  else begin
    x := n*lph;
    m := round(x);
    x := (x - m) + n*lpl;
    m := m + k*n;
    x := exp2(x);
    powpi2k := ldexpd(x,m);
  end;
end;


{---------------------------------------------------------------------------}
function powpi(n: longint): double;
  {-Return accurate powers of Pi, result = Pi^n}
begin
  powpi := powpi2k(0,n);
end;


{---------------------------------------------------------------------------}
function compound(x: double; n: longint): double;
  {-Return (1+x)^n; accurate version of Delphi/VP internal function}
var
  xx: dbl2;
begin
  {$ifdef FPC}
    xx.h := 0.0;   {Avoid silly warning}
  {$endif}
  ddto2d(1.0,x,xx);
  pow2di(xx,n,xx);
  compound := xx.h+xx.l;
end;


(*
{---------------------------------------------------------------------------}
function comprel_fast(x: double; n: double): double;
  {-Return the ((1+x)^n-1)/x}
begin
  {Faster comprel without dbl2, but there are large errors up to 512 eps. It}
  {is still much better than the Delphi implementation, see t_xdamat*.cmp   }
  if n=0 then comprel_fast := n
  else if abs(x) < sqrt_epsh/n then begin
    comprel_fast := (1.0 + 0.5*(n-1)*x)*n;
  end
  else comprel_fast := pow1pm1(x,n)/x;
end;
*)

{---------------------------------------------------------------------------}
function comprel(x,n: double): double;
  {-Return ((1+x)^n-1)/x; accurate version of Delphi/VP internal function}
var
  xx,tt: dbl2;
begin
  {Relative compound function: (compound(x,n)-1)/x, see math.annuity2}
  if IsNanOrInfD(x) then begin
    comprel := Nan_d;
    exit;
  end;
  if x=0.0 then comprel := n
  else begin
    {$ifdef FPC}
      xx.h := 0.0;   {Avoid silly warning}
    {$endif}
    ddto2d(1.0,x,xx);
    pow2di(xx,round(n),xx);
    tt.h := 1.0;
    tt.l := 0.0;
    sub2d(xx,tt,xx);
    tt.h := x;
    tt.l := 0.0;
    div2d(xx,tt,xx);
    comprel := xx.h+xx.l;
  end;
end;


{---------------------------------------------------------------------------}
function arccos(x: double): double;
  {-Return the inverse circular cosine of x, |x| <= 1}
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) then begin
    {Inf or NaN}
    arccos := Nan_D;
    exit;
  end;
  {basic formula arccos(x) = arctan(sqrt(1-x^2)/x))}
  if abs(x)=1.0 then begin
    if x<0.0 then arccos := Pi else arccos := 0.0;
  end
  else arccos := arctan2(sqrt((1.0-x)*(1.0+x)),x)
end;


{---------------------------------------------------------------------------}
function arccos1m(x: double): double;
  {-Return arccos(1-x), 0 <= x <= 2, accurate even for x near 0}
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) then begin
    {Inf or NaN}
    arccos1m := Nan_D;
    exit;
  end;
  if x>=0.5 then arccos1m := arccos(1.0-x)
  else begin
    {use arccos(z) = 2arcsin(sqrt((1-z)/2)), ref: }
    {http://functions.wolfram.com/01.13.06.0043.01}
    arccos1m := 2.0*arcsin(sqrt(0.5*x));
  end;
end;


{---------------------------------------------------------------------------}
function archav(x: double): double;
  {-Return the inverse haversine archav(x), 0 <= x <= 1}
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) then begin
    {Inf or NaN}
    archav := Nan_D;
    exit;
  end;
  {archav(x) = 2*arcsin(sqrt(x)) = 2*arctan2(sqrt(x), sqrt(1-x))}
  archav := 2.0*arctan2(sqrt(x), sqrt(1.0-x));
end;


{---------------------------------------------------------------------------}
function arcsin(x: double): double;
  {-Return the inverse circular sine of x, |x| <= 1}
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) then begin
    {Inf or NaN}
    arcsin := Nan_D;
    exit;
  end;
  {basic formula arcsin(x) = arctan(x/sqrt(1-x^2))}
  arcsin := arctan2(x, sqrt((1.0-x)*(1.0+x)))
end;


{---------------------------------------------------------------------------}
function sec(x: double): double;
  {-Return the circular secant of x, x mod Pi <> Pi/2}
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) then begin
    {Inf or NaN}
    sec := Nan_D;
    exit;
  end;
  sec := 1.0/cos(x);
end;


{---------------------------------------------------------------------------}
function csc(x: double): double;
  {-Return the circular cosecant of x, x mod Pi <> 0}
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) then begin
    {Inf or NaN}
    csc := Nan_D;
    exit;
  end;
  csc := 1.0/sin(x);
end;


{---------------------------------------------------------------------------}
function coth(x: double): double;
  {-Return the hyperbolic cotangent of x, x<>0}
begin
  if IsNanD(x) then coth := x
  else coth := 1.0/tanh(x);
end;


{---------------------------------------------------------------------------}
function sech(x: double): double;
  {-Return the hyperbolic secant of x}
begin
  if IsNanD(x) then sech := x
  else if abs(x) > ln_MaxDbl then sech := 0.0
  else sech := 1.0/cosh(x);
end;


{---------------------------------------------------------------------------}
function csch(x: double): double;
  {-Return the hyperbolic cosecant of x, x<>0}
begin
  if IsNanD(x) then csch := x
  else if abs(x) > ln_MaxDbl then csch := 0.0
  else csch := 1.0/sinh(x);
end;


{---------------------------------------------------------------------------}
function arcsinh(x: double): double;
  {-Return the inverse hyperbolic sine of x}
var
  t,z: double;
const
  t0 = 0.14e8;    {sqrt(t0^2 + 1) = t0}
  t1 = 0.89e308;  {t1 < MaxDouble/2 }
const
  CSN = 20;
  CSAH: array[0..CSN-1] of THexDblW = ({chebyshev((arcsinh(x)/x-1), x=-1..1, 1e-20);}
          ($8AEF,$E4C5,$68DE,$BFC0),  {-0.128200399117381863433721273592    }
          ($ACAC,$3DDC,$1C93,$BFAE),  {-0.588117611899517675652117571383e-1 }
          ($0708,$56CA,$5D1B,$3F73),  {+0.472746543221248156407252497561e-2 }
          ($D59D,$69F9,$2E99,$BF40),  {-0.493836316265361721013601747903e-3 }
          ($E3DF,$F7A9,$AC91,$3F0E),  {+0.585062070585574122874948352586e-4 }
          ($0176,$DF2F,$51A0,$BEDF),  {-0.746699832893136813547550692112e-5 }
          ($9DD1,$59D3,$CBFD,$3EB0),  {+0.100116935835581992659661920174e-5 }
          ($BF1A,$A6B4,$A938,$BE82),  {-0.139035438587083336086164721962e-6 }
          ($76E6,$8ACA,$48F3,$3E55),  {+0.198231694831727935473173603337e-7 }
          ($B201,$19A7,$C7A0,$BE28),  {-0.288474684178488436127472674074e-8 }
          ($3B2F,$8386,$531B,$3DFD),  {+0.426729654671599379534575435248e-9 }
          ($40B0,$AF1D,$95EA,$BDD1),  {-0.639760846543663578687526705064e-10}
          ($CD7A,$821C,$5425,$3DA5),  {+0.969916860890647041478747328871e-11}
          ($5755,$9D38,$1D44,$BD7A),  {-0.148442769720437708302434249871e-11}
          ($718F,$7AF0,$1DF8,$3D50),  {+0.229037379390274479880187613939e-12}
          ($3758,$08B5,$08D4,$BD24),  {-0.355883951327326451644472656709e-13}
          ($24D1,$EA4E,$0ED1,$3CF9),  {+0.556396940800567901824410323292e-14}
          ($3259,$1E84,$82FE,$BCCF),  {-0.874625095996246917448307610618e-15}
          ($A766,$5CA6,$E8ED,$3CA3),  {+0.138152488445267754259259259259e-15}
          ($C80F,$2587,$44AA,$BC79)); {-0.219166882829054218109516486403e-16}
        (*($5B9C,$92A9,$18D0,$3C50),  {+0.349046585251352023737824855452e-17}
          ($F017,$9989,$94D1,$BC24),  {-0.557857884196310306067662486759e-18}
          ($033B,$C6ED,$6648,$3BFA),  {+0.894451477596418179416166864814e-19}
          ($FBD5,$4037,$FB20,$BBD0),  {-0.143834333235934314272184000000e-19}
          ($C0AC,$0977,$E78A,$3BA5),  {+0.231922385806700541867783333333e-20} *)
var
  CSA: array[0..CSN-1] of double absolute CSAH;
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) then begin
    {Inf or NaN}
    arcsinh := x;
    exit;
  end;
  t := abs(x);
  if t<=1.0 then begin
    if t<sqrt_epsh then begin
      {arcsinh(x) = x*(1 - 1/6*x^2 + 3/40*x^4 + O(x^6))}
      arcsinh := x;
    end
    else begin
      z := 2.0*x*x - 1.0;
      t := CSEvalD(z, CSA, CSN);
      arcsinh := x + x*t;
    end;
  end
  else begin
    if t>=t0 then begin
      {skip sqrt() because sqrt(t^2+1) = t}
      if t<=t1 then z := ln(2.0*t)
      else z := ln(t)+ln2
    end
    else z := ln(t + sqrt(1.0+t*t));
    if x>0.0 then arcsinh := z
    else arcsinh := -z;
  end;
end;


{---------------------------------------------------------------------------}
function arccosh(x: double): double;
  {-Return the inverse hyperbolic cosine, x >= 1. Note: for x near 1 the}
  { function arccosh1p(x-1) should be used to reduce cancellation errors!}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsInfD(x) and (x>0) then arccosh := x else arccosh := Nan_D;
    exit;
  end;
  if x=1.0 then arccosh := 0
  else begin
    if x>1.5E8 then begin
      {skip sqrt() calculation because sqrt(x^2-1) = x}
      arccosh := ln(x)+ln2;
    end
    else if x>2.0 then begin
      {arccosh := ln(x+sqrt((x-1)*(x+1)))}
      arccosh := ln(2.0*x - 1.0/(x+sqrt(x*x-1.0)))
    end
    else begin
      {arccosh = ln1p(y+sqrt(2*y + y*y)), y=x-1}
      x := x - 1.0;
      arccosh := ln1p(x+sqrt(2.0*x+x*x));
    end;
  end;
end;


{---------------------------------------------------------------------------}
function arccosh1p(x: double): double;
  {-Return arccosh(1+x), x>=0, accurate even for x near 0}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsInfD(x) and (x>0) then arccosh1p := x else arccosh1p := Nan_D;
    exit;
  end;
  if x=0.0 then arccosh1p := 0
  else if x<1.0 then begin
    {arccosh(X) = ln(X + sqrt(X^2 - 1)), substituting X = 1+x gives}
    {arccosh1p(x) = ln1p(x+sqrt((X-1)*(X+1))) = ln1p(x+sqrt(x*(x+2)))}
    {             = ln1p(x+sqrt(2*x + x*x))}
    arccosh1p := ln1p(x+sqrt(2.0*x+x*x));
  end
  else arccosh1p := arccosh(1.0+x);
end;


{---------------------------------------------------------------------------}
function arctanh(x: double): double;
  {-Return the inverse hyperbolic tangent of x, |x| < 1}
var
  t: double;
const
  t0 = 1E-4;
begin
  if IsNanD(x) then begin
    arctanh := x;
    exit;
  end;
  t := abs(x);
  if t<=t0 then begin
    {arctanh(x) = x + 1/3*x^3 + 1/5*x^5 + 1/7*x^7 + O(x^9)}
    arctanh := x*(1.0 + sqr(x)/3.0);
  end
  else begin
    {arctanh(x) = 0.5*ln((1+x)/(1-x))     = 0.5*ln((1-x+2x)/(1-x)) }
    {           = 0.5*ln(1+2x/(1-x))      = 0.5*ln1p(2x/(1-x))     }
    {  or       = 0.5*(ln(1+x)-ln(1-x))   = 0.5*(ln1p(x)-ln1p(-x)) }
    {  or       = 0.5*ln(1+2x+2x^2/(1-x)) = 0.5*ln1p(2x+2x^2/(1-x))}
    if t<0.5 then begin
      t := 0.5*(ln1p(t)-ln1p(-t));
    end
    else t := 0.5*ln((1.0+t)/(1.0-t));
    if x>0.0 then arctanh := t
    else arctanh := -t;
  end;
end;


{---------------------------------------------------------------------------}
function arccot(x: double): double;
  {-Return the sign symmetric inverse circular cotangent; arccot(x) = arctan(1/x), x <> 0}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNanD(x) then arccot := Nan_D
    else arccot := 0;
    exit;
  end;
  if abs(x) > 1E-17 then arccot := arctan(1.0/x)
  else begin
    if x>=0.0 then arccot := Pi_2
    else arccot := -Pi_2;
  end;
end;


{---------------------------------------------------------------------------}
function arccotc(x: double): double;
  {-Return the continuous inverse circular cotangent; arccotc(x) = Pi/2 - arctan(x)}
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    if IsNanD(x) then arccotc := Nan_D
    else begin
      if x>0.0 then arccotc := 0
      else arccotc := Pi;
    end;
    exit;
  end;
  if x > 1E-17 then arccotc := arctan(1.0/x)
  else arccotc := Pi_2 - arctan(x);
end;


{---------------------------------------------------------------------------}
function arccsc(x: double): double;
  {-Return the inverse cosecant of x, |x| >= 1}
var
  y,z: double;
begin
  if IsNanD(x) then begin
    {Inf or NaN}
    arccsc := Nan_D;
    exit;
  end;
  {arccsc = arcsin(1.0/x) = arctan(1/sqrt(x^2-1))}
  y := abs(x);
  if y >= 1E8 then begin
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
function arcsec(x: double): double;
  {-Return the inverse secant of x, |x| >= 1}
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    arcsec := Nan_D;
    exit;
  end;
  if abs(x) >= 1e8 then begin
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
function arccoth(x: double): double;
  {-Return the inverse hyperbolic cotangent of x, |x| > 1}
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    arccoth := Nan_D;
    exit;
  end;
  t := abs(x);
  {compute t = arccoth(|x|)}
  if t<32.0 then t := 0.5*ln1p(2.0/(t-1.0))
  else t := -0.5*ln1p(-2.0/(t+1.0));
  {adjust sign}
  if x>0.0 then arccoth := t
  else arccoth := -t;
end;


{---------------------------------------------------------------------------}
function arcsech(x: double): double;
  {-Return the inverse hyperbolic secant of x, 0 < x <= 1}
var
  t: double;
begin
  if THexDblW(x)[3] and $7FF0=$7FF0 then begin
    {Inf or NaN}
    arcsech := Nan_D;
    exit;
  end;
  t := 1.0-x;
  if t=0 then arcsech := 0.0
  else if x<=1E-9 then begin
    {arcsech(x) = (ln(2)-ln(x)) -1/4*x^2 -3/32*x^4 + O(x^6)}
    arcsech := -ln(0.5*x);
  end
  else begin
    {arcsech(x) = arccosh(1/x), see arccosh branch for x<=2}
    t := t/x;
    arcsech := ln1p(t+sqrt(2.0*t+t*t));
  end;
end;


{---------------------------------------------------------------------------}
function arccsch(x: double): double;
  {-Return the inverse hyperbolic cosecant of x, x <> 0}
var
  t: double;
begin
  if IsNanD(x) then begin
    {Inf or NaN}
    arccsch := x;
    exit;
  end;
  t := abs(x);
  if t<=1.0 then begin
    if t<=5e-9 then begin
      {avoid overflow for 1/t^2}
      {arccsch(x) = (ln(2)-ln(x)) + 1/4*x^2 -3/32*x^4 + O(x^6)}
      t := -ln(0.5*t)
    end
    else t := ln(1.0/t + sqrt(1.0 + 1.0/sqr(t)));
    if x>0.0 then arccsch := t
    else arccsch := -t;
  end
  else arccsch := arcsinh(1.0/x);
end;


{---------------------------------------------------------------------------}
{------------- Degrees versions of trig / invtrig functions ----------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure trig_deg(x: double; var y,z: double; var n: integer; var m45: boolean);
  {-Reduce x in degrees mod 90; y=x mod 90, |y|<45. z=x/45, m45 if x is multiple of 45}
const
  XMAX = 5e13; {~2^53/180}
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
    y := floord(z);
    n := trunc(y - 16.0*floord(y/16.0));
    if odd(n) then begin
      inc(n);
      y := y + 1.0;
    end;
    n := (n shr 1) and 7;
    y := x - y*c45;
  end;
end;


{---------------------------------------------------------------------------}
procedure sincosd(x: double; var s,c: double);
  {-Return sin(x) and cos(x), x in degrees}
var
  y,ss,cc: double;
  n: integer;
  m45: boolean;
begin
  if IsNanOrInfD(x) then begin
    s := Nan_d;
    c := Nan_d;
    exit;
  end;
  trig_deg(x,y,ss,n,m45);
  ss := system.sin(y);
  cc := system.cos(y);
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
function sind(x: double): double;
  {-Return sin(x), x in degrees}
var
  y,z: double;
  n  : integer;
  m45: boolean;
begin
  if IsNanOrInfD(x) then begin
    sind := Nan_d;
    exit;
  end;
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
function cosd(x: double): double;
  {-Return cos(x), x in degrees}
var
  y,z: double;
  n  : integer;
  m45: boolean;
begin
  if IsNanOrInfD(x) then begin
    cosd := Nan_d;
    exit;
  end;
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
function tand(x: double): double;
  {-Return tan(x), x in degrees}
var
  y,z: double;
  n  : integer;
  m45: boolean;
begin
  if IsNanOrInfD(x) then begin
    tand := Nan_d;
    exit;
  end;
  trig_deg(x,y,z,n,m45);
  if m45 then begin
    z := abs(z);
    y := isign(x);
    case round(4.0*frac(0.25*z)) of
        0: tand := 0.0;
        1: tand := y;
        2: tand := PosInf_d;
      else tand := -y;          {case 3, but keep some stupid compilers happy}
    end;
  end
  else if odd(n) then tand := 0.0 - _cot(y*Pi_180)
  else tand := _tan(y*Pi_180);
end;


{---------------------------------------------------------------------------}
function tanPi(x: double): double;
  {-Return tan(Pi*x), result will be 0 for abs(x) >= 2^52}
var
  c,s: double;
begin
  if IsNanOrInfD(x) then begin
    tanPi := Nan_d;
    exit;
  end;
  sincospi(x,s,c);
  if c=0.0 then tanPi := PosInf_d
  else tanPi := s/c;
end;


{---------------------------------------------------------------------------}
function cotd(x: double): double;
  {-Return cot(x), x in degrees}
var
  y,z: double;
  n  : integer;
  m45: boolean;
begin
  if IsNanOrInfD(x) then begin
    cotd := Nan_d;
    exit;
  end;
  trig_deg(x,y,z,n,m45);
  if m45 then begin
    z := abs(z);
    y := isign(x);
    case round(4.0*frac(0.25*z)) of
        0: cotd := PosInf_d;
        1: cotd := y;
        2: cotd := 0.0;
      else cotd := -y;          {case 3, but keep some stupid compilers happy}
    end;
  end
  else if odd(n) then cotd := 0.0 - _tan(y*Pi_180)
  else cotd := _cot(y*Pi_180);
end;


{---------------------------------------------------------------------------}
function arctand(x: double): double;
  {-Return the inverse circular tangent of x, result in degrees}
begin
  {exact for double +-1, +-INF; 90 for x>2^65, -90 for x < -2^65}
  arctand := arctan(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
function arccosd(x: double): double;
  {-Return the inverse circular cosine of x, |x| <= 1, result in degrees}
begin
  arccosd := arccos(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
function arccotd(x: double): double;
  {-Return the sign symmetric inverse circular cotangent,}
  { arccotd(x) = arctand(1/x), x <> 0, result in degrees}
begin
  arccotd := arccot(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
function arccotcd(x: double): double;
  {-Return the continuous inverse circular cotangent;}
  { arccotcd(x) = 90 - arctand(x), result in degrees}
begin
  arccotcd := arccotc(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
function arcsind(x: double): double;
  {-Return the inverse circular sine of x, |x| <= 1, result in degrees}
begin
  arcsind := arcsin(x)/Pi_180;
end;


{---------------------------------------------------------------------------}
{---------------------- Elementary numerical functions ---------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function cbrt(x: double): double;
  {-Return the cube root of x}
var
  y: double;
begin
  if (THexDblW(x)[3] and $7FF0=$7FF0) or (x=0.0) or (abs(x)=1.0) then begin
    {x is 0, +1, -1, Inf, or NaN}
    cbrt := x;
  end
  else begin
    {calculate initial approximation}
    y := exp(ln(abs(x))/THREE);
    if x<0.0 then y := -y;
    {perform one Newton step}
    cbrt := y - (y - x/sqr(y))/THREE;
  end;
end;


{---------------------------------------------------------------------------}
function ceil(x: double): longint;
  {-Return the smallest integer >= x; |x|<=MaxLongint}
var
  i: longint;
begin
  x := modf(x,i);
  if x>0.0 then ceil := i+1
  else ceil := i;
end;


{---------------------------------------------------------------------------}
function ceild(x: double): double;
  {-Return the smallest integer >= x}
var
  t: double;
begin
  t := int(x);
  if (x<=0.0) or (x=t) then ceild := t
  else ceild := t + 1.0
end;


{---------------------------------------------------------------------------}
function floor(x: double): longint;
  {-Return the largest integer <= x; |x|<=MaxLongint}
var
  i: longint;
begin
  x := modf(x,i);
  if x<0.0 then floor := i-1
  else floor := i;
end;


{---------------------------------------------------------------------------}
function floord(x: double): double;
  {-Return the largest integer <= x}
var
  t: double;
begin
  t := int(x);
  if (x>=0.0) or (x=t) then floord := t
  else floord := t - 1.0;
end;


{---------------------------------------------------------------------------}
function hypot(x,y: double): double;
  {-Return sqrt(x*x + y*y)}
begin
  x := abs(x);
  y := abs(y);
  if x>y then hypot := x*sqrt(1.0+sqr(y/x))
  else if x>0.0 then hypot := y*sqrt(1.0+sqr(x/y)) {here y >= x > 0}
  else hypot := y;                                 {here x=0}
end;


{---------------------------------------------------------------------------}
function hypot3(x,y,z: double): double;
  {-Return sqrt(x*x + y*y + z*z)}
var
  ax,ay,az,r,s,t: double;
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
function modf(x: double; var ip: longint): double;
  {-Return frac(x) and trunc(x) in ip, |x|<=MaxLongint}
begin
  ip := trunc(x);
  modf := x-ip
end;


{---------------------------------------------------------------------------}
function nroot(x: double; n: integer): double;
  {-Return the nth root of x; n<>0, x >= 0 if n is even}
var
  y: double;
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
    nroot := copysignd(exp(ln(y)/n),x);
  end;
end;


{---------------------------------------------------------------------------}
function sqrt1pm1(x: double): double;
  {-Return sqrt(1+x)-1, accurate even for x near 0, x>=-1}
begin
  if IsNanOrInfD(x) then begin
    sqrt1pm1 := Nan_d;
    exit;
  end;
  if abs(x)>0.75 then sqrt1pm1 := sqrt(1.0+x)-1.0
  else begin
    {sqrt(1+x)-1 = (sqrt(1+x)-1)*(sqrt(1+x)+1)/(sqrt(1+x)+1) = x/(sqrt(1+x)-1)}
    sqrt1pm1 := x/(1.0+sqrt(1.0+x));
  end;
end;


{---------------------------------------------------------------------------}
function sqrt1pmx(x: double): double;
  {-Return sqrt(1+x^2)-x}
begin
  if IsNanOrInfD(x) then begin
    sqrt1pmx := Nan_d;
    exit;
  end;
  if x <= 0.0 then sqrt1pmx := hypot(1.0,x)-x
  else sqrt1pmx := 1.0/(hypot(1.0,x)+x);
end;


{---------------------------------------------------------------------}
{---------------------- Floating point functions ---------------------}
{---------------------------------------------------------------------}


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
    {+-INF, NaN, 0}
    m := d;
  end
  else begin
    if w < $0010 then begin
      d := d*double(H2_54);
      e := -54;
      w := THexDblW(d)[3] and $7FF0;
    end;
    inc(e, (w shr 4) - 1022);
    m := d;
    THexDblW(m)[3] := (THexDblW(m)[3] and $800F) or $3FE0;
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
function ilogb(x: double): longint;
  {-Return base 2 exponent of x. For finite x ilogb = floor(log2(|x|))}
  { otherwise -MaxLongint for x=0 and MaxLongint if x = +-INF or Nan. }
var
  e: integer;
  f: longint;
  m: double;
begin
  e := (THexDblW(x)[3] and $7FF0) shr 4;
  if e=$7FF then ilogb := MaxLongint
  else if x=0.0 then ilogb := -MaxLongint
  else if e<>0  then ilogb := e-$3FF
  else begin
    {e=0, x<>0: denormal use frexpd exponent}
    frexpd(x,m,f);
    ilogb := f-1;
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
function IsNaND(d: double): boolean;   {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is a NaN}
begin
  with TDblRec(d) do begin
    IsNaND := (hm and $7FF00000=$7FF00000) and ((hm and $000FFFFF<>0) or (lm<>0));
  end;
end;


{---------------------------------------------------------------------------}
function IsNaNorInfD(d: double): boolean;   {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if d is a NaN or infinite}
begin
  IsNaNorInfD := THexDblW(d)[3] and $7FF0=$7FF0;
end;


{---------------------------------------------------------------------------}
function IsInfS(s: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is +INF or -INF}
var
  L: longint absolute s;
begin
  IsInfS := L and $7FFFFFFF = $7F800000;
end;


{---------------------------------------------------------------------------}
function IsNaNS(s: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if s is a NaN}
var
  L: longint absolute s;
begin
  IsNaNS := (L and $7F800000 = $7F800000) and  (L and $7FFFFF <> 0);
end;


{---------------------------------------------------------------------------}
function IsNaNorInfS(s: single): boolean; {$ifdef HAS_INLINE} inline;{$endif}
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
  small: double = 1e-30;
var
  x,y: double;
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
        hm := longint($80000000);
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
function rint(x: double): double;
  {-Return the integral value nearest x for the current rounding mode}
const
  twopower52h: THexDblW = ($0000,$0000,$0000,$4330); {2^63}
var
  c: double absolute twopower52h;
begin
  if THexDblW(x)[3] and $7FF0 > $4340 then begin
    {x is either INF or NAN or has no fractional part}
    rint := x;
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
 0.16     03.12.09  we          used ldexp for scalbn
***************************************************************************)

{---------------------------------------------------------------------------}
{------------------  Cody and Waite argument reduction ---------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function rem_pio2_cw(x: double; var z: double): integer;
  {-Cody/Waite reduction of x: z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8}
var
  i: longint;
  y: double;
const
  HP1 : THexDblW = ($0000,$4000,$21fb,$3fe9);
  HP2 : THexDblW = ($0000,$0000,$442d,$3e64);
  HP3 : THexDblW = ($5170,$98cc,$4698,$3ce8);
var
  DP1 : double absolute HP1; {= 7.85398125648498535156E-1}
  DP2 : double absolute HP2; {= 3.77489470793079817668E-8}
  DP3 : double absolute HP3; {= 2.69515142907905952645E-15}
begin
  {This is my Pascal translation of the CW reduction given in}
  {sin.c from the Cephes Math Library Release 2.8: June, 2000}
  {Copyright 1985, 1995, 2000 by Stephen L. Moshier          }

  if x=0.0 then begin
    z := 0.0;
    rem_pio2_cw := 0;
  end
  else begin
    y := floord(x/Pi_4);
    i := trunc(y - 16.0*floord(y/16.0));
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


{Table of constants for 2/pi, 396 Hex digits (476 decimal) of 2/pi }
const
  ipio2: array[0..65] of longint = (
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
           $4D7327, $310606, $1556CA, $73A8C9, $60E27B, $C08C6B);
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
  z  := ldexpd(z,q0);             {actual value of z}
  z  := z - 8.0*floord(z*0.125);  {trim off integer >= 8}
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
      if carry<>0 then  z := z - ldexpd(1.0,q0);
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
    z := ldexpd(z,-q0);
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
  fw := ldexpd(1.0,q0);
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
function rem_pio2_ph(x: double; var z: double): integer;
  {-Payne/Hanek reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8}
var
  ax, ay: TDA02;
  e0: integer;
begin
  z := abs(x);
  if (THexDblW(x)[3] and $7FF0=$7FF0) or (x=0.0) then rem_pio2_ph := 0
  else begin
    e0 := ilogb(z)-23;
    z  := ldexpd(z,-e0);
    ax[0] := trunc(z); z := ldexpd(z-ax[0],24);
    ax[1] := trunc(z);
    ax[2] := trunc(ldexpd(z-ax[1],24));
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
function rem_pio2(x: double; var z: double): integer;
  {-Argument reduction of x:  z = x - n*Pi/2, |z| <= Pi/4, result = n mod 8.}
  { Uses Payne/Hanek if |x| >= ph_cutoff, Cody/Waite otherwise}
const
  tol: double = 2.384185791015625E-7;  {ph_cutoff*eps_d}
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
      {If x is near a multiple of Pi/2, the C/W relative error may be large.}
      {In this case redo the calculation with the Payne/Hanek algorithm.    }
      rem_pio2 := rem_pio2_ph(x,z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function rem_2pi(x: double): double;
  {-Return x mod 2*Pi}
var
  {$ifndef ExtendedSyntax_on} n: integer;{$endif}
  z: double;
const
  {Reduction constants: 2*Pi = Pi2H + Pi2M + Pi2L}
  H2PH : THexDblW = ($0000,$4000,$21fb,$4019);
  H2PM : THexDblW = ($0000,$0000,$442d,$3e94);
  H2PL : THexDblW = ($5170,$98cc,$4698,$3d18);
var
  Pi2H: double absolute H2PH; {= 6.283185005187988}
  Pi2M: double absolute H2PM; {= 3.019915766344639E-7}
  Pi2L: double absolute H2PL; {= 2.156121143263248E-14}

begin
  if IsNanOrInfD(x) then begin
    rem_2pi := Nan_d;
    exit;
  end;
  if abs(x) <= ph_cutoff then begin
    {Direct Cody/Waite style reduction. This is more efficient}
    {than calling rem_pio2_cw with additional adjustment.}
    z := floord(x/TwoPi);
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
  {All z with -2eps_d <= z < 0 will have z + TwoPi = TwoPi}
  if z<0.0 then z := ((z + Pi2H) + Pi2M) + Pi2L;
  rem_2pi := z;
end;


{---------------------------------------------------------------------------}
function rem_2pi_sym(x: double): double;
  {-Return x mod 2*Pi, -Pi <= result <= Pi}
var
  {$ifndef ExtendedSyntax_on} n: integer;{$endif}
  z: double;
begin
  {$ifndef ExtendedSyntax_on} n:={$endif} rem_pio2(0.25*x,z);
  rem_2pi_sym := 4.0*z;
end;


{---------------------------------------------------------------------------}
function rem_int2(x: double; var z: double): integer;
  {-Argument reduction of x: z*Pi = x*Pi - n*Pi/2, |z|<=1/4, result = n mod 8.}
  { Used for argument reduction in sin(Pi*x) and cos(Pi*x)}
var
  y: double;
  i: integer;
  w: word;
begin
  w := THexDblW(x)[3] and $7FF0;
  if (w=$7FF0) or (abs(x)<=0.25) then begin
    {Nan, Inf, or <= 1/4}
    rem_int2 := 0;
    z := x;
    exit;
  end;
  if frac(x)=0.0 then begin
    {Here x is an integer or abs(x) >= 2^52}
    z := 0.0;
    i := 0;
    {set i=2, if x is a odd}
    if (w < $4340) and (frac(0.5*x)<>0.0) then i:=2;
  end
  else begin
    {Here x is not an integer. First calculate x mod 2,}
    {this leaves Pi*x = Pi*(x mod 2) mod 2*Pi invariant}
    x := 0.5*x;
    x := 2.0*(x-floord(x));
    {then apply the Cody/Waite style range reduction}
    y := floord(4.0*x);
    i := trunc(y - 16.0*floord(y/16.0));
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
function sum2(const a: array of double; n: integer): double;
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
function dot2(const x,y: array of double; n: integer): double;
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
function rms(const a: array of double; n: integer): double;
  {-Calculate the RMS value sqrt(sum(a[i]^2, i=0..n-1)/n) of a double vector}
var
  scale, sumsq: double;
begin
  ssqd(a,n,scale,sumsq);
  rms := scale*sqrt(sumsq/n);
end;


{---------------------------------------------------------------------------}
procedure ssqd(const a: array of double; n: integer; var scale, sumsq: double);
  {-Calculate sum(a[i]^2, i=0..n-1) = scale^2*sumsq, scale>=0, sumsq>0}
var
  i: integer;
  t: double;
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
function sumsqr(const a: array of double; n: integer): double;
  {-Calculate sum(a[i]^2, i=0..n-1) of a double vector}
var
  scale, sumsq: double;
begin
  ssqd(a,n,scale,sumsq);
  sumsqr := sqr(scale)*sumsq;
end;


{---------------------------------------------------------------------------}
function norm2(const a: array of double; n: integer): double;
  {-Calculate the 2-norm = sqrt(sum(a[i]^2, i=0..n-1)) of a double vector}
var
  scale, sumsq: double;
begin
  ssqd(a,n,scale,sumsq);
  norm2 := scale*sqrt(sumsq);
end;


{---------------------------------------------------------------------------}
function PolEval(x: double; const a: array of double; n: integer): double;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}
var
  i: integer;
  s: double;
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
  uj,vj,r,s,t: double;
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
procedure PolEvalDeriv(x: double; const a: array of double; n: integer; var px,dp: double);
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
function PolEvalS(x: double; const a: array of single; n: integer): double;
  {-Evaluate polynomial; return a[0] + a[1]*x + ... + a[n-1]*x^(n-1)}
var
  i: integer;
  s: double;
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
function CSEvalD(x: double; const a: array of double; n: integer): double;
  {-Evaluate Chebyshev sum a[0]/2 + a[1]*T_1(x) +..+ a[n-1]*T_(n-1)(x) using Clenshaw algorithm}
var
  b0,b1,b2: double;
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
  CSEvalD := 0.5*(b0-b2);
end;


{---------------------------------------------------------------------------}
procedure CSEvalDDer(x: double; const a: array of double; n: integer; var px,dp: double);
  {-Evaluate Chebyshev sum p(x) = a[0]/2 + a[1]*T_1(x) +..+ a[n-1]*T_(n-1)(x),}
  { using Clenshaw algorithm, return pd = p(x) and dp = p'(x) }
var
  x2,b0,b1,b2,c0,c1,c2: double;
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
function mean(const a: array of double; n: integer): double;
  {-Compute accurate mean = sum(a[i], i=0..n-1)/n of a double vector}
begin
  if n<=0 then mean := 0.0
  else mean := sum2(a,n)/n;
end;


{---------------------------------------------------------------------------}
procedure mssqd(const a: array of double; n: integer; var mval, scale, sumsq: double);
  {-Calculate mean mval and ssqd sum((a[i]-mval)^2) of a double vector}
var
  i: integer;
  t: double;
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
procedure MeanAndStdDev(const a: array of double; n: integer; var mval, sdev: double);
  {-Accurate mean and sample standard deviation of a double vector}
var
  scale, sumsq: double;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then sdev := 0.0
  else sdev := scale*sqrt(sumsq/(n-1));
end;


{---------------------------------------------------------------------------}
function CoeffVar(const a: array of double; n: integer): double;
  {-Return the coefficient of variation (sdev/mean) of a double vector}
var
  mval, sdev: double;
begin
  MeanAndStdDev(a,n,mval,sdev);
  CoeffVar := SafeDiv(sdev,mval);
end;


{---------------------------------------------------------------------------}
procedure moment(const a: array of double; n: integer; var m1, m2, m3, m4, skew, kurt: double);
  {-Return the first 4 moments, skewness, and kurtosis of a double vector}
var
  q,d,scale,sumsq: double;
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
function ssdev(const a: array of double; n: integer): double;
  {-Return the sample standard deviation of a double vector}
var
  mval,scale,sumsq: double;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then ssdev := 0.0
  else ssdev := scale*sqrt(sumsq/(n-1));
end;


{---------------------------------------------------------------------------}
function psdev(const a: array of double; n: integer): double;
  {-Return the population standard deviation of a double vector}
var
  mval,scale,sumsq: double;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then psdev := 0.0
  else psdev := scale*sqrt(sumsq/n);
end;


{---------------------------------------------------------------------------}
function svar(const a: array of double; n: integer): double;
  {-Return the sample variance of a double vector}
var
  mval,scale,sumsq: double;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then svar := 0.0
  else svar := sqr(scale)*(sumsq/(n-1));
end;


{---------------------------------------------------------------------------}
function pvar(const a: array of double; n: integer): double;
  {-Return the population variance of a double vector}
var
  mval,scale,sumsq: double;
begin
  mssqd(a, n, mval, scale, sumsq);
  if n<2 then pvar := 0.0
  else pvar := sqr(scale)*(sumsq/n);
end;


{---------------------------------------------------------------------------}
function pcov(const x,y: array of double; n: integer): double;
  {-Calculate the population covariance = sum((x[i]-mx)(y[i]-my), i=0..n-1)/n of two double vectors}
var
  mx,my,s: double;
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
function tvar(const a: array of double; n: integer): double;
  {-Return the total variance of a double vector}
var
  mval,scale,sumsq: double;
begin
  mssqd(a, n, mval, scale, sumsq);
  tvar := sqr(scale)*(sumsq);
end;


{--------------------------------------------------------------------}
{------------------------- Other function ---------------------------}
{--------------------------------------------------------------------}

{---------------------------------------------------------------------------}
function DegToRad(x: double): double;
  {-Convert angle x from degrees to radians}
begin
  DegToRad := x*Pi_180;
end;


{---------------------------------------------------------------------------}
function RadToDeg(x: double): double;
  {-Convert angle x from radians to degrees}
begin
  RadToDeg := x/Pi_180;
end;


{---------------------------------------------------------------------------}
function angle2(x1,x2,y1,y2: double): double;
  {-Return the accurate angle between the vectors (x1,x2) and (y1,y2)}
var
  d1,d2: double;
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
function orient2d(x1,y1,x2,y2,x3,y3: double): double;
  {-Return the mathematical orientation of the three points (xi,yi): >0 if}
  { the order is counterclockwise, <0 if clockwise, =0 if they are collinear.}
  { Result is twice the signed area of the triangle defined by the points.}
var
  det,detl,detr: double;
const
  eps = 6.6613381478E-16;  {= (3 + 16*eps_d)*eps_d}
var
  v1,v2: array[0..5] of double;
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
  {but here I calculate orient2 as an accurate dot product using dot2. }

  v1[0] := x1;    v2[0] :=  y2;
  v1[1] := x1;    v2[1] := -y3;
  v1[2] := y1;    v2[2] := -x2;
  v1[3] := y1;    v2[3] :=  x3;
  v1[4] := x3;    v2[4] := -y2;
  v1[5] := y3;    v2[5] :=  x2;

  orient2d := dot2(v1,v2,6);
end;


{---------------------------------------------------------------------------}
function area_triangle(x1,y1,x2,y2,x3,y3: double): double;
  {-Return the area of the triangle defined by the points (xi,yi)}
begin
  area_triangle := 0.5*abs(orient2d(x1,y1,x2,y2,x3,y3));
end;


{---------------------------------------------------------------------------}
function in_triangle(x,y,x1,y1,x2,y2,x3,y3: double): boolean;
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
function in_triangle_ex(x,y,x1,y1,x2,y2,x3,y3: double): integer;
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


{$ifndef BIT16}
  type char8=ansichar;
{$else}
  type char8=char;
{$endif}

const
  hnib: array[0..15] of char8 = '0123456789ABCDEF';


{---------------------------------------------------------------------------}
function Dbl2Hex(d: double): string;
  {-Return d as a big-endian hex string}
var
  ad : THexDblA absolute d;
  i,j: integer;
  s  : string[16];
begin
  {$ifndef BIT16}
    {$ifdef FPC32Plus}
      s := '';  {Avoid silly warning}
    {$endif}
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
    {$ifdef FPC32Plus}
      t := '';  {Avoid silly warning}
    {$endif}
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
procedure Hex2Float(const hex: string; n: integer; var a: THexDblA; var code: integer);
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
  else if (n=7) or (n=3) then begin
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
procedure Hex2Dbl(const hex: string; var d: double; var code: integer);
  {-Convert big-endian hex string to double, leading $ is skipped, OK if code=0;}
  { hex must have 16 hex characters (17 if leading $), inverse of Dbl2Hex.}
var
  a: THexDblA;
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
  a: THexDblA;
  t: single absolute a;
begin
  Hex2Float(hex, 3, a, code);
  if code=0 then s := t;
end;


{---------------------------------------------------------------------------}
function fisEQd(x,y: double): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are bit-identical}
begin
  {FPC64 2.6.x with -O2/O3 will generate internal error for absolute xr.hm=yr.hm etc}
  fisEQd := (TDblRec(x).hm=TDblRec(y).hm) and (TDblRec(x).lm=TDblRec(y).lm);
end;


{---------------------------------------------------------------------------}
function fisNEd(x,y: double): boolean; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return true if x and y are not bit-identical}
begin
  {FPC64 2.6.x with -O2/O3  will generate internal error for absolute xr.hm<>yr.hm etc}
  fisNEd := (TDblRec(x).hm<>TDblRec(y).hm) or (TDblRec(x).lm<>TDblRec(y).lm);
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
function isign(x: double): integer;
  {-Return the sign of x, 0 if x=0 or NAN}
begin
  if IsNaNd(x) or (x=0.0) then isign := 0
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
function mins(x, y: single): single; {$ifdef HAS_INLINE} inline;{$endif}
  {-Return the minimum of two singles; x,y <> NAN}
begin
  if x<y then mins := x
  else mins := y;
end;


{---------------------------------------------------------------------------}
function RandG01: double;
  {-Random number from standard normal distribution (Mean=0, StdDev=1)}
var
  s,v1,v2: double;
{$ifdef J_OPT}
  {$J+}
{$endif}
const
  RGNxt01: double  = 0.0;
  RGValid: boolean = false;
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
function RandG(Mean, StdDev: double): double;
  {-Random number from Gaussian (normal) distribution with given mean}
  { and standard deviation |StdDev|}
begin
  RandG := Mean + StdDev*RandG01;
end;


{---------------------------------------------------------------------------}
function SafeDiv(x,y: double): double;
  {-Safe quotient x/y, +-Infinity if overflow, NaN if x=y=0}
var
  a: double;
begin
  a := abs(y);
  if a >= 1.0 then SafeDiv := x/y
  else begin
    if (a=0.0) and (x=0.0) then SafeDiv := Nan_d
    else begin
      if abs(x) > a*Maxdouble then begin
        if (THexDblW(x)[3] xor THexDblW(y)[3]) and $8000 = 0 then begin
          SafeDiv := PosInf_d;
        end
        else begin
          SafeDiv := NegInf_d;
        end;
      end
      else SafeDiv := Nan_d;
    end;
  end;
end;



{--------------------------------------------------------------------}
{-----------------  dbl2 (double-double) functions  -----------------}
{--------------------------------------------------------------------}

{$ifdef FPC2Plus}
  {$ifndef DEBUG}
    {$warn 5036 OFF}  {Local variable "xx" does not seem to be initialized}
  {$endif}
{$endif}

const
  CSD =  134217729.0; {2^27+1}  {Split constant for double}

{* IMPORTANT NOTE: The splitting of statements compared to the ext2 AMath }
{* routines is (at least) necessary for Delphi 64-bit compilers, they seem}
{* to optimize away some of the features, see DAMTools/dscrmt for similar }
{* observation. FPC (even with optimization) does not have these problems.}


{---------------------------------------------------------------------------}
procedure ddto2d(a,b: double; var x: dbl2); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = a + b using TwoSum algorithm}
var
  z: double;
begin
  {This the TwoSum algorithm, see Knuth [32], Sect. 4.2.2, Theorem A/B}
  {or Ogita et al. [8], Algorithm 3.1, or Graillat [46], Algorithm 2.1}
  {FastTwosum from older DAMATH versions is correct only if |a| >= |b|}
  x.h := a + b;
  z   := x.h - a;
  x.l := (a - (x.h - z)) + (b - z);
end;


{---------------------------------------------------------------------------}
procedure hhto2d(const a,b: THexDblW; var x: dbl2);
  {-Return x = a + b using xxto2x}
begin
  ddto2d(double(a), double(b), x);
end;


{---------------------------------------------------------------------------}
procedure dadd12(const a, b: double; var xh, xl: double); {$ifdef HAS_INLINE} inline;{$endif}
  {-Return [xh,xl] = a + b}
var
  z: double;
begin
  xh := a + b;
  z  := xh - a;
  xl := (a - (xh - z)) + (b - z);
end;


{---------------------------------------------------------------------------}
procedure dto2d(a: double; var x: dbl2);
  {-Return x = a}
begin
  x.h := a;
  x.l := 0.0;
end;


{---------------------------------------------------------------------------}
procedure dmul12(a,b: double; var xh,xl: double);
  {-Return x = a * b}
var
  t,a1,a2,b1,b2: double;
begin
  {G.W. Veltkamp's multiplication routine, see Dekker[64] mul12}
  {See also Linnainmaa[65]: exactmul2; and [8],[46]: TwoProduct}
  t  := a*CSD;
  a1 := a - t;
  a1 := a1 + t;
  a2 := a - a1;
  t  := b*CSD;
  b1 := b - t;
  b1 := b1 + t;
  b2 := b - b1;
  xh := a*b;
  xl := a1*b1;
  xl := xl - xh;
  t  := a1*b2;
  xl := xl + t;
  t  := a2*b1;
  xl := xl + t;
  t  := a2*b2;
  xl := xl + t;
end;


{---------------------------------------------------------------------------}
procedure sqr12d(a: double; var xh,xl: double);
  {-Return [xh,xl] = a^2}
var
  t,a1,a2: double;
begin
  {mul12 with b=a, avoids second split}
  t  := a*CSD;
  a1 := a - t;
  a1 := a1 + t;
  a2 := a - a1;
  xh := a*a;
  xl := a1*a1;
  xl := xl - xh;
  t  := a1*a2;
  xl := xl + t;
  xl := xl + t;
  t  := a2*a2;
  xl := xl + t;
end;


{---------------------------------------------------------------------------}
procedure mul2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a*b}
var
  u,v,zh,zl: double;
begin
  {Linnainmaa[65]: longmul}
  dmul12(a.h, b.h, zh, zl);
  u := a.h + a.l;
  u := u*b.l;
  v := a.l*b.h;
  u := u+v;
  zl  := u + zl;
  x.h := zh + zl;
  x.l := zh - x.h;
  x.l := x.l + zl;
end;


{---------------------------------------------------------------------------}
procedure mul21d(const a: dbl2; b: double; var x: dbl2);
  {-Return x = a*b}
var
  u,zh,zl: double;
begin
  {mul2 with b.h=b, b.l=0}
  dmul12(a.h, b, zh, zl);
  u   := a.l*b;
  zl  := u + zl;
  x.h := zh + zl;
  u   := zh - x.h;
  x.l := u + zl;
end;


{---------------------------------------------------------------------------}
procedure sqr2d(const a: dbl2; var x: dbl2);
  {-Return x = a^2}
var
  u,v,zh,zl: double;
begin
  {Simplified mul2}
  sqr12d(a.h, zh, zl);
  u := (a.h + a.l);
  u := u*a.l;
  v := a.l*a.h;
  u := u + v;
  zl  := u + zl;
  x.h := zh + zl;
  u   := zh - x.h;
  x.l := u + zl;
end;


{---------------------------------------------------------------------------}
procedure div2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a/b,  b<>0}
var
  u,qh,ql,zh,zl: double;
begin
  {Linnainmaa[65]: longdiv}
  zh := a.h/b.h;
  dmul12(b.h, zh, qh, ql);
  zl := a.h - qh;
  zl := zl - ql;
  zl := zl + a.l;
  u  := zh*b.l;
  zl := zl - u;
  u  := b.h+b.l;
  zl := zl/u;
  x.h := zh + zl;
  x.l := zh - x.h;
  x.l := x.l + zl;
end;


{---------------------------------------------------------------------------}
procedure div21d(const a: dbl2; b: double; var x: dbl2);
  {-Return x = a/b,  b<>0}
var
  qh,ql,zh,zl: double;
begin
  {div2 with b.l=0}
  zh := a.h/b;
  dmul12(b, zh, qh, ql);
  zl := a.h - qh;
  zl := zl - ql;
  zl := zl + a.l;
  zl := zl / b;
  x.h := zh + zl;
  x.l := zh - x.h;
  x.l := x.l + zl;
end;


{---------------------------------------------------------------------------}
procedure ddivrem(a,b: double; var q,r: double);
  {-Compute q,r with a = q*b + r,  assumes round to nearest}
var
  x,y: double;
begin
  q := a/b;
  dmul12(q,b,x,y);
  r := (a-x)-y;
end;


{---------------------------------------------------------------------------}
procedure inv2d(const b: dbl2; var x: dbl2);
  {-Return x = 1/b,  b<>0}
var
  u,qh,ql,zh,zl: double;
begin
  {div2 with a,h=1, a.l=0}
  zh := 1.0/b.h;
  dmul12(b.h, zh, qh, ql);
  zl := 1.0 - qh;
  zl := zl - ql;
  u  := zh*b.l;
  zl := zl - u;
  u  := b.h+b.l;
  zl := zl/u;
  x.h := zh + zl;
  x.l := zh - x.h;
  x.l := x.l + zl;
end;


{---------------------------------------------------------------------------}
procedure add2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a+b}
var
  u,zh,zl: double;
begin
  {Linnainmaa[65]: longadd}
  zh := a.h + b.h;
  zl := a.h - zh;
  u  := zl + zh;
  u  := a.h - u;
  u  := u + a.l;
  zl := zl + b.h;
  zl := zl + u;
  zl := zl + b.l;
  x.h := zh + zl;
  x.l := zh - x.h;
  x.l := x.l + zl;
end;


{---------------------------------------------------------------------------}
procedure add21d(const a: dbl2; b: double; var x: dbl2);
  {-Return x = a+b}
var
  u, zh,zl: double;
begin
  {add2 with b.l=0}
  zh := a.h + b;
  zl := a.h - zh;
  u  := zl + zh;
  u  := a.h - u;
  u  := u + a.l;
  zl := zl + b;
  zl := zl + u;
  x.h := zh + zl;
  x.l := zh - x.h;
  x.l := x.l + zl;
end;


{---------------------------------------------------------------------------}
procedure sub2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a-b}
var
  u,zh,zl: double;
begin
  {add2 with a + (-b)}
  zh := a.h - b.h;
  zl := a.h - zh;
  u  := zl + zh;
  u  := a.h - u;
  u  := u + a.l;
  zl := zl - b.h;
  zl := zl + u;
  zl := zl - b.l;
  x.h := zh + zl;
  x.l := zh - x.h;
  x.l := x.l + zl;
end;


{---------------------------------------------------------------------------}
procedure sqrt2d(const a: dbl2; var x: dbl2);
  {-Return x = sqrt(a),  a >= 0}
var
  qh,ql,zh,zl: double;
begin
  {Dekker[64] sqrt2, with mul12 replaced by sqr12}
  zh := sqrt(a.h);
  sqr12d(zh, qh, ql);
  zl := a.h - qh;
  zl := zl - ql;
  zl := zl + a.l;
  zl := 0.5*zl/zh;
  x.h := zh + zl;
  x.l := zh - x.h;
  x.l := x.l + zl;
end;


{---------------------------------------------------------------------------}
procedure pow2di(const a: dbl2; n: double; var y: dbl2);
  {-Return y = a^n, frac(n)=0, a<>0 if n<0}
var
  x: double;
  z: dbl2;
begin
  {$ifdef debug}
    if frac(n) <> 0.0 then begin
      writeln('pow2di:  frac(n) <>0');
    end;
  {$endif}
  x := int(n);
  if x<0.0 then begin
    {1/a^x instead of (1/a)^x maybe a bit more accurate but can overflow}
    inv2d(a,y);
    pow2di(y,-x,y);
  end
  else if x<=2.0 then begin
    if x=1.0 then y := a
    else if x=2.0 then sqr2d(a,y)
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
      if frac(x)<>0.0 then mul2d(y,z,y);
      x := int(x);
      if x=0.0 then exit
      else sqr2d(z,z);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure pi2d(var dpi: dbl2);
  {-Return x = Pi with double-double precision}
const
  pih: THexDblW = ($2D18,$5444,$21FB,$4009); {3.14159265358979E+0000}
  pil: THexDblW = ($5C07,$3314,$A626,$3CA1); {1.22464679914735E-0016}
begin
  dpi.l := double(pil);
  dpi.h := double(pih);
end;


{---------------------------------------------------------------------------}
function fma_d(a,b,c: double): double;
  {-Accurately compute a*b + c}
var
  x,y,z: double;
  ah,al,bh,bl: double;
begin

  {(x,y) = TwoProdct(a,b) [46] Algorithm 5}

  {Split a, [46] Algorithm 4}
  x  := a * b;
  z  := a*csd;
  ah := z - a;
  ah := z - ah;
  al := a - ah;

  {Split b}
  z  := b*csd;
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
  {(ah, al) = TwoSum(x,y) from [46] Algorithm 2}
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

  fma_d := x + y;

end;


{---------------------------------------------------------------------------}
procedure floor2d(const a: dbl2; var x: dbl2);
  {-Return x = floor(a)}
var
  t: double;
begin
  t := floord(a.h);
  if t<>a.h then begin
    x.h := t;
    x.l := 0.0;
  end
  else ddto2d(t,floord(a.l),x);
end;


{---------------------------------------------------------------------------}
procedure ceil2d(const a: dbl2; var x: dbl2);
  {-Return x = ceil(a)}
var
  t: double;
begin
  t := ceild(a.h);
  if t<>a.h then begin
    x.h := t;
    x.l := 0.0;
  end
  else ddto2d(t,ceild(a.l),x);
end;


{---------------------------------------------------------------------------}
procedure trunc2d(const a: dbl2; var x: dbl2);
  {-Return x = trunc(a)}
begin
  if a.h >= 0.0 then floor2d(a,x)
  else ceil2d(a,x);
end;


{---------------------------------------------------------------------------}
procedure ldexp2d(const a: dbl2; n: longint; var x: dbl2);
  {-Return x = a*2^n}
begin
  x.h := ldexpd(a.h, n);
  x.l := ldexpd(a.l, n);
end;


{---------------------------------------------------------------------------}
procedure abs2d(const a: dbl2; var x: dbl2);  {$ifdef HAS_INLINE} inline;{$endif}
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
procedure chs2d(const a: dbl2; var x: dbl2);   {$ifdef HAS_INLINE} inline;{$endif}
  {-Return x = -a}
begin
  x.h := -a.h;
  x.l := -a.l;
end;


{---------------------------------------------------------------------------}
function cmp2d(const a,b: dbl2): integer;
  {-Return sign(a-b)}
var
  t: dbl2;
begin
  sub2d(a,b,t);
  cmp2d := isign(t.h + t.l);
end;


{---------------------------------------------------------------------------}
procedure ln2_2d(var x: dbl2);
  {-Return x = ln(2) with double-double precision}
const
  ln2h: THexDblW = ($39EF,$FEFA,$2E42,$3FE6);  { 6.93147180559945E-0001}
  ln2l: THexDblW = ($803F,$3B39,$BC9E,$3C7A);  { 2.31904681384630E-0017}

begin
  x.l := double(ln2l);
  x.h := double(ln2h);
end;


{---------------------------------------------------------------------------}
procedure exp1_2d(var x: dbl2);
  {-Return x = exp(1) with double-double precision}
const
  e2h: THexDblW = ($5769,$8B14,$BF0A,$4005);  { 2.71828182845905E+0000}
  e2l: THexDblW = ($013A,$E2B1,$D57E,$3CA4);  { 1.44564689172925E-0016}
begin
  x.l := double(e2l);
  x.h := double(e2h);
end;


{---------------------------------------------------------------------------}
procedure exp2d(const a: dbl2; var x: dbl2);
  {-Return x = exp(a)}
var
  m: longint;
  r,s,t,p: dbl2;
  f: double;
  i: integer;
const
  IMAX  = 22;     {Max iteration index, IMAX! must be exact}
  eps   = 1e-34;  {eps_d^2/512}
  shift = 10;
  FAC   = 1 shl shift;
begin

  if a.h < ln_MinDbl then begin
    x.h := 0.0;
    x.l := 0.0;
    exit;
  end
  else if a.h > ln_MaxDbl then begin
    x.h := PosInf_d;
    x.l := 0.0;
    exit;
  end;

  {Range reduction}
  m := round(a.h/ln2);
  ln2_2d(r);
  mul21d(r,m,r);
  sub2d(a,r,r);
  r.h := r.h/FAC;
  r.l := r.l/FAC;

  {Taylor series}
  sqr2d(r,p);
  s.h := 0.5*p.h;
  s.l := 0.5*p.l;
  add2d(r,s,s);
  mul2d(p,r,p);
  f := 6.0;
  div21d(p,f,t);

  for i:=4 to IMAX do begin
    add2d(s,t,s);
    mul2d(p,r,p);
    f := f*i;
    div21d(p,f,t);
    if abs(t.h) < eps then break;
    {$ifdef debug}
      if i=IMAX then writeln('exp2d: i=IMAX');
    {$endif}
  end;
  add2d(s,t,s);
  {undo range reduction}
  for i:=1 to shift do begin
    sqr2d(s,t);
    s.h := 2*s.h;
    s.l := 2*s.l;
    add2d(s,t,s);
  end;
  add21d(s,1,s);
  ldexp2d(s, m, x);
end;


{---------------------------------------------------------------------------}
procedure ln2d(const a: dbl2; var x: dbl2);
  {-Return x = ln(a)}
var
  y,z: dbl2;
begin
  {initial approximation y = -ln(x)}
  y.h := -ln(a.h);
  y.l := 0;
  {One Newton step}
  exp2d(y,z);
  mul2d(z,a,z);
  add21d(z,-1.0,z);
  sub2d(z,y,x);
end;


{---------------------------------------------------------------------------}
procedure pow2d(const a,b: dbl2; var x: dbl2);
  {-Return x = a^b, a > 0}
var
  t: dbl2;
begin
  {x = exp(b*ln(a))}
  ln2d(a,t);
  mul2d(t,b,t);
  exp2d(t,x)
end;


{---------------------------------------------------------------------------}
procedure nroot2d(const a: dbl2; n: longint; var x: dbl2);
  {-Return x = a^(1/n), a > 0 if n is even}
var
  k: longint;
  z: double;
begin
  k := abs(n);
  z := a.h;
  if (k and 1 = 0) and (z<0.0) then begin
    x.h := nroot(z,n);
    x.l := 0;
    exit;
  end
  else abs2d(a,x);
  if k=2 then sqrt2d(x,x)
  else begin
    ln2d(x,x);
    div21d(x,k,x);
    exp2d(x,x);
  end;
  if z < 0 then chs2d(x,x);
  if k < 0 then inv2d(x,x);
end;


end.

