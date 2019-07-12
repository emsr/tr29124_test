unit mp_real;

{Multi precision real floating point arithmetic routines}

interface


{$ifdef VirtualPascal}
  {$X+} {needed for pchars/RESULT}
{$endif}

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
  {$X+} {needed for pchars}
{$endif}

uses
  BTypes, mp_types;

{$i mp_conf.inc}

(*************************************************************************

 DESCRIPTION   :  Multi precision real floating point arithmetic routines

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25S, FPC, VP, WDOSX

 EXTERNAL DATA :  (mp_types)

 MEMORY USAGE  :  heap

 DISPLAY MODE  :  ---

 REFERENCES    :  [3] D.E. Knuth: The Art of computer programming, Volume 2,
                      Seminumerical Algorithms, 3rd ed., 1998
                      http://www-cs-faculty.stanford.edu/~knuth/taocp.html
                  [8] Marcel Martin: NX - Numerics library of multiprecision
                      numbers for Delphi and Free Pascal, 2006-2009
                      www.ellipsa.eu/public/nx/index.html
                 [22] LiDIA - A Library for Computational Number Theory, 2006,
                      LiDIA-Group, Technische Universit„t Darmstadt, Fachbereich Informatik
                 [27] T. Papanikolaou: libF - Eine lange Gleitpunktarithmetik,
                      Diplomarbeit 1995, available online from
                      http://www.cdc.informatik.tu-darmstadt.de/reports/reports/papa.diplom.ps.gz
                 [31] J. Arndt, Matters Computational. Ideas, algorithms, source code.
                      http://www.jjj.de/fxt/#fxtbook
                 [35] R.P. Brent, P. Zimmermann: Modern Computer Arithmetic, Cambridge University Press, 2010.
                      A preliminary version (V0.5.9, Oct. 2010) of the book is available from
                      http://maths-people.anu.edu.au/~brent/pd/mca-cup-0.5.9.pdf
                      or http://arxiv.org/abs/1004.4710 (V0.5.1)
                 [37] R.M. Corless, G.H. Gonnet, D.E.G. Hare, D.J. Jeffrey, D.E. Knuth,
                      On the Lambert W Function, Adv. Comput. Math., 5 (1996), pp. 329-359.
                      http://www.apmaths.uwo.ca/~rcorless/frames/PAPERS/LambertW/LambertW.ps or
                      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.33.3583
                 [38] F. Johansson, Efficient implementation of the Hardy-Ramanujan-Rademacher
                      formula, 2012. Available from http://arxiv.org/abs/1205.5991v1
                 [41] [HMF]: M. Abramowitz, I.A. Stegun. Handbook of Mathematical
                      Functions. New York, 1970, http://www.math.sfu.ca/~cbm/aands/


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------

 0.0.01   01.11.07  W.Ehrhardt  Basic type definitions, mpf_init etc
 0.0.02   03.11.07  we          mpf_todouble/extended, s_mpf_normalize, mpf_set_int
 0.0.03   03.11.07  we          s_mpf_normalizep, mpf_initp, mpf_mul, mpf_expt_int
 0.0.04   03.11.07  we          mpf_exch, mpf_inv, mpf_div, mpf_is0
 0.0.05   04.11.07  we          mpf_abs, mpf_chs, mpf_mul_2k, s_mpf_incexp
 0.0.06   04.11.07  we          mpf_add, mpf_sub, s_mpf_addsub
 0.0.07   04.11.07  we          removed .flags, fix mpf_abs, mpf_chs
 0.0.08   04.11.07  we          mpf_cmp_mag, mpf_initp2
 0.0.09   05.11.07  we          mpf_get_default_prec, mpf_sqrt
 0.0.10   07.11.07  we          mpf_cmp, mpf_noinit, arg checks
 0.0.11   07.11.07  we          s_mpf_normalize: stick bit and postnormalization
 0.0.12   07.11.07  we          mpf_set_ext
 0.0.13   08.11.07  we          mpf_chg_prec, mpf_checksum
 0.0.14   08.11.07  we          mpf_int, mpf_frac, mpf_trunc
 0.0.15   09.11.07  we          mpf_add/sub for a=0, skip final sqr in mpf_expt_int
 0.0.16   10.11.07  we          mpf_is1, mpf_is1a, special cases in mpf_mul/div
 0.0.17   10.11.07  we          bugfix mpf_div
 0.0.18   10.11.07  we          mpf_div_mpi, mpf_mul_mpi, mpf_set_mpi
 0.0.19   10.11.07  we          s_mpf_addsub changed to s_mpf_inc, mpf_add/sub_mpi
 0.0.20   11.11.07  we          mpf_add_ext, mpf_sub_ext, mpf_div_ext, mpf_mul_ext
 0.0.21   11.11.07  we          mpf_read_decimal, mpf_read_radix
 0.0.22   12.11.07  we          Fix memory leak(s) if MPC_HaltOnError is not defined
 0.0.23   13.11.07  we          mpf_is_ge, mpf_is_lt
 0.0.24   14.11.07  we          s_mpf_to10_n, mpf_decimal, mpf_adecimal
 0.0.25   15.11.07  we          mpf_initp_multi_p, mpf_initp[x], mpf_clear[x],
 0.0.26   15.11.07  we          mpf_set1, overflow check in mpf_sqr
 0.0.27   16.11.07  we          s_mpf_inc: allow @x=@y, s_mpf_incf
 0.0.28   16.11.07  we          mpf_exp, bugfix mpf_read_radix for '-.xxx'
 0.0.29   16.11.07  we          s_mpf_incexp with BASM
 0.0.30   17.11.07  we          mpf_ln, use b with working precision in mpf_exp
 0.0.31   17.11.07  we          mpf_expt, mpf_iexpt
 0.0.32   18.11.07  we          exponent overflow checks via s_mpf_incexp
 0.0.33   19.11.07  we          mpf_pi_chud, mpf_pi_machin. pi via mp_read_unsigned_bin,
 0.0.34   20.11.07  we          mpf_set_pip uses AddrPiBytes (pi table linked only if needed)
 0.0.35   23.11.07  we          mpf_set_pip2k, mpf_set_pi2k, mpf_arctan
 0.0.36   25.11.07  we          mpf_trig, mpf_cos, mpf_sin, mpf_tan
 0.0.37   26.11.07  we          s_mpf_inc1
 0.0.38   26.11.07  we          rounding in mpf_set_pip2k, mpf_pi_chud/mpf_pi_machin outsourced
 0.0.39   26.11.07  we          mpf_is_eq, mpf_is_gt, mpf_is_le, mpf_is_ne
 0.0.40   26.11.07  we          mpf_random, s_mpf_to10_n (case a=1, removed c10)
 0.0.41   29.11.07  we          s_mpf_abs, s_mpf_chs, mpf_mul_d/int, mpf_div_d/int
 0.0.42   01.12.07  we          s_mpf_is0
 0.0.43   02.12.07  we          s_mpf_is_neg, bugfix mpf_div
 0.0.44   02.12.07  we          s_mpf_mod_pi2k (used in mpf_trig)
 0.0.45   09.12.07  we          mp_float to mp_types, mpf_ln with mpf_div_d/int
 0.0.46   09.12.07  we          use mpf_copyp in mp_trig to set ps^ and pc^
 0.0.47   09.12.07  we          improve mpf_sub_ext for a=0
 0.0.48   16.12.07  we          mpf_todecimal_n, mpf_is_eq_rel
 0.0.49   17.12.07  we          mpf_tohex_n, mpf_read_hex, (s_)mpf_toradix_n
 0.0.50   18.12.07  we          used s_mpf_normalizep in some routines for mpf_chg_prec
 0.0.51   18.12.07  we          (mpf_is_eq_rel_n)
 0.0.52   20.12.07  we          mpf_reldev
 0.0.53   21.12.07  we          mpf_trig: bugfix -eps, sum series for cos - 1
 0.0.54   21.12.07  we          improved mpf_exp

 1.4.00   01.01.08  we          check a<>0 in mpf_ln, s_mpf_add1
 1.4.01   01.01.08  we          s_mpf_dec1, mpf_expm1, mpf_ln1p
 1.4.02   02.01.08  we          mpf_cosh, mpf_sinh, mpf_tanh, mpf_atanh
 1.4.03   03.01.08  we          mpf_acosh, mpf_acosh1p, mpf_asinh
 1.4.04   03.01.08  we          mpf_arctan2, mpf_arccos, mpf_arcsin
 1.4.05   05.01.08  we          s_mpf_toradix_n: fix if rounded res = radix^ndd
 1.4.06   06.01.08  we          s_mpf_agm, mpf_agm, mpf_ccell1, mpf_ccell12
 1.4.07   10.01.08  we          mpf_set_mpi2k, mpf_ccell12: increased working precision

 1.5.00   16.01.08  we          Fix rounding/overflow in s_mpf_toradix_n
 1.5.01   16.01.08  we          mpf_writeln, mpf_write_decimal/radix, mpf_output_decimal/radix
 1.5.02   17.01.08  we          s_mpf_ldx, s_mpf_is_le0, s_mpf_is_ge0
 1.5.03   18.01.08  we          mpf_cosh/sinh: don't calculate inverse if it is to small
 1.5.04   19.01.08  we          mpf_acosh/asinh = ln(2a) for large a, mpf_sinh small fix
 1.5.05   20.01.08  we          mpf_ccell2, s_mpf_ccell12, s_mpf_incf with add32_ovr
 1.5.06   25.01.08  we          _set_ptab_const, mpf_set_ln2/10 functions
 1.5.07   25.01.08  we          mpf_10expt, mpf_2expt, mpf_log10, mpf_log2
 1.5.08   28.01.08  we          s_mpf_ldx returns MaxLongint if overflow, fix mpf_ln
 1.5.09   29.01.08  we          abs(b)=1 in mpf_mul/div_int/d
 1.5.10   30.01.08  we          mpf_read_radix: removed readxp quick hack
 1.5.11   31.01.08  we          {$x+} for VP and D1
 1.5.12   03.02.08  we          bugfix s_mpf_agm (don't use c in loop)

 1.6.00   11.06.08  we          MPXRange.Create('mpf_chg_prec: newprec out of range');

 1.7.00   23.08.08  we          Avoid FPC222 warning in mpf_decimal
 1.7.01   16.09.08  we          Accept 'b' or 'B' in binary exponents
 1.7.02   24.09.08  we          string replaced by mp_string

 1.9.00   02.12.08  we          Uses BTypes: char8, pchar8
 1.9.01   13.12.08  we          mpf_sumalt

 1.10.00  21.01.09  we          changes related to (s)mp_divrem
 1.10.01  21.02.09  we          mpf_sinhcosh
 1.10.02  23.02.09  we          mpf_read_radix: changed evaluation of fraction part
 1.10.03  24.02.09  we          mpf_round
 1.10.04  25.02.09  we          s_mpf_ldx(0)=-MaxLongint, mpf_ln1p for a=0,
 1.10.05  25.02.09  we          s_mpf_numbpart, mpf_numbpart
 1.10.06  25.02.09  we          mpf_add_int, mpf_sub_int
 1.10.07  27.02.09  we          s_mpf_frac, mpf_trig_ex

 1.11.00  23.03.09  we          mpf_sumaltf
 1.11.01  30.03.09  we          removed redefinition of pByte

 1.12.00  05.07.09  we          D12 fixes in mpf_read_radix and s_mpf_toradix_n

 1.13.00  14.08.09  we          small improvements in (s_)mpf_numbpart
 1.13.01  23.08.09  we          mpf_squad
 1.13.02  01.11.09  we          New names: mpf_arccosh,mpf_arccosh1p,mpf_arcsinh,mpf_arctanh,mpf_exp10,mpf_exp2
 1.13.03  01.11.09  we          mpf_cot,mpf_csc,mpf_sec,mpf_coth,mpf_csch,mpf_sech
 1.13.04  02.11.09  we          mpf_arccot,mpf_arccotc,mpf_arccsc,mpf_arcsec,mpf_arccoth,mpf_arccsch,mpf_arcsech
 1.13.05  02.11.09  we          changed mpf_arccosh domain to a >= 1
 1.13.06  16.11.09  we          mpf_init[x], x=2..5
 1.13.07  24.11.09  we          mpf_arccotc calls mpf_arccot for a>0
 1.13.08  24.11.09  we          inv_func/func_inv with 32 guard bits

 1.14.00  31.01.10  we          mpf_sqrt1pm1
 1.14.01  07.02.10  we          mpf_logbase
 1.14.02  09.02.10  we          improved guard bits calculation in mpf_expt
 1.14.03  09.02.10  we          mpf_cmp_ext, mpf_cmp_mag_ext

 1.15.00  05.05.10  we          basic s_mpf_toradix_alt
 1.15.01  06.05.10  we          mpf_(a)decimal_alt, mpf_todecimal_alt, mpf_toradix_alt
 1.15.02  07.05.10  we          mpf_write_radix_alt, mpf_write_decimal_alt
 1.15.03  07.05.10  we          mpf_output_radix_alt, mpf_output_decimal_alt
 1.15.04  08.05.10  we          s_mpf_toradix_alt: k = ndd-e-1
 1.15.05  08.05.10  we          mpf_round: fix quirk if @b = @a.mantissa
 1.15.07  09.05.10  we          mpf_todouble, mpf_toextended: pre-check max values
 1.15.08  13.05.10  we          mp_fract_sep
 1.15.09  13.05.10  we          improved s_mpf_toradix_alt: use s_mp_toradix_n for fractional part
 1.15.09  13.05.10  we          mpf_set_default_decprec
 1.15.10  14.05.10  we          improved s_mpf_toradix_alt: faster removal of trailing zeros
 1.15.11  15.05.10  we          s_mpf_toradix_alt: use IsPow2_w
 1.15.12  15.05.10  we          s_mpf_toradix_alt: improved for radix<>2^n

 1.16.00  12.06.10  we          Corrected some exception strings
 1.16.01  13.06.10  we          mpf_set_exp1, mpf_set_exp1p
 1.16.02  17.07.10  we          mpf_expt1pm1, mpf_exptm1

 1.18.00  17.02.11  we          mpf_lambertw
 1.18.01  23.06.11  we          mpf_sinc
 1.18.02  23.06.11  we          mpf_hav. mpf_archav
 1.18.03  24.06.11  we          mpf_gd, mpf_arcgd
 1.18.04  24.06.11  we          improved mpf_tanh
 1.18.05  24.06.11  we          mpf_coshm1

 1.19.00  12.11.11  we          Avoid some warnings for MAXDigits=32000 with MP_32BIT

 1.20.00  16.01.12  we          mpf_set_ext removed absolute
 1.20.01  16.01.12  we          mpf_set_dbl
 1.20.02  16.01.12  we          EXT64: redirect mpf_set_ext, mpf_toextended
 1.20.03  16.01.12  we          LambertW approximation with double
 1.20.04  12.02.12  we          s_mpf_lnagm
 1.20.05  13.02.12  we          improved mpf_ln1p

 1.21.00  30.05.12  we          faster s_mpf_numbpart
 1.21.01  16.07.12  we          s_mpf_rint

 1.24.00  15.12.12  we          Some word types changed to longint
 1.24.01  23.12.12  we          Removed MPC_Old_EleFun function names
 1.24.02  03.01.13  we          s_mpf_normalize: Fix type of up to TNInt

 1.25.00  23.01.13  we          improved mpf_todouble/mpf_toextended
 1.25.01  31.01.13  we          mpf_vers

 1.26.00  20.02.13  we          mpf_set_ext for denormals
 1.26.01  30.07.13  we          mpf_exp10m1, mpf_exp2m1
 1.26.02  31.07.13  we          mpf_log10p1, mpf_log2p1

 1.27.00  16.12.13  we          mpf_hypot
 1.27.01  09.02.14  we          mpf_divr_ext
 1.27.02  10.02.14  we          mpf_sincos
 1.27.03  11.02.14  we          mpf_exp10i
 1.27.04  13.02.14  we          s_mpf_is_gt0
 1.27.05  15.02.14  we          mpf_nroot

 1.30.00  16.08.14  we          delete superfluous second s_mpf_ldx(a) in mpf_exp

 1.32.00  09.04.15  we          BIT16-compatibility functions mpf_adecimal, mpf_adecimal_alt

 1.35.00  26.07.17  we          From fdamsupp: mpf_arccosd, mpf_arccotcd, mpf_arccotd,
                                mpf_arcsind, mpf_arctand, mpf_cbrt, mpf_cosd,
                                mpf_cospi, mpf_cotd, mpf_deg2rad, mpf_rad2deg,
                                mpf_rem_2pi, mpf_sind, mpf_sinpi, mpf_tand, mpf_tanpi
 1.35.01  27.07.17  we          Temp. var with increased precision in mpf_cbrt

 1.37.00  10.05.18  we          mpf_cell1
 1.37.01  10.05.18  we          mpf_cell2
 1.37.02  11.05.18  we          mpf_ccell1, mpf_ccell12, mpf_ccell2 for |k| > 1
 1.37.03  11.05.18  we          improved mpf_cell2

 1.38.00  25.06.18  we          mpf_??_ext changed to mpf_??_dbl
 1.38.01  26.06.18  we          separate functionw for mpf_add/sub/mul/div_ext
 1.38.02  27.06.18  we          const PINC for default precision elevation,
 1.38.03  28.06.18  we          some manual precision adjustments

 1.39.00  12.11.18  we          fix copy/paste bug in mpf_sinhcosh
 1.39.01  12.11.18  we          Error MP_PRECISION in s_mpf_mod_pi2k instead of assert
 1.39.02  16.11.18  we          fix subtle bug for mpf_arctan2(y,0)
**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2007-2018 Wolfgang Ehrhardt

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

type
  sumalt_term_func = function (k: longint; var num, den: longint): boolean;
  {sumalt_term_func must be of type (-1)^k*x(k), x(k)=num(k)/den(k) > 0}
  {result = false if x(k) cannot be evaluated, e.g. den > MaxLongint.}

  sumalt_fterm_func = function (k: longint; var term: mp_float): boolean;
  {sumalt_fterm_func must be of type (-1)^k*x(k), x(k) >= 0}
  {result = false if x(k) cannot be evaluated}


procedure mpf_abs(const a: mp_float; var b: mp_float);
  {-Absolute value, b = |a|}

procedure mpf_add(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a+b}

procedure mpf_add_dbl(const a: mp_float; b: double; var c: mp_float);
  {-Calculate c = a+b}

procedure mpf_add_int(const a: mp_float; b: longint; var c: mp_float);
  {-Calculate c = a+b}

procedure mpf_add_mpi(const a: mp_float; const b: mp_int; var c: mp_float);
  {-Calculate c = a+b}

function  mpf_adecimal(const a: mp_float; ndd: word): ansistring;
  {-Convert to decimal scientific representation with ndd digits, max 65000 digits for 32+ bit}

function  mpf_adecimal_alt(const a: mp_float; ndd: word): ansistring;
  {-Convert to decimal alternative representation with ndd digits, max 65000 digits for 32+ bit}

procedure mpf_agm(const a,b: mp_float; var c: mp_float);
  {-Calculate c = AGM(|a|,|b|)}

procedure mpf_arccos(const a: mp_float; var b: mp_float);
  {-Calculate b = arccos(a), |a| <= 1}

procedure mpf_arccosd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccos(a), b in degrees}

procedure mpf_arccosh(const a: mp_float; var b: mp_float);
  {-Calculate b = arccosh(a), a >= 1. Note: for a near 1 the function}
  { arccosh1p(a-1) should be used to reduce cancellation errors!}

procedure mpf_arccosh1p(const a: mp_float; var b: mp_float);
  {-Calculate b = arccosh(1+a), a>=0}

procedure mpf_arccot(const a: mp_float; var b: mp_float);
  {-Calculate the sign symmetric circular cotangent b = arccot(a) = arctan(1/a)}

procedure mpf_arccotc(const a: mp_float; var b: mp_float);
  {-Calculate the continuous circular cotangent b = arccotc(a) = Pi/2 - arctan(a)}

procedure mpf_arccotcd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccotc(a), b in degrees}

procedure mpf_arccotd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccot(a), b in degrees}

procedure mpf_arccoth(const a: mp_float; var b: mp_float);
  {-Calculate the inverse hyperbolic cotangent b = arccoth(a), |a| > 1}

procedure mpf_arccsc(const a: mp_float; var b: mp_float);
  {-Calculate the inverse circular cosecant b = arccsc(a), |a| >= 1}

procedure mpf_arccsch(const a: mp_float; var b: mp_float);
  {-Calculate the inverse hyperbolic cosecant b = arccsch(a)}

procedure mpf_arcgd(const a: mp_float; var b: mp_float);
  {-Calculate the inverse Gudermannian function b = arcgd(a) = arcsinh(tan(a))}

procedure mpf_archav(const a: mp_float; var b: mp_float);
  {-Calculate the inverse haversine b = archav(a) = 2*arcsin(sqrt(a))}

procedure mpf_arcsec(const a: mp_float; var b: mp_float);
  {-Calculate the inverse circular secant b = arcsec(a), |a| >= 1}

procedure mpf_arcsech(const a: mp_float; var b: mp_float);
  {-Calculate the inverse hyperbolic secant b = arcsech(a), 0 < a <= 1}

procedure mpf_arcsin(const a: mp_float; var b: mp_float);
  {-Calculate b = arcsin(a), |a| <= 1}

procedure mpf_arcsind(const a: mp_float; var b: mp_float);
  {-Calculate b = arcsin(a), b in degrees}

procedure mpf_arctan(const a: mp_float; var b: mp_float);
  {-Calculate b = arctan(a)}

procedure mpf_arctan2(const y,x: mp_float; var a: mp_float);
  {-Calculate a = arctan(y/x) with special treatment for zero x or y. a is}
  { the principal value of arg(x + i*y), i.e. -pi < arctan2(y,x) <= pi.}

procedure mpf_arctand(const a: mp_float; var b: mp_float);
  {-Calculate b = arctan(a), b in degrees}

procedure mpf_arcsinh(const a: mp_float; var b: mp_float);
  {-Calculate b = arcsinh(a)}

procedure mpf_arctanh(const a: mp_float; var b: mp_float);
  {-Calculate b = arctanh(a), |a| < 1}

procedure mpf_cbrt(const a: mp_float; var b: mp_float);
  {-Calculate b = a^(1/3) }

procedure mpf_cell1(const k: mp_float; var KK: mp_float);
  {-Calculate the complete elliptic integral of the first kind KK := K(k), |k|<>1, real part if |k|>1}

procedure mpf_cell2(const k: mp_float; var EK: mp_float);
  {-Calculate the complete elliptic integral of the 2nd kind EE := E(k), real part if |k|>1}

procedure mpf_ccell1(const k: mp_float; var CK: mp_float);
  {-Complementary complete elliptic integral of the first kind CK = CK(k)}
  { with k<>0 using AGM algorithm, real part if k<-1}

procedure mpf_ccell12(const k: mp_float; var CK, CE: mp_float);
  {-Complementary complete elliptic integrals of the 1st and 2nd kind using}
  { AGM algorithm; k<>0 and pCK <> pCE, with init checks, real parts of K'(k) if k<-1}

procedure mpf_ccell2(const k: mp_float; var CE: mp_float);
  {-Complementary complete elliptic integral of the 2nd kind CE = CE(k)}

function  mpf_checksum(const a: mp_float): longint;
  {-Return a checksum for a, -1 if mp_error<>MP_OKAY, -2 if not initialized}

procedure mpf_chg_prec(var a: mp_float; newprec: longint);
  {-Change bitprec of a to newprec}

procedure mpf_chs(const a: mp_float; var b: mp_float);
  {-Change sign, b = -a}

procedure mpf_clear(var a: mp_float);
  {-Clear an mp_float}

procedure mpf_clear2(var a,b: mp_float);
  {-Clear 2 mp_floats}

procedure mpf_clear3(var a,b,c: mp_float);
  {-Clear 3 mp_floats}

procedure mpf_clear4(var a,b,c,d: mp_float);
  {-Clear 4 mp_floats}

procedure mpf_clear5(var a,b,c,d,e: mp_float);
  {-Clear 5 mp_floats}

function  mpf_cmp(const a,b: mp_float): integer;
  {-Compare two mp_floats, return sign(a-b)}

function  mpf_cmp_dbl(const a: mp_float; b: double): integer;
  {-Compare a and b, return sign(a-b)}

function  mpf_cmp_mag(const a,b: mp_float): integer;
  {-Compare magnitude of two mp_floats, return sign(|a|-|b|)}

function  mpf_cmp_mag_dbl(const a: mp_float; b: double): integer;
  {-Compare magnitude of a and b, return sign(|a|-|b|)}

procedure mpf_copy(const a: mp_float; var b: mp_float);
  {-Copy a to b with b.bitprec=a.bitprec}

procedure mpf_copyp(const a: mp_float; var b: mp_float);
  {-Copy a to b, preserve b.bitprec}

procedure mpf_cos(const a: mp_float; var b: mp_float);
  {-Calculate b = cos(a)}

procedure mpf_cosd(const a: mp_float; var b: mp_float);
  {-Calculate b = cos(a), a in degrees}

procedure mpf_cospi(const a: mp_float; var b: mp_float);
  {-Calculate b = cos(Pi*a)}

procedure mpf_cosh(const a: mp_float; var b: mp_float);
  {-Calculate b = cosh(a), |a| < 2^31 * ln(2)}

procedure mpf_coshm1(const a: mp_float; var b: mp_float);
  {-Calculate b = cosh(a)-1, |a| < 2^31 * ln(2); special version for small a}

procedure mpf_cot(const a: mp_float; var b: mp_float);
  {-Calculate the circular cotangent b := cot(a), a mod Pi <> 0}

procedure mpf_cotd(const a: mp_float; var b: mp_float);
  {-Calculate b = cot(a), a in degrees}

procedure mpf_coth(const a: mp_float; var b: mp_float);
  {-Calculate the hyperbolic cotangent b = coth(a), a <> 0}

procedure mpf_csc(const a: mp_float; var b: mp_float);
  {-Calculate the circular cosecant b = csc(a), a mod Pi <> 0}

procedure mpf_csch(const a: mp_float; var b: mp_float);
  {-Calculate the hyperbolic cosecant b = csch(a), a  <> 0}

function  mpf_decimal(const a: mp_float; ndd: word): mp_string;
  {-Convert to decimal scientific representation with ndd digits, max 255 chars}

function  mpf_decimal_alt(const a: mp_float; ndd: word): mp_string;
  {-Convert to decimal alternative representation with ndd digits, max 255 chars}

procedure mpf_deg2rad(const a: mp_float; var b: mp_float);
  {-Convert a in degrees to b in radians}

procedure mpf_div(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a/b}

procedure mpf_div_d(const a: mp_float; b: mp_digit; var c: mp_float);
  {-Calculate c = a/b}

procedure mpf_div_dbl(const a: mp_float; b: double; var c: mp_float);
  {-Calculate c = a/b}

procedure mpf_div_int(const a: mp_float; b: longint; var c: mp_float);
  {-Calculate c = a/b}

procedure mpf_div_mpi(const a: mp_float; const b: mp_int; var c: mp_float);
  {-Calculate c = a/b}

procedure mpf_divr_dbl(a: double; const b: mp_float; var c: mp_float);
  {-Calculate c = a/b}

procedure mpf_exch(var a,b: mp_float);
  {-Exchange two mp_floats (including bitprec)}

procedure mpf_exp(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(a), a < 2^31 * ln(2)}

procedure mpf_exp10(const a: mp_float; var b: mp_float);
  {-Calculate b = 10^a}

procedure mpf_exp10i(n: longint; var a: mp_float);
  {-Calculate a = 10^n}

procedure mpf_exp10m1(const a: mp_float; var b: mp_float);
  {-Calculate b = 10^a-1}

procedure mpf_exp2(const a: mp_float; var b: mp_float);
  {-Calculate b = 2^a}

procedure mpf_exp2m1(const a: mp_float; var b: mp_float);
  {-Calculate b = 2^a-1}

procedure mpf_expm1(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(a)-1, a < 2^31 * ln(2); special version for small a}

procedure mpf_expt(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a^b, a>0}

procedure mpf_expt1pm1(const a,b: mp_float; var c: mp_float);
  {-Calculate c = (1+a)^b-1, a>-1}

procedure mpf_exptm1(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a^b-1, a>0}

procedure mpf_expt_int(const a: mp_float; b: longint; var c: mp_float);
  {-Calculate c = a^b}

procedure mpf_frac(const a: mp_float; var b: mp_float);
  {-Set b to the fractional part of a; frac(x)=x-int(x)}

procedure mpf_gd(const a: mp_float; var b: mp_float);
  {-Calculate the Gudermannian function b = gd(a) = arctan(sinh(a))}

function  mpf_get_default_prec: longint;
  {-Return current default (bit) precision, initial=240}

procedure mpf_hav(const a: mp_float; var b: mp_float);
  {-Calculate the haversine b = hav(a) = 0.5*(1 - cos(a))}

procedure mpf_hypot(const a,b: mp_float; var c: mp_float);
  {-Calculate c = sqrt(a^2 + b^2)}

procedure mpf_iexpt(a: longint; const b: mp_float; var c: mp_float);
  {-Calculate c = a^b, a>0}

procedure mpf_init(var a: mp_float);
  {-Initialize an mp_float with default precision}

procedure mpf_init2(var a,b: mp_float);
  {-Initialize two mp_floats with default precision}

procedure mpf_init3(var a,b,c: mp_float);
  {-Initialize 3 mp_floats with default precision}

procedure mpf_init4(var a,b,c,d: mp_float);
  {-Initialize 4 mp_floats with default precision}

procedure mpf_init5(var a,b,c,d,e: mp_float);
  {-Initialize 5 mp_floats with default precision}

procedure mpf_initp(var a: mp_float; prec: longint);
  {-Initialize an mp_float with bit precision prec}

procedure mpf_initp2(var a,b: mp_float; prec: longint);
  {-Initialize two mp_floats with bit precision prec}

procedure mpf_initp3(var a,b,c: mp_float; prec: longint);
  {-Initialize 3 mp_floats with bit precision prec}

procedure mpf_initp4(var a,b,c,d: mp_float; prec: longint);
  {-Initialize 4 mp_floats with bit precision prec}

procedure mpf_initp5(var a,b,c,d,e: mp_float; prec: longint);
  {-Initialize 5 mp_floats with bit precision prec}

procedure mpf_initp_multi_p(var pv: array of pmp_float; prec: longint);
  {-Initialize with bit precision prec a list of mp_floats given as a pointer}
  { vector; on error the already initialized mp_floats will be cleared}

procedure mpf_int(const a: mp_float; var b: mp_float);
  {-Set b to the integer part of a; i.e. is b rounded toward zero}

procedure mpf_inv(const a: mp_float; var b: mp_float);
  {-Calculate b = 1/a}

function  mpf_is0(const a: mp_float): boolean;
  {-Return true if a=0}

function  mpf_is1(const a: mp_float): boolean;
  {-Return true if a=1}

function  mpf_is1a(const a: mp_float): boolean;
  {-Return true if abs(a)=1}

function  mpf_is_eq(const a,b: mp_float): boolean;
  {-Return a = b}

function  mpf_is_eq_rel(const a,b: mp_float): boolean;
  {-Check if |a-b| <= r*2^(1-b.bitprec);  r=1 if b=0, r=|b| otherwise}

function  mpf_is_ge(const a,b: mp_float): boolean;
  {-Return a >= b}

function  mpf_is_gt(const a,b: mp_float): boolean;
  {-Return a > b}

function  mpf_is_le(const a,b: mp_float): boolean;
  {-Return a <= b}

function  mpf_is_lt(const a,b: mp_float): boolean;
  {-Return a < b}

function  mpf_is_ne(const a,b: mp_float): boolean;
  {-Return a <> b}

procedure mpf_lambertw(const a: mp_float; var b: mp_float);
  {-Calculate b = LambertW(a); principal branch, a >= -1/e}

procedure mpf_ln(const a: mp_float; var b: mp_float);
  {-Calculate b = ln(a), a>0}

procedure mpf_ln1p(const a: mp_float; var b: mp_float);
  {-Calculate b = ln(1+a), a > -1; special version for small a}

procedure mpf_log10(const a: mp_float; var b: mp_float);
  {-Calculate b = log10(a), a>0}

procedure mpf_log10p1(const a: mp_float; var b: mp_float);
  {-Calculate b = log10(1+a), a > -1}

procedure mpf_log2(const a: mp_float; var b: mp_float);
  {-Calculate b = log2(a), a>0}

procedure mpf_log2p1(const a: mp_float; var b: mp_float);
  {-Calculate b = log2(1+a),  a > -1}

procedure mpf_logbase(const b,x: mp_float; var y: mp_float);
  {-Calculate y = base b logarithm of x}

procedure mpf_mul(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a*b}

procedure mpf_mul_2k(const a: mp_float; k: longint; var b: mp_float);
  {-Calculate b = a*2^k}

procedure mpf_mul_d(const a: mp_float; b: mp_digit; var c: mp_float);
  {-Multiply by a digit}

procedure mpf_mul_dbl(const a: mp_float; b: double; var c: mp_float);
  {-Calculate c = a*b}

procedure mpf_mul_int(const a: mp_float; b: longint; var c: mp_float);
  {-Multiply by a 32 bit integer}

procedure mpf_mul_mpi(const a: mp_float; const b: mp_int; var c: mp_float);
  {-Calculate c = a*b}

function  mpf_not_init(const a: mp_float): boolean;
  {-Sanity check if a is initialized, does not catch all cases!}

procedure mpf_nroot(const a: mp_float; n: longint; var b: mp_float);
  {-Calculate the nth root of a: b = a^(1/n); n<>0, a >= 0  if n is even}

procedure mpf_numbpart(n: longint; var p: mp_int);
  {-Calculate the number of partitions of n with Hardy-Ramanujan-Rademacher formula}

procedure mpf_output_decimal(const a: mp_float; ndd: word);
  {-Write an mp_float to output using decimal scientific representation, ndd>0}
  { is the total number of digits (including one digit before the '.')}

procedure mpf_output_decimal_alt(const a: mp_float; ndd: word);
  {-Write an mp_float to output using decimal alternative representation, ndd>0 is}
  { the total number of digits, trailing '.' or '0' are suppressed. If a is too}
  { large or too small, scientific representation is used.}

procedure mpf_output_radix(const a: mp_float; radix,ndd: word);
  {-Write an mp_float to output using radix scientific representation, ndd>0}
  { is the total number of digits (including one digit before the '.')}

procedure mpf_output_radix_alt(const a: mp_float; radix,ndd: word);
  {-Write an mp_float to output using radix alternative representation, ndd>0 is}
  { the total number of digits, trailing '.' or '0' are suppressed. If a is too}
  { large or too small, scientific representation is used. NOTE: no radix prefix/suffix!}

procedure mpf_rad2deg(const a: mp_float; var b: mp_float);
  {-Convert a in radians to b in degrees}

procedure mpf_random(var a: mp_float);
  {-Set to a random number uniformly distributed in [0,1)}

procedure mpf_read_decimal(var a: mp_float; str: pchar8);
  {-Read a from ASCII float decimal string. str may contain a single '.'}
  { The exponent part is @[+|-]nnn (e or E can replace @). Integer or }
  { fractional part must be present.}

procedure mpf_read_hex(var a: mp_float; str: pchar8);
  {-Read a from ASCII float hexadecimal string. str may contain a single}
  { '.'. The exponent part is @[+|-]nnn (h or H can replace @). Integer }
  { or fractional part must be present.}

procedure mpf_read_radix(var a: mp_float; str: pchar8; radix: word);
  {-Read a from an ASCII float radix string. str may contain a single '.'.}
  { The exponent part <xp> is @[+|-]nnn, b/B, e/E, or h/H can replace @ if}
  { radix in [2,10,16]. The integer <ip> or fractional part <fp> must be }
  { present, nnn must be decimal, ie general format is <ip>.<fp>*radix^<xp>.}

function   mpf_reldev(const a,b: mp_float): double;
  {-Return abs((a-b)/b)*2^b.bitprec, special if b=0, or a-b=0}

procedure mpf_rem_2pi(const a: mp_float; var b: mp_float);
  {-Calculate b = a mod (2*Pi)}

procedure mpf_round(const a: mp_float; var b: mp_int);
  {-Round an mp_float to nearest mp_int}

procedure mpf_sec(const a: mp_float; var b: mp_float);
  {-Calculate the circular secant b = sec(a), a mod Pi <> Pi/2}

procedure mpf_sech(const a: mp_float; var b: mp_float);
  {-Calculate the hyperbolic secant b = sech(a)}

procedure mpf_set0(var a: mp_float);
  {-Set a=0, a.bitprec is preserved}

procedure mpf_set1(var a: mp_float);
  {-Set a=1, a.bitprec is preserved}

procedure mpf_set_default_decprec(dprec: word);
  {-Set default bit precision to dprec*log_2(10), i.e. dprec decimal digits}

procedure mpf_set_default_prec(prec: longint);
  {-Set new default (bit) precision}

procedure mpf_set_exp1(var a: mp_float);
  {-Set a to exp(1), preserve a.bitprec}

procedure mpf_set_exp1p(var a: mp_float; prec: longint);
  {-Set a to exp(1) with bit precision prec}

procedure mpf_set_dbl(var a: mp_float; b: double);
  {-Set a to a double. Error if b = NAN or INF}

procedure mpf_set_int(var a: mp_float; b: longint);
  {-Set a to a longint}

procedure mpf_set_ln10(var a: mp_float);
  {-Set a to ln(10), preserve a.bitprec}

procedure mpf_set_ln10p(var a: mp_float; prec: longint);
  {-Set a to ln(10) with bit precision prec}

procedure mpf_set_ln10p2k(var a: mp_float; k,prec: longint);
  {-Set a to ln(10)*2^k with bit precision prec}

procedure mpf_set_ln2(var a: mp_float);
  {-Set a to ln(2), preserve a.bitprec}

procedure mpf_set_ln2p(var a: mp_float; prec: longint);
  {-Set a to ln(2) with bit precision prec}

procedure mpf_set_ln2p2k(var a: mp_float; k,prec: longint);
  {-Set a to ln(2)*2^k with bit precision prec}

procedure mpf_set_mpi(var a: mp_float; const b: mp_int);
  {-Set a to an mp_int}

procedure mpf_set_mpi2k(var a: mp_float; const m: mp_int; e: longint);
  {-Set a to m*2^e (build a from mantissa and exponent)}

procedure mpf_set_pi(var a: mp_float);
  {-Set a to pi, preserve a.bitprec}

procedure mpf_set_pi2k(var a: mp_float; k: longint);
  {-Set a to pi*2^k, preserve a.bitprec}

procedure mpf_set_pip(var a: mp_float; prec: longint);
  {-Set a to pi with bit precision prec}

procedure mpf_set_pip2k(var a: mp_float; k,prec: longint);
  {-Set a to pi*2^k with bit precision prec}

procedure mpf_sin(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(a)}

procedure mpf_sinc(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(a)/a}

procedure mpf_sincos(const a: mp_float; var b,c: mp_float);
  {-Calculate b = sin(a), c = cos(a);  @b<>@c}

procedure mpf_sind(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(a), a in degrees}

procedure mpf_sinh(const a: mp_float; var b: mp_float);
  {-Calculate b = sinh(a), |a| < 2^31 * ln(2)}

procedure mpf_sinhcosh(const a: mp_float; var b,c: mp_float);
  {-Calculate b = sinh(a), c = cosh(a);  @b<>@c,  |a| < 2^31 * ln(2)}

procedure mpf_sinpi(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(Pi*a)}

procedure mpf_sqr(const a: mp_float; var b: mp_float);
  {-Calculate b = a*a}

procedure mpf_sqrt(const a: mp_float; var b: mp_float);
  {-Calculate b = sqrt(a)}

procedure mpf_sqrt1pm1(const a: mp_float; var b: mp_float);
  {-Calculate b = sqrt(1+a)-1 with increased accuracy for a near 0, a >= -1}

function  mpf_squad(const a,b,c: mp_float; var x1,y1,x2,y2: mp_float): integer;
  {-Solve the quadratic equation a*x^2 + b*x + c = 0. Result is the number}
  { of different solutions: 0 (if a=b=0), 1 (x1), or 2 (x1 and x2). If the}
  { result is = -2, x1 + i*y1 and x2 + i*y2 are the two complex solutions.}
  { x1,y1,x2,y2 should be different variables and not the same as a,b or c}

procedure mpf_sub(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a-b}

procedure mpf_sub_int(const a: mp_float; b: longint; var c: mp_float);
  {-Calculate c = a-b}

procedure mpf_sub_dbl(const a: mp_float; b: double; var c: mp_float);
  {-Calculate c = a-b}

procedure mpf_sub_mpi(const a: mp_float; const b: mp_int; var c: mp_float);
  {-Calculate c = a-b}

procedure mpf_sumalt(term: sumalt_term_func; n: longint; var s: mp_float; var Err: integer);
  {-Calculate s=sum(i=0..n, term(i)) of alternating series with sumalt algorithm.}
  { If n<4, n is adjusted to bitprec of s. Err<>0 if a term cannot be evaluated.}

procedure mpf_sumaltf(fterm: sumalt_fterm_func; n: longint; var s: mp_float; var Err: integer);
  {-Calculate s=sum(i=0..n, term(i)) of alternating series with sumalt algorithm.}
  { If n<4, n is adjusted to bitprec of s. Err<>0 if a term cannot be evaluated.}

procedure mpf_tan(const a: mp_float; var b: mp_float);
  {-Calculate b = tan(a)}

procedure mpf_tand(const a: mp_float; var b: mp_float);
  {-Calculate b = tan(a), a in degrees}

procedure mpf_tanh(const a: mp_float; var b: mp_float);
  {-Calculate b = tanh(a)}

procedure mpf_tanpi(const a: mp_float; var b: mp_float);
  {-Return b := tan(Pi*a)}

procedure mpf_todecimal_alt(const a: mp_float; ndd: word; str: pchar8; maxlen: word);
  {-Convert an mp_float to decimal alternative representation, ndd>0 is the total}
  { number of digits, trailing '.' or '0' are suppressed. If a is too large or}
  { too small, scientific representation is used.}

procedure mpf_todecimal_n(const a: mp_float; ndd: word; str: pchar8; maxlen: word);
  {-Convert an mp_float to decimal scientific representation, ndd>0 is the}
  { total number of digits (including one digit before the '.')}

function  mpf_todouble(const a: mp_float): double;
  {-Convert a to double, +-inf if too large}

procedure mpf_tohex_n(const a: mp_float; ndd: word; str: pchar8; maxlen: word);
  {-Convert an mp_float to hexadecimal scientific representation, ndd>0 is the}
  { total number of digits (including one digit before the '.')}

procedure mpf_toradix_alt(const a: mp_float; radix, ndd: word; str: pchar8; maxlen: longint);
  {-Convert an mp_float to radix alternative representation, ndd>0 is the total}
  { number of digits, trailing '.' or '0' are suppressed. If a is too large or}
  { too small, scientific representation is used. NOTE: no radix prefix/suffix!}

procedure mpf_toradix_n(const a: mp_float; radix, ndd: word; str: pchar8; maxlen: longint);
  {-Convert an mp_float to radix scientific representation, ndd>0 is the}
  { total number of digits (including one digit before the '.')}

procedure mpf_trig(const a: mp_float; pc,ps,pt: pmp_float);
  {-Calculate pc^=cos(a), ps^=sin(a), pt^=tan(a) if pointers are <> nil}

procedure mpf_trig_ex(const a: mp_float; mulpi: boolean; pc,ps,pt: pmp_float);
  {-Calculate pc^=cos(a), ps^=sin(a), pt^=tan(a) if pointers are <> nil.}
  { If mulpi, calculate pc^=cos(a*Pi), ps^=sin(a*Pi), pt^=tan(a*Pi).}

procedure mpf_trunc(const a: mp_float; var b: mp_int);
  {-Truncate mp_float to mp_int}

procedure mpf_vers(const a: mp_float; var b: mp_float);
  {-Calculate the versine b = vers(a) = 1-cos(a)}

procedure mpf_writeln(const msg: mp_string; const a: mp_float; ndd: word);
  {-Writeln a to output with leading msg}

procedure mpf_write_decimal(var tf: system.text; const a: mp_float; ndd: word);
  {-Write an mp_float to file tf using decimal scientific representation, ndd>0}
  { is the total number of digits (including one digit before the '.')}

procedure mpf_write_decimal_alt(var tf: system.text; const a: mp_float; ndd: word);
  {-Write an mp_float to file tf using decimal alternative representation, ndd>0 is}
  { the total number of digits, trailing '.' or '0' are suppressed. If a is too}
  { large or too small, scientific representation is used.}

procedure mpf_write_radix(var tf: system.text; const a: mp_float; radix,ndd: word);
  {-Write an mp_float to file tf using radix scientific representation, ndd>0}
  { is the total number of digits (including one digit before the '.')}

procedure mpf_write_radix_alt(var tf: system.text; const a: mp_float; radix,ndd: word);
  {-Write an mp_float to file tf using radix alternative representation, ndd>0 is}
  { the total number of digits, trailing '.' or '0' are suppressed. If a is too}
  { large or too small, scientific representation is used. NOTE: no radix prefix/suffix!}


{$ifndef EXT64}
procedure mpf_set_ext(var a: mp_float; x: extended);
  {-Set a to an extended. Error if x = NAN or INF}

function  mpf_toextended(const a: mp_float): extended;
  {-Convert a to extended, +-inf if too large}

procedure mpf_add_ext(const a: mp_float; b: extended; var c: mp_float);
  {-Calculate c = a+b}

procedure mpf_sub_ext(const a: mp_float; b: extended; var c: mp_float);
  {-Calculate c = a-b}

procedure mpf_mul_ext(const a: mp_float; b: extended; var c: mp_float);
  {-Calculate c = a*b}

procedure mpf_div_ext(const a: mp_float; b: extended; var c: mp_float);
  {-Calculate c = a/b}

procedure mpf_divr_ext(a: extended; const b: mp_float; var c: mp_float);
  {-Calculate c = a/b}

{$endif}


{#Z+}
{---------------------------------------------------------------------------}
{- 'Internal' functions, don't use them unless you know what you are doing -}
{---------------------------------------------------------------------------}
{#Z-}

procedure s_mpf_abs(var a: mp_float);
  {-Absolute value of an mp_float, no init check}

procedure s_mpf_add1(const a: mp_float; var b: mp_float);
  {-Calculate b = a+1; no init checks}

procedure s_mpf_agm(const a,b: mp_float; var c: mp_float; ps: pmp_float);
  {-Calculate c = AGM(|a|,|b|), if ps<>nil set ps^=sum(2^k*c_k^2); ps<>@c, no init check}

procedure s_mpf_ccell12(const k: mp_float; pCK, pCE: pmp_float);
  {-Complementary complete elliptic integrals of the 1st and 2nd kind using}
  { AGM algorithm; k<>0 and pCK <> pCE, with init checks, real parts of K'(k) if k<-1}

procedure s_mpf_chs(var a: mp_float);
  {-Change sign of an mp_float, no init check}

function  s_mpf_cmp_mag(const a,b: mp_float): integer;
  {-Compare magnitude of two mp_floats, return sign(|a|-|b|), no init check}

procedure s_mpf_dec1(var a: mp_float);
  {-Calculate a = a-1; no init check}

procedure s_mpf_inc(var x: mp_float; const y: mp_float);
  {-Calculate x = x+y}

procedure s_mpf_frac(const a: mp_float; var b: mp_float; podd: pBoolean);
  {-Set b to the fractional part of a; frac(x)=x-int(x), if podd<>nil set podd^=odd(trunc(a))}

function  s_mpf_incf(var x: mp_float; const y: mp_float): boolean;
  {-Calculate x = x+y; return true if x is changed}

procedure s_mpf_incexp(var a: mp_float; x: longint);
  {-Increment a.exponent by x, do nothing if a.mantissa=0}

procedure s_mpf_inc1(var a: mp_float);
  {-Calculate a = a+1; no init check}

function  s_mpf_is_ge0(const a: mp_float): boolean;
  {-Return true if a>=0, no init checks}

function  s_mpf_is_gt0(const a: mp_float): boolean;
  {-Return true if a>0, no init checks}

function  s_mpf_is_le0(const a: mp_float): boolean;
  {-Return true if a<=0, no init checks}

function  s_mpf_is_neg(const a: mp_float): boolean;
  {-Return true if a<0, no init checks}

function  s_mpf_is0(const a: mp_float): boolean;
  {-Return true if a=0, no init checks}

function  s_mpf_ldx(const a: mp_float): longint;
  {-Return ldx(a) = a.exponent+mp_bitsize(a.mantissa), no init checks, a must be}
  { normalized, ldx(a) = floor(log2(|a|))+1, -MaxLongint if a=0, MaxLongint if overflow}

procedure s_mpf_lnagm(const a: mp_float; var b: mp_float);
  {-Calculate b = ln(a) with AGM algorithm, a>0; no init checks}

procedure s_mpf_mod_pi2k(const a: mp_float; k: integer; var b: mp_float; var oddm: boolean);
  {-Calculate b = a mod pi*2^k, c in [0,pi*2^k); oddm: odd multiple of pi*2^k}
  { was used for reduction; no init check; extended precision used if necessary.}

procedure s_mpf_normalize(var a: mp_float);
  {-Normalize an mp_float}

procedure s_mpf_normalizep(var a: mp_float; prec: longint);
  {-Normalize an mp_float to bit precision prec}

procedure s_mpf_numbpart(n: longint; var p: mp_float);
  {-Calculate the number of partitions of n with Hardy-Ramanujan-Rademacher formula}

procedure s_mpf_toradix_n(const a: mp_float; radix, ndd: word; var maxlen: longint; var pstr: pchar8);
  {-Convert an mp_float to radix scientific representation, ndd>0 is the}
  { total number of digits (including one digit before the '.')}

procedure s_mpf_toradix_alt(const a: mp_float; radix, ndd: word; var maxlen: longint; var pstr: pchar8);
  {-Convert an mp_float to radix alternative representation, ndd>0 is the total}
  { number of digits, trailing '.' or '0' are suppressed. If a is too large or}
  { too small, scientific representation is used. NOTE: no radix prefix/suffix!}

{#Z+}
procedure _CheckBitPrec(const a: mp_float);
  {-Check a.bitprec; used for arg checks if init check is done otherwise}

procedure _set_ptab_const(var a: mp_float; e0,prec: longint; ptab: pointer);
  {-Set a to [ptab] with bit precision prec; e0=ldx(a)}

procedure s_mpf_expnewt(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(a), a < 2^31 * ln(2), !!!TEST ONLY!!!}

procedure s_mpf_rint(var a: mp_float);
  {-Set a = round(a)}
{#Z-}


implementation


uses
  {$ifdef BIT16}
    mp_rc16,
  {$else}
    mp_rc32,
  {$endif}
  mp_base;

const
  PINC = 16; {Default precision increment in some functions}

var
  mpf_default_prec : longint;  {current default bit precision, initial=240}


{---------------------------------------------------------------------------}
procedure _CheckBitPrec(const a: mp_float);
  {-Check a.bitprec; used for arg checks if init check is done otherwise}
begin
  with a do begin
    if (bitprec<MPF_MIN_PREC) or (bitprec>MPF_MAX_PREC) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('_CheckBitPrec');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure _set_ptab_const(var a: mp_float; e0,prec: longint; ptab: pointer);
  {-Set a to [ptab] with bit precision prec; e0=ldx(a)}
var
  nb: longint;
begin
  if mp_error<>MP_OKAY then exit;
  if prec<MPF_MIN_PREC then prec := MPF_MIN_PREC;
  if prec>MPF_MAX_PREC then prec := MPF_MAX_PREC;
  if prec>MAX_TCBITS   then prec := MAX_TCBITS;
  nb := (prec+7) div 8;
  {read an additional byte (if possible) for rounding}
  if nb<MAX_TCBITS div 8 then inc(nb);
  mp_read_unsigned_bin(a.mantissa, PByte(ptab)^,nb);
  a.exponent := e0-mp_bitsize(a.mantissa);
  a.bitprec  := prec;
  s_mpf_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpf_exp10(const a: mp_float; var b: mp_float);
  {-Calculate b = 10^a}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b);
  {$endif}
  if mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  mpf_initp(t, b.bitprec + 32);
  if mp_error=MP_OKAY then begin
    mpf_set_ln10(t);
    mpf_mul(t,a,t);
    mpf_exp(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_exp2(const a: mp_float; var b: mp_float);
  {-Calculate b = 2^a}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b);
  {$endif}
  if mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  mpf_initp(t, b.bitprec + 32);
  if mp_error=MP_OKAY then begin
    mpf_set_ln2(t);
    mpf_mul(t,a,t);
    mpf_exp(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_abs(const a: mp_float; var b: mp_float);
  {-Absolute value, b = |a|}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
    _CheckBitPrec(b);
  {$endif}
  mp_abs(a.mantissa, b.mantissa);
  if @a<>@b then begin
    b.exponent := a.exponent;
    if a.bitprec<>b.bitprec then s_mpf_normalize(b);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_arccosh(const a: mp_float; var b: mp_float);
  {-Calculate b = arccosh(a), a >= 1. Note: for a near 1 the function}
  { arccosh1p(a-1) should be used to reduce cancellation errors!}
var
  t: mp_float;
  l2: longint;
begin
  if mp_error<>MP_OKAY then exit;

  mpf_copyp(a,b);
  if mpf_is1(b) then begin
    mpf_set0(b);
    exit;
  end;

  l2 := s_mpf_ldx(b);
  if (l2 <= 0) or (s_mpf_is_neg(b)) then begin
    {if |a| < 1 then arccosh is undefined}
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_arccosh: a < 1');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  if l2<2 then begin
    {if a<2 use arccosh(a) = arccosh(1+(a-1)) = arccosh1p(a-1)}
    s_mpf_dec1(b);
    mpf_arccosh1p(b,b);
    exit;
  end;

  mpf_initp(t, b.bitprec + PINC);
  if mp_error<>MP_OKAY then exit;

  if l2 > succ(t.bitprec div 2) then begin
    {arccosh(a) = ln(2a) if a^2+1=a^2}
    s_mpf_incexp(b,1);
    mpf_ln(b,b);
  end
  else begin
    {arccosh(a) = ln(a + sqrt(a^2 - 1))}
    mpf_sqr(b,t);
    s_mpf_dec1(t);
    mpf_sqrt(t,t);
    s_mpf_inc(t,b);
    mpf_ln(t,b);
  end;
  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccosh1p(const a: mp_float; var b: mp_float);
  {-Calculate b = arccosh(1+a), a>=0}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_arccosh1p');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  {a must be positive}
  if s_mpf_is_neg(a) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_arccosh1p: a < 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mpf_initp(t, b.bitprec + PINC);
  if mp_error<>MP_OKAY then exit;

  {arccosh(x) = ln(x + sqrt(x^2 - 1)), substituting x=1+a gives}
  {arccosh1p(a) = ln1p(a+sqrt(a*(a+2)))}

  mpf_set_int(t,2);
  s_mpf_inc(t,a);
  mpf_mul(t,a,t);
  mpf_sqrt(t,t);
  s_mpf_inc(t,a);
  mpf_ln1p(t,b);

  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpf_add(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a+b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_add');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a) then begin
    mpf_copyp(b,c);
    exit;
  end;
  if s_mpf_is0(b) then begin
    mpf_copyp(a,c);
    exit;
  end;
  if @a=@b then begin
    mpf_mul_2k(a,1,c);
    exit;
  end;
  {make local copy of a}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_copyp(a,x);
    s_mpf_inc(x,b);
    mpf_exch(x,c);  {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_add_dbl(const a: mp_float; b: double; var c: mp_float);
  {-Calculate c = a+b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_add_dbl');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a) then begin
    mpf_set_dbl(c,b);
    exit;
  end;
  if b=0.0 then begin
    mpf_copyp(a,c);
    exit;
  end;
  {make local copy of b}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_dbl(x,b);
    s_mpf_inc(x,a);
    mpf_exch(x,c); {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_add_int(const a: mp_float; b: longint; var c: mp_float);
  {-Calculate c = a+b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_add_int');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a) then begin
    mpf_set_int(c,b);
    exit;
  end;
  if b=0.0 then begin
    mpf_copyp(a,c);
    exit;
  end;
  {make local copy of b}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_int(x,b);
    s_mpf_inc(x,a);
    mpf_exch(x,c); {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_add_mpi(const a: mp_float; const b: mp_int; var c: mp_float);
  {-Calculate c = a+b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mp_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_add_mpi');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a) then begin
    mpf_set_mpi(c,b);
    exit;
  end;
  if mp_is0(b) then begin
    mpf_copyp(a,c);
    exit;
  end;
  {make local copy of b}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_mpi(x,b);
    s_mpf_inc(x,a);
    mpf_exch(x,c); {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{$ifndef BIT16}
{---------------------------------------------------------------------------}
function mpf_adecimal(const a: mp_float; ndd: word): ansistring;
  {-Convert to decimal scientific representation with ndd digits, max 65000 digits for 32+ bit}
const
  LMAX=65000;
var
  l,ls: longint;
  pc: pchar8;
  {$ifndef RESULT}
    Result: ansistring;
  {$endif}
begin
  mpf_adecimal := '';
  {rough estimate > ndd + <1-> + 1<.> + 2<E-> + 12<expo>+ 1<#0>}
  {arg checks are done by mp_radix_size}
  if mp_error<>MP_OKAY then exit;
  if ndd+17<=LMAX then begin
    ls := ndd+17;
    SetLength(Result, ls);
    for l:=1 to length(Result) do Result[l] := #0;
    pc := @Result[1];
    s_mpf_toradix_n(a, 10, ndd, ls, pc);
    if mp_error=MP_OKAY then begin
      l := length(Result);
      {trim trailing #0}
      while (l>0) and (Result[l]=#0) do dec(l);
      if l<>length(Result) then SetLength(Result, l);
      {$ifndef RESULT}
        mpf_adecimal := Result;
      {$endif}
    end;
  end;
end;


{---------------------------------------------------------------------------}
function  mpf_adecimal_alt(const a: mp_float; ndd: word): ansistring;
  {-Convert to decimal alternative representation with ndd digits, max 65000 digits for 32+ bit}
const
  LMAX=$FF00;
var
  l,ls: longint;
  pc: pchar8;
  {$ifndef RESULT}
    Result: ansistring;
  {$endif}
begin
  mpf_adecimal_alt := '';
  {rough estimate > ndd + 1<-> + 1<.> + 2<E-> + 12<expo>+ 1<#0>}
  {arg checks are done by mp_radix_size}
  if mp_error<>MP_OKAY then exit;
  if ndd+17<=LMAX then begin
    ls := ndd+17;
    SetLength(Result, ls);
    for l:=1 to length(Result) do Result[l] := #0;
    pc := @Result[1];
    s_mpf_toradix_alt(a, 10, ndd, ls, pc);
    if mp_error=MP_OKAY then begin
      l := length(Result);
      {trim trailing #0}
      while (l>0) and (Result[l]=#0) do dec(l);
      if l<>length(Result) then SetLength(Result, l);
      {$ifndef RESULT}
        mpf_adecimal_alt := Result;
      {$endif}
    end;
  end;
end;
{$else}
{---------------------------------------------------------------------------}
function mpf_adecimal(const a: mp_float; ndd: word): ansistring;
  {-Convert to decimal scientific representation with ndd digits, compatibility function}
begin
  mpf_adecimal := mpf_decimal(a,ndd);
end;

{---------------------------------------------------------------------------}
function  mpf_adecimal_alt(const a: mp_float; ndd: word): ansistring;
  {-Convert to decimal alternative representation with ndd digits, compatibility function}
begin
  mpf_adecimal_alt := mpf_decimal_alt(a,ndd);
end;
{$endif}


{---------------------------------------------------------------------------}
procedure mpf_agm(const a,b: mp_float; var c: mp_float);
  {-Calculate c = AGM(|a|,|b|)}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_agm');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mpf_agm(a,b,c,nil);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccos(const a: mp_float; var b: mp_float);
  {-Calculate b = arccos(a), |a| <= 1}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_arccos');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set_pi2k(b,-1);
    exit;
  end;

  if mpf_is1a(a) then begin
    if s_mpf_is_neg(a) then mpf_set_pi(b) else mpf_set0(b);
    exit;
  end;

  if s_mpf_ldx(a) > 0 then begin
    {if |a| > 1 then arccos is undefined}
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_arccos: |a| > 1');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mpf_initp(t, b.bitprec + PINC);
  if mp_error<>MP_OKAY then exit;

  {arccos(a) := arctan2(sqrt(1 - a^2), a)}
  mpf_sqr(a,t);
  s_mpf_chs(t);
  s_mpf_inc1(t);
  mpf_sqrt(t,t);
  mpf_arctan2(t,a,b);

  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpf_arcsin(const a: mp_float; var b: mp_float);
  {-Calculate b = arcsin(a), |a| <= 1}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_arcsin');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  if mpf_is1a(a) then begin
    if s_mpf_is_neg(a) then begin
      mpf_set_pi2k(b,-1);
      s_mpf_chs(b);
    end
    else mpf_set_pi2k(b,-1);
    exit;
  end;

  if s_mpf_ldx(a) > 0 then begin
    {if |a| > 1 then arcsin is undefined}
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_arcsin: |a| > 1');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mpf_initp(t, b.bitprec + PINC);
  if mp_error<>MP_OKAY then exit;

  {arcsin(a) := arctan2(a, sqrt(1 - a^2))}
  mpf_sqr(a,t);
  s_mpf_chs(t);
  s_mpf_inc1(t);
  mpf_sqrt(t,t);
  mpf_arctan2(a,t,b);

  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpf_arctan(const a: mp_float; var b: mp_float);
  {-Calculate b = arctan(a)}
var
  x,y,t: mp_float;
  d: mp_int;
  i,f,bprec: longint;
  neg,inv: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_arctan');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  {if a<0 then arctan(a)=-arctan(|a|)}
  neg := s_mpf_is_neg(a);

  {arctan(+-1) = +-pi/4}
  if mpf_is1a(a) then begin
    mpf_set_pi2k(b,-2);
    if neg then s_mpf_chs(b);
    exit;
  end;

  {save result precision and use slightly larger working precision}
  bprec := b.bitprec;
  i := bprec + PINC;
  if i>MPF_MAX_PREC then i:=MPF_MAX_PREC;

  mp_init(d);
  if mp_error<>MP_OKAY then exit;

  mpf_initp3(x,y,t,i);
  if mp_error<>MP_OKAY then begin
    mp_clear(d);
    exit;
  end;

  {set working precision for b, and x=|a|}
  s_mpf_normalizep(b,i);
  mpf_abs(a,x);

  {if x > 1 then arctan(x)=pi/2-arctan(1/x)}
  inv := s_mpf_ldx(x) > 0;
  if inv then mpf_inv(x,x);

  mpf_set1(t);
  f := 0;
  {Range reduction: make x < 2^-10}
  while s_mpf_ldx(x) >= -10 do begin
    {arctan(x) = 2*arctan(x/(1+sqrt(x*x+1)))}
    mpf_sqr(x,y);
    s_mpf_inc1(y);
    mpf_sqrt(y,y);
    s_mpf_inc1(y);
    mpf_div(x,y,x);
    inc(f);
    if mp_error<>MP_OKAY then break;
  end;

  {set y = x^2, and b to the first term}
  mpf_sqr(x,y);
  mpf_copy(x,b);

  {arctan(x) = sum( (-1)^i*x^(2i+1)/(2i+1) )}
  i := 1;
  repeat
    inc(i,2);
    {Too many terms or error?}
    if (i<0) or (mp_error<>MP_OKAY) then break;
    mp_set_int(d,i);
    if i and 3 = 3 then mp_set_int(d,-i) else mp_set_int(d,i);
    mpf_mul(x,y,x);
    mpf_div_mpi(x,d,t);
  until not s_mpf_incf(b,t);

  {undo range reduction}
  if f>0 then s_mpf_incexp(b, f);
  if inv then begin
    mpf_set_pi2k(t,-1);
    mpf_sub(t,b,b);
  end;

  {a was < 0}
  if neg then s_mpf_chs(b);

  {normalize to original precision}
  s_mpf_normalizep(b,bprec);

  mpf_clear3(x,y,t);
  mp_clear(d);
end;


{---------------------------------------------------------------------------}
procedure mpf_arctan2(const y,x: mp_float; var a: mp_float);
  {-Calculate a = arctan(y/x) with special treatment for zero x or y. a is}
  { the principal value of arg(x + i*y), i.e. -pi < arctan2(y,x) <= pi.}
var
  xneg,yneg: boolean;
  t: mp_float;
begin
  {
  arctan2(0,0) =    0
  arctan2(z,z) =    pi/4  if z > 0
  arctan2(z,z) = -3*pi/4  if z < 0
  arctan2(y,0) =    0     if y > 0
  arctan2(y,0) =    pi    if y < 0
  arctan2(0,x) =    pi/2  if x > 0
  arctan2(0,x) =   -pi/2  if x > 0
  }

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(x) or mpf_not_init(y) or mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_arctan2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  xneg := s_mpf_is_neg(x);

  if @x=@y then begin
    {special case if x same var as y: arg along the line y=x}
    if s_mpf_is0(x) then mpf_set0(a)
    else begin
      {pi/4 for x,y > 0 and -3*pi/4 for x,y < 0}
      mpf_set_pi2k(a,-2);
      if xneg then mpf_mul_int(a,-3,a);
    end;
    exit;
  end;

  if s_mpf_is0(y) then begin
    if xneg then mpf_set_pi(a) else mpf_set0(a);
    exit;
  end;

  {1.39.02: Make copy of y sign before setting a to pi/2,}
  {         because a and y may be the same variable!}
  yneg := s_mpf_is_neg(y);

  if s_mpf_is0(x) then begin
    mpf_set_pi2k(a,-1);
    if yneg then s_mpf_chs(a);
    exit;
  end;

  mpf_div(y,x,a);
  mpf_arctan(a,a);

  if xneg then begin
    mpf_initp(t,a.bitprec);
    if mp_error=MP_OKAY then begin
      mpf_set_pi(t);
      if yneg then mpf_sub(a,t,a)
      else mpf_add(a,t,a);
      mpf_clear(t);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_arctand(const a: mp_float; var b: mp_float);
  {-Calculate b = arctan(a), b in degrees}
begin
  mpf_arctan(a,b);
  mpf_rad2deg(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccotd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccot(a), b in degrees}
begin
  mpf_arccot(a,b);
  mpf_rad2deg(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccotcd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccotc(a), b in degrees}
begin
  mpf_arccotc(a,b);
  mpf_rad2deg(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arcsind(const a: mp_float; var b: mp_float);
  {-Calculate b = arcsin(a), b in degrees}
begin
  mpf_arcsin(a,b);
  mpf_rad2deg(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccosd(const a: mp_float; var b: mp_float);
  {-Calculate b = arccos(a), b in degrees}
begin
  mpf_arccos(a,b);
  mpf_rad2deg(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arcsinh(const a: mp_float; var b: mp_float);
  {-Calculate b = arcsinh(a)}
var
  t: mp_float;
  l2: longint;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_arcsinh');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  l2 := s_mpf_ldx(a);
  if l2 < -(b.bitprec + PINC) div 2 then begin
    {arcsinh(a) = a*(1 - a^2/6 + O(a^4))}
    mpf_copyp(a,b);
    exit;
  end;

  mpf_initp(t, b.bitprec + PINC);
  if mp_error<>MP_OKAY then exit;

  neg := s_mpf_is_neg(a);
  mpf_abs(a,b);

  if l2>0 then begin
    {|a| >= 1.0}
    if l2 > succ(t.bitprec div 2) then begin
      {arcsinh(a) = ln(2a) if a^2+1=a^2}
      s_mpf_incexp(b,1);
      mpf_ln(b,b);
    end
    else begin
      {arcsinh(a) = ln(a + sqrt(a^2 + 1))}
      mpf_sqr(a,t);
      s_mpf_inc1(t);
      mpf_sqrt(t,t);
      s_mpf_inc(t,b);
      mpf_ln(t,b);
    end;
  end
  else begin
    {arcsinh(a) = ln1p(a*(1 + a/(1 + sqrt(1 + a^2))))}
    mpf_sqr(a,t);
    s_mpf_inc1(t);
    mpf_sqrt(t,t);
    s_mpf_inc1(t);
    mpf_div(b,t,t);
    s_mpf_inc1(t);
    mpf_mul(b,t,t);
    mpf_ln1p(t,b);
  end;

  if neg then s_mpf_chs(b);

  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpf_arctanh(const a: mp_float; var b: mp_float);
  {-Calculate b = arctanh(a), |a| < 1}
var
  t: mp_float;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_arctanh');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  if s_mpf_ldx(a) > 0 then begin
    {if |a| >= 1 then arctanh is undefined}
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_arctanh: |a| >= 1');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mpf_initp(t, b.bitprec + PINC);
  if mp_error<>MP_OKAY then exit;

  {arctanh(x) =  0.5*ln((1+x)/(1-x)) = -0.5*ln((1-x)/(1+x))}
  {           = -0.5*ln(1-2x/(1+x))  = -0.5*ln1p(-2x/(1+x))}

  neg := not s_mpf_is_neg(a);
  mpf_abs(a,b);
  s_mpf_add1(b,t);
  mpf_div(b,t,t);
  s_mpf_incexp(t,1);
  s_mpf_chs(t);
  mpf_ln1p(t,b);
  s_mpf_incexp(b,-1);
  if neg then s_mpf_chs(b);

  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpf_cbrt(const a: mp_float; var b: mp_float);
  {-Calculate b = a^(1/3) }
var
  an: boolean;
  x: mp_float;
begin
  an := s_mpf_is_neg(a);
  mpf_initp(x, b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_abs(a,x);
    mpf_ln(x,x);
    mpf_div_dbl(x,3.0,x);
    mpf_exp(x,b);
    mpf_clear(x);
  end;
  if an then s_mpf_chs(b);
end;


{---------------------------------------------------------------------------}
procedure mpf_cell1(const k: mp_float; var KK: mp_float);
  {-Calculate the complete elliptic integral of the first kind KK := K(k), |k|<>1, real part if |k|>1}
var
  s,t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(k) or mpf_not_init(KK) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_cell1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mpf_cmp_mag_dbl(k,1) = MP_EQ then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mpf_cell1: |k| = 1');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;
  if mpf_is0(k) then begin
    mpf_set_pi2k(KK,-1);
    exit;
  end;
  mpf_initp2(s,t,KK.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_set1(s);
    mpf_sub(s,k,s);
    mpf_add_dbl(k,1,t);
    s_mpf_agm(s,t,t,nil);
    mpf_set_pi2k(s,-1);
    mpf_div(s,t,KK);
    mpf_clear2(s,t);
  end;
end;


(*
{---------------------------------------------------------------------------}
procedure mpf_cell2(const k: mp_float; var EK: mp_float);
  {-Calculate the complete elliptic integral of the 2nd kind EE := E(k)}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(k) or mpf_not_init(EK) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_cell2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_initp(t, EK.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    {compute t = k' = sqrt(1-k^2)}
    mpf_sqr(k,t);
    s_mpf_dec1(t);
    s_mpf_chs(t);
    mpf_sqrt(t,t);
    if s_mpf_is_gt0(t) then s_mpf_ccell12(t, nil, @EK)
    else mpf_set1(EK);
    mpf_clear(t);
  end;
end;
*)

{---------------------------------------------------------------------------}
procedure mpf_cell2(const k: mp_float; var EK: mp_float);
  {-Calculate the complete elliptic integral of the 2nd kind EE := E(k), real part if |k|>1}
var
  s,t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(k) or mpf_not_init(EK) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_cell2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_initp2(s,t,EK.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_set1(s);
    mpf_sub(s,k,s);
    mpf_add_dbl(k,1,t);
    if s_mpf_is0(s) or s_mpf_is0(t) then begin
      mpf_set1(EK);
    end
    else begin
      s_mpf_agm(s,t,t,@s);
      mpf_mul_dbl(s, -0.25, s);
      s_mpf_inc1(s);
      mpf_div(s,t,s);
      mpf_set_pi2k(t,-1);
      mpf_mul(s,t,EK);
    end;
    mpf_clear2(s,t);
  end;
end;

{---------------------------------------------------------------------------}
procedure mpf_ccell1(const k: mp_float; var CK: mp_float);
  {-Complementary complete elliptic integral of the first kind CK = CK(k)}
  { with k<>0 using AGM algorithm, real part if k<-1}
var
  s: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(k) or mpf_not_init(CK) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_ccell1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(k) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mpf_ccell1: k = 0');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  mpf_initp(s,CK.bitprec + PINC);

  {calculate AGM(1,k)}
  mpf_set1(s);
  s_mpf_agm(s,k,s,nil);

  {CK = pi/2/AGM = pi/2/CK}
  mpf_set_pi2k(CK,-1);
  mpf_div(CK,s,CK);

  mpf_clear(s);
end;


{---------------------------------------------------------------------------}
procedure mpf_ccell12(const k: mp_float; var CK, CE: mp_float);
  {-Complementary complete elliptic integrals of the 1st and 2nd kind using}
  { AGM algorithm; k<>0 and pCK <> pCE, with init checks, real parts of K'(k) if k<-1}
begin
  s_mpf_ccell12(k, @CK, @CE);
end;


{---------------------------------------------------------------------------}
procedure mpf_ccell2(const k: mp_float; var CE: mp_float);
  {-Complementary complete elliptic integral of the 2nd kind CE = CE(k)}
begin
  s_mpf_ccell12(k, nil, @CE);
end;


{---------------------------------------------------------------------------}
function mpf_checksum(const a: mp_float): longint;
  {-Return a checksum for a, -1 if mp_error<>MP_OKAY, -2 if not initialized}
var
  adler: longint;
begin
  mpf_checksum := -1;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      mpf_checksum := -2;
      exit;
    end;
  {$endif}
  adler := mp_checksum(a.mantissa);
  with a do begin
    s_mp_checksum(adler,@exponent,sizeof(exponent));
    s_mp_checksum(adler,@bitprec, sizeof(bitprec));
  end;
  mpf_checksum := adler;
end;


{---------------------------------------------------------------------------}
procedure mpf_chg_prec(var a: mp_float; newprec: longint);
  {-Change bitprec of a to newprec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_chg_prec');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if (newprec<MPF_MIN_PREC) or (newprec>MPF_MAX_PREC) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mpf_chg_prec: newprec out of range');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}
  s_mpf_normalizep(a, newprec);
end;


{---------------------------------------------------------------------------}
procedure mpf_chs(const a: mp_float; var b: mp_float);
  {-Change sign, b = -a}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
    _CheckBitPrec(b);
  {$endif}
  mp_chs(a.mantissa, b.mantissa);
  if @a<>@b then begin
    b.exponent := a.exponent;
    if a.bitprec<>b.bitprec then s_mpf_normalize(b);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_clear(var a: mp_float);
  {-Clear an mp_float}
begin
  with a do begin
    bitprec := 0;
    exponent:= 0;
    mp_clear(mantissa);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_clear2(var a,b: mp_float);
  {-Clear 2 mp_floats}
begin
  mpf_clear(a);
  mpf_clear(b);
end;


{---------------------------------------------------------------------------}
procedure mpf_clear3(var a,b,c: mp_float);
  {-Clear 3 mp_floats}
begin
  mpf_clear2(a,b);
  mpf_clear(c);
end;


{---------------------------------------------------------------------------}
procedure mpf_clear4(var a,b,c,d: mp_float);
  {-Clear 4 mp_floats}
begin
  mpf_clear2(a,b);
  mpf_clear2(c,d);
end;


{---------------------------------------------------------------------------}
procedure mpf_clear5(var a,b,c,d,e: mp_float);
  {-Clear 5 mp_floats}
begin
  mpf_clear2(a,b);
  mpf_clear2(c,d);
  mpf_clear(e);
end;


{---------------------------------------------------------------------------}
function mpf_cmp_dbl(const a: mp_float; b: double): integer;
  {-Compare a and b, return sign(a-b)}
var
  x: mp_float;
begin
  mpf_cmp_dbl := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_cmp_dbl');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {make local copy of b}
  mpf_initp(x, a.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_dbl(x,b);
    mpf_cmp_dbl := mpf_cmp(a,x);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
function mpf_cmp(const a,b: mp_float): integer;
  {-Compare two mp_floats, return sign(a-b)}
begin
  {Value for last alternative, keep D6+ happy}
  mpf_cmp := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_cmp');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {compare based on sign}
  if a.mantissa.sign<>b.mantissa.sign then begin
    if a.mantissa.sign=MP_NEG then mpf_cmp := -1 else mpf_cmp := 1;
    exit;
  end;
  {compare magnitude}
  if a.mantissa.sign=MP_NEG then begin
    {if negative compare opposite direction}
    mpf_cmp := s_mpf_cmp_mag(b, a);
  end
  else mpf_cmp := s_mpf_cmp_mag(a, b);
end;


{---------------------------------------------------------------------------}
function mpf_cmp_mag(const a,b: mp_float): integer;
  {-Compare magnitude of two mp_floats, return sign(|a|-|b|)}
begin
  {Value for last alternative, keep D6+ happy}
  mpf_cmp_mag := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_cmp_mag');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_cmp_mag := s_mpf_cmp_mag(a,b);
end;


{---------------------------------------------------------------------------}
function mpf_cmp_mag_dbl(const a: mp_float; b: double): integer;
  {-Compare magnitude of a and b, return sign(|a|-|b|)}
var
  x: mp_float;
begin
  mpf_cmp_mag_dbl := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_cmp_mag_dbl');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {make local copy of b}
  mpf_initp(x, a.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_dbl(x,b);
    mpf_cmp_mag_dbl := mpf_cmp_mag(a,x);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_copy(const a: mp_float; var b: mp_float);
  {-Copy a to b with b.bitprec=a.bitprec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  if @a<>@b then begin
    b.exponent:= a.exponent;
    b.bitprec := a.bitprec;
    mp_copy(a.mantissa, b.mantissa);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_copyp(const a: mp_float; var b: mp_float);
  {-Copy a to b, preserve b.bitprec}
var
  pb: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b);
  {$endif}
  pb := b.bitprec;
  mpf_copy(a,b);
  if b.bitprec<>pb then s_mpf_normalizep(b, pb);
end;


{---------------------------------------------------------------------------}
procedure mpf_cos(const a: mp_float; var b: mp_float);
  {-Calculate b = cos(a)}
begin
  mpf_trig(a, @b, nil, nil);
end;


{---------------------------------------------------------------------------}
procedure mpf_cosd(const a: mp_float; var b: mp_float);
  {-Calculate b = cos(a), a in degrees}
begin
  mpf_deg2rad(a,b);
  mpf_cos(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_cosh(const a: mp_float; var b: mp_float);
  {-Calculate b = cosh(a), |a| < 2^31 * ln(2)}
var
  t: mp_float;
  pb: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_cosh');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;

  pb := b.bitprec;
  mpf_initp(t, pb + PINC);
  if mp_error<>MP_OKAY then exit;

  mpf_abs(a,t);
  b.bitprec := t.bitprec;

  {cosh(t) = (exp(|t|)+exp(-|t|))/2}
  mpf_exp(t,b);

  {don't calculate inverse if it is to small}
  if s_mpf_ldx(b) < 16 + b.bitprec div 2 then begin
    mpf_inv(b,t);
    mpf_add(b,t,b);
  end;
  s_mpf_incexp(b,-1);
  s_mpf_normalizep(b,pb);

  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpf_coshm1(const a: mp_float; var b: mp_float);
  {-Calculate b = cosh(a)-1, |a| < 2^31 * ln(2); special version for small a}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_coshm1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  mpf_initp(t, b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    if s_mpf_ldx(a)<2 then begin
      {if |a|<2 then coshm1 = 2.0*sqr(sinh(0.5*x))}
      mpf_mul_2k(a,-1,t);
      mpf_sinh(t,t);
      mpf_sqr(t,t);
      mpf_mul_2k(t,1,t);
    end
    else begin
      mpf_cosh(a,t);
      s_mpf_dec1(t);
    end;
    mpf_copyp(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
function mpf_decimal(const a: mp_float; ndd: word): mp_string;
  {-Convert to decimal scientific representation with ndd digits, max 255 chars}
var
  i: integer;
  d1,d2: __P2I;
  pc: pchar8;
  iostr: string[255];
begin
  {arg checks are done by mpf_toradix_n}
  mpf_decimal := '';
  if (mp_error<>MP_OKAY) or (ndd=0) then exit;
  pc := @iostr[1];
  mpf_toradix_n(a, 10, ndd, pc, 255);
  if mp_error<>MP_OKAY then exit;
  i := 0;
  d1 := __P2I(pc);
  d2 := __P2I(@iostr[1]);
  if (d2>d1) and (d2<d1+255) then i := integer(d2-d1);
  if (i>0) and (i<255) and (iostr[i]<>#0) and (iostr[i+1]=#0) then iostr[0] := char8(i)
  else begin
    iostr[0] := #255;
    for i:=1 to 255 do begin
      if iostr[i]=#0 then begin
        iostr[0] := char8(i-1);
        break;
      end;
    end;
  end;
  mpf_decimal := iostr;
end;


{---------------------------------------------------------------------------}
function mpf_decimal_alt(const a: mp_float; ndd: word): mp_string;
  {-Convert to decimal alternative representation with ndd digits, max 255 chars}
var
  i: integer;
  d1,d2: __P2I;
  pc: pchar8;
  iostr: string[255];
begin
  {arg checks are done by mpf_toradix_n}
  mpf_decimal_alt := '';
  if (mp_error<>MP_OKAY) or (ndd=0) then exit;
  pc := @iostr[1];
  mpf_toradix_alt(a, 10, ndd, pc, 255);
  if mp_error<>MP_OKAY then exit;
  i := 0;
  d1 := __P2I(pc);
  d2 := __P2I(@iostr[1]);
  if (d2>d1) and (d2<d1+255) then i := integer(d2-d1);
  if (i>0) and (i<255) and (iostr[i]<>#0) and (iostr[i+1]=#0) then iostr[0] := char8(i)
  else begin
    iostr[0] := #255;
    for i:=1 to 255 do begin
      if iostr[i]=#0 then begin
        iostr[0] := char8(i-1);
        break;
      end;
    end;
  end;
  mpf_decimal_alt := iostr;
end;


{---------------------------------------------------------------------------}
procedure mpf_div(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a/b}
var
  x: mp_float;
  diff: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_div');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(b) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_div: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  if s_mpf_is0(a) then begin
    mpf_set0(c);
    exit;
  end;

  if @a=@b then begin
    mpf_set1(c);
    exit;
  end;

  {make local copy of a with precision c.bitprec}
  mpf_initp(x, c.bitprec);
  if mp_error<>MP_OKAY then exit;
  mpf_copyp(a,x);

  diff :=  mp_bitsize(b.mantissa) - mp_bitsize(x.mantissa) + 1 + x.bitprec;
  s_mpf_incexp(x, -b.exponent);
  if diff>0 then begin
    s_mpf_incexp(x, -diff);
    mp_shl(x.mantissa, diff, x.mantissa);
  end;
  mp_div(x.mantissa, b.mantissa, x.mantissa);
  s_mpf_normalize(x);
  mpf_exch(x,c); {Note: x.bitprec=c.bitprec}
  mpf_clear(x);
end;


{---------------------------------------------------------------------------}
procedure mpf_deg2rad(const a: mp_float; var b: mp_float);
  {-Convert a in degrees to b in radians}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  mpf_initp(t,b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_set_pi(t);
    mpf_mul(t,a,t);
    mpf_div_d(t,180,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_div_d(const a: mp_float; b: mp_digit; var c: mp_float);
  {-Calculate c = a/b}
var
  diff: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_div_d');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if b=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_div_d: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  if s_mpf_is0(a) then begin
    mpf_set0(c);
    exit;
  end;

  mpf_copyp(a,c);
  if b=1 then exit;

  {same calculation as mpf_div with b.exponent=0, b.mantissa=b}
  diff :=  bitsize32(b) - mp_bitsize(c.mantissa) + 1 + c.bitprec;
  if diff>0 then begin
    s_mpf_incexp(c, -diff);
    mp_shl(c.mantissa, diff, c.mantissa);
  end;
  mp_div_d(c.mantissa, b, @c.mantissa, b);
  s_mpf_normalize(c);
end;


{---------------------------------------------------------------------------}
procedure mpf_div_dbl(const a: mp_float; b: double; var c: mp_float);
  {-Calculate c = a/b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_div_dbl');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if b=0.0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_div_dbl: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mpf_initp(x,c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_dbl(x,b);
    mpf_div(a,x,c);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_divr_dbl(a: double; const b: mp_float; var c: mp_float);
  {-Calculate c = a/b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_div_dbl');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(b) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_divr_dbl: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mpf_initp(x,c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_dbl(x,a);
    mpf_div(x,b,c);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_div_int(const a: mp_float; b: longint; var c: mp_float);
  {-Calculate c = a/b}
var
  diff: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_div_int');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if b=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_div_int: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  if s_mpf_is0(a) then begin
    mpf_set0(c);
    exit;
  end;
  mpf_copyp(a,c);
  if abs(b)=1 then begin
    if b<0 then s_mpf_chs(c);
    exit;
  end;
  {same calculation as mpf_div with b.exponent=0, b.mantissa=b}
  diff :=  bitsize32(abs(b)) - mp_bitsize(c.mantissa) + 1 + c.bitprec;
  if diff>0 then begin
    s_mpf_incexp(c, -diff);
    mp_shl(c.mantissa, diff, c.mantissa);
  end;
  mp_div_int(c.mantissa, b, @c.mantissa, diff);
  s_mpf_normalize(c);
end;


{---------------------------------------------------------------------------}
procedure mpf_div_mpi(const a: mp_float; const b: mp_int; var c: mp_float);
  {-Calculate c = a/b}
var
  diff: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mp_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_div_mpi');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mp_is0(b) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_div_mpi: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  if s_mpf_is0(a) then begin
    mpf_set0(c);
    exit;
  end;
  {same calculation as mpf_div with b.exponent=0, b.mantissa=b}
  mpf_copyp(a,c);
  diff :=  mp_bitsize(b) - mp_bitsize(c.mantissa) + 1 + c.bitprec;
  if diff>0 then begin
    s_mpf_incexp(c, -diff);
    mp_shl(c.mantissa, diff, c.mantissa);
  end;
  mp_div(c.mantissa, b, c.mantissa);
  s_mpf_normalize(c);
end;


{---------------------------------------------------------------------------}
procedure mpf_exch(var a,b: mp_float);
  {-Exchange two mp_floats (including bitprec)}
var
  t: mp_float;
begin
  if mp_error=MP_OKAY then begin
    t := a;
    a := b;
    b := t;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_exp(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(a), a < 2^31 * ln(2)}
var
  x: mp_float;
  c,f,n,l2,bprec: longint;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_exp');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;

  {Overflow for a > 2^31 * ln(2)}
  l2 := s_mpf_ldx(a);
  {Safe first check for over/underflow. Does not catch x=1.49E9 .. 2.14E9}
  if l2 > 30 then begin
    if s_mpf_is_neg(a) then begin
      mpf_set0(b);
      exit;
    end
    else if l2>31 then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXOverflow.Create('mpf_exp');
        {$else}
          RunError(MP_RTE_OVRFLOW);
        {$endif}
      {$else}
        set_mp_error(MP_OVERFLOW);
        exit;
      {$endif}
    end;
  end;

  {save result precision and use larger working precision}
  bprec := b.bitprec;
  {get approx. sqr count in reconstruction step}
  c := trunc(sqrt(0.5*b.bitprec));
  f := 32+l2+c;
  if f<32 then f := 32;
  n := bprec+f;
  if n>MPF_MAX_PREC then n:=MPF_MAX_PREC;

  mpf_initp(x,n);
  if mp_error<>MP_OKAY then exit;

  {if a<0, use exp(a)=1/exp(|a|)}
  neg := s_mpf_is_neg(a);
  mpf_abs(a,x);

  {set working precision for b, init series with 1}
  b.bitprec := n;
  mpf_set1(b);

  {Algorithm from T. Papanikolaou [27], section 4.4.2}
  with x do begin
    {calculate number of steps and perform range reduction if necessary}
    l2:= s_mpf_ldx(x);
    c := trunc(sqrt(0.5*bitprec)); {Note: c>= sqrt(37/2) > 4}
    f := 1 + l2 + c;
    if f>0 then begin
      {Range reduction to (0..2^-c)}
      s_mpf_incexp(x, -f);
      n := 1 + bitprec div c;
    end
    else begin
      {0>=f = 1+l2+c --> l2 <= -1-c < -5  --> abs(l2) > 5}
      n := 1 + bitprec div abs(l2);
    end;
  end;

  {Taylor sum from n downto 0, b is initialized with 1}
  while n>0 do begin
    mpf_mul(b,x,b);
    if n <= lv_digit_max then mpf_div_d(b,mp_digit(n),b)
    else mpf_div_int(b,n,b);
    s_mpf_inc1(b);
    dec(n);
  end;

  {if f>0, reconstruct result from range reduction}
  for n:=1 to f do mpf_sqr(b,b);

  {if a<0 uses exp(a)=1/exp(|a|)}
  if neg then mpf_inv(b,b);

  {normalize to original precision}
  s_mpf_normalizep(b,bprec);

  mpf_clear(x);
end;


{---------------------------------------------------------------------------}
procedure mpf_expm1(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(a)-1, a < 2^31 * ln(2); special version for small a}
var
  x,y: mp_float;
  n,l2,bprec: longint;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_expm1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  if s_mpf_ldx(a) >= 0 then begin
    {if |a| > 0.5 then use standard function b = exp(a)-1}
    mpf_exp(a,b);
    s_mpf_dec1(b);
    exit;
  end;

  {save result precision and use slightly larger working precision}
  bprec := b.bitprec;
  mpf_initp2(x,y,bprec + PINC);
  if mp_error<>MP_OKAY then exit;

  {if a<0, note sign and calculate expm1(|a|)}
  neg := s_mpf_is_neg(a);
  mpf_abs(a,x);

  {range reduction x < 1/2^16}
  l2 := s_mpf_ldx(x) + 16;
  if l2>0 then dec(x.exponent,l2);

  {sum Taylor series for exp(x)-1}
  mpf_copy(x,y);
  mpf_copy(x,b);
  n := 1;
  repeat
    inc(n);
    mpf_mul(y,x,y);
    mpf_div_int(y,n,y);
    if mp_error<>MP_OKAY then break;
  until not s_mpf_incf(b,y);

  {undo range reduction, expm1(2x)=expm1(x)*(2+expm1(x))}
  if l2>0 then begin
    mpf_set_int(x,2);
    for n:=1 to l2 do begin
      mpf_add(b,x,y);
      mpf_mul(b,y,b);
    end;
  end;

  {if a<0 uses expm1(-a)=-expm1(a)/(1+expm1(a))}
  if neg then begin
    mpf_chs(b,x);
    s_mpf_inc1(b);
    mpf_div(x,b,b);
  end;

  {normalize to original precision}
  s_mpf_normalizep(b,bprec);

  mpf_clear2(x,y);
end;


{---------------------------------------------------------------------------}
procedure mpf_expt(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a^b, a>0}
var
  t: mp_float;
  p: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_expt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is_le0(a) then p := PINC
  else begin
    p := abs(s_mpf_ldx(a));
    if p<PINC then p := PINC;
  end;
  mpf_initp(t, c.bitprec+p);
  if mp_error=MP_OKAY then begin
    mpf_ln(a,t);
    mpf_mul(t,b,t);
    mpf_exp(t,c);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_expt1pm1(const a,b: mp_float; var c: mp_float);
  {-Calculate c = (1+a)^b-1, a>-1}
var
  t: mp_float;
  p: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_expt1pm1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is_le0(a) then p:=0
  else begin
    p := abs(s_mpf_ldx(a));
    if p<32 then p := PINC;
  end;
  mpf_initp(t, c.bitprec+p);
  if mp_error=MP_OKAY then begin
    mpf_ln1p(a,t);
    mpf_mul(t,b,t);
    mpf_expm1(t,c);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_exp10m1(const a: mp_float; var b: mp_float);
  {-Calculate b = 10^a-1}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_exp2m1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_initp(t, b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_set_ln10(t);
    mpf_mul(t,a,t);
    mpf_expm1(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_exp10i(n: longint; var a: mp_float);
  {-Calculate a = 10^n}
begin
  mpf_set_int(a,10);
  mpf_expt_int(a,n,a);
end;


{---------------------------------------------------------------------------}
procedure mpf_exp2m1(const a: mp_float; var b: mp_float);
  {-Calculate b = 2^a-1}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_exp2m1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_initp(t, b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_set_ln2(t);
    mpf_mul(t,a,t);
    mpf_expm1(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_exptm1(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a^b-1, a>0}
var
  t: mp_float;
  p: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_exptm1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is_le0(a) then p:=0
  else begin
    p := abs(s_mpf_ldx(a));
    if p<32 then p := PINC;
  end;
  mpf_initp(t, c.bitprec+p);
  if mp_error=MP_OKAY then begin
    mpf_ln(a,t);
    mpf_mul(t,b,t);
    mpf_expm1(t,c);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_expt_int(const a: mp_float; b: longint; var c: mp_float);
  {-Calculate c = a^b}
var
  x: mp_float;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_expt_int');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if b<0 then begin
    if b = -MaxLongint then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mpf_expt_int: b = -MaxLongint');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end
    else begin
      {c = (a^-1)^(-b)}
      {will give error if a=0}
      mpf_inv(a,c);
      mpf_expt_int(c, -b, c);
      exit;
    end
  end;

  {easy outs}
  if s_mpf_is0(a) then begin
    mpf_set0(c);
    exit;
  end;
  if b=0 then begin
    mpf_set1(c);
    exit;
  end
  else if b=1 then begin
    mpf_copyp(a,c);
    exit;
  end;

  {make local copy of a}
  mpf_initp(x, c.bitprec);
  if mp_error<>MP_OKAY then exit;

  mpf_copyp(a,x);
  mpf_set1(c);

  {initially b>=2}
  while mp_error=MP_OKAY do begin
    if odd(b) then mpf_mul(c,x,c);
    b := b shr 1;
    if b=0 then break else mpf_sqr(x,x);
  end;

  mpf_clear(x);
end;


{---------------------------------------------------------------------------}
procedure mpf_frac(const a: mp_float; var b: mp_float);
  {-Set b to the fractional part of a; frac(x)=x-int(x)}
begin
  s_mpf_frac(a,b,nil);
end;


{---------------------------------------------------------------------------}
function mpf_get_default_prec: longint;
  {-Return current default (bit) precision, initial=240}
begin
  mpf_get_default_prec := mpf_default_prec;
end;


{---------------------------------------------------------------------------}
procedure mpf_hav(const a: mp_float; var b: mp_float);
  {-Calculate the haversine b = hav(a) = 0.5*(1 - cos(a))}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_hav');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  mpf_initp(t, b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    {hav(a) = sqr(sin(0.5*a))}
    mpf_mul_2k(a,-1,t);
    mpf_sin(t,t);
    mpf_sqr(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_hypot(const a,b: mp_float; var c: mp_float);
  {-Calculate c = sqrt(a^2 + b^2)}
var
  q: mp_float;
begin
  if mpf_is0(a) then begin
    mpf_abs(b,c);
  end
  else if mpf_is0(b) then begin
    mpf_abs(a,c);
  end
  else begin
    mpf_initp(q, c.bitprec+8);
    if mp_error=0 then begin
      case mpf_cmp_mag(a,b) of
          0:  begin
                {|a| = |b|}
                mpf_set_dbl(q,2.0);
                mpf_sqrt(q,q);
                mpf_mul(q,a,q);
              end;
          1:  begin
                {|a| > |b|,  c = |a|*sqrt(1 + (b/a)^2)}
                mpf_div(b,a,q);
                mpf_sqr(q,q);
                s_mpf_inc1(q);
                mpf_sqrt(q,q);
                mpf_mul(q,a,q);
              end;
        else  begin
                {|a| < |b|,  c = |b|*sqrt(1 + (b/a)^2)}
                mpf_div(a,b,q);
                mpf_sqr(q,q);
                s_mpf_inc1(q);
                mpf_sqrt(q,q);
                mpf_mul(q,b,q);
              end;
      end;
      mpf_abs(q,c);
      mpf_clear(q);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_iexpt(a: longint; const b: mp_float; var c: mp_float);
  {-Calculate c = a^b, a>0}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(c);
  {$endif}
  {Essentially a copy of mpf_expt in order to avoid use of temporary mp_floats}
  mpf_initp(t, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_int(t,a);
    mpf_ln(t,t);
    mpf_mul(t,b,t);
    mpf_exp(t,c);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_init(var a: mp_float);
  {-Initialize an mp_float with default precision}
begin
  mpf_initp(a, mpf_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_init2(var a,b: mp_float);
  {-Initialize two mp_floats with default precision}
begin
  mpf_initp2(a,b,mpf_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_init3(var a,b,c: mp_float);
  {-Initialize 3 mp_floats with default precision}
begin
  mpf_initp3(a,b,c,mpf_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_init4(var a,b,c,d: mp_float);
  {-Initialize 4 mp_floats with default precision}
begin
  mpf_initp4(a,b,c,d,mpf_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_init5(var a,b,c,d,e: mp_float);
  {-Initialize 5 mp_floats with default precision}
begin
  mpf_initp5(a,b,c,d,e,mpf_default_prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_initp(var a: mp_float; prec: longint);
  {-Initialize an mp_float with bit precision prec}
begin
  if mp_error<>MP_OKAY then exit;
  if prec<MPF_MIN_PREC then prec := MPF_MIN_PREC;
  if prec>MPF_MAX_PREC then prec := MPF_MAX_PREC;
  with a do begin
    bitprec := prec;
    exponent:= 0;
    mp_init_size(mantissa, (prec+pred(DIGIT_BIT)) div DIGIT_BIT);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_initp2(var a,b: mp_float; prec: longint);
  {-Initialize two mp_floats with bit precision prec}
var
  pa: array[0..1] of pmp_float;
begin
  pa[0] := @a;
  pa[1] := @b;
  mpf_initp_multi_p(pa, prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_initp3(var a,b,c: mp_float; prec: longint);
  {-Initialize 3 mp_floats with bit precision prec}
var
  pa: array[0..2] of pmp_float;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  mpf_initp_multi_p(pa, prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_initp4(var a,b,c,d: mp_float; prec: longint);
  {-Initialize 4 mp_floats with bit precision prec}
var
  pa: array[0..3] of pmp_float;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  mpf_initp_multi_p(pa, prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_initp5(var a,b,c,d,e: mp_float; prec: longint);
  {-Initialize 5 mp_floats with bit precision prec}
var
  pa: array[0..4] of pmp_float;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  pa[4] := @e;
  mpf_initp_multi_p(pa, prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_initp_multi_p(var pv: array of pmp_float; prec: longint);
  {-Initialize with bit precision prec a list of mp_floats given as a pointer}
  { vector; on error the already initialized mp_floats will be cleared}
var
  i,k: integer;
begin
  if mp_error<>MP_OKAY then exit;
  for i:=low(pv) to high(pv) do begin
    mpf_initp(pv[i]^, prec);
    if mp_error<>MP_OKAY then begin
      {error, clear all previous mp_floats}
      for k:=low(pv) to i-1 do mpf_clear(pv[k]^);
      break;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_int(const a: mp_float; var b: mp_float);
  {-Set b to the integer part of a; i.e. is b rounded toward zero}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_int');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_copyp(a,b);
  with b do begin
    if exponent < 0 then begin
      mp_shr(mantissa,-exponent,mantissa);
      exponent := 0;
    end;
  end;
  s_mpf_normalize(b);
end;


{---------------------------------------------------------------------------}
procedure mpf_inv(const a: mp_float; var b: mp_float);
  {-Calculate b = 1/a}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_inv');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_init(x);
  if mp_error=MP_OKAY then begin
    mpf_set1(x);
    mpf_div(x,a,b);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
function mpf_is0(const a: mp_float): boolean;
  {-Return true if a=0}
begin
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  mpf_is0 := (a.exponent=0) and mp_is0(a.mantissa);
end;


{---------------------------------------------------------------------------}
function mpf_is1(const a: mp_float): boolean;
  {-Return true if a=1}
var
  n: longint;
begin
  {init check in mp_is_pow2}
  mpf_is1 := (a.mantissa.sign=MP_ZPOS) and mp_is_pow2(a.mantissa,n) and ((a.exponent+n)=0);
end;


{---------------------------------------------------------------------------}
function mpf_is1a(const a: mp_float): boolean;
  {-Return true if abs(a)=1}
var
  n: longint;
begin
  {init check in mp_is_pow2}
  mpf_is1a := mp_is_pow2(a.mantissa,n) and ((a.exponent+n)=0);
end;


{---------------------------------------------------------------------------}
function mpf_is_eq(const a,b: mp_float): boolean;
  {-Return a = b}
begin
  mpf_is_eq := mpf_cmp(a,b)=MP_EQ;
end;


{---------------------------------------------------------------------------}
function mpf_is_eq_rel(const a,b: mp_float): boolean;
  {-Check if |a-b| <= r*2^(1-b.bitprec);  r=1 if b=0, r=|b| otherwise}
var
  c,d: mp_float;
begin
  mpf_is_eq_rel := false;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_is_eq_rel');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_initp2(c,d,b.bitprec);
  if mp_error<>MP_OKAY then exit;
  if s_mpf_is0(b) then mpf_set1(d) else mpf_abs(b,d);
  mpf_mul_2k(d,1-b.bitprec,d);
  mpf_sub(a,b,c);
  mpf_is_eq_rel := mpf_cmp_mag(c,d)<>MP_GT;
  mpf_clear2(c,d);
end;


{---------------------------------------------------------------------------}
function mpf_is_ge(const a,b: mp_float): boolean;
  {-Return a >= b}
begin
  mpf_is_ge := mpf_cmp(a,b)<>MP_LT;
end;


{---------------------------------------------------------------------------}
function mpf_is_gt(const a,b: mp_float): boolean;
  {-Return a > b}
begin
  mpf_is_gt := mpf_cmp(a,b)=MP_GT;
end;


{---------------------------------------------------------------------------}
function mpf_is_le(const a,b: mp_float): boolean;
  {-Return a <= b}
begin
  mpf_is_le := mpf_cmp(a,b)<>MP_GT;
end;


{---------------------------------------------------------------------------}
function mpf_is_lt(const a,b: mp_float): boolean;
  {-Return a < b}
begin
  mpf_is_lt := mpf_cmp(a,b)=MP_LT;
end;


{---------------------------------------------------------------------------}
function mpf_is_ne(const a,b: mp_float): boolean;
  {-Return a <> b}
begin
  mpf_is_ne := mpf_cmp(a,b)<>MP_EQ;
end;


{---------------------------------------------------------------------------}
procedure s_mpf_lnagm(const a: mp_float; var b: mp_float);
  {-Calculate b = ln(a) with AGM algorithm, a>0; no init checks}
var
  x,t: mp_float;
  n: longint;
begin
  if mp_error<>MP_OKAY then exit;

  {ln(1) = 0}
  if mpf_is1(a) then begin
    mpf_set0(b);
    exit;
  end;

  {Use slightly larger working precision, see below}
  n := b.bitprec + PINC;
  if n>MPF_MAX_PREC then n:=MPF_MAX_PREC;

  mpf_initp2(x,t,n);
  if mp_error<>MP_OKAY then exit;

  {Brent/Zimmermann [35], Ch.4.8.2: First AGM algorithm for the logarithm}

  {Write x = 2^n*a with x > 2^(b.bitprec/2). For these x the AGM code can}
  {be used, then ln(a) = ln(x/2^n) = ln(x) - n*ln(2) with a loss of about}
  {log2(b.bitprec) bits. This is compensated by using b.bitprec + 32.}

  {Compute x = 2^n*a}
  n := (n div 2) - s_mpf_ldx(a);
  if n<0 then n:=0;
  mpf_mul_2k(a,n,x);

  {If x is large enough ln(x) = Pi/2/agm(1,4/x)}
  mpf_set_int(t,4);
  mpf_div(t,x,x);
  mpf_set1(t);
  mpf_agm(t,x,x);
  mpf_set_pi2k(t,-1);
  mpf_div(t,x,x);
  if n>0 then begin
    {subtract n*ln(2)}
    mpf_set_ln2(t);
    mpf_mul_int(t,n,t);
    mpf_sub(x,t,x);
  end;
  mpf_copyp(x,b);
  mpf_clear2(x,t);
end;


{---------------------------------------------------------------------------}
procedure mpf_ln(const a: mp_float; var b: mp_float);
  {-Calculate b = ln(a), a>0}
var
  x,y,t: mp_float;
  i,f,bprec: longint;
  inv: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_ln');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {a must be positive}
  if s_mpf_is_le0(a) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_ln: a <= 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {ln(1) = 0}
  if mpf_is1(a) then begin
    mpf_set0(b);
    exit;
  end;

  {save result precision and use slightly larger working precision}
  bprec := b.bitprec;
  if bprec>=mpf_lna_cutoff then begin
    s_mpf_lnagm(a,b);
    exit;
  end;

  i := bprec + PINC;
  if i>MPF_MAX_PREC then i:=MPF_MAX_PREC;

  mpf_initp3(x,y,t,i);
  if mp_error<>MP_OKAY then exit;

  {set working precision for b}
  s_mpf_normalizep(b,i);

  {if a < 1 then x=1/a else x=a}

  inv := s_mpf_ldx(a) <= 0;
  if inv then mpf_inv(a,x) else mpf_copyp(a,x);

  if s_mpf_is0(x) then begin
    {it could be possible that 1/a = 0, avoid endless loop}
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_ln: 1/a = 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {Range reduction 1: make x < 2}
  f := 0;
  while s_mpf_ldx(x) > 1 do begin
    inc(f);
    mpf_sqrt(x,x);
    if mp_error<>MP_OKAY then break;
  end;

  {y=x-1}
  mpf_set1(t);
  mpf_sub(x,t,y);

  if s_mpf_ldx(y) >= -7 then begin
    {Range reduction 2: make  1 < x < 1+2^-8}
    for i:=1 to 8 do begin
      inc(f);
      mpf_sqrt(x,x);
    end;
    {recalc y=x-1}
    mpf_sub(x,t,y);
  end;

  {ln(x) = 2*sum[ ((x-1)/(x+1))^(2i+1)/(2i+1) ] }
  mpf_add(x,t,x);
  {t will be the term ((x-1)/(x+1))^(2i+1)/(2i+1)}
  mpf_div(y,x,t);
  {y = ((x-1)/(x+1))^2}
  mpf_sqr(t,y);
  {set b to first term, loop starts with i=2}
  mpf_copyp(t,b);

  i := 1;
  repeat
    inc(i,2);
    if (i<0) or (mp_error<>MP_OKAY) then break;
    mpf_mul(t,y,t);
    if i <= lv_digit_max then mpf_div_d(t,mp_digit(i),x)
    else mpf_div_int(t,i,x);
  until not s_mpf_incf(b,x);

  {Undo range reduction and multiply by 2}
  s_mpf_incexp(b, 1+f);

  {ln(1/a) = - ln(a)}
  if inv then s_mpf_chs(b);

  {normalize to original precision}
  s_mpf_normalizep(b,bprec);

  mpf_clear3(x,y,t);
end;


{---------------------------------------------------------------------------}
procedure s_ln1p_ama(const a: mp_float; var b: mp_float);
  {-Calculate b = ln(1+a), a>-1; AMath algorithm, no init checks}
var
  x,y: mp_float;
  i,bprec: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {save result precision and use slightly larger working precision}
  bprec := b.bitprec;
  i := bprec + PINC;
  if i>MPF_MAX_PREC then i:=MPF_MAX_PREC;

  mpf_initp2(x,y,i);
  if mp_error<>MP_OKAY then exit;

  {Use mpf version of AMath algorithm for ln1p:   y := 1.0 + x; }
  { if y=1.0 then ln1p := x else ln1p := ln(y) + (x-(y-1.0))/y; }
  {This is faster than s_mpf_ln1p, except if |a| is smaller than}
  {a certain limit and the s_mpf_ln1p series converges rapidly. }
  mpf_copy(a,x);
  s_mpf_add1(x,y);
  mpf_copy(y,b);
  s_mpf_dec1(b);
  if s_mpf_is0(b) then mpf_copy(x,b)
  else begin
    mpf_sub(x,b,b); {b = (x-(y-1))}
    mpf_div(b,y,b); {b = (x-(y-1))/y)}
    mpf_ln(y,x);
    mpf_add(x,b,b); {b = ln(y) + (x-(y-1))/y)}
  end;
  {normalize to original precision}
  s_mpf_normalizep(b,bprec);
  mpf_clear2(x,y);
end;


{---------------------------------------------------------------------------}
procedure mpf_ln1p(const a: mp_float; var b: mp_float);
  {-Calculate b = ln(1+a), a>-1; special version for small a}
var
  x,y,t: mp_float;
  i,f,bprec: longint;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_ln1p');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  i := s_mpf_ldx(a);
  if i >= 0 then begin
    {if |a| > 0.5 then use standard function b = ln(1+a)}
    s_mpf_add1(a,b);
    mpf_ln(b,b);
    exit;
  end;

  bprec := b.bitprec;
  if i > (-0.005)*bprec then begin
    s_ln1p_ama(a,b);
    exit;
  end;

  {save result precision and use slightly larger working precision}
  i := bprec + PINC;
  if i>MPF_MAX_PREC then i:=MPF_MAX_PREC;

  mpf_initp3(x,y,t,i);
  if mp_error<>MP_OKAY then exit;

  {set working precision for b}
  s_mpf_normalizep(b,i);
  mpf_copyp(a,x);

  {range reduction: make x < 2^-10, using ln1p(x)=2*ln1p(x/(1+sqrt(1+x)))}
  f := 0;
  while s_mpf_ldx(x) >= -10 do begin
    s_mpf_add1(x,y);
    mpf_sqrt(y,y);
    s_mpf_inc1(y);
    mpf_div(x,y,x);
    inc(f);
  end;

  {ln(1+x) = 2y*[1 + t/3 + t^2/5 + t^3/7 ...],  y=x/(2+x), t=y^2}
  {b accumulates the result starting with b=x/(2+x)}
  mpf_set_int(b,2);
  s_mpf_inc(b,x);
  mpf_div(x,b,b);

  {t = y^2}
  mpf_sqr(b,t);

  {x is the term y*z^i}
  mpf_copy(b,x);

  {add the higher terms until no change in sum}
  i := 1;
  repeat
    inc(i,2);
    if (i<0) or (mp_error<>MP_OKAY) then break;
    mpf_mul(x,t,x);
    if i <= lv_digit_max then mpf_div_d(x,mp_digit(i),y)
    else mpf_div_int(x,i,y);
  until not s_mpf_incf(b,y);

  {Undo range reduction and multiply by 2}
  s_mpf_incexp(b, 1+f);

  {normalize to original precision}
  s_mpf_normalizep(b,bprec);

  mpf_clear3(x,y,t);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_logxp1(const a: mp_float; var b: mp_float; log2: boolean);
  {-Calculate b = logx(1+a), x=2 if log2 otherwise x=10}
var
  x,y: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mpf_logxp1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_initp2(x,y,b.bitprec + PINC);
  mpf_ln1p(a,x);
  if log2 then mpf_set_ln2(y)
  else mpf_set_ln10(y);
  mpf_div(x,y,b);
  mpf_clear2(x,y);
end;


{---------------------------------------------------------------------------}
procedure mpf_log10(const a: mp_float; var b: mp_float);
  {-Calculate b = log10(a), a>0}
var
  t: mp_float;
begin
  mpf_ln(a,b);
  if mpf_is0(b) or (mp_error<>MP_OKAY) then exit;
  mpf_initp(t, b.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_ln10(t);
    mpf_div(b,t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_log10p1(const a: mp_float; var b: mp_float);
  {-Calculate b = log10(1+a),  a > -1}
begin
  s_mpf_logxp1(a,b,false);
end;


{---------------------------------------------------------------------------}
procedure mpf_log2(const a: mp_float; var b: mp_float);
  {-Calculate b = log2(a), a>0}
var
  t: mp_float;
begin
  mpf_ln(a,b);
  if mpf_is0(b) or (mp_error<>MP_OKAY) then exit;
  mpf_initp(t, b.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_ln2(t);
    mpf_div(b,t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_log2p1(const a: mp_float; var b: mp_float);
  {-Calculate b = log2(1+a),  a > -1}
begin
  s_mpf_logxp1(a,b,true);
end;


{---------------------------------------------------------------------------}
procedure mpf_logbase(const b,x: mp_float; var y: mp_float);
  {-Calculate y = base b logarithm of x}
var
  u,v: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(y) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_logbase');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_initp2(u,v,y.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_ln(b,u);
    mpf_ln(x,v);
    mpf_div(v,u,y);
    mpf_clear2(u,v);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_mul(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a*b}
var
  be: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(c);
  {$endif}
  if mpf_is0(a) or mpf_is0(b) then mpf_set0(c)
  else begin
    mp_mul(a.mantissa, b.mantissa, c.mantissa);
    {need temp var for b.exponent, if @c = @a or @b}
    be := b.exponent;
    c.exponent := a.exponent;
    s_mpf_incexp(c, be);
    s_mpf_normalize(c);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_mul_2k(const a: mp_float; k: longint; var b: mp_float);
  {-Calculate b = a*2^k}
begin
  {arg check in mpf_copyp}
  mpf_copyp(a,b);
  s_mpf_incexp(b, k);
end;


{---------------------------------------------------------------------------}
procedure mpf_mul_d(const a: mp_float; b: mp_digit; var c: mp_float);
  {-Multiply by a digit}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
    _CheckBitPrec(c);
  {$endif}
  if (b=0) or mpf_is0(a) then mpf_set0(c)
  else if b=1 then mpf_copyp(a,c)
  else begin
    mp_mul_d(a.mantissa, b, c.mantissa);
    c.exponent := a. exponent;
    s_mpf_normalize(c);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_mul_dbl(const a: mp_float; b: double; var c: mp_float);
  {-Calculate c = a*b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_mul_dbl');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a) or (b=0.0) then begin
    mpf_set0(c);
    exit;
  end;
  mpf_initp(x,c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_dbl(x,b);
    mpf_mul(a,x,c);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_mul_int(const a: mp_float; b: longint; var c: mp_float);
  {-Multiply by a 32 bit integer}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
    _CheckBitPrec(c);
  {$endif}
  if mpf_is0(a) or (b=0) then mpf_set0(c)
  else if abs(b)=1 then begin
    if b<0 then mpf_chs(a,c) else mpf_copyp(a,c);
  end
  else begin
    mp_mul_int(a.mantissa, b, c.mantissa);
    c.exponent := a. exponent;
    s_mpf_normalize(c);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_mul_mpi(const a: mp_float; const b: mp_int; var c: mp_float);
  {-Calculate c = a*b}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
    _CheckBitPrec(c);
  {$endif}
  if mpf_is0(a) or mp_is0(b) then mpf_set0(c)
  else begin
    mp_mul(a.mantissa, b, c.mantissa);
    c.exponent := a.exponent;
    s_mpf_normalize(c);
  end;
end;


{---------------------------------------------------------------------------}
function mpf_not_init(const a: mp_float): boolean;
  {-Sanity check if a is initialized, does not catch all cases!}
begin
  with a do begin
    mpf_not_init := mp_not_init(mantissa) or (bitprec<MPF_MIN_PREC) or (bitprec>MPF_MAX_PREC);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_nroot(const a: mp_float; n: longint; var b: mp_float);
  {-Calculate the nth root of a: b = a^(1/n); n<>0, a >= 0  if n is even}
var
  k: longint;
  s: word;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_nroot');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  k := abs(n);
  if k=1 then mpf_copyp(a,b)
  else if k=2 then mpf_sqrt(a,b)
  else begin
    s := a.mantissa.sign;
    if odd(k) then mpf_abs(a,b)
    else mpf_copyp(a,b);
    mpf_ln(b,b);
    mpf_div_int(b,k,b);
    mpf_exp(b,b);
    if s=MP_NEG then s_mpf_chs(b);
  end;
  if (mp_error=MP_OKAY) and (n<0) then mpf_inv(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_numbpart(n: longint; var p: mp_int);
  {-Calculate the number of partitions of n with Hardy-Ramanujan-Rademacher formula}
const
  NS = 16;
  small: array[0..NS] of byte = (1,1,2,3,5,7,11,15,22,30,42,56,77,101,135,176,231);
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(p) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_numbpart');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if n<=NS then begin
    if n>=0 then mp_set_int(p,small[n]) else mp_zero(p);
    exit;
  end;
  mpf_init(x);
  if mp_error=MP_OKAY then begin
    s_mpf_numbpart(n,x);
    mpf_round(x,p);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_output_decimal(const a: mp_float; ndd: word);
  {-Write an mp_float to output using decimal scientific representation, ndd>0}
  { is the total number of digits (including one digit before the '.')}
begin
  mpf_write_decimal(output, a, ndd);
end;


{---------------------------------------------------------------------------}
procedure mpf_output_decimal_alt(const a: mp_float; ndd: word);
  {-Write an mp_float to output using decimal alternative representation, ndd>0 is}
  { the total number of digits, trailing '.' or '0' are suppressed. If a is too}
  { large or too small, scientific representation is used.}
begin
  mpf_write_decimal_alt(output, a, ndd);
end;


{---------------------------------------------------------------------------}
procedure mpf_output_radix(const a: mp_float; radix,ndd: word);
  {-Write an mp_float to output using radix scientific representation, ndd>0}
  { is the total number of digits (including one digit before the '.')}
begin
  mpf_write_radix(output, a, radix, ndd);
end;


{---------------------------------------------------------------------------}
procedure mpf_output_radix_alt(const a: mp_float; radix,ndd: word);
  {-Write an mp_float to output using radix alternative representation, ndd>0 is}
  { the total number of digits, trailing '.' or '0' are suppressed. If a is too}
  { large or too small, scientific representation is used. NOTE: no radix prefix/suffix!}
begin
  mpf_write_radix_alt(output, a, radix, ndd);
end;


{---------------------------------------------------------------------------}
procedure mpf_rad2deg(const a: mp_float; var b: mp_float);
  {-Convert a in radians to b in degrees}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  mpf_initp(t,b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_set_pi(t);
    mpf_div(a,t,t);
    mpf_mul_d(t,180,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_random(var a: mp_float);
  {-Set to a random number uniformly distributed in [0,1)}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  with a do begin
    mp_rand_bits_ex(mantissa, bitprec, false);
    exponent := -bitprec;
  end;
  s_mpf_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpf_read_decimal(var a: mp_float; str: pchar8);
  {-Read a from ASCII float decimal string. str may contain a single '.'}
  { The exponent part is @[+|-]nnn (e or E can replace @). Integer or }
  { fractional part must be present.}
begin
  mpf_read_radix(a, str, 10);
end;


{---------------------------------------------------------------------------}
procedure mpf_read_hex(var a: mp_float; str: pchar8);
  {-Read a from ASCII float hexadecimal string. str may contain a single}
  { '.'. The exponent part is @[+|-]nnn (h or H can replace @). Integer }
  { or fractional part must be present.}
begin
  mpf_read_radix(a, str, 16);
end;


{---------------------------------------------------------------------------}
procedure mpf_read_radix(var a: mp_float; str: pchar8; radix: word);
  {-Read a from an ASCII float radix string. str may contain a single '.'.}
  { The exponent part <xp> is @[+|-]nnn, b/B, e/E, or h/H can replace @ if}
  { radix in [2,10,16]. The integer <ip> or fractional part <fp> must be }
  { present, nnn must be decimal, ie general format is <ip>.<fp>*radix^<xp>.}
var
  fs: pchar8;
  ip,fp: mp_int;
  x: mp_float;
  lf,xp,li: longint;
  {$ifdef VirtualPascal}
    ec: longint;
  {$else}
    ec: integer;
  {$endif}
  c: char8;
  neg: boolean;
  sep,dsep: string[5];

{$ifndef MPC_HaltOnError}
label leave;
{$endif}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_read_radix');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_init2(ip,fp);
  if mp_error<>MP_OKAY then exit;

  mpf_initp(x,a.bitprec);
  if mp_error<>MP_OKAY then begin
    mp_clear2(ip,fp);
    exit;
  end;

  if radix=2 then sep := '@bB'
  else if radix=10 then sep := '@eE'
  else if radix=16 then sep := '@hH'
  else sep := '@';
  dsep := mp_fract_sep{'.'}+sep;
  lf := 0;
  li := 0;
  xp := 0;

  {skip leading white space}
  repeat
    c := upcase(str^);
    if c=#0 then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mpf_read_radix: invalid syntax');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        goto leave;
      {$endif}
    end;
    if (c<>' ') and (c<>#9) then break;
    inc(str);
  until false;
  neg := false;

  if c='-' then begin
    neg := true;
    inc(str);
  end
  else if c='+' then inc(str);
  if pos(str^,dsep)>0 then mp_zero(ip)
  else begin
    {read integer part until '.@eE' or #0}
    fs := str;
    s_mp_read_radix(ip,str,radix,dsep,false);
    li := str-fs;
  end;
  if str^=mp_fract_sep{'.'} then begin
    {read fractional part}
    inc(str);
    if pos(str^,sep+#0)=0 then begin
      {remember starting position of floating part}
      fs := str;
      s_mp_read_radix(fp,str,radix,sep,false);
      {fractional part = fp/radix^lf}
      lf := str-fs;
    end;
  end;
  if str^<>#0 then begin
    {exponent part}
    inc(str);
    {$ifdef D12Plus}
      {Make the implicit type cast explicit to avoid warning}
      val(string(str),xp,ec);
    {$else}
      val(str,xp,ec);
    {$endif}
    if ec<>0 then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mpf_read_radix: invalid syntax');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        goto leave;
      {$endif}
    end;
  end;
  if lf>0 then begin
    {fractional part present, calculate a = ip + fp/radix^lf}
    (*
      {Old code gives 0.5 <> 5/10!}
      mpf_set_int(a,radix);
      mpf_expt_int(a,-lf,a);
      mpf_mul_mpi(a,fp,a);
      mpf_add_mpi(a,ip,a);
    *)
    mpf_set_mpi(a,fp);
    mpf_set_int(x,radix);
    mpf_expt_int(x,lf,x);
    mpf_div(a,x,a);
    mpf_add_mpi(a,ip,a);
  end
  else begin
    if li=0 then begin
      {no integer and fractional part}
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mpf_read_radix: invalid syntax');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        goto leave;
      {$endif}
    end
    else begin
      {no fractional part}
      mpf_set_mpi(a,ip);
    end;
  end;
  if xp<>0 then begin
    {multiply by radix^xp}
    mpf_set_int(x,radix);
    mpf_expt_int(x,xp,x);
    mpf_mul(x,a,a);
  end;
  if neg then s_mpf_chs(a);
{$ifndef MPC_HaltOnError}
leave:
{$endif}
  mp_clear2(ip,fp);
  mpf_clear(x);
end;


{---------------------------------------------------------------------------}
function mpf_reldev(const a,b: mp_float): double;
  {-Return abs((a-b)/b)*2^b.bitprec, special if b=0, or a-b=0}
var
  c: mp_float;
begin
  mpf_reldev := DblPosInf;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_reldev');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_initp(c,b.bitprec);
  if mp_error<>MP_OKAY then exit;
  mpf_sub(a,b,c);
  if s_mpf_is0(c) then mpf_reldev := 0.0
  else begin
    if not s_mpf_is0(b) then mpf_div(c,b,c);
    mpf_mul_2k(c,b.bitprec,c);
    mpf_reldev := abs(mpf_todouble(c));
  end;
  mpf_clear(c);
end;


{---------------------------------------------------------------------------}
procedure mpf_rem_2pi(const a: mp_float; var b: mp_float);
  {-Calculate b = a mod (2*Pi)}
var
  oddm: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  s_mpf_mod_pi2k(a,1,b,oddm);
end;


{---------------------------------------------------------------------------}
procedure mpf_round(const a: mp_float; var b: mp_int);
  {-Round an mp_float to nearest mp_int}
var
  msign: word;
  rbset: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_round');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  with a do begin
    if exponent>=0 then mp_shl(mantissa,exponent,b)
    else begin
      {Use variables to fix quirk if @b = @a.mantissa}
      {Get rounding bit of mantissa}
      rbset := mp_isbit(mantissa, -exponent-1);
      msign := mantissa.sign;
      mp_shr(mantissa,-exponent,b);
      if rbset then begin
        if msign=MP_NEG then mp_dec(b) else mp_inc(b);
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_set0(var a: mp_float);
  {-Set a=0, a.bitprec is preserved}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  mp_zero(a.mantissa);
  if mp_error=MP_OKAY then a.exponent:= 0;
end;


{---------------------------------------------------------------------------}
procedure mpf_set1(var a: mp_float);
  {-Set a=1, a.bitprec is preserved}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  with a do begin
    mp_set(mantissa,1);
    exponent := 0;
  end;
  s_mpf_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_default_prec(prec: longint);
  {-Set new default (bit) precision}
begin
  if prec<MPF_MIN_PREC then prec := MPF_MIN_PREC;
  if prec>MPF_MAX_PREC then prec := MPF_MAX_PREC;
  mpf_default_prec := prec;
end;


{---------------------------------------------------------------------------}
procedure mpf_set_default_decprec(dprec: word);
  {-Set default precision dprec*log_2(10), i.e. about dprec decimal digits}
begin
  mpf_set_default_prec(trunc(1.0 + 3.321928095*dprec));
end;


{---------------------------------------------------------------------------}
procedure mpf_set_exp1(var a: mp_float);
  {-Set a to exp(1), preserve a.bitprec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_set_exp1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {$ifdef MPC_E1Ln10Tab}
    _set_ptab_const(a, 2, a.bitprec, AddrEBytes);
  {$else}
    mpf_set_int(a,1);
    mpf_exp(a,a);
  {$endif}
end;


{---------------------------------------------------------------------------}
procedure mpf_set_exp1p(var a: mp_float; prec: longint);
  {-Set a to exp(1) with bit precision prec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_set_exp1p');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {$ifdef MPC_E1Ln10Tab}
    _set_ptab_const(a, 2, prec, AddrEBytes);
  {$else}
    a.bitprec := prec;
    mpf_set_int(a,1);
    mpf_exp(a,a);
  {$endif}
end;


{---------------------------------------------------------------------------}
procedure mpf_set_dbl(var a: mp_float; b: double);
  {-Set a to a double. Error if b = NAN or INF}
var
  i,e: integer;
  d: double;
  n: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_set_dbl');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  d := b;
  e := TMPHexDblW(d)[3] and $7FF0;
  n := TMPHexDblW(d)[3] and $8000 <> 0;

  if e=$7FF0 then begin
    {nan or inf}
    {$ifdef MPC_HaltOnError}
     {$ifdef MPC_UseExceptions}
       raise MPXBadArg.Create('mpf_set_dbl: NAN or INF');
     {$else}
       RunError(MP_RTE_BADARG);
     {$endif}
   {$else}
     set_mp_error(MP_BADARG);
     exit;
   {$endif}
  end;

  {zerofill sign and exp}
  TMPHexDblW(d)[3] := TMPHexDblW(d)[3] and $000F;
  if e<>0 then begin
    {Remove bias from exponent and adjust with mantissa length}
    e := (e shr 4) - (1022+53);
    {Insert hidden bit}
    TMPHexDblW(d)[3] := TMPHexDblW(d)[3] or $0010;
  end
  else begin
    {denormal, mpf exponent is const}
    e := -(1022+52);
  end;
  mp_set(a.mantissa, TMPHexDblA(d)[6]);
  for i:=5 downto 0 do begin
    mp_shl(a.mantissa, 8, a.mantissa);
    mp_add_d(a.mantissa, TMPHexDblA(d)[i], a.mantissa);
  end;
  if n then s_mpf_chs(a);
  a.exponent := e;
  s_mpf_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_int(var a: mp_float; b: longint);
  {-Set a to a longint}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  with a do begin
    mp_set_int(mantissa,b);
    exponent := 0;
  end;
  s_mpf_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_ln10(var a: mp_float);
  {-Set a to ln(10), preserve a.bitprec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  mpf_set_ln10p2k(a,0,a.bitprec);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_ln10p(var a: mp_float; prec: longint);
  {-Set a to ln(10) with bit precision prec}
begin
  mpf_set_ln10p2k(a,0,prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_ln10p2k(var a: mp_float; k,prec: longint);
  {-Set a to ln(10)*2^k with bit precision prec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_set_ln10p2k');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {$ifdef MPC_E1Ln10Tab}
    _set_ptab_const(a, 2, prec, AddrLn10Bytes);
  {$else}
    a.bitprec := prec;
    mpf_set_int(a,10);
    mpf_ln(a,a);
  {$endif}
  if k<>0 then s_mpf_incexp(a, k);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_ln2(var a: mp_float);
  {-Set a to ln(2), preserve a.bitprec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  mpf_set_ln2p2k(a,0,a.bitprec);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_ln2p(var a: mp_float; prec: longint);
  {-Set a to ln(2) with bit precision prec}
begin
  mpf_set_ln2p2k(a,0,prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_ln2p2k(var a: mp_float; k,prec: longint);
  {-Set a to ln(2)*2^k with bit precision prec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_set_ln2p2k');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  _set_ptab_const(a, 0, prec, AddrLn2Bytes);
  if k<>0 then s_mpf_incexp(a, k);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_mpi(var a: mp_float; const b: mp_int);
  {-Set a to an mp_int}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  with a do begin
    mp_copy(b, mantissa);
    exponent := 0;
  end;
  s_mpf_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_mpi2k(var a: mp_float; const m: mp_int; e: longint);
  {-Set a to m*2^e (build a from mantissa and exponent)}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  with a do begin
    mp_copy(m, mantissa);
    exponent := e;
  end;
  s_mpf_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_pi(var a: mp_float);
  {-Set a to pi, preserve a.bitprec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  mpf_set_pip2k(a,0,a.bitprec);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_pip(var a: mp_float; prec: longint);
  {-Set a to pi with bit precision prec}
begin
  mpf_set_pip2k(a,0,prec);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_pip2k(var a: mp_float; k,prec: longint);
  {-Set a to pi*2^k with bit precision prec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_set_pip2k');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  _set_ptab_const(a, 2, prec, AddrPiBytes);
  if k<>0 then s_mpf_incexp(a, k);
end;


{---------------------------------------------------------------------------}
procedure mpf_set_pi2k(var a: mp_float; k: longint);
  {-Set a to pi*2^k, preserve a.bitprec}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
  {$endif}
  mpf_set_pip2k(a,k,a.bitprec);
end;


{---------------------------------------------------------------------------}
procedure mpf_sin(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(a)}
begin
  mpf_trig(a, nil, @b, nil);
end;


{---------------------------------------------------------------------------}
procedure mpf_sinc(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(a)/a}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sinc');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mpf_is0(a) or (s_mpf_ldx(a) < -(b.bitprec + PINC) div 2) then begin
    {sinc(a) = 1 - a^2/6 + O(a^4))}
    mpf_set1(b);
    exit;
  end;
  mpf_initp(t, b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_sin(a,t);
    mpf_div(t,a,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_sincos(const a: mp_float; var b,c: mp_float);
  {-Calculate b = sin(a), c = cos(a);  @b<>@c}
begin
  mpf_trig(a,@c,@b,nil);
end;


{---------------------------------------------------------------------------}
procedure mpf_sind(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(a), a in degrees}
begin
  mpf_deg2rad(a,b);
  mpf_sin(b,b);
end;

{---------------------------------------------------------------------------}
procedure mpf_sinh(const a: mp_float; var b: mp_float);
  {-Calculate b = sinh(a), |a| < 2^31 * ln(2)}
var
  t: mp_float;
  pb: longint;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sinh');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  pb := b.bitprec;
  mpf_initp(t, pb + PINC);
  if mp_error<>MP_OKAY then exit;

  neg := s_mpf_is_neg(a);
  mpf_abs(a,t);
  b.bitprec := t.bitprec;

  if s_mpf_ldx(t) >= 0 then begin
    {t>=0.5, sinh(t) = (exp(t)-exp(-t))/2}
    mpf_exp(t,b);
    {don't calculate inverse if it is to small}
    if s_mpf_ldx(b) < 16 + b.bitprec div 2 then begin
      mpf_inv(b,t);
      mpf_sub(b,t,b);
    end;
  end
  else begin
    {t<0.5, sinh(t) = 0.5*expm1(t)*(1+exp(-t))}
    {sinh(t) = 0.5*y(1+1/(1+y)) with y = expm1(t)}
    mpf_expm1(t,b);
    s_mpf_add1(b,t);
    mpf_inv(t,t);
    s_mpf_inc1(t);
    mpf_mul(b,t,b);
  end;
  s_mpf_incexp(b,-1);
  s_mpf_normalizep(b,pb);
  if neg then s_mpf_chs(b);
  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpf_sinhcosh(const a: mp_float; var b,c: mp_float);
  {-Calculate b = sinh(a), c = cosh(a);  @b<>@c,  |a| < 2^31 * ln(2)}
var
  t: mp_float;
  pb,pc,prec: longint;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sinhcosh');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if @b=@c then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mpf_sinhcosh: @b=@c');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    mpf_set1(c);
    exit;
  end;

  pb := b.bitprec;
  pc := b.bitprec;
  if pb>pc then prec := pb else prec := pc;

  mpf_initp(t, prec + PINC);
  if mp_error<>MP_OKAY then exit;

  neg := s_mpf_is_neg(a);
  mpf_abs(a,t);

  b.bitprec := t.bitprec;
  c.bitprec := t.bitprec;

  if s_mpf_ldx(t) >= 0 then begin
    {t>=0.5, sinh(t) = (exp(t)-exp(-t))/2, cosh(t) = (exp(t)+exp(-t))/2}
    mpf_exp(t,b);
    {don't calculate inverse if it is to small}
    if s_mpf_ldx(b) < 16 + b.bitprec div 2 then begin
      mpf_inv(b,t);
      mpf_add(b,t,c);
      mpf_sub(b,t,b);
    end
    else mpf_copy(b,c);
  end
  else begin
    {t<0.5, sinh(t) = 0.5*expm1(t)*(1+exp(-t))}
    {sinh(t) = 0.5*y(1+1/(1+y)) with y = expm1(t)}
    {cosh(t) = 0.5*((1+y)+1/(1+y))}
    mpf_expm1(t,b);
    s_mpf_add1(b,c);
    {b=y, c=(1+y)}
    mpf_inv(c,t);
    mpf_add(c,t,c);
    s_mpf_inc1(t);
    mpf_mul(b,t,b);
  end;
  s_mpf_incexp(b,-1);
  s_mpf_incexp(c,-1);
  s_mpf_normalizep(b,pb);
  s_mpf_normalizep(c,pc);
  if neg then s_mpf_chs(b);
  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpf_sinpi(const a: mp_float; var b: mp_float);
  {-Calculate b = sin(Pi*a)}
begin
  mpf_trig_ex(a, true, nil, @b, nil);
end;


{---------------------------------------------------------------------------}
procedure mpf_sqr(const a: mp_float; var b: mp_float);
  {-Compute b = a*a}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    _CheckBitPrec(a);
    _CheckBitPrec(b);
  {$endif}
  mp_sqr(a.mantissa, b.mantissa);
  b.exponent := a.exponent;
  s_mpf_incexp(b,b.exponent);
  s_mpf_normalize(b);
end;


{---------------------------------------------------------------------------}
procedure mpf_sqrt(const a: mp_float; var b: mp_float);
  {-Calculate b = sqrt(a)}
var
  diff: longint;
begin
  if mp_error<>MP_OKAY then exit;

  {a must be positive}
  if s_mpf_is_neg(a) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_sqrt: a < 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {other arg check in mpf_copyp}
  {copy const a to working float b and normalize to b.bitprec}
  mpf_copyp(a,b);
  if s_mpf_is0(b) or mpf_is1(b) then exit;

  with b do begin
    {use bitprec+2 to allow more reliable rounding in normalize}
    diff := 2*(bitprec+2) - mp_bitsize(mantissa);
    if diff>=0 then begin
      {shift mantissa and make sure exponent is even}
      if odd(exponent xor diff) then inc(diff);
      s_mpf_incexp(b, -diff);
      mp_shl(mantissa, diff, mantissa);
    end
    else begin
      {This branch SHOULD not be taken because copyp normalizes to}
      {precision b.bitprec. Just paranoia, but make exponent even!}
      if odd(exponent) then begin
        s_mpf_incexp(b, 1);
        mp_shr(mantissa, 1, mantissa);
      end;
    end;
    mp_sqrt(mantissa, mantissa);
    exponent := exponent div 2;
  end;
  s_mpf_normalize(b);
end;


{---------------------------------------------------------------------------}
procedure mpf_sqrt1pm1(const a: mp_float; var b: mp_float);
  {-Calculate b = sqrt(1+a)-1 with increased accuracy for a near 0, a >= -1}
var
  t: mp_float;
begin
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sqrt1pm1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_initp(t, b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    s_mpf_add1(a,t);
    mpf_sqrt(t,t);
    if s_mpf_ldx(a) >= 0 then begin
      {a >= 0.5}
      s_mpf_dec1(t);
      mpf_copyp(t,b);
    end
    else begin
      {sqrt(1+x)-1 = x/(sqrt(1+x)+1)}
      s_mpf_add1(t,t);
      mpf_div(a,t,b);
    end;
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
function mpf_squad(const a,b,c: mp_float; var x1,y1,x2,y2: mp_float): integer;
  {-Solve the quadratic equation a*x^2 + b*x + c = 0. Result is the number}
  { of different solutions: 0 (if a=b=0), 1 (x1), or 2 (x1 and x2). If the}
  { result is = -2, x1 + i*y1 and x2 + i*y2 are the two complex solutions.}
  { x1,y1,x2,y2 should be different variables and not the same as a,b or c}
var
  b2,d,t: mp_float;
  p,pd : longint;
begin

  mpf_squad := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) or
       mpf_not_init(x1) or mpf_not_init(x2) or mpf_not_init(y1) or mpf_not_init(y2) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_squad');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Setup default return values}
  mpf_squad := 2;
  mpf_set0(x2);
  mpf_set0(y1);
  mpf_set0(y2);

  if mpf_is0(a) then begin
    {solve bx+c=0}
    if mpf_is0(b) then mpf_squad := 0
    else begin
      mpf_squad := 1;
      mpf_div(c,b,x1);
      s_mpf_chs(x1);
    end;
    exit;
  end;

  if mpf_is0(c) then begin
    {ax^2+bx = 0 = x(ax+b) with a<>0}
    {x1 := -b/a;}
    mpf_div(b,a,x1);
    s_mpf_chs(x1);
    if s_mpf_is_ge0(x1) then mpf_exch(x1,x2);
    exit;
  end;

  {Here a<>0 and c<>0. Set precision for discriminant}
  {calculation: maximum precision of b^2, 4ac, and x1}
  pd := 2*b.bitprec;
  p  := 2+a.bitprec+c.bitprec; if p>pd then pd := p;
  p  := x1.bitprec;            if p>pd then pd := p;

  mpf_initp(b2,b.bitprec+8);
  if mp_error<>MP_OKAY then exit;
  mpf_initp2(d,t,pd+8);
  if mp_error<>MP_OKAY then begin
    mpf_clear(b2);
    exit;
  end;

  {get b2 = -b/2}
  mpf_mul_2k(b,-1,b2);
  s_mpf_chs(b2);

  {d = discriminant(a,b2,c) = b2*b2 - a*c}
  mpf_sqr(b2,d);
  mpf_mul(a,c,t);
  mpf_sub(d,t,d);

  if mpf_is0(d) then begin
    {d=0, real double root}
    mpf_squad := 1;
    mpf_div(b2,a,x1);
    mpf_copy(x1,x2);
    exit;
  end
  else if s_mpf_is_neg(d) then begin
    {d<0: two complex roots}
    mpf_squad := -2;
    s_mpf_chs(d);
    mpf_sqrt(d,d);
    mpf_div(b2,a,x1);
    mpf_copy(x1,x2);
    mpf_div(d,a,y1);
    mpf_chs(y1,y2);
    if mpf_is_gt(y1,y2) then mpf_exch(y1,y2);
  end
  else begin
    {d>0: two real roots}
    mpf_sqrt(d,d);
    if s_mpf_is_ge0(b2) then mpf_add(b2,d,d)
    else mpf_sub(b2,d,d);
    mpf_div(d,a,x1);
    mpf_div(c,d,x2);
    if mpf_is_gt(x1,x2) then mpf_exch(x1,x2);
  end;
  mpf_clear3(b2,d,t);
end;


{---------------------------------------------------------------------------}
procedure mpf_sub(const a,b: mp_float; var c: mp_float);
  {-Calculate c = a-b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sub');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(b) then begin
    mpf_copyp(a,c);
    exit;
  end;
  if s_mpf_is0(a) then begin
    mpf_chs(b,c);
    exit;
  end;
  if @a=@b then begin
    mpf_set0(c);
    exit;
  end;
  {make local copy of a}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_chs(b,x);
    {x = a-b = -a + b}
    s_mpf_inc(x,a);
    mpf_exch(x,c);  {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_sub_dbl(const a: mp_float; b: double; var c: mp_float);
  {-Calculate c = a-b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sub_dbl');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if b=0.0 then begin
    mpf_copyp(a,c);
    exit;
  end;
  if s_mpf_is0(a) then begin
    {c = -b}
    mpf_set_dbl(c,-b);
    exit;
  end;
  {make local copy of a}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_dbl(x,-b);
    {c = (-b)+a}
    s_mpf_inc(x,a);
    mpf_exch(x,c);    {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_sub_int(const a: mp_float; b: longint; var c: mp_float);
  {-Calculate c = a-b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sub_int');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if b=0 then begin
    mpf_copyp(a,c);
    exit;
  end;
  if s_mpf_is0(a) then begin
    {c = -b}
    mpf_set_int(c,-b);
    exit;
  end;
  {make local copy of a}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_int(x,-b);
    {c = (-b)+a}
    s_mpf_inc(x,a);
    mpf_exch(x,c);    {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_sub_mpi(const a: mp_float; const b: mp_int; var c: mp_float);
  {-Calculate c = a-b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mp_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sub_mpi');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mp_is0(b) then begin
    mpf_copyp(a,c);
    exit;
  end;
  if s_mpf_is0(a) then begin
    mpf_set_mpi(c,b);
    s_mpf_chs(c);
    exit;
  end;
  {make local copy of a}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_mpi(x,b);
    s_mpf_chs(x);
    {c = (-b)+a}
    s_mpf_inc(x,a);
    mpf_exch(x,c);     {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_sumalt(term: sumalt_term_func; n: longint; var s: mp_float; var Err: integer);
  {-Calculate s=sum(i=0..n, term(i)) of alternating series with sumalt algorithm.}
  { If n<4, n is adjusted to bitprec of s. Err<>0 if a term cannot be evaluated.}
var
  b,c,t: mp_float;
  k,bprec,xn,xd: longint;
begin
  {H. Cohen, F. Rodriguez Villegas, D. Zagier: Convergence acceleration of   }
  {alternating series, Experimental Mathematics, vol.9, no.1, pp.3-12, 2000  }
  {available via http://projecteuclid.org/download/pdf_1/euclid.em/1046889587}

  {J. Arndt[31], Chapter: Numerical evaluation of power series}
  {Section: The magic sumalt algorithm for alternating series}

  Err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(s) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sumalt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {save result precision and use slightly larger working precision}
  bprec := s.bitprec;
  k := bprec + PINC;
  if k>MPF_MAX_PREC then k:=MPF_MAX_PREC;
  mpf_initp3(b,c,t,k);
  if mp_error<>MP_OKAY then exit;

  if n<4 then n := trunc(0.39321985067869744234*k); {ln(2)/ln(3 + sqrt(8)))}

  mpf_set1(b);
  mpf_mul_2k(b, 2*n-1, b);
  mpf_copy(b,c);

  {set working precision for s}
  mpf_set0(s);
  s_mpf_normalizep(s,k);

  for k:=n-1 downto 0 do begin
     {get x[k] = xn/xd}
     if not term(k,xn,xd) then begin
       Err := -1;
       break;
     end;
     {s := s + c * x[k]}
     mpf_div_int(c,xd,t);
     mpf_mul_int(t,xn,t);
     mpf_add(s,t,s);
     {b := b * ((2*k+1)*(k+1)) / (2*(n+k)*(n-k))}
     mpf_div_int(b, 2*(n+k)*(n-k), b);
     mpf_mul_int(b, (2*k+1)*(k+1), b);
     {c := c + b}
     mpf_add(c,b,c);
  end;
  if Err=0 then begin
    if mpf_is0(c) then Err := -2
    else begin
      mpf_div(s,c,s);
      {normalize to original precision}
      s_mpf_normalizep(s,bprec);
    end;
  end;
  mpf_clear3(b,c,t);
end;


{---------------------------------------------------------------------------}
procedure mpf_sumaltf(fterm: sumalt_fterm_func; n: longint; var s: mp_float; var Err: integer);
  {-Calculate s=sum(i=0..n, term(i)) of alternating series with sumalt algorithm.}
  { If n<4, n is adjusted to bitprec of s. Err<>0 if a term cannot be evaluated.}
var
  b,c,t: mp_float;
  k,bprec: longint;
begin
  {H. Cohen, F. Rodriguez Villegas, D. Zagier: Convergence acceleration of   }
  {alternating series, Experimental Mathematics, vol.9, no.1, pp.3-12, 2000  }
  {available via http://projecteuclid.org/download/pdf_1/euclid.em/1046889587}

  {J. Arndt[31], Chapter: Numerical evaluation of power series}
  {Section: The magic sumalt algorithm for alternating series}

  Err := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(s) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sumaltf');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {save result precision and use slightly larger working precision}
  bprec := s.bitprec;
  k := bprec + PINC;
  if k>MPF_MAX_PREC then k:=MPF_MAX_PREC;
  mpf_initp3(b,c,t,k);
  if mp_error<>MP_OKAY then exit;

  if n<4 then n := trunc(0.39321985067869744234*k); {ln(2)/ln(3 + sqrt(8)))}

  mpf_set1(b);
  mpf_mul_2k(b, 2*n-1, b);
  mpf_copy(b,c);

  {set working precision for s}
  mpf_set0(s);
  s_mpf_normalizep(s,k);

  for k:=n-1 downto 0 do begin
     {get k-th term x[k]}
     if not fterm(k,t) then begin
       Err := -1;
       break;
     end;
     {s := s + c * x[k]}
     mpf_mul(t,c,t);
     mpf_add(s,t,s);
     {b := b * ((2*k+1)*(k+1)) / (2*(n+k)*(n-k))}
     mpf_div_int(b, 2*(n+k)*(n-k), b);
     mpf_mul_int(b, (2*k+1)*(k+1), b);
     {c := c + b}
     mpf_add(c,b,c);
  end;
  if Err=0 then begin
    if mpf_is0(c) then Err := -2
    else begin
      mpf_div(s,c,s);
      {normalize to original precision}
      s_mpf_normalizep(s,bprec);
    end;
  end;
  mpf_clear3(b,c,t);
end;


{---------------------------------------------------------------------------}
procedure mpf_tan(const a: mp_float; var b: mp_float);
  {-Calculate b = tan(a)}
begin
  mpf_trig(a, nil, nil, @b);
end;


{---------------------------------------------------------------------------}
procedure mpf_tand(const a: mp_float; var b: mp_float);
  {-Calculate b = tan(a), a in degrees}
begin
  mpf_deg2rad(a,b);
  mpf_tan(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_tanh(const a: mp_float; var b: mp_float);
  {-Calculate b = tanh(a)}
var
  s,t: mp_float;
  pos: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_tanh');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  pos := s_mpf_is_ge0(a);
  if s_mpf_ldx(a) > 1.45*ln(b.bitprec+2.0) then begin
    {exp(-x) < 2^(-b.bitprec)}
    if pos then mpf_set1(b)
    else mpf_set_int(b,-1);
  end
  else begin
    mpf_initp2(s,t, b.bitprec + PINC);
    if mp_error<>MP_OKAY then exit;
    mpf_abs(a,t);
    {tanh(a) = -sign(a)*t/(t+2) with t = expm1(-2|a|)}
    s_mpf_chs(t);
    s_mpf_incexp(t,1);
    mpf_expm1(t,s);
    mpf_set_int(t,2);
    s_mpf_inc(t,s);
    mpf_div(s,t,s);
    if pos then s_mpf_chs(s);
    mpf_copyp(s,b);
    mpf_clear2(s,t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_tanpi(const a: mp_float; var b: mp_float);
  {-Return b := tan(Pi*a)}
begin
  mpf_trig_ex(a, true, nil, nil, @b);
end;


{---------------------------------------------------------------------------}
procedure mpf_todecimal_n(const a: mp_float; ndd: word; str: pchar8; maxlen: word);
  {-Convert an mp_float to decimal scientific representation, ndd>0 is the}
  { total number of digits (including one digit before the '.')}
begin
  mpf_toradix_n(a, 10, ndd, str, maxlen);
end;


{---------------------------------------------------------------------------}
procedure mpf_todecimal_alt(const a: mp_float; ndd: word; str: pchar8; maxlen: word);
  {-Convert an mp_float to decimal alternative representation, ndd>0 is the total}
  { number of digits, trailing '.' or '0' are suppressed. If a is too large or}
  { too small, scientific representation is used.}
begin
  mpf_toradix_alt(a, 10, ndd, str, maxlen);
end;


{---------------------------------------------------------------------------}
function mpf_todouble(const a: mp_float): double;
  {-Convert a to double, +-inf if too large}
const
  MaxDblHex: TMPHexDblW = ($ffff,$ffff,$ffff,$7fef);
var
  x: mp_float;
  p: longint;
begin
  mpf_todouble := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_todouble');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {make local copy of a}
  p := a.bitprec;
  if p<53 then p := 53;
  mpf_initp(x, p);
  if mp_error=MP_OKAY then begin
    mpf_set_dbl(x,double(MaxDblHex));
    if  mpf_cmp_mag(a,x) > 0 then begin
      if a.mantissa.sign=MP_NEG then mpf_todouble := DblNegInf
      else mpf_todouble := DblPosInf;
    end
    else begin
      mpf_copy(a,x);
      s_mpf_normalizep(x, 53);
      mpf_todouble := mp_todouble_ex(x.mantissa, x.exponent);
    end;
    mpf_clear(x);
  end;
end;


{$ifndef EXT64}
{---------------------------------------------------------------------------}
procedure mpf_set_ext(var a: mp_float; x: extended);
  {-Set a to an extended. Error if x = NAN or INF}
type
  txr = packed record
          mant: array[0..7] of byte;
          sexp: word;
        end;
var
  i: integer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_set_ext');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if txr(x).sexp and $7fff = $7fff then begin
    {nan or inf}
    {$ifdef MPC_HaltOnError}
     {$ifdef MPC_UseExceptions}
       raise MPXBadArg.Create('mpf_set_ext: NAN or INF');
     {$else}
       RunError(MP_RTE_BADARG);
     {$endif}
   {$else}
     set_mp_error(MP_BADARG);
     exit;
   {$endif}
  end;
  mp_set(a.mantissa, txr(x).mant[7]);
  for i:=6 downto 0 do begin
    mp_shl(a.mantissa, 8, a.mantissa);
    mp_add_d(a.mantissa, txr(x).mant[i], a.mantissa);
  end;
  if txr(x).sexp and $8000 <> 0 then s_mpf_chs(a);
  i := txr(x).sexp and $7fff;
  if i<>0 then a.exponent := longint(i) - ($3fff+63)
  else a.exponent := -16445;
  s_mpf_normalize(a);
end;


{---------------------------------------------------------------------------}
function mpf_toextended(const a: mp_float): extended;
  {-Convert a to double, +-inf if too large}
const
  MaxExtHex: TMPHexExtW = ($ffff,$ffff,$ffff,$ffff,$7ffe); {MaxExtended as Hex}
var
  x: mp_float;
  p: longint;
begin
  mpf_toextended := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_toextended');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {make local copy of a}
  p := a.bitprec;
  if p<64 then p := 64;
  mpf_initp(x, p);
  if mp_error=MP_OKAY then begin
    mpf_set_ext(x,extended(MaxExtHex));
    if  mpf_cmp_mag(a,x) > 0 then begin
      if a.mantissa.sign=MP_NEG then mpf_toextended := DblNegInf
      else mpf_toextended := DblPosInf;
    end
    else begin
      mpf_copy(a,x);
      s_mpf_normalizep(x, 64);
      mpf_toextended := mp_toextended_ex(x.mantissa, x.exponent);
    end;
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_add_ext(const a: mp_float; b: extended; var c: mp_float);
  {-Calculate c = a+b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_add_ext');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a) then begin
    mpf_set_ext(c,b);
    exit;
  end;
  if b=0.0 then begin
    mpf_copyp(a,c);
    exit;
  end;
  {make local copy of b}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_ext(x,b);
    s_mpf_inc(x,a);
    mpf_exch(x,c); {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_sub_ext(const a: mp_float; b: extended; var c: mp_float);
  {-Calculate c = a-b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_sub_ext');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if b=0.0 then begin
    mpf_copyp(a,c);
    exit;
  end;
  if s_mpf_is0(a) then begin
    {c = -b}
    mpf_set_ext(c,-b);
    exit;
  end;
  {make local copy of a}
  mpf_initp(x, c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_ext(x,-b);
    {c = (-b)+a}
    s_mpf_inc(x,a);
    mpf_exch(x,c);    {Note: x.bitprec=c.bitprec}
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_mul_ext(const a: mp_float; b: extended; var c: mp_float);
  {-Calculate c = a*b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_mul_ext');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(a) or (b=0.0) then begin
    mpf_set0(c);
    exit;
  end;
  mpf_initp(x,c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_ext(x,b);
    mpf_mul(a,x,c);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_div_ext(const a: mp_float; b: extended; var c: mp_float);
  {-Calculate c = a/b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_div_ext');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if b=0.0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_div_ext: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mpf_initp(x,c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_ext(x,b);
    mpf_div(a,x,c);
    mpf_clear(x);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_divr_ext(a: extended; const b: mp_float; var c: mp_float);
  {-Calculate c = a/b}
var
  x: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(b) or mpf_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_div_ext');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mpf_is0(b) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_divr_ext: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mpf_initp(x,c.bitprec);
  if mp_error=MP_OKAY then begin
    mpf_set_ext(x,a);
    mpf_div(x,b,c);
    mpf_clear(x);
  end;
end;



{$endif}


{---------------------------------------------------------------------------}
procedure mpf_tohex_n(const a: mp_float; ndd: word; str: pchar8; maxlen: word);
  {-Convert an mp_float to hexadecimal scientific representation, ndd>0 is the}
  { total number of digits (including one digit before the '.')}
begin
  mpf_toradix_n(a, 16, ndd, str, maxlen);
end;


{---------------------------------------------------------------------------}
procedure mpf_toradix_n(const a: mp_float; radix, ndd: word; str: pchar8; maxlen: longint);
  {-Convert an mp_float to radix scientific representation, ndd>0 is the}
  { total number of digits (including one digit before the '.')}
begin
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_toradix_n');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mpf_toradix_n(a, radix, ndd, maxlen, str);
end;


{---------------------------------------------------------------------------}
procedure mpf_toradix_alt(const a: mp_float; radix, ndd: word; str: pchar8; maxlen: longint);
  {-Convert an mp_float to radix alternative representation, ndd>0 is the total}
  { number of digits, trailing '.' or '0' are suppressed. If a is too large or}
  { too small, scientific representation is used. NOTE: no radix prefix/suffix!}
begin
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_toradix_alt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mpf_toradix_alt(a, radix, ndd, maxlen, str);
end;


{---------------------------------------------------------------------------}
procedure mpf_trig_ex(const a: mp_float; mulpi: boolean; pc,ps,pt: pmp_float);
  {-Calculate pc^=cos(a), ps^=sin(a), pt^=tan(a) if pointers are <> nil.}
  { If mulpi, calculate pc^=cos(a*Pi), ps^=sin(a*Pi), pt^=tan(a*Pi).}
var
  x,y,z: mp_float;
  prec: longint;
  lx,i: longint;
  negc, negst, sita, oddm: boolean;
label
  done;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a)
       or ((pc<>nil) and mpf_not_init(pc^))
       or ((ps<>nil) and mpf_not_init(ps^))
       or ((pt<>nil) and mpf_not_init(pt^))
    then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_trig_ex');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if (pc=nil) and (ps=nil) and (pt=nil) then exit;
  {handle trivial case}
  if s_mpf_is0(a) then begin
    if pc<>nil then mpf_set1(pc^);
    if ps<>nil then mpf_set0(ps^);
    if pt<>nil then mpf_set0(pt^);
    exit;
  end;

  {get maximum precision of output mp_floats}
  prec := 0;
  if (pc<>nil) and (pc^.bitprec>prec) then prec := pc^.bitprec;
  if (ps<>nil) and (ps^.bitprec>prec) then prec := ps^.bitprec;
  if (pt<>nil) and (pt^.bitprec>prec) then prec := pt^.bitprec;

  mpf_initp3(x,y,z,prec + (prec div 8));
  if mp_error<>MP_OKAY then exit;

  mpf_abs(a,x);
  negc  := false;
  negst := s_mpf_is_neg(a);
  sita  := (ps<>nil) or (pt<>nil);

  {Range reduction mod Pi}
  if mulpi then begin
    {use x=a*Pi, first reduce a mod 1 then multiply by Pi}
    s_mpf_frac(x,x,@oddm);
    if not s_mpf_is0(x) then begin
      mpf_set_pi(y);
      mpf_mul(x,y,x);
    end;
  end
  else begin
    {Reduce modulo pi. WARNING: there may be severe loss of precision}
    {in the mod pi reduction, if the "real" mantissa m of a=m*2^e has}
    {more than prec significant high bits!}
    s_mpf_mod_pi2k(x,0,x,oddm);
    if mp_error<>MP_OKAY then goto done;
  end;

  if oddm then begin
    {change sign because odd multiple of Pi has been used in range reduction}
    negc  := not negc;
    negst := not negst;
  end;

  lx := s_mpf_ldx(x);
  if s_mpf_is0(x) or (lx <= -(prec div 2)) then begin
    {x is zero or small enough to use only one term of Taylor series}
    if pc<>nil then begin
      mpf_set1(pc^);
      if negc then s_mpf_chs(pc^);
    end;
    if negst then s_mpf_chs(x);
    if ps<>nil then mpf_copyp(x,ps^);
    if pt<>nil then mpf_copyp(x,pt^);
    goto done;
  end;

  inc(lx,2);
  if lx>0 then begin
    {range reduction 2: make x < 1/4; "worst" x=0.499999....}
    s_mpf_incexp(x,-lx);
  end;

  {Taylor series h(x) = cos(x) - 1 = - x^2/2! + x^4/4! - x^6/6! + x^6/8! -}
  {Calculate h(x) instead of cos(x) to reduce rounding errors for x near 0}
  {z will accumulate the sum, y is the i-th term.}
  mpf_sqr(x,x);
  mpf_set0(z);
  mpf_set1(y);
  i := 1;
  repeat
    mpf_mul(y,x,y);
    if i<46340 then mpf_div_int(y,-i*succ(i),y)
    else mpf_div_dbl(y,-i*(1.0+i),y);
    inc(i,2);
  until not s_mpf_incf(z,y);

  for i:=1 to lx do begin
    {undo range reduction 2 using cos(x)-1 = 2(cos(x/2)^2 - 1)}
    {or h(x) = 2*h(x/2)*(h(x/2)+2), ie iterate z=2*z*(z+2)}
    mpf_set_int(y,2);
    if s_mpf_incf(y,z) then begin
      mpf_mul(z,y,z);
      inc(z.exponent);
    end
    else begin
      {z near 2, ie z=2*z*(z+2) near 4z}
      inc(z.exponent,2);
    end;
  end;

  {here z = h = cos - 1}
  if sita then begin
    {sin and/or tan are requested, first calculate sin = sqrt(1-cos^2) =}
    {sqrt((1-cos)*(1+cos)) = sqrt(|(h+2)*h|) to reduce rounding errors!}
    mpf_set_int(y,2);
    s_mpf_inc(y,z);
    mpf_mul(z,y,y);
    s_mpf_abs(y);
    mpf_sqrt(y,y);
    {then adjust sign if necessary}
    if negst then s_mpf_chs(y);
  end;

  {cos = h+1}
  s_mpf_inc1(z);
  {adjust sign if necessary}
  if negc then s_mpf_chs(z);

  {here y = sin, z = cos}
  if pc<>nil then mpf_copyp(z,pc^);
  if pt<>nil then mpf_div(y,z,pt^);
  if ps<>nil then mpf_copyp(y,ps^);

done:
  mpf_clear3(x,y,z);
end;


{---------------------------------------------------------------------------}
procedure mpf_trig(const a: mp_float; pc,ps,pt: pmp_float);
  {-Calculate pc^=cos(a), ps^=sin(a), pt^=tan(a) if pointers are <> nil}
begin
  mpf_trig_ex(a,false,pc,ps,pt);
end;


{---------------------------------------------------------------------------}
procedure mpf_trunc(const a: mp_float; var b: mp_int);
  {-Truncate mp_float to mp_int}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_trunc');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  with a do begin
    if exponent<0 then mp_shr(mantissa,-exponent,b)
    else if exponent>0 then mp_shl(mantissa,exponent,b)
    else mp_copy(mantissa,b);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_vers(const a: mp_float; var b: mp_float);
  {-Calculate the versine b = vers(a) = 1-cos(a)}
begin
  {Use vers(a) = 2*hav(a)}
  mpf_hav(a,b);
  if mp_error=MP_OKAY then mpf_mul_2k(b,1,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_writeln(const msg: mp_string; const a: mp_float; ndd: word);
  {-Writeln a to output with leading msg}
begin
  write(msg);
  mpf_write_decimal(output,a, ndd);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure mpf_write_decimal(var tf: system.text; const a: mp_float; ndd: word);
  {-Write an mp_float to file tf using decimal scientific representation, ndd>0}
  { is the total number of digits (including one digit before the '.')}
begin
  mpf_write_radix(tf, a, 10, ndd);
end;


{---------------------------------------------------------------------------}
procedure mpf_write_decimal_alt(var tf: system.text; const a: mp_float; ndd: word);
  {-Write an mp_float to file tf using decimal alternative representation, ndd>0 is}
  { the total number of digits, trailing '.' or '0' are suppressed. If a is too}
  { large or too small, scientific representation is used.}
begin
  mpf_write_radix_alt(tf, a, 10, ndd);
end;


{---------------------------------------------------------------------------}
procedure mpf_write_radix(var tf: system.text; const a: mp_float; radix,ndd: word);
  {-Write an mp_float to file tf using radix scientific representation, ndd>0}
  { is the total number of digits (including one digit before the '.')}
var
  ls,lw,la: longint;
  pc,pt: pchar8;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_write_radix');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  ls := ndd+20;
  if (ls<=$FF00) and (mp_error=MP_OKAY) then begin
    lw := ls and $ffff;
    pc := mp_getmem(lw);
    if pc<>nil then begin
      pt := pc;
      {save alloc count, lw is changed in s_mp_toradix_n} {1.0.14}
      la := lw;
      s_mpf_toradix_n(a,radix, ndd, lw, pt);
      write(tf, pc);
      mp_freemem(pointer(pc),la);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_write_radix_alt(var tf: system.text; const a: mp_float; radix,ndd: word);
  {-Write an mp_float to file tf using radix alternative representation, ndd>0 is}
  { the total number of digits, trailing '.' or '0' are suppressed. If a is too}
  { large or too small, scientific representation is used. NOTE: no radix prefix/suffix!}
var
  ls,lw,la: longint;
  pc,pt: pchar8;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_write_radix_alt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  ls := ndd+20;
  if (ls<=$FF00) and (mp_error=MP_OKAY) then begin
    lw := ls and $ffff;
    pc := mp_getmem(lw);
    if pc<>nil then begin
      pt := pc;
      {save alloc count, lw is changed in s_mp_toradix_n} {1.0.14}
      la := lw;
      s_mpf_toradix_alt(a,radix, ndd, lw, pt);
      write(tf, pc);
      mp_freemem(pointer(pc),la);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpf_abs(var a: mp_float);
  {-Absolute value of an mp_float, no init check}
begin
  with a.mantissa do begin
    if (mp_error=MP_OKAY) and (magic=MP_MAGIC) then sign := MP_ZPOS;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpf_add1(const a: mp_float; var b: mp_float);
  {-Calculate b = a+1; no init checks}
var
  bd,e1,e2,ed: longint;
  m: mp_int;
begin

  if s_mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;

  bd := s_mpf_ldx(a) - 1;

  {check if total bitsizes are too different}
  if abs(bd) > b.bitprec+1 then begin
    if bd<0 then begin
      {a is almost 0 compared to 1}
      mpf_set1(b);
    end
    else begin
      mpf_copyp(a,b);
    end;
    exit;
  end;
  if bd<0 then begin
    e1 := 0;
    e2 := a.exponent;
  end
  else begin
    e1 := a.exponent;
    e2 := 0;
  end;
  ed := e1-e2;

  mp_init(m);
  if mp_error<>MP_OKAY then exit;

  if ed>0 then begin
    if bd<0 then begin
      mp_2expt(m,ed);
      mp_add(a.mantissa, m, b.mantissa);
      b.exponent := e2;
    end
    else begin
      mp_shl(a.mantissa, ed, b.mantissa);
      s_mp_add_d(b.mantissa,1,b.mantissa);
      b.exponent := 0;
    end;
  end
  else if ed<0 then begin
    if bd<0 then begin
      mp_shl(a.mantissa, -ed, b.mantissa);
      s_mp_add_d(b.mantissa,1,b.mantissa);
      b.exponent := 0;
    end
    else begin
      mp_2expt(m,-ed);
      mp_add(a.mantissa, m, b.mantissa);
      b.exponent := e1;
    end;
  end
  else begin
    {ed=0}
    s_mp_add_d(a.mantissa,1,b.mantissa);
  end;
  mp_clear(m);
  s_mpf_normalize(b);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_agm(const a,b: mp_float; var c: mp_float; ps: pmp_float);
  {-Calculate c = AGM(|a|,|b|), if ps<>nil set ps^=sum(2^k*c_k^2); ps<>@c, no init check}
var
  g,m,t: mp_float;
  cnt, bd, lim, maxi: longint;
  cchg: boolean;
begin

  if mpf_is0(a) or mpf_is0(b) then begin
    mpf_set0(c);
    if ps<>nil then mpf_set0(ps^);
    exit;
  end;

  mpf_initp3(g,m,t,c.bitprec);
  {m ~ arith, g ~ geo,  g <= m}
  if mpf_cmp_mag(a,b)=MP_LT then begin
    mpf_abs(a,g);
    mpf_abs(b,m);
  end
  else begin
    mpf_abs(a,m);
    mpf_abs(b,g);
  end;
  if ps<>nil then mpf_set0(ps^);

  {limiting precision}
  lim := c.bitprec-1;

  {Upper bound for iteration cnt}
  bd  := abs(s_mpf_ldx(m)) + abs(s_mpf_ldx(g));
  maxi:= 4+round((ln(1.0+bd)+ln(lim))/ln(2));

  cnt := 0;
  repeat
    inc(cnt);
    mpf_sub(g,m,t);
    if s_mpf_is0(t) then break;
    bd := s_mpf_ldx(m) - s_mpf_ldx(t);
    mpf_mul(g,m,g);
    dec(t.exponent);
    {remember if m has changed, break later. If not, update sum with 2^cnt*t^2}
    cchg := s_mpf_incf(m,t);
    if ps<>nil then begin
      mpf_sqr(t,t);
      mpf_mul_2k(t,cnt,t);
      s_mpf_inc(ps^,t);
    end;
    {avoid sqrt calculation if done}
    if (not cchg) or (bd>=lim) or (cnt>maxi) then break;
    mpf_sqrt(g,g);
  until false;
  mpf_exch(m,c);
  mpf_clear3(g,m,t);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_ccell12(const k: mp_float; pCK, pCE: pmp_float);
  {-Complementary complete elliptic integrals of the 1st and 2nd kind using}
  { AGM algorithm; k<>0 and pCK <> pCE, with init checks, real parts of K'(k) if k<-1}
var
  c0,s,te,tk: mp_float;
  prec: longint;
begin

  {Ref: HMF[41] Abramowitz/Stegun; Handbook of Mathematical Functions, 1965}
  {     Ch. 17: Elliptic Integrals, Sec. 17.6: Arithmetic-Geometric Mean}

  { F(x,k) = integral(1/sqrt(1-t^2)/sqrt(1-k^2*t^2),t=0..x)}
  { CK = CK(k) = K'(k) = K(k') = F(1, sqrt(1-k^2))}

  { E(x,k) = integral(sqrt(1-k^2*t^2)/sqrt(1-t^2),t=0..x)}
  { CE = CE(k) = E'(k) = E(k') = E(1, sqrt(1-k^2))}

  { With the HMF notation k = m1^2, ie  m=1-k^2, cf the Ex3/4 in 17.8:}
  { CK=CK(1/9)=K(80/81) ~ 3.59145000, CE=CE(1/9)=E(80/81) ~ 1.019106047}

  { Maple V, Rel.4  notation is  EllipticCK(k) and  EllipticCE(k)}

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(k) or ((pCK<>nil) and mpf_not_init(pCK^))
       or ((pCE<>nil) and mpf_not_init(pCE^)) then
    begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mpf_ccell12');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Check k<>0 and  @CK <> @CE}
  if pCK=pCE then begin
    if pCK=nil then exit;
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('s_mpf_ccell12: @CK = @CE');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;
  if (s_mpf_is0(k)) and (pCK<>nil) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('s_mpf_ccell12: k = 0, pCK <> nil');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {get maximum precision of output mp_floats}
  if pCK<>nil then prec := pCK^.bitprec else prec := 0;
  if (pCE<>nil) and (prec < pCE^.bitprec) then prec := pCE^.bitprec;

  {do calculation with slightly increased working precision}
  mpf_initp4(c0,s,te,tk,prec + PINC);

  {remember c0^2=1-k^2 and call AGM(1,k)}
  mpf_set1(s);
  mpf_sqr(k,c0);
  mpf_sub(s,c0,c0);
  s_mpf_agm(s{=1},k,tk,@te);

  {CE = 1 - 0.5(c0^2 + CE)}
  s_mpf_inc(te,c0);
  s_mpf_incexp(te,-1);
  mpf_sub(s,te,te);

  {CK = pi/2/AGM = pi/2/CK}
  mpf_set_pi2k(s,-1);
  mpf_div(s,tk,tk);

  {E'(k) = K'(k)*(1 - 0.5*(c_0^2 + sum(2^k*c_k^2)))}
  if pCE<>nil then begin
    {CE(k) is negative for k<-1}
    mpf_mul(te,tk,te);
    if mpf_cmp_dbl(k,-1)=MP_LT then mpf_chs(te, pCE^)
    else mpf_abs(te,pCE^);
  end;
  if pCK<>nil then mpf_copyp(tk,pCK^);
  mpf_clear4(c0,s,te,tk);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_chs(var a: mp_float);
  {-Change sign of an mp_float, no init check}
begin
  with a.mantissa do begin
    if (mp_error=MP_OKAY) and (magic=MP_MAGIC) then begin
      if used>0 then sign := sign xor (MP_NEG xor MP_ZPOS)
      else sign := MP_ZPOS;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function s_mpf_cmp_mag(const a,b: mp_float): integer;
  {-Compare magnitude of two mp_floats, return sign(|a|-|b|), no init check}
var
  ba,bb: longint;
  t: mp_float;
  sig: word;
begin
  {assign default to keep D6+ happy, combine error test with a=b=0}
  s_mpf_cmp_mag := 0;
  if (mp_error<>MP_OKAY) or (s_mpf_is0(a) and s_mpf_is0(b)) then exit;

  if s_mpf_is0(a) then begin
    {a=0,  |b|>0}
    s_mpf_cmp_mag := -1;
    exit;
  end;
  if s_mpf_is0(b) then begin
    {|a|>0, b=0}
    s_mpf_cmp_mag := 1;
    exit;
  end;

  {if exponents are equal, compare mantissas}
  if a.exponent=b.exponent then begin
    s_mpf_cmp_mag := mp_cmp_mag(a.mantissa, b.mantissa);
    exit;
  end;

  {compare effective bitsizes}
  ba := s_mpf_ldx(a);
  bb := s_mpf_ldx(b);
  if ba>bb then begin
    s_mpf_cmp_mag := 1;
    exit;
  end
  else if ba<bb then begin
    s_mpf_cmp_mag := -1;
    exit;
  end;

  {effective bitsizes are equal but exponents are different}
  {subtract using the maximum bitprec of a and b}
  if a.bitprec < b.bitprec then begin
    mpf_initp(t, b.bitprec);
    if mp_error=MP_OKAY then begin
      mpf_copyp(a,t);
      sig := b.mantissa.sign;
      t.mantissa.sign := sig;
      mpf_sub(t,b,t);
    end
    else exit;
  end
  else begin
    mpf_initp(t, a.bitprec);
    if mp_error=MP_OKAY then begin
      mpf_copyp(b,t);
      sig := a.mantissa.sign;
      t.mantissa.sign := sig;
      mpf_sub(a,t,t);
    end
    else exit;
  end;

  {here t=sign(b)|a|-b or t=a - sign(a)|b|}
  if s_mpf_is0(t) then s_mpf_cmp_mag := 0
  else begin
    {if sig=negative reverse normal compare}
    if sig=MP_NEG then begin
      if t.mantissa.sign=MP_NEG then s_mpf_cmp_mag := 1
      else s_mpf_cmp_mag := -1;
    end
    else begin
      if t.mantissa.sign=MP_NEG then s_mpf_cmp_mag := -1
      else s_mpf_cmp_mag := 1;
    end;
  end;
  mpf_clear(t);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_dec1(var a: mp_float);
  {-Calculate a = a-1; no init check}
begin
  {a-1 = -(-a + 1)}
  s_mpf_chs(a);
  s_mpf_inc1(a);
  s_mpf_chs(a);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_frac(const a: mp_float; var b: mp_float; podd: pBoolean);
  {-Set b to the fractional part of a; frac(x)=x-int(x), if podd<>nil set podd^=odd(trunc(a))}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mpf_frac');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpf_copyp(a,b);
  with b do begin
    if exponent>=0 then begin
      {a is an integer}
      if podd<>nil then podd^ := (exponent=0) and mp_isodd(mantissa);
      mpf_set0(b);
    end
    else begin
      {exponent<0, test bit 0 of trunc(a) = mantissa shr (-exponent)}
      if podd<>nil then podd^ := mp_isbit(mantissa, -exponent);
      s_mp_mod_2k(mantissa,-exponent,mantissa);
      s_mpf_normalize(b);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpf_inc(var x: mp_float; const y: mp_float);
  {-Calculate x = x+y}
begin
  if s_mpf_incf(x,y) then ;
end;


{---------------------------------------------------------------------------}
procedure s_mpf_incexp(var a: mp_float; x: longint);
  {-Increment a.exponent by x, do nothing if a.mantissa=0}
var
  y : longint;
begin
  if (mp_error<>MP_OKAY) or mp_is0(a.mantissa) then exit;
  if add32_ovr(a.exponent,x,y) then begin
    if a.exponent<0 then begin
      {floating underflow, set a=0}
      mpf_set0(a);
    end
    else begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXOverflow.Create('s_mpf_incexp');
        {$else}
          RunError(MP_RTE_OVRFLOW);
        {$endif}
      {$else}
        set_mp_error(MP_OVERFLOW);
        exit;
      {$endif}
    end;
  end
  else a.exponent := y;
end;


{---------------------------------------------------------------------------}
function s_mpf_incf(var x: mp_float; const y: mp_float): boolean;
  {-Calculate x = x+y; return true if x is changed}
var
  m: mp_int;
  bd,e1,e2,ed: longint;
  ovr: boolean;
begin
  s_mpf_incf := true;


  if mp_error<>MP_OKAY then exit;
  if s_mpf_is0(y) then begin
    s_mpf_incf := false;
    exit;
  end;

  if s_mpf_is0(x) then begin
    mpf_copyp(y,x);
    exit;
  end;
  if @x=@y then begin
    {x is non zero, so increment exponent is OK}
    s_mpf_incexp(x,1);
    exit;
  end;

  {check if total bitsizes are very different}
  ovr := add32_ovr(s_mpf_ldx(x), -s_mpf_ldx(y), bd);
  if (not ovr) and (abs(bd) > x.bitprec+1) then begin
    if bd<0 then begin
      {x is almost 0 compared to y}
      mpf_copyp(y,x);
    end
    else s_mpf_incf := false;
    exit;
  end;

  {assign e1,e2 even if overflow to avoid compiler warnings}
  if bd<0 then begin
    e1 := y.exponent;
    e2 := x.exponent;
  end
  else begin
    e1 := x.exponent;
    e2 := y.exponent;
  end;

  ovr := ovr or add32_ovr(e1, -e2, ed);
  if ovr then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXOverflow.Create('s_mpf_incf');
      {$else}
        RunError(MP_RTE_OVRFLOW);
      {$endif}
    {$else}
      set_mp_error(MP_OVERFLOW);
      exit;
    {$endif}
  end;

  mp_init(m);
  if mp_error<>MP_OKAY then exit;

  if ed>0 then begin
    if bd<0 then begin
      mp_shl(y.mantissa, ed, m);
      mp_add(x.mantissa, m, x.mantissa);
    end
    else begin
      mp_shl(x.mantissa, ed, m);
      mp_add(y.mantissa, m, x.mantissa);
    end;
    x.exponent := e2;
  end
  else if ed<0 then begin
    if bd<0 then begin
      mp_shl(x.mantissa, -ed, m);
      mp_add(y.mantissa, m, x.mantissa);
    end
    else begin
      mp_shl(y.mantissa, -ed, m);
      mp_add(x.mantissa, m, x.mantissa);
    end;
    x.exponent := e1;
  end
  else begin
    {ed=0}
    mp_add(y.mantissa, x.mantissa, x.mantissa);
  end;

  mp_clear(m);
  s_mpf_normalize(x);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_inc1(var a: mp_float);
  {-Calculate a = a+1; no init check}
begin
  s_mpf_add1(a,a);
end;


{---------------------------------------------------------------------------}
function s_mpf_is0(const a: mp_float): boolean;
  {-Return true if a=0, no init checks}
begin
  s_mpf_is0 := (a.exponent=0) and (a.mantissa.used=0);
end;


{---------------------------------------------------------------------------}
function s_mpf_is_ge0(const a: mp_float): boolean;
  {-Return true if a>=0, no init checks}
begin
  s_mpf_is_ge0 := (a.mantissa.sign=MP_ZPOS);
end;


{---------------------------------------------------------------------------}
function s_mpf_is_gt0(const a: mp_float): boolean;
  {-Return true if a>0, no init checks}
begin
  s_mpf_is_gt0 := (a.mantissa.sign=MP_ZPOS) and (a.mantissa.used>0)
end;


{---------------------------------------------------------------------------}
function s_mpf_is_le0(const a: mp_float): boolean;
  {-Return true if a<=0, no init checks}
begin
  s_mpf_is_le0 := (a.mantissa.sign=MP_NEG) or ((a.exponent=0) and (a.mantissa.used=0));
end;


{---------------------------------------------------------------------------}
function s_mpf_is_neg(const a: mp_float): boolean;
  {-Return true if a<0, no init checks}
begin
  s_mpf_is_neg := a.mantissa.sign=MP_NEG;
end;


{---------------------------------------------------------------------------}
function s_mpf_ldx(const a: mp_float): longint;
  {-Return ldx(a) = a.exponent+mp_bitsize(a.mantissa), no init checks, a must be}
  { normalized, ldx(a) = floor(log2(|a|))+1, -MaxLongint if a=0, MaxLongint if overflow}
var
  ldx: longint;
begin
  {s_mpf_ldx := a.exponent + mp_bitsize(a.mantissa);}
  with a.mantissa do begin
    if used>0 then begin
      if add32_ovr(a.exponent, longint(used-1)*DIGIT_BIT + bitsize32(pdigits^[used-1]), ldx) then begin
        s_mpf_ldx := MaxLongint;
      end
      else s_mpf_ldx := ldx;
    end
    else s_mpf_ldx := -MaxLongint;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpf_numbpart(n: longint; var p: mp_float);
  {-Calculate the number of partitions of n with Hardy-Ramanujan-Rademacher formula}
var
  pfull,pq,q,sn: longint;
  bn,cn,x,y: mp_float;

  {$ifdef old_partv}
    {old version of code, can still be used for test with -dold_partv}
    {---------------------------------------------------------------------------}
    procedure sigma(h,k: longint; var a,b: longint);
      {-Calculate Dedekind sum sigma(h,k) = a+b/k}
    var
      aa,p,pp,r,s: longint;
    begin
      {sigma is computed with Knuth's algorithm [3], Chap 3.3.3, Ex.17}
      {specialized to c=0}
      a := 0;
      b := h;
      p := 1;
      pp:= 0;
      s := 1;
      while h>0 do begin
        aa := k div h;
        a  := a + (aa-3)*s;
        if h=1 then b := b +p*s;
        s := -s;
        r := k-aa*h;
        k := h;
        h := r;
        r := aa*p + pp;
        pp:= p;
        p := r;
      end;
    end;

    {-------------------------------------------------------------------}
    procedure Lq(n,q: longint; var L: mp_float);
    var
      h,a,b: longint;
      x,y: mp_float;
    begin
      if q=1 then mpf_set1(L)
      else begin
        mpf_initp2(x,y,L.bitprec);
        if mp_error<>MP_OKAY then exit;
        mpf_set0(L);
        {L(n,q) = sum(h=1..q-1 and gcd(h,q)=1,cos((g(h,q)-2*h*n)*Pi/q)))}
        for h:=1 to q-1 do begin
          if mp_error<>MP_OKAY then exit;
          if gcd32(h,q)=1 then begin
            {compute g(h,q)=q*s(h,q)=q/12*sigma(h,q)=q/12*(a+b/q)=(a*q+b)/12}
            if q<3 then mpf_set0(x)
            else begin
              sigma(h,q,a,b);
              mpf_set_int(x,a);
              mpf_mul_int(x,q,x);
              mpf_add_int(x,b,x);
              mpf_div_d(x,12,x);
            end;
            mpf_set_int(y,2*h);
            mpf_mul_int(y,n,y);
            mpf_sub(x,y,x);
            mpf_div_int(x,q,x);
            {x=cos(Pi*x)}
            mpf_trig_ex(x,true,@y,nil,nil);
            mpf_add(L,y,L);
          end;
        end;
        mpf_clear2(x,y);
      end;
    end;
  {$else}
    {-------------------------------------------------------------------}
    procedure Lq(n,k: longint; var s: mp_float);
      {-Compute A_k(n) = Lq(n,k) with Johansson method}
    var
      r,m,l: longint;
      x,y: mp_float;
    begin
      {This faster method for Lq(n,k) is based on the simple}
      {brute force Algorithm 1 from Fredrik Johansson [38]. }
      if k=1 then mpf_set1(s)
      else if k=2 then begin
        mpf_set1(s);
        if odd(n) then s_mpf_chs(s);
      end
      else begin
        mpf_initp2(x,y,s.bitprec);
        if mp_error<>MP_OKAY then exit;
        mpf_set0(s);
        m := n mod k;
        r := 2;
        for l:=0 to 2*k-1 do begin
          if mp_error<>MP_OKAY then exit;
          if m=0 then begin
            {x = (6*l+1)/6k}
            mpf_set_int(x,6*l+1);
            mpf_div_int(x,6*k,x);
            {y = cos(Pi*x)}
            mpf_trig_ex(x,true,@y,nil,nil);
            if odd(l) then mpf_sub(s,y,s)
            else mpf_add(s,y,s);
          end;
          m := m+r;
          if m>=k then m := m-k;
          r := r+3;
          if r>=k then r := r-k;
        end;
        mpf_set_int(x,k);
        mpf_div_d(x,3,x);
        mpf_sqrt(x,y);
        mpf_mul(s,y,s);
        mpf_clear2(x,y);
      end;
    end;
  {$endif}

  {-------------------------------------------------------------------}
  procedure Psiq(q: longint; const b, c: mp_float; var p: mp_float);
  var
    a,d: mp_float;
  begin
    {Psi(n,q) = (sqrt(q)/(2*sqrt(2)*b*Pi))*(a*cosh(d)-(sinh(d)/c))}
    {with a = sqrt(2/3)*Pi/q,  b = n-1/24,  c = sqrt(b),  d = a*c}
    mpf_initp2(a,d,p.bitprec);
    if mp_error<>MP_OKAY then exit;
    {a = sqrt(2/3)*Pi/q}
    mpf_set_pi(a);
    mpf_set_int(d,2);
    mpf_div_d(d,3,d);
    mpf_sqrt(d,d);
    mpf_mul(a,d,a);
    mpf_div_int(a,q,a);
    {d = a*c}
    mpf_mul(a,c,d);
    {p = sinh(d)/c,  d = a*cosh(d)}
    mpf_sinhcosh(d,p,d);
    mpf_mul(d,a,d);
    mpf_div(p,c,p);
    {p = a*cosh(d) - sinh(d)/c}
    mpf_sub(d,p,p);
    {d = sqrt(q/8)/(b*Pi)}
    mpf_set_pi(d);
    mpf_mul(d,b,d);
    mpf_set_int(a,q);
    s_mpf_incexp(a, -3);
    mpf_sqrt(a,a);
    mpf_div(a,d,d);
    {Final result = d*p}
    mpf_mul(d,p,p);
    mpf_clear2(a,d);
  end;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(p) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mpf_numbpart');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if n<0 then begin
    mpf_set0(p);
    exit;
  end
  else if n<2 then begin
    mpf_set1(p);
    exit;
  end;

  {c.f. the Pari implementation [12] based on Ralf Stephan's script and code}
  {numpart(n) = sum(q=1,5 + 0.24*sqrt(n),L(n,q)*Psi(n,q)))}

  {Get estimated bit prec based on first term of HRR series}
  {p(n) ~ exp(pi*sqrt(2n/3))/(4n*sqrt(3))}
  pfull := MPF_MIN_PREC + 128 + round(sqrt(2.0*n/3.0)*pi/ln(2.0));
  if pfull > MPF_MAX_PREC then begin
    {$ifdef MPC_ArgCheck}
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('s_mpf_numbpart');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    {$endif}
    mpf_set_int(p,-1);
    exit;
  end;

  p.bitprec := pfull;

  sn := round(5+0.24*sqrt(n));
  mpf_set0(p);
  mpf_initp4(bn,cn,x,y,p.bitprec);

  {Pre-compute bn,cn, because they are independent of q}
  {bn = n-1/24}
  mpf_set1(bn);
  mpf_div_int(bn,-24,bn);
  mpf_add_int(bn,n,bn);
  {cn = sqrt(bn)}
  mpf_sqrt(bn,cn);

  {sum backwards with increasing precision}
  for q:=sn downto 1 do begin
    pq := 8 + pfull div q;
    if pq<32 then pq := 32;
    x.bitprec := pq;
    y.bitprec := pq;
    Lq(n,q,x);
    if s_mpf_is0(x) or (s_mpf_ldx(x)<-pq) then continue;
    Psiq(q,bn,cn,y);
    mpf_mul(x,y,x);
    mpf_add(p,x,p);
  end;
  mpf_clear4(bn,cn,x,y);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_mod_pi2k(const a: mp_float; k: integer; var b: mp_float; var oddm: boolean);
  {-Calculate b = a mod pi*2^k, c in [0,pi*2^k); oddm: odd multiple of pi*2^k}
  { was used for reduction; no init check; extended precision used if necessary.}
var
  x,pi2k: mp_float;
  px: longint;
begin
  oddm := false;
  if mp_error<>MP_OKAY then exit;

  if s_mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;

  {Use extended precision greater than b.bitprec if necessary}
  {k should be small, ie about -2 .. 2, otherwise adjust px}
  px := s_mpf_ldx(a);
  if px < 0 then px := 0;
  px := px + b.bitprec + PINC;

  if px >= MPF_MAX_PREC then begin
    set_mp_error(MP_PRECISION);
    exit;
  end;

  mpf_initp2(x,pi2k,px);
  if mp_error<>MP_OKAY then exit;

  mpf_set_pi2k(pi2k,k);
  if mpf_cmp_mag(a,pi2k)=MP_LT then begin
    {abs(a) < pi2k, just copy a to x for sign adjustment}
    mpf_copyp(a,x);
  end
  else begin
    {x = a - int(a/pi2k)*pi2k}
    mpf_div(a,pi2k,x);
    mpf_int(x,x);
    oddm := mp_isbit(x.mantissa,-x.exponent);
    mpf_mul(x,pi2k,x);
    mpf_sub(a,x,x);
  end;

  if s_mpf_is_neg(x) then begin
    mpf_add(x,pi2k,x);
    oddm := not oddm;
  end;
  {$ifdef MPC_USE_Assert}
    {only if assert supported by compiler or debug}
    assert((x.mantissa.sign=MP_ZPOS) and mpf_is_lt(x,pi2k),' 0 <= a mod pi*2^k < pi*2^k');
  {$endif}
  s_mpf_normalizep(x,b.bitprec);
  mpf_exch(x,b);
  mpf_clear2(x,pi2k);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_normalize(var a: mp_float);
  {-Normalize an mp_float}
var
  bs,diff: longint;
  dm,d0: mp_digit;
  up: TNInt;
  mustround,sticky: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  with a do begin
    bs := mp_bitsize(mantissa);
    if bs > bitprec then begin
      diff := bs - bitprec;
      s_mpf_incexp(a, diff);
      {mustround true if rounding bit diff-1 is set}
      mustround := mp_isbit(mantissa, pred(diff));
      {sticky bit present if any bit below rounding bit is set}
      sticky    := mustround and (mp_cnt_lsb(mantissa) < pred(diff));
      mp_shr(mantissa, diff, mantissa);
      {mantissa is still non zero}
      if mustround then with mantissa do begin
        d0 := pdigits^[0];
        if odd(d0) then begin
          {must inc/dec an odd number, remember used and high digit}
          up := used;
          dm := pdigits^[up-1];
          {round to even with odd mantissa: inc if positive, dec}
          {if negative: 1.5 -> 2 not 1,  -1.5 -> -2 not -1}
          if sign=MP_ZPOS then mp_inc(mantissa) else mp_dec(mantissa);
          if (up<>used) or (dm<>pdigits^[up-1]) then begin
            {a.used or highest digit changed: check if bitsize>bitprec}
            if mp_bitsize(mantissa)>bitprec then begin
              {postnormalization}
              mp_shr(mantissa, 1, mantissa);
              s_mpf_incexp(a,1);
            end;
          end;
        end
        else begin
          {even mantissa, no postnormalization necessary}
          if sticky then pdigits^[0] := d0+1;
        end;
      end;
    end
    else if bs < bitprec then begin
      if bs=0 then mpf_set0(a)
      else begin
        diff := bitprec - bs;
        s_mpf_incexp(a, -diff);
        mp_shl(mantissa, diff, mantissa);
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpf_normalizep(var a: mp_float; prec: longint);
  {-Normalize an mp_float to bit precision prec}
begin
  if mp_error<>MP_OKAY then exit;
  if (prec>=MPF_MIN_PREC) and (prec<=MPF_MAX_PREC) then begin
    a.bitprec := prec;
    s_mpf_normalize(a);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpf_rint(var a: mp_float);
  {-Set a = round(a)}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mpf_rint');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {if exponent>=0 then a is already an 'integer'}
  if a.exponent<0 then begin
    mpf_round(a, a.mantissa);
    a.exponent := 0;
    s_mpf_normalize(a);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpf_toradix_n(const a: mp_float; radix, ndd: word; var maxlen: longint; var pstr: pchar8);
  {-Convert an mp_float to radix scientific representation, ndd>0 is the}
  { total number of digits (including one digit before the '.')}
const
  ovmsg: array[1..10] of char8 = '*overflow*';
var
  e: longint;
  x,y: mp_float;
  ps1,ps2: pchar8;
  i: word;
  d: mp_digit;
  t: string[20];
  xc: char8;

  procedure appchar(c: char8);
    {-Append char at pstr^}
  begin
    if maxlen>1 then begin
      pstr^ := c;
      inc(pstr);
      dec(maxlen);
    end;
  end;

begin

  if mp_error<>MP_OKAY then exit;
  if (radix < 2) or (radix > MAXRadix) or (maxlen <= ndd) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_toradix_n: radix out of range or maxlen too small');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  pstr^ := #0;
  if ndd=0 then exit;

  if s_mpf_is0(a) then begin
    {special case a=0}
    if maxlen>ndd+2+ord(mp_show_plus) then begin
      if mp_show_plus then appchar('+');
      appchar('0');
      appchar(mp_fract_sep{'.'});
      for i:=2 to ndd do appchar('0');
      pstr^ := #0;
      dec(maxlen);
    end;
    exit;
  end;

  {Step 0: get raw binary exponent e}
  e := s_mpf_ldx(a);
  if add32_ovr(a.exponent,mp_bitsize(a.mantissa),e) then begin
    {This is a conversion problem, bsx(a) overflows}
    for i:=1 to sizeof(ovmsg) do appchar(ovmsg[i]);
    pstr^ := #0;
    dec(maxlen);
    exit;
  end;

  mpf_initp2(x,y,a.bitprec + PINC);

  {Step 1: calculate abs(a) = x*radix^e with 1 <= x < radix}
  mpf_set_int(y,radix); {y=radix in Step 1}

  {Step 1a: get raw radix exponent e and raw x}
  if radix>2 then e := trunc(lograd[radix]*e);
  if e<>0 then begin
    mpf_expt_int(y,e,x);
    mpf_div(a,x,x);
    s_mpf_abs(x);
  end
  else mpf_abs(a,x);

  {Step 1b: make x<radix}
  while (mp_error=MP_OKAY) and mpf_is_ge(x,y) do begin
    mpf_div(x,y,x);
    inc(e);
  end;

  {Step 1c: make x>=1}
  while (mp_error=MP_OKAY) and (s_mpf_ldx(x) < 1) do begin
    mpf_mul_d(x,radix,x);
    dec(e);
  end;

  {Step 2: multiply by radix^ndd, this gives ndd+1 radix digits}
  if mpf_is1(x) then begin
    mp_set_pow(x.mantissa,radix,ndd);
  end
  else begin
    mp_set_pow(y.mantissa,radix,ndd);
    mpf_mul_mpi(x,y.mantissa,x);
    mpf_trunc(x,x.mantissa);
  end;

  {Step 3: divide by radix and round}
  mp_div_d(x.mantissa,radix,@x.mantissa,d);
  d := 2*d;
  if (d > radix) or ((d=radix) and mp_isodd(x.mantissa)) then begin
    mp_inc(x.mantissa);
    if mp_is_eq(x.mantissa,y.mantissa) then begin
      inc(e);
      {V1.4.05 - Fix if rounded x = radix^ndd: increment exponent AND }
      {divide by radix, otherwise there is an extra trailing '0'!}
      mp_div_d(x.mantissa, radix, @x.mantissa, d);
    end;
  end;

  {Step 4; convert exponent part to string}
  if e<>0 then begin
    case radix of
        2: xc := 'B';
       10: xc := 'E';
       16: xc := 'H';
      else xc := '@';
    end;
    str(e,t);
    {D12 warning fix: use char8 type cast}
    if e>0 then t := xc+char8('+')+t else t := xc+t;
  end
  else t := '';

  {Step 5: convert mantissa, fix and combine parts}
  if (mp_error=MP_OKAY) and (mp_radix_size(x.mantissa,radix)+2+length(t)<maxlen) then begin
    {insert sign char}
    if s_mpf_is_neg(a) then appchar('-')
    else if mp_show_plus then appchar('+');
    {insert '.', will be swapped later}
    ps1   := pstr;
    appchar(mp_fract_sep{'.'});
    ps2   := pstr;
    {convert mantissa}
    s_mp_toradix_n(x.mantissa, radix, false, maxlen, pstr);
    {swap '.' with first digit}
    if ps2^<>#0 then begin
      ps1^ := ps2^;
      ps2^ := mp_fract_sep{'.'};
    end;
    {append exponent}
    for i:=1 to length(t) do appchar(t[i]);
    pstr^ := #0;
    dec(maxlen);
  end;
  mpf_clear2(x,y);
end;


{---------------------------------------------------------------------------}
procedure s_mpf_toradix_alt(const a: mp_float; radix, ndd: word; var maxlen: longint; var pstr: pchar8);
  {-Convert an mp_float to radix alternative representation, ndd>0 is the total}
  { number of digits, trailing '.' or '0' are suppressed. If a is too large or}
  { too small, scientific representation is used. NOTE: no radix prefix/suffix!}
var
  x: mp_float;
  f,p,q: mp_int;
  k,e: longint;
  lra: double;
  prmap: ^TRadixCMap;
  d,rsqr: mp_digit;
  n: integer;
  psav: pchar8;
  reven: boolean;

  procedure appchar(c: char8);
    {-Append char at pstr^}
  begin
    if maxlen>1 then begin
      pstr^ := c;
      inc(pstr);
      dec(maxlen);
    end;
  end;

begin

  if mp_error<>MP_OKAY then exit;
  if (radix < 2) or (radix > MAXRadix) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('s_mpf_toradix_alt: radix out of range or maxlen too small');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  pstr^ := #0;
  if ndd=0 then exit;

  if s_mpf_is0(a) then begin
    {special case a=0}
    if maxlen>ndd+2+ord(mp_show_plus) then begin
      appchar('0');
      pstr^ := #0;
      dec(maxlen);
    end;
    exit;
  end;

  {Basic principle of the alternative output format: With a suitable k}
  {compute m = round(x*radix^k), p = m div radix^k, q = m div radix^k.}
  {Output the integer part p, and if q<>0 then output '.' and q. }

  {Get raw binary exponent e and approximate logbase(radix,a)}
  e := s_mpf_ldx(a);
  if add32_ovr(a.exponent,mp_bitsize(a.mantissa),e) then begin
    {set lra to force scientific format, (and avoid warning 'k not initialized')}
    lra := ndd;
    k := 0;
  end
  else begin
    if radix=2 then lra := e else lra:= e*lograd[radix];
    if lra<0 then e:=0 else e := trunc(lra);
    {Get suitable radix exponent k for conversion}
    k := ndd-e-1;
  end;
  if (lra>=ndd) or (lra>=maxlen) or ((lra+ndd)<1) or (k<0) then begin
    {Use scientific format if alternative format not applicable}
    s_mpf_toradix_n(a, radix, ndd, maxlen, pstr);
    exit;
  end;

  mpf_initp(x,a.bitprec + PINC);
  if mp_error<>0 then exit;
  mp_init3(f,p,q);
  if mp_error<>0 then begin
    mpf_clear(x);
    exit;
  end;

  {Multiply with radix^k, make positive, and round to nearest}
  mp_set_pow(f,radix,k);
  mpf_mul_mpi(a,f,x);
  s_mpf_abs(x);
  mpf_round(x,x.mantissa);

  {Divide by radix^k and process quotient and remainder}
  mp_divrem(x.mantissa,f,@p,@q);

  {Insert sign of input}
  if s_mpf_is_neg(a) then appchar('-')
  else if mp_show_plus then appchar('+');

  {Convert integer part}
  s_mp_toradix_n(p, radix, false, maxlen, pstr);

  {Convert fractional part only if non-zero.}
  if (maxlen>1) and (not mp_is0(q)) then begin
    if maxlen=2 then begin
      {Only fract_sep and #0 can be inserted}
      appchar(mp_fract_sep)
    end
    else begin
      if maxlen>k+1 then begin
        {save fract_sep position in psav}
        psav := pstr;
        {Add radix^k, this will result in '1' at fract_sep position}
        mp_add(q,f,q);

        {remove trailing zeros}
        if IsPow2_w(radix, n) then begin
          {use shifts if radix is power of two}
          e := mp_cnt_lsb(q);
          if e>0 then begin
            e := (e div n)*n;
            if e>0 then mp_shr(q, e, q);
          end;
        end
        else begin
          reven := radix and 1 = 0;
          {try to remove pairs of zeros first, i.e. divide by radix^2}
          if mp_radexp[radix]=1 then rsqr := radix else rsqr := sqr(radix);
          {first loop using max. power of radix that fits into an mp_digit}
          repeat
            {loop invariant: q<>0. Done if radix is even and q is odd}
            if reven and odd(q.pdigits^[0]) then break;
            mp_div_d(q,rsqr,@f,d);
            if d=0 then mp_exch(q,f)
            else begin
              {d = q mod rsqr is nonzero, do final division only if d=0 mod radix}
              if (d mod radix <> 0) or (reven and odd(q.pdigits^[0])) then break;
              mp_div_d(q,radix,@f,d);
              if d=0 then mp_exch(q,f)
              else break;
            end;
          until false;
        end;

        s_mp_toradix_n(q, radix, false, maxlen, pstr);
        {$ifdef MPC_USE_Assert}
          assert(psav^='1','psav^ = ''1'' in s_mpf_toradix_alt');
        {$endif}
        {replace the '1' at fract_sep position}
        psav^ := mp_fract_sep;
      end
      else begin
        {Here maxlen <= k+1 and the fractional part will be truncated. Use slow}
        {double division loop to get the correct truncated sequence of digits}
        appchar(mp_fract_sep);
        {local pointer to radix map}
        if mp_uppercase or (radix>36) then prmap := @mp_ucrmap else prmap := @mp_lcrmap;
        {Calculate the quotients q = q mod radix^(i-1) and the remainders}
        {p_i = q div radix^(i-1) with 0 <= p_i < radix, i=k...}
        while (maxlen>1) and (not mp_is0(q)) do begin
          mp_div_d(f,radix,@f,d);
          mp_divrem(q,f,@p,@q);
          if p.used>0 then begin
            d := p.pdigits^[0];
            {$ifdef MPC_USE_Assert}
              assert((p.used=1) and (d<radix),'d<radix in s_mpf_toradix_alt');
            {$endif}
          end;
          appchar(prmap^[d]);
        end;
      end
    end;
  end;
  pstr^ := #0;
  dec(maxlen);
  mpf_clear(x);
  mp_clear3(f,p,q);
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

type
  t_mpfunc = procedure(const a: mp_float; var b: mp_float);


{---------------------------------------------------------------------------}
procedure inv_func(f: t_mpfunc; const a: mp_float; var b: mp_float);
  {-Calculate b := f(1/a) with PINC guard bits}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
   {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b);
  {$endif}
  mpf_initp(t, b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_inv(a,t);
    f(t,t);
    mpf_copyp(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure func_inv(f: t_mpfunc; const a: mp_float; var b: mp_float);
  {-Calculate b := 1/f(a) with PINC guard bits}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
   {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b);
  {$endif}
  mpf_initp(t, b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    f(a,t);
    mpf_inv(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_cot(const a: mp_float; var b: mp_float);
  {-Calculate the circular cotangent b := cot(a), a mod Pi <> 0}
begin
  func_inv({$ifdef FPC_ProcVar}@{$endif}mpf_tan,a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_cotd(const a: mp_float; var b: mp_float);
  {-Calculate b = cot(a), a in degrees}
begin
  mpf_deg2rad(a,b);
  mpf_cot(b,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_sec(const a: mp_float; var b: mp_float);
  {-Calculate the circular secant b = sec(a), a mod Pi <> Pi/2}
begin
  func_inv({$ifdef FPC_ProcVar}@{$endif}mpf_cos,a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_csc(const a: mp_float; var b: mp_float);
  {-Calculate the circular cosecant b = csc(a), a mod Pi <> 0}
begin
  func_inv({$ifdef FPC_ProcVar}@{$endif}mpf_sin,a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_coth(const a: mp_float; var b: mp_float);
  {-Calculate the hyperbolic cotangent b = coth(a), a <> 0}
begin
  func_inv({$ifdef FPC_ProcVar}@{$endif}mpf_tanh,a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_sech(const a: mp_float; var b: mp_float);
  {-Calculate the hyperbolic secant b = sech(a)}
begin
  func_inv({$ifdef FPC_ProcVar}@{$endif}mpf_cosh,a,b);
end;



{---------------------------------------------------------------------------}
procedure mpf_csch(const a: mp_float; var b: mp_float);
  {-Calculate the hyperbolic cosecant b = csch(a), a <> 0}
begin
  func_inv({$ifdef FPC_ProcVar}@{$endif}mpf_sinh,a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_cospi(const a: mp_float; var b: mp_float);
  {-Calculate b = cos(Pi*a)}
begin
  mpf_trig_ex(a, true, @b, nil, nil);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccot(const a: mp_float; var b: mp_float);
  {-Calculate the sign symmetric inverse circular cotangent; b = arccot(a) = arctan(1/a), a <> 0}
var
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_arccot');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  neg := s_mpf_is_neg(a);
  if s_mpf_ldx(a) > -(b.bitprec + PINC) then begin
    mpf_inv(a,b);
    mpf_arctan(b,b);
  end
  else begin
    mpf_set_pi2k(b,-1);
    if neg then s_mpf_chs(b);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_arccotc(const a: mp_float; var b: mp_float);
  {-Calculate the continuous inverse circular cotangent; b = arccotc(a) = Pi/2 - arctan(a)}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
   {$ifdef MPC_ArgCheck}
    _CheckBitPrec(b);
  {$endif}
  if s_mpf_is_neg(a) then begin
    mpf_initp(t, b.bitprec);
    if mp_error=MP_OKAY then begin
      mpf_arctan(a,b);
      mpf_set_pi2k(t,-1);
      mpf_sub(t,b,b);
      mpf_clear(t);
    end;
  end
  else mpf_arccot(a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccsc(const a: mp_float; var b: mp_float);
  {-Calculate the inverse circular cosecant b = arccsc(a), |a| >= 1}
begin
  inv_func({$ifdef FPC_ProcVar}@{$endif}mpf_arcsin, a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arcsec(const a: mp_float; var b: mp_float);
  {-Calculate the inverse circular secant b = arcsec(a), |a| >= 1}
begin
  inv_func({$ifdef FPC_ProcVar}@{$endif}mpf_arccos,a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccoth(const a: mp_float; var b: mp_float);
  {-Calculate the inverse hyperbolic cotangent b = arccoth(a), |a| > 1}
begin
  inv_func({$ifdef FPC_ProcVar}@{$endif}mpf_arctanh,a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arccsch(const a: mp_float; var b: mp_float);
  {-Calculate the inverse hyperbolic cosecant b = arccsch(a)}
begin
  inv_func({$ifdef FPC_ProcVar}@{$endif}mpf_arcsinh,a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_arcsech(const a: mp_float; var b: mp_float);
  {-Calculate the inverse hyperbolic secant b = arcsech(a), 0 < a <= 1}
begin
  inv_func({$ifdef FPC_ProcVar}@{$endif}mpf_arccosh,a,b);
end;


{---------------------------------------------------------------------------}
procedure mpf_archav(const a: mp_float; var b: mp_float);
  {-Calculate the inverse haversine b = archav(a) = 2*arcsin(sqrt(a))}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_archav');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end
  else if mpf_is1(a) then begin
    mpf_set_pi(b);
    exit;
  end
  else if s_mpf_is_neg(a) or (s_mpf_ldx(a)>0) then begin
    {a must be in [0..1]}
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpf_archav: a < 0 or a > 1');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mpf_initp(t,b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_sqrt(a,t);
    mpf_arcsin(t,t);
    mpf_mul_2k(t,1,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_gd(const a: mp_float; var b: mp_float);
  {-Calculate the Gudermannian function b = gd(a) = arctan(sinh(a))}
var
  t: mp_float;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_gd');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mp_error<>MP_OKAY then exit;
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  mpf_initp(t,b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    if s_mpf_ldx(a) < 1.45*ln(t.bitprec) then begin
      {use gd(a) = arctan(sinh(a))}
      mpf_sinh(a,t);
      mpf_arctan(t,t);
      mpf_copyp(t,b);
    end
    else begin
      neg := s_mpf_is_neg(a);
      mpf_set_pi2k(b,-1);
      if neg then s_mpf_chs(b);
    end;
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_arcgd(const a: mp_float; var b: mp_float);
  {-Calculate the inverse Gudermannian function b = arcgd(a) = arcsinh(tan(a))}
var
  t: mp_float;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_arcgd');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  mpf_initp(t,b.bitprec + PINC);
  if mp_error=MP_OKAY then begin
    mpf_tan(a,t);
    mpf_arcsinh(t,t);
    mpf_copyp(t,b);
    mpf_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
function lw_app(x: double): double;
  {-Return initial approximation for Lambert W function (principal branch), x >= -1/e}
const
  em1h = 0.36785888671875;           {high part of 1/e}
  em1l = 0.205544526923215955238e-4; {low  part of 1/e}
  xlow =-0.36767578125;              {threshold for branch-point series}
const
  PC: array[0..11] of double = ( {Calculated with MPArith:}
           {Branch point series W(x) = sum(a[i]*y^i, i=0..), y=sqrt(x+1/e),}
           {a[i] = mu[i]*(sqrt(2e))^i, with mu[i] from [37] (4.23/4).}
           -1.0000000000000000000,
            2.3316439815971242034,
           -1.8121878856393634902,
            1.9366311144923597554,
           -2.3535512018816145169,
            3.0668589010506319128,
           -4.1753356002581771388,
            5.8580237298747741487,
           -8.4010322175239773709,
            12.250753501314460424,
           -18.100697012472442756,
            27.029044799010561650);
const
  p1 = 19.0/10.0;
  p2 = 17.0/60.0;
  q1 = 29.0/10.0;
  q2 = 101.0/60.0;
const
  eps_d: double = 2.2204460492503131E-16;
var
  q,w,z: double;
  i: integer;
begin
  if abs(x)<=5e-9 then begin
    {LambertW = x*(1-x+1.5*x^2+O(x^3), OK but main routine use mpf}
    w := x*(1.0-x);
  end
  else if x <= xlow then begin
    z := (x+em1h)+em1l;
    if abs(z) <= eps_d then begin
      {Tolerate some rounding error, this will give W(-exp(-1)) = -1}
      w := -1.0;
    end
    else if z<0 then w := -2
    else begin
      z := sqrt(z);
      w := PC[11];
      for i:=10 downto 0 do w := w*z + PC[i];
    end;
  end
  else begin
    {Initial approximations for iteration}
    if x<3.0 then begin
      {pade(LambertW(x), x, [3,2])}
      q := (p2*x+p1)*x +1.0;
      z := (q2*x+q1)*x +1.0;
      w := x*(q/z);
    end
    else begin
      {Asymptotic series, see [37] (4.19)}
      z := ln(x);
      q := ln(z);
      w := z - q + q/z;
      if q>3.0 then w := w + 0.5*q*(q-2.0)/sqr(z);
    end;
  end;
  lw_app := w;
end;


{---------------------------------------------------------------------------}
procedure mpf_lambertw(const a: mp_float; var b: mp_float);
  {-Calculate b = LambertW(a); principal branch, a >= -1/e}
var
  e,q,w,x,z: mp_float;
  s,t: double;
  i: integer;
  bsa,pe,pw,pz: longint;
const
  IMAX=200;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpf_lambertw');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mpf_is0(a) then begin
    mpf_set0(b);
    exit;
  end;
  mpf_initp5(e,q,w,x,z,b.bitprec + PINC);
  pw := w.bitprec;
  mpf_copyp(a,x);

  {calculate first approximation w and starting precision}
  bsa := s_mpf_ldx(a);
  pz  := -(bsa div 2);
  if pz<16 then pz := 16;

  if (bsa <= 1020) or s_mpf_is_neg(a) then begin
    {Argument is in double range}
    if bsa<-30 then begin
      {w := a - a^2 + 3/2*a^3}
      mpf_sqr(x,q);
      mpf_mul(q,x,w);
      mpf_mul_dbl(z,1.5,w);
      mpf_sub(w,q,w);
      mpf_add(w,a,w);
      t  := 0;
    end
    else begin
      {no mpf needed for starting value}
      s := mpf_todouble(a);
      t := lw_app(s);
      if t=-2.0 then begin
        {t value indicates s < -1/e}
        mpf_clear5(e,q,w,x,z);
        {$ifdef MPC_HaltOnError}
          {$ifdef MPC_UseExceptions}
            raise MPXRange.Create('mpf_lambertw: a < -1/e');
          {$else}
            RunError(MP_RTE_RANGE);
          {$endif}
        {$else}
          set_mp_error(MP_RANGE);
          exit;
        {$endif}
      end
      else if (t=-1.0) or (s<-0.36787944117) then begin
        {full precision from start if near branch point}
        pz := pw;
      end;
      mpf_set_dbl(w,t);
    end;
    if t=-1.0 then begin
      {x very close to -1/e, use the first terms of branch point series}
      {w = -1 + p - p^2/3 + 11/72*p^3, p = sqrt(2(e*x+1)), [37] (4.22)}
      {w := 2(e*x+1) = p^2}
      mpf_exp(w,w);
      mpf_div(a,w,w);
      s_mpf_inc1(w);
      mpf_mul_2k(w,1,w);
      if s_mpf_is_neg(w) then mpf_set_dbl(w,-1)
      else begin
        {z = p}
        mpf_sqrt(w,z);
        mpf_div_dbl(w,-3.0,w);
        mpf_add(w,z,w);
        mpf_sub_dbl(w,1.0,w);
      end;
    end;
  end
  else begin
    {a > 2^1000, use asymptotic form as initial approximation}
    {w = z - q + q/z + 0.5*q*(q-2)/z^2  with  z=ln(x), q=ln(z)}
    mpf_ln(x,z);
    mpf_ln(z,q);
    mpf_div(q,z,w);
    {e := 0.5*q*(q-2.0)/sqr(z) = (q/z)*(q-2.0)/z*0.5}
    mpf_sub_dbl(q,2.0,e);
    mpf_mul(e,w,e);
    mpf_div(e,z,e);
    mpf_mul_2k(e,-1,e);
    {w  = q/z + z - q + e}
    mpf_add(w,z,w);
    mpf_sub(w,q,w);
    mpf_add(w,e,w);
  end;

  for i:=1 to IMAX do begin
    {e := 1.0 + w;}
    s_mpf_add1(w,e);
    if mpf_is0(e) then break;
    {Increase working precision}
    if pz<pw then begin
      pz := 2*pz;
      if pz>pw then pz := pw;
      z.bitprec := pz;
    end;
    {z := ln(x/w) - w;}
    mpf_div(x,w,z);
    mpf_ln(z,z);
    mpf_sub(z,w,z);
    {q := 2.0*e*(e + z/1.5) - z;}
    mpf_div_dbl(z,1.5,q);
    mpf_add(q,e,q);
    mpf_mul(q,w,q);
    mpf_mul_2k(q,1,q);
    mpf_sub(q,z,q);
    {e := (z/e)*(q/(q-z));}
    mpf_div(z,e,e);
    mpf_mul(e,q,e);
    mpf_sub(q,z,q);
    mpf_div(e,q,e);
    {w := w*(1.0+e);}
    s_mpf_add1(e,q);
    mpf_mul(w,q,w);
    pe := s_mpf_ldx(e);
    if (pz>=pw) and (pe+pw <= 8) then break;
  end;

  mpf_copyp(w,b);
  mpf_clear5(e,q,w,x,z);
end;



{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

{---------------------------------------------------------------------------}
procedure s_mpf_expnewt(const a: mp_float; var b: mp_float);
  {-Calculate b = exp(a), a < 2^31 * ln(2), !!!TEST ONLY!!!}
var
  x,t: mp_float;
  i,bp0: longint;
const
  PIT  = 4000;       {approx. break-even precision standard / Newton exp}
  PITF = PIT div 2;
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mpf_not_init(a) or mpf_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mpf_expnewt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if s_mpf_is0(a) then begin
    mpf_set1(b);
    exit;
  end;
  bp0 := b.bitprec + PINC;

  if bp0 > MPF_MAX_PREC then bp0 := MPF_MAX_PREC;
  if bp0 <= PIT then begin
    mpf_exp(a,b);
    exit;
  end;
  i := bp0;
  while i>PITF do i := succ(i div 2);
  i := i + 4;
  mpf_initp2(x,t,i);
  if mp_error<>MP_OKAY then exit;
  {Get b = exp(a) as the root of f(x) = ln(x) - a = 0 using Newton}
  {iteration x_next = x - f(x)/f'(x) = x + x*(a - ln(x)), initial }
  {value x = exp(a) with initial bitprec less than PIT/2}
  mpf_exp(a,x);
  repeat
    if i<bp0 then begin
      i := 2*i;
      if i>bp0 then i := bp0;
    end;
    t.bitprec := i;
    s_mpf_normalizep(x,i);
    s_mpf_lnagm(x,t);
    mpf_sub(a,t,t);
    mpf_mul(t,x,t);
    mpf_add(t,x,x);
  until (i=bp0) or (mp_error<>MP_OKAY);
  mpf_copyp(x,b);
  mpf_clear2(x,t);
end;



{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

begin
  mpf_set_default_prec(240);
end.

