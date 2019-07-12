unit mp_base;

{Multi precision integer arithmetic basic routines}

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

 DESCRIPTION   :  Multi precision integer arithmetic basic routines

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA :  (mp_types)

 MEMORY USAGE  :  lots of heap

 DISPLAY MODE  :  ---

 REFERENCES    :  [1] LibTomMath V0.30+ by Tom St Denis
                  [2] MPI by M.J. Fromberger
                  [3] D.E. Knuth, The Art of computer programming:
                      Volume 1, Fundamental Algorithms, 3rd ed., 1997;
                      Volume 2, Seminumerical Algorithms, 3rd ed., 1998;
                      http://www-cs-faculty.stanford.edu/~knuth/taocp.html
                  [4] Forster, O.: Algorithmische Zahlentheorie, 1996
                  [5] (HAC) Menezes,A., von Oorschot,P., Vanstone, S: Handbook of
                      Applied Cryptography, 1996, www.cacr.math.uwaterloo.ca/hac
                  [8] Marcel Martin: NX - Numerics library of multiprecision
                      numbers for Delphi and Free Pascal, 2006-2009
                      www.ellipsa.eu/public/nx/index.html
                 [10] Crandall,R., C.Pomerance: Prime Numbers, A Computational
                      Perspective, 2nd ed., 2005
                 [15] The GNU Multiple Precision Arithmetic Library, http://gmplib.org/
                 [24] H. Cohen, A Course in Computational Algebraic Number Theory
                      4th printing, 2000
                 [29] V. Shoup, A Computational Introduction to Number Theory and
                      Algebra, Version 2, 2008, from http://shoup.net/ntb/
                 [30] J. v. zur Gathen, J. Gerhard, Modern computer algebra, 2nd ed., 2003
                      http://math-www.uni-paderborn.de/mca/
                 [33] C. Burnikel, J. Ziegler: Fast Recursive Division. MPI fr Informatik,
                      Forschungsbericht MPI-I-98-1-022 (1998); available via
                      http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.47.565
                 [34] P. Zimmermann, Karatsuba Square Root, INRIA Research Report RR-3805;
                      available from http://hal.inria.fr/inria-00072854/en/
                 [35] R.P. Brent, P. Zimmermann: Modern Computer Arithmetic, Cambridge University Press, 2010.
                      A preliminary version (V0.5.9, Oct. 2010) of the book is available from
                      http://maths-people.anu.edu.au/~brent/pd/mca-cup-0.5.9.pdf
                      or http://arxiv.org/abs/1004.4710 (V0.5.1)
                 [36] M. Bodrato, A.Zanoni, What About Toom-Cook Matrices Optimality?
                      available from http://bodrato.it/papers/#CIVV2006
                      See also M. Bodrato's page http://bodrato.it/software/


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.0.01   09.05.04  W.Ehrhardt  Initial version: BP7
 0.0.02   10.05.04  we          mp_copyprim
 0.0.03   11.05.04  we          mp_expand, mp_copy
 0.0.04   12.05.04  we          mp_newXY, mp_delete, LSign
 0.0.05   22.07.04  we          new design a la MPI/LibTomMath
                                mp_init, mp_clear, mp_exch, mp_zero
 0.0.06   23.07.04  we          mp_shrink, mp_grow, mp_init_size, mp_set,
                                mp_init_set, mp_clamp, mp_lshd
 0.0.07   23.07.04  we          cleanup/reordering, mp_iseven/odd/zero
 0.0.08   23.07.04  we          mp_copy
 0.0.09   23.07.04  we          debug/magic, mp_reverse, mp_toradix_n
 0.0.10   24.07.04  we          mp_div_d (only a=0,b=1 quick outs)
 0.0.11   24.07.04  we          mp_div_2d
 0.0.12   24.07.04  we          mp_rshd
 0.0.13   24.07.04  we          mp_mod_2d
 0.0.14   24.07.04  we          mp_mul_d
 0.0.15   24.07.04  we          mp_set_int
 0.0.16   24.07.04  we          mp_mul_2d
 0.0.17   24.07.04  we          fixed severe bug in mp_lshd
 0.0.18   25.07.04  we          debug code checking initialized mp_ints
 0.0.19   25.07.04  we          mp_read_radix
 0.0.20   25.07.04  we          mp_add_d
 0.0.21   26.07.04  we          mp_sub_d
 0.0.22   26.07.04  we          mp_div_2, mp_mul_2

 0.1.00   26.07.04  we          const parameters, single digit functions
                                complete, code clean up
 0.1.01   27.07.04  we          mp_cmp_d, mp_cmp_mag, mp_cmp
 0.1.02   27.07.04  we          InitCheck / RunError(210)
 0.1.03   27.07.04  we          s_mp_add/sub, mp_add/sub
 0.1.04   28.07.04  we          mp_rand, bugfix: mp_mul_d, s_mp_add
 0.1.05   28.07.04  we          $Q- in s_mp_sub/mp_sub_d
 0.1.06   29.07.04  we          s_mp_mul_digs {+tbd}, mp_mul {+tbd}
 0.1.07   29.07.04  we          mp_count_bits, mp_div
 0.1.08   30.07.04  we          mp_abs, mp_chs, mp_2expt
 0.1.09   30.07.04  we          mp_sqr {*tbd}, mp_expt, mp_expt_d
 0.1.10   01.08.04  we          mp_show_plus, mp_radix_size, mp_toradix
 0.1.11   01.08.04  we          mp_uppercase, MAXRadix, InitCheck->mp_argchk
 0.1.12   01.08.04  we          RunError with MP_RTE_xx constants
 0.1.13   01.08.04  we          mp_cmp_z, mp_cmp_int
 0.1.14   01.08.04  we          fix compiler warnings, formatting cleanup
 0.1.15   01.08.04  we          mp_sqrt, mp_mod
 0.1.16   02.08.04  we          mp_karatsuba_mul, mp_mul completed
 0.1.17   02.08.04  we          Bug? in Delphi mem mgt, fixed with mp_realloc in mp_grow
 0.1.18   02.08.04  we          Bug fixed 16 bit mp_realloc, mp_abs interfaced
 0.1.19   03.08.04  we          Fix Delphi32 with signs in Karatsuba
 0.1.20   03.08.04  we          BugFix mp_div (negative zero)
 0.1.21   03.08.04  we          small change in s_mp_mul_digs
 0.1.22   04.08.04  we          conditional define: use move in Karatsuba
 0.1.23   04.08.04  we          s_mp_sqr, mp_sqr_cutoff
 0.1.24   04.08.04  we          mp_expt_d with b: word
 0.1.25   04.08.04  we          mp_clear with fillchar(a,sizeof(a),0)
 0.1.26   04.08.04  we          {$i mp_conf.inc), MP_VERSION from mp_types
 0.1.27   04.08.04  we          mp_karatsuba_sqr, mp_conf.inc removed
 0.1.28   05.08.04  we          32Bit fix s_mp_sqr
 0.1.29   06.08.04  we          Bug search cleanups marked {*0.1.29}
 0.1.30   06.08.04  we          "Karatsuba" bug fixed in s_mp_add
 0.1.31   07.08.04  we          last digit optimization in mp_expt, 0 is even

 0.2.00   08.08.04  we          "release" version 0.2

 0.2.01   08.08.04  we          Debug: check mp_digits <= MP_DIGIT_MAX
 0.2.02   08.08.04  we          mp_unsigned_bin_size, mp_radix_size with size: longint
 0.2.03   09.08.04  we          mp_to_(un)signed_bin_n, mp_read_(un=signed_bin
 0.2.04   09.08.04  we          mp_mod_d, changed mp_div_d parameter list
 0.2.05   09.08.04  we          mp_rand renamed to mp_rand_radix, new mp_rand
 0.2.04   09.08.04  we          bugfix mp_div_2d
 0.2.05   09.08.04  we          optimize mp_init_copy
 0.2.06   09.08.04  we          LTM031 change applied to mp_2expt
 0.2.07   10.08.04  we          mp_init/clear_multi(_p)
 0.2.08   11.08.04  we          easy outs for mp_expt(_d)
 0.2.09   14.08.04  we          mp_init_set_int, mp_get_int
 0.2.10   18.08.04  we          optimized mp_mod_d
 0.2.11   24.08.04  we          mp_mul_int

 0.3.00   26.08.04  we          "release" version 0.3
 0.3.01   26.08.04  we          mp_error, most functions now procedures
 0.3.02   26.08.04  we          mp_radix_size: longint, remove D6+ warnings
 0.3.03   27.08.04  we          mp_result
 0.3.04   30.08.04  we          mp_radix_size with ln
 0.3.05   31.08.04  we          mp_is_<xy> functions, mp_memused, lograd table
 0.3.06   27.02.05  we          mp_freemem in mp_clear, mp_memused with mp_memstat
 0.3.07   16.05.05  we          FPC: $mode objfpc/$goto on
 0.3.08   09.08.05  we          Bugfix mp_toradix_n if a has only one digit
 0.3.09   15.08.05  we          mp_isbit
 0.3.10   17.08.05  we          severe bug fix mp_expt(a,b,c) for @b=@c

 0.4.00   20.08.05  we          use mp_set_error
 0.4.01   20.08.05  we          some mp_digit typecasts, radix/digits word in mp_rand_radix
 0.4.02   21.08.05  we          mp_conf.inc,  mp_memused to mp_base
 0.4.03   21.08.05  we          MPC_ArgCheck, MPC_Assert
 0.4.04   22.08.05  we          removed mp_argcheck, mp_errchk
 0.4.05   22.08.05  we          mp_init<i>, mp_clear<i>, i=2..5
 0.4.06   22.08.05  we          mp_2expt with longint b
 0.4.07   22.08.05  we          bugfix mp_init5
 0.4.08   23.08.05  we          MPC_HaltOnArgCheck
 0.4.09   23.08.05  we          IsPrime32, IsSPP32, IsSPP32A
 0.4.10   24.08.05  we          use MaxDigits related codes MP_MAXDIGITS, MP_RTE_OTHER
 0.4.11   24.08.05  we          fix IsPrime32 for BP7/{$R+,Q+}
 0.4.12   24.08.05  we          fix mp_init_set/_int
 0.4.13   26.08.05  we          uses mp_prng (mp_rand, mp_rand_radix)
 0.4.14   27.08.05  we          argcheck for digit in some mp_xxx_d functions
 0.4.15   27.08.05  we          functions mp_cmp_mag_d, mp_read_decimal
 0.4.16   28.08.05  we          usage of exceptions implemented
 0.4.17   29.08.05  we          mp_set_short, mp_isone, mp_is0, mp_is1
 0.4.18   07.09.05  we          mp_n_root, bugfix mp_mul_d
 0.4.19   08.09.05  we          some functions moved from mp_supp
 0.4.20   10.09.05  we          mp_reduce_2k functions
 0.4.21   10.09.05  we          mp_core routines integrated
 0.4.22   10.09.05  we          optimized mp_reduce_is_2k
 0.4.23   16.09.05  we          optimized mp_cmp_mag and s_mp_sub_d
 0.4.24   18.09.05  we          mp_let (alternative for mp_read_decimal)
 0.4.25   21.09.05  we          use mp_clearzero
 0.4.26   21.09.05  we          $argcheck pc<>pd in mp_div

 0.5.00   29.09.05  we          'internal' functions moved to end of interface
 0.5.01   30.09.05  we          IsPrime16, pbits16, pmask16
 0.5.02   02.10.05  we          changed 'SPP' to more conventional 'spsp'
 0.5.03   07.10.05  we          is_spsp32A: optimized BASM16
 0.5.04   08.10.05  we          more is_spsp32A optimization
 0.5.05   13.10.05  we          mp_isbit with BASM16
 0.5.06   14.10.05  we          BugFix is_spsp32A: reduce a := bases[k] mod N in
 0.5.07   15.10.05  we          BugFix is_spsp32A: no mod N, use MulMod32 or div
 0.5.08   15.10.05  we          BugFix N>=9080191 in IsPrime32
 0.5.09   20.11.05  we          ArgCheck in mp_isbit, new name: mp_montgomery_calcnorm
 0.5.10   21.11.05  we          mp_rand_bits, mp_rand: check digits<MAXDigits

 0.6.00   30.12.05  we          mp_count_bits/CountBits32 renamed to mp_bitsize/bitsize32
 0.6.01   30.12.05  we          MP_8BIT removed
 0.6.02   31.12.05  we          mp_set_int via array[0..3] of byte
 0.6.03   31.12.05  we          mp_mul_w, mp_set_w
 0.6.04   31.12.05  we          popcount16/32, mp_popcount
 0.6.05   10.01.06  we          Halt on BAD_ARG in mp_montgomery_setup
 0.6.06   28.01.06  we          s_mp_add_d; changes in mp_sub_d, s_mp_sqr, s_mp_chs

 0.7.00   19.03.06  we          mp_div_w
 0.7.01   04.08.06  we          Bugfix in mp_mod_2d, improve mp_reduce_2k
 0.7.02   08.08.06  we          mp_makeodd
 0.7.03   09.08.06  we          mp_clear6, mp_init6; removed mp_let, mp_isone
 0.7.04   10.08.06  we          mp_inc, mp_dec
 0.7.05   11.08.06  we          bugfix mp_set_short
 0.7.06   11.08.06  we          avoid FPC warnings: mp_set_int/mp_set_w
 0.7.07   11.08.06  we          fixed and improved mp_reduce
 0.7.08   12.08.06  we          s_mp_mul_high_digs: error if id<0
 0.7.09   13.08.06  we          mp_set_pow, mp_expt_int, mp_expt_d uses mp_expt_int
 0.7.10   15.08.06  we          mp_[x]_int with [x]: add,dec,div,inc,mod,sub
 0.7.11   26.08.06  we          rewrite is_spsp32A, new internal function _spsp32
 0.7.12   27.08.06  we          BIT32: bigalloc; FPC: use ReturnNilIfGrowHeapFails
 0.7.13   28.08.06  we          mp_n_root with longint parameter
 0.7.14   28.08.06  we          mp_clrbit, mp_setbit, mp_2expt uses mp_setbit
 0.7.15   30.08.06  we          mp_gr_mod, mp_gr_setup
 0.7.16   30.08.06  we          mp_shr, mp_shl, mp_div_2d(..,nil) replaced by mp_shr
 0.7.17   06.09.06  we          better initial approximation for mp_sqrt
 0.7.18   07.09.06  we          mp_mod_w

 0.8.00   18.09.06  we          mp_alloc interfaced

 0.9.00   26.12.06  we          some minor changes related to mp_clamp/assert
 0.9.01   27.12.06  we          mp_lshd2
 0.9.02   27.12.06  we          s_mp_div: Knuth's q calculation from Alg. D
 0.9.03   28.12.06  we          s_mp_div: separate treatment of single digit b
 0.9.04   29.12.06  we          $ifdef MPC_USE_Assert
 0.9.05   01.01.07  we          mp_prod_int, bigalloc renamed to IAlloc
 0.9.06   01.01.07  we          mp_mul_int optimized for small longints
 0.9.07   01.01.07  we          mp_mul optimized for single digit factors
 0.9.08   02.01.07  we          mp_toradix10_n, speed up 2.5 .. 4 for mp_toradix(..,10,..)
 0.9.09   03.01.07  we          changed mp_read_radix to accept uppercase and lowercase
 0.9.10   03.01.07  we          new mp_toradix_n (generalization of mp_toradix10_n)
 0.9.11   03.01.07  we          mp_toradix_n: local TRadixCMap
 0.9.12   04.01.07  we          bugfixed/improved mp_is_power_of_two, renamed to mp_is_pow2_d
 0.9.13   04.01.07  we          mp_is_pow2
 0.9.14   07.01.07  we          improved mp_read_radix

 1.0.00   11.04.07  we          mp_init_prim: use mp_precision if size=0
 1.0.01   11.04.07  we          mp_mul_d/w with factor w/d=0/1
 1.0.02   12.04.07  we          s_mp_read_radix
 1.0.03   12.04.07  we          Bugfix EstimateQDigit, improved normalization in s_mp_div
 1.0.04   13.04.07  we          Easy outs in mp_expt_int
 1.0.05   14.04.07  we          s_mp_toradix_n, off by 1 bugfix mp_radix_astr
 1.0.06   15.04.07  we          s_mp_write_radix
 1.0.07   16.04.07  we          mp_is1 without mp_cmp_d
 1.0.08   01.05.07  we          mp_todouble, ldexpd
 1.0.09   01.05.07  we          DblPosInf,DblNegInf,DblNaN; improved mp_todouble
 1.0.10   05.05.07  we          s_mp_read_radix with SignAllowed, bugfix mp_set_w/int
 1.0.11   05.05.07  we          rewrite mp_set_w; frexpd, bugfix mp_todouble
 1.0.12   07.05.07  we          0^0=1 in mp_expt
 1.0.13   11.05.07  we          Removed MulMod32, mp_shrink
 1.0.14   11.05.07  we          Bugfix s_mp_write_radix
 1.0.15   13.05.07  we          Corrected some exception strings
 1.0.16   13.05.07  we          MPAF prefix in assert strings

 1.1.00   27.06.07  we          mp_abs: Arg check done in mp_copy
 1.1.01   01.07.07  we          EstimateQDigit: removed = from q>=MP_MASK test
 1.1.02   10.07.07  we          mp_writeln, mp_clear[x]/mp_init[x] (x=7..9)
 1.1.03   15.07.07  we          mp_reduce: easy out if x<m, allow x<0
 1.1.04   21.07.07  we          improved mp_shl, easy out in mp_shr
 1.1.05   21.07.07  we          mp_sqrt: Initial double approximation based on highest mp_digit(s)
 1.1.06   26.07.07  we          isqrt32
 1.1.07   26.07.07  we          isqrt32 with FPU
 1.1.08   26.07.07  we          mp_sqrt: new recursive integer square root algorithm
 1.1.09   29.07.07  we          improved: mp_reduce, s_mp_mul_high_digs, s_mp_mul_digs
 1.1.10   04.08.07  we          renamed/new: s_mp_mod_2d, mp_mod_2d

 1.2.00   17.08.07  we          gcd32/gcd32u
 1.2.01   19.08.07  we          changed mp_mod_int to use mp_mod
 1.2.02   19.08.07  we          Bugfix mp_mod: don't adjust sign if result=0
 1.2.03   19.08.07  we          mp_mod_int without temporary mp_int, local BASM in loop
 1.2.04   26.08.07  we          gcd32/gcd32u call internal __gcd32
 1.2.05   02.09.07  we          new s_mp_mul_int used in mp_prod_int
 1.2.06   04.09.07  we          mp_rand_bits_ex
 1.2.07   05.09.07  we          s_mp_expt_dl, s_mp_expt_wl, changed mp_set_pow
 1.2.08   05.09.07  we          MP_32BIT/MP_16BIT versions of mp_set_int
 1.2.09   06.09.07  we          popcount16 and popcount32 return integer
 1.2.10   07.09.07  we          mp_checksum, s_mp_checksum
 1.2.11   10.09.07  we          use MP_INV_MASK to avoid warnings if DIGIT_BIT=31
 1.2.12   11.09.07  we          Bugfix in BIT16 version of _spsp32
 1.2.13   17.09.07  we          IAlloc, mp_alloc, mp_freemem inline $ifdef HAS_INLINE
 1.2.14   17.09.07  we          mp_init_multi/_p use mp_init_prim, mp_abs inline
 1.2.15   17.09.07  we          Bugfix s_mp_read_radix for DIGIT_BIT<10
 1.2.16   20.09.07  we          Removed inline (D9 bug with mp_freemem(pointer(pchar)...)
 1.2.17   21.09.07  we          mp_sqrt uses at most one local mp_int

 1.3.00   03.11.07  we          mp_toextended, frexpx, ldexpx
 1.3.01   05.11.07  we          mp_todouble_ex, mp_toextended_ex
 1.3.02   11.11.07  we          s_mp_read_radix: sep now a string
 1.3.03   14.11.07  we          mp_radix_astr: prefill result with #0
 1.3.04   19.11.07  we          merged routines from mp_supp
 1.3.05   22.11.07  we          complete rewrite of mp_to_unsigned_bin_n,
                                mp_read_unsigned_bin, and mp_rand
 1.3.06   26.11.07  we          add32_ovr
 1.3.07   26.11.07  we          fix mask/bit logic in mp_rand_bits_ex
 1.3.08   29.11.07  we          Removed some arg checks (done in mp_copy)
 1.3.09   09.12.07  we          s_mp_add_ovr -> add32_ovr, undo mp_supp merging
 1.3.10   17.12.07  we          moved lograd table to interface

 1.5.00   24.01.08  we          mp_and, mp_or, mp_xor
 1.5.01   31.01.08  we          {$x+} for VP and D1

 1.6.00   23.05.08  we          Optimized Argcheck in: mp_mod, mp_cmp_int, mp_todouble_ex, mp_toextended_ex
 1.6.01   24.05.08  we          mp_not_init_multi
 1.6.02   25.05.08  we          mp_hex/mp_ahex
 1.6.03   03.06.08  we          mp_is1a, bugfix mp_is1
 1.6.04   06.06.08  we          mp_is? routines return false if mp_error<>MP_OKAY
 1.6.05   08.06.08  we          improved carry propagation in s_mp_add_d,s_mp_sub_d
 1.6.06   10.06.08  we          fix mp_mod_d for a<0
 1.6.07   11.06.08  we          mp_init_prim: use size=4 if size=mp_allocprec=0

 1.7.00   23.08.08  we          Avoid FPC222 warning in isqrt32
 1.7.01   24.08.08  we          s_mp_toradix_n improved if radix is power of 2
 1.7.02   14.09.08  we          new function IsPow2_w
                                MP_16BIT: improved mp_div_d, mp_div_w, mp_mod, s_mp_div
 1.7.03   14.09.08  we          Fix FPC RTE 201 if R+ in mp_read_unsigned_bin
 1.7.04   15.09.08  we          mp_get/set_allocprec, mp_allocprec local,
                                mp_init_prim rounds up to multiple of mp_allocprec}
 1.7.05   17.09.08  we          mp_rand_ex, mp_is_longint; BIT16: (f)LeftShiftAdd
 1.7.06   17.09.08  we          mp_mod with mp_is_longint, BASM16 for 'shr DIGIT_BIT'
 1.7.07   21.09.08  we          mp_shr1
 1.7.08   24.09.08  we          mp_sign, removed mp_cmp_z, renamed ??_2d to ??_2k
 1.7.09   24.09.08  we          string replaced by mp_string
 1.7.10   26.09.08  we          mp_read_radix_str, mp_read_decimal_str
 1.7.11   27.09.08  we          invmod32
 1.7.12   02.10.08  we          improved mp_reduce

 1.8.00   04.10.08  we          BASM16 in s_mp_mul_(high)_digs
 1.8.01   04.10.08  we          BASM16 in s_mp_sqr
 1.8.02   04.10.08  we          BASM16 in mp_montgomery_reduce
 1.8.03   05.10.08  we          BASM16 in mp_shl, mp_shr
 1.8.04   05.10.08  we          BASM16 in mp_div_w and new s_mp_div_d, fixed IsPow2_w
 1.8.05   06.10.08  we          Simplified MP_32BIT power of two code in mp_div_w
 1.8.06   07.10.08  we          Check for power of two in mp_mod_int
 1.8.07   09.10.08  we          mp_set1
 1.8.08   11.10.08  we          Improved mp_reduce_2k_setup
 1.8.09   18.10.08  we          Improved mp_expt_int, mp_popcount, mp_gr_setup
 1.8.10   19.10.08  we          s_mp_mod_w
 1.8.11   23.10.08  we          mp_n_root with Halley's iteration or bisection method
 1.8.12   24.10.08  we          mp_n_root: check startup convergence of Halley steps
 1.8.13   25.10.08  we          s_mp_n_root2, mp_n_root2
 1.8.14   26.10.08  we          Fix check for exact root in Halley
 1.8.15   27.10.08  we          Halley: optimized to keep remainder
 1.8.16   31.10.08  we          Halley: bisection count = bitsize32(d)
 1.8.17   01.11.08  we          mp_mod_int: check if a is longint
 1.8.18   04.11.08  we          s_mp_n_root2: b=1 if n>=mp_bitsize(a)

 1.9.00   08.11.08  we          mp_rshd2, avoid some warnings
 1.9.01   08.11.08  we          Halley renamed to iroot, improved initial approximation
 1.9.02   30.11.08  we          mp_read_radix_arr
 1.9.03   02.12.08  we          Uses BTypes: char8, pchar8
 1.9.04   06.12.08  we          s_mp_is_le0
 1.9.05   26.12.08  we          IsPrime16/32,is_spsp32,is_spsp32A move to mp_prime
 1.9.06   02.01.09  we          s_mp_sqrtrem, mp_sqrtrem, renamed mp_sqrt to s_mp_sqrt
 1.9.07   03.01.09  we          mp_n_root2 with mp_sqrtrem
 1.9.08   06.01.09  we          improved mp_reduce_2k
 1.9.09   06.01.09  we          Replaced CHAR_BIT

 1.10.00  06.01.09  we          skip final sqr in mp_expt
 1.10.01  21.01.09  we          mp_shl1, renamed (s)mp_div to (s)mp_divrem, new mp_div
 1.10.02  25.01.09  we          improved mp_add/dec/inc/sub_int
 1.10.03  01.02.09  we          s_mp_ln, s_mp_log2, s_mp_set_ext
 1.10.04  04.02.09  we          s_mp_mod_is0
 1.10.05  16.02.09  we          mp_shlx, mp_shrx
 1.10.06  18.02.09  we          mp_divrem_newton

 1.11.00  07.03.09  we          s_mp_divrem renamed to s_mp_divrem_basecase, new s_mp_divrem with Burnikel/Ziegler
 1.11.01  08.03.09  we          improved bz_d3n2n, changed bz_divrem_pos argument list
 1.11.02  14.03.09  we          s_mp_toom3_mul, renamed s_mp_karatsuba_mul/sqr
 1.11.03  15.03.09  we          s_mp_toom3_mul: use s_mp_mod_2k, mp_shl
 1.11.04  15.03.09  we          s_mp_toom3_sqr, changed mp_sqr
 1.11.05  16.03.09  we          s_mp_toom3_mul: split with B ~ 2^(bitsize(max(a,b)/3)
 1.11.06  16.03.09  we          s_mp_fakeinit (used in s_mp_toom3_mul/sqr)
 1.11.07  17.03.09  we          complete rewrite of s_mp_karatsuba_mul/sqr
 1.11.08  20.03.09  we          mp_mul: separate handling of unbalanced factors
 1.11.09  21.03.09  we          mp_mul: fix sign for unbalanced part
 1.11.10  25.03.09  we          fix sign in mp_lshd2
 1.11.11  25.03.09  we          s_mp_toom3_mul uses Bodrato algorithm
 1.11.12  26.03.09  we          s_mp_fakeinit with mp_clamp
 1.11.13  26.03.09  we          s_mp_toom3_sqr uses Bodrato algorithm
 1.11.14  27.03.09  we          mp_clamp without init check, improved fake init in Toom-3
 1.11.15  27.03.09  we          removed MPC_UseToom3 (and mem check for BP7)
 1.11.16  29.03.09  we          improved s_mp_mul_high_digs
 1.11.17  30.03.09  we          mp_init_size2, removed mp_divrem_newton
 1.11.18  30.03.09  we          reduce temporary memory in s_mp_karatsuba_mul/sqr

 1.12.00  20.06.09  we          Fix mp_rand_bits_ex: add mp_clamp(a)

 1.14.00  13.02.10  we          MPC_MAXRadix64 adjustments

 1.16.00  05.06.10  we          mp_read_decimal_astr, mp_read_radix_astr
 1.16.01  25.07.10  we          const pa in mp_read_radix_arr

 1.19.00  03.11.11  we          Memory routines for MAXDigits=32000 with MP_32BIT

 1.20.00  12.01.12  we          Include files mp_bas64/32/16.inc
 1.20.01  13.01.12  we          rewrite (ld/fr)exp(x/d)
 1.20.02  13.01.12  we          adjust mp_toextended_ex for EXT64
 1.20.03  15.01.12  we          Full length for mp_radix_astr
 1.20.04  16.01.12  we          New add32_ovr in mp_bas64.inc
 1.20.05  17.01.12  we          Fix border cases in mp_lshd/2
 1.20.06  19.01.12  we          Subquadratic mp_radix_astr
 1.20.07  20.01.12  we          s_mp_write_radix uses mp_radix_astr for BIT32/64
 1.20.08  21.01.12  we          s_mp_radix_astr

 1.21.00  19.07.12  we          Improved mp_div_int for 0 < b < MP_DIGIT_MAX
 1.21.01  19.07.12  we          s_mp_mod_int, mp_mod_d calls s_mp_mod_int if MP_32BIT
 1.21.02  21.07.12  we          mp_product
 1.21.03  24.07.12  we          isqrt32 without FPU, new is_square32/ex
 1.21.04  25.07.12  we          exptmod32/jacobi32 from mp_numth
 1.21.05  30.07.12  we          Improved mp_set_pow, mp_expt_int: 0^0 = 1

 1.22.00  26.08.12  we          mp_mul: call mp_sqr if @a=@b
 1.22.01  29.08.12  we          kronjac32: used by kronecker32 and jacobi32
 1.22.02  02.09.12  we          mp_radix_astr etc with truncation for BIT16
 1.22.03  04.09.12  we          icbrt32

 1.23.00  23.09.12  we          improved BASM16 s_mp_div_d, toom3 with s_mp_div_d
 1.23.01  25.09.12  we          Fix: uses production code invmod32 in mp_bas16.inc
 1.23.02  01.10.12  we          xgcd32
 1.23.03  02.10.12  we          Avoid -0 in s_mp_divrem_basecase after mp_div_w
 1.23.04  03.10.12  we          mp_set_pow: special code if b is a power of two
 1.23.05  04.10.12  we          BASM16: improved mp_mul_d
 1.23.06  04.10.12  we          mp_div_int with s_mp_divrem_basecase
 1.23.07  05.10.12  we          Removed second s_mp_mod_2k in mp_mod_2k for a<0
 1.23.08  05.10.12  we          Reintroduce s_mp_shrink (mainly for BP7 real mode)
 1.23.09  12.10.12  we          Improved s_mp_n_root2: Use Newton if overflow probable, allow deeper recursion

 1.24.00  12.12.12  we          MPC_PurePascal
 1.24.01  15.12.12  we          Some word/integer types changed to longint or TNInt
 1.24.02  03.01.13  we          Fix some more words to TNInt
 1.24.03  04.01.13  we          Fix mp_set_int for MP_16bit and b=$80000000
 1.24.04  04.01.13  we          Renamed mp_expt_d to mp_expt_w
 1.24.05  05.01.13  we          Fix debug mode icbrt32 for VPC, D2/D3
 1.24.06  05.01.13  we          Fix FPC242/244 bug in debug mode s_mp_div_d

 1.30.00  27.09.14  we          mp_random
 1.30.01  30.09.14  we          Reintroduced mulmod32 function
 1.30.02  03.10.14  we          mlinsolve32

 1.36.00  12.09.17  we          Fix obscure 16-Bit (MSDOS) FP error in s_mp_n_root2

 1.37.00  10.05.18  we          Renamed mp_reverse to s_mp_reverse

 1.38.00  25.06.18  we          frexpx/ldxepx only if not EXT64
 1.38.01  26.06.18  we          s_mp_set_ext ->  s_mp_set_dbl, s_mp_ln/log2 return double

*************************************************************************)


(*-------------------------------------------------------------------------
  This code uses material/ideas from the following 3rd party libraries:
   - LibTomMath 0.30+ by Tom St Denis
   - MPI 1.8.6 by Michael J. Fromberger
  See the file '3rdparty.mpa' for the licenses.
----------------------------------------------------------------------------*)


(*-------------------------------------------------------------------------
 (C) Copyright 2004-2018 Wolfgang Ehrhardt

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


{#Z+}
{----------------------------------}
{32 bit and floating point routines}
{----------------------------------}
{#Z-}

function  add32_ovr(x,y: longint; var z: longint): boolean;
  {-Add z=x+y with overflow detection}

function  bitsize32(a: longint): integer;
  {-Return the number of bits in a (index of highest bit), 0 if no bit is set}

function  gcd32(A, B: longint): longint;
  {-Calculate GCD of two longints}

function  gcd32u(A, B: longint): longint;
  {-Calculate GCD of two longints (DWORD interpretation)}

function  exptmod32(a,b,c: longint): longint;
  {-Calculate a^b mod c if a>=0, b>=0, c>0; result=0 otherwise}

function  icbrt32(a: longint): longint;
  {-Return the integer cube root sign(a)*floor(|a|^(1/3))}

function  invmod32(a,b: longint): longint;
  {-Return a^-1 mod b, b>1. Result is 0 if gcd(a,b)<>1 or b<2}

function  isqrt32(a: longint): longint;
  {-Return floor(sqrt(abs(a))}

function  is_square32(a: longint): boolean;
  {-Test if a is square, i.e. test if a=sqr(isqrt32(a)), false if a<0}

function  is_square32ex(a: longint; var b: longint): boolean;
  {-Test if a is square, false if a<0. If yes, b = sqrt(a) else b is undefined}

function  jacobi32(a,b: longint): integer;
  {-Compute the Jacobi/Legendre symbol (a|b), b: odd and > 2}

function  kronecker32(a,b: longint): integer;
  {-Compute the Kronecker symbol (a|b)}

function  mlinsolve32(a,b,m: longint; var x0,d: longint): boolean;
  {-Solve ax=b mod m; m>1. If result=true then x0 >= 0, d=gcd(a,m) > 0}
  { and there are d solutions: x_k = x0 + k*(m/d) mod m, for k=0..d-1.}

function  mulmod32(a,b,n: longint): longint;
  {-Return a*b mod n, assumes n>0, a,b>=0}

function  popcount16(w: word): integer;
  {-Get population count = number of 1-bits in a word}

function  popcount32(l: longint): integer;
  {-Get population count = number of 1-bits in a longint}

procedure xgcd32(a,b: longint; var u1,u2,u3: longint);
  {-Extended gcd algorithm, calculate a*u1 + b*u2 = u3 = gcd(a,b); a,b <> -2^31}


function  DblPosInf: double;
  {-Return positive double infinity}

function  DblNegInf: double;
  {-Return negative double infinity}

function  DblNaN: double;
  {-Return double NaN (Not a Number)}

procedure frexpd(d: double; var m: double; var e: longint);
  {-Return m,e with d=m*2^e and 0.5 <= abs(m) < 1}

function  ldexpd(d: double; e: longint): double;
  {-Return d*2^e}

{$ifndef EXT64}
procedure frexpx(x: extended; var m: extended; var e: longint);
  {-Return m,e with x=m*2^e and 0.5 <= abs(m) < 1}

function  ldexpx(x: extended; e: longint): extended;
  {-Return x*2^e}
{$endif}


{#Z+}
{---------------------------------------}
{mp functions are sorted alphabetically }
{---------------------------------------}
{#Z-}

procedure mp_2expt(var a: mp_int; b: longint);
  {-Compute a = 2^b, b>=0,  error if b<0 or b>MP_MAXBIT}

procedure mp_abs(const a: mp_int; var b: mp_int);
  {-Absolute value, b = |a|}

procedure mp_add(const a,b: mp_int; var c: mp_int);
  {-High level addition (handles signs)}

procedure mp_add_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Single digit addition}

procedure mp_add_int(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a + b}

function  mp_adecimal(const a: mp_int): ansistring;
  {-Convert to decimal ansistring}

function  mp_ahex(const a: mp_int): ansistring;
  {-Convert to hex ansistring}

function  mp_radix_astr(const a: mp_int; radix: word): ansistring;
  {-Convert to radix representation ansistring (set mp_show_plus to show '+')}

procedure mp_read_decimal_astr(var a: mp_int; const s: ansistring);
  {-Read an mp_int from a decimal ansistring}

procedure mp_read_radix_astr(var a: mp_int; const s: ansistring; radix: word);
  {-Read an mp_int from an ansistring in given radix}

{$ifdef BIT32or64}
function  s_mp_radix_astr(const a: mp_int; radix: word; plus: boolean): ansistring;
  {-Convert to radix representation ansistring, plus to show '+'}
{$endif}

procedure mp_and(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a and b}

function  mp_bitsize(const a: mp_int): longint;
  {-Return the number of bits in a (index of highest bit), 0 if no bit is set}

function  mp_checksum(const a: mp_int): longint;
  {-Return a checksum for a, -1 if mp_error<>MP_OKAY, -2 if not initialized}

procedure mp_chs(const a: mp_int; var b: mp_int);
  {-Change sign, b = -a}

procedure mp_clear(var a: mp_int);
  {-Free an mp_int}

procedure mp_clear2(var a,b: mp_int);
  {-Clear 2 mp_ints}

procedure mp_clear3(var a,b,c: mp_int);
  {-Clear 3 mp_ints}

procedure mp_clear4(var a,b,c,d: mp_int);
  {-Clear 4 mp_ints}

procedure mp_clear5(var a,b,c,d,e: mp_int);
  {-Clear 5 mp_ints}

procedure mp_clear6(var a,b,c,d,e,f: mp_int);
  {-Clear 6 mp_ints}

procedure mp_clear7(var a,b,c,d,e,f,g: mp_int);
  {-Clear 7 mp_ints}

procedure mp_clear8(var a,b,c,d,e,f,g,h: mp_int);
  {-Clear 8 mp_ints}

procedure mp_clear9(var a,b,c,d,e,f,g,h,i: mp_int);
  {-Clear 9 mp_ints}

procedure mp_clear_multi(var vi: array of mp_int);
  {-Clear a vector of mp_ints}

procedure mp_clear_multi_p(const pv: array of pmp_int);
  {-Clear a list of mp_ints given as a ptr vector}

procedure mp_clrbit(var a: mp_int; n: longint);
  {-Clear bit n of a, no action if out of range, (1 = bit 0)}

function  mp_cmp(const a,b: mp_int): integer;
  {-Compare two mp_ints (signed), return sign(a-b)}

function  mp_cmp_d(const a: mp_int; b: mp_digit): integer;
  {-Compare a with an mp_digit, return sign(a-b)}

function  mp_cmp_int(const a: mp_int; b: longint): integer;
  {-Compare a with a longint, return sign(a-b)}

function  mp_cmp_mag(const a,b: mp_int): integer;
  {-Compare magnitude of two mp_ints (unsigned), return sign(|a|-|b|)}

function  mp_cmp_mag_d(const a: mp_int; b: mp_digit): integer;
  {-Compare |a| with a digit, return sign(|a|-b)}

function  mp_cnt_lsb(const a: mp_int): longint;
  {-Count the number of least significant bits which are zero}

procedure mp_copy(const a: mp_int; var b: mp_int);
  {-Copy an mp_int, b = a}

procedure mp_dec(var a: mp_int);
  {-Decrement an mp_int by 1}

procedure mp_dec_int(var a: mp_int; b: longint);
  {-Calculate a = a - b}

function  mp_decimal(const a: mp_int): mp_string;
  {-Convert to decimal, max 255 digits}

procedure mp_div(const a,b: mp_int; var c: mp_int);
  {-Integer signed division, c = a div b}

procedure mp_divrem(const a,b: mp_int; pc,pd: pmp_int);
  {-Integer signed division, pc^ = a div b, pd^ = a rem b; sign(pd^)=sign(a)}

procedure mp_div_2(const a: mp_int; var b: mp_int);
  {-Divide by 2, b = a/2}

procedure mp_div_2k(const a: mp_int; b: longint; var c: mp_int; pd: pmp_int);
  {-Divide by 2^b; quotient in c, optional remainder in pd^, sign(pd^)=sign(a)}

procedure mp_div_d(const a: mp_int; b: mp_digit; pc: pmp_int; var d: mp_digit);
  {-Single digit division, pc^ = a div b, d = a mod b}

procedure mp_div_int(const a: mp_int; b: longint; pc: pmp_int; var d: longint);
  {-Integer signed division, pc^ = a div b, d = a rem b; sign(d)=sign(a)}

procedure mp_div_w(const a: mp_int; b: word; pc: pmp_int; var r: word);
  {-Divide a by a single word b, pc^=sign(a)(|a| div b), r = |a| mod b}

procedure mp_exch(var a,b: mp_int);
  {-Exchange two mp_ints}

procedure mp_expt(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a^b, b>=0}

procedure mp_expt_w(const a: mp_int; b: word; var c: mp_int);
  {-Calculate c = a^b}

procedure mp_expt_int(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a^b, b>=0}

function  mp_gcd_int(const a: mp_int; b: longint): longint;
  {-Return gcd(a,b), b<>0}

function  mp_get_int(const a: mp_int): longint;
  {-Get the lower signed 31 bits of an mp_int}

procedure mp_gr_mod(var x: mp_int; const N, R: mp_int);
  {-Reduce x to x mod N, N > 1, using generalized reciprocal iteration.}
  { Result is >= 0. R is from mp_gr_setup. Compared to the similar}
  { Barrett reduction the restricted range 0<x<N^2 is not required.}

procedure mp_gr_setup(var RN: mp_int; const N: mp_int);
  {-Calculate the generalized reciprocal for N>0}

function  mp_hex(const a: mp_int): mp_string;
  {-Convert to hex string, max 255 digits}

procedure mp_inc(var a: mp_int);
  {-Increment an mp_int by 1}

procedure mp_inc_int(var a: mp_int; b: longint);
  {-Calculate a = a + b}

procedure mp_init(var a: mp_int);
  {-Initialize an mp_int}

procedure mp_init2(var a,b: mp_int);
  {-Initialize 2 mp_ints}

procedure mp_init3(var a,b,c: mp_int);
  {-Initialize 3 mp_ints}

procedure mp_init4(var a,b,c,d: mp_int);
  {-Initialize 4 mp_ints}

procedure mp_init5(var a,b,c,d,e: mp_int);
  {-Initialize 5 mp_ints}

procedure mp_init6(var a,b,c,d,e,f: mp_int);
  {-Initialize 6 mp_ints}

procedure mp_init7(var a,b,c,d,e,f,g: mp_int);
  {-Initialize 7 mp_ints}

procedure mp_init8(var a,b,c,d,e,f,g,h: mp_int);
  {-Initialize 8 mp_ints}

procedure mp_init9(var a,b,c,d,e,f,g,h,i: mp_int);
  {-Initialize 9 mp_ints}

procedure mp_init_copy(var a: mp_int; const b: mp_int);
  {-Create a, then copy b into it}

procedure mp_init_multi(var vi: array of mp_int);
  {-Initialize a vector of mp_ints}

procedure mp_init_multi_p(var pv: array of pmp_int);
  {-Initialize a list of mp_ints given as a ptr vector}

procedure mp_init_set(var a: mp_int; b: mp_digit);
  {-Initialize and set a digit}

procedure mp_init_set_int(var a: mp_int; b: longint);
  {-Initialize and set a to a longint}

procedure mp_init_size(var a: mp_int; size: longint);
  {-Initialize a to size digits, rounded up to multiple of mp_allocprec}

procedure mp_init_size2(var a,b: mp_int; size: longint);
  {-Initialize a and b to size digits, rounded up to multiple of mp_allocprec}

function  mp_isbit(const a: mp_int; n: longint): boolean;
  {-Test if bit n of a is set, (1 = bit 0)}

function  mp_iseven(const a: mp_int): boolean;
  {-Initialized and even}

function  mp_isodd(const a: mp_int): boolean;
  {-Initialized and odd}

function  mp_is_eq(const a,b: mp_int): boolean;
  {-Return a = b}

function  mp_is_ge(const a,b: mp_int): boolean;
  {-Return a >= b}

function  mp_is_gt(const a,b: mp_int): boolean;
  {-Return a > b}

function  mp_is_le(const a,b: mp_int): boolean;
  {-Return a <= b}

function  mp_is_lt(const a,b: mp_int): boolean;
  {-Return a < b}

function  mp_is_ne(const a,b: mp_int): boolean;
  {-Return a <> b}

function  mp_is_longint(const a: mp_int; var b: longint): boolean;
  {-Test if a fits into longint, if true set b := a}

function  mp_is_pow2(const a: mp_int; var n: longint): boolean;
  {-Check if |a| is a power of 2, if true, return n with |a|=2^n}

function  mp_is_pow2_d(d: mp_digit; var n: integer): boolean;
  {-Check if d is power of 2, if true, return n with d=2^n}

function  IsPow2_w(w: word; var n: integer): boolean;
  {-Check if w is power of 2, if true, return n with w=2^n}

function  mp_iszero(const a: mp_int): boolean;
  {-Initialized and zero}

function  mp_is0(const a: mp_int): boolean;
  {-Initialized and = 0}

function  mp_is1(const a: mp_int): boolean;
  {-Initialized and a = 1}

function  mp_is1a(const a: mp_int): boolean;
  {-Initialized and abs(a) = 1}

procedure mp_lshd(var a: mp_int; b: longint);
  {-Shift left a certain amount of digits}

procedure mp_lshd2(const a: mp_int; var b: mp_int; cnt: longint);
  {-Set b to a shifted left by cnt digits}

procedure mp_makeodd(const a: mp_int; var b: mp_int; var s: longint);
  {-Return b,s with a = 2^s*b if a<>0, b=0,s=-1 otherwise}

procedure mp_mod(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a mod b, 0 <= c < b}

procedure mp_mod_2k(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a mod 2^b, 0 <= c < 2^b}

procedure mp_mod_d(const a: mp_int; b: mp_digit; var c: mp_digit);
  {-Calculate c = a mod b, 0 <= c < b (digit version)}

procedure mp_mod_int(const a: mp_int; b: longint; var c: longint);
  {-Calculate c = a mod b}

procedure mp_mod_w(const a: mp_int; b: word; var r: word);
  {-Calculate r = a mod b for a single word b}

procedure mp_montgomery_calcnorm(var R: mp_int; const m: mp_int);
  {-Calculate R = B^n mod m, n=number of digits in m, B=2^DIGIT_BIT}

procedure mp_montgomery_reduce(var x: mp_int; const n: mp_int; rho: mp_digit);
  {-Calculate x = xR^-1 (mod n) via Montgomery reduction}

procedure mp_montgomery_setup(const n: mp_int; var rho: mp_digit);
  {-Calculate rho = -1/n mod B for Montgomery reduction, B=2^DIGIT_BIT}

procedure mp_mul(const a,b: mp_int; var c: mp_int);
  {-High level multiplication, c = a*b}

procedure mp_mul_2(const a: mp_int; var b: mp_int);
  {-Multiply by 2, b = 2*a}

procedure mp_mul_2k(const a: mp_int; b: longint; var c: mp_int);
  {-Shift left a, c = a*2^b; c=a if b<=0 [synonym for mp_shl]}

procedure mp_mul_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Multiply by a digit}

procedure mp_mul_int(const a: mp_int; b: longint; var c: mp_int);
  {-Multiply by a 32 bit integer}

procedure mp_mul_w(const a: mp_int; b: word; var c: mp_int);
  {-Multiply by a word}

function  mp_not_init(const a: mp_int): boolean;
  {-Sanity check if a is initialized, does not catch all cases!}

function  mp_not_init_multi(const a: array of mp_int): boolean;
  {-Sanity check if all elements of a are initialized, does not catch all cases!}

procedure mp_n_root(const a: mp_int; n: longint; var b: mp_int);
  {-Calculate the n'th root of a, a must >=0 if n is even; b=0 if n<1}

procedure mp_n_root2(const a: mp_int; n: longint; var b: mp_int; pr: pmp_int);
  {-Calculate the n'th root of a, pr^=a-b^n, a must be >=0 if n is even; b=0,pr^=0 if n<1}

procedure mp_output_decimal(const a: mp_int);
  {-Write decimal representation to output}

procedure mp_output_radix(const a: mp_int; radix: word);
  {-Write radix representation to output}

procedure mp_or(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a or b}

function  mp_popcount(const a: mp_int): longint;
  {-Get population count = number of 1-bits in a}

procedure mp_product(var a: mp_int; n,m: longint);
  {-Return the product a = n*(n+1)*...(m-1)*m, a=1 if m<n}

procedure mp_prod_int(var a: mp_int; const b: array of longint; n: longint);
  {-Calculate a = product of first n elements of longint array b}

function  mp_radix_size(const a: mp_int; radix: word): longint;
  {-Return size of ASCII representation (incl. sign and #0)}

function  mp_radix_str(const a: mp_int; radix: word): mp_string;
  {-Convert to radix representation, max 255 digits}

procedure mp_rand(var a: mp_int; digits: longint);
  {-Make a pseudo-random mp_int of a given digit size}

procedure mp_random(const a: mp_int; var b: mp_int);
  {-Make a pseudo-random mp_int b with  0 <= b < |a|}

procedure mp_rand_ex(var a: mp_int; digits: longint; sethi: boolean);
  {-Make a pseudo-random mp_int of a given digit size, if not sethi}
  { then a[digits-1] may be zero (and a.used will be decremented)}

procedure mp_rand_bits(var a: mp_int; bits: longint);
  {-Make a pseudo-random mp_int of a given bit size}

procedure mp_rand_bits_ex(var a: mp_int; bits: longint; sethi: boolean);
  {-Make pseudo-random a with bitsize <= bits, if sethi highest bit is set}

procedure mp_rand_radix(var a: mp_int; radix: word; digits: longint);
  {-Make a pseudo-random mp_int of order radix^digits}

procedure mp_read_decimal(var a: mp_int; str: pchar8);
  {-Read an mp_int from a decimal ASCII pchar}

procedure mp_read_decimal_str(var a: mp_int; const s: mp_string);
  {-Read an mp_int from a decimal ASCII string[255]}

procedure mp_read_radix(var a: mp_int; str: pchar8; radix: word);
  {-Read an mp_int from a ASCII pchar in given radix}

procedure mp_read_radix_arr(var a: mp_int; const pa: array of pchar8; radix: word);
  {-Read an mp_int from concatenated pchar array pa in given radix,}
  { max 65000 chars. Mainly used for 16 bit compatibility with max.}
  { length of string literals = 255 and line length = 127 chars    }

procedure mp_read_radix_str(var a: mp_int; const s: mp_string; radix: word);
  {-Read an mp_int from a ASCII string[255] in given radix}

procedure mp_read_signed_bin(var a: mp_int; const b; numbytes: longint);
  {-Read signed bin, big endian, first byte is 0=positive or 1=negative}

procedure mp_read_unsigned_bin(var a: mp_int; const b; numbytes: longint);
  {-Reads a unsigned mp_int, assumes the msb is stored first [big endian]}

procedure mp_reduce(var x: mp_int; const m, mu: mp_int);
  {-Reduce x mod m via Barrett, assumes x<m^2, mu is from mp_reduce_setup}

procedure mp_reduce_setup(var mu: mp_int; const a: mp_int);
  {-Pre-calculate the value required for Barrett reduction}

procedure mp_reduce_2k(var a: mp_int; const n: mp_int; d: mp_digit);
  {-Reduce a mod n where n is of the form 2^p-d, @a<>@n}

procedure mp_reduce_2k_setup(const a: mp_int; var d: mp_digit);
  {-Determine setup value d for unrestricted diminished radix reduction, a>=0}

function  mp_reduce_is_2k(const a: mp_int): boolean;
  {-Determine if mp_reduce_2k can be used}

function  mp_result: integer;
  {-Return and reset mp_error}


procedure mp_rshd(var a: mp_int; b: longint);
  {-Shift right a certain amount of digits}

procedure mp_rshd2(const a: mp_int; var b: mp_int; cnt: longint);
  {-Set b to a shifted right by cnt digits}

procedure mp_set(var a: mp_int; b: mp_digit);
  {-Set a to digit b}

procedure mp_set1(var a: mp_int);
  {-Set a=1}

procedure mp_set_int(var a: mp_int; b: longint);
  {-Set a to a longint}

procedure mp_set_pow(var a: mp_int; b,c: longint);
  {-Set a to b^c, a=0 for c<0}

procedure mp_set_short(var a: mp_int; b: shortint);
  {-Set a to a shortint}

procedure mp_set_w(var a: mp_int; w: word);
  {-Set a to a word}

procedure mp_setbit(var a: mp_int; n: longint);
  {-Set bit n of a, error if n<0 or n>MP_MAXBIT (1 = bit 0)}

procedure mp_shl(const a: mp_int; b: longint; var c: mp_int);
  {-Shift left a, c = a*2^b; c=a if b<=0 }

procedure mp_shlx(const a: mp_int; b: longint; var c: mp_int);
  {-Shift left a b bits if b>=0, shift right |b| if b<0}

procedure mp_shl1(var a: mp_int);
  {-Shift left a by 1}

procedure mp_shr(const a: mp_int; b: longint; var c: mp_int);
  {-Shift right a, c = a/2^b; c=a if b<=0}

procedure mp_shrx(const a: mp_int; b: longint; var c: mp_int);
  {-Shift right a b bits if b>=0, shift left |b| if b<0}

procedure mp_shr1(var a: mp_int);
  {-Divide a by 2, a = a/2}

function  mp_sign(const a: mp_int): integer;
  {-Return sign(a): -1 if a<0, 0 if a=0, +1 if a>0}

function  mp_signed_bin_size(const a: mp_int): longint;
  {-Get the size in bytes for an signed equivalent}

procedure mp_sqr(const a: mp_int; var b: mp_int);
  {-Compute b = a*a}

procedure mp_sqrt(const a: mp_int; var b: mp_int);
  {-Compute b = floor(sqrt(a)), a >=0 using Karatsuba or recursive Newton}

procedure mp_sqrtrem(const n: mp_int; var s,r: mp_int);
  {-Compute Karatsuba square root s and remainder r of n >= 0, n = s^2 + r}

procedure mp_sub(const a,b: mp_int; var c: mp_int);
  {-High level subtraction (handles signs)}

procedure mp_sub_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Single digit subtraction}

procedure mp_sub_int(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a - b}

function  mp_to_signed_bin_n(const a: mp_int; var b; n: longint): longint;
  {-Store in signed big-endian format, max n bytes; return no. of bytes stored}

function  mp_to_unsigned_bin_n(const a: mp_int; var b; n: longint): longint;
  {-Store in unsigned big-endian format, max n bytes; return no. of bytes stored}

function  mp_todouble(const a: mp_int): double;
  {-Convert a to double, +-inf if too large}

function  mp_todouble_ex(const a: mp_int; x: longint): double;
  {-Convert a*2^x to double, +-inf if too large}

{$ifndef EXT64}
function  mp_toextended(const a: mp_int): extended;
  {-Convert a to extended, +-inf if too large}

function  mp_toextended_ex(const a: mp_int; x: longint): extended;
  {-Convert a*2^x to extended, +-inf if too large}
{$endif}

procedure mp_toradix(const a: mp_int; str: pchar8; radix: word);
  {-Store mp_int as a ASCII string in a given radix, better use mp_toradix_n}

procedure mp_toradix_n(const a: mp_int; str: pchar8; radix: word; maxlen: longint);
  {-Store an mp_int as an ASCII string with a given radix (2..MAXRadix)}

function  mp_unsigned_bin_size(const a: mp_int): longint;
  {-Get the size in bytes for an unsigned equivalent}

procedure mp_write_decimal(var tf: system.text; const a: mp_int);
  {-Write decimal representation to file tf}

procedure mp_write_radix(var tf: system.text; const a: mp_int; radix: word);
  {-Write radix representation to file tf}

procedure mp_writeln(const msg: mp_string; const a: mp_int);
  {-Writeln a to stdout with leading msg}

procedure mp_xor(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a xor b}

procedure mp_zero(var a: mp_int);
  {-Set a to zero}

{#Z+}
{---------------------------------------------------------------------------}
{- 'Internal' functions, don't use them unless you know what you are doing -}
{---------------------------------------------------------------------------}
{#Z-}

function  IAlloc(lsize: longint): pointer;
  {-Allocate heap > 64K, return nil if error, no diagnostics}

function  mp_alloc(size: longint): pointer;
  {-Allocate and zero heap, return nil if error}

procedure mp_clamp(var a: mp_int);
  {-Trim unused digits}

procedure mp_freemem(var p: pointer; size: longint);
  {-Deallocate heap if p<>nil, p will be set to nil}

function  mp_getmem(size: longint): pointer;
  {-Allocate heap, return nil if error}

procedure mp_grow(var a: mp_int; size: longint);
  {-Grow an mp_int to a given size (new part is zerofilled)}

procedure mp_init_prim(var a: mp_int; size: longint);
  {-Initialize a to size digits, rounded up to multiple of mp_allocprec}

function  mp_realloc(p: pointer; oldsize, newsize: longint): pointer;
  {-Reallocate heap to new size, if newsize>oldsize the new allocated space is zerofilled}

function  mp_get_allocprec: integer;
  {-Return current value of mp_allocprec}

procedure mp_set_allocprec(prec: integer);
  {-Set new alloc prec 8..64, will be rounded up to power of 2}

procedure s_mp_add(const a,b: mp_int; var c: mp_int);
  {-Low level addition c=a+b, based on HAC pp.594, algorithm 14.7}

procedure s_mp_add_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Single digit addition, no init check, b<>0}

procedure s_mp_checksum(var adler: longint; Msg: pointer; Len: longint);
  {-Update Adler32 checksum with Msg data; init with adler=1}

procedure s_mp_chs(var a: mp_int);
  {-Change sign of an mp_int, no init check}

procedure s_mp_divrem_basecase(const a,b: mp_int; pc,pd: pmp_int);
  {-Integer signed division using Knuth's basecase algorithm D;    }
  { pc^ = a div b, pd^ = a rem b; sign(pd^)=sign(a); no init check.}

procedure s_mp_divrem(const a,b: mp_int; pc,pd: pmp_int);
  {-Integer signed division using recursive Burnikel-Ziegler algorithm:}
  { pc^ = a div b, pd^ = a rem b; sign(pd^)=sign(a); no init check.}
  { Knuth's algorithm D is used for bitsizes < mp_bz_cutoff.}

procedure s_mp_div_d(const a: mp_int; b: mp_digit; pc: pmp_int; var d: mp_digit);
  {-Single digit division, pc^ = a div b, d = a mod b, no init check}

procedure s_mp_expt_dl(a: mp_digit; b: longint; var c: mp_int);
  {-Calculate c = a^b, return 0 for b<0}

procedure s_mp_expt_wl(a: word; b: longint; var c: mp_int);
  {-Calculate c = a^b, return 0 for b<0}

procedure s_mp_fakeinit(const a: mp_int; i0,i1: longint; var b: mp_int);
  {-Make a positive fake mp_int b = a[i0..i1], b=0 if i0>=a.used, or i0>i1.}
  { Internal use only, no init check of a. DANGER: b uses the memory of a, }
  { therefore use b only as and on CONST parameter, DO NOT clear or GROW b!}

function  s_mp_is_le0(const a: mp_int): boolean;
  {-Return true if a<=0, no init check}

procedure s_mp_karatsuba_mul(const a,b: mp_int; var c: mp_int);
  {-Calculate c = |a| * |b| using Karatsuba multiplication}

procedure s_mp_karatsuba_sqr(const a: mp_int; var b: mp_int);
  {-Karatsuba squaring, compute b = a*a using three half size squarings}

function  s_mp_ln(const a: mp_int): double;
  {-Calculate ln(a), a>0. Result=0 for a<=0}

function  s_mp_log2(const a: mp_int): double;
  {-Calculate log2(a), a>0. Result=0 for a<=0}

procedure s_mp_mod_2k(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a mod 2^b, -(|a| mod 2^b) if a < 0}

function  s_mp_mod_is0(const a,b: mp_int): boolean;
  {-Test if a mod b = 0, ie if a is a multiple of b}

procedure s_mp_mod_w(const a: mp_int; b: word; var r: word);
  {-Calculate r = |a| mod b for a single word b}

procedure s_mp_mul_digs(const a,b: mp_int; var c: mp_int; digs: longint);
  {-Multiply |a| * |b| and only compute up to digs digits of result}

procedure s_mp_mul_high_digs(const a,b: mp_int; var c: mp_int; digs: longint);
  {-Multiply |a| * |b| and does not compute the lower digs digits}

procedure s_mp_mul_int(const a: mp_int; b: longint; var c,tmp: mp_int);
  {-Multiply by a 32 bit integer, c=a*b, tmp is an initialized temporary}

function  s_mp_n_root2(const a: mp_int; n: longint; var b: mp_int; pr: pmp_int): boolean;
  {-Calculate the n'th root of a, n>=2, pr^ = a-b^n, a must be >=0 }
  { if n is even; return true, if b is an exact root; no init check}

procedure s_mp_read_radix(var a: mp_int; var str: pchar8; radix: word; const sep: mp_string; SignAllowed: boolean);
  {-Read an ASCII pchar in a given radix into a, breaks on sep and #0.}
  { No init check for a, used for mp_read_radix and mpr_read_radix}

procedure s_mp_reverse(var s; len: longint);
  {-Reverse an array of char, used for radix code}

procedure s_mp_set_dbl(var a: mp_int; x: double; toinf: boolean);
  {-Set a to an double; if toinf, 'round' |a| outward. Error if x=NAN or INF}

procedure s_mp_shrink(var a: mp_int);
  {-Shrink RAM allocated for an mp_int}

procedure s_mp_sqr(const a: mp_int; var b: mp_int);
  {-Low level squaring, b = a*a, HAC pp.596-597, algorithm 14.16}

procedure s_mp_sqrt(const a: mp_int; var b: mp_int);
  {-Compute b = floor(sqrt(a)), a >=0 using recursive integer Newton square root, no init check}

procedure s_mp_sqrtrem(const n: mp_int; var s,r: mp_int);
  {-Compute Karatsuba square root s and remainder r of n >= 0, n = s^2 + r, no init check}

procedure s_mp_sub(const a,b: mp_int; var c: mp_int);
  {-Low level subtraction (assumes |a| > |b|), HAC pp.595 algorithm 14.9}

procedure s_mp_sub_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Single digit subtraction, no init check, b<>0}

procedure s_mp_toom3_mul(const a,b: mp_int; var c: mp_int);
  {-Calculate c = |a| * |b| using Toom-3 multiplication}

procedure s_mp_toom3_sqr(const a: mp_int; var b: mp_int);
  {-Compute b = a*a using Toom-3 squaring}

procedure s_mp_toradix_n(const a: mp_int; radix: word; plus: boolean; var maxlen: longint; var str: pchar8);
  {-Convert an mp_int to ASCII for a given radix, plus=show '+', no init check}

procedure s_mp_write_radix(var tf: system.text; const a: mp_int; radix: word; plus: boolean);
  {-Write radix representation to file tf}

{#Z+}
const
  lograd: array[3..MAXRadix] of double = {ln(2)/ln(radix)}
            (0.630929753571457437, 0.500000000000000000,
             0.430676558073393051, 0.386852807234541587,
             0.356207187108022176, 0.333333333333333333,
             0.315464876785728719, 0.301029995663981195,
             0.289064826317887859, 0.278942945651129843,
             0.270238154427319741, 0.262649535037193548,
             0.255958024809815489, 0.250000000000000000,
             0.244650542118226030, 0.239812466568131445,
             0.235408913366638236, 0.231378213159759174,
             0.227670248696952998, 0.224243824217575439,
             0.221064729457503746, 0.218104291985531559,
             0.215338279036696525, 0.212746053553363154,
             0.210309917857152479, 0.208014597676509458,
             0.205846832460434457, 0.203795047090506190,
             0.201849086582099851, 0.200000000000000000,
             0.198239863170560532, 0.196561632232822608,
             0.194959021893786308, 0.193426403617270793
           {$ifdef MPC_MAXRadix64},
             0.191958720006560152, 0.190551412426773387,
             0.189200359516870037, 0.187901824709107578,
             0.186652411238943373, 0.185449023415368901,
             0.184288833148706181, 0.183169250913633614,
             0.182087900469938248, 0.181042596780040217,
             0.180031326656692640, 0.179052231751041368,
             0.178103593554011088, 0.177183820135557908,
             0.176291434388882096, 0.175425063581954521,
             0.174583430048044930, 0.173765342871439997,
             0.172969690445077101, 0.172195433794098117,
             0.171441600573913449, 0.170707279663720130,
             0.169991616286914032, 0.169293807598781433,
             0.168613098689501110, 0.167948778957041942,
             0.167300178810174135, 0.166666666666666667);
           {$else}
             );
           {$endif}
{#Z-}


implementation


uses
  mp_prng;  {functions for PRNG generation}

{$ifdef BIT16}
  type TArrayOfByte = packed array[0..$F000] of byte; {helper type}
{$else}
  type TArrayOfByte = packed array[0..MaxLongint-1] of byte; {helper type}
{$endif}

var
  mp_allocprec: integer;
  mp_allocmask: longint;

const
  PopCnt: array[byte] of byte = {number of 1-bits in a byte}
           (0,1,1,2,1,2,2,3,1,2,2,3,2,3,3,4,1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,
            1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
            1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
            2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
            1,2,2,3,2,3,3,4,2,3,3,4,3,4,4,5,2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,
            2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
            2,3,3,4,3,4,4,5,3,4,4,5,4,5,5,6,3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,
            3,4,4,5,4,5,5,6,4,5,5,6,5,6,6,7,4,5,5,6,5,6,6,7,5,6,6,7,6,7,7,8);

const
  _DPInfA: TMPHexDblW = (0,0,0,$7FF0); {Bit pattern positive double infinity}
  _DNInfA: TMPHexDblW = (0,0,0,$FFF0); {Bit pattern negative double infinity}
  _DNaNA : TMPHexDblW = (0,0,0,$FFFF); {Bit pattern double NaN(Not a Number)}


{---------------------------------------------------------------------------}
function  DblPosInf: double;
  {-Return positive double infinity}
begin
  DblPosInf := double(_DPInfA);
end;


{---------------------------------------------------------------------------}
function  DblNegInf: double;
 {-Return negative double infinity}
begin
  DblNegInf:= double(_DNInfA);
end;


{---------------------------------------------------------------------------}
function  DblNaN: double;
 {-Return double NaN (Not a Number)}
begin
  DblNaN:= double(_DNanA);
end;


{---------------------------------------------------------------------------}
procedure frexpd(d: double; var m: double; var e: longint);
  {-Return the mantissa m and exponent e of d with d = m*2^e, 0.5 <= abs(m) < 1;}
  { if d is 0, +-INF, NaN or denormal, return m=d, e=0}
begin
  e := TMPHexDblW(d)[3] and $7FF0;
  {First check is INF or NAN, then if d is zero/denormal}
  if (e=$7FF0) or (e=0) or (d=0.0) then e := 0
  else begin
    e := (e shr 4) - $3FE;
    TMPHexDblW(d)[3] := (TMPHexDblW(d)[3] and $800F) or $3FE0;
  end;
  m := d;
end;


{---------------------------------------------------------------------------}
function ldexpd(d: double; e: longint): double;
  {-Return d*2^e}
var
  i: integer;
const
  H2_54: TMPHexDblW = ($0000,$0000,$0000,$4350);  {2^54}
begin
  {if +-INF, NaN, 0 or if e=0 return d}
  i := (TMPHexDblW(d)[3] and $7FF0) shr 4;
  if (i=$7FF) or (e=0) or (d=0.0) then ldexpd := d
  else if i=0 then begin
    {Denormal: result = d*2^54*2^(e-54)}
    ldexpd := ldexpd(d*double(H2_54), e-54);
  end
  else begin
    e := e+i;
    if e>$7FE then begin
      {overflow}
      if d>0.0 then ldexpd := double(_DPInfA)
      else ldexpd := double(_DNInfA);;
    end
    else if e<1 then begin
      {underflow or denormal}
      if e<-53 then ldexpd := 0.0
      else begin
        {Denormal: result = d*2^(e+54)/2^54}
        inc(e,54);
        TMPHexDblW(d)[3] := (TMPHexDblW(d)[3] and $800F) or (e shl 4 and $7FF0);
        ldexpd := d/double(H2_54);
      end;
    end
    else begin
      TMPHexDblW(d)[3] := (TMPHexDblW(d)[3] and $800F) or (e shl 4 and $7FF0);
      ldexpd := d;
    end;
  end;
end;


{$ifndef EXT64}
{---------------------------------------------------------------------------}
procedure frexpx(x: extended; var m: extended; var e: longint);
  {-Return the mantissa m and exponent e of x with x = m*2^e, 0.5 <= abs(m) < 1;}
  { if x is 0, +-INF, NaN or denormal, return m=x, e=0}
begin
  e  := TMPHexExtW(x)[4] and $7FFF;
  {First check is INF or NAN, then if x is zero/denormal}
  if (e=$7FFF) or (e=0) or (x=0.0) then e := 0
  else begin
    dec(e,$3FFE);
    TMPHexExtW(x)[4] := (TMPHexExtW(x)[4] and $8000) or $3FFE;
  end;
  m := x;
end;


{---------------------------------------------------------------------------}
function ldexpx(x: extended; e: longint): extended;
  {-Return x*2^e}
var
  i: integer;
const
  H2_64: TMPHexExtW = ($0000,$0000,$0000,$8000,$403f);  {2^64}
begin
  {if +-INF, NaN, 0 or if e=0 return x}
  i := TMPHexExtW(x)[4] and $7FFF;
  if (i=$7FFF) or (e=0) or (x=0.0) then ldexpx := x
  else if i=0 then begin
    {Denormal: result = x*2^64*2^(e-64)}
    ldexpx := ldexpx(x*extended(H2_64), e-64);
  end
  else begin
    e := e+i;
    if e>$7FFE then begin
      {overflow}
      if x>0.0 then ldexpx := double(_DPInfA)
      else ldexpx := double(_DNInfA);
    end
    else if e<1 then begin
      {underflow or denormal}
      if e<-63 then ldexpx := 0.0
      else begin
        {Denormal: result = x*2^(e+64)/2^64}
        inc(e,64);
        TMPHexExtW(x)[4] := (TMPHexExtW(x)[4] and $8000) or (e and $7FFF);
        ldexpx := x/extended(H2_64);
      end;
    end
    else begin
      TMPHexExtW(x)[4] := (TMPHexExtW(x)[4] and $8000) or (e and $7FFF);
      ldexpx := x;
    end;
  end;
end;
{$endif}


{---------------------------------------------------------------------------}
{------------ begin basic compiler dependent core routines -----------------}
{------------ should be grouped together on top of mp_base -----------------}
{---------------------------------------------------------------------------}

{$ifdef MPC_PurePascal}
  {$i mp_bas64.inc}
{$else}
  {$ifdef BIT16}
    {$i mp_bas16.inc}
  {$else}
    {$i mp_bas32.inc}
  {$endif}
{$endif}

{$ifdef MPC_FPU_ISQRT}
{---------------------------------------------------------------------------}
function isqrt32(a: longint): longint;
  {-Return floor(sqrt(abs(a))}
begin
  isqrt32 := trunc(sqrt(abs(a)));
end;
{$else}
{---------------------------------------------------------------------------}
function isqrt32(a: longint): longint;
  {-Return floor(sqrt(abs(a))}
var
  n,r,b,t: longint;
begin
  {Based on Wilco Dijkstra's C code from comp.sys.arm}
  {http://groups.google.com/group/comp.sys.arm/msg/2f9fe6ab123f6a6d?dmode=source&output=gplain}
  {WE changes: local n, abs(arg), partial unroll}
  r := 0;
  n := abs(a);
  if n >= $10000 then b := $40000000
  else b := $4000;
  repeat
    t := r+b;
    if n >= t then begin
      dec(n, t);
      r := t + b;
    end;
    r := r shr 1;
    b := b shr 2;
    t := r+b;
    if n >= t then begin
      dec(n, t);
      r := t + b;
    end;
    r := r shr 1;
    b := b shr 2;
    t := r+b;
    if n >= t then begin
      dec(n, t);
      r := t + b;
    end;
    r := r shr 1;
    b := b shr 2;
    t := r+b;
    if n >= t then begin
      dec(n, t);
      r := t + b;
    end;
    r := r shr 1;
    b := b shr 2;
    t := r+b;
    if n >= t then begin
      dec(n, t);
      r := t + b;
    end;
    r := r shr 1;
    b := b shr 2;
    t := r+b;
    if n >= t then begin
      dec(n, t);
      r := t + b;
    end;
    r := r shr 1;
    b := b shr 2;
    t := r+b;
    if n >= t then begin
      dec(n, t);
      r := t + b;
    end;
    r := r shr 1;
    b := b shr 2;
    t := r+b;
    if n >= t then begin
      dec(n, t);
      r := t + b;
    end;
    r := r shr 1;
    b := b shr 2;
  until b=0;
  isqrt32 := r;
end;
{$endif}


const {bit i is cleared if i mod 128 is a square}
  ba_128: array[0..15] of byte = ($ec,$fd,$fc,$fd,$ed,$fd,$fd,$fd,$ec,$fd,$fd,$fd,$ed,$fd,$fd,$fd);

const {set of byte size squares}
  sqrset = [0,1,4,9,16,25,36,49,64,81,100,121,144,169,196,225];

{---------------------------------------------------------------------------}
function is_square32(a: longint): boolean;
  {-Test if a is square, i.e. test if a=sqr(isqrt32(a)), false if a<0}
var
  i: integer;
begin
  is_square32 := false;
  if a<256 then begin
    if a<0 then exit;
    is_square32 := a in sqrset;
  end
  else begin
    {First use the mp_is_square2 mod 128 test with a 82.03% rejection rate}
    i := a and 127;
    if ba_128[i shr 3] and (1 shl (i and 7)) <> 0 then exit;
    is_square32 := sqr(isqrt32(a))=a;
  end;
end;


{---------------------------------------------------------------------------}
function is_square32ex(a: longint; var b: longint): boolean;
  {-Test if a is square, false if a<0. If yes, b = sqrt(a) else b is undefined}
var
  i: integer;
begin
  is_square32ex := false;
  if a<256 then begin
    if a<0 then exit;
    if a in sqrset then begin
      is_square32ex := true;
      b := isqrt32(a);
    end;
  end
  else begin
    {First use the mp_is_square2 mod 128 test with a 82.03% rejection rate}
    i := a and 127;
    if ba_128[i shr 3] and (1 shl (i and 7)) <> 0 then exit;
    b := isqrt32(a);
    is_square32ex := sqr(b)=a;
  end;
end;


{$ifdef MPC_PurePascal}
{---------------------------------------------------------------------------}
function exptmod32(a,b,c: longint): longint;
  {-Calculate a^b mod c if a>=0, b>=0, c>0; result=0 otherwise}
var
  p: longint;
begin
  exptmod32 := 0;
  if (c<2) or (a<=0) or (b<0) then exit;
  if a>=c then begin
    a := a mod c;
    if a=0 then exit;
  end;
  if b<2 then begin
    if b=0 then exptmod32 := 1 else exptmod32 := a;
    exit;
  end;
  p := 1;
  while b<>0 do begin
    if odd(b) then p := int64(p)*a mod c;
    a := int64(a)*a mod c;
    b := b shr 1;
  end;
  exptmod32 := p;
end;
{$else}

{---------------------------------------------------------------------------}
function exptmod32(a,b,c: longint): longint;
  {-Calculate a^b mod c if a>=0, b>=0, c>0; result=0 otherwise}
begin
  exptmod32 := 0;
  if (c<2) or (a<=0) or (b<0) then exit;
  if a>=c then begin
    a := a mod c; if a=0 then exit;
  end;
  if b<2 then begin
    if b=0 then exptmod32 := 1 else exptmod32 := a;
    exit;
  end;
  {$ifdef BIT32}
    asm
           push  edi
           push  esi
           push  ebx
           mov   edi,$40000000
           mov   ebx,[b]
      @@1: test  edi,ebx
           jnz   @@2
           shr   edi,1
           jmp   @@1
      @@2: shr   edi,1
           mov   esi,[a]
           mov   ecx,[c]
           mov   edx,esi
      @@3: mov   eax,edx
           mul   edx
           div   ecx
           test  edi,ebx
           jz    @@4
           mov   eax,esi
           mul   edx
           div   ecx
      @@4: shr   edi,1
           jnz   @@3
         {$ifdef FPC}
           mov   [a],edx
         {$else}
           mov   [@result],edx
         {$endif}
           pop   ebx
           pop   esi
           pop   edi
    end;
  {$else}
    asm
                   mov  di,$4000
           db $66; shl  di,16
           db $66; mov  bx,word ptr [b]
      @@1: db $66; test di,bx
           db $66; jnz  @@2
           db $66; shr  di,1
                   jmp  @@1
      @@2: db $66; shr  di,1
           db $66; mov  si,word ptr [a]
           db $66; mov  cx,word ptr [c]
           db $66; mov  dx,si
      @@3: db $66; mov  ax,dx
           db $66; mul  dx
           db $66; div  cx
           db $66; test di,bx
                   jz   @@4
           db $66; mov  ax,si
           db $66; mul  dx
           db $66; div  cx
      @@4: db $66; shr  di,1
                   jnz  @@3
           db $66; mov  word ptr [@result],dx
    end;
  {$endif}
  {$ifdef FPC}
    exptmod32 := a;
  {$endif}
end;

{$endif}


{$ifndef BIT16}

{$ifdef HAS_CARD32}
  {---------------------------------------------------------------------------}
  function icbrt32(a: longint): longint;
    {-Return the integer cube root sign(a)*floor(|a|^(1/3))}
  var
    s: integer;
    x,y,y2,b: cardinal; {must be declared unsigned for test x >= b}
  begin
    {Ref: http://www.hackersdelight.org/hdcodetxt/icbrt.c.txt}
    x  := abs(a);
    y  := 0;
    s  := 30;
    y2 := 0;
    while s>=0 do begin
      y2 := y2 shl 2;
      y  := y+y;
      b  := (3*(y2 + y) + 1) shl s;
      if x >= b then begin
        x  := x - b;
        y2 := y2 + y + y + 1;
        y  := y + 1;
      end;
      dec(s,3);
    end;
    if a<0 then icbrt32 := -y
    else icbrt32 := y;
  end;
{$else}

  {For VPC, D2/D3; cardinal is only 31 bit}
  {---------------------------------------------------------------------------}
  function icbrt32(a: longint): longint;
    {-Return the integer cube root sign(a)*floor(|a|^(1/3))}
  var
    x: longint;
  begin
    {Ref: http://www.hackersdelight.org/hdcodetxt/icbrt.c.txt}

    {Note: There will be RTEs for large a if Range!? checking (not obverflow)}
    {is done, e.g. in debug mode. Here a 32-bit BASM port of the Pascal code:}
    x := abs(a);
    if x>7 then begin
      {y=ebx, y2=esi, s=cx, b=edi, x=eax}
      asm
            push edi
            push esi
            push ebx
            mov  eax,[x]
            xor  ebx,ebx       {y  := 0}
            xor  esi,esi       {y2 := 0}
            mov  cx,30         {s  := 30}
      @@1:  shl  esi,2         {y2 := y2 shl 2}
            add  ebx,ebx
            mov  edx,ebx
            add  edx,esi       {edx:= y2+y}
            mov  edi,edx
            add  edi,edi
            add  edi,edx       {edi:= 3*(y2+y)}
            inc  edi
            shl  edi,cl        {b  := 3*(y2 + y) + 1) shl s}
            cmp  edi,eax
            ja   @@2
            sub  eax,edi       {x  := x - b}
            add  esi,ebx
            inc  ebx           {y  := y + 1}
            add  esi,ebx       {y2 := y2 + y + (y+1)}
      @@2:  sub  cx,3
            jnl  @@1
            mov  [x],ebx
            pop  ebx
            pop  esi
            pop  edi
      end;
      if a<0 then icbrt32 := -x
      else icbrt32 := x;
    end
    else begin
      { -7 <= a <= 7, ie icbrt = -1,0,1}
      if a>0 then icbrt32 := 1
      else if a<0 then icbrt32 := -1
      else icbrt32 := 0;
    end;
  end;
{$endif}

{$else}
{---------------------------------------------------------------------------}
function icbrt32(a: longint): longint;
  {-Return the integer cube root sign(a)*floor(|a|^(1/3))}
var
  x: longint;
begin
  {Ref: http://www.hackersdelight.org/hdcodetxt/icbrt.c.txt}

  {Note: Floating point implementation icbrt = trunc(exp(ln(abs(a)/3)) }
  {does not work without various small corrections for different ranges}
  {of a. The fixed code is not faster even for 16 bit, therefore the 32}
  {bit Pascal code is translated to 16 bit BASM.}
  x := abs(a);
  if x>7 then begin
    {y=ebx, y2=esi, s=cx, b=edi, x=eax}
    asm
         db $66;   mov  ax,word[x]
         db $66;   xor  bx,bx        {y  := 0}
         db $66;   xor  si,si        {y2 := 0}
                   mov  cx,30        {s  := 30}
    @@1: db $66;   shl  si,2         {y2 := y2 shl 2}
         db $66;   add  bx,bx
         db $66;   mov  dx,bx
         db $66;   add  dx,si        {edx:= y2+y}
         db $66;   mov  di,dx
         db $66;   add  di,di
         db $66;   add  di,dx        {edi:= 3*(y2+y)}
         db $66;   inc  di
         db $66;   shl  di,cl        {b  := 3*(y2 + y) + 1) shl s}
         db $66;   cmp  di,ax
                   ja   @@2
         db $66;   sub  ax,di        {x  := x - b}
         db $66;   add  si,bx
         db $66;   inc  bx           {y  := y + 1}
         db $66;   add  si,bx        {y2 := y2 + y + (y+1)}
    @@2:           sub  cx,3
                   jnl  @@1
         db $66;   mov  word[x],bx
    end;
    if a<0 then icbrt32 := -x
    else icbrt32 := x;
  end
  else begin
    { -7 <= a <= 7, ie icbrt = -1,0,1}
    if a>0 then icbrt32 := 1
    else if a<0 then icbrt32 := -1
    else icbrt32 := 0;
  end;
end;
{$endif}

{------------------------------------------------------------------------}
{------------------------------------------------------------------------}
{-----------  end of basic compiler dependent core routines   -----------}
{------------------------------------------------------------------------}
{------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function kronjac32(a,b: longint; kron: boolean): integer;
  {-Compute the Kronecker/Jacobi/Legendre symbol (a|b)}
var
  res,m8: integer;
  t: longint;
begin
  {Uses classic Jacobi method using only mod see e.g. [10] Alg. 2.3.5.}
  kronjac32 := 0;

  {initialize accumulated result}
  res := 1;

  if kron then begin
    {Compute Kronecker symbol, reduce to Jacobi with}
    if b=0 then begin
      if abs(a)=1 then kronjac32 := 1;
      exit;
    end;

    {here b<>0, make b positive}
    if b<0 then begin
      {(a|b) = (a|-1)*(a|-b)}
      if a<0 then res := -res;
      b := abs(b);
    end;

    {if b even, reduce to odd case}
    if b and 1 = 0 then begin
      {(a|2)=0 if a is even}
      if a and 1 = 0 then exit;
      {divide out powers of 4}
      while b and 3 = 0 do b := b shr 2;
      {if b is even, divide by 2 and adjust result}
      if b and 1 = 0 then begin
        b  := b shr 1;
        m8 := a and 7;
        if (m8=3) or (m8=5) then res := -res;
      end;
    end;
  end
  else begin
    {Compute Jacobi/Legendre symbol: Check range of b}
    if (b<3) or (b and 1 = 0) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('kronjac32: b<3 or even');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  end;

  {Here b is positive and odd, and we actually calculate a Jacobi symbol.}
  {if a<0 use property of Jacobi symbol to make it positive. Don't rely  }
  {on Pascal mod function to return positive result.}
  if a<0 then begin
    {(-a|b) = (a|b)*(-1|b); change sign if b=3 (mod 4)}
    if b and 3 = 3 then res := -res;
    a := -a;
  end;

  {do initial reduction to force a < b}
  if a>=b then a := a mod b;

  {loop invariant: a < b}
  while a<>0 do begin
    {divide out powers of 4}
    while a and 3 = 0 do a := a shr 2;
    {if a is even, divide by 2 and adjust result}
    if a and 1 = 0 then begin
      a  := a shr 1;
      m8 := b and 7;
      { (2|b) = -1  if b = +-3 (mod 8)}
      if (m8=3) or (m8=5) then res := -res;
    end;
    if a=1 then begin
      kronjac32 := res;
      exit;
    end;
    {swap variables, reduce and adjust result using quadratic reciprocity}
    if (a and b and 3) = 3 then res := -res;
    t := b;
    b := a;
    a := t mod a;
  end;
  if b=1 then kronjac32 := res;
end;


{---------------------------------------------------------------------------}
function jacobi32(a,b: longint): integer;
  {-Compute the Jacobi/Legendre symbol (a|b), b: odd and > 2}
begin
  jacobi32 := kronjac32(a,b,false);
end;


{---------------------------------------------------------------------------}
function kronecker32(a,b: longint): integer;
  {-Compute the Kronecker symbol (a|b)}
begin
  kronecker32 := kronjac32(a,b,true);
end;


{---------------------------------------------------------------------------}
procedure xgcd32(a,b: longint; var u1,u2,u3: longint);
  {-Extended gcd algorithm, calculate a*u1 + b*u2 = u3 = gcd(a,b); a,b <> -2^31}
var
  v1,v2,v3,t1,t2,t3,q: longint;
  alt0, blt0, bgta: boolean;
begin

  {Algorithm and notation is from Knuth [3] Alg.X. Normalisation}
  {is used to avoid overflows, see e.g. Shoup [29], Theorem 4.3.}

  {Note that the optimization suggested in Knuth by suppressing }
  {v2 and t2 and using u2 = (u3 - a*u1) div b will overflow even}
  {for relative small values like a = 765432 and b = 654321!    }

  {Normalise a >= b >= 0}
  alt0 := a<0;  if alt0 then a := -a;
  blt0 := b<0;  if blt0 then b := -b;
  bgta := b>a;
  if bgta then begin
    q := a;
    a := b;
    b := q;
  end;

  {Initialize}
  u1 := 1;
  u2 := 0;
  u3 := a;
  v1 := 0;
  v2 := 1;
  v3 := b;

  {Euclidian loop}
  while v3<>0 do begin
    q  := u3 div v3;
    t1 := u1 - q*v1;
    t2 := u2 - q*v2;
    t3 := u3 - q*v3;
    u1 := v1;
    u2 := v2;
    u3 := v3;
    v1 := t1;
    v2 := t2;
    v3 := t3;
  end;

  {Undo normalisation}
  if bgta then begin
    q  := u1;
    u1 := u2;
    u2 := q;
  end;
  if alt0 then u1 := -u1;
  if blt0 then u2 := -u2;
end;


{---------------------------------------------------------------------------}
function mlinsolve32(a,b,m: longint; var x0,d: longint): boolean;
  {-Solve ax=b mod m; m>1. If result=true then x0 >= 0, d=gcd(a,m) > 0}
  { and there are d solutions: x_k = x0 + k*(m/d) mod m, for k=0..d-1.}
var
  y: longint;
begin
  mlinsolve32 := false;
  if m>0 then begin
    {normalize a,b because mod may return negative results}
    if (a<0) or (a >= m) then begin
      a := a mod m;
      if a<0 then a := a + m;
    end;
    if (b<0) or (b >= m) then begin
      b := b mod m;
      if b<0 then b := b + m;
    end;
    xgcd32(a,m,x0,y,d);
    if b mod d = 0 then begin
      mlinsolve32 := true;
      while x0<0 do x0 := x0+m;
      x0 := mulmod32(x0,b div d,m);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_set_allocprec(prec: integer);
  {-Set new alloc prec 8..64, will be rounded up to power of 2}
begin
  if prec<=8 then mp_allocprec := 8
  else if prec<=16 then mp_allocprec := 16
  else if prec<=32 then mp_allocprec := 32
  else mp_allocprec := 64;
  mp_allocmask := pred(mp_allocprec);
end;


{---------------------------------------------------------------------------}
function mp_get_allocprec: integer;
  {-Return current value of mp_allocprec}
begin
  mp_get_allocprec := mp_allocprec;
end;


{---------------------------------------------------------------------------}
procedure mp_freemem(var p: pointer; size: longint);
  {-Deallocate heap if p<>nil, p will be set to nil}
begin
  if p<>nil then begin
    freemem(p, size);
    p := nil;
    {$ifdef MPC_Diagnostic}
      dec(mp_memstat.MemDiff, size);
    {$endif}
  end;
end;


{---------------------------------------------------------------------------}
function mp_getmem(size: longint): pointer;
  {-Allocate heap, return nil if error}
var
  p: pointer;
begin
  p := IAlloc(size);
  {$ifdef MPC_Diagnostic}
    if p<>nil then begin
      inc(mp_memstat.MemDiff, size);
      inc(mp_memstat.ACntLow);
      if mp_memstat.ACntLow=0 then inc(mp_memstat.ACntHigh);
    end;
  {$endif}
  mp_getmem := p;
end;


{--------------------------------------------------------------------------}
function mp_alloc(size: longint): pointer;
  {-Allocate and zero heap, return nil if error}
var
  p: pointer;
begin
  p := mp_getmem(size);
  mp_alloc := p;
  if p<>nil then fillchar(p^, size, 0);
end;


{---------------------------------------------------------------------------}
procedure mp_init_prim(var a: mp_int; size: longint);
  {-Initialize a to size digits, rounded up to multiple of mp_allocprec}
begin
  if mp_error<>MP_OKAY then exit;
  if size=0 then size := mp_allocprec
  else size := (size+mp_allocmask) and (not mp_allocmask);
  with a do begin
    if size>MAXDigits then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXMaxDigits.Create('mp_init_prim');
        {$else}
          RunError(MP_RTE_OTHER);
        {$endif}
      {$else}
        set_mp_error(MP_MAXDIGITS);
        fillchar(a, sizeof(a), 0);
        exit;
      {$endif}
    end;
    {allocate memory required and clear it}
    pdigits := mp_alloc(sizeof(mp_digit)*size);
    if pdigits=nil then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXMemory.Create('mp_init_prim: alloc');
        {$else}
          RunError(MP_RTE_MEM);
        {$endif}
      {$else}
        set_mp_error(MP_MEM);
        fillchar(a, sizeof(a), 0);
        exit;
      {$endif}
    end;
    {set the used to zero, allocated digits to the}
    {default precision and sign to positive       }
    used  := 0;
    alloc := size;
    sign  := MP_ZPOS;
    magic := MP_MAGIC;
    {$ifdef MPC_Diagnostic}
      inc(mp_memstat.InitDiff);
    {$endif}
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_clamp(var a: mp_int);
  {-Trim unused digits}
begin
  if mp_error<>MP_OKAY then exit;
  {Decrease used while the most significant digit is zero. This is    }
  {used to ensure that leading zero digits are trimmed and the leading}
  {"used" digit will be non-zero. Typically very fast. Also fixes the }
  {sign if there are no more leading digits.                          }

  {No arg check since mp_clamp will be called after working with a    }
  with a do begin
    if pdigits<>nil then begin
      while (used>0) and (pdigits^[pred(used)]=0) do dec(used);
    end;
    {reset the sign flag if used=0}
    if used=0 then sign := MP_ZPOS;
  end;
end;


{---------------------------------------------------------------------------}
function IsPow2_w(w: word; var n: integer): boolean;
  {-Check if w is power of 2, if true, return n with w=2^n}
var
  i: integer;
begin
  IsPow2_w := false;
  {First do a quick check to get the result}
  if (w=0) or (w and pred(w) <> 0) then exit;
  {now get bit index in highest digit}
  for i:= 0 to 15 do begin
    if w=(word(1) shl i) then begin
      IsPow2_w := true;
      n := i;
      exit;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function popcount16(w: word): integer;
  {-Get population count = number of 1-bits in a word}
var
  ba: packed array[0..1] of byte absolute w;
begin
  popcount16 := PopCnt[ba[0]]+PopCnt[ba[1]];
end;


{---------------------------------------------------------------------------}
function popcount32(l: longint): integer;
  {-Get population count = number of 1-bits in a longint}
var
  ba: packed array[0..3] of byte absolute l;
begin
  popcount32 := PopCnt[ba[0]]+PopCnt[ba[1]]+PopCnt[ba[2]]+PopCnt[ba[3]];
end;


{---------------------------------------------------------------------------}
procedure mp_2expt(var a: mp_int; b: longint);
  {-Compute a = 2^b, b>=0,  error if b<0 or b>MP_MAXBIT}
begin
  {all checking is done in the used routines}
  {clear all bits}
  mp_zero(a);
  {set the specified bit}
  mp_setbit(a,b);
end;


{---------------------------------------------------------------------------}
procedure mp_abs(const a: mp_int; var b: mp_int);
  {-Absolute value, b = |a|}
begin
  {Arg check in mp_copy}
  mp_copy(a, b);
  if mp_error=MP_OKAY then begin
    {force the sign of b to positive}
    b.sign := MP_ZPOS;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_add(const a,b: mp_int; var c: mp_int);
  {-High level addition (handles signs)}
var
  cmp: integer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_add');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {easy outs} {*0.1.29}
  if a.used=0 then begin
    mp_copy(b,c);
    exit;
  end;
  if b.used=0 then begin
    mp_copy(a,c);
    exit;
  end;

  {handle two cases, not four}
  if a.sign=b.sign then begin
    {both positive or both negative}
    {add their magnitudes, copy the sign}
    c.sign := a.sign;
    s_mp_add(a, b, c);
  end
  else begin

    {one positive, the other negative                 }
    {subtract the one with the greater magnitude from }
    {the one of the lesser magnitude.  The result gets}
    {the sign of the one with the greater magnitude.  }

    {Easy out, and make inputs to s_mp_sub unequal    } {*0.1.29}
    cmp := mp_cmp_mag(a, b);
    if cmp=MP_EQ then mp_zero(c)
    else if cmp=MP_LT then begin
      c.sign := b.sign;
      s_mp_sub(b, a, c);
    end
    else begin
      c.sign := a.sign;
      s_mp_sub(a, b, c);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_add_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Single digit addition}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_add_d');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if b>MP_DIGIT_MAX then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_add_d: b>MP_DIGIT_MAX');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}
  {V0.6.06}
  if b=0 then mp_copy(a,c) else s_mp_add_d(a,b,c);
end;


{---------------------------------------------------------------------------}
procedure mp_add_int(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a + b}
var
  t: mp_int;
  d: mp_digit;
begin
  if abs(b)<=lv_digit_max then begin
    {$ifdef MPC_ArgCheck}
      if mp_not_init(a) or mp_not_init(c) then begin
        {$ifdef MPC_HaltOnArgCheck}
          {$ifdef MPC_UseExceptions}
            raise MPXNotInit.Create('mp_add_int');
          {$else}
            RunError(MP_RTE_NOTINIT);
          {$endif}
        {$else}
          set_mp_error(MP_NOTINIT);
          exit;
        {$endif}
      end;
    {$endif}
    d := abs(b);
    if b>0 then s_mp_add_d(a,d,c)
    else if b<0 then s_mp_sub_d(a,d,c)
    else mp_copy(a,c);
  end
  else begin
    mp_init_set_int(t,b);
    {Debug checking in called functions}
    if mp_error=MP_OKAY then begin
      mp_add(a,t,c);
      mp_clear(t);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function mp_adecimal(const a: mp_int): ansistring;
  {-Convert to decimal ansistring}
begin
  mp_adecimal := mp_radix_astr(a, 10);
end;


{---------------------------------------------------------------------------}
function mp_ahex(const a: mp_int): ansistring;
  {-Convert to hex ansistring}
begin
  mp_ahex := mp_radix_astr(a, 16);
end;


{---------------------------------------------------------------------------}
procedure mp_and(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a and b}
var
  pa,pb,pc: pmp_digit;
  olduse, min: longint;
  i: TNInt;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_and');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if a.used < b.used then min := a.used else min := b.used;
  if min=0 then begin
    mp_zero(c);
    exit;
  end;

  if c.alloc < min then begin
    mp_grow(c, min);
    if mp_error<>MP_OKAY then exit;
  end;

  pa  := pmp_digit(a.pdigits);
  pb  := pmp_digit(b.pdigits);
  pc  := pmp_digit(c.pdigits);

  {get old used digit count and set new one}
  olduse := c.used;
  c.used := min;
  c.sign := a.sign and b.sign;

  {process digits 0 .. min-1}
  for i:=0 to pred(min) do begin
    pc^ := pa^ and pb^;
    inc(pa);
    inc(pb);
    inc(pc);
  end;

  {clear c above min if necessary}
  if olduse>min then begin
    for i:=min to pred(olduse) do begin
      pc^ := 0;
      inc(pc);
    end;
  end;
  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
function mp_bitsize(const a: mp_int): longint;
  {-Return the number of bits in a (index of highest bit), 0 if no bit is set}
begin
  mp_bitsize := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_bitsize');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {bits used = (used-1)*DIGIT_BIT + bits in highest digit}
  if a.used>0 then begin
    mp_bitsize := longint(a.used-1)*DIGIT_BIT + bitsize32(a.pdigits^[a.used-1]);
  end;
end;


{---------------------------------------------------------------------------}
function mp_checksum(const a: mp_int): longint;
  {-Return a checksum for a, -1 if mp_error<>MP_OKAY, -2 if not initialized}
var
  adler: longint;
begin
  mp_checksum := -1;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      mp_checksum := -2;
      exit;
    end;
  {$endif}
  adler := 1;
  with a do begin
    s_mp_checksum(adler,@used, sizeof(used));
    s_mp_checksum(adler,@sign, sizeof(sign));
    s_mp_checksum(adler,pdigits, longint(used)*sizeof(mp_digit));
  end;
  mp_checksum := adler;
end;


{---------------------------------------------------------------------------}
procedure mp_chs(const a: mp_int; var b: mp_int);
  {-Change sign, b = -a}
begin
  {Arg check in mp_copy}
  mp_copy(a, b);
  if mp_error<>MP_OKAY then exit;
  if b.used>0 then b.sign := a.sign xor (MP_NEG xor MP_ZPOS);
end;


{---------------------------------------------------------------------------}
procedure mp_clear(var a: mp_int);
  {-Free an mp_int}
begin
  {mp_clear will be executed even if mp_error<>MP_OKAY!}
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_clear');
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
    {do nothing if a hasn't been freed previously}
    if pdigits<>nil then begin
      {first zero the digits}
      if mp_clearzero then fillchar(pdigits^, sizeof(mp_digit) * used, 0);
      {free ram}
      mp_freemem(pointer(pdigits), sizeof(mp_digit) * alloc);
      {reset members to make debugging easier}
      fillchar(a,sizeof(a),0);
      {$ifdef MPC_Diagnostic}
        dec(mp_memstat.InitDiff);
      {$endif}
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_clear2(var a,b: mp_int);
  {-Clear 2 mp_ints}
begin
  mp_clear(a);
  mp_clear(b);
end;


{---------------------------------------------------------------------------}
procedure mp_clear3(var a,b,c: mp_int);
  {-Clear 3 mp_ints}
begin
  mp_clear(a);
  mp_clear(b);
  mp_clear(c);
end;


{---------------------------------------------------------------------------}
procedure mp_clear4(var a,b,c,d: mp_int);
  {-Clear 4 mp_ints}
begin
  mp_clear(a);
  mp_clear(b);
  mp_clear(c);
  mp_clear(d);
end;


{---------------------------------------------------------------------------}
procedure mp_clear5(var a,b,c,d,e: mp_int);
  {-Clear 5 mp_ints}
begin
  mp_clear(a);
  mp_clear(b);
  mp_clear(c);
  mp_clear(d);
  mp_clear(e);
end;


{---------------------------------------------------------------------------}
procedure mp_clear6(var a,b,c,d,e,f: mp_int);
  {-Clear 6 mp_ints}
begin
  mp_clear3(a,b,c);
  mp_clear3(d,e,f);
end;


{---------------------------------------------------------------------------}
procedure mp_clear7(var a,b,c,d,e,f,g: mp_int);
  {-Clear 7 mp_ints}
begin
  mp_clear4(a,b,c,d);
  mp_clear3(e,f,g);
end;


{---------------------------------------------------------------------------}
procedure mp_clear8(var a,b,c,d,e,f,g,h: mp_int);
  {-Clear 8 mp_ints}
begin
  mp_clear4(a,b,c,d);
  mp_clear4(e,f,g,h);
end;


{---------------------------------------------------------------------------}
procedure mp_clear9(var a,b,c,d,e,f,g,h,i: mp_int);
  {-Clear 9 mp_ints}
begin
  mp_clear4(a,b,c,d);
  mp_clear5(e,f,g,h,i);
end;


{---------------------------------------------------------------------------}
procedure mp_clear_multi(var vi: array of mp_int);
  {-Clear a vector of mp_ints}
var
  i: integer;
begin
  for i:=low(vi) to high(vi) do mp_clear(vi[i]);
end;


{---------------------------------------------------------------------------}
procedure mp_clear_multi_p(const pv: array of pmp_int);
  {-Clear a list of mp_ints given as a ptr vector}
var
  i: integer;
begin
  for i:=low(pv) to high(pv) do mp_clear(pv[i]^);
end;


{---------------------------------------------------------------------------}
procedure mp_clrbit(var a: mp_int; n: longint);
  {-Clear bit n of a, no action if out of range, (1 = bit 0)}
var
  d,k: TNInt;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_clrbit');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if (n<0) or (n>MP_MAXBIT) then exit;
  {$ifdef BASM16}
    asm
       mov ax, word ptr [n]
       mov dx, word ptr [n+2]
       mov cx, MP_DIGIT_BIT
       div cx
       mov [k],ax
       mov [d],dx
    end;
  {$else}
    k := n div MP_DIGIT_BIT;
    d := mp_digit(n mod MP_DIGIT_BIT);
  {$endif}
  with a do begin
    if (pdigits<>nil) and (used>0) and (k<used) then begin
      pdigits^[k] := pdigits^[k] and ((not (mp_digit(1) shl d)) and MP_MASK);
      if k=used-1 then mp_clamp(a);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function mp_cmp(const a,b: mp_int): integer;
  {-Compare two mp_ints (signed), return sign(a-b)}
begin
  mp_cmp := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cmp');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {compare based on signs}
  if a.sign<>b.sign then begin
    if a.sign=MP_NEG then mp_cmp := -1 else mp_cmp := 1;
    exit;
  end;
  {compare digits}
  if a.sign=MP_NEG then begin
    {if negative, compare opposite direction}
    mp_cmp := mp_cmp_mag(b, a);
  end
  else mp_cmp := mp_cmp_mag(a, b);
end;


{---------------------------------------------------------------------------}
function mp_cmp_d(const a: mp_int; b: mp_digit): integer;
  {-Compare a with an mp_digit, return sign(a-b)}
var
  a0: mp_digit;
begin
  mp_cmp_d := -1;   {default result, used if a<0}
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cmp_d');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if b>MP_DIGIT_MAX then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_cmp_d: b>MP_DIGIT_MAX');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}
  if (mp_error<>MP_OKAY) or (a.sign=MP_NEG) then exit;
  if a.used>1 then mp_cmp_d := 1 {compare based on magnitude}
  else begin
    {compare the only digit of a to b}
    if a.used=1 then a0 := a.pdigits^[0] else a0:=0;
    if a0>b then mp_cmp_d := 1
    else if a0<b then mp_cmp_d := -1
    else mp_cmp_d := 0;
  end;
end;


{---------------------------------------------------------------------------}
function mp_cmp_int(const a: mp_int; b: longint): integer;
  {-Compare a with a longint, return sign(a-b)}
var
  t: mp_int;
begin
  mp_cmp_int := 0;
  if mp_error<>MP_OKAY then exit;
  {MPC_ArgCheck in mp_cmp}
  mp_init_set_int(t, b);
  if mp_error=MP_OKAY then begin
    mp_cmp_int := mp_cmp(a, t);
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
function mp_cmp_mag(const a,b: mp_int): integer;
  {-Compare magnitude of two mp_ints (unsigned), return sign(|a|-|b|)}
var
  pa,pb: pmp_digit;
  i: TNInt;
begin
  {Value for last alternative, keep D6+ happy}
  mp_cmp_mag := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cmp_mag');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {compare based on # of non-zero digits}
  if a.used > b.used then begin
    mp_cmp_mag := 1;
    exit;
  end
  else if a.used < b.used then begin
    mp_cmp_mag := -1;
    exit;
  end
  else if a.used>0 then begin
    {a.used=b.used, if a.used=0 then both are 0}
    pa := @a.pdigits^[a.used-1];
    pb := @b.pdigits^[b.used-1];
    for i:=1 to a.used do begin
      {look for different highest digit}
      if pa^ <> pb^ then begin
        if pa^ < pb^ then mp_cmp_mag := -1
        else mp_cmp_mag := 1;
        exit;
      end;
      dec(pa);
      dec(pb);
    end;
    {all digits are equal, def. result 0 applies}
  end;
end;


{---------------------------------------------------------------------------}
function mp_cmp_mag_d(const a: mp_int; b: mp_digit): integer;
  {-Compare |a| with a digit, return sign(|a|-b)}
var
  a0: mp_digit;
begin
  mp_cmp_mag_d := 1;   {default result, used if a.used>1}
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cmp_mag_d');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if b>MP_DIGIT_MAX then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_cmp_mag_d: b>MP_DIGIT_MAX');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}
  if (mp_error<>MP_OKAY) or (a.used>1) then exit;
  {compare the only digit of a to b}
  if a.used=1 then a0 := a.pdigits^[0] else a0:=0;
  if a0>b then mp_cmp_mag_d := 1
  else if a0<b then mp_cmp_mag_d := -1
  else mp_cmp_mag_d := 0;
end;


{---------------------------------------------------------------------------}
function mp_cnt_lsb(const a: mp_int): longint;
  {-Count the number of least significant bits which are zero}
var
  x: longint;
  q, qq: mp_digit;
const
  lnz: array[0..15] of word = (4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0);
begin
  mp_cnt_lsb := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_cnt_lsb');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {easy out}
  if mp_iszero(a) then exit;

  {scan lower digits until non-zero}
  x := 0;
  while (x<a.used) and (a.pdigits^[x]=0) do inc(x);

  q := a.pdigits^[x];
  x := x*DIGIT_BIT;

  {now scan this digit until a 1 is found}
  if q and 1 = 0 then begin
    repeat
      qq := q and 15;
      inc(x,lnz[qq]);
      q := q shr 4;
    until qq<>0;
  end;
  mp_cnt_lsb := x;
end;


{---------------------------------------------------------------------------}
procedure mp_copy(const a: mp_int; var b: mp_int);
  {-Copy an mp_int, b = a}
var
  i: TNInt;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_copy');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {if dst=src do nothing}
  if @a=@b then exit;
  {grow dest}
  if b.alloc < a.used then begin
    mp_grow(b, a.used);
    if mp_error<>MP_OKAY then exit;
  end;
  move(a.pdigits^, b.pdigits^, a.used*sizeof(mp_digit));
  if b.used>a.used then for i:=a.used to b.used-1 do b.pdigits^[i]:=0;
  b.used := a.used;
  b.sign := a.sign;
end;


{---------------------------------------------------------------------------}
procedure mp_dec(var a: mp_int);
  {-Decrement an mp_int by 1}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_dec');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mp_sub_d(a,1,a);
end;


{---------------------------------------------------------------------------}
function mp_decimal(const a: mp_int): mp_string;
  {-Convert to decimal, max 255 digits}
begin
  mp_decimal := mp_radix_str(a, 10);
end;


{---------------------------------------------------------------------------}
procedure mp_dec_int(var a: mp_int; b: longint);
  {-Calculate a = a - b}
begin
  mp_sub_int(a,b,a);
end;


{---------------------------------------------------------------------------}
procedure mp_div(const a,b: mp_int; var c: mp_int);
  {-Integer signed division, c = a div b}
begin
  mp_divrem(a,b,@c,nil);
end;


{---------------------------------------------------------------------------}
procedure mp_divrem(const a,b: mp_int; pc,pd: pmp_int);
  {-Integer signed division, pc^ = a div b, pd^ = a rem b; sign(pd^)=sign(a)}
begin
  if pc=nil then begin
    if pd<>nil then mp_mod(a,b,pd^); {xx}
    exit;
  end;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or ((pc<>nil) and mp_not_init(pc^)) or ((pd<>nil) and mp_not_init(pd^)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_divrem');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if pc=pd then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_divrem: pc=pd');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
    if (@a=pd) and (@b=pc) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_divrem: (@a=pd) and (@b=pc)');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  {$endif}
  {Check b=0}
  if b.used=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_divrem: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  {call Burnikel/Ziegler or Knuth's basecase algorithm D}
  s_mp_divrem(a,b,pc,pd);
end;


{---------------------------------------------------------------------------}
procedure mp_div_2k(const a: mp_int; b: longint; var c: mp_int; pd: pmp_int);
  {-Divide by 2^b; quotient in c, optional remainder in pd^, sign(pd^)=sign(a)}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(c) or ((pd<>nil) and mp_not_init(pd^)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_div_2k');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if a=0 or b<=0 then we do no work}
  if (b<=0) or (a.used=0) then begin
    mp_copy(a, c);
    if pd<>nil then mp_zero(pd^);
    exit;
  end;

  mp_init(t);
  if mp_error<>MP_OKAY then exit;
  {get the remainder}
  if pd<>nil then begin
    s_mp_mod_2k(a, b, t);
    if mp_error<>MP_OKAY then begin
      mp_clear(t);
      exit;
    end;
  end;
  {get the quotient via mp_shr}
  mp_shr(a,b,c);
  if pd<>nil then mp_exch(t, pd^);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_div_2(const a: mp_int; var b: mp_int);
  {-Divide by 2, b = a/2}
var
  oldused: longint;
  i: TNInt;
  r,rr: mp_digit;
  pa,pb: pmp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_div_2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if a=0 then just copy}
  if a.used=0 then begin
    mp_copy(a, b);
    exit;
  end;

  {grow dest}
  if b.alloc < a.used then begin
    mp_grow(b, a.used);
    if mp_error<>MP_OKAY then exit;
  end;

  oldused := b.used;
  b.used  := a.used;
  pa := @a.pdigits^[b.used-1];
  pb := @b.pdigits^[b.used-1];
  {carry}
  r  := 0;
  for i:=b.used-1 downto 0 do begin
    {get the carry for the next iteration}
    rr := pa^ and 1;
    {shift the current digit, add in carry and store}
    pb^ := (pa^ shr 1) or (r shl (DIGIT_BIT - 1));
    {forward carry to next iteration}
    r := rr;
    dec(pa);
    dec(pb);
  end;

  {zero excess digits}
  pb := @b.pdigits^[b.used];
  for i:=b.used to oldused-1 do begin
    pb^ := 0;
    inc(pb);
  end;
  b.sign := a.sign;
  mp_clamp(b);
end;


{---------------------------------------------------------------------------}
procedure mp_shr1(var a: mp_int);
  {-Divide a by 2, a = a/2}
var
  i: TNInt;
  r,rr: mp_digit;
  pa: pmp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_shr1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if a.used=0 then exit;

  pa := @a.pdigits^[pred(a.used)];
  {carry}
  r  := 0;
  for i:=pred(a.used) downto 0 do begin
    {get the carry for the next iteration}
    rr := pa^ and 1;
    {shift the current digit, add in carry and store}
    pa^ := (pa^ shr 1) or (r shl (DIGIT_BIT - 1));
    {forward carry to next iteration}
    r := rr;
    dec(pa);
  end;
  mp_clamp(a);
end;


{---------------------------------------------------------------------------}
procedure s_mp_div_d(const a: mp_int; b: mp_digit; pc: pmp_int; var d: mp_digit);
  {-Single digit division, pc^=sign(a)(|a| div b), r = |a| mod b, no init check}
var
  q: mp_int;
  pa,pq: pmp_digit;
  i2: integer;
  ix: TNInt;
  {$ifndef BASM16}
    w: mp_word;
    t: mp_digit;
  {$endif}
begin
  {cannot divide by zero}
  if b=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('s_mp_div_d: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {*tbd: div 3?}

  {quick outs}
  if (b=1) or (a.used=0) then begin
    d := 0;
    if pc<>nil then mp_copy(a, pc^);
    exit;
  end;

  {power of two?}
  if mp_is_pow2_d(b, i2) then begin
    d := a.pdigits^[0] and pred(b);
    if pc<>nil then mp_shr(a, i2, pc^);
    exit;
  end;

  {no easy answer [c'est la vie]. Just division}
  mp_init_size(q, a.used);
  if mp_error <> MP_OKAY then exit;

  q.used := a.used;
  q.sign := a.sign;
  {a<>0 -> a.used>0}

  ix := pred(a.used);

  pa := pmp_digit(@a.pdigits^[ix]);
  pq := pmp_digit(@q.pdigits^[ix]);

  {$ifdef BASM16}
    asm
           mov   cx,[ix]
           push  ds
           lds   si,[pa]
           les   di,[pq]
           sub   ax,ax              {w=0}
           mov   bx,[b]
           std
           mov   dx,ax
      @@1: shr   dx,16-DIGIT_BIT    {before shr: ax=dx=w}
           shl   ax,DIGIT_BIT
           add   ax,[si]            {w   := (w shl DIGIT_BIT) or pa^}
           sub   si,2
           div   bx
           stosw                    {pa^ := w div b}
           mov   ax,dx              {w   := w mod b}
           dec   cx
           jns   @@1
           cld
           lds   si,[d]
           mov   [si],ax
           pop   ds
    end;
  {$else}
    w := 0;
    for ix:= pred(a.used) downto 0 do begin
      w := (w shl DIGIT_BIT) or pa^;
      if w >= b then begin
        t := w div b;
        w := w - mp_word(t) * mp_word(b); {w := w mod b}
        {The following is buggy for FPC242,244 with debug!!}
        {dec(w, mp_word(t) * mp_word(b));}
      end
      else t:=0;
      pq^ := t;
      dec(pa);
      dec(pq);
    end;
    d := mp_digit(w);
  {$endif}
  if pc<>nil then begin
    {*0.9.00: clamp only if needed}
    mp_clamp(q);
    mp_exch(q, pc^);
  end;
  mp_clear(q);
end;


{---------------------------------------------------------------------------}
procedure mp_div_d(const a: mp_int; b: mp_digit; pc: pmp_int; var d: mp_digit);
  {-Single digit division, pc^=sign(a)(|a| div b), r = |a| mod b}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or ((pc<>nil) and mp_not_init(pc^)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_div_d');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if b>MP_DIGIT_MAX then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_div_d: b>MP_DIGIT_MAX');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}
  s_mp_div_d(a,b,pc,d);
end;


{---------------------------------------------------------------------------}
procedure mp_div_int(const a: mp_int; b: longint; pc: pmp_int; var d: longint);
  {-Integer signed division, pc^ = a div b, d = a rem b; sign(d)=sign(a)}
var
  t: mp_int;
  r: mp_digit;
begin
  if (b>0) and (b<=lv_digit_max) then begin
    {Although mp_divrem finally reaches special code for small b in   }
    {s_mp_divrem_basecase, avoid that overhead and mp_init if possible}
    mp_div_d(a,b,pc,r);
    d := r;
    if a.SIGN=MP_NEG then d := -d;
  end
  else begin
    mp_init_set_int(t,b);
    {Debug checking in called functions}
    if mp_error=MP_OKAY then begin
      s_mp_divrem_basecase(a,t,pc,@t);
      {mp_divrem(a,t,pc,@t);}
      d := mp_get_int(t);
      mp_clear(t);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_div_w(const a: mp_int; b: word; pc: pmp_int; var r: word);
  {-Divide a by a single word b, pc^=sign(a)(|a| div b), r = |a| mod b}
var
  q: mp_int;
  pa,pq: pmp_digit;
  ix: TNInt;
  i2: integer;
  {$ifndef BASM16}
    w: mp_word;
    t: mp_digit;
  {$endif}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or ((pc<>nil) and mp_not_init(pc^)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_div_w');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {cannot divide by zero}
  if b=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_div_w: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {check for trivial cases}
  if (b=1) or (a.used=0) then begin
    r := 0;
    if pc<>nil then mp_copy(a, pc^);
    exit;
  end;

  {power of two?}
  if IsPow2_w(b, i2) then begin
    {$ifdef MP_32BIT}
      {mp_digit is at least 16 bit, so a.pdigits^[0] has all bits for masking}
      r := a.pdigits^[0] and pred(b);
    {$else}
      if (a.used=1) or (i2<=DIGIT_BIT) then r := a.pdigits^[0] and pred(b)
      else begin
        {DIGIT_BIT in 7..15, so lowest two digit have the 16 bits needed}
        r := ((mp_word(a.pdigits^[1]) shl DIGIT_BIT) or a.pdigits^[0]) and pred(b);
      end;
    {$endif}
    if pc<>nil then mp_shr(a, i2, pc^);
    exit;
  end;

  mp_init_size(q, a.used);
  if mp_error <> MP_OKAY then exit;

  q.used := a.used;
  q.sign := a.sign;

  ix := pred(a.used);
  pa := pmp_digit(@a.pdigits^[ix]);
  pq := pmp_digit(@q.pdigits^[ix]);

  {$ifdef BASM16}
    asm
           mov   cx,[ix]
           push  ds
           lds   si,[pa]
           les   di,[pq]
           sub   ax,ax              {w=0}
           mov   bx,[b]
           std
      @@1: mov   dx,ax
           shr   dx,16-DIGIT_BIT
           shl   ax,DIGIT_BIT
           add   ax,[si]            {w   := (w shl DIGIT_BIT) or pa^}
           sub   si,2
           div   bx
           stosw                    {pa^ := w div b}
           mov   ax,dx              {w   := w mod b}
           dec   cx
           jns   @@1
           cld
           lds   si,[r]
           mov   [si],ax
           pop   ds
    end;
  {$else}
    w := 0;
    {w<b is loop invariant}
    for ix:=a.used-1 downto 0 do begin
      w := (w shl DIGIT_BIT) or pa^;
      {w < b*2^DIGIT_BIT}
      if w >= b then begin
        {t < b*2^DIGIT_BIT div b = 2^DIGIT_BIT}
        {w<b and t<2^DIGIT_BIT}
        t := w div b;
        dec(w, mp_word(t) * mp_word(b)); {w := w mod b}
      end
      else begin
        t:=0;
        {w<b and t<2^DIGIT_BIT}
      end;
      pq^ := t;
      dec(pa);
      dec(pq);
    end;
    {still w<b}
    r := word(w);
  {$endif}

  if pc<>nil then begin
    {*0.9.00: clamp only if needed}
    mp_clamp(q);
    mp_exch(q, pc^);
  end;
  mp_clear(q);
end;


{---------------------------------------------------------------------------}
procedure mp_exch(var a,b: mp_int);
  {-Exchange two mp_ints}
var
  t: mp_int;
begin
  if mp_error=MP_OKAY then begin
    t := a;
    a := b;
    b := t;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_expt(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a^b, b>=0}
label
  __R,__X;
var
  r,x: mp_int;
  d: mp_digit;
  dig, dh: TNInt;
  bit: integer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_expt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {easy outs}
  if b.sign=MP_NEG then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_expt: b < 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  if mp_is0(a) then begin
    if mp_is0(b) then mp_set1(c) else mp_zero(c);
    exit;
  end;

  mp_init_copy(x,a);
  if mp_error<>MP_OKAY then exit;

  mp_init_set(r,1);
  if mp_error<>MP_OKAY then goto __X;

  {if b=0 then c=1, do nothing}
  if b.used>0 then begin
    {Loop over low-order digits in ascending order }
    {*we: assume b is clamped, special treatment of}
    {highest digit to avoid unnecessary squaring   }
    dh := b.used-1;
    for dig:=0 to dh-1 do begin
      d := b.pdigits^[dig];
      for bit:=1 to DIGIT_BIT do begin
        if odd(d) then begin
          mp_mul(r,x,r);
          if mp_error<>MP_OKAY then goto __R;
        end;
        mp_sqr(x,x);
        if mp_error<>MP_OKAY then goto __R;
        d := d shr 1;
      end;
    end;
    {use only minimum nonzero bits in last digit}
    d := b.pdigits^[dh];
    while d<>0 do begin
      if odd(d) then mp_mul(r,x,r);
      d := d shr 1;
      {skip final sqr}
      if d<>0 then mp_sqr(x,x);
      if mp_error<>MP_OKAY then goto __R;
    end;
  end;
  mp_copy(r,c);

__R: mp_clear(r);
__X: mp_clear(x);
end;


{---------------------------------------------------------------------------}
procedure mp_expt_w(const a: mp_int; b: word; var c: mp_int);
  {-Calculate c = a^b}
begin
  mp_expt_int(a,b,c);
end;


{---------------------------------------------------------------------------}
procedure mp_expt_int(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a^b, b>=0}
var
  x: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_expt_int');
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
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_expt_int: b < 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {easy outs}
  if b=0 then begin
    mp_set1(c);
    exit;
  end
  else if mp_iszero(a) then begin
    mp_zero(c);
    exit;
  end
  else if b=1 then begin
    mp_copy(a,c);
    exit;
  end
  else if b=2 then begin
    mp_sqr(a,c);
    exit;
  end;

  {make local copy of a}
  mp_init_copy(x,a);
  if mp_error<>MP_OKAY then exit;

  mp_set1(c);
  {b >2, therefore loop is executed at least once}
  while mp_error=MP_OKAY do begin
    if odd(b) then mp_mul(c,x,c);
    b := b shr 1;
    if b=0 then break;
    mp_sqr(x,x);
  end;
  mp_clear(x);
end;


{---------------------------------------------------------------------------}
function mp_gcd_int(const a: mp_int; b: longint): longint;
  {-Return gcd(a,b), b<>0}
var
  c: longint;
begin
  mp_mod_int(a,b,c);
  if mp_error=MP_OKAY then mp_gcd_int := gcd32(b,c)
  else mp_gcd_int := 0;
end;


{---------------------------------------------------------------------------}
function mp_get_int(const a: mp_int): longint;
  {-Get the lower signed 31 bits of an mp_int}
var
  i: integer;
  res: longint;
begin
  mp_get_int:=0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_get_int');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {get number of digits of the lsb we have to read}
  i := (31+DIGIT_BIT-1) div DIGIT_BIT;
  if i>a.used then i:=a.used;
  dec(i);
  if i<0 then exit;

  {get most significant digit of result}
  res := a.pdigits^[i];

  while i>0 do begin
    dec(i);
    {$ifdef BIT16}
      LeftShiftAdd(res,a.pdigits^[i]);
    {$else}
      res := (res shl DIGIT_BIT) or longint(a.pdigits^[i]);
    {$endif}
  end;

  res := res and $7FFFFFFF;
  if a.sign=MP_NEG then res := -res;
  mp_get_int := res;
end;


{---------------------------------------------------------------------------}
procedure mp_gr_mod(var x: mp_int; const N, R: mp_int);
  {-Reduce x to x mod N, N > 1, using generalized reciprocal iteration.}
  { Result is >= 0. R is from mp_gr_setup. Compared to the similar}
  { Barrett reduction the restricted range 0<x<N^2 is not required.}
var
  d: mp_int;
  s: longint;
  cmp: integer;
  sigx: word;
begin
  {$ifdef MPC_ArgCheck}
    if mp_not_init(x) or mp_not_init(N) or mp_not_init(R) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_gr_mod');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mp_cmp_d(N,2)=MP_LT then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_gr_mod: N < 2');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  if s_mp_is_le0(R) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_gr_mod: R <= 0');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {Ref: Crandall/Pomerance Alg. 9.2.10 and Marcel Martin's IReciprocalMod}
  cmp := mp_cmp_mag(x,N);
  if cmp=MP_LT then begin
    {abs(x)<N: no reduction just make result positiv}
    if x.sign=MP_NEG then mp_add(x,N,x);
  end
  else if cmp=MP_EQ then begin
    {abs(x)=N, result is zero}
    mp_zero(x);
  end
  else begin
    {the reduction case: abs(x) > N}
    mp_init(d);
    if mp_error=MP_OKAY then begin
      {remember sign(x) and make x positive}
      sigx := x.sign;
      x.sign := MP_ZPOS;
      s := mp_bitsize(R)-1;
      {loop to overcome the Barrett restriction x < N^2-1}
      repeat
        {use Bosselaers' form which is better if div is not needed}
        {x mod N = x - N*floor(R*floor(x/2^(s-1))/2^(s+1))}
        mp_shr(x,s-1,d);
        mp_mul(d,R,d);
        mp_shr(d,s+1,d);
        mp_mul(d,N,d);
        mp_sub(x,d,x);
        if mp_is_lt(x,N) then break else mp_sub(x,N,x);
        if mp_is_lt(x,N) then break else mp_sub(x,N,x);
      until mp_is_lt(x,N) or (mp_error<>MP_OKAY);
      {if original x was negative, return (-x) residue class}
      if (x.used>0) and (sigx=MP_NEG) then mp_sub(N,x,x);
      mp_clear(d);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_gr_setup(var RN: mp_int; const N: mp_int);
  {-Calculate the generalized reciprocal for N>1, @RN<>@N}
var
  b: longint;
  d: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(N) or mp_not_init(RN) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_gr_setup');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Check N>1 and @N<>@RN}
  if mp_cmp_d(N,2)=MP_LT then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_gr_setup: N < 2');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;
  if @N=@RN then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_gr_setup: @N=@RN');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {RN := 4^bitsize(N-1)/N, c.f. Crandall/Pomerance [10], Def. 9.2.7}
  {avoid calculation of N-1: bitsize(N)=bitsize(N-1) if N<>2^k}
  b := mp_bitsize(N);
  {if N is power of 2 then decrement b}
  d := N.pdigits^[pred(N.used)];
  if (d and (d-1) = 0) and (mp_popcount(N)=1) then dec(b);
  mp_2expt(RN,b shl 1);
  {C/P optimize further by using discrete Newton iteration without}
  {division (only shifts etc). Here a simple div is used instead.}
  mp_div(RN,N,RN);
end;


{---------------------------------------------------------------------------}
procedure mp_grow(var a: mp_int; size: longint);
  {-Grow an mp_int to a given size (new part is zerofilled)}
var
  np: PDigitArray;
  xs: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_grow');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {if the alloc size is smaller alloc more ram}
  if a.alloc<size then begin
    if size>MAXDigits then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXMaxDigits.Create('mp_grow');
        {$else}
          RunError(MP_RTE_OTHER);
        {$endif}
      {$else}
        set_mp_error(MP_MAXDIGITS);
        exit;
      {$endif}
    end
    else begin
      {ensure there are always at least mp_allocprec digits extra on top}
      xs := mp_allocprec + (mp_allocprec - (size and mp_allocmask));
      inc(size, xs);
      {Note: new part is zerofilled}
      np := mp_realloc(a.pdigits, sizeof(mp_digit)*a.alloc, sizeof(mp_digit)*size);
      if np=nil then begin
        {$ifdef MPC_HaltOnError}
          {$ifdef MPC_UseExceptions}
            raise MPXMemory.Create('mp_grow: realloc');
          {$else}
            RunError(MP_RTE_MEM);
          {$endif}
        {$else}
          set_mp_error(MP_MEM);
          exit;
        {$endif}
      end
      else begin
        a.pdigits := np;
        a.alloc := size;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_shrink(var a: mp_int);
  {-Shrink RAM allocated for an mp_int}
var
  tmp: pointer;
  ns : longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_shrink');
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
    if (pdigits<>nil) and (alloc>used) then begin
      if used=0 then ns := mp_allocprec
      else ns := used;
      tmp := mp_realloc(pdigits, sizeof(mp_digit)*alloc, sizeof(mp_digit)*ns);
      if tmp=nil then begin
        {$ifdef MPC_HaltOnError}
          {$ifdef MPC_UseExceptions}
            raise MPXMemory.Create('s_mp_shrink: realloc');
          {$else}
            RunError(MP_RTE_MEM);
          {$endif}
        {$else}
          set_mp_error(MP_MEM);
          exit;
        {$endif}
      end;
      pdigits := tmp;
      alloc   := ns;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function mp_hex(const a: mp_int): mp_string;
  {-Convert to hex string, max 255 digits}
begin
  mp_hex := mp_radix_str(a, 16);
end;


{---------------------------------------------------------------------------}
procedure mp_inc(var a: mp_int);
  {-Increment an mp_int by 1}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_inc');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mp_add_d(a,1,a);
end;


{---------------------------------------------------------------------------}
procedure mp_inc_int(var a: mp_int; b: longint);
  {-Calculate a = a + b}
begin
  mp_add_int(a,b,a);
end;


{---------------------------------------------------------------------------}
procedure mp_init(var a: mp_int);
  {-Initialize an mp_int}
begin
  mp_init_prim(a, 0);
end;


{---------------------------------------------------------------------------}
procedure mp_init2(var a,b: mp_int);
  {-Initialize 2 mp_ints}
var
  pa: array[0..1] of pmp_int;
begin
  pa[0] := @a;
  pa[1] := @b;
  mp_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mp_init3(var a,b,c: mp_int);
  {-Initialize 3 mp_ints}
var
  pa: array[0..2] of pmp_int;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  mp_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mp_init4(var a,b,c,d: mp_int);
  {-Initialize 4 mp_ints}
var
  pa: array[0..3] of pmp_int;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  mp_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mp_init5(var a,b,c,d,e: mp_int);
  {-Initialize 5 mp_ints}
var
  pa: array[0..4] of pmp_int;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  pa[4] := @e;
  mp_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mp_init6(var a,b,c,d,e,f: mp_int);
  {-Initialize 6 mp_ints}
var
  pa: array[0..5] of pmp_int;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  pa[4] := @e;
  pa[5] := @f;
  mp_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mp_init7(var a,b,c,d,e,f,g: mp_int);
  {-Initialize 7 mp_ints}
var
  pa: array[0..6] of pmp_int;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  pa[4] := @e;
  pa[5] := @f;
  pa[6] := @g;
  mp_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mp_init8(var a,b,c,d,e,f,g,h: mp_int);
  {-Initialize 8 mp_ints}
var
  pa: array[0..7] of pmp_int;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  pa[4] := @e;
  pa[5] := @f;
  pa[6] := @g;
  pa[7] := @h;
  mp_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mp_init9(var a,b,c,d,e,f,g,h,i: mp_int);
  {-Initialize 9 mp_ints}
var
  pa: array[0..8] of pmp_int;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  pa[4] := @e;
  pa[5] := @f;
  pa[6] := @g;
  pa[7] := @h;
  pa[8] := @i;
  mp_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mp_init_copy(var a: mp_int; const b: mp_int);
  {-Create a, then copy b into it}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_init_copy');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_init_size(a,b.used);
  mp_copy(b, a);
end;


{---------------------------------------------------------------------------}
procedure mp_init_multi(var vi: array of mp_int);
  {-Initialize a vector of mp_ints}
  { on error the already initialized mp_ints will be cleared}
var
  i,k: integer;
begin
  if mp_error<>MP_OKAY then exit;
  for i:=low(vi) to high(vi) do begin
    mp_init_prim(vi[i], 0);
    if mp_error<>MP_OKAY then begin
      {error, clear all previous mp_ints}
      for k:=low(vi) to i-1 do mp_clear(vi[k]);
      break;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_init_multi_p(var pv: array of pmp_int);
  {-Initialize a list of mp_ints given as a ptr vector}
  { on error the already initialized mp_ints will be cleared}
var
  i,k: integer;
begin
  if mp_error<>MP_OKAY then exit;
  for i:=low(pv) to high(pv) do begin
    mp_init_prim(pv[i]^, 0);
    if mp_error<>MP_OKAY then begin
      {error, clear all previous mp_ints}
      for k:=low(pv) to i-1 do mp_clear(pv[k]^);
      break;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_init_set(var a: mp_int; b: mp_digit);
  {-Initialize and set a digit}
begin
  if mp_error<>MP_OKAY then exit;
  mp_init(a);
  if mp_error=MP_OKAY then begin
    mp_set(a,b);
    if mp_error<>MP_OKAY then mp_clear(a);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_init_set_int(var a: mp_int; b: longint);
  {-Initialize and set a to a longint}
begin
  if mp_error<>MP_OKAY then exit;
  mp_init(a);
  if mp_error=MP_OKAY then begin
    mp_set_int(a,b);
    if mp_error<>MP_OKAY then mp_clear(a);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_init_size(var a: mp_int; size: longint);
  {-Initialize a to size digits, rounded up to multiple of mp_allocprec}
begin
  if mp_error<>MP_OKAY then exit;
  {pad size so there are always extra digits}
  if size>MAXDigits then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMaxDigits.Create('mp_init_size');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$else}
      set_mp_error(MP_MAXDIGITS);
      exit;
    {$endif}
  end;
  mp_init_prim(a, size);
end;


{---------------------------------------------------------------------------}
procedure mp_init_size2(var a,b: mp_int; size: longint);
  {-Initialize a and b to size digits, rounded up to multiple of mp_allocprec}
begin
  if mp_error<>MP_OKAY then exit;
  {pad size so there are always extra digits}
  if size>MAXDigits then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMaxDigits.Create('mp_init_size2');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$else}
      set_mp_error(MP_MAXDIGITS);
      exit;
    {$endif}
  end;
  mp_init_prim(a, size);
  if mp_error=MP_OKAY then begin
    mp_init_prim(b, size);
    if mp_error<>MP_OKAY then mp_clear(a);
  end;
end;


{---------------------------------------------------------------------------}
function mp_isbit(const a: mp_int; n: longint): boolean;
  {-Test if bit n of a is set, (1 = bit 0)}
var
  d,k: TNInt;
begin
  mp_isbit := false;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_isbit');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if (n<0) or (n>MP_MAXBIT) then exit;
  {$ifdef BASM16}
    asm
       mov ax, word ptr [n]
       mov dx, word ptr [n+2]
       mov cx, MP_DIGIT_BIT
       div cx
       mov [k],ax
       mov [d],dx
    end;
  {$else}
    k := n div MP_DIGIT_BIT;
    d := mp_digit(n mod MP_DIGIT_BIT);
  {$endif}
  with a do begin
    mp_isbit := (pdigits<>nil) and (used>0) and (k<used) and odd(pdigits^[k] shr d);
  end;
end;


{---------------------------------------------------------------------------}
function mp_iseven(const a: mp_int): boolean;
  {-Initialized and even}
begin
  with a do mp_iseven := (mp_error=MP_OKAY) and (magic=MP_MAGIC) and (pdigits<>nil) and ((used=0) or (not odd(pdigits^[0])));
end;


{---------------------------------------------------------------------------}
function mp_isodd(const a: mp_int): boolean;
  {-Initialized and odd}
begin
  with a do mp_isodd := (mp_error=MP_OKAY) and (magic=MP_MAGIC) and (pdigits<>nil) and (used>0) and odd(pdigits^[0]);
end;


{---------------------------------------------------------------------------}
function mp_iszero(const a: mp_int): boolean;
  {-Initialized and zero}
begin
  with a do mp_iszero := (mp_error=MP_OKAY) and (magic=MP_MAGIC) and (pdigits<>nil) and (used=0);
end;


{---------------------------------------------------------------------------}
function mp_is0(const a: mp_int): boolean;
  {-Initialized and = 0}
begin
  with a do mp_is0 := (mp_error=MP_OKAY) and (magic=MP_MAGIC) and (pdigits<>nil) and (used=0);
end;


{---------------------------------------------------------------------------}
function mp_is1(const a: mp_int): boolean;
  {-Initialized and a = 1}
begin
  with a do begin
    mp_is1 := (mp_error=MP_OKAY) and (magic=MP_MAGIC) and (sign=MP_ZPOS) and (pdigits<>nil) and (used=1) and (pdigits^[0]=1);
  end;
end;


{---------------------------------------------------------------------------}
function mp_is1a(const a: mp_int): boolean;
  {-Initialized and abs(a) = 1}
begin
  with a do mp_is1a := (mp_error=MP_OKAY) and (magic=MP_MAGIC) and (pdigits<>nil) and (used=1) and (pdigits^[0]=1);
end;


{---------------------------------------------------------------------------}
function mp_is_eq(const a,b: mp_int): boolean;
  {-Return a = b}
begin
  mp_is_eq := (mp_error=MP_OKAY) and (mp_cmp(a,b)=MP_EQ);
end;


{---------------------------------------------------------------------------}
function mp_is_ge(const a,b: mp_int): boolean;
  {-Return a >= b}
begin
  mp_is_ge := (mp_error=MP_OKAY) and (mp_cmp(a,b)<>MP_LT);
end;


{---------------------------------------------------------------------------}
function mp_is_gt(const a,b: mp_int): boolean;
  {-Return a > b}
begin
  mp_is_gt := (mp_error=MP_OKAY) and (mp_cmp(a,b)=MP_GT);
end;


{---------------------------------------------------------------------------}
function mp_is_le(const a,b: mp_int): boolean;
  {-Return a <= b}
begin
  mp_is_le := (mp_error=MP_OKAY) and (mp_cmp(a,b)<>MP_GT);
end;


{---------------------------------------------------------------------------}
function mp_is_lt(const a,b: mp_int): boolean;
  {-Return a < b}
begin
  mp_is_lt := (mp_error=MP_OKAY) and (mp_cmp(a,b)=MP_LT);
end;


{---------------------------------------------------------------------------}
function mp_is_ne(const a,b: mp_int): boolean;
  {-Return a <> b}
begin
  mp_is_ne := (mp_error=MP_OKAY) and (mp_cmp(a,b)<>MP_EQ);
end;


{---------------------------------------------------------------------------}
function mp_is_longint(const a: mp_int; var b: longint): boolean;
  {-Test if a fits into longint, if true set b := a}
const
  DC_MAXLONG = (30+DIGIT_BIT) div DIGIT_BIT;
var
  i,u: TNInt;
begin
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_longint');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        mp_is_longint := false;
        exit;
      {$endif}
    end;
  {$endif}
  u := a.used;
  if u=0 then b := 0
  else begin
    if u=1 then b := longint(a.pdigits^[0])
    else begin
      if (u>DC_MAXLONG) or (u=DC_MAXLONG) and (mp_bitsize(a)>31) then begin
        mp_is_longint := false;
        exit;
      end;
      b := 0;
      for i:=pred(u) downto 0 do begin
        {$ifdef BASM16}
          LeftShiftAdd(b,a.pdigits^[i]);
        {$else}
          b := (b shl DIGIT_BIT) or longint(a.pdigits^[i]);
        {$endif}
      end;
    end;
    if a.sign=MP_NEG then b := -b;
  end;
  mp_is_longint := true;
end;


{---------------------------------------------------------------------------}
function mp_is_pow2(const a: mp_int; var n: longint): boolean;
  {-Check if |a| is a power of 2, if true, return n with |a|=2^n}
var
  i: integer;
  d: mp_digit;
  p: pmp_digit;
  k,m: TNInt;
begin
  mp_is_pow2 := false;
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_is_pow2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {zero is no power of 2}
  if a.used=0 then exit;

  {m: index of highest non zero digit}
  m := pred(a.used);

  {if highest digit d is zero or more than 1 bit set, return false}
  d := a.pdigits^[m];
  if (d=0) or (d and (d-1) <> 0) then exit;

  {check for zero digits digits 0..m-1}
  if m>0 then begin
    p := pmp_digit(a.pdigits);
    for k:=0 to m-1 do begin
      if p^<>0 then exit;
      inc(p);
    end;
  end;

  {get bit index in highest digit}
  for i:= 0 to DIGIT_BIT-1 do begin
    if d=(mp_digit(1) shl i) then begin
      n := i+longint(m)*DIGIT_BIT;
      mp_is_pow2 := true;
      exit;
    end;
  end;
  {if we arrive here something is fishy, but result is false}
end;


{---------------------------------------------------------------------------}
function mp_is_pow2_d(d: mp_digit; var n: integer): boolean;
  {-Check if d is power of 2, if true, return n with d=2^n}
var
  i: integer;
begin
  mp_is_pow2_d := false;
  {First do a quick check to get the result}
  if (d=0) or (d and pred(d) <> 0) then exit;
  {now get bit index in highest digit}
  for i:= 0 to DIGIT_BIT-1 do begin
    if d=(mp_digit(1) shl i) then begin
      mp_is_pow2_d := true;
      n := i;
      exit;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_lshd2(const a: mp_int; var b: mp_int; cnt: longint);
  {-Set b to a shifted left by cnt digits}
var
  ns: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_lshd2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {if a=0 or cnt is less than zero then b := a}
  if (cnt<=0) or (a.used=0) then begin
    mp_copy(a,b);
    exit;
  end;

  ns := longint(cnt)+a.used;
  if ns > MaxDigits then begin {V1.20.04: changed <= to <}
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMaxDigits.Create('mp_lshd2');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$else}
      set_mp_error(MP_MAXDIGITS);
      exit;
    {$endif}
  end;

  b.sign := a.sign;
  {grow b to fit the new digits}
  if b.alloc < ns then begin
    mp_grow(b, ns);
    if mp_error<>MP_OKAY then exit;
  end
  else if b.used>ns then begin
    {clear high digits of b above new and old .used}
    fillchar(b.pdigits^[ns], (b.used-ns)*sizeof(mp_digit), 0);
  end;
  {shift used digits, used>0!}
  move(a.pdigits^[0], b.pdigits^[cnt], a.used*sizeof(mp_digit));
  {clear low cnt digits, cnt>0!}
  fillchar(b.pdigits^, cnt*sizeof(mp_digit), 0);
  b.used := ns;
  mp_clamp(b);
end;


{---------------------------------------------------------------------------}
procedure mp_lshd(var a: mp_int; b: longint);
  {-Shift left a certain amount of digits}
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_lshd');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {if b is less than zero return}
  {if a=0, return} {*0.1.29}
  if (b<=0) or (a.used=0) then exit;

  with a do begin
    if longint(b)+used > MaxDigits then begin   {V1.20.04: changed <= to <}
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXMaxDigits.Create('mp_lshd');
        {$else}
          RunError(MP_RTE_OTHER);
        {$endif}
      {$else}
        set_mp_error(MP_MAXDIGITS);
        exit;
      {$endif}
    end;
    {grow to fit the new digits}
    if alloc < used+b then begin
      mp_grow(a, used+b);
      if mp_error<>MP_OKAY then exit;
    end;
    {shift used digits, used>0!}
    move(pdigits^[0], pdigits^[b], used*sizeof(mp_digit));
    {clear low b digits, b>0!}
    fillchar(pdigits^, b*sizeof(mp_digit), 0);
    inc(used, b);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_makeodd(const a: mp_int; var b: mp_int; var s: longint);
  {-Return b,s with a = 2^s*b if a<>0, b=0,s=-1 otherwise}
begin
  {count the number of least significant bits which are zero}
  s := mp_cnt_lsb(a);
  if mp_error=MP_OKAY then begin
    {now divide a by 2^s}
    if (s=0) and (a.used=0) then begin
      mp_zero(b);
      s := -1;
    end
    else mp_shr(a, s, b);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_mod(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a mod b, 0 <= c < b}
var
  t: mp_int;
  s: longint;
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_mod');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mp_is_longint(b,s) then begin
    {b fits into a longint and we do not need a quotient here, so use mp_mod_int}
    {$ifdef MP_32BIT}
      mp_mod_int(a,s,s);
      mp_set_int(c,s);
      exit;
    {$else}
      {For MP_16BIT s_mp_div is faster for small b, check if b is not small}
      if {$ifndef BIT16}(s and $7FFF0000 <> 0){$else}b.used>1{$endif} then begin
        mp_mod_int(a,s,s);
        mp_set_int(c,s);
        exit;
      end;
    {$endif}
  end;

  mp_init(t);
  if mp_error=MP_OKAY then begin
    {call Burnikel/Ziegler or Knuth's basecase algorithm D}
    s_mp_divrem(a, b, nil, @t);
    if mp_error=MP_OKAY then begin
      {don't adjust sign if result=0}
      if (t.used>0) and (t.sign<>b.sign) then mp_add(b, t, c)
      else mp_exch(t, c);
    end;
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_mod_2k(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a mod 2^b, 0 <= c < 2^b}
var
  t: mp_int;
begin
  {checks done in s_mp_mod_2k}
  s_mp_mod_2k(a,b,c);
  if c.sign=MP_NEG then begin
    mp_init(t);
    mp_2expt(t,b);
    mp_add(c,t,c);
    mp_clear(t);
  end
end;


{---------------------------------------------------------------------------}
procedure s_mp_mod_int(const a: mp_int; b: longint; var c: longint);
  {-Calculate r = a mod b for a single longint b, no init check}
var
  w,q: longint;
  ix: TNInt;
  d: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  q := abs(b);
  {cannot divide by zero}
  if q=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('s_mp_mod_int: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  {check for trivial cases}
  if (q=1) or (a.used=0) then begin
    c := 0;
    exit;
  end;

  if q and (q-1) = 0 then begin
    {b is power of 2}
    w := abs(mp_get_int(a)) and (q-1);
  end
  else begin
    if mp_is_longint(a,w) then w := abs(w) mod q
    else begin
      w := 0;
      {note: a.used>0}
      for ix:=a.used-1 downto 0 do begin
        d := a.pdigits^[ix];
        {$ifdef MPC_PurePascal}
          w := ((int64(w) shl DIGIT_BIT) or d) mod q;
        {$else}
          {$ifdef MP_32BIT}
              asm
                mov    eax,[w]
                sub    edx,edx
                shld   edx,eax,DIGIT_BIT
                shl    eax,DIGIT_BIT
                or     eax,[d]
                div    [q]
                mov    [w],edx
              end;
          {$else}
            {$ifdef BIT16}
              asm
                db $66; mov  ax,word ptr [w]
                db $66; mov  dx,ax
                db $66; shl  ax,DIGIT_BIT
                db $66; shr  dx,32-DIGIT_BIT
                        or   ax,[d]
                db $66; div  word ptr [q]
                db $66; mov  word ptr [w],dx
              end;
            {$else}
              asm
                movzx  ecx,[d]
                mov    eax,[w]
                mov    edx,eax
                shl    eax,DIGIT_BIT
                shr    edx,32-DIGIT_BIT
                or     eax,ecx
                div    [q]
                mov    [w],edx
              end;
            {$endif}
          {$endif}
        {$endif}
      end;
    end;
  end;
  {adjust sign for non-zero result}
  if w<>0 then begin
    {w is positive. First if a<0 make remainder negative}
    if a.sign=MP_NEG then w := -w;
    {if remainder and modulus have different sign, add modulus}
    if w xor b < 0 then w := w+b;
  end;
  c := w;
end;


{---------------------------------------------------------------------------}
procedure mp_mod_d(const a: mp_int; b: mp_digit; var c: mp_digit);
  {-Calculate c = a mod b, 0 <= c < b (digit version)}
{$ifndef MP_32BIT}
var
  w: mp_word;
  ix: TNInt;
{$endif}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a)  then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_mod_d');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if b>MP_DIGIT_MAX then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_mod_d: b>MP_DIGIT_MAX');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}

  {$ifdef MP_32BIT}
    {use s_mp_mod_int because the Delphi 'optimizer' produces very}
    {slow code for the standard version, FPC is much better here. }
    s_mp_mod_int(a, longint(b), longint(c));
  {$else}
    if b=0 then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_mod_d: b=0');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;

    if a.used=0 then c := 0
    else if (b and pred(b)) = 0 then begin
      {b=power of two (incl b=1), simple and no need to get shift count}
      c := a.pdigits^[0] and pred(b);
    end
    else begin
      w := 0;
      {a.used>0}
      for ix:=a.used-1 downto 0 do begin
        {$ifdef BIT16}
          LeftShiftAdd(w, a.pdigits^[ix]);
        {$else}
          w := (w shl DIGIT_BIT) or a.pdigits^[ix];
        {$endif}
        if w>=b then w := w mod b;
      end;
      c := mp_digit(w);
    end;

    {Here 0 <= c < b. Return b-c if a < 0}
    if (c<>0) and (a.sign=MP_NEG) then c := b-c;
  {$endif}

end;


{---------------------------------------------------------------------------}
procedure mp_mod_int(const a: mp_int; b: longint; var c: longint);
  {-Calculate r = a mod b for a single longint b}
begin
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_mod_int');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mp_mod_int(a,b,c);
end;


{---------------------------------------------------------------------------}
function s_mp_mod_is0(const a,b: mp_int): boolean;
  {-Test if a mod b = 0, ie if a is a multiple of b}
var
  t: mp_int;
begin
  s_mp_mod_is0 := false;
  if mp_error<>MP_OKAY then exit;
  mp_init(t);
  if mp_error=MP_OKAY then begin
    mp_mod(a,b,t);
    s_mp_mod_is0 := mp_is0(t);
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_mod_w(const a: mp_int; b: word; var r: word);
  {-Calculate r = |a| mod b for a single word b}
var
  w: mp_word;
  ix: TNInt;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_mod_aw');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {cannot divide by zero}
  if b=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_mod_aw: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {check for trivial cases}
  if (b=1) or (a.used=0) then begin
    r := 0;
    exit;
  end;

  w := 0;
  {note: a.used>0}
  for ix:=a.used-1 downto 0 do begin
    {$ifdef BIT16}
      w := fLeftShiftAdd(w, a.pdigits^[ix]) mod b;
    {$else}
      w := (w shl DIGIT_BIT) or a.pdigits^[ix];
      w := w mod b;
    {$endif}
  end;
  r := word(w);
end;


{---------------------------------------------------------------------------}
procedure mp_mod_w(const a: mp_int; b: word; var r: word);
  {-Calculate r = a mod b for a single word b}
begin
  s_mp_mod_w(a,b,r);
  if (mp_error=MP_OKAY) and (a.sign=MP_NEG) and (r>0) and (r<=b) then r := b-r;
end;


{---------------------------------------------------------------------------}
procedure mp_montgomery_calcnorm(var R: mp_int; const m: mp_int);
  {-Calculate R = B^n mod m, n=number of digits in m, B=2^DIGIT_BIT}
var
  x, bits: integer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(R) or mp_not_init(m) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_montgomery_calcnorm');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {The method is slightly modified to shift B unconditionally up to just under}
  {the leading bit of m.  This saves a lot of multiple precision shifting}

  {how many bits of last digit does m use}
  bits := mp_bitsize(m) mod DIGIT_BIT;
  if m.used>1 then mp_2expt(R, (m.used-1)*DIGIT_BIT + bits - 1)
  else begin
    mp_set(R, 1);
    bits := 1;
  end;
  if mp_error<>MP_OKAY then exit;
  {now compute R * B mod m}
  for x:= bits to DIGIT_BIT do begin
    mp_mul_2(R, R);
    if mp_cmp_mag(R, m) <> MP_LT then s_mp_sub(R, m, R);
    if mp_error<>MP_OKAY then exit;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_montgomery_reduce(var x: mp_int; const n: mp_int; rho: mp_digit);
  {-Calculate x = xR^-1 (mod n) via Montgomery reduction}
var
  ix, iy: TNInt;
  tmpn, tmpx: pmp_digit;
  digs: longint;
{$ifndef BASM16}
  mu,t,u: mp_digit;
  r: mp_word;
{$endif}
begin
  if mp_error<>MP_OKAY then exit;

  {$ifdef MPC_ArgCheck}
    if mp_not_init(x) or mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_montgomery_reduce');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if n.used<1 then exit;
  digs := n.used * 2 + 1;
  {grow the input as required}
  if x.alloc < digs then mp_grow(x, digs);
  if mp_error<>MP_OKAY then exit;
  x.used := digs;
  for ix:=0 to n.used-1 do begin
    {alias for digits of the modulus}
    tmpn := pmp_digit(n.pdigits);

    {alias for the digits of x [the input]}
    tmpx := @x.pdigits^[ix];

    {$ifdef BASM16}
      iy := n.used;
      asm
           push  ds
           lds   si,[tmpn]
           les   di,[tmpx]
           mov   ax,[rho]
           mul   word ptr es:[di]
           and   ax,MP_MASK
           mov   bx,ax                       {bx=mu=tmpx^*rho and MP_MASK}
           sub   cx,cx
           cld
      @@1: lodsw                             {ax=mp_word(tmpn^)}
           mul   bx                          {mp_word(mu)*mp_word(tmpn^)}
           add   ax,cx
           adc   dx,0                        {mp_word(mu)*mp_word(tmpn^) + mp_word(u)}
           add   ax,es:[di]
           adc   dx,0                        {r := mp_word(mu)*mp_word(tmpn^) + mp_word(u) + mp_word(tmpx^)}
           mov   cx,ax
           and   ax,MP_MASK
           stosw                             {tmpx^ := mp_digit(r and MP_MASK)}
           db    $0F,$AC,$D1,DIGIT_BIT       {shrd cx,dx,DIGIT_BIT}
           dec   [iy]
           jnz   @@1
           jcxz  @@3
      @@2: add   cx,es:[di]
           mov   ax,MP_MASK
           and   ax,cx
           stosw
           shr   cx,DIGIT_BIT
           jnz   @@2
      @@3: pop   ds
      end;
    {$else}
      {mu = ai * rho mod B}
      {The value of rho must be precalculated via montgomery_setup() such that
       it equals -1/n mod B this allows the following inner loop to reduce the
       input one digit at a time}
      mu := mp_word(x.pdigits^[ix]) * mp_word(rho) and MP_MASK;
      {set the carry to zero}
      u := 0;
      {Multiply and add in place}
      for iy:=0 to n.used-1 do begin
        {compute product and sum}
        r := mp_word(mu) * mp_word(tmpn^) + mp_word(u) + mp_word(tmpx^);
        {get carry}
        u := mp_digit(r shr DIGIT_BIT);
        {fix digit}
        tmpx^ := mp_digit(r and MP_MASK);
        inc(tmpn);
        inc(tmpx);
      end;
      {At this point the ix'th digit of x should be zero}
      {propagate carries upwards as required}
      while u<>0 do begin
        t := tmpx^ + u;
        u := t shr DIGIT_BIT;
        tmpx^ := t and MP_MASK;
        inc(tmpx);
      end;
    {$endif}

  end;

  {at this point the n.used'th least significant digits of x are all zero
   which means we can shift x to the right by n.used digits and the
   residue is unchanged.}

  {x = x/B^n.used}
  mp_clamp(x);
  mp_rshd(x, n.used);

  {if x >= n then x = x - n}
  if mp_cmp_mag(x, n) <> MP_LT then s_mp_sub(x, n, x);
end;


{---------------------------------------------------------------------------}
procedure mp_montgomery_setup(const n: mp_int; var rho: mp_digit);
  {-Calculate rho = -1/n mod B for Montgomery reduction, B=2^DIGIT_BIT}
var
  x, b: mp_digit;
begin
  {Fast inversion mod 2^k, based on the fact that
    XA = 1 (mod 2^n)  =>  (X(2-XA)) A = 1 (mod 2^2n)
                      =>  2*X*A - X*X*A*A = 1
                      =>  2*(1) - (1)     = 1 }
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_montgomery_setup');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  b := n.pdigits^[0];
  if b and 1 = 0 then begin
    {n must be odd}
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_montgomery_setup: n is even');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {temporarily turn range checks off}
  {$R-}
  x := (((b + 2) and 4) shl 1) + b;  {here x*a==1 mod 2^4}
  x := x * (2 - b * x);              {here x*a==1 mod 2^8}
  {FPC V2 gives "unreachable code warnings" for < 32bit, that's OK}
  if sizeof(mp_digit) > 1 then begin
    x := x * (2 - b * x);            {here x*a==1 mod 2^16}
    {$ifdef MP_32BIT}
      if sizeof(mp_digit) > 2 then begin
        x := x * (2 - b * x);          {here x*a==1 mod 2^32}
      end;
    {$endif}
  end;
  {$ifdef RangeChecks_on} {$R+} {$endif}
  {rho = -1/n mod B }
  rho := ((mp_word(1) shl mp_word(MP_DIGIT_BIT)) - x) and MP_MASK;
end;


{---------------------------------------------------------------------------}
procedure mp_mul(const a,b: mp_int; var c: mp_int);
  {-High level multiplication, c = a*b}
var
  s,t,ai,bb: mp_int;
  pa,pb: pmp_int;
  n,minu,maxu,i0,i1,i: longint;
  sig: word;
begin
  if mp_error<>MP_OKAY then exit;
  if @a=@b then begin
    mp_sqr(a,c);
    exit;
  end;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_mul');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {get result sign, calculated product sign is not used}
  if a.sign=b.sign then sig:=MP_ZPOS else sig:=MP_NEG;

  {get size/addr of arguments. pa points to larger, pb to smaller factor}
  if a.used >= b.used then begin
    pa := @a;
    pb := @b;
    maxu := a.used;
    minu := b.used;
  end
  else begin
    pb := @a;
    pa := @b;
    minu := a.used;
    maxu := b.used;
  end;

  if (minu < mp_mul_cutoff) or (maxu < 2*minu) then begin
    {small or approximately balanced factors}
    if minu=0 then begin
      mp_zero(c);
      exit;
    end;
    if minu >= mp_t3m_cutoff then s_mp_toom3_mul(a, b, c)
    else if minu >= mp_mul_cutoff then s_mp_karatsuba_mul(a, b, c)
    else begin
      {optimization for single digit factors, note: a,b <>0}
      if a.used=1 then begin
        {sign adjust is done below}
        mp_mul_d(b,a.pdigits^[0],c);
      end
      else if b.used=1 then begin
        {sign adjust is done below}
        mp_mul_d(a,b.pdigits^[0],c);
      end
      else s_mp_mul_digs(a, b, c, a.used + b.used + 1);
    end;
  end
  else begin
    {Large and unbalanced factors: partition larger factor 'a' in chunks of}
    {size equal to smaller factor a = a[n]*X^n .. a[1]*X + a[0] and compute}
    {a*b as (..((a[n]*b*X + a[n-1])*X .. )*X + a[0], X=2^(b.used*DIGIT_BIT)}
    mp_init2(s,t);
    n := maxu div minu;
    i0 := n*minu;
    i1 := i0+pred(minu);
    {make local positive copy of b}
    bb := pb^;
    bb.sign := MP_ZPOS;
    {make fake chunk a[n]}
    s_mp_fakeinit(pa^,i0,i1,ai);
    mp_mul(ai,bb,s);
    for i:=pred(n) downto 0 do begin
      dec(i0,minu);
      dec(i1,minu);
      {ai represents a[i]}
      s_mp_fakeinit(pa^,i0,i1,ai);
      mp_mul(ai,bb,t);
      {multiply by X}
      mp_lshd(s,minu);
      mp_add(s,t,s);
    end;
    mp_exch(s,c);
    mp_clear2(s,t);
  end;
  {adjust result sign}
  if c.used=0 then c.sign := MP_ZPOS else c.sign := sig;
end;


{---------------------------------------------------------------------------}
procedure mp_mul_int(const a: mp_int; b: longint; var c: mp_int);
  {-Multiply by a 32 bit integer}
var
  ba : longint;
  tmp: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  ba := abs(b);

  {$ifdef MP_32BIT}
    if ba and longint(MP_INV_MASK) = 0 then begin
      {arg check in mp_mul_d}
      mp_mul_d(a,mp_digit(ba),c);
      if b<0 then s_mp_chs(c);
      exit;
    end;
  {$else}
    if ba<=$FFFF then begin
      {arg check in mp_mul_w}
      mp_mul_w(a,word(ba),c);
      if b<0 then s_mp_chs(c);
      exit;
    end;
  {$endif}

  {arg check in mp_mul}
  mp_init_set_int(tmp, b);
  if mp_error=MP_OKAY then begin
    mp_mul(a,tmp,c);
    mp_clear(tmp);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_mul_2(const a: mp_int; var b: mp_int);
  {-Multiply by 2, b = 2*a}
var
  oldused: longint;
  i: TNInt;
  r,rr: mp_digit;
  pa,pb: pmp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_mul_2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if a=0 then just copy}
  if a.used=0 then begin
    mp_copy(a, b);
    exit;
  end;

  {grow dest}
  if b.alloc < a.used+1 then begin
    mp_grow(b, a.used+1);
    if mp_error<>MP_OKAY then exit;
  end;

  oldused := b.used;
  b.used := a.used;
  pa := pmp_digit(a.pdigits);
  pb := pmp_digit(b.pdigits);
  {carry}
  r  := 0;
  for i:=0 to a.used-1 do begin
    {get what will be the *next* carry bit from the MSB of the current digit}
    rr := pa^ shr (DIGIT_BIT - 1);
    {now shift up this digit, add in the carry [from the previous]}
    pb^ := ((pa^ shl 1) or r) and MP_MASK;
    {copy carry that would be from the source digit into the next iteration}
    r := rr;
    inc(pa);
    inc(pb);
  end;

  {new leading digit?}
  if r<>0 then begin
    {add a MSB which is always 1 at this point}
    pb^ := 1;
    inc(b.used);
    inc(pb);
  end;

  {now zero any excess digits on the destination that we did not write to}
  for i:=b.used to oldused-1 do begin
    pb^ := 0;
    inc(pb);
  end;
  b.sign := a.sign;
end;


{---------------------------------------------------------------------------}
procedure mp_mul_2k(const a: mp_int; b: longint; var c: mp_int);
  {-Shift left a, c = a*2^b; c=a if b<=0 [synonym for mp_shl]}
begin
  mp_shl(a,b,c);
end;


{---------------------------------------------------------------------------}
procedure mp_mul_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Multiply by a digit}
var
  pa,pc: pmp_digit;
{$ifndef BASM16}
  u,r,bw: mp_word;
{$endif}
  n: TNInt;
  olduse: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_mul_d');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if b>MP_DIGIT_MAX then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_mul_d: b>MP_DIGIT_MAX');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}

  {if a=0 or b=0, set c=0 and exit}
  if (a.used=0) or (b=0) then begin
    mp_zero(c);
    exit;
  end;

  {trivial case b=1, just copy a}
  if b=1 then begin
    mp_copy(a,c);
    exit;
  end;

  {make sure c is big enough to hold a*b}
  if c.alloc < a.used + 1 then begin
    mp_grow(c, a.used+1);
    if mp_error<>MP_OKAY then exit;
  end;

  {get the original destinations used count}
  olduse := c.used;
  {set the sign}
  c.sign := a.sign;
  {setup pointers}
  pa := pmp_digit(a.pdigits);
  pc := pmp_digit(c.pdigits);

{$ifdef BASM16}
  n := a.used;
  asm
             push ds
             lds  si,[pa]
             les  di,[pc]
    db $66;  sub  cx,cx             {ecx=u=0}
    db $66;  sub  bx,bx             {ebx=bw=b}
             mov  bx,[b]

    @@1:
    db $66;  db   $0F,$B7,$04       {movzx  eax,word ptr [si]}
             add  si,2
    db $66;  mul  bx
    db $66;  add  cx,ax             {r = pa^ * bw + u}
             mov  ax,cx
             and  ax,MP_MASK
             mov  es:[di],ax        {pc^ := r and MP_MASK}
             add  di,2
    db $66;  shr  cx,DIGIT_BIT      {u := r shr DIGIT_BIT}
             dec  [n]
             jnz  @@1
             mov  es:[di],cx        {pc^ := mp_digit(u)}
             pop  ds
  end;
{$else}
  {zero carry}
  u  := 0;
  bw := b;

  for n:=0 to a.used-1 do begin
    {compute product and carry sum for this term}
    r := pa^ * bw + u;
    inc(pa);
    {mask off higher bits to get a single digit}
    pc^ := r and MP_MASK;
    inc(pc);
    {send carry into next iteration}
    u := r shr DIGIT_BIT;
  end;

  {store final carry [if any]}
  pc^ := mp_digit(u);
{$endif}

  {set used count}
  c.used := a.used + 1;

  {now zero digits above the top}
  if c.used<olduse then begin
    pc := @c.pdigits^[c.used];
    for n:=c.used to olduse-1 do begin
      pc^ := 0;
      inc(pc);
    end;
  end;

  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
procedure mp_mul_w(const a: mp_int; b: word; var c: mp_int);
  {-Multiply by a word}
var
  pa,pc: pmp_digit;
  u,r,bw: mp_word;
  n: TNInt;
  olduse: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_mul_w');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if a=0 or b=0, set c=0 and exit}
  if (a.used=0) or (b=0) then begin
    mp_zero(c);
    exit;
  end;

  {trivial case b=1, just copy a}
  if b=1 then begin
    mp_copy(a,c);
    exit;
  end;

  {a word needs up to 2 mp_digits, an mp_word can hold an mp_digit*word}

  {make sure c is big enough to hold a*b}
  if c.alloc < a.used+2 then begin
    mp_grow(c, a.used+2);
    if mp_error<>MP_OKAY then exit;
  end;

  {get the original destinations used count}
  olduse := c.used;
  {set the sign}
  c.sign := a.sign;
  {setup pointers}
  pa := pmp_digit(a.pdigits);
  pc := pmp_digit(c.pdigits);

  {zero carry}
  u  := 0;
  bw := b;

  for n:=0 to a.used-1 do begin
    {compute product and carry sum for this term}
    r := pa^ * bw + u;
    inc(pa);
    {mask off higher bits to get a single digit}
    pc^ := r and MP_MASK;
    inc(pc);
    {send carry into next iteration}
    {$ifdef BASM16}
      asm
        db  $66;  mov ax, word ptr [r]
        db  $66;  shr ax, DIGIT_BIT
        db  $66;  mov word ptr [u], ax
      end;
    {$else}
      u := r shr DIGIT_BIT;
    {$endif}
  end;

  {store final carries [if any]}
  pc^ := mp_digit(u and MP_MASK);
  inc(pc);
  pc^ := mp_digit((u shr DIGIT_BIT));

  {set used count}
  c.used := a.used + 2;

  {now zero digits above the top}
  if c.used<olduse then begin
    pc := @c.pdigits^[c.used];
    for n:=c.used to olduse-1 do begin
      pc^ := 0;
      inc(pc);
    end;
  end;

  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
function mp_not_init(const a: mp_int): boolean;
  {-Sanity check if a is initialized, does not catch all cases!}
begin
  mp_not_init := (a.magic<>MP_MAGIC) or (a.pdigits=nil) or (a.used > a.alloc);
end;


{---------------------------------------------------------------------------}
function mp_not_init_multi(const a: array of mp_int): boolean;
  {-Sanity check if all elements of a are initialized, does not catch all cases!}
var
  i: integer;
begin
  mp_not_init_multi := true;
  if mp_error<>MP_OKAY then exit;
  for i:=0 to high(a) do begin
    if mp_not_init(a[i]) then exit;
  end;
  mp_not_init_multi := false;
end;



(*********************************************************************

Transformation of the standard Halley formula for nth root calculation.
From [http://en.wikipedia.org/wiki/Halley's_method] to the form from
[http://yacas.sourceforge.net/Algochapter4.html#c4s11]

x~ = x - 2*f(x)*f'(x)/[2*f'(x)^2 - f''(x)*f(x)]

x~ = x - f(x)/[f'(x) - f''(x)/f'(x)*f(x)/2]

f(x)=x^n-a, f'(x)=n*x^(n-1), f''(x)=n*(n-1)*x^(n-2), f''(x)/f'(x)=(n-1)/x

x~ = x - f(x)/[n*x^(n-1) - (n-1)/2x*(x^n-a)]

   = x - 2x*(x^n-a)/[2n*x^n - (n-1)*(x^n-a)]

   = x - 2x*(x^n-a)/(2n*x^n - n*x^n + n*a + x^n - a)

   = x - 2x*(x^n-a)/(n*x^n + n*a + x^n -a)

   = x - 2x*(x^n-a)/[n*(x^n+a) -(a-x^n)]

   = x*(1 - 2*(x^n-a)/[n*(x^n+a) - (a-x^n)])

   = x*([n*(x^n+a) - (a-x^n) - 2*(x^n-a)]/[n*(x^n+a) - (a-x^n)])

   = x*([n*x^n + n*a - a + x^n - 2x^n + 2a]/[n*(x^n+a) - (a-x^n)])

   = x*([n*x^n + n*a - x^n + a]/[n*(x^n+a) - (a-x^n)])

         n*(a+x^n) + (a-x^n)
x~ = x * -------------------
         n*(a+x^n) - (a-x^n)

*********************************************************************)


{---------------------------------------------------------------------------}
function s_mp_n_root2(const a: mp_int; n: longint; var b: mp_int; pr: pmp_int): boolean;
  {-Calculate the n'th root of a, n>=2, pr^ = a-b^n, a must be >=0 }
  { if n is even; return true, if b is an exact root; no init check}
var
  level: integer; {recursion level}
  bsa: longint;   {bitsize(a)}
  siga: word;     {sign of a}
  UseNewton: boolean;
  nsfac: double;
const
  powtab: array[1..31] of mp_digit = {Pari/GP: vector(31,n,floor(128*2^(n/32.0)))}
            (130, 133, 136, 139, 142, 145, 148, 152, 155, 158, 162,
             165, 169, 173, 177, 181, 184, 189, 193, 197, 201, 206,
             210, 215, 219, 224, 229, 234, 239, 245, 250);

  {---------------------------------------------------}
  function iroot(var z: mp_int; d: longint): boolean;
    {-Internal root proc, z := z^(1/d), z>0, d>1; return true if known exact}
  const
    ProgCheck = 10;
  var
    w,x,y,t: mp_int;
    iter, ichk, r, bsz: longint;
    done: boolean;
    s: double;
    {meaningful alias names for bisection steps}
    xl: mp_int absolute w;  {left  }
    xc: mp_int absolute t;  {center}
    xr: mp_int absolute x;  {right }
  begin
    iroot := false;
    {Create initial approximation}
    mp_init(x);
    if mp_error<>MP_OKAY then exit;

    bsz := mp_bitsize(z);

    r := bsz div trunc(1.5*d*DIGIT_BIT);
    if r>0 then begin
      {Initial approximation x0 based on high order bits, cf. Keith Matthews, }
      {Computing m-th roots, The College Mathematics Journal 19 (1988) 174-176}
      {http://www.numbertheory.org/keith/mthroot.html}
      mp_rshd2(z,x,r*d);
      inc(level);
      if iroot(x,d) then ;
      dec(level);
      mp_inc(x);
      mp_lshd(x,r);
    end
    else begin
      {Initial approximation for base of recursion or small z,d values.}
      {Uses integer part and five bits of fractional part of ln(root). }
      s := s_mp_log2(z)/d;
      mp_2expt(x,trunc(s));
      iter := 1+trunc(frac(s)*32);
      if iter=32 then mp_shl1(x)
      else begin
        mp_mul_d(x,powtab[iter],x);
        mp_shr(x,7,x);
      end;
    end;

    {Create remaining local mp_ints}
    mp_init3(t,w,y);
    if mp_error<>MP_OKAY then begin
      mp_clear(x);
      exit;
    end;

    iter := 0;
    ichk := ProgCheck;  {iter value for progress check}
    done := false;

    while (not done) and (mp_error=MP_OKAY) do begin
      mp_exch(x,y);             {y = x0 (exch is faster than copy)}
      mp_expt_int(y,d-1,t);     {t = x0^(d-1)}
      mp_mul(y,t,x);            {x = x0^d}
      mp_sub(z,x,w);            {w = z-x0^d}
      if mp_is0(w) then begin
        {exact root accidently found}
        if (level=0) and (pr<>nil) then mp_zero(pr^);
        mp_exch(z,y);
        mp_clear4(t,w,x,y);
        iroot := true;
        exit;
      end;
      if (iter=0) and (w.sign=MP_ZPOS) then begin
        {initial approximation (now in y) was too low, increment and retry}
        if r>0 then begin
          mp_2expt(t,r*DIGIT_BIT-1);
          mp_add(y,t,x);
        end
        else begin
          mp_mul_d(y,130,x);
          mp_shr(x,7,x);
          if mp_is_le(x,y) then begin
            mp_mul_d(y,5,x);
            mp_shr(x,2,x);
            if mp_is_le(x,y) then begin
              {give up, use simple working approximation}
              mp_2expt(x, (bsz + pred(d)) div d);
            end;
          end;
        end;
        continue;
      end;
      inc(iter);
      if iter=ichk then begin
        {Check startup convergence of Halley steps, if too slow}
        {do some bisection steps and continue with Halley}
        if mp_bitsize(w)+ProgCheck > bsz then begin
          inc(ichk,bitsize32(d));
          mp_2expt(xr, (bsz + pred(d)) div d);
          mp_2expt(xl, (bsz - {pred(d)}-1) div d);
          while (mp_error=MP_OKAY) and mp_is_lt(xl,xr) do begin
            {Compute center xc = (xl+xr)/2, and y =xc^d}
            mp_add(xl,xr,xc);
            mp_shr1(xc);
            mp_expt_int(xc, d, y);
            {Exchange xc with xl or xr depending on sign}
            {of y-z or exit immediately if xc^d=z}
            case mp_cmp_mag(y,z) of
              -1: s_mp_add_d(xc,1,xl);
               0: begin
                    {exact root found}
                    if (pr<>nil) and (level=0) then mp_zero(pr^);
                    mp_exch(z,xc);
                    mp_clear4(t,w,x,y);
                    iroot := true;
                    exit;
                  end;
               1: mp_exch(xr,xc);
            end;
            inc(iter);
            if iter>ichk then break;
          end;
          {Here x=xr will be the next x0 for Halley or Newton, hopefully }
          {a better approximation. Note: We must use xr although xl might}
          {be the better approximation because x0^d must be >= z.}
          ichk := -1; {no more checks}
          continue;
        end;
      end;

      if UseNewton then begin
        {Use Newton if Halley is supposed to produce overflows}
        {x1 = ((d-1)*x0^d + z)/(d*x0^(d-1) = ((d-1)*x+z)/(d*t)}
        {   = (d*x + z-w)/(d*t) = (d*x + w)/(d*t) = (x+w/d)/t }
        mp_div_int(w,d,@w,r);
        mp_add(x,w,x);
        mp_div(x,t,x);
      end
      else begin
        { Halley's iteration in main loop: }
        {           d*(z+x0^d) + (z-x0^d)  }
        { x1 = x0 * ---------------------  }
        {           d*(z+x0^d) - (z-x0^d)  }
        mp_add(z,x,x);      {x = z+x0^d}
        mp_mul_int(x,d,x);  {x = d*(z+x0^d)}
        mp_add(x,w,t);      {t = d*(z+x0^d) + (z-x0^d)}
        mp_sub(x,w,x);      {x = d*(z+x0^d) - (z-x0^d)}
        mp_mul(y,t,t);
        mp_div(t,x,x);      {x1= x0*t/x}
      end;
      done := mp_cmp_mag(x,y)>=0; {done if x1 >= x0}
    end;

    {Root is in y, but z - y^d <> 0 otherwise we would not be here.}
    {Exception is mp_error<>0, but then result should be false too!}
    if (level=0) and (pr<>nil) then begin
      if UseNewton then begin
        {z-y^d not directly available, only t = y^(d-1): pr^ = z-y*t}
        mp_mul(y,t,t);
        mp_sub(z,t,pr^);
      end
      else begin
        {Halley: w = z-y^d}
        mp_copy(w,pr^);
      end;
    end;

    mp_exch(z,y);
    mp_clear4(t,w,x,y);
  end;

begin
  s_mp_n_root2 := false;
  if mp_error<>MP_OKAY then exit;

  {save sign of a}
  siga := a.sign;
  if (n and 1 = 0) and (siga=MP_NEG) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('s_mp_n_root2: a<0 and n even');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  {calculate root of b=abs(a)}
  mp_abs(a,b);
  if mp_is0(b) or mp_is1(b) then begin
    s_mp_n_root2 := true;
    if pr<>nil then mp_zero(pr^);
  end
  else begin
    bsa := mp_bitsize(a);
    if n >= bsa then begin
      {a < 2^n, therefore b := a^(1/n) < 2,  ie b = 1, pr^ = a-b^n = a-1}
      if pr<>nil then mp_sub_d(a,1,pr^);
      mp_set1(b);
    end
    else begin
      {$ifdef MSDOS}
        nsfac := 0.01*MaxAvail/a.used;
        {Next statement fixes an obscure BP7 FP RTE 207 which occurs}
        {if output is redirected to a file!? Not with DPMI!}
        if nsfac > 0.95 then nsfac := 0.95;
      {$else}
       nsfac := 0.95;
      {$endif}
      UseNewton := bsa >= trunc(MP_MAXBIT*nsfac*n/(n+1));
      level := 0;
      s_mp_n_root2 := iroot(b,n);
    end;
  end;
  {copy original sign}
  b.sign := siga;
  if (pr<>nil) and (pr^.used<>0) then pr^.sign := siga;
end;


{---------------------------------------------------------------------------}
procedure mp_n_root2(const a: mp_int; n: longint; var b: mp_int; pr: pmp_int);
  {-Calculate the n'th root of a, pr^=a-b^n, a must be >=0 if n is even; b=0,pr^=0 if n<1}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or ((pr<>nil) and mp_not_init(pr^)) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_n_root2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if (n<1) or mp_iszero(a) then begin
    mp_zero(b);
    if pr<>nil then mp_zero(pr^);
  end
  else if n=1 then begin
    mp_copy(a,b);
    if pr<>nil then mp_zero(pr^);
  end
  else if n=2 then begin
    if pr=nil then mp_sqrt(a,b)
    else begin
      mp_sqrtrem(a,b,pr^);
    end;
  end
  else s_mp_n_root2(a,n,b,pr);
end;


{---------------------------------------------------------------------------}
procedure mp_n_root(const a: mp_int; n: longint; var b: mp_int);
  {-Calculate the n'th root of a, a must be >=0 if n is even; b=0 if n<1}
begin
  mp_n_root2(a,n,b,nil);
end;


{---------------------------------------------------------------------------}
procedure mp_output_decimal(const a: mp_int);
  {-Write decimal representation to output}
begin
  mp_write_decimal(output, a);
end;


{---------------------------------------------------------------------------}
procedure mp_output_radix(const a: mp_int; radix: word);
  {-Write radix representation to output}
begin
  mp_write_radix(output, a, radix);
end;


{---------------------------------------------------------------------------}
procedure mp_or(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a or b}
var
  pa,pb,pc,px: pmp_digit;
  olduse, min, max: longint;
  i: TNInt;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_or');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if a.used > b.used then max := a.used else max := b.used;

  {easy outs}
  if max=0 then begin
    mp_zero(c);
    exit;
  end;
  if a.used=0 then begin
    mp_copy(b,c);
    exit;
  end;
  if b.used=0 then begin
    mp_copy(a,c);
    exit;
  end;

  {must grow before assigning pointer otherwise realloc can}
  {change the digits memory, ie px points to the old digits}
  if c.alloc < max then begin
    mp_grow(c, max);
    if mp_error<>MP_OKAY then exit;
  end;

  if a.used > b.used then begin
    min := b.used;
    px  := @a.pdigits^[min];
  end
  else begin
    min := a.used;
    px  := @b.pdigits^[min];
  end;

  pa  := pmp_digit(a.pdigits);
  pb  := pmp_digit(b.pdigits);
  pc  := pmp_digit(c.pdigits);

  {get old used digit count and set new one}
  olduse := c.used;
  c.used := max;
  c.sign := a.sign or b.sign;

  {process digits 0 .. min-1}
  for i:=0 to pred(min) do begin
    pc^ := pa^ or pb^;
    inc(pa);
    inc(pb);
    inc(pc);
  end;

  {copy digits min .. max-1}
  if min<max then begin
    for i:=min to pred(max) do begin
      pc^ := px^;
      inc(px);
      inc(pc);
    end;
  end;

  {clear c above max if necessary}
  if olduse>max then begin
    for i:=max to pred(olduse) do begin
      pc^ := 0;
      inc(pc);
    end;
  end;
  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
function mp_popcount(const a: mp_int): longint;
  {-Get population count = number of 1-bits in a}
var
  i: TNInt;
  bw: longint;
  pa: pmp_digit;
begin
  mp_popcount := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_popcount');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if a.used>0 then begin
    bw := 0;
    pa := pmp_digit(a.pdigits);
    for i:=0 to a.used-1 do begin
      {$ifdef MP_32BIT}
        inc(bw,popcount32(pa^ and MP_DIGIT_MAX));
      {$else}
        inc(bw,popcount16(pa^ and MP_DIGIT_MAX));
      {$endif}
      inc(pa);
    end;
    mp_popcount := bw;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_prod_int(var a: mp_int; const b: array of longint; n: longint);
  {-Calculate a = product of first n elements of longint array b}
var
  thresh: longint;            {threshold for recursive calls}
  level:  integer;            {recursion level}
  q: array[0..32] of mp_int;  {stack for recursion}

  procedure rprod(var p: mp_int; L,H: longint);
    {-Calculate p=b[L]*b[L+1]*...b[H]}
  var
    m: longint;
  begin
    if mp_error<>MP_OKAY then exit;
    {increase recursion level}
    inc(level);
    if H-L<=thresh then begin
      if L>H then mp_set1(p)
      else begin
        mp_set_int(p,b[L]);
        {skip loop if p is zero}
        while (L<H) and (p.used>0) and (mp_error=MP_OKAY) do begin
          inc(L);
          s_mp_mul_int(p,b[L],p,q[0]);
        end;
      end;
    end
    else begin
      {recursive calls}
      m := (L + H) div 2;
      rprod(p,L,m);
      {calculate second half only if p<>0}
      if (p.used>0) and (mp_error=MP_OKAY) then begin
        rprod(q[level],m+1,H);
        mp_mul(p,q[level],p);
      end;
    end;
    dec(level);
  end;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_prod_int');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if n<2 then begin
    {Easy out: empty product or single element}
    if n=1 then mp_set_int(a,b[low(b)])
    else mp_set1(a);
    exit;
  end;
  {value of thresh is not critical}
  thresh := mp_mul_cutoff*DIGIT_BIT div 31;
  {initialize recursion level}
  level := 0;
  {initialize recursion stack}
  mp_init_multi(q);
  rprod(a,low(b),low(b)+n-1);
  mp_clear_multi(q);
end;


{---------------------------------------------------------------------------}
procedure mp_product(var a: mp_int; n,m: longint);
  {-Return the product a = n*(n+1)*...(m-1)*m, a=1 if m<n}
const
  thresh = 16;                {threshold for recursive calls}
var
  level:  integer;            {recursion level}
  q: array[0..32] of mp_int;  {stack for recursion}

  procedure rprod(var p: mp_int; L,R: longint);
    {-Calculate p=L*(L+1]*...(R-1)*R}
  var
    m: longint;
  begin
    {This is a simplified version of mp_prod_int}
    if mp_error<>MP_OKAY then exit;
    {increase recursion level}
    inc(level);
    if R-L<=thresh then begin
      if L>R then mp_set1(p)
      else begin
        mp_set_int(p,L);
        {skip loop if p is zero}
        while (L<R) and (p.used>0) and (mp_error=MP_OKAY) do begin
          inc(L);
          s_mp_mul_int(p,L,p,q[0]);
        end;
      end;
    end
    else begin
      {recursive calls}
      m := (L + R) div 2;
      rprod(p,L,m);
      {calculate second half only if p<>0}
      if (p.used>0) and (mp_error=MP_OKAY) then begin
        rprod(q[level],m+1,R);
        mp_mul(p,q[level],p);
      end;
    end;
    dec(level);
  end;

begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_product');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if m<n then begin
    mp_set1(a);
    exit;
  end
  else if m=n then begin
    mp_set_int(a,n);
    exit;
  end;

  {here n < m}
  if (n<=0) and (m>=0) then begin
    {one factor is zero}
    mp_zero(a);
    exit;
  end;

  {initialize recursion level}
  level := 0;
  {initialize recursion stack}
  mp_init_multi(q);
  rprod(a,n,m);
  mp_clear_multi(q);
end;


{---------------------------------------------------------------------------}
function mp_radix_size(const a: mp_int; radix: word): longint;
  {-Return size of ASCII representation (incl. sign and #0)}
var
  digs: longint;
  sz: integer;
begin
  mp_radix_size := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_radix_size');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {check range of the radix}
  if (radix < 2) or (radix > MAXRadix) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_radix_size: radix out of range');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {sz: number of chars for sign plus the NULL byte that would be required.}
  if (a.sign=MP_NEG) or mp_show_plus then sz:=2 else sz:=1;

  {get bits}
  digs := mp_bitsize(a);
  if radix<>2 then begin
    {radix<>2, adjust digs}
    digs := 1+trunc(digs*lograd[radix]);
  end;
  mp_radix_size := digs + sz;
end;


{---------------------------------------------------------------------------}
function mp_radix_str(const a: mp_int; radix: word): mp_string;
  {-Convert to radix representation, max 255 digits}
var
  i: integer;
  ls: longint;
  iostr: string[255];
begin
  {arg checks are done by mp_radix_size}
  mp_radix_str := '';
  ls := mp_radix_size(a, radix);
  if mp_error<>MP_OKAY then exit;
  if ls<=255 then mp_toradix_n(a, @iostr[1], radix, 255);
  if (ls>255) or (mp_error<>MP_OKAY) then exit;
  iostr[0] := #255;
  for i:=1 to 255 do begin
    if iostr[i]=#0 then begin
      iostr[0] := char8(i-1);
      break;
    end;
  end;
  mp_radix_str := iostr;
end;


{$ifdef BIT32or64}
{---------------------------------------------------------------------------}
function s_mp_radix_astr(const a: mp_int; radix: word; plus: boolean): ansistring;
  {-Convert to radix representation ansistring, plus to show '+'}
var
  i,im,k,ls,ps,rsq: longint;
  t: mp_int;
  w: longint;
  str: pchar8;
  {$ifndef RESULT}
    Result: ansistring;
  {$endif}
  rmap: TRadixCMap;
  pad: boolean;
const
  NRP=32;
var
  rp2: array[0..NRP] of mp_int;

  {--------------------------------------------------------------------}
  procedure fillpad(cnt: longint);
    {-Insert cnt padding zeros}
  begin
    while cnt>0 do begin
      Result[ps] := '0';
      inc(ps);
      dec(cnt);
    end;
  end;

  {--------------------------------------------------------------------}
  procedure convert(const x: mp_int; level: integer);
    {-Recursive fast integer output}
  var
    q,r: mp_int;
    s,dd: longint;
  begin
    if mp_error<>MP_OKAY then exit;
    if level=0 then begin
      if (not mp_is_longint(x,s)) or (s>=rsq) then begin
        {$ifdef MPC_HaltOnError}
          {$ifdef MPC_UseExceptions}
            raise MPXRange.Create('s_mp_radix_astr: internal error');
          {$else}
            RunError(MP_RTE_RANGE);
          {$endif}
        {$else}
          set_mp_error(MP_RANGE);
          exit;
        {$endif}
      end;
      if s=0 then begin
        if pad then begin
          Result[ps] := '0';
          inc(ps);
          Result[ps] := '0';
          inc(ps);
        end;
        exit;
      end;
      dd := s div radix;
      {s := s mod radix}
      dec(s,dd*radix);
      if dd=0 then begin
        if pad then begin
          Result[ps] := '0';
          inc(ps);
        end;
      end
      else begin
        pad := true;
        Result[ps] := rmap[dd];
        inc(ps);
      end;
      if s=0 then begin
        if pad then begin
          Result[ps] := '0';
          inc(ps);
        end;
      end
      else begin
        pad := true;
        Result[ps] := rmap[s];
        inc(ps);
      end;
    end
    else begin
      mp_init2(q,r);
      if mp_error=MP_OKAY then begin
        {Get quotient q and remainder r with x = q*radix^(2*level) + r}
        s_mp_divrem(x,rp2[level],@q,@r);
        dec(level);
        {convert quotient}
        if q.used>0 then convert(q,level)
        else if pad then fillpad(2 shl level);
        mp_clear(q);
        {convert remainder}
        if r.used>0 then convert(r,level)
        else if pad then fillpad(2 shl level);
        mp_clear(r);
      end;
    end;
  end;

begin
  {Fast subquadratic radix conversion function. Uses ideas from Algorithm 1.26}
  {FastIntegerOutput in Brent/Zimmermann [35], Ch. 1.7 Base conversion; and the}
  {GMP [15] documentation, Ch. Algorithms/Radix Conversion/Binary to Radix.}
  s_mp_radix_astr := '';
  {$ifndef RESULT}
    Result := '';  {Some FPC versions give warnings}
  {$endif}

  {Arg checks are done by mp_radix_size}
  ls := mp_radix_size(a, radix)+16;
  if mp_error<>MP_OKAY then exit;

  {easy out if a=0}
  if a.used=0 then begin
    if plus then s_mp_radix_astr := '+0'
    else s_mp_radix_astr := '0';
    exit;
  end;

  if ls<=65000 then begin
    {check if old algorithm should be used, first check tests if radix is power of two}
    if (radix and pred(radix)=0) or (ls<80) or (a.used<=5*mp_radexp[radix]) then begin
      setlength(result,ls);
      k := length(Result);
      if (k>0) and (k<=65000) then begin
        fillchar(Result[1],k,#0);
        str := pchar8(@Result[1]);
        w := k;
        s_mp_toradix_n(a,radix,plus,w,str);
        if mp_error=MP_OKAY then begin
          {trim trailing #0}
          while (k>0) and (Result[k]=#0) do dec(k);
          if k<>length(Result) then SetLength(Result, k);
          {$ifndef RESULT}
            s_mp_radix_astr := Result;
          {$endif}
        end;
        exit;
      end
    end;
  end;

  mp_init_multi(rp2);
  if mp_error<>MP_OKAY then exit;

  mp_init(t);
  if mp_error<>MP_OKAY then begin
    mp_clear_multi(rp2);
    exit;
  end;

  mp_abs(a, t);
  if mp_error=MP_OKAY then begin
    ps := 1;
    setlength(result,ls);
    if a.sign=MP_NEG then begin
      result[ps] := '-';
      inc(ps);
    end
    else if plus then begin
      result[ps] := '+';
      inc(ps);
    end;

    {get local radix map}
    if mp_uppercase or (radix>36) then rmap := mp_ucrmap else rmap := mp_lcrmap;

    {Pre-compute power table radix^(2i), also get starting level im,}
    {i.e the smallest value with radix^(2*im) >= sqrt(t)}
    k := mp_bitsize(t);
    mp_set(rp2[0],radix);
    for i:=1 to NRP do begin
      if 2*mp_bitsize(rp2[i-1]) > k+1 then begin
        im := i-1;
        break;
      end;
      mp_sqr(rp2[i-1], rp2[i]);
      im := i;
    end;

    rsq := sqr(radix);
    {pad: flag to pad with zeros. Initially false for omitting leading zeros}
    pad := false;
    {do the recursive conversion}
    convert(t,im);
    if mp_error=MP_OKAY then begin
      {adjust string length}
      dec(ps);
      setlength(result,ps);
    end
    else Result := '';
  end;
  mp_clear_multi(rp2);
  mp_clear(t);
  {$ifndef RESULT}
    s_mp_radix_astr := Result;
  {$endif}
end;


{---------------------------------------------------------------------------}
function mp_radix_astr(const a: mp_int; radix: word): ansistring;
  {-Convert to radix representation ansistring (set mp_show_plus to show '+')}
begin
  mp_radix_astr := s_mp_radix_astr(a,radix,mp_show_plus);
end;

{$else}

{---------------------------------------------------------------------------}
function mp_radix_astr(const a: mp_int; radix: word): ansistring;
  {-Convert to radix representation 'ansistring' (set mp_show_plus to show '+')}
var
  ls: longint;
  x,y: mp_int;
  s: string;
  spsav: boolean;
const
  CMAX = 248;
begin
  {Arg checks are done by mp_radix_size}
  ls := mp_radix_size(a, radix)+16;
  if ls < 255 then mp_radix_astr := mp_radix_str(a,radix)
  else begin
    mp_init2(x,y);
    if mp_error=MP_OKAY then begin
      {result to long, return <sign>...<abs(a) mod radix^CMAX>}
      mp_set_pow(x,radix,CMAX);
      mp_abs(a,y);
      mp_mod(y,x,y);
      {save mp_show_plus status}
      spsav := mp_show_plus;
      mp_show_plus := false;
      s := mp_radix_str(y,radix);
      mp_show_plus := spsav;
      {left fill with 0}
      while length(s)<CMAX do s := '0'+s;
      {add sign and ellipsis}
      s := '...'+s;
      if a.sign=MP_NEG then s := '-'+s
      else if mp_show_plus then s := '+'+s;
      mp_radix_astr := s;
      mp_clear2(x,y);
    end;
  end;
end;

{$endif}


{---------------------------------------------------------------------------}
procedure mp_rand_ex(var a: mp_int; digits: longint; sethi: boolean);
  {-Make a pseudo-random mp_int of a given digit size, if not sethi}
  { then a[digits-1] may be zero (and a.used will be decremented)}
var
  d: mp_digit;
  i: TNInt;
begin
  {InitCheck in mp_zero}
  if mp_error<>MP_OKAY then exit;

  if digits>=MAXDigits then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXMaxDigits.Create('mp_rand_ex');
      {$else}
        RunError(MP_RTE_OTHER);
      {$endif}
    {$else}
      set_mp_error(MP_MAXDIGITS);
      exit;
    {$endif}
  end;

  mp_zero(a);
  if digits=0 then exit;

  if a.alloc<digits then begin
    mp_grow(a,digits);
    if mp_error<>MP_OKAY then exit;
  end;

  with a do begin
    used := digits;
    if sethi then begin
      {make highest digit non-zero}
      repeat
        d := mp_random_digit;
        if mp_error<>MP_OKAY then break;
      until d<>0;
    end
    else d := mp_random_digit;
    pdigits^[digits-1] := d;
    if digits>1 then begin
      for i:=0 to digits-2 do pdigits^[i] := mp_random_digit;
    end;
  end;
  if not sethi then mp_clamp(a);
end;


{---------------------------------------------------------------------------}
procedure mp_rand(var a: mp_int; digits: longint);
  {-Make a pseudo-random mp_int of a given digit size}
begin
  mp_rand_ex(a, digits, true);
end;


{---------------------------------------------------------------------------}
procedure mp_random(const a: mp_int; var b: mp_int);
  {-Make a pseudo-random mp_int b with  0 <= b < |a|}
var
  bs: longint;
begin
  bs := mp_bitsize(a);
  if mp_error<>MP_OKAY then exit;
  if bs<2 then mp_zero(b)
  else begin
    mp_rand_bits_ex(b, bs, false);
    while (mp_error=MP_OKAY) and not (mp_cmp_mag(b,a)=MP_LT) do begin
      if a.sign=MP_ZPOS then mp_sub(b,a,b)
      else mp_add(b,a,b);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rand_bits_ex(var a: mp_int; bits: longint; sethi: boolean);
  {-Make pseudo-random a with bitsize <= bits, if sethi highest bit is set}
var
  digits: longint;
  rem: integer;
  m1,m2: mp_digit;
begin
  {InitCheck in mp_rand...}
  if bits<=0 then begin
    mp_zero(a);
    exit;
  end;
  {make random digits, here bits>0 and digits>0}
  digits := (bits+(DIGIT_BIT-1)) div DIGIT_BIT;
  rem := bits mod DIGIT_BIT;
  mp_rand(a, digits);
  {$ifdef MPC_USE_Assert}
    {only if assert supported by compiler or debug}
    assert(a.used=digits, MPAF+'a.used=digits in mp_rand_bits_ex');
  {$endif}
  {Mask highest digit and set highest bit}
  if rem=0 then begin
    m2 := mp_digit(1) shl (DIGIT_BIT-1);
    m1 := MP_MASK;
  end
  else begin
    m2 := mp_digit(1) shl (rem-1);
    m1 := m2 or mp_digit(m2-1);
  end;
  {highest bit is untouched if sethi=false}
  if not sethi then m2:=0;
  with a do if used>0 then pdigits^[used-1] := (pdigits^[used-1] and m1) or m2;
  if m2=0 then mp_clamp(a);
end;


{---------------------------------------------------------------------------}
procedure mp_rand_bits(var a: mp_int; bits: longint);
  {-Make a pseudo-random mp_int of a given bit size}
  { Note: bits are not completely random because highest bit is always set}
begin
  mp_rand_bits_ex(a,bits,true);
end;


{---------------------------------------------------------------------------}
procedure mp_rand_radix(var a: mp_int; radix: word; digits: longint);
  {-Make a pseudo-random mp_int of order radix^digits}
var
  d: mp_digit;
begin
  {InitCheck in mp_zero}
  if mp_error<>MP_OKAY then exit;

  {make sure the radix is ok}
  if (radix < 2) or (radix > MAXRadix) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_rand_radix: radix out of range');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_zero(a);
  if digits=0 then exit;

  {get non zero digit}
  repeat
    d := mp_random_radix(radix) and MP_MASK;
    if mp_error<>MP_OKAY then break;
  until d<>0;
  mp_add_d(a, d, a);

  while (digits>1) and (mp_error=MP_OKAY) do begin
    dec(digits);
    mp_mul_d(a, mp_digit(radix), a);
    mp_add_d(a, mp_random_radix(radix) and MP_MASK, a);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_read_decimal(var a: mp_int; str: pchar8);
  {-Read an mp_int from a decimal ASCII pchar}
begin
  {arg check in mp_read_radix}
  mp_read_radix(a,str,10);
end;


{---------------------------------------------------------------------------}
procedure mp_read_decimal_str(var a: mp_int; const s: mp_string);
  {-Read an mp_int from a ASCII string[255] in given radix}
begin
  mp_read_radix_str(a,s,10);
end;


{---------------------------------------------------------------------------}
procedure s_mp_read_radix(var a: mp_int; var str: pchar8; radix: word; const sep: mp_string; SignAllowed: boolean);
  {-Read an ASCII pchar in a given radix into a, breaks on sep and #0.}
  { No init check for a, used for mp_read_radix and mpr_read_radix}
var
  p: integer;
  c: char8;
  neg,OK: boolean;
  d,rp: mp_digit;
  ri,rn: integer;
begin
  if mp_error<>MP_OKAY then exit;

  {make sure the radix is ok}
  if (radix < 2) or (radix > MAXRadix) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('s_mp_read_radix: radix out of range');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {set the integer to the default of zero}
  mp_zero(a);
  neg := false;
  OK  := false;

  {get number of ASCII digits and radix power that fit into an mp_digit}
  rn := mp_radexp[radix];
  rp := mp_radpow[radix];
  ri := 0;
  d  := 0;

  {skip leading white space}
  repeat
    {$ifdef MPC_MAXRadix64}
      c := str^;
      if radix<=36 then c := upcase(c);
    {$else}
      c := upcase(str^);
    {$endif}
    if c=#0 then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('s_mp_read_radix: invalid char in str');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
    inc(str);
  until (c<>' ') and (c<>#9);

  {process first char}
  p := pos(c,mp_ucrmap)-1;
  if (p>=0) and (p<radix) then begin
    {found digit}
    OK := true;
    ri := 1;
    d  := mp_digit(p);
  end
  else begin
    if ((c='-') or (c='+')) and not SignAllowed then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('s_mp_read_radix: invalid number syntax (sign not allowed)');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
    if c='-' then neg := true
    else if c<>'+' then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('s_mp_read_radix: invalid number syntax');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  end;

  {process remaining chars}
  repeat
    if ri<rn then begin
      {$ifdef MPC_MAXRadix64}
        c := str^;
        if radix<=36 then c := upcase(c);
      {$else}
        c := upcase(str^);
      {$endif}
      if (c=#0) or (pos(c,sep)>0) then break;
      inc(str);
      {if c=' ' then continue;}         {???? space and/or group separator ????}
      p := pos(c,mp_ucrmap)-1;
      if (p<0) or (p>=radix) then begin
        {$ifdef MPC_HaltOnError}
          {$ifdef MPC_UseExceptions}
            raise MPXBadArg.Create('s_mp_read_radix: invalid char in str');
          {$else}
            RunError(MP_RTE_BADARG);
          {$endif}
        {$else}
          set_mp_error(MP_BADARG);
          exit;
        {$endif}
      end;
      {found digit}
      OK := true;

      {accumulate next digits in d as long as possible}
      inc(ri);
      if ri=1 then d:= mp_digit(p)
      else d := d*mp_digit(radix) + mp_digit(p);
    end;

    if ri=rn then begin
      {update a with accumulated digits and reset radix power counter}
      mp_mul_d(a,rp, a);
      mp_add_d(a, d, a);
      if mp_error<>MP_OKAY then exit;
      ri := 0;
    end;
  until false;

  if not OK then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('s_mp_read_radix: invalid number syntax');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  if ri>0 then begin
    {process partial accumulated digit}
    rp := radix;
    {calculate radix power, this could have been done in main loop}
    {but there would have been many unnecessary multiplications}
    while ri>1 do begin
      dec(ri);
      rp := rp*radix;
    end;
    mp_mul_d(a,rp, a);
    mp_add_d(a, d, a);
    if mp_error<>MP_OKAY then exit;
  end;

  if neg and (a.used<>0) then a.sign := MP_NEG;
end;


{---------------------------------------------------------------------------}
procedure mp_read_radix(var a: mp_int; str: pchar8; radix: word);
  {-Read an mp_int from a ASCII pchar in given radix}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_read_radix');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mp_read_radix(a,str,radix,#0,true);
end;


{---------------------------------------------------------------------------}
procedure mp_read_radix_arr(var a: mp_int; const pa: array of pchar8; radix: word);
  {-Read an mp_int from concatenated pchar array pa in given radix,}
  { max 65000 chars. Mainly used for 16 bit compatibility with max.}
  { length of string literals = 255 and line length = 127 chars    }
const
  LMAX=$FF00;
type
  tcbuf = array[0..LMAX] of char8;
var
  pb: ^tcbuf;
  pc: pchar8;
  i : integer;
  k : word;
begin
  {arg check in mp_zero}
  mp_zero(a);
  pb := mp_getmem(sizeof(tcbuf));
  if pb<>nil then begin
    k := 0;
    for i:=low(pa) to high(pa) do begin
      pc := pointer(pa[i]);
      while (k<LMAX) and (pc^<>#0) do begin
        pb^[k] := pc^;
        inc(pc);
        inc(k);
      end;
      if k=LMAX then break;
    end;
    pb^[k] := #0;
    mp_read_radix(a, pchar8(pb), radix);
    mp_freemem(pointer(pb),sizeof(tcbuf));
  end;
end;


{$ifdef BIT32or64}
{---------------------------------------------------------------------------}
procedure mp_read_radix_astr(var a: mp_int; const s: ansistring; radix: word);
  {-Read an mp_int from an ansistring in given radix}
begin
  {Note: pAnsiChar(s) is terminated with #0 even if s=''}
  mp_read_radix(a, pAnsiChar(s), radix);
end;


{---------------------------------------------------------------------------}
procedure mp_read_decimal_astr(var a: mp_int; const s: ansistring);
  {-Read an mp_int from a decimal ansistring}
begin
  mp_read_radix_astr(a, s, 10);
end;

{$else}

{---------------------------------------------------------------------------}
procedure mp_read_radix_astr(var a: mp_int; const s: ansistring; radix: word);
  {-Read an mp_int from an ansistring in given radix}
begin
  mp_read_radix_str(a,s,radix);
end;


{---------------------------------------------------------------------------}
procedure mp_read_decimal_astr(var a: mp_int; const s: ansistring);
  {-Read an mp_int from a decimal ansistring}
begin
  mp_read_radix_str(a,s,10);
end;
{$endif}


{---------------------------------------------------------------------------}
procedure mp_read_radix_str(var a: mp_int; const s: mp_string; radix: word);
  {-Read a string [ASCII] in a given radix into a}
var
  buf: array[0..255] of char8;
  ls: integer;
begin
  {checks in mp_read_radix}
  ls := length(s);
  move(s[1],buf,ls);
  buf[ls] := #0;
  mp_read_radix(a,buf,radix);
end;


{---------------------------------------------------------------------------}
procedure mp_read_signed_bin(var a: mp_int; const b; numbytes: longint);
  {-Read signed bin, big endian, first byte is 0=positive or 1=negative}
var
  ba: TArrayOfByte absolute b;
begin
  if mp_error<>MP_OKAY then exit;

  {read magnitude}
  if numbytes=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_read_signed_bin: numbytes=0');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end
  else mp_read_unsigned_bin(a, ba[1], numbytes - 1);
  if mp_error=MP_OKAY then begin
    {first byte is 0 for positive, non-zero for negative}
    if (ba[0]<>0) and (a.used>0) then a.sign := MP_NEG;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_read_unsigned_bin(var a: mp_int; const b; numbytes: longint);
  {-Reads a unsigned mp_int, assumes the msb is stored first [big endian]}
var
  ba: TArrayOfByte absolute b;
  bbuf: mp_word; {buffer for sequential bits}
  n: longint;
  i: TNint;
  lim,bits: integer;
  numbits: longint;
const
  mask: array[0..7] of word = ($0000,$0001,$0003,$0007,$000F,$001F,$003F,$007F);
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_read_unsigned_bin');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_zero(a);
  if numbytes=0 then exit;

  {Calculate number of trailing highest bits, these must go into a  }
  {separate digit, i.e. big endianness is digit-based not bit-based.}
  numbits := longint(numbytes)*8;
  lim := numbits mod DIGIT_BIT;
  if lim=0 then lim := DIGIT_BIT;

  {index of highest digit}
  n := (numbits-1) div DIGIT_BIT;

  {single allocate enough digits}
  if a.alloc<n+2 then begin
    mp_grow(a,n+2);
    if mp_error<>MP_OKAY then exit;
  end;
  a.used := n+1;

  {initialize bit buffer and bit count}
  bits := 0;
  bbuf := 0;

  with a do begin
    for i:=0 to numbytes-1 do begin
      {accumulate next 8 bits}
      bbuf := (bbuf shl 8) or ba[i];
      inc(bits,8);
      if bits >= lim then begin
        dec(bits, lim);
        pdigits^[n] := mp_digit(bbuf shr bits);
        bbuf := bbuf and mask[bits];
        {$ifopt R+}
          {FPC produces RTE 201 for dec(0)}
          if n>0 then dec(n);
        {$else}
          dec(n);
        {$endif}
        {after first (partial) digit set limit bit count to DIGIT_BIT}
        lim := DIGIT_BIT;
      end;
    end;
  end;
  mp_clamp(a);
end;


{---------------------------------------------------------------------------}
procedure mp_reduce(var x: mp_int; const m, mu: mp_int);
  {-Reduce x mod m via Barrett, assumes x<m^2, mu is from mp_reduce_setup}
  { x should be positive otherwise slow mp_mod is used}
var
  um: TNInt;
  q : mp_int;
const
  limit: longint = longint(MP_DIG1 shl pred(DIGIT_BIT));

begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(x) or mp_not_init(m) or mp_not_init(mu) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_reduce');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if s_mp_is_le0(m) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_reduce: m <= 0');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
    if s_mp_is_le0(mu) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_reduce: mu <= 0');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}

  um := m.used;

  {x must be positive otherwise infinite loop, *V0.7.07}
  {allow x<0 instead of exception, *V1.1.03 }
  if x.sign=MP_NEG then begin
    {m is positive, use add if possible}
    if x.used<um then mp_add(x,m,x) else mp_mod(x,m,x);
    exit;
  end;

  {here x>=0; easy out if x<m, includes case x=0 since um>0}
  if x.used<um then exit;

  {here x.used >= um; use mp_mod if x.used > 2*um}
  if x.used-um>um then begin
    mp_mod(x,m,x);
    exit;
  end;

  {q = x}
  mp_init_copy(q, x);
  if mp_error<>MP_OKAY then exit;

  {from HAC algorithm 14.42}
  {q1 = x / b^(k-1)}
  mp_rshd(q, um-1);

  {done if q1=0, *we:0.7.07}
  if q.used>0 then begin
    if (um>limit) or ((um >= mp_mul_cutoff) and (q.used >= mp_mul_cutoff)) then begin
      mp_mul(q, mu, q)
    end
    else begin
      {*V1.1.09: use digs=um-1 to avoid large number of back subtracts below}
      s_mp_mul_high_digs(q, mu, q, um-1);
    end;

    {q3 = q2 / b^(k+1) */}
    mp_rshd(q, um + 1);

    {x = x mod b^(k+1), quick (no division)}
    s_mp_mod_2k(x, DIGIT_BIT * (um + 1), x);

    {q = q * m mod b^(k+1), quick (no division)}
    s_mp_mul_digs(q, m, q, um + 1);

    {x = x - q}
    mp_sub(x,q,x);

    {If x < 0, add b^(k+1) to it}
    if x.sign=MP_NEG then begin
      mp_set(q, 1);
      mp_lshd(q, um + 1);
      mp_add(x, q, x);
    end;

    {Back off if it's too big, use mp_cmp_mag because x and m are positive}
    while (mp_error=MP_OKAY) and (mp_cmp_mag(x,m)<>MP_LT) do begin
      s_mp_sub(x, m, x);
    end;
  end;
  mp_clear(q);
end;


{---------------------------------------------------------------------------}
procedure mp_reduce_setup(var mu: mp_int; const a: mp_int);
  {-Pre-calculate the value required for Barrett reduction.}
  { For a given modulus a it calculates the value required in mu}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(mu) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_reduce_setup');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {Setup local variable. This is paranoid: @a=@mu not useful, but...}
  mp_init_size(t, 2*(a.used+1));
  if mp_error=MP_OKAY then begin
    mp_2expt(t, a.used * 2 * DIGIT_BIT);
    mp_div(t, a, mu);
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_reduce_2k(var a: mp_int; const n: mp_int; d: mp_digit);
  {-Reduce a mod n where n is of the form 2^p-d, @a<>@n}
var
  q: mp_int;
  p: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(n) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_reduce_2k');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if @a=@n then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_reduce_2k: @a=@n');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;

  {$endif}
  if s_mp_is_le0(n) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_reduce_2k: n<=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {if a<0 do a recursive call with -a}
  if a.sign=MP_NEG then begin
    s_mp_chs(a);
    mp_reduce_2k(a,n,d);
    if a.used>0 then begin
      {change sign mod n}
      mp_sub(n,a,a);
    end;
    exit;
  end;

  {done if a<n}
  if a.used<n.used then exit;

  p := mp_bitsize(n);
  mp_init(q);
  if mp_error<>MP_OKAY then exit; {*0.7.01}

  while mp_error=MP_OKAY do begin
    {q = a/2^p, a = a mod 2^p}
    mp_div_2k(a, p, q, @a);
    {Skip trivial case q=0} {*0.7.01}
    if q.used > 0 then begin
      if d<>1 then begin
        {q = q * d}
        mp_mul_d(q, d, q);
      end;
      {a = a + q}
      s_mp_add(a, q, a);
    end;
    if mp_cmp_mag(a, n)=MP_LT then break;
    s_mp_sub(a, n, a);
  end;
  mp_clear(q);
end;


{---------------------------------------------------------------------------}
procedure mp_reduce_2k_setup(const a: mp_int; var d: mp_digit);
  {-Determine setup value d for unrestricted diminished radix reduction, a>=0}
var
  bs: integer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_reduce_2k_setup');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mp_is_le0(a) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_reduce_2k_setup: a<=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  {a is of the form 2^k - digit}
  d := a.pdigits^[0];
  if a.used>1 then
    {if more than one digit, d = 2^DIGIT_BIT-a[0] = (MP_DIGIT_MAX-a[0])+1}
    d := (MP_DIGIT_MAX-d)+1
  else begin
    {d := 2^bitsize32(a[0]) - a[0]}
    bs := bitsize32(d);
    if bs<DIGIT_BIT then begin
      {$ifdef FPC}
        d := mp_digit(longint(1 shl bs)-longint(d));
      {$else}
        d := (1 shl bs)-d;
      {$endif}
    end
    else d := (MP_DIGIT_MAX-d)+1;
  end;
end;


{---------------------------------------------------------------------------}
function mp_reduce_is_2k(const a: mp_int): boolean;
  {-Determine if mp_reduce_2k can be used}
var
  i: TNInt;
  da: mp_digit;
  pa: pmp_digit;
begin
  mp_reduce_is_2k := false;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_reduce_is_2k');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mp_is_le0(a) then exit;
  {if a.used=1 then skip loop}
  if a.used>1 then begin
    pa := @a.pdigits^[1];
    {check all digits except lowest and highest, note a.used>=2!}
    for i:=1 to a.used-2 do begin
      {if pa^ and MP_MASK <> MP_MASK then exit;}
      if pa^ <> MP_MASK then exit;
      inc(pa);
    end;
    {check highest digit: if we shift all the low}
    {bits away the remaining part should be zero }
    da := pa^;
    while odd(da) do da := da shr 1;
    if da<>0 then exit;
  end;
  mp_reduce_is_2k := mp_error=MP_OKAY;
end;


{---------------------------------------------------------------------------}
function mp_result: integer;
  {-Return and reset mp_error}
begin
  mp_result := mp_error;
  set_mp_error(MP_OKAY);
end;


{---------------------------------------------------------------------------}
procedure s_mp_reverse(var s; len: longint);
  {-Reverse an array of char, used for radix code}
var
  sa: packed array[0..$FFF0] of char8 absolute s;
  i: longint;
  t: char8;
begin
  if len>1 then begin
    i := 0;
    dec(len);
    while i<len do begin
      t := sa[i];
      sa[i] := sa[len];
      sa[len] := t;
      inc(i);
      dec(len);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rshd(var a: mp_int; b: longint);
  {-Shift right a certain amount of digits}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_rshd');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}

    end;
  {$endif}
  {if b <= 0 then ignore}
  if b<=0 then exit;
  with a do begin
    {if b > used then simply zero it and return}
    if used<=b then begin
      mp_zero(a);
      exit;
    end;
    {shift digits}
    move(pdigits^[b], pdigits^[0], (used-b)*sizeof(mp_digit));
    {zero the top digits}
    fillchar(pdigits^[used-b],b*sizeof(mp_digit),0);
    dec(used, b);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_rshd2(const a: mp_int; var b: mp_int; cnt: longint);
  {-Set b to a shifted right by cnt digits}
begin
  {Check is mp_copy}
  mp_copy(a,b);
  mp_rshd(b,cnt);
end;


{---------------------------------------------------------------------------}
procedure mp_set(var a: mp_int; b: mp_digit);
  {-Set a to digit b}
begin
  if mp_error<>MP_OKAY then exit;
  {mp_zero does ArgCheck!}
  mp_zero(a);
  b := b and MP_MASK;
  if (b<>0) and (mp_error=MP_OKAY) then with a do begin
    if (pdigits=nil) or (alloc=0) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_set');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end
    else begin
      used := 1;
      pdigits^[0] := b;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_set1(var a: mp_int);
  {-Set a=1}
begin
  if mp_error<>MP_OKAY then exit;
  {mp_zero does ArgCheck!}
  mp_zero(a);
  if mp_error=MP_OKAY then with a do begin
    if (pdigits=nil) or (alloc=0) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_set1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end
    else begin
      used := 1;
      pdigits^[0] := 1;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_setbit(var a: mp_int; n: longint);
  {-Set bit n of a, error if n<0 or n>MP_MAXBIT (1 = bit 0)}
var
  d,k: TNInt;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_setbit');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if (n<0) or (n>MP_MAXBIT) then begin
    {$ifdef MPC_HaltOnArgCheck}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_setbit: (n<0) or (n>MP_MAXBIT)');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  {$ifdef BASM16}
    asm
       mov ax, word ptr [n]
       mov dx, word ptr [n+2]
       mov cx, MP_DIGIT_BIT
       div cx
       mov [k],ax
       mov [d],dx
    end;
  {$else}
    k := n div MP_DIGIT_BIT;
    d := mp_digit(n mod MP_DIGIT_BIT);
  {$endif}
  if k>=a.used then begin
    {greater than current maximum used mp_digit}
    {grow a to accommodate the single bit}
    mp_grow(a,k+1);
    if mp_error<>MP_OKAY then exit;
    {set the used count of where the bit will go}
    a.used := k+1;
  end;
  a.pdigits^[k] := a.pdigits^[k] or (mp_digit(1) shl d);
end;


{$ifdef MP_32BIT}
{---------------------------------------------------------------------------}
procedure mp_set_int(var a: mp_int; b: longint);
  {-Set a to a longint}
var
  hd : mp_digit;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {mp_zero does ArgCheck!}
  if b=0 then begin
    mp_zero(a);
    exit;
  end;
  if b>0 then neg := false
  else begin
    neg := true;
    b := -b;
  end;
  {mp_digit has at least 16 bits, ie a longint max two mp_digits}
  hd := b shr DIGIT_BIT;
  if hd=0 then mp_set(a, mp_digit(b))
  else begin
    mp_set(a, hd);
    mp_lshd(a,1);
    a.pdigits^[0] := mp_digit(b and MP_MASK);
  end;
  if neg then a.sign := MP_NEG;
end;

{$else}

{---------------------------------------------------------------------------}
procedure mp_set_int(var a: mp_int; b: longint);
  {-Set a to a longint}
var
  ba: packed array[0..3] of byte absolute b;
  j: integer;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {mp_zero does ArgCheck!}
  if b=0 then begin
    mp_zero(a);
    exit;
  end;
  if b>0 then neg := false
  else begin
    neg := true;
    b := -b;
  end;
  if (b>0) and (b<=MP_DIGIT_MAX) then mp_set(a,b)
  else begin
    j := 3;
    while ba[j]=0 do dec(j);
    mp_set(a,ba[j]);
    while j>0 do begin
      dec(j);
      mp_shl(a, 8, a);
      if mp_error<>MP_OKAY then exit;
      inc(a.pdigits^[0],mp_digit(ba[j]));
    end;
  end;
  if neg then a.sign := MP_NEG;
end;
{$endif}


{---------------------------------------------------------------------------}
procedure mp_set_pow(var a: mp_int; b,c: longint);
  {-Set a to b^c, a=0 for c<0}
var
  ba: longint;
  i: integer;
begin
  if mp_error<>MP_OKAY then exit;

  {arg check in called routines}
  ba := abs(b);
  if (ba=1) or (c=0) then begin
    mp_set1(a);
    if (b<0) and odd(c) then s_mp_chs(a);
    exit;
  end
  else if (b=0) or (c<0) then begin
    mp_zero(a);
    exit;
  end;
  if (ba<0) or (ba and (ba-1)=0) then begin
    {b is a power of two, ba<0 if b = $80000000!}
    i:=0;
    while ba and 1 = 0 do begin
      inc(i);
      ba := ba shr 1;
    end;
    mp_2expt(a, c*i);
    if (b<0) and odd(c) then s_mp_chs(a);
    exit;
  end;
  {$ifdef MP_32BIT}
    if ba and longint(MP_INV_MASK) = 0 then begin
      s_mp_expt_dl(mp_digit(ba),c,a);
      if (b<0) and odd(c) then s_mp_chs(a);
      exit;
    end;
  {$else}
    if ba<=$FFFF then begin
      s_mp_expt_wl(word(ba),c,a);
      if (b<0) and odd(c) then s_mp_chs(a);
      exit;
    end;
  {$endif}
  mp_set_int(a,b);
  if c<>1 then mp_expt_int(a,c,a);
end;


{---------------------------------------------------------------------------}
procedure mp_set_short(var a: mp_int; b: shortint);
  {-Set a to a shortint}
begin
  if mp_error<>MP_OKAY then exit;
  if b=-128 then mp_set_int(a,b)
  else begin
    {mp_zero does ArgCheck!}
    mp_zero(a);
    if (b<>0) and (mp_error=MP_OKAY) then with a do begin
      if (pdigits=nil) or (alloc=0) then begin
        {$ifdef MPC_HaltOnError}
          {$ifdef MPC_UseExceptions}
            raise MPXNotInit.Create('mp_set_short');
          {$else}
            RunError(MP_RTE_NOTINIT);
          {$endif}
        {$else}
          set_mp_error(MP_NOTINIT);
          exit;
        {$endif}
      end
      else begin
        used := 1;
        pdigits^[0] := abs(b) and MP_MASK;
        if b<0 then a.sign := MP_NEG;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_set_w(var a: mp_int; w: word);
  {-Set a to a word}
var
  hd: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {mp_set does ArgCheck!}
  hd := w shr DIGIT_BIT;
  if hd=0 then mp_set(a, w)
  else begin
    mp_set(a, hd);
    mp_lshd(a,1);
    a.pdigits^[0] := mp_digit(w and MP_MASK);
    mp_clamp(a);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_shl(const a: mp_int; b: longint; var c: mp_int);
  {-Shift left by a certain bit count}
var
  shift, mask, r, d: mp_digit;
  {$ifndef BASM16}
    rr: mp_digit;
  {$endif}
  pc: pmp_digit;
  i: TNInt;
  nu,bd: longint;
begin
  {init check in mp_copy, (does nothing if @a=@c)}
  mp_copy(a, c);
  if mp_error<>MP_OKAY then exit;
  if (a.used=0) or (b<=0) then exit;

  if b>MP_MAXBIT then begin
    {$ifdef MPC_HaltOnArgCheck}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_shl: b>MP_MAXBIT');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  bd := b div DIGIT_BIT;

  nu := c.used + bd + 1;
  if c.alloc < nu then begin
    mp_grow(c, nu);
    if mp_error<>MP_OKAY then exit;
  end;

  {shift by as many digits in the bit count}
  if bd>0 then begin
    mp_lshd(c, bd);
    if mp_error<>MP_OKAY then exit;
  end;

  {shift any bit count < DIGIT_BIT}
  d := b mod DIGIT_BIT;
  if (d<>0) and (bd<c.used) then begin
    {bitmask for carries}
    mask := (1 shl d) - 1;
    {shift for msbs}
    shift := DIGIT_BIT - d;
    pc := pmp_digit(@c.pdigits^[bd]);

    {$ifdef BASM16}
      i := c.used-bd;
      asm
             les  di,[pc]
             mov  dl,byte ptr [d]
             mov  dh,byte ptr [shift]
             sub  bx,bx               {bx=r}
             cld
        @@1: mov  ax,es:[di]
             mov  si,ax
             mov  cl,dh
             shr  ax,cl
             and  ax,[mask]
             xchg ax,bx               {bx=(pc^ shr shift) and mask, ax=r}
             mov  cl,dl
             shl  si,cl
             or   ax,si
             and  ax,MP_MASK
             stosw                    {pc^ := ((pc^ shl d) or r) and MP_MASK}
             dec  [i]
             jnz  @@1
             mov  [r],bx              {store r,pc for final carry adjust}
             mov  word ptr [pc],di
      end;
    {$else}
      {carry}
      r := 0;
      {don't waste time on the zero digits we have just created}
      {with mp_lshd, ie don't touch digits 0..bd-1}
      for i:=bd to c.used-1 do begin
        {get the higher bits of the current word}
        rr := (pc^ shr shift) and mask;
        {shift the current word and OR in the carry}
        pc^ := ((pc^ shl d) or r) and MP_MASK;
        inc(pc);
        {set the carry to the carry bits of the current word}
        r := rr;
      end;
    {$endif}
    {set final carry}
    if r<>0 then begin
      pc^ := r;
      inc(c.used);
    end;
  end;
  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
procedure mp_shlx(const a: mp_int; b: longint; var c: mp_int);
  {-Shift left a b bits if b>=0, shift right |b| if b<0}
begin
  if b>=0 then mp_shl(a,b,c)
  else mp_shr(a,-b,c)
end;


{---------------------------------------------------------------------------}
procedure mp_shl1(var a: mp_int);
  {-Shift left a by 1}
var
  shift, r: mp_digit;
  {$ifndef BASM16}
    rr: mp_digit;
  {$endif}
  pc: pmp_digit;
  i: TNInt;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_shl1');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if a.used=0 then exit;

  if a.alloc<a.used+1 then begin
    if a.pdigits^[pred(a.used)] > (MP_DIGIT_MAX shr 1) then begin
      mp_grow(a, a.used + 1);
      if mp_error<>MP_OKAY then exit;
    end;
  end;

  {bitmask for carries}
  shift := DIGIT_BIT - 1;
  pc := pmp_digit(a.pdigits);

  {$ifdef BASM16}
    i := a.used;
    asm
           les  di,[pc]
           mov  dh,byte ptr [shift]
           sub  bx,bx               {bx=r}
           cld
      @@1: mov  ax,es:[di]
           mov  si,ax
           mov  cl,dh
           shr  ax,cl
           and  ax,1
           xchg ax,bx               {bx=(pc^ shr shift) and mask, ax=r}
           shl  si,1
           or   ax,si
           and  ax,MP_MASK
           stosw                    {pc^ := ((pc^ shl d) or r) and MP_MASK}
           dec  [i]
           jnz  @@1
           mov  [r],bx              {store r,pc for final carry adjust}
           mov  word ptr [pc],di
    end;
  {$else}
    {carry}
    r := 0;
    {don't waste time on the zero digits we have just created}
    {with mp_lshd, ie don't touch digits 0..bd-1}
    for i:=0 to a.used-1 do begin
      {get the higher bits of the current word}
      rr := (pc^ shr shift) and 1;
      {shift the current word and OR in the carry}
      pc^ := ((pc^ shl 1) or r) and MP_MASK;
      inc(pc);
      {set the carry to the carry bits of the current word}
      r := rr;
    end;
  {$endif}
  {set final carry}
  if r<>0 then begin
    pc^ := r;
    inc(a.used);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_shr(const a: mp_int; b: longint; var c: mp_int);
  {-Shift right a, c = a/2^b; c=a if b<=0}
var
  d, mask, shift: mp_digit;
  {$ifndef BASM16}
    r, rr: mp_digit;
  {$endif}
  pc: pmp_digit;
  x: TNInt;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_shr');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if a=0 or b<=0 then we do no work}
  if (b<=0) or (a.used=0) then begin
    mp_copy(a, c);
    exit;
  end;

  if b>MP_MAXBIT then begin
    mp_zero(c);
    exit;
  end;

  {copy}
  mp_copy(a, c);
  if mp_error<>MP_OKAY then exit;

  {shift by as many digits in the bit count}
  if b>=DIGIT_BIT then mp_rshd(c, b div DIGIT_BIT);

  {shift any bit count < DIGIT_BIT}
  d := b mod DIGIT_BIT;
  if (d<>0) and (c.used>0) then begin  {0.2.06: c.used>0}
    mask := (1 shl d) - 1;
    {shift for lsb}
    shift := DIGIT_BIT - d;
    pc := @c.pdigits^[c.used-1];

    {$ifdef BASM16}
      x := c.used;
      asm
           les  di,[pc]
           mov  dl,byte ptr [d]
           mov  dh,byte ptr [shift]
           sub  bx,bx                    {bx=r}
           std
      @@1: mov  ax,es:[di]
           mov  si,ax
           and  ax,[mask]                {ax=rr=pc^ and mask}
           mov  cl,dh
           shl  bx,cl
           xchg ax,bx                    {ax=r shl shift, r=rr}
           mov  cl,dl
           shr  si,cl
           or   ax,si
           stosw                         {pc := (pc^ shr d) or (r shl shift)}
           dec  [x]
           jnz  @@1
           cld
      end;
    {$else}
      {carry}
      r := 0;
      for x:=c.used-1 downto 0 do begin
        {get the lower bits of this word in a temp}
        rr := pc^ and mask;
        {shift the current word and mix in the carry bits from the previous word}
        pc^ := (pc^ shr d) or (r shl shift);
        dec(pc);
        {set the carry to the carry bits of the current word found above}
        r := rr;
      end;
    {$endif}
  end;
  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
procedure mp_shrx(const a: mp_int; b: longint; var c: mp_int);
  {-Shift right a b bits if b>=0, shift left |b| if b<0}
begin
  if b>=0 then mp_shr(a,b,c)
  else mp_shl(a,-b,c)
end;


{---------------------------------------------------------------------------}
function mp_sign(const a: mp_int): integer;
  {-Return sign(a): -1 if a<0, 0 if a=0, +1 if a>0}
begin
  mp_sign := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sign');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if a.used>0 then begin
    if a.sign=MP_NEG then mp_sign := -1 else mp_sign := +1;
  end;
end;


{---------------------------------------------------------------------------}
function mp_signed_bin_size(const a: mp_int): longint;
  {-Get the size in bytes for an signed equivalent}
begin
  {Size should always be < about 32000 = max allocated bytes}
  mp_signed_bin_size := mp_unsigned_bin_size(a) + 1;
end;


{---------------------------------------------------------------------------}
procedure mp_sqr(const a: mp_int; var b: mp_int);
  {-Compute b = a*a}
begin
  if mp_error<>MP_OKAY then exit;
  {Debug checking in called functions}
  if a.used >= mp_t3s_cutoff then s_mp_toom3_sqr(a, b)
  else if a.used >= mp_sqr_cutoff then s_mp_karatsuba_sqr(a, b)
  else s_mp_sqr(a, b);
  b.sign := MP_ZPOS;
end;


{---------------------------------------------------------------------------}
procedure rec_sqrt(const a: mp_int; var b: mp_int);
  {-Compute b=floor(sqrt(a)), a>0, @a <> @b. Internal: no init. checks etc}
var
  t: mp_int;
  s: longint;
begin
  s := mp_bitsize(a);
  if mp_error<>MP_OKAY then exit;
  if s<=31 then begin
    s := mp_get_int(a);
    mp_set_int(b, isqrt32(s));
    exit;
  end;
  mp_init(t);
  if mp_error<>MP_OKAY then exit;
  {Init. approx. = (1+rec_sqrt(a div 2^(2s))*2^s, s = ceil(bitsize(a)/4)}
  s := (s + 3) div 4;
  mp_shr(a,2*s,b);
  rec_sqrt(b,t);
  mp_inc(t);
  mp_shl(t,s,t);
  {Newton iteration; init. approx t >= floor(sqrt(a)), see below}
  repeat
    if mp_error<>MP_OKAY then break;
    mp_exch(t,b);
    {t := (t + (a div t)) div 2}
    mp_div(a,b,t);
    mp_add(b,t,t);
    mp_div_2(t,t);
  until mp_cmp(t,b)<>MP_LT;
  mp_clear(t);
end;

(*--------------------------------------------------------------------------
 The idea for the recursive integer Newton square root algorithm is simple:
 It uses the standard integer Newton algorithm where the initial approximation
 is derived from the integer root of the higher bits of the argument a, if the
 bit size is not too small. The crucial point is to prove that the recursive
 initial approximation is at least floor(sqrt(a)). This is the convergence
 condition for the standard Newton loop.

 PROPOSITION: rec_sqrt(a) computes floor(sqrt(a)) for all integers a >= 1

 PROOF: As rec_sqrt finally uses a standard Newton iteration for integer
 square roots it remains to show that the initial approximation for the
 iteration is >= floor(sqrt(a)), i.e. we must show

 t = t(a) = (1+rec_sqrt(a div 2^(2s))*2^s >= floor(sqrt(a))
 with s = ceil(bitsize(a) div 4)

 The proof is by complete induction:

 1) for a=1, rec_sqrt(1) = 1 = floor(sqrt(a))

 2) Let a>1. If bitsize(a) is less than 32 rec_sqrt is calculated with
    isqrt32. So let a be big enough for recursion. Let R(a) = rec_sqrt(a),
    a' = a div 2^(2s) and s = ceil(bitsize(a) div 4).

 We have the following four facts:

   (1) R(a') = floor(sqrt(a'))        induction hypothesis for a' < a

   (2) R(a')^2 <= a' < (R(a')+1)^2    integer sqrt property

   (3) (R(a')+1)^2 >= a'+1            from (2) because all terms are integers

   (4) (x div y)*y + x mod y = x,     0 <= x mod y < y;  for x>=0, y>0 integer

 t^2  = 2^(2s)*(1+R(a'))^2
     >= 2^(2s)*(a'+1)                   fact (3)
      = 2^(2s)*(a div 2^(2s) + 1)       definition of a'
      = (a div 2^(2s))*2^(2s) + 2^(2s)
      = a - a mod 2^(2s) + 2^(2s)       fact (4)
      = a + (2^(2s) - a mod 2^(2s))
      > a                               since 2^(2s) > a mod 2^(2s)

 and therefore t >= floor(sqrt(a)). QED
----------------------------------------------------------------------------*)


{---------------------------------------------------------------------------}
procedure s_mp_sqrt(const a: mp_int; var b: mp_int);
  {-Compute b = floor(sqrt(a)), a >=0 using recursive integer Newton square root, no init check}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {must be positive}
  if a.sign=MP_NEG then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('s_mp_sqrt: a < 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {easy out for a=0/1}
  if mp_iszero(a) or mp_is1(a) then begin
    mp_copy(a,b);
    exit;
  end;

  if @a=@b then begin
    {Must use temporary local mp_int}
    mp_init(t);
    if mp_error=MP_OKAY then begin
      rec_sqrt(a,t);
      {all the work is already done by rec_sqrt, store result}
      mp_exch(t,b);
      mp_clear(t);
    end;
  end
  else rec_sqrt(a,b);
end;


{---------------------------------------------------------------------------}
procedure s_mp_sqrtrem(const n: mp_int; var s,r: mp_int);
  {-Compute Karatsuba square root s and remainder r of n >= 0, n = s^2 + r, no init check}
var
  x,y: longint;
  q,u: mp_int;
begin

  {Primary reference: P. Zimmermann [34], Karatsuba Square Root.}
  {See also Modern Computer Arithmetic [35], Algorithm 1.12 SqrtRem}

  if mp_error<>MP_OKAY then exit;

  if n.used<mp_sqrt_cutoff then begin
    if mp_is_longint(n,y) then begin
      x := isqrt32(y);
      mp_set_int(s, x);
      mp_set_int(r, y-sqr(x));
    end
    else begin
      {use recursive integer Newton square root algorithm}
      {and compute the remainder manually}
      s_mp_sqrt(n,s);
      mp_sqr(s,r);
      mp_sub(n,r,r);
    end;
  end
  else begin
    mp_init2(q,u);
    { Algorithm SqrtRem(n = a_3*b^3 + a_2*b^2 + a_1*b + a_0)  }
    { Input: 0 <= a_i < b with a_3 >= b/4                     }
    { Output: (s, r) such that s^2 <= n = s^2 + r < (s + 1)^2 }
    {   (s, r) = SqrtRem(a_3*b + a_2)                         }
    {   (q, u) = DivRem(r*b + a_1, 2*s)                       }
    {        s = s*b + q                                      }
    {        r = u*b + a_0 - q^2                              }
    {        if r < 0 then                                    }
    {          r = r + 2*s - 1                                }
    {          s = s - 1                                      }
    {        end                                              }
    { return (s, r)                                           }

    {Here we set x = bitsize(n) div 4, b=2^x}
    x := mp_bitsize(n) shr 2;

    {a_3*b + a_2 = q = n div b^2 = n shr 2x}
    mp_shr(n,x+x,q);

    {(s, r) = SqrtRem(a_3*b + a_2)}
    s_mp_sqrtrem(q,s,r);

    {q = a_1 = (n shr x) mod 2^x}
    mp_shr(n,x,q);
    mp_mod_2k(q,x,q);

    {r*b + a_1}
    mp_shl(r,x,r);
    mp_add(r,q,r);
    mp_shl(s,1,q);

    {(q, u) = DivRem(r*b + a_1, 2*s)}
    mp_divrem(r,q,@q,@u);

    {s = s*b + q}
    mp_shl(s,x,s);
    mp_add(s,q,s);

    {r = u*b - q^2 + a_0,  a_0 = n mod 2^x}
    mp_shl(u,x,r);
    mp_sqr(q,q);
    mp_sub(r,q,r);
    mp_mod_2k(n,x,q);
    mp_add(r,q,r);

    if r.sign=MP_NEG then begin
      {r = r + 2*s - 1}
      mp_add(r,s,r);
      mp_add(r,s,r);
      mp_dec(r);
      {s = s - 1}
      mp_dec(s);
    end;
    mp_clear2(q,u);
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_sqrtrem(const n: mp_int; var s,r: mp_int);
  {-Compute Karatsuba square root s and remainder r of n >= 0, n = s^2 + r}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(n) or mp_not_init(s) or mp_not_init(r) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sqrtrem');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {must be positive}
  if n.sign=MP_NEG then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_sqrtrem: n < 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  if @s=@r then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXBadArg.Create('mp_sqrtrem: @s=@r');
      {$else}
        RunError(MP_RTE_BADARG);
      {$endif}
    {$else}
      set_mp_error(MP_BADARG);
      exit;
    {$endif}
  end;

  {easy outs for n=0/1}
  if mp_iszero(n) or mp_is1(n) then begin
    mp_copy(n,s);
    mp_zero(r);
    exit;
  end;

  if (@n=@s) or (@n=@r) then begin
    mp_init_copy(t,n);
    if mp_error=MP_OKAY then begin
      s_mp_sqrtrem(t,s,r);
      mp_clear(t);
    end
  end
  else s_mp_sqrtrem(n,s,r);
end;


{---------------------------------------------------------------------------}
procedure mp_sqrt(const a: mp_int; var b: mp_int);
  {-Compute b = floor(sqrt(a)), a >=0 using Karatsuba or recursive Newton}
var
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sqrt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {must be positive}
  if a.sign=MP_NEG then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mp_sqrt: a < 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  if a.used>=mp_sqrt_cutoff then begin
    {Must use temporary local mp_int}
    mp_init(t);
    if mp_error=MP_OKAY then begin
      mp_sqrtrem(a,b,t);
      mp_clear(t);
    end;
  end
  else s_mp_sqrt(a,b);
end;


{---------------------------------------------------------------------------}
procedure mp_sub(const a,b: mp_int; var c: mp_int);
  {-High level subtraction (handles signs)}
var
  cmp: integer;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sub');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if a.sign<>b.sign then begin
    {subtract a negative from a positive, OR}
    {subtract a positive from a negative.   }
    {In either case, ADD their magnitudes,  }
    {and use the sign of the first number.  }
    c.sign := a.sign;
    s_mp_add(a, b, c);
  end
  else begin
    {subtract a positive from a positive, OR }
    {subtract a negative from a negative.    }
    {First, take the difference between their}
    {magnitudes, then... }
    {make s_mp_sub arguments unequal in magnitude}
    {and easy out our equal} {*0.1.29}
    cmp := mp_cmp_mag(a, b);
    if cmp=MP_EQ then mp_zero(c)
    else if cmp=MP_GT then begin
      {Copy the sign from the first}
      c.sign := a.sign;
      {The first has a larger or equal magnitude}
      s_mp_sub(a, b, c);
    end
    else begin
      {The result has the *opposite* sign from the first number}
      c.sign := a.sign xor (MP_NEG xor MP_ZPOS);
      {The second has a larger magnitude}
      s_mp_sub(b, a, c);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_sub_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Single digit subtraction}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_sub_d');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if b>MP_DIGIT_MAX then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mp_sub_d: b>MP_DIGIT_MAX');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;
  {$endif}
  if b=0 then mp_copy(a,c) else s_mp_sub_d(a,b,c);
end;


{---------------------------------------------------------------------------}
procedure mp_sub_int(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a - b}
var
  t: mp_int;
  d: mp_digit;
begin
  if abs(b)<=lv_digit_max then begin
    {$ifdef MPC_ArgCheck}
      if mp_not_init(a) or mp_not_init(c) then begin
        {$ifdef MPC_HaltOnArgCheck}
          {$ifdef MPC_UseExceptions}
            raise MPXNotInit.Create('mp_sub_int');
          {$else}
            RunError(MP_RTE_NOTINIT);
          {$endif}
        {$else}
          set_mp_error(MP_NOTINIT);
          exit;
        {$endif}
      end;
    {$endif}
    d := abs(b);
    if b<0 then s_mp_add_d(a,d,c)
    else if b>0 then s_mp_sub_d(a,d,c)
    else mp_copy(a,c);
  end
  else begin
    mp_init_set_int(t,b);
    {Debug checking in called functions}
    if mp_error=MP_OKAY then begin
      mp_sub(a,t,c);
      mp_clear(t);
    end;
  end;
end;


{---------------------------------------------------------------------------}
function mp_to_signed_bin_n(const a: mp_int; var b; n: longint): longint;
  {-Store in signed big-endian format, max n bytes; return no. of bytes stored}
var
  ba: packed array[0..$F000] of byte absolute b;
begin
  mp_to_signed_bin_n := 0;
  if mp_error<>MP_OKAY then exit;
  if n>0 then begin
    {Checking and processing done mp_to_unsigned_bin_n}
    mp_to_signed_bin_n := 1+mp_to_unsigned_bin_n(a, ba[1], n-1);
    if mp_error=MP_OKAY then begin
      {store sign byte}
      if a.sign=MP_ZPOS then ba[0]:=0 else ba[0]:=1;
    end;
  end;
end;


{---------------------------------------------------------------------------}
function mp_to_unsigned_bin_n(const a: mp_int; var b; n: longint): longint;
  {-Store in unsigned big-endian format, max n bytes; return no. of bytes stored}
var
  ba: TArrayOfByte absolute b;
  i,x: longint;
  bits: integer;
  bbuf: mp_word;
label
  maxn;
begin
  mp_to_unsigned_bin_n := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_to_unsigned_bin_n');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mp_iszero(a) or (n=0) then exit;

  {initialize bit buffer, bit counter, and byte counter}
  bbuf := 0;
  bits := 0;
  x    := 0;

  for i:=0 to a.used-1 do begin
    {accumulate next digit into bit buffer}
    bbuf := bbuf or mp_word(a.pdigits^[i]) shl bits;
    inc(bits, DIGIT_BIT);
    while bits>=8 do begin
      if x<n then begin
        ba[x] := byte(bbuf and $FF);
        inc(x);
      end
      else goto maxn;
      dec(bits,8);
      bbuf := bbuf shr 8;
    end;
  end;

  {here bbuf must be checked! bits may be non zero,}
  {even if all the remaining highest bits are zero!}
  while (bbuf<>0) and (x<n) do begin
    ba[x] := byte(bbuf and $FF);
    inc(x);
    bbuf := bbuf shr 8;
  end;

maxn:

  {Clamp high order zero bytes}
  while (x>0) and (ba[x-1]=0) do dec(x);
  if x>1 then s_mp_reverse(b, x);

  mp_to_unsigned_bin_n := x;
end;


{---------------------------------------------------------------------------}
function mp_todouble(const a: mp_int): double;
  {-Convert a to double, +-inf if too large}
begin
  mp_todouble := mp_todouble_ex(a,0);
end;


{---------------------------------------------------------------------------}
function mp_todouble_ex(const a: mp_int; x: longint): double;
  {-Convert a*2^x to double, +-inf if too large}
var
  bs: longint;
  d: double;
  i,i0: TNInt;
const
  umin = (54 div DIGIT_BIT) + 2; {minimum a.used for skipping low digits}
begin
  mp_todouble_ex := 0.0;
  if (mp_error<>MP_OKAY) or mp_is0(a) then exit;

  {if bitsize too large return +-inf}
  bs := mp_bitsize(a)+x;
  if bs>1024 then begin
    if a.sign=MP_NEG then mp_todouble_ex := DblNegInf
    else mp_todouble_ex := DblPosInf;
  end
  else if bs >= -1074 then begin
    {loop through significant digits of a}
    d := 0.0;
    {get index of lowest digit to accumulate}
    if a.used>umin then i0 := a.used-umin else i0:=0;
    for i:=a.used-1 downto i0 do d := ldexpd(d,DIGIT_BIT) + a.pdigits^[i];
    d := ldexpd(d,i0*DIGIT_BIT+x);
    if a.sign=MP_NEG then d := -d;
    mp_todouble_ex := d;
  end;
end;


{$ifndef EXT64}
{---------------------------------------------------------------------------}
function mp_toextended(const a: mp_int): extended;
  {-Convert a to extended, +-inf if too large}
begin
  mp_toextended := mp_toextended_ex(a,0);
end;


{---------------------------------------------------------------------------}
function mp_toextended_ex(const a: mp_int; x: longint): extended;
  {-Convert a*2^x to extended, +-inf if too large}
var
  bs: longint;
  d: extended;
  i,i0: TNInt;
const
  umin = (64 div DIGIT_BIT) + 2; {minimum a.used for skipping low digits}
  {$ifdef EXT64}
    MAXPOS=1024;
    MAXNEG=-1074;
  {$else}
    MAXPOS=16384;
    MAXNEG=-16445;
  {$endif}
begin
  mp_toextended_ex := 0.0;
  if (mp_error<>MP_OKAY) or mp_is0(a) then exit;
  {if bitsize too large return +-inf}
  bs := mp_bitsize(a)+x;
  if bs>MAXPOS then begin
    if a.sign=MP_NEG then mp_toextended_ex := DblNegInf
    else mp_toextended_ex := DblPosInf;
  end
  else if bs >= MAXNEG then begin
    {loop through significant digits of a}
    d := 0.0;
    {get index of lowest digit to accumulate}
    if a.used>umin then i0 := a.used-umin else i0:=0;
    for i:=a.used-1 downto i0 do d := ldexpx(d,DIGIT_BIT) + a.pdigits^[i];
    d := ldexpx(d,i0*DIGIT_BIT+x);
    if a.sign=MP_NEG then d := -d;
    mp_toextended_ex := d;
  end;
end;
{$endif}

{---------------------------------------------------------------------------}
procedure mp_toradix(const a: mp_int; str: pchar8; radix: word);
  {-Store mp_int as a ASCII string in a given radix, better use mp_toradix_n}
begin
  {NO maxlen check, possible BUFFER OVERFLOW!! }
  mp_toradix_n(a, str, radix, 65000);
end;


{---------------------------------------------------------------------------}
procedure s_mp_toradix_n(const a: mp_int; radix: word; plus: boolean; var maxlen: longint; var str: pchar8);
  {-Convert an mp_int to ASCII for a given radix, plus=show '+', no init check,}
  { after return str points to the final #0 char}
var
  digs: longint;
  s0: pchar8;
  t: mp_int;
  d,rp,rp1: mp_digit;
  i,ix,ri: integer;
  last,p2: boolean;
  prmap: ^TRadixCMap;
begin
  if mp_error<>MP_OKAY then exit;
  if (a.sign=MP_NEG) or plus then digs:=3 else digs:=2;

  {check range of radix/maxlen, here digs = minimum maxlen}
  if (radix < 2) or (radix > MAXRadix) or (maxlen < digs) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('s_mp_toradix_n: radix out of range or maxlen too small');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {easy out if a=0}
  if a.used=0 then begin
    if plus then begin
      str^ := '+';
      inc(str);
    end;
    str^ := '0';
    inc(str);
    str^ := #0;
    exit;
  end;

  mp_init(t);
  if mp_error<>MP_OKAY then exit;

  mp_abs(a, t);
  if mp_error=MP_OKAY then begin
    if a.sign=MP_NEG then begin
      {if a<0 then output '-', dec maxlen}
      str^:='-';
      inc(str);
      dec(maxlen);
    end
    else if plus then begin
      {output '+', dec maxlen}
      str^ := '+';
      inc(str);
      dec(maxlen);
    end;

    {use local pointer to radix map}
    if mp_uppercase or (radix>36) then prmap := @mp_ucrmap else prmap := @mp_lcrmap;

    {get number of ASCII digits and radix power that fit into an mp_digit}
    ri := mp_radexp[radix];
    rp := mp_radpow[radix];

    {remember first digits position}
    s0 := str;

    {initialize digit counter}
    digs := 0;

    {power of two?}
    p2  := mp_is_pow2_d(rp, ix);
    rp1 := pred(rp);
    while (mp_error=MP_OKAY) and (t.used<>0) and (maxlen>1) do begin
      {radix division loop: divide by rp and get ri ASCII digits}
      if p2 then begin
        d := t.pdigits^[0] and rp1;
        mp_shr(t,ix,t);
      end
      else s_mp_div_d(t, rp, @t, d);
      {special flag: no trailing '0' in last chunk}
      last := t.used=0;
      for i:=1 to ri do begin
        str^ := prmap^[d mod radix];
        inc(str);
        d := d div radix;
        inc(digs);
        dec(maxlen);
        {skip if last no-zero digit reached or maxlen-1}
        if (last and (d=0)) or (maxlen<=1) then break;
      end;
    end;

    if mp_error=MP_OKAY then begin
      {reverse the digits part of the string}
      s_mp_reverse(s0^, digs);
      {append a #0 so the string is properly terminated}
      str^ := #0;
    end;
  end;
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_toradix_n(const a: mp_int; str: pchar8; radix: word; maxlen: longint);
  {-Convert an mp_int to an ASCII string for a given radix (2..MAXRadix)}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_toradix_n');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mp_toradix_n(a,radix,mp_show_plus,maxlen,str);
end;


{---------------------------------------------------------------------------}
function mp_unsigned_bin_size(const a: mp_int): longint;
  {-Get the size in bytes for an unsigned equivalent}
var
  size: longint;
begin
  mp_unsigned_bin_size := 0;
  if mp_error<>MP_OKAY then exit;
  {Arg check in mp_bitsize}
  size := (mp_bitsize(a)+7) div 8;
  {Size should always be < about 32000 = max allocated bytes}
  {$ifdef MPC_USE_Assert}
    {only if assert supported by compiler or debug}
    assert(size < $FFFF, MPAF+'size < $FFFF in mp_unsigned_bin_size');
  {$endif}
  mp_unsigned_bin_size := size and $FFFF;
end;


{---------------------------------------------------------------------------}
procedure mp_write_decimal(var tf: system.text; const a: mp_int);
  {-Write decimal representation to file tf}
begin
  mp_write_radix(tf, a, 10);
end;


{---------------------------------------------------------------------------}
procedure mp_writeln(const msg: mp_string; const a: mp_int);
  {-Writeln a to stdout with leading msg}
begin
  write(msg);
  mp_write_decimal(output,a);
  writeln;
end;


{---------------------------------------------------------------------------}
procedure mp_xor(const a,b: mp_int; var c: mp_int);
  {-Calculate c = a xor b}
var
  pa,pb,pc,px: pmp_digit;
  olduse, min, max: longint;
  i: TNInt;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_xor');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if a.used > b.used then max := a.used else max := b.used;

  {easy outs}
  if max=0 then begin
    mp_zero(c);
    exit;
  end;
  if a.used=0 then begin
    mp_copy(b,c);
    exit;
  end;
  if b.used=0 then begin
    mp_copy(a,c);
    exit;
  end;

  {must grow before assigning pointer otherwise realloc can}
  {change the digits memory, ie px points to the old digits}

  if c.alloc < max then begin
    mp_grow(c, max);
    if mp_error<>MP_OKAY then exit;
  end;

  if a.used > b.used then begin
    min := b.used;
    px  := @a.pdigits^[min];
  end
  else begin
    min := a.used;
    px  := @b.pdigits^[min];
  end;

  pa  := pmp_digit(a.pdigits);
  pb  := pmp_digit(b.pdigits);
  pc  := pmp_digit(c.pdigits);

  {get old used digit count and set new one}
  olduse := c.used;
  c.used := max;
  c.sign := a.sign xor b.sign;

  {process digits 0 .. min-1}
  for i:=0 to pred(min) do begin
    pc^ := pa^ xor pb^;
    inc(pa);
    inc(pb);
    inc(pc);
  end;

  {copy digits min .. max-1}
  if min<max then begin
    for i:=min to pred(max) do begin
      pc^ := px^;
      inc(px);
      inc(pc);
    end;
  end;

  {clear c above max if necessary}
  if olduse>max then begin
    for i:=max to pred(olduse) do begin
      pc^ := 0;
      inc(pc);
    end;
  end;
  mp_clamp(c);
end;


{$ifdef BIT16}
{---------------------------------------------------------------------------}
procedure s_mp_write_radix(var tf: system.text; const a: mp_int; radix: word; plus: boolean);
  {-Write radix representation to file tf}
var
  lw,la,ls: longint;
  pc,pt: pchar8;
begin
  {arg checks are done by  mp_radix_size}
  ls := mp_radix_size(a, radix);
  if (ls<=$FF00) and (mp_error=MP_OKAY) then begin
    lw := ls and $ffff;
    pc := mp_getmem(lw);
    if pc<>nil then begin
      pt := pc;
      {save alloc count, lw is changed in s_mp_toradix_n} {1.0.14}
      la := lw;
      s_mp_toradix_n(a,radix, plus, lw, pt);
      write(tf, pc);
      mp_freemem(pointer(pc),la);
    end;
  end;
end;
{$else}
{---------------------------------------------------------------------------}
procedure s_mp_write_radix(var tf: system.text; const a: mp_int; radix: word; plus: boolean);
  {-Write radix representation to file tf}
begin
  if mp_error=MP_OKAY then write(tf,s_mp_radix_astr(a,radix,plus));
end;
{$endif}


{---------------------------------------------------------------------------}
procedure mp_write_radix(var tf: system.text; const a: mp_int; radix: word);
  {-Write radix representation to file tf}
begin
  s_mp_write_radix(tf,a,radix,mp_show_plus);
end;


{---------------------------------------------------------------------------}
procedure mp_zero(var a: mp_int);
  {-Set a to zero}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_zero');
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
    sign := MP_ZPOS;
    used := 0;
    if (alloc>0) and (pdigits<>nil) then begin
      fillchar(pdigits^, sizeof(mp_digit) * alloc, 0);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_add(const a,b: mp_int; var c: mp_int);
  {-Low level addition c=a+b, based on HAC pp.594, algorithm 14.7}
var
  pa,pb,pc,px: pmp_digit;
  olduse, min, max: longint;
  i: TNInt;
  u,t: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_add');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {find sizes, we let |a| <= |b| which means we have to sort}
  {them. "px" will point to the input with the most digits  }

  if a.used > b.used then max := a.used else max := b.used;

  {initialize result}

  {must grow before assigning pointer otherwise realloc }
  {can change the digits memory, ie px points to the old}
  {digits, and (max+1)th digit is stored into nirvana   }

  if c.alloc < max+1 then begin
    mp_grow(c, max + 1);
    if mp_error<>MP_OKAY then exit;
  end;

  if a.used > b.used then begin
    min := b.used;
    px  := @a.pdigits^[min];
  end
  else begin
    min := a.used;
    px  := @b.pdigits^[min];
  end;

  pa  := pmp_digit(a.pdigits);
  pb  := pmp_digit(b.pdigits);
  pc  := pmp_digit(c.pdigits);

  {get old used digit count and set new one}
  olduse := c.used;
  c.used := max + 1;

  {zero the carry}
  u := 0;
  i := 0;
  while i<min do begin
    {Compute the sum at one digit, t[i] = a[i] + b[i] + u}
    t := pa^ + pb^ + u;
    {U = carry bit of t[i]}
    u := t shr mp_digit(DIGIT_BIT);
    {take away carry bit from t[i]}
    pc^ := t and MP_MASK;
    inc(pa);
    inc(pb);
    inc(pc);
    inc(i);
  end;
  {now copy higher words if any, that is in a+b}
  {if a or b has more digits add those in}
  if min<max then begin
    while i<max do begin
      {t[i] = x[i] + u}
      t := px^ + u;
      {u = carry bit of t[i]}
      u := t shr mp_digit(DIGIT_BIT);
      {take away carry bit from t[i]}
      pc^ := t and MP_MASK;
      inc(px);
      inc(pc);
      inc(i);
    end;
  end;

  {add carry}
  pc^ := u;
  inc(pc);
  inc(i);

  {clear digits above oldused}
  while i<olduse do begin
    pc^ := 0;
    inc(pc);
    inc(i);
  end;
  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
procedure s_mp_add_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Single digit addition, no init check, b<>0}
var
  pa, pc: pmp_digit;
  p: mp_int;
  ix: TNInt;
  oldused: longint;
  t, mu: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;

  {grow c as required}
  if c.alloc < a.used + 1 then begin
    mp_grow(c, a.used + 1);
    if mp_error<>MP_OKAY then exit;
  end;
  pa := pmp_digit(a.pdigits);

  {if a is negative and |a| >= b, call c = |a| - b}
  if (a.sign=MP_NEG) and ((a.used>1) or (pa^ >= b)) then begin
    {make positive copy of a, pdigits are from a but unchanged!}
    p := a;
    p.sign := MP_ZPOS;
    {c = |a| - b}
    s_mp_sub_d(p, b, c);
    c.sign := MP_NEG;
    mp_clamp(c);
    exit;
  end;

  {old number of used digits in c}
  oldused := c.used;
  {sign always positive}
  c.sign := MP_ZPOS;
  pc := pmp_digit(c.pdigits);

  {if a is positive}
  if a.sign=MP_ZPOS then begin
    {add digit, after this we're propagating the carry}
    t := pa^ + b;
    inc(pa);
    mu  := t shr DIGIT_BIT;
    pc^ := t and MP_MASK;
    inc(pc);
    ix := 1;
    {now handle rest of the digits, first loop while carry is non-zero}
    while (mu<>0) and (ix<a.used) do begin
      t := pa^ + mu;
      inc(pa);
      mu  := t shr DIGIT_BIT;
      pc^ := t and MP_MASK;
      inc(pc);
      inc(ix);
    end;
    {copy remaining without carry}
    while ix<a.used do begin
      pc^ := pa^;
      inc(pa);
      inc(pc);
      inc(ix);
    end;
    {set final carry}
    inc(ix);
    pc^ := mu;
    {setup size}
    c.used := a.used + 1;
  end
  else begin
    {a was negative and |a| < b}
    c.used := 1;
    {the result is a single digit}
    if a.used=1 then pc^ := b-pa^ else pc^ := b;
    inc(pc);
    {setup count so the clearing of oldused can fall through correctly}
    ix := 1;
  end;

  {now zero up to oldused}
  while ix<oldused do begin
    pc^ := 0;
    inc(pc);
    inc(ix);
  end;

  mp_clamp(c);
end;


{$ifdef BASM16}
{---------------------------------------------------------------------------}
procedure s_mp_checksum(var adler: longint; Msg: pointer; Len: longint);
  {-Update Adler32 checksum with Msg data; init with adler=1}
const
  BASE = 65521; {max. prime < 65536 }
  NMAX =  5552; {max. n with 255n(n+1)/2 + (n+1)(BASE-1) < 2^32}
type
  LH    = packed record
            L,H: word;
          end;
var
  s1,s2: longint;
  n: integer;
begin
  s1 := LH(adler).L;
  s2 := LH(adler).H;
  while Len > 0 do begin
    if Len<NMAX then n := Len else n := NMAX;
    {BASM increases speed from about 52 cyc/byte to about 3.7 cyc/byte}
    asm
                    mov  cx,[n]
            db $66; mov  ax,word ptr [s1]
            db $66; mov  di,word ptr [s2]
                    les  si,[msg]
      @@1:  db $66, $26, $0f, $b6, $1c      {movzx ebx,es:[si]}
                    inc  si
            db $66; add  ax,bx              {inc(s1, PByte(Msg)^)}
            db $66; add  di,ax              {inc(s2, s1}
                    dec  cx
                    jnz  @@1
            db $66; sub  cx,cx
                    mov  cx,BASE
            db $66; sub  dx,dx
            db $66; div  cx
            db $66; mov  word ptr [s1],dx   {s1 := s1 mod BASE}
            db $66; sub  dx,dx
            db $66; mov  ax,di
            db $66; div  cx
            db $66; mov  word ptr [s2],dx   {s2 := s2 mod BASE}
                    mov  word ptr [msg],si  {save offset for next chunk}
    end;
    dec(len, n);
  end;
  LH(adler).L := word(s1);
  LH(adler).H := word(s2);
end;
{$else}
{---------------------------------------------------------------------------}
procedure s_mp_checksum(var adler: longint; Msg: pointer; Len: longint);
  {-Update Adler32 checksum with Msg data; init with adler=1}
const
  BASE = 65521; {max. prime < 65536 }
  NMAX =  3854; {max. n with 255n(n+1)/2 + (n+1)(BASE-1) < 2^31}
type
  LH   = packed record
           L,H: word;
         end;
var
  s1,s2: longint;
  i,n: integer;
begin
  s1 := LH(adler).L;
  s2 := LH(adler).H;
  while Len > 0 do begin
    if Len<NMAX then n := Len else n := NMAX;
    for i:=1 to n do begin
      inc(s1, PByte(Msg)^);
      inc(Ptr2Inc(Msg));
      inc(s2, s1);
    end;
    s1 := s1 mod BASE;
    s2 := s2 mod BASE;
    dec(len, n);
  end;
  LH(adler).L := word(s1);
  LH(adler).H := word(s2);
end;
{$endif}


{---------------------------------------------------------------------------}
procedure s_mp_chs(var a: mp_int);
  {-Change sign of an mp_int, no init check}
begin
  if (mp_error=MP_OKAY) and (a.magic=MP_MAGIC) then with a do begin
    if used>0 then sign := sign xor (MP_NEG xor MP_ZPOS)
    else sign := MP_ZPOS;
  end;
end;


{---------------------------------------------------------------------------}
{------------------------  division routines -------------------------------}
{---------------------------------------------------------------------------}


{---------------------------------------------------------------------------}
function EstimateQDigit(x2,x1,x0,y1,y0: mp_digit): mp_digit;
  {-Calculate q as in Knuth's Algorithm D, step D3}
var
  q,w: mp_word;
begin
  {$ifdef BASM16}
    asm
       db $66;  sub ax,ax
       db $66;  sub dx,dx
       db $66;  sub cx,cx
                mov ax,[x2]
                mov cx,[y1]
       db $66;  shl ax,DIGIT_BIT
                or  ax,[x1]
       db $66;  div cx
       db $66;  mov word ptr [q], ax
       db $66;  mov word ptr [w], dx
    end;
    while (q>MP_MASK) or (q*y0 > fLeftShiftAdd(w,x0)) do begin
      dec(q);
      inc(w,y1);
      if w>=MP_MASK then break;
    end;
  {$else}
    w := (mp_word(x2) shl DIGIT_BIT) or x1;
    q := w div y1;
    w := w mod y1;
    while (q>MP_MASK) or (q*y0 > ((w shl DIGIT_BIT) or x0)) do begin
      dec(q);
      inc(w,y1);
      if w>=MP_MASK then break;
    end;
  {$endif}
  if q>MP_MASK then begin
    {Bugfix in V1.0.03}
    q := MP_MASK;
  end;
  EstimateQDigit := q;
end;


{---------------------------------------------------------------------------}
procedure s_mp_divrem_basecase(const a,b: mp_int; pc,pd: pmp_int);
  {-Integer signed division using Knuth's basecase algorithm D;    }
  { pc^ = a div b, pd^ = a rem b; sign(pd^)=sign(a); no init check.}
label
  __X, __Y, __Z, __Q;
var
  q, x, y, z: mp_int;
  n, t, i, it: TNInt;
  norm: integer;
  asign, bsign: word;
  y0,y1: mp_digit;
{$ifdef MP_16BIT}
  w: word;
{$endif}
begin
  if mp_error<>MP_OKAY then exit;

  { c*b + d == a [e.g. a/b, c=quotient, d=remainder]}
  { HAC[5], algorithm 14.20 and Knuth[3] Algorithm D}

  {if |a| < |b| then c=0, d=a}
  if mp_cmp_mag(a,b)=MP_LT then begin
    if pd<>nil then mp_copy(a,pd^);
    if pc<>nil then mp_zero(pc^);
    exit;
  end;

  {remember signs}
  asign := a.sign;
  bsign := b.sign;

  if b.used<2 then begin
    {single digit b, possibly signed}
    {although b should be <> 0, let mp_div_d handle b=0 correctly}
    if b.used=0 then y0:=0 else y0:=b.pdigits^[0];
    s_mp_div_d(a,y0,pc,y0);
    if pd<>nil then begin
      {sign of remainder = sign(pc^) = sign(a)}
      mp_set(pd^,y0);
      if y0<>0 then pd^.sign := asign;
    end;
    if pc<>nil then with pc^ do begin
      if used<>0 then begin
        if asign=bsign then sign := MP_ZPOS else sign := MP_NEG;
      end;
    end;
    exit;
  end;

{$ifdef MP_16BIT}
  {This bridges the gap between DIGIT_BIT and 16, which reduces}
  {to bitsize=16 with the standard configuration but the gap is}
  {of course larger for smaller values of DIGIT_BIT}
  if (b.used=2) and (bitsize32(b.pdigits^[1])+DIGIT_BIT<17) then begin
    w := (word(b.pdigits^[1]) shl DIGIT_BIT) or b.pdigits^[0];
    mp_div_w(a,w,pc,w);
    if pd<>nil then begin
      {sign of remainder = sign(pc^) = sign(a)}
      mp_set_w(pd^,w);
      if w<>0 then pd^.sign := asign;
    end;
    if pc<>nil then with pc^ do begin
      if used<>0 then begin
        if asign=bsign then sign := MP_ZPOS else sign := MP_NEG;
      end;
    end;
    exit;
  end;
{$endif}

  {Here b.used>1 and a.used>1}

  {step 1: initialize q}
  mp_init_size(q, a.used + 2); if mp_error<>MP_OKAY then exit;
  q.used := a.used + 2;

  mp_init(z);         if mp_error<>MP_OKAY then goto __Q;
  mp_init_copy(x, a); if mp_error<>MP_OKAY then goto __Z;
  mp_init_copy(y, b); if mp_error<>MP_OKAY then goto __X;

  {make local copies positive}
  x.sign := MP_ZPOS;
  y.sign := MP_ZPOS;

  {step 2: normalize both x and y, ensure that y[t] >= b/2, [b=2^DIGIT_BIT]}
  {Note: LTM norm value is off by one in many cases. Fixed in V1.0.03}

  norm := mp_bitsize(y) mod DIGIT_BIT;
  if norm > 0 then begin
    norm := DIGIT_BIT - norm;
    mp_shl(x, norm, x);
    mp_shl(y, norm, y);
  end;

  n := x.used - 1;
  t := y.used - 1;

  {$ifdef MPC_USE_Assert}
    {only if assert supported by compiler or debug}
    assert((n>0) and (t>0), MPAF+'(n>0) and (t>0) in s_mp_div');
    assert(y.pdigits^[t] >= 1 shl (DIGIT_BIT-1), MPAF+'y.pdigits^[t] >= 1 shl (DIGIT_BIT-1)');
  {$endif}

  {while (x >= y*b^n-t) do [ q[n-t] += 1; x -= y*b^(n-t)] }
  mp_lshd2(y, z, n-t);  {z = y*b^(n-t)}
  while mp_cmp(x, z)<>MP_LT do begin
    inc(q.pdigits^[n-t]);
    mp_sub(x, z, x);
    if mp_error<>MP_OKAY then goto __Y;
  end;

  {get upper two digits of y for quotient digit estimation, t>1!}
  y0 := y.pdigits^[t-1];
  y1 := y.pdigits^[t];

  {step 3: for i from n downto (t + 1)}
  for i:=n downto t+1 do begin
    if i>x.used then continue;
    it := i-t-1;
    {step 3.1 and 3.2: estimate and adjust next quotient digit}
    q.pdigits^[it] := EstimateQDigit(x.pdigits^[i],x.pdigits^[i-1],x.pdigits^[i-2],y1,y0);
    {step 3.3: x = x - q[i-t-1] * y * b^(i-t-1)}
    mp_mul_d(y, q.pdigits^[it], z);
    mp_lshd(z, it);
    mp_sub(x, z, x);
    {if x < 0 then [ x = x + y*b^(i-t-1); q[i-t-1] -= 1; ]}
    if x.sign=MP_NEG then begin
      mp_lshd2(y, z, it);
      mp_add(x, z, x);
      {q[i-t-1] -= 1}
      q.pdigits^[it] := (q.pdigits^[it]+MP_MASK) and MP_MASK;
    end;
    if mp_error<>MP_OKAY then goto __Y;
  end; {for i}

  if pc<>nil then begin
    {q is the positive quotient, make it negative if non-zero and sign(a)<>sign(b)}
    mp_clamp(q);
    if (asign<>bsign) and (q.used>0) then q.sign := MP_NEG;
    mp_exch(q, pc^);
  end;

  if pd<>nil then begin
    {x is the remainder which we have to sign-adjust and normalize}
    if x.used>0 then x.sign := asign;
    if norm<>0 then mp_shr(x, norm, x);
    mp_exch(x, pd^);
  end;

__Y: mp_clear(y);
__X: mp_clear(x);
__Z: mp_clear(z);
__Q: mp_clear(q);

end;


{---------------------------------------------------------------------------}
procedure bz_d2n1n(const a,b: mp_int; n: longint; var q,r: mp_int); forward;

{---------------------------------------------------------------------------}
procedure bz_d3n2n(const a12,a3,b,b1,b2: mp_int; n: longint; sb: boolean; var q,r: mp_int);
  {-Divide 3n-bit (a12,a3) by 2n-bit b=(b1,b2), using a indirect recursion via}
  { bz_d2n1n; a,b positive, b<>0, a<2^n*b, internal, no init check; @r=@a12 is}
  { allowed. sb is true if (b1,b2)=2b, ie in the odd(n) case of bz_d2n1n, used}
  { for correction if a - (q*b+r) is negative.}
var
  t: mp_int;
begin
  mp_init(t);
  if mp_error=MP_OKAY then begin
    mp_shr(a12,n,t);
    if mp_cmp_mag(t,b1) < 0 then begin
      {tricky: using a12 directly in bz_d2n1n does not work}
      mp_copy(a12,t);
      bz_d2n1n(t,b1,n,q,r);
    end
    else begin
      mp_shl(b1,n,t);
      mp_sub(a12,t,r);
      mp_add(r,b1,r);
      mp_2expt(q, n);
      mp_dec(q);
    end;
    mp_mul(q,b2,t);
    mp_shl(r,n,r);
    mp_or(r,a3,r);
    mp_sub(r,t,r);
    if r.sign=MP_NEG then begin
      {add back, use b or 2b depending on sb = (shift b flag)}
      if sb then mp_shl(b,1,t)
      else mp_copy(b,t);
      repeat
        if mp_error<>MP_OKAY then break;
        mp_dec(q);
        mp_add(r,t,r);
      until (mp_error<>MP_OKAY) or (r.sign=MP_ZPOS);
    end;
    mp_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure bz_d2n1n(const a,b: mp_int; n: longint; var q,r: mp_int);
  {-Divide 2n-bit a by n-bit b, using a indirect recursion via bz_d3n2n}
  { a,b positive, b<>0, a<2^n*b, internal, no init check}
var
  b1,b2,q1,t: mp_int;
  nh: longint;
begin
  if b.used<=mp_bz_cutoff then begin
    {use Knuth's basecase algorithm D}
    s_mp_divrem_basecase(a,b,@q,@r);
    exit;
  end;
  mp_init4(b1,b2,q1,t);
  if mp_error=MP_OKAY then begin
    if n and 1 = 0 then begin
      {n is even, divide (a1,a2,a3) by (b1,b2)}
      nh := n div 2;
      mp_div_2k(b,nh,b1,@b2);
      mp_shr(a,n,t);
      mp_shr(a,nh,q);
      mp_mod_2k(q,nh,q);
      bz_d3n2n(t,q,b,b1,b2,nh,false,q1,r);
      mp_mod_2k(a,nh,t);
      bz_d3n2n(r,t,b,b1,b2,nh,false,q,r);
      mp_shl(q1,nh,q1);
      mp_or(q1,q,q);
    end
    else begin
      {n is odd, the original Burnikel-Ziegler uses baseline division in this}
      {case. We follow the Python idea and divide 2*(a1,a2,a3) by 2*(b1,b2), }
      {but the shifts are NOT applied to a/b, they are virtual or with temps.}
      inc(n);
      nh := n div 2;
      mp_shl(b,1,t);
      mp_div_2k(t,nh,b1,@b2);
      mp_shl(a,1,t);
      mp_shr(t,n,q1);
      mp_shr(t,nh,q);
      mp_mod_2k(q,nh,q);
      bz_d3n2n(q1,q,b,b1,b2,nh,true,q1,r);
      mp_mod_2k(t,nh,t);
      bz_d3n2n(r,t,b,b1,b2,nh,true,q,r);
      mp_shl(q1,nh,q1);
      mp_or(q1,q,q);
      {q is OK, but remainder must be halved}
      mp_shr1(r);
    end;
  end;
  mp_clear4(b1,b2,q1,t);
end;


{---------------------------------------------------------------------------}
procedure bz_divrem_pos(const a,b: mp_int; calcq: boolean; var q,r: mp_int);
  {-Divide a by b using Burnikel-Ziegler algorithm, a,b>0 e; calcq=true if q}
  { shall be calculated. q,r: no overlap with a,b; internal use, no checks.}
var
  i,m,n: longint;
  s,t: mp_int;
begin

  {Ref: C. Burnikel, J. Ziegler [33]. See also fast_div.py}
  {by M. Dickinson from http://bugs.python.org/issue3451  }

  n := mp_bitsize(b);
  m := mp_bitsize(a) div n;
  mp_init2(s,t);
  if mp_error<>MP_OKAY then exit;

  if calcq then mp_zero(q);
  mp_shr(a,n*m,r);
  mp_mod_2k(r,n,r);
  for i:=m-1 downto 0 do begin
    mp_shr(a,n*i,t);
    mp_mod_2k(t,n,t);
    mp_shl(r,n,r);
    mp_or(r,t,t);
    bz_d2n1n(t,b,n,s,r);
    if calcq then begin
      {accumulate quotient only if requested}
      mp_shl(q,n,q);
      mp_or(q,s,q);
    end;
  end;
  mp_clear2(s,t);
end;


{---------------------------------------------------------------------------}
procedure s_mp_divrem(const a,b: mp_int; pc,pd: pmp_int);
  {-Integer signed division using recursive Burnikel-Ziegler algorithm:}
  { pc^ = a div b, pd^ = a rem b; sign(pd^)=sign(a); no init check.}
  { Knuth's algorithm D is used for bitsizes < mp_bz_cutoff.}
var
  aa,bb,q,r: mp_int;
  asign, bsign: word;
begin

  if (a.used<=b.used) or (b.used<=mp_bz_cutoff) or ((a.used-b.used)<=mp_bz_cutoff) then begin
    {use Knuth's basecase algorithm D}
    s_mp_divrem_basecase(a,b,pc,pd);
    exit;
  end;

  mp_init2(q,r);
  if mp_error<>MP_OKAY then exit;

  {remember signs}
  asign := a.sign;
  bsign := b.sign;

  {make local positive copies of const parameters}
  aa := a; aa.sign := MP_ZPOS;
  bb := b; bb.sign := MP_ZPOS;

  {call Burnikel-Ziegler algorithm for positive integers}
  bz_divrem_pos(aa,bb,pc<>nil,q,r);

  if pc<>nil then begin
    {q is the positive quotient, make it negative if non-zero and sign(a)<>sign(b)}
    if (asign<>bsign) and (q.used>0) then q.sign := MP_NEG;
    mp_exch(q, pc^);
  end;

  if pd<>nil then begin
    {r is the remainder which we have to sign-adjust and normalize}
    if r.used>0 then r.sign := asign;
    mp_exch(r, pd^);
  end;
  mp_clear2(q,r);
end;


{---------------------------------------------------------------------------}
procedure s_mp_expt_dl(a: mp_digit; b: longint; var c: mp_int);
  {-Calculate c = a^b, return 0 for b<0}
var
  m: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {Init check in set routines}
  if (a=1) or (b=0) then begin
    mp_set1(c);
    exit;
  end
  else if (a=0) or (b<0) then begin
    mp_zero(c);
    exit;
  end
  else if a=2 then begin
    mp_2expt(c, b);
    exit;
  end
  else begin
    mp_set(c,a);
    if b>1 then begin
      m := $40000000;
      while m and b = 0 do m := m shr 1;
      m := m shr 1;
      repeat
        mp_sqr(c,c);
        if b and m <> 0 then mp_mul_d(c,a,c);
        m := m shr 1;
      until m=0;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_expt_wl(a: word; b: longint; var c: mp_int);
  {-Calculate c = a^b, return 0 for b<0}
var
  m: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {Init check in set routines}
  if (a=1) or (b=0) then begin
    mp_set1(c);
    exit;
  end
  else if (a=0) or (b<0) then begin
    mp_zero(c);
    exit;
  end
  else if a=2 then begin
    mp_2expt(c, b);
    exit;
  end
  else begin
    mp_set_w(c,a);
    if b>1 then begin
      m := $40000000;
      while m and b = 0 do m := m shr 1;
      m := m shr 1;
      repeat
        mp_sqr(c,c);
        if b and m <> 0 then mp_mul_w(c,a,c);
        m := m shr 1;
      until m=0;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mp_fakeinit(const a: mp_int; i0,i1: longint; var b: mp_int);
  {-Make a positive fake mp_int b = a[i0..i1], b=0 if i0>=a.used, or i0>i1.}
  { Internal use only, no init check of a. DANGER: b uses the memory of a, }
  { therefore use b only as and on CONST parameter, DO NOT clear or GROW b!}
begin
  if mp_error<>MP_OKAY then exit;
  with b do begin
    if (i0>=a.used) or (i0>i1) then begin
      pdigits := a.pdigits;
      used    := 0;
    end
    else begin
      {here a.used>0, and i0<a.used}
      if i1>=a.used then i1 := a.used-1;
      pdigits := PDigitArray(@a.pdigits^[i0]);
      used := 1+i1-i0;
    end;
    alloc  := used;
    sign   := MP_ZPOS;
    magic  := MP_MAGIC;
    mp_clamp(b);
  end;
end;


{---------------------------------------------------------------------------}
function s_mp_is_le0(const a: mp_int): boolean;
  {-Return true if a<=0, no init check}
begin
  s_mp_is_le0 := (a.used=0) or (a.sign=MP_NEG);
end;


{---------------------------------------------------------------------------}
procedure s_mp_karatsuba_mul(const a,b: mp_int; var c: mp_int);
  {-Calculate c = |a| * |b| using Karatsuba multiplication}
var
  n,m: TNInt;
  a0,a1,b0,b1: mp_int;
  t1,t2: mp_int;
begin
  if mp_error<>MP_OKAY then exit;

  {c = |a| * |b| with Karatsuba multiplication using 3 half size}
  {multiplications. With B ~ 2^(mp_bitsize(max(|a|,|b|)/2) write}
  { |a| = a1*B + a0 and |b| = b1*B + b0. Then |a * b| =         }
  { a1b1 * B^2 + ((a1 + a0)(b1 + b0) - (a0b0 + a1b1))*B + a0b0  }
  {This divide-and-conquer algorithm leads to the O(n^lg(3)) or }
  {O(n^1.584) asymptotics (for balanced arguments with n digits)}

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_karatsuba_mul');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Contrary to LTM use maximum of digits}
  if b.used>a.used then n:=b.used else n:=a.used;
  if n=0 then begin
    mp_zero(c);
    exit;
  end;

  n  := (n+1) div 2;
  {B = 2^(n*DIGIT_BIT)}

  mp_init_size2(t1,t2,2*n+4);
  if mp_error<>MP_OKAY then exit;

  {Since @a=@c or @b=@c is allowed, remember size of product}
  m := a.used+b.used;

  {|a| = a1 * B + a0}
  s_mp_fakeinit(a,  0,     pred(n),a0);
  s_mp_fakeinit(a,  n,pred(a.used),a1);

  {|b| = b1 * B + b0}
  s_mp_fakeinit(b,  0,     pred(n),b0);
  s_mp_fakeinit(b,  n,pred(b.used),b1);

  {c = a1*b1*B^2 + ((a1+a0)*(b1+b0) - a1*b1 - a0*b0)*B + a0*b0}
  {  = c2*B^2 + c1*B + c0}

  s_mp_add(a1, a0, t1);      {t1 = a1 + a0}
  s_mp_add(b1, b0, t2);      {t2 = b1 + b0}
  mp_mul(t1, t2, t1);        {t1 = (a1 + a0) * (b1 + b0)}

  {now calc the products c0 = a0*b0 and c2 = a1*b1}
  mp_mul(a0, b0, t2);        {t2 = a0*b0}
  mp_mul(a1, b1, c);         { c = a1*b1}


  {compute c1. s_mp_sub can be used because (a1+a0)*(b1+b0) > a1*b1 + a0*b0}
  s_mp_sub(t1, t2, t1);
  s_mp_sub(t1,  c, t1);      {t1 = (a1+a0)*(b1+b0) - a1*b1 - a0*b0}

  {now grow c if necessary and 'multiply' by powers of B}
  if c.alloc < m then mp_grow(c,m);
  mp_lshd(c, n);
  s_mp_add(c,t1,c);          {c = c2*B + c1}
  mp_lshd(c, n);
  s_mp_add(c,t2,c);          {c = (c2*B + c1)*B + c0}

  mp_clear2(t1,t2);
end;


{---------------------------------------------------------------------------}
procedure s_mp_karatsuba_sqr(const a: mp_int; var b: mp_int);
  {-Karatsuba squaring, compute b = a*a using three half size squarings}
var
  n,m: TNInt;
  a0,a1: mp_int;
  t0,t1: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {See comments of s_mp_karatsuba_mul for details. It is essentially }
  {the same algorithm but merely tuned to perform recursive squarings}

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_karatsuba_sqr');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if a.used=0 then begin
    mp_zero(b);
    exit;
  end;

  n := succ(a.used) div 2;
  {B = 2^(n*DIGIT_BIT)}

  mp_init_size2(t0,t1,2*n+4);
  if mp_error<>MP_OKAY then exit;

  {Since @a=@b is allowed, remember size of product}
  m := a.used shl 1;

  {|a| = a1 * B + a0}
  s_mp_fakeinit(a,  0,     pred(n),a0);
  s_mp_fakeinit(a,  n,pred(a.used),a1);

  {b = a1^2*B^2 + ((a1+a0)^2 - a1^2 - a0^2)*B + a0^2}
  {  = b2*B^2 + b1*B + b0}

  s_mp_add(a1, a0, t1);      {t1 = a1 + a0}

  {now recursively compute the half size squares}
  mp_sqr(t1, t1);            {t1 = (a1 + a0)^2}
  mp_sqr(a0, t0);            {t2 = a0^2}
  mp_sqr(a1, b);             { b = a1^2}

  {compute middle term ('coefficient' of B)}
  s_mp_sub(t1, t0, t1);
  s_mp_sub(t1,  b, t1);      {t1 = (a1+a0)^2 - a0^2 - a1^2}

  {now grow b if necessary and 'multiply' by powers of B}
  if b.alloc < m then mp_grow(b,m);
  mp_lshd(b, n);
  s_mp_add(b,t1,b);          {b = b2*B + b1}
  mp_lshd(b, n);
  s_mp_add(b,t0,b);          {b = (b2*B + b1)*B + b0}

  mp_clear2(t0,t1);
end;


{---------------------------------------------------------------------------}
function s_mp_ln(const a: mp_int): double;
  {-Calculate ln(a), a>0. Result=0 for a<=0}
begin
  s_mp_ln := s_mp_log2(a)*ln(2.0);
end;


{---------------------------------------------------------------------------}
function s_mp_log2(const a: mp_int): double;
  {-Calculate log2(a), a>0. Result=0 for a<=0}
var
  s: double;
  e: integer;
  n,k: TNInt;
const
  MaxExp = MP_DIGIT_BIT+64;
begin
  s_mp_log2 := 0;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_log2');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if s_mp_is_le0(a) then exit;
  { The calculations of log2(a) uses the significant mp_digits of a.}
  { With B=2^MP_DIGIT_BIT, n=a.used-1, a_i = a.pdigits^[i] we have: }
  { a = a_n*B^n + a_(n-1)*B^(n-1) + ... a_1*B + a_0                 }
  {   = B^n*(a_n + a_(n-1)*B^-1 + a_(n-2)*B^-2 + ... =: B^n*C       }
  { log2(a) = n*log2(B) + log2(C) = n*MP_DIGIT_BIT + ln(C)/ln(2)    }
  n := pred(a.used);
  e := 0;
  s := 0.0;
  for k:=n downto 0 do begin
    s := s + ldexpd(a.pdigits^[k], -e);
    e := e + MP_DIGIT_BIT;
    if e > MaxExp then break;
  end;
  s_mp_log2 := n*longint(MP_DIGIT_BIT) + ln(s)/ln(2.0);
end;


{---------------------------------------------------------------------------}
procedure s_mp_mod_2k(const a: mp_int; b: longint; var c: mp_int);
  {-Calculate c = a mod 2^b, -(|a| mod 2^b) if a < 0}
var
  n,m,ndig,nbit: TNInt;
  dp: pmp_digit;
begin

  {Arg check in mp_zero/mp_copy}

  {if b is <= 0 then zero the int}
  if b<=0 then begin
    mp_zero(c);
    exit;
  end;

  mp_copy(a, c);
  if mp_error<>MP_OKAY then exit;

  {if the modulus is larger than the value then return}
  if b >= longint(a.used) * DIGIT_BIT then exit;  {0.7.01: >=}
  {The bug was fixed in LTM 0.33 but unfortunately Tom did not document it}

  ndig := b div DIGIT_BIT;
  nbit := b mod DIGIT_BIT;
  m    := succ(ndig);

  {Flush all the bits above 2^d in its digit}
  dp := @c.pdigits^[ndig];
  dp^ := dp^ and mp_digit((1 shl nbit) - 1);

  if m<c.used then begin
    {*0.9.00}
    {Flush all digits above the one with 2^b in it}
    for n:=m to c.used-1 do begin
      inc(dp);
      dp^ := 0;
    end;
    c.used := m;
  end;

  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
procedure s_mp_mul_digs(const a,b: mp_int; var c: mp_int; digs: longint);
  {-Multiply |a|*|b| and only compute up to digs digits of result}
  { HAC pp. 595, algorithm 14.12  Modified so you can control how   }
  { many digits of output are created.                              }
var
  r: mp_word;
  u, tmpx: mp_digit;
  pt, py: pmp_digit;
  my, ix, iy: TNInt;
  minu: longint;
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_mul_digs');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if a.used<b.used then minu:=a.used else minu:=b.used;

  {Product is zero?}
  if (minu=0) or (digs=0) then begin
    mp_zero(c);
    exit;
  end;

  {*Note: Comba canceled because MP_WARRAY too small for 8/16 bit}
  mp_init_size(t, digs);
  if mp_error<>MP_OKAY then exit;

  t.used := digs;

  {compute the digits of the product directly}
  for ix:=0 to a.used-1 do begin
    {copy of the digit from a used within the nested loop}
    tmpx := a.pdigits^[ix];
    {skip loop if multiplier is zero, *V1.1.09. This is OK}
    {because init_size has already zero-filled  t.pdigits^}
    if (tmpx<>0) and (digs>ix) then begin
      {limit ourselves to making digs digits of output}
      my := digs-ix;
      if b.used<my then my:=b.used;
      {an alias for the destination shifted ix places}
      pt := @t.pdigits^[ix];
      {an alias for the digits of b }
      py := pmp_digit(b.pdigits);
      {$ifdef BASM16}
        iy := my;
        asm
             push  ds
             lds   si,[py]
             les   di,[pt]
             mov   bx,word ptr [tmpx]
             sub   cx,cx
             cld
        @@1: lodsw
             mul   bx
             add   ax,cx
             adc   dx,0
             add   ax,es:[di]
             adc   dx,0
             mov   cx,ax
             and   ax,MP_MASK
             stosw
             db    $0F,$AC,$D1,DIGIT_BIT        {shrd cx,dx,DIGIT_BIT}
             dec   [iy]
             jnz   @@1
             pop   ds
             mov   [u],cx
             mov   word ptr [pt],di
        end;
      {$else}
        {set the carry to zero}
        u := 0;
        {compute the columns of the output and propagate the carry}
        for iy:=0 to my-1 do begin
          {compute the column as a mp_word}
          r := mp_word(pt^) + mp_word(tmpx)*mp_word(py^) + mp_word(u);
          inc(py);
          {the new column is the lower part of the result}
          pt^ := mp_digit(r and mp_word(MP_MASK));
          inc(pt);
          {get the carry word from the result}
          u := mp_digit(r shr mp_word(DIGIT_BIT));
        end;
      {$endif}
      {set carry if it is placed below digs}
      {we: use my because iy may undefined}
      if ix+my < digs then pt^ := u;
    end;
  end;
  mp_clamp(t);
  mp_exch(t, c);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure s_mp_mul_high_digs(const a,b: mp_int; var c: mp_int; digs: longint);
  {-Multiply |a| * |b| and does not compute the lower digs digits}
  { meant to get the higher part of the product}
var
  pt, py: pmp_digit;
  u, tmpx: mp_digit;
  r: mp_word;
  pa, pb: longint;
  ix, iy, id: TNInt;
  t: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_mul_high_digs');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  mp_init_size(t, a.used + b.used + 1);
  if mp_error<>MP_OKAY then exit;
  {Since t is zerofilled we skip the inner loops for zero multiplier}

  t.used := a.used + b.used + 1;
  pa := a.used;
  pb := b.used;

  for ix:=0 to pa-1 do begin
    {left hand side of A[ix] * B[iy]}
    tmpx := a.pdigits^[ix];
    {skip loop if multiplier is zero, *V1.1.09}
    if tmpx<>0 then begin
      {*V1.1.09: Test and adjust indices to allow digs downto zero}
      id := digs-ix; if id<0 then id := 0;
      if pb>id then begin
        {alias to the address of where the digits will be stored }
        pt := @t.pdigits^[id+ix];
        {alias for where to read the right hand side from}
        py := @b.pdigits^[id];
        {clear the carry}
        {$ifdef BASM16}
          iy := pb-id;
          asm
               push  ds
               lds   si,[py]
               les   di,[pt]
               mov   bx,word ptr [tmpx]
               sub   cx,cx                       {u:=0}
               cld
          @@1: lodsw
               mul   bx
               add   ax,cx
               adc   dx,0
               add   ax,es:[di]
               adc   dx,0
               mov   cx,ax
               and   ax,MP_MASK
               stosw
               db    $0F,$AC,$D1,DIGIT_BIT       {shrd cx,dx,DIGIT_BIT}
               dec   [iy]
               jnz   @@1
               mov   es:[di],cx
               pop   ds
          end;
        {$else}
          u := 0;
          for iy:=id to pb-1 do begin
            {calculate the double precision result}
            r := mp_word(pt^) + mp_word(tmpx)*mp_word(py^) + mp_word(u);
            inc(py);
            {get the lower part}
            pt^ := mp_digit(r and mp_word(MP_MASK));
            inc(pt);
            {carry the carry}
            u := mp_digit(r shr mp_word(DIGIT_BIT));
          end;
          pt^ := u;
        {$endif}
      end;
    end;
  end;

  mp_clamp(t);
  mp_exch(t, c);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure s_mp_mul_int(const a: mp_int; b: longint; var c,tmp: mp_int);
  {-Multiply by a 32 bit integer, c=a*b, tmp is an initialized temporary}
var
  ba: longint;
begin
  if mp_error<>MP_OKAY then exit;
  ba := abs(b);
  {$ifdef MP_32BIT}
    if ba and longint(MP_INV_MASK) = 0 then begin
      {arg check in mp_mul_d}
      mp_mul_d(a,mp_digit(ba),c);
      if b<0 then s_mp_chs(c);
      exit;
    end;
  {$else}
    if ba<=$FFFF then begin
      {arg check in mp_mul_w}
      mp_mul_w(a,word(ba),c);
      if b<0 then s_mp_chs(c);
      exit;
    end;
  {$endif}

  {arg check in mp_mul}
  mp_set_int(tmp, b);
  mp_mul(a,tmp,c);
end;


{---------------------------------------------------------------------------}
procedure s_mp_set_dbl(var a: mp_int; x: double; toinf: boolean);
  {-Set a to an double; if toinf, 'round' |a| outward. Error if x=NAN or INF}
type
  txr = packed record
          mant: array[0..7] of byte;
          sexp: word;
        end;
var
  xr: txr absolute x;
  i: integer;
  expo: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_set_dbl');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if xr.sexp and $7fff = $7fff then begin
    {NAN or INF}
    {$ifdef MPC_HaltOnError}
     {$ifdef MPC_UseExceptions}
       raise MPXBadArg.Create('s_mp_set_dbl: NAN or INF');
     {$else}
       RunError(MP_RTE_BADARG);
     {$endif}
   {$else}
     set_mp_error(MP_BADARG);
     exit;
   {$endif}
  end;
  if abs(x)<1.0 then begin
    if toinf and (x<>0.0) then mp_set1(a)
    else mp_zero(a);
  end
  else begin
    {convert 64 bit mantissa}
    mp_set(a, xr.mant[7]);
    for i:=6 downto 0 do begin
      mp_shl(a, 8, a);
      mp_add_d(a, xr.mant[i], a);
    end;
    {shift with exponent}
    expo := longint(xr.sexp and $7fff) -($3fff+63);
    if expo<0 then begin
      mp_shr(a,-expo,a);
      if toinf then mp_inc(a);
    end
    else if expo>=0 then begin
      if toinf then mp_inc(a);
      mp_shl(a,expo,a);
    end;
  end;
  if x<0.0 then s_mp_chs(a);
end;


{---------------------------------------------------------------------------}
procedure s_mp_sqr(const a: mp_int; var b: mp_int);
  {-Low level squaring, b = a*a, HAC pp.596-597, algorithm 14.16}
var
  t: mp_int;
  ix, iy: TNInt;
  pa: longint;
  r: mp_word;
  tmpx,u: mp_digit;
  pt,py: pmp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_sqr');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  pa := a.used;
  if pa<1 then begin
    mp_zero(b);
    exit;
  end;

  mp_init_size(t, 2*pa+1);
  if mp_error<>MP_OKAY then exit;

  {default used is maximum possible size}
  t.used := 2*pa + 1;

  for ix:=0 to pa-1 do begin
    {first calculate the digit at 2*ix}
    tmpx := a.pdigits^[ix];
    {calculate double precision result}
    {*we: sqr(int64) produces failures!?}
    r := mp_word(tmpx)*mp_word(tmpx) + t.pdigits^[2*ix];

    {store lower part in result}
    t.pdigits^[2*ix] := mp_digit(r and MP_MASK);
    {get the carry}
    {$ifdef BASM16}
      asm
        db  $66;  mov ax, word ptr [r]
        db  $66;  shr ax, DIGIT_BIT
                  mov [u], ax
      end;
    {$else}
      u := mp_digit(r shr DIGIT_BIT);
    {$endif}

    {left hand side of A[ix] * A[iy]}
    {alias for where to store the results}
    pt := pmp_digit(@t.pdigits^[2*ix+1]);
    py := pmp_digit(@a.pdigits^[ix+1]);

    {$ifdef BASM16}
      iy := pa-ix-1;
      if iy>0 then begin
        asm
             push  ds
             lds   si,[py]
             les   di,[pt]
             mov   bx,word ptr [tmpx]
             mov   cx,[u]
             cld
        @@1: lodsw                           {ax=py^}
             mul   bx                        {r=tmpx*py^}
             add   ax,ax
             adc   dx,dx                     {2r}
             add   ax,cx
             adc   dx,0                      {2r+u}
             add   ax,es:[di]
             adc   dx,0                      {r := 2r+u+pt^}
             mov   cx,ax
             and   ax,MP_MASK
             stosw                           {pt^ := r and MP_MASK}
             db    $0F,$AC,$D1,DIGIT_BIT     {shrd cx,dx,DIGIT_BIT} {u := r shr DIGIT_BIT}
             dec   [iy]
             jnz   @@1
             jcxz  @@3
        @@2: sub   dx,dx
             add   cx,es:[di]
             adc   dx,dx                     {r = pt^ + u}
             mov   ax,cx
             and   ax,MP_MASK
             stosw                           {pt^ := r and MP_MASK}
             db    $0F,$AC,$D1,DIGIT_BIT     {shrd cx,dx,DIGIT_BIT} {u := r shr DIGIT_BIT}
             jnz   @@2
        @@3: pop   ds
             mov   [u],cx
             mov   word ptr [pt],di
        end;
      end;
      {propagate upwards}
      while u<>0 do begin
        r := mp_word(pt^) + mp_word(u);
        pt^ := mp_digit(r and MP_MASK);
        inc(pt);
        asm
          db  $66;  mov ax, word ptr [r]
          db  $66;  shr ax, DIGIT_BIT
                    mov [u], ax
        end;
      end;
    {$else}
      for iy:=1 to pa-ix-1 do begin
        {first calculate the product}
        r := mp_word(tmpx)*mp_word(py^);
        {now calculate the double precision result, note we use}
        {addition instead of *2 since it's easier to optimize}
        r  :=  mp_word(pt^)+ r + r + mp_word(u);
        { store lower part}
        pt^ := mp_digit(r and MP_MASK);
        inc(py);
        inc(pt);
        u := mp_digit(r shr DIGIT_BIT);
      end;
      {propagate upwards}
      while u<>0 do begin
        r := mp_word(pt^) + mp_word(u);
        pt^ := mp_digit(r and MP_MASK);
        inc(pt);
        u := mp_digit(r shr DIGIT_BIT);
      end;
    {$endif}
  end;

  mp_clamp(t);
  mp_exch(t, b);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure s_mp_sub(const a,b: mp_int; var c: mp_int);
  {-Low level subtraction (assumes |a| >= |b|), HAC pp.595 algorithm 14.9}
var
  pa,pb,pc: pmp_digit;
  olduse, min, max: longint;
  i: TNInt;
  u,t: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_sub');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {find sizes}
  min := b.used;
  max := a.used;

  {$ifdef MPC_USE_Assert}
    {only if assert supported by compiler or debug}
    {Check violation of entry condition |a| >= |b|}
    assert(mp_cmp_mag(a,b)<>MP_LT,MPAF+'|a| >= |b| in s_mp_sub');
  {$endif}

  {initialize result}
  if c.alloc < max then begin
    mp_grow(c, max);
    if mp_error<>MP_OKAY then exit;
  end;

  {get old used digit count and set new one}
  olduse := c.used;
  c.used := max;

  pa  := pmp_digit(a.pdigits);
  pb  := pmp_digit(b.pdigits);
  pc  := pmp_digit(c.pdigits);

  {set carry to zero}
  u := 0;
  i := 0;

  {$Q-,R-} {Temporarily turn off range/overflow checking}
  while i<min do begin
    {t[i] = a[i] - b[i] - u}
    t := pa^ - pb^ - u;              {May cause overflow of mp_digit -> $Q-}
    {u = carry bit (sign bit) of t[i]}
    {Note this saves performing an AND operation since          }
    {if a carry does occur it will propagate all the way to the }
    {MSB.  As a result a single shift is enough to get the carry}

    {take away carry bit from t[i]}
    u := t shr mp_digit((8*sizeof(mp_digit) - 1));
    pc^ := t and MP_MASK;
    inc(pa);
    inc(pb);
    inc(pc);
    inc(i);
  end;

  {now copy higher words if any, e.g. if a has more digits than b}
  while i<max do begin
    {t[i] = a[i] - u}
    t := pa^ - u;                    {May cause overflow of mp_digit -> $Q-}
    {u = carry bit of t[i]}
    u := t shr mp_digit((8*sizeof(mp_digit)-1));
    {take away carry bit from t[i]}
    pc^ := t and MP_MASK;
    inc(pa);
    inc(pc);
    inc(i);
  end;

  {restore range/overflow check status}
  {$ifdef RangeChecks_on}    {$R+} {$endif}
  {$ifdef OverflowChecks_on} {$Q+} {$endif}

  {clear digits above oldused}
  while i<olduse do begin
    pc^ := 0;
    inc(pc);
    inc(i);
  end;
  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
procedure s_mp_sub_d(const a: mp_int; b: mp_digit; var c: mp_int);
  {-Single digit subtraction, no init check, b<>0}
var
  pa, pc: pmp_digit;
  p: mp_int;
  ix: TNInt;
  oldused: longint;
  t, mu: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;
  {grow c as required}
  if c.alloc < a.used + 1 then begin
    mp_grow(c, a.used + 1);
    if mp_error<>MP_OKAY then exit;
  end;

  {if a is negative just do an unsigned addition}
  if a.sign=MP_NEG then begin
    {make positive copy of a, pdigits is from a but unchanged!}
    p := a;
    p.sign := MP_ZPOS;
    s_mp_add_d(p, b, c);
    mp_clamp(c); {V0.6.06}
    c.sign := MP_NEG;
    exit;
  end;

  pa := pmp_digit(a.pdigits);
  pc := pmp_digit(c.pdigits);
  oldused := c.used;

  {if a <= b simply fix the single digit}
  if ((a.used=1) and (pa^ <= b)) or (a.used=0) then begin
   if a.used=1 then pc^ := b - pa^ else pc^ := b;
   inc(pc);
   ix := 1;
   {negative/1digit}
   c.sign := MP_NEG;
   c.used := 1;
  end
  else begin
    {positive/size}
    c.sign := MP_ZPOS;
    c.used := a.used;
    {t := pa^ - .. and t shr .. can produce range errors}
    {so turn off range/overflow check temporarily}
    {$Q-,R-}
    mu := b;
    ix := 0;
    {note: a.used>0, loop is executed at least once}
    while (mu<>0) and (ix<a.used) do begin
      t  := pa^ - mu;
      inc(pa);
      mu  := t shr (8*sizeof(mp_digit) - 1);
      pc^ := t and MP_MASK;
      inc(pc);
      inc(ix);
    end;
    {restore range/overflow check status}
    {$ifdef RangeChecks_on}    {$R+} {$endif}
    {$ifdef OverflowChecks_on} {$Q+} {$endif}
    while ix<a.used do begin
      pc^ := pa^;
      inc(pa);
      inc(pc);
      inc(ix);
    end;
  end;

  {zero excess digits}
  while ix<oldused do begin
    pc^ := 0;
    inc(pc);
    inc(ix);
  end;
  mp_clamp(c);
end;


{---------------------------------------------------------------------------}
procedure s_mp_toom3_mul(const a,b: mp_int; var c: mp_int);
  {-Calculate c = |a| * |b| using Toom-3 multiplication}
var
  w0,w1,w2,w3,w4,a0,a1,a2,b0,b1,b2: mp_int;
  n: TNInt;
  d: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;

  {Multiplication using the Toom-Cook 3-way algorithm, i.e. splitting  }
  {|a| and |b| in three parts, |a| = a2*B^2 + a1*B + a0, similar for b,}
  {with B ~ 2^(bitsize(max(|a|,|b|)/3). Toom-3 is more complicated than}
  {Karatsuba (~Toom-2) but has smaller asymptotic complexity O(n^1.465)}
  {More info: [1]/doc Ch. 5.3.5, [3] Ch. 4.3.3, and [35] Ch 1.3.3. This}
  {implementation uses M. Bodrato's [36] optimal Toom-3 code. See also }
  {http://bodrato.it/toom-cook/ and http://en.wikipedia.org/wiki/Toom3 }

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) or mp_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_toom3_mul');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if (a.used=0) or (b.used=0) then begin
    mp_zero(c);
    exit;
  end;

  {Use splitting parameter B ~ 2^(bitsize(max(|a|,|b|)/3). Note that TSD [1]}
  {uses min(|a|,|b|), but this is inferior and often slower than Karatsuba.}

  if a.used > b.used then n := a.used else n := b.used;
  n  := (n+2) div 3;
  mp_init5(w0,w1,w2,w3,w4);

  {B = 2^(n*DIGIT_BIT)}

  {|a| = a2 * B^2 + a1 * B + a0}
  s_mp_fakeinit(a,  0,     pred(n),a0);
  s_mp_fakeinit(a,  n,   pred(n+n),a1);
  s_mp_fakeinit(a,n+n,pred(a.used),a2);
  {|b| = b2 * B^2 + b1 * B + b0}
  s_mp_fakeinit(b,  0,     pred(n),b0);
  s_mp_fakeinit(b,  n,   pred(n+n),b1);
  s_mp_fakeinit(b,n+n,pred(b.used),b2);

  {w0 = a0+a2}
  s_mp_add(a0, a2, w0);
  {w4 = b0+b2}
  s_mp_add(b0, b2, w4);
  {w3 = w0-a1 = a2-a1+a0}
  mp_sub(w0, a1, w3);
  {w0 = w0+a1 = a2+a1+a0}
  s_mp_add(w0, a1, w0);
  {w2 = w4-b1 = b2-b1+b0}
  mp_sub(w4, b1, w2);
  {w4 = w4+b1 = b2+b1+b0}
  s_mp_add(w4, b1, w4);
  {w1 = w3*w2}
  mp_mul(w3, w2, w1);
  {w2 = w0*w4}
  mp_mul(w0, w4, w2);
  {w0 = 2(w0+a2)-a0 = 4a2+2a1+a0}
  s_mp_add(w0, a2, w0);
  mp_shl1(w0);
  s_mp_sub(w0, a0, w0);
  {w4 = 2(w4+b2)-b0 = 4b2+2b1+b0}
  s_mp_add(w4, b2, w4);
  mp_shl1(w4);
  s_mp_sub(w4, b0, w4);
  {w3 = w0*w4}
  mp_mul(w0, w4, w3);
  {w0 = a0*b0}
  mp_mul(a0, b0, w0);
  {w4 = a2*b2}
  mp_mul(a2, b2, w4);

  {now solve the linear system

   | 1  0  0  0  0 |   |x4|    |w4|
   | 16 8  4  2  1 |   |x3|    |w3|
   | 1  1  1  1  1 | * |x2|  = |w2|
   | 1 -1  1 -1  1 |   |x1|    |w1|
   | 0  0  0  0  1 |   |x0|    |w0|

   for x[i] using 9 subs, 2 shifts, 1 small div.
   x[i] will overlay w[i], obviously x0=w0, x4=w4.
  }

  {w3 = (w3-w1)/3}
  mp_sub(w3, w1, w3);
  s_mp_div_d(w3, 3, @w3, d);
  {w1 = (w2-w1)/2}
  mp_sub(w2, w1, w1);
  mp_shr1(w1);
  {w2 = w2-w0}
  s_mp_sub(w2, w0, w2);
  {w3 = (w3-w2)/2 - 2*w4}
  s_mp_sub(w3, w2, w3);
  mp_shr1(w3);
  s_mp_sub(w3, w4, w3);
  s_mp_sub(w3, w4, w3);
  {w2 = w2-w1-w4}
  s_mp_sub(w2, w4, w2);
  s_mp_sub(w2, w1, w2);
  {w1 = w1-w3}
  s_mp_sub(w1, w3, w1);

  {now grow c if necessary and 'multiply' by powers of B}
  if c.alloc < a.used + b.used then mp_grow(c,a.used + b.used);
  {c = w4*B^4 + w3*B^2 + w2*B^2 + w1*B + w0}
  {  = (((w4*B + w3)*B + w2)*B + w1)*B + w0}
  mp_lshd2(w4,c,n);  s_mp_add(c,w3,c);
  mp_lshd(c,n);      s_mp_add(c,w2,c);
  mp_lshd(c,n);      s_mp_add(c,w1,c);
  mp_lshd(c,n);      s_mp_add(c,w0,c);

  mp_clear5(w0,w1,w2,w3,w4);
end;


{---------------------------------------------------------------------------}
procedure s_mp_toom3_sqr(const a: mp_int; var b: mp_int);
  {-Compute b = a*a using Toom-3 squaring}
var
  w0, w1, w2, w3, a0, a1, a2: mp_int;
  n: TNInt;
  d: mp_digit;
begin
  if mp_error<>MP_OKAY then exit;

  {See s_mp_toom3_mul for algorithm and references}

  {$ifdef MPC_ArgCheck}
    if mp_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mp_toom3_sqr');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  n := a.used;
  if n=0 then begin
    mp_zero(b);
    exit;
  end;

  {Note b is used as w4}
  mp_init4(w0,w1,w2,w3);
  n  := (n+2) div 3;

  {B = 2^(n*DIGIT_BIT)}
  {|a| = a2 * B^2 + a1 * B + a0}
  s_mp_fakeinit(a,  0,     pred(n),a0);
  s_mp_fakeinit(a,  n,   pred(n+n),a1);
  s_mp_fakeinit(a,n+n,pred(a.used),a2);

  {w3 = a0+a2}
  s_mp_add(a0,a2,w3);
  {w1 = w3-a1 = a2-a1+a0}
  mp_sub(w3, a1, w1);
  {w3 = w3+a1 = a2+a1+a0}
  s_mp_add(w3, a1, w3);
  {w2 = w3^2}
  mp_sqr(w3, w2);
  {w3 = 2(w3+a2)-a0 = 4a2+2a1+a0}
  s_mp_add(w3, a2, w3);
  mp_shl1(w3);
  s_mp_sub(w3, a0, w3);
  {w1 = w1^2}
  mp_sqr(w1, w1);
  {w3 = w3^2}
  mp_sqr(w3, w3);
  {w0 = a0^2}
  mp_sqr(a0, w0);
  {w4 = a2^2}
  mp_sqr(a2, b);

  {As for multiplication solve the linear system

   | 1  0  0  0  0 |   |x4|    |w4|
   | 16 8  4  2  1 |   |x3|    |w3|
   | 1  1  1  1  1 | * |x2|  = |w2|
   | 1 -1  1 -1  1 |   |x1|    |w1|
   | 0  0  0  0  1 |   |x0|    |w0|

   for x[i] using 9 subs, 2 shifts, 1 small div.
   x[i] will overlay w[i], obviously x0=w0, x4=w4.
  }

  {w3 = (w3-w1)/3}
  mp_sub(w3, w1, w3);
  s_mp_div_d(w3, 3, @w3, d);
  {w1 = (w2-w1)/2}
  mp_sub(w2, w1, w1);
  mp_shr1(w1);
  {w2 = w2-w0}
  s_mp_sub(w2, w0, w2);
  {w3 = (w3-w2)/2 - 2*w4}
  s_mp_sub(w3, w2, w3);
  mp_shr1(w3);
  s_mp_sub(w3, b, w3);
  s_mp_sub(w3, b, w3);
  {w2 = w2-w1-w4}
  s_mp_sub(w2,  b, w2);
  s_mp_sub(w2, w1, w2);
  {w1 = w1-w3}
  s_mp_sub(w1, w3, w1);

  {now grow b if necessary and 'multiply' by powers of B}
  if b.alloc < (a.used shl 1) then mp_grow(b,a.used shl 1);
  {b = w4*B^4 + w3*B^2 + w2*B^2 + w1*B + w0}
  {  = (((w4*B + w3)*B + w2)*B + w1)*B + w0}
  mp_lshd(b,n);   s_mp_add(b,w3,b);
  mp_lshd(b,n);   s_mp_add(b,w2,b);
  mp_lshd(b,n);   s_mp_add(b,w1,b);
  mp_lshd(b,n);   s_mp_add(b,w0,b);

  mp_clear4(w0,w1,w2,w3);
end;


begin
  mp_set_allocprec(16);    {minimum allocation unit/mask}
end.
