unit mp_ratio;

{Multi precision rational arithmetic routines}

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

 DESCRIPTION   :  Multi precision rational arithmetic routines

 REQUIREMENTS  :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA :  (mp_types)

 MEMORY USAGE  :  heap

 DISPLAY MODE  :  ---

 REFERENCES    :  [ 3] Knuth,D.E.: TAOCP Vol 2, Seminumerical Algorithms
                  [14] IMATH library by M.J. Fromberger
                  [22] LiDIA - A Library for Computational Number Theory, 2006, LiDIA-Group

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 1.0.00   11.04.07  W.Ehrhardt  Initial version: mpr_init, mpr_clear
 1.0.01   11.04.07  we          s_mpr_normalize, mpr_set_int, mpr_(a)decimal
 1.0.02   11.04.07  we          mpr_copy, mpr_abs, mpr_chs, mpr_zero
 1.0.03   11.04.07  we          mpr_add, mpr_sub, mpr_init_xx, mpr_clear_xx
 1.0.04   11.04.07  we          mpr_inv, mpr_mul, mpr_div, mpr_exch
 1.0.05   12.04.07  we          mpr_read_radix
 1.0.06   12.04.07  we          mpr_cmp, mpr_cmp_mag, functions mpr_is_??
 1.0.07   13.04.07  we          mpr_expt, xx.used tests in mpr_is_eq
 1.0.08   13.04.07  we          GCDs in mpr_mul, mpr_div, s_mpr_add_sub
 1.0.09   14.04.07  we          improved mpr_is_ne/ge/le
 1.0.10   14.04.07  we          mpr_toradix_n, mpr_radix_str, mpr_write_radix etc
 1.0.11   15.04.07  we          improved s_mpr_add_sub
 1.0.12   15.04.07  we          improved mpr_write_radix
 1.0.13   15.04.07  we          mpr_set1, mpr_?op?_mpi
 1.0.14   16.04.07  we          speedup s_mpr_add_sub (about 20%)
 1.0.15   22.04.07  we          mpr_read_decimal
 1.0.16   22.04.07  we          mpr_ceil, mpr_floor, mpr_frac, mpr_trunc
 1.0.17   29.04.07  we          mpr_tofloat_n, mpr_tofloat_str
 1.0.18   29.04.07  we          fix rounding in mpr_tofloat_n
 1.0.19   30.04.07  we          s_mpr_tofloat_n, mpr_tofloat_astr
 1.0.20   30.04.07  we          mpr_round, mpr_is0, mpr_is_mp
 1.0.21   01.05.07  we          mpr_todouble
 1.0.22   06.05.07  we          mpr_read_float_radix/decimal, bugfix s_mpr_cmp_mag
 1.0.23   06.05.07  we          mpr_read_double
 1.0.24   11.05.07  we          bugfix mpr_radix_size
 1.0.25   13.05.07  we          Corrected some exception strings
 1.0.26   13.05.07  we          MPAF prefix in assert strings

 1.2.00   07.09.07  we          mpr_checksum
 1.2.01   17.09.07  we          "uses mp_base, mp_numth" in implementation

 1.3.00   04.11.07  we          mpr_init with mp_allocprec
 1.3.01   12.11.07  we          Fix memory leak(s) if MPC_HaltOnError is not defined
 1.3.02   14.11.07  we          mpr_radix_astr: prefill result with #0

 1.5.00   31.01.08  we          {$x+} for VP and D1

 1.6.00   06.06.08  we          function mp_gcd1 is used

 1.7.00   14.09.08  we          s_mpr_normalize no gcd if |num|=1 or |den|=1

 1.8.00   09.10.08  we          use mp_set1

 1.9.00   02.12.08  we          Uses BTypes: char8, pchar8

 1.10.00  21.01.09  we          changes related to (s)mp_divrem

 1.11.01  01.04.09  we          mpr_harmonic, improved mpr_todouble

 1.14.00  13.02.10  we          MPC_MAXRadix64 adjustments

 1.15.00  13.05.10  we          mp_fract_sep
 1.15.01  20.05.10  we          mpr_div_int, mpr_mul_int
 1.15.02  21.05.10  we          mpr_div_2, mpr_mul_2
 1.15.04  22.05.10  we          mpr_div_2k, mpr_mul_2k

 1.20.00  21.01.12  we          mpr_radix_astr/mpr_adecimal without length restriction

 1.24.00  15.12.12  we          Some word types changed to longint
 1.24.01  17.12.12  we          s_mpr_tofloat_n: @r changed to nil in divrem
 1.24.02  03.01.13  we          Fix s_mpr_cmp_mag ds.. variables to longint

 1.37.00  14.06.18  we          mp_reverse changed to s_mp_reverse
 1.37.01  18.06.18  we          Fixed s_mpr_cmp_mag if a=0

**************************************************************************)


(*-------------------------------------------------------------------------
  This code uses material/ideas from the following 3rd party libraries:
   IMATH 1.1+ by Michael J. Fromberger
  See the file '3rdparty.mpa' for the licenses.
----------------------------------------------------------------------------*)

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

procedure mpr_abs(const a: mp_rat; var b: mp_rat);
  {-Absolute value, b = |a|}

procedure mpr_add(const a,b: mp_rat; var c: mp_rat);
  {-Add two mp_rats: c = a+b}

procedure mpr_add_mpi(const a: mp_rat; const b: mp_int; var c: mp_rat);
  {-Add mp_int to mp_rat: c = a+b}

{$ifdef BIT32or64}
function  mpr_adecimal(const a: mp_rat): ansistring;
  {-Convert to decimal ansistring}
{$endif}

procedure mpr_ceil(const a: mp_rat; var b: mp_int);
  {-Return b := ceil(a)}

function  mpr_checksum(const a: mp_rat): longint;
  {-Return a checksum for a, -1 if mp_error<>MP_OKAY, -2 if not initialized}

procedure mpr_chs(const a: mp_rat; var b: mp_rat);
  {-Change sign, b = -a}

function  mpr_cmp(const a,b: mp_rat): integer;
  {-Compare two mp_rats (signed), return sign(a-b)}

function  mpr_cmp_mag(const a,b: mp_rat): integer;
  {-Compare magnitude of two mp_rats (unsigned), return sign(|a|-|b|)}

procedure mpr_clear(var a: mp_rat);
  {-Free an mp_rat}

procedure mpr_clear2(var a,b: mp_rat);
  {-Clear 2 mp_rats}

procedure mpr_clear3(var a,b,c: mp_rat);
  {-Clear 3 mp_rats}

procedure mpr_clear4(var a,b,c,d: mp_rat);
  {-Clear 4 mp_rats}

procedure mpr_clear5(var a,b,c,d,e: mp_rat);
  {-Clear 5 mp_rats}

procedure mpr_clear6(var a,b,c,d,e,f: mp_rat);
  {-Clear 6 mp_rats}

procedure mpr_clear_multi(var vi: array of mp_rat);
  {-Clear a vector of mp_rats}

procedure mpr_clear_multi_p(const pv: array of pmp_rat);
  {-Clear a list of mp_rats given as a ptr vector}

procedure mpr_copy(const a: mp_rat; var b: mp_rat);
  {-Copy an mp_rat, b: = a}

function  mpr_decimal(const a: mp_rat): mp_string;
  {-Convert a to decimal, max 255 chars}

procedure mpr_div(const a,b: mp_rat; var c: mp_rat);
  {-Divide two mp_rats: c = a/b, b<>0}

procedure mpr_div_2(const a: mp_rat; var b: mp_rat);
  {-Divide mp_rat by 2: b = a/2}

procedure mpr_div_2k(const a: mp_rat; k: longint; var b: mp_rat);
  {-Divide mp_rat by 2^k: b = a/2^k}

procedure mpr_div_int(const a: mp_rat; b: longint; var c: mp_rat);
  {-Divide mp_rat by longint: c = a/b, b<>0}

procedure mpr_div_mpi(const a: mp_rat; const b: mp_int; var c: mp_rat);
  {-Divide mp_rat by mp_int: c = a/b, b<>0}

procedure mpr_exch(var a,b: mp_rat);
  {-Exchange two mp_rats}

procedure mpr_expt(const a: mp_rat; b: longint; var c: mp_rat);
  {-Calculate c = a^b, a<>0 for b<0, 0^0=1}

procedure mpr_floor(const a: mp_rat; var b: mp_int);
  {-Return b := floor(a)}

procedure mpr_frac(const a: mp_rat; var b: mp_rat);
  {-Return b := frac(a) = a - trunc(a)}

procedure mpr_harmonic(n: longint; var hn: mp_rat);
  {-Compute the harmonic number hn = sum(1/i, i=1..n) with binary splitting}

procedure mpr_init(var a: mp_rat);
  {-Initialize an mp_rat}

procedure mpr_init2(var a,b: mp_rat);
  {-Initialize 2 mp_rats}

procedure mpr_init3(var a,b,c: mp_rat);
  {-Initialize 3 mp_rats}

procedure mpr_init4(var a,b,c,d: mp_rat);
  {-Initialize 4 mp_rats}

procedure mpr_init5(var a,b,c,d,e: mp_rat);
  {-Initialize 5 mp_rats}

procedure mpr_init6(var a,b,c,d,e,f: mp_rat);
  {-Initialize 6 mp_rats}

procedure mpr_init_copy(var a: mp_rat; const b: mp_rat);
  {-Create a, then copy b into it}

procedure mpr_init_multi(var vi: array of mp_rat);
  {-Initialize a vector of mp_rats.}
  { On error the already initialized mp_rats will be cleared}

procedure mpr_init_multi_p(var pv: array of pmp_rat);
  {-Initialize a list of mp_rats given as a ptr vector.}
  { On error the already initialized mp_rats will be cleared}

procedure mpr_init_size(var a: mp_rat; nsize, dsize: longint);
  {-Initialize an mp_rat to given number of digits}

procedure mpr_inv(const a: mp_rat; var b: mp_rat);
  {-Invert an mp_rat, b = 1/a, a<>0}

function  mpr_is_eq(const a,b: mp_rat): boolean;
  {-Return a = b}

function  mpr_is_ge(const a,b: mp_rat): boolean;
  {-Return a >= b}

function  mpr_is_gt(const a,b: mp_rat): boolean;
  {-Return a > b}

function  mpr_is_le(const a,b: mp_rat): boolean;
  {-Return a <= b}

function  mpr_is_lt(const a,b: mp_rat): boolean;
  {-Return a < b}

function  mpr_is_mp(const a: mp_rat): boolean;
  {-Return true if a is initialized and a.den=1}

function  mpr_is_ne(const a,b: mp_rat): boolean;
  {-Return a <> b}

function  mpr_is0(const a: mp_rat): boolean;
  {-Return a=0}

procedure mpr_mul(const a,b: mp_rat; var c: mp_rat);
  {-Multiply two mp_rats: c = a*b}

procedure mpr_mul_2(const a: mp_rat; var b: mp_rat);
  {-Multiply mp_rat by 2: b = 2*a}

procedure mpr_mul_2k(const a: mp_rat; k: longint; var b: mp_rat);
  {-Multiply mp_rat by 2^k: b = a*2^k}

procedure mpr_mul_int(const a: mp_rat; b: longint; var c: mp_rat);
  {-Multiply mp_rat by longint: c = a*b}

procedure mpr_mul_mpi(const a: mp_rat; const b: mp_int; var c: mp_rat);
  {-Multiply mp_rat and mp_int: c = a*b}

function  mpr_not_init(const a: mp_rat): boolean;
  {-Sanity check if a is initialized, does not catch all cases!}

procedure mpr_output_decimal(const a: mp_rat);
  {-Write decimal representation to output}

procedure mpr_output_radix(const a: mp_rat; radix: word);
  {-Write radix representation to output}

{$ifdef BIT32or64}
function  mpr_radix_astr(const a: mp_rat; radix: word): ansistring;
  {-Convert to radix representation ansistring}
{$endif}

function  mpr_radix_size(const a: mp_rat; radix: word): longint;
  {-Return size of ASCII representation (incl. sign and #0)}

function  mpr_radix_str(const a: mp_rat; radix: word): mp_string;
  {-Convert to radix representation, max 255 digits}

procedure mpr_read_decimal(var a: mp_rat; str: pchar8);
  {-Read a ASCII decimal string radix into a. str may contain a single '/'}

procedure mpr_read_double(var a: mp_rat; d: double);
  {-Convert d to an mp_rat}

procedure mpr_read_float_decimal(var a: mp_rat; str: pchar8);
  {-Read a ASCII float decimal string into a. str may contain a single '.'}

procedure mpr_read_float_radix(var a: mp_rat; str: pchar8; radix: word);
  {-Read a ASCII float radix string into a. str may contain a single '.'}

procedure mpr_read_radix(var a: mp_rat; str: pchar8; radix: word);
  {-Read a ASCII radix string into a. str may contain a single '/'}

procedure mpr_round(const a: mp_rat; var b: mp_int);
  {-Return b := round(a), round(-a)=-round(a), round(n+0.5)=n+1}

procedure mpr_set(var a: mp_rat; const n, d: mp_int);
  {-Set a to n/d, d<>0}

procedure mpr_set1(var a: mp_rat; const n: mp_int);
  {-Set a to n/1}

procedure mpr_set_int(var a: mp_rat; n, d: longint);
  {-Set a to n/d, d<>0}

procedure mpr_sub(const a,b: mp_rat; var c: mp_rat);
  {-Subtract two mp_rats: c = a-b}

procedure mpr_sub_mpi(const a: mp_rat; const b: mp_int; var c: mp_rat);
  {-Subtract mp_int from mp_rat: c = a-b}

function  mpr_todouble(const a: mp_rat): double;
  {-Convert a to double, +-inf if too large, 0 if too small}

{$ifdef BIT32or64}
function  mpr_tofloat_astr(const a: mp_rat; radix, prec: word): ansistring;
  {-Convert to float representation with prec, max 65000 digits}
{$endif}

procedure mpr_tofloat_n(const a: mp_rat; str: pchar8; radix, prec, maxlen: word);
  {-Convert to float format for a given radix, prec digits after '.'}

function  mpr_tofloat_str(const a: mp_rat; radix, prec: word): mp_string;
  {-Convert to float representation, prec digits after '.', max 255 chars}

procedure mpr_toradix_n(const a: mp_rat; str: pchar8; radix: word; maxlen: longint);
  {-Convert an mp_rat to an ASCII string for a given radix (2..MAXRadix)}

procedure mpr_trunc(const a: mp_rat; var b: mp_int);
  {-Return b := trunc(a)}

procedure mpr_write_decimal(var tf: system.text; const a: mp_rat);
  {-Write decimal representation to file tf}

procedure mpr_write_radix(var tf: system.text; const a: mp_rat; radix: word);
  {-Write radix representation to file tf}

procedure mpr_zero(var a: mp_rat);
  {-Set a to zero}



{#Z+}
{---------------------------------------------------------------------------}
{- 'Internal' functions, don't use them unless you know what you are doing -}
{---------------------------------------------------------------------------}
{#Z-}

procedure s_mpr_add_sub(const a,b: mp_rat; var c: mp_rat; sub: boolean);
  {-Add or subtract two mp_rats}

procedure s_mpr_normalize(var a: mp_rat);
  {-Normalize a, assumes a is initialized}

procedure s_mpr_normsign(var a: mp_rat);
  {-Make denominator positive, assumes a is initialized}

procedure s_mpr_qr(const num, den: mp_int; var q,r: mp_int);
  {-Calculate q := abs(num) div abs(den); r := abs(num) mod abs(den), no init check}

procedure s_mpr_tofloat_n(const a: mp_rat; var str: pchar8; radix, prec: word; maxlen: longint);
  {-Convert to float format for a given radix, prec digits after '.', assumes a is initialized}

procedure szf_toradix_n(const a: mp_int; radix,prec,maxlen: word; var str: pchar8);
  {-Convert a positive mp_int to ASCII for a given radix, no init check,}
  { prec digits, leading zeros, on return str points to the final #0 char}


implementation

uses
  mp_base, mp_modul;


{---------------------------------------------------------------------------}
procedure s_mpr_normsign(var a: mp_rat);
  {-Make denominator positive, assumes a is initialized}
begin
  if mp_error<>MP_OKAY then exit;
  with a do begin
    {make denominator positive}
    if den.sign=MP_NEG then begin
      den.sign := MP_ZPOS;
      if num.used>0 then num.sign := num.sign xor (MP_NEG xor MP_ZPOS);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpr_normalize(var a: mp_rat);
  {-Normalize a, assumes a is initialized}
var
  gcd: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  with a do begin
    {sanity check for zero denominator}
    if mp_is0(den) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('s_mpr_normalize: a.den=0');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
        exit;
      {$endif}
    end;

    {normalize 0/? to 0/1}
    if mp_is0(num) then begin
      mp_set1(den);
      exit;
    end;

    if not (mp_is1a(num) or mp_is1a(den)) then begin
      {If the gcd of the numerator and denominator is <> 1 then divide it out}
      mp_init(gcd);
      if mp_error<>MP_OKAY then exit;
      if not mp_gcd1(num, den, gcd) then begin
        mp_div(num, gcd, num);
        mp_div(den, gcd, den);
      end;
      mp_clear(gcd);
    end;

    {make denominator positive}
    if den.sign=MP_NEG then begin
      den.sign := MP_ZPOS;
      if num.used>0 then num.sign := num.sign xor (MP_NEG xor MP_ZPOS);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure s_mpr_add_sub(const a,b: mp_rat; var c: mp_rat; sub: boolean);
  {-Add or subtract two mp_rats}
var
  t1,t2,g: mp_int;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(b) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('s_mpr_add_sup');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if b is zero just copy a}
  if b.num.used=0 then begin
    mpr_copy(a,c);
    exit;
  end;

  {easy out if a is zero}
  if a.num.used=0 then begin
    {if we subtract c := -b else c := b}
    if sub then mpr_chs(b,c) else mpr_copy(b,c);
    exit;
  end;

  if mp_is_eq(a.den, b.den) then begin
    {denominators are equal}
    {c.num =  a.num * b.den + a.den * b.num}
    {c.den =  a.den * b.den                }
    if sub then mp_sub(a.num, b.num, c.num)
           else mp_add(a.num, b.num, c.num);
    mp_copy(a.den, c.den);
    s_mpr_normalize(c);
    exit;
  end;

  mp_init3(t1,t2,g);
  if mp_error<>MP_OKAY then exit;
  if mp_gcd1(a.den, b.den, g) then begin
    {gcd of denominators is 1}
    mp_mul(a.num, b.den, t1);
    mp_mul(b.num, a.den, t2);
    if sub then mp_sub(t1,t2,c.num) else mp_add(t1,t2,c.num);
    if c.num.used=0 then mp_set1(c.den)
    else mp_mul(a.den,b.den,c.den);
  end
  else begin
    {Speedup by calculating second gcd here before calculating c.den. The}
    {result is a product of smaller factors in c.den compared to  first  }
    {multiply greater factors and divide the product during normalization}

    {Formulas:
       g = gcd(a.den, b.den)
       n = a.num*(b.den/g) +/- b.num*(a.den/g)
       h = gcd(n,g)
       c.num = n/h;
       c.den = (a.den/g)*(b.den/h)

     Implementation: h overwrites g; t1,t2 used to accumulate n
       t1 = a.den/g
       t2 = b.den/g
       c.num = a.num*(b.den/g) +/- b.num*(a.den/g) = a.num*t2 +/- b.num*t1
       g = gcd(c.num,g)
       c.num = c.num/g
       c.den = t1*(b.den/h)
    }
    mp_div(a.den, g, t1);
    mp_div(b.den, g, t2);
    mp_mul(a.num, t2, t2);
    mp_mul(b.num, t1, c.num);
    if sub then mp_sub(t2, c.num, c.num) else mp_add(t2, c.num, c.num);
    if mp_gcd1(c.num,g,g) then mp_mul(b.den, t1, c.den)
    else begin
      mp_div(c.num, g, c.num);
      mp_div(b.den, g, t2);
      mp_mul(t1,t2,c.den);
    end;
  end;
  mp_clear3(t1,t2,g);
end;


{---------------------------------------------------------------------------}
function s_mpr_cmp_mag(const a,b: mp_rat): integer;
  {-Compare magnitude of two mp_rats (unsigned), return sign(|a|-|b|), no init check}
var
  dsad,dsbd,dsx,dsy: longint;
  bsad,bsbd,bsx,bsy: longint;
  x,y: mp_int;
begin
  {Value for last alternative, keep D6+ happy}
  s_mpr_cmp_mag := 0;
  if mp_error<>MP_OKAY then exit;

  if a.num.used=0 then begin
    {a is zero; if b<>0 then return -1 else 0}
    if b.num.used<>0 then s_mpr_cmp_mag := -1;
    exit;
  end;

  if b.num.used=0 then begin
    {b is zero, a <> 0, return +1}
    s_mpr_cmp_mag := +1;
    exit;
  end;

  dsad := a.den.used;
  dsbd := b.den.used;

  {$ifdef MPC_USE_Assert}
    {only if assert supported by compiler or debug}
    assert((dsad>0) and (dsbd>0), MPAF+'(a.den.used>0) and (b.den.used>0) in s_mpr_cmp_mag');
  {$endif}

  {In the general case we have to compare x=a.num*b.den with y=b.num*a.den}
  {Before actually calculating the products some checks are made based on}
  {the digit and bit sizes}

  {Check based on digit sizes}
  dsx := a.num.used+dsbd;
  dsy := b.num.used+dsad;
  if dsx > dsy+1 then begin
    s_mpr_cmp_mag := 1;
    exit;
  end;
  if dsy > dsx+1 then begin
    s_mpr_cmp_mag := -1;
    exit;
  end;

  {Check based on bit sizes}
  bsad := mp_bitsize(a.den);
  bsbd := mp_bitsize(b.den);
  bsx  := mp_bitsize(a.num)+bsbd;
  bsy  := mp_bitsize(b.num)+bsad;
  if bsx > bsy+1 then begin
    s_mpr_cmp_mag := 1;
    exit;
  end;
  if bsy > bsx+1 then begin
    s_mpr_cmp_mag := -1;
    exit;
  end;

  {If the denominators have the same bit size check whether they are equal}
  if (bsad=bsbd) and mp_is_eq(a.den, b.den) then begin
    {If the denominators are equal, we just compare the numerators}
    s_mpr_cmp_mag := mp_cmp_mag(a.num, b.num);
    exit;
  end;

  {The hard case: cross multiply and compare}
  mp_init2(x,y);
  if mp_error=MP_OKAY then begin
    {We have to compare x=a.num*b.den with y=b.num*a.den}
    mp_mul(a.num,b.den,x);
    mp_mul(b.num,a.den,y);
    s_mpr_cmp_mag := mp_cmp_mag(x, y);
    mp_clear2(x,y);
  end;

end;


{---------------------------------------------------------------------------}
procedure s_mpr_qr(const num, den: mp_int; var q,r: mp_int);
  {-Calculate q := abs(num) div abs(den); r := abs(num) mod abs(den), no init check}
begin
  mp_abs(num,q);
  mp_abs(den,r);
  mp_divrem(q,r,@q,@r);
end;


{---------------------------------------------------------------------------}
procedure s_mpr_tofloat_n(const a: mp_rat; var str: pchar8; radix, prec: word; maxlen: longint);
  {-Convert to float format for a given radix, prec digits after '.', assumes a is initialized}
var
  n,r,f,p: mp_int;
  intplus: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {check range of radix/maxlen}
  if (radix < 2) or (radix > MAXRadix) or (prec>=maxlen) or (maxlen=0) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('s_mpr_tofloat_n: radix/prec out of range or maxlen too small');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  str^ := #0;

  mp_init4(n,r,f,p);
  if mp_error<>MP_OKAY then exit;
  intplus := mp_show_plus;

  {calculate n,r with |a| = n + r/a.den}
  s_mpr_qr(a.num,a.den,n,r);

  {if Maxlen is too small return empty string}
  if maxlen >= mp_radix_size(n, radix)+prec then begin

    {We have to calculate the frac(a) approximation first in order}
    {to decide if the string starts with '-' for a<0, trunc(a)=0!}

    {Multiply by power of radix: p := r*radix^prec}
    mp_set_pow(p,radix,prec);
    mp_mul(p,r,r);

    if mp_roundfloat and not mp_is1(a.den) then begin
      {Rounding: calculate f=(2*r0*radix^prec + a.den)/(2*a.den)}
      mp_shl(a.den,1,f);
      mp_shl1(r);
      mp_add(r,a.den,r);
      mp_divrem(r,f,@f,nil);
      if mp_is_ge(f,p) then begin
        {Here f=(2*r0*radix^prec+a.den)/(2*a.den) is >= radix^prec, so the}
        {rounding propagates to the left of the radix point! Increment the}
        {integer part, and subtract radix*prec from f; f remains positive.}
        mp_inc(n);
        mp_sub(f,p,f);
      end;
    end
    else begin
      {truncate, no special action}
      mp_divrem(r,a.den,@f,nil);
    end;

    if a.num.sign=MP_NEG then begin
      if n.used=0 then begin
        {if a<0, trunc(a)=0, f<>0 then insert '-'}
        if (maxlen>0) and (f.used<>0) then begin
          str^ := '-';
          inc(str);
          dec(maxlen);
          {do not show '+' for n=0}
          intplus := false;
        end;
      end
      else n.sign := MP_NEG;
    end;
    s_mp_toradix_n(n,radix,intplus,maxlen,str);
    if mp_error<>MP_OKAY then exit;
    if (maxlen>2) and (prec>0) then begin
      str^ := '.';
      inc(str);
      dec(maxlen);
      szf_toradix_n(f,radix,prec,maxlen,str);
    end;
  end;
  mp_clear4(n,r,f,p);
end;


{---------------------------------------------------------------------------}
procedure szf_toradix_n(const a: mp_int; radix,prec,maxlen: word; var str: pchar8);
  {-Convert a positive mp_int to ASCII for a given radix, no init check,}
  { prec digits, leading zeros, on return str points to the final #0 char}
var
  digs: word;
  s0: pchar8;
  t: mp_int;
  d,rp: mp_digit;
  i,ri: integer;
  prmap: ^TRadixCMap;
begin
  if mp_error<>MP_OKAY then exit;
  digs:=1+prec; {prec digits and #0}
  if a.sign=MP_NEG then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('szf_toradix_n: a < 0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  {check range of radix/maxlen, here digs = minimum maxlen}
  if (radix < 2) or (radix > MAXRadix) or (maxlen < digs) then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('szf_toradix_n: radix/prec out of range or maxlen too small');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_init_copy(t,a);
  if mp_error<>MP_OKAY then exit;

  {use local pointer to radix map}
  if mp_uppercase or (radix>36) then prmap := @mp_ucrmap else prmap := @mp_lcrmap;

  {get number of ASCII digits and radix power that fit into an mp_digit}
  ri := mp_radexp[radix];
  rp := mp_radpow[radix];

  {remember first digits position}
  s0 := str;

  {initialize digit counter}
  digs := 0;

  while (mp_error=MP_OKAY) and (prec>0) do begin
    {radix division loop: divide by rp and get ri ASCII digits}
    mp_div_d(t, rp, @t, d);
    {special flag: no trailing '0' in last chunk}
    for i:=1 to ri do begin
      str^ := prmap^[d mod radix];
      inc(str);
      d := d div radix;
      inc(digs);
      dec(prec);
      if prec=0 then break;
    end;
  end;

  if mp_error=MP_OKAY then begin
    {reverse the digits part of the string}
    s_mp_reverse(s0^, digs);
    {append a #0 so the string is properly terminated}
    str^ := #0;
  end;
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpr_add(const a,b: mp_rat; var c: mp_rat);
  {-Add two mp_rats: c = a+b}
begin
  s_mpr_add_sub(a,b,c,false);
end;


{---------------------------------------------------------------------------}
procedure mpr_add_mpi(const a: mp_rat; const b: mp_int; var c: mp_rat);
  {-Add mp_int to mp_rat: c = a+b}
var
  t: mp_rat;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mp_not_init(b) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_add_mpi');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpr_init(t);
  if mp_error=MP_OKAY then begin
    mpr_set1(t,b);
    mpr_add(a,t,c);
    mpr_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_abs(const a: mp_rat; var b: mp_rat);
  {-Absolute value, b = |a|}
begin
  mpr_copy(a, b);
  if mp_error=MP_OKAY then with b do begin
    {force the sign of b to positive}
    num.sign := MP_ZPOS;
    den.sign := MP_ZPOS;
  end;
end;


{$ifdef BIT32or64}
{---------------------------------------------------------------------------}
function mpr_adecimal(const a: mp_rat): ansistring;
  {-Convert to decimal ansistring}
begin
  mpr_adecimal := mpr_radix_astr(a, 10);
end;
{$endif}


{---------------------------------------------------------------------------}
procedure mpr_ceil(const a: mp_rat; var b: mp_int);
  {-Return b := ceil(a)}
var
  r: mp_int;
  asig: integer;
begin
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_ceil');
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
  mp_init(r);
  if mp_error<>MP_OKAY then exit;
  asig := a.num.sign;
  s_mpr_qr(a.num,a.den,b,r);
  if (asig=MP_ZPOS) and (r.used<>0) then mp_inc(b);
  if asig=MP_NEG then b.sign := MP_NEG;
  mp_clear(r);
end;


{---------------------------------------------------------------------------}
function mpr_checksum(const a: mp_rat): longint;
  {-Return a checksum for a, -1 if mp_error<>MP_OKAY, -2 if not initialized}
var
  adler: longint;
begin
  mpr_checksum := -1;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      mpr_checksum := -2;
      exit;
    end;
  {$endif}
  adler := 1;
  with a.num do begin
    s_mp_checksum(adler,@used, sizeof(used));
    s_mp_checksum(adler,@sign, sizeof(sign));
    s_mp_checksum(adler,pdigits, longint(used)*sizeof(mp_digit));
  end;
  with a.den do begin
    s_mp_checksum(adler,@used, sizeof(used));
    s_mp_checksum(adler,@sign, sizeof(sign));
    s_mp_checksum(adler,pdigits, longint(used)*sizeof(mp_digit));
  end;
  mpr_checksum := adler;
end;


{---------------------------------------------------------------------------}
procedure mpr_chs(const a: mp_rat; var b: mp_rat);
  {-Change sign, b = -a}
begin
  if mp_error<>MP_OKAY then exit;
  mpr_copy(a, b);
  if mp_error=MP_OKAY then with b do begin
    if num.used>0 then num.sign := num.sign xor (MP_NEG xor MP_ZPOS);
    s_mpr_normsign(b);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_clear(var a: mp_rat);
  {-Free an mp_rat}
begin
  with a do mp_clear2(num, den);
end;


{---------------------------------------------------------------------------}
procedure mpr_clear2(var a,b: mp_rat);
  {-Clear 2 mp_rats}
begin
  mpr_clear(a);
  mpr_clear(b);
end;


{---------------------------------------------------------------------------}
procedure mpr_clear3(var a,b,c: mp_rat);
  {-Clear 3 mp_rats}
begin
  mpr_clear2(a,b);
  mpr_clear(c);
end;


{---------------------------------------------------------------------------}
procedure mpr_clear4(var a,b,c,d: mp_rat);
  {-Clear 4 mp_rats}
begin
  mpr_clear2(a,b);
  mpr_clear2(c,d);
end;


{---------------------------------------------------------------------------}
procedure mpr_clear5(var a,b,c,d,e: mp_rat);
  {-Clear 5 mp_rats}
begin
  mpr_clear2(a,b);
  mpr_clear2(c,d);
  mpr_clear(e);
end;


{---------------------------------------------------------------------------}
procedure mpr_clear6(var a,b,c,d,e,f: mp_rat);
  {-Clear 6 mp_rats}
begin
  mpr_clear2(a,b);
  mpr_clear2(c,d);
  mpr_clear2(e,f);
end;


{---------------------------------------------------------------------------}
procedure mpr_clear_multi(var vi: array of mp_rat);
  {-Clear a vector of mp_rats}
var
  i: integer;
begin
  for i:=low(vi) to high(vi) do mpr_clear(vi[i]);
end;


{---------------------------------------------------------------------------}
procedure mpr_clear_multi_p(const pv: array of pmp_rat);
  {-Clear a list of mp_rats given as a ptr vector}
var
  i: integer;
begin
  for i:=low(pv) to high(pv) do mpr_clear(pv[i]^);
end;


{---------------------------------------------------------------------------}
function mpr_cmp(const a,b: mp_rat): integer;
  {-Compare two mp_rats (signed), return sign(a-b)}
begin
  {Value for last alternative, keep D6+ happy}
  mpr_cmp := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_cmp');
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
  if a.num.sign<>b.num.sign then begin
    if a.num.sign=MP_NEG then mpr_cmp := -1 else mpr_cmp := 1;
    exit;
  end;
  {compare magnitude}
  if a.num.sign=MP_NEG then begin
    {if negative compare opposite direction}
    mpr_cmp := s_mpr_cmp_mag(b, a);
  end
  else mpr_cmp := s_mpr_cmp_mag(a, b);
end;


{---------------------------------------------------------------------------}
function mpr_cmp_mag(const a,b: mp_rat): integer;
  {-Compare magnitude of two mp_rats (unsigned), return sign(|a|-|b|)}
begin
  {Value for last alternative, keep D6+ happy}
  mpr_cmp_mag := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_cmp_mag');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpr_cmp_mag := s_mpr_cmp_mag(a, b);
end;


{---------------------------------------------------------------------------}
procedure mpr_copy(const a: mp_rat; var b: mp_rat);
  {-Copy an mp_rat, b: = a}
begin
  if mp_error<>MP_OKAY then exit;
  {ArgCheck in mp_copy}
  mp_copy(a.num, b.num);
  mp_copy(a.den, b.den);
end;


{---------------------------------------------------------------------------}
procedure mpr_div(const a,b: mp_rat; var c: mp_rat);
  {-Divide two mp_rats: c = a/b, b<>0}
var
  g1,g2,t1,t2,t3: mp_int;
  is11,is21: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(b) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_div');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if b.num.used=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpr_div: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {if a is zero, so is c}
  if a.num.used=0 then begin
    mpr_zero(c);
    exit;
  end;

  if @a=@b then begin
    {a and b same varibale <>0, so c=1}
    mpr_set_int(c,1,1);
    exit;
  end;

  {(a.num/a.den)/(b.num/d.den) = (a.num*b.den)/(a.den*b.num), so we}
  {use same technique as in mpr_mul, but here are two differences: }
  {1. a third temporary must be used before we can store c.num and }
  {2. c.den can be negative and the results must be sign-normalized}

  mp_init5(g1,g2,t1,t2,t3);
  if mp_error=MP_OKAY then begin
    is11 := mp_gcd1(a.num, b.num, g1);
    is21 := mp_gcd1(a.den, b.den, g2);
    {c.num = (a.num/g1)*(b.den/g2)}
    if is11 then mp_copy(a.num, t1) else mp_div(a.num, g1, t1);
    if is21 then mp_copy(b.den, t2) else mp_div(b.den, g2, t2);
    mp_mul(t1,t2,t3);
    {c.den = (b.num/g1)*(a.den/g2)}
    if is11 then mp_copy(b.num, t1) else mp_div(b.num, g1, t1);
    if is21 then mp_copy(a.den, t2) else mp_div(a.den, g2, t2);
    mp_mul(t1,t2,c.den);
    mp_exch(t3,c.num);
    s_mpr_normsign(c);
    mp_clear5(g1,g2,t1,t2,t3);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_div_2(const a: mp_rat; var b: mp_rat);
  {-Divide mp_rat by 2: b = a/2}
begin
  mpr_copy(a,b);
  if (mp_error<>MP_OKAY) or mpr_is0(b) then exit;
  if mp_iseven(b.num) then mp_shr1(b.num)
  else mp_shl1(b.den);
end;


{---------------------------------------------------------------------------}
procedure mpr_div_2k(const a: mp_rat; k: longint; var b: mp_rat);
  {-Divide mp_rat by 2^k: b = a/2^k}
var
  s: longint;
begin
  if k<0 then mpr_mul_2k(a,abs(k),b)
  else if k=1 then mpr_div_2(a,b)
  else begin
    mpr_copy(a,b);
    if (mp_error<>MP_OKAY) or mpr_is0(b) or (k=0) then exit;
    with b do begin
      if mp_isodd(num) then mp_shl(den,k,den)
      else begin
        s := mp_cnt_lsb(num);
        if s >= k then mp_shr(num,k,num)
        else begin
          mp_shr(num,s,num);
          mp_shl(den,k-s,den);
        end;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_div_int(const a: mp_rat; b: longint; var c: mp_rat);
  {-Divide mp_rat by longint: c = a/b, b<>0}
var
  g: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_div_int');
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
        raise MPXRange.Create('mpr_div_int: b=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  {if a is zero, so is c}
  if a.num.used=0 then begin
    mpr_zero(c);
    exit;
  end;

  g := mp_gcd_int(a.num,b);
  if g=1 then begin
    mp_copy(a.num, c.num);
    mp_mul_int(a.den, b, c.den);
  end
  else begin
    mp_mul_int(a.den, b div g, c.den);
    mp_div_int(a.num, g, @c.num, b);
  end;
  s_mpr_normsign(c);
end;


{---------------------------------------------------------------------------}
procedure mpr_div_mpi(const a: mp_rat; const b: mp_int; var c: mp_rat);
  {-Divide mp_rat by mp_int: c = a/b, b<>0}
var
  t: mp_rat;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mp_not_init(b) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_div_mpi');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpr_init(t);
  if mp_error=MP_OKAY then begin
    mpr_set1(t,b);
    mpr_div(a,t,c);
    mpr_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
function mpr_decimal(const a: mp_rat): mp_string;
  {-Convert r to decimal, max 255 chars}
begin
  mpr_decimal := mpr_radix_str(a, 10);
end;


{---------------------------------------------------------------------------}
procedure mpr_exch(var a,b: mp_rat);
  {-Exchange two mp_rats}
var
  t: mp_rat;
begin
  if mp_error<>MP_OKAY then exit;
  t := a;
  a := b;
  b := t;
end;


{---------------------------------------------------------------------------}
procedure mpr_expt(const a: mp_rat; b: longint; var c: mp_rat);
  {-Calculate c = a^b, a<>0 for b<0, 0^0=1}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_expt');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {Handle a=0}
  if a.den.used=0 then begin
    if b<0 then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXRange.Create('mpr_expt: a=0, b<0');
        {$else}
          RunError(MP_RTE_RANGE);
        {$endif}
      {$else}
        set_mp_error(MP_RANGE);
      {$endif}
    end
    else if b>0 then mpr_zero(c)
    else mpr_set_int(c,1,1);
    exit;
  end;

  {Setup for exponentiation and handle easy cases}
  if b<0 then begin
    mpr_inv(a,c);
    b := abs(b);
  end
  else begin
    if b=0 then mpr_set_int(c,1,1) else mpr_copy(a,c);
  end;
  if b<=1 then exit;

  {here we have to calculate c^b, b>1. Since rats are always stored}
  {in lowest terms, calculate c.num^b/c.den^b without normalization}
  mp_expt_int(c.num, b, c.num);
  mp_expt_int(c.den, b, c.den);
end;


{---------------------------------------------------------------------------}
procedure mpr_floor(const a: mp_rat; var b: mp_int);
  {-Return b := floor(a)}
var
  r: mp_int;
  asig: integer;
begin
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_floor');
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
  mp_init(r);
  if mp_error<>MP_OKAY then exit;
  asig := a.num.sign;
  s_mpr_qr(a.num,a.den,b,r);
  if (asig=MP_NEG) and (b.used<>0) then b.sign := MP_NEG;
  if (asig=MP_NEG) and (r.used<>0) then mp_dec(b);
  mp_clear(r);
end;


{---------------------------------------------------------------------------}
procedure mpr_frac(const a: mp_rat; var b: mp_rat);
  {-Return b := frac(a) = a - trunc(a)}
var
  t: mp_int;
begin
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_frac');
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
  mp_init(t);
  if mp_error<>MP_OKAY then exit;
  mpr_trunc(a,t);
  mpr_sub_mpi(a,t,b);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mpr_harmonic(n: longint; var hn: mp_rat);
  {-Compute the harmonic number hn = sum(1/i, i=1..n) with binary splitting}

  procedure bs_harm(a,b: longint; var r: mp_rat);
    {-Binary splitting calculation of sum(1/i, a<=i<b)}
  var
    l: mp_rat;
    m: longint;
  begin
    {See Alg. 4 in Fredrik Johansson's 'How (not) to compute harmonic numbers'}
    {http://fredrik-j.blogspot.com/2009/02/how-not-to-compute-harmonic-numbers.html}
    if b-a=1 then begin
      mpr_set_int(r,1,a);
    end
    else begin
      mpr_init(l);
      if mp_error=MP_OKAY then begin
        m := (a+b) div 2;
        bs_harm(a,m,l);
        bs_harm(m,b,r);
        mpr_add(l,r,r);
        mpr_clear(l);
      end;
    end;
  end;

begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(hn) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_harmonic');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if n<1 then mpr_zero(hn)
  else bs_harm(1,n+1,hn);
end;

{---------------------------------------------------------------------------}
procedure mpr_init(var a: mp_rat);
  {-Initialize an mp_rat}
begin
  mpr_init_size(a,0,0);
end;


{---------------------------------------------------------------------------}
procedure mpr_init2(var a,b: mp_rat);
  {-Initialize 2 mp_rats}
var
  pa: array[0..1] of pmp_rat;
begin
  pa[0] := @a;
  pa[1] := @b;
  mpr_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mpr_init3(var a,b,c: mp_rat);
  {-Initialize 3 mp_rats}
var
  pa: array[0..2] of pmp_rat;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  mpr_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mpr_init4(var a,b,c,d: mp_rat);
  {-Initialize 4 mp_rats}
var
  pa: array[0..3] of pmp_rat;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  mpr_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mpr_init5(var a,b,c,d,e: mp_rat);
  {-Initialize 5 mp_rats}
var
  pa: array[0..4] of pmp_rat;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  pa[4] := @e;
  mpr_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mpr_init6(var a,b,c,d,e,f: mp_rat);
  {-Initialize 6 mp_rats}
var
  pa: array[0..5] of pmp_rat;
begin
  pa[0] := @a;
  pa[1] := @b;
  pa[2] := @c;
  pa[3] := @d;
  pa[4] := @e;
  pa[5] := @f;
  mpr_init_multi_p(pa);
end;


{---------------------------------------------------------------------------}
procedure mpr_init_copy(var a: mp_rat; const b: mp_rat);
  {-Create a, then copy b into it}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_init_copy');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_init_copy(a.num, b.num);  if mp_error<>MP_OKAY then exit;
  mp_init_copy(a.den, b.den);  if mp_error<>MP_OKAY then mp_clear(a.num);
end;


{---------------------------------------------------------------------------}
procedure mpr_init_multi(var vi: array of mp_rat);
  {-Initialize a vector of mp_rats}
  { on error the already initialized mp_rats will be cleared}
var
  i,k: integer;
begin
  if mp_error<>MP_OKAY then exit;
  for i:=low(vi) to high(vi) do begin
    mpr_init(vi[i]);
    if mp_error<>MP_OKAY then begin
      {error, clear all previous mp_rats}
      for k:=low(vi) to i-1 do mpr_clear(vi[k]);
      break;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_init_multi_p(var pv: array of pmp_rat);
  {-Initialize a list of mp_rats given as a ptr vector}
  { on error the already initialized mp_rats will be cleared}
var
  i,k: integer;
begin
  if mp_error<>MP_OKAY then exit;
  for i:=low(pv) to high(pv) do begin
    mpr_init(pv[i]^);
    if mp_error<>MP_OKAY then begin
      {error, clear all previous mp_rats}
      for k:=low(pv) to i-1 do mpr_clear(pv[k]^);
      break;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_init_size(var a: mp_rat; nsize, dsize: longint);
  {-Initialize an mp_rat to given number of digits}
begin
  if mp_error<>MP_OKAY then exit;
  with a do begin
    mp_init_size(num, nsize);
    if mp_error=MP_OKAY then begin
      mp_init_size(den, dsize);
      if mp_error=MP_OKAY then begin
        mp_set1(den);
        if mp_error<>MP_OKAY then mpr_clear(a);
      end
      else mp_clear(num);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_inv(const a: mp_rat; var b: mp_rat);
  {-Invert an mp_rat, b = 1/a, a<>0}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_inv');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if a.den.used=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpr_inv: a=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  if @a=@b then begin
    {a and b same variable, exchange numerator and denominator}
    mp_exch(b.num,b.den);
  end
  else begin
    mp_copy(a.num, b.den);
    mp_copy(a.den, b.num);
  end;
  {a should be normalized, so adjust only the signs of b}
  s_mpr_normsign(b);
end;


{---------------------------------------------------------------------------}
function mpr_is0(const a: mp_rat): boolean;
  {-Return a=0}
begin
  mpr_is0 := (a.num.magic=MP_MAGIC) and (a.num.pdigits<>nil) and (a.num.used=0)
         and (a.den.magic=MP_MAGIC) and (a.den.pdigits<>nil) and (a.den.used>0);
end;


{---------------------------------------------------------------------------}
function mpr_is_eq(const a,b: mp_rat): boolean;
  {-Return a = b}
begin
  {No need to use mpr_cmp, normalized rats are equal if the numerators and}
  {denominators are equal, short cuts for a<>b are done with xx.used tests}
  mpr_is_eq := (a.num.used=b.num.used) and (a.den.used=b.den.used)
               and mp_is_eq(a.num, b.num) and mp_is_eq(a.den, b.den)
end;


{---------------------------------------------------------------------------}
function mpr_is_ge(const a,b: mp_rat): boolean;
  {-Return a >= b}
begin
  mpr_is_ge := mpr_is_eq(a,b) or mpr_is_gt(a,b);
end;


{---------------------------------------------------------------------------}
function mpr_is_gt(const a,b: mp_rat): boolean;
  {-Return a > b}
begin
  mpr_is_gt := mpr_cmp(a,b)=MP_GT;
end;


{---------------------------------------------------------------------------}
function mpr_is_le(const a,b: mp_rat): boolean;
  {-Return a <= b}
begin
  mpr_is_le := mpr_is_eq(a,b) or mpr_is_lt(a,b);
end;


{---------------------------------------------------------------------------}
function mpr_is_lt(const a,b: mp_rat): boolean;
  {-Return a < b}
begin
  mpr_is_lt := mpr_cmp(a,b)=MP_LT;
end;


{---------------------------------------------------------------------------}
function mpr_is_mp(const a: mp_rat): boolean;
  {-Return true if a is initialized and a.den=1}
begin
  mpr_is_mp := (a.num.magic=MP_MAGIC) and (a.num.pdigits<>nil) and mp_is1(a.den);
end;


{---------------------------------------------------------------------------}
function mpr_is_ne(const a,b: mp_rat): boolean;
  {-Return a <> b}
begin
  mpr_is_ne := mp_is_ne(a.num, b.num) or mp_is_ne(a.den, b.den);
end;


{---------------------------------------------------------------------------}
procedure mpr_mul(const a,b: mp_rat; var c: mp_rat);
  {-Multiply two mp_rats: c = a*b}
var
  g1,g2,t1,t2: mp_int;
  is11,is21: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(b) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_mul');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  {if a or b is zero, so is z}
  if (a.num.used=0) or (b.num.used=0) then begin
    mpr_zero(c);
    exit;
  end;

  {if a==b just do two squares}
  if @a=@b then begin
    mp_sqr(a.num, c.num);
    mp_sqr(a.den, c.den);
    exit;
  end;

  mp_init4(g1,g2,t1,t2);
  if mp_error=MP_OKAY then begin
    is11 := mp_gcd1(a.num, b.den, g1);
    is21 := mp_gcd1(a.den, b.num, g2);
    {c.num = (a.num/g1)*(b.num/g2)}
    if is11 then mp_copy(a.num, t1) else mp_div(a.num, g1, t1);
    if is21 then mp_copy(b.num, t2) else mp_div(b.num, g2, t2);
    mp_mul(t1,t2,c.num);
    {c.den = (b.den/g1)*(a.den/g2)}
    if is11 then mp_copy(b.den, t1) else mp_div(b.den, g1, t1);
    if is21 then mp_copy(a.den, t2) else mp_div(a.den, g2, t2);
    mp_mul(t1,t2,c.den);
    mp_clear4(g1,g2,t1,t2);
  end;

end;


{---------------------------------------------------------------------------}
procedure mpr_mul_2(const a: mp_rat; var b: mp_rat);
  {-Multiply mp_rat by 2: b = 2*a}
begin
  mpr_copy(a,b);
  if (mp_error<>MP_OKAY) or mpr_is0(b) then exit;
  if mp_iseven(b.den) then mp_shr1(b.den)
  else mp_shl1(b.num);
end;


{---------------------------------------------------------------------------}
procedure mpr_mul_2k(const a: mp_rat; k: longint; var b: mp_rat);
  {-Multiply mp_rat by 2^k: b = a*2^k}
var
  s: longint;
begin
  if k<0 then mpr_div_2k(a,abs(k),b)
  else if k=1 then mpr_mul_2(a,b)
  else begin
    mpr_copy(a,b);
    if (mp_error<>MP_OKAY) or mpr_is0(b) or (k=0) then exit;
    with b do begin
      if mp_isodd(den) then mp_shl(num,k,num)
      else begin
        s := mp_cnt_lsb(den);
        if s >= k then mp_shr(den,k,den)
        else begin
          mp_shl(num,k-s,num);
          mp_shr(den,s,den);
        end;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_mul_int(const a: mp_rat; b: longint; var c: mp_rat);
  {-Multiply mp_rat by longint: c = a*b}
var
  g: longint;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_mul_int');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if (b=0) or (a.num.used=0) then begin
    mpr_zero(c);
    exit;
  end
  else if abs(b)=1 then begin
    if b>0 then mpr_copy(a,c)
    else mpr_chs(a,c);
    exit;
  end;

  g := mp_gcd_int(a.den,b);
  if g=1 then begin
    mp_mul_int(a.num, b, c.num);
    mp_copy(a.den, c.den);
  end
  else begin
    mp_mul_int(a.num, b div g, c.num);
    mp_div_int(a.den, g, @c.den, b);
  end;
  s_mpr_normsign(c);
end;


{---------------------------------------------------------------------------}
procedure mpr_mul_mpi(const a: mp_rat; const b: mp_int; var c: mp_rat);
  {-Multiply mp_rat and mp_int: c = a*b}
var
  t: mp_rat;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mp_not_init(b) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_mul_mpi');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpr_init(t);
  if mp_error=MP_OKAY then begin
    mpr_set1(t,b);
    mpr_mul(a,t,c);
    mpr_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_output_decimal(const a: mp_rat);
  {-Write decimal representation to output}
begin
  mpr_write_decimal(output, a);
end;


{---------------------------------------------------------------------------}
procedure mpr_output_radix(const a: mp_rat; radix: word);
  {-Write radix representation to output}
begin
  mpr_write_radix(output, a, radix);
end;


{---------------------------------------------------------------------------}
function mpr_not_init(const a: mp_rat): boolean;
  {-Sanity check if a is initialized, does not catch all cases!}
begin
  mpr_not_init := mp_not_init(a.num) or mp_not_init(a.den);
end;


{$ifdef BIT32or64}
{---------------------------------------------------------------------------}
function mpr_radix_astr(const a: mp_rat; radix: word): ansistring;
  {-Convert to radix representation ansistring}
begin
  mpr_radix_astr := '';
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_write_radix');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  if mp_is1(a.den) then mpr_radix_astr := mp_radix_astr(a.num,radix)
  else mpr_radix_astr := mp_radix_astr(a.num,radix) + '/' + s_mp_radix_astr(a.den,radix,false);
end;
{$endif}


{---------------------------------------------------------------------------}
function mpr_radix_size(const a: mp_rat; radix: word): longint;
  {-Return size of ASCII representation (incl. sign and #0)}
var
  digs: longint;
begin
  mpr_radix_size := 0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_radix_size');
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
        raise MPXRange.Create('mpr_radix_size: radix out of range');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  digs := mp_radix_size(a.num, radix);
  if not mp_is1(a.den) then begin
    {add size of den, counts '/' and #0 cancel}
    digs := digs + mp_radix_size(a.den, radix);
    {den is positive but no '+' is generated}
    if mp_show_plus then dec(digs);
  end;
  mpr_radix_size := digs;
end;


{---------------------------------------------------------------------------}
function mpr_radix_str(const a: mp_rat; radix: word): mp_string;
  {-Convert to radix representation, max 255 digits}
var
  i: integer;
  ls: longint;
  iostr: string[255];
begin
  {arg checks are done by mp_radix_size}
  mpr_radix_str := '';
  ls := mpr_radix_size(a, radix);
  if mp_error<>MP_OKAY then exit;
  if ls<=255 then mpr_toradix_n(a, @iostr[1], radix, 255);
  if (ls>255) or (mp_error<>MP_OKAY) then exit;
  iostr[0] := #255;
  for i:=1 to 255 do begin
    if iostr[i]=#0 then begin
      iostr[0] := char8(i-1);
      break;
    end;
  end;
  mpr_radix_str := iostr;
end;


{---------------------------------------------------------------------------}
procedure mpr_read_decimal(var a: mp_rat; str: pchar8);
  {-Read a ASCII decimal string into a. str may contain a single '/'}
begin
  {the denominator can have sign characters}
  mpr_read_radix(a,str,10);
end;


{---------------------------------------------------------------------------}
procedure mpr_read_double(var a: mp_rat; d: double);
  {-Convert d to an mp_rat}
var
  m: double;
  e: longint;
  i: integer;
  dig: mp_digit;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {mpr_zero does ArgCheck!}
  mpr_zero(a);
  if d=0.0 then exit;
  neg := d<0.0;
  d := abs(d);
  frexpd(d,m,e);
  for i:=0 to succ(54 div DIGIT_BIT) do begin
    if m=0.0 then break;
    mp_lshd(a.num,1);
    m := ldexpd(m,DIGIT_BIT);
    dig := trunc(m) and MP_MASK;
    {This should reduce the number of mantissa bits by DIGIT_BIT, the}
    {for loop is used to avoid an infinite loop due to rounding errors}
    m   := m - dig;
    dec(e,DIGIT_BIT);
    mp_add_d(a.num,dig,a.num);
  end;
  if e<0 then mp_shl(a.den,-e,a.den)
  else mp_shl(a.num,e,a.num);
  s_mpr_normalize(a);
  if neg then mpr_chs(a,a);
end;


{---------------------------------------------------------------------------}
procedure mpr_read_float_decimal(var a: mp_rat; str: pchar8);
  {-Read a ASCII float decimal string into a. str may contain a single '.'}
begin
  mpr_read_float_radix(a, str, 10);
end;


{---------------------------------------------------------------------------}
procedure mpr_read_float_radix(var a: mp_rat; str: pchar8; radix: word);
  {-Read a ASCII float radix string into a. str may contain a single '.'}
var
  fs: pchar8;
  t: mp_int;
  c: char8;
  neg: boolean;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_read_float_radix');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_init(t);
  if mp_error<>MP_OKAY then exit;


  {skip leading white space}
  repeat
    c := upcase(str^);
    if c=#0 then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mpr_read_float_radix: invalid syntax');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        mp_clear(t);
        set_mp_error(MP_BADARG);
        exit;
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
  if c=mp_fract_sep{'.'} then mp_zero(t)
  else begin
    {read integer part until '.' or #0}
    s_mp_read_radix(t,str,radix,mp_fract_sep{'.'},false);
  end;
  if str^=mp_fract_sep{'.'} then begin
    {read floating part into numerator}
    inc(str);
    {remember starting position of floating part}
    fs := str;
    s_mp_read_radix(a.num,str,radix,#0,false);
    if str<>fs then begin
      {set denominator to radix^(number of floating digits)}
      mp_set_pow(a.den,radix,str-fs);
      {normalize and add integer part}
      s_mpr_normalize(a);
      mpr_add_mpi(a,t,a);
    end
    else mpr_set1(a,t);
  end
  else mpr_set1(a,t);
  if neg then mpr_chs(a,a);
  mp_clear(t);
end;



{---------------------------------------------------------------------------}
procedure mpr_read_radix(var a: mp_rat; str: pchar8; radix: word);
  {-Read a ASCII radix string into a. str may contain a single '/'}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_read_radix');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mp_read_radix(a.num,str,radix,'/',true);
  if str^='/' then begin
    inc(str);
    s_mp_read_radix(a.den,str,radix,#0,false);
  end
  else mp_set1(a.den);
  s_mpr_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpr_round(const a: mp_rat; var b: mp_int);
  {-Return b := round(a), round(-a)=-round(a), round(n+0.5)=n+1}
var
  n,d: mp_int;
  asig: integer;
begin
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_round');
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
  mp_init2(n,d);
  if mp_error<>MP_OKAY then exit;
  asig := a.num.sign;
  {round(n/d) = trunc((2|n|+d)/(2d))}
  mp_abs(a.num,n);
  mp_shl1(n);
  mp_add(n,a.den,n);
  mp_shl(a.den,1,d);
  s_mpr_qr(n,d,b,d);
  if asig=MP_NEG then b.sign := MP_NEG;
  mp_clear2(n,d);
end;


{---------------------------------------------------------------------------}
procedure mpr_set(var a: mp_rat; const n, d: mp_int);
  {-Set a to n/d, d<>0}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mp_not_init(n) or mp_not_init(d) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_set');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
    if (@a.num=@n) or (@a.num=@d) or (@a.den=@n) or (@a.den=@d) then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mpr_set: invalid addr');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        exit;
      {$endif}
    end;
  {$endif}
  if d.used=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpr_set: d=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mp_copy(n, a.num);
  mp_copy(d, a.den);
  if mp_error<>MP_OKAY then exit;
  s_mpr_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpr_set_int(var a: mp_rat; n, d: longint);
  {-Set a to n/d, d<>0}
begin
  if mp_error<>MP_OKAY then exit;
  {mp_set_int does MPC_ArgCheck}
  {cannot divide by zero}
  if d=0 then begin
    {$ifdef MPC_HaltOnError}
      {$ifdef MPC_UseExceptions}
        raise MPXRange.Create('mpr_set_int: d=0');
      {$else}
        RunError(MP_RTE_RANGE);
      {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;
  mp_set_int(a.num, n);
  mp_set_int(a.den, d);
  s_mpr_normalize(a);
end;


{---------------------------------------------------------------------------}
procedure mpr_set1(var a: mp_rat; const n: mp_int);
  {-Set a to n/1}
begin
  if mp_error<>MP_OKAY then exit;
  {Argcheck by mp_copy and mp_set}
  mp_copy(n, a.num);
  mp_set(a.den, 1);
  {no need to call normalize}
end;


{---------------------------------------------------------------------------}
procedure mpr_sub(const a,b: mp_rat; var c: mp_rat);
  {-Subtract two mp_rats: c = a-b}
begin
  s_mpr_add_sub(a,b,c,true);
end;


{---------------------------------------------------------------------------}
procedure mpr_sub_mpi(const a: mp_rat; const b: mp_int; var c: mp_rat);
  {-Subtract mp_int from mp_rat: c = a-b}
var
  t: mp_rat;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mp_not_init(b) or mpr_not_init(c) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_sub_mpi');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mpr_init(t);
  if mp_error=MP_OKAY then begin
    mpr_set1(t,b);
    mpr_sub(a,t,c);
    mpr_clear(t);
  end;
end;


{---------------------------------------------------------------------------}
function mpr_todouble(const a: mp_rat): double;
  {-Convert a to double, +-inf if too large, 0 if too small}
var
  bsn,bsd,scn,scd: longint; {bitsizes and shift counts}
  an,ad: mp_int;
begin
  mpr_todouble := 0.0;
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_todouble');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}

  if mpr_is0(a) then exit
  else if mpr_is_mp(a) then begin
    {a is an mp_int use mp_todouble}
    mpr_todouble := mp_todouble(a.num);
    exit;
  end;

  mp_init2(an,ad);
  if mp_error=MP_OKAY then begin
    {get approximation for log2(a)}
    bsn := mp_bitsize(a.num);
    bsd := mp_bitsize(a.den);
    if bsn-bsd>1024 then begin
      {quotient too large return +-inf}
      if a.num.sign=MP_NEG then mpr_todouble := DblNegInf
      else mpr_todouble := DblPosInf;
    end
    else if bsn-bsd < -1024 then begin
      {quotient too small return zero}
      mpr_todouble := 0.0;
    end
    else begin
      {shift num and den into range 0..2^64}
      if bsn>64 then scn := bsn-64 else scn := 0;
      if bsd>64 then scd := bsd-64 else scd := 0;
      mp_shr(a.num, scn, an);
      mp_shr(a.den, scd, ad);
      {result := ((a.num*2^-scn)/(a.den*2^-scd))*2^(scn-scd)}
      mpr_todouble := ldexpd(mp_todouble(an)/mp_todouble(ad),scn - scd);
    end;
    mp_clear2(an,ad);
  end;
end;


{$ifdef BIT32or64}
{---------------------------------------------------------------------------}
function mpr_tofloat_astr(const a: mp_rat; radix, prec: word): ansistring;
  {-Convert to float representation with prec, max 65000 digits}
const
  LMAX=65000;
var
  l: longint;
  s0,s1: pchar8;
  {$ifndef RESULT}
    Result: ansistring;
  {$endif}
begin
  mpr_tofloat_astr := '';
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_tofloat_astr');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  SetLength(Result, LMAX);
  s0 := pchar8(@Result[1]);
  s1 := s0;
  s_mpr_tofloat_n(a, s1, radix, prec, LMAX);
  if mp_error=MP_OKAY then begin
    l := (s1-s0)+1;
    {trim trailing #0}
    while (l>0) and (Result[l]=#0) do dec(l);
    SetLength(Result, l);
    {$ifndef RESULT}
      mpr_tofloat_astr := Result;
    {$endif}
  end
  else mpr_tofloat_astr := '';
end;
{$endif}


{---------------------------------------------------------------------------}
procedure mpr_tofloat_n(const a: mp_rat; str: pchar8; radix, prec, maxlen: word);
  {-Convert to float format for a given radix, prec digits after '.'}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_tofloat_n');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mpr_tofloat_n(a, str, radix, prec, maxlen);
end;


{---------------------------------------------------------------------------}
function mpr_tofloat_str(const a: mp_rat; radix, prec: word): mp_string;
  {-Convert to float representation with prec, max 255 digits}
var
  i: integer;
  iostr: string[255];
begin
  {arg checks are done by mpr_tofloat_n}
  mpr_tofloat_str := '';
  if mp_error<>MP_OKAY then exit;
  mpr_tofloat_n(a, @iostr[1], radix, prec, 255);
  if mp_error<>MP_OKAY then exit;
  iostr[0] := #255;
  for i:=1 to 255 do begin
    if iostr[i]=#0 then begin
      iostr[0] := char8(i-1);
      break;
    end;
  end;
  mpr_tofloat_str := iostr;
end;


{---------------------------------------------------------------------------}
procedure mpr_toradix_n(const a: mp_rat; str: pchar8; radix: word; maxlen: longint);
  {-Convert an mp_rat to an ASCII string for a given radix (2..MAXRadix)}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_toradix_n');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  s_mp_toradix_n(a.num,radix,mp_show_plus,maxlen,str);
  if mp_error<>MP_OKAY then exit;
  if (not mp_is1(a.den)) and (maxlen>2) then begin
    str^ := '/';
    inc(str);
    dec(maxlen);
    s_mp_toradix_n(a.den,radix,false,maxlen,str);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_trunc(const a: mp_rat; var b: mp_int);
  {-Return b := trunc(a)}
var
  r: mp_int;
  asig: integer;
begin
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) or mp_not_init(b) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_trunc');
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
  mp_init(r);
  if mp_error<>MP_OKAY then exit;
  asig := a.num.sign;
  s_mpr_qr(a.num,a.den,b,r);
  if asig=MP_NEG then b.sign := MP_NEG;
  mp_clear(r);
end;


{---------------------------------------------------------------------------}
procedure mpr_write_decimal(var tf: system.text; const a: mp_rat);
  {-Write decimal representation to file tf}
begin
  mpr_write_radix(tf, a, 10);
end;


{---------------------------------------------------------------------------}
procedure mpr_write_radix(var tf: system.text; const a: mp_rat; radix: word);
  {-Write radix representation to file tf}
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mpr_not_init(a) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mpr_write_radix');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  mp_write_radix(tf,a.num,radix);
  if not mp_is1(a.den) then begin
    write(tf,'/');
    s_mp_write_radix(tf,a.den,radix,false);
  end;
end;


{---------------------------------------------------------------------------}
procedure mpr_zero(var a: mp_rat);
  {-Set a to zero}
begin
  {arg checks are done by mpr_zero/set}
  if mp_error<>MP_OKAY then exit;
  mp_zero(a.num);
  mp_set1(a.den);
end;


{---------------------------------------------------------------------------}
{---------------------------------------------------------------------------}

end.

