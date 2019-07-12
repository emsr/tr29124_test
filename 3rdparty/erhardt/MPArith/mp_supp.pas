unit mp_supp;

{Supplemental routines for MP integers}

interface

{$i STD.INC}

uses
  mp_types;

{$i mp_conf.inc}

(*************************************************************************

 DESCRIPTION     :  Supplemental routines for MP integers

 REQUIREMENTS    :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA   :  (mp_types)

 MEMORY USAGE    :  heap

 DISPLAY MODE    :  ---

 REFERENCES      :  [2] MPI, M.J. Fromberger


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
          Aug 2004  W.Ehrhardt  Initial raw version
 0.3.01   26.08.04  we          with mp_... routines
 0.3.02   14.01.05  we          with mp_arctan
 0.3.03   28.03.05  we          with {$ifdef BIT16} for {$x+}
 0.3.04   10.08.05  we          Bugfix memory leak in mp_write_radix
 0.3.05   17.08.05  we          mp_dumpu/a routines

 0.4.00   20.08.05  we          use mp_set_error
 0.4.01   20.08.05  we          mp_output_radix, Radix: word in mp_arctan
 0.4.02   21.08.05  we          mp_conf.inc, MPC_Diagnostic: mp_memused from mp_base
 0.4.03   21.08.05  we          MPC_ArgCheck
 0.4.04   22.08.05  we          removed mp_argcheck, mp_errchk
 0.4.04   22.08.05  we          mp_dump_memused
 0.4.05   24.08.05  we          improved arg checks
 0.4.06   24.08.05  we          ansistring: mp_radix_astr, mp_adecimal
 0.4.07   27.08.05  we          improved mp_arctan
 0.4.08   28.08.05  we          usage of exceptions implemented
 0.4.09   08.09.05  we          some functions moved to mp_base

 0.7.00   19.03.06  we          mp_arctanw
 0.7.01   11.08.06  we          Avoid FPC warnings: mp_arctanw
 0.7.02   12.08.06  we          mp_dump_diagctr: dump diagnostic counters

 1.1.00   27.06.07  we          mp_reset_diagctr

 1.3.00   22.12.07  we          mpf_dump_me

 1.4.00   05.01.08  we          mp_atanhw

 1.7.00   24.09.08  we          string replaced by mp_string
 1.7.01   26.09.08  we          DIGIT_BIT dependent formatting for mp_dumpa/u

 1.10.00  21.01.09  we          changes related to (s)mp_divrem

 1.24.00  03.01.13  we          Fix some words to TNInt
 1.24.01  04.01.13  we          Renamed mp_expt_d to mp_expt_w

**************************************************************************)

(*-------------------------------------------------------------------------
  This code uses material/ideas from the following 3rd party libraries:
   - MPI 1.8.6 by Michael J. Fromberger
  See the file '3rdparty.mpa' for the licenses.
----------------------------------------------------------------------------*)


(*-------------------------------------------------------------------------
 (C) Copyright 2004-2013 Wolfgang Ehrhardt

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


procedure mp_arctan(mul, x, Radix: mp_digit; prec: word; var sum: mp_int);
  {-Compute sum := mul * arctan(1/x), to prec Radix digits}

procedure mp_arctanw(mul, x, Radix, prec: word; var sum: mp_int);
  {-Compute sum := mul*arctan(1/x), to prec Radix digits, word version.}
  { This version is less general for 32 bit mp_digits but extends}
  { the ranges of prec and x for 16 bit mp_digits and allows one }
  { Machin formula with good performance for all mp_digits sizes.}

procedure mp_atanhw(mul, x, Radix, prec: word; var sum: mp_int);
  {-Compute sum := mul*atanh(1/x), to prec Radix digits}

procedure mp_dumpu(const a: mp_int; const hdr: mp_string);
  {-Dump used fields of a}

procedure mp_dumpa(const a: mp_int; const hdr: mp_string);
  {-Dump all fields of a}

procedure mp_dump_diagctr;
  {-Dump diagnostic counters}

procedure mp_reset_diagctr;
  {-Reset diagnostic counters to 0}

procedure mp_dump_meminfo;
  {-Write mp_memstat to output (if MPC_Diagnostic is defined)}

procedure mp_dump_memused;
  {-Write mp_memused to output (if MPC_Diagnostic is defined)}

function  mp_memused: longint;
  {-Return total allocated memory, MaxLongint if MPC_Diagnostic not defined}

procedure mpf_dump_me(const a: mp_float);
  {-Dump mantissa and exponent of a}


implementation

uses
  mp_base;


{---------------------------------------------------------------------------}
procedure mp_arctan(mul, x, Radix: mp_digit; prec: word; var sum: mp_int);
  {-Compute sum := mul * arctan(1/x), to prec Radix digits}
var
  t,v:  mp_int;
  q,rd: mp_digit;
  sign: integer;
  DMax: mp_word;
const
  EXTRA = 10;
begin
  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(sum) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_arctan');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {Algorithm from MPI/Utils pi.c}
  sign := 1;
  q := 1;
  DMax := MP_DIGIT_MAX;
  {Check x^2 <= MP_DIGIT_MAX, Radix <= MaxRadix}
  if (sqr(mp_word(x)) > DMax) or (Radix>MAXRadix) then begin
    {$ifdef MPC_HaltOnError}
     {$ifdef MPC_UseExceptions}
       raise MPXRange.Create('mp_arctan: x too large');
     {$else}
       RunError(MP_RTE_RANGE);
     {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_init(t);
  mp_set(t, Radix);
  mp_init(v);

  {Increase precision to guard against rounding errors}
  mp_expt_w(t, prec + EXTRA, t);
  mp_mul_d(t,mul,t);
  mp_mul_d(t,x,t);

  {                1     1       1       1         }
  { arctan(1/x) =  - - ----- + ----- - ----- + ... }
  {                x   3 x^3   5 x^5   7 x^7       }

  x := x*x;
  mp_zero(sum);

  {prec must be <= trunc(MP_DIGIT_MAX*ln(x)/ln(Radix))-EXTRA}
  {otherwise the q will become greater than MP_DIGIT_MAX    }

  repeat
    mp_div_d(t, x, @t, rd);
    if (sign<0) and (rd<>0) then mp_inc(t);

    mp_div_d(t, q, @v, rd);
    if (sign<0) and (rd<>0) then mp_inc(v);

    if sign>0 then mp_add(sum, v, sum) else mp_sub(sum, v, sum);

    sign := -sign;
    q := q+2;
    if q>MP_DIGIT_MAX then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_arctan: prec too large');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        break;
      {$endif}
    end;
  until mp_iszero(t);
  {Chop off extra precision}
  {sum := trunc(sum/Radix^EXTRA + 0.5}
  mp_set(v,Radix);
  mp_expt_w(v,EXTRA,v);
  mp_div_2(v,t);
  mp_add(sum,t,sum);
  mp_div(sum, v, sum);

  mp_clear(v);
  mp_clear(t);
end;


{---------------------------------------------------------------------------}
procedure mp_arctanw(mul, x, Radix, prec: word; var sum: mp_int);
  {-Compute sum := mul*arctan(1/x), to prec Radix digits, word version.}
  { This version is less general for 32 bit mp_digits but extends}
  { the ranges of prec and x for 16 bit mp_digits and allows one }
  { Machin formula with good performance for all mp_digits sizes.}
var
  t,v:  mp_int;
  q,rd: word;
  sign: integer;
const
  EXTRA = 10;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(sum) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_arctanw');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {Algorithm from MPI/Utils pi.c}
  sign := 1;
  q := 1;
  {Check x^2 <= $FFFF, Radix <= MaxRadix}
  if (x > $FF) or (Radix>MAXRadix) then begin
    {$ifdef MPC_HaltOnError}
     {$ifdef MPC_UseExceptions}
       raise MPXRange.Create('mp_arctanw: x or Radix too large');
     {$else}
       RunError(MP_RTE_RANGE);
     {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_init2(t,v);
  if mp_error<>MP_OKAY then exit;

  {Increase precision to guard against rounding errors}
  mp_set_w(t, Radix);
  mp_expt_w(t, prec + EXTRA, t);
  mp_mul_w(t,mul,t);
  mp_mul_w(t,x,t);

  {                1     1       1       1         }
  { arctan(1/x) =  - - ----- + ----- - ----- + ... }
  {                x   3 x^3   5 x^5   7 x^7       }

  x := x*x;
  mp_zero(sum);

  {prec must be <= trunc(2^16*ln(x)/ln(Radix))-EXTRA}
  {otherwise the q will become greater than 2^16    }

  repeat
    mp_div_w(t, x, @t, rd);
    if (sign<0) and (rd<>0) then mp_inc(t);

    mp_div_w(t, q, @v, rd);
    if (sign<0) and (rd<>0) then mp_inc(v);

    if sign>0 then mp_add(sum, v, sum) else mp_sub(sum, v, sum);

    sign := -sign;
    q := q+2;
    if q>$FFF0 then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_arctanw: prec too large');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        break;
      {$endif}
    end;
  until mp_iszero(t) or (mp_error<>MP_OKAY);
  {Chop off extra precision}
  {sum := trunc(sum/Radix^EXTRA + 0.5}
  mp_set(v,Radix);
  mp_expt_w(v,EXTRA,v);
  mp_div_2(v,t);
  mp_add(sum,t,sum);
  mp_div(sum, v, sum);

  mp_clear2(v,t);
end;


{---------------------------------------------------------------------------}
procedure mp_atanhw(mul, x, Radix, prec: word; var sum: mp_int);
  {-Compute sum := mul*atanh(1/x), to prec Radix digits}
var
  t,v:  mp_int;
  q,rd: word;
const
  EXTRA = 10;
begin

  if mp_error<>MP_OKAY then exit;
  {$ifdef MPC_ArgCheck}
    if mp_not_init(sum) then begin
      {$ifdef MPC_HaltOnArgCheck}
        {$ifdef MPC_UseExceptions}
          raise MPXNotInit.Create('mp_atanh');
        {$else}
          RunError(MP_RTE_NOTINIT);
        {$endif}
      {$else}
        set_mp_error(MP_NOTINIT);
        exit;
      {$endif}
    end;
  {$endif}
  {Algorithm from MPI/Utils pi.c}
  q := 1;
  {Check x^2 <= $FFFF, Radix <= MaxRadix}
  if (x > $FF) or (Radix>MAXRadix) then begin
    {$ifdef MPC_HaltOnError}
     {$ifdef MPC_UseExceptions}
       raise MPXRange.Create('mp_atanhw: x or Radix too large');
     {$else}
       RunError(MP_RTE_RANGE);
     {$endif}
    {$else}
      set_mp_error(MP_RANGE);
      exit;
    {$endif}
  end;

  mp_init2(t,v);
  if mp_error<>MP_OKAY then exit;

  {Increase precision to guard against rounding errors}
  mp_set_w(t, Radix);
  mp_expt_w(t, prec + EXTRA, t);
  mp_mul_w(t,mul,t);
  mp_mul_w(t,x,t);

  {                1     1       1       1         }
  { atanh(1/x) =   - + ----- + ----- + ----- + ... }
  {                x   3 x^3   5 x^5   7 x^7       }

  x := x*x;
  mp_zero(sum);

  {prec must be <= trunc(2^16*ln(x)/ln(Radix))-EXTRA}
  {otherwise the q will become greater than 2^16    }

  repeat
    mp_div_w(t, x, @t, rd);
    mp_div_w(t, q, @v, rd);
    mp_add(sum, v, sum);
    q := q+2;
    if q>$FFF0 then begin
      {$ifdef MPC_HaltOnError}
        {$ifdef MPC_UseExceptions}
          raise MPXBadArg.Create('mp_atanhw: prec too large');
        {$else}
          RunError(MP_RTE_BADARG);
        {$endif}
      {$else}
        set_mp_error(MP_BADARG);
        break;
      {$endif}
    end;
  until mp_iszero(t) or (mp_error<>MP_OKAY);
  {Chop off extra precision}
  {sum := trunc(sum/Radix^EXTRA + 0.5}
  mp_set(v,Radix);
  mp_expt_w(v,EXTRA,v);
  mp_div_2(v,t);
  mp_add(sum,t,sum);
  mp_div(sum, v, sum);

  mp_clear2(v,t);
end;


{---------------------------------------------------------------------------}
procedure mp_dump_diagctr;
  {-Dump diagnostic counters}
{$ifdef MPC_Diagnostic}
var
  i: integer;
{$endif}
begin
  {$ifdef MPC_Diagnostic}
    write('mp_diagctr ');
    for i:=0 to MAX_DiagCTR do write(i:1,':',mp_diagctr[i],'  ');
    writeln;
  {$endif}
end;


{---------------------------------------------------------------------------}
procedure mp_reset_diagctr;
  {-Reset diagnostic counters to 0}
begin
  {$ifdef MPC_Diagnostic}
    fillchar(mp_diagctr, sizeof(mp_diagctr),0);
  {$endif}
end;


{---------------------------------------------------------------------------}
procedure mp_dump_meminfo;
  {-Write mp_memstat to output (if MPC_Diagnostic is defined)}
begin
  {$ifdef MPC_Diagnostic}
    writeln('Memory status:');
    with mp_memstat do begin
      writeln(' MemDiff : ', MemDiff);
      writeln(' InitDiff: ', InitDiff);
      writeln(' ACntHigh: ', ACntHigh);
      writeln(' ACntLow : ', ACntLow);
    end;
  {$endif}
end;


{---------------------------------------------------------------------------}
procedure mp_dump_memused;
  {-Write mp_memused to output (if MPC_Diagnostic is defined)}
begin
  {$ifdef MPC_Diagnostic}
    writeln(' MemDiff : ', mp_memstat.MemDiff);
  {$endif}
end;


{---------------------------------------------------------------------------}
function mp_memused: longint;
  {-Return total allocated memory, MaxLongint if MPC_Diagnostic not defined}
begin
  {$ifdef MPC_Diagnostic}
    mp_memused := mp_memstat.MemDiff;
  {$else}
    mp_memused := MaxLongint;
  {$endif}
end;


{---------------------------------------------------------------------------}
procedure mp_dumpu(const a: mp_int; const hdr: mp_string);
  {-Dump used fields of a}
var
  i: TNInt;
const
  width = 1 + (DIGIT_BIT*146+484) div 485;   {146/485 appr. ln(2)/ln(10)}
  num   = 70 div width;
begin
  writeln('*** Dump of ', hdr);
  if mp_not_init(a) then begin
    writeln(' *(not initialized)*');
    exit;
  end;
  with a do begin
    writeln('alloc: ', alloc, '  used: ', used, '  sign: ', sign);
    write  ('pdigits: ');
    if used=0 then writeln('---')
    else begin
      for i:=0 to pred(used) do begin
        write(pdigits^[i]:width);
        if i=pred(used) then writeln
        else if i mod num = pred(num) then begin
          writeln;
          write('         ');
        end;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mp_dumpa(const a: mp_int; const hdr: mp_string);
  {-Dump all fields of a}
var
  i: TNInt;
const
  width = 1 + (DIGIT_BIT*146+484) div 485;   {146/485 appr. ln(2)/ln(10)}
  num   = 70 div width;
begin
  writeln('*** Dump of ', hdr);
  if mp_not_init(a) then begin
    writeln(' *(not initialized)*');
    exit;
  end;
  with a do begin
    writeln('alloc: ', alloc, '  used: ', used, '  sign: ', sign);
    write('pdigits: ');
    if alloc=0 then writeln('---')
    else begin
      for i:=0 to pred(alloc) do begin
        write(pdigits^[i]:width);
        if i=pred(alloc) then writeln
        else if i mod num = pred(num) then begin
          writeln;
          write('         ');
        end;
      end;
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mpf_dump_me(const a: mp_float);
  {-Dump mantissa and exponent of a}
begin
  writeln('<M:',mp_decimal(a.mantissa), ', X:', a.exponent,'>');
end;



end.

