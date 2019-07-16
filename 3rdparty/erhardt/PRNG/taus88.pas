unit Taus88;

(*************************************************************************

 DESCRIPTION   :  Taus88 based pseudo random number generator
                  Period about 2**88

 REQUIREMENTS  :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  [1] P. L'Ecuyer, "Maximally Equidistributed Combined Tausworthe
                      Generators", Mathematics of Computation 65, 213 (1996), 203-213
                  [2] http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme.ps
                      Online version of [1] with corrections


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     23.04.05  W.Ehrhardt  Initial BP7 version, prng_next, self-test
 0.11     24.04.05  we          Remaining functions
 0.12     24.04.05  we          {$N+} only for BIT16, long/word/double
 0.13     24.04.05  we          BASM16
 0.14     24.04.05  we          prng_double: first add then divide
 0.15     11.05.05  we          Bugfix and consts M,A
 0.16     29.05.05  we          renamed unit to PRNG, .sn to .nr
 0.17     30.05.05  we          new selftest values, constant ICNT
 0.18     02.08.05  we          Bugfix: inc s3 in prng_init3
 0.19     02.08.05  we          Changed prng_ to taus88_
 0.20     05.11.08  we          taus88_dword function
 0.21     02.12.08  we          BTypes/Ptr2Inc
 0.22     07.01.09  we          Uses BTypes moved to implementation
 0.23     14.06.12  we          Fix bug in _read for trailing max 3 bytes
 0.24     21.12.17  we          taus88_rangel/w
 0.25     25.12.17  we          rangel for old non-basm compilers
**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2005-2017 Wolfgang Ehrhardt

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

interface

{$i std.inc}

{$ifdef BIT16}
  {$N+}
{$endif}

type
  taus88_ctx = record
                 s1,s2,s3: longint;  {state variables }
                 nr      : longint;  {next rand = s1 xor s2 xor s3}
               end;

procedure taus88_init(var ctx: taus88_ctx; seed: longint);
  {-Init context from seed}

procedure taus88_init0(var ctx: taus88_ctx);
  {-Init context from randseed}

procedure taus88_init3(var ctx: taus88_ctx; seed1,seed2,seed3: longint);
  {-Init all context with separate seeds, the initial seeds}
  { seed1, seed2, seed3 should be >= 2, 8, 16 resp. or < 0  }

procedure taus88_next(var ctx: taus88_ctx);
  {-Next step of PRNG}

procedure taus88_read(var ctx: taus88_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}

function  taus88_long(var ctx: taus88_ctx): longint;
  {-Next random positive longint}

function  taus88_dword(var ctx: taus88_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}

function  taus88_word(var ctx: taus88_ctx): word;
  {-Next random word}

function  taus88_double(var ctx: taus88_ctx): double;
  {-Next random double [0..1) with 32 bit precision}

function  taus88_double53(var ctx: taus88_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}

function  taus88_rangew(var ctx: taus88_ctx; range: word): word;
  {-Next random word in range 0..range-1}

function  taus88_rangel(var ctx: taus88_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}

function  taus88_selftest: boolean;
  {-Simple self-test of PRNG}


implementation

uses
  BTypes, rndrange;

const
   TC1 = longint($FFFFFFFE);  {4294967294}
   TC2 = longint($FFFFFFF8);  {4294967288}
   TC3 = longint($FFFFFFF0);  {4294967280}

  ICNT = 6;                   {Init count to reach valid state}

{$ifdef BASM16}

{---------------------------------------------------------------------------}
procedure taus88_next(var ctx: taus88_ctx);
  {-Next step of PRNG}
begin
  asm
            push ds
            lds  si, [ctx]
    {s1 := ((s1 and TC1) shl 12) xor (((s1 shl 13) xor s1) shr 19);}
    db $66; mov  ax, word ptr taus88_ctx[si].s1
    db $66; mov  dx, ax
    db $66; shl  dx, 13
    db $66; xor  dx, ax
    db $66; shr  dx, 19
            and  ax, TC1  {use only low word, high word = FFFF}
    db $66; shl  ax, 12
    db $66; xor  ax, dx
    db $66; mov  bx, ax
    db $66; mov  word ptr taus88_ctx[si].s1,ax

    {s2 := ((s2 and TC2) shl  4) xor (((s2 shl  2) xor s2) shr 25);}
    db $66; mov  ax, word ptr taus88_ctx[si].s2
    db $66; mov  dx, ax
    db $66; shl  dx, 2
    db $66; xor  dx, ax
    db $66; shr  dx, 25
            and  ax, TC2  {use only low word, high word = FFFF}
    db $66; shl  ax, 4
    db $66; xor  ax, dx
    db $66; xor  bx, ax
    db $66; mov  word ptr taus88_ctx[si].s2,ax

    {s3 := ((s3 and TC3) shl 17) xor (((s3 shl  3) xor s3) shr 11);}
    db $66; mov  ax, word ptr taus88_ctx[si].s3
    db $66; mov  dx, ax
    db $66; shl  dx, 3
    db $66; xor  dx, ax
    db $66; shr  dx, 11
            and  ax, TC3  {use only low word, high word = FFFF}
    db $66; shl  ax, 17
    db $66; xor  ax, dx
    db $66; xor  bx, ax
    db $66; mov  word ptr taus88_ctx[si].s3,ax

    {nr := s1 xor s2 xor s3;}
    db $66; mov  word ptr taus88_ctx[si].nr,bx
            pop  ds
  end;
end;

{$else}

{---------------------------------------------------------------------------}
procedure taus88_next(var ctx: taus88_ctx);
  {-Next step of PRNG}
begin
  with ctx do begin
    s1 := ((s1 and TC1) shl 12) xor (((s1 shl 13) xor s1) shr 19);
    s2 := ((s2 and TC2) shl  4) xor (((s2 shl  2) xor s2) shr 25);
    s3 := ((s3 and TC3) shl 17) xor (((s3 shl  3) xor s3) shr 11);
    nr := s1 xor s2 xor s3;
  end;
end;

{$endif}



{---------------------------------------------------------------------------}
procedure taus88_init3(var ctx: taus88_ctx; seed1,seed2,seed3: longint);
  {-Init all context with separate seeds, the initial seeds}
  { seed1, seed2, seed3 should be >= 2, 8, 16 resp. or < 0  }
var
  i: integer;
begin
  with ctx do begin
    s1 := seed1;
    s2 := seed2;
    s3 := seed3;
    {Make sure most significant bits are set, c.f. [2], correction on p.10}
    if s1 and TC1=0 then inc(s1,2);
    if s2 and TC2=0 then inc(s2,8);
    if s3 and TC3=0 then inc(s3,16);   {Bugfix: inc s3!}
    {Make sure ctx has valid state, c.f. [2], p.3}
    for i:=1 to ICNT do taus88_next(ctx);
  end;
end;


{---------------------------------------------------------------------------}
procedure taus88_init(var ctx: taus88_ctx; seed: longint);
  {-Init context from seed}
const
  M=69069;
  A=1;
var
  seed2,seed3: longint;
begin
  seed2 := M*seed+A;
  seed3 := M*seed2+A;
  taus88_init3(ctx,seed,seed2,seed3);
end;


{---------------------------------------------------------------------------}
procedure taus88_init0(var ctx: taus88_ctx);
  {-Init context from randseed}
begin
  taus88_init(ctx, randseed);
end;


{---------------------------------------------------------------------------}
function taus88_long(var ctx: taus88_ctx): longint;
  {-Next random positive longint}
begin
  taus88_next(ctx);
  taus88_long := ctx.nr shr 1;
end;


{---------------------------------------------------------------------------}
function taus88_dword(var ctx: taus88_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}
begin
  taus88_next(ctx);
  {$ifdef HAS_CARD32}
    taus88_dword := cardinal(ctx.nr);
  {$else}
    taus88_dword := ctx.nr;
  {$endif}
end;


{---------------------------------------------------------------------------}
function taus88_word(var ctx: taus88_ctx): word;
  {-Next random word}
type
  TwoWords = packed record
               L,H: word
             end;
begin
  taus88_next(ctx);
  taus88_word := TwoWords(ctx.nr).H;
end;


{---------------------------------------------------------------------------}
function taus88_double(var ctx: taus88_ctx): double;
  {-Next random double [0..1) with 32 bit precision}
begin
  taus88_next(ctx);
  taus88_double := (ctx.nr + 2147483648.0) / 4294967296.0;
end;


{---------------------------------------------------------------------------}
function taus88_double53(var ctx: taus88_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}
var
  hb,lb: longint;
begin
  taus88_next(ctx);
  hb := ctx.nr shr 5;
  taus88_next(ctx);
  lb := ctx.nr shr 6;
  taus88_double53 := (hb*67108864.0+lb)/9007199254740992.0;
end;


{---------------------------------------------------------------------------}
procedure taus88_read(var ctx: taus88_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}
type
  plong = ^longint;
begin
  while len>3 do begin
    taus88_next(ctx);
    plong(dest)^ := ctx.nr;
    inc(Ptr2Inc(dest),4);
    dec(len, 4);
  end;
  if len>0 then begin
    taus88_next(ctx);
    move(ctx.nr, dest^, len and 3);
  end;
end;


{---------------------------------------------------------------------------}
function taus88_rangew(var ctx: taus88_ctx; range: word): word;
  {-Next random word in range 0..range-1}
begin
  taus88_next(ctx);
  taus88_rangew := rangew(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function taus88_rangel(var ctx: taus88_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}
begin
  taus88_next(ctx);
  taus88_rangel := rangel(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function taus88_selftest: boolean;
  {-Simple self-test of PRNG}
var
  ctx: taus88_ctx;
  i: integer;
begin
  {Test values from t_taus88.c based on [2]}
  taus88_init3(ctx,2,8,16);
  {  500: $4f87c056
    1000: $fed1530b
    5000: $91d993a0
   10000: $438f5e83
   20000: $f4645720
   30000: $d5ccc345 }
  for i:=1 to 500-ICNT do taus88_next(ctx);
  taus88_selftest := ctx.nr=longint($4f87c056);
end;

end.
