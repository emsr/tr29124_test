unit taus113;

(*************************************************************************

 DESCRIPTION   :  Taus113 based pseudo random number generator
                  Period about 2**113

 REQUIREMENTS  :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  [1] P. L'Ecuyer, "Tables of Maximally-Equidistributed
                      Combined LFSR Generators", Mathematics of Computation,
                      Vol. 68, 225 (1999), 261-269.
                  [2] http://www.iro.umontreal.ca/~lecuyer/myftp/papers/tausme2.ps
                      Online version of [1]


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     02.08.05  W.Ehrhardt  Initial BP7 version from taus88
 0.11     02.08.05  we          References, bugfix
 0.12     05.11.08  we          taus113_dword function
 0.13     04.12.08  we          BTypes/Ptr2Inc
 0.14     14.06.12  we          Fix bug in _read for trailing max 3 bytes
 0.15     22.12.17  we          taus113_rangel/w
 0.16     25.12.17  we          rangel for old non-basm compilers
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
  taus113_ctx = record
               nr         : longint;  {next rand = s1 xor s2 xor s3 xor s4}
               s1,s2,s3,s4: longint;  {state variables }
             end;

procedure taus113_init(var ctx: taus113_ctx; seed: longint);
  {-Init context from seed}

procedure taus113_init0(var ctx: taus113_ctx);
  {-Init context from randseed}

procedure taus113_init4(var ctx: taus113_ctx; seed1,seed2,seed3,seed4: longint);
  {-Init all context with separate seeds, the initial seeds}
  { seed1, seed2, seed3, seed4 should be >= 2, 8, 16, 128 resp. or < 0  }

procedure taus113_next(var ctx: taus113_ctx);
  {-Next step of PRNG}

procedure taus113_read(var ctx: taus113_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}

function  taus113_long(var ctx: taus113_ctx): longint;
  {-Next random positive longint}

function  taus113_dword(var ctx: taus113_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}

function  taus113_word(var ctx: taus113_ctx): word;
  {-Next random word}

function  taus113_double(var ctx: taus113_ctx): double;
  {-Next random double [0..1) with 32 bit precision}

function  taus113_double53(var ctx: taus113_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}

function  taus113_rangew(var ctx: taus113_ctx; range: word): word;
  {-Next random word in range 0..range-1}

function  taus113_rangel(var ctx: taus113_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}

function  taus113_selftest: boolean;
  {-Simple self-test of PRNG}


implementation


uses
  BTypes, rndrange;

const
   TC1 = longint($FFFFFFFE);  {4294967294}
   TC2 = longint($FFFFFFF8);  {4294967288}
   TC3 = longint($FFFFFFF0);  {4294967280}
   TC4 = longint($FFFFFF80);  {4294967168}
  ICNT = 10;                  {Init count to reach valid state}

{$ifdef BASM16}

{---------------------------------------------------------------------------}
procedure taus113_next(var ctx: taus113_ctx);
  {-Next step of PRNG}
begin
  asm
            push ds
            lds  si, [ctx]

    {s1 := ((s1 and TC1) shl 18) xor (((s1 shl  6) xor s1) shr 13);}
    db $66; mov  ax, word ptr taus113_ctx[si].s1
    db $66; mov  dx, ax
    db $66; shl  dx, 6
    db $66; xor  dx, ax
    db $66; shr  dx, 13
            and  ax, TC1  {use only low word, high word = FFFF}
    db $66; shl  ax, 18
    db $66; xor  ax, dx
    db $66; mov  bx, ax
    db $66; mov  word ptr taus113_ctx[si].s1,ax

    {s2 := ((s2 and TC2) shl  2) xor (((s2 shl  2) xor s2) shr 27);}
    db $66; mov  ax, word ptr taus113_ctx[si].s2
    db $66; mov  dx, ax
    db $66; shl  dx, 2
    db $66; xor  dx, ax
    db $66; shr  dx, 27
            and  ax, TC2  {use only low word, high word = FFFF}
    db $66; shl  ax, 2
    db $66; xor  ax, dx
    db $66; xor  bx, ax
    db $66; mov  word ptr taus113_ctx[si].s2,ax

    {s3 := ((s3 and TC3) shl  7) xor (((s3 shl 13) xor s3) shr 21);}
    db $66; mov  ax, word ptr taus113_ctx[si].s3
    db $66; mov  dx, ax
    db $66; shl  dx, 13
    db $66; xor  dx, ax
    db $66; shr  dx, 21
            and  ax, TC3  {use only low word, high word = FFFF}
    db $66; shl  ax, 7
    db $66; xor  ax, dx
    db $66; xor  bx, ax
    db $66; mov  word ptr taus113_ctx[si].s3,ax

    {s4 := ((s4 and TC4) shl 13) xor (((s4 shl  3) xor s4) shr 12);}
    db $66; mov  ax, word ptr taus113_ctx[si].s4
    db $66; mov  dx, ax
    db $66; shl  dx, 3
    db $66; xor  dx, ax
    db $66; shr  dx, 12
            and  ax, TC4  {use only low word, high word = FFFF}
    db $66; shl  ax, 13
    db $66; xor  ax, dx
    db $66; xor  bx, ax
    db $66; mov  word ptr taus113_ctx[si].s4,ax

    {nr := s1 xor s2 xor s3 xor s4;}
    db $66; mov  word ptr taus113_ctx[si].nr,bx
            pop  ds
  end;
end;

{$else}

{---------------------------------------------------------------------------}
procedure taus113_next(var ctx: taus113_ctx);
  {-Next step of PRNG}
begin
  with ctx do begin
    s1 := ((s1 and TC1) shl 18) xor (((s1 shl  6) xor s1) shr 13);
    s2 := ((s2 and TC2) shl  2) xor (((s2 shl  2) xor s2) shr 27);
    s3 := ((s3 and TC3) shl  7) xor (((s3 shl 13) xor s3) shr 21);
    s4 := ((s4 and TC4) shl 13) xor (((s4 shl  3) xor s4) shr 12);
    nr := s1 xor s2 xor s3 xor s4;
  end;
end;

{$endif}



{---------------------------------------------------------------------------}
procedure taus113_init4(var ctx: taus113_ctx; seed1,seed2,seed3,seed4: longint);
  {-Init all context with separate seeds, the initial seeds}
  { seed1, seed2, seed3, seed4 should be >= 2, 8, 16, 128 resp. or < 0  }
var
  i: integer;
begin
  with ctx do begin
    s1 := seed1;
    s2 := seed2;
    s3 := seed3;
    s4 := seed4;
    {Make sure most significant bits are set}
    if s1 and TC1=0 then inc(s1,2);
    if s2 and TC2=0 then inc(s2,8);
    if s3 and TC3=0 then inc(s3,16);
    if s4 and TC4=0 then inc(s4,128);
    {Make sure ctx has valid state}
    for i:=1 to ICNT do taus113_next(ctx);
  end;
end;


{---------------------------------------------------------------------------}
procedure taus113_init(var ctx: taus113_ctx; seed: longint);
  {-Init context from seed}
const
  M=69069;
  A=1;
var
  seed2,seed3,seed4: longint;
begin
  seed2 := M*seed+A;
  seed3 := M*seed2+A;
  seed4 := M*seed3+A;
  taus113_init4(ctx,seed,seed2,seed3,seed4);
end;


{---------------------------------------------------------------------------}
procedure taus113_init0(var ctx: taus113_ctx);
  {-Init context from randseed}
begin
  taus113_init(ctx, randseed);
end;


{---------------------------------------------------------------------------}
function taus113_long(var ctx: taus113_ctx): longint;
  {-Next random positive longint}
begin
  taus113_next(ctx);
  taus113_long := ctx.nr shr 1;
end;


{---------------------------------------------------------------------------}
function taus113_dword(var ctx: taus113_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}
begin
  taus113_next(ctx);
  {$ifdef HAS_CARD32}
    taus113_dword := cardinal(ctx.nr);
  {$else}
    taus113_dword := ctx.nr;
  {$endif}
end;


{---------------------------------------------------------------------------}
function taus113_word(var ctx: taus113_ctx): word;
  {-Next random word}
type
  TwoWords = packed record
               L,H: word
             end;
begin
  taus113_next(ctx);
  taus113_word := TwoWords(ctx.nr).H;
end;


{---------------------------------------------------------------------------}
function taus113_double(var ctx: taus113_ctx): double;
  {-Next random double [0..1) with 32 bit precision}
begin
  taus113_next(ctx);
  taus113_double := (ctx.nr + 2147483648.0) / 4294967296.0;
end;


{---------------------------------------------------------------------------}
function taus113_double53(var ctx: taus113_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}
var
  hb,lb: longint;
begin
  taus113_next(ctx);
  hb := ctx.nr shr 5;
  taus113_next(ctx);
  lb := ctx.nr shr 6;
  taus113_double53 := (hb*67108864.0+lb)/9007199254740992.0;
end;


{---------------------------------------------------------------------------}
procedure taus113_read(var ctx: taus113_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}
type
  plong = ^longint;
begin
  while len>3 do begin
    taus113_next(ctx);
    plong(dest)^ := ctx.nr;
    inc(Ptr2Inc(dest),4);
    dec(len, 4);
  end;
  if len>0 then begin
    taus113_next(ctx);
    move(ctx.nr, dest^, len and 3);
  end;
end;


{---------------------------------------------------------------------------}
function taus113_rangew(var ctx: taus113_ctx; range: word): word;
  {-Next random word in range 0..range-1}
begin
  taus113_next(ctx);
  taus113_rangew := rangew(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function taus113_rangel(var ctx: taus113_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}
begin
  taus113_next(ctx);
  taus113_rangel := rangel(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function taus113_selftest: boolean;
  {-Simple self-test of PRNG}
var
  ctx: taus113_ctx;
  i: integer;
begin
  {Test values from t_taus113.c based on [2]}
  taus113_init4(ctx,2,8,16,128);
  { First: 00180820
      500: b2209ae9
     1000: 37d83989
     2500: 3bfcadbd
    10000: 26620794
    20000: 7aaa2b26
    30000: a98802f6}
  for i:=1 to 500-ICNT do taus113_next(ctx);
  taus113_selftest := ctx.nr=longint($b2209ae9);
end;

end.
