unit well1024;

(*************************************************************************

 DESCRIPTION   :  well1024a based pseudo random number generator
                  Period about 2**1024

 REQUIREMENTS  :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  [1] http://www.iro.umontreal.ca/~panneton/well/WELL1024a.c
                  [2] F. Panneton, P. L'Ecuyer and M. Matsumoto, "Improved
                      Long-Period Generators Based on Linear Recurrences Modulo 2",
                      http://www.iro.umontreal.ca/~lecuyer/myftp/papers/lfsr04.pdf

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     10.02.14  W.Ehrhardt  Initial BP7 version, prng_next
 0.11     11.02.14  we          Remaining functions
 0.12     12.02.14  we          Selftest
 0.13     04.12.17  we          BIT16 inline code
 0.14     04.12.17  we          BASM16 code
 0.15     04.12.17  we          SHR8
 0.16     22.12.17  we          well1024a_rangel/w
 0.17     25.12.17  we          rangel for old non-basm compilers
**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2014-2017 Wolfgang Ehrhardt

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
  TW1024a_State = array[0..31] of longint;

type
  well1024a_ctx = record
                    nr     : longint;       {next random }
                    state  : TW1024a_State; {state vector}
                    state_i: longint;       {state index }
                  end;


procedure well1024a_init(var ctx: well1024a_ctx; seed: longint);
  {-Init context from seed}

procedure well1024a_init0(var ctx: well1024a_ctx);
  {-Init context from randseed}

procedure well1024a_next(var ctx: well1024a_ctx);
  {-Next step of PRNG}

procedure well1024a_read(var ctx: well1024a_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}

function  well1024a_long(var ctx: well1024a_ctx): longint;
  {-Next random positive longint}

function  well1024a_dword(var ctx: well1024a_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}

function  well1024a_word(var ctx: well1024a_ctx): word;
  {-Next random word}

function  well1024a_double(var ctx: well1024a_ctx): double;
  {-Next random double [0..1) with 32 bit precision}

function  well1024a_double53(var ctx: well1024a_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}

function  well1024a_rangew(var ctx: well1024a_ctx; range: word): word;
  {-Next random word in range 0..range-1}

function  well1024a_rangel(var ctx: well1024a_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}

function  well1024a_selftest: boolean;
  {-Simple self-test of PRNG}


implementation


uses
  BTypes, rndrange;


{$ifndef BIT16}

{---------------------------------------------------------------------------}
procedure well1024a_next(var ctx: well1024a_ctx);
  {-Next step of PRNG}
var
  z0,z1,z2,i: longint;
begin
  with ctx do begin
    z1 := state[(state_i + 3)  and 31];
    z1 := state[state_i] xor z1 xor (z1 shr 8);
    z2 := state[(state_i + 24)  and 31];
    z0 := state[(state_i + 10)  and 31];
    z2 := z2 xor (z2 shl 19) xor z0 xor (z0 shl 14);
    i  := (state_i + 31) and 31;
    z0 := state[i];
    state[state_i] := z1 xor z2;
    state[i] := z0 xor (z0 shl 11) xor z1 xor (z1 shl 7) xor z2 xor (z2 shl 13);
    state_i  := i;
    nr := state[i];
  end;
end;

{$else}

{$ifdef BASM16}

{---------------------------------------------------------------------------}
procedure well1024a_next(var ctx: well1024a_ctx);
  {-Next step of PRNG}
var
  z0,z1,z2,z3: longint;
  i: integer;
begin
  with ctx do begin
    z1 := state[state_i];
    z0 := state[(state_i + 10)  and 31];
    z2 := state[(state_i + 24)  and 31];
    z3 := state[(state_i +  3)  and 31];
    asm
      {z1 := z1 xor z3 xor (z3 shr 8);}
      db $66;   mov ax,word[z3]
      db $66;   mov dx,ax
      db $66;   shr ax,8
      db $66;   xor ax,dx
      db $66;   xor word[z1],ax
      {z2 := z2 xor LShift(z2,19) xor z0 xor LShift(z0,14);}
      db $66;   mov ax,word[z2]
      db $66;   shl ax,19
      db $66;   mov dx,word[z0]
      db $66;   xor ax,dx
      db $66;   shl dx,14
      db $66;   xor ax,dx
      db $66;   xor word[z2],ax
    end;
    i  := (state_i + 31) and 31;
    z0 := state[i];
    state[state_i] := z1 xor z2;
    asm
      {z3 := z0 xor LShift(z0,11) xor z1 xor LShift(z1,7) xor z2 xor LShift(z2,13);}
      db $66;   mov ax,word[z0]
      db $66;   mov dx,ax
      db $66;   shl ax,11
      db $66;   xor dx,ax      {z0 xor LShift(z0,11)}

      db $66;   mov ax,word[z1]
      db $66;   mov cx,ax
      db $66;   shl ax,7
      db $66;   xor ax,cx      {z1 xor LShift(z1,7)}
      db $66;   xor dx,ax

      db $66;   mov ax,word[z2]
      db $66;   mov cx,ax
      db $66;   shl cx,13
      db $66;   xor ax,cx      {z2 xor LShift(z2,13)}

      db $66;   xor ax,dx
      db $66;   mov word[z3],ax
    end;
    state[i] := z3;
    state_i  := i;
    nr := state[i];
  end;
end;

{$else}

{** TP 5X}
{---------------------------------------------------------------------------}
function LShift(X: longint; c: word): longint;
  {-Shift X left by c bits}
inline(
  $59/           {   pop    cx    }
  $58/           {   pop    ax    }
  $5A/           {   pop    dx    }
  $83/$E1/$1F/   {   and    cx,31 }
  $74/$06/       {   je     0792  }
  $D1/$E0/       {L: shl    ax,1  }
  $D1/$D2/       {   rcl    dx,1  }
  $E2/$FA);      {   loop   L     }
                 {X:              }

{---------------------------------------------------------------------------}
function SHR8(X: longint): longint;
  {-Shift X left by 8 bits}
inline(
  $58/           {   pop    ax    }
  $5A/           {   pop    dx    }
  $8A/$C4/       {   mov    al,ah }
  $8A/$E2/       {   mov    ah,dl }
  $8A/$D6/       {   mov    dl,dh }
  $2A/$F6);       {   sub    dh,dh }


{---------------------------------------------------------------------------}
procedure well1024a_next(var ctx: well1024a_ctx);
  {-Next step of PRNG}
var
  z0,z1,z2: longint;
  i: integer;
begin
  with ctx do begin
    z1 := state[(state_i + 3)  and 31];
    z1 := state[state_i] xor z1 xor SHR8(z1);
    z2 := state[(state_i + 24)  and 31];
    z0 := state[(state_i + 10)  and 31];
    z2 := z2 xor LShift(z2,19) xor z0 xor LShift(z0,14);
    i  := (state_i + 31) and 31;
    z0 := state[i];
    state[state_i] := z1 xor z2;
    state[i] := z0 xor LShift(z0,11) xor z1 xor LShift(z1,7) xor z2 xor LShift(z2,13);
    state_i  := i;
    nr := state[i];
  end;
end;
{$endif}

{$endif}


{---------------------------------------------------------------------------}
procedure well1024a_init(var ctx: well1024a_ctx; seed: longint);
  {-Init context from seed}
const
  M=69069;
  A=1;
var
  i: integer;
begin
  for i:=0 to 31 do begin
    ctx.state[i] := seed;
    seed := M*seed+A;
  end;
  ctx.state_i := 0;
end;


{---------------------------------------------------------------------------}
procedure well1024a_init0(var ctx: well1024a_ctx);
  {-Init context from randseed}
begin
  well1024a_init(ctx, randseed);
end;


{---------------------------------------------------------------------------}
function well1024a_long(var ctx: well1024a_ctx): longint;
  {-Next random positive longint}
begin
  well1024a_next(ctx);
  well1024a_long := ctx.nr shr 1;
end;


{---------------------------------------------------------------------------}
function well1024a_dword(var ctx: well1024a_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}
begin
  well1024a_next(ctx);
  {$ifdef HAS_CARD32}
    well1024a_dword := cardinal(ctx.nr);
  {$else}
    well1024a_dword := ctx.nr;
  {$endif}
end;


{---------------------------------------------------------------------------}
function well1024a_word(var ctx: well1024a_ctx): word;
  {-Next random word}
type
  TwoWords = packed record
               L,H: word
             end;
begin
  well1024a_next(ctx);
  well1024a_word := TwoWords(ctx.nr).H;
end;


{---------------------------------------------------------------------------}
function well1024a_double(var ctx: well1024a_ctx): double;
  {-Next random double [0..1) with 32 bit precision}
begin
  well1024a_next(ctx);
  well1024a_double := (ctx.nr + 2147483648.0) / 4294967296.0;
end;


{---------------------------------------------------------------------------}
function well1024a_double53(var ctx: well1024a_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}
var
  hb,lb: longint;
begin
  well1024a_next(ctx);
  hb := ctx.nr shr 5;
  well1024a_next(ctx);
  lb := ctx.nr shr 6;
  well1024a_double53 := (hb*67108864.0+lb)/9007199254740992.0;
end;


{---------------------------------------------------------------------------}
procedure well1024a_read(var ctx: well1024a_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}
type
  plong = ^longint;
begin
  while len>3 do begin
    well1024a_next(ctx);
    plong(dest)^ := ctx.nr;
    inc(Ptr2Inc(dest),4);
    dec(len, 4);
  end;
  if len>0 then begin
    well1024a_next(ctx);
    move(ctx.nr, dest^, len and 3);
  end;
end;


{---------------------------------------------------------------------------}
function well1024a_rangew(var ctx: well1024a_ctx; range: word): word;
  {-Next random word in range 0..range-1}
begin
  well1024a_next(ctx);
  well1024a_rangew := rangew(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function well1024a_rangel(var ctx: well1024a_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}
begin
  well1024a_next(ctx);
  well1024a_rangel := rangel(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function well1024a_selftest: boolean;
  {-Simple self-test of PRNG}
var
  ctx: well1024a_ctx;
  i: integer;
begin
  well1024a_init(ctx,123456789);
  for i:=1 to 256 do well1024a_next(ctx);
  well1024a_selftest := ctx.nr = $23A29784;
end;

end.
