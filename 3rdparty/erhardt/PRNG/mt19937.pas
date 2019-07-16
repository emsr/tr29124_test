unit MT19937;

(*************************************************************************

 DESCRIPTION   :  MT19937 based pseudo random number generator
                    Period about 2**19937

 REQUIREMENTS  :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    : [1] Matsumoto, M. and Nishimura, T. "Mersenne Twister:
                     A 623-Dimensionally Equidistributed Uniform Pseudo-Random
                     Number Generator", ACM Transactions on Modeling and Computer
                     Simulation, Vol. 8, (1998) 3-30.

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     18.07.05  W.Ehrhardt  Initial BP7 port of mt19937ar.c
 0.11     18.07.05  we          TP5-7
 0.12     24.07.05  we          Tempering with BASM16
 0.13     05.11.08  we          mt19937_dword function
 0.14     04.12.08  we          BTypes/Ptr2Inc
 0.15     14.06.12  we          Fix bug in _read for trailing max 3 bytes
 0.16     22.12.17  we          mt19937_rangel/w
 0.17     25.12.17  we          rangel for old non-basm compilers
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


(*---------------------------------------------------------------------------
 This Pascal code is based on mt19937ar.c, a C-program  with initialization
 improved 2002/1/26.  Coded by Takuji Nishimura and Makoto Matsumoto.

 Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
 All rights reserved.

 Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions
 are met:

   1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

   3. The names of its contributors may not be used to endorse or promote
      products derived from this software without specific prior written
      permission.

 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
---------------------------------------------------------------------------*)

interface

{$i std.inc}

{$ifdef BIT16}
  {$N+}
{$endif}

uses
  BTypes;

const
  N19937 = 624;
  M19937 = 397;

type

  mt19937_state = array[0..N19937-1] of longint;

  mt19937_ctx   = record
                    mt : mt19937_state; {state vector}
                    mti: integer;       {state index }
                    nr : longint;       {next random }
                  end;



procedure mt19937_init (var ctx: mt19937_ctx; seed: longint);
  {-Init context from seed}

procedure mt19937_init0(var ctx: mt19937_ctx);
  {-Init context from randseed}

{$ifdef CONST}
procedure mt19937_inita(var ctx: mt19937_ctx; const key: array of longint; klen: integer);
  {-Init all context variables with separate seeds, klen: number of seeds}
{$else}
procedure mt19937_inita(var ctx: mt19937_ctx; var KArr; klen: integer);
  {-Init all context variables with separate seeds, klen: number of seeds}
{$endif}

procedure mt19937_next(var ctx: mt19937_ctx);
  {-Next step of PRNG}

procedure mt19937_read(var ctx: mt19937_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}

function  mt19937_long(var ctx: mt19937_ctx): longint;
  {-Next random positive longint}

function  mt19937_dword(var ctx: mt19937_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}

function  mt19937_word(var ctx: mt19937_ctx): word;
  {-Next random word}

function  mt19937_double(var ctx: mt19937_ctx): double;
  {-Next random double [0..1) with 32 bit precision}

function  mt19937_double53(var ctx: mt19937_ctx): double;
  {-Next random double in [0..1) with full double 53 bit precision}

function  mt19937_rangew(var ctx: mt19937_ctx; range: word): word;
  {-Next random word in range 0..range-1}

function  mt19937_rangel(var ctx: mt19937_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}

function  mt19937_selftest: boolean;
  {-Simple self-test of MT19937 PRNG}


implementation

uses
  rndrange;

const
  mag01: array[0..1] of longint = (0, longint($9908b0df));

const
  Upper_Mask = longint($80000000);
  Lower_mask = longint($7fffffff);


{---------------------------------------------------------------------------}
procedure mt19937_init(var ctx: mt19937_ctx; seed: longint);
  {-Init context from seed}
begin
  with ctx do begin
    mt[0]:= seed;
    mti  := 1;
    while mti<N19937 do begin
      {See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier.}
      {In the previous versions, MSBs of the seed affect  }
      {only MSBs of the array mt[].                       }
      {2002/01/09 modified by Makoto Matsumoto            }
      mt[mti] := (1812433253*(mt[mti-1] xor (mt[mti-1] shr 30))+mti);
      inc(mti);
    end;
  end;
end;


{---------------------------------------------------------------------------}
procedure mt19937_init0(var ctx: mt19937_ctx);
  {-Init context from randseed}
begin
  mt19937_init(ctx, randseed);
end;


{---------------------------------------------------------------------------}
{$ifdef CONST}
procedure mt19937_inita(var ctx: mt19937_ctx; const key: array of longint; klen: integer);
  {-Init all context variables with separate seeds, klen: number of seeds}
{$else}
procedure mt19937_inita(var ctx: mt19937_ctx; var KArr; klen: integer);
  {-Init all context variables with separate seeds, klen: number of seeds}
var
  key: mt19937_state absolute KArr; {T5-6 do not have open arrrays}
{$endif}
var
  i,j,k: integer;
begin
  {$ifdef CONST}
    if klen>high(key)+1 then klen := high(key)+1;
  {$endif}
  mt19937_init(ctx, 19650218);
  with ctx do begin
    i := 1;
    j := 0;
    k := klen;
    if N19937>k then k := N19937;
    while k>0 do begin
      mt[i] := (mt[i] xor ((mt[i-1] xor (mt[i-1] shr 30))*1664525))+key[j]+j;
      inc(i);
      inc(j);
      if i>=N19937 then begin
        mt[0] := mt[N19937-1];
        i := 1;
      end;
      if j>=klen then j:=0;
      dec(k);
    end;
    k := N19937-1;
    while k>0 do begin
      mt[i] := (mt[i] xor ((mt[i-1] xor (mt[i-1] shr 30))*1566083941))-i;
      inc(i);
      if i>=N19937 then begin
        mt[0] := mt[N19937-1];
        i := 1;
      end;
      dec(k);
    end;
    {MSB is 1; assuring non-zero initial array}
    mt[0]:= longint($80000000);
  end;
end;


{---------------------------------------------------------------------------}
procedure mt19937_next(var ctx: mt19937_ctx);
  {-Next step of PRNG}
var
  kk: integer;
  y : longint;
begin
  with ctx do begin
    if mti >= N19937 then begin
      {generate N19937 numbers at one time}
      kk := 0;
      while kk<N19937-M19937 do begin
        y := (mt[kk] and Upper_Mask) or (mt[kk+1] and Lower_mask);
        mt[kk] := mt[kk+M19937] xor (y shr 1) xor mag01[y and 1];
        inc(kk);
      end;
      while kk<N19937-1 do begin
        y := (mt[kk] and Upper_Mask) or (mt[kk+1] and Lower_mask);
        mt[kk] := mt[kk+(M19937-N19937)] xor (y shr 1) xor mag01[y and 1];
        inc(kk);
      end;
      y := (mt[N19937-1] and Upper_Mask) or (mt[0] and Lower_mask);
      mt[N19937-1] := mt[M19937-1] xor (y shr 1) xor mag01[y and 1];
      mti := 0;
    end;
    y := mt[mti];
    inc(mti);
    {Tempering}
    {$ifdef BASM16}
      asm
                les di,[ctx]
        db $66; mov dx, word ptr [y]
        db $66; mov ax,dx
        db $66; shr ax,11
        db $66; xor dx,ax               {y := y xor (y shr 11);}

        db $66; mov ax,dx
        db $66; shl ax,7
        db $66; and ax,$5680;dw $9d2c;
        db $66; xor dx,ax               {y := y xor ((y shl 7) and $9d2c5680}

        db $66; mov ax,dx
        db $66; shl ax,15
        db $66; and ax,$0000;dw $efc6;
        db $66; xor dx,ax               {y := y xor ((y shl 15) and $efc60000}

        db $66; mov ax,dx
        db $66; shr ax,18
        db $66; xor ax,dx               {nr:= y xor (y shr 18)}
        db $66; mov word ptr mt19937_ctx(es:[di]).nr,ax
      end;
    {$else}
      y := y xor (y shr 11);
      y := y xor ((y shl  7) and longint($9d2c5680));
      y := y xor ((y shl 15) and longint($efc60000));
      nr:= y xor (y shr 18);
    {$endif}
  end;
end;


{---------------------------------------------------------------------------}
function mt19937_long(var ctx: mt19937_ctx): longint;
  {-Next random positive longint}
begin
  mt19937_next(ctx);
  mt19937_long := ctx.nr shr 1;
end;


{---------------------------------------------------------------------------}
function mt19937_dword(var ctx: mt19937_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}
begin
  mt19937_next(ctx);
  {$ifdef HAS_CARD32}
    mt19937_dword := cardinal(ctx.nr);
  {$else}
    mt19937_dword := ctx.nr;
  {$endif}
end;


{---------------------------------------------------------------------------}
function mt19937_word(var ctx: mt19937_ctx): word;
  {-Next random word}
type
  TwoWords = packed record
               L,H: word
             end;
begin
  mt19937_next(ctx);
  mt19937_word := TwoWords(ctx.nr).H;
end;


{---------------------------------------------------------------------------}
function mt19937_double(var ctx: mt19937_ctx): double;
  {-Next random double [0..1) with 32 bit precision}
begin
  mt19937_next(ctx);
  mt19937_double := (ctx.nr + 2147483648.0) / 4294967296.0;
end;


{---------------------------------------------------------------------------}
function mt19937_double53(var ctx: mt19937_ctx): double;
  {-Next random double in [0..1) with full double 53 bit precision}
var
  hb,lb: longint;
begin
  mt19937_next(ctx);
  hb := ctx.nr shr 5;
  mt19937_next(ctx);
  lb := ctx.nr shr 6;
  mt19937_double53 := (hb*67108864.0+lb)/9007199254740992.0;
end;


{---------------------------------------------------------------------------}
procedure mt19937_read(var ctx: mt19937_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}
type
  plong = ^longint;
begin
  while len>3 do begin
    mt19937_next(ctx);
    plong(dest)^ := ctx.nr;
    inc(Ptr2Inc(dest),4);
    dec(len, 4);
  end;
  if len>0 then begin
    mt19937_next(ctx);
    move(ctx.nr, dest^, len and 3);
  end;
end;


{---------------------------------------------------------------------------}
function mt19937_rangew(var ctx: mt19937_ctx; range: word): word;
  {-Next random word in range 0..range-1}
begin
  mt19937_next(ctx);
  mt19937_rangew := rangew(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function mt19937_rangel(var ctx: mt19937_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}
begin
  mt19937_next(ctx);
  mt19937_rangel := rangel(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function mt19937_selftest: boolean;
  {-Simple self-test of mt19937 PRNG}
var
  ctx: mt19937_ctx;
  i: integer;
const
  STV : array[0..3] of longint = ($123,$234,$345,$456);
  snr = longint($CE3BCD2E);
begin
  mt19937_inita(ctx, STV, 4);
  for i:=1 to 1000 do mt19937_next(ctx);
  mt19937_selftest := ctx.nr=snr;
end;

end.
