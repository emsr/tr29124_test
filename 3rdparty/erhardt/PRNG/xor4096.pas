unit xor4096;

(*************************************************************************

 DESCRIPTION   :  xor4096 based pseudo random number generator
                    Period at least 2**4096-1.

 REQUIREMENTS  :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  [1] xorgens.c version 3.04, R. P. Brent, 20060628 available
                      http://wwwmaths.anu.edu.au/~brent/ftp/random/xorgens304.zip
                  [2] R. P. Brent "Some long-period random number generators using
                      shifts and xors",  Preprint: 28 June 2006, available
                      as rpb224.pdf in xorgens304.zip from [1]

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     09.04.07  W.Ehrhardt  Initial BP7 port of xorgens.c
 0.11     09.04.07  we          xor4096_init calls xor4096_next
 0.12     09.04.07  we          BASM16 version of xor4096_next
 0.13     10.04.07  we          Improved BASM16 version
 0.14     10.04.07  we          Improved Pascal xor4096_next
 0.15     10.04.07  we          Improved BIT16  xor4096_next
 0.16     11.04.07  we          BASM16: removed local Weyl longint
 0.17     17.11.07  we          Corrected typos in selftest description
 0.18     05.11.08  we          xor4096_dword function
 0.19     04.12.08  we          BTypes/Ptr2Inc
 0.20     14.06.12  we          Fix bug in _read for trailing max 3 bytes
 0.21     22.12.17  we          xor4096_rangel/w
 0.22     25.12.17  we          rangel for old non-basm compilers
 **************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2007-2017 Wolfgang Ehrhardt

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
  xor4096_ctx   = record
                    x : array[0..127] of longint; {state vector}
                    i : integer;                  {state index }
                    w : longint;                  {Weyl state  }
                    nr: longint;                  {next random }
                  end;


procedure xor4096_init (var ctx: xor4096_ctx; seed: longint);
  {-Init context from seed}

procedure xor4096_init0(var ctx: xor4096_ctx);
  {-Init context from randseed}

procedure xor4096_next(var ctx: xor4096_ctx);
  {-Next step of PRNG}

procedure xor4096_read(var ctx: xor4096_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}

function  xor4096_long(var ctx: xor4096_ctx): longint;
  {-Next random positive longint}

function  xor4096_dword(var ctx: xor4096_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}

function  xor4096_word(var ctx: xor4096_ctx): word;
  {-Next random word}

function  xor4096_double(var ctx: xor4096_ctx): double;
  {-Next random double [0..1) with 32 bit precision}

function  xor4096_double53(var ctx: xor4096_ctx): double;
  {-Next random double in [0..1) with full double 53 bit precision}

function  xor4096_rangew(var ctx: xor4096_ctx; range: word): word;
  {-Next random word in range 0..range-1}

function  xor4096_rangel(var ctx: xor4096_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}

function  xor4096_selftest: boolean;
  {-Simple self-test of xor4096 PRNG}


implementation

uses
  BTypes, rndrange;


const
  Weyl = $61c88647; {Weyl generator increment}

type
  LH = packed record
         L,H: word
       end;

{---------------------------------------------------------------------------}
procedure xor4096_init(var ctx: xor4096_ctx; seed: longint);
  {-Init context from seed}
var
  k: integer;
  v: longint;
begin
  v := seed;
  {starting value must be non-zero}
  if v=0 then v := not v;
  with ctx do begin
    {Avoid correlations for close seeds, recurrence has period 2^32-1}
    for k:=0 to 31 do begin
      v := v xor (v shl 10);  v := v xor (v shr 15);
      v := v xor (v shl  4);  v := v xor (v shr 13);
    end;
    {Initialize circular array}
    w := v;
    for k:=0 to 127 do begin
      v := v xor (v shl 10);  v := v xor (v shr 15);
      v := v xor (v shl  4);  v := v xor (v shr 13);
      w := w + Weyl;
      x[k] := v + w;
    end;
    {Discard first 512 results}
    i := 127;
    {Because xor4096_next changes w, the current Weyl state is saved and}
    {restored. The original c code doubles the source except that w is  }
    {unchanged. As the init code is normally called only once, the small}
    {speed penalty is no problem here and optimization is done in next. }
    v := w;
    for k:=0 to 511 do xor4096_next(ctx);
    {Restore Weyl state}
    w := v;
  end;
end;


{---------------------------------------------------------------------------}
procedure xor4096_init0(var ctx: xor4096_ctx);
  {-Init context from randseed}
begin
  xor4096_init(ctx, randseed);
end;


{$ifdef BASM16}
{---------------------------------------------------------------------------}
procedure xor4096_next(var ctx: xor4096_ctx);
  {-Next step of PRNG}
begin
  asm
            {save ds and point ds:si to ctx record}
            push ds
            lds  si, [ctx]
            {i := succ(i) and 127;}
            mov  bx, [si].xor4096_ctx.i
            inc  bx
            and  bx,127
            mov  [si].xor4096_ctx.i,bx
            {t := x[i]}
            shl  bx,2
    db $66; mov  cx, word ptr [si].xor4096_ctx.x[bx]
            {t := t xor (t shl a);}
    db $66; mov  ax,cx
    db $66; shl  ax,17
    db $66; xor  cx,ax
            {t := t xor (t shr b);}
    db $66; mov  ax,cx
    db $66; shr  ax,12
    db $66; xor  cx,ax
            {save index i in bx}
            mov  ax,bx
            {v := x[(i+33) and 127];}
            {Note by is already scaled to longint access}
            add  bx,4*33
            and  bx,511
    db $66; mov  dx, word ptr [si].xor4096_ctx.x[bx]
            {restore index i in bx}
            mov  bx,ax
            {v := v xor (v shl c);}
    db $66; mov  ax,dx
    db $66; shl  ax,13
    db $66; xor  dx,ax
            {v := v xor (v shr b);}
    db $66; mov  ax,dx
    db $66; shr  ax,15
    db $66; xor  dx,ax
            {Update circular array}
            {v := v xor t;}
    db $66; xor  dx,cx
            {x[i] := v;}
    db $66; mov  word ptr [si].xor4096_ctx.x[bx],dx
            {Update Weyl generator}
            {w := w + Weyl;}
    db $66; mov  ax, word ptr [si].xor4096_ctx.w
    db $66; add  ax, (Weyl and $ffff); dw (weyl shr 16);
    db $66; mov  word ptr [si].xor4096_ctx.w,ax
            {ax := (w xor (w shr 16)}
    db $66; mov  cx,ax
    db $66; shr  ax,16
    db $66; xor  ax,cx
            {Store next random result}
            {nr := v + (w xor (w shr 16)}
    db $66; add  ax,dx
    db $66; mov  word ptr [si].xor4096_ctx.nr,ax

            {Restore Pascal's ds}
            pop  ds
  end;
end;

{$else}

{---------------------------------------------------------------------------}
procedure xor4096_next(var ctx: xor4096_ctx);
  {-Next step of PRNG}
var
  v,t: longint;
const
  a=17; b=12; c=13; d=15;
begin
  with ctx do begin
    {(I + L^a)(I + R^b) */}
    i := succ(i) and 127;
    t := x[i];
    t := t xor (t shl a);
    t := t xor (t shr b);
    {(I + L^c)(I + R^d) and xor with t}
    v := x[(i+33) and 127];
    v := v xor (v shl c);
    v := (v xor (v shr d)) xor t;
    {Update circular array}
    x[i] := v;
    {Update Weyl generator}
    w := w + Weyl;
    {Store next random result}
    {$ifdef BIT16}
      nr := v + (w xor LH(w).H);
    {$else}
      nr := v + (w xor (w shr 16));
    {$endif}
  end;
end;
{$endif}


{---------------------------------------------------------------------------}
function xor4096_long(var ctx: xor4096_ctx): longint;
  {-Next random positive longint}
begin
  xor4096_next(ctx);
  xor4096_long := ctx.nr shr 1;
end;


{---------------------------------------------------------------------------}
function xor4096_dword(var ctx: xor4096_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}
begin
  xor4096_next(ctx);
  {$ifdef HAS_CARD32}
    xor4096_dword := cardinal(ctx.nr);
  {$else}
    xor4096_dword := ctx.nr;
  {$endif}
end;


{---------------------------------------------------------------------------}
function xor4096_word(var ctx: xor4096_ctx): word;
  {-Next random word}
begin
  xor4096_next(ctx);
  xor4096_word := LH(ctx.nr).H;
end;


{---------------------------------------------------------------------------}
function xor4096_double(var ctx: xor4096_ctx): double;
  {-Next random double [0..1) with 32 bit precision}
begin
  xor4096_next(ctx);
  xor4096_double := (ctx.nr + 2147483648.0) / 4294967296.0;
end;


{---------------------------------------------------------------------------}
function xor4096_double53(var ctx: xor4096_ctx): double;
  {-Next random double in [0..1) with full double 53 bit precision}
var
  hb,lb: longint;
begin
  xor4096_next(ctx);
  hb := ctx.nr shr 5;
  xor4096_next(ctx);
  lb := ctx.nr shr 6;
  xor4096_double53 := (hb*67108864.0+lb)/9007199254740992.0;
end;


{---------------------------------------------------------------------------}
procedure xor4096_read(var ctx: xor4096_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}
type
  plong = ^longint;
begin
  while len>3 do begin
    xor4096_next(ctx);
    plong(dest)^ := ctx.nr;
    inc(Ptr2Inc(dest),4);
    dec(len, 4);
  end;
  if len>0 then begin
    xor4096_next(ctx);
    move(ctx.nr, dest^, len and 3);
  end;
end;


{---------------------------------------------------------------------------}
function xor4096_rangew(var ctx: xor4096_ctx; range: word): word;
  {-Next random word in range 0..range-1}
begin
  xor4096_next(ctx);
  xor4096_rangew := rangew(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function xor4096_rangel(var ctx: xor4096_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}
begin
  xor4096_next(ctx);
  xor4096_rangel := rangel(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function xor4096_selftest: boolean;
  {-Simple self-test of xor4096 PRNG}
var
  i: integer;
  ctx: xor4096_ctx;
const
  t10000 = 820417911; {10000th random number, calculated with t_xor4096.c}
begin
  xor4096_init(ctx, 1);
  for i:=1 to 10000 do xor4096_next(ctx);
  xor4096_selftest := ctx.nr=t10000;
end;

end.
