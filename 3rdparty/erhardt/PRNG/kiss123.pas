unit kiss123;

(*************************************************************************

 DESCRIPTION   :  KISS123 based pseudo random number generator
                  Period about 2**123

 REQUIREMENTS  :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  [1] G. Marsaglia, "Random numbers for C: End, at last?",
                      Usenet post (Jan. 21, 1999) to sci.math, sci.crypt ...
                      Message-ID: <36A5BB98.F2561DFF@stat.fsu.edu>


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     04.02.07  W.Ehrhardt  Initial Pascal version (taus113 layaut)
 0.11     04.02.07  we          BASM16
 0.12     04.02.07  we          Improved BASM16
 0.13     04.02.07  we          Improved BASM16 (calc s3 after s4)
 0.14     04.02.07  we          Simplified selftest
 0.15     10.04.07  we          30% speedup for BIT16 kiss123_next
 0.16     05.11.08  we          kiss123_dword function
 0.17     05.11.08  we          restore lost line in 32 bit kiss123_next
 0.18     04.12.08  we          BTypes/Ptr2Inc
 0.19     14.06.12  we          Fix bug in _read for trailing max 3 bytes
 0.20     22.12.17  we          kiss123_rangel/w
 0.21     25.12.17  we          rangel for old non-basm compilers
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
  kiss123_ctx = record
               nr         : longint;  {next 32 bit rand}
               s1,s2,s3,s4: longint;  {state variables }
             end;

procedure kiss123_init(var ctx: kiss123_ctx; seed: longint);
  {-Init context from seed}

procedure kiss123_init0(var ctx: kiss123_ctx);
  {-Init context from randseed}

procedure kiss123_init4(var ctx: kiss123_ctx; seed1,seed2,seed3,seed4: longint);
  {-Init all context with separate seeds, the initial seeds}
  { seed2, seed3, seed4 should be <> 0 }

procedure kiss123_next(var ctx: kiss123_ctx);
  {-Next step of PRNG}

procedure kiss123_read(var ctx: kiss123_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}

function  kiss123_long(var ctx: kiss123_ctx): longint;
  {-Next random positive longint}

function  kiss123_dword(var ctx: kiss123_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}

function  kiss123_word(var ctx: kiss123_ctx): word;
  {-Next random word}

function  kiss123_double(var ctx: kiss123_ctx): double;
  {-Next random double [0..1) with 32 bit precision}

function  kiss123_double53(var ctx: kiss123_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}

function  kiss123_rangew(var ctx: kiss123_ctx; range: word): word;
  {-Next random word in range 0..range-1}

function  kiss123_rangel(var ctx: kiss123_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}


function  kiss123_selftest: boolean;
  {-Simple self-test of PRNG}


implementation

uses
  BTypes, rndrange;


type
  LH = packed record
         L,H: word
       end;

{G. Marsaglia's original definitions:
#define znew  (z=36969*(z&65535)+(z>>16))
#define wnew  (w=18000*(w&65535)+(w>>16))
#define MWC   ((znew<<16)+wnew )
#define SHR3  (jsr^=(jsr<<17), jsr^=(jsr>>13), jsr^=(jsr<<5))
#define CONG  (jcong=69069*jcong+1234567)
#define KISS  ((MWC^CONG)+SHR3)}

{Mapping to ctx: s1=jcong, s2=jsr, s3=z, s4=w}

{$ifdef BASM16}
{---------------------------------------------------------------------------}
procedure kiss123_next(var ctx: kiss123_ctx);
  {-Next step of PRNG}
const
  C1234567: longint = 1234567;
  C69069  : longint = 69069;
begin
  asm
    db $66; mov  di, word ptr [C1234567]
    db $66; mov  ax, word ptr [C69069]
            push ds
            lds  si, [ctx]
            {s1 := 69069*s1 + 1234567;}
    db $66; mul  word ptr kiss123_ctx[si].s1
    db $66; add  di,ax
    db $66; mov  word ptr kiss123_ctx[si].s1,di
            {nr := s2 xor (s2 shl 17);}
    db $66; mov  ax,word ptr kiss123_ctx[si].s2
    db $66; mov  cx,ax
    db $66; shl  ax,17
    db $66; xor  cx,ax
            {nr := nr xor (nr shr 13);}
    db $66; mov  ax,cx
    db $66; shr  ax,13
    db $66; xor  cx,ax
            {s2 := nr xor (nr shl 5);}
    db $66; mov  ax,cx
    db $66; shl  ax,5
    db $66; xor  cx,ax
    db $66; mov  word ptr kiss123_ctx[si].s2,cx
            {s4 := 18000 * (s4 and $FFFF) + (s4 shr 16);}
            mov  ax,18000
            mul  word ptr kiss123_ctx[si].s4
            add  ax,word ptr kiss123_ctx[si].s4+2
            adc  dx,0
            mov  word ptr kiss123_ctx[si].s4,ax
            mov  word ptr kiss123_ctx[si].s4+2,dx
            {s3 := 36969 * (s3 and $FFFF) + (s3 shr 16);}
            mov  ax,36969
            mul  word ptr kiss123_ctx[si].s3
            add  ax,word ptr kiss123_ctx[si].s3+2
            adc  dx,0
            mov  word ptr kiss123_ctx[si].s3,ax
            mov  word ptr kiss123_ctx[si].s3+2,dx
            {nr := (s3 shl 16) + s4;}
            {low word of s3 is still in ax}
    db $66; shl  ax,16
    db $66; add  ax,word ptr kiss123_ctx[si].s4
            {nr := s2 + (s1 xor nr);}
    db $66; xor  ax,di
    db $66; add  ax,cx
    db $66; mov  word ptr kiss123_ctx[si].nr,ax
            pop  ds
  end;
end;

{$else}

{$ifdef BIT16}

{---------------------------------------------------------------------------}
function mulw(A,B: word): longint;
  {-reverse byte order in longint}
inline(
  $58/          { pop  ax }
  $5A/          { pop  dx }
  $F7/$E2);     { mul  dx }

{---------------------------------------------------------------------------}
procedure kiss123_next(var ctx: kiss123_ctx);
  {-Next step of PRNG}
begin
  with ctx do begin
    s1 := 69069*s1 + 1234567;
    nr := s2 xor (s2 shl 17);
    nr := nr xor (nr shr 13);
    s2 := nr xor (nr shl 5);
    s3 := mulw(36969, LH(s3).L) + LH(s3).H;
    s4 := mulw(18000, LH(s4).L) + LH(s4).H;
    nr := (s3 shl 16) + s4;
    nr := s2 + (s1 xor nr);
  end;
end;

{$else}

{---------------------------------------------------------------------------}
procedure kiss123_next(var ctx: kiss123_ctx);
  {-Next step of PRNG}
begin
  with ctx do begin
    s1 := 69069*s1 + 1234567;
    nr := s2 xor (s2 shl 17);
    nr := nr xor (nr shr 13);
    s2 := nr xor (nr shl 5);
    s3 := 36969 * (s3 and $FFFF) + (s3 shr 16);
    s4 := 18000 * (s4 and $FFFF) + (s4 shr 16);
    nr := (s3 shl 16) + s4;
    nr := s2 + (s1 xor nr);
  end;
end;

{$endif}


{$endif}



{---------------------------------------------------------------------------}
procedure kiss123_init4(var ctx: kiss123_ctx; seed1,seed2,seed3,seed4: longint);
  {-Init all context with separate seeds, the initial seeds}
  { seed2, seed3, seed4 should be <> 0 }
begin
  with ctx do begin
    s1 := seed1;
    s2 := seed2;
    s3 := seed3;
    s4 := seed4;
    if s2=0 then s2:=1;
    if s3=0 then s3:=1;
    if s4=0 then s4:=1;
  end;
end;


{---------------------------------------------------------------------------}
procedure kiss123_init(var ctx: kiss123_ctx; seed: longint);
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
  kiss123_init4(ctx,seed,seed2,seed3,seed4);
end;


{---------------------------------------------------------------------------}
procedure kiss123_init0(var ctx: kiss123_ctx);
  {-Init context from randseed}
begin
  kiss123_init(ctx, randseed);
end;


{---------------------------------------------------------------------------}
function kiss123_long(var ctx: kiss123_ctx): longint;
  {-Next random positive longint}
begin
  kiss123_next(ctx);
  kiss123_long := ctx.nr shr 1;
end;


{---------------------------------------------------------------------------}
function kiss123_dword(var ctx: kiss123_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}
begin
  kiss123_next(ctx);
  {$ifdef HAS_CARD32}
    kiss123_dword := cardinal(ctx.nr);
  {$else}
    kiss123_dword := ctx.nr;
  {$endif}
end;


{---------------------------------------------------------------------------}
function kiss123_word(var ctx: kiss123_ctx): word;
  {-Next random word}
begin
  kiss123_next(ctx);
  kiss123_word := LH(ctx.nr).H;
end;


{---------------------------------------------------------------------------}
function kiss123_double(var ctx: kiss123_ctx): double;
  {-Next random double [0..1) with 32 bit precision}
begin
  kiss123_next(ctx);
  kiss123_double := (ctx.nr + 2147483648.0) / 4294967296.0;
end;


{---------------------------------------------------------------------------}
function kiss123_double53(var ctx: kiss123_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}
var
  hb,lb: longint;
begin
  kiss123_next(ctx);
  hb := ctx.nr shr 5;
  kiss123_next(ctx);
  lb := ctx.nr shr 6;
  kiss123_double53 := (hb*67108864.0+lb)/9007199254740992.0;
end;


{---------------------------------------------------------------------------}
procedure kiss123_read(var ctx: kiss123_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}
type
  plong = ^longint;
begin
  while len>3 do begin
    kiss123_next(ctx);
    plong(dest)^ := ctx.nr;
    inc(Ptr2Inc(dest),4);
    dec(len, 4);
  end;
  if len>0 then begin
    kiss123_next(ctx);
    move(ctx.nr, dest^, len and 3);
  end;
end;


{---------------------------------------------------------------------------}
function kiss123_rangew(var ctx: kiss123_ctx; range: word): word;
  {-Next random word in range 0..range-1}
begin
  kiss123_next(ctx);
  kiss123_rangew := rangew(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function kiss123_rangel(var ctx: kiss123_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}
begin
  kiss123_next(ctx);
  kiss123_rangel := rangel(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function kiss123_selftest: boolean;
  {-Simple self-test of PRNG}
var
  ctx: kiss123_ctx;
  i: longint;
begin
  {Values from G. Marsaglia's table setup for test cases}
  kiss123_init4(ctx,12345,34221,12345,65435);
  {256 steps in table setup + 1000000 standard kiss steps}
  for i:=1 to 1000256 do kiss123_next(ctx);
  kiss123_selftest := ctx.nr=longint(1372460312);
end;

end.
