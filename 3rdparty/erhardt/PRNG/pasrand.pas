unit pasrand;

(*************************************************************************

 DESCRIPTION   :  Pascal/Delphi compatible LC generator
                  Period about 2**32

 REQUIREMENTS  :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  ---


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     26.12.17  W.Ehrhardt  Initial BP7 version from taus88
 0.11     26.12.17  we          Selftest from 500 Pascal/Delphi randoms
 0.12     27.12.17  we          BASM16 in pasrand_next
**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2017 Wolfgang Ehrhardt

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
  pasrand_ctx = record
                  nr: longint;  {next rand}
                end;

procedure pasrand_init(var ctx: pasrand_ctx; seed: longint);
  {-Init context from seed}

procedure pasrand_init0(var ctx: pasrand_ctx);
  {-Init context from randseed}

procedure pasrand_next(var ctx: pasrand_ctx);
  {-Next step of PRNG}

procedure pasrand_read(var ctx: pasrand_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}

function  pasrand_long(var ctx: pasrand_ctx): longint;
  {-Next random positive longint}

function  pasrand_dword(var ctx: pasrand_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}

function  pasrand_word(var ctx: pasrand_ctx): word;
  {-Next random word}

function  pasrand_double(var ctx: pasrand_ctx): double;
  {-Next random double [0..1) with 32 bit precision}

function  pasrand_double53(var ctx: pasrand_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}

function  pasrand_rangew(var ctx: pasrand_ctx; range: word): word;
  {-Next random word in range 0..range-1}

function  pasrand_rangel(var ctx: pasrand_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}

function  pasrand_selftest: boolean;
  {-Simple self-test of PRNG}


implementation

uses
  BTypes, rndrange;


{---------------------------------------------------------------------------}
procedure pasrand_next(var ctx: pasrand_ctx);
  {-Next step of PRNG}
begin
  {TP and Delphi compatible LCG}
{$ifdef BASM16}
  asm
    db $66,$BA,$05,$84,$08,$08    {mov edx,$8088405}
             les  si,[ctx]
    db $66;  mov  ax, word ptr es:[si].pasrand_ctx.nr
    db $66;  imul dx
    db $66;  inc  ax
    db $66;  mov  word ptr es:[si].pasrand_ctx.nr,ax
  end;
{$else}
  with ctx do nr := $8088405*nr + 1;
{$endif}

end;


{---------------------------------------------------------------------------}
procedure pasrand_init(var ctx: pasrand_ctx; seed: longint);
  {-Init context from seed}
begin
  ctx.nr := seed;
end;


{---------------------------------------------------------------------------}
procedure pasrand_init0(var ctx: pasrand_ctx);
  {-Init context from randseed}
begin
  pasrand_init(ctx, randseed);
end;


{---------------------------------------------------------------------------}
function pasrand_long(var ctx: pasrand_ctx): longint;
  {-Next random positive longint}
begin
  pasrand_next(ctx);
  pasrand_long := ctx.nr shr 1;
end;


{---------------------------------------------------------------------------}
function pasrand_dword(var ctx: pasrand_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}
begin
  pasrand_next(ctx);
  {$ifdef HAS_CARD32}
    pasrand_dword := cardinal(ctx.nr);
  {$else}
    pasrand_dword := ctx.nr;
  {$endif}
end;


{---------------------------------------------------------------------------}
function pasrand_word(var ctx: pasrand_ctx): word;
  {-Next random word}
type
  TwoWords = packed record
               L,H: word
             end;
begin
  pasrand_next(ctx);
  pasrand_word := TwoWords(ctx.nr).H;
end;


{---------------------------------------------------------------------------}
function pasrand_double(var ctx: pasrand_ctx): double;
  {-Next random double [0..1) with 32 bit precision}
begin
  pasrand_next(ctx);
  pasrand_double := (ctx.nr + 2147483648.0) / 4294967296.0;
end;


{---------------------------------------------------------------------------}
function pasrand_double53(var ctx: pasrand_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}
var
  hb,lb: longint;
begin
  pasrand_next(ctx);
  hb := ctx.nr shr 5;
  pasrand_next(ctx);
  lb := ctx.nr shr 6;
  pasrand_double53 := (hb*67108864.0+lb)/9007199254740992.0;
end;


{---------------------------------------------------------------------------}
procedure pasrand_read(var ctx: pasrand_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}
type
  plong = ^longint;
begin
  while len>3 do begin
    pasrand_next(ctx);
    plong(dest)^ := ctx.nr;
    inc(Ptr2Inc(dest),4);
    dec(len, 4);
  end;
  if len>0 then begin
    pasrand_next(ctx);
    move(ctx.nr, dest^, len and 3);
  end;
end;


{---------------------------------------------------------------------------}
function pasrand_rangew(var ctx: pasrand_ctx; range: word): word;
  {-Next random word in range 0..range-1}
begin
  pasrand_next(ctx);
  pasrand_rangew := rangew(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function pasrand_rangel(var ctx: pasrand_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}
begin
  pasrand_next(ctx);
  pasrand_rangel := rangel(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function pasrand_selftest: boolean;
  {-Simple self-test of PRNG}
var
  ctx: pasrand_ctx;
  i: integer;
begin
  pasrand_init(ctx,12345);
  for i:=1 to 500 do pasrand_next(ctx);
  {Test values from BP7..D18}
  pasrand_selftest := ctx.nr=longint($529B5315);
end;

end.
