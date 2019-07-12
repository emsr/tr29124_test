unit salsar;

(*************************************************************************

 DESCRIPTION   :  o Cryptographic pseudo random number generator based on
                    Salsa20 encryption function by D.J. Bernstein
                  o Period about 2**70 bytes
                  o uses salsa unit

 REQUIREMENTS  :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  [1] Original version for ECRYPT Stream Cipher Project
                      http://www.ecrypt.eu.org/stream/ciphers/salsa20/salsa20.zip
                      http://www.ecrypt.eu.org/stream/ciphers/salsa20/salsa20source.zip
                  [2] salsa20,12,8-ref.c version 20060209, D.J. Bernstein, public domain
                      http://cr.yp.to/streamciphers/submissions.tar.gz
                  [3] Snuffle 2005: the Salsa20 encryption function:
                      http://cr.yp.to/snuffle.html

 REMARKS       :  1. The RECOMMENDED init procedure is salsar_inita, this
                     covers the full key/IV range
                  2. salsar_init uses only max 32 independent seed bits !!!!
                  3. salsar_init0 may use even less than 32 bits        !!!!
                     (depends on compiler initialisation of randseed)   !!!!
                  4. Default number of rounds is 8, can be changed with salsar_set_rounds

 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     23.04.06  W.Ehrhardt  Initial BP7 version based on aesr
 0.11     23.04.06  we          Selftest for SR_Rounds=8,12,20
 0.12     24.04.06  we          Use 128 bit keys
 0.13     24.04.06  we          salsar_set_rounds
 0.14     24.04.06  we          salsar_get_rounds
 0.15     05.11.08  we          salsar_dword function
 0.16     14.06.12  we          Fix bug in _read for trailing max 3 bytes
 0.17     22.12.17  we          salsar_rangel/w
 0.18     25.12.17  we          rangel for old non-basm compilers
**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2006-2017 Wolfgang Ehrhardt

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

uses salsa20;

const
  SR_KSB_SIZE = 32;   {key stream buffer size in longints}
  SSARR_CNT   =  6;   {Number of longints in seed array, = 10 for 256 bit keys}

type
  salsar_sarr = array[0..SSARR_CNT-1] of longint;   {Seed array: 128 bit key, 64 bit IV}
  salsar_kstr = array[0..SR_KSB_SIZE-1] of longint; {key stream buffer type}

  salsar_ctx  = record
                cc: salsa_ctx;   {stream cipher context}
                ks: salsar_kstr; {key streean buffer   }
                nr: longint;     {next random result   }
                ib: integer;     {index into buffer }
              end;

procedure salsar_inita(var ctx: salsar_ctx; {$ifdef CONST} const {$else} var {$endif} SArr: salsar_sarr);
  {-Init all context variables with separate seeds}

procedure salsar_init(var ctx: salsar_ctx; seed: longint);
  {-Init context from seed}

procedure salsar_init0(var ctx: salsar_ctx);
  {-Init context from randseed}

procedure salsar_next(var ctx: salsar_ctx);
  {-Next step of PRNG}

procedure salsar_read(var ctx: salsar_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}

function  salsar_long(var ctx: salsar_ctx): longint;
  {-Next random positive longint}

function  salsar_dword(var ctx: salsar_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}

function  salsar_word(var ctx: salsar_ctx): word;
  {-Next random word}

function  salsar_double(var ctx: salsar_ctx): double;
  {-Next random double [0..1) with 32 bit precision}

function  salsar_double53(var ctx: salsar_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}

function  salsar_rangew(var ctx: salsar_ctx; range: word): word;
  {-Next random word in range 0..range-1}

function  salsar_rangel(var ctx: salsar_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}

function  salsar_get_rounds: word;
  {-Return the current used number of salsa rounds (8,12,20)}

procedure salsar_set_rounds(rounds: word);
  {-Set salsa rounds to 8,12,20; default = 8}

function  salsar_selftest: boolean;
  {-Simple self-test of PRNG}


implementation

uses
  BTypes, rndrange;

var
  SR_Rounds: word;


{---------------------------------------------------------------------------}
procedure salsar_next(var ctx: salsar_ctx);
  {-Next step of PRNG}
begin
  {call basic Salsa20 keystream generator}
  with ctx do begin
    if ib >= SR_KSB_SIZE then begin
      salsa_keystream_bytes(cc, @ks, sizeof(ks));
      ib := 0;
    end;
    nr := ks[ib];
    inc(ib);
  end;
end;


{---------------------------------------------------------------------------}
procedure salsar_inita(var ctx: salsar_ctx; {$ifdef CONST} const {$else} var {$endif} SArr: salsar_sarr);
  {-Init all context variables with separate seeds}
begin
  with ctx do begin
    salsa_xkeysetup(cc, @Sarr, 128, SR_Rounds);
    salsa_ivsetup(cc, @Sarr[SSARR_CNT-2]);
    ib := SR_KSB_SIZE;
  end;
end;


{---------------------------------------------------------------------------}
procedure salsar_init(var ctx: salsar_ctx; seed: longint);
  {-Init context from seed}
const
  M=69069;
  A=1;
var
  SArr: salsar_sarr;
  i: integer;
begin
  SArr[0] := seed;
  {Use simple LCG for next seeds}
  for i:=1 to SSARR_CNT-1 do SArr[i] := M*SArr[i-1] + A;
  salsar_inita(ctx,SArr);
end;


{---------------------------------------------------------------------------}
procedure salsar_init0(var ctx: salsar_ctx);
  {-Init context from randseed}
begin
  salsar_init(ctx, randseed);
end;


{---------------------------------------------------------------------------}
function salsar_long(var ctx: salsar_ctx): longint;
  {-Next random positive longint}
begin
  salsar_next(ctx);
  {make positive, highest bit=0}
  salsar_long := ctx.nr shr 1;
end;


{---------------------------------------------------------------------------}
function salsar_dword(var ctx: salsar_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}
begin
  salsar_next(ctx);
  {$ifdef HAS_CARD32}
    salsar_dword := cardinal(ctx.nr);
  {$else}
    salsar_dword := ctx.nr;
  {$endif}
end;


{---------------------------------------------------------------------------}
function salsar_word(var ctx: salsar_ctx): word;
  {-Next random word}
type
  TwoWords = packed record
               L,H: word
             end;
begin
  salsar_next(ctx);
  salsar_word := TwoWords(ctx.nr).H;
end;


{---------------------------------------------------------------------------}
function salsar_double(var ctx: salsar_ctx): double;
  {-Next random double [0..1) with 32 bit precision}
begin
  salsar_next(ctx);
  salsar_double := (ctx.nr + 2147483648.0) / 4294967296.0;
end;


{---------------------------------------------------------------------------}
function salsar_double53(var ctx: salsar_ctx): double;
  {-Next random double in [0..1) with 53 bit precision}
var
  hb,lb: longint;
begin
  salsar_next(ctx);
  hb := ctx.nr shr 5;
  salsar_next(ctx);
  lb := ctx.nr shr 6;
  salsar_double53 := (hb*67108864.0+lb)/9007199254740992.0;
end;


{---------------------------------------------------------------------------}
procedure salsar_read(var ctx: salsar_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}
begin
  with ctx do begin
    {use "native" salsa__keystream_bytes}
    salsa_keystream_bytes(cc, dest, len);
    {Reset index into buffer}
    ib := SR_KSB_SIZE;
  end;
end;


{---------------------------------------------------------------------------}
function salsar_rangew(var ctx: salsar_ctx; range: word): word;
  {-Next random word in range 0..range-1}
begin
  salsar_next(ctx);
  salsar_rangew := rangew(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function salsar_rangel(var ctx: salsar_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}
begin
  salsar_next(ctx);
  salsar_rangel := rangel(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
procedure salsar_set_rounds(rounds: word);
  {-Set salsa rounds to 8,12,20; default = 8}
begin
  if (rounds=12) or (rounds=20) then SR_Rounds := rounds
  else SR_Rounds := 8;
end;


{---------------------------------------------------------------------------}
function salsar_get_rounds: word;
  {-Return the current used number of salsa rounds (8,12,20)}
begin
  salsar_get_rounds := SR_Rounds;
end;


{---------------------------------------------------------------------------}
function salsar_selftest: boolean;
  {-Simple self-test of PRNG}
var
  SArr: salsar_sarr;
  ctx : salsar_ctx;
  i,k : integer;
(* 256
const
  tv: array[0..5] of longint =
        { set 1, vector# 0     set 5, vector# 0}
        (longint($d1a8069a),  longint($046c0b0b),     { 8 rounds}
         longint($5fcbf076),  longint($92945bf4),     {12 rounds}
         longint($7784bc64),  longint($22aeb9fb));    {20 rounds}
*)
const
  tv: array[0..5] of longint =
        { set 1, vector# 0     set 5, vector# 0}
        (longint($2002774d),  longint($c00878d2),     { 8 rounds}
         longint($a46bcc0d),  longint($bfebbf90),     {12 rounds}
         longint($0397be63),  longint($9d8e1946));    {20 rounds}

begin
  salsar_selftest := false;
  for k:=0 to 1 do begin
    {Use test data from sets 1/5, vector# 0}
    fillchar(SArr, sizeof(SArr), 0);
    SArr[(SSARR_CNT-2)*k] := 128;      {either key[0]:=$80 or IV[0] := $80}
    salsar_inita(ctx,SArr);
    {Get 128 x 32 Bits}
    for i:=1 to 128 do salsar_next(ctx);
    {Check against stream bytes [508..511]}
    if SR_Rounds=8 then i:=0
    else if SR_Rounds=12 then i:=2
    else i:=4;
    if ctx.nr<>tv[k+i] then exit;
  end;
  salsar_selftest := true;
end;

begin
  SR_Rounds := 8;
end.
