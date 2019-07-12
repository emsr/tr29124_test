unit TT800;

(*************************************************************************

 DESCRIPTION     :  TT800 based pseudo random number generator
                    Period about 2**800

 REQUIREMENTS    :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA   :  ---

 MEMORY USAGE    :  ---

 DISPLAY MODE    :  ---

 REFERENCES      :  [1] Makoto Matsumoto and Yoshiharu Kurita, "Twisted GFSR
                        Generators II", ACM Transactions on Modelling and Computer
                        Simulation, Vol. 4, No. 3, 1994, pages 254-266.


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     25.12.99  W.Ehrhardt  Initial BP7 port of TT800.C
 0.20     29.05.05  we          uses ctx like Taus88 based PRNG
 0.21     29.05.05  we          prefix tt800_
 0.22     30.05.05  we          selftest with 1000 calls
 0.23     02.08.05  we          Tempering with BASM16
 0.24     05.11.08  we          tt800_dword function
 0.25     04.12.08  we          BTypes/Ptr2Inc
 0.26     14.06.12  we          Fix bug in _read for trailing max 3 bytes
 0.27     21.12.17  we          tt800_rangel/w
 0.28     25.12.17  we          rangel for old non-basm compilers
**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 1999-2017 Wolfgang Ehrhardt

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

{Original comment}

(* A C-program for TT800 : July 8th 1996 Version                    *)
(* by M. Matsumoto, email: matumoto@math.keio.ac.jp                 *)
(* genrand() generate one pseudorandom number with double precision *)
(* which is uniformly distributed on [0,1]-interval                 *)
(* for each call.  One may choose any initial 25 seeds              *)
(* except all zeros.                                                *)

(* See: ACM Transactions on Modelling and Computer Simulation,      *)
(* Vol. 4, No. 3, 1994, pages 254-266.                              *)


interface

{$i std.inc}

{$ifdef BIT16}
  {$N+}
{$endif}

uses
  BTypes;

const
  N800 = 25;
  M800 =  7;

type
  tt800_state = array[0..N800-1] of longint;

  tt800_ctx   = record
                  mt: tt800_state;  {state vector}
                  nr: longint;      {next random }
                  si: integer;      {state index }
                end;

procedure tt800_init (var ctx: tt800_ctx; seed: longint);
  {-Init context from seed}

procedure tt800_init0(var ctx: tt800_ctx);
  {-Init context from randseed}

procedure tt800_inita(var ctx: tt800_ctx; {$ifdef CONST} const {$else} var {$endif} SArr: tt800_state);
  {-Init all context variables with separate seeds}

procedure tt800_next(var ctx: tt800_ctx);
  {-Next step of PRNG}

procedure tt800_read(var ctx: tt800_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}

function  tt800_long(var ctx: tt800_ctx): longint;
  {-Next random positive longint}

function  tt800_dword(var ctx: tt800_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}

function  tt800_word(var ctx: tt800_ctx): word;
  {-Next random word}

function  tt800_double(var ctx: tt800_ctx): double;
  {-Next random double [0..1) with 32 bit precision}

function  tt800_double53(var ctx: tt800_ctx): double;
  {-Next random double in [0..1) with full double 53 bit precision}

function  tt800_rangew(var ctx: tt800_ctx; range: word): word;
  {-Next random word in range 0..range-1}

function  tt800_rangel(var ctx: tt800_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}


function  tt800_selftest: boolean;
  {-Simple self-test of TT800 PRNG}


implementation

uses
  rndrange;

const
  mag01: array[0..1] of longint = (0, longint($8ebfd028));


{---------------------------------------------------------------------------}
procedure tt800_inita(var ctx: tt800_ctx; {$ifdef CONST} const {$else} var {$endif} SArr: tt800_state);
  {-Init all context variables with separate seeds}
var
  i: integer;
  t: longint;
begin
  with ctx do begin
    t := 0;
    for i:=0 to N800-1 do begin
      mt[i] := SArr[i];
      t := t or SArr[i];
    end;
    {Force state to be non-zero}
    if t=0 then mt[0] := 1;
    nr := 0;
    si := 0;
  end;
end;


{---------------------------------------------------------------------------}
procedure tt800_next(var ctx: tt800_ctx);
  {-Next step of PRNG}
var
  kk: integer;
  y : longint;
begin
  with ctx do begin
    if si=N800 then begin
      for kk:=0 to N800-1 do begin
        mt[kk] := mt[(kk+M800) mod N800] xor (mt[kk] shr 1) xor mag01[mt[kk] and 1];
      end;
      si := 0;
    end;
    {Tempering}
    y := mt[si];
    inc(si);
    {$ifdef BASM16}
      asm
                les di,[ctx]
        db $66; mov dx, word ptr [y]

        db $66; mov ax,dx
        db $66; shl ax,7
        db $66; and ax,$2500;dw $2b5b;
        db $66; xor dx,ax               {y := y xor ((y shl 7) and $2b5b2500}

        db $66; mov ax,dx
        db $66; shl ax,15
        db $66; and ax,$0000;dw $db8b;
        db $66; xor dx,ax               {y := y xor ((y shl 15) and $db8b0000}

        db $66; mov ax,dx
        db $66; shr ax,16
        db $66; xor ax,dx               {nr:= y xor (y shr 16)}

        db $66; mov word ptr tt800_ctx(es:[di]).nr,ax
      end;
    {$else}
      y := y xor ((y shl  7) and longint($2b5b2500));
      y := y xor ((y shl 15) and longint($db8b0000));
      nr:= y xor (y shr 16);
    {$endif}
  end;
end;


{---------------------------------------------------------------------------}
procedure tt800_init(var ctx: tt800_ctx; seed: longint);
  {-Init context from seed}
const
  M=69069;
  A=1;
var
  SArr: tt800_state;
  i: integer;
begin
  SArr[0] := seed;
  for i:=1 to N800-1 do begin
    SArr[i] := M*SArr[i-1]+A;
  end;
  tt800_inita(ctx,Sarr);
end;


{---------------------------------------------------------------------------}
procedure tt800_init0(var ctx: tt800_ctx);
  {-Init context from randseed}
begin
  tt800_init(ctx, randseed);
end;


{---------------------------------------------------------------------------}
function tt800_long(var ctx: tt800_ctx): longint;
  {-Next random positive longint}
begin
  tt800_next(ctx);
  tt800_long := ctx.nr shr 1;
end;


{---------------------------------------------------------------------------}
function tt800_dword(var ctx: tt800_ctx): {$ifdef HAS_CARD32}cardinal{$else}longint{$endif};
  {-Next 32 bit random dword (cardinal or longint)}
begin
  tt800_next(ctx);
  {$ifdef HAS_CARD32}
    tt800_dword := cardinal(ctx.nr);
  {$else}
    tt800_dword := ctx.nr;
  {$endif}
end;


{---------------------------------------------------------------------------}
function tt800_word(var ctx: tt800_ctx): word;
  {-Next random word}
type
  TwoWords = packed record
               L,H: word
             end;
begin
  tt800_next(ctx);
  tt800_word := TwoWords(ctx.nr).H;
end;


{---------------------------------------------------------------------------}
function tt800_double(var ctx: tt800_ctx): double;
  {-Next random double [0..1) with 32 bit precision}
begin
  tt800_next(ctx);
  tt800_double := (ctx.nr + 2147483648.0) / 4294967296.0;
end;


{---------------------------------------------------------------------------}
function tt800_double53(var ctx: tt800_ctx): double;
  {-Next random double in [0..1) with full double 53 bit precision}
var
  hb,lb: longint;
begin
  tt800_next(ctx);
  hb := ctx.nr shr 5;
  tt800_next(ctx);
  lb := ctx.nr shr 6;
  tt800_double53 := (hb*67108864.0+lb)/9007199254740992.0;
end;


{---------------------------------------------------------------------------}
procedure tt800_read(var ctx: tt800_ctx; dest: pointer; len: longint);
  {-Read len bytes from the PRNG to dest}
type
  plong = ^longint;
begin
  while len>3 do begin
    tt800_next(ctx);
    plong(dest)^ := ctx.nr;
    inc(Ptr2Inc(dest),4);
    dec(len, 4);
  end;
  if len>0 then begin
    tt800_next(ctx);
    move(ctx.nr, dest^, len and 3);
  end;
end;


{---------------------------------------------------------------------------}
function tt800_rangew(var ctx: tt800_ctx; range: word): word;
  {-Next random word in range 0..range-1}
begin
  tt800_next(ctx);
  tt800_rangew := rangew(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function tt800_rangel(var ctx: tt800_ctx; range: longint): longint;
  {-Next random longint in range 0..range-1}
begin
  tt800_next(ctx);
  tt800_rangel := rangel(ctx.nr, range);
end;


{---------------------------------------------------------------------------}
function tt800_selftest: boolean;
  {-Simple self-test of TT800 PRNG}
var
  ctx: tt800_ctx;
  i: integer;

{$ifdef StrictLong}
  {$warnings off}
  {$R-} {avoid D9 errors!}
{$endif}

{seed for selftest uses original code state vector}
const
  OSeed: tt800_state=($95f24dab, $0b685215, $e76ccae7, $af3ec239, $715fad23,
                      $24a590ad, $69e4b5ef, $bf456141, $96bc1b7b, $a7bdf825,
                      $c1de75b7, $8858a9c9, $2da87693, $b657f9dd, $ffdc8a9f,
                      $8121da71, $8b823ecb, $885d05f5, $4e20cd47, $5a9ad5d9,
                      $512c0c03, $ea857ccd, $4cc1d30f, $8891a8a1, $a6b7aadb);
{$ifdef StrictLong}
  {$warnings on}
  {$ifdef RangeChecks_on}
    {$R+}
  {$endif}
{$endif}

(*  First 50 unsigned long numbers from t_tt800l.c

    3169973338   2724982910    347012937   1735893326   2282497071
    3975116866     62755666    500522132    129776071   1978109378
    4040131704   3800592193   3057303977   1468369496    370579849
    3630178833     51910867    819270944    476180518    190380673
    1370447020   1620916304    663482756   1354889312   4000276916
     868393086   1441698743   1086138563   1899869374   3717419747
    2455034041   2617437696   1595651084   4148285605   1860328467
     928897371    263340857   4091726170   2359987311   1669697327
    1882626857   1635656338    897501559   3233276032    373770970
    2950632840   2706386845   3294066568   3819538748   1902519841

 1000: $1dd4585f
 5000: $c22ff751
10000: $aa4465c3
20000: $3d055cbd
30000: $0196d2ae
*)


begin
  tt800_inita(ctx, OSeed);
  for i:=1 to 1000 do tt800_next(ctx);
  tt800_selftest := ctx.nr=$1dd4585f;
end;

end.
