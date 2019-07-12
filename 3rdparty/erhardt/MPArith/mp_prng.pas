unit mp_prng;

{MP interface to (C)PRNG, functions for PRNG generation}

interface

{$i STD.INC}

uses mp_types;

{$i mp_conf.inc}

(*************************************************************************

 DESCRIPTION     :  MP interface to (C)PRNG, functions for PRNG generation,
                    generators are plugged in via $ifdefs/mp_conf

 REQUIREMENTS    :  BP7, D1-D7/D9-D10/D12/D17-D18, FPC, VP

 EXTERNAL DATA   :  (mp_types)

 MEMORY USAGE    :  heap

 DISPLAY MODE    :  ---

 REFERENCES      :  ---


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------

 0.1.00   26.08.05  WEhrhardt   initial BP7 random only version
 0.1.01   26.08.05  we          bugfix mp_random_radix
 0.1.02   26.08.05  we          fix for BIT32 and FPC
 0.1.03   01.09.05  we          ISAAC routines
 0.1.04   14.09.05  we          mp_random_seed
 0.1.05   12.11.05  we          mp_random_byte
 0.1.06   13.11.05  we          mp_random_read
 0.1.07   20.11.05  we          fix mp_random_digit for Delphi/DIGIT_BIT=31
 0.1.08   22.11.05  we          mp_random_randomize
 0.1.09   31.12.05  we          mp_random_word
 0.1.10   05.08.06  we          MPC_UseTSC implemented
 0.1.11   12.08.06  we          Added missing mp_random_word if no ISAAC
 0.1.12   09.09.06  we          Added some type casts type(random())
 0.1.13   11.03.07  we          Generators via inc files, taus88
 0.1.14   11.03.07  we          Removed RTL random - FPC 2.0.4: incompatible
                                randseed handling, others: non-random low bits
 0.1.15   10.07.07  we          mp_random_seed1, removed small bias in mp_random_radix
 0.1.16   20.08.07  we          changed mp_random_int, new mp_random_long

 1.14.00  13.02.10  we          MPC_MAXRadix64 adjustment

**************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2005-2015 Wolfgang Ehrhardt

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


function  mp_random_byte: byte;
  {-Return a random byte}

function  mp_random_digit: mp_digit;
  {-Return a random mp_digit}

function  mp_random_int: longint;
  {-Return a random signed longint}

function  mp_random_long: longint;
  {-Return a random positive longint}

function  mp_random_radix(radix: word): word;
  {-Return a random word in the range 0..radix-1}

procedure mp_random_randomize;
  {-Initialize PRNG via randomize/TSC}

procedure mp_random_read(dest: pointer; len: word);
  {-Read len bytes from the PRNG to dest}

procedure mp_random_seed(const seed: array of longint);
  {-Initialize PRNG with array of longint}

procedure mp_random_seed1(seed: longint);
  {-Initialize PRNG with a longint}

function  mp_random_word: word;
  {-Return a random word}


implementation

{$ifdef MPC_UseISAAC}
  {$i mp_isaac.inc}
{$else}
  {$i mp_taus.inc}
{$endif}


{$ifdef MPC_ACCEPTBIAS}  {included in mp_conf.inc V0.42+}

{---------------------------------------------------------------------------}
function  mp_random_radix(radix: word): word;
  {-Return a random word in the range 0..radix-1}
begin
  if radix<2 then mp_random_radix := 0
  else mp_random_radix := mp_random_word mod radix;
end;

{$else}

const
  rul: array[2..MaxRadix] of word =  {radix*($FFFF div radix)}
         (65534, 65535, 65532, 65535, 65532, 65534, 65528, 65529,
          65530, 65527, 65532, 65533, 65534, 65535, 65520, 65535,
          65520, 65531, 65520, 65520, 65516, 65527, 65520, 65525,
          65520, 65529, 65520, 65511, 65520, 65534, 65504, 65505,
          65518, 65520, 65520
        {$ifdef MPC_MAXRadix64},
          65527, 65512, 65520, 65520, 65518, 65520, 65532, 65516,
          65520, 65504, 65518, 65520, 65513, 65500, 65535, 65520,
          65508, 65502, 65505, 65520, 65493, 65482, 65490, 65520,
          65514, 65534, 65520, 65472);
        {$else}
          );
        {$endif}

var
  lastrad: word;  {last used radix}
  lastrul: word;  {lastrad*($FFFF div lastrad)}

{---------------------------------------------------------------------------}
function  mp_random_radix(radix: word): word;
  {-Return a random word in the range 0..radix-1}
var
  w,u: word;
begin
  if radix<2 then mp_random_radix := 0
  else begin
    {mp_random_word mod radix produces a small theoretical bias for}
    {uniform mp_random_word. Repeat calls to mp_random_word until  }
    {the result is less than radix*($FFFF div radix) than take mod.}
    if radix<=MaxRadix then u := rul[radix]
    else if radix=lastrad then u := lastrul
    else begin
      u :=  radix*($FFFF div radix);
      lastrul := u;
      lastrad := radix;
    end;
    repeat
      w := mp_random_word;
    until w<u;
    mp_random_radix := w mod radix;
  end;
end;

{$endif}


{---------------------------------------------------------------------------}
procedure mp_random_seed1(seed: longint);
  {-Initialize PRNG with a longint}
begin
  mp_random_seed(seed);
end;


begin
  lastrad := 0;
  {$ifdef MPC_Randomize}
    {randomize generator at startup}
    mp_random_randomize;
  {$else}
    {start generator with known value}
    mp_random_seed1(0);
  {$endif}
end.

