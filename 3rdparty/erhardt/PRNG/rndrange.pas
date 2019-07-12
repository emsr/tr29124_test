unit rndrange;

(*************************************************************************

 DESCRIPTION   :  Transform 32-bit random number to 0..range-1 like random(range)}

 REQUIREMENTS  :  TP5-7, D1-D7/D9-D10/D12/D17-D18/D25S, FPC, VP

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  ---


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     21.12.17  W.Ehrhardt  Initial BP7 version
 0.11     21.12.17  we          Other compilers
 0.12     21.12.17  we          single line for HAS_UINT64
 0.13     21.12.17  we          inline for HAS_INLINE
 0.14     22.12.17  we          abs(range) in rangel
 0.15     25.12.17  we          rangel/w for old non-basm compilers
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

function rangew(rnd: longint; range: word): word; {$ifdef HAS_INLINE}inline;{$endif}
  {-Transform 32-bit rnd to 0..range-1 like random(range)}

function rangel(rnd, range: longint): longint; {$ifdef HAS_INLINE}inline;{$endif}
  {-Transform 32-bit rnd to 0..range-1 like random(range)}


implementation

uses
  BTypes;

{$ifdef HAS_UINT64}
  {--------------------------------------------------}
  function rangew(rnd: longint; range: word): word;{$ifdef HAS_INLINE}inline;{$endif}
  begin
    rangew := (uint32(rnd)*uint64(range) shr 32) and $ffff;
  end;

  {--------------------------------------------------}
  function rangel(rnd, range: longint): longint;{$ifdef HAS_INLINE}inline;{$endif}
  begin
    rangel := longint(uint32(rnd)*uint64(abs(range)) shr 32);
  end;

{$else}
  {$ifdef BIT32}
    {--------------------------------------------------}
    function rangew(rnd: longint; range: word): word;
    begin
      asm
        mov   eax,[rnd]
        movzx edx,[range]
        mul   edx
        mov   @result,dx
      end
    end;
    {--------------------------------------------------}
    function rangel(rnd, range: longint): longint;
    begin
      asm
        mov   eax,[range]
        cdq
        xor   eax,edx
        sub   eax,edx
        {here eax = abs(range)}
        mul   [rnd]
        mov   @result,edx
      end
    end;
  {$else}
    {$ifdef BASM16}
      {--------------------------------------------------}
      function rangew(rnd: longint; range: word): word;
      begin
        asm
          db $66; mov ax,word ptr [rnd]
          db $66; xor dx,dx
                  mov dx,word ptr [range]
          db $66; mul dx
                  mov @result,dx
        end
      end;
      {--------------------------------------------------}
      function rangel(rnd, range: longint): longint;
      begin
        asm
          db $66; mov ax,word ptr [range]
          db $66,$99; {cdq}
          db $66; xor ax,dx
          db $66; sub ax,dx
          {eax = abs(range)}
          db $66; mul word ptr [rnd]
          db $66; mov word ptr[@result], dx
        end
      end;
    {$else}
      {--------------------------------------------------}
      function mullh(x,y: longint): longint;
        {-Multiply x*|y| and return the high 32-bit of the 64-bit result}
      var
        z: array[0..1] of longint;
        w: array[0..7] of byte absolute z;
        u: array[0..3] of byte absolute x;
        v: array[0..3] of byte absolute y;
        t,k: word;
        i,j: integer;
      const
        m=4;
        n=4;
      begin
        {Schoolbook method: Knuth's algorithm M}
        z[0] := 0;
        z[1] := 0;
        y := abs(y);
        for j:=0 to n-1 do begin
          if v[j]=0 then w[j+m] := 0
          else begin
            k := 0;
            for i:=0 to m-1 do begin
              t := k + word(u[i])*word(v[j]) + w[i+j];
              w[i+j] := byte(t);
              k := t div 256;
            end;
            w[j+m] := k;
          end;
        end;
        mullh := z[1];
      end;
      {--------------------------------------------------}
      function rangew(rnd: longint; range: word): word;
      begin
        rangew := word(mullh(rnd, word(range)));
      end;
      {--------------------------------------------------}
      function rangel(rnd, range: longint): longint;
      begin
        rangel := mullh(rnd, range);
      end;
    {$endif}
  {$endif}
{$endif}

end.
