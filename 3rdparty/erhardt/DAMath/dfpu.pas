unit DFPU;

{Simple FPU unit for DAMath: Rounding and Precision, Get/SetExceptionMask}

interface

{$i STD.INC}

{$ifdef BIT16}
{$N+,F+}
{$endif}


(*************************************************************************

 DESCRIPTION   :  Simple FPU unit for DAMath: Rounding and Precision, Get/SetExceptionMask

 REQUIREMENTS  :  BP7, D2-D7/D9-D10/D12/D17-D18/D25, FPC, VP, WDOSX

 EXTERNAL DATA :  ---

 MEMORY USAGE  :  ---

 DISPLAY MODE  :  ---

 REFERENCES    :  ---


 Version  Date      Author      Modification
 -------  --------  -------     ------------------------------------------
 0.10     28.07.17  W.Ehrhardt  Initial version
 0.11     28.07.17  we          UNIT_SCOPE
 0.12     29.07.17  we          look-up arrays
 0.13     29.07.17  we          look-up arrays dont work for FPC 64-bit
 0.14     29.07.17  we          Name constants for rounding/precision modes
 0.15     29.07.17  we          Fix GetPrecisionMode
 0.16     29.07.17  we          sizeof asserts in debug mode
 0.17     29.07.17  we          Changed type of names to string[10]
 0.18     22.08.18  we          Get/SetExceptionMask

***************************************************************************)


(*-------------------------------------------------------------------------
 (C) Copyright 2017-2018 Wolfgang Ehrhardt

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

{$define fpu_usemath}

{$ifdef BIT16}
  {$undef fpu_usemath}
{$else}
  {$ifdef VirtualPascal}
    {$undef fpu_usemath}
  {$endif}

  {$ifdef FPC}
    {$ifdef VER1}
      {$undef fpu_usemath}
    {$endif}
  {$endif}

  {$ifdef Delphi}
    {$ifndef CONDITIONALEXPRESSIONS}
      {$undef fpu_usemath}
    {$endif}
  {$endif}
{$endif}


type
  TFPURoundingMode  = (rmNearest, rmDown, rmUp, rmTruncate);

type
  TFPUPrecisionMode = (pmSingle, pmReserved, pmDouble, pmExtended);

const
  RoundingNames:  array[TFPURoundingMode]  of string[10] = ('rmNearest', 'rmDown', 'rmUp', 'rmTruncate');
  PrecisionNames: array[TFPUPrecisionMode] of string[10] = ('pmSingle', 'pmReserved', 'pmDouble', 'pmExtended');

function  GetRoundMode: TFPURoundingMode;
  {-Return the current rounding mode}

function  SetRoundMode(NewRoundMode: TFPURoundingMode): TFPURoundingMode;
  {-Set new rounding mode and return the old mode}

function  GetPrecisionMode: TFPUPrecisionMode;
  {-Return the current precision control mode}

function  SetPrecisionMode(NewPrecision: TFPUPrecisionMode): TFPUPrecisionMode;
  {-Set new precision control mode and return the old precision}


const
  {Bitmasks of single exceptions}
  bexInvalidOp    = $01;
  bexDenormalized = $02;
  bexZeroDivide   = $04;
  bexOverflow     = $08;
  bexUnderflow    = $10;
  bexPrecision    = $20;

  bexAllArithmeticExceptions = $3F;


procedure GetExceptionMask(var Mask: byte);
  {-Return the current exception mask}

procedure SetExceptionMask(NewMask: byte);
  {-Set new exception mask}


implementation


{$ifdef fpu_usemath}

{$ifdef UNIT_SCOPE}
uses
  system.math;

{---------------------------------------------------------------------------}
function  GetRoundMode: TFPURoundingMode;
  {-Return the current rounding mode}
begin
  GetRoundMode := TFPURoundingMode(system.math.GetRoundMode);
end;


{---------------------------------------------------------------------------}
function  SetRoundMode(NewRoundMode: TFPURoundingMode): TFPURoundingMode;
  {-Set new rounding mode and return the old mode}
begin
  SetRoundMode := TFPURoundingMode(system.math.SetRoundMode(system.math.TFPURoundingMode(NewRoundMode)));
end;


{---------------------------------------------------------------------------}
function  GetPrecisionMode: TFPUPrecisionMode;
  {-Return the current precision control mode}
begin
  GetPrecisionMode := TFPUPrecisionMode(system.math.GetPrecisionMode);
end;


{---------------------------------------------------------------------------}
function  SetPrecisionMode(NewPrecision: TFPUPrecisionMode): TFPUPrecisionMode;
  {-Set new precision control mode and return the old precision}
begin
  SetPrecisionMode := TFPUPrecisionMode(system.math.SetPrecisionMode(system.math.TFPUPrecisionMode(NewPrecision)));
end;


{---------------------------------------------------------------------------}
procedure GetExceptionMask(var Mask: byte);
  {-Return the current exception mask}
var
  MaskSet: TFPUExceptionMask;
begin
  MaskSet := System.Math.GetExceptionMask;
  Mask := 0;
  if exInvalidOp    in MaskSet then Mask := Mask or bexInvalidOp;
  if exDenormalized in MaskSet then Mask := Mask or bexDenormalized;
  if exZeroDivide   in MaskSet then Mask := Mask or bexZeroDivide;
  if exOverflow     in MaskSet then Mask := Mask or bexOverflow;
  if exUnderflow    in MaskSet then Mask := Mask or bexUnderflow;
  if exPrecision    in MaskSet then Mask := Mask or bexPrecision;
end;

procedure SetExceptionMask(NewMask: byte);
  {-Set new exception mask}
var
  MaskSet: TFPUExceptionMask;
begin
  MaskSet := [];
  if NewMask and bexInvalidOp    <> 0 then MaskSet := MaskSet + [exInvalidOp];
  if NewMask and bexDenormalized <> 0 then MaskSet := MaskSet + [exDenormalized];
  if NewMask and bexZeroDivide   <> 0 then MaskSet := MaskSet + [exZeroDivide];
  if NewMask and bexOverflow     <> 0 then MaskSet := MaskSet + [exOverflow];
  if NewMask and bexUnderflow    <> 0 then MaskSet := MaskSet + [exUnderflow];
  if NewMask and bexPrecision    <> 0 then MaskSet := MaskSet + [exPrecision];
  System.Math.SetExceptionMask(MaskSet);
end;


{$else}

uses
  math;

{---------------------------------------------------------------------------}
function  GetRoundMode: TFPURoundingMode;
  {-Return the current rounding mode}
begin
  GetRoundMode := TFPURoundingMode(math.GetRoundMode);
end;


{---------------------------------------------------------------------------}
function  SetRoundMode(NewRoundMode: TFPURoundingMode): TFPURoundingMode;
  {-Set new rounding mode and return the old mode}
begin
  SetRoundMode := TFPURoundingMode(math.SetRoundMode(math.TFPURoundingMode(NewRoundMode)));
end;


{---------------------------------------------------------------------------}
function  GetPrecisionMode: TFPUPrecisionMode;
  {-Return the current precision control mode}
begin
  GetPrecisionMode := TFPUPrecisionMode(math.GetPrecisionMode);
end;


{---------------------------------------------------------------------------}
function  SetPrecisionMode(NewPrecision: TFPUPrecisionMode): TFPUPrecisionMode;
  {-Set new precision control mode and return the old precision}
begin
  SetPrecisionMode := TFPUPrecisionMode(math.SetPrecisionMode(math.TFPUPrecisionMode(NewPrecision)));
end;


{---------------------------------------------------------------------------}
procedure GetExceptionMask(var Mask: byte);
  {-Return the current exception mask}
var
  MaskSet: TFPUExceptionMask;
begin
  MaskSet := Math.GetExceptionMask;
  Mask := 0;
  if exInvalidOp    in MaskSet then Mask := Mask or bexInvalidOp;
  if exDenormalized in MaskSet then Mask := Mask or bexDenormalized;
  if exZeroDivide   in MaskSet then Mask := Mask or bexZeroDivide;
  if exOverflow     in MaskSet then Mask := Mask or bexOverflow;
  if exUnderflow    in MaskSet then Mask := Mask or bexUnderflow;
  if exPrecision    in MaskSet then Mask := Mask or bexPrecision;
end;

procedure SetExceptionMask(NewMask: byte);
  {-Set new exception mask}
var
  MaskSet: TFPUExceptionMask;
begin
  MaskSet := [];
  if NewMask and bexInvalidOp    <> 0 then MaskSet := MaskSet + [exInvalidOp];
  if NewMask and bexDenormalized <> 0 then MaskSet := MaskSet + [exDenormalized];
  if NewMask and bexZeroDivide   <> 0 then MaskSet := MaskSet + [exZeroDivide];
  if NewMask and bexOverflow     <> 0 then MaskSet := MaskSet + [exOverflow];
  if NewMask and bexUnderflow    <> 0 then MaskSet := MaskSet + [exUnderflow];
  if NewMask and bexPrecision    <> 0 then MaskSet := MaskSet + [exPrecision];
  Math.SetExceptionMask(MaskSet);
end;

{$endif}

{$else}

{Do not use math}

{$ifdef BIT16}
  {$define need_get87}
  {$define need_set87}
{$else}
  {$ifdef VirtualPascal}
    {$define need_get87}
    {$define need_set87}
  {$endif}

  {$ifdef FPC}
    {$ifdef VER1}
      {$define need_get87}
      {$define need_set87}
    {$endif}
  {$endif}

  {$ifdef Delphi}
    {$ifndef CONDITIONALEXPRESSIONS}
      {$ifdef VER90 }
        {$define need_set87}
      {$endif}
      {Delphi 3+ has Set8087CW}
      {Delphi 6+ has Get8087CW}
      {$define need_get87}
    {$endif}
  {$endif}
{$endif}


{$ifdef need_get87}
  {----------------------------------------------------------------}
  function Get8087CW: word;
    {-Return the FPU control word}
  var
    cw: word;
  begin
    asm
      fstcw [cw]
      fwait
    end;
    Get8087CW := cw;
  end;
{$endif}

{$ifdef need_set87}
  {----------------------------------------------------------------}
  procedure Set8087CW(cw: word);
    {-Set new FPU control word}
  begin
    asm
      fnclex
      fldcw  [cw]
      fwait
    end;
  end;
{$endif}


{---------------------------------------------------------------------------}
function GetRoundMode: TFPURoundingMode;
  {-Return the current rounding mode}
begin
  GetRoundMode := TFPURoundingMode((Get8087CW shr 10) and 3);
end;


{---------------------------------------------------------------------------}
function SetRoundMode(NewRoundMode: TFPURoundingMode): TFPURoundingMode;
  {-Set new rounding mode and return the old mode}
var
  CW: word;
begin
  CW := Get8087CW;
  Set8087CW((CW and $F3FF) or (ord(NewRoundMode) shl 10));
  SetRoundMode := TFPURoundingMode((CW shr 10) and 3);
end;


{---------------------------------------------------------------------------}
function GetPrecisionMode: TFPUPrecisionMode;
  {-Return the current precision control mode}
begin
  GetPrecisionMode := TFPUPrecisionMode((Get8087CW shr 8) and 3);
end;


{---------------------------------------------------------------------------}
function SetPrecisionMode(NewPrecision: TFPUPrecisionMode): TFPUPrecisionMode;
  {-Set new precision control mode and return the old precision}
var
  CW: word;
begin
  CW := Get8087CW;
  Set8087CW((CW and $FCFF) or (ord(NewPrecision) shl 8));
  SetPrecisionMode := TFPUPrecisionMode((CW shr 8) and 3);
end;


{---------------------------------------------------------------------------}
procedure GetExceptionMask(var Mask: byte);
  {-Return the current exception mask}
begin
  Mask := Get8087CW and $3F;
end;


{---------------------------------------------------------------------------}
procedure SetExceptionMask(NewMask: byte);
  {-Set new exception mask}
var
  CW: word;
begin
  CW := Get8087CW;
  {$ifdef Delphi}
    {$ifndef D6Plus}
     asm
       fnclex  {don't raise pending exceptions enabled by the new mask}
     end;
    {$endif}
  {$endif}
  Set8087CW((CW and $FFC0) or (NewMask and $3F));
end;
{$endif}


{$ifdef fpu_usemath}
{$ifdef debug}
  {$ASSERTIONS ON}
{$endif}
begin
  {$ifdef UNIT_SCOPE}
    assert(sizeof(TFPURoundingMode) = sizeof(system.math.TFPURoundingMode), 'sizeof TFPURoundingMode');
    assert(sizeof(TFPUPrecisionMode)= sizeof(system.math.TFPUPrecisionMode),'sizeof TFPUPrecisionMode');
  {$else}
    assert(sizeof(TFPURoundingMode) = sizeof(math.TFPURoundingMode), 'sizeof TFPURoundingMode');
    assert(sizeof(TFPUPrecisionMode)= sizeof(math.TFPUPrecisionMode),'sizeof TFPUPrecisionMode');
  {$endif}
{$endif}

end.
