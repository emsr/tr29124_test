unit pfdu_crt;

{Use of pfdu with win/crt and keypress check,  (c) WE Sep.2005}
{win/crt is linked via uses, progress is assigned in unit initialization}

interface

{$i STD.INC}

uses
  {$ifdef WINCRT} WinCRT, {$else} CRT, {$endif}
  mp_pfu, pfdu;


implementation


{---------------------------------------------------------------------------}
procedure progress(checkonly: boolean; cnt,maxcnt: longint; var cancel: boolean); {$ifdef BIT16} far; {$endif}
  {-Standard progress function}
begin
  if (cnt > 0) and not checkonly then write('.');
  cancel := abort or (keypressed and (readkey=#27));
  if cancel then begin
    abort := true;
  end;
end;


begin
  {$ifdef FPC}
    mp_set_progress(@progress);
  {$else}
    mp_set_progress(progress);
  {$endif}
end.
