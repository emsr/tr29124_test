{Simple test for xor4096 unit, we Dec.2017}

program t_rnd_17;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

{$ifdef BIT16}
  {$N+}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  xor4096;
var
  rng: xor4096_ctx;
  i: integer;
const
  NR = 14*8;
begin
  writeln('Test for xor4096 unit: functions rangel/w     (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', xor4096_selftest);
  xor4096_init(rng,1);
  for i:=1 to NR do begin
    write(xor4096_rangew(rng,100):5);
    if (i mod 14 = 0) or (i=NR) then writeln;
  end;
  for i:=1 to NR do begin
    write(xor4096_rangel(rng,1234567):10);
    if (i mod 7 = 0) or (i=NR) then writeln;
  end;
end.
