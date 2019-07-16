{Simple test for well1024a unit, we Dec.2017}

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
  well1024;
var
  rng: well1024a_ctx;
  i: integer;
const
  NR = 14*8;
begin
  writeln('Test for well1024 unit: functions rangel/w     (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', well1024a_selftest);
  well1024a_init(rng,1);
  for i:=1 to NR do begin
    write(well1024a_rangew(rng,100):5);
    if (i mod 14 = 0) or (i=NR) then writeln;
  end;
  for i:=1 to NR do begin
    write(well1024a_rangel(rng,1234567):10);
    if (i mod 7 = 0) or (i=NR) then writeln;
  end;
end.
