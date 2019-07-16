{Simple test for salsar unit, we Dec.2017}

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
  salsar;
var
  rng: salsar_ctx;
  i: integer;
const
  NR = 14*8;
begin
  writeln('Test for salsar unit: functions rangel/w     (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', salsar_selftest);
  salsar_init(rng,1);
  for i:=1 to NR do begin
    write(salsar_rangew(rng,100):5);
    if (i mod 14 = 0) or (i=NR) then writeln;
  end;
  for i:=1 to NR do begin
    write(salsar_rangel(rng,1234567):10);
    if (i mod 7 = 0) or (i=NR) then writeln;
  end;
end.
