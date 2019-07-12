{Simple test for pasrand unit, we 12.2017}

program t_rnd_07;

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
  pasrand;
var
  rng: pasrand_ctx;
  i: integer;
const
  NR = 14*8;
begin
  writeln('Test for pasrand unit: functions rangel/w     (c) 2017 W.Ehrhardt');
  writeln('Selftest: ', pasrand_selftest);
  pasrand_init(rng,1);
  for i:=1 to NR do begin
    write(pasrand_rangew(rng,100):5);
    if (i mod 14 = 0) or (i=NR) then writeln;
  end;
  for i:=1 to NR do begin
    write(pasrand_rangel(rng,1234567):10);
    if (i mod 7 = 0) or (i=NR) then writeln;
  end;
end.
