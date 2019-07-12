{Simple test for pasrand unit, we 12.2017}

program t_rnd_02;

{$i STD.INC}

{$ifdef BIT16}
  {$N+}
{$endif}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  pasrand;

const
  NMAX = 100000000;
var
  n: longint;
  rng: pasrand_ctx;
begin
  writeln('Test for pasrand unit: ', NMAX, ' pasrand_next calls     (c) 2071 W.Ehrhardt');
  pasrand_init(rng,1);
  for n:=1 to NMAX do begin
    pasrand_next(rng);
  end;
  writeln('Done');
end.
