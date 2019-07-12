{Simple test for well1024a unit, we 12.2017}

program t_rnd_a2;

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
  well1024;

const
  NMAX = 100000000;
var
  n: longint;
  rng: well1024a_ctx;
begin
  writeln('Test for well1024a unit: ', NMAX, ' well1024a_next calls     (c) 2017 W.Ehrhardt');
  well1024a_init(rng,1);
  for n:=1 to NMAX do begin
    well1024a_next(rng);
  end;
  writeln('Done');
end.
