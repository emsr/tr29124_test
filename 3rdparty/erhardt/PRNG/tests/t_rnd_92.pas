{Simple test for xor4096 unit, we Apr.2007}

program t_rnd_92;

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
  xor4096;

const
  NMAX = 100000000;
var
  n: longint;
  rng: xor4096_ctx;
begin
  writeln('Test for xor4096 unit: ', NMAX, ' xor4096_next calls     (c) 2007 W.Ehrhardt');
  xor4096_init(rng,1);
  for n:=1 to NMAX do begin
    xor4096_next(rng);
  end;
  writeln('Done');
end.
