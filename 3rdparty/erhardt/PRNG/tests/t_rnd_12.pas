{Simple test for taus88 unit, we Apr.2005}

program t_rnd_12;

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
  taus88;

const
  NMAX = 100000000;
var
  n: longint;
  rng: taus88_ctx;
begin
  writeln('Test for taus88 unit: ', NMAX, ' taus88_next calls     (c) 2005 W.Ehrhardt');
  taus88_init(rng,1);
  for n:=1 to NMAX do begin
    taus88_next(rng);
  end;
  writeln('Done');
end.
