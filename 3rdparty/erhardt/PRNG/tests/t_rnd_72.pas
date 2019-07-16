{Simple test for aesr unit, we Aug.2005}

program t_rnd_72;

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
  aesr;

const
  NMAX = 100000000;
var
  n: longint;
  rng: aesr_ctx;
begin
  writeln('Test for aesr unit: ', NMAX, ' aesr_next calls     (c) 2005 W.Ehrhardt');
  aesr_init(rng,1);
  for n:=1 to NMAX do begin
    aesr_next(rng);
  end;
  writeln('Done');
end.
