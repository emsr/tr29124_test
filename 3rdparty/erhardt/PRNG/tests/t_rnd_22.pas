{Simple test for tt800 unit, we May 2005}

program t_rnd_22;

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
  tt800;

const
  NMAX = 100000000;
var
  n: longint;
  rng: tt800_ctx;
begin
  writeln('Test for tt800 unit: ', NMAX, ' tt800_next calls     (c) 2005 W.Ehrhardt');
  tt800_init(rng,1);
  for n:=1 to NMAX do begin
    tt800_next(rng);
  end;
  writeln('Done');
end.
