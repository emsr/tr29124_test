{Simple test for kiss123 unit, we Feb.2007}

program t_rnd_62;

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
  kiss123;

const
  NMAX = 100000000;
var
  n: longint;
  rng: kiss123_ctx;
begin
  writeln('Test for kiss123 unit: ', NMAX, ' kiss123_next calls     (c) 2007 W.Ehrhardt');
  kiss123_init(rng,1);
  for n:=1 to NMAX do begin
    kiss123_next(rng);
  end;
  writeln('Done');
end.
