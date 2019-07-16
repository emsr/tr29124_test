{Simple test for kiss123 unit, we Feb.2007}

program t_rnd_60;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  kiss123;

begin
  writeln('Simple test for kiss123 unit     (c) 2007 W.Ehrhardt');
  writeln('PRNG selftest: ', kiss123_selftest);
end.
