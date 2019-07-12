{Simple test for xor4096 unit, we Apr.2007}

program t_rnd_90;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  xor4096;

begin
  writeln('Simple test for xor4096 unit     (c) 2007 W.Ehrhardt');
  writeln('PRNG selftest: ', xor4096_selftest);
end.
