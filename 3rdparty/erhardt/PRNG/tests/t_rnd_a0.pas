{Simple test for well1024 unit, we 12.2017}

program t_rnd_a0;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  well1024;

begin
  writeln('Simple test for well1024 unit     (c) 2017 W.Ehrhardt');
  writeln('PRNG selftest: ', well1024a_selftest);
end.
