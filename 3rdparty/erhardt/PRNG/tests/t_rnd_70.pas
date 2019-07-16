{Simple test for aesr unit, we Aug.2005}

program t_rnd_70;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  aesr;

begin
  writeln('Simple test for aesr unit     (c) 2005 W.Ehrhardt');
  writeln('aesr self test: ', aesr_selftest);
end.
