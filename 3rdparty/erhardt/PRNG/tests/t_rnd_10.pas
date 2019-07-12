{Simple test for taus88 unit, we Apr.2005}

program t_rnd_10;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  Taus88;

begin
  writeln('Simple test for taus88 unit     (c) 2005 W.Ehrhardt');
  writeln('Taus88 selftest: ', taus88_selftest);
end.
