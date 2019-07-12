{Simple test for tt800 unit, we May 2005}

program t_rnd_20;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  tt800;

begin
  writeln('Simple test for tt800 unit     (c) 2005 W.Ehrhardt');
  writeln('TT800 selftest: ',tt800_selftest);
end.
