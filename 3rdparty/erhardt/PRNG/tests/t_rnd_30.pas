{Simple test for ISAAC unit, we May 2005}

program t_rnd_30;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  isaac;
begin
  writeln('Simple test for ISAAC unit     (c) 2005 W.Ehrhardt');
  writeln('ISAAC selftest: ',isaac_selftest);
end.
