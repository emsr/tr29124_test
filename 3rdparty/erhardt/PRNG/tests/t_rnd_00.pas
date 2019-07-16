{Simple test for pasrand unit, we 12.2017}

program t_rnd_00;

{$i STD.INC}

{$ifdef APPCONS}
  {$apptype console}
{$endif}

uses
  {$ifdef WINCRT}
     wincrt,
  {$endif}
  pasrand;

begin
  writeln('Simple test for pasrand unit     (c) 2071 W.Ehrhardt');
  writeln('pasrand selftest: ', pasrand_selftest);
end.
